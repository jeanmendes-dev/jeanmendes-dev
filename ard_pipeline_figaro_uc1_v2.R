# ============================================================
# CLINICAL ADaM → ML-READY ARD PIPELINE
# STUDY: SHP647UC301-FIGARO-UC1
# ============================================================

# ===============================
# LIBRARIES
# ===============================
library(data.table)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

# ===============================
# CONFIG
# ===============================
DATA_PATH    <- "/domino/datasets/local/clinical-trial-data/SHP647UC301-FIGARO-UC1/data"

OUTPUT_ARD      <- "/mnt/figaro_uc1_ml_dataset_v2.csv"
OUTPUT_COVERAGE <- "/mnt/figaro_uc1_ml_coverage_report_v2.csv"

# ---------------------------------------------------------------
# Datasets confirmados no diagnóstico FIGARO UC1
#
# Particularidades vs estudos UNIFI / GALAXI:
#   1. Sem sufixo .sas7bdat nos nomes de arquivo
#   2. 9 datasets — estudo de fase inicial, n=12 sujeitos
#   3. Sem PARAMCD — identificadores por domínio ADaM
#      (LBTESTCD, VSTESTCD, EGTESTCD, MITESTCD, MOTESTCD)
#   4. Sem AVISITN nativo — usa VISITNUM + VISIT
#   5. Flag de análise principal: SAFFL
#      Flags de população: ITTFL, SAFFL, PPROTFL, ENRLFL (ADSL)
#   6. Identificação de braço: ARM / ARMCD (sem TRT01P/TRT02P)
# ---------------------------------------------------------------
EXPECTED_DATASETS <- c(
  "adae",   # Adverse Events          — 28 linhas,  sem VISITNUM, sem PARAMCD
  "adcm",   # Concomitant Medications — 211 linhas, sem VISITNUM, sem PARAMCD
  "adeg",   # ECG                     — 336 linhas, tem VISITNUM
  "adex",   # Exposure                — 33 linhas,  tem VISITNUM
  "adlb",   # Laboratory              — 3796 linhas, tem VISITNUM  ← âncora spine
  "admi",   # Microbiology            — 850 linhas, tem VISITNUM
  "admo",   # Morphology/Histology    — 1829 linhas, tem VISITNUM
  "adsl",   # Subject-Level           — 12 linhas,  sem VISITNUM  ← estático puro
  "advs"    # Vital Signs             — 321 linhas, tem VISITNUM
)

# ---------------------------------------------------------------
# NAMESPACE DE AVISITN
#
#   1 – 99   → VISITNUM real do protocolo
#              10 = SCREENING
#              20 = BASELINE/WEEK 0
#              40 = WEEK 4
#              50 = WEEK 8
#              60 = WEEK 12
#              71 = WEEK 16/ET PART 1
#              72 = WEEK 16/ET PART 2
#              73 = WEEK 16/ET PART 3
#              90 = WEEK 32
#
#   101+     → eventos de adae/adcm sem visita protocolar
#              PSEUDO_1 = AVISITN 101, PSEUDO_2 = 102 ...
#              Offset = 100 (maior VISITNUM do protocolo = 90)
#
# ADSL não entra na spine — é estático puro.
# As colunas STC_ADSL_* chegam ao ARD via join por USUBJID.
# ---------------------------------------------------------------
PSEUDO_OFFSET <- 100L

# ===============================
# LOGGING
# ===============================
log_warn <- function(msg) warning(paste0("[WARNING] ", msg))
log_info <- function(msg) message(paste0("[INFO] ", msg))

# ===============================
# 1. LOAD DATASETS
# ===============================
load_datasets <- function(path) {
  files <- list.files(path = path, pattern = "\\.csv$", full.names = TRUE)

  datasets <- lapply(files, function(f) {
    tryCatch({
      as.data.frame(fread(f, sep = "auto"))
    }, error = function(e) {
      log_warn(paste("Falha ao ler:", f))
      return(NULL)
    })
  })

  names(datasets) <- tolower(tools::file_path_sans_ext(basename(files)))
  datasets <- datasets[!sapply(datasets, is.null)]

  missing_ds <- setdiff(EXPECTED_DATASETS, names(datasets))
  if (length(missing_ds) > 0) {
    log_warn(paste("Datasets esperados não encontrados:",
                   paste(missing_ds, collapse = ", ")))
  }

  log_info(paste("Carregados", length(datasets), "datasets"))
  return(datasets)
}

# ===============================
# 2. NORMALISE VISIT KEY
# ===============================
# Cria AVISITN a partir de VISITNUM (datasets protocolares) ou
# sequência com offset (adae/adcm sem VISITNUM), garantindo
# namespaces separados que nunca colidem.
normalise_visit_key <- function(df, dataset_name) {

  if (!("USUBJID" %in% names(df))) return(df)

  # Caso 1: AVISITN já existe — não esperado no FIGARO UC1, mas defensivo
  if ("AVISITN" %in% names(df)) return(df)

  # Caso 2: VISITNUM disponível — promover para AVISITN (range 1–99)
  if ("VISITNUM" %in% names(df)) {
    df <- df %>%
      rename(AVISITN = VISITNUM) %>%
      mutate(AVISIT = if ("VISIT" %in% names(.)) VISIT
                      else paste0("VISIT_", AVISITN))
    return(df)
  }

  # Caso 3: sem VISITNUM (adae, adcm) — offset 100+, ordem por data
  log_warn(paste("Criando pseudo AVISITN (offset", PSEUDO_OFFSET, ") para:", dataset_name))

  date_col <- intersect(c("ADT", "ASTDT", "ASTDTC"), names(df))[1]

  if (!is.na(date_col)) {
    df <- df %>%
      arrange(USUBJID, .data[[date_col]]) %>%
      group_by(USUBJID) %>%
      mutate(AVISITN = PSEUDO_OFFSET + row_number(),
             AVISIT  = paste0("PSEUDO_", row_number())) %>%
      ungroup()
  } else {
    df <- df %>%
      arrange(USUBJID) %>%
      group_by(USUBJID) %>%
      mutate(AVISITN = PSEUDO_OFFSET + row_number(),
             AVISIT  = paste0("PSEUDO_", row_number())) %>%
      ungroup()
  }

  return(df)
}

# ===============================
# 3. BUILD SPINE
# ===============================
# Âncora primária: adlb (3796 linhas — maior dataset com visitas).
# admo (1829) e admi (850) complementam histologia e microbiologia.
# adeg, adex, advs contribuem visitas adicionais.
# adae e adcm contribuem os AVISITNs pseudo (101+).
# adsl excluído da spine — estático, join por USUBJID em merge_ard.
#
# Deduplicação: para cada USUBJID×AVISITN mantém o label real
# (não-PSEUDO) se existir, senão mantém o PSEUDO_*.
build_spine <- function(datasets) {

  spine_parts <- list()

  anchor <- "adlb"
  if (anchor %in% names(datasets)) {
    spine_parts[[anchor]] <- datasets[[anchor]] %>%
      select(any_of(c("USUBJID", "AVISIT", "AVISITN"))) %>%
      distinct()
    log_info(paste("Usando", anchor, "como âncora da spine"))
  } else {
    log_warn("Âncora adlb não encontrada — spine construída de todos os datasets")
  }

  # Complemento: demais datasets com AVISITN, exceto adsl
  for (nm in names(datasets)) {
    if (nm %in% c(anchor, "adsl")) next
    df <- datasets[[nm]]
    if (all(c("USUBJID", "AVISITN") %in% names(df))) {
      spine_parts[[nm]] <- df %>%
        select(any_of(c("USUBJID", "AVISIT", "AVISITN"))) %>%
        distinct()
    }
  }

  spine_raw <- bind_rows(spine_parts) %>% distinct()

  # Deduplicar: para cada USUBJID×AVISITN, preferir label real a PSEUDO_*
  spine <- spine_raw %>%
    group_by(USUBJID, AVISITN) %>%
    summarise(
      AVISIT = {
        real <- AVISIT[!grepl("^PSEUDO_", AVISIT) & !is.na(AVISIT) & AVISIT != ""]
        if (length(real) > 0) real[1] else AVISIT[!is.na(AVISIT) & AVISIT != ""][1]
      },
      .groups = "drop"
    ) %>%
    arrange(USUBJID, AVISITN)

  # Verificações de sanidade
  dec_n <- sum(spine$AVISITN != as.integer(spine$AVISITN), na.rm = TRUE)
  if (dec_n > 0) log_warn(paste("AVISITNs decimais na spine:", dec_n))

  blank_n <- sum(is.na(spine$AVISIT) | spine$AVISIT == "", na.rm = TRUE)
  if (blank_n > 0) log_warn(paste("AVISITs em branco na spine:", blank_n))

  log_info(paste("Spine criada com", nrow(spine), "linhas e",
                 n_distinct(spine$USUBJID), "sujeitos únicos"))
  return(spine)
}

# ===============================
# 4. CLASSIFICATION (vectorized)
# ===============================
classify_all_variables <- function(df, vars) {

  dt <- as.data.table(df[, c("USUBJID", vars), drop = FALSE])

  variability <- dt[, lapply(.SD, function(x) uniqueN(x, na.rm = TRUE)),
                    by = USUBJID, .SDcols = vars]

  is_longitudinal <- sapply(vars,
    function(v) any(variability[[v]] > 1, na.rm = TRUE))

  return(ifelse(is_longitudinal, "LONGITUDINAL", "STATIC"))
}

# ===============================
# 5. SMART AGGREGATION
# ===============================
aggregate_value <- function(x) {

  if (is.numeric(x))                   return(mean(x, na.rm = TRUE))
  if (is.character(x) || is.factor(x)) return(first(na.omit(x)))
  return(first(x))
}

# ===============================
# 6. EXTRACT VARIABLES
# ===============================
# GARANTIA DE GRAIN: cada objeto extraído tem exactamente uma linha
# por USUBJID (estático) ou USUBJID×AVISITN (longitudinal).
# A agregação acontece ANTES de qualquer join — elimina a fonte dos
# produtos cartesianos que geravam AVISITNs decimais (10.1, 10.2…).
extract_dataset <- function(df, dataset_name) {

  vars <- setdiff(names(df), c("USUBJID", "AVISIT", "AVISITN"))

  extracted  <- list()
  used_vars  <- c()

  # Prefixo limpo — sem sufixo .sas7bdat no FIGARO UC1
  clean_name <- toupper(dataset_name)

  class_map   <- classify_all_variables(df, vars)
  long_vars   <- vars[class_map == "LONGITUDINAL"]
  static_vars <- vars[class_map == "STATIC"]

  log_info(paste0("  ", dataset_name, ": ", length(long_vars),
                  " longitudinais / ", length(static_vars), " estáticas"))

  # --- LONGITUDINAIS: uma linha por USUBJID × AVISITN ---
  if (length(long_vars) > 0) {

    dt_long  <- as.data.table(df[, c("USUBJID", "AVISITN", long_vars), drop = FALSE])
    agg_long <- dt_long[, lapply(.SD, aggregate_value),
                        by = .(USUBJID, AVISITN), .SDcols = long_vars]

    n_dup <- nrow(agg_long) - nrow(unique(agg_long[, .(USUBJID, AVISITN)]))
    if (n_dup > 0) {
      log_warn(paste0("  ", dataset_name,
                      ": ", n_dup, " duplicatas USUBJID×AVISITN após agregação"))
    }

    for (v in long_vars) {
      colname <- paste0("AVAL_", clean_name, "_", v)
      tmp     <- agg_long[, .(USUBJID, AVISITN, value = get(v))]
      setnames(tmp, "value", colname)
      extracted[[colname]] <- as.data.frame(tmp)
      used_vars <- c(used_vars, v)
    }
  }

  # --- ESTÁTICAS: uma linha por USUBJID ---
  if (length(static_vars) > 0) {

    dt_static  <- as.data.table(df[, c("USUBJID", static_vars), drop = FALSE])
    agg_static <- dt_static[, lapply(.SD, aggregate_value),
                             by = USUBJID, .SDcols = static_vars]

    for (v in static_vars) {
      colname <- paste0("STC_", clean_name, "_", v)
      tmp     <- agg_static[, .(USUBJID, value = get(v))]
      setnames(tmp, "value", colname)
      extracted[[colname]] <- as.data.frame(tmp)
      used_vars <- c(used_vars, v)
    }
  }

  return(list(
    extracted = extracted,
    used_vars = used_vars,
    all_vars  = vars
  ))
}

# ===============================
# 7. MERGE ARD
# ===============================
merge_ard <- function(spine, extracted_all) {

  ard <- spine

  for (ds in extracted_all) {
    for (obj in ds$extracted) {

      if (all(c("USUBJID", "AVISITN") %in% names(obj))) {

        # Guardar rail de grain: re-agregar antes do join
        obj_dt <- as.data.table(obj)
        obj_dt <- obj_dt[, lapply(.SD, aggregate_value),
                         by = .(USUBJID, AVISITN)]
        obj    <- as.data.frame(obj_dt)

        ard <- left_join(ard, obj, by = c("USUBJID", "AVISITN"))

      } else {
        ard <- left_join(ard, obj, by = "USUBJID")
      }
    }
  }

  # Verificar grain final
  dup_check <- ard %>% count(USUBJID, AVISITN) %>% filter(n > 1)
  if (nrow(dup_check) > 0) {
    log_warn(paste("Duplicatas USUBJID×AVISITN no ARD final:",
                   nrow(dup_check), "combinações"))
  } else {
    log_info("Grain OK — sem duplicatas USUBJID×AVISITN")
  }

  # Verificar AVISITNs decimais residuais
  dec_n <- sum(ard$AVISITN != as.integer(ard$AVISITN), na.rm = TRUE)
  if (dec_n > 0) {
    log_warn(paste("AVISITNs decimais no ARD final:", dec_n, "linhas"))
  } else {
    log_info("AVISITN OK — sem valores decimais")
  }

  # Verificar AVISITs em branco
  blank_n <- sum(is.na(ard$AVISIT) | ard$AVISIT == "", na.rm = TRUE)
  if (blank_n > 0) {
    log_warn(paste("AVISITs em branco no ARD final:", blank_n, "linhas"))
  } else {
    log_info("AVISIT OK — sem células em branco")
  }

  # Ordenar colunas: USUBJID | AVISITN | AVISIT | restante
  other_cols <- setdiff(names(ard), c("USUBJID", "AVISITN", "AVISIT"))
  ard <- ard %>% select(USUBJID, AVISITN, AVISIT, all_of(other_cols))

  log_info(paste("ARD final:", nrow(ard), "linhas ×", ncol(ard), "colunas"))
  return(ard)
}

# ===============================
# 8. COVERAGE REPORT
# ===============================
build_coverage <- function(extracted_all, dataset_names) {

  report <- map2_df(extracted_all, dataset_names, function(ds, name) {

    total     <- length(ds$all_vars)
    extracted <- length(ds$used_vars)
    missing   <- setdiff(ds$all_vars, ds$used_vars)

    tibble(
      dataset           = name,
      total_vars        = total,
      extracted_vars    = extracted,
      missing_vars      = total - extracted,
      extraction_rate   = round(100 * extracted / total, 2),
      missing_var_names = paste(missing, collapse = ", ")
    )
  })

  print(report)
  write.csv(report, OUTPUT_COVERAGE, row.names = FALSE)
  log_info(paste("Coverage report salvo em:", OUTPUT_COVERAGE))

  return(report)
}

# ===============================
# 9. MAIN PIPELINE
# ===============================
run_pipeline <- function() {

  log_info("=== FIGARO UC1 ARD PIPELINE v2 INICIADO ===")

  # 1. Carregar datasets
  datasets <- load_datasets(DATA_PATH)

  # 2. Normalizar chave de visita (VISITNUM → AVISITN; offset 100 para pseudos)
  datasets <- imap(datasets, normalise_visit_key)

  # 3. Spine baseada em adlb + complemento (adsl excluído)
  spine <- build_spine(datasets)

  # 4. Extrair variáveis de cada dataset
  extracted_all <- list()

  for (name in names(datasets)) {
    log_info(paste("Processando dataset:", name))
    tryCatch({
      res <- extract_dataset(datasets[[name]], name)
      extracted_all[[name]] <- res
    }, error = function(e) {
      log_warn(paste("Erro ao processar", name, ":", conditionMessage(e)))
    })
  }

  # 5. Montar ARD (grain garantido + ordem de colunas correcta)
  ard <- merge_ard(spine, extracted_all)

  # 6. Salvar ARD
  write.csv(ard, OUTPUT_ARD, row.names = FALSE)
  log_info(paste("ARD salvo em:", OUTPUT_ARD))

  # 7. Coverage report
  build_coverage(extracted_all, names(datasets))

  log_info("=== PIPELINE CONCLUÍDO COM SUCESSO ===")

  invisible(ard)
}

# ===============================
# RUN
# ===============================
run_pipeline()
