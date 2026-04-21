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

OUTPUT_ARD      <- "/mnt/figaro_uc1_ml_dataset_v1.csv"
OUTPUT_COVERAGE <- "/mnt/figaro_uc1_ml_coverage_report_v1.csv"

# ---------------------------------------------------------------
# Datasets confirmados no diagnóstico FIGARO UC1
#
# ATENÇÃO — convenções distintas de todos os estudos anteriores:
#
#   1. SEM sufixo .sas7bdat nos nomes de arquivo
#      (UNIFI load-1899/1901, GALAXI-3 tinham .sas7bdat;
#       GALAXI-1/2 e FIGARO UC1 não têm)
#
#   2. Apenas 9 datasets (vs 27-41 nos estudos anteriores)
#      Estudo de fase inicial / piloto com n=12 sujeitos
#
#   3. SEM PARAMCD em vários datasets — usam colunas de domínio
#      ADaM padrão (LBTESTCD, VSTESTCD, EGTESTCD, etc.)
#      em vez de PARAMCD/PARAM
#
#   4. Nenhum dataset exclusivo de UC (sem ADMAYO, ADUCEIS,
#      ADIBDQ, ADSF36, ADCORT, ADLBEF, ADHEMA, ADCHEM)
#
#   5. FLAG DE ANÁLISE PRINCIPAL: SAFFL
#      Sem SAFW8FL, SAF2FL, RANDFL, RAND2FL, ENTERMFL, etc.
#      Flags de população: ITTFL, SAFFL, PPROTFL, ENRLFL (via ADSL)
#
#   6. IDENTIFICAÇÃO DE BRAÇO: ARM / ARMCD (sem TRT01P/TRT02P)
# ---------------------------------------------------------------
EXPECTED_DATASETS <- c(
  "adae",   # Adverse Events          — 28 linhas,  sem AVISITN, sem PARAMCD
  "adcm",   # Concomitant Medications — 211 linhas, sem AVISITN, sem PARAMCD
  "adeg",   # ECG                     — 336 linhas, tem ADT/VISITNUM, sem AVISITN
  "adex",   # Exposure                — 33 linhas,  tem VISITNUM, sem AVISITN
  "adlb",   # Laboratory              — 3796 linhas, tem ADT/VISITNUM, sem AVISITN
  "admi",   # Microbiology            — 850 linhas, tem ADT/VISITNUM, sem AVISITN
  "admo",   # Morphology (histology)  — 1829 linhas, tem ADT/VISITNUM, sem AVISITN
  "adsl",   # Subject-Level           — 12 linhas,  sem AVISITN, sem PARAMCD
  "advs"    # Vital Signs             — 321 linhas, tem ADT/VISITNUM, sem AVISITN
)

# ---------------------------------------------------------------
# Datasets SEM AVISITN no diagnóstico FIGARO UC1
#
# Nenhum dataset tem coluna AVISITN — o estudo usa VISITNUM + VISIT
# como identificador de visita, sem o campo numérico AVISITN
# padrão dos estudos UNIFI/GALAXI.
#
# Estratégia: usar VISITNUM como proxy de AVISITN onde disponível;
# fallback para sequência temporal por ADT quando VISITNUM ausente;
# fallback final para row_number() por USUBJID.
# ---------------------------------------------------------------
DATASETS_COM_VISITNUM <- c("adeg", "adex", "adlb", "admi", "admo", "advs")
DATASETS_SEM_VISITNUM <- c("adae", "adcm", "adsl")

# ---------------------------------------------------------------
# Datasets com identificadores de teste por domínio (sem PARAMCD)
#   adeg  — EGTESTCD / EGTEST  (ECG)
#   adlb  — LBTESTCD / LBTEST  (Laboratory)
#   admi  — MITESTCD / MITEST  (Microbiology)
#   admo  — MOTESTCD / MOTEST  (Morphology)
#   advs  — VSTESTCD / VSTEST  (Vital Signs)
# ---------------------------------------------------------------

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
# FIGARO UC1 não tem AVISITN — usa VISITNUM como chave de visita.
#
# NAMESPACE DE AVISITN:
#   1–99   → VISITNUM real do protocolo (adlb, advs, adeg, adex, admi, admo)
#            Ex: 10=SCREENING, 20=BASELINE/WEEK 0, 40=WEEK 4, 50=WEEK 8…
#            O maior VISITNUM observado é 90 (WEEK 32).
#
#   101+   → eventos sem visita protocolar (adae, adcm)
#            PSEUDO_1 → AVISITN=101, PSEUDO_2 → AVISITN=102…
#            Offset = MAX_VISITNUM + 1 = 100 + seq
#            Garante que PSEUDO_N e AVISITN=N nunca colidam.
#
# ADSL não entra na spine — é dataset estático puro (1 linha/sujeito).
# As colunas STC_ADSL_* chegam ao ARD via join por USUBJID, sem precisar
# de AVISITN. Não há AVISITN=0 artificial.

PSEUDO_OFFSET <- 100L   # maior VISITNUM do protocolo; pseudos começam em 101

normalise_visit_key <- function(df, dataset_name) {

  if (!("USUBJID" %in% names(df))) return(df)

  # Caso 1: AVISITN já existe (não esperado no FIGARO UC1, mas defensivo)
  if ("AVISITN" %in% names(df)) return(df)

  # Caso 2: VISITNUM disponível — promover para AVISITN (range protocolar)
  if ("VISITNUM" %in% names(df)) {
    df <- df %>%
      rename(AVISITN = VISITNUM) %>%
      mutate(AVISIT = if ("VISIT" %in% names(.)) VISIT
                      else paste0("VISIT_", AVISITN))
    return(df)
  }

  # Caso 3: sem VISITNUM — sequência com offset para isolar namespace
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
# Âncora primária: adlb (3.796 linhas — maior dataset com visitas)
# Cobre laboratório longitudinal, visitas mais densas do estudo.
# admo (1.829 linhas) e admi (850) complementam visitas de histologia
# e microbiologia.
#
# ADSL não entra na spine: é dataset estático puro — suas colunas
# chegam ao ARD via join por USUBJID. Injetar AVISITN=0 artificial
# criaria uma linha extra sem dados de visita.
#
# DEDUPLICAÇÃO DE AVISIT:
# Para cada USUBJID × AVISITN, mantém o label real (não-PSEUDO) se
# existir, senão mantém o PSEUDO_*. Garante AVISITN sempre inteiro
# e AVISIT sempre preenchido.
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

  # Complemento: demais datasets com AVISITN (após normalização)
  # ADSL excluído intencionalmente — é estático, não tem visitas
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

  # Deduplicar: para cada USUBJID × AVISITN, manter label real se existir,
  # senão manter o PSEUDO_*. Elimina duplicatas residuais e garante
  # AVISITN sempre inteiro e AVISIT sempre preenchido.
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
  decimal_check <- spine %>% filter(AVISITN != as.integer(AVISITN))
  if (nrow(decimal_check) > 0) {
    log_warn(paste("AVISITNs decimais ainda presentes:", nrow(decimal_check), "linhas"))
  }

  blank_check <- spine %>% filter(is.na(AVISIT) | AVISIT == "")
  if (nrow(blank_check) > 0) {
    log_warn(paste("AVISITs em branco ainda presentes:", nrow(blank_check), "linhas"))
  }

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

  if (is.numeric(x))                 return(mean(x, na.rm = TRUE))
  if (is.character(x) || is.factor(x)) return(first(na.omit(x)))
  return(first(x))
}

# ===============================
# 6. EXTRACT VARIABLES
# ===============================
extract_dataset <- function(df, dataset_name) {

  vars <- setdiff(names(df), c("USUBJID", "AVISIT", "AVISITN"))

  extracted  <- list()
  used_vars  <- c()

  # Nome limpo para prefixo — FIGARO UC1 não tem sufixo .sas7bdat
  clean_name <- toupper(dataset_name)

  class_map <- classify_all_variables(df, vars)

  long_vars   <- vars[class_map == "LONGITUDINAL"]
  static_vars <- vars[class_map == "STATIC"]

  log_info(paste0("  ", dataset_name, ": ", length(long_vars),
                  " longitudinais / ", length(static_vars), " estáticas"))

  # --- LONGITUDINAIS: agregação por USUBJID × AVISITN ---
  if (length(long_vars) > 0) {

    dt_long <- as.data.table(df[, c("USUBJID", "AVISITN", long_vars),
                                drop = FALSE])

    agg_long <- dt_long[, lapply(.SD, aggregate_value),
                        by = .(USUBJID, AVISITN), .SDcols = long_vars]

    for (v in long_vars) {
      colname <- paste0("AVAL_", clean_name, "_", v)
      tmp     <- agg_long[, .(USUBJID, AVISITN, value = get(v))]
      setnames(tmp, "value", colname)
      extracted[[colname]] <- as.data.frame(tmp)
      used_vars <- c(used_vars, v)
    }
  }

  # --- ESTÁTICAS: agregação por USUBJID ---
  if (length(static_vars) > 0) {

    dt_static <- as.data.table(df[, c("USUBJID", static_vars), drop = FALSE])

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
        ard <- left_join(ard, obj, by = c("USUBJID", "AVISITN"))
      } else {
        ard <- left_join(ard, obj, by = "USUBJID")
      }
    }
  }

  dup_check <- ard %>%
    count(USUBJID, AVISITN) %>%
    filter(n > 1)

  if (nrow(dup_check) > 0) {
    log_warn(paste("Duplicatas USUBJID × AVISITN detectadas:",
                   nrow(dup_check), "combinações"))
  }

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
  fwrite(report, OUTPUT_COVERAGE)
  log_info(paste("Coverage report salvo em:", OUTPUT_COVERAGE))

  return(report)
}

# ===============================
# 9. MAIN PIPELINE
# ===============================
run_pipeline <- function() {

  log_info("=== FIGARO UC1 ARD PIPELINE INICIADO ===")

  # 1. Carregar
  datasets <- load_datasets(DATA_PATH)

  # 2. Normalizar chave de visita (VISITNUM → AVISITN; fallback temporal)
  datasets <- imap(datasets, normalise_visit_key)

  # 3. Spine baseada em adlb + adsl
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

  # 5. Montar ARD
  ard <- merge_ard(spine, extracted_all)

  # 6. Salvar ARD
  fwrite(ard, OUTPUT_ARD)
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
