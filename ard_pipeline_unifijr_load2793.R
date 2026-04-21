# ============================================================
# CLINICAL ADaM → ML-READY ARD PIPELINE
# STUDY: CNTO1275PUC3001-UNIFI-JR / load-2793
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
DATA_PATH    <- "/domino/datasets/local/clinical-trial-data/CNTO1275PUC3001-UNIFI-JR/load-2793/Data/_csv"

OUTPUT_ARD      <- "/mnt/unifijr_load2793_ml_dataset_v1.csv"
OUTPUT_COVERAGE <- "/mnt/unifijr_load2793_ml_coverage_report_v1.csv"

# ---------------------------------------------------------------
# Datasets confirmados no diagnóstico UNIFI JR (load-2793)
#
# UNIFI JR é o estudo PEDIÁTRICO de UC (CNTO1275PUC3001).
# Apenas 7 datasets — muito mais enxuto que o UNIFI adulto
# (33 datasets no load-1899).
# ---------------------------------------------------------------
EXPECTED_DATASETS <- c(
  "adbdc",     # Baseline disease characteristics — 1.531 linhas, sem AVISITN
  "adeff",     # Efficacy endpoints               — 6.251 linhas, tem AVISITN
  "adhist",    # Histology                        — 11.945 linhas, tem AVISITN  ← âncora spine
  "adlbef",    # Lab efficacy (CRP, calprotectina)— 6.487 linhas, tem AVISITN
  "admayo",    # Mayo scores                      — 58.353 linhas, tem AVISITN
  "adpucai",   # PUCAI score (pediátrico)         — 9.044 linhas, tem AVISITN
  "adsl",      # Subject-level                    — 168 linhas,   sem AVISITN
  "advpucai"   # PUCAI virtual/diary              — 11.571 linhas, tem AVISITN
)

# ---------------------------------------------------------------
# Datasets SEM AVISITN no diagnóstico UNIFI JR
#   adbdc — sem AVISITN; tem ASTDT (baseline cross-sectional)
#   adsl  — sem AVISITN; 1 linha por sujeito (168 sujeitos)
#           → estático puro, não entra na spine
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Particularidades do UNIFI JR vs UNIFI adulto (load-1899/1901)
#
# 1. ESTUDO PEDIÁTRICO — n=168 sujeitos vs 1.331 no adulto
#    Indicação: Colite Ulcerosa Pediátrica (pUC)
#
# 2. APENAS 7 DATASETS
#    Sem: adae, adchem, adcm, adcort, addeath, adds, addv,
#         adecon, adeq5d, adex, adhema, adho, adibdq, adis,
#         admh, adpc, adsaf, adsf36, adsg, adtfi, adtfm,
#         adtte, aduceis, advs, advscort
#
# 3. DATASETS EXCLUSIVOS DO UNIFI JR (pediátrico)
#    adpucai   — PUCAI score (Pediatric Ulcerative Colitis
#                Activity Index): ABPAIN, RBLEED, STLCONS,
#                NUMSTL, NOCTSTL, ACTLEV, PUCAITS
#    advpucai  — PUCAI virtual/diary com duplo período
#                (sufixos 1/2 nos PARAMCDs)
#    adhist    — histologia com Geboes/RHI scores:
#                ARCHCHG, CHINFIN, EOS, NEUTLP, EROULC,
#                ULCER, RHI, GBTOT, RHITOT, NITOT...
#
# 4. FLAGS DE ANÁLISE — UNIFI JR
#    Mantém: RANDFL, SAFFL
#    Novos vs UNIFI adulto:
#      SAFRFL, SAFCRFL      — randomised/crossover safety
#      FASFL, FASRFL, FASCRFL, FASRESFL — Full Analysis Set
#                                          variants
#    Ausentes vs adulto: SAFW8FL, SAF2FL, COMPLFL, ENTERMFL,
#                        ENEXTFL, RAND2FL
#
# 5. ESTRUTURA DE TRATAMENTO — 3 PERÍODOS
#    TRT01P/TRT01A, TRT02P/TRT02A (com TR02PG1/TR02AG1),
#    TRT03P/TRT03A — mesma estrutura do UNIFI adulto
#
# 6. ÂNCORA DA SPINE — ADMAYO
#    admayo: 58.353 linhas — maior dataset, cobre endpoints
#    centrais UC: MAYO, MMAYO, PMAYO, SFSCORE, RBSCORE,
#    ENSCORE, PGSCORE — equivalente ao adulto.
#    Complementada com adsl para cobertura dos 168 sujeitos.
#
# 7. SUFIXO DE FICHEIRO — SEM .sas7bdat
#    (igual ao GALAXI-1/2 e FIGARO UC1;
#     diferente do UNIFI adulto load-1899/1901 e GALAXI-3)
# ---------------------------------------------------------------

ANALYSIS_FLAGS_UNIFIJR <- c(
  "RANDFL", "SAFFL", "SAFRFL", "SAFCRFL",
  "FASFL", "FASRFL", "FASCRFL", "FASRESFL"
)

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

  # Avisar sobre datasets esperados mas não encontrados
  missing_ds <- setdiff(EXPECTED_DATASETS, names(datasets))
  if (length(missing_ds) > 0) {
    log_warn(paste("Datasets esperados não encontrados:",
                   paste(missing_ds, collapse = ", ")))
  }

  log_info(paste("Carregados", length(datasets), "datasets"))
  return(datasets)
}

# ===============================
# 2. ADD PSEUDO VISIT
# ===============================
# Apenas adbdc não tem AVISITN (adsl é excluído da spine).
# Para adbdc, cria sequência temporal por ASTDT.
add_pseudo_visit <- function(df, dataset_name) {

  if (!("USUBJID" %in% names(df))) return(df)

  if (!("AVISITN" %in% names(df))) {

    log_warn(paste("Criando pseudo AVISITN para:", dataset_name))

    date_col <- intersect(
      c("ADT", "ASTDT", "ASTDTC", "DVSTDT", "DTHDTC"),
      names(df)
    )[1]

    if (!is.na(date_col)) {
      df <- df %>%
        arrange(USUBJID, .data[[date_col]]) %>%
        group_by(USUBJID) %>%
        mutate(
          AVISITN = row_number(),
          AVISIT  = paste0("PSEUDO_", AVISITN)
        ) %>%
        ungroup()
    } else {
      df <- df %>%
        arrange(USUBJID) %>%
        group_by(USUBJID) %>%
        mutate(
          AVISITN = row_number(),
          AVISIT  = paste0("PSEUDO_", AVISITN)
        ) %>%
        ungroup()
    }
  }

  return(df)
}

# ===============================
# 3. BUILD SPINE
# ===============================
# Âncora primária: admayo (58.353 linhas)
# Endpoints centrais pUC: MAYO, MMAYO, PMAYO, SFSCORE,
# RBSCORE, ENSCORE, PGSCORE — inclui variantes LV/LT/DIÁRIO.
# Complementada com adsl para cobertura total dos 168 sujeitos.
# adsl excluído da spine — estático, join por USUBJID em merge_ard.
build_spine <- function(datasets) {

  spine_parts <- list()

  anchor <- "admayo"
  if (anchor %in% names(datasets)) {
    spine_parts[[anchor]] <- datasets[[anchor]] %>%
      select(any_of(c("USUBJID", "AVISIT", "AVISITN"))) %>%
      distinct()
    log_info(paste("Usando", anchor, "como âncora da spine"))
  } else {
    log_warn("Âncora admayo não encontrada — spine construída de todos os datasets")
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

  spine <- bind_rows(spine_parts) %>% distinct()

  # Garantir cobertura total de sujeitos via adsl (168 sujeitos)
  if ("adsl" %in% names(datasets)) {
    subj_adsl <- datasets[["adsl"]] %>%
      select(USUBJID) %>%
      distinct() %>%
      mutate(AVISITN = 0L, AVISIT = "BASELINE")

    spine <- bind_rows(spine, subj_adsl) %>% distinct()
  }

  log_info(paste("Spine criada com", nrow(spine), "linhas e",
                 n_distinct(spine$USUBJID), "sujeitos únicos"))
  return(spine)
}

# ===============================
# 4. CLASSIFICATION (vectorized)
# ===============================
# Classifica TODAS as variáveis de um dataset em uma única passagem.
# admayo (58k linhas) e adhist (12k) beneficiam da abordagem batch.
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
# GARANTIA DE GRAIN: cada objeto extraído tem exactamente
# uma linha por USUBJID (estático) ou USUBJID×AVISITN (longitudinal).
# Prefixos de coluna sem sufixo .sas7bdat (UNIFI JR não usa sufixo).
extract_dataset <- function(df, dataset_name) {

  vars <- setdiff(names(df), c("USUBJID", "AVISIT", "AVISITN"))

  extracted  <- list()
  used_vars  <- c()

  # Sem sufixo .sas7bdat no UNIFI JR
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

        # Guardar rail de grain antes do join
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

  dup_check <- ard %>%
    count(USUBJID, AVISITN) %>%
    filter(n > 1)

  if (nrow(dup_check) > 0) {
    log_warn(paste("Duplicatas USUBJID × AVISITN detectadas:",
                   nrow(dup_check), "combinações"))
  } else {
    log_info("Grain OK — sem duplicatas USUBJID × AVISITN")
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

  log_info("=== UNIFI JR load-2793 ARD PIPELINE INICIADO ===")

  # 1. Carregar
  datasets <- load_datasets(DATA_PATH)

  # 2. Pseudo-visita para datasets sem AVISITN (adbdc)
  datasets <- imap(datasets, add_pseudo_visit)

  # 3. Spine baseada em admayo + adsl
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
