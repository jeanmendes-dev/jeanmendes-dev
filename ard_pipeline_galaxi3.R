# ============================================================
# CLINICAL ADaM → ML-READY ARD PIPELINE
# STUDY: CNTO1959CRD3001-GALAXI-GAL3-WK48
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
DATA_PATH    <- "/domino/datasets/local/clinical-trial-data/CNTO1959CRD3001-GALAXI-GAL3-WK48/load-2491/Data/csv"

OUTPUT_ARD      <- "/mnt/galaxi3_ml_dataset_v1.csv"
OUTPUT_COVERAGE <- "/mnt/galaxi3_ml_coverage_report_v1.csv"

# ---------------------------------------------------------------
# Datasets confirmados no diagnóstico GALAXI-3
# Nota: os arquivos têm sufixo ".sas7bdat" no nome — os CSVs
# gerados mantêm esse sufixo; o pipeline normaliza para lowercase.
# ---------------------------------------------------------------
EXPECTED_DATASETS <- c(
  "adae.sas7bdat",
  "adbdc.sas7bdat",
  "adcdai.sas7bdat",
  "addisp.sas7bdat",
  "adex.sas7bdat",
  "adice.sas7bdat",
  "adlbef.sas7bdat",
  "adrxfail.sas7bdat",
  "adrxhist.sas7bdat",
  "adsaf.sas7bdat",
  "adsl.sas7bdat",
  "advscdai.sas7bdat"
)

# ---------------------------------------------------------------
# Datasets SEM AVISITN no diagnóstico (recebem pseudo-visita)
#   adae.sas7bdat    — sem AVISITN, tem ASTDT/ASTDY
#   adbdc.sas7bdat   — sem AVISITN, tem ASTDT
#   addisp.sas7bdat  — sem AVISITN, tem ASTDT
#   adice.sas7bdat   — sem AVISITN, tem ADT/ADY
#   adrxfail.sas7bdat — sem AVISITN (cross-sectional)
#   adrxhist.sas7bdat — sem AVISITN (cross-sectional)
#   adsl.sas7bdat    — sem AVISITN (1 linha por sujeito)
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Flags de análise específicos do GALAXI-3
# (diferem do GAL1/GAL2 que usam FASFL/SAFFL/USMFL)
# ---------------------------------------------------------------
ANALYSIS_FLAGS_G3 <- c("RANDFL", "ATASFL", "ATSAFFL", "PASFL", "PSAFFL", "MPASFL")

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
  # Normaliza esperados para comparação (sem extensão)
  expected_norm <- tolower(tools::file_path_sans_ext(EXPECTED_DATASETS))
  missing_ds <- setdiff(expected_norm, names(datasets))
  if (length(missing_ds) > 0) {
    log_warn(paste("Datasets esperados não encontrados:", paste(missing_ds, collapse = ", ")))
  }

  log_info(paste("Carregados", length(datasets), "datasets"))
  return(datasets)
}

# ===============================
# 2. ADD PSEUDO VISIT (CRÍTICO)
# ===============================
add_pseudo_visit <- function(df, dataset_name) {

  if (!("USUBJID" %in% names(df))) return(df)

  if (!("AVISITN" %in% names(df))) {

    log_warn(paste("Criando pseudo AVISITN para:", dataset_name))

    # Ordernar por data disponível para preservar sequência temporal
    date_col <- intersect(c("ADT", "ASTDT", "ASTDTC"), names(df))[1]

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
# Âncora primária: advscdai.sas7bdat (309.550 linhas — maior dataset,
# endpoints centrais: CDAI, AP, SF, PRO2, CREM*, ACDAI* — GAL3 inclui
# PARAMCDs extras como ACDAIPG1, ACREMG1, CREMG3 vs GAL1/GAL2).
# Complementada com adsl.sas7bdat para garantir todos os 525 sujeitos.
build_spine <- function(datasets) {

  spine_parts <- list()

  anchor <- "advscdai.sas7bdat"
  if (anchor %in% names(datasets)) {
    spine_parts[[anchor]] <- datasets[[anchor]] %>%
      select(any_of(c("USUBJID", "AVISIT", "AVISITN"))) %>%
      distinct()
    log_info(paste("Usando", anchor, "como âncora da spine"))
  } else {
    log_warn("Âncora advscdai.sas7bdat não encontrada — spine construída de todos os datasets")
  }

  # Complemento: demais datasets com AVISITN real
  for (nm in names(datasets)) {
    if (nm == anchor) next
    df <- datasets[[nm]]
    if (all(c("USUBJID", "AVISITN") %in% names(df))) {
      spine_parts[[nm]] <- df %>%
        select(any_of(c("USUBJID", "AVISIT", "AVISITN"))) %>%
        distinct()
    }
  }

  spine <- bind_rows(spine_parts) %>% distinct()

  # Garantir cobertura total de sujeitos via ADSL
  adsl_key <- "adsl.sas7bdat"
  if (adsl_key %in% names(datasets)) {
    subj_adsl <- datasets[[adsl_key]] %>%
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
# 4. CLASSIFICATION
# ===============================
classify_variable <- function(df, var) {

  tmp <- df %>%
    select(USUBJID, all_of(var)) %>%
    distinct()

  variability <- tmp %>%
    group_by(USUBJID) %>%
    summarise(n = n_distinct(.data[[var]]), .groups = "drop")

  if (any(variability$n > 1, na.rm = TRUE)) {
    return("LONGITUDINAL")
  } else {
    return("STATIC")
  }
}

# ===============================
# 5. SMART AGGREGATION
# ===============================
aggregate_value <- function(x) {

  if (is.numeric(x)) {
    return(mean(x, na.rm = TRUE))
  }

  if (is.character(x) || is.factor(x)) {
    return(first(na.omit(x)))
  }

  return(first(x))
}

# ===============================
# 6. EXTRACT VARIABLES
# ===============================
extract_dataset <- function(df, dataset_name) {

  vars <- setdiff(names(df), c("USUBJID", "AVISIT", "AVISITN"))

  extracted <- list()
  used_vars  <- c()

  # Nome limpo para prefixo de coluna (remove sufixo .sas7bdat)
  clean_name <- toupper(gsub("\\.sas7bdat$", "", dataset_name))

  for (v in vars) {

    class_type <- classify_variable(df, v)

    if (class_type == "LONGITUDINAL") {

      tmp <- df %>%
        select(USUBJID, AVISITN, value = all_of(v)) %>%
        group_by(USUBJID, AVISITN) %>%
        summarise(value = aggregate_value(value), .groups = "drop")

      colname <- paste0("AVAL_", clean_name, "_", v)
      names(tmp)[3] <- colname

      extracted[[colname]] <- tmp
      used_vars <- c(used_vars, v)

    } else {

      tmp <- df %>%
        select(USUBJID, value = all_of(v)) %>%
        group_by(USUBJID) %>%
        summarise(value = aggregate_value(value), .groups = "drop")

      colname <- paste0("STC_", clean_name, "_", v)
      names(tmp)[2] <- colname

      extracted[[colname]] <- tmp
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
    log_warn(paste("Duplicatas USUBJID × AVISITN detectadas:", nrow(dup_check), "combinações"))
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

  log_info("=== GALAXI-3 ARD PIPELINE INICIADO ===")

  # 1. Carregar
  datasets <- load_datasets(DATA_PATH)

  # 2. Pseudo-visita para datasets sem AVISITN
  datasets <- imap(datasets, add_pseudo_visit)

  # 3. Spine baseada em advscdai.sas7bdat + adsl.sas7bdat
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
