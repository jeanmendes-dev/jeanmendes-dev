# ============================================================
# CLINICAL ADaM → ML-READY ARD PIPELINE
# STUDY: CNTO1275UCO3001-UNIFI / load-1901
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
DATA_PATH    <- "/domino/datasets/local/clinical-trial-data/CNTO1275UCO3001-UNIFI/load-1901/data"

OUTPUT_ARD      <- "/mnt/unifi_load1901_ml_dataset_v1.csv"
OUTPUT_COVERAGE <- "/mnt/unifi_load1901_ml_coverage_report_v1.csv"

# ---------------------------------------------------------------
# Datasets confirmados no diagnóstico UNIFI load-1901
# ---------------------------------------------------------------
EXPECTED_DATASETS <- c(
  "adae.sas7bdat",
  "adbdc.sas7bdat",
  "adchem.sas7bdat",
  "adcm.sas7bdat",
  "adcort.sas7bdat",
  "addeath.sas7bdat",
  "adds.sas7bdat",
  "addv.sas7bdat",
  "adecon.sas7bdat",
  "adeq5d.sas7bdat",
  "adex.sas7bdat",
  "adhema.sas7bdat",
  "adho.sas7bdat",
  "adibdq.sas7bdat",
  "adis.sas7bdat",
  "adlbef.sas7bdat",
  "admalig.sas7bdat",   # NOVO vs load-1899: dataset de malignidades
  "admayo2.sas7bdat",   # NOVO vs load-1899: substitui admayo.sas7bdat
  "adpc.sas7bdat",
  "adsaf.sas7bdat",
  "adsf36.sas7bdat",
  "adsg.sas7bdat",
  "adsl.sas7bdat",
  "adtfi.sas7bdat",
  "adtfm.sas7bdat",
  "aduceis.sas7bdat",
  "advscort.sas7bdat"
)

# ---------------------------------------------------------------
# Ausentes no load-1901 vs load-1899
#   adeff.sas7bdat    — não presente neste load
#   adhist.sas7bdat   — não presente neste load
#   adhist1.sas7bdat  — não presente neste load
#   adhist2.sas7bdat  — não presente neste load
#   admh.sas7bdat     — não presente neste load
#   adtte.sas7bdat    — não presente neste load
#   advs.sas7bdat     — não presente neste load
#
# Novos no load-1901 vs load-1899
#   admalig.sas7bdat  — malignidades (67 linhas; sem AVISITN)
#   admayo2.sas7bdat  — substitui admayo: 1.079.244 linhas, PARAMCDs
#                       de remissão/resposta Mayo com sufixos 2/4
#                       (RMGL2, RMGL4, RESPE2, RESPE4, EHL2, EHL4…)
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Datasets SEM AVISITN confirmados no diagnóstico load-1901
#   adae.sas7bdat     — sem AVISITN; tem ASTDT/ASTDY
#   adbdc.sas7bdat    — sem AVISITN; tem ASTDT (baseline cross-sectional)
#   adcm.sas7bdat     — sem AVISITN; cross-sectional de medicações
#   adcort.sas7bdat   — sem AVISITN; tem ADT/ADY (corticosteroid daily)
#   addeath.sas7bdat  — sem AVISITN; 3 linhas apenas
#   adds.sas7bdat     — sem AVISITN; disposição de sujeitos
#   addv.sas7bdat     — sem AVISITN; tem DVSTDT
#   adho.sas7bdat     — sem AVISITN; tem ASTDT
#   admalig.sas7bdat  — sem AVISITN; tem ASTDT (67 linhas)
#   adsg.sas7bdat     — sem AVISITN; tem ASTDT
#   adsl.sas7bdat     — sem AVISITN; 1 linha por sujeito (1.331 sujeitos)
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Particularidades do load-1901 vs load-1899
#
# 1. ÂNCORA DA SPINE — ADMAYO2
#    admayo2.sas7bdat: 1.079.244 linhas — maior dataset deste load,
#    PARAMCDs de remissão/resposta Mayo com sufixos temporais 2/4
#    (ex: RMGL2, RMGL4, RESPE2, RESPE4, EHL2, EHL4, MMRESP2/4).
#    Substitui admayo.sas7bdat do load-1899.
#
# 2. ESTRUTURA DE TRATAMENTO EXPANDIDA
#    load-1901 tem até TRT04P/TRT04A e subgrupos adicionais:
#    T03ACG2N até T03ACG9N no ADAE (vs T03ACG1N apenas no 1899).
#    Flags extras: TRTLTEFL, PDAFL, ADAFL, DAFL.
#
# 3. FLAGS DE ANÁLISE — LOAD-1901
#    Mantém: SAFFL, SAFW8FL, SAF2FL, COMPLFL, COMP20FL,
#            ENTERMFL, ENEXTFL, RANDFL, RAND2FL
#    Adiciona: TRTLTEFL, PDAFL, ADAFL (ausentes no load-1899)
#
# 4. ADLBEF MUITO MAIOR
#    load-1899: 152.994 linhas | load-1901: 860.930 linhas
#    Novos PARAMCDs com sufixos TF1/TF2 (ex: CALPT1, CALPT2,
#    CRPLTF1, CRPLTF2, LTFT1, LTFT2).
#
# 5. ADCORT MUITO MAIOR
#    load-1899: 306.112 linhas | load-1901: 1.276.653 linhas
#    Inclui ANL02FL, ANL03FL e TF2ADT ausentes no 1899.
#
# 6. ADECON MUITO MAIOR
#    load-1899: 71.860 linhas | load-1901: 437.025 linhas
#    PARAMCDs com sufixos TF1/TF2 para dois períodos de follow-up.
#
# 7. ADVSCORT MUITO MAIOR
#    load-1899: 39.960 linhas | load-1901: 806.204 linhas
#    Inclui CSUSEFL, TFFL1-4 e PARAMCDs extras CORTLTF1-4.
# ---------------------------------------------------------------

ANALYSIS_FLAGS_1901 <- c(
  "SAFFL", "SAFW8FL", "SAF2FL", "COMPLFL", "COMP20FL",
  "ENTERMFL", "ENEXTFL", "RANDFL", "RAND2FL",
  "TRTLTEFL", "PDAFL", "ADAFL"   # novos no load-1901
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
  expected_norm <- tolower(tools::file_path_sans_ext(EXPECTED_DATASETS))
  missing_ds    <- setdiff(expected_norm, names(datasets))
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
add_pseudo_visit <- function(df, dataset_name) {

  if (!("USUBJID" %in% names(df))) return(df)

  if (!("AVISITN" %in% names(df))) {

    log_warn(paste("Criando pseudo AVISITN para:", dataset_name))

    # Prioridade de data para ordenação temporal
    date_col <- intersect(
      c("ADT", "ASTDT", "DVSTDT", "ASTDTC", "DTHDTC"),
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
# Âncora primária: admayo2.sas7bdat (1.079.244 linhas)
# PARAMCDs centrais UC com sufixos de período 2/4:
# RMGL2, RMGL4, RESPE2, RESPE4, EHL2, EHL4, MMRESP2, MMRESP4.
# Substitui admayo.sas7bdat (790.256 linhas) do load-1899.
# Complementada com adsl.sas7bdat para os 1.331 sujeitos.
build_spine <- function(datasets) {

  spine_parts <- list()

  anchor <- "admayo2.sas7bdat"
  if (anchor %in% names(datasets)) {
    spine_parts[[anchor]] <- datasets[[anchor]] %>%
      select(any_of(c("USUBJID", "AVISIT", "AVISITN"))) %>%
      distinct()
    log_info(paste("Usando", anchor, "como âncora da spine"))
  } else {
    log_warn("Âncora admayo2.sas7bdat não encontrada — spine construída de todos os datasets")
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

  # Garantir cobertura total de sujeitos via ADSL (1.331 sujeitos)
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
# 4. CLASSIFICATION (vectorized)
# ===============================
# Classifica TODAS as variáveis de um dataset em uma única passagem,
# evitando N chamadas de group_by separadas sobre datasets grandes
# (datasets como adlbef 860k linhas, adcort 1.27M linhas e
#  admayo2 1.08M linhas travariam com abordagem coluna-a-coluna).
classify_all_variables <- function(df, vars) {

  dt <- as.data.table(df[, c("USUBJID", vars), drop = FALSE])

  # n_distinct por USUBJID para cada variável — uma única passagem
  variability <- dt[, lapply(.SD, function(x) uniqueN(x, na.rm = TRUE)),
                    by = USUBJID, .SDcols = vars]

  # Variável é LONGITUDINAL se qualquer sujeito tem > 1 valor distinto
  is_longitudinal <- sapply(vars, function(v) any(variability[[v]] > 1, na.rm = TRUE))

  return(ifelse(is_longitudinal, "LONGITUDINAL", "STATIC"))
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

  extracted  <- list()
  used_vars  <- c()

  # Nome limpo para prefixo: remove sufixo .sas7bdat, uppercase
  clean_name <- toupper(gsub("\\.sas7bdat$", "", dataset_name))

  # Classificar todas as variáveis de uma vez (evita loop de group_by)
  class_map <- classify_all_variables(df, vars)

  long_vars   <- vars[class_map == "LONGITUDINAL"]
  static_vars <- vars[class_map == "STATIC"]

  log_info(paste0("  ", dataset_name, ": ", length(long_vars),
                  " longitudinais / ", length(static_vars), " estáticas"))

  # --- LONGITUDINAIS: agregação por USUBJID × AVISITN ---
  if (length(long_vars) > 0) {

    dt_long <- as.data.table(df[, c("USUBJID", "AVISITN", long_vars), drop = FALSE])

    # Agregar todas as colunas longitudinais em uma só operação
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

  log_info("=== UNIFI load-1901 ARD PIPELINE INICIADO ===")

  # 1. Carregar
  datasets <- load_datasets(DATA_PATH)

  # 2. Pseudo-visita para datasets sem AVISITN
  datasets <- imap(datasets, add_pseudo_visit)

  # 3. Spine baseada em admayo2.sas7bdat + adsl.sas7bdat
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
