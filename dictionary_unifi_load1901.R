# ============================================================
# ADVANCED ML FEATURE DICTIONARY (CLINICAL-AWARE)
# STUDY: CNTO1275UCO3001-UNIFI / load-1901
# ============================================================

library(data.table)
library(dplyr)
library(stringr)
library(purrr)

# ===============================
# CONFIG
# ===============================
INPUT_ARD   <- "/mnt/unifi_load1901_ml_dataset_v1.csv"
DATA_PATH   <- "/domino/datasets/local/clinical-trial-data/CNTO1275UCO3001-UNIFI/load-1901/data"
OUTPUT_DICT <- "/mnt/unifi_load1901_feature_dictionary.csv"

# ---------------------------------------------------------------
# Diferenças load-1901 vs load-1899 relevantes para o dicionário
#
# Datasets presentes (27 vs 33 do 1899):
#   NOVO  — admalig  (malignidades; sem PARAMCD/PARAM)
#   NOVO  — admayo2  (âncora Mayo; substitui admayo do 1899)
#   AUSENTE — adeff, adhist, adhist1, adhist2, admh, adtte, advs
#
# Prefixos de coluna afetados:
#   AVAL_ADMAYO2_* / STC_ADMAYO2_* em vez de AVAL_ADMAYO_* / STC_ADMAYO_*
#
# Flags extras nas colunas: TRTLTEFL, PDAFL, ADAFL, DAFL
# TRT04P / TRT04A e subgrupos T03ACGxN até T03ACG9N
# ---------------------------------------------------------------

# ===============================
# LOAD ARD
# ===============================
df   <- fread(INPUT_ARD)
cols <- names(df)

# ===============================
# LOAD SOURCE DATASETS (for PARAM lookup)
# ===============================
files    <- list.files(DATA_PATH, pattern = "\\.csv$", full.names = TRUE)
datasets <- lapply(files, function(f) fread(f))
names(datasets) <- tolower(tools::file_path_sans_ext(basename(files)))

# ===============================
# PARSE VARIABLE
# ===============================
# Extrai prefixo (AVAL/STC), nome do dataset de origem e variável original
# a partir do padrão: <PREFIXO>_<DATASET>_<VARIAVEL>
parse_var <- function(var_name) {

  if (!str_detect(var_name, "^(AVAL|STC)_")) {
    return(tibble(
      variable     = var_name,
      type         = "ID",
      dataset      = NA_character_,
      original_var = var_name
    ))
  }

  parts <- str_split(var_name, "_", simplify = TRUE)

  tibble(
    variable     = var_name,
    type         = parts[1],
    dataset      = parts[2],
    original_var = paste(parts[-c(1, 2)], collapse = "_")
  )
}

meta <- map_df(cols, parse_var)

# ===============================
# BUILD PARAM DICTIONARY
# ===============================
# Percorre todos os datasets do load-1901 e extrai pares PARAMCD → PARAM
param_lookup <- map_df(names(datasets), function(ds_name) {

  df_ds <- datasets[[ds_name]]

  if (!all(c("PARAMCD", "PARAM") %in% names(df_ds))) return(NULL)

  df_ds %>%
    as.data.frame() %>%
    distinct(PARAMCD, PARAM) %>%
    mutate(dataset = ds_name)
})

# ===============================
# DOMAIN-SPECIFIC DEFINITIONS
# ===============================
define_variable_meaning <- function(type, dataset, var) {

  var_upper <- toupper(var)

  # ---- STANDARD VARIABLES ----
  standard_dict <- list(
    # Demographics / subject
    AGE      = "Age of subject",
    SEX      = "Sex of subject",
    RACE     = "Race of subject",
    COUNTRY  = "Country of subject",
    ETHNIC   = "Ethnicity of subject",
    REGION   = "Geographic region of subject",
    SITEID   = "Study site identifier",
    SUBJID   = "Subject identifier",

    # Treatment
    TRT01A   = "Actual treatment arm (period 1)",
    TRT01P   = "Planned treatment arm (period 1)",
    TRT02A   = "Actual treatment arm (period 2)",
    TRT02P   = "Planned treatment arm (period 2)",
    TRT03A   = "Actual treatment arm (period 3)",
    TRT03P   = "Planned treatment arm (period 3)",
    TRT04A   = "Actual treatment arm (period 4) [load-1901]",
    TRT04P   = "Planned treatment arm (period 4) [load-1901]",
    TR01AG1  = "Actual treatment arm group 1 (period 1)",
    TR01PG1  = "Planned treatment arm group 1 (period 1)",
    TR02AG1  = "Actual treatment arm group 1 (period 2)",
    TR02PG1  = "Planned treatment arm group 1 (period 2)",
    TR02AG2  = "Actual treatment arm group 2 (period 2)",
    TR02PG2  = "Planned treatment arm group 2 (period 2)",

    # Analysis flags
    SAFFL    = "Safety population flag",
    SAFW8FL  = "Safety (Week 8 induction) population flag",
    SAF2FL   = "Safety (maintenance) population flag",
    RANDFL   = "Randomised population flag",
    RAND2FL  = "Re-randomised (maintenance) population flag",
    COMPLFL  = "Completer population flag",
    COMP20FL = "Completer (20-week) population flag",
    ENTERMFL = "Entered maintenance flag",
    ENEXTFL  = "Extended maintenance flag",
    TRTLTEFL = "Treatment long-term extension flag [load-1901]",
    PDAFL    = "Placebo/dose adaptation flag [load-1901]",
    ADAFL    = "Adaptive dosing flag [load-1901]",
    DAFL     = "Dose adaptation flag [load-1901]",

    # Adverse events
    AESEV    = "Severity of adverse event",
    ASEV     = "Severity of adverse event (analysed)",
    AESER    = "Serious adverse event indicator",
    ASER     = "Serious adverse event indicator (analysed)",
    AETERM   = "Reported adverse event term",
    AEDECOD  = "Adverse event MedDRA preferred term",
    AEBODSYS = "Adverse event body system / SOC",
    AEREL    = "Adverse event relationship to study drug",
    AEOUT    = "Adverse event outcome",
    AEACN    = "Action taken with study drug",
    TRTEMFL  = "Treatment-emergent adverse event flag",

    # Labs
    AVAL     = "Analysis value",
    BASE     = "Baseline value",
    CHG      = "Change from baseline",
    PCHG     = "Percent change from baseline",
    ANRLO    = "Analysis normal range lower limit",
    ANRHI    = "Analysis normal range upper limit",
    ANRIND   = "Analysis normal range indicator",
    BNRIND   = "Baseline normal range indicator",
    LBTOXGR  = "Laboratory toxicity grade",
    ABLFL    = "Baseline record flag",

    # Malignancies (admalig — new in load-1901)
    NEOTYPE  = "Neoplasm type",
    MALTYPE  = "Malignancy type",
    INSITUFL = "In situ malignancy flag",
    MALFL    = "Confirmed malignancy flag",
    SEERAFL  = "SEER malignancy registry flag",

    # Dates
    ASTDT    = "Analysis start date",
    AENDT    = "Analysis end date",
    ADT      = "Analysis date",
    RANDDT   = "Randomisation date",
    TRTSDT   = "Treatment start date",
    TRTEDT   = "Treatment end date"
  )

  if (var_upper %in% names(standard_dict)) {
    base_meaning <- standard_dict[[var_upper]]
  } else {
    base_meaning <- str_replace_all(var, "_", " ")
  }

  # ---- ADD TYPE CONTEXT ----
  if (type == "AVAL") {
    return(paste(base_meaning, "- longitudinal measurement across visits"))
  }

  if (type == "STC") {
    return(paste(base_meaning, "- subject-level (static) characteristic"))
  }

  return(base_meaning)
}

# ===============================
# ENRICH WITH PARAM
# ===============================
# Para variáveis originadas de colunas PARAMCD, recupera o label PARAM completo.
# Nota: dataset nos prefixos está em uppercase (ex: ADMAYO2), param_lookup em lowercase.
get_param_meaning <- function(dataset_prefix, var) {

  if (is.na(dataset_prefix) || is.na(var)) return(NA_character_)

  # Normaliza para lowercase com sufixo .sas7bdat para bater com param_lookup
  ds_key <- paste0(tolower(dataset_prefix), ".sas7bdat")
  ds     <- param_lookup %>% filter(dataset == ds_key)

  if (nrow(ds) == 0) return(NA_character_)

  match <- ds %>% filter(PARAMCD == var)

  if (nrow(match) > 0) return(match$PARAM[1])

  return(NA_character_)
}

# ===============================
# GENERATE MEANING
# ===============================
meta <- meta %>%
  rowwise() %>%
  mutate(
    param_desc = get_param_meaning(dataset, original_var),
    meaning = case_when(
      !is.na(param_desc) & type == "AVAL" ~
        paste(param_desc, "- longitudinal clinical measurement"),

      !is.na(param_desc) & type == "STC" ~
        paste(param_desc, "- baseline/static clinical characteristic"),

      TRUE ~ define_variable_meaning(type, dataset, original_var)
    )
  ) %>%
  ungroup()

# ===============================
# DATA TYPE
# ===============================
detect_data_type <- function(vec) {
  if (is.numeric(vec))   return("numeric")
  if (is.character(vec)) return("categorical")
  if (is.logical(vec))   return("binary")
  return(class(vec)[1])
}

meta$data_type <- map_chr(df, detect_data_type)[meta$variable]

# ===============================
# STATS
# ===============================
get_stats <- function(vec) {
  if (is.numeric(vec)) {
    return(c(
      missing_pct = round(mean(is.na(vec)) * 100, 2),
      unique_n    = n_distinct(vec),
      mean        = round(mean(vec, na.rm = TRUE), 4),
      sd          = round(sd(vec, na.rm = TRUE), 4)
    ))
  } else {
    return(c(
      missing_pct = round(mean(is.na(vec)) * 100, 2),
      unique_n    = n_distinct(vec),
      mean        = NA_real_,
      sd          = NA_real_
    ))
  }
}

stats <- as.data.frame(t(sapply(df, get_stats)))

meta <- bind_cols(meta, stats)

# ===============================
# FINAL DICTIONARY
# ===============================
dictionary <- meta %>%
  select(
    variable,
    dataset,
    type,
    original_var,
    data_type,
    missing_pct,
    unique_n,
    mean,
    sd,
    meaning
  ) %>%
  arrange(dataset, type, variable)

fwrite(dictionary, OUTPUT_DICT)

cat("\n✅ Feature dictionary (load-1901) saved at:", OUTPUT_DICT, "\n")
cat("   Variables documented:", nrow(dictionary), "\n")
