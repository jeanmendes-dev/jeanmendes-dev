# ============================================================
# ADVANCED ML FEATURE DICTIONARY (CLINICAL-AWARE)
# ============================================================

library(data.table)
library(dplyr)
library(stringr)
library(purrr)

# ===============================
# CONFIG
# ===============================
INPUT_ARD <- "/mnt/galaxi2_ml_dataset_v2.csv"
DATA_PATH <- "/domino/datasets/local/clinical-trial-data/CNTO1959CRD3001-GALAXI-GAL2-WK48/load-2490/Data/csv"
OUTPUT_DICT <- "/mnt/galaxi2_feature_dictionary_v2.csv"

# ===============================
# LOAD DATA
# ===============================
df <- fread(INPUT_ARD)
cols <- names(df)

# ===============================
# LOAD SOURCE DATASETS (for PARAM)
# ===============================
files <- list.files(DATA_PATH, pattern = "\\.csv$", full.names = TRUE)

datasets <- lapply(files, function(f) fread(f))
names(datasets) <- tolower(tools::file_path_sans_ext(basename(files)))

# ===============================
# PARSE VARIABLE
# ===============================
parse_var <- function(var_name) {
  
  if (!str_detect(var_name, "^(AVAL|STC)_")) {
    return(tibble(
      variable = var_name,
      type = "ID",
      dataset = NA,
      original_var = var_name
    ))
  }
  
  parts <- str_split(var_name, "_", simplify = TRUE)
  
  tibble(
    variable = var_name,
    type = parts[1],
    dataset = parts[2],
    original_var = paste(parts[-c(1,2)], collapse = "_")
  )
}

meta <- map_df(cols, parse_var)

# ===============================
# BUILD PARAM DICTIONARY
# ===============================
param_lookup <- map_df(names(datasets), function(ds_name) {
  
  df_ds <- datasets[[ds_name]]
  
  if (!all(c("PARAMCD", "PARAM") %in% names(df_ds))) return(NULL)
  
  df_ds %>%
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
    AGE = "Age of subject",
    SEX = "Sex of subject",
    RACE = "Race of subject",
    COUNTRY = "Country of subject",
    TRT01A = "Actual treatment arm",
    TRT01P = "Planned treatment arm",
    AESEV = "Severity of adverse event",
    AESER = "Serious adverse event flag",
    AETERM = "Reported adverse event term"
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
get_param_meaning <- function(dataset, var) {
  
  ds <- param_lookup %>% filter(dataset == dataset)
  
  if (nrow(ds) == 0) return(NA)
  
  match <- ds %>% filter(PARAMCD == var)
  
  if (nrow(match) > 0) {
    return(match$PARAM[1])
  }
  
  return(NA)
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
  if (is.numeric(vec)) return("numeric")
  if (is.character(vec)) return("categorical")
  if (is.logical(vec)) return("binary")
  return(class(vec)[1])
}

meta$data_type <- map_chr(df, detect_data_type)[meta$variable]

# ===============================
# STATS
# ===============================
get_stats <- function(vec) {
  if (is.numeric(vec)) {
    return(c(
      missing_pct = mean(is.na(vec)) * 100,
      unique_n = n_distinct(vec),
      mean = mean(vec, na.rm = TRUE),
      sd = sd(vec, na.rm = TRUE)
    ))
  } else {
    return(c(
      missing_pct = mean(is.na(vec)) * 100,
      unique_n = n_distinct(vec),
      mean = NA,
      sd = NA
    ))
  }
}

stats <- as.data.frame(t(sapply(df, get_stats)))

meta <- bind_cols(meta, stats)

# ===============================
# FINAL
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
  )

fwrite(dictionary, OUTPUT_DICT)

cat("\n✅ Feature dictionary (ENRICHED) saved at:", OUTPUT_DICT, "\n")