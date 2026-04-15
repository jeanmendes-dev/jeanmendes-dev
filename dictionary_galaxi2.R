# ============================================================
# ML FEATURE DICTIONARY GENERATOR
# ============================================================

library(data.table)
library(dplyr)
library(stringr)
library(purrr)

# ===============================
# CONFIG
# ===============================
INPUT_ARD <- "/mnt/galaxi2_ml_dataset_v2.csv"
OUTPUT_DICT <- "/mnt/galaxi2_ml_dictionary.csv"

# ===============================
# LOG
# ===============================
log_info <- function(msg) message(paste0("[INFO] ", msg))

# ===============================
# 1. LOAD DATA
# ===============================
log_info("Loading ARD dataset...")

df <- fread(INPUT_ARD)

cols <- names(df)

# ===============================
# 2. PARSE VARIABLE NAME
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
# 3. GENERATE MEANING (AUTO)
# ===============================
generate_meaning <- function(type, dataset, var) {
  
  # Clean var name
  clean_var <- str_replace_all(var, "_", " ")
  
  # Rules
  if (type == "AVAL") {
    return(paste("Valor longitudinal de", clean_var, "no dataset", toupper(dataset)))
  }
  
  if (type == "STC") {
    return(paste("Característica estática do paciente:", clean_var, "(dataset", toupper(dataset), ")"))
  }
  
  return(paste("Variável identificadora:", clean_var))
}

meta <- meta %>%
  mutate(
    meaning = pmap_chr(
      list(type, dataset, original_var),
      generate_meaning
    )
  )

# ===============================
# 4. DETECT VARIABLE TYPE (DATA)
# ===============================
detect_data_type <- function(vec) {
  if (is.numeric(vec)) return("numeric")
  if (is.character(vec)) return("categorical")
  if (is.logical(vec)) return("binary")
  return(class(vec)[1])
}

data_types <- map_chr(df, detect_data_type)

meta$data_type <- data_types[meta$variable]

# ===============================
# 5. ADD BASIC STATS (ML GOLD)
# ===============================
log_info("Calculating statistics...")

get_stats <- function(vec) {
  
  if (is.numeric(vec)) {
    return(tibble(
      missing_pct = mean(is.na(vec)) * 100,
      unique_n = n_distinct(vec),
      mean = mean(vec, na.rm = TRUE),
      sd = sd(vec, na.rm = TRUE)
    ))
  } else {
    return(tibble(
      missing_pct = mean(is.na(vec)) * 100,
      unique_n = n_distinct(vec),
      mean = NA,
      sd = NA
    ))
  }
}

stats <- map_df(df, get_stats)

meta <- bind_cols(meta, stats)

# ===============================
# 6. FINAL OUTPUT
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

# ===============================
# SAVE
# ===============================
fwrite(dictionary, OUTPUT_DICT)

log_info(paste("Feature dictionary saved to:", OUTPUT_DICT))