# ============================================================
# CLINICAL ADaM → ML-READY ARD PIPELINE
# ============================================================

# ===============================
# LIBRARIES
# ===============================
library(data.table)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(readr)

# ===============================
# CONFIG
# ===============================
DATA_PATH <- "/domino/datasets/local/clinical-trial-data/CNTO1959CRD3001-GALAXI-GAL2-WK48/load-2490/Data/csv"
OUTPUT_ARD <- "/mnt/galaxi2_ml_dataset.csv"
OUTPUT_COVERAGE <- "/mnt/galaxi2_ml_coverage_report.csv"

# ===============================
# LOGGING UTILITY
# ===============================
log_warn <- function(msg) {
  warning(paste0("[WARNING] ", msg))
}

log_info <- function(msg) {
  message(paste0("[INFO] ", msg))
}

# ===============================
# 1. LOAD ALL DATASETS
# ===============================
load_datasets <- function(path) {
  files <- list.files(path = path, pattern = "\\.csv$", full.names = TRUE)
  
  datasets <- lapply(files, function(f) {
    dt <- tryCatch(
      fread(f, sep = "auto"),
      error = function(e) {
        log_warn(paste("Failed to read:", f))
        return(NULL)
      }
    )
    return(as.data.frame(dt))
  })
  
  names(datasets) <- tolower(tools::file_path_sans_ext(basename(files)))
  
  datasets <- datasets[!sapply(datasets, is.null)]
  
  log_info(paste("Loaded", length(datasets), "datasets"))
  return(datasets)
}

# ===============================
# 2. BUILD MASTER SPINE
# ===============================
build_spine <- function(datasets) {
  
  visit_datasets <- datasets[sapply(datasets, function(df) {
    all(c("USUBJID", "AVISITN") %in% names(df))
  })]
  
  if (length(visit_datasets) == 0) {
    stop("No datasets with USUBJID and AVISITN found")
  }
  
  spine <- bind_rows(lapply(visit_datasets, function(df) {
    df %>%
      select(any_of(c("USUBJID", "AVISIT", "AVISITN"))) %>%
      distinct()
  })) %>%
    distinct()
  
  log_info(paste("Spine created with", nrow(spine), "rows"))
  
  return(spine)
}

# ===============================
# 3. VARIABLE CLASSIFICATION
# ===============================
classify_variable <- function(df, var) {
  
  if (!("USUBJID" %in% names(df))) return("IGNORE")
  
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
# 4. EXTRACT VARIABLES
# ===============================
extract_dataset <- function(df, dataset_name) {
  
  vars <- setdiff(names(df), c("USUBJID", "AVISIT", "AVISITN"))
  
  extracted <- list()
  used_vars <- c()
  
  for (v in vars) {
    
    class_type <- classify_variable(df, v)
    
    if (class_type == "IGNORE") next
    
    if (class_type == "LONGITUDINAL") {
      
      if (!all(c("USUBJID", "AVISITN") %in% names(df))) next
      
      tmp <- df %>%
        select(USUBJID, AVISITN, value = all_of(v)) %>%
        group_by(USUBJID, AVISITN) %>%
        summarise(value = first(value), .groups = "drop")
      
      colname <- paste0("AVAL_", dataset_name, "_", v)
      names(tmp)[3] <- colname
      
      extracted[[colname]] <- tmp
      used_vars <- c(used_vars, v)
      
    } else if (class_type == "STATIC") {
      
      tmp <- df %>%
        select(USUBJID, value = all_of(v)) %>%
        group_by(USUBJID) %>%
        summarise(value = first(value), .groups = "drop")
      
      colname <- paste0("STC_", dataset_name, "_", v)
      names(tmp)[2] <- colname
      
      extracted[[colname]] <- tmp
      used_vars <- c(used_vars, v)
    }
  }
  
  return(list(
    extracted = extracted,
    used_vars = used_vars,
    all_vars = vars
  ))
}

# ===============================
# 5. MERGE INTO ARD
# ===============================
merge_ard <- function(spine, extracted_all) {
  
  ard <- spine
  
  for (ds in extracted_all) {
    for (obj in ds$extracted) {
      
      if (all(c("USUBJID", "AVISITN") %in% names(obj))) {
        ard <- left_join(ard, obj, by = c("USUBJID", "AVISITN"))
      } else if ("USUBJID" %in% names(obj)) {
        ard <- left_join(ard, obj, by = "USUBJID")
      }
    }
  }
  
  # Check duplicates
  dup_check <- ard %>%
    count(USUBJID, AVISITN) %>%
    filter(n > 1)
  
  if (nrow(dup_check) > 0) {
    log_warn("Duplicate USUBJID × AVISITN detected")
  }
  
  return(ard)
}

# ===============================
# 6. COVERAGE REPORT
# ===============================
build_coverage <- function(extracted_all, dataset_names) {
  
  report <- map2_df(extracted_all, dataset_names, function(ds, name) {
    
    total <- length(ds$all_vars)
    extracted <- length(ds$used_vars)
    missing <- setdiff(ds$all_vars, ds$used_vars)
    
    tibble(
      dataset = name,
      total_vars = total,
      extracted_vars = extracted,
      missing_vars = total - extracted,
      extraction_rate = round(100 * extracted / total, 2),
      missing_var_names = paste(missing, collapse = ", ")
    )
  })
  
  print(report)
  write_csv(report, OUTPUT_COVERAGE)
  
  return(report)
}

# ===============================
# 7. MAIN PIPELINE
# ===============================
run_pipeline <- function() {
  
  datasets <- load_datasets(DATA_PATH)
  
  spine <- build_spine(datasets)
  
  extracted_all <- list()
  
  for (name in names(datasets)) {
    log_info(paste("Processing dataset:", name))
    
    res <- extract_dataset(datasets[[name]], name)
    extracted_all[[name]] <- res
  }
  
  ard <- merge_ard(spine, extracted_all)
  
  write_csv(ard, OUTPUT_ARD)
  log_info(paste("ARD saved to:", OUTPUT_ARD))
  
  coverage <- build_coverage(extracted_all, names(datasets))
  
  log_info("Pipeline completed successfully")
}

# ===============================
# RUN
# ===============================
run_pipeline()