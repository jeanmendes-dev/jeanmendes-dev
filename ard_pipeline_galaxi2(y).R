# ============================================================
# CLINICAL ADaM → ML-READY ARD PIPELINE (v2)
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
DATA_PATH <- "/domino/datasets/local/clinical-trial-data/CNTO1959CRD3001-GALAXI-GAL2-WK48/load-2490/Data/csv"

OUTPUT_ARD <- "/mnt/galaxi2_ml_dataset_v2.csv"
OUTPUT_COVERAGE <- "/mnt/galaxi2_ml_coverage_report_v2.csv"

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
      log_warn(paste("Failed to read:", f))
      return(NULL)
    })
  })
  
  names(datasets) <- tolower(tools::file_path_sans_ext(basename(files)))
  datasets <- datasets[!sapply(datasets, is.null)]
  
  log_info(paste("Loaded", length(datasets), "datasets"))
  return(datasets)
}

# ===============================
# 2. ADD PSEUDO VISIT (CRÍTICO)
# ===============================
add_pseudo_visit <- function(df, dataset_name) {
  
  if (!("USUBJID" %in% names(df))) return(df)
  
  if (!("AVISITN" %in% names(df))) {
    
    log_warn(paste("Creating pseudo AVISITN for:", dataset_name))
    
    df <- df %>%
      arrange(USUBJID) %>%
      group_by(USUBJID) %>%
      mutate(
        AVISITN = row_number(),
        AVISIT = paste0("PSEUDO_", AVISITN)
      ) %>%
      ungroup()
  }
  
  return(df)
}

# ===============================
# 3. BUILD SPINE
# ===============================
build_spine <- function(datasets) {
  
  spine <- bind_rows(lapply(datasets, function(df) {
    if (all(c("USUBJID", "AVISITN") %in% names(df))) {
      df %>%
        select(any_of(c("USUBJID", "AVISIT", "AVISITN"))) %>%
        distinct()
    }
  })) %>%
    distinct()
  
  log_info(paste("Spine created with", nrow(spine), "rows"))
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
  used_vars <- c()
  
  for (v in vars) {
    
    class_type <- classify_variable(df, v)
    
    if (class_type == "LONGITUDINAL") {
      
      tmp <- df %>%
        select(USUBJID, AVISITN, value = all_of(v)) %>%
        group_by(USUBJID, AVISITN) %>%
        summarise(value = aggregate_value(value), .groups = "drop")
      
      colname <- paste0("AVAL_", dataset_name, "_", v)
      names(tmp)[3] <- colname
      
      extracted[[colname]] <- tmp
      used_vars <- c(used_vars, v)
      
    } else {
      
      tmp <- df %>%
        select(USUBJID, value = all_of(v)) %>%
        group_by(USUBJID) %>%
        summarise(value = aggregate_value(value), .groups = "drop")
      
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
    log_warn("Duplicate USUBJID × AVISITN detected")
  }
  
  return(ard)
}

# ===============================
# 8. COVERAGE REPORT
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
  fwrite(report, OUTPUT_COVERAGE)
  
  return(report)
}

# ===============================
# 9. MAIN PIPELINE
# ===============================
run_pipeline <- function() {
  
  datasets <- load_datasets(DATA_PATH)
  
  # 🔥 APPLY PSEUDO VISIT
  datasets <- imap(datasets, add_pseudo_visit)
  
  spine <- build_spine(datasets)
  
  extracted_all <- list()
  
  for (name in names(datasets)) {
    log_info(paste("Processing dataset:", name))
    
    res <- extract_dataset(datasets[[name]], name)
    extracted_all[[name]] <- res
  }
  
  ard <- merge_ard(spine, extracted_all)
  
  fwrite(ard, OUTPUT_ARD)
  log_info(paste("ARD saved to:", OUTPUT_ARD))
  
  build_coverage(extracted_all, names(datasets))
  
  log_info("Pipeline completed successfully")
}

# ===============================
# RUN
# ===============================
run_pipeline()