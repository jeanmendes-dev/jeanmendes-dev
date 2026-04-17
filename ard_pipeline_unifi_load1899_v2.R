# ============================================================
# CLINICAL ADaM -> ML-READY ARD PIPELINE
# STUDY: CNTO1275UCO3001-UNIFI / load-1899
# OPTIMIZED + TYPE-STABLE VERSION
# ============================================================

library(data.table)
library(dplyr)
library(stringr)
library(purrr)

# ===============================
# CONFIG
# ===============================
DATA_PATH <- "/domino/datasets/local/clinical-trial-data/CNTO1275UCO3001-UNIFI/load-1899/data"

OUTPUT_ARD <- "/mnt/unifi_load1899_ml_dataset.csv"
OUTPUT_COVERAGE <- "/mnt/unifi_load1899_ml_coverage_report.csv"

# ===============================
# LOGGING
# ===============================
log_info <- function(msg) message(paste0("[INFO] ", msg))
log_warn <- function(msg) warning(paste0("[WARNING] ", msg))

# ===============================
# 1. LOAD DATASETS
# ===============================
load_datasets <- function(path) {
  files <- list.files(path = path, pattern = "\\.csv$", full.names = TRUE)
  
  datasets <- lapply(files, function(f) {
    tryCatch({
      fread(f, sep = "auto", showProgress = FALSE)
    }, error = function(e) {
      log_warn(paste("Failed to read:", f, "|", conditionMessage(e)))
      return(NULL)
    })
  })
  
  names(datasets) <- tolower(tools::file_path_sans_ext(basename(files)))
  datasets <- datasets[!sapply(datasets, is.null)]
  
  log_info(paste("Loaded", length(datasets), "datasets"))
  return(datasets)
}

# ===============================
# 2. ADD PSEUDO VISIT
# ===============================
add_pseudo_visit <- function(dt, dataset_name) {
  
  if (!("USUBJID" %in% names(dt))) return(dt)
  
  if (!("AVISITN" %in% names(dt))) {
    
    log_warn(paste("Creating pseudo AVISITN for:", dataset_name))
    
    date_candidates <- c("ADT", "ASTDT", "DVSTDT", "AENDT", "DTHDTC", "AESDT")
    date_col <- intersect(date_candidates, names(dt))[1]
    
    if (!is.na(date_col)) {
      setorderv(dt, c("USUBJID", date_col), na.last = TRUE)
    } else {
      setorderv(dt, "USUBJID", na.last = TRUE)
    }
    
    dt[, AVISITN := seq_len(.N), by = USUBJID]
    dt[, AVISIT := paste0("PSEUDO_", AVISITN)]
  }
  
  return(dt)
}

# ===============================
# 3. BUILD SPINE
# ===============================
build_spine <- function(datasets) {
  
  spine_parts <- list()
  
  if ("admayo.sas7bdat" %in% names(datasets)) {
    spine_parts[["admayo"]] <- unique(
      datasets[["admayo.sas7bdat"]][, .(USUBJID, AVISIT, AVISITN)]
    )
    log_info("Using admayo.sas7bdat as primary spine anchor")
  }
  
  for (nm in names(datasets)) {
    dt <- datasets[[nm]]
    if (all(c("USUBJID", "AVISITN") %in% names(dt))) {
      spine_parts[[nm]] <- unique(dt[, .(USUBJID, AVISIT, AVISITN)])
    }
  }
  
  spine <- unique(rbindlist(spine_parts, fill = TRUE, use.names = TRUE))
  
  if ("adsl.sas7bdat" %in% names(datasets)) {
    adsl_subj <- unique(datasets[["adsl.sas7bdat"]][, .(USUBJID)])
    adsl_subj[, `:=`(AVISITN = 0L, AVISIT = "BASELINE")]
    spine <- unique(rbindlist(list(spine, adsl_subj), fill = TRUE, use.names = TRUE))
  }
  
  setkey(spine, USUBJID, AVISITN)
  
  log_info(paste(
    "Spine created with", nrow(spine), "rows and",
    uniqueN(spine$USUBJID), "unique subjects"
  ))
  
  return(spine)
}

# ===============================
# 4. VARIABLE CLASSIFICATION
# ===============================
classify_all_variables <- function(dt, vars) {
  
  if (length(vars) == 0) {
    return(setNames(character(0), character(0)))
  }
  
  tmp <- dt[, c("USUBJID", vars), with = FALSE]
  
  variability <- tmp[, lapply(.SD, function(x) uniqueN(x, na.rm = TRUE)),
                     by = USUBJID, .SDcols = vars]
  
  result <- sapply(vars, function(v) {
    if (any(variability[[v]] > 1, na.rm = TRUE)) "LONGITUDINAL" else "STATIC"
  })
  
  return(result)
}

# ===============================
# 5. TYPE-STABLE AGGREGATION
# ===============================
aggregate_value <- function(x) {
  
  if (is.integer(x)) {
    x_non_na <- x[!is.na(x)]
    if (length(x_non_na) == 0) return(as.integer(NA))
    return(as.integer(x_non_na[1]))
  }
  
  if (is.double(x) || is.numeric(x)) {
    if (all(is.na(x))) return(NA_real_)
    return(as.numeric(mean(x, na.rm = TRUE)))
  }
  
  if (is.logical(x)) {
    x_non_na <- x[!is.na(x)]
    if (length(x_non_na) == 0) return(NA)
    return(as.logical(x_non_na[1]))
  }
  
  if (is.character(x)) {
    x_non_na <- x[!is.na(x)]
    if (length(x_non_na) == 0) return(NA_character_)
    return(as.character(x_non_na[1]))
  }
  
  if (is.factor(x)) {
    x_non_na <- x[!is.na(x)]
    if (length(x_non_na) == 0) return(NA_character_)
    return(as.character(x_non_na[1]))
  }
  
  x_non_na <- x[!is.na(x)]
  if (length(x_non_na) == 0) return(NA_character_)
  return(as.character(x_non_na[1]))
}

# ===============================
# 6. EXTRACT DATASET AS BLOCKS
# ===============================
extract_dataset_blocks <- function(dt, dataset_name) {
  
  vars <- setdiff(names(dt), c("USUBJID", "AVISIT", "AVISITN"))
  
  class_map <- classify_all_variables(dt, vars)
  
  long_vars <- names(class_map)[class_map == "LONGITUDINAL"]
  static_vars <- names(class_map)[class_map == "STATIC"]
  
  clean_name <- toupper(gsub("\\.sas7bdat$", "", dataset_name))
  
  log_info(paste0(
    dataset_name, ": ",
    length(long_vars), " longitudinal / ",
    length(static_vars), " static"
  ))
  
  long_block <- NULL
  static_block <- NULL
  
  if (length(long_vars) > 0) {
    long_block <- dt[, lapply(.SD, aggregate_value),
                     by = .(USUBJID, AVISITN),
                     .SDcols = long_vars]
    
    new_long_names <- c("USUBJID", "AVISITN",
                        paste0("AVAL_", clean_name, "_", long_vars))
    setnames(long_block, names(long_block), new_long_names)
    setkey(long_block, USUBJID, AVISITN)
  }
  
  if (length(static_vars) > 0) {
    static_block <- dt[, lapply(.SD, aggregate_value),
                       by = USUBJID,
                       .SDcols = static_vars]
    
    new_static_names <- c("USUBJID",
                          paste0("STC_", clean_name, "_", static_vars))
    setnames(static_block, names(static_block), new_static_names)
    setkey(static_block, USUBJID)
  }
  
  return(list(
    long_block = long_block,
    static_block = static_block,
    used_vars = c(long_vars, static_vars),
    all_vars = vars
  ))
}

# ===============================
# 7. MERGE ARD
# ===============================
merge_ard <- function(spine, extracted_all) {
  
  ard <- copy(spine)
  
  for (nm in names(extracted_all)) {
    ds <- extracted_all[[nm]]
    
    if (!is.null(ds$long_block)) {
      ard <- ds$long_block[ard, on = .(USUBJID, AVISITN)]
    }
    
    if (!is.null(ds$static_block)) {
      ard <- ds$static_block[ard, on = .(USUBJID)]
    }
    
    log_info(paste("Merged dataset block:", nm))
  }
  
  dup_check <- ard[, .N, by = .(USUBJID, AVISITN)][N > 1]
  
  if (nrow(dup_check) > 0) {
    log_warn(paste("Duplicate USUBJID x AVISITN found:", nrow(dup_check)))
  }
  
  log_info(paste("ARD final:", nrow(ard), "rows x", ncol(ard), "columns"))
  return(ard)
}

# ===============================
# 8. COVERAGE REPORT
# ===============================
build_coverage <- function(extracted_all) {
  
  report <- rbindlist(lapply(names(extracted_all), function(nm) {
    ds <- extracted_all[[nm]]
    total <- length(ds$all_vars)
    extracted <- length(ds$used_vars)
    missing <- setdiff(ds$all_vars, ds$used_vars)
    
    data.table(
      dataset = nm,
      total_vars = total,
      extracted_vars = extracted,
      missing_vars = total - extracted,
      extraction_rate = round(100 * extracted / total, 2),
      missing_var_names = paste(missing, collapse = ", ")
    )
  }), fill = TRUE)
  
  print(report)
  fwrite(report, OUTPUT_COVERAGE)
  log_info(paste("Coverage report saved to:", OUTPUT_COVERAGE))
  
  return(report)
}

# ===============================
# 9. MAIN
# ===============================
run_pipeline <- function() {
  
  log_info("=== UNIFI load-1899 ARD PIPELINE STARTED ===")
  
  datasets <- load_datasets(DATA_PATH)
  datasets <- imap(datasets, add_pseudo_visit)
  
  spine <- build_spine(datasets)
  
  extracted_all <- list()
  
  for (nm in names(datasets)) {
    log_info(paste("Processing dataset:", nm))
    extracted_all[[nm]] <- extract_dataset_blocks(datasets[[nm]], nm)
  }
  
  log_info("Starting final merge...")
  ard <- merge_ard(spine, extracted_all)
  
  log_info("Writing final CSV...")
  fwrite(ard, OUTPUT_ARD, showProgress = TRUE)
  log_info(paste("ARD saved to:", OUTPUT_ARD))
  
  build_coverage(extracted_all)
  
  log_info("=== PIPELINE COMPLETED SUCCESSFULLY ===")
  invisible(ard)
}

run_pipeline()