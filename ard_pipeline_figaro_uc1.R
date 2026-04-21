# =============================================================================
# ARD PIPELINE — FIGARO-UC1 (SHP647UC301)
# Analysis-Ready Dataset for Machine Learning
# Grain: USUBJID × AVISITN
# Generated from: diagnosis_figaro_uc1.docx
# =============================================================================
# ARCHITECTURE:
#   Section 1 — Configuration & helpers
#   Section 2 — Data ingestion
#   Section 3 — Visit spine construction
#   Section 4 — Longitudinal feature extraction  (AVAL_xxx)
#   Section 5 — Static feature extraction        (STC_xxx)
#   Section 6 — ARD assembly (left-join merge)
#   Section 7 — Output
#   Section 8 — Coverage report
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# =============================================================================
# SECTION 1 — CONFIGURATION
# =============================================================================

## ---- Paths ------------------------------------------------------------------
DATA_DIR  <- "/domino/datasets/local/clinical-trial-data/SHP647UC301-FIGARO-UC1/data"
OUT_DIR   <- "."          # write outputs here; adjust as needed
ARD_FILE  <- file.path(OUT_DIR, "ard_figaro_uc1.csv")
RPT_FILE  <- file.path(OUT_DIR, "coverage_report_figaro_uc1.txt")

## ---- Dataset registry -------------------------------------------------------
# Authoritative source: diagnosis_figaro_uc1.docx
# Each entry declares the dataset type and all columns that exist in the file.
# type = "longitudinal"  → keyed on USUBJID + visit identifiers (VISITNUM/VISIT)
# type = "static"        → one row per USUBJID (ADSL)
# type = "event"         → event-level, no VISITNUM; will be aggregated to static counts/flags

DATASET_REGISTRY <- list(

  adsl = list(
    type = "static",
    cols = c("STUDYID","USUBJID","SUBJID","SITEID","AGE","AGEU","SEX","RACE",
             "ETHNIC","COUNTRY","ARMCD","ARM","ACTARMCD","ACTARM","RFSTDTC",
             "RFENDTC","RFXSTDTC","RFXENDTC","RFICDTC","DTHFL","DTHDTC",
             "INVNAM","BRTHDTC","TRTSDT","TRTEDT","RANDDT","TRTDUR",
             "ITTFL","SAFFL","PPROTFL","ENRLFL"),
    # Columns to exclude from the ARD (pure admin / redundant identifiers)
    exclude = c("STUDYID","SUBJID","SITEID","INVNAM","AGEU")
  ),

  adlb = list(
    type        = "longitudinal",
    cols        = c("STUDYID","USUBJID","ARM","ARMCD","LBSEQ","LBTESTCD",
                    "LBTEST","LBCAT","AVAL","AVALU","BASE","CHG","PCHG",
                    "ABLFL","LBORRES","LBORRESU","LBNRIND","LBSTNRLO",
                    "LBSTNRHI","LBDTC","ADT","ADY","VISITNUM","VISIT",
                    "SAFFL","LBBLFL","LBFAST","LBSPEC"),
    param_col   = "LBTESTCD",   # pivot: one column per PARAMCD value
    value_col   = "AVAL",
    visit_col   = "VISITNUM",
    avisit_col  = "VISIT",
    exclude     = c("STUDYID","LBSEQ","ARM","ARMCD","AVALU","LBDTC","ADT",
                    "ABLFL","LBORRES","LBORRESU","LBNRIND","LBSTNRLO",
                    "LBSTNRHI","SAFFL","LBBLFL","LBFAST","LBSPEC","ADY")
  ),

  advs = list(
    type        = "longitudinal",
    cols        = c("STUDYID","USUBJID","ARM","ARMCD","VSSEQ","VSTESTCD",
                    "VSTEST","AVAL","AVALU","BASE","CHG","PCHG","ABLFL",
                    "VSORRES","VSORRESU","VSDTC","ADT","ADY","VISITNUM",
                    "VISIT","SAFFL"),
    param_col   = "VSTESTCD",
    value_col   = "AVAL",
    visit_col   = "VISITNUM",
    avisit_col  = "VISIT",
    exclude     = c("STUDYID","VSSEQ","ARM","ARMCD","AVALU","VSDTC","ADT",
                    "ABLFL","VSORRES","VSORRESU","SAFFL","ADY","BASE",
                    "CHG","PCHG")
  ),

  adeg = list(
    type        = "longitudinal",
    cols        = c("STUDYID","USUBJID","ARM","ARMCD","EGSEQ","EGTESTCD",
                    "EGTEST","EGCAT","AVAL","AVALU","BASE","CHG","PCHG",
                    "ABLFL","EGORRES","EGORRESU","EGSTAT","EGDTC","ADT",
                    "ADY","VISITNUM","VISIT","SAFFL"),
    param_col   = "EGTESTCD",
    value_col   = "AVAL",
    visit_col   = "VISITNUM",
    avisit_col  = "VISIT",
    exclude     = c("STUDYID","EGSEQ","ARM","ARMCD","AVALU","EGDTC","ADT",
                    "ABLFL","EGORRES","EGORRESU","EGSTAT","SAFFL","ADY",
                    "BASE","CHG","PCHG","EGCAT")
  ),

  admi = list(
    type        = "longitudinal",
    cols        = c("STUDYID","USUBJID","ARM","ARMCD","MISEQ","MITESTCD",
                    "MITEST","MICAT","MISPEC","MILOC","MIORRES","AVAL",
                    "AVALU","BASE","CHG","ABLFL","MIDTC","ADT","ADY",
                    "VISITNUM","VISIT","SAFFL","MISTAT","MIREASND","MISTRESC"),
    param_col   = "MITESTCD",
    value_col   = "AVAL",
    visit_col   = "VISITNUM",
    avisit_col  = "VISIT",
    exclude     = c("STUDYID","MISEQ","ARM","ARMCD","AVALU","MIDTC","ADT",
                    "ABLFL","MIORRES","MISTAT","MIREASND","MISTRESC",
                    "SAFFL","ADY","BASE","CHG","MICAT","MISPEC","MITEST")
  ),

  admo = list(
    type        = "longitudinal",
    cols        = c("STUDYID","USUBJID","ARM","ARMCD","MOSEQ","MOTESTCD",
                    "MOTEST","MOCAT","MOLOC","MOORRES","AVAL","AVALU",
                    "BASE","CHG","ABLFL","MODTC","ADT","ADY","VISITNUM",
                    "VISIT","SAFFL","MOSTRESC"),
    param_col   = "MOTESTCD",
    value_col   = "AVAL",
    visit_col   = "VISITNUM",
    avisit_col  = "VISIT",
    exclude     = c("STUDYID","MOSEQ","ARM","ARMCD","AVALU","MODTC","ADT",
                    "ABLFL","MOORRES","MOSTRESC","SAFFL","ADY","BASE",
                    "CHG","MOCAT","MOLOC","MOTEST")
  ),

  adex = list(
    type        = "longitudinal",
    cols        = c("STUDYID","USUBJID","ARM","ARMCD","EXSEQ","EXTRT",
                    "EXDOSE","EXDOSU","EXDOSFRM","EXDOSFRQ","EXROUTE",
                    "AVAL","AVALU","ADURU","EXSTDTC","EXENDTC","ASTDT",
                    "AENDT","ASTDY","AENDY","VISITNUM","VISIT","SAFFL"),
    param_col   = NULL,         # no PARAMCD; AVAL is a single exposure measure
    value_col   = "AVAL",
    visit_col   = "VISITNUM",
    avisit_col  = "VISIT",
    exclude     = c("STUDYID","EXSEQ","ARM","ARMCD","AVALU","EXSTDTC",
                    "EXENDTC","ASTDT","AENDT","ASTDY","AENDY","SAFFL",
                    "ADURU","EXDOSU","EXDOSFRM","EXDOSFRQ","EXROUTE")
  ),

  # Event datasets: no VISITNUM; aggregated to per-subject counts/flags (→ STC_)
  adae = list(
    type = "event",
    cols = c("STUDYID","USUBJID","ARM","ARMCD","AESEQ","AETERM","AEDECOD",
             "AEBODSYS","AESOC","AESEV","AESER","AEREL","AEOUT","AEACN",
             "AESTDTC","AEENDTC","ASTDT","AENDT","ASTDY","AENDY","TRTEMFL",
             "ANL01FL","SAFFL","AECONTRT","AESPID","AEHLGT","AEHLGTCD",
             "AELLT","AELLTCD","AEPTCD","AEHLT","AEHLTCD"),
    # Aggregations to materialise as STC_ columns
    aggregations = list(
      N_AE_TOTAL    = ~ n(),
      N_AE_SERIOUS  = ~ sum(AESER  == "Y", na.rm = TRUE),
      N_AE_RELATED  = ~ sum(AEREL  == "Y", na.rm = TRUE),
      FLG_AE_TRTEM  = ~ as.integer(any(TRTEMFL == "Y", na.rm = TRUE))
    )
  ),

  adcm = list(
    type = "event",
    cols = c("STUDYID","USUBJID","ARM","ARMCD","CMSEQ","CMTRT","CMDECOD",
             "CMCLAS","CMCLASCD","CMDOSE","CMDOSU","CMDOSFRQ","CMROUTE",
             "CMSTDTC","CMENDTC","ASTDT","AENDT","ASTDY","AENDY","TRTEMFL",
             "SAFFL","CMINDC","CMCAT"),
    aggregations = list(
      N_CM_TOTAL   = ~ n(),
      FLG_CM_TRTEM = ~ as.integer(any(TRTEMFL == "Y", na.rm = TRUE))
    )
  )
)

# =============================================================================
# SECTION 2 — DATA INGESTION
# =============================================================================

#' Load all datasets listed in the registry from CSV files on disk.
#' Returns a named list of data.frames; warns if a file is missing.
load_datasets <- function(registry, data_dir) {
  ds_list <- list()
  for (nm in names(registry)) {
    path <- file.path(data_dir, paste0(toupper(nm), ".csv"))
    if (!file.exists(path)) {
      warning(sprintf("[LOAD] File not found: %s — dataset '%s' will be skipped.", path, nm))
      next
    }
    df <- as.data.frame(data.table::fread(path, sep = "auto", na.strings = c("", "NA", ".")))
    # Validate that all expected columns are present
    expected <- registry[[nm]]$cols
    missing_cols <- setdiff(expected, names(df))
    if (length(missing_cols) > 0) {
      warning(sprintf("[LOAD] %s: expected columns not found in file: %s",
                      toupper(nm), paste(missing_cols, collapse = ", ")))
    }
    ds_list[[nm]] <- df
    message(sprintf("[LOAD] %-8s — %d rows × %d cols", toupper(nm), nrow(df), ncol(df)))
  }
  ds_list
}

# =============================================================================
# SECTION 3 — VISIT SPINE CONSTRUCTION
# =============================================================================

#' Build the master USUBJID × AVISITN spine from ADSL + all longitudinal datasets.
#' ADSL provides the subject universe; longitudinal datasets add the visit grid.
build_spine <- function(datasets, registry) {

  # Subject universe from ADSL
  adsl <- datasets[["adsl"]]
  all_subjects <- unique(adsl[["USUBJID"]])
  message(sprintf("[SPINE] %d unique subjects from ADSL.", length(all_subjects)))

  # Collect all USUBJID × VISITNUM × VISIT combos from longitudinal datasets
  visit_rows <- list()
  for (nm in names(registry)) {
    cfg <- registry[[nm]]
    if (cfg$type != "longitudinal") next
    df  <- datasets[[nm]]
    if (is.null(df)) next
    vc  <- cfg$visit_col      # e.g. "VISITNUM"
    avc <- cfg$avisit_col     # e.g. "VISIT"
    if (!all(c("USUBJID", vc, avc) %in% names(df))) next

    sub <- df[, c("USUBJID", vc, avc), drop = FALSE]
    sub <- sub[!is.na(sub[[vc]]), ]
    names(sub) <- c("USUBJID", "AVISITN", "AVISIT")
    visit_rows[[nm]] <- sub
  }

  spine <- bind_rows(visit_rows) %>%
    distinct(USUBJID, AVISITN, AVISIT) %>%
    filter(USUBJID %in% all_subjects) %>%
    arrange(USUBJID, AVISITN)

  # Warn about subjects in spine but not in ADSL (should not happen)
  orphans <- setdiff(spine$USUBJID, all_subjects)
  if (length(orphans) > 0)
    warning(sprintf("[SPINE] %d subjects in longitudinal data not found in ADSL: %s",
                    length(orphans), paste(orphans, collapse = ", ")))

  message(sprintf("[SPINE] Final spine: %d rows (%d subjects × visits).",
                  nrow(spine), n_distinct(spine$USUBJID)))
  spine
}

# =============================================================================
# SECTION 4 — LONGITUDINAL FEATURE EXTRACTION  (→ AVAL_xxx)
# =============================================================================

#' Extract longitudinal features from a single parameterised dataset
#' (e.g. ADLB, ADVS, ADEG, ADMI, ADMO).
#' Pivots PARAMCD × AVAL wide, producing one AVAL_<PARAMCD> column per parameter.
#'
#' For non-parameterised longitudinal datasets (e.g. ADEX, param_col = NULL),
#' the single AVAL column is kept as-is.
extract_longitudinal <- function(df, cfg, ds_name) {

  vc  <- cfg$visit_col
  avc <- cfg$avisit_col
  pc  <- cfg$param_col
  val <- cfg$value_col

  required <- c("USUBJID", vc, val)
  if (!is.null(pc)) required <- c(required, pc)
  missing_req <- setdiff(required, names(df))
  if (length(missing_req) > 0) {
    warning(sprintf("[LONG] %s: required columns missing — %s. Skipping.",
                    ds_name, paste(missing_req, collapse = ", ")))
    return(NULL)
  }

  # --- Duplicate guard -------------------------------------------------------
  if (!is.null(pc)) {
    key_cols <- c("USUBJID", vc, pc)
  } else {
    key_cols <- c("USUBJID", vc)
  }
  dupes <- df %>%
    group_by(across(all_of(key_cols))) %>%
    filter(n() > 1) %>%
    nrow()
  if (dupes > 0) {
    warning(sprintf("[LONG] %s: %d duplicate rows on key {%s}. Taking first value per key.",
                    ds_name, dupes, paste(key_cols, collapse = " × ")))
    df <- df %>%
      group_by(across(all_of(key_cols))) %>%
      slice(1) %>%
      ungroup()
  }

  if (!is.null(pc)) {
    # Pivot wide: one column per PARAMCD level
    wide <- df %>%
      select(all_of(c("USUBJID", vc, avc, pc, val))) %>%
      pivot_wider(
        id_cols     = c("USUBJID", !!sym(vc), !!sym(avc)),
        names_from  = all_of(pc),
        values_from = all_of(val),
        values_fn   = function(x) {
          if (length(x) > 1) {
            warning(sprintf("[LONG] %s: multiple values for same key after dedup — taking mean.", ds_name))
            mean(x, na.rm = TRUE)
          } else x
        }
      )
    # Rename value columns: AVAL_<PARAMCD>
    param_cols <- setdiff(names(wide), c("USUBJID", vc, avc))
    names(wide)[names(wide) %in% param_cols] <-
      paste0("AVAL_", str_replace_all(param_cols, "[^A-Za-z0-9]", "_"))

  } else {
    # No param column: keep AVAL as-is, rename to AVAL_<DS_NAME>
    wide <- df %>%
      select(all_of(c("USUBJID", vc, avc, val))) %>%
      rename(!!paste0("AVAL_", toupper(ds_name)) := all_of(val))
  }

  # Standardise visit columns to AVISITN / AVISIT
  wide <- wide %>%
    rename(AVISITN = all_of(vc), AVISIT = all_of(avc))

  message(sprintf("[LONG] %-8s — %d AVAL_ columns extracted.",
                  toupper(ds_name), sum(str_starts(names(wide), "AVAL_"))))
  wide
}

# =============================================================================
# SECTION 5 — STATIC FEATURE EXTRACTION  (→ STC_xxx)
# =============================================================================

#' Extract static features from ADSL.
#' All non-excluded columns become STC_<COL> columns (one value per USUBJID).
extract_adsl_static <- function(df, cfg) {

  exclude <- c("USUBJID", cfg$exclude)
  keep    <- setdiff(names(df), exclude)
  keep    <- intersect(keep, cfg$cols)   # only columns declared in registry

  static <- df %>%
    select(USUBJID, all_of(keep)) %>%
    distinct(USUBJID, .keep_all = TRUE)

  # Check for subjects with multiple rows (should never happen in ADSL)
  dupes <- static %>% group_by(USUBJID) %>% filter(n() > 1) %>% nrow()
  if (dupes > 0)
    warning(sprintf("[STATIC] ADSL: %d duplicate USUBJIDs found after dedup.", dupes))

  # Rename to STC_ prefix
  feature_cols <- setdiff(names(static), "USUBJID")
  names(static)[names(static) %in% feature_cols] <-
    paste0("STC_", feature_cols)

  message(sprintf("[STATIC] ADSL     — %d STC_ columns extracted.", length(feature_cols)))
  static
}


# =============================================================================
# SECTION 6 — ARD ASSEMBLY
# =============================================================================

#' Merge all feature tables onto the spine using left joins.
#' The spine guarantees no subject-visit rows are dropped.
assemble_ard <- function(spine, longitudinal_list, static_list) {

  ard <- spine   # start from the full USUBJID × AVISITN spine

  # --- Join longitudinal features (AVAL_xxx) ---------------------------------
  for (nm in names(longitudinal_list)) {
    feat <- longitudinal_list[[nm]]
    if (is.null(feat)) next
    before <- nrow(ard)
    ard <- left_join(ard, feat, by = c("USUBJID", "AVISITN", "AVISIT"))
    after <- nrow(ard)
    if (after != before)
      warning(sprintf("[MERGE] Joining %s increased row count %d → %d. Check for VISITNUM duplicates.",
                      toupper(nm), before, after))
    message(sprintf("[MERGE] %-8s joined. ARD now %d cols.", toupper(nm), ncol(ard)))
  }

  # --- Join static features (STC_xxx) ----------------------------------------
  for (nm in names(static_list)) {
    feat <- static_list[[nm]]
    if (is.null(feat)) next
    before <- nrow(ard)
    ard <- left_join(ard, feat, by = "USUBJID")
    after <- nrow(ard)
    if (after != before)
      warning(sprintf("[MERGE] Joining %s (static) changed row count %d → %d. Check for duplicate USUBJIDs.",
                      toupper(nm), before, after))
    message(sprintf("[MERGE] %-8s (static) joined. ARD now %d cols.", toupper(nm), ncol(ard)))
  }

  # Final sort
  ard <- ard %>% arrange(USUBJID, AVISITN)

  message(sprintf("[ARD] Final ARD: %d rows × %d columns.", nrow(ard), ncol(ard)))
  ard
}

# =============================================================================
# SECTION 7 — OUTPUT
# =============================================================================

#' Write ARD to CSV using write.csv to guarantee comma separator (locale-safe).
write_ard <- function(ard, path) {
  write.csv(ard, file = path, row.names = FALSE, na = "")
  message(sprintf("[OUT] ARD written to: %s", path))
}

# =============================================================================
# SECTION 8 — COVERAGE REPORT
# =============================================================================

#' For every dataset in the registry, report:
#'   - Total columns available in the raw file
#'   - Columns extracted into the ARD (as AVAL_ or STC_)
#'   - Columns not extracted (gaps)
#'   - Extraction rate %
coverage_report <- function(ard, datasets, registry, save_path = NULL) {

  sep <- paste(rep("=", 70), collapse = "")
  lines <- c(sep,
             "  COVERAGE REPORT — FIGARO-UC1 ARD PIPELINE",
             sprintf("  Generated: %s", Sys.time()),
             sep)

  ard_cols <- names(ard)

  for (nm in names(registry)) {
    cfg <- registry[[nm]]
    df  <- datasets[[nm]]
    if (is.null(df)) {
      lines <- c(lines, sprintf("\n[%s]  *** FILE NOT LOADED — SKIPPED ***", toupper(nm)))
      next
    }

    all_cols      <- cfg$cols
    excluded_cols <- c(cfg$exclude, "STUDYID", "USUBJID")   # always excluded from features
    candidate_cols <- setdiff(all_cols, excluded_cols)

    if (cfg$type == "longitudinal") {
      # Count distinct PARAMCD values that became AVAL_ columns in the ARD
      pc <- cfg$param_col
      if (!is.null(pc) && pc %in% names(df)) {
        params          <- unique(df[[pc]])
        params          <- params[!is.na(params)]
        expected_avals  <- paste0("AVAL_", str_replace_all(params, "[^A-Za-z0-9]", "_"))
        extracted       <- intersect(expected_avals, ard_cols)
        not_extracted   <- setdiff(expected_avals, ard_cols)
        total_available <- length(params)
        total_extracted <- length(extracted)
      } else {
        # Single-value longitudinal (ADEX)
        col_nm         <- paste0("AVAL_", toupper(nm))
        total_available <- 1
        total_extracted <- if (col_nm %in% ard_cols) 1L else 0L
        not_extracted   <- if (total_extracted == 0) col_nm else character(0)
      }
      label <- "AVAL_"

    } else if (cfg$type == "static") {
      expected_stcs  <- paste0("STC_", candidate_cols)
      extracted      <- intersect(expected_stcs, ard_cols)
      not_extracted  <- setdiff(expected_stcs, ard_cols)
      total_available <- length(candidate_cols)
      total_extracted <- length(extracted)
      label <- "STC_"

    } else {  # event
      agg_nms        <- paste0("STC_", names(cfg$aggregations))
      extracted      <- intersect(agg_nms, ard_cols)
      not_extracted  <- setdiff(agg_nms, ard_cols)
      total_available <- length(agg_nms)
      total_extracted <- length(extracted)
      label <- "STC_ (aggregated)"
    }

    rate <- if (total_available > 0)
      round(100 * total_extracted / total_available, 1) else NA_real_

    lines <- c(lines,
               "",
               sprintf("┌─ %-8s  [%s]", toupper(nm), cfg$type),
               sprintf("│  Raw columns available : %d", length(all_cols)),
               sprintf("│  Feature candidates    : %d  (after excluding admin cols)", total_available),
               sprintf("│  Extracted (%s)       : %d", label, total_extracted),
               sprintf("│  Not extracted (gaps)  : %d", length(not_extracted)),
               sprintf("│  Extraction rate       : %.1f%%", rate))

    if (length(not_extracted) > 0) {
      gap_lines <- strwrap(paste(not_extracted, collapse = ", "), width = 60,
                           prefix = "│    ", initial = "│  Gaps: ")
      lines <- c(lines, gap_lines)
    }
    lines <- c(lines, "└" %+% paste(rep("─", 69), collapse = ""))
  }

  lines <- c(lines, "", sep,
             sprintf("  ARD SUMMARY: %d rows × %d columns", nrow(ard), ncol(ard)),
             sprintf("  AVAL_ columns: %d", sum(str_starts(ard_cols, "AVAL_"))),
             sprintf("  STC_  columns: %d", sum(str_starts(ard_cols, "STC_"))),
             sep)

  cat(paste(lines, collapse = "\n"), "\n")

  if (!is.null(save_path)) {
    writeLines(lines, save_path)
    message(sprintf("[REPORT] Coverage report saved to: %s", save_path))
  }

  invisible(lines)
}

# String concatenation helper (used in coverage report)
`%+%` <- function(a, b) paste0(a, b)

# =============================================================================
# MAIN — EXECUTION SEQUENCE
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("  ARD PIPELINE — FIGARO-UC1  |  ", Sys.time())
message(paste(rep("=", 60), collapse = ""), "\n")

## 1. Load all datasets --------------------------------------------------------
message("--- [1/6] Loading datasets ---")
datasets <- load_datasets(DATASET_REGISTRY, DATA_DIR)

## 2. Build visit spine --------------------------------------------------------
message("\n--- [2/6] Building USUBJID × AVISITN spine ---")
spine <- build_spine(datasets, DATASET_REGISTRY)

## 3. Extract longitudinal features --------------------------------------------
message("\n--- [3/6] Extracting longitudinal features (AVAL_) ---")
longitudinal_features <- list()

for (nm in names(DATASET_REGISTRY)) {
  cfg <- DATASET_REGISTRY[[nm]]
  if (cfg$type != "longitudinal") next
  df  <- datasets[[nm]]
  if (is.null(df)) next
  longitudinal_features[[nm]] <- extract_longitudinal(df, cfg, nm)
}

## 4. Extract static features --------------------------------------------------
message("\n--- [4/6] Extracting static features (STC_) ---")
static_features <- list()

# ADSL → STC_
if (!is.null(datasets[["adsl"]])) {
  static_features[["adsl"]] <- extract_adsl_static(datasets[["adsl"]], DATASET_REGISTRY[["adsl"]])
}

# Event datasets → aggregated STC_
# Each event dataset is aggregated to one row per USUBJID using the
# aggregation specs defined in the registry. We evaluate each formula
# expression inside a summarise() call, catching errors individually
# so a bad column in one aggregation does not abort the whole dataset.

aggregate_event_dataset <- function(df, cfg, ds_name) {
  col_nms  <- names(cfg$aggregations)
  agg_list <- vector("list", length(col_nms))
  names(agg_list) <- col_nms

  for (col_nm in col_nms) {
    expr_rhs <- rlang::f_rhs(cfg$aggregations[[col_nm]])
    agg_list[[col_nm]] <- tryCatch(
      rlang::expr(!!expr_rhs),
      error = function(e) {
        warning(sprintf("[EVENT] %s / %s — expression error: %s", ds_name, col_nm, e$message))
        rlang::expr(NA_real_)
      }
    )
  }

  result <- df %>%
    group_by(USUBJID) %>%
    summarise(!!!agg_list, .groups = "drop")

  # Rename feature columns with STC_ prefix
  feat_cols <- setdiff(names(result), "USUBJID")
  names(result)[names(result) %in% feat_cols] <- paste0("STC_", feat_cols)

  message(sprintf("[EVENT] %-8s — %d STC_ columns aggregated.",
                  toupper(ds_name), length(feat_cols)))
  result
}

for (nm in names(DATASET_REGISTRY)) {
  cfg <- DATASET_REGISTRY[[nm]]
  if (cfg$type != "event") next
  df  <- datasets[[nm]]
  if (is.null(df)) next
  static_features[[nm]] <- aggregate_event_dataset(df, cfg, nm)
}

## 5. Assemble ARD -------------------------------------------------------------
message("\n--- [5/6] Assembling ARD (left-join merge) ---")
ard <- assemble_ard(spine, longitudinal_features, static_features)

## 6. Write output & coverage report ------------------------------------------
message("\n--- [6/6] Writing output and coverage report ---")
write_ard(ard, ARD_FILE)
coverage_report(ard, datasets, DATASET_REGISTRY, save_path = RPT_FILE)

message("\n✓ Pipeline complete.")
