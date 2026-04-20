# ── AJUSTA ESTES PATHS ──────────────────────────────────────────────────────
adsl_path <- "/domino/datasets/local/clinical-trial-data/SHP647UC301-FIGARO-UC1/data/adam_csv/ADSL.csv"
dm_path   <- "/domino/datasets/local/clinical-trial-data/SHP647UC301-FIGARO-UC1/data/sdtm_csv/dm.csv"
# ────────────────────────────────────────────────────────────────────────────

adsl <- read.csv(adsl_path, stringsAsFactors = FALSE)
dm   <- read.csv(dm_path,   stringsAsFactors = FALSE)
names(dm) <- toupper(names(dm))

cat("Sujeitos no ADSL       :", nrow(adsl), "\n")
cat("USUBJIDs únicos        :", n_distinct(adsl$USUBJID), "\n")
cat("Duplicatas             :", nrow(adsl) - n_distinct(adsl$USUBJID), "\n")
cat("Sujeitos no DM         :", nrow(dm), "\n")
cat("Em ADSL mas não em DM  :", sum(!adsl$USUBJID %in% dm$USUBJID), "\n")
cat("Em DM mas não em ADSL  :", sum(!dm$USUBJID %in% adsl$USUBJID), "\n")