# Verifica contagem de sujeitos
adsl <- read.csv("ADSL.csv")

cat("Sujeitos no ADSL       :", nrow(adsl), "\n")
cat("USUBJIDs únicos        :", n_distinct(adsl$USUBJID), "\n")
cat("Duplicatas             :", nrow(adsl) - n_distinct(adsl$USUBJID), "\n")

# Cruza com DM original
dm <- read.csv("dm.csv")
names(dm) <- toupper(names(dm))
cat("Sujeitos no DM         :", nrow(dm), "\n")
cat("Em ADSL mas não em DM  :", sum(!adsl$USUBJID %in% dm$USUBJID), "\n")
cat("Em DM mas não em ADSL  :", sum(!dm$USUBJID %in% adsl$USUBJID), "\n")