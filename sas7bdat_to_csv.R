library(haven)
library(readr)
library(fs)
library(dplyr)

input_dir  <- "/domino/datasets/local/clinical-trial-data/SHP647UC301-FIGARO-UC1/data/sas7bdat"
output_dir <- "/domino/datasets/local/clinical-trial-data/SHP647UC301-FIGARO-UC1/data/csv"

dir_create(output_dir)

sas_files <- dir_ls(input_dir, regexp = "\\.sas7bdat$", type = "file")

log_conv <- list()

for (f in sas_files) {
  nome_base <- tolower(path_ext_remove(path_file(f)))
  out_file  <- path(output_dir, paste0(nome_base, ".csv"))
  
  status <- "OK"
  nrows_ <- NA_integer_
  ncols_ <- NA_integer_
  msg_   <- ""
  
  tryCatch({
    df <- read_sas(f)
    nrows_ <- nrow(df)
    ncols_ <- ncol(df)
    names(df) <- toupper(names(df))
    write_csv(df, out_file, na = "")
  }, error = function(e) {
    status <<- "ERROR"
    msg_ <<- conditionMessage(e)
  })
  
  log_conv[[length(log_conv) + 1]] <- data.frame(
    source_file = path_file(f),
    output_file = path_file(out_file),
    status = status,
    nrows = nrows_,
    ncols = ncols_,
    message = msg_,
    stringsAsFactors = FALSE
  )
}

log_df <- bind_rows(log_conv)
write_csv(log_df, path(output_dir, "_conversion_log.csv"))

cat("Conversão concluída.\n")
print(log_df)