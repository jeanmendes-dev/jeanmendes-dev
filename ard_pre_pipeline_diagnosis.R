# DIAGNÓSTICO — rode antes do pipeline
csv_files <- list.files(
  path = "/domino/datasets/local/clinical-trial-data/CNTO1959CRD3001-GALAXI-GAL3-WK48/load-2491/Data/csv",
  pattern = "\\.csv$", full.names = TRUE
)
datasets_auto <- lapply(csv_files, function(f) as.data.frame(data.table::fread(f, sep = "auto")))
names(datasets_auto) <- tolower(tools::file_path_sans_ext(basename(csv_files)))

for (nome in names(datasets_auto)) {
  cat("\n========", toupper(nome), "========\n")
  cat("Colunas:", paste(names(datasets_auto[[nome]]), collapse = ", "), "\n")
  cat("Linhas:", nrow(datasets_auto[[nome]]), "\n")
  if ("PARAMCD" %in% names(datasets_auto[[nome]])) {
    cat("PARAMCDs únicos:", paste(unique(datasets_auto[[nome]]$PARAMCD), collapse = ", "), "\n")
  }
}