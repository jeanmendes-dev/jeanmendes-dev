library(data.table)
library(stringr)

# ============================================================
# CONFIG
# ============================================================
base_path <- "/domino/datasets/local/clinical-trial-data/SHP647UC301-FIGARO-UC1/data/csv"

# Termos relacionados ao MES / Mayo Endoscopic Score
mes_terms <- c(
  "MES",
  "MAYO",
  "MAYO ENDOSCOPIC",
  "MAYO ENDOSCOPIC SCORE",
  "ENDOSCOPIC SCORE",
  "MUCOSAL APPEARANCE",
  "ENDOSCOPY",
  "UCEIS"
)

# Padrões de colunas onde geralmente podem existir descrições clínicas
candidate_col_patterns <- c(
  "TEST", "TESTCD", "CAT", "SCAT", "OBJ", "TERM", "DESC", "NAM", "LABEL", "METHOD"
)

# Número de linhas para amostra inicial
sample_n <- 5000

# ============================================================
# FUNÇÕES AUXILIARES
# ============================================================
safe_read <- function(file, nrows = sample_n) {
  tryCatch(
    fread(
      file,
      nrows = nrows,
      sep = ",",
      header = TRUE,
      fill = TRUE,
      quote = "\"",
      encoding = "UTF-8"
    ),
    error = function(e) {
      message("Erro ao ler: ", basename(file), " -> ", e$message)
      return(NULL)
    }
  )
}

clean_chr <- function(x) {
  toupper(trimws(as.character(x)))
}

find_term_hits <- function(vec, terms) {
  vec2 <- clean_chr(vec)
  hits <- vec2[
    Reduce(`|`, lapply(terms, function(term) str_detect(vec2, fixed(toupper(term)))))
  ]
  hits <- unique(hits)
  hits[!is.na(hits) & hits != ""]
}

detect_candidate_cols <- function(df) {
  nm <- names(df)
  nm[grepl(paste(candidate_col_patterns, collapse = "|"), nm, ignore.case = TRUE)]
}

detect_score_cols <- function(df) {
  nm <- names(df)
  nm[grepl("(STRESN|STRESC|ORRES|ORRESU|AVAL|SCORE)", nm, ignore.case = TRUE)]
}

domain_from_file <- function(file) {
  toupper(gsub("\\.csv$", "", basename(file), ignore.case = TRUE))
}

# ============================================================
# LISTAR ARQUIVOS
# ============================================================
files <- list.files(base_path, pattern = "\\.csv$", full.names = TRUE)

if (length(files) == 0) {
  stop("Nenhum arquivo CSV encontrado em: ", base_path)
}

# ============================================================
# DISCOVERY
# ============================================================
summary_results <- list()
content_results <- list()

for (f in files) {
  cat("Scanning:", basename(f), "\n")
  
  df <- safe_read(f)
  if (is.null(df) || ncol(df) == 0) next
  
  nm <- names(df)
  domain <- domain_from_file(f)
  
  has_usubjid <- "USUBJID" %in% toupper(nm)
  has_visit   <- "VISIT" %in% toupper(nm)
  has_visitnum <- "VISITNUM" %in% toupper(nm)
  
  visit_col <- nm[toupper(nm) == "VISIT"][1]
  visitnum_col <- nm[toupper(nm) == "VISITNUM"][1]
  
  candidate_cols <- detect_candidate_cols(df)
  score_cols <- detect_score_cols(df)
  
  # Procurar MES em nomes de colunas
  mes_colname_hits <- nm[
    Reduce(`|`, lapply(mes_terms, function(term) {
      str_detect(toupper(nm), fixed(toupper(term)))
    }))
  ]
  
  # Procurar MES em conteúdo de colunas candidatas
  file_content_hits <- list()
  
  for (col in candidate_cols) {
    hits <- find_term_hits(df[[col]], mes_terms)
    if (length(hits) > 0) {
      file_content_hits[[col]] <- hits
    }
  }
  
  # Resumo por arquivo
  summary_results[[length(summary_results) + 1]] <- data.table(
    file = basename(f),
    domain = domain,
    n_rows_sampled = nrow(df),
    n_cols = ncol(df),
    has_usubjid = has_usubjid,
    has_visit = has_visit,
    has_visitnum = has_visitnum,
    visit_col = ifelse(is.na(visit_col), "", visit_col),
    visitnum_col = ifelse(is.na(visitnum_col), "", visitnum_col),
    candidate_cols = paste(candidate_cols, collapse = " | "),
    score_cols = paste(score_cols, collapse = " | "),
    mes_in_colnames = paste(mes_colname_hits, collapse = " | "),
    mes_hit_count_in_content = length(file_content_hits)
  )
  
  # Detalhes dos hits encontrados no conteúdo
  if (length(file_content_hits) > 0) {
    for (col in names(file_content_hits)) {
      vals <- file_content_hits[[col]]
      content_results[[length(content_results) + 1]] <- data.table(
        file = basename(f),
        domain = domain,
        column_name = col,
        matched_value = vals
      )
    }
  }
}

# ============================================================
# CONSOLIDAR
# ============================================================
discovery_summary <- rbindlist(summary_results, fill = TRUE)
discovery_summary <- discovery_summary[order(
  -has_usubjid, -has_visit, -has_visitnum, -mes_hit_count_in_content, file
)]

if (length(content_results) > 0) {
  discovery_hits <- unique(rbindlist(content_results, fill = TRUE))
} else {
  discovery_hits <- data.table(
    file = character(),
    domain = character(),
    column_name = character(),
    matched_value = character()
  )
}

# ============================================================
# MOSTRAR RESULTADOS
# ============================================================
cat("\n==================== DISCOVERY SUMMARY ====================\n")
print(discovery_summary)

cat("\n==================== DISCOVERY HITS ====================\n")
print(discovery_hits)

# ============================================================
# SALVAR RELATÓRIOS
# ============================================================
fwrite(discovery_summary, file.path(base_path, "mes_discovery_summary.csv"))
fwrite(discovery_hits, file.path(base_path, "mes_discovery_hits.csv"))

cat("\nArquivos gerados com sucesso:\n")
cat("- mes_discovery_summary.csv\n")
cat("- mes_discovery_hits.csv\n")