library(data.table)

arquivo <- "/domino/datasets/local/clinical-trial-data/CNT01275UC03001-UNIFI/load-1899/data/unifi_load1899_ml_dataset.csv"

param_interesse <- "STRICT"

# ============================================================
# 1) Ler só o cabeçalho / amostra para descobrir as colunas
# ============================================================
amostra <- fread(arquivo, nrows = 1000)

# Colunas que terminam com _PARAMCD
paramcd_cols <- grep("_PARAMCD$", names(amostra), value = TRUE)

cat("Colunas *_PARAMCD encontradas:\n")
print(paramcd_cols)

# ============================================================
# 2) Ver em quais colunas *_PARAMCD o valor procurado aparece
# ============================================================
colunas_com_match <- paramcd_cols[
  sapply(paramcd_cols, function(col) {
    any(toupper(trimws(amostra[[col]])) == toupper(trimws(param_interesse)), na.rm = TRUE)
  })
]

cat("\nColunas onde o PARAMCD apareceu na amostra:\n")
print(colunas_com_match)

# ============================================================
# 3) Se não aparecer na amostra, ler apenas colunas *_PARAMCD
#    do arquivo inteiro para procurar com mais confiança
# ============================================================
if (length(colunas_com_match) == 0) {
  
  cat("\nNao apareceu na amostra. Procurando no arquivo inteiro apenas nas colunas *_PARAMCD...\n")
  
  dt_paramcd <- fread(
    arquivo,
    select = unique(c("USUBJID", "AVISIT", "AVISITN", paramcd_cols))
  )
  
  colunas_com_match <- paramcd_cols[
    sapply(paramcd_cols, function(col) {
      any(toupper(trimws(dt_paramcd[[col]])) == toupper(trimws(param_interesse)), na.rm = TRUE)
    })
  ]
  
  cat("\nColunas onde o PARAMCD apareceu no arquivo:\n")
  print(colunas_com_match)
}

# ============================================================
# 4) Para cada coluna encontrada, montar colunas associadas
# ============================================================
if (length(colunas_com_match) == 0) {
  cat("\nNenhuma coluna *_PARAMCD contem o valor:", param_interesse, "\n")
  
} else {
  
  resultados <- list()
  
  for (paramcd_col in colunas_com_match) {
    
    prefixo <- sub("_PARAMCD$", "", paramcd_col)
    
    param_col <- paste0(prefixo, "_PARAM")
    aval_col  <- paste0(prefixo, "_AVAL")
    avalc_col <- paste0(prefixo, "_AVALC")
    
    cols_desejadas <- c("USUBJID", "AVISIT", "AVISITN", paramcd_col)
    
    if (param_col %in% names(amostra)) cols_desejadas <- c(cols_desejadas, param_col)
    if (aval_col  %in% names(amostra)) cols_desejadas <- c(cols_desejadas, aval_col)
    if (avalc_col %in% names(amostra)) cols_desejadas <- c(cols_desejadas, avalc_col)
    
    # Ler só as colunas necessárias
    dt <- fread(arquivo, select = unique(cols_desejadas))
    
    # Filtrar o PARAMCD de interesse
    dt_filtrado <- dt[
      toupper(trimws(get(paramcd_col))) == toupper(trimws(param_interesse))
    ]
    
    # Renomear para facilitar leitura
    setnames(
      dt_filtrado,
      old = intersect(c(paramcd_col, param_col, aval_col, avalc_col), names(dt_filtrado)),
      new = c(
        "MATCH_PARAMCD_COL",
        if (param_col %in% names(dt_filtrado)) "MATCH_PARAM" else NULL,
        if (aval_col  %in% names(dt_filtrado)) "MATCH_AVAL" else NULL,
        if (avalc_col %in% names(dt_filtrado)) "MATCH_AVALC" else NULL
      )
    )
    
    # Adicionar informação de origem
    dt_filtrado[, ORIGIN_PARAMCD_COLUMN := paramcd_col]
    dt_filtrado[, ORIGIN_PREFIX := prefixo]
    
    # Reordenar colunas
    ordem <- c(
      "ORIGIN_PREFIX",
      "ORIGIN_PARAMCD_COLUMN",
      "USUBJID",
      "AVISIT",
      "AVISITN",
      "MATCH_PARAMCD_COL",
      intersect(c("MATCH_PARAM", "MATCH_AVAL", "MATCH_AVALC"), names(dt_filtrado))
    )
    
    dt_filtrado <- dt_filtrado[, ..ordem]
    
    resultados[[paramcd_col]] <- dt_filtrado
  }
  
  resultado_final <- rbindlist(resultados, fill = TRUE)
  
  cat("\nResultado final:\n")
  print(head(resultado_final, 50))
  
  # Opcional: salvar
  # fwrite(resultado_final, paste0("busca_paramcd_", param_interesse, ".csv"))
}