# ============================================================
# VISUALIZAR CSV GRANDE + FILTRAR PARAMCD DE INTERESSE
# ============================================================

library(DBI)
library(duckdb)

# ---- caminho do arquivo ----
arquivo <- "/domino/datasets/local/clinical-trial-data/CNT01275UC03001-UNIFI/load-1899/data/unifi_load1899_ml_dataset.csv"

# ---- conectar DuckDB em memória ----
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")

# ============================================================
# 1) VER UMA AMOSTRA DO ARQUIVO
# ============================================================

sample_df <- dbGetQuery(con, sprintf("
  SELECT *
  FROM read_csv_auto('%s', SAMPLE_SIZE=-1)
  LIMIT 20
", arquivo))

print(sample_df)

# Ver nomes das colunas
colunas <- dbGetQuery(con, sprintf("
  SELECT *
  FROM read_csv_auto('%s', SAMPLE_SIZE=-1)
  LIMIT 1
", arquivo))

print(names(colunas))

# ============================================================
# 2) FILTRO POR PARAMCD DE INTERESSE
#    Exemplo: STRICT em AVAL_ADBDC_PARAMCD
# ============================================================

param_interesse <- "STRICT"

resultado <- dbGetQuery(con, sprintf("
  SELECT
    USUBJID,
    AVISIT,
    AVISITN,
    AVAL_ADBDC_PARAMCD,
    AVAL_ADBDC_PARAM,
    AVAL_ADBDC_AVAL,
    AVAL_ADBDC_AVALC
  FROM read_csv_auto('%s', SAMPLE_SIZE=-1)
  WHERE UPPER(TRIM(AVAL_ADBDC_PARAMCD)) = UPPER(TRIM('%s'))
  LIMIT 200
", arquivo, param_interesse))

print(resultado)

# ============================================================
# 3) SE QUISER VER QUAIS PARAMCD EXISTEM NO ADBDC
# ============================================================

paramcd_disponiveis <- dbGetQuery(con, sprintf("
  SELECT DISTINCT AVAL_ADBDC_PARAMCD
  FROM read_csv_auto('%s', SAMPLE_SIZE=-1)
  WHERE AVAL_ADBDC_PARAMCD IS NOT NULL
  ORDER BY AVAL_ADBDC_PARAMCD
", arquivo))

print(paramcd_disponiveis)
