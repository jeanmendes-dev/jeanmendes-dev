# =============================================================================
# ADMAYO – Mayo Endoscopic Score (MES)
# Estudo : SHP647UC301 / FIGARO-UC1
# Input  : CSVs SDTM (MO, MI, QS*, SUPPMO, SUPPMI) + ADSL.csv já gerado
# Output : ADMAYO.csv  →  grain: USUBJID × VISIT (wide)
#          ADMAYO_LONG.csv → grain: USUBJID × VISIT × TESTCD (long, opcional)
#
# Lógica de detecção:
#   1. Varre MO procurando TESTCDs/CATs que contenham "MAYO" ou "MES"
#   2. Varre MI  procurando scores histológicos (Geboes, Robarts, NHI)
#   3. Varre QS  (qser, qscc, qscd) procurando scores já calculados
#   4. Pivota tudo para wide USUBJID × VISIT
#   5. Deriva flags clínicos: remissão, mucosal healing, resposta
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# ── CONFIGURAÇÃO ──────────────────────────────────────────────────────────────
sdtm_dir <- "/domino/datasets/local/clinical-trial-data/SHP647UC301-FIGARO-UC1/data/sdtm_csv"
adam_dir <- "/domino/datasets/local/clinical-trial-data/SHP647UC301-FIGARO-UC1/data/adam_csv"

dir.create(adam_dir, showWarnings = FALSE, recursive = TRUE)

# ── UTILITÁRIOS ───────────────────────────────────────────────────────────────
iso_date <- function(x) as.Date(substr(as.character(x), 1, 10))

calc_ady <- function(adt, trtsdt) {
  d <- as.numeric(adt - trtsdt)
  ifelse(is.na(d), NA_real_, ifelse(d >= 0, d + 1, d))
}

read_csv_upper <- function(fname) {
  path <- file.path(sdtm_dir, fname)
  if (!file.exists(path)) { message("  [SKIP] ", fname); return(NULL) }
  df <- read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  names(df) <- toupper(names(df))
  message(sprintf("  [READ] %-16s  %d × %d", fname, nrow(df), ncol(df)))
  df
}

pivot_supp <- function(supp_df, seq_var) {
  if (is.null(supp_df) || nrow(supp_df) == 0) return(NULL)
  needed <- c("USUBJID", "IDVARVAL", "QNAM", "QVAL")
  if (!all(needed %in% names(supp_df))) return(NULL)
  supp_df %>%
    select(USUBJID, IDVARVAL, QNAM, QVAL) %>%
    pivot_wider(names_from = QNAM, values_from = QVAL,
                values_fn = ~ .x[1]) %>%
    rename(!!seq_var := IDVARVAL) %>%
    mutate(across(everything(), as.character))
}

merge_supp <- function(parent, supp_wide, seq_var) {
  if (is.null(supp_wide)) return(parent)
  new_cols <- setdiff(names(supp_wide), c("USUBJID", seq_var))
  parent %>%
    mutate(across(all_of(seq_var), as.character)) %>%
    left_join(
      supp_wide %>% mutate(across(everything(), as.character)),
      by = c("USUBJID", seq_var)
    )
}

# ── LEITURA DO ADSL ───────────────────────────────────────────────────────────
message("\n========== LENDO ADSL ==========")
adsl_path <- file.path(adam_dir, "ADSL.csv")
if (!file.exists(adsl_path)) {
  stop("ADSL.csv não encontrado em adam_dir. Execute sdtm_to_adam.R primeiro.")
}
adsl <- read.csv(adsl_path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
names(adsl) <- toupper(names(adsl))
adsl_key <- adsl %>%
  select(USUBJID, TRTSDT, TRTEDT, ARM, ARMCD, SAFFL, ITTFL,
         any_of(c("PPROTFL", "RANDDT"))) %>%
  mutate(TRTSDT = as.Date(TRTSDT), TRTEDT = as.Date(TRTEDT))

message(sprintf("  ADSL: %d sujeitos", nrow(adsl_key)))

# ── LEITURA DE MO + SUPPMO ────────────────────────────────────────────────────
message("\n========== LENDO MO + SUPPMO ==========")
mo      <- read_csv_upper("mo.csv")
suppmo  <- read_csv_upper("suppmo.csv")

if (!is.null(suppmo)) {
  suppmo_wide <- pivot_supp(suppmo, "MOSEQ")
  if (!is.null(mo) && !is.null(suppmo_wide))
    mo <- merge_supp(mo, suppmo_wide, "MOSEQ")
}

# ── LEITURA DE MI + SUPPMI ────────────────────────────────────────────────────
message("\n========== LENDO MI + SUPPMI ==========")
mi      <- read_csv_upper("mi.csv")
suppmi  <- read_csv_upper("suppmi.csv")

if (!is.null(suppmi)) {
  suppmi_wide <- pivot_supp(suppmi, "MISEQ")
  if (!is.null(mi) && !is.null(suppmi_wide))
    mi <- merge_supp(mi, suppmi_wide, "MISEQ")
}

# ── LEITURA DE QS RELEVANTES ──────────────────────────────────────────────────
message("\n========== LENDO QS (endoscopia/mayo) ==========")
# qser = Endoscopic Remission; qscc = CCIS; qscd = CDAI scores
qser  <- read_csv_upper("qser.csv")
qscc  <- read_csv_upper("qscc.csv")
qscd  <- read_csv_upper("qscd.csv")

# =============================================================================
# BLOCO DE DIAGNÓSTICO — imprime estrutura real dos dados
# =============================================================================
message("\n========== DIAGNÓSTICO: estrutura dos dados ==========")

diag_domain <- function(df, name, test_col, result_col, loc_col = NULL, cat_col = NULL) {
  if (is.null(df) || nrow(df) == 0) {
    message(sprintf("  %s: vazio ou não encontrado", name)); return(invisible(NULL))
  }
  message(sprintf("\n  --- %s (%d obs) ---", name, nrow(df)))
  message("  Colunas: ", paste(names(df), collapse = ", "))

  if (test_col %in% names(df)) {
    tbl <- df %>%
      count(!!sym(test_col), name = "n") %>%
      arrange(desc(n))
    message(sprintf("  Valores únicos em %s (%d):", test_col, nrow(tbl)))
    print(as.data.frame(tbl), row.names = FALSE)
  }
  if (!is.null(cat_col) && cat_col %in% names(df)) {
    message(sprintf("  Valores únicos em %s:", cat_col))
    print(as.data.frame(df %>% distinct(!!sym(cat_col))), row.names = FALSE)
  }
  if (!is.null(loc_col) && loc_col %in% names(df)) {
    message(sprintf("  Valores únicos em %s:", loc_col))
    print(as.data.frame(df %>% distinct(!!sym(loc_col))), row.names = FALSE)
  }
  if (result_col %in% names(df)) {
    message(sprintf("  Distribuição de %s:", result_col))
    print(as.data.frame(df %>% count(!!sym(result_col))), row.names = FALSE)
  }
  # Visitas disponíveis
  visit_col <- intersect(c("VISIT", "AVISIT"), names(df))[1]
  if (!is.na(visit_col)) {
    message(sprintf("  Visitas (%s):", visit_col))
    print(as.data.frame(df %>% distinct(!!sym(visit_col))), row.names = FALSE)
  }
}

diag_domain(mo,   "MO",   "MOTESTCD", "MOORRES",  "MOLOC",  "MOCAT")
diag_domain(mi,   "MI",   "MITESTCD", "MIORRES",  "MILOC",  "MICAT")
diag_domain(qser, "QSER", "QSTESTCD", "QSORRES",  NULL,     "QSCAT")
diag_domain(qscc, "QSCC", "QSTESTCD", "QSORRES",  NULL,     "QSCAT")

# =============================================================================
# DETECÇÃO AUTOMÁTICA DE TESTCDs RELACIONADOS A MAYO/MES
# =============================================================================
message("\n========== DETECÇÃO DE TESTCDs MAYO/MES ==========")

# Padrões conhecidos no FIGARO-UC1 / SHP647 e estudos UC similares
MES_PATTERNS <- c(
  # MO domain — Mayo Endoscopic Score por segmento e global
  "MES",        # prefixo genérico
  "MAYO",       # alternativa
  "MESOVRL",    # overall MES
  "MESCSIG",    # sigmoid/rectum
  "MESCTRS",    # transverse
  "MESCDES",    # descending
  "MESCASC",    # ascending/cecum
  "MESCRC",     # cecum/right colon (alternativa)
  "MESCEC",     # cecum
  "MESRTCL",    # right colon
  "MESTCL",     # transverse colon
  "MESLCL",     # left colon
  "MESRCT",     # rectum
  "ENDOSCOPY",
  "MUCOS"       # mucosal appearance
)

# Padrões histológicos (MI)
HISTO_PATTERNS <- c(
  "GEBOES", "GBS",          # Geboes Score
  "ROBARTS", "RHI",         # Robarts Histopathology Index
  "NHI", "NANCY",           # Nancy Histological Index
  "HISTO", "HISTOL",
  "CRYPTO",                 # Cryptitis
  "NEUINF",                 # Neutrophil infiltration
  "EROSION", "ULCER"
)

detect_testcds <- function(df, test_col, patterns, domain_name) {
  if (is.null(df) || !(test_col %in% names(df))) return(character(0))
  all_tests <- unique(toupper(trimws(df[[test_col]])))
  # Match por substring
  matched <- all_tests[
    sapply(all_tests, function(t)
      any(sapply(patterns, function(p) grepl(p, t, ignore.case = TRUE)))
    )
  ]
  if (length(matched) > 0)
    message(sprintf("  [%s] TESTCDs detectados: %s", domain_name,
                    paste(matched, collapse = ", ")))
  else
    message(sprintf("  [%s] nenhum TESTCD matched — verifique DIAGNÓSTICO acima", domain_name))
  matched
}

mes_testcds_mo    <- detect_testcds(mo,   "MOTESTCD", MES_PATTERNS,    "MO")
mes_testcds_mi    <- detect_testcds(mi,   "MITESTCD", HISTO_PATTERNS,  "MI")
mes_testcds_qser  <- detect_testcds(qser, "QSTESTCD", MES_PATTERNS,    "QSER")
mes_testcds_qscc  <- detect_testcds(qscc, "QSTESTCD", MES_PATTERNS,    "QSCC")

# =============================================================================
# EXTRAÇÃO E PADRONIZAÇÃO DOS SCORES
# =============================================================================
message("\n========== EXTRAINDO SCORES ==========")

# Coluna de avaliador (central reader vs local) — preferir central se existir
get_eval_filter <- function(df, eval_col = "MOEVAL") {
  if (!(eval_col %in% names(df))) return(df)
  evals <- unique(toupper(trimws(df[[eval_col]])))
  central_vals <- evals[grepl("CENTRAL|INDEPEND|BLINDED|CORE LAB", evals)]
  if (length(central_vals) > 0) {
    message(sprintf("  Filtrando avaliador CENTRAL: %s", paste(central_vals, collapse=", ")))
    df %>% filter(toupper(trimws(!!sym(eval_col))) %in% central_vals)
  } else {
    df
  }
}

# ── MO: extrai linhas MES ────────────────────────────────────────────────────
mo_mes <- NULL
if (!is.null(mo) && length(mes_testcds_mo) > 0) {
  mo_mes <- mo %>%
    filter(toupper(trimws(MOTESTCD)) %in% mes_testcds_mo) %>%
    get_eval_filter("MOEVAL") %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ADT    = iso_date(MODTC),
      ADY    = calc_ady(ADT, TRTSDT),
      AVAL   = suppressWarnings(as.numeric(MOORRES)),
      AVALC  = as.character(MOORRES),
      ABLFL  = ifelse(!is.na(ADT) & !is.na(TRTSDT) & ADT < TRTSDT, "Y", "N"),
      DOMAIN = "MO",
      # Normaliza nome da visita
      AVISIT = toupper(trimws(VISIT)),
      AVISITN = suppressWarnings(as.numeric(VISITNUM))
    ) %>%
    select(USUBJID, ARM, ARMCD, SAFFL, ITTFL,
           AVISIT, AVISITN, ADT, ADY,
           TESTCD = MOTESTCD, TEST = MOTEST,
           any_of(c("LOC" = "MOLOC", "CAT" = "MOCAT",
                    "EVAL" = "MOEVAL", "STRESC" = "MOSTRESC")),
           AVAL, AVALC, ABLFL, DOMAIN)
  message(sprintf("  MO MES rows: %d", nrow(mo_mes)))
}

# ── MI: extrai linhas histológicas ───────────────────────────────────────────
mi_histo <- NULL
if (!is.null(mi) && length(mes_testcds_mi) > 0) {
  mi_histo <- mi %>%
    filter(toupper(trimws(MITESTCD)) %in% mes_testcds_mi) %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ADT    = iso_date(MIDTC),
      ADY    = calc_ady(ADT, TRTSDT),
      AVAL   = suppressWarnings(as.numeric(MIORRES)),
      AVALC  = as.character(MIORRES),
      ABLFL  = ifelse(!is.na(ADT) & !is.na(TRTSDT) & ADT < TRTSDT, "Y", "N"),
      DOMAIN = "MI",
      AVISIT = toupper(trimws(VISIT)),
      AVISITN = suppressWarnings(as.numeric(VISITNUM))
    ) %>%
    select(USUBJID, ARM, ARMCD, SAFFL, ITTFL,
           AVISIT, AVISITN, ADT, ADY,
           TESTCD = MITESTCD, TEST = MITEST,
           any_of(c("LOC" = "MILOC", "CAT" = "MICAT",
                    "STRESC" = "MISTRESC")),
           AVAL, AVALC, ABLFL, DOMAIN)
  message(sprintf("  MI Histology rows: %d", nrow(mi_histo)))
}

# ── QS: extrai scores MES já calculados ──────────────────────────────────────
qs_mes <- NULL
qs_sources <- list(
  list(df = qser, label = "QSER"),
  list(df = qscc, label = "QSCC")
)
qs_mes_list <- map(qs_sources, function(src) {
  df <- src$df; label <- src$label
  if (is.null(df) || nrow(df) == 0) return(NULL)
  testcds <- detect_testcds(df, "QSTESTCD", MES_PATTERNS, label)
  if (length(testcds) == 0) return(NULL)
  df %>%
    filter(toupper(trimws(QSTESTCD)) %in% testcds) %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ADT    = iso_date(QSDTC),
      ADY    = calc_ady(ADT, TRTSDT),
      AVAL   = suppressWarnings(as.numeric(QSORRES)),
      AVALC  = as.character(QSORRES),
      ABLFL  = ifelse(!is.na(ADT) & !is.na(TRTSDT) & ADT < TRTSDT, "Y", "N"),
      DOMAIN = label,
      AVISIT = toupper(trimws(VISIT)),
      AVISITN = suppressWarnings(as.numeric(VISITNUM))
    ) %>%
    select(USUBJID, ARM, ARMCD, SAFFL, ITTFL,
           AVISIT, AVISITN, ADT, ADY,
           TESTCD = QSTESTCD, TEST = QSTEST,
           any_of(c("CAT" = "QSCAT")),
           AVAL, AVALC, ABLFL, DOMAIN)
}) %>% compact()

if (length(qs_mes_list) > 0) qs_mes <- bind_rows(qs_mes_list)

# =============================================================================
# COMBINA TODAS AS FONTES
# =============================================================================
all_scores <- bind_rows(
  mo_mes,
  mi_histo,
  qs_mes
)

if (nrow(all_scores) == 0) {
  stop(paste0(
    "\nNenhum score MES/histológico detectado automaticamente.\n",
    "Verifique o DIAGNÓSTICO impresso acima e ajuste MES_PATTERNS ou HISTO_PATTERNS.\n",
    "Exemplo: adicione o TESTCD real encontrado no diagnóstico:\n",
    "  mes_testcds_mo <- c('SEU_TESTCD_AQUI')\n"
  ))
}

message(sprintf("\n  Total de registros combinados: %d", nrow(all_scores)))

# =============================================================================
# BASELINE POR TESTCD
# =============================================================================
all_scores <- all_scores %>%
  group_by(USUBJID, TESTCD) %>%
  mutate(
    # Baseline = último valor com ABLFL=="Y" antes de TRTSDT
    BASE = {
      bl <- AVAL[ABLFL == "Y" & !is.na(AVAL)]
      dt <- ADT[ABLFL == "Y" & !is.na(AVAL)]
      if (length(bl) > 0) bl[which.max(dt)] else NA_real_
    },
    CHG  = AVAL - BASE,
    PCHG = ifelse(!is.na(BASE) & BASE != 0, (CHG / BASE) * 100, NA_real_)
  ) %>%
  ungroup()

# =============================================================================
# FLAGS CLÍNICOS (baseados em MES overall — MESOVRL ou equivalente)
# =============================================================================
# Detecta qual TESTCD representa o score GLOBAL (não por segmento)
overall_pattern <- c("OVRL", "OVERALL", "TOTAL", "GLOBAL", "MES$", "MAYO$",
                     "MESOVRL", "MAYOSCORE")

overall_testcd <- all_scores %>%
  filter(DOMAIN == "MO") %>%
  pull(TESTCD) %>%
  unique() %>%
  { .[sapply(., function(t) any(grepl(paste(overall_pattern, collapse="|"), t, ignore.case=TRUE)))] }

# Fallback: se não detectar "overall", usa o TESTCD de maior frequência em MO
if (length(overall_testcd) == 0 && !is.null(mo_mes)) {
  overall_testcd <- mo_mes %>%
    count(TESTCD) %>%
    slice_max(n, n = 1) %>%
    pull(TESTCD)
  message(sprintf("  [WARN] TESTCD overall não detectado — usando '%s' como fallback",
                  overall_testcd))
} else {
  message(sprintf("  [OK] TESTCD overall detectado: %s",
                  paste(overall_testcd, collapse = ", ")))
}

# Adiciona flags ao dataset long
all_scores <- all_scores %>%
  mutate(
    # Flags baseados no score overall
    IS_OVERALL = TESTCD %in% overall_testcd,
    # Remissão endoscópica: MES ≤ 1
    MESREMFL   = ifelse(IS_OVERALL & !is.na(AVAL), ifelse(AVAL <= 1, "Y", "N"), NA_character_),
    # Mucosal healing: MES = 0
    MESHEALFL  = ifelse(IS_OVERALL & !is.na(AVAL), ifelse(AVAL == 0, "Y", "N"), NA_character_),
    # Deep remission: MES = 0 ou 1 (sem ulceração friável)
    MESDEEPFL  = ifelse(IS_OVERALL & !is.na(AVAL), ifelse(AVAL <= 1, "Y", "N"), NA_character_),
    # Resposta endoscópica: redução ≥ 1 ponto vs baseline
    MESRESPFL  = ifelse(IS_OVERALL & !is.na(CHG), ifelse(CHG <= -1, "Y", "N"), NA_character_),
    # MES worsening: aumento vs baseline
    MESWORSFL  = ifelse(IS_OVERALL & !is.na(CHG), ifelse(CHG > 0, "Y", "N"), NA_character_)
  ) %>%
  select(-IS_OVERALL)

# =============================================================================
# FORMATO LONG → ADMAYO_LONG
# =============================================================================
message("\n========== GERANDO ADMAYO_LONG ==========")

admayo_long <- all_scores %>%
  arrange(USUBJID, AVISITN, TESTCD) %>%
  select(
    USUBJID, ARM, ARMCD, SAFFL, ITTFL,
    AVISIT, AVISITN, ADT, ADY,
    DOMAIN, TESTCD, TEST,
    any_of(c("LOC", "CAT", "EVAL")),
    AVAL, AVALC, BASE, CHG, PCHG, ABLFL,
    MESREMFL, MESHEALFL, MESDEEPFL, MESRESPFL, MESWORSFL
  )

write.csv(admayo_long,
          file.path(adam_dir, "ADMAYO_LONG.csv"),
          row.names = FALSE, na = "")
message(sprintf("  [SAVE] ADMAYO_LONG.csv  %d obs × %d vars",
                nrow(admayo_long), ncol(admayo_long)))

# =============================================================================
# PIVOT WIDE: USUBJID × VISIT (formato para Foundation Model)
# =============================================================================
message("\n========== GERANDO ADMAYO (wide USUBJID × VISIT) ==========")

# Cria nome de coluna limpo: TESTCD_metric
# ex: MESOVRL_AVAL, MESOVRL_BASE, MESOVRL_CHG, MESOVRL_REMFL
make_wide_cols <- function(df) {
  df %>%
    # Para cada TESTCD: AVAL, BASE, CHG, flags
    select(USUBJID, AVISIT, AVISITN, ADT, ADY,
           TESTCD,
           AVAL, BASE, CHG, PCHG, ABLFL,
           MESREMFL, MESHEALFL, MESDEEPFL, MESRESPFL, MESWORSFL) %>%
    # Garante um registro por USUBJID × VISIT × TESTCD (pega o primeiro se duplicado)
    group_by(USUBJID, AVISIT, AVISITN, ADT, ADY, TESTCD) %>%
    slice(1) %>%
    ungroup() %>%
    pivot_wider(
      id_cols     = c(USUBJID, AVISIT, AVISITN, ADT, ADY),
      names_from  = TESTCD,
      values_from = c(AVAL, BASE, CHG, PCHG, ABLFL,
                      MESREMFL, MESHEALFL, MESDEEPFL, MESRESPFL, MESWORSFL),
      names_glue  = "{TESTCD}_{.value}"
    )
}

# Wide do MO (scores endoscópicos)
wide_mo <- if (!is.null(mo_mes) && nrow(mo_mes) > 0)
  make_wide_cols(mo_mes) else NULL

# Wide do MI (histologia)
wide_mi <- if (!is.null(mi_histo) && nrow(mi_histo) > 0)
  make_wide_cols(
    mi_histo %>% rename(
      MESREMFL = MESREMFL, MESHEALFL = MESHEALFL,
      MESDEEPFL = MESDEEPFL, MESRESPFL = MESRESPFL, MESWORSFL = MESWORSFL
    ) %>%
    mutate(across(any_of(c("MESREMFL","MESHEALFL","MESDEEPFL","MESRESPFL","MESWORSFL")),
                  ~ NA_character_))
  ) else NULL

# Junta MO e MI pelo grain USUBJID × VISIT
admayo_wide_list <- compact(list(wide_mo, wide_mi))

if (length(admayo_wide_list) == 0) {
  stop("Nenhuma tabela wide gerada. Verifique os dados.")
}

admayo_wide <- reduce(admayo_wide_list,
                      full_join,
                      by = c("USUBJID", "AVISIT", "AVISITN", "ADT", "ADY"))

# Junta variáveis de sujeito do ADSL
admayo <- adsl_key %>%
  select(USUBJID, ARM, ARMCD, SAFFL, ITTFL,
         any_of(c("PPROTFL", "TRTSDT", "TRTEDT"))) %>%
  right_join(admayo_wide, by = "USUBJID") %>%
  arrange(USUBJID, AVISITN)

write.csv(admayo,
          file.path(adam_dir, "ADMAYO.csv"),
          row.names = FALSE, na = "")

message(sprintf("  [SAVE] ADMAYO.csv  %d obs × %d vars",
                nrow(admayo), ncol(admayo)))

# =============================================================================
# SUMÁRIO FINAL
# =============================================================================
message("\n========================================")
message("  ADMAYO gerado com sucesso")
message("========================================")
message(sprintf("  Sujeitos únicos   : %d", n_distinct(admayo$USUBJID)))
message(sprintf("  Visitas únicas    : %d", n_distinct(admayo$AVISIT)))
message(sprintf("  Colunas wide      : %d", ncol(admayo)))
message(sprintf("  Fontes usadas     : %s",
                paste(unique(all_scores$DOMAIN), collapse = ", ")))
message(sprintf("  TESTCDs overall   : %s",
                paste(overall_testcd, collapse = ", ")))

# Tabela de cobertura por visita
message("\n  Cobertura por visita (N sujeitos com AVAL overall não-NA):")
if (length(overall_testcd) > 0) {
  overall_aval_col <- paste0(overall_testcd[1], "_AVAL")
  if (overall_aval_col %in% names(admayo)) {
    cov <- admayo %>%
      group_by(AVISIT, AVISITN) %>%
      summarise(
        N_total = n(),
        N_aval  = sum(!is.na(!!sym(overall_aval_col))),
        N_remis = sum(!!sym(paste0(overall_testcd[1], "_MESREMFL")) == "Y",
                      na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(AVISITN)
    print(as.data.frame(cov), row.names = FALSE)
  }
}

message("\n  Arquivos salvos:")
message(sprintf("    %s", file.path(adam_dir, "ADMAYO.csv")))
message(sprintf("    %s", file.path(adam_dir, "ADMAYO_LONG.csv")))
message("========================================\n")

# =============================================================================
# NOTA PARA O UTILIZADOR
# =============================================================================
# Se o diagnóstico acima mostrar TESTCDs que não foram detectados
# automaticamente, adicione-os manualmente assim:
#
#   mes_testcds_mo <- c(mes_testcds_mo, "SEU_TESTCD")
#
# e re-execute a partir do bloco "EXTRAÇÃO E PADRONIZAÇÃO DOS SCORES".
#
# Definições clínicas usadas (ajuste conforme SAP do estudo):
#   MESREMFL  = "Y" se MES overall ≤ 1  (remissão endoscópica)
#   MESHEALFL = "Y" se MES overall = 0  (mucosal healing)
#   MESDEEPFL = "Y" se MES overall ≤ 1  (remissão profunda — combinar com histologia)
#   MESRESPFL = "Y" se CHG overall ≤ -1 (resposta endoscópica)
#   MESWORSFL = "Y" se CHG overall > 0  (piora)
# =============================================================================
