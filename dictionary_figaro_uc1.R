# ============================================================
# ADVANCED ML FEATURE DICTIONARY (CLINICAL-AWARE)
# STUDY: SHP647UC301-FIGARO-UC1
# ============================================================

library(data.table)
library(dplyr)
library(stringr)
library(purrr)

# ===============================
# CONFIG
# ===============================
INPUT_ARD   <- "/mnt/figaro_uc1_ml_dataset_v2.csv"
DATA_PATH   <- "/domino/datasets/local/clinical-trial-data/SHP647UC301-FIGARO-UC1/data"
OUTPUT_DICT <- "/mnt/figaro_uc1_feature_dictionary.csv"

# ---------------------------------------------------------------
# MES GT — Modified Endoscopic Score / Geboes Transform
#
# O time está interessado nos valores MES GT, que no ARD do
# FIGARO UC1 estão distribuídos em dois datasets:
#
#   ADMI (Microbiology / Histology biopsy):
#     MITESTCD = G0  |  MITEST = "Epithelial Damage"
#     Categoria clínica: Geboes score — histologia colónica global
#     Colunas-chave no ARD:
#       AVAL_ADMI_AVAL    → score Geboes na visita (range 0–1.28 obs.)
#       AVAL_ADMI_BASE    → score Geboes no baseline
#       AVAL_ADMI_CHG     → variação absoluta vs baseline
#       AVAL_ADMI_ABLFL   → flag de baseline (Y = visita de referência)
#       AVAL_ADMI_MIORRES → valor original reportado
#
#   ADMO (Morphology / Endoscopy):
#     MOTESTCD = ILER  |  MOTEST = "Exploration results: Ileum"
#     Categoria clínica: exploração endoscópica do íleum
#     Colunas-chave no ARD:
#       AVAL_ADMO_AVAL    → score endoscópico na visita (range 0–8.6 obs.)
#       AVAL_ADMO_BASE    → score endoscópico no baseline
#       AVAL_ADMO_CHG     → variação absoluta vs baseline
#       AVAL_ADMO_ABLFL   → flag de baseline
#       AVAL_ADMO_MOORRES → resultado original da exploração
#
# Ambas as colunas têm ~94% de missing (n=12 sujeitos; biopsia e
# endoscopia só em visitas específicas: SCREENING e WEEK 16).
# ---------------------------------------------------------------

# ===============================
# LOAD ARD
# ===============================
df   <- fread(INPUT_ARD)
cols <- names(df)

# ===============================
# LOAD SOURCE DATASETS (para enriquecimento de labels)
# FIGARO UC1: sem sufixo .sas7bdat, sem PARAMCD/PARAM —
# usa LBTESTCD/LBTEST, VSTESTCD/VSTEST, EGTESTCD/EGTEST,
# MITESTCD/MITEST, MOTESTCD/MOTEST por domínio.
# ===============================
files    <- list.files(DATA_PATH, pattern = "\\.csv$", full.names = TRUE)
datasets <- lapply(files, function(f) {
  tryCatch(fread(f), error = function(e) NULL)
})
names(datasets) <- tolower(tools::file_path_sans_ext(basename(files)))
datasets <- datasets[!sapply(datasets, is.null)]

# ---------------------------------------------------------------
# Lookup de labels por domínio — cada domínio ADaM do FIGARO UC1
# usa um par de colunas código/label diferente (sem PARAMCD/PARAM)
# ---------------------------------------------------------------
domain_testcode_map <- list(
  adlb = c("LBTESTCD", "LBTEST"),
  adeg = c("EGTESTCD", "EGTEST"),
  advs = c("VSTESTCD", "VSTEST"),
  admi = c("MITESTCD", "MITEST"),
  admo = c("MOTESTCD", "MOTEST")
)

# Constrói lookup unificado: dataset + testcd → label
test_lookup <- map_df(names(domain_testcode_map), function(ds_name) {
  if (!(ds_name %in% names(datasets))) return(NULL)
  df_ds  <- datasets[[ds_name]]
  cols_p <- domain_testcode_map[[ds_name]]
  if (!all(cols_p %in% names(df_ds))) return(NULL)
  df_ds %>%
    as.data.frame() %>%
    select(TESTCD = all_of(cols_p[1]), TESTLABEL = all_of(cols_p[2])) %>%
    distinct() %>%
    mutate(dataset = ds_name)
})

# ===============================
# PARSE VARIABLE
# ===============================
parse_var <- function(var_name) {

  if (!str_detect(var_name, "^(AVAL|STC)_")) {
    return(tibble(
      variable     = var_name,
      type         = "ID",
      dataset      = NA_character_,
      original_var = var_name
    ))
  }

  parts <- str_split(var_name, "_", simplify = TRUE)
  tibble(
    variable     = var_name,
    type         = parts[1],
    dataset      = parts[2],
    original_var = paste(parts[-c(1, 2)], collapse = "_")
  )
}

meta <- map_df(cols, parse_var)

# ===============================
# DOMAIN-SPECIFIC DEFINITIONS
# ===============================
define_variable_meaning <- function(type, dataset, var) {

  var_upper <- toupper(var)
  ds_lower  <- tolower(dataset)

  # ---- STANDARD VARIABLES ----
  standard_dict <- list(
    # Grain
    USUBJID  = "Unique subject identifier",
    AVISITN  = "Visit number (protocol visits: 10=SCREENING, 20=BASELINE/WEEK 0, 40=WEEK 4, 50=WEEK 8, 60=WEEK 12, 71-73=WEEK 16, 90=WEEK 32; 101+=pseudo events)",
    AVISIT   = "Visit label",

    # Treatment / demographics
    ARM      = "Treatment arm (label)",
    ARMCD    = "Treatment arm code",
    ACTARM   = "Actual treatment arm (label)",
    ACTARMCD = "Actual treatment arm code",
    AGE      = "Age of subject",
    AGEU     = "Age units",
    SEX      = "Sex of subject",
    RACE     = "Race of subject",
    ETHNIC   = "Ethnicity of subject",
    COUNTRY  = "Country of subject",
    SITEID   = "Study site identifier",
    SUBJID   = "Subject identifier",
    INVNAM   = "Investigator name",

    # Population flags
    SAFFL    = "Safety population flag",
    ITTFL    = "Intent-to-treat population flag",
    PPROTFL  = "Per-protocol population flag",
    ENRLFL   = "Enrolled population flag",
    ANL01FL  = "Analysis flag 01",
    TRTEMFL  = "Treatment-emergent flag",
    ABLFL    = "Baseline record flag",

    # Dates / timing
    ASTDT    = "Analysis start date",
    AENDT    = "Analysis end date",
    ADT      = "Analysis date",
    ADY      = "Analysis relative day",
    ASTDY    = "Analysis start relative day",
    AENDY    = "Analysis end relative day",
    TRTSDT   = "Treatment start date",
    TRTEDT   = "Treatment end date",
    RANDDT   = "Randomisation date",
    TRTDUR   = "Treatment duration",
    RFSTDTC  = "Reference start date/time (ISO)",
    RFENDTC  = "Reference end date/time (ISO)",
    RFXSTDTC = "Reference start date/time of exposure (ISO)",
    RFXENDTC = "Reference end date/time of exposure (ISO)",
    RFICDTC  = "Date/time of informed consent (ISO)",
    BRTHDTC  = "Date of birth (ISO)",
    DTHFL    = "Death flag",
    DTHDTC   = "Date of death (ISO)",

    # Analysis values
    AVAL     = "Analysis value (numeric)",
    AVALC    = "Analysis value (character)",
    AVALU    = "Analysis value units",
    BASE     = "Baseline value",
    CHG      = "Change from baseline",
    PCHG     = "Percent change from baseline",

    # Adverse events
    AESEQ    = "Adverse event sequence number",
    AETERM   = "Adverse event reported term",
    AEDECOD  = "Adverse event MedDRA preferred term",
    AEBODSYS = "Adverse event body system",
    AESOC    = "Adverse event system organ class",
    AESEV    = "Adverse event severity",
    AESER    = "Serious adverse event flag",
    AEREL    = "Adverse event relatedness to study drug",
    AEOUT    = "Adverse event outcome",
    AEACN    = "Action taken with study drug",
    AECONTRT = "Concomitant treatment given for AE flag",
    AESTDTC  = "Adverse event start date/time (ISO)",
    AEENDTC  = "Adverse event end date/time (ISO)",
    AESPID   = "Adverse event sponsor ID",
    AEHLGT   = "Adverse event high level group term",
    AEHLGTCD = "Adverse event high level group term code",
    AEHLT    = "Adverse event high level term",
    AEHLTCD  = "Adverse event high level term code",
    AELLT    = "Adverse event lowest level term",
    AELLTCD  = "Adverse event lowest level term code",
    AEPTCD   = "Adverse event preferred term code",

    # Concomitant medications
    CMSEQ    = "Medication sequence number",
    CMTRT    = "Reported medication name",
    CMDECOD  = "Standardised medication name",
    CMCLAS   = "Medication class",
    CMCLASCD = "Medication class code",
    CMDOSE   = "Medication dose",
    CMDOSU   = "Medication dose units",
    CMDOSFRQ = "Medication dose frequency",
    CMROUTE  = "Medication route of administration",
    CMSTDTC  = "Medication start date/time (ISO)",
    CMENDTC  = "Medication end date/time (ISO)",
    CMINDC   = "Medication indication",
    CMCAT    = "Medication category",

    # Lab
    LBSEQ    = "Lab test sequence number",
    LBTESTCD = "Lab test code",
    LBTEST   = "Lab test name",
    LBCAT    = "Lab test category",
    LBORRES  = "Lab original result",
    LBORRESU = "Lab original result units",
    LBNRIND  = "Lab normal range indicator",
    LBSTNRLO = "Lab standard normal range lower limit",
    LBSTNRHI = "Lab standard normal range upper limit",
    LBDTC    = "Lab date/time (ISO)",
    LBBLFL   = "Lab baseline flag",
    LBFAST   = "Fasting flag",
    LBSPEC   = "Lab specimen type",

    # ECG
    EGSEQ    = "ECG sequence number",
    EGTESTCD = "ECG test code",
    EGTEST   = "ECG test name",
    EGCAT    = "ECG category",
    EGORRES  = "ECG original result",
    EGORRESU = "ECG original result units",
    EGSTAT   = "ECG completion status",
    EGDTC    = "ECG date/time (ISO)",

    # Exposure
    EXSEQ    = "Exposure sequence number",
    EXTRT    = "Exposure treatment name",
    EXDOSE   = "Exposure dose",
    EXDOSU   = "Exposure dose units",
    EXDOSFRM = "Exposure dose form",
    EXDOSFRQ = "Exposure dose frequency",
    EXROUTE  = "Exposure route of administration",
    EXSTDTC  = "Exposure start date/time (ISO)",
    EXENDTC  = "Exposure end date/time (ISO)",
    ADURU    = "Exposure duration units",

    # Vital signs
    VSSEQ    = "Vital sign sequence number",
    VSTESTCD = "Vital sign test code",
    VSTEST   = "Vital sign test name",
    VSORRES  = "Vital sign original result",
    VSORRESU = "Vital sign original result units",
    VSDTC    = "Vital sign date/time (ISO)",

    # *** MES GT — HISTOLOGY (ADMI) ***
    MISEQ    = "Histology/biopsy sequence number",
    MITESTCD = "Histology test code [MES GT: G0 = Geboes Epithelial Damage score]",
    MITEST   = "Histology test name [MES GT: 'Epithelial Damage' = Geboes colonic global score]",
    MICAT    = "Histology test category [MES GT: 'Colonic Global']",
    MISPEC   = "Histology specimen type",
    MILOC    = "Histology specimen location",
    MIORRES  = "Histology original result value",
    MISTRESC = "Histology standardised result (character)",
    MIDTC    = "Histology date/time (ISO)",
    MISTAT   = "Histology completion status",
    MIREASND = "Histology reason not done",

    # *** MES GT — ENDOSCOPY / MORPHOLOGY (ADMO) ***
    MOSEQ    = "Endoscopy/morphology sequence number",
    MOTESTCD = "Endoscopy test code [MES GT: ILER = Ileum exploration result]",
    MOTEST   = "Endoscopy test name [MES GT: 'Exploration results: Ileum']",
    MOCAT    = "Endoscopy test category",
    MOLOC    = "Endoscopy location",
    MOORRES  = "Endoscopy original result",
    MOSTRESC = "Endoscopy standardised result (character)",
    MODTC    = "Endoscopy date/time (ISO)"
  )

  if (var_upper %in% names(standard_dict)) {
    base_meaning <- standard_dict[[var_upper]]
  } else {
    base_meaning <- str_replace_all(var, "_", " ")
  }

  if (type == "AVAL") return(paste(base_meaning, "— longitudinal measurement across visits"))
  if (type == "STC")  return(paste(base_meaning, "— subject-level (static) characteristic"))
  return(base_meaning)
}

# ===============================
# GET DOMAIN TEST LABEL
# ===============================
# Para colunas de domínio (AVAL), tenta recuperar o label do
# test code do dataset de origem (ex: LBTESTCD → LBTEST label)
get_test_label <- function(dataset_prefix, var) {

  if (is.na(dataset_prefix) || is.na(var)) return(NA_character_)

  ds_key <- tolower(dataset_prefix)

  # Mapear variável para o código de teste esperado
  testcd_var <- switch(ds_key,
    adlb = "LBTESTCD", adeg = "EGTESTCD", advs = "VSTESTCD",
    admi = "MITESTCD", admo = "MOTESTCD", NULL
  )

  if (is.null(testcd_var) || var != testcd_var) return(NA_character_)

  ds <- test_lookup %>% filter(dataset == ds_key)
  if (nrow(ds) == 0) return(NA_character_)

  # Retornar todos os valores únicos do test label para este dataset
  labels <- unique(ds$TESTLABEL)
  if (length(labels) == 0) return(NA_character_)
  return(paste(labels, collapse = " | "))
}

# ===============================
# FLAG MES GT COLUMNS
# ===============================
# Identifica colunas directamente relacionadas com o endpoint MES GT
is_mes_gt_column <- function(dataset_prefix, var) {
  ds  <- toupper(dataset_prefix)
  v   <- toupper(var)

  # ADMI: Geboes score (histologia)
  if (ds == "ADMI" && v %in% c("AVAL", "BASE", "CHG", "PCHG", "ABLFL",
                                 "MIORRES", "MISTRESC", "MITESTCD", "MITEST",
                                 "MICAT", "MIDTC", "ADT", "ADY")) return(TRUE)

  # ADMO: endoscopia íleum
  if (ds == "ADMO" && v %in% c("AVAL", "BASE", "CHG", "PCHG", "ABLFL",
                                 "MOORRES", "MOSTRESC", "MOTESTCD", "MOTEST",
                                 "MOCAT", "MODTC", "ADT", "ADY")) return(TRUE)
  return(FALSE)
}

# ===============================
# GENERATE MEANING
# ===============================
meta <- meta %>%
  rowwise() %>%
  mutate(
    test_label = get_test_label(dataset, original_var),
    mes_gt_flag = is_mes_gt_column(dataset, original_var),
    meaning = case_when(
      # Colunas MES GT — enriquecer com contexto clínico explícito
      mes_gt_flag & !is.na(test_label) & type == "AVAL" ~
        paste0("[MES GT] ", define_variable_meaning(type, dataset, original_var),
               " | test: ", test_label),
      mes_gt_flag & type == "AVAL" ~
        paste0("[MES GT] ", define_variable_meaning(type, dataset, original_var)),
      mes_gt_flag & type == "STC" ~
        paste0("[MES GT] ", define_variable_meaning(type, dataset, original_var)),

      # Demais colunas com label de teste disponível
      !is.na(test_label) & type == "AVAL" ~
        paste0(define_variable_meaning(type, dataset, original_var),
               " | test: ", test_label),

      TRUE ~ define_variable_meaning(type, dataset, original_var)
    )
  ) %>%
  ungroup()

# ===============================
# DATA TYPE
# ===============================
detect_data_type <- function(vec) {
  if (is.numeric(vec))   return("numeric")
  if (is.character(vec)) return("categorical")
  if (is.logical(vec))   return("binary")
  return(class(vec)[1])
}

meta$data_type <- map_chr(df, detect_data_type)[meta$variable]

# ===============================
# STATS
# ===============================
get_stats <- function(vec) {
  if (is.numeric(vec)) {
    return(c(
      missing_pct = round(mean(is.na(vec)) * 100, 2),
      unique_n    = n_distinct(vec),
      mean        = round(mean(vec, na.rm = TRUE), 4),
      sd          = round(sd(vec, na.rm = TRUE), 4)
    ))
  } else {
    return(c(
      missing_pct = round(mean(is.na(vec)) * 100, 2),
      unique_n    = n_distinct(vec),
      mean        = NA_real_,
      sd          = NA_real_
    ))
  }
}

stats <- as.data.frame(t(sapply(df, get_stats)))

meta <- bind_cols(meta, stats)

# ===============================
# FINAL DICTIONARY
# ===============================
dictionary <- meta %>%
  select(
    variable,
    dataset,
    type,
    original_var,
    mes_gt_flag,     # ← coluna de destaque para o time
    data_type,
    missing_pct,
    unique_n,
    mean,
    sd,
    meaning
  ) %>%
  # Ordenar: MES GT primeiro, depois por dataset e tipo
  arrange(desc(mes_gt_flag), dataset, type, variable)

write.csv(dictionary, OUTPUT_DICT, row.names = FALSE)

# ===============================
# SUMMARY REPORT
# ===============================
cat("\n✅ Feature dictionary (FIGARO UC1) saved at:", OUTPUT_DICT, "\n")
cat("   Variables documented:", nrow(dictionary), "\n")
cat("   MES GT flagged columns:", sum(dictionary$mes_gt_flag), "\n\n")

cat("=== MES GT COLUMN REFERENCE FOR THE TEAM ===\n\n")

cat("HISTOLOGY — Geboes score (ADMI, MITESTCD=G0)\n")
cat("  Score value (per visit):  AVAL_ADMI_AVAL\n")
cat("  Baseline score:           AVAL_ADMI_BASE\n")
cat("  Change from baseline:     AVAL_ADMI_CHG\n")
cat("  Baseline flag:            AVAL_ADMI_ABLFL\n")
cat("  Original reported value:  AVAL_ADMI_MIORRES\n")
cat("  Test code:                AVAL_ADMI_MITESTCD  (G0 = Epithelial Damage)\n")
cat("  Visit available:          SCREENING (AVISITN=10), WEEK 16/ET PART 2 (AVISITN=72)\n\n")

cat("ENDOSCOPY — Ileum exploration (ADMO, MOTESTCD=ILER)\n")
cat("  Score value (per visit):  AVAL_ADMO_AVAL\n")
cat("  Baseline score:           AVAL_ADMO_BASE\n")
cat("  Change from baseline:     AVAL_ADMO_CHG\n")
cat("  Baseline flag:            AVAL_ADMO_ABLFL\n")
cat("  Original result:          AVAL_ADMO_MOORRES\n")
cat("  Test code:                AVAL_ADMO_MOTESTCD  (ILER = Exploration results: Ileum)\n")
cat("  Visit available:          SCREENING (AVISITN=10), WEEK 16/ET PART 2 (AVISITN=72)\n\n")

cat("NOTE: Both columns have ~94% missing — biopsies/endoscopies only at\n")
cat("      specific protocol visits (Screening and Week 16/ET Part 2).\n")
