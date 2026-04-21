# ============================================================
# ADVANCED ML FEATURE DICTIONARY (CLINICAL-AWARE)
# STUDY: CNTO1275PUC3001-UNIFI-JR / load-2793
# ============================================================

library(data.table)
library(dplyr)
library(stringr)
library(purrr)

# ===============================
# CONFIG
# ===============================
INPUT_ARD   <- "/mnt/unifijr_load2793_ml_dataset_v1.csv"
DATA_PATH   <- "/domino/datasets/local/clinical-trial-data/CNTO1275PUC3001-UNIFI-JR/load-2793/Data/_csv"
OUTPUT_DICT <- "/mnt/unifijr_load2793_feature_dictionary.csv"

# ---------------------------------------------------------------
# Diferenças load-2793 (UNIFI JR) vs load-1899 (UNIFI adulto)
# relevantes para o dicionário:
#
# 1. SEM sufixo .sas7bdat — prefixos de coluna são AVAL_ADMAYO_*,
#    AVAL_ADHIST_*, etc. (vs AVAL_ADMAYO_* do adulto que tinha
#    AVAL_ADMAYO.SAS7BDAT_* antes da limpeza de sufixo)
#
# 2. Datasets presentes (7 vs 33 do adulto):
#    adbdc, adeff, adhist, adlbef, admayo, adpucai, adsl, advpucai
#
# 3. Datasets EXCLUSIVOS do UNIFI JR (pediátrico):
#    adpucai   — PUCAI score pediátrico (ABPAIN, RBLEED, PUCAITS…)
#    advpucai  — PUCAI virtual/diary com duplo período (sufixos 1/2)
#    adhist    — histologia rica (RHI, Geboes, GBTOT, NITOT, 40+ PARAMCDs)
#
# 4. Flags de análise distintos:
#    UNIFI JR usa FASFL/FASRFL/FASCRFL/FASRESFL/SAFRFL/SAFCRFL
#    em vez de SAFW8FL/SAF2FL/COMPLFL/ENTERMFL/ENEXTFL/RAND2FL do adulto
#
# 5. Bug corrigido do load-1899:
#    get_param_meaning: filter(dataset == dataset) comparava a coluna
#    consigo própria (sempre TRUE). Corrigido usando variável local
#    explícita ds_key para a comparação.
# ---------------------------------------------------------------

# ===============================
# LOAD ARD
# ===============================
df   <- fread(INPUT_ARD)
cols <- names(df)

# ===============================
# LOAD SOURCE DATASETS (for PARAM lookup)
# ===============================
files    <- list.files(DATA_PATH, pattern = "\\.csv$", full.names = TRUE)
datasets <- lapply(files, function(f) {
  tryCatch(fread(f), error = function(e) NULL)
})
names(datasets) <- tolower(tools::file_path_sans_ext(basename(files)))
datasets <- datasets[!sapply(datasets, is.null)]

# ===============================
# BUILD PARAM DICTIONARY
# ===============================
# Percorre todos os datasets e extrai pares PARAMCD → PARAM.
# UNIFI JR não tem sufixo .sas7bdat — nomes ficam diretos
# (ex: "admayo", "adhist", "adpucai").
param_lookup <- map_df(names(datasets), function(ds_name) {

  df_ds <- datasets[[ds_name]]

  if (!all(c("PARAMCD", "PARAM") %in% names(df_ds))) return(NULL)

  df_ds %>%
    as.data.frame() %>%
    distinct(PARAMCD, PARAM) %>%
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

  standard_dict <- list(

    # --- Grain ---
    USUBJID   = "Unique subject identifier",
    AVISITN   = "Visit number (numeric)",
    AVISIT    = "Visit label",

    # --- Demographics ---
    AGE       = "Age of subject",
    AGEU      = "Age units",
    AGEGR1    = "Age group 1",
    AGEGR1N   = "Age group 1 (numeric)",
    SEX       = "Sex of subject",
    RACE      = "Race of subject",
    ETHNIC    = "Ethnicity of subject",
    COUNTRY   = "Country of subject",
    REGION1   = "Geographic region 1",
    REGION2   = "Geographic region 2",
    SITEID    = "Study site identifier",
    SUBJID    = "Subject identifier",
    WEIGHTBL  = "Baseline body weight (kg)",
    HEIGHTBL  = "Baseline height (cm)",
    BMIBL     = "Baseline BMI",
    BSABL     = "Baseline body surface area",
    WGTGR1    = "Weight group 1",
    WGTGR1N   = "Weight group 1 (numeric)",
    WGTGR2    = "Weight group 2",
    WGTGR2N   = "Weight group 2 (numeric)",

    # --- Treatment ---
    ARM       = "Randomised treatment arm (label)",
    ARMCD     = "Randomised treatment arm code",
    ACTARM    = "Actual treatment arm (label)",
    ACTARMCD  = "Actual treatment arm code",
    TRT01A    = "Actual treatment period 1",
    TRT01AN   = "Actual treatment period 1 (numeric)",
    TRT01P    = "Planned treatment period 1",
    TRT01PN   = "Planned treatment period 1 (numeric)",
    TRT02A    = "Actual treatment period 2",
    TRT02AN   = "Actual treatment period 2 (numeric)",
    TRT02P    = "Planned treatment period 2",
    TRT02PN   = "Planned treatment period 2 (numeric)",
    TR02AG1   = "Actual treatment period 2 subgroup 1",
    TR02AG1N  = "Actual treatment period 2 subgroup 1 (numeric)",
    TR02PG1   = "Planned treatment period 2 subgroup 1",
    TR02PG1N  = "Planned treatment period 2 subgroup 1 (numeric)",
    TRT03A    = "Actual treatment period 3",
    TRT03AN   = "Actual treatment period 3 (numeric)",
    TRT03P    = "Planned treatment period 3",
    TRT03PN   = "Planned treatment period 3 (numeric)",

    # --- Analysis flags (UNIFI JR specific) ---
    RANDFL    = "Randomised population flag",
    ENRLFL    = "Enrolled population flag",
    SAFFL     = "Safety population flag",
    SAFRFL    = "Safety randomised subpopulation flag",
    SAFCRFL   = "Safety crossover randomised subpopulation flag",
    FASFL     = "Full analysis set flag",
    FASRFL    = "Full analysis set randomised subpopulation flag",
    FASCRFL   = "Full analysis set crossover randomised flag",
    FASRESFL  = "Full analysis set re-screened flag",
    PKFL      = "PK population flag",
    PKRFL     = "PK randomised subpopulation flag",
    PKCRFL    = "PK crossover randomised subpopulation flag",
    PKRESFL   = "PK re-screened flag",
    IASFL     = "Immunogenicity analysis set flag",
    IASRFL    = "Immunogenicity analysis set randomised flag",
    IASCRFL   = "Immunogenicity analysis set crossover flag",
    IASRESFL  = "Immunogenicity analysis set re-screened flag",
    RESCRNFL  = "Re-screening flag",
    CORBFL    = "Corticosteroid baseline flag",
    IMBFL     = "Immunomodulator baseline flag",
    BIO5HA    = "Biologic 5-hour baseline flag",
    BI5HACRF  = "Biologic 5-hour baseline crossover flag",
    CRSPM0    = "Clinical response at month 0",
    DCRSFL    = "Discontinuation clinical response flag",
    RSPIWRS   = "Response at induction week (randomised)",
    DAFL      = "Dose adaptation flag",
    SUBFL     = "Sub-study flag",
    UNBLFL    = "Unblinded flag",
    DCTFL     = "Discontinuation from treatment flag",
    DSTFL     = "Discontinuation from study flag",

    # --- Analysis values ---
    AVAL      = "Analysis value (numeric)",
    AVALC     = "Analysis value (character)",
    AVALCAT1  = "Analysis value category 1",
    BASE      = "Baseline value",
    BASE2     = "Secondary baseline value",
    CHG       = "Change from baseline",
    CHG2      = "Change from secondary baseline",
    PCHG      = "Percent change from baseline",
    DTYPE     = "Derivation type",
    ABLFL     = "Baseline record flag",
    ABLFL2    = "Secondary baseline record flag",
    APOBLFL   = "Post-baseline record flag",
    ANL01FL   = "Analysis flag 01",
    ANL02FL   = "Analysis flag 02",
    ANL03FL   = "Analysis flag 03",
    ANL04FL   = "Analysis flag 04",
    ANL05FL   = "Analysis flag 05",
    ANL06FL   = "Analysis flag 06",
    ICE       = "ICE (intercurrent event) flag",

    # --- Dates ---
    ADT       = "Analysis date",
    ADTM      = "Analysis date/time",
    ADTF      = "Analysis date imputation flag",
    ADY       = "Analysis relative day",
    ASTDT     = "Analysis start date",
    ASTDTC    = "Analysis start date/time (ISO)",
    ARFSTDT   = "Reference start date (induction)",
    ARFSTDTM  = "Reference start date/time (induction)",
    TRTSDT    = "Treatment start date",
    TRTSDTM   = "Treatment start date/time",
    TRTEDT    = "Treatment end date",
    TRTEDTM   = "Treatment end date/time",
    TRTEDY    = "Treatment end relative day",
    RFICDT    = "Informed consent date",
    RANDDT    = "Randomisation date",
    RANDDTM   = "Randomisation date/time",
    EOSDT     = "End of study date",
    EOSDTM    = "End of study date/time",
    EOTDT     = "End of treatment date",
    EOTDTM    = "End of treatment date/time",
    DTHDT     = "Death date",
    DTHDY     = "Death relative day",
    DCSDT     = "Discontinuation from study date",
    DSTDT     = "Discontinuation from study date",
    DSTDY     = "Discontinuation from study relative day",
    DCTDT     = "Discontinuation from treatment date",
    DCTDY     = "Discontinuation from treatment relative day",
    WK40DT    = "Week 40 date",
    WK44DT    = "Week 44 date",

    # --- Study period dates (AP01-03) ---
    AP01SDT   = "Analysis period 1 start date",
    AP01SDTM  = "Analysis period 1 start date/time",
    AP01EDT   = "Analysis period 1 end date",
    AP01EDTM  = "Analysis period 1 end date/time",
    AP02SDT   = "Analysis period 2 start date",
    AP02SDTM  = "Analysis period 2 start date/time",
    AP02EDT   = "Analysis period 2 end date",
    AP02EDTM  = "Analysis period 2 end date/time",
    AP03SDT   = "Analysis period 3 start date",
    AP03SDTM  = "Analysis period 3 start date/time",
    AP03EDT   = "Analysis period 3 end date",
    AP03EDTM  = "Analysis period 3 end date/time",

    # --- PUCAI (Pediatric UC Activity Index) — EXCLUSIVE pUC ---
    ABPAIN    = "PUCAI — abdominal pain score",
    RBLEED    = "PUCAI — rectal bleeding score",
    STLCONS   = "PUCAI — stool consistency score",
    NUMSTL    = "PUCAI — number of stools per day",
    NOCTSTL   = "PUCAI — nocturnal stools score",
    ACTLEV    = "PUCAI — activity level score",
    PUCAITS   = "PUCAI — total score (0–85)",
    PUCREMI1  = "PUCAI — remission flag period 1",
    PUCREMI2  = "PUCAI — remission flag period 2",
    PUCLCHG1  = "PUCAI — last change from baseline period 1",
    PUCLCHG2  = "PUCAI — last change from baseline period 2",
    CORFRE90  = "Corticosteroid-free at 90 days flag",
    SSTRTSDT  = "Sub-study treatment start date",
    SSVISIT   = "Sub-study visit label",
    SSVISITN  = "Sub-study visit number",
    SSRFL     = "Sub-study randomised flag",
    SSICDT    = "Sub-study informed consent date",
    SSTRTVIS  = "Sub-study treatment visit",
    SSTRVISN  = "Sub-study treatment visit number",
    MISOBSFL  = "Missing observation flag",
    ICEDSDT   = "ICE discontinuation start date",

    # --- Mayo Score (UC) ---
    SFSCORE   = "Mayo stool frequency subscore",
    RBSCORE   = "Mayo rectal bleeding subscore",
    ENSCORE   = "Mayo endoscopy subscore",
    PGSCORE   = "Mayo physician global assessment subscore",
    MAYO      = "Mayo score total",
    MMAYO     = "Modified Mayo score",
    PMAYO     = "Partial Mayo score",
    ABSSTOOL  = "Absolute stool count",
    ENFSCORE  = "Endoscopy friability subscore",
    ENCSCORE  = "Endoscopy colour subscore",
    ENLSCORE  = "Endoscopy loss of vascular pattern subscore",
    NUMMAYO   = "Number of Mayo diary days",

    # --- Histology (ADHIST — rich in UNIFI JR) ---
    ARCHCHG   = "Histology — architectural change score",
    CHINFIN   = "Histology — chronic inflammatory infiltrate",
    EOS       = "Histology — eosinophils in lamina propria",
    NEUTLP    = "Histology — neutrophils in lamina propria",
    NEUTEP    = "Histology — neutrophils in epithelium",
    CRYPTDS   = "Histology — crypt destruction score",
    EROULC    = "Histology — erosion or ulceration",
    ULCER     = "Histology — ulcer score",
    AINFCIN   = "Histology — acute inflammatory cells in crypt",
    OVLGRD    = "Histology — overall grade",
    RHI       = "Histology — Robarts Histopathology Index (RHI)",
    GBTOT     = "Histology — Geboes score total",
    RHITOT    = "Histology — RHI total",
    NITOT     = "Histology — Nancy Index total",
    GBTOTP    = "Histology — Geboes score total (post-baseline)",
    RHITOTP   = "Histology — RHI total (post-baseline)",
    NITOTP    = "Histology — Nancy Index total (post-baseline)",
    HEMIGS    = "Histology — histological remission Geboes score",
    HEMIRHI   = "Histology — histological remission RHI",
    HADISGS   = "Histology — histological disease Geboes score",
    HADISRHI  = "Histology — histological disease RHI",
    HADISNI   = "Histology — histological disease Nancy Index",
    CHINFINI  = "Histology — chronic inflammatory cells (induction)",
    HISTFLG1  = "Histology analysis flag 1",
    HISTFLG2  = "Histology analysis flag 2",
    SRCDOM    = "Source domain",
    SRCVAR    = "Source variable",
    VISIT     = "Visit label (source)",
    VISITNUM  = "Visit number (source)",
    VISITD    = "Visit description",
    VISNUMD   = "Visit number description",
    QSDTC     = "Questionnaire date/time (ISO)",
    MODTC     = "Morphology date/time (ISO)",

    # --- Lab efficacy (ADLBEF) ---
    CALPRO    = "Calprotectin",
    CRPABFL   = "CRP analysis baseline flag",
    CALABFL   = "Calprotectin analysis baseline flag",
    LACABFL   = "Lactoferrin analysis baseline flag",
    CRPICE    = "CRP ICE flag",
    CALICE    = "Calprotectin ICE flag",
    LTFICE    = "Lactoferrin ICE flag",
    CALNOR    = "Calprotectin normalisation",
    CRPNOR    = "CRP normalisation",
    LTFNOR    = "Lactoferrin normalisation",
    LTF       = "Lactoferrin",

    # --- Baseline disease characteristics (ADBDC) ---
    AGEDIAG   = "Age at diagnosis",
    DISDURN   = "Disease duration",
    EXTDIS    = "Extent of disease",
    CRPBL     = "CRP at baseline",
    CALPRBL   = "Calprotectin at baseline",
    LTFBL     = "Lactoferrin at baseline",
    MAYOBL    = "Mayo score at baseline",
    MMAYOBL   = "Modified Mayo score at baseline",
    ENDOBL    = "Endoscopy subscore at baseline",
    PUCAIBL   = "PUCAI at baseline",
    IMPACTBL  = "IMPACT at baseline",
    ALCOHOL   = "Alcohol use flag",
    TOBACCO   = "Tobacco use flag",
    OPIOID    = "Opioid use flag",
    PARAMN    = "Parameter numeric code",
    ASTDTC    = "Analysis start date/time (ISO)",

    # --- Efficacy (ADEFF) ---
    CRMI8P1   = "Clinical remission induction 8 weeks period 1",
    SREM      = "Sustained remission",
    ENIMP     = "Endoscopic improvement",
    CORTFCRM  = "Corticosteroid-free remission",
    REMS1P    = "Remission sub-period 1",
    CRSP      = "Clinical response",
    CRSPDA    = "Clinical response dose adaptation",
    CRSPA     = "Clinical response additional",
    SREMA     = "Sustained remission additional",
    SSVISIT   = "Sub-study visit label",
    SSVISITN  = "Sub-study visit number (numeric)",
    ARFSTDT   = "Analysis reference start date (induction)",
    ARFSTDTM  = "Analysis reference start date/time",
    CORBFL    = "Corticosteroid baseline flag",
    CRSPM0    = "Clinical response month 0",
    DCRSFL    = "Discontinuation response flag",
    RSPIWRS   = "Response at induction week (randomised)",
    CRSPM0D   = "Clinical response month 0 dose adaptation",
    DDCRSFL   = "Dose adaptation discontinuation flag",
    DAFL      = "Dose adaptation flag",
    SUBFL     = "Sub-study flag",

    # --- Discontinuation ---
    DCTREAS   = "Reason for treatment discontinuation",
    DCTREASP  = "Reason for treatment discontinuation (verbatim)",
    DCTSSFL   = "Sub-study discontinuation flag",
    DCTSSREA  = "Sub-study discontinuation reason",
    DCTSSRSP  = "Sub-study discontinuation reason (verbatim)",
    DCSREAS   = "Reason for study discontinuation",
    DCSREASP  = "Reason for study discontinuation (verbatim)",

    # --- Unblinding ---
    UNBLFL    = "Unblinded flag",
    UNBLTDT   = "Unblinding date",
    UNBLDY    = "Unblinding relative day",
    UNBREAS   = "Unblinding reason",

    # --- Other baseline ---
    B5ASA     = "5-ASA use at baseline",
    BIOFAIL   = "Biologic failure flag",
    BMIBLG1   = "BMI group 1 at baseline",
    BMIBLG1N  = "BMI group 1 at baseline (numeric)",
    RANUM     = "Randomisation number",
    LTSTDT    = "Last treatment start date",
    LTENDT    = "Last treatment end date",
    LTVISIT   = "Last treatment visit",
    LTSTDY    = "Last treatment start relative day",
    LTENDY    = "Last treatment end relative day",
    LDOSE     = "Last dose",
    LDOSEU    = "Last dose units",
    LSVISIT   = "Last study visit",
    LOR01FL   = "Loss of response 01 flag",
    LOR02FL   = "Loss of response 02 flag",
    LOWEXFL   = "Low exposure flag",
    DOSADDT   = "Dose adaptation date",
    DOSATRDT  = "Dose adaptation re-treatment date",
    DOSATVIS  = "Dose adaptation visit",
    RANDDTM   = "Randomisation date/time",
    UNBLTDT   = "Unblinding date",
    TR01SDT   = "Treatment period 1 start date",
    TR01SDTM  = "Treatment period 1 start date/time",
    TR01EDT   = "Treatment period 1 end date",
    TR01EDTM  = "Treatment period 1 end date/time",
    TR02SDT   = "Treatment period 2 start date",
    TR02SDTM  = "Treatment period 2 start date/time",
    TR02EDT   = "Treatment period 2 end date",
    TR03SDT   = "Treatment period 3 start date",
    TR03SDTM  = "Treatment period 3 start date/time",
    TR03EDT   = "Treatment period 3 end date",
    TR03EDTM  = "Treatment period 3 end date/time"
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
# ENRICH WITH PARAM
# ===============================
# Correcção do bug do load-1899: filter(dataset == dataset) comparava
# a coluna com ela própria. Aqui usamos ds_key como variável local.
get_param_meaning <- function(dataset_prefix, var) {

  if (is.na(dataset_prefix) || is.na(var)) return(NA_character_)

  # UNIFI JR sem sufixo — normalizar para lowercase directo
  ds_key <- tolower(dataset_prefix)

  ds <- param_lookup %>% filter(dataset == ds_key)

  if (nrow(ds) == 0) return(NA_character_)

  match <- ds %>% filter(PARAMCD == var)

  if (nrow(match) > 0) return(match$PARAM[1])

  return(NA_character_)
}

# ===============================
# GENERATE MEANING
# ===============================
meta <- meta %>%
  rowwise() %>%
  mutate(
    param_desc = get_param_meaning(dataset, original_var),
    meaning = case_when(
      !is.na(param_desc) & type == "AVAL" ~
        paste(param_desc, "— longitudinal clinical measurement"),

      !is.na(param_desc) & type == "STC" ~
        paste(param_desc, "— baseline/static clinical characteristic"),

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
    data_type,
    missing_pct,
    unique_n,
    mean,
    sd,
    meaning
  ) %>%
  arrange(dataset, type, variable)

write.csv(dictionary, OUTPUT_DICT, row.names = FALSE)

cat("\n✅ Feature dictionary (UNIFI JR load-2793) saved at:", OUTPUT_DICT, "\n")
cat("   Variables documented:", nrow(dictionary), "\n")
