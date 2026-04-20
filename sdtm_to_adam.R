# =============================================================================
# SDTM → ADaM Conversion Script
# Estudo: SHP647UC301 / FIGARO-UC1
# Input:  CSV (SDTM) — 83 arquivos (40 domínios + 43 SUPPs)
# Output: CSV (ADaM)  — ADSL, ADAE, ADLB, ADVS, ADEX, ADCM, ADEG,
#                        ADPR, ADHO, ADMI, ADMO, ADQS
# Autor:  gerado automaticamente
# =============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(purrr)

# =============================================================================
# 0. CONFIGURAÇÃO DE PATHS
# =============================================================================

# !! AJUSTE ESTES DOIS PATHS PARA O SEU AMBIENTE !!
sdtm_dir <- "/domino/datasets/local/clinical-trial-data/SHP647UC301-FIGARO-UC1/data/sdtm_csv"
adam_dir <- "/domino/datasets/local/clinical-trial-data/SHP647UC301-FIGARO-UC1/data/adam_csv"

dir.create(adam_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. FUNÇÕES UTILITÁRIAS
# =============================================================================

# Lê um CSV SDTM com colunas em maiúsculas; retorna NULL se não existir
read_sdtm <- function(domain) {
  path <- file.path(sdtm_dir, paste0(tolower(domain), ".csv"))
  if (!file.exists(path)) {
    message(sprintf("  [SKIP] %s.csv não encontrado", domain))
    return(NULL)
  }
  df <- read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  names(df) <- toupper(names(df))
  message(sprintf("  [READ] %-12s  %d obs × %d vars", paste0(domain, ".csv"), nrow(df), ncol(df)))
  df
}

# Lê e pivota um SUPP; retorna NULL se não existir ou estiver vazio
read_supp <- function(domain) {
  df <- read_sdtm(paste0("supp", tolower(domain)))
  if (is.null(df) || nrow(df) == 0) return(NULL)
  # Garante colunas obrigatórias
  needed <- c("USUBJID", "IDVAR", "IDVARVAL", "QNAM", "QVAL")
  if (!all(needed %in% names(df))) return(NULL)
  df %>%
    select(USUBJID, IDVAR, IDVARVAL, QNAM, QVAL) %>%
    pivot_wider(names_from = QNAM, values_from = QVAL,
                values_fn = ~ .x[1])
}

# Merge SUPP pivotado de volta ao domínio pai
merge_supp <- function(parent, supp_wide, seq_var) {
  if (is.null(supp_wide)) return(parent)
  # Remove colunas duplicadas que já existam no parent
  new_cols <- setdiff(names(supp_wide), c("USUBJID", "IDVAR", "IDVARVAL"))
  supp_clean <- supp_wide %>%
    select(USUBJID, IDVARVAL, all_of(new_cols)) %>%
    rename(!!seq_var := IDVARVAL) %>%
    mutate(across(everything(), as.character))
  parent %>%
    mutate(across(all_of(seq_var), as.character)) %>%
    left_join(supp_clean, by = c("USUBJID", seq_var))
}

# Converte ISO8601 date string para Date (aceita YYYY-MM-DD ou YYYY-MM-DDThh:mm)
iso_date <- function(x) as.Date(substr(as.character(x), 1, 10))

# Calcula dia relativo ao estudo (ADY): positivo = on/após TRTSDT; negativo = antes
calc_ady <- function(adt, trtsdt) {
  d <- as.numeric(adt - trtsdt)
  ifelse(is.na(d), NA_real_, ifelse(d >= 0, d + 1, d))
}

# Salva ADaM como CSV e reporta
write_adam <- function(df, name) {
  path <- file.path(adam_dir, paste0(toupper(name), ".csv"))
  write.csv(df, path, row.names = FALSE, na = "")
  message(sprintf("  [SAVE] %-10s  %d obs × %d vars → %s",
                  toupper(name), nrow(df), ncol(df), basename(path)))
  invisible(df)
}

# =============================================================================
# 2. LEITURA DOS DOMÍNIOS SDTM
# =============================================================================
message("\n========== LENDO DOMÍNIOS SDTM ==========")

dm     <- read_sdtm("dm")
dmmp   <- read_sdtm("dmmp")
ds     <- read_sdtm("ds")
ex     <- read_sdtm("ex")
ae     <- read_sdtm("ae")
cm     <- read_sdtm("cm")
lb     <- read_sdtm("lb")
vs     <- read_sdtm("vs")
eg     <- read_sdtm("eg")
pr     <- read_sdtm("pr")
ho     <- read_sdtm("ho")
mi     <- read_sdtm("mi")
mo     <- read_sdtm("mo")
mh     <- read_sdtm("mh")
ec     <- read_sdtm("ec")
su     <- read_sdtm("su")
sv     <- read_sdtm("sv")
se     <- read_sdtm("se")
ss     <- read_sdtm("ss")
ie     <- read_sdtm("ie")
dv     <- read_sdtm("dv")
pc     <- read_sdtm("pc")
pe     <- read_sdtm("pe")
rp     <- read_sdtm("rp")
is_df  <- read_sdtm("is")
xn     <- read_sdtm("xn")
facd   <- read_sdtm("facd")
faae   <- read_sdtm("faae")

# Domínios QS (questionários)
qsai   <- read_sdtm("qsai")   # IBDQ / Adequacy of Information
qsca   <- read_sdtm("qsca")   # CAIS
qscc   <- read_sdtm("qscc")   # CCIS
qscd   <- read_sdtm("qscd")   # CDAI / CD scores
qscr   <- read_sdtm("qscr")   # CRP-related QS
qseq   <- read_sdtm("qseq")   # EQ-5D
qser   <- read_sdtm("qser")   # ERIQ / Endoscopy
qsib   <- read_sdtm("qsib")   # IBDQ
qspg   <- read_sdtm("qspg")   # PGA
qss2   <- read_sdtm("qss2")   # SFS2 / Stool Frequency
qsta   <- read_sdtm("qsta")   # (vazio na checklist)
qsts   <- read_sdtm("qsts")   # TSQM
qswp   <- read_sdtm("qswp")   # WPAI

# =============================================================================
# 3. LEITURA E MERGE DOS SUPPs
# =============================================================================
message("\n========== LENDO E PIVOTANDO SUPPs ==========")

supp_dm    <- read_supp("dm")
supp_dmmp  <- read_supp("dmmp")
supp_ds    <- read_supp("ds")
supp_ex    <- read_supp("ex")
supp_ae    <- read_supp("ae")
supp_cm    <- read_supp("cm")
supp_lb    <- read_supp("lb")
supp_vs    <- read_supp("vs")
supp_eg    <- read_supp("eg")
supp_pr    <- read_supp("pr")
supp_ho    <- read_supp("ho")
supp_mi    <- read_supp("mi")
supp_mo    <- read_supp("mo")
supp_mh    <- read_supp("mh")
supp_ec    <- read_supp("ec")
supp_is    <- read_supp("is")
supp_se    <- read_supp("se")
supp_ss    <- read_supp("ss")
supp_su    <- read_supp("su")
supp_sv    <- read_supp("sv")
supp_rp    <- read_supp("rp")
supp_xn    <- read_supp("xn")
supp_facd  <- read_supp("facd")
supp_pe    <- read_supp("pe")
supp_pc    <- read_supp("pc")

# SUPP QS individuais
supp_qsai  <- read_supp("qsai")
supp_qsca  <- read_supp("qsca")
supp_qscc  <- read_supp("qscc")
supp_qscd  <- read_supp("qscd")
supp_qscr  <- read_supp("qscr")
supp_qseq  <- read_supp("qseq")
supp_qser  <- read_supp("qser")
supp_qsib  <- read_supp("qsib")
supp_qspg  <- read_supp("qspg")
supp_qss2  <- read_supp("qss2")
supp_qsts  <- read_supp("qsts")
supp_qswp  <- read_supp("qswp")

# Merge SUPPs nos respectivos domínios
if (!is.null(dm) && !is.null(supp_dm))
  dm   <- merge_supp(dm,   supp_dm,   "DMSEQ")
if (!is.null(ex) && !is.null(supp_ex))
  ex   <- merge_supp(ex,   supp_ex,   "EXSEQ")
if (!is.null(ae) && !is.null(supp_ae))
  ae   <- merge_supp(ae,   supp_ae,   "AESEQ")
if (!is.null(cm) && !is.null(supp_cm))
  cm   <- merge_supp(cm,   supp_cm,   "CMSEQ")
if (!is.null(lb) && !is.null(supp_lb))
  lb   <- merge_supp(lb,   supp_lb,   "LBSEQ")
if (!is.null(vs) && !is.null(supp_vs))
  vs   <- merge_supp(vs,   supp_vs,   "VSSEQ")
if (!is.null(eg) && !is.null(supp_eg))
  eg   <- merge_supp(eg,   supp_eg,   "EGSEQ")
if (!is.null(pr) && !is.null(supp_pr))
  pr   <- merge_supp(pr,   supp_pr,   "PRSEQ")
if (!is.null(ho) && !is.null(supp_ho))
  ho   <- merge_supp(ho,   supp_ho,   "HOSEQ")
if (!is.null(mi) && !is.null(supp_mi))
  mi   <- merge_supp(mi,   supp_mi,   "MISEQ")
if (!is.null(mo) && !is.null(supp_mo))
  mo   <- merge_supp(mo,   supp_mo,   "MOSEQ")
if (!is.null(ds) && !is.null(supp_ds))
  ds   <- merge_supp(ds,   supp_ds,   "DSSEQ")

# =============================================================================
# 4. ADSL – Subject-Level Analysis Dataset
# =============================================================================
message("\n========== ADSL ==========")

if (!is.null(dm)) {

  adsl <- dm %>%
    select(
      STUDYID, USUBJID, SUBJID, SITEID,
      AGE, AGEU, SEX, RACE, ETHNIC, COUNTRY,
      ARMCD, ARM, ACTARMCD, ACTARM,
      RFSTDTC, RFENDTC,
      any_of(c("RFXSTDTC", "RFXENDTC", "RFICDTC",
                "DTHFL", "DTHDTC",
                "INVID", "INVNAM", "BRTHDTC"))
    ) %>%
    mutate(
      TRTSDT  = iso_date(coalesce(
        if ("RFXSTDTC" %in% names(.)) RFXSTDTC else NA_character_,
        RFSTDTC)),
      TRTEDT  = iso_date(coalesce(
        if ("RFXENDTC" %in% names(.)) RFXENDTC else NA_character_,
        RFENDTC)),
      RANDDT  = iso_date(RFSTDTC),
      TRTDUR  = as.numeric(TRTEDT - TRTSDT) + 1,
      # Populações — refinadas logo abaixo
      ITTFL   = "Y",
      SAFFL   = "N",
      PPROTFL = "N",
      ENRLFL  = "Y"
    )

  # SAFFL: pelo menos uma administração registrada em EX
  if (!is.null(ex) && nrow(ex) > 0) {
    expostos <- ex %>%
      filter(!is.na(EXSTDTC) & trimws(EXSTDTC) != "") %>%
      distinct(USUBJID) %>%
      mutate(SAFFL = "Y")
    adsl <- adsl %>%
      rows_update(expostos, by = "USUBJID", unmatched = "ignore")
  }

  # PPROTFL: completou o estudo
  if (!is.null(ds) && nrow(ds) > 0) {
    completers <- ds %>%
      filter(toupper(trimws(DSDECOD)) == "COMPLETED") %>%
      distinct(USUBJID) %>%
      mutate(PPROTFL = "Y")
    adsl <- adsl %>%
      rows_update(completers, by = "USUBJID", unmatched = "ignore")
  }

  # Adiciona informações de DMMP (dados médicos adicionais) se disponível
  if (!is.null(dmmp) && nrow(dmmp) > 0) {
    dmmp_wide <- dmmp %>%
      select(USUBJID, any_of(c("DMMPTESTCD", "DMMPTEST", "DMMPVAL", "DMMPORRES"))) %>%
      distinct()
    adsl <- adsl %>%
      left_join(dmmp_wide %>% distinct(USUBJID), by = "USUBJID")
  }

  write_adam(adsl, "adsl")
} else {
  message("  [ERRO] DM não disponível — ADSL não gerado. Todos os outros ADaMs serão afetados.")
  adsl <- data.frame()
}

# Atalho para join de variáveis de ADSL (frequentemente reutilizadas)
adsl_key <- if (nrow(adsl) > 0) {
  adsl %>% select(USUBJID, TRTSDT, TRTEDT, ARM, ARMCD, SAFFL, ITTFL, PPROTFL)
} else NULL

# =============================================================================
# 5. ADAE – Adverse Events
# =============================================================================
message("\n========== ADAE ==========")

if (!is.null(ae) && nrow(ae) > 0 && !is.null(adsl_key)) {

  adae <- ae %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ASTDT   = iso_date(AESTDTC),
      AENDT   = iso_date(AEENDTC),
      ASTDY   = calc_ady(ASTDT, TRTSDT),
      AENDY   = calc_ady(AENDT, TRTSDT),
      # Treatment-Emergent: início ≥ TRTSDT
      TRTEMFL = ifelse(!is.na(ASTDT) & !is.na(TRTSDT) & ASTDT >= TRTSDT, "Y", "N"),
      AESER   = toupper(trimws(AESER)),
      ANL01FL = "Y"
    ) %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      AESEQ, AETERM, AEDECOD, AEBODSYS, AESOC,
      AESEV, AESER, AEREL, AEOUT, AEACN,
      AESTDTC, AEENDTC, ASTDT, AENDT, ASTDY, AENDY,
      TRTEMFL, ANL01FL, SAFFL,
      any_of(c("AETOXGR", "AECONTRT", "AEONGO",
                # colunas pivotadas do SUPPAE
                "AERELNST", "AESPID", "AEHLGT", "AEHLGTCD",
                "AELLT", "AELLTCD", "AEPTCD", "AESOCD",
                "AEHLT", "AEHLTCD"))
    )

  write_adam(adae, "adae")
}

# =============================================================================
# 6. ADLB – Laboratory Data
# =============================================================================
message("\n========== ADLB ==========")

if (!is.null(lb) && nrow(lb) > 0 && !is.null(adsl_key)) {

  adlb <- lb %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ADT   = iso_date(LBDTC),
      ADY   = calc_ady(ADT, TRTSDT),
      AVAL  = suppressWarnings(as.numeric(LBORRES)),
      AVALU = LBORRESU,
      ABLFL = ifelse(!is.na(ADT) & !is.na(TRTSDT) & ADT < TRTSDT, "Y", "N")
    ) %>%
    group_by(USUBJID, LBTESTCD) %>%
    mutate(
      # Baseline = último valor pré-tratamento
      BASE = {
        bl_vals <- AVAL[ABLFL == "Y" & !is.na(AVAL)]
        bl_dts  <- ADT[ABLFL == "Y" & !is.na(AVAL)]
        if (length(bl_vals) > 0) bl_vals[which.max(bl_dts)] else NA_real_
      },
      CHG  = AVAL - BASE,
      PCHG = ifelse(!is.na(BASE) & BASE != 0, (CHG / BASE) * 100, NA_real_)
    ) %>%
    ungroup() %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      LBSEQ, LBTESTCD, LBTEST, LBCAT,
      AVAL, AVALU, BASE, CHG, PCHG, ABLFL,
      LBORRES, LBORRESU, LBNRIND, LBSTNRLO, LBSTNRHI,
      LBDTC, ADT, ADY,
      VISITNUM, VISIT, SAFFL,
      any_of(c("LBBLFL", "LBFAST", "LBSPEC", "LBMETHOD"))
    )

  write_adam(adlb, "adlb")
}

# =============================================================================
# 7. ADVS – Vital Signs
# =============================================================================
message("\n========== ADVS ==========")

if (!is.null(vs) && nrow(vs) > 0 && !is.null(adsl_key)) {

  advs <- vs %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ADT   = iso_date(VSDTC),
      ADY   = calc_ady(ADT, TRTSDT),
      AVAL  = suppressWarnings(as.numeric(VSORRES)),
      AVALU = VSORRESU,
      ABLFL = ifelse(!is.na(ADT) & !is.na(TRTSDT) & ADT < TRTSDT, "Y", "N")
    ) %>%
    group_by(USUBJID, VSTESTCD) %>%
    mutate(
      BASE = {
        bl_vals <- AVAL[ABLFL == "Y" & !is.na(AVAL)]
        bl_dts  <- ADT[ABLFL == "Y" & !is.na(AVAL)]
        if (length(bl_vals) > 0) bl_vals[which.max(bl_dts)] else NA_real_
      },
      CHG  = AVAL - BASE,
      PCHG = ifelse(!is.na(BASE) & BASE != 0, (CHG / BASE) * 100, NA_real_)
    ) %>%
    ungroup() %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      VSSEQ, VSTESTCD, VSTEST, any_of("VSCAT"),
      AVAL, AVALU, BASE, CHG, PCHG, ABLFL,
      VSORRES, VSORRESU, any_of("VSNRIND"),
      VSDTC, ADT, ADY, VISITNUM, VISIT, SAFFL
    )

  write_adam(advs, "advs")
}

# =============================================================================
# 8. ADEG – ECG
# =============================================================================
message("\n========== ADEG ==========")

if (!is.null(eg) && nrow(eg) > 0 && !is.null(adsl_key)) {

  adeg <- eg %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ADT   = iso_date(EGDTC),
      ADY   = calc_ady(ADT, TRTSDT),
      AVAL  = suppressWarnings(as.numeric(EGORRES)),
      AVALU = EGORRESU,
      ABLFL = ifelse(!is.na(ADT) & !is.na(TRTSDT) & ADT < TRTSDT, "Y", "N")
    ) %>%
    group_by(USUBJID, EGTESTCD) %>%
    mutate(
      BASE = {
        bl_vals <- AVAL[ABLFL == "Y" & !is.na(AVAL)]
        bl_dts  <- ADT[ABLFL == "Y" & !is.na(AVAL)]
        if (length(bl_vals) > 0) bl_vals[which.max(bl_dts)] else NA_real_
      },
      CHG  = AVAL - BASE,
      PCHG = ifelse(!is.na(BASE) & BASE != 0, (CHG / BASE) * 100, NA_real_)
    ) %>%
    ungroup() %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      EGSEQ, EGTESTCD, EGTEST, any_of("EGCAT"),
      AVAL, AVALU, BASE, CHG, PCHG, ABLFL,
      EGORRES, EGORRESU, any_of(c("EGNRIND", "EGSTAT")),
      EGDTC, ADT, ADY, VISITNUM, VISIT, SAFFL
    )

  write_adam(adeg, "adeg")
}

# =============================================================================
# 9. ADEX – Exposure
# =============================================================================
message("\n========== ADEX ==========")

if (!is.null(ex) && nrow(ex) > 0 && !is.null(adsl_key)) {

  adex <- ex %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ASTDT  = iso_date(EXSTDTC),
      AENDT  = iso_date(EXENDTC),
      ASTDY  = calc_ady(ASTDT, TRTSDT),
      AENDY  = calc_ady(AENDT, TRTSDT),
      AVAL   = suppressWarnings(as.numeric(EXDOSE)),
      AVALU  = EXDOSU,
      ADURU  = as.numeric(AENDT - ASTDT) + 1
    ) %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      EXSEQ, EXTRT, EXDOSE, EXDOSU,
      any_of(c("EXDOSFRM", "EXDOSFRQ", "EXROUTE")),
      AVAL, AVALU, ADURU,
      EXSTDTC, EXENDTC, ASTDT, AENDT, ASTDY, AENDY,
      VISITNUM, VISIT, SAFFL
    )

  write_adam(adex, "adex")
}

# =============================================================================
# 10. ADCM – Concomitant Medications
# =============================================================================
message("\n========== ADCM ==========")

if (!is.null(cm) && nrow(cm) > 0 && !is.null(adsl_key)) {

  adcm <- cm %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ASTDT   = iso_date(CMSTDTC),
      AENDT   = iso_date(CMENDTC),
      ASTDY   = calc_ady(ASTDT, TRTSDT),
      AENDY   = calc_ady(AENDT, TRTSDT),
      TRTEMFL = ifelse(!is.na(ASTDT) & !is.na(TRTSDT) & ASTDT >= TRTSDT, "Y", "N")
    ) %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      CMSEQ, CMTRT, CMDECOD,
      any_of(c("CMCLAS", "CMCLASCD")),
      any_of(c("CMDOSE", "CMDOSU", "CMDOSFRQ", "CMROUTE")),
      CMSTDTC, CMENDTC, ASTDT, AENDT, ASTDY, AENDY,
      TRTEMFL, SAFFL,
      any_of(c("CMINDC", "CMCAT"))
    )

  write_adam(adcm, "adcm")
}

# =============================================================================
# 11. ADPR – Procedures
# =============================================================================
message("\n========== ADPR ==========")

if (!is.null(pr) && nrow(pr) > 0 && !is.null(adsl_key)) {

  adpr <- pr %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ASTDT   = iso_date(PRSTDTC),
      AENDT   = iso_date(PRENDTC),
      ASTDY   = calc_ady(ASTDT, TRTSDT),
      TRTEMFL = ifelse(!is.na(ASTDT) & !is.na(TRTSDT) & ASTDT >= TRTSDT, "Y", "N")
    ) %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      PRSEQ, PRTRT, PRDECOD,
      any_of(c("PRCAT", "PRSCAT", "PRPRESP", "PROCCUR")),
      PRSTDTC, any_of("PRENDTC"), ASTDT, AENDT, ASTDY,
      TRTEMFL, SAFFL,
      any_of(c("VISITNUM", "VISIT"))
    )

  write_adam(adpr, "adpr")
}

# =============================================================================
# 12. ADHO – Hospitalizations
# =============================================================================
message("\n========== ADHO ==========")

if (!is.null(ho) && nrow(ho) > 0 && !is.null(adsl_key)) {

  adho <- ho %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ASTDT  = iso_date(HOSTDTC),
      AENDT  = iso_date(HOENDTC),
      ASTDY  = calc_ady(ASTDT, TRTSDT),
      AENDY  = calc_ady(AENDT, TRTSDT),
      ADURU  = as.numeric(AENDT - ASTDT) + 1,
      TRTEMFL = ifelse(!is.na(ASTDT) & !is.na(TRTSDT) & ASTDT >= TRTSDT, "Y", "N")
    ) %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      HOSEQ,
      any_of(c("HOTERM", "HODECOD", "HOCAT", "HOINDC")),
      HOSTDTC, HOENDTC, ASTDT, AENDT, ASTDY, AENDY, ADURU,
      TRTEMFL, SAFFL
    )

  write_adam(adho, "adho")
}

# =============================================================================
# 13. ADMI – Microscopy (biópsia)
# =============================================================================
message("\n========== ADMI ==========")

if (!is.null(mi) && nrow(mi) > 0 && !is.null(adsl_key)) {

  admi <- mi %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ADT   = iso_date(MIDTC),
      ADY   = calc_ady(ADT, TRTSDT),
      AVAL  = suppressWarnings(as.numeric(MIORRES)),
      AVALU = if ("MIORRESU" %in% names(.)) MIORRESU else NA_character_,
      ABLFL = ifelse(!is.na(ADT) & !is.na(TRTSDT) & ADT < TRTSDT, "Y", "N")
    ) %>%
    group_by(USUBJID, MITESTCD) %>%
    mutate(
      BASE = {
        bl_vals <- AVAL[ABLFL == "Y" & !is.na(AVAL)]
        bl_dts  <- ADT[ABLFL == "Y" & !is.na(AVAL)]
        if (length(bl_vals) > 0) bl_vals[which.max(bl_dts)] else NA_real_
      },
      CHG  = AVAL - BASE
    ) %>%
    ungroup() %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      MISEQ, MITESTCD, MITEST,
      any_of(c("MICAT", "MISCAT", "MISPEC", "MILOC")),
      MIORRES, any_of("MIORRESU"), AVAL, AVALU, BASE, CHG, ABLFL,
      MIDTC, ADT, ADY, VISITNUM, VISIT, SAFFL,
      any_of(c("MISTAT", "MIREASND", "MISTRESC"))
    )

  write_adam(admi, "admi")
}

# =============================================================================
# 14. ADMO – Morphology / Endoscopy findings
# =============================================================================
message("\n========== ADMO ==========")

if (!is.null(mo) && nrow(mo) > 0 && !is.null(adsl_key)) {

  admo <- mo %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ADT   = iso_date(MODTC),
      ADY   = calc_ady(ADT, TRTSDT),
      AVAL  = suppressWarnings(as.numeric(MOORRES)),
      AVALU = if ("MOORRESU" %in% names(.)) MOORRESU else NA_character_,
      ABLFL = ifelse(!is.na(ADT) & !is.na(TRTSDT) & ADT < TRTSDT, "Y", "N")
    ) %>%
    group_by(USUBJID, MOTESTCD) %>%
    mutate(
      BASE = {
        bl_vals <- AVAL[ABLFL == "Y" & !is.na(AVAL)]
        bl_dts  <- ADT[ABLFL == "Y" & !is.na(AVAL)]
        if (length(bl_vals) > 0) bl_vals[which.max(bl_dts)] else NA_real_
      },
      CHG  = AVAL - BASE
    ) %>%
    ungroup() %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      MOSEQ, MOTESTCD, MOTEST,
      any_of(c("MOCAT", "MOSCAT", "MOLOC", "MOSPEC")),
      MOORRES, any_of("MOORRESU"), AVAL, AVALU, BASE, CHG, ABLFL,
      MODTC, ADT, ADY, VISITNUM, VISIT, SAFFL,
      any_of(c("MOSTAT", "MOREASND", "MOSTRESC", "MORESCAT"))
    )

  write_adam(admo, "admo")
}

# =============================================================================
# 15. ADQS – Questionnaires (todos os QS consolidados)
# =============================================================================
message("\n========== ADQS ==========")

# Lista de todos os QS disponíveis com seus SUPPs
qs_list <- list(
  list(df = qsai,  supp = supp_qsai,  label = "QSAI"),
  list(df = qsca,  supp = supp_qsca,  label = "QSCA"),
  list(df = qscc,  supp = supp_qscc,  label = "QSCC"),
  list(df = qscd,  supp = supp_qscd,  label = "QSCD"),
  list(df = qscr,  supp = supp_qscr,  label = "QSCR"),
  list(df = qseq,  supp = supp_qseq,  label = "QSEQ"),
  list(df = qser,  supp = supp_qser,  label = "QSER"),
  list(df = qsib,  supp = supp_qsib,  label = "QSIB"),
  list(df = qspg,  supp = supp_qspg,  label = "QSPG"),
  list(df = qss2,  supp = supp_qss2,  label = "QSS2"),
  list(df = qsts,  supp = supp_qsts,  label = "QSTS"),
  list(df = qswp,  supp = supp_qswp,  label = "QSWP")
)

process_qs <- function(qs_df, qs_supp, label) {
  if (is.null(qs_df) || nrow(qs_df) == 0) return(NULL)
  if (!is.null(qs_supp) && nrow(qs_supp) > 0)
    qs_df <- merge_supp(qs_df, qs_supp, "QSSEQ")
  qs_df %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ADT   = iso_date(QSDTC),
      ADY   = calc_ady(ADT, TRTSDT),
      AVAL  = suppressWarnings(as.numeric(QSORRES)),
      AVALC = as.character(QSORRES),
      ABLFL = ifelse(!is.na(ADT) & !is.na(TRTSDT) & ADT < TRTSDT, "Y", "N"),
      QSDOM = label
    ) %>%
    group_by(USUBJID, QSTESTCD) %>%
    mutate(
      BASE = {
        bl_vals <- AVAL[ABLFL == "Y" & !is.na(AVAL)]
        bl_dts  <- ADT[ABLFL == "Y" & !is.na(AVAL)]
        if (length(bl_vals) > 0) bl_vals[which.max(bl_dts)] else NA_real_
      },
      CHG  = AVAL - BASE,
      PCHG = ifelse(!is.na(BASE) & BASE != 0, (CHG / BASE) * 100, NA_real_)
    ) %>%
    ungroup() %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD, QSDOM,
      QSSEQ, QSTESTCD, QSTEST, any_of(c("QSCAT", "QSSCAT")),
      QSORRES, AVAL, AVALC, BASE, CHG, PCHG, ABLFL,
      QSDTC, ADT, ADY,
      any_of(c("VISITNUM", "VISIT")), SAFFL
    )
}

qs_combined <- map(qs_list, ~ process_qs(.x$df, .x$supp, .x$label)) %>%
  compact() %>%
  bind_rows()

if (nrow(qs_combined) > 0) {
  write_adam(qs_combined, "adqs")
} else {
  message("  [SKIP] ADQS — nenhum dado QS disponível")
}

# =============================================================================
# 16. ADMH – Medical History
# =============================================================================
message("\n========== ADMH ==========")

if (!is.null(mh) && nrow(mh) > 0 && !is.null(adsl_key)) {

  # Merge SUPPMH se disponível
  if (!is.null(supp_mh) && nrow(supp_mh) > 0)
    mh <- merge_supp(mh, supp_mh, "MHSEQ")

  admh <- mh %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ASTDT = iso_date(MHSTDTC),
      AENDT = iso_date(MHENDTC),
      ASTDY = calc_ady(ASTDT, TRTSDT)
    ) %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      MHSEQ, MHTERM, MHDECOD, MHBODSYS,
      any_of(c("MHSOC", "MHCAT", "MHSCAT", "MHPRESP", "MHOCCUR")),
      MHSTDTC, any_of("MHENDTC"), ASTDT, AENDT, ASTDY,
      SAFFL
    )

  write_adam(admh, "admh")
}

# =============================================================================
# 17. ADSU – Substance Use
# =============================================================================
message("\n========== ADSU ==========")

if (!is.null(su) && nrow(su) > 0 && !is.null(adsl_key)) {

  if (!is.null(supp_su) && nrow(supp_su) > 0)
    su <- merge_supp(su, supp_su, "SUSEQ")

  adsu <- su %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ASTDT = iso_date(SUSTDTC),
      AENDT = iso_date(SUENDTC),
      ASTDY = calc_ady(ASTDT, TRTSDT),
      AVAL  = suppressWarnings(as.numeric(SUORRES)),
      AVALU = if ("SUORRESU" %in% names(.)) SUORRESU else NA_character_
    ) %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      SUSEQ, SUTRT, SUDECOD, any_of(c("SUCAT", "SUSCAT")),
      SUORRES, AVAL, AVALU,
      SUSTDTC, any_of("SUENDTC"), ASTDT, AENDT, ASTDY,
      SAFFL
    )

  write_adam(adsu, "adsu")
}

# =============================================================================
# 18. ADPC – Pharmacokinetics Concentrations
# =============================================================================
message("\n========== ADPC ==========")

if (!is.null(pc) && nrow(pc) > 0 && !is.null(adsl_key)) {

  if (!is.null(supp_pc) && nrow(supp_pc) > 0)
    pc <- merge_supp(pc, supp_pc, "PCSEQ")

  adpc <- pc %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ADT   = iso_date(PCDTC),
      ADY   = calc_ady(ADT, TRTSDT),
      AVAL  = suppressWarnings(as.numeric(PCORRES)),
      AVALU = if ("PCORRESU" %in% names(.)) PCORRESU else NA_character_
    ) %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      PCSEQ, PCTESTCD, PCTEST, any_of(c("PCCAT", "PCSPEC")),
      PCORRES, AVAL, AVALU,
      PCDTC, ADT, ADY,
      any_of(c("VISITNUM", "VISIT", "PCRFTDTC", "PCELTM")),
      SAFFL
    )

  write_adam(adpc, "adpc")
}

# =============================================================================
# 19. ADPE – Physical Examination
# =============================================================================
message("\n========== ADPE ==========")

if (!is.null(pe) && nrow(pe) > 0 && !is.null(adsl_key)) {

  if (!is.null(supp_pe) && nrow(supp_pe) > 0)
    pe <- merge_supp(pe, supp_pe, "PESEQ")

  adpe <- pe %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ADT   = iso_date(PEDTC),
      ADY   = calc_ady(ADT, TRTSDT),
      ABLFL = ifelse(!is.na(ADT) & !is.na(TRTSDT) & ADT < TRTSDT, "Y", "N")
    ) %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      PESEQ, PETESTCD, PETEST, any_of(c("PECAT", "PEORRES", "PESTRESC")),
      PEDTC, ADT, ADY, ABLFL,
      any_of(c("VISITNUM", "VISIT")), SAFFL
    )

  write_adam(adpe, "adpe")
}

# =============================================================================
# 20. ADSV – Subject Visits
# =============================================================================
message("\n========== ADSV ==========")

if (!is.null(sv) && nrow(sv) > 0 && !is.null(adsl_key)) {

  if (!is.null(supp_sv) && nrow(supp_sv) > 0)
    sv <- merge_supp(sv, supp_sv, "SVSEQ")

  adsv <- sv %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      SVSTDT = iso_date(SVSTDTC),
      SVENDT = iso_date(SVENDTC),
      SVSTDY = calc_ady(SVSTDT, TRTSDT)
    ) %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      SVSEQ, VISITNUM, VISIT,
      SVSTDTC, any_of("SVENDTC"), SVSTDT, SVENDT, SVSTDY,
      SAFFL
    )

  write_adam(adsv, "adsv")
}

# =============================================================================
# 21. ADDV – Protocol Deviations
# =============================================================================
message("\n========== ADDV ==========")

if (!is.null(dv) && nrow(dv) > 0 && !is.null(adsl_key)) {

  if (!is.null(supp_dv) && nrow(supp_dv) > 0) {
    supp_dv <- read_supp("dv")
    dv <- merge_supp(dv, supp_dv, "DVSEQ")
  }

  addv <- dv %>%
    left_join(adsl_key, by = "USUBJID") %>%
    mutate(
      ASTDT = iso_date(DVSTDTC),
      ASTDY = calc_ady(ASTDT, TRTSDT)
    ) %>%
    select(
      STUDYID, USUBJID, ARM, ARMCD,
      DVSEQ, DVTERM, DVDECOD, any_of(c("DVCAT", "DVSCAT")),
      DVSTDTC, ASTDT, ASTDY,
      SAFFL
    )

  write_adam(addv, "addv")
}

# =============================================================================
# SUMÁRIO FINAL
# =============================================================================
message("\n========================================")
message("         CONVERSÃO CONCLUÍDA")
message("========================================")

adam_files <- list.files(adam_dir, pattern = "\\.csv$", full.names = TRUE)
total_obs  <- 0

for (fp in adam_files) {
  df <- tryCatch(read.csv(fp, nrows = 1, check.names = FALSE), error = function(e) NULL)
  n_lines <- as.integer(system(paste("wc -l <", shQuote(fp)), intern = TRUE))
  n_obs   <- max(0, n_lines - 1)
  n_vars  <- if (!is.null(df)) ncol(df) else NA
  total_obs <- total_obs + n_obs
  message(sprintf("  %-12s  %6d obs  %3d vars", basename(fp), n_obs, n_vars))
}

message(sprintf("\nTotal de observações geradas: %d", total_obs))
message(sprintf("Arquivos ADaM salvos em: %s", adam_dir))
message("========================================\n")
