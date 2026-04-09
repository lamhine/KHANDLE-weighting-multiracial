#------------------------------------------------------------------
# Script: 03_clean_KHANDLE_othervars.R
# Purpose: Clean and process KHANDLE data other than race
# Original author: Taylor Mobley (adapted by Tracy Lam-Hine)
# Repo: github.com/lamhine/KHANDLE-weighting-multiracial
# Last updated: 2026-01-12
# Notes:
#   - Supports both t10 and t16 analysis setups
#   - Loads race_mece_brfss5 from 01_clean_KHANDLE_race.R output
#   - Avoids hard stops for non-core variables (defensive / portable)
#------------------------------------------------------------------

source("config/config.R")

#------------------------------------------------------------------
# Load KHANDLE data WITH cleaned race already attached
#------------------------------------------------------------------
khandle_cleaned_race <- readRDS(file.path(processed_data_dir, "02_khandle_cleaned_race.rds"))

#------------------------------------------------------------------
# Utilities: to_num(), get_col(), na_if_in(), recode_binary01_khandle()
# provided by R/utils.R via config
#------------------------------------------------------------------

#------------------------------------------------------------------
# Decide which analysis setup we are in (t10 vs t16)
#   If you have a flag in config (recommended), use it; otherwise default to "t16".
#   This only affects *which subset* we save as "core" vs "full".
#------------------------------------------------------------------

# analysis_tag is set by config/config.R

#------------------------------------------------------------------
# 1) Variable selection
#   Keep a *core required* set, plus a wide optional set.
#   We do NOT stop for missing optional vars.
#------------------------------------------------------------------

vars_keep_required <- c(
  "STUDYID",
  "W1_INTERVIEW_AGE",
  "W1_D_GENDER",
  "race_mece_brfss5"
)

vars_keep_optional <- c(
  # SES / demographics
  "W1_D_EDUCATION",
  "W1_MARITAL_STATUS",
  "W1_INCOME_RANGE",
  "W1_INCOME_WORRY",
  "W1_COUNTRY_BORN",
  "W1_LADDER1",
  "W1_MILITARY",
  "W1_LANGUAGE",
  
  # Retirement / employment
  "W1_EMP_RETIRED_AGE",
  "W1_EMP_NOTWORKING",
  "W1_EMP_FULLTIME",
  "W1_EMP_PARTTIME",
  "W1_EMP_DISABLE",
  "W1_EMP_HOMEMAKER",
  "W1_EMP_WORKING",
  "W1_EMP_OTHER",
  "W1_EMP_DK",
  "W1_EMP_REFUSED",
  
  # Growing up conditions
  "W1_GROWINGUP_FINANCE",
  "W1_GROWINGUP_GOHUNGRY",
  "W1_GROWINGUP_HOUSING",
  
  # Cognition / depression
  "W1_D_SENAS_AVG_COGNITION",
  "W1_D_SENAS_AVG_COGNITION_Z",
  "W1_D_SENAS_EXEC_Z",
  "W1_D_SENAS_SEM_Z",
  "W1_D_SENAS_VRMEM_Z",
  "W1_NIHTLBX_DEPR_RAW",
  "W1_NIHTLBX_DEPR_THETA",
  
  # Health / behaviors / impairments
  "W1_HEALTH",
  "W1_PAIN",
  "W1_SMK",
  "W1_SMK_NOW",
  "W1_SENSIMP_VISION",
  "W1_SENSIMP_HEARING",
  
  # ADL/IADL
  "W1_DAILY_LIVING_ADL1", "W1_DAILY_LIVING_ADL2", "W1_DAILY_LIVING_ADL3",
  "W1_DAILY_LIVING_ADL4", "W1_DAILY_LIVING_ADL5", "W1_DAILY_LIVING_ADL6",
  "W1_DAILY_LIVING_ADL7", "W1_DAILY_LIVING_ADL8", "W1_DAILY_LIVING_ADL9",
  "W1_DAILY_LIVING_IADL1", "W1_DAILY_LIVING_IADL2", "W1_DAILY_LIVING_IADL3",
  
  # eCog
  "W1_ECOG_CONCERNED_THINKING",
  
  # Physical activity
  "W1_PA_HVY_WRK", "W1_PA_LT_EX", "W1_PA_LT_HSE", "W1_PA_LT_WRK",
  "W1_PA_VIG_EX", "W1_PA_VIG_HSE",
  
  # Cog impairment probabilities
  "W1_COGIMP_PROB_ADJ_BL_DX", "W1_COGIMP_PROB_ADJ_BL_NODX",
  "W1_COGIMP_PROB_FIN_DX", "W1_COGIMP_PROB_FIN_NODX"
)

missing_required <- setdiff(vars_keep_required, names(khandle_cleaned_race))
if (length(missing_required) > 0) {
  stop("Missing REQUIRED variables in khandle_cleaned_race: ", paste(missing_required, collapse = ", "))
}

vars_keep <- c(vars_keep_required, intersect(vars_keep_optional, names(khandle_cleaned_race)))

khandle_clean <- khandle_cleaned_race[, vars_keep, drop = FALSE]
khandle_clean <- as.data.frame(khandle_clean)
names(khandle_clean) <- tolower(names(khandle_clean))

#------------------------------------------------------------------
# 2) Harmonized / cleaned baseline variables
#------------------------------------------------------------------

# Age
khandle_clean$age <- to_num(khandle_clean$w1_interview_age)

# Sex (male)
khandle_clean$male <- dplyr::case_when(
  to_num(khandle_clean$w1_d_gender) == 1 ~ 1,
  to_num(khandle_clean$w1_d_gender) == 2 ~ 0,
  TRUE ~ NA_real_
)

# Education (W1_D_EDUCATION): collapse 9/10 -> NA (as in your script)
khandle_clean$education <- NA_real_
if ("w1_d_education" %in% names(khandle_clean)) {
  khandle_clean$education <- to_num(khandle_clean$w1_d_education)
  khandle_clean$education[khandle_clean$education %in% c(9, 10)] <- NA_real_
}

# Retired age (89.99 -> NA)
khandle_clean$retired_age <- NA_real_
if ("w1_emp_retired_age" %in% names(khandle_clean)) {
  ra <- to_num(khandle_clean$w1_emp_retired_age)
  ra[ra == 89.99] <- NA_real_
  khandle_clean$retired_age <- ra
}

# Keep SENAS + depression as-is (typically already NA’d upstream)
# (No changes needed; leaving raw fields intact is safest for MI.)

#------------------------------------------------------------------
# 3) Copy “raw” items into harmonized names
#    Then recode 88/99 (and only those) to NA for those survey-coded items.
#------------------------------------------------------------------

# Map of newname -> oldname (lowercase)
map_vars <- c(
  marital_status      = "w1_marital_status",
  income_range        = "w1_income_range",
  income_worry        = "w1_income_worry",
  ladder              = "w1_ladder1",
  country_born        = "w1_country_born",
  military            = "w1_military",
  language            = "w1_language",
  
  growingup_finance   = "w1_growingup_finance",
  growingup_gohungry  = "w1_growingup_gohungry",
  growingup_housing   = "w1_growingup_housing",
  
  health              = "w1_health",
  pain                = "w1_pain",
  smk                 = "w1_smk",
  smk_now             = "w1_smk_now",
  
  sensimp_vision      = "w1_sensimp_vision",
  sensimp_hearing     = "w1_sensimp_hearing",
  
  adl1                = "w1_daily_living_adl1",
  adl2                = "w1_daily_living_adl2",
  adl3                = "w1_daily_living_adl3",
  adl4                = "w1_daily_living_adl4",
  adl5                = "w1_daily_living_adl5",
  adl6                = "w1_daily_living_adl6",
  adl7                = "w1_daily_living_adl7",
  adl8                = "w1_daily_living_adl8",
  adl9                = "w1_daily_living_adl9",
  iadl1               = "w1_daily_living_iadl1",
  iadl2               = "w1_daily_living_iadl2",
  iadl3               = "w1_daily_living_iadl3",
  
  concerned_thinking  = "w1_ecog_concerned_thinking",
  
  pa_hvy_wrk          = "w1_pa_hvy_wrk",
  pa_lt_ex            = "w1_pa_lt_ex",
  pa_lt_hse           = "w1_pa_lt_hse",
  pa_lt_wrk           = "w1_pa_lt_wrk",
  pa_vig_ex           = "w1_pa_vig_ex",
  pa_vig_hse          = "w1_pa_vig_hse"
)

for (nm in names(map_vars)) {
  src <- map_vars[[nm]]
  if (src %in% names(khandle_clean)) {
    khandle_clean[[nm]] <- khandle_clean[[src]]
  } else {
    khandle_clean[[nm]] <- NA
  }
}

# Recode 88/99 -> NA for those survey-coded fields (and only those)
recode_88_99 <- names(map_vars)
recode_88_99 <- recode_88_99[recode_88_99 %in% names(khandle_clean)]

for (v in recode_88_99) {
  x <- to_num(khandle_clean[[v]])
  x[x %in% c(88, 99)] <- NA_real_
  khandle_clean[[v]] <- x
}

#------------------------------------------------------------------
# 4) Restrict to modeling race groups + create race_summary
#     Modeling groups: Asian, Black, Hispanic, White, Multiracial
#     (Explicitly excludes other categories such as American Indian or Alaska Native)
#------------------------------------------------------------------

# Restrict analytic sample to the race groups used for modeling (ALLOWED_RACE5 from R/race_constants.R)
khandle_clean <- subset(
  khandle_clean,
  !is.na(race_mece_brfss5) & race_mece_brfss5 %in% ALLOWED_RACE5
)

# Numeric race code used in downstream models
# 1=Asian, 2=Black, 3=Hispanic, 4=White, 5=Multiracial
khandle_clean$race_summary <- dplyr::case_when(
  khandle_clean$race_mece_brfss5 == "Asian" ~ 1,
  khandle_clean$race_mece_brfss5 == "Black" ~ 2,
  khandle_clean$race_mece_brfss5 == "Hispanic" ~ 3,
  khandle_clean$race_mece_brfss5 == "White" ~ 4,
  khandle_clean$race_mece_brfss5 == "Multiracial" ~ 5,
  TRUE ~ NA_real_
)

#------------------------------------------------------------------
# 5) Calculated variables
#------------------------------------------------------------------

# Retired indicator
# NOTE: your prior logic used NOTWORKING/FULLTIME/PARTTIME in a way that can’t all mean “retired”.
# Here is a safer approach:
#   - If "working" indicators exist, define working first.
#   - Define retired as NOT working, AND either has a retired_age reported OR has explicit notworking==1.
# This is conservative and reduces accidental 1s.
notworking <- to_num(get_col(khandle_clean, "w1_emp_notworking"))
fulltime   <- to_num(get_col(khandle_clean, "w1_emp_fulltime"))
parttime   <- to_num(get_col(khandle_clean, "w1_emp_parttime"))
working_i  <- to_num(get_col(khandle_clean, "w1_emp_working"))

# Working: any explicit working flag OR FT/PT
khandle_clean$working <- dplyr::case_when(
  working_i == 1 | fulltime == 1 | parttime == 1 ~ 1,
  working_i == 0 & fulltime == 0 & parttime == 0 ~ 0,
  TRUE ~ NA_real_
)

# Retired: retired_age observed OR notworking==1, but not currently working
khandle_clean$retired <- dplyr::case_when(
  khandle_clean$working == 1 ~ 0,
  khandle_clean$working == 0 & (!is.na(khandle_clean$retired_age) | notworking == 1) ~ 1,
  khandle_clean$working == 0 & (is.na(khandle_clean$retired_age) & notworking == 0) ~ 0,
  TRUE ~ NA_real_
)

# Smoking status (mirror BRFSS logic)
# smk: ever smoker; smk_now: current smoker
khandle_clean$smoke_status <- dplyr::case_when(
  khandle_clean$smk == 1 & is.na(khandle_clean$smk_now) ~ 1,
  khandle_clean$smk == 1 & khandle_clean$smk_now == 1 ~ 2,
  khandle_clean$smk == 1 ~ 1,
  khandle_clean$smk == 0 ~ 0,
  TRUE ~ NA_real_
)

# smk_now variant where NA's are set to 0 when never-smoker
khandle_clean$smk_now_v2 <- dplyr::case_when(
  khandle_clean$smk == 0 & is.na(khandle_clean$smk_now) ~ 0,
  TRUE ~ khandle_clean$smk_now
)

# Born in US (only define when country_born observed)
khandle_clean$born_us <- dplyr::case_when(
  is.na(khandle_clean$country_born) ~ NA_real_,
  khandle_clean$country_born == 1 ~ 1,
  TRUE ~ 0
)

# English interview (only define when language observed)
# Your prior: 1033=English, else 0. Keep same.
khandle_clean$english_interview <- dplyr::case_when(
  is.na(khandle_clean$language) ~ NA_real_,
  to_num(khandle_clean$language) == 1033 ~ 1,
  TRUE ~ 0
)

# eCog concerned thinking (binary)
khandle_clean$concerned_thinking_v2 <- dplyr::case_when(
  is.na(khandle_clean$concerned_thinking) ~ NA_real_,
  khandle_clean$concerned_thinking == 1 ~ 1,
  TRUE ~ 0
)

#------------------------------------------------------------------
# 6) Create “core” vs “full” outputs (t10 vs t16 friendly)
#   - FULL keeps everything you pulled + derived variables.
#   - CORE is a minimal set used in most weighting models, and is what I’d recommend
#     you point the downstream pipeline at unless you confirm you need extras.
#------------------------------------------------------------------

core_vars <- c(
  "studyid",
  "age",
  "male",
  "race_mece_brfss5",
  "race_summary",
  "education",
  "marital_status",
  "income_range",
  "income_worry",
  "born_us",
  "english_interview",
  "military",
  "ladder",
  "working",
  "retired",
  "smk",
  "smk_now_v2",
  "smoke_status",
  "health"
)

core_vars <- core_vars[core_vars %in% names(khandle_clean)]
khandle_core <- khandle_clean[, core_vars, drop = FALSE]

#------------------------------------------------------------------
# 7) Save artifacts
#------------------------------------------------------------------

# Full dataset (superset; safest for MI + future model changes)
saveRDS(khandle_clean, file = file.path(processed_data_dir, "03_khandle_clean_full.rds"))

# Core dataset (recommended for the actual weighting model unless you need extras)
saveRDS(khandle_core, file = file.path(processed_data_dir, "03_khandle_clean_core.rds"))

# Backward-compatible name: if prior pipeline expects "KHANDLE_clean.rds",
# point it at CORE (most consistent with how the BRFSS side is used).
saveRDS(khandle_core, file = file.path(processed_data_dir, "KHANDLE_clean.rds"))

# Optional CSV (core only; smaller + less fragile)
khandle_core_for_csv <- khandle_core
khandle_core_for_csv[is.na(khandle_core_for_csv)] <- "."
write.csv(khandle_core_for_csv, file.path(processed_data_dir, "KHANDLE_clean.csv"), row.names = FALSE)

message("03_clean_KHANDLE_othervars.R completed successfully. Saved: ",
        "\n - 03_khandle_clean_full.rds",
        "\n - 03_khandle_clean_core.rds",
        "\n - KHANDLE_clean.rds (core)",
        "\n - KHANDLE_clean.csv (core)")