#------------------------------------------------------------------
# Script: 02_harmonize_KHANDLE_BRFSS.R
# Purpose: Harmonize imputed KHANDLE + BRFSS into one combined dataset
# Original author: Eleanor Hayes-Larson (modified by Tracy Lam-Hine)
# Repo: github.com/lamhine/KHANDLE-weighting-multiracial
# Last updated: 2026-01-13
# Notes:
#   - Compatible with both t10_replication and t16_update tracks
#   - Assumes 01_multiple_imputation.R outputs:
#       processed_data_dir/khandle_pmm.rds
#       processed_data_dir/brfss_pmm.rds
#------------------------------------------------------------------

source("config/config.R")

# Helpers: to_num(), assert_has_vars(), race_to_code() provided by R/utils.R and R/race_constants.R via config

#------------------------------------------------------------------
# Load imputed long datasets (created by 01_multiple_imputation.R)
#------------------------------------------------------------------

khandle_pmm <- readRDS(file.path(processed_data_dir, "khandle_pmm.rds"))
brfss_pmm   <- readRDS(file.path(processed_data_dir, "brfss_pmm.rds"))

#------------------------------------------------------------------
# Detect analysis setup from BRFSS years (align harmonization to actual data)
#------------------------------------------------------------------

analysis_setup <- detect_analysis_setup(brfss_pmm)
message("02_harmonize_KHANDLE_BRFSS.R detected analysis_setup = ", analysis_setup)

#------------------------------------------------------------------
# KHANDLE: required columns from khandle_pmm (new naming)
#------------------------------------------------------------------

k_req <- c(
  "studyid", "imp",
  "w1_interview_age", "male", "race_cat",
  "marital_status", "education", "income_range", "military",
  "health", "smoke_status",
  # Optional in some tracks, but you’ve been carrying them in full:
  "sensimp_vision", "adl1", "adl2", "adl9", "pa_lt_ex", "pa_vig_ex",
  "w1_cogimp_prob_fin_dx"
)
assert_has_vars(khandle_pmm, intersect(k_req, names(khandle_pmm)), "KHANDLE PMM")

khandle_harm <- khandle_pmm %>%
  transmute(
    id_h   = as.character(studyid),
    imp_h  = to_num(imp),
    brfss_h = 0,
    
    age_h  = to_num(w1_interview_age),
    male_h = to_num(male),
    
    # race_cat is a label in PMM
    race_summary_h = race_to_code(as.character(race_cat)),
    
    # KHANDLE english_interview is not part of your current minimal MI input;
    # keep as NA (so the combined dataset has the field consistently)
    english_interview_h = NA_real_,
    
    marital_h = dplyr::case_when(
      to_num(marital_status) %in% 1:2 ~ 1,
      to_num(marital_status) %in% 3:6 ~ 0,
      TRUE ~ NA_real_
    ),
    
    education_h = factor(dplyr::case_when(
      to_num(education) %in% c(0, 1) ~ 1,
      to_num(education) == 2         ~ 2,
      to_num(education) == 7         ~ 3,
      to_num(education) %in% 3:4     ~ 4,
      to_num(education) %in% 5:6     ~ 5,
      TRUE ~ NA_real_
    )),
    
    military_h = to_num(military),
    
    goodhealth_h = dplyr::case_when(
      to_num(health) %in% 1:3 ~ 1,
      to_num(health) %in% 4:5 ~ 0,
      TRUE ~ NA_real_
    ),
    
    blind_h = dplyr::case_when(
      to_num(sensimp_vision) %in% 1:4 ~ 0,
      to_num(sensimp_vision) %in% c(5, 6) ~ 1,
      TRUE ~ NA_real_
    ),
    
    smoke_status_h = dplyr::case_when(
      is.na(to_num(smoke_status)) ~ NA_real_,
      to_num(smoke_status) == 2 ~ 1,
      to_num(smoke_status) %in% c(0, 1) ~ 0,
      TRUE ~ NA_real_
    ),
    
    adl_walking_h  = dplyr::case_when(
      is.na(to_num(adl1)) & is.na(to_num(adl2)) ~ NA_real_,
      (to_num(adl1) %in% 3:5) | (to_num(adl2) %in% 3:5) ~ 1,
      TRUE ~ 0
    ),
    
    adl_dressing_h = dplyr::case_when(
      is.na(to_num(adl9)) ~ NA_real_,
      to_num(adl9) %in% 3:5 ~ 1,
      TRUE ~ 0
    ),
    
    exercise_h = dplyr::case_when(
      is.na(to_num(pa_lt_ex)) & is.na(to_num(pa_vig_ex)) ~ NA_real_,
      (to_num(pa_lt_ex) %in% 1:3) | (to_num(pa_vig_ex) %in% 1:3) ~ 1,
      TRUE ~ 0
    ),
    
    # KHANDLE income per-person (marital approximates household size)
    income_pp_h = dplyr::case_when(
      to_num(income_range) == 1  ~ 10000,
      to_num(income_range) == 2  ~ 15000,
      to_num(income_range) == 3  ~ 20000,
      to_num(income_range) == 4  ~ 25000,
      to_num(income_range) == 5  ~ 35000,
      to_num(income_range) == 6  ~ 45000,
      to_num(income_range) == 7  ~ 55000,
      to_num(income_range) == 8  ~ 65000,
      to_num(income_range) == 9  ~ 75000,
      to_num(income_range) == 10 ~ 100000,
      to_num(income_range) == 11 ~ 125000,
      to_num(income_range) == 12 ~ 150000,
      to_num(income_range) == 13 ~ 175000,
      TRUE ~ NA_real_
    ) / sqrt(dplyr::if_else(to_num(marital_status) %in% 1:2, 2, 1)),
    
    brfss_sampwt_h = 1,
    county_code = NA,
    
    cogimp_prob_fin_dx = to_num(w1_cogimp_prob_fin_dx)
  )

#------------------------------------------------------------------
# BRFSS harmonization
#------------------------------------------------------------------

# living wage lookup is optional; if county missing, income_pp_h becomes NA anyway
lw_path <- file.path(raw_data_dir, "livingwage_MIT.csv")
livingwage <- NULL
if (file.exists(lw_path)) {
  livingwage <- read.csv(lw_path) %>%
    transmute(
      county_code = as.character(FIPS.code),
      Living.wage = suppressWarnings(as.numeric(Living.wage))
    )
}

# Required BRFSS PMM fields (based on your new MI script)
b_req <- c(
  ".id", "imp",
  "year", "age", "male", "race_cat",
  "marital_status", "education", "income", "military",
  "vision_impair", "english_interview", "health", "smoke_status",
  "diff_walk", "diff_dress", "phys_activity",
  "hhsize_v2", "county_code", "cdc_finalwt"
)
assert_has_vars(brfss_pmm, b_req, "BRFSS PMM")

brfss_imp2 <- brfss_pmm %>%
  mutate(
    age = dplyr::if_else(to_num(age) <= 90, to_num(age), 90),
    county_code = as.character(county_code)
  )

if (!is.null(livingwage)) {
  brfss_imp2 <- brfss_imp2 %>% left_join(livingwage, by = "county_code")
} else {
  brfss_imp2 <- brfss_imp2 %>% mutate(Living.wage = NA_real_)
}

# Income mapping differs between older (t10) vs newer (t16) BRFSS codebooks
# We keep both mappings but select based on analysis_setup.
income_to_cont_t10 <- function(income) {
  dplyr::case_when(
    income == 1 ~ 10000,
    income == 2 ~ 15000,
    income == 3 ~ 20000,
    income == 4 ~ 25000,
    income == 5 ~ 35000,
    income == 6 ~ 50000,
    income == 7 ~ 75000,
    income == 8 ~ 100000,
    income == 9 ~ 125000,
    income == 10 ~ 175000,
    TRUE ~ NA_real_
  )
}

income_to_cont_t16 <- function(income) {
  # If t16 uses the same coding in your cleaned file, keep it aligned.
  # If you later confirm different bins for 2019–2023, update here centrally.
  income_to_cont_t10(income)
}

brfss_imp2 <- brfss_imp2 %>%
  mutate(
    income_num = to_num(income),
    income_cont = if (analysis_setup == "t10") income_to_cont_t10(income_num) else income_to_cont_t16(income_num),
    
    # MIT living wage correction only if Living.wage present
    income_corrected = dplyr::case_when(
      is.na(income_cont) ~ NA_real_,
      is.na(Living.wage) ~ income_cont,
      TRUE ~ income_cont / (Living.wage / 16.5310111)
    ),
    
    income_pp_h = dplyr::case_when(
      is.na(income_corrected) ~ NA_real_,
      is.na(to_num(hhsize_v2)) ~ NA_real_,
      to_num(hhsize_v2) <= 0 ~ NA_real_,
      TRUE ~ income_corrected / sqrt(to_num(hhsize_v2))
    )
  )

brfss_harm <- brfss_imp2 %>%
  transmute(
    id_h  = as.character(.id),
    imp_h = to_num(imp),
    brfss_h = 1,
    
    age_h  = to_num(age),
    male_h = to_num(male),
    
    race_summary_h = race_to_code(as.character(race_cat)),
    
    english_interview_h = to_num(english_interview),
    
    marital_h = dplyr::case_when(
      to_num(marital_status) %in% c(1, 6) ~ 1,
      to_num(marital_status) %in% 2:5 ~ 0,
      TRUE ~ NA_real_
    ),
    
    education_h = factor(dplyr::case_when(
      to_num(education) %in% c(1, 2) ~ 1,
      to_num(education) == 3 ~ 2,
      to_num(education) == 4 ~ 3,
      to_num(education) == 5 ~ 4,
      to_num(education) == 6 ~ 5,
      TRUE ~ NA_real_
    )),
    
    military_h = to_num(military),
    
    goodhealth_h = dplyr::case_when(
      to_num(health) %in% 1:3 ~ 1,
      to_num(health) %in% 4:5 ~ 0,
      TRUE ~ NA_real_
    ),
    
    blind_h = to_num(vision_impair),
    
    smoke_status_h = dplyr::case_when(
      is.na(to_num(smoke_status)) ~ NA_real_,
      to_num(smoke_status) == 2 ~ 1,
      to_num(smoke_status) %in% c(0, 1) ~ 0,
      TRUE ~ NA_real_
    ),
    
    adl_walking_h  = to_num(diff_walk),
    adl_dressing_h = to_num(diff_dress),
    exercise_h     = to_num(phys_activity),
    
    income_pp_h = income_pp_h,
    
    brfss_sampwt_h = to_num(cdc_finalwt),
    county_code = county_code,
    
    cogimp_prob_fin_dx = NA_real_
  )

#------------------------------------------------------------------
# Combine + derived SES variables
#------------------------------------------------------------------

khandle_brfss_harmonized <- bind_rows(khandle_harm, brfss_harm)

# Weighted BRFSS income quantiles (only from BRFSS; requires non-missing income_pp_h and weights)
wtd_quarts <- Hmisc::wtd.quantile(
  brfss_harm$income_pp_h,
  weights = brfss_harm$brfss_sampwt_h,
  probs = c(0.25, 0.5, 0.75),
  na.rm = TRUE
)

khandle_brfss_harmonized <- khandle_brfss_harmonized %>%
  mutate(
    income_gt50k_pp_h = dplyr::case_when(
      is.na(income_pp_h) ~ NA_real_,
      income_pp_h > 50000 ~ 1,
      TRUE ~ 0
    ),
    income_gt65k_pp_h = dplyr::case_when(
      is.na(income_pp_h) ~ NA_real_,
      income_pp_h > 65000 ~ 1,
      TRUE ~ 0
    ),
    income_gtmed_pp_h = dplyr::case_when(
      is.na(income_pp_h) ~ NA_real_,
      income_pp_h > wtd_quarts[2] ~ 1,
      TRUE ~ 0
    ),
    
    income_quartile = factor(dplyr::case_when(
      is.na(income_pp_h) ~ NA_real_,
      income_pp_h <= wtd_quarts[1] ~ 1,
      income_pp_h <= wtd_quarts[2] ~ 2,
      income_pp_h <= wtd_quarts[3] ~ 3,
      TRUE ~ 4
    )),
    
    education2_h = factor(dplyr::case_when(
      education_h %in% c(1, 2, 3) ~ 1,
      education_h == 4 ~ 2,
      education_h == 5 ~ 3,
      TRUE ~ NA_real_
    )),
    
    education3_h = factor(dplyr::case_when(
      education_h %in% c(1, 2) ~ 1,
      education_h == 3 ~ 2,
      education_h == 4 ~ 3,
      education_h == 5 ~ 4,
      TRUE ~ NA_real_
    ))
  )

#------------------------------------------------------------------
# Quick checks
#------------------------------------------------------------------

cat("\nRace by source (rows):\n")
print(table(khandle_brfss_harmonized$race_summary_h,
            khandle_brfss_harmonized$brfss_h,
            useNA = "ifany"))

cat("\nBRFSS county missingness (imp 1):\n")
print(
  khandle_brfss_harmonized %>%
    filter(brfss_h == 1, imp_h == 1) %>%
    summarise(
      n = n(),
      n_missing_county = sum(is.na(county_code)),
      pct_missing_county = mean(is.na(county_code))
    )
)

cat("\nIncome_pp_h missingness by source (imp 1):\n")
print(
  khandle_brfss_harmonized %>%
    filter(imp_h == 1) %>%
    group_by(brfss_h) %>%
    summarise(
      n = n(),
      pct_missing_income_pp = mean(is.na(income_pp_h)),
      .groups = "drop"
    )
)

#------------------------------------------------------------------
# Save final dataset
#------------------------------------------------------------------

saveRDS(khandle_brfss_harmonized, file.path(processed_data_dir, "02_khandle_brfss_harmonized.rds"))
message("Saved: ", file.path(processed_data_dir, "02_khandle_brfss_harmonized.rds"))
message("02_harmonize_KHANDLE_BRFSS.R completed successfully.")
