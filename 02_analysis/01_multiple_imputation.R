#------------------------------------------------------------------
# Script: 01_multiple_imputation.R
# Purpose: Multiple imputation for KHANDLE + BRFSS (t10 + t16 compatible)
# Repo: github.com/lamhine/KHANDLE-weighting-multiracial
# Last updated: 2026-01-13
#------------------------------------------------------------------

source("config/config.R")

#------------------------------------------------------------------
# Helpers
# to_num(), get_col(), assert_has_vars() provided by R/utils.R via config
#------------------------------------------------------------------

print_table <- function(x, name, max_levels = 30) {
  cat("\n---", name, "---\n")
  ux <- unique(x)
  if (length(ux) > max_levels) {
    cat("Too many levels (", length(ux), "). Showing head of values + NA count.\n", sep = "")
    print(head(sort(ux), 20))
    cat("NA count:", sum(is.na(x)), "\n")
  } else {
    print(table(x, useNA = "ifany"))
  }
}

pct_missing <- function(x) 100 * sum(is.na(x)) / length(x)

preflight_report <- function(dat, vars, dat_name) {
  cat("\n====================================================\n")
  cat("Preflight:", dat_name, "\n")
  cat("N rows:", nrow(dat), "\n")
  cat("====================================================\n")
  
  miss <- data.frame(
    var = vars,
    pct_miss = sapply(vars, function(v) pct_missing(dat[[v]]))
  )
  miss <- miss[order(miss$pct_miss), ]
  print(miss)
  
  cat("\nValue tables (including NA):\n")
  for (v in vars) {
    print_table(dat[[v]], v)
  }
  
  invisible(miss)
}

disable_if_single_level <- function(dat, vars, method) {
  for (v in vars) {
    if (v %in% names(dat) && !identical(method[[v]], "")) {
      x <- dat[[v]]
      if (is.factor(x)) {
        obs <- droplevels(x[!is.na(x)])
        if (length(obs) == 0 || nlevels(obs) < 2) {
          method[[v]] <- ""
        }
      }
    }
  }
  method
}

#------------------------------------------------------------------
# Global MI settings
#------------------------------------------------------------------

m_imps <- 20
maxit  <- 10
seed   <- 20260112
default_method <- c("pmm", "pmm", "pmm", "pmm")

#------------------------------------------------------------------
# Load cleaned inputs
#------------------------------------------------------------------

khandle_clean <- readRDS(file.path(processed_data_dir, "03_khandle_clean_full.rds"))
brfss_clean   <- readRDS(file.path(processed_data_dir, "02_brfss_clean.rds"))

#------------------------------------------------------------------
# Detect analysis setup (t10 vs t16) from BRFSS years
#------------------------------------------------------------------

analysis_setup <- detect_analysis_setup(brfss_clean)
message("01_multiple_imputation.R detected analysis_setup = ", analysis_setup)

#------------------------------------------------------------------
# KHANDLE: variables to impute (harmonizable only, minimal model)
#------------------------------------------------------------------

k_required <- c(
  "studyid",
  "age", "male", "race_mece_brfss5",
  "marital_status", "education", "income_range", "military",
  "health", "smoke_status"
)

k_optional <- c(
  "sensimp_vision",
  "adl1", "adl2", "adl9",
  "pa_lt_ex", "pa_vig_ex",
  "w1_cogimp_prob_fin_dx", "english_interview"
)

assert_has_vars(khandle_clean, k_required, "KHANDLE clean full")

k_in <- khandle_clean %>%
  transmute(
    studyid = as.character(studyid),
    
    w1_interview_age = to_num(age),
    male = to_num(male),
    race_cat = as.character(race_mece_brfss5),
    
    marital_status = to_num(marital_status),
    education = to_num(education),
    income_range = to_num(income_range),
    military = to_num(military),
    
    sensimp_vision = to_num(get_col(., "sensimp_vision")),
    english_interview = to_num(get_col(., "english_interview")),
    health = to_num(health),
    smoke_status = to_num(smoke_status),
    
    adl1 = to_num(get_col(., "adl1")),
    adl2 = to_num(get_col(., "adl2")),
    adl9 = to_num(get_col(., "adl9")),
    
    pa_lt_ex  = to_num(get_col(., "pa_lt_ex")),
    pa_vig_ex = to_num(get_col(., "pa_vig_ex")),
    
    w1_cogimp_prob_fin_dx = to_num(get_col(., "w1_cogimp_prob_fin_dx"))
  )

k_in <- k_in %>%
  mutate(
    marital_status = dplyr::na_if(marital_status, 77),
    sensimp_vision = dplyr::na_if(sensimp_vision, 77),
    adl1 = dplyr::na_if(adl1, 77),
    adl2 = dplyr::na_if(adl2, 77),
    adl9 = dplyr::na_if(adl9, 77)
  )

all_na_cols <- names(k_in)[vapply(k_in, function(x) all(is.na(x)), logical(1))]
if (length(all_na_cols) > 0) {
  message("Dropping all-NA KHANDLE columns: ", paste(all_na_cols, collapse = ", "))
  k_in <- k_in[, setdiff(names(k_in), all_na_cols), drop = FALSE]
}

bad_smk <- setdiff(
  sort(unique(k_in$smoke_status[!is.na(k_in$smoke_status)])),
  c(0, 1, 2)
)
if (length(bad_smk) > 0) {
  stop("KHANDLE smoke_status has unexpected values: ", paste(bad_smk, collapse = ", "))
}

k_preflight_vars <- names(k_in)
preflight_report(k_in, k_preflight_vars, "KHANDLE in.data")

k_continuous <- c("w1_interview_age", "w1_cogimp_prob_fin_dx")
k_binary     <- intersect(c("male", "english_interview", "military"), names(k_in))
k_ordinal    <- c(
  "education", "income_range", "sensimp_vision", "health",
  "pa_lt_ex", "pa_vig_ex", "adl1", "adl2", "adl9",
  "smoke_status"
)
k_categorical <- c("race_cat", "marital_status")

drop_all_na <- function(vars, dat) {
  vars[vars %in% names(dat) & !sapply(vars, function(v) all(is.na(dat[[v]])))]
}
k_continuous  <- drop_all_na(k_continuous, k_in)
k_binary      <- drop_all_na(k_binary, k_in)
k_ordinal     <- drop_all_na(k_ordinal, k_in)
k_categorical <- drop_all_na(k_categorical, k_in)

all_k_typed <- c(k_continuous, k_binary, k_ordinal, k_categorical, "studyid")
stopifnot(length(setdiff(names(k_in), all_k_typed)) == 0)

for (v in k_binary)      k_in[[v]] <- factor(k_in[[v]])
for (v in k_ordinal)     k_in[[v]] <- factor(k_in[[v]], ordered = TRUE)
for (v in k_categorical) k_in[[v]] <- factor(k_in[[v]], ordered = FALSE)

set.seed(seed)
k_ini <- mice(k_in, maxit = 0, defaultMethod = default_method, seed = seed)

if (!is.null(k_ini$loggedEvents) && nrow(k_ini$loggedEvents) > 0) {
  message("\nKHANDLE mice init loggedEvents:")
  print(k_ini$loggedEvents)
}

k_meth <- k_ini$method
k_pred <- k_ini$predictorMatrix

k_meth["studyid"] <- ""
k_pred[, "studyid"] <- 0

k_meth[k_continuous]  <- "pmm"
k_meth[k_binary]      <- "logreg"
k_meth[k_ordinal]     <- "polr"
k_meth[k_categorical] <- "polyreg"

k_meth["race_cat"] <- ""
k_pred["race_cat", ] <- 0

vars_to_check <- c(k_continuous, k_binary, k_ordinal, k_categorical)
k_meth <- disable_if_single_level(k_in, vars_to_check, k_meth)

set.seed(seed)
k_mids <- mice(
  k_in,
  m = m_imps,
  maxit = maxit,
  predictorMatrix = k_pred,
  method = k_meth,
  defaultMethod = default_method,
  seed = seed,
  printFlag = TRUE
)

khandle_pmm <- complete(k_mids, action = "long", include = FALSE) %>%
  as_tibble() %>%
  rename(imp = .imp)

#------------------------------------------------------------------
# BRFSS: variables to impute (harmonizable only, minimal model)
# Decision: DO IMPUTE hhsize_v2 for t16 (do not reconstruct).
#------------------------------------------------------------------

if ("hhsize_std" %in% names(brfss_clean) && !("hhsize_v2" %in% names(brfss_clean))) {
  brfss_clean$hhsize_v2 <- brfss_clean$hhsize_std
}

b_harmonizable <- c(
  "year", "age", "male", "raceth",
  "marital_status", "education_ca", "income", "military",
  "vision_impair", "english_interview", "health", "smoke_status",
  "diff_walk", "diff_dress", "phys_activity", "hhsize_v2",
  "county_code", "cdc_finalwt"
)

assert_has_vars(brfss_clean, b_harmonizable, "BRFSS clean")

b_in <- brfss_clean %>%
  transmute(
    year = to_num(year),
    
    age = to_num(age),
    male = to_num(male),
    race_cat = as.character(raceth),
    
    marital_status = to_num(marital_status),
    education = to_num(education_ca),
    income = to_num(income),
    military = to_num(military),
    
    vision_impair = to_num(vision_impair),
    english_interview = to_num(english_interview),
    health = to_num(health),
    smoke_status = to_num(smoke_status),
    
    diff_walk = to_num(diff_walk),
    diff_dress = to_num(diff_dress),
    phys_activity = to_num(phys_activity),
    
    hhsize_v2 = to_num(hhsize_v2),
    
    county_code = county_code,
    cdc_finalwt = cdc_finalwt
  ) %>%
  mutate(
    # guardrails before imputation
    hhsize_v2 = dplyr::case_when(
      is.na(hhsize_v2) ~ NA_real_,
      hhsize_v2 <= 0 ~ NA_real_,
      hhsize_v2 > 20 ~ 20,
      TRUE ~ hhsize_v2
    )
  )

b_preflight_vars <- names(b_in)
preflight_report(b_in, b_preflight_vars, "BRFSS in.data")

# Types
b_continuous <- c("age", "cdc_finalwt")

b_binary <- c(
  "male", "english_interview", "military", "vision_impair",
  "diff_walk", "diff_dress", "phys_activity"
)

# Treat hhsize as ordered categorical so it imputes cleanly and stays plausible
# (you can widen levels later if needed)
b_ordinal <- c("education", "income", "health", "smoke_status", "hhsize_v2")

b_categorical <- c("race_cat", "marital_status", "county_code")  # year is predictor-only

all_b_typed <- c(b_continuous, b_binary, b_ordinal, b_categorical, "year")
stopifnot(length(setdiff(names(b_in), all_b_typed)) == 0)

for (v in b_binary)      b_in[[v]] <- factor(b_in[[v]])
for (v in b_ordinal)     b_in[[v]] <- factor(b_in[[v]], ordered = TRUE)
for (v in b_categorical) b_in[[v]] <- factor(b_in[[v]], ordered = FALSE)
b_in[["year"]] <- factor(b_in[["year"]], ordered = FALSE)

set.seed(seed)
b_ini <- mice(b_in, maxit = 0, defaultMethod = default_method, seed = seed)

b_meth <- b_ini$method
b_pred <- b_ini$predictorMatrix

b_meth[b_continuous]  <- "pmm"
b_meth[b_binary]      <- "logreg"
b_meth[b_ordinal]     <- "polr"
b_meth[b_categorical] <- "polyreg"

# Predictor-only (never impute), but keep as predictors:
# - year should predict others
# - race_cat should predict others
# Never impute:
no_impute_b <- c("race_cat", "county_code", "cdc_finalwt", "year")

# If you still want english_interview fixed for t16, keep this line.
# Otherwise, remove it.
if (analysis_setup == "t16") {
  no_impute_b <- c(no_impute_b, "english_interview")
}

for (v in no_impute_b) b_meth[v] <- ""
for (v in no_impute_b) b_pred[v, ] <- 0

# Block county_code and cdc_finalwt as predictors (keeps things stable)
b_pred[, "county_code"] <- 0
b_pred[, "cdc_finalwt"] <- 0

# race_cat predictor-only (keep column), do not predict it
b_pred["race_cat", ] <- 0

# year predictor-only (keep column), do not predict it
b_pred["year", ] <- 0

vars_to_check <- c(b_continuous, b_binary, b_ordinal, b_categorical, "year")
b_meth <- disable_if_single_level(b_in, vars_to_check, b_meth)

set.seed(seed)
b_mids <- mice(
  b_in,
  m = m_imps,
  maxit = maxit,
  predictorMatrix = b_pred,
  method = b_meth,
  defaultMethod = default_method,
  seed = seed,
  printFlag = TRUE
)

brfss_pmm <- complete(b_mids, action = "long", include = FALSE) %>%
  as_tibble() %>%
  rename(imp = .imp)

#------------------------------------------------------------------
# Save artifacts
#------------------------------------------------------------------

saveRDS(khandle_pmm, file.path(processed_data_dir, "khandle_pmm.rds"))
saveRDS(brfss_pmm,   file.path(processed_data_dir, "brfss_pmm.rds"))

saveRDS(k_mids, file.path(processed_data_dir, "khandle_mids.rds"))
saveRDS(b_mids, file.path(processed_data_dir, "brfss_mids.rds"))

message("Saved: ", file.path(processed_data_dir, "khandle_pmm.rds"))
message("Saved: ", file.path(processed_data_dir, "brfss_pmm.rds"))
message("Saved: ", file.path(processed_data_dir, "khandle_mids.rds"))
message("Saved: ", file.path(processed_data_dir, "brfss_mids.rds"))
message("01_multiple_imputation.R completed successfully.")