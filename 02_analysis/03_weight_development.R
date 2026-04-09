#------------------------------------------------------------------
# Script: 03_develop_weights.R
# Purpose: Develop KHANDLE generalizability weights using harmonized KHANDLE+BRFSS
# Original author: Eleanor Hayes-Larson (modified by Tracy Lam-Hine)
# Repo: github.com/lamhine/KHANDLE-weighting-multiracial
# Last updated: 2026-04-08
# Notes:
#   - t10 + t16 compatible
#   - Uses readRDS/saveRDS (no load/save)
#   - Fits weight model within each imputation
#   - Produces diagnostics + balance (overall and by race)
#   - t16-specific: excludes english_interview_h (NA in KHANDLE)
#   - FIX: predict on full data within each imp so p_hat aligns with nrow(d)
#   - Helpers from R/utils.R, R/race_constants.R, R/balance_helpers.R via config
#------------------------------------------------------------------

source("config/config.R")

#------------------------------------------------------------------
# Load harmonized data
#------------------------------------------------------------------

dat <- readRDS(file.path(processed_data_dir, "02_khandle_brfss_harmonized.rds"))

#------------------------------------------------------------------
# Detect analysis_setup (t10 vs t16) from BRFSS years in harmonized
#------------------------------------------------------------------

if (!("imp_h" %in% names(dat))) stop("Missing imp_h in harmonized dataset.")
if (!("brfss_h" %in% names(dat))) stop("Missing brfss_h in harmonized dataset.")

analysis_setup <- detect_analysis_setup(dat)
message("03_develop_weights.R detected analysis_setup = ", analysis_setup)

#------------------------------------------------------------------
# Mission-critical checks: race coding + factor ordering
#------------------------------------------------------------------

if (!("race_summary_h" %in% names(dat))) stop("Missing race_summary_h in harmonized dataset.")
bad_codes <- setdiff(sort(unique(dat$race_summary_h)), c(1, 2, 3, 4, 5, NA))
bad_codes <- bad_codes[!is.na(bad_codes)]
if (length(bad_codes) > 0) stop("race_summary_h has unexpected codes: ", paste(bad_codes, collapse = ", "))

dat <- add_race_factor(dat)

assert_allowed_values(dat$race_summary_f, ALLOWED_RACE5, "race_summary_f")
print(with(dat, table(race_summary_f, brfss_h, useNA = "ifany")))

#------------------------------------------------------------------
# Define covariate set for balance plots
#------------------------------------------------------------------

covars_balance <- c(
  "age_h", "male_h",
  "education3_h",
  "income_gtmed_pp_h",
  "adl_walking_h", "adl_dressing_h",
  "english_interview_h",
  "goodhealth_h",
  "blind_h",
  "exercise_h",
  "marital_h",
  "military_h",
  "smoke_status_h"
)

if (analysis_setup == "t16") {
  covars_balance <- setdiff(covars_balance, "english_interview_h")
}

missing_covs <- setdiff(covars_balance, names(dat))
if (length(missing_covs) > 0) stop("Missing expected covariates: ", paste(missing_covs, collapse = ", "))

#------------------------------------------------------------------
# Weighting model specification
#
# NOTE on education balance: KHANDLE is ~50% college-educated (edu=4)
# vs ~16% in BRFSS. The logistic model correctly identifies this
# (OR ~9.9 for edu=4), but raw weights of 9-1,236 require aggressive
# trimming (1st/99th percentile). After trimmed IOSW, the overall
# education gap shrinks from 33.6pp to ~1.3pp. However, within-race
# SMDs remain large because KHANDLE composition differs dramatically
# by race (e.g., 69% of KHANDLE Asians are college-educated).
# See 03a_weight_diagnostics.R for full analysis.
#------------------------------------------------------------------

dat <- dat %>% mutate(khandle_h = 1 - brfss_h)

if (analysis_setup == "t10") {
  weight_formula <- khandle_h ~
    race_summary_f + male_h + age_h + education3_h +
    income_gtmed_pp_h + adl_walking_h + english_interview_h + goodhealth_h
} else {
  weight_formula <- khandle_h ~
    race_summary_f + male_h + age_h + education3_h +
    income_gtmed_pp_h + adl_walking_h + goodhealth_h
}

#------------------------------------------------------------------
# Prep common coercions
#------------------------------------------------------------------

dat <- dat %>%
  mutate(
    male_h = as01_num(male_h),
    english_interview_h = as01_num(english_interview_h),
    military_h = as01_num(military_h),
    blind_h = as01_num(blind_h),
    adl_walking_h = as01_num(adl_walking_h),
    adl_dressing_h = as01_num(adl_dressing_h),
    exercise_h = as01_num(exercise_h),
    goodhealth_h = as01_num(goodhealth_h),
    marital_h = as01_num(marital_h)
  )

if ("smoke_status_h" %in% names(dat)) {
  dat$smoke_status_h <- to_num(dat$smoke_status_h)
}

#------------------------------------------------------------------
# Step A: stacked diagnostics (optional; uses the pooled data)
#------------------------------------------------------------------

p_khandle_stack <- mean(dat$khandle_h, na.rm = TRUE)
message("p_khandle (stacked): ", signif(p_khandle_stack, 4))

covbal_before <- compute_covbal_by_race(
  dat %>% mutate(w_all = 1),
  vars = covars_balance,
  weightvar = "w_all"
)

covbal_overall_before <- compute_covbal_overall(
  dat %>% mutate(w_all = 1),
  vars = covars_balance,
  weightvar = "w_all"
)

fit_stack <- glm(
  formula = weight_formula,
  data = dat,
  family = binomial(link = "logit"),
  weights = dat$brfss_sampwt_h
)

p_hat_stack <- predict(fit_stack, type = "response")
print(glm_diag(fit_stack, p_hat_stack, label = "stacked"))

w_raw_stack <- (1 - p_hat_stack) / p_hat_stack
sw_stack <- w_raw_stack * (p_khandle_stack / (1 - p_khandle_stack))
sw_stack[dat$brfss_h == 1] <- 1

sw_stack_trim <- sw_stack
kh_idx_stack <- which(dat$brfss_h == 0 & !is.na(sw_stack))
if (length(kh_idx_stack) > 0) {
  sw_stack_trim[kh_idx_stack] <- trim_weights(sw_stack[kh_idx_stack], lower = 0.01, upper = 0.99)
}

covbal_after <- compute_covbal_by_race(
  dat %>% mutate(w_all = sw_stack_trim),
  vars = covars_balance,
  weightvar = "w_all"
)

covbal_overall_after <- compute_covbal_overall(
  dat %>% mutate(w_all = sw_stack_trim),
  vars = covars_balance,
  weightvar = "w_all"
)

#------------------------------------------------------------------
# Save covariate balance figures (stacked diagnostics)
#------------------------------------------------------------------

fig_dir <- file.path(results_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

p_cov_race <- plot_covbal_paper(
  before_df = covbal_before,
  after_df  = covbal_after,
  title = paste0("Covariate balance by race/ethnicity (", analysis_tag, ")"),
  y_lim = c(-1, 1.25)
)

p_cov_overall <- plot_covbal_overall_paper(
  before_df = covbal_overall_before,
  after_df  = covbal_overall_after,
  title = paste0("Overall covariate balance (", analysis_tag, ")"),
  y_lim = c(-1, 1)
)

if (!is.null(p_cov_race)) {
  save_fig(p_cov_race, paste0("03_covbal_by_race_", analysis_tag), width = 12, height = 6, fig_dir = fig_dir)
  print(p_cov_race)
}

if (!is.null(p_cov_overall)) {
  save_fig(p_cov_overall, paste0("03_covbal_overall_", analysis_tag), width = 12, height = 6, fig_dir = fig_dir)
  print(p_cov_overall)
}

message("Saved balance figures to: ", fig_dir)

overall_balance_summary <- summarize_balance_improvement(
  before_df = covbal_overall_before,
  after_df  = covbal_overall_after
)
print(overall_balance_summary$metrics)

#------------------------------------------------------------------
# Step B: fit weights within each imputation + store diagnostics
#   FIX: predict on full d, initialize sw as NA, compute only for ok_pred
#------------------------------------------------------------------

m <- sort(unique(dat$imp_h))
m <- m[!is.na(m)]
if (length(m) == 0) stop("No imp_h found.")

m_chr <- as.character(m)

weights_list      <- setNames(vector("list", length(m_chr)), m_chr)
diag_list         <- setNames(vector("list", length(m_chr)), m_chr)
covbal_unw        <- setNames(vector("list", length(m_chr)), m_chr)
covbal_wt         <- setNames(vector("list", length(m_chr)), m_chr)
covbal_unw_overall <- setNames(vector("list", length(m_chr)), m_chr)
covbal_wt_overall  <- setNames(vector("list", length(m_chr)), m_chr)

for (j in m) {

  key <- as.character(j)
  d <- dat %>% filter(imp_h == j)

  tab_j <- with(d, table(race_summary_f, brfss_h, useNA = "ifany"))

  p_khandle <- with(
    d,
    sum(khandle_h * brfss_sampwt_h, na.rm = TRUE) / sum(brfss_sampwt_h, na.rm = TRUE)
  )

  covbal_unw[[key]] <- compute_covbal_by_race(
    d %>% mutate(w_all = 1),
    vars = covars_balance,
    weightvar = "w_all"
  )

  covbal_unw_overall[[key]] <- compute_covbal_overall(
    d %>% mutate(w_all = 1),
    vars = covars_balance,
    weightvar = "w_all"
  )

  fit <- glm(
    formula = weight_formula,
    data = d,
    family = binomial(link = "logit"),
    weights = d$brfss_sampwt_h
  )

  # Predict on full d (keeps length = nrow(d); rows not used in fit become NA)
  p_hat <- predict(fit, newdata = d, type = "response")

  # Initialize weights as NA_real_ so missing predictions don't get bogus weights
  sw <- rep(NA_real_, nrow(d))

  # Only compute for rows with valid predictions
  ok_pred <- !is.na(p_hat) & p_hat > 0 & p_hat < 1
  sw[ok_pred] <- ((1 - p_hat[ok_pred]) / p_hat[ok_pred]) * (p_khandle / (1 - p_khandle))

  # By definition, BRFSS rows have generalizability weight 1
  sw[d$brfss_h == 1] <- 1

  # Trim KHANDLE only (and only where sw is not NA)
  sw_trim <- sw
  kh_idx <- which(d$brfss_h == 0 & !is.na(sw))
  if (length(kh_idx) > 0) {
    sw_trim[kh_idx] <- trim_weights(sw[kh_idx], lower = 0.01, upper = 0.99)
  }

  d <- d %>%
    mutate(
      p_hat = p_hat,
      sw_final = sw_trim
    )

  covbal_wt[[key]] <- compute_covbal_by_race(
    d %>% mutate(w_all = sw_final),
    vars = covars_balance,
    weightvar = "w_all"
  )

  covbal_wt_overall[[key]] <- compute_covbal_overall(
    d %>% mutate(w_all = sw_final),
    vars = covars_balance,
    weightvar = "w_all"
  )

  di <- glm_diag(fit, p_hat, label = paste0("imp_", j))
  wd <- weight_diag(d, w_var = "sw_final", group_var = "race_summary_f")

  diag_list[[key]] <- list(
    imp = j,
    race_source_tab = tab_j,
    glm = di,
    weight_overall = wd$overall,
    weight_by_race = wd$by_race
  )

  weights_list[[key]] <- d
}

#------------------------------------------------------------------
# Combine and save
#------------------------------------------------------------------

weights_all <- bind_rows(weights_list)

saveRDS(weights_all, file.path(processed_data_dir, "03_khandle_brfss_harmonized_weights.rds"))
saveRDS(diag_list,   file.path(processed_data_dir, "03_weight_model_diagnostics.rds"))
saveRDS(covbal_unw,  file.path(processed_data_dir, "03_covbal_unweighted_by_imp.rds"))
saveRDS(covbal_wt,   file.path(processed_data_dir, "03_covbal_weighted_by_imp.rds"))
saveRDS(covbal_unw_overall, file.path(processed_data_dir, "03_covbal_unweighted_overall_by_imp.rds"))
saveRDS(covbal_wt_overall,  file.path(processed_data_dir, "03_covbal_weighted_overall_by_imp.rds"))

message("Saved weighted stack + diagnostics to: ", processed_data_dir)
message("03_develop_weights.R completed successfully.")
