#------------------------------------------------------------------
# Script: 05_gcomp_analysis.R
# Purpose: G-computation (outcome modeling) to transport KHANDLE outcome model to BRFSS
# Repo: github.com/lamhine/KHANDLE-weighting-multiracial
# Last updated: 2026-04-08
# Notes:
#   - t10 + t16 compatible (auto-detect from BRFSS years in harmonized file)
#   - Uses readRDS/saveRDS only
#   - Outcome model fit in KHANDLE; predicted into BRFSS
#   - Includes Multiracial category
#   - Produces (MI-averaged point estimates):
#       * prediction diagnostics
#       * crude prevalence: KHANDLE observed; BRFSS OM-predicted (weighted)
#       * PR/PD vs White standardized to:
#           - all covariates in OM ("std_all")
#           - age/sex only ("std_agesex")
#   - Outputs tagged with analysis_tag to avoid overwriting across tracks
#   - Helpers from R/utils.R, R/race_constants.R, R/analysis_helpers.R via config
#------------------------------------------------------------------

source("config/config.R")

message("05_gcomp_analysis.R using analysis_tag = ", analysis_tag)

#------------------------------------------------------------------
# Inputs
#------------------------------------------------------------------

harm_path <- file.path(processed_data_dir, "02_khandle_brfss_harmonized.rds")
dat <- readRDS(harm_path)

#------------------------------------------------------------------
# Setup
#------------------------------------------------------------------

analysis_setup <- detect_analysis_setup(dat)
message("05_gcomp_analysis.R detected analysis_setup = ", analysis_setup)

# Outcome variable in KHANDLE
outcome_var <- "cogimp_prob_fin_dx"
if (!(outcome_var %in% names(dat))) stop("Missing outcome: ", outcome_var)

dat <- add_race_factor(dat)
dat <- dat %>% mutate(agecat_h = make_agecat(age_h))

#------------------------------------------------------------------
# Standard populations (age/sex)
#------------------------------------------------------------------

if (!("brfss_sampwt_h" %in% names(dat))) stop("Missing brfss_sampwt_h in harmonized data.")

std_pops <- compute_std_populations(dat)
brfss_std <- std_pops$brfss_std
khandle_std <- std_pops$khandle_std

#------------------------------------------------------------------
# Outcome model specification
#   - Fit in KHANDLE only
#   - t16-specific: exclude english_interview_h terms (NA in KHANDLE)
#   - Linear probability model (lm), matching original approach
#------------------------------------------------------------------

if (!("english_interview_h" %in% names(dat))) {
  analysis_setup <- "t16"
  message("05_gcomp_analysis.R: english_interview_h not found; forcing analysis_setup = t16 for OM.")
}

if (analysis_setup == "t10") {
  om_formula <- stats::as.formula(
    paste0(
      outcome_var,
      " ~ race_summary_f + male_h + age_h + education3_h + income_gtmed_pp_h + ",
      "adl_walking_h + english_interview_h + goodhealth_h + ",
      "race_summary_f:male_h + race_summary_f:age_h + race_summary_f:education3_h + ",
      "race_summary_f:income_gtmed_pp_h + race_summary_f:adl_walking_h + ",
      "race_summary_f:goodhealth_h"
    )
  )
} else {
  om_formula <- stats::as.formula(
    paste0(
      outcome_var,
      " ~ race_summary_f + male_h + age_h + education3_h + income_gtmed_pp_h + ",
      "adl_walking_h + goodhealth_h + ",
      "race_summary_f:male_h + race_summary_f:age_h + race_summary_f:education3_h + ",
      "race_summary_f:income_gtmed_pp_h + race_summary_f:adl_walking_h + ",
      "race_summary_f:goodhealth_h"
    )
  )
}

#' Weighted mean (safe, for internal use)
weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) == 0) return(NA_real_)
  sum(w[ok] * x[ok]) / sum(w[ok])
}

#------------------------------------------------------------------
# One-imputation estimation
#------------------------------------------------------------------

estimate_one_imp <- function(d_imp) {
  # Split
  k <- d_imp %>% filter(brfss_h == 0)
  b <- d_imp %>% filter(brfss_h == 1)

  # Fit OM in KHANDLE only
  fit <- stats::lm(om_formula, data = k)

  # Predict for everyone (KHANDLE + BRFSS)
  pred_all <- stats::predict(fit, newdata = d_imp)
  d_imp <- d_imp %>% mutate(pred_p = pred_all)

  # Prediction diagnostics (LPM can go <0 or >1)
  pred_diag <- d_imp %>%
    group_by(brfss_h) %>%
    summarise(
      n = n(),
      prop_lt0 = mean(pred_p < 0, na.rm = TRUE),
      prop_gt1 = mean(pred_p > 1, na.rm = TRUE),
      .groups = "drop"
    )

  # "Crude" prevalence
  crude_overall <- bind_rows(
    b %>%
      summarise(
        sample = "BRFSS",
        group = "Overall",
        est = weighted_mean_safe(stats::predict(fit, newdata = b), b$brfss_sampwt_h)
      ),
    k %>%
      summarise(
        sample = "KHANDLE",
        group = "Overall",
        est = mean(.data[[outcome_var]], na.rm = TRUE)
      )
  )

  crude_byrace <- bind_rows(
    b %>%
      mutate(pred_p = stats::predict(fit, newdata = b)) %>%
      group_by(race_summary_f) %>%
      summarise(
        sample = "BRFSS",
        group = as.character(race_summary_f),
        est = weighted_mean_safe(pred_p, brfss_sampwt_h),
        .groups = "drop"
      ),
    k %>%
      group_by(race_summary_f) %>%
      summarise(
        sample = "KHANDLE",
        group = as.character(race_summary_f),
        est = mean(.data[[outcome_var]], na.rm = TRUE),
        .groups = "drop"
      )
  )

  # Standardize to ALL covariates in OM via replicate-everyone-to-each-race
  base_rows <- d_imp %>% select(-race_summary_f)

  base_rep <- base_rows[rep(seq_len(nrow(base_rows)), times = length(ALLOWED_RACE5)), , drop = FALSE]
  race_rep <- rep(ALLOWED_RACE5, each = nrow(base_rows))

  reps <- base_rep %>%
    mutate(race_summary_f = factor(race_rep, levels = ALLOWED_RACE5))

  reps <- reps %>%
    mutate(pred_p = stats::predict(fit, newdata = reps))

  reps <- reps %>%
    mutate(w_std = if_else(brfss_h == 1, brfss_sampwt_h, 1))

  std_all <- reps %>%
    group_by(brfss_h, race_summary_f) %>%
    summarise(
      est = weighted_mean_safe(pred_p, w_std),
      .groups = "drop"
    ) %>%
    mutate(
      sample = if_else(brfss_h == 1, "BRFSS", "KHANDLE"),
      group = as.character(race_summary_f)
    ) %>%
    select(sample, group, est)

  std_all_PRPD <- std_all %>%
    group_by(sample) %>%
    mutate(
      ref = est[group == REF_RACE][1],
      PR_std_all = est / ref,
      PD_std_all = est - ref
    ) %>%
    ungroup() %>%
    select(sample, group, PR_std_all, PD_std_all)

  # Standardize to age/sex only
  std_agesex_k <- k %>%
    filter(!is.na(agecat_h), !is.na(male_h)) %>%
    group_by(race_summary_f, agecat_h, male_h) %>%
    summarise(est_cell = mean(.data[[outcome_var]], na.rm = TRUE), .groups = "drop") %>%
    left_join(khandle_std, by = c("agecat_h", "male_h")) %>%
    group_by(race_summary_f) %>%
    summarise(est = sum(est_cell * prop_agesex, na.rm = TRUE), .groups = "drop") %>%
    transmute(sample = "KHANDLE", group = as.character(race_summary_f), est)

  std_agesex_b <- b %>%
    mutate(pred_p = stats::predict(fit, newdata = b)) %>%
    filter(!is.na(agecat_h), !is.na(male_h)) %>%
    group_by(race_summary_f, agecat_h, male_h) %>%
    summarise(est_cell = mean(pred_p, na.rm = TRUE), .groups = "drop") %>%
    left_join(brfss_std, by = c("agecat_h", "male_h")) %>%
    group_by(race_summary_f) %>%
    summarise(est = sum(est_cell * prop_agesex, na.rm = TRUE), .groups = "drop") %>%
    transmute(sample = "BRFSS", group = as.character(race_summary_f), est)

  std_agesex <- bind_rows(std_agesex_k, std_agesex_b)

  std_agesex_PRPD <- std_agesex %>%
    group_by(sample) %>%
    mutate(
      ref = est[group == REF_RACE][1],
      PR_std_agesex = est / ref,
      PD_std_agesex = est - ref
    ) %>%
    ungroup() %>%
    select(sample, group, PR_std_agesex, PD_std_agesex)

  list(
    pred_diag = pred_diag,
    crude_overall = crude_overall,
    crude_byrace = crude_byrace,
    std_all_PRPD = std_all_PRPD,
    std_agesex_PRPD = std_agesex_PRPD
  )
}

#------------------------------------------------------------------
# Run over imputations
#------------------------------------------------------------------

imps <- sort(unique(dat$imp_h))
imps <- imps[!is.na(imps)]
if (length(imps) == 0) stop("No imp_h found.")

res_by_imp <- vector("list", length(imps))
names(res_by_imp) <- as.character(imps)

for (j in imps) {
  d_imp <- dat %>% filter(imp_h == j)
  res_by_imp[[as.character(j)]] <- estimate_one_imp(d_imp)
}

#------------------------------------------------------------------
# Combine across imputations (simple averaging of point estimates)
#------------------------------------------------------------------

bind_component <- function(name) {
  out <- lapply(names(res_by_imp), function(k) {
    res_by_imp[[k]][[name]] %>% mutate(imp_h = as.integer(k))
  })
  bind_rows(out)
}

pred_diag_all      <- bind_component("pred_diag")
crude_overall_all  <- bind_component("crude_overall")
crude_byrace_all   <- bind_component("crude_byrace")
std_all_all        <- bind_component("std_all_PRPD")
std_agesex_all     <- bind_component("std_agesex_PRPD")

avg_by <- function(df, group_vars, value_vars) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(across(all_of(value_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

results_pred_diag <- avg_by(pred_diag_all, c("brfss_h"), c("prop_lt0", "prop_gt1"))

results_crude_overall <- avg_by(crude_overall_all, c("sample", "group"), c("est"))
results_crude_byrace  <- avg_by(crude_byrace_all,  c("sample", "group"), c("est"))

results_std_all_PRPD <- avg_by(std_all_all, c("sample", "group"), c("PR_std_all", "PD_std_all"))
results_std_agesex_PRPD <- avg_by(std_agesex_all, c("sample", "group"), c("PR_std_agesex", "PD_std_agesex"))

#------------------------------------------------------------------
# Save outputs (RDS) -- tagged to avoid overwriting across tracks
#------------------------------------------------------------------

out_tag <- analysis_tag

saveRDS(results_pred_diag, file.path(processed_data_dir, paste0("05_gcomp_pred_diag_", out_tag, ".rds")))
saveRDS(results_crude_overall, file.path(processed_data_dir, paste0("05_gcomp_crude_overall_", out_tag, ".rds")))
saveRDS(results_crude_byrace,  file.path(processed_data_dir, paste0("05_gcomp_crude_byrace_", out_tag, ".rds")))
saveRDS(results_std_all_PRPD,  file.path(processed_data_dir, paste0("05_gcomp_std_all_PRPD_", out_tag, ".rds")))
saveRDS(results_std_agesex_PRPD, file.path(processed_data_dir, paste0("05_gcomp_std_agesex_PRPD_", out_tag, ".rds")))

message("Saved g-computation results to: ", processed_data_dir)
message("Outputs tagged with: ", out_tag)
message("05_gcomp_analysis.R completed successfully.")
