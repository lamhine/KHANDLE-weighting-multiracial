#------------------------------------------------------------------
# Script: 04_iosw_analysis.R
# Purpose: Estimate naive vs IOSW-weighted KHANDLE prevalence and disparities
# Repo: github.com/lamhine/KHANDLE-weighting-multiracial
# Last updated: 2026-04-08
# Notes:
#   - t10 + t16 compatible (auto-detect from BRFSS years in harmonized file)
#   - Uses RDS inputs from 03_develop_weights.R
#   - Uses sw_final created per imputation (no refitting weights here)
#   - Produces:
#       (1) Naive KHANDLE prevalence (overall, by race)
#       (2) Naive KHANDLE standardized to KHANDLE age/sex
#       (3) IOSW-weighted KHANDLE prevalence (overall, by race)
#       (4) IOSW-weighted KHANDLE standardized to BRFSS age/sex
#       (5) PR and PD vs White (standardized versions)
#   - Inference:
#       * Rubin's rules for (overall/by-race) naive + weighted means (approx SEs)
#       * Standardized estimates: point estimates MI-averaged; SE not computed here
#   - Helpers from R/utils.R, R/race_constants.R, R/analysis_helpers.R via config
#------------------------------------------------------------------

source("config/config.R")

#------------------------------------------------------------------
# Inputs
#------------------------------------------------------------------

harm_path <- file.path(processed_data_dir, "02_khandle_brfss_harmonized.rds")
wts_path  <- file.path(processed_data_dir, "03_khandle_brfss_harmonized_weights.rds")

dat_harm <- readRDS(harm_path)
dat_wts  <- readRDS(wts_path)

message("04_iosw_analysis.R using analysis_tag = ", analysis_tag)

dat_harm <- add_race_factor(dat_harm)
dat_wts  <- add_race_factor(dat_wts)

#------------------------------------------------------------------
# Basic checks
#------------------------------------------------------------------

req_vars_harm <- c("brfss_h", "imp_h", "age_h", "male_h", "race_summary_h", "brfss_sampwt_h")
req_vars_wts  <- c("brfss_h", "khandle_h", "imp_h", "age_h", "male_h", "race_summary_h", "sw_final")

assert_has_vars(dat_harm, req_vars_harm, "harmonized data")
assert_has_vars(dat_wts, req_vars_wts, "weighted stack")

# Outcome variable (t10 vs t16 can differ; set here)
outcome_var <- "cogimp_prob_fin_dx"
if (!(outcome_var %in% names(dat_wts))) stop("Missing outcome in dat_wts: ", outcome_var)

#------------------------------------------------------------------
# Age categories + standard populations
#------------------------------------------------------------------

dat_harm <- dat_harm %>% mutate(agecat_h = make_agecat(age_h))
dat_wts  <- dat_wts  %>% mutate(agecat_h = make_agecat(age_h))

std_pops <- compute_std_populations(dat_harm)
brfss_std <- std_pops$brfss_std
khandle_std <- std_pops$khandle_std

#------------------------------------------------------------------
# Core estimators per imputation
#------------------------------------------------------------------

races <- levels(dat_wts$race_summary_f)

estimate_one_imp <- function(d_wts_imp) {
  # Keep KHANDLE only for estimates
  k <- d_wts_imp %>% filter(brfss_h == 0)

  # 1) Naive KHANDLE overall + by race
  overall_unw <- umean_se(k[[outcome_var]])

  byrace_unw <- k %>%
    group_by(race_summary_f) %>%
    summarise(
      est = umean_se(.data[[outcome_var]])["est"],
      se  = umean_se(.data[[outcome_var]])["se"],
      .groups = "drop"
    )

  # 2) Naive standardized to KHANDLE age/sex
  cells_unw <- k %>%
    filter(!is.na(agecat_h), !is.na(male_h)) %>%
    group_by(race_summary_f, agecat_h, male_h) %>%
    summarise(
      est = umean_se(.data[[outcome_var]])["est"],
      .groups = "drop"
    )

  std_unw <- standardize_from_cells(cells_unw, khandle_std)

  # 3) IOSW-weighted KHANDLE overall + by race
  overall_wtd <- wmean_se(k[[outcome_var]], k$sw_final)

  byrace_wtd <- k %>%
    group_by(race_summary_f) %>%
    summarise(
      est = wmean_se(.data[[outcome_var]], sw_final)["est"],
      se  = wmean_se(.data[[outcome_var]], sw_final)["se"],
      .groups = "drop"
    )

  # 4) IOSW-weighted standardized to BRFSS age/sex
  cells_wtd <- k %>%
    filter(!is.na(agecat_h), !is.na(male_h)) %>%
    group_by(race_summary_f, agecat_h, male_h) %>%
    summarise(
      est = wmean_se(.data[[outcome_var]], sw_final)["est"],
      .groups = "drop"
    )

  std_wtd <- standardize_from_cells(cells_wtd, brfss_std)

  # PR/PD vs White (standardized)
  prpd_unw <- std_unw %>%
    rename(est_unw_std = est_std) %>%
    left_join(
      std_unw %>% filter(race_summary_f == REF_RACE) %>% transmute(ref_unw_std = est_std),
      by = character()
    ) %>%
    mutate(
      PR_unw_std = est_unw_std / ref_unw_std,
      PD_unw_std = est_unw_std - ref_unw_std
    )

  prpd_wtd <- std_wtd %>%
    rename(est_wtd_std = est_std) %>%
    left_join(
      std_wtd %>% filter(race_summary_f == REF_RACE) %>% transmute(ref_wtd_std = est_std),
      by = character()
    ) %>%
    mutate(
      PR_wtd_std = est_wtd_std / ref_wtd_std,
      PD_wtd_std = est_wtd_std - ref_wtd_std
    )

  list(
    overall_unw = overall_unw,
    overall_wtd = overall_wtd,
    byrace_unw  = byrace_unw,
    byrace_wtd  = byrace_wtd,
    std_unw     = std_unw,
    std_wtd     = std_wtd,
    prpd_unw    = prpd_unw,
    prpd_wtd    = prpd_wtd
  )
}

#------------------------------------------------------------------
# Run per-imputation estimates
#------------------------------------------------------------------

imps <- sort(unique(dat_wts$imp_h))
imps <- imps[!is.na(imps)]
if (length(imps) == 0) stop("No imp_h found in dat_wts.")

res_by_imp <- vector("list", length(imps))
names(res_by_imp) <- as.character(imps)

for (j in imps) {
  d_wts_imp <- dat_wts %>% filter(imp_h == j)

  if (nrow(d_wts_imp) == 0) next
  if (!("sw_final" %in% names(d_wts_imp))) stop("Missing sw_final in dat_wts imp=", j)

  res_by_imp[[as.character(j)]] <- estimate_one_imp(d_wts_imp)
}

#------------------------------------------------------------------
# Combine across imputations
#------------------------------------------------------------------

# Overall
overall_unw_mi <- rubin_combine(
  est = sapply(res_by_imp, function(x) x$overall_unw[["est"]]),
  se  = sapply(res_by_imp, function(x) x$overall_unw[["se"]])
) %>% mutate(stat = "Overall_unweighted")

overall_wtd_mi <- rubin_combine(
  est = sapply(res_by_imp, function(x) x$overall_wtd[["est"]]),
  se  = sapply(res_by_imp, function(x) x$overall_wtd[["se"]])
) %>% mutate(stat = "Overall_weighted")

results_overall <- bind_rows(overall_unw_mi, overall_wtd_mi)

# By race (unweighted + weighted)
combine_byrace <- function(component = c("byrace_unw", "byrace_wtd")) {
  component <- match.arg(component)

  out <- lapply(races, function(r) {
    est <- sapply(res_by_imp, function(x) {
      tmp <- x[[component]] %>% filter(race_summary_f == r)
      if (nrow(tmp) == 0) return(NA_real_)
      tmp$est[1]
    })
    se <- sapply(res_by_imp, function(x) {
      tmp <- x[[component]] %>% filter(race_summary_f == r)
      if (nrow(tmp) == 0) return(NA_real_)
      tmp$se[1]
    })
    rubin_combine(est, se) %>% mutate(race = r)
  })

  bind_rows(out)
}

byrace_unw_mi <- combine_byrace("byrace_unw") %>%
  mutate(stat = "By_race_unweighted") %>%
  select(stat, race, est, se, m)

byrace_wtd_mi <- combine_byrace("byrace_wtd") %>%
  mutate(stat = "By_race_weighted") %>%
  select(stat, race, est, se, m)

results_byrace <- bind_rows(byrace_unw_mi, byrace_wtd_mi)

# Standardized prevalence (unweighted to KHANDLE; weighted to BRFSS)
combine_std <- function(component = c("std_unw", "std_wtd"), stat_name) {
  component <- match.arg(component)

  out <- lapply(races, function(r) {
    est <- sapply(res_by_imp, function(x) {
      tmp <- x[[component]] %>% filter(race_summary_f == r)
      if (nrow(tmp) == 0) return(NA_real_)
      tmp$est_std[1]
    })
    tibble(
      m = sum(is.finite(est)),
      est = mean(est, na.rm = TRUE),
      se = NA_real_,
      race = r,
      stat = stat_name
    )
  })

  bind_rows(out)
}

std_unw_mi <- combine_std("std_unw", stat_name = "Std_to_KHANDLE_unweighted")
std_wtd_mi <- combine_std("std_wtd", stat_name = "Std_to_BRFSS_weighted")

results_std <- bind_rows(
  std_unw_mi %>% select(stat, race, est, m),
  std_wtd_mi %>% select(stat, race, est, m)
)

# PR/PD vs White (standardized; MI-averaged)
combine_prpd <- function(component = c("prpd_unw", "prpd_wtd"),
                         pr_col, pd_col,
                         stat_prefix) {
  component <- match.arg(component)

  out <- lapply(races, function(r) {
    pr <- sapply(res_by_imp, function(x) {
      tmp <- x[[component]] %>% filter(race_summary_f == r)
      if (nrow(tmp) == 0) return(NA_real_)
      tmp[[pr_col]][1]
    })
    pd <- sapply(res_by_imp, function(x) {
      tmp <- x[[component]] %>% filter(race_summary_f == r)
      if (nrow(tmp) == 0) return(NA_real_)
      tmp[[pd_col]][1]
    })

    tibble(
      race = r,
      PR = mean(pr, na.rm = TRUE),
      PD = mean(pd, na.rm = TRUE),
      m_PR = sum(is.finite(pr)),
      m_PD = sum(is.finite(pd)),
      stat = stat_prefix
    )
  })

  bind_rows(out)
}

prpd_unw_mi <- combine_prpd(
  "prpd_unw",
  pr_col = "PR_unw_std",
  pd_col = "PD_unw_std",
  stat_prefix = "PR_PD_unweighted_std"
)

prpd_wtd_mi <- combine_prpd(
  "prpd_wtd",
  pr_col = "PR_wtd_std",
  pd_col = "PD_wtd_std",
  stat_prefix = "PR_PD_weighted_std"
)

results_prpd <- bind_rows(prpd_unw_mi, prpd_wtd_mi) %>%
  mutate(is_reference = (race == REF_RACE))

#------------------------------------------------------------------
# Track-specific output names (avoid overwriting between t10 and t16)
#------------------------------------------------------------------

out_tag <- analysis_tag

saveRDS(results_overall, file.path(processed_data_dir, paste0("04_iosw_results_overall_", out_tag, ".rds")))
saveRDS(results_byrace,  file.path(processed_data_dir, paste0("04_iosw_results_by_race_", out_tag, ".rds")))
saveRDS(results_std,     file.path(processed_data_dir, paste0("04_iosw_results_standardized_", out_tag, ".rds")))
saveRDS(results_prpd,    file.path(processed_data_dir, paste0("04_iosw_results_pr_pd_", out_tag, ".rds")))

message("Saved IOSW analysis outputs to: ", processed_data_dir)
message("Outputs tagged with: ", out_tag)
message("04_iosw_analysis.R completed successfully.")
