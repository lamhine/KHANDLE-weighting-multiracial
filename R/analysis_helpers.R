#------------------------------------------------------------------
# R/analysis_helpers.R
# Purpose: Shared analysis/estimation functions for the pipeline
#------------------------------------------------------------------

#' Detect analysis setup (t10 vs t16) from BRFSS years in a data frame
#' @param df Data frame with brfss_h and year columns
#' @return Character string "t10" or "t16"
detect_analysis_setup <- function(df) {
  if (!("brfss_h" %in% names(df))) {
    # If no brfss_h column, try to detect from year alone
    if (!("year" %in% names(df))) {
      message("detect_analysis_setup: no brfss_h or year found; defaulting to t16")
      return("t16")
    }
    yr_min <- suppressWarnings(min(to_num(df$year), na.rm = TRUE))
  } else {
    if (!("year" %in% names(df))) {
      message("detect_analysis_setup: year not found; defaulting to t16")
      return("t16")
    }
    yr_min <- suppressWarnings(min(to_num(df$year[df$brfss_h == 1]), na.rm = TRUE))
  }
  ifelse(is.finite(yr_min) && yr_min <= 2018, "t10", "t16")
}

#' Weighted mean + SE (Horvitz-Thompson-style approximation)
#'
#' This treats rows as approximately iid, which is appropriate for KHANDLE
#' IOSW estimates since KHANDLE is a cohort study (not a complex survey).
#' The BRFSS survey design affects weight estimation, but that uncertainty
#' is captured by Rubin's rules across multiple imputations (weights are
#' re-estimated per imputation). Within-imputation SEs here approximate
#' the sampling variability of the weighted KHANDLE mean.
wmean_se <- function(y, w) {
  ok <- is.finite(y) & is.finite(w) & w > 0
  y <- y[ok]
  w <- w[ok]
  if (length(y) < 2) return(c(est = NA_real_, se = NA_real_))
  est <- sum(w * y) / sum(w)
  se  <- sqrt(sum((w^2) * (y - est)^2) / (sum(w)^2))
  c(est = est, se = se)
}

#' Unweighted mean + SE
umean_se <- function(y) {
  y <- y[is.finite(y)]
  if (length(y) < 2) return(c(est = NA_real_, se = NA_real_))
  est <- mean(y)
  se  <- stats::sd(y) / sqrt(length(y))
  c(est = est, se = se)
}

#' Combine estimates across imputations using Rubin's rules
rubin_combine <- function(est, se) {
  ok <- is.finite(est) & is.finite(se)
  est <- est[ok]
  se  <- se[ok]
  m <- length(est)
  if (m == 0) return(tibble::tibble(m = 0, est = NA_real_, se = NA_real_))

  qbar <- mean(est)
  ubar <- mean(se^2)
  b    <- stats::var(est)
  if (!is.finite(b)) b <- 0
  tvar <- ubar + (1 + 1/m) * b
  tibble::tibble(m = m, est = qbar, se = sqrt(tvar))
}

#' Compute standard age/sex populations from harmonized data (imp=1)
#' @return List with $brfss_std and $khandle_std data frames
compute_std_populations <- function(dat_harm) {
  assert_has_vars(dat_harm, c("brfss_h", "imp_h", "age_h", "male_h", "brfss_sampwt_h"), "dat_harm")

  dat_harm <- dat_harm %>%
    dplyr::mutate(agecat_h = make_agecat(age_h))

  brfss_std <- dat_harm %>%
    dplyr::filter(brfss_h == 1, imp_h == 1) %>%
    dplyr::filter(!is.na(agecat_h), !is.na(male_h)) %>%
    dplyr::group_by(agecat_h, male_h) %>%
    dplyr::summarise(n = sum(brfss_sampwt_h, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(prop_agesex = n / sum(n))

  khandle_std <- dat_harm %>%
    dplyr::filter(brfss_h == 0, imp_h == 1) %>%
    dplyr::filter(!is.na(agecat_h), !is.na(male_h)) %>%
    dplyr::group_by(agecat_h, male_h) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(prop_agesex = n / sum(n))

  if (nrow(brfss_std) == 0) stop("BRFSS standard population is empty (check imp_h==1 BRFSS rows).")
  if (nrow(khandle_std) == 0) stop("KHANDLE standard population is empty (check imp_h==1 KHANDLE rows).")

  list(brfss_std = brfss_std, khandle_std = khandle_std)
}

#' Standardize cell-level estimates to a standard population
standardize_from_cells <- function(cell_df, std_df) {
  cell_df %>%
    dplyr::left_join(std_df, by = c("agecat_h", "male_h")) %>%
    dplyr::group_by(race_summary_f) %>%
    dplyr::summarise(
      est_std = sum(est * prop_agesex, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Standardize cell-level estimates with SE propagation
standardize_from_cells_with_se <- function(cells_df, std_df) {
  cells_df %>%
    dplyr::left_join(std_df, by = c("agecat_h", "male_h")) %>%
    dplyr::group_by(race_summary_f) %>%
    dplyr::summarise(
      est_std = sum(est * prop_agesex, na.rm = TRUE),
      var_std = sum((prop_agesex^2) * (se^2), na.rm = TRUE),
      se_std  = sqrt(var_std),
      .groups = "drop"
    ) %>%
    dplyr::select(race_summary_f, est_std, se_std)
}
