#------------------------------------------------------------------
# R/balance_helpers.R
# Purpose: Covariate balance, weight diagnostics, and plotting helpers
# Extracted from 03_weight_development.R
#------------------------------------------------------------------

#' Plain-English labels for harmonized covariates shown in Figure 1.
#' Order matches Hayes-Larson 2022 Table 1.
COVAR_LABELS <- c(
  "age_h"             = "Age (years)",
  "male_h"            = "Male sex",
  "education3_h=1"    = "Education: less than high school",
  "education3_h=2"    = "Education: high school diploma/GED",
  "education3_h=3"    = "Education: trade school/some college",
  "education3_h=4"    = "Education: college graduate or higher",
  "income_gtmed_pp_h" = "Per-capita income above BRFSS median",
  "marital_h"         = "Married/living with partner",
  "military_h"        = "Military service",
  "goodhealth_h"      = "Self-rated health good or better",
  "adl_walking_h"     = "Difficulty walking/climbing stairs",
  "blind_h"           = "Blind or vision impairment",
  "adl_dressing_h"    = "Difficulty dressing",
  "smoke_status_h"    = "Current smoker",
  "exercise_h"        = "Exercised in last month"
)

#' Map raw covariate names (incl. factor level suffixes like "education3_h=1")
#' to plain-English labels. Unknown names pass through unchanged.
pretty_covar_labels <- function(var) {
  out <- COVAR_LABELS[var]
  ifelse(is.na(out), var, out)
}

#' Return a factor with display order matching Hayes-Larson 2022 Table 1.
#' Unknown labels are appended at the end. Reverses for ggplot so that the
#' first variable above appears at the top of the y-axis.
pretty_covar_factor <- function(var) {
  labels <- pretty_covar_labels(var)
  known  <- unname(COVAR_LABELS)
  unknown <- setdiff(unique(labels), known)
  factor(labels, levels = rev(c(known, unknown)))
}

#' Coerce to 0/1 numeric (handles factors/characters)
as01_num <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  x_num <- suppressWarnings(as.numeric(x))
  dplyr::case_when(
    is.na(x_num) ~ NA_real_,
    x_num %in% c(0, 1) ~ x_num,
    TRUE ~ NA_real_
  )
}

#' Assert that values are within an allowed set
assert_allowed_values <- function(x, allowed, varname) {
  vals <- sort(unique(as.character(x)))
  vals <- vals[!is.na(vals)]
  bad <- setdiff(vals, allowed)
  if (length(bad) > 0) {
    stop(
      varname, " has unexpected value(s): ", paste(bad, collapse = ", "),
      "\nAllowed: ", paste(allowed, collapse = ", ")
    )
  }
  invisible(TRUE)
}

#' Effective sample size from weights
effective_n <- function(w) {
  w <- w[!is.na(w) & w > 0]
  if (length(w) == 0) return(NA_real_)
  (sum(w)^2) / sum(w^2)
}

#' GLM diagnostics summary
glm_diag <- function(fit, p, label) {
  tibble::tibble(
    label = label,
    converged = isTRUE(fit$converged),
    any_na_coef = any(is.na(stats::coef(fit))),
    min_p = suppressWarnings(min(p, na.rm = TRUE)),
    p01   = suppressWarnings(stats::quantile(p, 0.01, na.rm = TRUE)),
    p50   = suppressWarnings(stats::quantile(p, 0.50, na.rm = TRUE)),
    p99   = suppressWarnings(stats::quantile(p, 0.99, na.rm = TRUE)),
    max_p = suppressWarnings(max(p, na.rm = TRUE))
  )
}

#' Weight distribution diagnostics
weight_diag <- function(df, w_var, group_var = "race_summary_f") {
  stopifnot(w_var %in% names(df), "brfss_h" %in% names(df), group_var %in% names(df))

  kh <- df %>% dplyr::filter(brfss_h == 0)
  w  <- kh[[w_var]]
  w  <- w[!is.na(w)]

  if (length(w) == 0) {
    return(list(
      overall = tibble::tibble(n = 0),
      by_race = tibble::tibble()
    ))
  }

  qs <- stats::quantile(w, probs = c(0, .5, .9, .95, .99, 1), na.rm = TRUE, names = FALSE)
  overall <- tibble::tibble(
    n = length(w),
    min = qs[1], p50 = qs[2], p90 = qs[3], p95 = qs[4], p99 = qs[5], max = qs[6],
    mean = mean(w), sd = stats::sd(w), cv = stats::sd(w) / mean(w),
    ess = effective_n(w)
  )

  by_race <- df %>%
    dplyr::filter(brfss_h == 0) %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(
      n = sum(!is.na(.data[[w_var]])),
      min = min(.data[[w_var]], na.rm = TRUE),
      p50 = stats::quantile(.data[[w_var]], .50, na.rm = TRUE),
      p95 = stats::quantile(.data[[w_var]], .95, na.rm = TRUE),
      p99 = stats::quantile(.data[[w_var]], .99, na.rm = TRUE),
      max = max(.data[[w_var]], na.rm = TRUE),
      mean = mean(.data[[w_var]], na.rm = TRUE),
      sd = stats::sd(.data[[w_var]], na.rm = TRUE),
      cv = stats::sd(.data[[w_var]], na.rm = TRUE) / mean(.data[[w_var]], na.rm = TRUE),
      ess = effective_n(.data[[w_var]]),
      .groups = "drop"
    ) %>%
    dplyr::rename(race = 1)

  list(overall = overall, by_race = by_race)
}

#' Trim weights at quantile boundaries
trim_weights <- function(w, lower = 0.01, upper = 0.99) {
  stopifnot(is.numeric(w))
  ok <- !is.na(w)
  if (sum(ok) == 0) return(w)
  q <- stats::quantile(w[ok], probs = c(lower, upper), na.rm = TRUE, names = FALSE)
  w2 <- w
  w2[ok] <- pmin(pmax(w[ok], q[1]), q[2])
  w2
}

#------------------------------------------------------------------
# Covariate balance (weighted SMD)
#------------------------------------------------------------------

w_mean <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(w[ok] * x[ok]) / sum(w[ok])
}

w_var <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (sum(ok) < 2) return(NA_real_)
  mu <- w_mean(x[ok], w[ok])
  sum(w[ok] * (x[ok] - mu)^2) / sum(w[ok])
}

w_sd <- function(x, w) sqrt(w_var(x, w))

smd_numeric <- function(x, treat, w_total, treat_is_brfss = 1) {
  ok <- !is.na(x) & !is.na(treat) & !is.na(w_total) & w_total > 0
  if (!any(ok)) return(NA_real_)

  x0 <- x[ok & treat == 0]
  w0 <- w_total[ok & treat == 0]
  x1 <- x[ok & treat == treat_is_brfss]
  w1 <- w_total[ok & treat == treat_is_brfss]

  if (length(x0) < 2 || length(x1) < 2) return(NA_real_)

  mu0 <- w_mean(x0, w0)
  mu1 <- w_mean(x1, w1)
  sd1 <- w_sd(x1, w1)

  if (is.na(sd1) || sd1 == 0) return(NA_real_)
  (mu0 - mu1) / sd1
}

expand_var_to_columns <- function(d, v) {
  x <- d[[v]]

  if (is.factor(x) || is.character(x)) {
    x <- as.factor(x)
    levs <- levels(x)
    out <- vector("list", length(levs))
    names(out) <- paste0(v, "=", levs)
    for (k in seq_along(levs)) out[[k]] <- as.numeric(x == levs[k])
    return(out)
  }

  list(setNames(list(to_num(x)), v))[[1]]
}

safe_covariates_for_stratum <- function(d, vars, treat = "brfss_h", min_n = 10) {
  keep <- c()

  if (length(unique(d[[treat]][!is.na(d[[treat]])])) < 2) return(keep)

  for (v in vars) {
    if (!v %in% names(d)) next
    x <- d[[v]]

    n0 <- sum(!is.na(x) & d[[treat]] == 0)
    n1 <- sum(!is.na(x) & d[[treat]] == 1)
    if (n0 < min_n || n1 < min_n) next

    if (is.factor(x) || is.character(x)) {
      if (length(unique(as.character(x[!is.na(x)]))) < 2) next
    } else {
      if (length(unique(to_num(x[!is.na(x)]))) < 2) next
    }

    keep <- c(keep, v)
  }

  keep
}

compute_covbal_by_race <- function(data, vars, weightvar,
                                   race_var = "race_summary_f",
                                   treat = "brfss_h",
                                   min_n = 10) {

  races <- levels(data[[race_var]])
  out <- vector("list", length(races))
  names(out) <- races

  for (r in races) {
    d <- data[data[[race_var]] == r, , drop = FALSE]
    if (length(unique(d[[treat]][!is.na(d[[treat]])])) < 2) next

    v_keep <- safe_covariates_for_stratum(d, vars, treat = treat, min_n = min_n)
    if (length(v_keep) == 0) next

    w_total <- d$brfss_sampwt_h * d[[weightvar]]

    res_rows <- list()

    for (v in v_keep) {
      cols <- expand_var_to_columns(d, v)

      if (is.list(cols) && !is.null(names(cols))) {
        for (nm in names(cols)) {
          smd <- smd_numeric(cols[[nm]], d[[treat]], w_total, treat_is_brfss = 1)
          res_rows[[length(res_rows) + 1]] <- tibble::tibble(var = nm, std.eff.sz = smd)
        }
      } else {
        smd <- smd_numeric(cols, d[[treat]], w_total, treat_is_brfss = 1)
        res_rows[[length(res_rows) + 1]] <- tibble::tibble(var = v, std.eff.sz = smd)
      }
    }

    res <- dplyr::bind_rows(res_rows)
    if (nrow(res) == 0) next
    res$race <- r
    out[[r]] <- res
  }

  out_df <- dplyr::bind_rows(out)
  if (nrow(out_df) == 0) return(out_df)

  out_df %>% dplyr::select(race, var, std.eff.sz)
}

compute_covbal_overall <- function(data, vars, weightvar, treat = "brfss_h") {
  w_total <- data$brfss_sampwt_h * data[[weightvar]]
  res_rows <- list()

  for (v in vars) {
    cols <- expand_var_to_columns(data, v)

    if (is.list(cols) && !is.null(names(cols))) {
      for (nm in names(cols)) {
        smd <- smd_numeric(cols[[nm]], data[[treat]], w_total, treat_is_brfss = 1)
        res_rows[[length(res_rows) + 1]] <- tibble::tibble(var = nm, std.eff.sz = smd)
      }
    } else {
      smd <- smd_numeric(cols, data[[treat]], w_total, treat_is_brfss = 1)
      res_rows[[length(res_rows) + 1]] <- tibble::tibble(var = v, std.eff.sz = smd)
    }
  }

  dplyr::bind_rows(res_rows)
}

#------------------------------------------------------------------
# Balance plotting
#------------------------------------------------------------------

plot_covbal_paper <- function(before_df, after_df,
                              title = NULL,
                              y_lim = c(-1, 1.25)) {

  df <- dplyr::bind_rows(
    before_df %>% dplyr::mutate(stage = "Unweighted"),
    after_df  %>% dplyr::mutate(stage = "Weighted")
  ) %>%
    dplyr::filter(is.finite(std.eff.sz)) %>%
    dplyr::mutate(abs_smd = abs(std.eff.sz))

  if (nrow(df) == 0) return(NULL)

  var_order <- df %>%
    dplyr::filter(stage == "Unweighted") %>%
    dplyr::group_by(var) %>%
    dplyr::summarise(m = mean(abs_smd, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(m) %>%
    dplyr::pull(var)

  df <- df %>%
    dplyr::mutate(var = factor(var, levels = var_order))

  ggplot2::ggplot(df, ggplot2::aes(x = std.eff.sz, y = var, colour = race)) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_point(size = 2) +
    ggplot2::coord_cartesian(xlim = y_lim) +
    ggplot2::facet_wrap(~ stage, nrow = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = title,
      x = "Standardized mean difference (KHANDLE - BRFSS) / SD[BRFSS]",
      colour = "Race/ethnicity"
    )
}

plot_covbal_overall_paper <- function(before_df, after_df,
                                      title = NULL,
                                      y_lim = c(-1, 1)) {

  df <- dplyr::bind_rows(
    before_df %>% dplyr::mutate(stage = "Unweighted"),
    after_df  %>% dplyr::mutate(stage = "Weighted")
  ) %>%
    dplyr::filter(is.finite(std.eff.sz)) %>%
    dplyr::mutate(abs_smd = abs(std.eff.sz))

  if (nrow(df) == 0) return(NULL)

  var_order <- df %>%
    dplyr::filter(stage == "Unweighted") %>%
    dplyr::group_by(var) %>%
    dplyr::summarise(m = mean(abs_smd, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(m) %>%
    dplyr::pull(var)

  df <- df %>%
    dplyr::mutate(var = factor(var, levels = var_order))

  ggplot2::ggplot(df, ggplot2::aes(x = std.eff.sz, y = var)) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_point(size = 2) +
    ggplot2::coord_cartesian(xlim = y_lim) +
    ggplot2::facet_wrap(~ stage, nrow = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = title,
      x = "Standardized mean difference (KHANDLE - BRFSS) / SD[BRFSS]"
    )
}

summarize_balance_improvement <- function(before_df, after_df) {
  out <- before_df %>%
    dplyr::select(var, smd_before = std.eff.sz) %>%
    dplyr::left_join(after_df %>% dplyr::select(var, smd_after = std.eff.sz), by = "var") %>%
    dplyr::mutate(
      abs_before = abs(smd_before),
      abs_after  = abs(smd_after),
      improved   = abs_after < abs_before
    )

  metrics <- tibble::tibble(
    mean_abs_smd_before = mean(out$abs_before, na.rm = TRUE),
    mean_abs_smd_after  = mean(out$abs_after,  na.rm = TRUE),
    p90_abs_smd_before  = stats::quantile(out$abs_before, 0.90, na.rm = TRUE, names = FALSE),
    p90_abs_smd_after   = stats::quantile(out$abs_after,  0.90, na.rm = TRUE, names = FALSE),
    prop_improved       = mean(out$improved, na.rm = TRUE),
    n_covariates        = sum(!is.na(out$abs_before) & !is.na(out$abs_after))
  )

  list(detail = out, metrics = metrics)
}

#------------------------------------------------------------------
# Figure saving
#------------------------------------------------------------------

save_fig <- function(p, name_stub, width = 10, height = 6, fig_dir = NULL) {
  stopifnot(inherits(p, "ggplot"))
  if (is.null(fig_dir)) fig_dir <- file.path(here::here("03_results"), "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file.path(fig_dir, paste0(name_stub, ".png")),
                  plot = p, width = width, height = height, units = "in", dpi = 300)
  ggplot2::ggsave(file.path(fig_dir, paste0(name_stub, ".pdf")),
                  plot = p, width = width, height = height, units = "in")
}
