#------------------------------------------------------------------
# Script: 06_results_summary.R
# Purpose: Summarize IOSW results + diagnostics + produce figures/tables
# Repo: github.com/lamhine/KHANDLE-weighting-multiracial
# Last updated: 2026-04-08
#
# Outputs (clear labeling):
#   TABLE A + FIGURE 2 (main): Prevalence (Unweighted KHANDLE vs Generalized-to-BRFSS), NOT standardized
#   FIGURE 1 (main): Overall covariate balance overlay (Unweighted vs Weighted), NOT race-stratified
#   TABLE C + FIGURE 3 (main): PR/PD vs White (standardized), clustered series (Generalized vs Unweighted) WITH error bars
#   TABLE B + FIGURE S1 (supplement): Age/sex standardized prevalence to BRFSS across tags
#
# Notes:
#   - Helpers from R/utils.R, R/race_constants.R, R/analysis_helpers.R, R/balance_helpers.R via config
#   - ALLOWED_RACE5, REF_RACE, make_agecat, wmean_se, umean_se, rubin_combine,
#     standardize_from_cells_with_se all provided by shared files
#------------------------------------------------------------------

source("config/config.R")

suppressPackageStartupMessages({
  library(tidyr)
})

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package 'openxlsx' is required for .xlsx outputs. Install via install.packages('openxlsx').")
}

message("06_results_summary.R using analysis_tag = ", analysis_tag)

#------------------------------------------------------------------
# Helpers (script-specific)
#------------------------------------------------------------------

ci95 <- function(est, se, z = 1.96) {
  tibble(ci_lo = est - z * se, ci_hi = est + z * se)
}

write_xlsx <- function(x, path) {
  openxlsx::write.xlsx(x, file = path, overwrite = TRUE)
}

first_existing <- function(paths, required = TRUE, label = NULL) {
  ok <- paths[file.exists(paths)]
  if (length(ok) > 0) return(ok[1])
  if (required) {
    if (is.null(label)) label <- "file"
    stop("Could not find required ", label, " among:\n", paste(paths, collapse = "\n"))
  }
  NA_character_
}

find_by_pattern <- function(dir, pattern, required = TRUE, label = NULL) {
  hits <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(hits) > 0) return(hits[1])
  if (required) {
    if (is.null(label)) label <- "file"
    stop("Could not find required ", label, " in ", dir, " with pattern: ", pattern)
  }
  NA_character_
}

discover_tags <- function(processed_parent_dir,
                          preferred = c("t10_replication", "t16_update"),
                          fallback_tag = NULL) {
  dirs <- list.dirs(processed_parent_dir, full.names = FALSE, recursive = FALSE)
  tags <- preferred[preferred %in% dirs]
  if (length(tags) == 0 && !is.null(fallback_tag) && fallback_tag %in% dirs) tags <- fallback_tag
  if (length(tags) == 0) stop("No analysis_tag folders found under processed_data_parent_dir.")
  tags
}

read_iosw_outputs <- function(tag, processed_parent_dir) {
  p <- file.path(processed_parent_dir, tag)

  files <- c(
    overall = paste0("04_iosw_results_overall_", tag, ".rds"),
    byrace  = paste0("04_iosw_results_by_race_", tag, ".rds"),
    std     = paste0("04_iosw_results_standardized_", tag, ".rds"),
    prpd    = paste0("04_iosw_results_pr_pd_", tag, ".rds")
  )

  paths <- setNames(file.path(p, files), names(files))

  if (!all(file.exists(paths))) {
    missing <- names(paths)[!file.exists(paths)]
    stop(
      "Missing IOSW output(s) in ", p, ":\n",
      paste0("  - ", missing, ": ", paths[missing], collapse = "\n")
    )
  }

  list(
    overall = readRDS(paths["overall"]),
    byrace  = readRDS(paths["byrace"]),
    std     = readRDS(paths["std"]),
    prpd    = readRDS(paths["prpd"]),
    dir     = p
  )
}

#------------------------------------------------------------------
# Output folders (clear separation)
#------------------------------------------------------------------

tags <- discover_tags(
  processed_parent_dir = processed_data_parent_dir,
  preferred = c("t10_replication", "t16_update"),
  fallback_tag = analysis_tag
)

message("06_results_summary.R will summarize tags: ", paste(tags, collapse = ", "))

results_parent_dir <- results_dir

fig_dir  <- file.path(results_parent_dir, "figures")
tab_dir  <- file.path(results_parent_dir, "tables")
supp_dir <- file.path(results_parent_dir, "supplement")

dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

tagset_label <- paste(tags, collapse = "__")

#------------------------------------------------------------------
# FIGURE 1 (main): Overall covariate balance overlay (unweighted vs weighted), rotated
#
# NOTE: Uses true overall SMD (compute_covbal_overall) rather than
# averaging within-race SMDs, which can be misleading when KHANDLE
# composition differs dramatically by race (e.g., education).
# Falls back to by-race averaging if overall files don't exist.
#------------------------------------------------------------------

read_covbal_overall_list <- function(tag_dir, prefix_candidates = c("03", "04")) {
  # Prefer true overall SMD files; fall back to by-race files
  candidates_unw_overall <- file.path(tag_dir, paste0(prefix_candidates, "_covbal_unweighted_overall_by_imp.rds"))
  candidates_wt_overall  <- file.path(tag_dir, paste0(prefix_candidates, "_covbal_weighted_overall_by_imp.rds"))

  unw_path <- first_existing(candidates_unw_overall, required = FALSE)
  wt_path  <- first_existing(candidates_wt_overall,  required = FALSE)

  if (!is.na(unw_path) && !is.na(wt_path)) {
    return(list(unw = readRDS(unw_path), wt = readRDS(wt_path), type = "overall"))
  }

  # Fallback: by-race files (old approach)
  candidates_unw <- file.path(tag_dir, paste0(prefix_candidates, "_covbal_unweighted_by_imp.rds"))
  candidates_wt  <- file.path(tag_dir, paste0(prefix_candidates, "_covbal_weighted_by_imp.rds"))

  unw_path <- first_existing(candidates_unw, required = FALSE)
  wt_path  <- first_existing(candidates_wt,  required = FALSE)

  if (is.na(unw_path) || is.na(wt_path)) return(NULL)

  message("NOTE: Using by-race covbal files (overall files not found). ",
          "Re-run 03_weight_development.R to generate true overall SMDs.")
  list(unw = readRDS(unw_path), wt = readRDS(wt_path), type = "by_race")
}

bind_covbal_list <- function(covbal_list, wtd_label) {
  out <- lapply(names(covbal_list), function(k) {
    df <- covbal_list[[k]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>% mutate(imp = as.integer(k), wtd = wtd_label)
  })
  bind_rows(out)
}

collapse_covbal_to_overall <- function(df, type = "overall") {
  if (!("var" %in% names(df))) stop("Covbal df missing 'var'.")
  if (!("std.eff.sz" %in% names(df))) stop("Covbal df missing 'std.eff.sz'.")

  if (type == "overall" || !("race" %in% names(df))) {
    # True overall SMD (already computed at the population level)
    df %>%
      transmute(var, imp, wtd, std_eff_sz = std.eff.sz)
  } else {
    # Fallback: average within-race SMDs (less accurate)
    df %>%
      group_by(var, imp, wtd) %>%
      summarise(std_eff_sz = mean(std.eff.sz, na.rm = TRUE), .groups = "drop")
  }
}

for (tag in tags) {

  tag_dir <- file.path(processed_data_parent_dir, tag)
  covbal <- read_covbal_overall_list(tag_dir)

  if (is.null(covbal)) {
    message("Covariate balance stacks not found for tag ", tag, "; skipping Figure 1.")
    next
  }

  cov_unw_df <- bind_covbal_list(covbal$unw, "Unweighted")
  cov_wt_df  <- bind_covbal_list(covbal$wt,  "Weighted")

  plotdat <- bind_rows(cov_unw_df, cov_wt_df) %>%
    filter(!is.na(std.eff.sz)) %>%
    collapse_covbal_to_overall(type = covbal$type) %>%
    group_by(var, wtd) %>%
    summarise(mean_smd = mean(std_eff_sz, na.rm = TRUE), .groups = "drop") %>%
    mutate(wtd = factor(wtd, levels = c("Unweighted", "Weighted")))

  fig1 <- ggplot(plotdat, aes(x = mean_smd, y = var, colour = wtd)) +
    geom_vline(xintercept = 0) +
    geom_point(size = 2.4, position = position_dodge(width = 0.35)) +
    theme_bw() +
    labs(
      title = paste0("Figure 1. Overall covariate balance before and after generalizability weighting (", tag, ")"),
      x = "Standardized mean difference (KHANDLE - BRFSS) / SD[BRFSS]",
      y = NULL,
      colour = NULL
    ) +
    theme(legend.position = "bottom") +
    coord_cartesian(xlim = c(-1, 1.25))

  fig1_path <- file.path(fig_dir, paste0("Figure1_covariate_balance_overall_", tag, ".png"))
  ggsave(fig1_path, plot = fig1, width = 9, height = 6.5, dpi = 300)

  message("Saved Figure 1 to:\n  ", fig1_path)
}

#------------------------------------------------------------------
# TABLE A + FIGURE 2 (main): Prevalence (unweighted KHANDLE vs generalized-to-BRFSS), NOT standardized
#------------------------------------------------------------------

make_prev_unw_wtd_tbl <- function(res, tag) {

  allowed_order <- c("Overall", "Asian", "Black", "Hispanic", "Multiracial", "White")

  overall_tbl <- res$overall %>%
    filter(stat %in% c("Overall_unweighted", "Overall_weighted")) %>%
    transmute(
      analysis_tag = tag,
      group = "Overall",
      series = case_when(
        stat == "Overall_unweighted" ~ "Unweighted KHANDLE",
        stat == "Overall_weighted"   ~ "Generalized to BRFSS (weighted)",
        TRUE ~ NA_character_
      ),
      est = est,
      se  = se
    )

  byrace_tbl <- res$byrace %>%
    filter(stat %in% c("By_race_unweighted", "By_race_weighted")) %>%
    transmute(
      analysis_tag = tag,
      group = race,
      series = case_when(
        stat == "By_race_unweighted" ~ "Unweighted KHANDLE",
        stat == "By_race_weighted"   ~ "Generalized to BRFSS (weighted)",
        TRUE ~ NA_character_
      ),
      est = est,
      se  = se
    )

  bind_rows(overall_tbl, byrace_tbl) %>%
    mutate(
      group = factor(group, levels = allowed_order),
      series = factor(series, levels = c("Unweighted KHANDLE", "Generalized to BRFSS (weighted)"))
    ) %>%
    bind_cols(ci95(.$est, .$se)) %>%
    mutate(
      est_pct   = 100 * est,
      ci_lo_pct = 100 * ci_lo,
      ci_hi_pct = 100 * ci_hi
    ) %>%
    arrange(group, series)
}

plot_prev_unw_wtd <- function(prev_tbl, title = NULL) {
  ggplot(prev_tbl, aes(x = group, y = est_pct, colour = series)) +
    geom_point(position = position_dodge(width = 0.5), size = 2.5) +
    geom_errorbar(
      aes(ymin = ci_lo_pct, ymax = ci_hi_pct),
      position = position_dodge(width = 0.5),
      width = 0
    ) +
    theme_bw() +
    labs(
      title = title,
      x = NULL,
      y = "Prevalence (%)",
      colour = NULL
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

for (tag in tags) {

  res <- read_iosw_outputs(tag, processed_data_parent_dir)
  prev_tbl_unstd <- make_prev_unw_wtd_tbl(res, tag)

  tabA_path <- file.path(tab_dir, paste0("TableA_prevalence_unweighted_vs_generalized_", tag, ".xlsx"))
  write_xlsx(
    prev_tbl_unstd %>%
      mutate(
        est_pct = round(est_pct, 1),
        ci_lo_pct = round(ci_lo_pct, 1),
        ci_hi_pct = round(ci_hi_pct, 1)
      ),
    tabA_path
  )

  fig2 <- plot_prev_unw_wtd(
    prev_tbl_unstd,
    title = paste0("Figure 2. Cognitive impairment prevalence: unweighted KHANDLE vs generalized to CA-BRFSS (", tag, ")")
  )

  fig2_path <- file.path(fig_dir, paste0("Figure2_prevalence_unweighted_vs_generalized_", tag, ".png"))
  ggsave(fig2_path, plot = fig2, width = 9, height = 5.5, dpi = 300)

  message("Saved Table A to:\n  ", tabA_path)
  message("Saved Figure 2 to:\n  ", fig2_path)
}


#------------------------------------------------------------------
# TABLE C + FIGURE 3 (main): PR/PD vs White (standardized), clustered series WITH error bars
#------------------------------------------------------------------

outcome_var <- "cogimp_prob_fin_dx"

# Compute per-imputation PR/PD AND their SEs
compute_prpd_by_imp <- function(dat_wts, dat_harm, outcome_var) {

  dat_wts <- dat_wts %>%
    mutate(
      race_summary_f = factor(race_summary_h, levels = 1:5, labels = ALLOWED_RACE5),
      agecat_h = make_agecat(age_h)
    )

  dat_harm <- dat_harm %>%
    mutate(agecat_h = make_agecat(age_h))

  std_pops <- compute_std_populations(dat_harm)
  brfss_std <- std_pops$brfss_std
  khandle_std <- std_pops$khandle_std

  imps <- sort(unique(dat_wts$imp_h))
  imps <- imps[!is.na(imps)]

  out <- lapply(imps, function(j) {

    d <- dat_wts %>% filter(imp_h == j, brfss_h == 0)
    if (nrow(d) == 0) return(NULL)

    # UNWEIGHTED cells -> KHANDLE standard
    cells_unw <- d %>%
      filter(!is.na(agecat_h), !is.na(male_h)) %>%
      group_by(race_summary_f, agecat_h, male_h) %>%
      summarise(
        est = umean_se(.data[[outcome_var]])["est"],
        se  = umean_se(.data[[outcome_var]])["se"],
        .groups = "drop"
      )

    std_unw <- standardize_from_cells_with_se(cells_unw, khandle_std)

    # WEIGHTED cells -> BRFSS standard
    cells_wt <- d %>%
      filter(!is.na(agecat_h), !is.na(male_h)) %>%
      group_by(race_summary_f, agecat_h, male_h) %>%
      summarise(
        est = wmean_se(.data[[outcome_var]], sw_final)["est"],
        se  = wmean_se(.data[[outcome_var]], sw_final)["se"],
        .groups = "drop"
      )

    std_wt <- standardize_from_cells_with_se(cells_wt, brfss_std)

    # Reference values
    ref_unw <- std_unw %>% filter(race_summary_f == REF_RACE)
    ref_wt  <- std_wt  %>% filter(race_summary_f == REF_RACE)

    if (nrow(ref_unw) == 0 || nrow(ref_wt) == 0) return(NULL)

    # Delta method for PR/PD SEs
    make_prpd <- function(std_df, ref_row, series_label) {

      p_ref  <- ref_row$est_std[1]
      se_ref <- ref_row$se_std[1]

      std_df %>%
        mutate(
          imp_h  = j,
          series = series_label,

          PR = est_std / p_ref,
          PD = est_std - p_ref,

          se_PD = sqrt(se_std^2 + se_ref^2),
          se_PR = sqrt((se_std^2) / (p_ref^2) + ((est_std^2) * (se_ref^2)) / (p_ref^4))
        ) %>%
        select(imp_h, race_summary_f, series, PR, PD, se_PR, se_PD)
    }

    prpd_unw <- make_prpd(std_unw, ref_unw, "Unweighted KHANDLE")
    prpd_wt  <- make_prpd(std_wt,  ref_wt,  "Generalized to BRFSS (weighted)")

    bind_rows(prpd_unw, prpd_wt)
  })

  bind_rows(out)
}

summarise_prpd <- function(prpd_long, tag) {

  prpd_long %>%
    mutate(
      analysis_tag = tag,
      race = as.character(race_summary_f),
      comparison = paste0(race, " vs ", REF_RACE),
      series = factor(series, levels = c("Unweighted KHANDLE", "Generalized to BRFSS (weighted)"))
    ) %>%
    filter(race != REF_RACE) %>%
    pivot_longer(
      cols = c(PR, PD, se_PR, se_PD),
      names_to = "name",
      values_to = "value"
    ) %>%
    mutate(
      measure = case_when(
        name %in% c("PR", "se_PR") ~ "PR",
        name %in% c("PD", "se_PD") ~ "PD",
        TRUE ~ NA_character_
      ),
      field = case_when(
        name %in% c("PR", "PD") ~ "est",
        name %in% c("se_PR", "se_PD") ~ "se",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(measure), !is.na(field)) %>%
    select(-name) %>%
    pivot_wider(names_from = field, values_from = value) %>%
    group_by(analysis_tag, series, race, comparison, measure) %>%
    summarise(
      rubin_combine(est = est, se = se),
      .groups = "drop"
    ) %>%
    bind_cols(ci95(.$est, .$se)) %>%
    mutate(
      est_plot = if_else(measure == "PD", 100 * est, est),
      lo_plot  = if_else(measure == "PD", 100 * ci_lo, ci_lo),
      hi_plot  = if_else(measure == "PD", 100 * ci_hi, ci_hi)
    )
}

plot_prpd <- function(df, tag) {

  desired <- c("Asian vs White", "Black vs White", "Hispanic vs White", "Multiracial vs White")

  df <- df %>%
    mutate(
      comparison = factor(comparison, levels = desired),
      series = factor(series, levels = c("Unweighted KHANDLE", "Generalized to BRFSS (weighted)")),
      measure = factor(measure, levels = c("PR", "PD"))
    )

  ref_lines <- tibble(
    measure = factor(c("PR", "PD"), levels = c("PR", "PD")),
    y0 = c(1, 0)
  )

  ggplot(df, aes(x = comparison, y = est_plot, colour = series)) +
    geom_hline(data = ref_lines, aes(yintercept = y0), inherit.aes = FALSE) +
    geom_pointrange(
      aes(ymin = lo_plot, ymax = hi_plot),
      position = position_dodge(width = 0.5),
      linewidth = 0.6
    ) +
    facet_wrap(~ measure, scales = "free_y", nrow = 1) +
    theme_bw() +
    labs(
      title = paste0("Figure 3. Age- and sex-adjusted PR/PD vs White (", tag, ")"),
      x = NULL,
      y = NULL,
      colour = NULL
    ) +
    theme(legend.position = "bottom")
}

for (tag in tags) {

  res <- read_iosw_outputs(tag, processed_data_parent_dir)

  wts_path <- find_by_pattern(
    dir = res$dir,
    pattern = "khandle_brfss_harmonized_weights.*\\.rds$",
    required = TRUE,
    label = "harmonized weights stack"
  )

  harm_path <- find_by_pattern(
    dir = res$dir,
    pattern = "02_khandle_brfss_harmonized\\.rds$|khandle_brfss_harmonized\\.rds$",
    required = TRUE,
    label = "harmonized data"
  )

  dat_wts  <- readRDS(wts_path)
  dat_harm <- readRDS(harm_path)

  if (!(outcome_var %in% names(dat_wts))) stop("Missing outcome in weights stack for tag ", tag, ": ", outcome_var)

  req_harm <- c("brfss_h", "imp_h", "age_h", "male_h", "brfss_sampwt_h")
  if (length(setdiff(req_harm, names(dat_harm))) > 0) {
    stop("Harmonized data missing required vars for PR/PD (", tag, "): ",
         paste(setdiff(req_harm, names(dat_harm)), collapse = ", "))
  }

  req_wts <- c("brfss_h", "imp_h", "age_h", "male_h", "race_summary_h", "sw_final")
  if (length(setdiff(req_wts, names(dat_wts))) > 0) {
    stop("Weights stack missing required vars for PR/PD (", tag, "): ",
         paste(setdiff(req_wts, names(dat_wts)), collapse = ", "))
  }

  prpd_by_imp <- compute_prpd_by_imp(dat_wts, dat_harm, outcome_var)

  tabC <- summarise_prpd(prpd_by_imp, tag)
  tabC_path <- file.path(tab_dir, paste0("TableC_pr_pd_vs_white_", tag, ".xlsx"))
  write_xlsx(tabC, tabC_path)

  fig3 <- plot_prpd(tabC, tag)
  fig3_path <- file.path(fig_dir, paste0("Figure3_pr_pd_vs_white_", tag, ".png"))
  ggsave(fig3_path, plot = fig3, width = 11, height = 5.5, dpi = 300)

  message("Saved Table C to:\n  ", tabC_path)
  message("Saved Figure 3 to:\n  ", fig3_path)
}

#------------------------------------------------------------------
# TABLE B + FIGURE S1 (supplement): Age/sex standardized prevalence to BRFSS across tags
#------------------------------------------------------------------

extract_abs_std_brfss <- function(res, tag) {
  out <- res$std %>%
    filter(stat == "Std_to_BRFSS_weighted") %>%
    transmute(analysis_tag = tag, group = race, est = est)

  # Borrow SE from by-race weighted (approx; same as prior approach)
  se_tbl <- res$byrace %>%
    filter(stat == "By_race_weighted") %>%
    select(race, se)

  out %>%
    left_join(se_tbl, by = c("group" = "race")) %>%
    bind_cols(ci95(.$est, .$se)) %>%
    mutate(
      est_pct = 100 * est,
      ci_lo_pct = 100 * ci_lo,
      ci_hi_pct = 100 * ci_hi
    )
}

extract_abs_overall_std_brfss <- function(res, tag) {

  std_by_race <- res$std %>%
    filter(stat == "Std_to_BRFSS_weighted") %>%
    transmute(race, est_r = est)

  wts_path <- find_by_pattern(
    dir = res$dir,
    pattern = "khandle_brfss_harmonized_weights.*\\.rds$",
    required = TRUE,
    label = "harmonized weights stack"
  )

  dat_wts <- readRDS(wts_path)

  dat_wts <- dat_wts %>%
    mutate(race_summary_f = factor(race_summary_h, levels = 1:5, labels = ALLOWED_RACE5))

  kh_race_prop <- dat_wts %>%
    filter(brfss_h == 0, imp_h == 1) %>%
    group_by(race_summary_f) %>%
    summarise(w = sum(sw_final, na.rm = TRUE), .groups = "drop") %>%
    mutate(prop = w / sum(w)) %>%
    transmute(race = as.character(race_summary_f), prop)

  overall_est <- std_by_race %>%
    left_join(kh_race_prop, by = "race") %>%
    summarise(est = sum(est_r * prop, na.rm = TRUE)) %>%
    pull(est)

  overall_se <- res$overall %>%
    filter(stat == "Overall_weighted") %>%
    summarise(se = se[1]) %>%
    pull(se)

  tibble(analysis_tag = tag, group = "Overall", est = overall_est, se = overall_se) %>%
    bind_cols(ci95(.$est, .$se)) %>%
    mutate(
      est_pct = 100 * est,
      ci_lo_pct = 100 * ci_lo,
      ci_hi_pct = 100 * ci_hi
    )
}

prev_std_list <- lapply(tags, function(tag) {
  res <- read_iosw_outputs(tag, processed_data_parent_dir)
  bind_rows(
    extract_abs_overall_std_brfss(res, tag),
    extract_abs_std_brfss(res, tag)
  )
})

prev_std_tbl <- bind_rows(prev_std_list) %>%
  mutate(
    group = factor(group, levels = c("Overall", "Asian", "Black", "Hispanic", "Multiracial", "White")),
    analysis_tag = factor(analysis_tag, levels = tags)
  ) %>%
  arrange(analysis_tag, group)

tabB_path <- file.path(supp_dir, paste0("TableB_prevalence_age_sex_standardized_to_BRFSS_", tagset_label, ".xlsx"))
write_xlsx(prev_std_tbl, tabB_path)

figS1 <- ggplot(prev_std_tbl, aes(x = group, y = est_pct)) +
  geom_pointrange(aes(ymin = ci_lo_pct, ymax = ci_hi_pct), linewidth = 0.5) +
  facet_wrap(~ analysis_tag, ncol = 1) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Prevalence (%) (age/sex standardized to BRFSS)",
    title = "Figure S1. Cognitive impairment prevalence (age/sex standardized to BRFSS)"
  )

figS1_path <- file.path(supp_dir, paste0("FigureS1_prevalence_age_sex_standardized_to_BRFSS_", tagset_label, ".png"))
ggsave(figS1_path, plot = figS1, width = 9, height = 6, dpi = 300)

message("Saved Table B to:\n  ", tabB_path)
message("Saved Figure S1 to:\n  ", figS1_path)

#------------------------------------------------------------------
# Done
#------------------------------------------------------------------

message("06_results_summary.R completed successfully.")
message("Figures saved to: ", fig_dir)
message("Tables saved to:  ", tab_dir)
message("Supplement saved to: ", supp_dir)
