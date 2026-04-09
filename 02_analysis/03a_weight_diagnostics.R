#------------------------------------------------------------------
# Script: 03a_weight_diagnostics.R
# Purpose: Diagnose education balance degradation after IOSW weighting
# Run AFTER: 03_weight_development.R
# Notes:
#   - Investigates why education3_h=4 (college+) SMD worsens after weighting
#   - Examines weight distributions by education level
#   - Compares trimmed vs untrimmed effects
#   - Helps decide whether to modify the weight model
#------------------------------------------------------------------

source("config/config.R")

#------------------------------------------------------------------
# Load data
#------------------------------------------------------------------

dat_wts <- readRDS(file.path(processed_data_dir, "03_khandle_brfss_harmonized_weights.rds"))
dat_harm <- readRDS(file.path(processed_data_dir, "02_khandle_brfss_harmonized.rds"))

dat_wts <- add_race_factor(dat_wts)
dat_harm <- add_race_factor(dat_harm)

# Use imputation 1 for diagnostics
d1 <- dat_wts %>% filter(imp_h == 1)
k1 <- d1 %>% filter(brfss_h == 0)
b1 <- d1 %>% filter(brfss_h == 1)

cat("\n====================================================================\n")
cat("DIAGNOSTIC 1: Education distribution in KHANDLE vs BRFSS (imp=1)\n")
cat("====================================================================\n")

cat("\n--- KHANDLE unweighted counts ---\n")
print(table(k1$education3_h, useNA = "ifany"))
cat("\n--- KHANDLE unweighted proportions ---\n")
print(prop.table(table(k1$education3_h, useNA = "ifany")))

cat("\n--- BRFSS unweighted counts ---\n")
print(table(b1$education3_h, useNA = "ifany"))
cat("\n--- BRFSS weighted proportions (using brfss_sampwt_h) ---\n")
brfss_edu_wtd <- tapply(b1$brfss_sampwt_h, b1$education3_h, sum, na.rm = TRUE)
print(brfss_edu_wtd / sum(brfss_edu_wtd))

cat("\n--- Gap: KHANDLE prop - BRFSS weighted prop ---\n")
kh_prop <- prop.table(table(k1$education3_h))
br_prop <- brfss_edu_wtd / sum(brfss_edu_wtd)
gap <- as.numeric(kh_prop) - as.numeric(br_prop[names(kh_prop)])
names(gap) <- names(kh_prop)
print(round(gap, 3))

cat("\n====================================================================\n")
cat("DIAGNOSTIC 2: Education x Race cross-tab (imp=1)\n")
cat("====================================================================\n")

cat("\n--- KHANDLE: education3_h by race (counts) ---\n")
print(table(k1$race_summary_f, k1$education3_h, useNA = "ifany"))

cat("\n--- KHANDLE: education3_h by race (row proportions) ---\n")
print(round(prop.table(table(k1$race_summary_f, k1$education3_h), margin = 1), 3))

cat("\n--- BRFSS: education3_h by race (counts) ---\n")
print(table(b1$race_summary_f, b1$education3_h, useNA = "ifany"))

cat("\n====================================================================\n")
cat("DIAGNOSTIC 3: Weight distribution by education level (KHANDLE only, imp=1)\n")
cat("====================================================================\n")

cat("\n--- sw_final summary by education3_h ---\n")
wt_by_edu <- k1 %>%
  group_by(education3_h) %>%
  summarise(
    n = n(),
    mean_wt = mean(sw_final, na.rm = TRUE),
    sd_wt = sd(sw_final, na.rm = TRUE),
    min_wt = min(sw_final, na.rm = TRUE),
    p25_wt = quantile(sw_final, 0.25, na.rm = TRUE),
    median_wt = median(sw_final, na.rm = TRUE),
    p75_wt = quantile(sw_final, 0.75, na.rm = TRUE),
    p95_wt = quantile(sw_final, 0.95, na.rm = TRUE),
    max_wt = max(sw_final, na.rm = TRUE),
    .groups = "drop"
  )
print(wt_by_edu)

cat("\n--- sw_final summary by education3_h x race ---\n")
wt_by_edu_race <- k1 %>%
  group_by(race_summary_f, education3_h) %>%
  summarise(
    n = n(),
    mean_wt = round(mean(sw_final, na.rm = TRUE), 3),
    median_wt = round(median(sw_final, na.rm = TRUE), 3),
    max_wt = round(max(sw_final, na.rm = TRUE), 3),
    .groups = "drop"
  )
print(wt_by_edu_race, n = Inf)

cat("\n====================================================================\n")
cat("DIAGNOSTIC 4: Effective reweighted education distribution (KHANDLE)\n")
cat("====================================================================\n")

cat("\n--- After IOSW weighting, what does KHANDLE education look like? ---\n")
k1_edu_wtd <- tapply(k1$sw_final, k1$education3_h, sum, na.rm = TRUE)
cat("KHANDLE IOSW-weighted proportions:\n")
print(round(k1_edu_wtd / sum(k1_edu_wtd), 3))
cat("\nBRFSS survey-weighted proportions (target):\n")
print(round(br_prop, 3))
cat("\nRemaining gap after IOSW:\n")
remaining_gap <- as.numeric(k1_edu_wtd / sum(k1_edu_wtd)) - as.numeric(br_prop[names(kh_prop)])
names(remaining_gap) <- names(kh_prop)
print(round(remaining_gap, 3))

cat("\n====================================================================\n")
cat("DIAGNOSTIC 5: Effect of weight trimming on education balance\n")
cat("====================================================================\n")

# Check if p_hat is available; if not, re-estimate from the model
if ("p_hat" %in% names(k1)) {
  cat("\n--- p_hat distribution by education3_h (KHANDLE, imp=1) ---\n")
  phat_by_edu <- k1 %>%
    group_by(education3_h) %>%
    summarise(
      n = n(),
      mean_phat = round(mean(p_hat, na.rm = TRUE), 4),
      min_phat = round(min(p_hat, na.rm = TRUE), 4),
      max_phat = round(max(p_hat, na.rm = TRUE), 4),
      .groups = "drop"
    )
  print(phat_by_edu)
  cat("\nInterpretation: p_hat is P(KHANDLE | covariates). Higher p_hat means the\n")
  cat("model thinks that covariate profile is more 'KHANDLE-like'. These people\n")
  cat("get LOWER weights (downweighted). If education3_h=4 has high p_hat,\n")
  cat("the model is correctly identifying them as over-represented in KHANDLE,\n")
  cat("but trimming may prevent the weights from fully correcting this.\n")
}

# Compare untrimmed vs trimmed weights for education3_h=4
cat("\n--- Untrimmed vs trimmed raw weight for education3_h=4 ---\n")
if ("p_hat" %in% names(k1)) {
  k1_edu4 <- k1 %>% filter(education3_h == 4)
  p_khandle_approx <- nrow(k1) / nrow(d1)
  raw_w <- ((1 - k1_edu4$p_hat) / k1_edu4$p_hat) * (p_khandle_approx / (1 - p_khandle_approx))

  cat("N with education3_h=4:", nrow(k1_edu4), "\n")
  cat("Raw (untrimmed) weight summary:\n")
  print(summary(raw_w))
  cat("\nTrimmed weight summary (sw_final):\n")
  print(summary(k1_edu4$sw_final))
  cat("\nProportion of edu=4 weights that were trimmed:\n")
  trimmed <- sum(abs(raw_w - k1_edu4$sw_final) > 0.001, na.rm = TRUE)
  cat(trimmed, "of", nrow(k1_edu4), "(", round(100 * trimmed / nrow(k1_edu4), 1), "%)\n")
}

cat("\n====================================================================\n")
cat("DIAGNOSTIC 6: Logistic model coefficients for education3_h\n")
cat("====================================================================\n")

# Refit the weight model on imp=1 to see coefficients
d1 <- d1 %>% mutate(khandle_h = 1 - brfss_h)

analysis_setup <- detect_analysis_setup(d1)

if (analysis_setup == "t10") {
  weight_formula <- khandle_h ~
    race_summary_f + male_h + age_h + education3_h +
    income_gtmed_pp_h + adl_walking_h + english_interview_h + goodhealth_h
} else {
  weight_formula <- khandle_h ~
    race_summary_f + male_h + age_h + education3_h +
    income_gtmed_pp_h + adl_walking_h + goodhealth_h
}

d1 <- d1 %>%
  mutate(
    male_h = as01_num(male_h),
    english_interview_h = as01_num(english_interview_h),
    adl_walking_h = as01_num(adl_walking_h),
    goodhealth_h = as01_num(goodhealth_h)
  )

fit1 <- glm(
  formula = weight_formula,
  data = d1,
  family = binomial(link = "logit"),
  weights = d1$brfss_sampwt_h
)

cat("\n--- Logistic model summary (imp=1) ---\n")
print(summary(fit1))

cat("\n--- Odds ratios for education terms ---\n")
edu_coefs <- coef(fit1)[grep("education3_h", names(coef(fit1)))]
cat("Coefficients (log-odds):\n")
print(round(edu_coefs, 3))
cat("\nOdds ratios:\n")
print(round(exp(edu_coefs), 3))

cat("\n====================================================================\n")
cat("SUMMARY\n")
cat("====================================================================\n")
cat("\nKey question: Is the education3_h=4 imbalance caused by:\n")
cat("  (a) The logistic model not being flexible enough?\n")
cat("  (b) Weight trimming clipping the correction?\n")
cat("  (c) Fundamental over-representation that can't be fixed without\n")
cat("      extreme weights that would inflate variance?\n")
cat("\nLook at:\n")
cat("  - Diagnostic 1: How big is the education gap?\n")
cat("  - Diagnostic 4: How much does IOSW reduce it?\n")
cat("  - Diagnostic 5: Are weights being heavily trimmed for edu=4?\n")
cat("  - If trimming is the cause, consider relaxing to 0.5%/99.5%\n")
cat("  - If the model is the cause, consider GBM weights via twang\n")
cat("\n03a_weight_diagnostics.R completed.\n")
