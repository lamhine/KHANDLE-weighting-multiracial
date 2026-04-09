#------------------------------------------------------------------
# Script: 01_clean_BRFSS.R
# Purpose: Clean BRFSS microdata after 00_pre-cleaning_BRFSS.R
# Original author: Taylor Mobley
# Notes: Code provided by Elizabeth Hayes-Larson and adapted by Tracy Lam-Hine
# Repo: github.com/lamhine/KHANDLE-weighting-multiracial
# Last updated: 2026-01-13
#------------------------------------------------------------------

# Load configuration
source("config/config.R")

#------------------------------------------------------------------
# Load pooled BRFSS data (output from 00_pre-cleaning_BRFSS.R)
#------------------------------------------------------------------

brfss_pooled <- readRDS(file.path(processed_data_dir, "01_brfss_pooled.rds"))

#------------------------------------------------------------------
# Utilities: to_num() and recode_binary01_brfss() provided by R/utils.R via config
#------------------------------------------------------------------

#------------------------------------------------------------------
# Keep only variables needed for cleaning/analysis
#------------------------------------------------------------------

vars_keep_required <- c(
  "AGE", "YEAR",
  "male", "raceth",
  "_LLCPWT", "cdc_finalwt",
  "_STSTR", "_PSU"
)

vars_keep_optional <- c(
  # version/track indicators
  "_LCPWTV1", "_LCPWTV2", "_LCPWTV3",
  
  # demographics / SES
  "MARITAL",
  "EDUCA", "EDUCAC",
  "INCOM02", "INCOM03", "INCOME3",
  "INCOMCDC", "INCOMECDC",
  
  # household size inputs
  "HHSIZE", "HHADULT", "CHILDREN",
  
  # employment
  "EMPLOY2", "EMPLOY1",
  
  # county
  "COUNTY1", "CTYCODE2", "CPCOUNTY",
  
  # housing tenure
  "OWNHOME", "RENTHOM1",
  
  # language
  "SPANEN2", "english_interview",
  
  # health insurance
  "HAVEPLN3", "_HLTHPLN", "_HLTHPL1",
  
  # military/veteran
  "MILITAR2", "VETERAN3",
  
  # optional environment/housing context
  "INDOORS", "HOUSTYPE", "OUTWORK",
  
  # health items
  "BLIND",
  "DEAF",
  "GENHLTH",
  "MENTHLTH", "PHYSHLTH",
  
  # functioning
  "DIFFWALK",
  "DIFDRES2", "DIFFDRESS",
  "DIFFERND", "DIFFALON",
  "REMEM2", "DECIDE",
  
  # outcomes/conditions
  "CIM_THNK", "CIMEMLOS", "CIMEMLO1",
  "STROKE2", "CVDSTRK3",
  "ANGINA", "CVDCRHD4",
  "ARTHRITD", "HAVARTH4",
  "KIDNEY", "CHCKDNY2",
  "DIABCOR3", "DIABETE4",
  "DEPRESS1", "ADDEPEV3",
  "EXERANY1", "EXERANY2",
  
  # smoking
  "SMOKE100",
  "SMKEVDA2", "SMOKDAY2",
  
  # anthropometrics
  "HEIGHT", "HEIGHT3",
  "WEIGHT", "WEIGH2"
)

missing_required <- setdiff(vars_keep_required, names(brfss_pooled))
if (length(missing_required) > 0) {
  stop(
    "01_clean_BRFSS.R: These REQUIRED variables are missing from 01_brfss_pooled.rds:\n",
    paste(missing_required, collapse = ", ")
  )
}

vars_keep <- c(vars_keep_required, intersect(vars_keep_optional, names(brfss_pooled)))

clean_data <- brfss_pooled[, vars_keep]
names(clean_data) <- tolower(names(clean_data))

#------------------------------------------------------------------
# Standardize VERSION from _LCPWTV1/_LCPWTV2/_LCPWTV3
#   Prefer V3 > V2 > V1
#------------------------------------------------------------------

clean_data$version <- NA_real_
if ("_lcpwtv1" %in% names(clean_data)) clean_data$version <- to_num(clean_data$`_lcpwtv1`)
if ("_lcpwtv2" %in% names(clean_data)) clean_data$version <- dplyr::coalesce(to_num(clean_data$`_lcpwtv2`), clean_data$version)
if ("_lcpwtv3" %in% names(clean_data)) clean_data$version <- dplyr::coalesce(to_num(clean_data$`_lcpwtv3`), clean_data$version)

#------------------------------------------------------------------
# Income variables
#   - CA income: INCOM02 (2014-2020), INCOM03 (2021), INCOME3 (2022-2023)
#   - CDC income: INCOMCDC vs INCOMECDC
#------------------------------------------------------------------

clean_data$income_cdc <- NA_real_
if ("incomcdc" %in% names(clean_data)) clean_data$income_cdc <- to_num(clean_data$incomcdc)
if ("incomecdc" %in% names(clean_data)) clean_data$income_cdc <- dplyr::coalesce(to_num(clean_data$incomecdc), clean_data$income_cdc)

clean_data$income_raw <- NA_real_
if ("incom02" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2020
  clean_data$income_raw[idx] <- to_num(clean_data$incom02[idx])
}
if ("incom03" %in% names(clean_data)) {
  idx <- clean_data$year == 2021
  clean_data$income_raw[idx] <- to_num(clean_data$incom03[idx])
}
if ("income3" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$income_raw[idx] <- to_num(clean_data$income3[idx])
}

clean_data$income <- dplyr::case_when(
  clean_data$income_raw %in% c(77, 99) ~ NA_real_,
  TRUE ~ clean_data$income_raw
)

#------------------------------------------------------------------
# Education variables
#------------------------------------------------------------------

clean_data$educdc_v2 <- NA_real_
if ("educac" %in% names(clean_data)) {
  clean_data$educdc_v2 <- dplyr::if_else(to_num(clean_data$educac) == 9, NA_real_, to_num(clean_data$educac))
}

clean_data$educa_v2 <- NA_real_
if ("educa" %in% names(clean_data)) {
  e <- to_num(clean_data$educa)
  clean_data$educa_v2 <- dplyr::case_when(
    e == 88 ~ 0,
    e %in% c(9, 77, 99) ~ NA_real_,
    TRUE ~ e
  )
}

clean_data$educa_v2_plus1 <- clean_data$educa_v2 + 1

clean_data$educa_v2_collapse <- dplyr::case_when(
  clean_data$educa_v2_plus1 %in% 5:7 ~ 5,
  clean_data$educa_v2_plus1 %in% c(8, 9) ~ 6,
  TRUE ~ clean_data$educa_v2_plus1
)

clean_data$education    <- clean_data$educdc_v2
clean_data$education_ca <- clean_data$educa_v2_collapse

#------------------------------------------------------------------
# Race summary
#------------------------------------------------------------------

clean_data$race_summary <- dplyr::case_when(
  is.na(clean_data$raceth) ~ NA_real_,
  clean_data$raceth == "Asian" ~ 1,
  clean_data$raceth == "Black" ~ 2,
  clean_data$raceth == "Hispanic" ~ 3,
  clean_data$raceth == "White" ~ 4,
  clean_data$raceth == "Multiracial" ~ 7,
  TRUE ~ 8
)

#------------------------------------------------------------------
# Marital status (MARITAL all years)
#------------------------------------------------------------------

clean_data$marital_status <- NA_real_
if ("marital" %in% names(clean_data)) {
  m <- to_num(clean_data$marital)
  clean_data$marital_status <- dplyr::case_when(
    m %in% c(7, 9, 77, 99) ~ NA_real_,
    TRUE ~ m
  )
}

#------------------------------------------------------------------
# Interview language
#------------------------------------------------------------------

if (!("english_interview" %in% names(clean_data))) clean_data$english_interview <- NA_real_

if (all(is.na(clean_data$english_interview)) && "spanen2" %in% names(clean_data)) {
  s <- to_num(clean_data$spanen2)
  clean_data$english_interview <- dplyr::case_when(
    clean_data$year %in% c(2015, 2016, 2017, 2018) & s == 1 ~ 1,
    clean_data$year == 2014 & s == 2 ~ 1,
    clean_data$year %in% c(2015, 2016, 2017, 2018) & s == 2 ~ 0,
    clean_data$year == 2014 & s == 1 ~ 0,
    TRUE ~ NA_real_
  )
}

#------------------------------------------------------------------
# Employment: EMPLOY2 (2014-2021), EMPLOY1 (2022-2023)
#------------------------------------------------------------------

clean_data$employ_std <- NA_real_
if ("employ2" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$employ_std[idx] <- to_num(clean_data$employ2[idx])
}
if ("employ1" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$employ_std[idx] <- to_num(clean_data$employ1[idx])
}

clean_data$retired <- dplyr::case_when(
  clean_data$employ_std == 7 ~ 1,
  clean_data$employ_std %in% c(9, 77, 99) ~ NA_real_,
  is.na(clean_data$employ_std) ~ NA_real_,
  TRUE ~ 0
)

clean_data$retired_v2 <- dplyr::case_when(
  clean_data$year == 2015 & clean_data$employ_std == 77 ~ 1,
  TRUE ~ clean_data$retired
)

clean_data$working <- dplyr::case_when(
  clean_data$employ_std %in% c(1, 2) ~ 1,
  clean_data$year == 2015 & clean_data$employ_std == 77 ~ 0,
  clean_data$employ_std %in% c(9, 77, 99) ~ NA_real_,
  is.na(clean_data$employ_std) ~ NA_real_,
  TRUE ~ 0
)

#------------------------------------------------------------------
# County FIPS
#   - 2014–2022: COUNTY1
#   - 2023: CPCOUNTY
#   (CTYCODE2 is kept in vars_keep_optional for compatibility, but not used here)
#------------------------------------------------------------------

clean_data$county_raw <- NA_real_

if ("county1" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2022
  clean_data$county_raw[idx] <- to_num(clean_data$county1[idx])
}

if ("cpcounty" %in% names(clean_data)) {
  idx <- clean_data$year == 2023
  clean_data$county_raw[idx] <- to_num(clean_data$cpcounty[idx])
}

clean_data$county_code <- dplyr::case_when(
  # match prior behavior: treat these as missing where they occur
  clean_data$county_raw %in% c(777, 888, 999) ~ NA_real_,
  TRUE ~ clean_data$county_raw
)

#------------------------------------------------------------------
# Home ownership: OWNHOME (2014-2021), RENTHOM1 (2022-2023)
#------------------------------------------------------------------

clean_data$ownhome_std <- NA_real_
if ("ownhome" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$ownhome_std[idx] <- to_num(clean_data$ownhome[idx])
}
if ("renthom1" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$ownhome_std[idx] <- to_num(clean_data$renthom1[idx])
}

clean_data$ownhome_v2 <- dplyr::case_when(
  clean_data$ownhome_std %in% c(7, 9, 77, 99) ~ NA_real_,
  TRUE ~ clean_data$ownhome_std
)

clean_data$ownhome_v3 <- dplyr::case_when(
  clean_data$ownhome_std %in% c(7, 9, 77, 99) ~ NA_real_,
  clean_data$ownhome_std == 1 ~ 1,
  clean_data$ownhome_std %in% c(2, 3) ~ 0,
  TRUE ~ NA_real_
)

#------------------------------------------------------------------
# Household size
#   Use HHSIZE only when provided by BRFSS.
#   Do NOT reconstruct from HHADULT + CHILDREN.
#------------------------------------------------------------------

clean_data$hhsize_std <- NA_real_

if ("hhsize" %in% names(clean_data)) {
  h <- to_num(clean_data$hhsize)
  h[h %in% c(77, 99)] <- NA_real_
  h[h < 1 | h > 20] <- NA_real_
  clean_data$hhsize_std <- h
}

#------------------------------------------------------------------
# Health variables
#------------------------------------------------------------------

# General health: GENHLTH
clean_data$health <- NA_real_
if ("genhlth" %in% names(clean_data)) {
  g <- to_num(clean_data$genhlth)
  g[g %in% c(7, 9, 77, 99)] <- NA_real_
  clean_data$health <- g
}

# Mental health days
clean_data$mental_health <- NA_real_
if ("menthlth" %in% names(clean_data)) {
  m <- to_num(clean_data$menthlth)
  m[m %in% c(77, 99)] <- NA_real_
  m[m == 88] <- 0
  clean_data$mental_health <- m
}

# Physical health days
clean_data$physical_health <- NA_real_
if ("physhlth" %in% names(clean_data)) {
  p <- to_num(clean_data$physhlth)
  p[p %in% c(77, 99)] <- NA_real_
  p[p == 88] <- 0
  clean_data$physical_health <- p
}

# Height / weight (year-specific names) -- safe when columns are missing
clean_data$height_std <- NA_real_
clean_data$weight_std <- NA_real_

if ("height" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$height_std[idx] <- to_num(clean_data$height[idx])
}
if ("height3" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$height_std[idx] <- to_num(clean_data$height3[idx])
}
if ("weight" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$weight_std[idx] <- to_num(clean_data$weight[idx])
}
if ("weigh2" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$weight_std[idx] <- to_num(clean_data$weigh2[idx])
}

clean_data$height_v2 <- dplyr::case_when(
  clean_data$height_std %in% c(7777, 9999) ~ NA_real_,
  TRUE ~ clean_data$height_std
)
clean_data$weight_v2 <- dplyr::case_when(
  clean_data$weight_std %in% c(7777, 9999) ~ NA_real_,
  TRUE ~ clean_data$weight_std
)

# Diabetes (year-specific names) -- safe
clean_data$diabetes_std <- NA_real_
if ("diabcor3" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$diabetes_std[idx] <- to_num(clean_data$diabcor3[idx])
}
if ("diabete4" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$diabetes_std[idx] <- to_num(clean_data$diabete4[idx])
}

clean_data$diabetes_v2 <- dplyr::case_when(
  clean_data$diabetes_std %in% c(7, 9, 77, 99) ~ NA_real_,
  TRUE ~ clean_data$diabetes_std
)
clean_data$diabetes_v3 <- dplyr::case_when(
  clean_data$diabetes_std %in% c(7, 9, 77, 99) ~ NA_real_,
  clean_data$diabetes_std == 1 ~ 1,
  clean_data$diabetes_std == 2 ~ 0,
  TRUE ~ NA_real_
)

# Smoking (SMOKE100 all years; current smoking differs)
clean_data$smoke <- NA_real_
if ("smoke100" %in% names(clean_data)) {
  s <- to_num(clean_data$smoke100)
  s[s %in% c(7, 9, 77, 99)] <- NA_real_
  clean_data$smoke <- dplyr::case_when(
    s == 1 ~ 1,
    s == 2 ~ 0,
    TRUE ~ NA_real_
  )
}

clean_data$smoke_now_std <- NA_real_
if ("smkevda2" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$smoke_now_std[idx] <- to_num(clean_data$smkevda2[idx])
}
if ("smokday2" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$smoke_now_std[idx] <- to_num(clean_data$smokday2[idx])
}

clean_data$smoke_now <- dplyr::case_when(
  clean_data$smoke_now_std %in% c(7, 9, 99) ~ NA_real_,
  TRUE ~ clean_data$smoke_now_std
)

clean_data$smoke_status <- dplyr::case_when(
  to_num(clean_data$smoke100) == 1 & is.na(clean_data$smoke_now_std) ~ 1,
  to_num(clean_data$smoke100) == 1 & clean_data$smoke_now_std == 1 ~ 2,
  to_num(clean_data$smoke100) == 1 ~ 1,
  to_num(clean_data$smoke100) == 2 ~ 0,
  TRUE ~ NA_real_
)

clean_data$smoke_now_v2 <- dplyr::case_when(
  clean_data$smoke == 0 & is.na(clean_data$smoke_now) ~ 0,
  TRUE ~ clean_data$smoke_now
)

#------------------------------------------------------------------
# YES/NO-format variables: safe harmonization + recode to 1/0
#------------------------------------------------------------------

# Military/Veteran
clean_data$military_raw <- NA_real_
if ("militar2" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$military_raw[idx] <- to_num(clean_data$militar2[idx])
}
if ("veteran3" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$military_raw[idx] <- to_num(clean_data$veteran3[idx])
}

clean_data$military <- dplyr::case_when(
  clean_data$year %in% 2014:2021 ~ recode_binary01_brfss(clean_data$military_raw, na_codes = c(3, 66, 7, 77, 8, 88, 9, 99, 96)),
  clean_data$year %in% 2022:2023 ~ recode_binary01_brfss(clean_data$military_raw, na_codes = c(7)),
  TRUE ~ NA_real_
)

# Insurance
clean_data$insured_raw <- NA_real_
if ("havepln3" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2020
  clean_data$insured_raw[idx] <- to_num(clean_data$havepln3[idx])
}
if ("_hlthpln" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2021:2022
  clean_data$insured_raw[idx] <- to_num(clean_data$`_hlthpln`[idx])
}
if ("_hlthpl1" %in% names(clean_data)) {
  idx <- clean_data$year == 2023
  clean_data$insured_raw[idx] <- to_num(clean_data$`_hlthpl1`[idx])
}

clean_data$insured <- dplyr::case_when(
  clean_data$year %in% 2014:2020 ~ recode_binary01_brfss(clean_data$insured_raw, na_codes = c(7, 77, 9, 99)),
  clean_data$year %in% 2021:2023 ~ recode_binary01_brfss(clean_data$insured_raw, na_codes = c(9)),
  TRUE ~ NA_real_
)

# Physical activity
clean_data$phys_activity_raw <- NA_real_
if ("exerany1" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$phys_activity_raw[idx] <- to_num(clean_data$exerany1[idx])
}
if ("exerany2" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$phys_activity_raw[idx] <- to_num(clean_data$exerany2[idx])
}
clean_data$phys_activity <- recode_binary01_brfss(clean_data$phys_activity_raw, na_codes = c(7, 9, 77, 99))

# Vision impairment (blind)
clean_data$vision_impair_raw <- NA_real_
if ("blind" %in% names(clean_data)) clean_data$vision_impair_raw <- to_num(clean_data$blind)
clean_data$vision_impair <- recode_binary01_brfss(clean_data$vision_impair_raw, na_codes = c(7, 9, 77, 99))

# Difficulty walking
clean_data$diff_walk_raw <- NA_real_
if ("diffwalk" %in% names(clean_data)) clean_data$diff_walk_raw <- to_num(clean_data$diffwalk)
clean_data$diff_walk <- recode_binary01_brfss(clean_data$diff_walk_raw, na_codes = c(7, 9, 77, 99))

# Difficulty dressing
clean_data$diff_dress_raw <- NA_real_
if ("difdres2" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$diff_dress_raw[idx] <- to_num(clean_data$difdres2[idx])
}
if ("diffdress" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$diff_dress_raw[idx] <- to_num(clean_data$diffdress[idx])
}
clean_data$diff_dress <- recode_binary01_brfss(clean_data$diff_dress_raw, na_codes = c(7, 9, 77, 99))

# Difficulty errands
clean_data$diff_errands_raw <- NA_real_
if ("differnd" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$diff_errands_raw[idx] <- to_num(clean_data$differnd[idx])
}
if ("diffalon" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$diff_errands_raw[idx] <- to_num(clean_data$diffalon[idx])
}
clean_data$diff_errands <- recode_binary01_brfss(clean_data$diff_errands_raw, na_codes = c(7, 9, 77, 99))

# Difficulty remember
clean_data$diff_remember_raw <- NA_real_
if ("remem2" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$diff_remember_raw[idx] <- to_num(clean_data$remem2[idx])
}
if ("decide" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$diff_remember_raw[idx] <- to_num(clean_data$decide[idx])
}
clean_data$diff_remember <- recode_binary01_brfss(clean_data$diff_remember_raw, na_codes = c(7, 9, 77, 99))

# Stroke
clean_data$stroke_raw <- NA_real_
if ("stroke2" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$stroke_raw[idx] <- to_num(clean_data$stroke2[idx])
}
if ("cvdstrk3" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$stroke_raw[idx] <- to_num(clean_data$cvdstrk3[idx])
}
clean_data$stroke_v2 <- recode_binary01_brfss(clean_data$stroke_raw, na_codes = c(7, 9, 77, 99))

# Angina / CHD
clean_data$angina_raw <- NA_real_
if ("angina" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$angina_raw[idx] <- to_num(clean_data$angina[idx])
}
if ("cvdcrhd4" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$angina_raw[idx] <- to_num(clean_data$cvdcrhd4[idx])
}
clean_data$angina_v2 <- recode_binary01_brfss(clean_data$angina_raw, na_codes = c(7, 9, 77, 99))

# Arthritis
clean_data$arthritis_raw <- NA_real_
if ("arthritd" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$arthritis_raw[idx] <- to_num(clean_data$arthritd[idx])
}
if ("havarth4" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$arthritis_raw[idx] <- to_num(clean_data$havarth4[idx])
}
clean_data$arthritis_v2 <- recode_binary01_brfss(clean_data$arthritis_raw, na_codes = c(7, 9, 77, 99))

# Kidney
clean_data$kidney_raw <- NA_real_
if ("kidney" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$kidney_raw[idx] <- to_num(clean_data$kidney[idx])
}
if ("chckdny2" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$kidney_raw[idx] <- to_num(clean_data$chckdny2[idx])
}
clean_data$kidney_v2 <- recode_binary01_brfss(clean_data$kidney_raw, na_codes = c(7, 9, 77, 99))

# Depression
clean_data$depress_raw <- NA_real_
if ("depress1" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2014:2021
  clean_data$depress_raw[idx] <- to_num(clean_data$depress1[idx])
}
if ("addepev3" %in% names(clean_data)) {
  idx <- clean_data$year %in% 2022:2023
  clean_data$depress_raw[idx] <- to_num(clean_data$addepev3[idx])
}
clean_data$depress_v2 <- recode_binary01_brfss(clean_data$depress_raw, na_codes = c(7, 9, 77, 99))

#------------------------------------------------------------------
# Optional: Cognitive module item (keep raw harmonized)
#   CIM_THNK (2015, 2020), CIMEMLOS (2022), CIMEMLO1 (2023)
#------------------------------------------------------------------

clean_data$cim_thnk_std <- NA_real_
if ("cim_thnk" %in% names(clean_data)) {
  idx <- clean_data$year %in% c(2015, 2020)
  clean_data$cim_thnk_std[idx] <- to_num(clean_data$cim_thnk[idx])
}
if ("cimemlos" %in% names(clean_data)) {
  idx <- clean_data$year == 2022
  clean_data$cim_thnk_std[idx] <- to_num(clean_data$cimemlos[idx])
}
if ("cimemlo1" %in% names(clean_data)) {
  idx <- clean_data$year == 2023
  clean_data$cim_thnk_std[idx] <- to_num(clean_data$cimemlo1[idx])
}

#------------------------------------------------------------------
# Save cleaned dataset to processed_data_dir
#------------------------------------------------------------------

saveRDS(clean_data, file = file.path(processed_data_dir, "02_brfss_clean.rds"))
message("Saved cleaned BRFSS dataset: ", file.path(processed_data_dir, "02_brfss_clean.rds"))
message("01_clean_BRFSS.R completed successfully.")