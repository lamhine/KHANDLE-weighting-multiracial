#------------------------------------------------------------------
# R/race_constants.R
# Purpose: Race/ethnicity coding constants and helpers
#------------------------------------------------------------------

# Canonical 5-category race/ethnicity labels used throughout the pipeline
ALLOWED_RACE5 <- c("Asian", "Black", "Hispanic", "White", "Multiracial")

# Reference group for PR/PD comparisons
REF_RACE <- "White"

#' Map race/ethnicity labels to numeric codes
#' 1=Asian, 2=Black, 3=Hispanic, 4=White, 5=Multiracial
race_to_code <- function(x_chr) {
  dplyr::case_when(
    is.na(x_chr) ~ NA_real_,
    x_chr == "Asian" ~ 1,
    x_chr == "Black" ~ 2,
    x_chr == "Hispanic" ~ 3,
    x_chr == "White" ~ 4,
    x_chr == "Multiracial" ~ 5,
    TRUE ~ NA_real_
  )
}

#' Add race_summary_f factor column from race_summary_h numeric codes
add_race_factor <- function(df) {
  if (!("race_summary_h" %in% names(df))) stop("Missing race_summary_h.")
  df %>%
    dplyr::mutate(
      race_summary_f = factor(
        race_summary_h,
        levels = 1:5,
        labels = ALLOWED_RACE5
      )
    )
}

#' Create age categories for standardization: 65-<75, 75-<85, 85+
make_agecat <- function(age) {
  dplyr::case_when(
    !is.finite(age) ~ NA_integer_,
    age >= 65 & age < 75 ~ 1L,
    age >= 75 & age < 85 ~ 2L,
    age >= 85 ~ 3L,
    TRUE ~ NA_integer_
  )
}
