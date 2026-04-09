#------------------------------------------------------------------
# R/utils.R
# Purpose: General-purpose utility functions shared across the pipeline
#------------------------------------------------------------------

#' Convert to numeric, suppressing warnings (handles factors/labelled)
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

#' Safe column getter: returns vector of NAs if column is missing
get_col <- function(df, v) {
  if (v %in% names(df)) df[[v]] else rep(NA, nrow(df))
}

#' Replace specific codes with NA (vectorized)
na_if_in <- function(x, codes) {
  x <- to_num(x)
  x[x %in% codes] <- NA_real_
  x
}

#' Assert that required variables exist in a data frame
assert_has_vars <- function(dat, vars, dat_name) {
  miss <- setdiff(vars, names(dat))
  if (length(miss) > 0) {
    stop(dat_name, " is missing required variable(s): ", paste(miss, collapse = ", "))
  }
  invisible(TRUE)
}

#' Recode BRFSS yes/no-style vars to 1/0
#' BRFSS convention: 1 = Yes, 2 = No. Values in na_codes become NA.
recode_binary01_brfss <- function(x, na_codes = numeric()) {
  x <- to_num(x)
  x[x %in% na_codes] <- NA_real_
  dplyr::case_when(
    x == 1 ~ 1,
    x == 2 ~ 0,
    TRUE ~ NA_real_
  )
}

#' Recode binary variables with flexible value mapping (KHANDLE convention)
recode_binary01_khandle <- function(x, yes = 1, no = 0, no_alt = 2, na_codes = c(88, 99)) {
  x <- to_num(x)
  x[x %in% na_codes] <- NA_real_
  dplyr::case_when(
    x == yes ~ 1,
    x == no  ~ 0,
    x == no_alt ~ 0,
    TRUE ~ NA_real_
  )
}
