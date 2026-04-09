#------------------------------------------------------------------
# Script: 00_pre-cleaning_BRFSS.R
# Purpose: Prepping CA BRFSS for data cleaning script (dual workflow)
# Original author: Tracy Lam-Hine
# Notes: Adapted from SAS code written and provided by Taylor Mobley
# Last updated: 2026-01-12
#------------------------------------------------------------------

source("config/config.R")

#####################################################
# File map
#####################################################

brfss_files <- list(
  t10_replication = list(
    "2014" = "brfs_14pr1.sas7bdat",
    "2015" = "brfs_15pr1.sas7bdat",
    "2016" = "brfs_16pr1.sas7bdat",
    "2017" = "brfs_17pr2.sas7bdat",
    "2018" = "brfs_18pr1.sas7bdat"
  ),
  t16_update = list(
    "2019" = "brfs_19p4.sas7bdat",
    "2020" = "brfs_20pr2.sas7bdat",
    "2021" = "brfs_21p.sas7bdat",
    "2022" = "BRFS_22core.sas7bdat",
    "2023" = "brfs_23core_r.sas7bdat"
  )
)

#####################################################
# Determine which tags to run
# - Prefer analysis_tags if config defines it
# - Otherwise fall back to single analysis_tag
#####################################################

if (exists("analysis_tags")) {
  tags_to_run <- as.character(analysis_tags)
} else if (exists("analysis_tag")) {
  tags_to_run <- as.character(analysis_tag)
} else {
  stop("Neither analysis_tags nor analysis_tag found in config/config.R")
}

# Keep only tags that exist in brfss_files (prevents silent typos)
tags_to_run <- tags_to_run[tags_to_run %in% names(brfss_files)]
if (length(tags_to_run) == 0) {
  stop("No valid analysis tag(s) to run. Expected one of: ",
       paste(names(brfss_files), collapse = ", "))
}

#####################################################
# Minimal utilities (NOT for variable selection)
#####################################################

# to_num() is provided by R/utils.R via config

harmonize_types_for_bind <- function(dfs) {
  all_names <- sort(unique(unlist(lapply(dfs, names))))
  class_map <- lapply(all_names, function(v) {
    classes <- unique(unlist(lapply(dfs, function(d) {
      if (v %in% names(d)) class(d[[v]])[1] else NA_character_
    })))
    classes[!is.na(classes)]
  })
  names(class_map) <- all_names
  
  conflicting <- names(class_map)[vapply(class_map, function(x) length(unique(x)) > 1, logical(1))]
  if (length(conflicting) == 0) return(dfs)
  
  message("Harmonizing types for bind_rows (coercing to character): ",
          paste(conflicting, collapse = ", "))
  
  lapply(dfs, function(d) {
    for (v in conflicting) {
      if (v %in% names(d)) d[[v]] <- as.character(d[[v]])
    }
    d
  })
}

#####################################################
# Main loop over analysis tags
#####################################################

for (tag in tags_to_run) {
  
  message("------------------------------------------------------------")
  message("Running BRFSS pre-cleaning for analysis_tag: ", tag)
  
  file_map <- brfss_files[[tag]]
  if (is.null(file_map) || length(file_map) == 0) stop("No BRFSS file map found for: ", tag)
  
  out_dir <- file.path(processed_data_parent_dir, tag)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  #####################################################
  # PART 1: Read each year + create standardized fields
  #####################################################
  
  brfss_list <- purrr::imap(file_map, function(fname, year_chr) {
    
    in_path <- file.path(raw_data_dir, fname)
    if (!file.exists(in_path)) stop("Missing BRFSS raw file: ", in_path)
    
    df <- haven::read_sas(in_path)
    
    # Standardize names to uppercase
    names(df) <- toupper(names(df))
    
    # Add YEAR early (makes debugging easier)
    df <- df %>% dplyr::mutate(YEAR = as.integer(year_chr))
    y <- df$YEAR[1]
    
    # ----------------------------
    # REQUIRED VARIABLES (hard-coded)
    # ----------------------------
    
    if (!("AGE" %in% names(df))) stop("Expected AGE not found for YEAR=", y)
    
    # Sex: SEX (2014-15), SEX1 (2016-17), SEX2 (2018-21), SEXVAR (2022-23)
    sex_var <- dplyr::case_when(
      y %in% c(2014, 2015) ~ "SEX",
      y %in% c(2016, 2017) ~ "SEX1",
      y %in% 2018:2021     ~ "SEX2",
      y %in% 2022:2023     ~ "SEXVAR",
      TRUE ~ NA_character_
    )
    if (is.na(sex_var) || !(sex_var %in% names(df))) {
      stop("Expected sex var (", sex_var, ") not found for YEAR=", y)
    }
    
    if (!("_HISPANC" %in% names(df))) stop("Expected _HISPANC not found for YEAR=", y)
    
    # Race: _MRACE1 (all years except 2022), _MRACE2 (2022)
    race_var <- if (y == 2022) "_MRACE2" else "_MRACE1"
    if (!(race_var %in% names(df))) {
      stop("Expected race var (", race_var, ") not found for YEAR=", y)
    }
    
    if (!("_LLCPWT" %in% names(df))) stop("Expected _LLCPWT not found for YEAR=", y)
    
    # ----------------------------
    # Create standardized variables
    # ----------------------------
    
    df <- df %>%
      dplyr::mutate(
        AGE_NUM = to_num(.data[["AGE"]]),
        
        male = dplyr::case_when(
          to_num(.data[[sex_var]]) == 1 ~ 1,
          to_num(.data[[sex_var]]) == 2 ~ 0,
          TRUE ~ NA_real_
        ),
        
        HISP_FLAG = dplyr::case_when(
          to_num(.data[["_HISPANC"]]) == 1 ~ 1,
          to_num(.data[["_HISPANC"]]) == 2 ~ 0,
          TRUE ~ NA_real_
        ),
        
        # NOTE: For 2022, _MRACE2 uses a different scheme:
        #   - 6 = Multiracial (and "Other" effectively drops out)
        # For all other years:
        #   - 7 = Multiracial
        #   - 6 = Other
        MRACE_CODE = to_num(.data[[race_var]])
      )
    
    # english_interview:
    # - used ONLY for t10_replication (2014–2018)
    # - set to NA for t16_update
    if (tag == "t10_replication") {
      
      if (!("SPANEN2" %in% names(df))) stop("Expected SPANEN2 not found for YEAR=", y)
      s <- to_num(df[["SPANEN2"]])
      
      df <- df %>%
        dplyr::mutate(
          english_interview = dplyr::case_when(
            y %in% c(2015, 2016, 2017, 2018) & s == 1 ~ 1,
            y == 2014 & s == 2 ~ 1,
            y %in% c(2015, 2016, 2017, 2018) & s == 2 ~ 0,
            y == 2014 & s == 1 ~ 0,
            TRUE ~ NA_real_
          )
        )
      
    } else {
      
      df <- df %>% dplyr::mutate(english_interview = NA_real_)
      
    }
    
    # Keep DATE stable if present
    if ("DATE" %in% names(df)) df[["DATE"]] <- as.character(df[["DATE"]])
    
    df
  })
  
  # Harmonize types across years (necessary for bind_rows stability)
  brfss_list <- harmonize_types_for_bind(brfss_list)
  
  #####################################################
  # PART 1B: Scale CDC weights by pooled-year proportion
  #####################################################
  
  year_sample_sizes <- purrr::map_dfr(
    brfss_list,
    ~ tibble::tibble(YEAR = unique(.x$YEAR), n = nrow(.x))
  ) %>%
    dplyr::mutate(proportion = n / sum(n))
  
  brfss_list_weighted <- purrr::map2(
    brfss_list,
    year_sample_sizes$proportion,
    function(df, prop) {
      df %>%
        dplyr::mutate(cdc_finalwt = to_num(.data[["_LLCPWT"]]) * prop)
    }
  )
  
  brfss_pooled <- dplyr::bind_rows(brfss_list_weighted)
  
  #####################################################
  # PART 2: Apply KHANDLE weighting restriction criteria
  #####################################################
  
  brfss_pooled <- brfss_pooled %>%
    dplyr::filter(AGE_NUM >= 65)
  
  brfss_pooled <- brfss_pooled %>%
    dplyr::mutate(
      raceth = dplyr::case_when(
        HISP_FLAG == 1 ~ "Hispanic",
        
        HISP_FLAG == 0 & MRACE_CODE == 1 ~ "White",
        HISP_FLAG == 0 & MRACE_CODE == 2 ~ "Black",
        HISP_FLAG == 0 & MRACE_CODE == 4 ~ "Asian",
        
        HISP_FLAG == 0 & YEAR == 2022 & MRACE_CODE == 6 ~ "Multiracial",
        HISP_FLAG == 0 & YEAR != 2022 & MRACE_CODE == 7 ~ "Multiracial",
        
        MRACE_CODE %in% c(77, 88, 99) ~ NA_character_,
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(raceth %in% c("Hispanic", "White", "Black", "Asian", "Multiracial"))
  
  #####################################################
  # PART 3: Data checks (light + relevant)
  #####################################################
  
  brfss_pooled %>%
    dplyr::group_by(YEAR) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      missing_age = sum(is.na(AGE_NUM)),
      missing_male = sum(is.na(male)),
      missing_raceth = sum(is.na(raceth)),
      missing_cdc_final = sum(is.na(cdc_finalwt)),
      missing_english_interview = sum(is.na(english_interview)),
      .groups = "drop"
    ) %>%
    print()
  
  brfss_pooled %>%
    dplyr::summarise(
      cdc_mean = mean(cdc_finalwt, na.rm = TRUE),
      cdc_sd   = stats::sd(cdc_finalwt, na.rm = TRUE)
    ) %>%
    print()
  
  if (tag == "t10_replication") {
    brfss_pooled %>%
      dplyr::count(YEAR, english_interview) %>%
      dplyr::arrange(YEAR, english_interview) %>%
      print(n = Inf)
  }
  
  if (tag == "t16_update") {
    brfss_pooled %>%
      dplyr::filter(YEAR == 2022) %>%
      dplyr::count(MRACE_CODE, raceth) %>%
      dplyr::arrange(MRACE_CODE, raceth) %>%
      print(n = Inf)
  }
  
  #####################################################
  # SAVE FILES (namespaced by analysis tag)
  #####################################################
  
  save_path <- file.path(out_dir, "01_brfss_pooled.rds")
  saveRDS(brfss_pooled, save_path)
  message("Saved pooled analytic BRFSS dataset, 65+, restricted raceth: ", save_path)
}