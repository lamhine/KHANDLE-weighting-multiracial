#------------------------------------------------------------------
# Script: 02_clean_KHANDLE_race.R
# Purpose: Clean KHANDLE race data (single-run via config analysis_tag)
# Author: Tracy Lam-Hine
# Repo: github.com/lamhine/KHANDLE-weighting-multiracial
# Last updated: 2026-01-12
# Notes:
#   - Runs ONE analysis at a time using `analysis_tag` from config/config.R
#       * t10_replication: no free-text available (structured-only)
#       * t16_update: uses free-text + lookup + manual overrides (adds domains only)
#   - Multiracial = 2+ non-Hispanic race domains (Hispanic overrides)
#   - Saves to:
#       processed_data_parent_dir/{analysis_tag}/02_khandle_cleaned_race.rds
#------------------------------------------------------------------

source("config/config.R")

#------------------------------------------------------------------
# Require analysis_tag from config
#------------------------------------------------------------------

stopifnot(exists("analysis_tag"))
stopifnot(is.character(analysis_tag), length(analysis_tag) == 1)
stopifnot(analysis_tag %in% c("t10_replication", "t16_update"))

#------------------------------------------------------------------
# File map
#------------------------------------------------------------------

khandle_files <- c(
  t10_replication = "khandle_all_waves_20210218.sas7bdat",
  t16_update      = "khandle_all_waves_20250110.sas7bdat"
)

khandle_file <- khandle_files[[analysis_tag]]
if (is.null(khandle_file)) stop("No KHANDLE file mapped for analysis_tag: ", analysis_tag)

in_path <- file.path(raw_data_dir, khandle_file)
if (!file.exists(in_path)) stop("KHANDLE raw file not found: ", in_path)

# IMPORTANT: write into the PARENT processed-data folder
out_dir <- file.path(processed_data_parent_dir, analysis_tag)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("------------------------------------------------------------")
message("Running analysis_tag: ", analysis_tag)
message("Reading: ", in_path)

#------------------------------------------------------------------
# Helper: normalize free-text
#------------------------------------------------------------------

normalize_text <- function(x) {
  x %>%
    tidyr::replace_na("") %>%
    stringr::str_squish() %>%
    stringr::str_to_lower()
}

#------------------------------------------------------------------
# Load free-text lookup (used only for t16_update when free-text exists)
#------------------------------------------------------------------

if (!file.exists("free_text_lookup_table.csv")) {
  stop("Missing required file in project root: free_text_lookup_table.csv")
}

text_lookup_raw <- readr::read_csv(
  "free_text_lookup_table.csv",
  show_col_types = FALSE
) %>%
  janitor::clean_names()

# Common Excel artifact: a leading row-number column named "...1"
if (ncol(text_lookup_raw) >= 2 && names(text_lookup_raw)[1] == "...1") {
  text_lookup_raw <- text_lookup_raw %>% dplyr::select(-1)
}

# First remaining column is assumed to be the raw key
key_col <- names(text_lookup_raw)[1]

text_lookup <- text_lookup_raw %>%
  dplyr::rename(text_raw_lookup = !!rlang::sym(key_col)) %>%
  dplyr::mutate(
    text_norm = normalize_text(text_raw_lookup),
    
    add_aian  = as.integer(add_aian),
    add_asian = as.integer(add_asian),
    add_black = as.integer(add_black),
    add_hisp  = as.integer(add_hispanic),
    add_white = as.integer(add_white),
    add_other = as.integer(add_other),
    flag      = as.logical(flag)
  ) %>%
  dplyr::transmute(
    text_norm,
    add_aian  = dplyr::coalesce(add_aian, 0L),
    add_asian = dplyr::coalesce(add_asian, 0L),
    add_black = dplyr::coalesce(add_black, 0L),
    add_hisp  = dplyr::coalesce(add_hisp, 0L),
    add_white = dplyr::coalesce(add_white, 0L),
    add_other = dplyr::coalesce(add_other, 0L),
    flag      = dplyr::coalesce(flag, FALSE)
  ) %>%
  dplyr::filter(text_norm != "") %>%
  dplyr::distinct(text_norm, .keep_all = TRUE)

#------------------------------------------------------------------
# Manual overrides for flagged study IDs (only relevant when free-text exists)
#------------------------------------------------------------------

manual_overrides <- tibble::tribble(
  ~STUDYID,  ~add_aian_manual, ~add_asian_manual, ~add_black_manual, ~add_hisp_manual, ~add_white_manual, ~add_other_manual, ~note,
  "008721",  0L,              0L,               0L,               1L,              1L,               0L, "Mix of both -> Hispanic + White",
  "039651",  0L,              0L,               0L,               0L,              1L,               0L, "Part Indian and white -> add White (yields Black+White given structured)",
  "001441",  0L,              0L,               0L,               0L,              0L,               0L, "Latin American and indian -> keep as structured",
  "020781",  0L,              0L,               0L,               0L,              0L,               0L, "Indian -> keep as structured (South Asian checked)",
  "032641",  0L,              0L,               0L,               0L,              0L,               0L, "Indian + Mexican checked -> keep as Hispanic",
  "049721",  0L,              0L,               0L,               0L,              0L,               0L, "Asian+White structured -> already multiracial",
  "208651",  1L,              0L,               0L,               0L,              1L,               0L, "American Indian and Basque -> reinforce lookup (optional)"
)

#------------------------------------------------------------------
# 1) Define race-domain flags from STRUCTURED W1_ETHNICITY_* vars
#    - Robust to missing W1_ETHNICITY_77 (absent in T10)
#    - Robust to missing W1_ETHNICITY_21_TEXT (absent in T10)
#------------------------------------------------------------------

make_structured_domains <- function(df) {
  
  needed_min <- c(
    "STUDYID",
    paste0("W1_ETHNICITY_", c(1:21, 88, 99)),  # do NOT require 77
    "W1_D_RACE_SUMMARY"
  )
  
  missing_min <- setdiff(needed_min, names(df))
  if (length(missing_min) > 0) {
    stop("Missing required variables in khandle_raw: ", paste(missing_min, collapse = ", "))
  }
  
  has_text <- "W1_ETHNICITY_21_TEXT" %in% names(df)
  has_77   <- "W1_ETHNICITY_77" %in% names(df)
  
  if (!has_77) {
    df <- df %>% dplyr::mutate(W1_ETHNICITY_77 = 0L)
  }
  
  df %>%
    dplyr::mutate(
      has_black_struct = (W1_ETHNICITY_1 == 1 | W1_ETHNICITY_2 == 1 | W1_ETHNICITY_3 == 1),
      
      has_hisp_struct  = (W1_ETHNICITY_4 == 1 | W1_ETHNICITY_5 == 1 | W1_ETHNICITY_6 == 1 |
                            W1_ETHNICITY_7 == 1 | W1_ETHNICITY_8 == 1),
      
      has_asian_struct = (W1_ETHNICITY_9 == 1 | W1_ETHNICITY_10 == 1 | W1_ETHNICITY_11 == 1 |
                            W1_ETHNICITY_12 == 1 | W1_ETHNICITY_13 == 1 | W1_ETHNICITY_14 == 1 |
                            W1_ETHNICITY_15 == 1 | W1_ETHNICITY_16 == 1),
      
      has_aian_struct  = (W1_ETHNICITY_17 == 1),
      
      has_white_struct = (W1_ETHNICITY_18 == 1 | W1_ETHNICITY_19 == 1 | W1_ETHNICITY_20 == 1),
      
      # "Other race domain" is only turned on by lookup/manual additions in this project
      has_other_struct = FALSE,
      
      has_refuse_struct = (W1_ETHNICITY_88 == 1),
      has_dk_struct     = (W1_ETHNICITY_99 == 1),
      
      text_raw  = if (has_text) as.character(W1_ETHNICITY_21_TEXT) else NA_character_,
      text_norm = if (has_text) normalize_text(text_raw) else NA_character_
    )
}

#------------------------------------------------------------------
# 2) Apply free-text lookup + overrides to build FINAL domain flags
#    - Only adds domains (never removes)
#    - If no free-text, returns structured-only domains
#------------------------------------------------------------------

apply_text_additions <- function(df_struct, text_lookup, manual_overrides) {
  
  has_text_norm <- "text_norm" %in% names(df_struct) &&
    any(!is.na(df_struct$text_norm) & df_struct$text_norm != "")
  
  if (!has_text_norm) {
    return(
      df_struct %>%
        dplyr::mutate(
          flag = FALSE,
          has_black = has_black_struct,
          has_hisp  = has_hisp_struct,
          has_asian = has_asian_struct,
          has_aian  = has_aian_struct,
          has_white = has_white_struct,
          has_other = has_other_struct
        )
    )
  }
  
  df_struct %>%
    dplyr::mutate(STUDYID = as.character(STUDYID)) %>%
    dplyr::left_join(text_lookup, by = "text_norm") %>%
    dplyr::left_join(manual_overrides, by = "STUDYID") %>%
    dplyr::mutate(
      flag = dplyr::coalesce(flag, FALSE),
      
      add_aian_lu  = add_aian,
      add_asian_lu = add_asian,
      add_black_lu = add_black,
      add_hisp_lu  = add_hisp,
      add_white_lu = add_white,
      add_other_lu = add_other,
      
      add_aian_manual  = dplyr::coalesce(add_aian_manual, 0L),
      add_asian_manual = dplyr::coalesce(add_asian_manual, 0L),
      add_black_manual = dplyr::coalesce(add_black_manual, 0L),
      add_hisp_manual  = dplyr::coalesce(add_hisp_manual, 0L),
      add_white_manual = dplyr::coalesce(add_white_manual, 0L),
      add_other_manual = dplyr::coalesce(add_other_manual, 0L),
      
      # If flagged, suppress lookup contributions (use manual only)
      add_aian_final  = dplyr::if_else(flag, 0L, dplyr::coalesce(add_aian_lu, 0L))  + as.integer(add_aian_manual),
      add_asian_final = dplyr::if_else(flag, 0L, dplyr::coalesce(add_asian_lu, 0L)) + as.integer(add_asian_manual),
      add_black_final = dplyr::if_else(flag, 0L, dplyr::coalesce(add_black_lu, 0L)) + as.integer(add_black_manual),
      add_hisp_final  = dplyr::if_else(flag, 0L, dplyr::coalesce(add_hisp_lu, 0L))  + as.integer(add_hisp_manual),
      add_white_final = dplyr::if_else(flag, 0L, dplyr::coalesce(add_white_lu, 0L)) + as.integer(add_white_manual),
      add_other_final = dplyr::if_else(flag, 0L, dplyr::coalesce(add_other_lu, 0L)) + as.integer(add_other_manual),
      
      has_black = has_black_struct | (add_black_final > 0),
      has_hisp  = has_hisp_struct  | (add_hisp_final  > 0),
      has_asian = has_asian_struct | (add_asian_final > 0),
      has_aian  = has_aian_struct  | (add_aian_final  > 0),
      has_white = has_white_struct | (add_white_final > 0),
      
      # "Other domain" only becomes TRUE via lookup/manual
      has_other = has_other_struct | (add_other_final > 0)
    )
}

#------------------------------------------------------------------
# 3) Create BRFSS-aligned MECE categorical variables
#    - Hispanic override; multiracial excludes Hispanic
#    - Patch BRFSS-5 edge cases using W1_D_RACE_SUMMARY
#------------------------------------------------------------------

make_mece_race <- function(df) {
  
  # Only include ethnicity indicator vars that actually exist
  eth_ind <- paste0("W1_ETHNICITY_", c(1:21, 77, 88, 99))
  eth_ind <- eth_ind[eth_ind %in% names(df)]
  
  # Build edge trigger
  edge_trigger <- rep(FALSE, nrow(df))
  if ("W1_ETHNICITY_17" %in% names(df)) edge_trigger <- edge_trigger | (df$W1_ETHNICITY_17 == 1)
  if ("W1_ETHNICITY_21" %in% names(df)) edge_trigger <- edge_trigger | (df$W1_ETHNICITY_21 == 1)
  if ("W1_ETHNICITY_88" %in% names(df)) edge_trigger <- edge_trigger | (df$W1_ETHNICITY_88 == 1)
  if ("W1_ETHNICITY_99" %in% names(df)) edge_trigger <- edge_trigger | (df$W1_ETHNICITY_99 == 1)
  
  df %>%
    mutate(
      # (initial) count non-Hispanic race domains
      n_nonhisp_races = (has_black + has_asian + has_aian + has_white + has_other),
      
      race_mece = case_when(
        has_hisp ~ "Hispanic",
        !has_hisp & n_nonhisp_races >= 2 ~ "Multiracial",
        !has_hisp & n_nonhisp_races == 1 & has_white ~ "White",
        !has_hisp & n_nonhisp_races == 1 & has_black ~ "Black",
        !has_hisp & n_nonhisp_races == 1 & has_asian ~ "Asian",
        !has_hisp & n_nonhisp_races == 1 & has_aian  ~ "American Indian or Alaska Native",
        !has_hisp & n_nonhisp_races == 1 & has_other ~ "Other",
        TRUE ~ NA_character_
      ),
      
      race_mece_brfss5 = case_when(
        race_mece %in% c("Hispanic", "White", "Black", "Asian", "Multiracial") ~ race_mece,
        TRUE ~ NA_character_
      ),
      
      # How many ethnicity boxes endorsed (only among columns that exist)
      n_eth_endorsed = if (length(eth_ind) > 0) {
        rowSums(across(all_of(eth_ind)) == 1, na.rm = TRUE)
      } else {
        NA_real_
      },
      
      edge_eligible = is.na(race_mece_brfss5) &
        !is.na(n_eth_endorsed) & (n_eth_endorsed == 1) &
        edge_trigger,
      
      # Map KHANDLE summary to your BRFSS labels
      w1_summary_mapped = case_when(
        W1_D_RACE_SUMMARY == "LatinX" ~ "Hispanic",
        W1_D_RACE_SUMMARY == "White" ~ "White",
        W1_D_RACE_SUMMARY == "Black" ~ "Black",
        W1_D_RACE_SUMMARY == "Asian" ~ "Asian",
        W1_D_RACE_SUMMARY == "Native American" ~ "American Indian or Alaska Native",
        TRUE ~ NA_character_
      ),
      
      # For edge-eligible rows only, turn on the corresponding domain flag
      has_hisp  = if_else(edge_eligible & w1_summary_mapped == "Hispanic", TRUE, has_hisp),
      has_white = if_else(edge_eligible & w1_summary_mapped == "White", TRUE, has_white),
      has_black = if_else(edge_eligible & w1_summary_mapped == "Black", TRUE, has_black),
      has_asian = if_else(edge_eligible & w1_summary_mapped == "Asian", TRUE, has_asian),
      has_aian  = if_else(edge_eligible & w1_summary_mapped == "American Indian or Alaska Native", TRUE, has_aian),
      
      # Recompute domain count AFTER the possible patch
      n_nonhisp_races = as.integer(has_black) + as.integer(has_asian) + as.integer(has_aian) +
        as.integer(has_white) + as.integer(has_other),
      
      # Now define source + patched BRFSS5 using the mapped value
      race_mece_source = case_when(
        !is.na(race_mece_brfss5) ~ "self_report",
        edge_eligible & !is.na(w1_summary_mapped) ~ "khandle_summary",
        TRUE ~ "missing"
      ),
      
      race_mece_brfss5 = case_when(
        !is.na(race_mece_brfss5) ~ race_mece_brfss5,
        edge_eligible ~ w1_summary_mapped,
        TRUE ~ NA_character_
      ),
      
      race_mece = case_when(
        !is.na(race_mece) ~ race_mece,
        edge_eligible ~ w1_summary_mapped,
        TRUE ~ NA_character_
      ),
      
      # IMPORTANT: recompute count at the very end, using final has_* flags
      n_nonhisp_races = (has_black + has_asian + has_aian + has_white + has_other)
    ) %>%
    select(-n_eth_endorsed, -edge_eligible, -w1_summary_mapped)
}

#------------------------------------------------------------------
# Run
#------------------------------------------------------------------

khandle_raw <- haven::read_sas(in_path)

# 1) Structured domains
khandle_domains_struct <- make_structured_domains(khandle_raw)

# 2) Apply free-text additions only for t16_update (and only if text exists)
if (analysis_tag == "t16_update") {
  khandle_domains_aug <- apply_text_additions(
    df_struct        = khandle_domains_struct,
    text_lookup      = text_lookup,
    manual_overrides = manual_overrides
  )
} else {
  khandle_domains_aug <- khandle_domains_struct %>%
    dplyr::mutate(
      flag = FALSE,
      has_black = has_black_struct,
      has_hisp  = has_hisp_struct,
      has_asian = has_asian_struct,
      has_aian  = has_aian_struct,
      has_white = has_white_struct,
      has_other = has_other_struct
    )
}

# 3) MECE coding
khandle_race_final <- make_mece_race(khandle_domains_aug)

khandle_race_final <- khandle_race_final %>%
  mutate(
    has_hisp  = if_else(race_mece_source == "khandle_summary" & race_mece_brfss5 == "Hispanic", TRUE, has_hisp),
    has_white = if_else(race_mece_source == "khandle_summary" & race_mece_brfss5 == "White", TRUE, has_white),
    has_black = if_else(race_mece_source == "khandle_summary" & race_mece_brfss5 == "Black", TRUE, has_black),
    has_asian = if_else(race_mece_source == "khandle_summary" & race_mece_brfss5 == "Asian", TRUE, has_asian),
    has_aian  = if_else(race_mece_source == "khandle_summary" & race_mece_brfss5 == "American Indian or Alaska Native", TRUE, has_aian)
  )

# 4) Diagnostics
message("Diagnostics: race_mece distribution")
print(khandle_race_final %>% dplyr::count(race_mece, sort = TRUE))

message("Diagnostics: race_mece_brfss5 distribution")
print(khandle_race_final %>% dplyr::count(race_mece_brfss5, sort = TRUE))

message("Diagnostics: race_mece_source distribution")
print(khandle_race_final %>% dplyr::count(race_mece_source, sort = TRUE))

if (analysis_tag == "t16_update") {
  message("Diagnostics: free-text lookup match rate (among those with non-empty text)")
  
  diag_lookup <- khandle_race_final %>%
    dplyr::mutate(
      has_text_nonempty = !is.na(text_norm) & text_norm != "",
      matched_lookup = has_text_nonempty & (text_norm %in% text_lookup$text_norm)
    ) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_text_nonempty = sum(has_text_nonempty),
      n_matched_among_text = sum(matched_lookup),
      prop_matched_among_text = dplyr::if_else(n_text_nonempty > 0, n_matched_among_text / n_text_nonempty, NA_real_)
    )
  
  print(diag_lookup)
  
  message("Diagnostics: flagged rows (if any)")
  flagged <- khandle_race_final %>%
    dplyr::filter(flag) %>%
    dplyr::select(
      STUDYID, text_raw, text_norm,
      has_hisp_struct, has_asian_struct, has_black_struct, has_white_struct, has_aian_struct, has_other_struct,
      dplyr::starts_with("add_"), n_nonhisp_races, race_mece, race_mece_brfss5
    ) %>%
    dplyr::arrange(STUDYID)
  
  print(flagged)
}

# 5) Save compact race output for join
khandle_race_out <- khandle_race_final %>%
  dplyr::transmute(
    STUDYID = as.character(STUDYID),
    text_raw,
    has_hisp, has_black, has_asian, has_white, has_aian, has_other,
    n_nonhisp_races,
    race_mece,
    race_mece_brfss5,
    race_mece_source
  )

# 6) Attach to khandle_raw and save full dataset
khandle_cleaned_race <- khandle_raw %>%
  dplyr::mutate(STUDYID = as.character(STUDYID)) %>%
  dplyr::left_join(khandle_race_out, by = "STUDYID")

save_path <- file.path(out_dir, "02_khandle_cleaned_race.rds")
saveRDS(khandle_cleaned_race, file = save_path)

message("Saved: ", save_path)
message("02_clean_KHANDLE_race.R completed successfully.")