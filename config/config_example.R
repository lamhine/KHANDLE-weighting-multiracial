# config_example.R
# Purpose: Define user-specific settings for KHANDLE replication project
# NOTE: Copy this file to 'config.R' and update the paths below.
#       config.R is gitignored to keep personal paths out of version control.

# ------------------------------------------------------------------
# Set secure Box paths (external to GitHub repo)
# ------------------------------------------------------------------

raw_data_dir <- "/path/to/your/Box/KHANDLE/01_raw"

# Parent folder containing BOTH pipelines:
#   - t16_update
#   - t10_replication
processed_data_parent_dir <- "/path/to/your/Box/KHANDLE/02_processed"

# Choose which analysis most scripts should run by default (t10_replication or t16_update).
analysis_tag <- "t16_update"

# Active processed data directory used by scripts that read/write a single pipeline
processed_data_dir <- file.path(processed_data_parent_dir, analysis_tag)

# Results directory
results_parent_dir <- here::here("03_results")
results_dir <- results_parent_dir

# ------------------------------------------------------------------
# Validate directory paths
# ------------------------------------------------------------------

if (!dir.exists(raw_data_dir)) {
  stop("ERROR: The raw data directory does not exist. Update 'config.R' with the correct path.")
}

if (!dir.exists(processed_data_parent_dir)) {
  dir.create(processed_data_parent_dir, recursive = TRUE)
  message("Created missing processed data parent directory: ", processed_data_parent_dir)
}

if (!dir.exists(processed_data_dir)) {
  dir.create(processed_data_dir, recursive = TRUE)
  message("Created missing processed data directory: ", processed_data_dir)
}

if (!dir.exists(results_parent_dir)) {
  dir.create(results_parent_dir, recursive = TRUE)
  message("Created missing results parent directory: ", results_parent_dir)
}

# Confirm paths
message("Using raw data directory: ", raw_data_dir)
message("Using processed data parent directory: ", processed_data_parent_dir)
message("Using processed data directory: ", processed_data_dir)
message("Using results directory: ", results_dir)
message("Active analysis_tag: ", analysis_tag)

# ------------------------------------------------------------------
# Load required packages
# ------------------------------------------------------------------

if (!require("pacman")) {
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}

pacman::p_load(
  rstudioapi,
  here,
  tidyverse,
  stringr,
  haven,
  rlang,
  data.table,
  survey,
  srvyr,
  janitor,
  summarytools,
  tableone,
  mice,
  openxlsx,
  twang,
  Hmisc,
  spatstat,
  boot,
  cowplot,
  ggpubr
)

# ------------------------------------------------------------------
# Global options
# ------------------------------------------------------------------

set.seed(12345)                  # reproducibility
options(scipen = 999)            # Avoid scientific notation
options(survey.lonely.psu = "adjust")  # Handle single PSU strata

# ------------------------------------------------------------------
# Source shared utility files
# ------------------------------------------------------------------

source(here::here("R/utils.R"))
source(here::here("R/race_constants.R"))
source(here::here("R/analysis_helpers.R"))
source(here::here("R/balance_helpers.R"))
