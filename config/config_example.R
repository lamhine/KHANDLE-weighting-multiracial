# config_example.R
# Purpose: Define user-specific settings for KHANDLE replication project
# NOTE: Save this file as 'config.R' and update the paths below to match your environment

# ------------------------------------------------------------------
# Set secure PHI-compliant paths (external to GitHub repo)
# ------------------------------------------------------------------

# Replace these paths with your own Box or local storage directories
raw_data_dir <- "/path/to/your/Box/KHANDLE/01_raw"
processed_data_dir <- "/path/to/your/Box/KHANDLE/02_processed"

# Results directory is assumed to live inside the GitHub repo
results_dir <- here::here("02_results")

# ------------------------------------------------------------------
# Validate directory paths
# ------------------------------------------------------------------

if (!dir.exists(raw_data_dir)) {
  stop("ERROR: The raw data directory does not exist. Update 'config.R' with the correct path.")
}
if (!dir.exists(processed_data_dir)) {
  dir.create(processed_data_dir, recursive = TRUE)
  message("Created missing processed data directory: ", processed_data_dir)
}
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
  message("Created missing results directory: ", results_dir)
}

# Confirm paths
message("Using raw data directory: ", raw_data_dir)
message("Using processed data directory: ", processed_data_dir)
message("Using results directory: ", results_dir)

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
  haven,
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

options(scipen = 999)                    # Avoid scientific notation
options(survey.lonely.psu = "adjust")   # Handle single PSU strata in survey analysis