#!/usr/bin/env Rscript
# ---------------------------
# Calculates means for mineralization and qPCR data by main factors
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------

################## Setup
# General setup scripts
source("code/setup/setup.R")
library(gt)
library(gtsummary)

# Read in data
n2o_data <- read.csv("data/prepped_data/n2o_data.csv") %>%
  mutate(
    Day = as.factor(Day)
  )
co2_data <- read.csv("data/prepped_data/co2_data.csv") %>%
  mutate(
    Day = as.factor(Day)
  )
nh4_no3_min_data <- read.csv("data/prepped_data/mineralization_data.csv") %>%
  mutate(
    Day = as.factor(Day)
  )

################ Summary statistics of mineralization by AmoA
# Helper function
calc_summary_stats <- function(n_df, n_vars, voi) {
  n_df %>%
    select(
      any_of(c(n_vars)),
      {{ voi }},
    ) %>%
    tbl_summary(
      by = {{ voi }},
      statistic = list(all_continuous() ~ "{mean}")
    ) %>%
    add_p() %>%
    as_tibble() %>%
    clean_names()
}

calc_summary_stats(
  co2_data,
  co2_vars,
  Crop
)

# Generate summary statistics for each main factor
# TODO: Is there a way to do this with gtsummary::tbl_strata?
crop_summary_statistics <- calc_summary_stats(Crop, co2_data)
treatment_summary_statistics <- calc_summary_stats(Treatment)
addition_summary_statistics <- calc_summary_stats(Addition)

# Combine all summary statistics, removing redundant ID columns
all_summaries <- bind_cols(
  addition_summary_statistics,
  crop_summary_statistics %>% select(-characteristic),
  treatment_summary_statistics %>% select(-characteristic)
)


write.table(
  all_summaries,
  "results/stats/qpcr_mineralization_summary_statistics.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
