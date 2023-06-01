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
mineralization_and_qpcr_data <- read.csv("data/prepped_data/mineralization_and_qpcr_data.csv")

################ Summary statistics of mineralization by AmoA
# Helper function
calc_summary_stats <- function(voi) {
  mineralization_and_qpcr_data %>%
    select(
      any_of(qpcr_variables),
      {{ voi }},
    ) %>%
    tbl_summary(
      by = {{ voi }},
      statistic = list(all_continuous() ~ "{mean}")
    ) %>%
    add_p() %>%
    as_tibble() %>%
    clean_names()
  #  %>%
  # mutate(
  #   characteristic = case_when(
  #     str_detect(characteristic, "ave_") ~ str_replace(characteristic, "ave_", "amoA "),
  #     characteristic == "f1r2_ave" ~ "F1R2",
  #     TRUE ~ characteristic
  #   )
  # )
}

# Generate summary statistics for each main factor
# TODO: Is there a way to do this with gtsummary::tbl_strata?
crop_summary_statistics <- calc_summary_stats(Crop)
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
