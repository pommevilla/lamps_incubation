#!/usr/bin/env Rscript
# ---------------------------
# Prints out summary tables
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------

################## Setup
# General setup scripts
source("code/setup.R")
library(gt)
library(gtsummary)

# Read in data
mineralization_and_amoa_data <- read.csv("data/prepped_data/mineralization_and_qpcr_data.csv")

################ Summary statistics of mineralization by AmoA

# crop_summary_statistics <-

calc_summary_stats <- function(voi) {
  mineralization_and_amoa_data %>%
    select(
      contains(c("ave", "rel", "abs")),
      {{ voi }},
    ) %>%
    tbl_summary(
      by = {{ voi }},
      statistic = list(all_continuous() ~ "{mean}")
    ) %>%
    add_p() %>%
    as_tibble() %>%
    clean_names() %>%
    mutate(
      characteristic = case_when(
        str_detect(characteristic, "ave_") ~ str_replace(characteristic, "ave_", "amoA "),
        characteristic == "f1r2_ave" ~ "F1R2",
        TRUE ~ characteristic
      )
    )
}

crop_summary_statistics <- calc_summary_stats(crop)
treatment_summary_statistics <- calc_summary_stats(treatment)
addition_summary_statistics <- calc_summary_stats(addition)


all_summaries <- bind_cols(
  addition_summary_statistics,
  crop_summary_statistics,
  treatment_summary_statistics
) %>%
  select(-c(6, 10))


write.table(
  all_summaries,
  "results/stats/qpcr_mineralization_summary_statistics.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
