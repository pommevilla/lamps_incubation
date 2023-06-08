#!/usr/bin/env Rscript
# ---------------------------
# Preps qPCR data for downstream analysis
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup/setup.R")

amoa_data <- read.csv(here("data/qpcr", "Incubation_Biomark-qPCR_all_20230525.csv")) %>%
  mutate(
    Crop = if_else(Crop == "Mxg", "Miscanthus", "Corn"),
    ave_sum = ave_012 + ave_025 + ave_039,
  ) %>%
  rename(
    Day = DNA_DOE,
  ) %>%
  # We drop this from here since it'll show up in the norB dataset.
  select(-norB.001)

# We keep Crop, Treatment, Addition, and Day to check that these columns
# match between datasets.
norb_data <- read.csv(
  here::here("data/qpcr", "overall_norB.csv")
) %>%
  select(
    Sample.Name, Crop, Treatment, Addition,
    Day = DNA_DOE,
    norB.001, norB.006, norB.039, norB.sum, cnorB
  ) %>%
  mutate(
    Crop = if_else(
      Crop == "Mxg", "Miscanthus", "Corn"
    ),
    across(
      contains("norB"),
      log,
      .names = "log_{.col}"
    )
  )

# After doing checks, we can drop the extra columns from above.
qpcr_data <- left_join(
  amoa_data,
  norb_data %>%
    select(-c(Crop, Treatment, Addition, Day)),
  by = "Sample.Name"
) %>%
  mutate(Treatment = case_when(
    Treatment == "100N" ~ "112N",
    Treatment == "300N" ~ "336N",
    TRUE ~ Treatment
  ))

write.csv(
  qpcr_data,
  here("data/prepped_data", "qpcr_data.csv"),
  row.names = FALSE,
  quote = FALSE
)
