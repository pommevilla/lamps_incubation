#!/usr/bin/env Rscript
# ---------------------------
# Preps qPCR data for downstream analysis
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup/setup.R")

qpcr_data <- read.csv(here("data/qpcr", "Incubation_Biomark-qPCR_all_20230525.csv")) %>%
  mutate(Treatment = case_when(
    Treatment == "100N" ~ "112N",
    Treatment == "300N" ~ "336N",
    TRUE ~ Treatment
  )) %>%
  mutate(
    Crop = if_else(
      Crop == "Mxg", "Miscanthus", "Corn"
    ),
    log_norB = log10(norB.001),
    ave_sum = ave_012 + ave_025 + ave_039,
    log_sum = log(ave_sum)
  ) %>%
  rename(
    Day = DNA_DOE,
    ave_norB = norB.001
  )

write.csv(
  qpcr_data,
  here("data/prepped_data", "qpcr_data.csv"),
  row.names = FALSE,
  quote = FALSE
)
