#!/usr/bin/env Rscript
# ---------------------------
# Calculates net mineralization rates for the samples
# Outputs a CSV for mineralization calculations and another
# for mineralization joined with qPCR data.
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup.R")

################## Read in data
# The different calculations for ammonium and nitrate are contained in the first sheet,
# while the day of extraction is contained in the Sample Info sheet.
# We read in both and join them together so we can calculate stuff.
mineralization_data <- read_xlsx("data/FILE_4247_sjh.xlsx", sheet = "Sheet1") %>%
  clean_names() %>%
  mutate(
    crop = str_to_lower(crop),
    treatment = as.factor(case_when(
      treatment == "100N" ~ "112N",
      treatment == "300N" ~ "336N",
      TRUE ~ treatment
    ))
  ) %>%
  separate(initial_name, c("sample", "timepoint"), sep = -1, remove = FALSE)

################## Calculate Net N mineralization
# Net N mineralization is calculated via taking the difference in the sum of
# nitrate and ammonium between the two days,, divided by the length of time
# between the two measurements. For example, to calculate the net M
# mineralization between days 5 and 0, we would use:
# [(Day 5-NO3N_mg.kg + Day 5-NH4N.mg.kg) - (Day 0-NO3N_mg.kg + Day 0-NH4N.mg.kg)] / (5-0)
mineralization_data <- mineralization_data %>%
  group_by(sample) %>%
  mutate(
    net_n_1 = no3n_mg_kg_1 + nh4n_mg_kg_1,
    net_n_2 = no3n_mg_kg_2 + nh4n_mg_kg_2
  ) %>%
  mutate(
    delta = day - lag(day, default = 0),
    net_mineralization_1 = if_else(
      delta == 0,
      0,
      (net_n_1 - lag(net_n_1)) / delta
    ),
    net_mineralization_2 = if_else(
      delta == 0,
      0,
      (net_n_2 - lag(net_n_2)) / delta
    )
  ) %>%
  ungroup()


################## Write out net mineralization data
write.csv(
  mineralization_data,
  here("data/prepped_data", "mineralization_data.csv"),
  row.names = FALSE,
  quote = FALSE
)

################## Join with qPCR data
# We'll join the qPCR data with the mineralization data and save it separately
# since there are fewer qPCR samples.
qpcr_data <- read.csv("data/Incubation_Biomark-qPCR_all_20230328.csv") %>%
  clean_names() %>%
  mutate(across(
    c(crop, addition),
    str_to_lower
  ))

mineralization_and_qpcr_data <- left_join(
  qpcr_data %>%
    select(-c(order, crop, treatment, addition, dna_doe, nh4n_mgl, no3n_mgl)),
  mineralization_data,
  by = c(
    "sample_name" = "initial_name"
  )
)

write.csv(
  mineralization_and_qpcr_data,
  here("data/prepped_data", "mineralization_and_qpcr_data.csv"),
  row.names = FALSE,
  quote = FALSE
)
