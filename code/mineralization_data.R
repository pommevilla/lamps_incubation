#!/usr/bin/env Rscript
# ---------------------------
# Calculates net mineralization rates for the samples
# Outputs a CSV for mineralization calculations and another
# for mineralization joined with qPCR data.
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup/setup.R")

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
  separate(rev_name, c("sample", "timepoint"), sep = -1, remove = FALSE)

################## Calculate Net N mineralization
# Net N mineralization is calculated via taking the difference in the sum of
# nitrate and ammonium between two days, divided by the length of time
# between the two measurements. For example, to calculate the net M
# mineralization between days 32 and 5, we would use:
# [(Day 32-NO3N_mg.kg + Day 32-NH4N.mg.kg) - (Day 5-NO3N_mg.kg + Day 5-NH4N.mg.kg)] / (32-5)
# We calculate two versions - one relative to the previous measurement (as above) called
# net_min_rate_rel, and one relative to the first measurement called net_min_rate_abs.
# We do similar calculations for net nitrification
mineralization_data <- mineralization_data %>%
  group_by(sample) %>%
  mutate(
    net_n_min = no3n_mg_kg_1 + nh4n_mg_kg_1,
  ) %>%
  mutate(
    delta = day - lag(day, default = 0),
    net_min_rate_rel = if_else(
      delta == 0,
      0,
      (net_n_min - lag(net_n_min)) / delta
    ),
    net_min_rate_abs = if_else(
      delta == 0,
      0,
      (net_n_min - first(net_n_min)) / delta
    ),
    net_nitr_rate_rel = if_else(
      delta == 0,
      0,
      (no3n_mg_kg_1 - lag(no3n_mg_kg_1)) / delta
    ),
    net_nitr_rate_abs = if_else(
      delta == 0,
      0,
      (no3n_mg_kg_1 - first(no3n_mg_kg_1)) / delta
    ),
    crop = if_else(
      crop == "corn", "Corn", "Miscanthus"
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
