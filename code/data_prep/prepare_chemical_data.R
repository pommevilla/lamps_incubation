#!/usr/bin/env Rscript
# ---------------------------
# Prepares all the chemical data for plotting/analysis
# Also calculates the net abs/rel mineralization and nitrification rates.
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup/setup.R")

################## N2O data
n2o_data <- read.csv("data/chemical/N2O.csv") %>%
    mutate(
        addition = case_when(
            addition == "Cntrl" ~ "Control",
            addition == "+C" ~ "Carbon",
            TRUE ~ "Nitrogen"
        )
    ) %>%
    mutate(
        crop = if_else(
            crop == "Mxg", "Miscanthus", crop
        )
    ) %>%
    rename(
        Crop = crop,
        Treatment = fert,
        Addition = addition,
        Day = doe
    )

write.csv(
    n2o_data,
    here::here("data/prepped_data", "n2o_data.csv"),
    row.names = FALSE,
    quote = FALSE
)

################## CO2 data
tgas <- read.csv("data/chemical/tgas1.csv") %>%
    select(
        sample,
        Day = doe_T1,
        sample_id = sampling_id_T1,
        Crop = crop_T1,
        Treatment = fert_T1,
        Addition = addition_T1,
        CO2_flux_ug_g_d,
        cum_CO2_flux_ug_g,
    ) %>%
    mutate(
        Crop = if_else(
            Crop == "Mxg", "Miscanthus", Crop
        ),
        Addition = case_when(
            Addition == "Cntrl" ~ "Control",
            Addition == "+C" ~ "Carbon",
            TRUE ~ "Nitrogen"
        )
    )

write.csv(
    tgas,
    here::here("data/prepped_data", "co2_data.csv"),
    row.names = FALSE,
    quote = FALSE
)
