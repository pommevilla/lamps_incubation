#!/usr/bin/env Rscript
# ---------------------------
# Statistical analysis of priming effect data
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
# ################## Setup
source("code/setup/setup.R")

pe_data <- read.csv("results/priming_effects.csv")
tukey_results <- read.csv("results/stats/n_tukey_results.long.csv")
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

last_days <- c(144, 146)
cumulative_emissions_pe <- pe_data %>%
    filter(Day %in% last_days)

############## End of experiment summaries
summarise_at_day <- function(n_df, voi, sum_day, vars) {
    sum_df <- n_df %>%
        filter(Day == sum_day) %>%
        group_by(!!sym(voi)) %>%
        summarise(across(
            any_of(vars),
            list(
                mean = mean,
                sd = sd
            ),
            .names = "{.col}_{.fn}"
        )) %>%
        pivot_longer(-1) %>%
        rename(level = 1) %>%
        mutate(factor = voi, .before = 1)

    return(sum_df)
}

summarise_variables_by_treatment <- function(voi) {
    co2_means <- summarise_at_day(co2_data, {{ voi }}, 146, co2_vars)
    n2o_means <- summarise_at_day(n2o_data, {{ voi }}, 144, n2o_vars)
    inorganic_n_means <- summarise_at_day(nh4_no3_min_data, {{ voi }}, 144, mineralization_variables)

    bind_rows(
        co2_means, n2o_means, inorganic_n_means
    )
}

eoe_summaries <- bind_rows(
    summarise_variables_by_treatment("Crop"),
    summarise_variables_by_treatment("Addition"),
    summarise_variables_by_treatment("Treatment")
)

write.csv(
    eoe_summaries,
    here::here("results/eoe_summaries/all_eoe_summaries.csv"),
    row.names = FALSE,
    quote = FALSE
)

summarise_addition_at_day <- function(n_df, voi, sum_day, vars) {
    sum_df <- n_df %>%
        filter(Day == sum_day) %>%
        group_by(Addition, !!sym(voi)) %>%
        summarise(across(
            any_of(vars),
            list(
                mean = mean,
                sd = sd
            ),
            .names = "{.col}_{.fn}"
        )) %>%
        ungroup() %>%
        unite("level", Addition, !!sym(voi), sep = " ") %>%
        pivot_longer(-1) %>%
        mutate(factor = voi, .before = 1)

    return(sum_df)
}

summarise_addition_in_treatment_at_day <- function(voi) {
    co2_means <- summarise_addition_at_day(co2_data, {{ voi }}, 146, co2_vars)
    n2o_means <- summarise_addition_at_day(n2o_data, {{ voi }}, 144, n2o_vars)
    inorganic_n_means <- summarise_addition_at_day(nh4_no3_min_data, {{ voi }}, 144, mineralization_variables)

    bind_rows(
        co2_means, n2o_means, inorganic_n_means
    )
}

addition_eoe_summaries <- bind_rows(
    summarise_addition_in_treatment_at_day("Crop"),
    summarise_addition_in_treatment_at_day("Treatment")
)

write.csv(
    addition_eoe_summaries,
    here::here("results/eoe_summaries/addition_eoe_summaries.csv"),
    row.names = FALSE,
    quote = FALSE
)


############## Summarise addition effects by crop in N rate
summarise_crop_addition_at_day <- function(n_df, sum_day, vars) {
    sum_df <- n_df %>%
        filter(Day == sum_day) %>%
        group_by(Crop, Addition, Treatment) %>%
        summarise(across(
            any_of(vars),
            list(
                mean = mean,
                sd = sd
            ),
            .names = "{.col}_{.fn}"
        )) %>%
        ungroup() %>%
        unite("level", Addition, Treatment, sep = " ") %>%
        pivot_longer(-(1:2))
    # %>%
    # mutate(Crop = plant_type, .before = 1)

    return(sum_df)
}

summarise_crop_addition_in_treatment_at_day <- function() {
    co2_means <- summarise_crop_addition_at_day(co2_data, 146, co2_vars)
    n2o_means <- summarise_crop_addition_at_day(n2o_data, 144, n2o_vars)
    inorganic_n_means <- summarise_crop_addition_at_day(nh4_no3_min_data, 144, mineralization_variables)

    bind_rows(
        co2_means, n2o_means, inorganic_n_means
    )
}

crop_addition_treatment_eoe_summaries <- summarise_crop_addition_in_treatment_at_day()

write.csv(
    crop_addition_treatment_eoe_summaries,
    here::here("results/eoe_summaries/crop_addition_treatment_eoe_summaries.csv"),
    row.names = FALSE,
    quote = FALSE
)


############## Tukey HSD results exploring
tukey_results %>%
    select(term, contrast, estimate, adj.p.value, n_var) %>%
    filter(term == "Treatment:Day") %>%
    filter(n_var == "cum_CO2_flux_ug_g") %>%
    # filter(n_var == "cum_N2O_flux_ug_g") %>%
    # filter(contrast == "Nitrogen:146-Carbon:146")
    filter(str_detect(contrast, "146"))

############## Priming Effects
cumulative_emissions_pe %>%
    filter(var == "cum_N2O_flux_ug_g") %>%
    filter(factor == "Crop") %>%
    select(level, Addition, mean_n, pe_mean_n, pe_perc_diff)


tukey_results %>%
    select(term, contrast, estimate, adj.p.value, n_var) %>%
    filter(term == "Addition:Treatment:Day") %>%
    separate(contrast, into = c("contrast_1", "contrast_2"), sep = "-") %>%
    filter(str_detect(contrast_1, "146") & str_detect(contrast_2, "146")) %>%
    filter(str_detect(contrast_1, "0N") & str_detect(contrast_2, "0N"))


tukey_results %>%
    filter(term == "Crop:Addition:Treatment:Day") %>%
    separate(contrast, into = c("contrast_1", "contrast_2"), sep = "-") %>%
    filter(n_var == "cum_nh4") %>%
    # filter(str_detect(contrast_1, "146") & str_detect(contrast_2, "146"))
    filter(str_detect(contrast_1, "144") & str_detect(contrast_2, "144")) %>%
    filter(str_detect(contrast_1, "0N") & str_detect(contrast_2, "0N")) %>%
    # filter(str_detect(contrast_1, "Corn") & str_detect(contrast_2, "Corn"))
    filter(str_detect(contrast_1, "Miscanthus") & str_detect(contrast_2, "Miscanthus"))
########################
n2o_data %>%
    filter(Day == 144) %>%
    group_by(Crop, Treatment, Addition) %>%
    summarise(across(
        any_of(contains("N2O")),
        list(mean = mean)
    ))


######################## Repeat for mineralization rates at days 4 and 15

all_inorganic_n_summaries <- nh4_no3_min_data %>%
    group_by(Crop, Addition, Treatment, Day) %>%
    summarise(across(
        any_of(mineralization_variables),
        list(
            mean = mean,
            sd = sd
        ),
        .names = "{.col}_{.fn}"
    )) %>%
    ungroup() %>%
    pivot_longer(-(1:4))

write.csv(
    all_inorganic_n_summaries,
    here::here("results/eoe_summaries/all_inorganic_pe_summaries.csv"),
    row.names = FALSE,
    quote = FALSE
)
