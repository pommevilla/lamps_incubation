#!/usr/bin/env Rscript
# ---------------------------
# Calculates priming effects.
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup/setup.R")

## Read in data
mineralization_data <- read.csv("data/prepped_data/mineralization_data.csv")
n2o_data <- read.csv("data/prepped_data/n2o_data.csv")
co2_data <- read.csv("data/prepped_data/co2_data.csv")

co2_data %>% head()

# We'll calculate PE effects based on the end of the experiment,
# and for the first four sampling days.
pe_days <- c(5, 10, 16, 24, 146)
co2_data %>%
    # filter(Day %in% pe_days) %>%
    filter(Day == 146) %>%
    group_by(Crop, Day, Addition) %>%
    summarise(
        mean_n = mean(cum_CO2_flux_ug_g),
        sd = sd(cum_CO2_flux_ug_g),
        se = sd(cum_CO2_flux_ug_g) / sqrt(n())
    ) %>%
    # mutate(
    #     pe = mean_n - lag(mean_n),
    # )
    mutate(across(
        mean_n:se,
        ~ .x - nth(.x, 2),
        .names = "pe_{.col}"
    ))
ggplot(aes(Addition, mean_n, fill = Addition, shape = Addition)) +
    geom_errorbar(
        aes(ymin = mean_n - se, ymax = mean_n + se),
        width = 0.5,
        position = position_dodge(width = 0.1)
    ) +
    geom_point(size = 5) +
    facet_grid(~Crop) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_fill_manual(values = addition_colors)

# Calculates priming effects by comparing emissions to the emissions of the control variable.
# Also calculates the percent difference between the +C/+N treatments and the control.
calculate_pes <- function(this_df, foi, voi) {
    this_var <- deparse(substitute(voi))

    this_df %>%
        group_by(!!sym(foi), Day, Addition) %>%
        summarise(
            mean_n = mean({{ voi }}),
            sd = sd({{ voi }}),
            se = sd({{ voi }}) / sqrt(n())
        ) %>%
        mutate(across(
            mean_n:se,
            ~ .x - nth(.x, 2),
            .names = "pe_{.col}"
        )) %>%
        mutate(
            pe_perc_diff = pe_mean_n / nth(mean_n, 2) * 100
        ) %>%
        rename(
            level = 1
        ) %>%
        mutate(
            factor = foi,
            .before = level
        ) %>%
        mutate(
            var = this_var
        ) %>%
        ungroup() %>%
        mutate(across(
            where(is.numeric),
            ~ round(., 3)
        ))
}

co2_priming_effects <- bind_rows(
    calculate_pes(co2_data, "Treatment", cum_CO2_flux_ug_g),
    calculate_pes(co2_data, "Crop", cum_CO2_flux_ug_g),
)

n2o_priming_effects <- bind_rows(
    calculate_pes(n2o_data, "Treatment", cum_N2O_flux_ug_g),
    calculate_pes(n2o_data, "Crop", cum_N2O_flux_ug_g),
)

mineralization_priming_effects <- bind_rows(
    calculate_pes(mineralization_data, "Treatment", net_min_rate_rel),
    calculate_pes(mineralization_data, "Treatment", net_min_rate_abs),
    calculate_pes(mineralization_data, "Treatment", net_nitr_rate_rel),
    calculate_pes(mineralization_data, "Treatment", net_nitr_rate_abs),
    calculate_pes(mineralization_data, "Crop", net_min_rate_rel),
    calculate_pes(mineralization_data, "Crop", net_min_rate_abs),
    calculate_pes(mineralization_data, "Crop", net_nitr_rate_rel),
    calculate_pes(mineralization_data, "Crop", net_nitr_rate_abs),
)

priming_effects <- bind_rows(
    co2_priming_effects,
    n2o_priming_effects,
    mineralization_priming_effects
)

write.csv(
    priming_effects,
    here::here("results/priming_effects.csv"),
    row.names = FALSE,
    quote = FALSE
)
