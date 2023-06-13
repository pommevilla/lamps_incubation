source("code/setup/setup.R")

lms_rs <- read.csv("results/stats/qpcr_n_correlations_lms.csv")

lms_rs %>%
    filter(treatment == "Crop") %>%
    filter(n_var == "N2ON_flux_ug_g_d") %>%
    filter(p_value < 0.05)

aov(cum_CO2_flux_ug_g ~ Treatment * Crop * Day * Addition, data = co2_data) %>%
    summary()

aov(cum_N2O_flux_ug_g ~ Treatment * Crop * Day * Addition, data = n2o_data) %>%
    summary()

tukey_results <- tukey_results %>%
    separate(contrast, into = c("contrast_1", "contrast_2"), sep = "-")

tukey_results %>%
    # filter(term == "Crop:Addition:Treatment:Day") %>%
    filter(term == "Crop:Day") %>%
    filter(n_var == "cum_N2O_flux_ug_g") %>%
    # filter(n_var == "cum_no3") %>%
    # filter(n_var == "cum_CO2_flux_ug_g") %>%
    # filter(adj.p.value < 0.05) %>%
    filter(str_detect(contrast_1, "144") & str_detect(contrast_2, "144"))
# filter(str_detect(contrast_1, "146") & str_detect(contrast_2, "146")) %>%
# filter(str_detect(contrast_1, "0N") & str_detect(contrast_2, "0N"))  %>%
# filter(str_detect(contrast_1, "Control") & str_detect(contrast_2, "Control")) %>%
# filter(str_detect(contrast_1, "Miscanthus") & str_detect(contrast_2, "Miscanthus"))



n2o_data %>%
    filter(Day == 144) %>%
    group_by(Treatment, Crop) %>%
    summarise(
        mean_n = mean(cum_N2O_flux_ug_g),
    ) %>%
    ungroup() %>%
    mutate(
        perc_diff = calc_percentage(last(mean_n), mean_n)
    )

co2_data %>%
    filter(Day == 146) %>%
    group_by(Addition) %>%
    summarise(
        mean_n = mean(cum_CO2_flux_ug_g),
    ) %>%
    ungroup() %>%
    mutate(
        perc_diff = calc_percentage(mean_n, nth(mean_n, 2))
    )
