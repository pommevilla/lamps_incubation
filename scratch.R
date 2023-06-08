source("code/setup/setup.R")

lms_rs <- read.csv("results/stats/qpcr_n_correlations_lms.csv")

lms_rs %>%
    filter(treatment == "Crop") %>%
    filter(n_var == "N2ON_flux_ug_g_d") %>%
    filter(p_value < 0.05)
