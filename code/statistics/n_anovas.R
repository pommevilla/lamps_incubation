#!/usr/bin/env Rscript
# ---------------------------
# ANOVAS for CO2, N2O, NH4, NO3, mineralization, and nitrification
# against main factors.
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
# ################## Setup
source("code/setup/setup.R")

######## Read in data
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

##### Helper functions
# Runs ANOVAS of Crop, Addition, Treatment, and Day against the LHS variables,
# which are a list of the variables in the dataframe to test against.
# n_df is the dataframe used in the anova
run_anovas <- function(lhs_vars, n_df, rhs_formula = "~ Crop * Addition * Treatment * Day") {
    formulae <- lapply(
        lhs_vars,
        function(x) as.formula(paste(x, rhs_formula))
    )

    anova_res <- lapply(formulae, function(x) broom::tidy(aov(x, data = n_df)))
    names(anova_res) <- format(formulae)
    names(anova_res) <- str_extract(names(anova_res), "^[^\\s~]+")

    anova_res <- lapply(
        seq_along(anova_res),
        function(i) anova_res[[i]] %>% mutate(n_var = names(anova_res)[[i]])
    ) %>%
        bind_rows() %>%
        filter(term != "Residuals") %>%
        mutate(
            term = str_replace(term, "day", "Day"),
            term = str_replace(term, "crop", "Crop"),
            term = str_replace(term, "addition", "Addition"),
            term = str_replace(term, "treatment", "Treatment")
        ) %>%
        mutate(sig = get_p_sig(p.value)) %>%
        mutate(p.value = round(p.value, 3)) %>%
        mutate(order = str_count(term, ":")) %>%
        mutate(term = fct_reorder(term, order))

    return(anova_res)
}

################## Run ANOVAS
# We'll run ANOVAs for each N2O, CO2, Nh4, NO3, mineralization, and nitrification.
n2o_anovas <- run_anovas(c("N2ON_flux_ug_g_d", "cum_N2O_flux_ug_g"), n2o_data)
co2_anovas <- run_anovas(c("CO2_flux_ug_g_d", "cum_CO2_flux_ug_g"), co2_data)
nh4_no3_min_anovas <- run_anovas(
    c(
        "net_min_rate_rel", "net_min_rate_abs",
        "net_nitr_rate_rel", "net_nitr_rate_abs",
        "no3n_mg_kg_1", "cum_no3",
        "nh4n_mg_kg_1", "cum_nh4"
    ),
    nh4_no3_min_data
)

# Combining all the results together into long format
anova_results_long <- bind_rows(
    n2o_anovas, co2_anovas, nh4_no3_min_anovas
)

write.csv(
    anova_results_long,
    here::here("results/stats", "n_anova_results.long.csv"),
    row.names = FALSE,
    quote = FALSE
)

# Pivoting to a wide format for presentations
anova_results_wide <- anova_results_long %>%
    select(n_var, term, p.value) %>%
    pivot_wider(names_from = n_var, values_from = p.value)

write.csv(
    anova_results_wide,
    here::here("results/stats", "n_anova_results.wide.csv"),
    row.names = FALSE,
    quote = FALSE
)

########## Tukey HSD
# Helper function to run Tukey HSD tests on N data
get_tukey_results <- function(n_vars, n_df) {
    results <- data.frame(
        term = character(),
        contrast = character(),
        null.value = numeric(),
        estimate = numeric(),
        conf.low = numeric(),
        conf.high = numeric(),
        adj.p.value = numeric(),
        n_var = character()
    )

    for (n_var in n_vars) {
        anova_model <- aov(
            as.formula(paste(n_var, "~ Crop * Addition * Treatment * Day")),
            data = n_df
        )

        these_tukey_results <- TukeyHSD(anova_model) %>%
            broom::tidy() %>%
            mutate(n_var = n_var)

        results <- bind_rows(results, these_tukey_results)
    }

    return(results)
}

n2o_vars <- c("N2ON_flux_ug_g_d", "cum_N2O_flux_ug_g")
tester <- get_tukey_results(n2o_vars, n2o_data)

co2_vars <- c("CO2_flux_ug_g_d", "cum_CO2_flux_ug_g")
tester2 <- get_tukey_results(co2_vars, co2_data)

inorganic_n_vars <- c(
    "no3n_mg_kg_1", "cum_no3", "nh4n_mg_kg_1", "cum_nh4",
    "net_min_rate_rel", "net_min_rate_abs",
    "net_nitr_rate_rel", "net_nitr_rate_abs"
)

tester3 <- get_tukey_results(inorganic_n_vars, nh4_no3_min_data)

tukey_results_long <- bind_rows(tester, tester2, tester3)

write.csv(
    tukey_results_long,
    here::here("results/stats/n_tukey_results.long.csv"),
    row.names = FALSE,
    quote = FALSE
)
