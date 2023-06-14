#!/usr/bin/env Rscript
# ---------------------------
# ANOVAS for CO2, N2O, NH4, NO3, mineralization, and nitrification
# against main factors.
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
# ################## Setup
source("code/setup/setup.R")

######## Read in data
qpcr_data <- read.csv("data/prepped_data/qpcr_data.csv") %>%
    mutate(
        Day = as.factor(Day),
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
        mutate(sig = get_p_sig(p.value)) %>%
        mutate(across(
            sumsq:p.value,
            ~ round(., 3)
        )) %>%
        mutate(order = str_count(term, ":")) %>%
        mutate(term = fct_reorder(term, order))

    return(anova_res)
}

################## Run ANOVAS
anova_results_long <- run_anovas(qpcr_variables, qpcr_data)

write.csv(
    anova_results_long,
    here::here("results/stats/qpcr_anova_results_long.csv"),
    row.names = FALSE,
    quote = FALSE
)

# Pivoting to a wide format for presentations
anova_results_wide <- anova_results_long %>%
    select(n_var, term, p.value) %>%
    pivot_wider(names_from = n_var, values_from = p.value)

write.csv(
    anova_results_wide,
    here::here("results/stats", "qpcr_anova_results_wide.csv"),
    row.names = FALSE,
    quote = FALSE
)

# Visually inspecting for normality
# Untransformed values appear to mostly be normally distributed,
qpcr_data %>%
    pivot_longer(
        any_of(qpcr_variables)
    ) %>%
    mutate(
        nor = if_else(
            str_detect(name, "log"),
            "logged",
            "non_logged",
        )
    ) %>%
    ggplot(aes(value)) +
    geom_histogram() +
    facet_grid(nor ~ name, scales = "free")

#################### Tukey results
get_tukey_results <- function(qpcr_vars, n_df) {
    results <- data.frame(
        term = character(),
        contrast = character(),
        null.value = numeric(),
        estimate = numeric(),
        conf.low = numeric(),
        conf.high = numeric(),
        adj.p.value = numeric(),
        gene_var = character()
    )

    for (qpcr_var in qpcr_vars) {
        anova_model <- aov(
            as.formula(paste(qpcr_var, "~ Crop * Addition * Treatment * Day")),
            data = n_df
        )

        these_tukey_results <- TukeyHSD(anova_model) %>%
            broom::tidy() %>%
            mutate(qpcr_var = qpcr_var)

        results <- bind_rows(results, these_tukey_results)
    }

    results <- results %>%
        mutate(adj.p.value = round(adj.p.value, 3)) %>%
        select(term, contrast, estimate, adj.p.value, qpcr_var)

    return(results)
}

qpcr_tukey <- get_tukey_results(qpcr_variables, qpcr_data) %>%
    separate(contrast, into = c("contrast_1", "contrast_2"), sep = "-")

write.csv(
    qpcr_tukey,
    here::here("results/stats/qpcr_tukey_results.csv"),
    row.names = FALSE,
    quote = FALSE
)

################## Exploring Tukey results
qpcr_tukey %>%
    head()

aov(norB.001 ~ Crop * Addition * Treatment * Day, data = qpcr_data) %>%
    summary()
