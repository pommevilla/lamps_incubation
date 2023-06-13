#!/usr/bin/env Rscript
# ---------------------------
# Train linear models and calculate correlations for qPCR data
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
# ################## Setup
source("code/setup/setup.R")

######## Read in data
qpcr_min_data <- read.csv("data/prepped_data/mineralization_and_qpcr_data.csv") %>%
    mutate(
        Day = as.factor(Day)
    )

# Variables for the for loops
n_variables <- c(
    mineralization_variables,
    n2o_vars
)

grouping_variables <- c("Treatment", "Crop", "Addition")

# Empty data frame to store the results
result_df <- data.frame(
    treatment = character(),
    level = character(),
    qpcr_var = character(),
    n_var = character(),
    statistic = character(),
    estimate = numeric(),
    p_value = numeric()
)

# Outer loop groups the data frame
for (group_var in grouping_variables) {
    grouped_data <- qpcr_min_data %>%
        group_by(!!sym(group_var))

    print(paste("Grouping by", group_var))

    # These inner loops will run linear regressions and Pearson/Spearman
    # correlations of the qpcr_var against the n_var for each factor
    # level of the grouping variable.
    for (qpcr_var in qpcr_variables) {
        print(paste("Running", qpcr_var))
        for (n_var in n_variables) {
            # This will run a linear regression of qpcr_var ~ n_var for each
            # factor level and record the r^2 and p-value in the result df.
            # print("Running linear regression...")
            model_formula <- as.formula(paste(qpcr_var, "~", n_var))
            lm_result <- grouped_data %>%
                do(broom::tidy(lm(model_formula, data = .)))


            for (factor_level in unique(grouped_data[[group_var]])) {
                this_lm_result <- lm_result %>%
                    filter(!!sym(group_var) == factor_level) %>%
                    filter(term == n_var)

                result_df <- result_df %>%
                    add_row(
                        treatment = group_var,
                        level = factor_level,
                        qpcr_var = qpcr_var,
                        n_var = n_var,
                        statistic = "lm_r_squared",
                        estimate = this_lm_result$estimate,
                        p_value = this_lm_result$p.value,
                    )
            }

            # print("Running correlations...")
            for (cor_type in c("pearson", "spearman")) {
                cors <- grouped_data %>%
                    summarise(out = as_tibble(
                        cor.test(
                            !!sym(qpcr_var),
                            !!sym(n_var),
                            method = cor_type
                        )[c("estimate", "p.value")]
                    )) %>%
                    unnest(out)

                for (factor_level in unique(grouped_data[[group_var]])) {
                    this_result <- cors %>%
                        filter(!!sym(group_var) == factor_level)

                    result_df <- result_df %>%
                        add_row(
                            treatment = group_var,
                            level = factor_level,
                            qpcr_var = qpcr_var,
                            n_var = n_var,
                            statistic = cor_type,
                            estimate = this_result$estimate,
                            p_value = this_result$p.value,
                        )
                }
            }
        }
    }
}

# Sanity check:
#    20 qpcr variables
#  *  10 n variables
#  *  (2 + 3 + 3) grouping variables
#  *  3 statistics
#  = 4800 rows
result_df <- result_df %>%
    mutate(across(
        estimate:p_value,
        ~ round(., 3)
    ))

write.csv(
    result_df,
    "results/stats/qpcr_n_correlations_lms_all.csv",
    row.names = FALSE,
    quote = FALSE
)

write.csv(
    result_df %>% filter(!str_detect(qpcr_var, "log")),
    "results/stats/qpcr_n_correlations_lms_ave.csv",
    row.names = FALSE,
    quote = FALSE
)

write.csv(
    result_df %>% filter(str_detect(qpcr_var, "log")),
    "results/stats/qpcr_n_correlations_lms_log.csv",
    row.names = FALSE,
    quote = FALSE
)
