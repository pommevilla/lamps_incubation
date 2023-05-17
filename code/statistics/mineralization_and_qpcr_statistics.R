#!/usr/bin/env Rscript
# ---------------------------
# ANOVA analysis of mineralization against main factors.
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------

################## Setup
# General setup scripts
source("code/setup.R")

# Read in data
mineralization_and_amoa_data <- read.csv("data/prepped_data/mineralization_and_qpcr_data.csv")

# Function to plot tidied ANOVA results
plot_tidied_anova_results <- function(tidied_anova_results) {
  tidied_anova_results %>%
    ggplot(aes(flux_var, term, fill = sig)) +
    geom_tile(color = "white", size = 1) +
    scale_fill_viridis_d(option = "magma", direction = -1) +
    theme(
      panel.border = element_blank(),
      axis.line = element_blank(),
      strip.text = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_markdown()
    ) +
    labs(
      x = "",
      y = "",
      fill = "Significance"
    )
}

################ ANOVA
mineralization_variables <- c(
  "net_min_rate_rel", "net_min_rate_abs",
  "net_nitr_rate_rel", "net_nitr_rate_abs",
  "no3n_mg_kg_1", "nh4n_mg_kg_1"
)

qpcr_variables <- c(
  "log_012", "ave_012", "log_025", "ave_025",
  "log_039", "ave_039", "f1r2_log", "f1r2_ave"
)

# Helper function to run a formula string against the mineralization data
run_anovas <- function(lhs_vars, rhs_formula) {
  formulae <- lapply(
    lhs_vars,
    function(x) as.formula(paste(x, rhs_formula))
  )

  anova_res <- lapply(formulae, function(x) broom::tidy(aov(x, data = mineralization_and_amoa_data)))
  names(anova_res) <- format(formulae)
  names(anova_res) <- str_extract(names(anova_res), "^[^\\s~]+")

  anova_res <- lapply(
    seq_along(anova_res),
    function(i) anova_res[[i]] %>% mutate(flux_var = names(anova_res)[[i]])
  ) %>%
    bind_rows() %>%
    filter(term != "Residuals") %>%
    mutate(sig = case_when(
      p.value < 0.05 & p.value > 0.01 ~ "*",
      p.value < 0.01 & p.value > 0.001 ~ "**",
      p.value < 0.001 ~ "***",
      TRUE ~ "NS"
    )) %>%
    mutate(
      flux_var = str_replace(flux_var, "f1r2_ave", "F1R2"),
      flux_var = str_replace(flux_var, "ave_012", "AmoA 12"),
      flux_var = str_replace(flux_var, "ave_025", "AmoA 25"),
      flux_var = str_replace(flux_var, "ave_039", "AmoA 39"),
      term = str_replace(term, "day", "Day"),
      term = str_replace(term, "crop", "Crop"),
      term = str_replace(term, "addition", "Addition"),
      term = str_replace(term, "treatment", "Treatment")
    ) %>%
    mutate(order = str_count(term, ":")) %>%
    mutate(term = fct_reorder(term, order))

  return(anova_res)
}

################## ANOVA of mineralization variables against main factors
factors_min_anovas <- run_anovas(
  lhs_vars = mineralization_variables,
  rhs_formula = "~ crop * addition * treatment * day"
)


# Write out results
write.csv(
  factors_min_anovas,
  here::here("results/stats", "factors_mineralization_anovas.csv"),
  quote = FALSE,
  row.names = FALSE
)

# Writing out results in a wide format
factors_min_anovas %>%
  select(term, p.value, flux_var) %>%
  pivot_wider(names_from = "flux_var", values_from = "p.value") %>%
  write.csv(
    here::here("results/stats", "factors_mineralization_anovas_wide.csv"),
    quote = FALSE,
    row.names = FALSE
  )

# Plotting
factors_min_anovas <- factors_min_anovas %>%
  mutate(
    flux_var = case_when(
      flux_var == "net_min_rate_rel" ~ "Net rel.<br>min. rate",
      flux_var == "net_min_rate_abs" ~ "Net abs.<br>min. rate",
      flux_var == "net_nitr_rate_rel" ~ "Net rel.<br>nitr. rate",
      flux_var == "net_nitr_rate_abs" ~ "Net abs.<br>nitr. rate",
      flux_var == "no3n_mg_kg_1" ~ paste0(nitrate_label, "<br>", flux_units),
      flux_var == "nh4n_mg_kg_1" ~ paste0(ammonia_label, "<br>", flux_units)
    )
  )

# Abridged version for presentations
factors_min_anovas %>%
  filter(!flux_var %in% c("Net abs.<br>min. rate", "Net rel.<br>nitr. rate", "Net abs.<br>nitr. rate")) %>%
  plot_tidied_anova_results()

ggsave(
  here::here("figures/stats/factors_mineralization_anovas_abridged.png"),
)

plot_tidied_anova_results(factors_min_anovas)

ggsave(
  here::here("figures/stats/factors_mineralization_anovas.png"),
  width = 8,
  height = 7,
  unit = "in"
)

# Tukey HSD
TukeyHSD(factors_min_anovas)

tidied_factors_day_anova_results <- broom::tidy(factors_day_anova_results) %>%
  mutate(sig = p.value < 0.05) %>%
  mutate(p.value = round(p.value, 2)) %>%
  filter(term != "Residuals") %>%
  mutate(term = str_replace(term, "crop", "Crop")) %>%
  mutate(term = str_replace(term, "addition", "Addition")) %>%
  mutate(term = str_replace(term, "treatment", "Treatment")) %>%
  mutate(term = str_replace(term, "dna_doe", "DOE")) %>%
  mutate(order = str_count(term, ":"))

write.csv(
  tidied_qpcr_day_anova_results,
  here::here("results", "mineralization_factors_anova_results.csv"),
  quote = FALSE,
  row.names = FALSE
)


################## ANOVA of qPCR variables against main factors
factors_qpcr_anovas <- run_anovas(
  lhs_vars = qpcr_variables,
  rhs_formula = "~ crop * addition * treatment * day"
)

# Writing out results
write.csv(
  factors_qpcr_anovas,
  here::here("results/stats", "factors_qpcr_anovas.csv"),
  quote = FALSE,
  row.names = FALSE
)

# Writing out results in a wide format
factors_qpcr_anovas %>%
  select(term, p.value, flux_var) %>%
  pivot_wider(names_from = "flux_var", values_from = "p.value") %>%
  write.csv(
    here::here("results/stats", "factors_qpcr_anovas_wide.csv"),
    quote = FALSE,
    row.names = FALSE
  )

# Plotting the qPCR anova results without the log GCN
plot_tidied_anova_results(
  factors_qpcr_anovas %>%
    filter(!str_detect(flux_var, "log"))
)

ggsave(
  here::here("figures/stats/factors_qpcr_anovas.png"),
  width = 8,
  height = 7,
  unit = "in"
)

################ Relating AmoAs to mineralization variables

log_mlr_model <- lm(
  net_mineralization_1 ~ f1r2_log + log_012 + log_025 + log_039,
  data = mineralization_and_amoa_data
)

summary(log_mlr_model)
tidy_log_mlr_model <- broom::tidy(log_mlr_model) %>%
  select(term, p.value) %>%
  filter(term != "(Intercept)")

ave_mlr_model <- lm(
  net_mineralization_1 ~ f1r2_ave + ave_012 + ave_025 + ave_039,
  data = mineralization_and_amoa_data
)

summary(ave_mlr_model)
tidy_ave_mlr_model <- broom::tidy(ave_mlr_model) %>%
  select(term, p.value) %>%
  filter(term != "(Intercept)") %>%
  mutate()

multiple_linear_regression_tests <- bind_rows(
  tidy_ave_mlr_model,
  tidy_log_mlr_model
) %>%
  mutate(p.value = round(p.value, 3)) %>%
  mutate(significant = p.value < 0.05)

write.csv(
  multiple_linear_regression_tests,
  here::here("results", "mineralization_multiple_linear_regression_tests.csv"),
  quote = FALSE,
  row.names = FALSE
)

# Linear regressions of the qPCR against nitrate, ammonium, and net m
# We want to see how each individual
combined_equations <- c(
  " ~ ave_012 * ave_025 * ave_039 * f1r2_ave * day",
  " ~ log_012 * log_025 * log_039 * f1r2_log * day"
)

results <- data.frame(
  dep_var = character(),
  indep_var = character(),
  r_squared = numeric(),
  p_value = numeric()
)

# Iterate over each dependent variable and independent variable combination
for (dep_var in mineralization_variables) {
  # Training LMs on one variable at a time
  for (ind_var in qpcr_variables) {
    model <- lm(mineralization_and_amoa_data[[dep_var]] ~ mineralization_and_amoa_data[[ind_var]])

    r_squared <- summary(model)$r.squared
    p_value <- summary(model)$coefficients[2, 4]

    results <- rbind(results, data.frame(
      dependent_var = dep_var,
      independent_var = ind_var,
      r_squared = r_squared,
      p_value = p_value
    ))
  }

  # Training multiple linear regression on all variables

  # Do it again for all variables
  for (eq in combined_equations) {
    this_formula <- paste0(dep_var, eq)
    model <- lm(this_formula, data = mineralization_and_amoa_data)

    r_squared <- summary(model)$adj.r.squared
    p_value <- summary(model)$coefficients[2, 4]

    results <- rbind(results, data.frame(
      dependent_var = dep_var,
      independent_var = eq,
      R_squared = r_squared,
      p_value = p_value
    ))
  }
}



# View results
results %>%
  mutate(significant = if_else(
    p_value < 0.05, "*", ""
  )) %>%
  mutate(
    r_squared = round(r_squared, 3),
    r_squared = paste0(r_squared, significant)
  ) %>%
  select(dependent_var, independent_var, r_squared) %>%
  pivot_wider(names_from = dependent_var, values_from = r_squared)


######## Correlations between qPCR and mineralization variables
cor_results <- data.frame(
  variable1 = character(),
  variable2 = character(),
  correlation_coefficient = numeric(),
  p_value = numeric()
)

# Perform cor.test for each combination
for (qpcr_var in qpcr_variables) {
  for (min_var in mineralization_variables) {
    cor_test <- cor.test(
      mineralization_and_amoa_data[[qpcr_var]],
      mineralization_and_amoa_data[[min_var]],
      method = "spearman"
    )

    this_result <- data.frame(
      variable1 = qpcr_var,
      variable2 = min_var,
      correlation_coefficient = cor_test$estimate,
      p_value = cor_test$p.value
    )

    cor_results <- rbind(cor_results, this_result)
  }
}

cor_results <- cor_results %>%
  mutate(
    correlation_coefficient = round(correlation_coefficient, 3),
    p_value = round(p_value, 3)
  ) %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )) %>%
  mutate(
    cor_sig = paste0(correlation_coefficient, significance)
  )

write.csv(
  cor_results,
  file = "results/stats/amoa_qpcr_correlations.csv",
  row.names = FALSE,
  quote = FALSE
)

cor_results %>%
  filter(!str_detect(variable1, "log")) %>%
  select(variable1, variable2, cor_sig) %>%
  pivot_wider(id_cols = variable2, names_from = variable1, values_from = cor_sig) %>%
  view()
mutate(
  variable2 = case_when(
    variable2 == "net_min_rate_rel" ~ "Net rel.<br>min. rate",
    variable2 == "net_min_rate_abs" ~ "Net abs.<br>min. rate",
    variable2 == "net_nitr_rate_rel" ~ "Net rel.<br>nitr. rate",
    variable2 == "net_nitr_rate_abs" ~ "Net abs.<br>nitr. rate",
    variable2 == "no3n_mg_kg_1" ~ paste0(nitrate_label, "<br>(mg/kg)"),
    variable2 == "nh4n_mg_kg_1" ~ paste0(ammonia_label, "<br>(mg/kg)")
  ),
  variable1 = case_when(
    variable1 == "f1r2_ave" ~ "F1R2",
    TRUE ~ str_replace(variable1, "ave_", "amoA ")
  )
) %>%
  ggplot(aes(variable1, variable2, fill = correlation_coefficient)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = round(correlation_coefficient, 2)), size = 4) +
  theme(
    text = element_text(size = 16),
    axis.line = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_markdown()
  ) +
  # scale_fill_viridis_c(option = "magma") +
  scale_fill_gradientn(colors = met.brewer("Hiroshige", direction = -1)) +
  coord_fixed() +
  labs(
    x = "qPCR",
    y = "",
    fill = "Pearson"
  )


cor_results %>%
  select(variable1, variable2, correlation_coefficient) %>%
  pivot_wider(
    id_cols = variable1,
    names_from = variable2,
    values_from = correlation_coefficient
  )
