#!/usr/bin/env Rscript
################## Setup
library(tidyverse)
library(readxl)
library(ggtext)

source("code/setup.R")

################## Read in data
# The different calculations for ammonium and nitrate are contained in the first sheet,
# while the day of extraction is contained in the Sample Info sheet.
# We read in both and join them together so we can calculate stuff.
mineralization_data <- read_xlsx("data/FILE_4247_sjh.xlsx", sheet = "Sheet1")
# dna_doe_info <- read_xlsx("data/FILE_4247_sjh.xlsx", sheet = "Sample Info") %>% 
#   select(`Sample Name`, `Sample Order`, `Date Order`, DNA_DOE)

mineralization_data <- mineralization_data %>% 
  # left_join(dna_doe_info, by = c("Initial_name" = "Sample Name")) %>% 
  clean_names() %>% 
  mutate(
    crop = str_to_lower(crop),
    treatment = as.factor(case_when(
      treatment == "100N" ~ "112N",
      treatment == "300N" ~ "336N",
      TRUE ~ treatment
    ))
  ) %>% 
  separate(initial_name, c("sample", "timepoint"), sep = -1, remove = FALSE)

################## Calculate Net N mineralization
# Net N mineralization is calculated via taking the difference in the sum of
# nitrate and ammonium between the two days,, divided by the length of time between the
# two measurements. For example, to calculate the net M mineralizaiton between 
# days 5 and 0, we would use:
# [(Day 5-NO3N_mg.kg + Day 5-NH4N.mg.kg) - (Day 0-NO3N_mg.kg + Day 0-NH4N.mg.kg)] / (5-0)		
# We have two different calculations for ammonium and nitrate, so we'll calculate both

mineralization_data <- mineralization_data %>%
  group_by(sample) %>% 
  mutate(
    net_n_1 = no3n_mg_kg_1 + nh4n_mg_kg_1,
    net_n_2 = no3n_mg_kg_2 + nh4n_mg_kg_2
  ) %>% 
  mutate(
    # net_mineralization_1 = (net_n_1 - lag(net_n_1, default = net_n_1[1]) / (dna_doe - lag(dna_doe)))
    delta = day - lag(day, default = 0),
    net_mineralization_1 = if_else(
      delta == 0,
      0,
      (net_n_1 - lag(net_n_1)) / delta
    ), 
    net_mineralization_2 = if_else(
      delta == 0,
      0,
      (net_n_2 - lag(net_n_2)) / delta
    )
  ) %>% 
  ungroup()
  

################## Plots

# This function will take in a grouping variable, a palette, and a y variable
# and plot that y variable against time. 
plot_mineralization <- function(voi, palette, net_min_var, label) {
  this_dodge <- position_dodge(width = 0.1)

  
  data_df <- mineralization_data %>%
    group_by(dna_doe, {{ voi }}) %>%
    summarize(
      mean_n = mean({{ net_min_var }}),
      sd = sd({{ net_min_var }}),
      se = sd({{ net_min_var}}) / sqrt(n())
    ) %>% 
    ungroup()
  
  p <- data_df %>%
    ggplot(aes(dna_doe, mean_n, fill = {{ voi }}, group = {{ voi }}, shape = {{ voi }})) +
    # geom_hline(yintercept = 23, color = "gray", linetype = "dashed") +
    geom_errorbar(
      aes(ymin = mean_n - se, ymax = mean_n + se),
      width = 1,
      position = this_dodge
    ) +
    geom_line() +
    geom_point(size = 4, position = this_dodge) +
    scale_fill_manual(values = palette) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_x_continuous(breaks = unique(mineralization_data$dna_doe)) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      plot.title = element_text(hjust = 0.5),
      aspect.ratio = 1
    ) +
    labs(
      y = "Net minieralization",
      x = "Day of experiment",
      fill = label,
      shape = label
    )
  
  return(p)
}


# Plot of net mineralization against factors
net_mineralization_plot <- plot_mineralization(crop, crop_colors, net_mineralization_1, label = "Crop") + 
  plot_mineralization(addition, addition_colors, net_mineralization_1, label = "Addition") +
  plot_mineralization(treatment, fertilization_colors, net_mineralization_1, label = "Fertilization Level") 

ggsave(
  "net_mineralization.tiff",
  width = 1383,
  height = 422,
  units = "px",
  device = "tiff"
)

# The next plot will show nitrate, ammonium, and net mineralization broken down by
# factors. The rows will be the factors and the columns will be nitrate, ammonium, and 
# net m. We start by building up the rows:
crop_row_plots <- (plot_mineralization(crop, crop_colors, no3n_mg_kg_1, label = "Crop") + 
  theme(legend.position = "none") +
  labs(
    x = "",
    y = "NO3",
    title = "NO3"
  )
) +
(
  plot_mineralization(crop, crop_colors, nh4n_mg_kg_1, label = "Crop") +
    theme(legend.position = "none") +
    labs(
      x = "",
      y = "NH4",
      title = "NH4"
    )
) +
(
  plot_mineralization(crop, crop_colors, net_mineralization_1, label = "Crop") +
    labs(
      x = "",
      y = "Net mineralization",
      title = "Net mineralization"
    )
) 

addition_row_plots <- (plot_mineralization(addition, addition_colors, no3n_mg_kg_1, label = "Addition") + 
                     theme(legend.position = "none") +
                     labs(
                       x = "",
                       y = "NO3",
                     )
) +
  (
    plot_mineralization(addition, addition_colors, nh4n_mg_kg_1, label = "Addition") +
      theme(legend.position = "none") +
      labs(
        x = "",
        y = "NH4",
      )
  ) +
  (
    plot_mineralization(addition, addition_colors, net_mineralization_1, label = "Addition") +
      labs(
        x = "",
        y = "Net mineralization",
      )
  ) 


treatment_row_plots <- (plot_mineralization(treatment, fertilization_colors, no3n_mg_kg_1, label = "N rate") + 
                         theme(legend.position = "none") +
                         labs(
                           x = "DOE",
                           y = "NO3",
                         )
) +
  (
    plot_mineralization(treatment, fertilization_colors, nh4n_mg_kg_1, label = "N rate") +
      theme(legend.position = "none") +
      labs(
        x = "DOE",
        y = "NH4",
      )

  ) +
  (
    plot_mineralization(treatment, fertilization_colors, net_mineralization_1, label = "N rate") +
      labs(
        x = "DOE",
        y = "Net mineralization",
      )
  ) 

# Combining the row plots together...
crop_row_plots / addition_row_plots / treatment_row_plots

# ...and saving them.
ggsave(
  here("figures", "nitrate_ammonium_net_m_by_factors.png"),
  # width = 1383,
  # height = 422,
  # units = "px",
  # device = "tiff"
)


################## Relating mineralization to qPCR data
# Combining the mineralization data with the AmoA data by day

latest_incubation_data <- read.csv("data/Incubation_Biomark-qPCR_all_20230328.csv") %>%
  clean_names() %>% 
  mutate(treatment = str_sub(treatment, end = -2)) %>% 
  mutate(across(
    c(crop, addition),
    str_to_lower
  )) 

mineralization_and_amoa_data <- left_join(
  latest_incubation_data,
  mineralization_data,
  by = c(
    "sample_name" = "initial_name"
    )
)

mineralization_and_amoa_data %>% 
  filter(addition.x == "control") %>% 
  select(net_mineralization_1, f1r2_ave, ave_012, ave_025, ave_039, sum_012025039) %>% 
  cor()


################ ANOVA
# ANOVA results for net mineralization as explained by the different AmoAs
# and DOE. This tests if there's a relationship between the factors and net 
# mineralization.
qpcr_day_anova_results <- aov(
  net_mineralization_1 ~ f1r2_ave * ave_012 * ave_025 * ave_039 * dna_doe.x, 
  data = mineralization_and_amoa_data
)


tidied_qpcr_day_anova_results <- broom::tidy(qpcr_day_anova_results) %>% 
  mutate(sig = p.value < 0.05) %>% 
  mutate(p.value = round(p.value, 2)) %>% 
  filter(term != "Residuals") %>% 
  mutate(term = str_replace(term, "dna_doe.x", "DOE")) %>% 
  mutate(term = str_replace(term, "f1r2_ave", "F1R2")) %>% 
  mutate(term = str_replace(term, "ave_012", "AmoA 12")) %>% 
  mutate(term = str_replace(term, "ave_025", "AmoA 25")) %>% 
  mutate(term = str_replace(term, "ave_039", "AmoA 39")) %>% 
  mutate(order = str_count(term, ":")) 

write.csv(
  tidied_qpcr_day_anova_results,
  here::here("results", "mineralization_qpcr_anova_results.csv"),
  quote = FALSE,
  row.names = FALSE
)

# Looking at net mineralization against main factors
factors_day_anova_results <- aov(
  net_mineralization_1 ~ crop * addition * treatment * dna_doe, 
  data = mineralization_data %>% mutate(dna_doe = as.factor(dna_doe))
)

TukeyHSD(factors_day_anova_results)

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

plot_tidied_anova_mineralization_results <- function(mineralization_anova_results) {
  mineralization_anova_results %>% 
    ggplot(aes(1, term, fill = sig)) +
    geom_tile(color = "white", size = 1) +
    facet_wrap(~ order, scales = "free", ncol = 1) +
    scale_fill_viridis_d(option = "magma", end = 0.75) +
    theme(
      aspect.ratio = 1,
      axis.text.x = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      strip.text = element_blank(),
      axis.ticks = element_blank(),
    ) +
    labs(
      x = "",
      y = "",
      fill = "Significance"
    )
}

factors_min_anova_plot <- plot_tidied_anova_mineralization_results(tidied_factors_day_anova_results)
qpcr_min_anova_plot <- plot_tidied_anova_mineralization_results(tidied_qpcr_day_anova_results) +
  theme(
    legend.position = "none"
  ) +
  scale_y_discrete(position = "right")

factors_min_anova_plot +
  qpcr_min_anova_plot

################ Tukey Post Hoc
# We're doing post-hoc analysis on the factor ANOVA.
tukey_results <- TukeyHSD(factors_day_anova_results)

tukey_results$dna_doe %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  select(rowname, "p.adj" = `p adj`) %>%  
  # mutate(
  #   sig = case_when(
  #     p.adj < 0.05 & p.adj > 0.01 ~ "*",
  #     p.adj < 0.01 & p.adj > 0.001 ~ "**",
  #     p.adj < 0.001 ~ "***",
  #     TRUE ~ "NS"
  #   )
  # ) 
  filter(p.adj < 0.05)

 ################ Multiple linear regression of gene copies against net m

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
  mutate

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
nitrogen_variables <- c("no3n_mg_kg_1", "nh4n_mg_kg_1", "net_mineralization_1")
qpcr_variables <- c(
  "log_012", "ave_012", "log_025", "ave_025", 
  "log_039", "ave_039", "f1r2_log", "f1r2_ave"
)

combined_equations <- c(
  " ~ ave_012 + ave_025 + ave_039 + f1r2_ave",
  " ~ log_012 + log_025 + log_039 + f1r2_log"
)

results <- data.frame(
  dep_var = character(), 
  indep_var = character(),
  r_squared = numeric(),
  p_value = numeric()
)

# Iterate over each dependent variable and independent variable combination
for (dep_var in nitrogen_variables) {
  # Training LMs on one variable at a time
  for (ind_var in qpcr_variables) {
    model <- lm(mineralization_and_amoa_data[[dep_var]] ~ mineralization_and_amoa_data[[ind_var]])
    
    R_squared <- summary(model)$r.squared
    p_value <- summary(model)$coefficients[2, 4]
    
    results <- rbind(results, data.frame(dependent_var = dep_var,
                                         independent_var = ind_var,
                                         R_squared = R_squared,
                                         p_value = p_value))
  }
  
  # Training multiple linear regression on all variables
  
  # Do it again for all variables 
  for (eq in combined_equations) {
    this_formula <- paste0(dep_var, eq)
    model <- lm(this_formula, data = mineralization_and_amoa_data)
    
    R_squared <- summary(model)$adj.r.squared
    p_value <- summary(model)$coefficients[2, 4]
    
    results <- rbind(results, data.frame(dependent_var = dep_var,
                                         independent_var = eq,
                                         R_squared = R_squared,
                                         p_value = p_value))
  }
  
}



# View results
results %>% 
  mutate(significant = if_else(
    p_value < 0.05, "*", ""
  )) %>% 
  mutate(
    R_squared = round(R_squared, 3),
    R_squared = paste0(R_squared, significant)
  ) %>% 
  select(dependent_var, independent_var, R_squared) %>% 
  pivot_wider(names_from = dependent_var, values_from = R_squared)

# Correlations
mineralization_and_amoa_data %>% 
  select(any_of(c(nitrogen_variables, qpcr_variables))) %>% 
  cor() %>% 
  reshape2::melt() %>% 
  filter(
    Var1 %in% qpcr_variables,
    Var2 %in% nitrogen_variables
  ) %>% 
  ggplot(aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white", size = 1) +
  theme(
    text = element_text(size = 16),
    axis.line = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_viridis_c(option = "magma") +
  coord_fixed() +
  labs(
    x = "qPCR",
    y = "",
    fill = "Pearson"
  )


################ Summary statistics of mineralization by AmoA
library(gt)
library(gtsummary)

crop_summary_statistics <- mineralization_and_amoa_data %>% 
  select(
    net_mineralization_1, 
    contains(c("ave", "log")), 
    crop = crop.x
  ) %>% 
  tbl_summary(
    by = crop,
    statistic = list(all_continuous() ~ "{mean}")
  ) %>% 
  add_p() %>% 
  as_tibble() %>% 
  clean_names()

treatment_summary_statistics <- mineralization_and_amoa_data %>% 
  select(
    net_mineralization_1, 
    contains(c("ave", "log")), 
    treatment = treatment.x,
  ) %>% 
  tbl_summary(
    by = treatment,
    statistic = list(all_continuous() ~ "{mean}")
  ) %>% 
  add_p() %>% 
  as_tibble() %>% 
  clean_names()

addition_summary_statistics <- mineralization_and_amoa_data %>% 
  select(
    net_mineralization_1, 
    contains(c("ave", "log")), 
    addition = addition.x,
  ) %>% 
  tbl_summary(
    by = addition,
    statistic = list(all_continuous() ~ "{mean}")
  ) %>% 
  add_p() %>% 
  as_tibble() %>% 
  clean_names()

all_summaries <- bind_cols(
  addition_summary_statistics, 
  crop_summary_statistics, 
  treatment_summary_statistics
) %>% 
  select(-c(6, 10))


write.table(
  all_summaries,
  here::here("results", "qpcr_mineralization_summary_statistics.tsv"),
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)
