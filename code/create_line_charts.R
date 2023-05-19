#!/usr/bin/env Rscript
source("code/setup/setup.R")

plot_amoas <- function(voi, palette, label, raw = FALSE) {
  amoa_dodge <- position_dodge(width = 0.1)

  if (raw == TRUE) {
    data_df <- data.raw
  } else {
    data_df <- data.priming
  }

  data_df <- data_df %>%
    select(doe, {{ voi }}, contains("amoA")) %>%
    pivot_longer(cols = contains("amoA")) %>%
    group_by(doe, {{ voi }}) %>%
    summarise(
      mean = mean(value),
      sd = sd(value),
      se = sd(value) / sqrt(n())
    ) %>%
    mutate(addition = str_to_title({{ voi }})) %>%
    group_by({{ voi }}) %>%
    mutate(cumulative = cumsum(mean))

  p <- data_df %>%
    ggplot(aes(factor(doe), mean, fill = {{ voi }}, group = {{ voi }}, shape = {{ voi }})) +
    geom_hline(yintercept = 23, color = "gray", linetype = "dashed") +
    geom_errorbar(
      aes(ymin = mean - se, ymax = mean + se),
      width = 0.5,
      position = amoa_dodge
    ) +
    geom_line(position = amoa_dodge) +
    scale_y_reverse(breaks = seq(14, 24, 2), limits = c(23, 13)) +
    geom_point(size = 4, position = amoa_dodge) +
    scale_fill_manual(values = palette) +
    scale_shape_manual(values = c(21, 22, 23)) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black")
    ) +
    labs(
      y = "CT",
      x = "Day of experiment",
      fill = label,
      shape = label
    )

  return(p)
}

plot_n2o_lines <- function(LOI, palette, label, show_se = FALSE, cumulative = FALSE, show_legend = FALSE) {
  df <- n2o %>%
    group_by({{ LOI }}, doe) %>%
    summarize(
      mean = mean(N2ON_flux_ug_g_d),
      sd = sd(N2ON_flux_ug_g_d),
      se = sd / sqrt(n())
    ) %>%
    group_by({{ LOI }}) %>%
    mutate(cumulative = cumsum(mean))

  if (cumulative == TRUE) {
    p <- ggplot(df, aes(doe, cumulative, fill = {{ LOI }}, group = {{ LOI }}, shape = {{ LOI }}))
  } else {
    p <- ggplot(df, aes(doe, mean, fill = {{ LOI }}, group = {{ LOI }}, shape = {{ LOI }}))
  }

  p <- p +
    geom_line() +
    geom_point(size = 4) +
    geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
    # scale_y_continuous(breaks = c(-0.05, 0, 0.05, 0.1, 0.15)) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_color_manual(
      values = palette,
      labels = c("Carbon", "Control", "Nitrogen")
    ) +
    scale_fill_manual(
      values = palette,
      labels = c("Carbon", "Control", "Nitrogen")
    ) +
    scale_x_continuous(breaks = day_breaks, labels = as.character(day_breaks), expand = c(0, 0)) +
    labs(
      x = "Day of Experiment",
      y = "n2o flux",
      color = label,
      fill = label
    )


  if (show_se == TRUE) {
    p <- p + geom_errorbar(aes(ymin = mean - se, ymax = mean + se))
  }

  if (show_legend == FALSE) {
    p <- p + theme(legend.position = "none")
  }

  return(p)
}

plot_co2_lines <- function(LOI, VOI, palette, label, show_se = FALSE, cumulative = FALSE, show_legend = FALSE) {
  df <- co2 %>%
    group_by({{ LOI }}, doe_T1) %>%
    summarize(
      mean = mean({{ VOI }}),
      sd = sd({{ VOI }}),
      se = sd / sqrt(n())
    ) %>%
    group_by({{ LOI }}) %>%
    mutate(cumulative = cumsum(mean))

  if (cumulative == TRUE) {
    p <- ggplot(df, aes(doe_T1, cumulative, fill = {{ LOI }}, group = {{ LOI }}, shape = {{ LOI }}))
  } else {
    p <- ggplot(df, aes(doe_T1, mean, fill = {{ LOI }}, group = {{ LOI }}, shape = {{ LOI }}))
  }

  p <- p +
    geom_line() +
    geom_point(size = 4) +
    scale_fill_manual(values = palette) +
    scale_shape_manual(values = c(21, 22, 23)) +
    geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
    scale_x_continuous(breaks = day_breaks, labels = as.character(day_breaks), expand = c(0, 0)) +
    labs(
      x = "Day of Experiment",
      y = "co2 flux",
      color = label,
      fill = label
    )

  if (show_se == TRUE) {
    p <- p + geom_errorbar(aes(ymin = mean - se, ymax = mean + se))
  }

  if (show_legend == FALSE) {
    p <- p + theme(legend.position = "none")
  }

  return(p)
}

plot_co2_lines.addition <- function(cumulative = FALSE) {
  data.df <- tgas %>%
    group_by(addition_T0, doe_T1) %>%
    summarize(
      mean = mean(CO2_flux_ug_g_d),
      sd = sd(CO2_flux_ug_g_d),
      se = sd / sqrt(n())
    ) %>%
    mutate(cumulative = cumsum(mean))

  if (cumulative) {
    p <- data.df %>% ggplot(aes(doe_T1, cumulative, fill = addition_T0, shape = addition_T0))
  } else {
    p <- data.df %>% ggplot(aes(doe_T1, mean, fill = addition_T0, shape = addition_T0))
  }

  p <- p +
    geom_line() +
    geom_point(size = 4) +
    scale_fill_manual(values = addition_colors) +
    scale_shape_manual(values = c(21, 22, 23)) +
    geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
    scale_x_continuous(breaks = day_breaks, labels = as.character(day_breaks), expand = c(0, 0)) +
    labs(
      x = "Day of Experiment",
      y = "co2 flux",
      color = "Addition",
      fill = "Addition"
    )

  return(p)
}

data.priming <- read.csv("data/priming_amoA_deltaCt.csv", header = T) %>%
  rename(sample_id = X) %>%
  mutate(addition = str_to_title(addition)) %>%
  mutate(across(c(addition, fert_level, crop), as.factor))

data.priming.long <- data.priming %>%
  pivot_longer(cols = starts_with("amoA"), names_to = "amoA", values_to = "deltaCT")

data.priming.long$sample_id <- fct_reorder(data.priming.long$sample_id, parse_number(data.priming.long$sample_id))

df <- data.priming[, -1]
rownames(df) <- data.priming[, 1]

metadata <- df %>%
  select(fert_level:field_rep) %>%
  mutate(across(everything(), as.factor))

amoa_counts <- df %>%
  select(starts_with("amoA"))

n2o <- readr::read_csv("data/N2O.csv") %>%
  select(-1) %>%
  mutate(crop_within_block = interaction(crop, rep)) %>%
  mutate(field_sample = interaction(crop, fert, rep)) %>%
  mutate(addition = case_when(
    addition == "Cntrl" ~ "Control",
    addition == "+C" ~ "Carbon",
    TRUE ~ "Nitrogen"
  )) %>%
  mutate(addition = as.factor(addition)) %>%
  mutate(crop = as.factor(str_to_lower(crop))) %>%
  filter(fert != "112N") %>%
  mutate(fert = as.factor(str_sub(fert, end = -2)))

co2 <- read.csv("data/priming_calc.csv") %>%
  filter(fert_T1 != "112N") %>%
  mutate(fert_T1 = as.factor(str_sub(fert_T1, end = -2))) %>%
  mutate(crop_T1 = as.factor(str_to_lower(crop_T1)))

tgas <- read.csv("data/tgas1.csv") %>%
  mutate(addition_T0 = case_when(
    addition_T0 == "Cntrl" ~ "Control",
    addition_T0 == "+C" ~ "Carbon",
    TRUE ~ "Nitrogen"
  ))

amoa_abundance_addition <- plot_amoas(addition, addition_colors, "Addition")
amoa_abundance_fert <- plot_amoas(fert_level, fertilization_colors, "Fertilization Level")
amoa_abundance_crop <- plot_amoas(crop, crop_colors, "Crop")

n20_flux_crop <- plot_n2o_lines(crop, crop_colors, "Crop", show_se = FALSE)
n20_flux_fert <- plot_n2o_lines(fert, fertilization_colors, "Fertilization Level", show_se = FALSE)
n20_flux_addition <- plot_n2o_lines(addition, addition_colors, "Addition", show_se = FALSE)

co2_flux_crop <- plot_co2_lines(crop_T1, Cntrl_CO2_flux_ug_g_d_soil, crop_colors, "Crop", show_se = FALSE)
co2_flux_fert <- plot_co2_lines(fert_T1, Cntrl_CO2_flux_ug_g_d_soil, fertilization_colors, "Fertilization Level", show_se = FALSE)
co2_flux_addition <- plot_co2_lines.addition(cumulative = FALSE)

fert_plot <- (n20_flux_fert + theme(legend.position = "none")) +
  (co2_flux_fert + theme(legend.position = "none")) +
  amoa_abundance_fert

ggsave("figures/fert_plot.png", fert_plot, height = 4, width = 11)

crop_plot <- (n20_flux_crop + theme(legend.position = "none")) +
  (co2_flux_crop + theme(legend.position = "none")) +
  amoa_abundance_crop

ggsave("figures/crop_plot.png", crop_plot, height = 4, width = 11)

addition_plot <- (n20_flux_addition + theme(legend.position = "none")) +
  (co2_flux_addition + theme(legend.position = "none")) +
  amoa_abundance_addition

ggsave("figures/addition_plot.png", addition_plot, height = 4, width = 11)

n20_flux_crop.cumulative <- plot_n2o_lines(crop, crop_colors, "Crop", show_se = FALSE, cumulative = TRUE)
n20_flux_fert.cumulative <- plot_n2o_lines(fert, fertilization_colors, "Fertilization Level", show_se = FALSE, cumulative = TRUE)
n20_flux_addition.cumulative <- plot_n2o_lines(addition, addition_colors, "Addition", show_se = FALSE, cumulative = TRUE)

co2_flux_crop.cumulative <- plot_co2_lines(crop_T1, Cntrl_CO2_flux_ug_g_d_soil, crop_colors, "Crop", show_se = FALSE, cumulative = TRUE)
co2_flux_fert.cumulative <- plot_co2_lines(fert_T1, Cntrl_CO2_flux_ug_g_d_soil, fertilization_colors, "Fertilization Level", show_se = FALSE, cumulative = TRUE)
co2_flux_addition.cumulative <- plot_co2_lines.addition(cumulative = TRUE)

fert_cumulative_plot <- (n20_flux_fert.cumulative + theme(legend.position = "none")) +
  (co2_flux_fert.cumulative)

fert_cumulative_plot

ggsave("figures/fert_plot_cumulative.png", fert_cumulative_plot, height = 4, width = 11)

crop_cumulative_plot <- (n20_flux_crop.cumulative + theme(legend.position = "none")) +
  (co2_flux_crop.cumulative)

ggsave("figures/crop_plot_cumulative.png", crop_cumulative_plot, height = 4, width = 11)

addition_cumulative_plot <- (n20_flux_addition.cumulative + theme(legend.position = "none")) +
  (co2_flux_addition.cumulative)

co2_flux_addition.cumulative + n20_flux_addition.cumulative

ggsave("figures/addition_plot_cumulative.png", addition_cumulative_plot, height = 4, width = 11)
