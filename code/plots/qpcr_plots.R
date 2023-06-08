#!/usr/bin/env Rscript
# ---------------------------
# Plots net mineralization figures
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup/setup.R")

# Read in data
qpcr_data <- read.csv("data/prepped_data/qpcr_data.csv")

################### Plots
# Plots qPCR data by factor over days.
#   voi: variable of interest. The name of a factor (unquoted) to plot. Can be Crop, Treatment, or Addition.
#   palette: the palette used for the colors. These are defined in setup.R. Can be crop_colors, fertilization_colors, or addition_colors.
#   plot_vars: list of variables to plot. Can be a list of any column names (quoted) in qpcr_data.
#     The variables will be reordered according to this vector. There are some predefined vectors in setup.R.
#   label: the label for the legend.
#   free_y: whether or not to free the y-axis for each facet. Defaults to FALSE.
plot_qpcr <- function(voi, palette, plot_vars, label, free_y = FALSE) {
  this_dodge <- position_dodge(width = 0.1)

  data_df <- qpcr_data %>%
    pivot_longer(any_of(plot_vars)) %>%
    group_by(Day, name, {{ voi }}) %>%
    summarize(
      mean_n = mean(value),
      sd = sd(value),
      se = sd(value) / sqrt(n())
    ) %>%
    ungroup() %>%
    mutate(
      name = factor(name, levels = plot_vars, labels = make_nice_qpcr_names(plot_vars)),
    )


  p <- data_df %>%
    ggplot(aes(Day, mean_n, fill = {{ voi }}, group = {{ voi }}, shape = {{ voi }})) +
    geom_errorbar(
      aes(ymin = mean_n - se, ymax = mean_n + se),
      width = 1,
      position = this_dodge
    ) +
    geom_line() +
    geom_point(size = 2, position = this_dodge) +
    scale_fill_manual(values = palette) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_x_continuous(breaks = qpcr_day_breaks) +
    scale_y_continuous(labels = scales::scientific) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      aspect.ratio = 1,
      strip.text = element_markdown(color = "black"),
      strip.background = element_blank(),
      axis.title.y = element_markdown()
    ) +
    labs(
      y = "",
      x = "",
    )

  if (free_y) {
    p <- p + facet_wrap(~name, scales = "free_y", nrow = 1)
  } else {
    p <- p + facet_wrap(~name, nrow = 1)
  }

  return(p)
}

######### amoA plots
# non-free y
amoa_crop_plot <- plot_qpcr(Crop, crop_colors, ave_amoa_qpcr_variables, "Crop") +
  labs(y = gcn_unit)
amoa_addition_plot <- plot_qpcr(Addition, addition_colors, ave_amoa_qpcr_variables, "Addition") +
  theme(strip.text = element_blank()) +
  labs(y = gcn_unit)
amoa_treatment_plot <- plot_qpcr(Treatment, fertilization_colors, ave_amoa_qpcr_variables, "N Rate") +
  theme(strip.text = element_blank()) +
  labs(x = "Day", y = gcn_unit)

(amoa_crop_plot / amoa_addition_plot / amoa_treatment_plot)

ggsave(
  "figures/qpcr/amoa_qpcr_by_factors.png",
  width = 3000,
  height = 1900,
  units = "px"
)

# Free y
amoa_crop_plot <- plot_qpcr(Crop, crop_colors, ave_amoa_qpcr_variables, "Crop", free_y = TRUE) +
  labs(y = gcn_unit)
amoa_addition_plot <- plot_qpcr(Addition, addition_colors, ave_amoa_qpcr_variables, "Addition", free_y = TRUE) +
  theme(strip.text = element_blank()) +
  labs(y = gcn_unit)
amoa_treatment_plot <- plot_qpcr(Treatment, fertilization_colors, ave_amoa_qpcr_variables, "N Rate", free_y = TRUE) +
  theme(strip.text = element_blank()) +
  labs(x = "Day", y = gcn_unit)

(amoa_crop_plot / amoa_addition_plot / amoa_treatment_plot)

ggsave(
  "figures/qpcr/amoa_qpcr_by_factors_free_y.png",
  width = 4000,
  height = 2100,
  units = "px"
)

######### norB plots
# non-free y
norb_crop_plot <- plot_qpcr(Crop, crop_colors, ave_norb_qpcr_variables, "Crop") +
  labs(y = gcn_unit)
norb_addition_plot <- plot_qpcr(Addition, addition_colors, ave_norb_qpcr_variables, "Addition") +
  theme(strip.text = element_blank()) +
  labs(y = gcn_unit)
norb_treatment_plot <- plot_qpcr(Treatment, fertilization_colors, ave_norb_qpcr_variables, "N Rate") +
  theme(strip.text = element_blank()) +
  labs(x = "Day", y = gcn_unit)

(norb_crop_plot / norb_addition_plot / norb_treatment_plot)

ggsave(
  "figures/qpcr/norb_qpcr_by_factors.png",
  width = 3000,
  height = 1900,
  units = "px"
)

# Free y
norb_crop_plot <- plot_qpcr(Crop, crop_colors, ave_norb_qpcr_variables, "Crop", free_y = TRUE) +
  labs(y = gcn_unit)
norb_addition_plot <- plot_qpcr(Addition, addition_colors, ave_norb_qpcr_variables, "Addition", free_y = TRUE) +
  theme(strip.text = element_blank()) +
  labs(y = gcn_unit)
norb_treatment_plot <- plot_qpcr(Treatment, fertilization_colors, ave_norb_qpcr_variables, "N Rate", free_y = TRUE) +
  theme(strip.text = element_blank()) +
  labs(x = "Day", y = gcn_unit)

(norb_crop_plot / norb_addition_plot / norb_treatment_plot)

ggsave(
  "figures/qpcr/norb_qpcr_by_factors_free_y.png",
  width = 4000,
  height = 2200,
  units = "px"
)

#################### qPCR barchart
# How much more are we detecting with our primers?
# The `how_much_more` column is how many more times abundant each primer is
# relative to F1R2, the classic literature primer.

# Plots how much more abundant each of the plot_vars is relative to the
# baseline_qpcr variable
calc_rel_abundance <- function(qpcr_df, plot_vars, baseline_qpcr, gene) {
  avg_qpcr_numbers <- qpcr_df %>%
    select(any_of(plot_vars)) %>%
    pivot_longer(everything()) %>%
    group_by(name) %>%
    summarise(
      mean = mean(value),
      sd = sd(value),
      se = sd / sqrt(n())
    )

  baseline_abundance <- avg_qpcr_numbers %>%
    filter(name == baseline_qpcr) %>%
    pull(mean)

  avg_qpcr_numbers <- avg_qpcr_numbers %>%
    mutate(
      how_much_more = round(mean / baseline_abundance, 1),
      how_much_more = paste0(how_much_more, "x"),
      rel_abund = mean / sum(mean)
    ) %>%
    mutate(
      gene = gene,
      baseline = if_else(name == baseline_qpcr, "baseline", "not baseline"),
      name = factor(name, levels = plot_vars, labels = make_nice_qpcr_names(plot_vars))
    )

  return(avg_qpcr_numbers)
}

norb_rel_abundances <- calc_rel_abundance(qpcr_data, ave_norb_qpcr_variables, "cnorB", "norB")
amoa_rel_abundances <- calc_rel_abundance(qpcr_data, ave_amoa_qpcr_variables, "F1R2_ave", "amoA")


combined_abundances <- bind_rows(
  amoa_rel_abundances,
  norb_rel_abundances
)

plot_rel_abundances <- function(abundance_df) {
  abundance_df %>%
    ggplot(aes(name, mean, fill = baseline)) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    geom_col(width = 0.5) +
    geom_label(
      inherit.aes = FALSE,
      data = abundance_df %>% filter(baseline != "baseline"),
      aes(x = name, y = mean, label = how_much_more),
      nudge_y = -0.5e8,
      size = 2
    ) +
    labs(
      x = "",
      y = "",
      fill = ""
    ) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      axis.title.y = element_markdown(),
      strip.text = element_markdown(color = "black"),
      strip.background = element_blank(),
    ) +
    scale_y_continuous(
      labels = scales::scientific,
      expand = expansion(0, 0.4)
    ) +
    scale_fill_manual(
      values = c("baseline" = "grey", "not baseline" = "black"),
      labels = c("baseline" = "Existing primers", "not baseline" = "MFP primers")
    ) +
    facet_grid(~gene, scales = "free_x")
}

plot_rel_abundances(combined_abundances) +
  labs(y = paste0(gcn_unit, "<br>(abundance rel. to classic primers  in label)"))

ggsave(
  here::here("figures/qpcr/qpcr_rel_abundances_barchart.png"),
  height = 1900,
  width = 3000,
  units = "px"
)


########## Plotting abundances by addition within N-rate/crop type
plot_addition_qpcr_in_group <- function(voi, plot_vars, free_y = FALSE) {
  this_dodge <- position_dodge(width = 0.1)

  data_df <- qpcr_data %>%
    pivot_longer(any_of(plot_vars)) %>%
    group_by(Addition, Day, name, {{ voi }}) %>%
    summarize(
      mean_n = mean(value),
      sd = sd(value),
      se = sd(value) / sqrt(n())
    ) %>%
    mutate(
      name = factor(name, levels = plot_vars),
    ) %>%
    ungroup()

  p <- data_df %>%
    ggplot(aes(Day, mean_n, fill = Addition, group = Addition, shape = Addition)) +
    geom_errorbar(
      aes(ymin = mean_n - se, ymax = mean_n + se),
      width = 1,
      position = this_dodge
    ) +
    geom_line() +
    geom_point(size = 4, position = this_dodge) +
    scale_fill_manual(values = addition_colors) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_x_continuous(breaks = qpcr_day_breaks) +
    scale_y_continuous(labels = scales::scientific) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      aspect.ratio = 1,
      strip.text = element_markdown(color = "black"),
      strip.background = element_blank(),
      axis.title.y = element_markdown()
    ) +
    labs(
      y = "",
      x = "",
      fill = "Addition",
      shape = "Addition"
    )

  return(p)
}

plot_addition_qpcr_in_group(Treatment, ave_norb_qpcr_variables, free_y = TRUE) +
  facet_grid(Treatment ~ name)

ggsave(
  here::here("figures/qpcr/norb_addition_in_nrate.png"),
  width = 4500,
  height = 3000,
  units = "px"
)

plot_addition_qpcr_in_group(Treatment, ave_amoa_qpcr_variables, free_y = TRUE) +
  facet_grid(Treatment ~ name)

ggsave(
  here::here("figures/qpcr/amoa_addition_in_nrate.png"),
  width = 4500,
  height = 3000,
  units = "px"
)

plot_addition_qpcr_in_group(Crop, ave_norb_qpcr_variables, free_y = TRUE) +
  facet_grid(Crop ~ name)

ggsave(
  here::here("figures/qpcr/norb_addition_in_crop.png"),
  width = 4500,
  height = 3000,
  units = "px"
)

plot_addition_qpcr_in_group(Crop, ave_amoa_qpcr_variables, free_y = TRUE) +
  facet_grid(Crop ~ name)

ggsave(
  here::here("figures/qpcr/amoa_addition_in_crop.png"),
  width = 4500,
  height = 3000,
  units = "px"
)
