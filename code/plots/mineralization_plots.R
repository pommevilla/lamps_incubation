#!/usr/bin/env Rscript
# ---------------------------
# Plots net mineralization figures
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup.R")

# Read in data
mineralization_data <- read.csv(here("data/prepped_data", "mineralization_data.csv"))

################## Plots
# This function will take in a grouping variable, a palette, and a y variable
# and plot that y variable against time. The optional zero line is to add a horizontal
# line at x = 0.
plot_mineralization <- function(voi, palette, net_min_var, label, zero_line = FALSE) {
  this_dodge <- position_dodge(width = 0.1)

  data_df <- mineralization_data %>%
    group_by(day, {{ voi }}) %>%
    summarize(
      mean_n = mean({{ net_min_var }}),
      sd = sd({{ net_min_var }}),
      se = sd({{ net_min_var }}) / sqrt(n())
    ) %>%
    ungroup()

  p <- data_df %>%
    ggplot(aes(day, mean_n, fill = {{ voi }}, group = {{ voi }}, shape = {{ voi }}))

  if (zero_line) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "#848884")
  }

  p <- p +
    geom_errorbar(
      aes(ymin = mean_n - se, ymax = mean_n + se),
      width = 1,
      position = this_dodge
    ) +
    geom_line() +
    geom_point(size = 4, position = this_dodge) +
    scale_fill_manual(values = palette) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_x_continuous(breaks = day_breaks) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      aspect.ratio = 1
    ) +
    labs(
      x = "",
      y = "",
      fill = label,
      shape = label
    )

  return(p)
}

# Plot of net mineralization against factors
net_mineralization_plot <- plot_mineralization(crop, crop_colors, net_min_rate_rel, label = "Crop") +
  plot_mineralization(addition, addition_colors, net_min_rate_rel, label = "Addition") +
  plot_mineralization(treatment, fertilization_colors, net_min_rate_rel, label = "Fertilization Level")

net_mineralization_plot

ggsave(
  "figures/net_mineralization.png"
)

ggsave(
  "net_mineralization.tiff",
  width = 1383,
  height = 422,
  units = "px",
  device = "tiff"
)

# The next plot will show nitrate, ammonium, and net mineralization broken down by
# factors. The rows will be the factors and the columns will be nitrate, ammonium, and
# net m rate.
# Helper function to plot individual rows nicely.
# TODO: Figure out why purrr::map doesn't want to cooperate with this entire function
plot_nnn_row <- function(voi, palette, label, titles = FALSE, xaxis = FALSE) {
  p1 <- plot_mineralization({{ voi }}, palette, no3n_mg_kg_1, label = label, zero_line = TRUE) +
    theme(
      legend.position = "none",
      plot.title = element_markdown(),
      plot.subtitle = element_markdown()
    )

  p2 <- plot_mineralization({{ voi }}, palette, nh4n_mg_kg_1, label = label, zero_line = TRUE) +
    theme(
      legend.position = "none",
      plot.title = element_markdown(),
      plot.subtitle = element_markdown()
    )

  p3 <- plot_mineralization({{ voi }}, palette, net_min_rate_rel, label = label, zero_line = TRUE)


  if (titles) {
    p1 <- p1 +
      labs(
        title = nitrate_label,
        subtitle = flux_units
      )
    p2 <- p2 +
      labs(
        title = ammonia_label,
        subtitle = flux_units
      )
    p3 <- p3 +
      labs(
        title = "Net mineralization rate (rel.)",
        subtitle = paste(flux_units, per_day_unit)
      ) +
      theme(
        plot.subtitle = element_markdown()
      )
  }

  if (xaxis) {
    p1 <- p1 + labs(x = "Day")
    p2 <- p2 + labs(x = "Day")
    p3 <- p3 + labs(x = "Day")
  }

  p <- p1 + p2 + p3 + plot_layout(ncol = 3)

  return(p)
}

# Putting rows together...
wrap_plots(
  plot_nnn_row(crop, crop_colors, "Crop", titles = TRUE),
  plot_nnn_row(addition, addition_colors, "Addition"),
  plot_nnn_row(treatment, fertilization_colors, "N Rate", xaxis = TRUE),
  nrow = 3
)

# ...and saving them.
ggsave(
  "figures/mineralization/nitrate_ammonium_net_m_by_factors.png",
  width = 3000,
  height = 2900,
  units = "px",
)

## Absolute and relative net mineralization and nitrification
plot_armn_row <- function(voi, palette, label, titles = FALSE, xaxis = FALSE) {
  p1 <- plot_mineralization({{ voi }}, palette, net_min_rate_rel, label = label, zero_line = TRUE) +
    theme(
      legend.position = "none",
      plot.title = element_markdown(),
      plot.subtitle = element_markdown()
    )

  p2 <- plot_mineralization({{ voi }}, palette, net_min_rate_abs, label = label, zero_line = TRUE) +
    theme(
      legend.position = "none",
      plot.subtitle = element_markdown()
    )

  p3 <- plot_mineralization({{ voi }}, palette, net_nitr_rate_rel, label = label, zero_line = TRUE) +
    theme(
      legend.position = "none",
      plot.subtitle = element_markdown()
    )

  p4 <- plot_mineralization({{ voi }}, palette, net_nitr_rate_abs, label = label, zero_line = TRUE) +
    theme(
      plot.subtitle = element_markdown()
    )

  # This sucks
  # TODO: Figure out why purrr::map doesn't like this
  if (titles) {
    p1 <- p1 + labs(
      title = "Net mineralization rate (rel.)",
      subtitle = paste(flux_units, per_day_unit)
    )
    p2 <- p2 + labs(
      title = "Net mineralization rate (abs.)",
      subtitle = paste(flux_units, per_day_unit)
    )
    p3 <- p3 + labs(
      title = "Net nitrification rate (rel.)",
      subtitle = paste(flux_units, per_day_unit)
    )
    p4 <- p4 + labs(
      title = "Net nitrification rate (abs.)",
      subtitle = paste(flux_units, per_day_unit)
    )
  }

  # So does this
  if (xaxis) {
    p1 <- p1 + labs(x = "Day")
    p2 <- p2 + labs(x = "Day")
    p3 <- p3 + labs(x = "Day")
    p4 <- p4 + labs(x = "Day")
  }

  p <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)

  return(p)
}

wrap_plots(
  plot_armn_row(crop, crop_colors, "Crop", titles = TRUE),
  plot_armn_row(addition, addition_colors, "Addition"),
  plot_armn_row(treatment, fertilization_colors, "N rate", xaxis = TRUE),
  nrow = 3
)

ggsave(
  "figures/mineralization/net_min_nitr_abs_rel.png",
  width = 4300,
  height = 3000,
  units = "px",
)
