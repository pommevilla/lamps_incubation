#!/usr/bin/env Rscript
# ---------------------------
# Plots net mineralization figures
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup.R")

ammonium_label <- "$NH_4$"

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

net_mineralization_plot +
  theme(
    text = element_text(family = "Times New Roman")
  )

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
    plot_mineralization(crop, crop_colors, net_mineralization_1, label = "Crop", zero_line = TRUE) +
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
      theme(
        legend.position = "none",
        axis.title.x = element_markdown()
      ) +
      labs(
        x = "",
        y = "NH4",
      )
  ) +
  (
    plot_mineralization(addition, addition_colors, net_mineralization_1, label = "Addition", zero_line = TRUE) +
      labs(
        x = "",
        y = "Net mineralization",
      )
  )


treatment_row_plots <- (plot_mineralization(treatment, fertilization_colors, no3n_mg_kg_1, label = "N rate") +
  theme(legend.position = "none") +
  labs(
    x = "Day",
    y = "NO3",
  )
) +
  (
    plot_mineralization(treatment, fertilization_colors, nh4n_mg_kg_1, label = "N rate") +
      theme(legend.position = "none") +
      labs(
        x = "Day",
        y = "NH4",
      )

  ) +
  (
    plot_mineralization(treatment, fertilization_colors, net_mineralization_1, label = "N rate", zero_line = TRUE) +
      labs(
        x = "Day",
        y = "Net mineralization",
      )
  )

# Combining the row plots together...
crop_row_plots / addition_row_plots / treatment_row_plots

# ...and saving them.
ggsave(
  "figures/nitrate_ammonium_net_m_by_factors.png",
  width = 4000,
  height = 4000,
  units = "px",
)
