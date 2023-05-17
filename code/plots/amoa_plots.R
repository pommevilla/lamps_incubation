#!/usr/bin/env Rscript
# ---------------------------
# Plots net mineralization figures
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup.R")

# Read in data
qpcr_data <- read.csv(here("data", "Incubation_Biomark-qPCR_all_20230328.csv")) %>%
  mutate(Treatment = case_when(
    Treatment == "100N" ~ "112N",
    Treatment == "300N" ~ "336N",
    TRUE ~ Treatment
  ))


################### Plots
# This function will take in a grouping variable and a palette and will
# plot all the amoAs by that factor.
plot_amoa_qpcr <- function(voi, palette, label, free_y = FALSE) {
  this_dodge <- position_dodge(width = 0.1)

  data_df <- qpcr_data %>%
    pivot_longer(contains("ave")) %>%
    group_by(DNA_DOE, name, {{ voi }}) %>%
    summarize(
      mean_n = mean(value),
      sd = sd(value),
      se = sd(value) / sqrt(n())
    ) %>%
    mutate(
      name = if_else(
        str_detect(name, "ave_"),
        str_replace(name, "ave_", "amoA "),
        "F1R2"
      )
    ) %>%
    ungroup()

  p <- data_df %>%
    ggplot(aes(DNA_DOE, mean_n, fill = {{ voi }}, group = {{ voi }}, shape = {{ voi }})) +
    geom_errorbar(
      aes(ymin = mean_n - se, ymax = mean_n + se),
      width = 1,
      position = this_dodge
    ) +
    geom_line() +
    geom_point(size = 4, position = this_dodge) +
    scale_fill_manual(values = palette) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_x_continuous(breaks = qpcr_day_breaks) +
    scale_y_continuous(labels = scales::scientific) +
    # scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      plot.title = element_text(hjust = 0.5),
      aspect.ratio = 1,
      strip.text = element_markdown(color = "black"),
      strip.background = element_blank(),
      axis.title.y = element_markdown()
    ) +
    labs(
      y = "",
      x = "",
      fill = label,
      shape = label
    )

  if (free_y) {
    p <- p + facet_wrap(~name, scales = "free_y", nrow = 1)
  } else {
    p <- p + facet_wrap(~name, nrow = 1)
  }

  return(p)
}

plot_amoa_qpcr(Addition, addition_colors, "Addition")

# Plot of net mineralization against factors
qpcr_crop_plot <- plot_amoa_qpcr(Crop, crop_colors, "Crop")
qpcr_addition_plot <- plot_amoa_qpcr(Addition, addition_colors, "Addition") +
  theme(strip.text = element_blank()) +
  labs(y = gcn_unit)
qpcr_treatment_plot <- plot_amoa_qpcr(Treatment, fertilization_colors, "N Rate") +
  theme(strip.text = element_blank()) +
  labs(x = "Day")

(qpcr_crop_plot / qpcr_addition_plot / qpcr_treatment_plot)

ggsave(
  "figures/amoa_qpcr/qpcr_by_factors.png",
  width = 3000,
  height = 1900,
  units = "px"
)

qpcr_crop_plot <- plot_amoa_qpcr(Crop, crop_colors, "Crop", free_y = TRUE)
qpcr_addition_plot <- plot_amoa_qpcr(Addition, addition_colors, "Addition", free_y = TRUE) +
  theme(strip.text = element_blank()) +
  labs(y = gcn_unit)
qpcr_treatment_plot <- plot_amoa_qpcr(Treatment, fertilization_colors, "N Rate", free_y = TRUE) +
  theme(strip.text = element_blank()) +
  labs(x = "Day")

(qpcr_crop_plot / qpcr_addition_plot / qpcr_treatment_plot)

ggsave(
  "figures/amoa_qpcr/qpcr_by_factors_free_y.png",
  width = 3000,
  height = 1900,
  units = "px"
)

# ggsave(
#   "figures/qpcr_by_factors.tiff",
#   width = 1383,
#   height = 1383,
#   units = "px",
#   device = "tiff"
# )
