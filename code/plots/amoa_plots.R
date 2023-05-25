#!/usr/bin/env Rscript
# ---------------------------
# Plots net mineralization figures
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup/setup.R")

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



#################### qPCR barchart
# Across all treatments and days, what is the average abundance of each amoA?
# The `how_much_more` column is how many more times abundant each primer is
# relative to F1R2, the classic literature primer.
# avg_qpcr_numbers <-

avg_qpcr_numbers <- qpcr_data %>%
  select(contains("ave")) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    se = sd / sqrt(n())
  ) %>%
  mutate(
    name = case_when(
      name == "F1R2_ave" ~ "F1R2",
      TRUE ~ str_replace(name, "ave_", "amoA ")
    )
  ) %>%
  mutate(
    how_much_more = round(mean / first(mean), 1),
    how_much_more = paste0(how_much_more, "x")
  ) %>%
  mutate(
    rel_abund = mean / sum(mean)
  )

avg_qpcr_numbers %>%
  ggplot(aes(name, mean)) +
  geom_col(width = 0.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  geom_label(
    data = avg_qpcr_numbers %>% filter(name != "F1R2"),
    aes(x = name, y = mean, label = how_much_more),
    nudge_y = -0.25e8
  ) +
  labs(
    x = "Primer set",
    y = gcn_unit,
  ) +
  theme(
    panel.border = element_blank(),
    axis.title.y = element_markdown(),
  ) +
  scale_y_continuous(
    labels = scales::scientific,
    expand = expansion(0, 0.2)
  )

ggsave(
  here::here("figures/amoa_qpcr/amoA_qpcr_barchart.png"),
)
