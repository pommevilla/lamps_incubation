#!/usr/bin/env Rscript
# ---------------------------
# Scatterplots of amoA abundance against variables
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup/setup.R")

# Read in data
mineralization_and_qpcr_data <- read.csv(here("data/prepped_data", "mineralization_and_qpcr_data.csv")) %>%
  rename(
    Crop = crop,
    Treatment = treatment,
    Addition = addition
  )


################## Plots
# This function will take in a grouping variable, a palette, and a y variable
# and plot that y variable against time. The optional zero line is to add a horizontal
# line at x = 0.

plot_amoa_min_scatter <- function(amoa_oi, min_voi, foi, palette, show_legend = FALSE) {
  # Get info R^2 and significance for linear model
  this_formula <- as.formula(paste(deparse(substitute(min_voi)), "~", deparse(substitute(amoa_oi))))
  this_model <- lm(this_formula, data = mineralization_and_qpcr_data)
  r_squared <- summary(this_model)$r.squared %>%
    round(., 3)
  p_value <- summary(this_model)$coefficients[2, 4] %>%
    round(., 3)
  p_sig <- get_p_sig(p_value)
  r_squared_label <- paste0("R<sup>2</sup> = ", r_squared, p_sig)

  # Get Pearson's R and significance
  this_cor <- cor.test(
    mineralization_and_qpcr_data %>% pull({{ amoa_oi }}),
    mineralization_and_qpcr_data %>% pull({{ min_voi }})
  )
  pearson_r <- this_cor$estimate %>%
    round(., 3)
  p_sig <- get_p_sig(this_cor$p.value)
  pearson_r_label <- paste0("ρ = ", pearson_r, p_sig)

  # Plotting
  p <- mineralization_and_qpcr_data %>%
    ggplot(aes({{ amoa_oi }}, {{ min_voi }}, fill = {{ foi }}, shape = {{ foi }})) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#848884") +
    geom_point(size = 4) +
    scale_fill_manual(values = palette) +
    scale_shape_manual(values = c(21, 22, 23)) +
    geom_abline(
      intercept = coef(this_model)[1],
      slope = coef(this_model)[2]
    ) +
    annotate(
      geom = "richtext",
      x = Inf, y = Inf, hjust = 1, vjust = 1,
      label = paste0(r_squared_label, "<br>", pearson_r_label),
    ) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      plot.title = element_markdown(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      aspect.ratio = 1,
      axis.title.y = element_markdown()
    ) +
    labs(
      x = ""
    )


  if (show_legend == FALSE) {
    p <- p + theme(legend.position = "none")
  }

  return(p)
}


############### amoa_12
amoa_12_crop_row <- wrap_plots(
  plot_amoa_min_scatter(ave_012, net_min_rate_rel, Crop, crop_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
      title = "Net mineralization rate (rel.)"
    ),
  plot_amoa_min_scatter(ave_012, net_min_rate_abs, Crop, crop_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
      title = "Net mineralization rate (abs.)"
    ),
  plot_amoa_min_scatter(ave_012, no3n_mg_kg_1, Crop, crop_colors) +
    labs(
      y = flux_units,
      title = nitrate_label
    ),
  plot_amoa_min_scatter(ave_012, nh4n_mg_kg_1, Crop, crop_colors, show_legend = TRUE) +
    labs(
      y = flux_units,
      title = ammonia_label
    ),
  nrow = 1
)


amoa_12_addition_row <- wrap_plots(
  plot_amoa_min_scatter(ave_012, net_min_rate_rel, Addition, addition_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(ave_012, net_min_rate_abs, Addition, addition_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(ave_012, no3n_mg_kg_1, Addition, addition_colors) +
    labs(
      y = flux_units,
    ),
  plot_amoa_min_scatter(ave_012, nh4n_mg_kg_1, Addition, addition_colors, show_legend = TRUE) +
    labs(
      y = flux_units,
    ),
  nrow = 1
)

amoa_12_fert_row <- wrap_plots(
  plot_amoa_min_scatter(ave_012, net_min_rate_rel, Treatment, fertilization_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(ave_012, net_min_rate_abs, Treatment, fertilization_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(ave_012, no3n_mg_kg_1, Treatment, fertilization_colors) +
    labs(
      y = flux_units,
    ),
  plot_amoa_min_scatter(ave_012, nh4n_mg_kg_1, Treatment, fertilization_colors, show_legend = TRUE) +
    labs(
      y = flux_units,
    ),
  nrow = 1
)

amoa_12_crop_row / amoa_12_addition_row / amoa_12_fert_row + plot_annotation(
  caption = "amoA 12 (gene copies g<sup>-1</sup>)",
  theme = theme(plot.caption = element_markdown(hjust = 0.5, size = 10))
)

ggsave(
  "figures/scatterplots/amoa12_vs_min.png",
  width = 4500,
  height = 3000,
  units = "px"
)


############### amoa_25
amoa_25_crop_row <- wrap_plots(
  plot_amoa_min_scatter(ave_025, net_min_rate_rel, Crop, crop_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
      title = "Net mineralization rate (rel.)"
    ),
  plot_amoa_min_scatter(ave_025, net_min_rate_abs, Crop, crop_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
      title = "Net mineralization rate (abs.)"
    ),
  plot_amoa_min_scatter(ave_025, no3n_mg_kg_1, Crop, crop_colors) +
    labs(
      y = flux_units,
      title = nitrate_label
    ),
  plot_amoa_min_scatter(ave_025, nh4n_mg_kg_1, Crop, crop_colors, show_legend = TRUE) +
    labs(
      y = flux_units,
      title = ammonia_label
    ),
  nrow = 1
)


amoa_25_addition_row <- wrap_plots(
  plot_amoa_min_scatter(ave_025, net_min_rate_rel, Addition, addition_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(ave_025, net_min_rate_abs, Addition, addition_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(ave_025, no3n_mg_kg_1, Addition, addition_colors) +
    labs(
      y = flux_units,
    ),
  plot_amoa_min_scatter(ave_025, nh4n_mg_kg_1, Addition, addition_colors, show_legend = TRUE) +
    labs(
      y = flux_units,
    ),
  nrow = 1
)

amoa_25_fert_row <- wrap_plots(
  plot_amoa_min_scatter(ave_025, net_min_rate_rel, Treatment, fertilization_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(ave_025, net_min_rate_abs, Treatment, fertilization_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(ave_025, no3n_mg_kg_1, Treatment, fertilization_colors) +
    labs(
      y = flux_units,
    ),
  plot_amoa_min_scatter(ave_025, nh4n_mg_kg_1, Treatment, fertilization_colors, show_legend = TRUE) +
    labs(
      y = flux_units,
    ),
  nrow = 1
)

amoa_25_crop_row / amoa_25_addition_row / amoa_25_fert_row + plot_annotation(
  caption = "amoA 25 (gene copies g<sup>-1</sup>)",
  theme = theme(plot.caption = element_markdown(hjust = 0.5, size = 10))
)

ggsave(
  "figures/scatterplots/amoa25_vs_min.png",
  width = 4500,
  height = 3000,
  units = "px"
)

############### amoa_39
amoa_39_crop_row <- wrap_plots(
  plot_amoa_min_scatter(ave_039, net_min_rate_rel, Crop, crop_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
      title = "Net mineralization rate (rel.)"
    ),
  plot_amoa_min_scatter(ave_039, net_min_rate_abs, Crop, crop_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
      title = "Net mineralization rate (abs.)"
    ),
  plot_amoa_min_scatter(ave_039, no3n_mg_kg_1, Crop, crop_colors) +
    labs(
      y = flux_units,
      title = nitrate_label
    ),
  plot_amoa_min_scatter(ave_039, nh4n_mg_kg_1, Crop, crop_colors, show_legend = TRUE) +
    labs(
      y = flux_units,
      title = ammonia_label
    ),
  nrow = 1
)


amoa_39_addition_row <- wrap_plots(
  plot_amoa_min_scatter(ave_039, net_min_rate_rel, Addition, addition_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(ave_039, net_min_rate_abs, Addition, addition_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(ave_039, no3n_mg_kg_1, Addition, addition_colors) +
    labs(
      y = flux_units,
    ),
  plot_amoa_min_scatter(ave_039, nh4n_mg_kg_1, Addition, addition_colors, show_legend = TRUE) +
    labs(
      y = flux_units,
    ),
  nrow = 1
)

amoa_39_fert_row <- wrap_plots(
  plot_amoa_min_scatter(ave_039, net_min_rate_rel, Treatment, fertilization_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(ave_039, net_min_rate_abs, Treatment, fertilization_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(ave_039, no3n_mg_kg_1, Treatment, fertilization_colors) +
    labs(
      y = flux_units,
    ),
  plot_amoa_min_scatter(ave_039, nh4n_mg_kg_1, Treatment, fertilization_colors, show_legend = TRUE) +
    labs(
      y = flux_units,
    ),
  nrow = 1
)

amoa_39_crop_row / amoa_39_addition_row / amoa_39_fert_row + plot_annotation(
  caption = "amoA 39 (gene copies g<sup>-1</sup>)",
  theme = theme(plot.caption = element_markdown(hjust = 0.5, size = 10))
)

ggsave(
  "figures/scatterplots/amoa39_vs_min.png",
  width = 4500,
  height = 3000,
  units = "px"
)


############### F1R2
f1r2_crop_row <- wrap_plots(
  plot_amoa_min_scatter(f1r2_ave, net_min_rate_rel, Crop, crop_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
      title = "Net mineralization rate (rel.)"
    ),
  plot_amoa_min_scatter(f1r2_ave, net_min_rate_abs, Crop, crop_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
      title = "Net mineralization rate (abs.)"
    ),
  plot_amoa_min_scatter(f1r2_ave, no3n_mg_kg_1, Crop, crop_colors) +
    labs(
      y = flux_units,
      title = nitrate_label
    ),
  plot_amoa_min_scatter(f1r2_ave, nh4n_mg_kg_1, Crop, crop_colors, show_legend = TRUE) +
    labs(
      y = flux_units,
      title = ammonia_label
    ),
  nrow = 1
)


f1r2_addition_row <- wrap_plots(
  plot_amoa_min_scatter(f1r2_ave, net_min_rate_rel, Addition, addition_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(f1r2_ave, net_min_rate_abs, Addition, addition_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(f1r2_ave, no3n_mg_kg_1, Addition, addition_colors) +
    labs(
      y = flux_units,
    ),
  plot_amoa_min_scatter(f1r2_ave, nh4n_mg_kg_1, Addition, addition_colors, show_legend = TRUE) +
    labs(
      y = flux_units,
    ),
  nrow = 1
)

f1r2_fert_row <- wrap_plots(
  plot_amoa_min_scatter(f1r2_ave, net_min_rate_rel, Treatment, fertilization_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(f1r2_ave, net_min_rate_abs, Treatment, fertilization_colors) +
    labs(
      y = paste(flux_units, per_day_unit),
    ),
  plot_amoa_min_scatter(f1r2_ave, no3n_mg_kg_1, Treatment, fertilization_colors) +
    labs(
      y = flux_units,
    ),
  plot_amoa_min_scatter(f1r2_ave, nh4n_mg_kg_1, Treatment, fertilization_colors, show_legend = TRUE) +
    labs(
      y = flux_units,
    ),
  nrow = 1
)

f1r2_crop_row / f1r2_addition_row / f1r2_fert_row + plot_annotation(
  caption = "F1R2 (gene copies g<sup>-1</sup>)",
  theme = theme(plot.caption = element_markdown(hjust = 0.5, size = 10))
)

ggsave(
  "figures/scatterplots/F1R2_vs_min.png",
  width = 4500,
  height = 3000,
  units = "px"
)
