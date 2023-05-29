#!/usr/bin/env Rscript
# ---------------------------
# Plots NH4, NO3, N2O, and CO2 by factors over time
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup/setup.R")

######## Read in data
n2o_data <- read.csv("data/prepped_data/n2o_data.csv")
co2_data <- read.csv("data/prepped_data/co2_data.csv")
# nh4_no3_mineralization_data <- read.csv("data/prepped_data/mineralization_data.csv") %>%
#     rename(
#         Day = day,
#         Crop = crop,
#         Treatment = treatment,
#         Addition = addition
#     )

######## Helper functions
plot_lines_over_time <- function(input_df, foi, voi, palette,
                                 plot_title = "", x_label = "", y_label = "",
                                 show_legend = TRUE) {
    p <- input_df %>%
        group_by({{ foi }}, Day) %>%
        summarize(
            mean = mean({{ voi }}),
            sd = sd({{ voi }}),
            se = sd / sqrt(n())
        ) %>%
        ungroup() %>%
        ggplot(aes(Day, mean, fill = {{ foi }}, shape = {{ foi }})) +
        geom_line() +
        geom_errorbar(
            aes(ymin = mean - se, ymax = mean + se),
            width = 4,
        ) +
        geom_point(size = 4) +
        scale_fill_manual(values = palette) +
        scale_shape_manual(values = c(21, 22, 23)) +
        labs(
            x = x_label,
            y = y_label,
            title = plot_title,
        ) +
        theme(
            aspect.ratio = 1
        )

    if (!show_legend) {
        p <- p + theme(legend.position = "none")
    }

    return(p)
}

plot_factor_line <- function(foi, palette) {
    plot_lines_over_time(n2o_data, {{ foi }}, N2ON_flux_ug_g_d, palette, plot_title = n2o_label, show_legend = FALSE) +
        plot_lines_over_time(co2_data, {{ foi }}, CO2_flux_ug_g_d, palette, plot_title = co2_label, show_legend = FALSE) +
        plot_lines_over_time(n2o_data, {{ foi }}, cum_N2O_flux_ug_g, palette, plot_title = paste("Cumulative", n2o_label), show_legend = FALSE) +
        plot_lines_over_time(co2_data, {{ foi }}, cum_CO2_flux_ug_g, palette, plot_title = paste("Cumulative", co2_label)) +
        plot_layout(nrow = 1)
}


######## Plotting


crop_line <- plot_factor_line(Crop, crop_colors)
addition_line <- plot_factor_line(Addition, addition_colors)
treatment_line <- plot_factor_line(Treatment, fertilization_colors)

crop_line / addition_line / treatment_line

ggsave(
    here::here("figures/n_figures", "n_line_plots.png"),
    width = 4500,
    height = 3000,
    units = "px"
)
