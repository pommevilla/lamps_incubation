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
nh4_no3_min_data <- read.csv("data/prepped_data/mineralization_data.csv")


######## Helper functions
# This plots a single factor over time.
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
            aspect.ratio = 1,
            axis.title.y = element_markdown()
        )

    if (!show_legend) {
        p <- p + theme(legend.position = "none")
    }

    return(p)
}

# This plots a row of plots for a given factor and palette.
# There are separate functions for the cumulative version because...
# TODO: Add arguments for x and y labels
plot_factor_line <- function(foi, palette, show_titles = FALSE) {
    if (!show_titles) {
        plot_lines_over_time(
            n2o_data, {{ foi }}, N2ON_flux_ug_g_d, palette,
            plot_title = "", show_legend = FALSE, y_label = flux_units
        ) +
            plot_lines_over_time(
                co2_data, {{ foi }}, CO2_flux_ug_g_d, palette,
                plot_title = "", show_legend = FALSE
            ) +
            plot_lines_over_time(
                nh4_no3_min_data, {{ foi }}, no3n_mg_kg_1, palette,
                plot_title = "", show_legend = FALSE
            ) +
            plot_lines_over_time(
                nh4_no3_min_data, {{ foi }}, nh4n_mg_kg_1, palette,
                plot_title = ""
            ) +
            plot_layout(nrow = 1)
    } else {
        plot_lines_over_time(
            n2o_data, {{ foi }}, N2ON_flux_ug_g_d, palette,
            plot_title = n2o_label, show_legend = FALSE, y_label = flux_units
        ) +
            plot_lines_over_time(
                co2_data, {{ foi }}, CO2_flux_ug_g_d, palette,
                plot_title = co2_label, show_legend = FALSE
            ) +
            plot_lines_over_time(
                nh4_no3_min_data, {{ foi }}, no3n_mg_kg_1, palette,
                plot_title = nitrate_label, show_legend = FALSE
            ) +
            plot_lines_over_time(
                nh4_no3_min_data, {{ foi }}, nh4n_mg_kg_1, palette,
                plot_title = ammonia_label
            ) +
            plot_layout(nrow = 1)
    }
}

plot_cumulative_factor_line <- function(foi, palette, show_titles = FALSE) {
    if (!show_titles) {
        plot_lines_over_time(
            n2o_data, {{ foi }}, cum_N2O_flux_ug_g, palette,
            plot_title = "", show_legend = FALSE, y_label = flux_units
        ) +
            plot_lines_over_time(
                co2_data, {{ foi }}, cum_CO2_flux_ug_g, palette,
                plot_title = "", show_legend = FALSE
            ) +
            plot_lines_over_time(
                nh4_no3_min_data, {{ foi }}, cum_no3, palette,
                plot_title = "", show_legend = FALSE
            ) +
            plot_lines_over_time(
                nh4_no3_min_data, {{ foi }}, cum_nh4, palette,
                plot_title = ""
            ) +
            plot_layout(nrow = 1)
    } else {
        plot_lines_over_time(
            n2o_data, {{ foi }}, cum_N2O_flux_ug_g, palette,
            plot_title = n2o_label, show_legend = FALSE, y_label = flux_units
        ) +
            plot_lines_over_time(
                co2_data, {{ foi }}, cum_CO2_flux_ug_g, palette,
                plot_title = co2_label, show_legend = FALSE
            ) +
            plot_lines_over_time(
                nh4_no3_min_data, {{ foi }}, cum_no3, palette,
                plot_title = nitrate_label, show_legend = FALSE
            ) +
            plot_lines_over_time(
                nh4_no3_min_data, {{ foi }}, cum_nh4, palette,
                plot_title = ammonia_label
            ) +
            plot_layout(nrow = 1)
    }

    # plot_lines_over_time(n2o_data, {{ foi }}, cum_N2O_flux_ug_g, palette, plot_title = n2o_label, show_legend = FALSE) +
    #     plot_lines_over_time(co2_data, {{ foi }}, cum_CO2_flux_ug_g, palette, plot_title = co2_label, show_legend = FALSE) +
    #     plot_lines_over_time(nh4_no3_min_data, {{ foi }}, cum_no3, palette, plot_title = nitrate_label, show_legend = FALSE) +
    #     plot_lines_over_time(nh4_no3_min_data, {{ foi }}, cum_nh4, palette, plot_title = ammonia_label) +
    #     plot_layout(nrow = 1)
}

######## Plotting
# We'll start by plotting all the main factors separately. This gives us an overall picture how each factors affects
# response variables overall.
crop_line <- plot_factor_line(Crop, crop_colors, show_titles = TRUE)
addition_line <- plot_factor_line(Addition, addition_colors)
treatment_line <- plot_factor_line(Treatment, fertilization_colors)

crop_line / addition_line / treatment_line

ggsave(
    here::here("figures/n_figures/n_line_plots.png"),
    width = 4500,
    height = 3000,
    units = "px"
)

plot_cumulative_factor_line(Crop, crop_colors, show_titles = TRUE) / plot_cumulative_factor_line(Addition, addition_colors) / plot_cumulative_factor_line(Treatment, fertilization_colors)

ggsave(
    here::here("figures/n_figures/n_line_plots_cumulative.png"),
    width = 4500,
    height = 3000,
    units = "px"
)

######## Plotting Addition by other factors
# To investigate how different soils react to different management practices, we'll plot addition by the other factors.
n2o_data %>%
    group_by(Addition, Day, Crop) %>%
    summarize(
        mean = mean(N2ON_flux_ug_g_d),
        sd = sd(N2ON_flux_ug_g_d),
        se = sd / sqrt(n())
    ) %>%
    ungroup() %>%
    ggplot(aes(Day, mean, fill = Addition, shape = Addition)) +
    geom_line() +
    geom_errorbar(
        aes(ymin = mean - se, ymax = mean + se),
        width = 4,
    ) +
    geom_point(size = 4) +
    scale_fill_manual(values = addition_colors) +
    scale_shape_manual(values = c(21, 22, 23)) +
    labs(
        x = "x_label",
        y = "y_label",
        title = ,
    ) +
    theme(
        aspect.ratio = 1
    ) +
    facet_wrap(~Crop)

plot_lines_over_time(
    n2o_data %>% group_by(Crop),
    Addition,
    N2ON_flux_ug_g_d,
    addition_colors,
) +
    facet_grid(~Crop)

plot_addition_lines_over_time <- function(input_df, foi, voi,
                                          plot_title = "", x_label = "", y_label = "",
                                          show_legend = TRUE,
                                          show_strip_title = FALSE) {
    p <- input_df %>%
        group_by(Addition, {{ foi }}, Day) %>%
        summarize(
            mean = mean({{ voi }}),
            sd = sd({{ voi }}),
            se = sd / sqrt(n())
        ) %>%
        ungroup() %>%
        ggplot(aes(Day, mean, fill = Addition, shape = Addition)) +
        geom_line() +
        geom_errorbar(
            aes(ymin = mean - se, ymax = mean + se),
            width = 4,
        ) +
        geom_point(size = 3) +
        scale_fill_manual(values = addition_colors) +
        scale_shape_manual(values = c(21, 22, 23)) +
        labs(
            x = x_label,
            y = y_label,
            title = plot_title,
        ) +
        theme(
            aspect.ratio = 1,
            strip.background = element_rect(color = "black", fill = "#C6C6C6"),
            axis.title.y = element_markdown(size = 10),
            strip.text = element_blank()
        ) +
        facet_wrap(formula(paste("~", enquo(foi))))

    if (!show_legend) {
        p <- p + theme(legend.position = "none")
    }

    if (show_strip_title) {
        p <- p + theme(strip.text = element_markdown(color = "black", size = 14))
    }

    return(p)
}

# This and the next one are the same functions, but I CBA to refactor.
plot_addition_factor_grid <- function(foi) {
    plot_addition_lines_over_time(
        n2o_data, {{ foi }}, N2ON_flux_ug_g_d,
        y_label = make_flux_label(n2o_label),
        show_legend = FALSE, show_strip_title = TRUE
    ) +
        plot_addition_lines_over_time(
            co2_data, {{ foi }}, CO2_flux_ug_g_d,
            y_label = make_flux_label(co2_label),
            show_legend = FALSE
        ) +
        plot_addition_lines_over_time(
            nh4_no3_min_data, {{ foi }}, no3n_mg_kg_1,
            y_label = make_flux_label(nitrate_label),
            show_legend = FALSE
        ) +
        plot_addition_lines_over_time(
            nh4_no3_min_data, {{ foi }}, nh4n_mg_kg_1,
            y_label = make_flux_label(ammonia_label)
        ) +
        plot_layout(ncol = 1, guides = "collect")
}

plot_addition_factor_grid(Treatment)



ggsave(
    here::here("figures/n_figures/n_line_addition_nrate.png"),
    width = 2700,
    height = 3500,
    units = "px"
)

plot_addition_factor_grid(Crop)

ggsave(
    here::here("figures/n_figures/n_line_addition_crop.png"),
    width = 2000,
    height = 3500,
    units = "px"
)

plot_addition_cum_factor_grid <- function(foi) {
    plot_addition_lines_over_time(
        n2o_data, {{ foi }}, cum_N2O_flux_ug_g,
        y_label = make_flux_label(n2o_label, cumulative = TRUE),
        show_legend = FALSE, show_strip_title = TRUE
    ) +
        plot_addition_lines_over_time(
            co2_data, {{ foi }}, cum_CO2_flux_ug_g,
            y_label = make_flux_label(co2_label, cumulative = TRUE),
            show_legend = FALSE
        ) +
        plot_addition_lines_over_time(
            nh4_no3_min_data, {{ foi }}, cum_no3,
            y_label = make_flux_label(nitrate_label, cumulative = TRUE),
            show_legend = FALSE
        ) +
        plot_addition_lines_over_time(
            nh4_no3_min_data, {{ foi }}, cum_nh4,
            y_label = make_flux_label(ammonia_label, cumulative = TRUE)
        ) +
        plot_layout(ncol = 1, guides = "collect")
}

plot_addition_cum_factor_grid(Treatment)

ggsave(
    here::here("figures/n_figures/n_line_addition_nrate_cumulative.png"),
    width = 2700,
    height = 4000,
    units = "px"
)

plot_addition_cum_factor_grid(Crop)

ggsave(
    here::here("figures/n_figures/n_line_addition_crop_cumulative.png"),
    width = 2000,
    height = 4000,
    units = "px"
)

plot_addition_factor_line(Treatment)
n2o_addition_crop_line <- plot_addition_lines_over_time(
    n2o_data,
    Crop,
    cum_N2O_flux_ug_g,
    y_label = n2o_label
)
