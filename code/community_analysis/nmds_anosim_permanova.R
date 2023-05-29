#!/usr/bin/env Rscript
# ---------------------------
# PERMANOVA, NMDS, ANOSIM analysis of communities
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup/setup.R")

# Loads the create_phyloseq function
source("code/setup/create_phyloseq.R")

# Read in data
mineralization_and_qpcr_data <- read.csv(here("data/prepped_data", "mineralization_and_qpcr_data.csv")) %>%
    mutate(Day = as.factor(Day))

# Separating out the metadata can be helpful for later
all_metadata <- mineralization_and_qpcr_data %>%
    select(sample_name, Day, Crop, Treatment, Addition) %>%
    mutate(Day = as.factor(Day)) %>%
    column_to_rownames("sample_name") %>%
    sample_data()

all_dfs <- list(
    all_ave = mineralization_and_qpcr_data %>%
        select(contains("ave")),
    mfp_ave = mineralization_and_qpcr_data %>%
        select(starts_with("ave")),
    all_log = mineralization_and_qpcr_data %>%
        select(contains("log")),
    mfp_log = mineralization_and_qpcr_data %>%
        select(starts_with("log"))
)


# Separating out the community matrix (ys) and the metadata (xs)
xs <- mineralization_and_qpcr_data %>%
    select(Day, Crop, Treatment, Addition)
ys <- mineralization_and_qpcr_data %>%
    select(contains(c("ave", "log")))

################## ANOSIM
# Run ANOSIM on log/ave mfp/all primer sets.
# Answers the question: are groups significantly different by treatment?
# R^2 measures how much of the variation is explained by the grouping factor,
# with higher being better. A negative R^2 implies that within-group differences
# are higher than without-group differences.
anosim_results <- data.frame(
    primers = character(),
    qpcr_values = character(),
    factor = character(),
    r = numeric(),
    p_value = numeric()
)


for (i in seq_along(all_dfs)) {
    this_info <- names(all_dfs)[i] %>%
        str_split("_") %>%
        unlist()
    print(paste0("ANOSIM: Looking at ", this_info[1], " primers with ", this_info[2], " qPCR values"))

    # Creates the dissimilarity matrix of distances
    ys_dis_m <- vegdist(
        all_dfs[[i]],
    )

    # Testing: are groups significantly different by treatment?
    for (this_factor_name in colnames(xs)) {
        this_anosim <- anosim(ys_dis_m, xs[[this_factor_name]])

        anosim_results <- rbind(anosim_results, data.frame(
            primers = this_info[1],
            qpcr_values = this_info[2],
            factor = this_factor_name,
            r = round(this_anosim$statistic, 3),
            p_value = round(this_anosim$signif, 3)
        ))
    }
}

# Write out ANOSIM results
write.csv(
    anosim_results,
    here::here("results", "community_analysis", "anosim_results.csv"),
    row.names = FALSE,
    quote = FALSE
)


################## NMDS
# Helper function to plot NMDS
plot_nmds <- function(nmds_scores, voi, palette, by_day = TRUE) {
    if (by_day) {
        p <- nmds_scores %>%
            ggplot(aes(NMDS1, NMDS2, shape = {{ voi }}, fill = {{ voi }})) +
            facet_wrap(~Day)
    } else {
        p <- nmds_scores %>%
            ggplot(aes(NMDS1, NMDS2, shape = Day, fill = {{ voi }})) +
            guides(fill = guide_legend(override.aes = list(shape = 21)))
    }

    p <- p +
        geom_point(size = 4) +
        scale_shape_manual(values = c(21, 22, 23)) +
        scale_fill_manual(values = palette) +
        theme(
            aspect.ratio = 1,
            strip.text = element_text(color = "black", size = 12),
            strip.background = element_rect(color = "NA", fill = "white"),
            panel.border = element_rect(color = "black", fill = "NA"),
            axis.line = element_blank()
        )

    return(p)
}

nmds_stresses <- data.frame(
    primers = character(),
    qpcr_values = character(),
    stress = numeric()
)

for (i in seq_along(all_dfs)) {
    this_info <- names(all_dfs)[i] %>%
        str_split("_") %>%
        unlist()
    print(paste0("NMDS: Looking at ", this_info[1], " primers with ", this_info[2], " qPCR values"))

    # Creates the dissimilarity matrix of distances
    ys_nmds <- metaMDS(
        all_dfs[[i]],
        distance = "bray"
    )

    this_stress <- round(ys_nmds$stress, 3)

    this_stress_caption <- paste0(
        "Stress: ",
        this_stress
    )

    ys_nmds_scores <- as.data.frame(scores(ys_nmds, display = "sites")) %>%
        merge(
            xs %>% mutate(Day = as.factor(Day)),
            by = "row.names"
        )

    nmds_stresses <- rbind(nmds_stresses, data.frame(
        primers = this_info[1],
        qpcr_values = this_info[2],
        stress = this_stress
    ))

    # NMDS plot of everything together
    this_plot_title <- paste0(this_info[1], " primers, ", this_info[2], " qPCR values")
    (
        plot_nmds(ys_nmds_scores, Day, c(crop_colors, "#FFCF00"), by_day = FALSE) +
            plot_nmds(ys_nmds_scores, Crop, crop_colors, by_day = FALSE) +
            labs(y = "")
    ) /
        (
            plot_nmds(ys_nmds_scores, Treatment, fertilization_colors, by_day = FALSE) +
                plot_nmds(ys_nmds_scores, Addition, addition_colors, by_day = FALSE) +
                labs(y = "")
        ) +
        plot_annotation(
            title = this_plot_title,
            subtitle = this_stress_caption,
            theme = theme(
                plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5),
            )
        )


    this_filename <- paste0("nmds_", names(all_dfs)[i], ".png")
    ggsave(
        here::here("figures/community_analysis/nmds_plots", this_filename),
        width = 4500,
        height = 3000,
        units = "px"
    )

    # NMDS plots by day
    (
        plot_nmds(ys_nmds_scores, Crop, crop_colors) +
            labs(x = "")
    ) / (
        plot_nmds(ys_nmsd_scores, Addition, addition_colors) +
            labs(x = "") +
            theme(strip.text = element_blank())
    ) / (
        plot_nmds(ys_nmds_scores, Treatment, fertilization_colors) +
            theme(strip.text = element_blank())
    ) +
        plot_annotation(
            title = this_plot_title,
            subtitle = this_stress_caption,
            theme = theme(
                plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5),
            )
        )

    this_filename <- paste0("nmds_", names(all_dfs)[i], "_by_day.png")
    ggsave(
        here::here("figures/community_analysis/nmds_plots", this_filename),
        width = 4500,
        height = 3000,
        units = "px"
    )
}

write.csv(
    nmds_stresses,
    here::here("results", "community_analysis", "nmds_stresses.csv"),
    row.names = FALSE,
    quote = FALSE
)

################## PERMANOVA
# Are there significant differences in community structure between treatments?
# Compared to dbRDA, this test is more concerned with detecting overall differences in
# community structure between groups, not quantifying those differences

# This is for all the days lumped together.
adonis_results <- data.frame(
    term = character(),
    df = numeric(),
    SumOfSqs = numeric(),
    R2 = numeric(),
    statistic = numeric(),
    p.value = numeric(),
    primers = character(),
    qpcr_values = character()
)

for (i in seq_along(all_dfs)) {
    this_df <- all_dfs[[i]]
    this_info <- names(all_dfs)[i] %>%
        str_split("_") %>%
        unlist()

    print(paste0("PERMANOVA: Looking at ", this_info[1], " primers with ", this_info[2], " qPCR values"))

    this_permanova <- adonis2(
        this_df ~ Treatment * Addition * Crop * Day,
        data = xs
    )

    this_permanova_tidied <- this_permanova %>%
        broom::tidy() %>%
        mutate(primers = this_info[1], qpcr_values = this_info[2]) %>%
        filter(term != "Total")

    adonis_results <- rbind(adonis_results, this_permanova_tidied)
}


write.csv(
    adonis_results %>%
        mutate(across(
            c(SumOfSqs, R2, statistic, p.value),
            ~ round(., 3)
        )),
    here::here("results", "community_analysis", "adonis_results.csv"),
    row.names = FALSE,
    quote = FALSE
)
