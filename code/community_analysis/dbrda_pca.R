#!/usr/bin/env Rscript
# ---------------------------
# Dimensional analysis of amoA communities
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup/setup.R")

# Loads the create_phyloseq function
source("code/setup/create_phyloseq.R")

# Read in data
mineralization_and_qpcr_data <- read.csv(here("data/prepped_data", "mineralization_and_qpcr_data.csv")) %>%
    rename(
        Crop = crop,
        Day = day,
        Treatment = treatment,
        Addition = addition
    ) %>%
    mutate(Day = as.factor(Day))

# Separating out the community matrix (ys) and the metadata (xs)
xs <- mineralization_and_qpcr_data %>%
    select(Day, Crop, Treatment, Addition)
ys <- mineralization_and_qpcr_data %>%
    select(contains(c("ave", "log")))

# Helper functions

# Plots a single CAP plot
plot_cap_with_arrows <- function(qpcr_cap, voi, palette, variances) {
    cap_samples <- scores(qpcr_cap, tidy = TRUE) %>%
        filter(score == "sites") %>%
        bind_cols(xs)

    cap_arrows <- scores(qpcr_cap, tidy = TRUE) %>%
        filter(score == "factorbiplot")

    x_label <- paste0("CAP1 (", variances[2], "%)")
    y_label <- paste0("CAP2 (", variances[3], "%)")

    cap_samples %>%
        ggplot(aes(CAP1, CAP2, fill = {{ voi }}, shape = Day)) +
        geom_point(size = 4) +
        geom_segment(
            inherit.aes = FALSE,
            data = cap_arrows,
            aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
        ) +
        geom_label(
            inherit.aes = FALSE,
            data = cap_arrows,
            aes(x = CAP1 * 1.1, y = CAP2 * 1.1, label = label),
        ) +
        scale_shape_manual(values = c(21, 22, 23)) +
        scale_fill_manual(values = palette) +
        theme(
            axis.line = element_blank(),
            panel.border = element_rect(color = "black")
        ) +
        guides(fill = guide_legend(override.aes = list(shape = 21))) +
        labs(
            x = x_label,
            y = y_label
        )
}

# Plots all CAP plots for a particular configuration
plot_all_caps <- function(qpcr_cap, arrow_p_values, cap_variances) {
    day_plot <- plot_cap_with_arrows(qpcr_cap, Day, qpcr_day_colors, cap_variances) +
        labs(
            x = "",
        ) +
        annotate(
            geom = "richtext",
            x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
            label = arrow_p_values %>%
                filter(factor == "Day") %>%
                pull(label)
        )

    crop_plot <- plot_cap_with_arrows(qpcr_cap, Crop, crop_colors, cap_variances) +
        labs(
            x = "",
            y = ""
        ) + annotate(
            geom = "richtext",
            x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
            label = arrow_p_values %>%
                filter(factor == "Crop") %>%
                pull(label)
        )

    addition_plot <- plot_cap_with_arrows(qpcr_cap, Addition, addition_colors, cap_variances) +
        annotate(
            geom = "richtext",
            x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
            label = arrow_p_values %>%
                filter(factor == "Addition") %>%
                pull(label)
        )

    treatment_plot <- plot_cap_with_arrows(qpcr_cap, Treatment, fertilization_colors, cap_variances) +
        labs(
            y = ""
        ) +
        annotate(
            geom = "richtext",
            x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
            label = arrow_p_values %>%
                filter(factor == "Treatment") %>%
                pull(label)
        )

    return(
        (day_plot + crop_plot) / (addition_plot + treatment_plot)
    )
}


################## dbRDA

# Creating phyloseq objets for each configuration of all/MFP primers and log/avg qPCR values


all_dfs <- list(
    all_ave = ys %>%
        select(contains("ave")),
    mfp_ave = ys %>%
        select(starts_with("ave")),
    all_log = ys %>%
        select(contains("log")) %>%
        select(-log_sum),
    mfp_log = ys %>%
        select(starts_with("log")) %>%
        select(-log_sum)
)

dbrda_r2 <- data.frame(
    primers = character(),
    qpcr_values = character(),
    r2 = numeric(),
    adj_r2 = numeric()
)

dbrda_anova_results <- data.frame(
    primers = character(),
    qpcr_values = character(),
    term = character(),
    df = numeric(),
    SumOfSqs = numeric(),
    statistics = numeric(),
    p.value = numeric(),
    primers = character(),
    qpcr_values = character()
)

##### Big nasty for loop
for (i in seq_along(all_dfs)) {
    this_info <- names(all_dfs)[i] %>%
        str_split("_") %>%
        unlist()
    this_df <- all_dfs[[i]]

    print(paste0("dbRDA: Looking at ", this_info[1], " primers with ", this_info[2], " qPCR values"))
    print(this_df %>% colnames())

    this_cap <- capscale(
        this_df ~ Day + Treatment + Addition + Crop,
        xs,
        distance = "bray"
    )

    # This is the proportion of variance explained by each CAP axis
    this_cap_var_explained <- summary(this_cap)$cont %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        select(rowname, importance.CAP1, importance.CAP2) %>%
        filter(rowname == "Proportion Explained") %>%
        mutate(across(
            importance.CAP1:importance.CAP2,
            ~ round(., 3) * 100
        ))


    # Calculating the arrows...
    this_envfit <- envfit(
        this_cap,
        xs,
        perm = 9999
    )

    # ...and getting the pvalues
    arrows_p_values <- this_envfit$factors$pvals %>%
        as.data.frame() %>%
        rownames_to_column("factor") %>%
        rename(pvalue = 2) %>%
        mutate(
            pvalue = round(pvalue, 3),
            label = paste0(
                "p = ", pvalue, get_p_sig(pvalue)
            )
        )

    # This is the Rsquared of the model against the community structure.
    this_r2 <- RsquareAdj(this_cap)$r.squared %>% round(., 3)
    this_adj_r2 <- RsquareAdj(this_cap)$adj.r.squared %>% round(., 3)

    this_plot_title <- paste0(
        "CAP Analysis of ", this_info[1], " primers with ", this_info[2], " qPCR values"
    )
    this_plot_subtitle <- paste0(
        "R2 = ", this_r2, ", Adj R2 = ", this_adj_r2
    )


    plot_all_caps(this_cap, arrows_p_values, this_cap_var_explained) +
        plot_annotation(
            title = this_plot_title,
            subtitle = this_plot_subtitle,
            theme = theme(
                plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5)
            )
        )

    ggsave(
        filename = here("figures/capscale", paste0(this_info[1], "_", this_info[2], "_cap.png")),
        width = 4500,
        height = 3000,
        units = "px"
    )

    # Marginal ANOVA of terms in the model against community composition
    these_anova_results <- anova(this_cap, by = "margin") %>%
        broom::tidy() %>%
        filter(term != "Residual") %>%
        mutate(across(
            df:p.value,
            ~ round(., 3)
        )) %>%
        mutate(
            primers = this_info[1],
            qpcr_values = this_info[2]
        )


    # Updating R2 and ANOVA dataframes
    dbrda_r2 <- dbrda_r2 %>%
        add_row(
            primers = this_info[1],
            qpcr_values = this_info[2],
            r2 = this_r2,
            adj_r2 = this_adj_r2
        )

    dbrda_anova_results <- rbind(
        dbrda_anova_results,
        these_anova_results
    )
}

### Write out results
write.csv(
    dbrda_r2,
    here::here("results/community_analysis", "dbrda_r2.csv"),
    row.names = FALSE,
    quote = FALSE
)

write.csv(
    dbrda_anova_results,
    here::here("results/community_analysis", "dbrda_anova_results.csv"),
    row.names = FALSE,
    quote = FALSE
)

################## Scratch



plot_all_caps(qpcr_cap, tester, tester_variances) +
    plot_annotation()

RsquareAdj(qpcr_cap)
anova(qpcr_cap, by = "terms")

plot(qpcr_cap)
plot(qpcr_cap, display = "bp")
env_arrows <- envfit(qpcr_cap, xs, perm = 999)
envfit(qpcr_cap, xs)

tester <- env_arrows$factors$pvals %>%
    as.data.frame() %>%
    rownames_to_column("factor") %>%
    rename(pvalue = 2) %>%
    mutate(
        pvalue = round(pvalue, 3),
        label = paste0(
            "p = ", pvalue, get_p_sig(pvalue)
        )
    )

plot(env_arrows, col = "red")

cap_summary <- summary(qpcr_cap)
tester_variances <- cap_summary$cont %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    select(rowname, importance.CAP1, importance.CAP2) %>%
    filter(rowname == "Proportion Explained") %>%
    mutate(across(
        importance.CAP1:importance.CAP2,
        ~ round(., 3) * 100
    ))
