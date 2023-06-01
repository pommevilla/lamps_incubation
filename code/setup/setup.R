#!/usr/bin/env Rscript
# ---------------------------
# Various setup stuff for the project
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
library(tidyverse)
library(lubridate)
library(readxl)
library(janitor)
library(MetBrewer)
library(patchwork)
library(here)
library(ggtext)


####### Plotting setup
theme_set(
  theme_light() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      text = element_text(family = "Times New Roman"),
      plot.title = element_markdown(hjust = 0.5),
      plot.subtitle = element_markdown(hjust = 0.5)
    )
)

nice_expansion <- expansion(add = 0, mult = c(0, 0.1))
nice_dodge <- position_dodge(width = 0.1)

# Factor colors
addition_colors <- met.brewer("Benedictus", 3)
names(addition_colors) <- c("Control", "Carbon", "Nitrogen")

fertilization_colors <- met.brewer("VanGogh3", 3, type = "continuous")
names(fertilization_colors) <- c("0N", "112N", "336N")

crop_colors <- met.brewer("Java", 2, direction = -1)
names(crop_colors) <- c("Corn", "Miscanthus")

qpcr_day_colors <- met.brewer("Greek", 3)
names(qpcr_day_colors) <- c("4", "30", "86")

# To make a consistent x-axis for days
day_breaks <- c(0, 4, 15, 30, 59, 86, 113, 144)
qpcr_day_breaks <- c(5, 32, 87)

# Units; requires ggtext::element_markdown
flux_units <- "mg N kg<sup>-1</sup>"
per_day_unit <- "day<sup>-1</sup>"
gcn_unit <- "gene copies g<sup>-1</sup>"

# Labels; requires ggtext::element_markdown
ammonia_label <- "NH<sub>4</sub><sup>+</sup>-N"
nitrate_label <- "NO<sub>3</sub><sup>-</sup>-N"
co2_label <- "CO<sub>2</sub>"
n2o_label <- "N<sub>2</sub>O"


# data.priming <- read.csv("data/priming_amoA_deltaCt.csv", header = T) %>%
#   rename(sample_id = X) %>%
#   mutate(fert_level = as.factor(case_when(
#     fert_level == 0 ~ "0N",
#     fert_level == 336 ~ "336N",
#     TRUE ~ "112N"
#   )))

# data.raw <- read.csv("data/priming_amoA_rawCt.csv", header = T) %>%
#   rename(sample_id = X) %>%
#   mutate(fert_level = as.factor(fert_level))

# data.priming.long <- data.priming %>%
#   pivot_longer(cols = starts_with("amoA"), names_to = "amoA", values_to = "deltaCT")

# data.raw.long <- data.raw %>%
#   pivot_longer(cols = starts_with("amoA"), names_to = "amoA", values_to = "CT")

# data.priming.long$sample_id <- fct_reorder(data.priming.long$sample_id, parse_number(data.priming.long$sample_id))

# df <- data.priming[, -1]
# rownames(df) <- data.priming[, 1]

# metadata <- df %>%
#   select(fert_level:field_rep) %>%
#   mutate(across(everything(), as.factor))


# amoa_counts <- df %>%
#   select(starts_with("amoA"))

# n2o <- readr::read_csv("data/N2O.csv") %>%
#   select(-1) %>%
#   mutate(crop_within_block = interaction(crop, rep)) %>%
#   mutate(field_sample = interaction(crop, fert, rep)) %>%
#   mutate(addition = case_when(
#     addition == "Cntrl" ~ "Control",
#     addition == "+C" ~ "+ Carbon",
#     TRUE ~ "+ Nitrogen"
#   ))

# co2 <- read.csv("data/priming_calc.csv")

# tgas <- read.csv("data/tgas1.csv")

# variable names
mineralization_variables <- c(
  "net_min_rate_rel", "net_min_rate_abs",
  "net_nitr_rate_rel", "net_nitr_rate_abs",
  "no3n_mg_kg_1", "nh4n_mg_kg_1"
)

qpcr_variables <- c(
  "log_012", "ave_012",
  "log_025", "ave_025",
  "log_039", "ave_039",
  "F1R2_log", "F1R2_ave",
  "ave_sum", "log_sum",
  "norB.001", "log_norB.001",
  "norB.006", "log_norB.006",
  "cnorB", "log_cnorB",
  "norb_sum", "log_norb_sum"
)

# These are variable lists used for plotting in qpcr_plots.R.
# They're in the order desired for plotting in the plot functions
ave_amoa_qpcr_variables <- c(
  "ave_012", "ave_025", "ave_039", "F1R2_ave", "ave_sum"
)

log_amoa_qpcr_variables <- c(
  "log_012", "log_025", "log_039", "F1R2_log", "log_sum"
)

ave_norb_qpcr_variables <- c(
  "norB.001", "norB.006", "norb_sum", "cnorB"
)

log_norb_qpcr_variables <- c(
  "log_norB.001", "log_norB.006", "log_norb_sum", "log_cnorB"
)


# Helper functions
get_p_sig <- function(p_value) {
  sig_marker <- case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )

  return(sig_marker)
}
