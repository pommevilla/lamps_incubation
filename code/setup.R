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

theme_set(
  theme_light() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5)
    )
)

nice_expansion <- expansion(add = 0, mult = c(0, 0.1))

data.priming <- read.csv("data/priming_amoA_deltaCt.csv", header = T) %>%
  rename(sample_id = X) %>%
  mutate(fert_level = as.factor(case_when(
    fert_level == 0 ~ "0N",
    fert_level == 336 ~ "336N",
    TRUE ~ "112N"
  )))

data.raw <- read.csv("data/priming_amoA_rawCt.csv", header = T) %>%
  rename(sample_id = X) %>%
  mutate(fert_level = as.factor(fert_level))

data.priming.long <- data.priming %>%
  pivot_longer(cols = starts_with("amoA"), names_to = "amoA", values_to = "deltaCT")

data.raw.long <- data.raw %>%
  pivot_longer(cols = starts_with("amoA"), names_to = "amoA", values_to = "CT")

data.priming.long$sample_id <- fct_reorder(data.priming.long$sample_id, parse_number(data.priming.long$sample_id))

df <- data.priming[, -1]
rownames(df) <- data.priming[, 1]

metadata <- df %>%
  select(fert_level:field_rep) %>%
  mutate(across(everything(), as.factor))


amoa_counts <- df %>%
  select(starts_with("amoA"))

n2o <- readr::read_csv("data/N2O.csv") %>%
  select(-1) %>%
  mutate(crop_within_block = interaction(crop, rep)) %>%
  mutate(field_sample = interaction(crop, fert, rep)) %>%
  mutate(addition = case_when(
    addition == "Cntrl" ~ "Control",
    addition == "+C" ~ "+ Carbon",
    TRUE ~ "+ Nitrogen"
  ))

co2 <- read.csv("data/priming_calc.csv")

tgas <- read.csv("data/tgas1.csv")

# Colors
addition_colors <- met.brewer("Benedictus", 3)
names(addition_colors) <- levels(data.priming$addition)

fertilization_colors <- met.brewer("VanGogh3", 3, type = "continuous")
names(fertilization_colors) <- c("0N", "112N", "336N")

crop_colors <- met.brewer("Java", 2, direction = -1)
names(crop_colors) <- levels(data.priming$crop)

# To make a consistent x-axis for days
day_breaks <- c(0, 4, 15, 30, 59, 86, 113, 144)
qpcr_day_breaks <- c(5, 32, 87)
