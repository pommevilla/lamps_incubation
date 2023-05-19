#!/usr/bin/env Rscript
# ---------------------------
# Dimensional analysis of amoA communities
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Setup
source("code/setup.R")

# ordinate
library(vegan)
library(phyloseq)

# Read in data
mineralization_and_qpcr_data <- read.csv(here("data/prepped_data", "mineralization_and_qpcr_data.csv")) %>%
  rename(
    Day = day,
    Crop = crop,
    Treatment = treatment,
    Addition = addition
  )

################## Creating the phyloseq object
otus <- mineralization_and_qpcr_data %>%
  select(sample_name, contains(c("ave"))) %>%
  column_to_rownames("sample_name") %>%
  t() %>%
  otu_table(., taxa_are_rows = TRUE)

metadata <- mineralization_and_qpcr_data %>%
  select(sample_name, Day, Crop, Treatment, Addition) %>%
  column_to_rownames("sample_name") %>%
  sample_data()

taxonomy <- matrix(
  sample(letters, nrow(otus) * 7, replace = TRUE),
  nrow = nrow(otus),
  ncol = 7,
  dimnames = list(
    rownames(otus),
    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  )
) %>%
  tax_table()

physeq <- phyloseq(
  otus,
  taxonomy,
  metadata
)

################## Meangingless bar plots
plot_bar(
  physeq,
  fill = "Family"
)

################## distance-based redundancy analysis (dbRDA)
dba <- ordinate(
  physeq = physeq,
  method = "CAP",
  formula = ~ Day + Crop + Treatment + Addition,
  distance = "bray"
)

anova(dba, by = "margin")
RsquareAdj(dba)
