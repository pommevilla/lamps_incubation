#!/usr/bin/env Rscript
# ---------------------------
# Sets up the phyloseq object for the amoAs used in community analysis
# Assumes that code/setup/setup.R has been run and that
# mineralization_and_qpcr_data has been loaded
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
################## Load packages
library(vegan)
library(phyloseq)

################## Creating the phyloseq object
otus <- mineralization_and_qpcr_data %>%
    select(sample_name, contains(c("ave", "log"))) %>%
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

rm(otus, metadata, taxonomy)
