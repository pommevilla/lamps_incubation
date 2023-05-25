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

################## Helper function to create phyloseq objects
# qpcr_df is dataframe with all of the rows and columns from the qPCR data
# logged is for whether or the analysis should be for the logged or average qPCR counts
# include_f1r2 is for whether or not to include the existing primers
create_phyloseq <- function(qpcr_df, logged = FALSE, include_f1r2 = TRUE) {
    if (logged) {
        otu_df <- qpcr_df %>%
            select(sample_name, contains("log"))
    } else {
        otu_df <- qpcr_df %>%
            select(sample_name, contains("ave"))
    }

    if (!include_f1r2) {
        otu_df <- otu_df %>%
            select(-contains("f1r2"))
    }

    otus <- otu_df %>%
        select(-any_of("log_sum")) %>%
        column_to_rownames("sample_name") %>%
        t() %>%
        otu_table(., taxa_are_rows = TRUE)

    metadata <- qpcr_df %>%
        select(sample_name, Day, Crop, Treatment, Addition) %>%
        mutate(Day = as.factor(Day)) %>%
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

    return(physeq)
}

################## Metadata for qPCR data
