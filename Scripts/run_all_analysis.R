#!/usr/bin/env Rscript

# Main script to run all analysis steps using modular functions
# Usage: Rscript Scripts/run_all_analysis.R

# Load required packages
required_packages <- c("dplyr", "tidyverse", "vegan", "ggplot2", "RColorBrewer", "pheatmap", "stringr", "tidyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}
# change to your own project root directory
setwd("/Users/jiayi/Desktop/metagenomic_pipeline/QIIME")

# Source function files
source("Scripts/table_combine_functions.R")
source("Scripts/taxonomic_analysis_functions.R")
source("Scripts/mob_analysis_functions.R")

# Step 1: Table combine
cat("\n[Step 1] Merging feature table and taxonomy...\n")
combine_result <- combine_feature_taxonomy(
  ft_path = "Results/denoise_mode/exported-table/feature-table.tsv",
  taxonomy_path = "Results/denoise_mode/exported-taxonomy/taxonomy.tsv",
  output_dir = "Results/denoise_mode"
)
feature_table <- combine_result$feature_table
taxonomy <- combine_result$taxonomy
merged_ranks <- combine_result$merged_ranks
cat("  Output: taxonomy-abundance-table.tsv, final-table-with-ranks.tsv\n")

# Step 2: Diversity analysis
cat("\n[Step 2] Diversity analysis...\n")
metadata <- read.table("Data/metadata/metadata.tsv", header = TRUE, sep = "\t", row.names = 1, comment.char = "")
metadata <- metadata[colnames(feature_table), ]
tax_table <- t(sapply(taxonomy$Taxon, parse_taxonomy2))
rownames(tax_table) <- rownames(taxonomy)

# Filter out only Chloroplast
cat("Filtering out Chloroplast...\n")
remove_idx <- (
  grepl("Chloroplast", tax_table[, "Order"], ignore.case = TRUE) |
  grepl("Chloroplast", tax_table[, "Family"], ignore.case = TRUE) |
  grepl("Chloroplast", tax_table[, "Genus"], ignore.case = TRUE)
)
feature_table_filtered <- feature_table[!remove_idx, ]
tax_table_filtered <- tax_table[!remove_idx, ]

# Subset by domain
kingdoms <- tax_table_filtered[, "Kingdom"]
bacteria_idx <- grepl("Bacteria", kingdoms, ignore.case = TRUE)
feature_table_bacteria <- feature_table_filtered[bacteria_idx, ]
tax_table_bacteria <- tax_table_filtered[bacteria_idx, ]
bacteria_archaea_idx <- grepl("Bacteria|Archaea", kingdoms, ignore.case = TRUE)
feature_table_bacteria_archaea <- feature_table_filtered[bacteria_archaea_idx, ]
tax_table_bacteria_archaea <- tax_table_filtered[bacteria_archaea_idx, ]

# Create output directories
dir.create("Results/denoise_mode/bacteria_only", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/denoise_mode/bacteria_archaea", showWarnings = FALSE, recursive = TRUE)

# Run analysis for both datasets
run_analysis(feature_table_bacteria, tax_table_bacteria, metadata, "Results/denoise_mode/bacteria_only/bacteria_only", "Bacteria Only")
run_analysis(feature_table_bacteria_archaea, tax_table_bacteria_archaea, metadata, "Results/denoise_mode/bacteria_archaea/bacteria_archaea", "Bacteria + Archaea")

# Step 3: MOB analysis
cat("\n[Step 3] MOB analysis...\n")
input_file <- "Results/denoise_mode/final-table-with-ranks.tsv"
header <- scan(input_file, what = "character", nlines = 1, sep = "\t", quiet = TRUE)
data <- read.table(input_file, sep = "\t", skip = 1, header = FALSE, stringsAsFactors = FALSE)
if (ncol(data) == length(header) + 1) {
  colnames(data) <- c("FeatureID", header)
} else {
  colnames(data) <- header
}
metadata_cols <- c("FeatureID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
sample_cols <- setdiff(colnames(data), metadata_cols)
sample_cols <- setdiff(sample_cols, c("A1", "A2"))
long_data <- data %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "SampleID", values_to = "Count") %>%
  filter(Count > 0) %>%
  mutate(
    Domain = trimws(str_replace_all(Domain, "[a-z]__", "")),
    Phylum = trimws(str_replace_all(Phylum, "[a-z]__", "")),
    Class = trimws(str_replace_all(Class, "[a-z]__", "")),
    Order = trimws(str_replace_all(Order, "[a-z]__", "")),
    Family = trimws(str_replace_all(Family, "[a-z]__", "")),
    Genus = trimws(str_replace_all(Genus, "[a-z]__", ""))
  )
long_data_bacteria <- long_data %>% filter(grepl("Bacteria", Domain, ignore.case = TRUE))
long_data_bacteria_archaea <- long_data %>% filter(grepl("Bacteria|Archaea", Domain, ignore.case = TRUE))
run_mob_analysis(long_data_bacteria, "Results/denoise_mode/bacteria_only/bacteria_only", "Bacteria Only")
run_mob_analysis(long_data_bacteria_archaea, "Results/denoise_mode/bacteria_archaea/bacteria_archaea", "Bacteria + Archaea")
cat("  MOB analysis complete.\n")

cat("\nAll analysis steps finished.\n")
