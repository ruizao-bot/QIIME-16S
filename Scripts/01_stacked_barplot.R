# ==============================================================================
# Family-Level Taxonomic Composition Analysis (Stacked Bar Plot)
# ==============================================================================
# Purpose: Visualize taxonomic composition at Family level across samples
# Input: final-table-with-ranks.tsv (QIIME2 feature table with taxonomy)
# Output: Family-level stacked bar plots (PDF and PNG)
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(RColorBrewer)


# Config variables (inlined)
BASE_DIR <- "/Users/jiayi/Desktop/metagenomic_pipeline/QIIME"
INPUT_FILE <- "Results/denoise_mode/final-table-with-ranks.tsv"
OUTPUT_DIR <- "Results/denoise_mode"
EXCLUDE_SAMPLES <- c("A1", "A2")
TOP_N_FAMILIES <- 75
TOP_N_DOMINANT <- 4

setwd(BASE_DIR)

# ------------------------------------------------------------------------------
# Data Loading and Preprocessing
# ------------------------------------------------------------------------------
cat("Loading data from:", INPUT_FILE, "\n")

# Read data with proper header handling
header <- scan(INPUT_FILE, what = "character", nlines = 1, sep = "\t", quiet = TRUE)
data <- read.table(INPUT_FILE, sep = "\t", skip = 1, header = FALSE, stringsAsFactors = FALSE)

# Fix column names
if (ncol(data) == length(header) + 1) {
  colnames(data) <- c("FeatureID", header)
} else {
  colnames(data) <- header
}

# Define metadata and sample columns
metadata_cols <- c("FeatureID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
sample_cols <- setdiff(colnames(data), metadata_cols)
sample_cols <- setdiff(sample_cols, EXCLUDE_SAMPLES)

cat("Loaded", nrow(data), "features across", length(sample_cols), "samples\n")

# Convert to long format and clean taxonomy
long_data <- data %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "SampleID", values_to = "Count") %>%
  filter(Count > 0) %>%
  mutate(
    # Clean Taxonomy: Remove prefixes (e.g., "f__", "g__") and trim whitespace
    Phylum = trimws(str_replace_all(Phylum, "[a-z]__", "")),
    Class = trimws(str_replace_all(Class, "[a-z]__", "")),
    Order = trimws(str_replace_all(Order, "[a-z]__", "")),
    Family = trimws(str_replace_all(Family, "[a-z]__", "")),
    Genus = trimws(str_replace_all(Genus, "[a-z]__", ""))
  )

# ------------------------------------------------------------------------------
# Family-Level Relative Abundance Calculation
# ------------------------------------------------------------------------------
cat("Calculating family-level relative abundances...\n")

# Calculate total abundance per sample
long_data_grouped <- long_data %>%
  group_by(SampleID) %>%
  mutate(TotalAbundance = sum(Count)) %>%
  ungroup()

# Calculate relative abundance for each Family
family_abundance <- long_data_grouped %>%
  group_by(SampleID, Family) %>%
  summarise(
    Count = sum(Count),
    TotalAbundance = first(TotalAbundance),
    .groups = "drop"
  ) %>%
  mutate(RelativeAbundance = Count / TotalAbundance * 100)

# Add grouping information (Surrounding: Bark, Soil, Wood)
family_abundance <- family_abundance %>%
  mutate(Surrounding = case_when(
    grepl("B$", SampleID) ~ "Bark",
    grepl("S$", SampleID) ~ "Soil",
    TRUE ~ "Wood"
  ))

# ------------------------------------------------------------------------------
# Identify Top Families and Prepare Plot Data
# ------------------------------------------------------------------------------
cat("Selecting top", TOP_N_FAMILIES, "families...\n")

# Identify top families by total relative abundance
top_families <- family_abundance %>%
  group_by(Family) %>%
  summarise(TotalRelAbundance = sum(RelativeAbundance)) %>%
  arrange(desc(TotalRelAbundance)) %>%
  slice(1:TOP_N_FAMILIES) %>%  
  pull(Family)

# Group low-abundance families as "Others"
family_plot_data <- family_abundance %>%
  mutate(Family_Group = ifelse(Family %in% top_families, Family, "Others")) %>%
  group_by(SampleID, Surrounding, Family_Group) %>%
  summarise(RelativeAbundance = sum(RelativeAbundance), .groups = "drop")

# ------------------------------------------------------------------------------
# Create Color Palette
# ------------------------------------------------------------------------------
n_families <- length(unique(family_plot_data$Family_Group))
family_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_families)
names(family_colors) <- sort(unique(family_plot_data$Family_Group))
family_colors["Others"] <- "grey80"  # Set Others to grey

# ------------------------------------------------------------------------------
# Generate Stacked Bar Plot
# ------------------------------------------------------------------------------
cat("Generating stacked bar plot...\n")

p_family <- ggplot(family_plot_data, aes(x = SampleID, y = RelativeAbundance, fill = Family_Group)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_grid(~ Surrounding, scales = "free_x", space = "free_x") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  scale_fill_manual(values = family_colors) +
  labs(x = "Sample", y = "Relative Abundance (%)", fill = "Family", 
       title = "Taxonomic Composition at Family Level") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# ------------------------------------------------------------------------------
# Save Outputs
# ------------------------------------------------------------------------------
# Display plot
print(p_family)

# Save plot
output_pdf <- file.path(OUTPUT_DIR, "stacked_barplot_family_level.pdf")
output_png <- file.path(OUTPUT_DIR, "stacked_barplot_family_level.png")

ggsave(output_pdf, p_family, width = 12, height = 7)
ggsave(output_png, p_family, width = 12, height = 7, dpi = 300)

cat("\n=== Family-level analysis complete ===\n")
cat("Stacked bar plot saved to:\n")
cat("  -", output_pdf, "\n")
cat("  -", output_png, "\n")

# Save processed data for downstream analyses
save(long_data, file = file.path(OUTPUT_DIR, "processed_long_data.RData"))
cat("  - Processed data saved to:", file.path(OUTPUT_DIR, "processed_long_data.RData"), "\n")
