# ==============================================================================
# MOB (Methanotrophic Bacteria) Analysis and Visualization Script
# ==============================================================================
# Purpose: 
#   1. Visualize taxonomic composition at Family level across different surroundings
#   2. Analyze MOB distribution and relative abundance
#   3. Identify dominant/rare MOB groups and key species by site
# 
# Input: final-table-with-ranks.tsv (QIIME2 feature table with taxonomy)
# Output: 
#   - Stacked bar plots (Family level and MOB genus level)
#   - MOB relative abundance tables
#   - MOB classification and key species tables
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP: Load Libraries and Set Paths
# ------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(RColorBrewer)

setwd("/Users/jiayi/Desktop/metagenomic_pipeline/QIIME")

# Set input and output file paths
input_file <- "Results/denoise_mode/final-table-with-ranks.tsv"
output_plot <- "Results/denoise_mode/stacked_barplot_MOB.pdf"

# ------------------------------------------------------------------------------
# 2. DATA LOADING AND PREPROCESSING
# ------------------------------------------------------------------------------
# Read data with proper header handling
header <- scan(input_file, what = "character", nlines = 1, sep = "\t", quiet = TRUE)
data <- read.table(input_file, sep = "\t", skip = 1, header = FALSE, stringsAsFactors = FALSE)

# Fix column names
if (ncol(data) == length(header) + 1) {
  colnames(data) <- c("FeatureID", header)
} else {
  colnames(data) <- header
}

# Define metadata and sample columns
metadata_cols <- c("FeatureID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
sample_cols <- setdiff(colnames(data), metadata_cols)
sample_cols <- setdiff(sample_cols, c("A1", "A2")) # Filter out specific samples if needed

# Convert to long format and clean taxonomy
long_data <- data %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "SampleID", values_to = "Count") %>%
  filter(Count > 0) %>% # Remove zero counts for efficiency
  mutate(
    # Clean Taxonomy: Remove prefixes (e.g., "f__", "g__") and trim whitespace
    Phylum = trimws(str_replace_all(Phylum, "[a-z]__", "")),
    Class = trimws(str_replace_all(Class, "[a-z]__", "")),
    Order = trimws(str_replace_all(Order, "[a-z]__", "")),
    Family = trimws(str_replace_all(Family, "[a-z]__", "")),
    Genus = trimws(str_replace_all(Genus, "[a-z]__", ""))
  )


# ------------------------------------------------------------------------------
# 3. FAMILY-LEVEL TAXONOMIC ANALYSIS
# ------------------------------------------------------------------------------
# Calculate total abundance per sample
long_data_grouped <- long_data %>%
  group_by(SampleID) %>%
  mutate(TotalAbundance = sum(Count)) %>%
  ungroup()

# Calculate relative abundance for each Family in each sample
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

# Identify top 75 families across all samples by total relative abundance
top_families <- family_abundance %>%
  group_by(Family) %>%
  summarise(TotalRelAbundance = sum(RelativeAbundance)) %>%
  arrange(desc(TotalRelAbundance)) %>%
  slice(1:75) %>%  
  pull(Family)

# Group low-abundance families as "Others"
family_plot_data <- family_abundance %>%
  mutate(Family_Group = ifelse(Family %in% top_families, Family, "Others")) %>%
  group_by(SampleID, Surrounding, Family_Group) %>%
  summarise(RelativeAbundance = sum(RelativeAbundance), .groups = "drop")

# Create color palette for families (assign grey to "Others")
n_families <- length(unique(family_plot_data$Family_Group))
family_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_families)
names(family_colors) <- sort(unique(family_plot_data$Family_Group))
family_colors["Others"] <- "grey80"

# Generate Family-level stacked bar plot
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

# Display and save Family-level plot
print(p_family)
ggsave("Results/denoise_mode/stacked_barplot_family_level.pdf", p_family, width = 12, height = 7)
ggsave("Results/denoise_mode/stacked_barplot_family_level.png", p_family, width = 12, height = 7, dpi = 300)

print("Family-level stacked bar plot saved to:")
print("  - denoise_mode/stacked_barplot_family_level.pdf")
print("  - denoise_mode/stacked_barplot_family_level.png")

# ------------------------------------------------------------------------------
# 4. MOB (METHANOTROPHIC BACTERIA) IDENTIFICATION AND FILTERING
# ------------------------------------------------------------------------------
# Define target MOB genera for analysis
target_genera <- c(
  "Methylacidiphilum", 
  "Crenothrix", 
  "Methylibium", 
  "Methylobacillus", 
  "Methylobacterium", 
  "Methylocella", 
  "Methylomonas", 
  "Methylonatrum", 
  "Methylophaga", 
  "Methylopila", 
  "Methylosinus",
  "Methylovirgula"
)

# Create regex pattern for genus matching
target_pattern <- paste(target_genera, collapse = "|")

# Filter data for MOB taxa
filtered_data <- long_data %>%
  mutate(
    # Extract matching genus name
    MatchedGenus = str_extract(Genus, regex(target_pattern, ignore_case = TRUE)),
    
    # Handle special cases: Unclassified MOB families
    IsUnclassifiedMethylophilaceae = grepl("Methylophilaceae", Family, ignore.case = TRUE) & 
                                     (Genus == "" | Genus == "uncultured" | Genus == "Unassigned"),
    IsUnclassifiedMethylacidiphilaceae = grepl("Methylacidiphilaceae", Family, ignore.case = TRUE) & 
                                     (Genus == "" | Genus == "uncultured" | Genus == "Unassigned")
  ) %>%
  filter(!is.na(MatchedGenus) | IsUnclassifiedMethylophilaceae | IsUnclassifiedMethylacidiphilaceae) %>%
  mutate(
    Taxon = case_when(
      !is.na(MatchedGenus) ~ str_to_title(MatchedGenus),
      IsUnclassifiedMethylophilaceae ~ "Unclassified Methylophilaceae",
      IsUnclassifiedMethylacidiphilaceae ~ "Unclassified Methylacidiphilaceae"
    ),
    # Add "Candidatus" prefix for Methylacidiphilum
    Taxon = ifelse(Taxon == "Methylacidiphilum", "Candidatus Methylacidiphilum", Taxon)
  )

# Validate that MOB taxa were found
if (nrow(filtered_data) == 0) {
  print("DEBUG: Top 20 Families found in your data (to help debug):")
  print(head(sort(table(long_data$Family), decreasing=TRUE), 20))
  stop("No target MOB species found! Check the debug output above to see what Families are actually in your data.")
}

# ------------------------------------------------------------------------------
# 5. MOB RELATIVE ABUNDANCE CALCULATION AND VISUALIZATION
# ------------------------------------------------------------------------------
# Calculate relative abundance WITHIN MOB (normalized to MOB total per sample)
plot_data <- filtered_data %>%
  group_by(SampleID, Taxon) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  group_by(SampleID) %>%
  mutate(RelativeAbundance = Count / sum(Count) * 100) %>%
  ungroup()

# Add grouping information (Surrounding: Bark, Soil, Wood)
plot_data <- plot_data %>%
  mutate(Group = case_when(
    grepl("B$", SampleID) ~ "Bark",
    grepl("S$", SampleID) ~ "Soil",
    TRUE ~ "Wood"
  ))

# Define custom colors for MOB taxa
mob_colors <- c(
  "Candidatus Methylacidiphilum" = "#CD5C5C", # IndianRed
  "Crenothrix" = "black",
  "Methylibium" = "#FF8C00", # DarkOrange
  "Methylobacillus" = "#98FB98", # PaleGreen
  "Methylobacterium" = "#7FFF00", # Chartreuse
  "Methylocella" = "#8B4513", # SaddleBrown
  "Methylomonas" = "#E6E6FA", # Lavender
  "Methylonatrum" = "#483D8B", # DarkSlateBlue
  "Methylophaga" = "#800080", # Purple
  "Methylopila" = "#C0C0C0", # Silver/Grey
  "Methylosinus" = "#DA70D6", # Orchid
  "Unclassified Methylophilaceae" = "#AFEEEE"  # PaleTurquoise
)

# Set factor levels for consistent legend ordering
plot_data$Taxon <- factor(plot_data$Taxon, levels = sort(unique(plot_data$Taxon)))

# Create MOB stacked bar plot
p <- ggplot(plot_data, aes(x = SampleID, y = RelativeAbundance, fill = Taxon)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = mob_colors) +
  labs(x = "Tree", y = "Relative abundance (%)", fill = "MOB - 16S") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(size = 12),
    legend.position = "right",
    legend.text = element_text(size = 10, face = "italic")
  )

# Display and save MOB plot
print(p)
ggsave(output_plot, p, width = 10, height = 6)
ggsave(gsub(".pdf", ".png", output_plot), p, width = 10, height = 6)

print(paste("MOB plotting complete. Saved to:", output_plot))

# ------------------------------------------------------------------------------
# 6. MOB TOTAL RELATIVE ABUNDANCE PER SAMPLE
# ------------------------------------------------------------------------------
# Extract species identifier from SampleID (first letter)
plot_data <- plot_data %>%
  mutate(Species = str_extract(SampleID, "^[A-Z]"))

# Calculate total MOB abundance per sample
mob_totals <- filtered_data %>%
  group_by(SampleID) %>%
  summarise(MOB_Total = sum(Count)) %>%
  ungroup()

# Get total abundance for all taxa per sample
all_totals <- long_data %>%
  group_by(SampleID) %>%
  summarise(Total_Abundance = sum(Count)) %>%
  ungroup()

# Calculate MOB relative abundance as percentage of total community
mob_rel_abundance <- mob_totals %>%
  left_join(all_totals, by = "SampleID") %>%
  mutate(
    MOB_RelAbundance = MOB_Total / Total_Abundance * 100,
    Species = str_extract(SampleID, "^[A-Z]"),
    Surrounding = case_when(
      grepl("B$", SampleID) ~ "Bark",
      grepl("S$", SampleID) ~ "Soil",
      TRUE ~ "Wood"
    )
  )

# Format per-sample results
mob_per_sample <- mob_rel_abundance %>%
  select(SampleID, Species, Surrounding, MOB_RelAbundance) %>%
  rename(`MOB_RelAbundance (%)` = MOB_RelAbundance)

# Display and save per-sample MOB abundance
print("=== MOB Relative Abundance Per Sample ===")
print(mob_per_sample)
write.csv(mob_per_sample, "Results/denoise_mode/MOB_relative_abundance_per_sample.csv", row.names = FALSE)

print("MOB relative abundance per sample saved to:")
print("  - denoise_mode/MOB_relative_abundance_per_sample.csv")

# ------------------------------------------------------------------------------
# 7. DOMINANT vs RARE MOB GROUP CLASSIFICATION
# ------------------------------------------------------------------------------
# Calculate total relative abundance for each MOB genus across all samples
mob_genus_totals <- plot_data %>%
  group_by(Taxon) %>%
  summarise(Total_RelAbundance = sum(RelativeAbundance)) %>%
  arrange(desc(Total_RelAbundance)) %>%
  ungroup()

# Define top 4 as dominant/key groups
key_groups <- mob_genus_totals %>%
  slice(1:4) %>%
  pull(Taxon)

# Define remaining as rare groups
rare_groups <- mob_genus_totals %>%
  slice(5:n()) %>%
  pull(Taxon)

# Create classification table
mob_classification <- mob_genus_totals %>%
  mutate(
    Group_Type = case_when(
      Taxon %in% key_groups ~ "Dominant/Key Group",
      TRUE ~ "Rare Group"
    ),
    Rank = row_number()
  )

# Display and save classification
print("=== MOB Group Classification ===")
print(mob_classification)

print("\nDominant/Key MOB Groups (Top 4):")
print(key_groups)

print("\nRare MOB Groups:")
print(rare_groups)

write.csv(mob_classification, "Results/denoise_mode/MOB_group_classification.csv", row.names = FALSE)

print("\nMOB group classification saved to:")
print("  - denoise_mode/MOB_group_classification.csv")

# ------------------------------------------------------------------------------
# 8. KEY SPECIES IDENTIFICATION BY SITE (SURROUNDING)
# ------------------------------------------------------------------------------
# Calculate total relative abundance for each genus in each surrounding
mob_by_site <- plot_data %>%
  group_by(Group, Taxon) %>%
  summarise(Total_RelAbundance = sum(RelativeAbundance), .groups = "drop") %>%
  arrange(Group, desc(Total_RelAbundance))

# Get top 4 key species for each site
key_species_by_site <- mob_by_site %>%
  group_by(Group) %>%
  slice(1:3) %>%
  mutate(Rank = row_number()) %>%
  ungroup()

# Display and save site-specific key species
print("\n=== Key MOB Species by Site (Top 4 per Surrounding) ===")
print(key_species_by_site)

write.csv(key_species_by_site, "Results/denoise_mode/MOB_key_species_by_site.csv", row.names = FALSE)

print("\nKey species by site saved to:")
print("  - denoise_mode/MOB_key_species_by_site.csv")

# ==============================================================================
# ANALYSIS COMPLETE
# ==============================================================================
print("\n=== All analyses completed successfully ===")
print("Output files generated:")
print("  1. denoise_mode/stacked_barplot_family_level.pdf")
print("  2. denoise_mode/stacked_barplot_family_level.png")
print("  3. denoise_mode/stacked_barplot_MOB.pdf")
print("  4. denoise_mode/stacked_barplot_MOB.png")
print("  5. denoise_mode/MOB_relative_abundance_per_sample.csv")
print("  6. denoise_mode/MOB_group_classification.csv")
print("  7. denoise_mode/MOB_key_species_by_site.csv")
