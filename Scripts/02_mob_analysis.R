# ==============================================================================
# MOB (Methanotrophic Bacteria) Comprehensive Analysis
# ==============================================================================
# Purpose: 
#   1. Identify and filter MOB taxa
#   2. Visualize MOB distribution with stacked bar plots
#   3. Calculate MOB total relative abundance per sample
#   4. Classify dominant vs rare MOB groups
#   5. Identify key species by site (Bark/Soil/Wood)
#
# Input: processed_long_data.RData (from 01_stacked_barplot.R)
# Output: 
#   - MOB stacked bar plots (PDF/PNG)
#   - MOB relative abundance tables (CSV)
#   - MOB classification results (CSV)
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
OUTPUT_DIR <- "Results/denoise_mode"
TARGET_MOB_GENERA <- c(
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
MOB_COLORS <- c(
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
  "Unclassified Methylophilaceae" = "#AFEEEE",  # PaleTurquoise
  "Unclassified Methylacidiphilaceae" = "#FFB6C1" # LightPink
)
TOP_N_DOMINANT <- 4

setwd(BASE_DIR)

# ------------------------------------------------------------------------------
# Load Preprocessed Data
# ------------------------------------------------------------------------------
cat("Loading preprocessed data...\n")
load(file.path(OUTPUT_DIR, "processed_long_data.RData"))
cat("Loaded", nrow(long_data), "records\n")

# ------------------------------------------------------------------------------
# MOB Identification and Filtering
# ------------------------------------------------------------------------------
cat("\n=== MOB Identification ===\n")
cat("Searching for", length(TARGET_MOB_GENERA), "target genera...\n")

# Create regex pattern for genus matching
target_pattern <- paste(TARGET_MOB_GENERA, collapse = "|")

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
  cat("\nERROR: No MOB taxa found!\n")
  cat("Top 20 Families in your data:\n")
  print(head(sort(table(long_data$Family), decreasing=TRUE), 20))
  stop("No target MOB species found! Check the families above.")
}

cat("Found", length(unique(filtered_data$Taxon)), "MOB taxa:\n")
print(unique(filtered_data$Taxon))

# ------------------------------------------------------------------------------
# MOB Relative Abundance Calculation (Within MOB Community)
# ------------------------------------------------------------------------------
cat("\n=== Calculating MOB relative abundances ===\n")

# Calculate relative abundance WITHIN MOB (normalized to MOB total per sample)
plot_data <- filtered_data %>%
  group_by(SampleID, Taxon) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  group_by(SampleID) %>%
  mutate(RelativeAbundance = Count / sum(Count) * 100) %>%
  ungroup()

# Add grouping information (Surrounding: Bark, Soil, Wood)
plot_data <- plot_data %>%
  mutate(
    Group = case_when(
      grepl("B$", SampleID) ~ "Bark",
      grepl("S$", SampleID) ~ "Soil",
      TRUE ~ "Wood"
    ),
    Species = str_extract(SampleID, "^[A-Z]")
  )

# Set factor levels for consistent legend ordering
plot_data$Taxon <- factor(plot_data$Taxon, levels = sort(unique(plot_data$Taxon)))

# ------------------------------------------------------------------------------
# Create MOB Stacked Bar Plot
# ------------------------------------------------------------------------------
cat("Generating MOB stacked bar plot...\n")

p_mob <- ggplot(plot_data, aes(x = SampleID, y = RelativeAbundance, fill = Taxon)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = MOB_COLORS) +
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

# Display plot
print(p_mob)

# Save plot
output_pdf <- file.path(OUTPUT_DIR, "stacked_barplot_MOB.pdf")
output_png <- file.path(OUTPUT_DIR, "stacked_barplot_MOB.png")

ggsave(output_pdf, p_mob, width = 10, height = 6)
ggsave(output_png, p_mob, width = 10, height = 6)

cat("MOB stacked bar plot saved to:\n")
cat("  -", output_pdf, "\n")
cat("  -", output_png, "\n")

# ------------------------------------------------------------------------------
# MOB Total Relative Abundance Per Sample (vs Total Community)
# ------------------------------------------------------------------------------
cat("\n=== Calculating MOB total relative abundance per sample ===\n")

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

# Display summary statistics
cat("MOB Relative Abundance Summary:\n")
cat("  Mean:", round(mean(mob_per_sample$`MOB_RelAbundance (%)`), 2), "%\n")
cat("  Median:", round(median(mob_per_sample$`MOB_RelAbundance (%)`), 2), "%\n")
cat("  Range:", round(min(mob_per_sample$`MOB_RelAbundance (%)`), 2), "-", 
    round(max(mob_per_sample$`MOB_RelAbundance (%)`), 2), "%\n")

# Save results
output_csv1 <- file.path(OUTPUT_DIR, "MOB_relative_abundance_per_sample.csv")
write.csv(mob_per_sample, output_csv1, row.names = FALSE)
cat("Saved to:", output_csv1, "\n")

# ------------------------------------------------------------------------------
# Dominant vs Rare MOB Group Classification
# ------------------------------------------------------------------------------
cat("\n=== Classifying dominant and rare MOB groups ===\n")

# Calculate total relative abundance for each MOB genus across all samples
mob_genus_totals <- plot_data %>%
  group_by(Taxon) %>%
  summarise(Total_RelAbundance = sum(RelativeAbundance)) %>%
  arrange(desc(Total_RelAbundance)) %>%
  ungroup()

# Define top N as dominant/key groups
key_groups <- mob_genus_totals %>%
  slice(1:TOP_N_DOMINANT) %>%
  pull(Taxon)

# Define remaining as rare groups
rare_groups <- mob_genus_totals %>%
  slice((TOP_N_DOMINANT + 1):n()) %>%
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

# Display results
cat("\nDominant/Key MOB Groups (Top", TOP_N_DOMINANT, "):\n")
for (i in 1:length(key_groups)) {
  cat("  ", i, ".", as.character(key_groups[i]), "\n")
}

if (length(rare_groups) > 0) {
  cat("\nRare MOB Groups (", length(rare_groups), "):\n")
  for (rg in rare_groups) {
    cat("  -", as.character(rg), "\n")
  }
}

# Save classification
output_csv2 <- file.path(OUTPUT_DIR, "MOB_group_classification.csv")
write.csv(mob_classification, output_csv2, row.names = FALSE)
cat("Saved to:", output_csv2, "\n")

# ------------------------------------------------------------------------------
# Key Species Identification by Site (Surrounding)
# ------------------------------------------------------------------------------
cat("\n=== Identifying key MOB species by site ===\n")

# Calculate total relative abundance for each genus in each surrounding
mob_by_site <- plot_data %>%
  group_by(Surrounding, Taxon) %>%
  summarise(Total_RelAbundance = sum(RelativeAbundance), .groups = "drop") %>%
  arrange(Surrounding, desc(Total_RelAbundance))

# Get top 4 key species for each site
key_species_by_site <- mob_by_site %>%
  group_by(Surrounding) %>%
  slice(1:TOP_N_DOMINANT) %>%
  mutate(Rank = row_number()) %>%
  ungroup()

# Display results by surrounding
for (site in c("Bark", "Soil", "Wood")) {
  cat("\n", site, "- Top", TOP_N_DOMINANT, "MOB species:\n")
  site_data <- key_species_by_site %>% filter(Surrounding == site)
  for (i in 1:nrow(site_data)) {
    cat("  ", site_data$Rank[i], ".", as.character(site_data$Taxon[i]), 
        "(", round(site_data$Total_RelAbundance[i], 2), "%)\n")
  }
}

# Save results
output_csv3 <- file.path(OUTPUT_DIR, "MOB_key_species_by_site.csv")
write.csv(key_species_by_site, output_csv3, row.names = FALSE)
cat("\nSaved to:", output_csv3, "\n")

# ------------------------------------------------------------------------------
# Analysis Summary
# ------------------------------------------------------------------------------
cat("\n" , rep("=", 70), "\n", sep = "")
cat("MOB ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")
cat("\nOutput files generated:\n")
cat("  1.", output_pdf, "\n")
cat("  2.", output_png, "\n")
cat("  3.", output_csv1, "\n")
cat("  4.", output_csv2, "\n")
cat("  5.", output_csv3, "\n")
cat("\n")
