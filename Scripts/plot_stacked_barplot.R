library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(RColorBrewer)
setwd("/Users/jiayi/Desktop/metagenomic_pipeline/QIIME")
# Set input and output file paths
input_file <- "Results/denoise_mode05/final-table-with-ranks.tsv"
output_plot <- "Results/denoise_mode05/stacked_barplot_MOB.pdf"

# Read data
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
sample_cols <- setdiff(sample_cols, c("A1", "A2")) # Filter out A1, A2

# Convert to long format
long_data <- data %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "SampleID", values_to = "Count") %>%
  filter(Count > 0) %>% # Remove zero counts to speed up
  mutate(
    # 1. Clean Taxonomy: Remove "f__", "g__", etc. and trim whitespace
    Family = str_replace_all(Family, "[a-z]__", ""),
    Genus = str_replace_all(Genus, "[a-z]__", ""),
    Family = trimws(Family),
    Genus = trimws(Genus)
  )


# --- DEFINING TARGETS ---
# List of Genera to search for (Updated based on your data)
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

# Create a regex pattern for matching
target_pattern <- paste(target_genera, collapse = "|")

# Filter and label data
filtered_data <- long_data %>%
  mutate(
    # Extract the matching genus name
    MatchedGenus = str_extract(Genus, regex(target_pattern, ignore_case = TRUE)),
    
    # Handle special cases for Unclassified Families
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
    # Add "Candidatus" prefix if needed
    Taxon = ifelse(Taxon == "Methylacidiphilum", "Candidatus Methylacidiphilum", Taxon)
  )

# Check if we found anything
if (nrow(filtered_data) == 0) {
  print("DEBUG: Top 20 Families found in your data (to help debug):")
  print(head(sort(table(long_data$Family), decreasing=TRUE), 20))
  stop("No target species found! Check the debug output above to see what Families are actually in your data.")
}

# --- RE-CALCULATE RELATIVE ABUNDANCE WITHIN MOB ---
# The reference plot sums to 100%, meaning it shows the proportion WITHIN the Methylotrophs
plot_data <- filtered_data %>%
  group_by(SampleID, Taxon) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  group_by(SampleID) %>%
  mutate(RelativeAbundance = Count / sum(Count) * 100) %>%
  ungroup()

# Add grouping info
plot_data <- plot_data %>%
  mutate(Group = case_when(
    grepl("B$", SampleID) ~ "Bark",
    grepl("S$", SampleID) ~ "Soil",
    TRUE ~ "Wood"
  ))

# Define Colors (Matched to your reference image + New Taxa)
mob_colors <- c(
  "Candidatus Methylacidiphilum" = "#CD5C5C", # IndianRed
  "Unclassified Methylacidiphilaceae" = "#CD5C5C", # Same as above (related)
  "Crenothrix" = "black",
  "Methylibium" = "#FF8C00", # DarkOrange
  "Methylobacillus" = "#98FB98", # PaleGreen
  "Methylobacterium" = "#7FFF00", # Chartreuse
  "Methylocella" = "#8B4513", # SaddleBrown
  "Methylomonas" = "#E6E6FA", # Lavender/White-ish
  "Methylonatrum" = "#483D8B", # DarkSlateBlue
  "Methylophaga" = "#800080", # Purple
  "Methylopila" = "#C0C0C0", # Silver/Grey
  "Methylosinus" = "#DA70D6", # Orchid
  "Unclassified Methylophilaceae" = "#AFEEEE", # PaleTurquoise
  "Methylovirgula" = "#FFD700" # Gold (New)
)

# Ensure factor levels match the legend order (Alphabetical usually)
plot_data$Taxon <- factor(plot_data$Taxon, levels = sort(unique(plot_data$Taxon)))

# Plotting
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
p
# Save plots
ggsave(output_plot, p, width = 10, height = 6)
ggsave(gsub(".pdf", ".png", output_plot), p, width = 10, height = 6)

print(paste("Plotting complete. Saved to:", output_plot))

