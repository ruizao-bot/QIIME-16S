#!/usr/bin/env Rscript

# Minimal downstream analysis: genus barplot, heatmap, NMDS/PCoA

# Load required packages
required_packages <- c("tidyverse", "vegan", "ggplot2", "RColorBrewer", "pheatmap")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

setwd("/Users/jiayi/Desktop/metagenomic_pipeline/QIIME")

# Import data
feature_table <- read.table(
  "Results/denoise_mode/exported-table/feature-table.tsv",
  header = TRUE,
  sep = "\t",
  skip = 1,
  row.names = 1,
  comment.char = ""
)
feature_table <- feature_table[, !(colnames(feature_table) %in% c("A1", "A2")), drop = FALSE]

taxonomy <- read.table(
  "Results/denoise_mode/exported-taxonomy/taxonomy.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  quote = ""
)

metadata <- read.table(
  "Data/metadata/metadata.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  comment.char = ""
)
metadata <- metadata[colnames(feature_table), ]

# Parse taxonomy
parse_taxonomy <- function(tax_string) {
  tax_clean <- gsub("\\s*\\([^\\)]+\\)", "", tax_string)
  levels <- strsplit(tax_clean, ";")[[1]]
  levels <- trimws(levels)
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  result <- rep(NA, 7)
  names(result) <- tax_levels
  result[1:min(length(levels), 7)] <- levels[1:min(length(levels), 7)]
  return(result)
}
tax_table <- t(sapply(taxonomy$Taxon, parse_taxonomy))
rownames(tax_table) <- rownames(taxonomy)

# Relative abundance
rel_abundance <- sweep(feature_table, 2, colSums(feature_table), "/") * 100

# Aggregate by genus
aggregate_by_level <- function(feat_table, tax_table, level) {
  taxa <- tax_table[, level]
  taxa[is.na(taxa)] <- "Unclassified"
  agg_table <- aggregate(feat_table, by = list(Taxa = taxa), FUN = sum)
  rownames(agg_table) <- agg_table$Taxa
  agg_table$Taxa <- NULL
  return(as.data.frame(agg_table))
}
genus_abund <- aggregate_by_level(rel_abundance, tax_table, "Genus")

# Genus barplot
plot_stacked_barplot <- function(abund_table, metadata, top_n = 15, group_var = "sample_type", title = "Genus-level Composition") {
  mean_abund <- rowMeans(abund_table)
  top_taxa <- names(sort(mean_abund, decreasing = TRUE)[1:top_n])
  plot_data <- abund_table[top_taxa, , drop = FALSE]
  others <- colSums(abund_table[!rownames(abund_table) %in% top_taxa, , drop = FALSE])
  plot_data <- rbind(plot_data, Others = others)
  plot_data_long <- plot_data %>%
    tibble::rownames_to_column("Taxa") %>%
    pivot_longer(-Taxa, names_to = "Sample", values_to = "Abundance") %>%
    left_join(metadata %>% tibble::rownames_to_column("Sample"), by = "Sample")
  n_taxa <- nrow(plot_data)
  if (n_taxa <= 12) {
    colors <- c(brewer.pal(n_taxa - 1, "Set3"), "grey80")[1:n_taxa]
  } else {
    colors <- c(colorRampPalette(brewer.pal(12, "Set3"))(n_taxa - 1), "grey80")
  }
  p <- ggplot(plot_data_long, aes(x = Sample, y = Abundance, fill = Taxa)) +
    geom_bar(stat = "identity", color = "black", size = 0.1) +
    scale_fill_manual(values = colors) +
    facet_grid(~ get(group_var), scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "right",
          strip.background = element_rect(fill = "lightblue")) +
    labs(title = title, y = "Relative Abundance (%)", x = "Sample") +
    guides(fill = guide_legend(ncol = 1))
  return(p)
}
p_genus <- plot_stacked_barplot(genus_abund, metadata, top_n = 15, group_var = "sample_type", title = "Genus-level Composition")
print(p_genus)
ggsave("Results/denoise_mode/genus_barplot.pdf", p_genus, width = 14, height = 7)
cat("Genus-level stacked bar plot saved\n\n")

# Heatmap
plot_heatmap <- function(abund_table, metadata, top_n = 30, cluster_samples = TRUE, group_var = "sample_type") {
  mean_abund <- rowMeans(abund_table)
  top_taxa <- names(sort(mean_abund, decreasing = TRUE)[1:top_n])
  heat_data <- abund_table[top_taxa, ]
  heat_data_log <- log10(heat_data + 0.01)
  annotation_col <- metadata[, group_var, drop = FALSE]
  pheatmap::pheatmap(
    heat_data_log,
    cluster_rows = TRUE,
    cluster_cols = cluster_samples,
    annotation_col = annotation_col,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    fontsize = 8,
    fontsize_row = 7,
    fontsize_col = 8,
    main = paste("Top", top_n, "Taxa Heatmap (log10 scale)"),
    filename = "Results/denoise_mode/taxa_heatmap.pdf",
    width = 10,
    height = 12
  )
}
plot_heatmap(genus_abund, metadata, top_n = 30, group_var = "sample_type")
cat("Heatmap saved\n\n")

# NMDS/PCoA ordination
bc_dist <- vegdist(t(feature_table), method = "bray")
set.seed(123)
nmds <- metaMDS(bc_dist, k = 2, trymax = 100)
nmds_scores <- as.data.frame(scores(nmds)$sites)
nmds_scores$Sample <- rownames(nmds_scores)
nmds_scores <- cbind(nmds_scores, metadata)
p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = sample_type, shape = species)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = sample_type), level = 0.95) +
  theme_bw() +
  labs(title = paste("NMDS Ordination (Stress =", round(nmds$stress, 3), ")"),
       color = "Sample Type", shape = "Species") +
  scale_color_brewer(palette = "Set1")
print(p_nmds)
ggsave("Results/denoise_mode/nmds_ordination.pdf", p_nmds, width = 8, height = 6)
cat("NMDS ordination complete\n")
cat("Stress value:", round(nmds$stress, 3), "\n\n")
