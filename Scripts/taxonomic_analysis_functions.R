# taxonomic_analysis_functions.R
# Functions for diversity analysis and visualization

parse_taxonomy2 <- function(tax_string) {
  tax_clean <- gsub("\\s*\\([^\\)]+\\)", "", tax_string)
  levels <- strsplit(tax_clean, ";")[[1]]
  levels <- trimws(levels)
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  result <- rep(NA, 7)
  names(result) <- tax_levels
  result[1:min(length(levels), 7)] <- levels[1:min(length(levels), 7)]
  return(result)
}

run_analysis <- function(ft, tt, meta, output_prefix, dataset_name) {
  # Relative abundance calculation
  rel_abundance <- sweep(ft, 2, colSums(ft), "/") * 100

  # Aggregate by genus (or other taxonomic level)
  aggregate_by_level <- function(feat_table, tax_table, level) {
    common_features <- intersect(rownames(feat_table), rownames(tax_table))
    feat_table <- feat_table[common_features, , drop = FALSE]
    tax_table <- tax_table[common_features, , drop = FALSE]
    taxa <- tax_table[, level]
    taxa[is.na(taxa)] <- "Unclassified"
    agg_table <- aggregate(feat_table, by = list(Taxa = taxa), FUN = sum)
    rownames(agg_table) <- agg_table$Taxa
    agg_table$Taxa <- NULL
    return(as.data.frame(agg_table))
  }
  genus_abund <- aggregate_by_level(rel_abundance, tt, "Genus")

  # Heatmap plotting
  plot_heatmap <- function(abund_table, metadata, top_n = 30, cluster_samples = TRUE, group_var = "sample_type") {
    mean_abund <- rowMeans(abund_table)
    top_taxa <- names(sort(mean_abund, decreasing = TRUE)[1:top_n])
    heat_data <- abund_table[top_taxa, ]
    heat_data_log <- log10(heat_data + 0.01)
    annotation_col <- metadata[, group_var, drop = FALSE]
    sample_types <- unique(metadata[, group_var])
    annotation_colors <- list()
    if (length(sample_types) == 3) {
      annotation_colors[[group_var]] <- c("green3", "purple3", "gold")
      names(annotation_colors[[group_var]]) <- sort(sample_types)
    }
    pheatmap::pheatmap(
      heat_data_log,
      cluster_rows = TRUE,
      cluster_cols = cluster_samples,
      annotation_col = annotation_col,
      annotation_colors = annotation_colors,
      color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
      fontsize = 8,
      fontsize_row = 7,
      fontsize_col = 8,
      main = paste("Top", top_n, "Taxa Heatmap -", dataset_name, "\n[log10(Relative Abundance %)]"),
      filename = paste0(output_prefix, "_taxa_heatmap.pdf"),
      width = 10,
      height = 12
    )
  }
  plot_heatmap(genus_abund, meta, top_n = 30, group_var = "sample_type")
  cat("Heatmap saved\n")

  # Alpha diversity calculation and boxplot
  alpha_div <- data.frame(
    Sample = colnames(ft),
    Observed = colSums(ft > 0),
    Shannon = vegan::diversity(t(ft), index = "shannon"),
    Simpson = vegan::diversity(t(ft), index = "simpson")
  )
  alpha_div <- cbind(alpha_div, meta[alpha_div$Sample, ])
  alpha_long <- alpha_div %>%
    tidyr::pivot_longer(cols = c(Observed, Shannon, Simpson), names_to = "Metric", values_to = "Value")
  p_alpha <- ggplot2::ggplot(alpha_long, ggplot2::aes(x = sample_type, y = Value, fill = sample_type)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(~ Metric, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = paste("Alpha Diversity -", dataset_name), x = "Sample Type", y = "Value") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  print(p_alpha)
  ggplot2::ggsave(paste0(output_prefix, "_alpha_diversity_boxplot.pdf"), p_alpha, width = 10, height = 4)
  cat("Alpha diversity boxplot saved\n")

  # PCoA ordination using Unweighted UniFrac from QIIME2
  cat("\n=== PCoA with Unweighted UniFrac ===\n")
  unweighted_pcoa_file <- "Results/denoise_mode/diversity/exported-unweighted_unifrac-pcoa/ordination.txt"
  unweighted_dist_file <- "Results/denoise_mode/diversity/exported-unweighted_unifrac/distance-matrix.tsv"
  lines <- readLines(unweighted_pcoa_file)
  prop_idx <- grep("^Proportion", lines)
  site_idx <- grep("^Site", lines)
  
  # Fix: Split the proportion line by tab before converting to numeric
  if (length(prop_idx) > 0) {
    prop_line <- lines[prop_idx[1] + 1]
    prop_vals <- strsplit(prop_line, "\t")[[1]]
    # Remove empty strings and convert to numeric
    prop_vals <- prop_vals[prop_vals != ""]
    percent_var <- as.numeric(prop_vals) * 100
    percent_var <- round(percent_var, 2)
    cat("Parsed percent variance:", paste(percent_var[1:2], collapse=", "), "\n")
  } else {
    percent_var <- c(NA, NA)
    cat("Warning: 'Proportion explained' not found in PCoA file\n")
  }
  
  coord_lines <- lines[(site_idx + 1):length(lines)]
  # Only keep lines that look like sample coordinates (start with a sample name, not empty, not Biplot etc)
  # Use '^[A-Za-z0-9]' instead of '\w' for R regex compatibility
  coord_lines <- coord_lines[grepl('^[A-Za-z0-9]', coord_lines) & coord_lines != ""]
  # Split and check matrix
  coord_matrix <- do.call(rbind, strsplit(coord_lines, "\t"))
  # Check for valid numeric columns (PC1, PC2)
  suppressWarnings({
    valid_idx <- !is.na(as.numeric(coord_matrix[,2])) & !is.na(as.numeric(coord_matrix[,3]))
  })
  coord_matrix <- coord_matrix[valid_idx, ]
  # Optional: print for debugging
  if (any(!valid_idx)) {
    cat("Warning: Some coordinate lines were invalid and removed.\n")
  }
  pcoa_scores <- data.frame(
    Sample = coord_matrix[, 1],
    PC1 = as.numeric(coord_matrix[, 2]),
    PC2 = as.numeric(coord_matrix[, 3]),
    stringsAsFactors = FALSE
  )
  rownames(pcoa_scores) <- pcoa_scores$Sample
  common_samples <- intersect(pcoa_scores$Sample, rownames(meta))
  cat("Debug: common_samples\n")
  print(common_samples)
  pcoa_scores <- pcoa_scores[common_samples, ]
  cat("Debug: pcoa_scores after subsetting\n")
  print(head(pcoa_scores))
  pcoa_scores <- cbind(pcoa_scores, meta[pcoa_scores$Sample, ])
  cat("Debug: pcoa_scores after cbind meta\n")
  print(head(pcoa_scores))
  unifrac_dist_data <- read.table(unweighted_dist_file, sep = "\t", header = TRUE, row.names = 1)
  unifrac_dist_data <- unifrac_dist_data[common_samples, common_samples]
  pcoa_dist <- as.dist(unifrac_dist_data)
  cat("Using QIIME2 Unweighted UniFrac results\n")
  cat("\n=== Plotting PCoA ===\n")
  p_pcoa <- ggplot2::ggplot(pcoa_scores, ggplot2::aes(x = PC1, y = PC2, color = sample_type, shape = species)) +
    ggplot2::geom_point(size = 4, alpha = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste("PCoA - Unweighted UniFrac -", dataset_name),
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)"),
      color = "Sample Type",
      shape = "Species"
    ) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
    )
  print(p_pcoa)
  ggplot2::ggsave(paste0(output_prefix, "_pcoa_unweighted_unifrac.pdf"), p_pcoa, width = 10, height = 7)
  cat("PCoA plot saved to:", paste0(output_prefix, "_pcoa_unweighted_unifrac.pdf\n"))

  # PERMANOVA analysis
  cat("\n=== PERMANOVA Analysis ===\n")
  dist_samples <- labels(pcoa_dist)
  meta_subset <- meta[dist_samples, ]
  if (length(unique(meta_subset$species)) > 1) {
    permanova_species <- vegan::adonis2(pcoa_dist ~ species, data = meta_subset, permutations = 999)
    cat("\n1. Tree Species Effect:\n")
    print(permanova_species)
  }
  if (length(unique(meta_subset$sample_type)) > 1) {
    permanova_tissue <- vegan::adonis2(pcoa_dist ~ sample_type, data = meta_subset, permutations = 999)
    cat("\n2. Tissue Effect:\n")
    print(permanova_tissue)
  }
  if (length(unique(meta_subset$species)) > 1 && length(unique(meta_subset$sample_type)) > 1) {
    permanova_combined <- vegan::adonis2(pcoa_dist ~ species + sample_type, data = meta_subset, permutations = 999)
    cat("\n3. Combined Effect:\n")
    print(permanova_combined)
    permanova_interaction <- vegan::adonis2(pcoa_dist ~ species * sample_type, data = meta_subset, permutations = 999)
    cat("\n4. Interaction Effect:\n")
    print(permanova_interaction)
  }
  sink(paste0(output_prefix, "_PERMANOVA_results.txt"))
  cat("PERMANOVA Test Results -", dataset_name, "\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Distance metric: Unweighted UniFrac\n")
  cat("Permutations: 999\n\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  if (length(unique(meta$species)) > 1) {
    cat("1. Tree Species Effect (dm ~ tree_species):\n")
    print(permanova_species)
    cat("\n")
  }
  if (length(unique(meta$sample_type)) > 1) {
    cat("2. Tissue Effect (dm ~ tissue):\n")
    print(permanova_tissue)
    cat("\n")
  }
  if (length(unique(meta$species)) > 1 && length(unique(meta$sample_type)) > 1) {
    cat("3. Combined Effect (dm ~ tree_species + tissue):\n")
    print(permanova_combined)
    cat("\n")
    cat("4. Interaction Effect (dm ~ tree_species * tissue):\n")
    print(permanova_interaction)
    cat("\n")
  }
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("\nInterpretation:\n")
  cat("- Pr(>F) < 0.001: *** (highly significant)\n")
  cat("- Pr(>F) < 0.01:  **  (significant)\n")
  cat("- Pr(>F) < 0.05:  *   (marginally significant)\n")
  cat("- Pr(>F) >= 0.05: not significant\n\n")
  cat("PERMANOVA tests whether community composition differs significantly between groups.\n")
  cat("R² represents the proportion of variance explained by each factor.\n")
  sink()
  cat("PERMANOVA results saved to file\n")
  cat("\n", dataset_name, "analysis complete\n")
}
