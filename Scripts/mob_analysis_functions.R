# mob_analysis_functions.R
# Functions for MOB (Methanotrophic Bacteria) analysis

run_mob_analysis <- function(data_subset, output_prefix, dataset_name) {
  cat("\n=== Processing", dataset_name, "===\n")
  # MOB identification and filtering
  target_genera <- c(
    "Methylacidiphilum", "Crenothrix", "Methylibium", "Methylobacillus", 
    "Methylobacterium", "Methylocella", "Methylomonas", "Methylonatrum", 
    "Methylophaga", "Methylopila", "Methylosinus", "Methylovirgula"
  )
  target_pattern <- paste(target_genera, collapse = "|")
  filtered_data <- data_subset %>%
    mutate(
      MatchedGenus = str_extract(Genus, regex(target_pattern, ignore_case = TRUE)),
      IsUnclassifiedMethylophilaceae = grepl("Methylophilaceae", Family, ignore.case = TRUE) & 
                                       (Genus == "" | Genus == "uncultured" | Genus == "Unassigned"),
      IsUnclassifiedMethylacidiphilaceae = grepl("Methylacidiphilaceae", Family, ignore.case = TRUE) & 
                                       (Genus == "" | Genus == "uncultured" | Genus == "Unassigned")
    ) %>%
    filter(!is.na(MatchedGenus) | IsUnclassifiedMethylophilaceae | IsUnclassifiedMethylacidiphilaceae) %>%
    mutate(
      Taxon = dplyr::case_when(
        !is.na(MatchedGenus) ~ stringr::str_to_title(MatchedGenus),
        IsUnclassifiedMethylophilaceae ~ "Unclassified Methylophilaceae",
        IsUnclassifiedMethylacidiphilaceae ~ "Unclassified Methylacidiphilaceae"
      ),
      Taxon = ifelse(Taxon == "Methylacidiphilum", "Candidatus Methylacidiphilum", Taxon)
    )
  if (nrow(filtered_data) == 0) {
    cat("WARNING: No MOB species found in", dataset_name, "\n")
    return(NULL)
  }
  # MOB relative abundance calculation and visualization
  plot_data <- filtered_data %>%
    dplyr::group_by(SampleID, Taxon) %>%
    dplyr::summarise(Count = sum(Count), .groups = "drop") %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(RelativeAbundance = Count / sum(Count) * 100) %>%
    dplyr::ungroup()
  plot_data <- plot_data %>%
    dplyr::mutate(Group = dplyr::case_when(
      grepl("B$", SampleID) ~ "Bark",
      grepl("S$", SampleID) ~ "Soil",
      TRUE ~ "Wood"
    ))
  mob_colors <- c(
    "Candidatus Methylacidiphilum" = "#CD5C5C",
    "Crenothrix" = "black",
    "Methylibium" = "#FF8C00",
    "Methylobacillus" = "#98FB98",
    "Methylobacterium" = "#7FFF00",
    "Methylocella" = "#8B4513",
    "Methylomonas" = "#E6E6FA",
    "Methylonatrum" = "#483D8B",
    "Methylophaga" = "#800080",
    "Methylopila" = "#C0C0C0",
    "Methylosinus" = "#DA70D6",
    "Unclassified Methylophilaceae" = "#AFEEEE"
  )
  plot_data$Taxon <- factor(plot_data$Taxon, levels = sort(unique(plot_data$Taxon)))
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = SampleID, y = RelativeAbundance, fill = Taxon)) +
    ggplot2::geom_bar(stat = "identity", width = 0.8) +
    ggplot2::facet_grid(~ Group, scales = "free_x", space = "free_x") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = mob_colors) +
    ggplot2::labs(x = "Tree", y = "Relative abundance (%)", fill = "MOB - 16S", 
         title = paste("MOB Distribution -", dataset_name)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5), 
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", colour = "black"),
      strip.text = ggplot2::element_text(size = 12),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 10, face = "italic")
    )
  print(p)
  ggplot2::ggsave(paste0(output_prefix, "_stacked_barplot_MOB.pdf"), p, width = 10, height = 6)
  ggplot2::ggsave(paste0(output_prefix, "_stacked_barplot_MOB.png"), p, width = 10, height = 6)
  # MOB total relative abundance per sample
  mob_totals <- filtered_data %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(MOB_Total = sum(Count)) %>%
    dplyr::ungroup()
  all_totals <- data_subset %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(Total_Abundance = sum(Count)) %>%
    dplyr::ungroup()
  mob_rel_abundance <- mob_totals %>%
    dplyr::left_join(all_totals, by = "SampleID") %>%
    dplyr::mutate(
      MOB_RelAbundance = MOB_Total / Total_Abundance * 100,
      Species = stringr::str_extract(SampleID, "^[A-Z]"),
      Surrounding = dplyr::case_when(
        grepl("B$", SampleID) ~ "Bark",
        grepl("S$", SampleID) ~ "Soil",
        TRUE ~ "Wood"
      )
    )
  mob_per_sample <- mob_rel_abundance %>%
    dplyr::select(SampleID, Species, Surrounding, MOB_RelAbundance) %>%
    dplyr::rename(`MOB_RelAbundance (%)` = MOB_RelAbundance)
  print("=== MOB Relative Abundance Per Sample ===")
  print(mob_per_sample)
  write.csv(mob_per_sample, paste0(output_prefix, "_MOB_relative_abundance_per_sample.csv"), row.names = FALSE)
  # Dominant vs rare MOB group classification
  mob_genus_totals <- plot_data %>%
    dplyr::group_by(Taxon) %>%
    dplyr::summarise(Total_RelAbundance = sum(RelativeAbundance)) %>%
    dplyr::arrange(dplyr::desc(Total_RelAbundance)) %>%
    dplyr::ungroup()
  key_groups <- mob_genus_totals %>%
    dplyr::slice(1:4) %>%
    dplyr::pull(Taxon)
  mob_classification <- mob_genus_totals %>%
    dplyr::mutate(
      Group_Type = dplyr::case_when(
        Taxon %in% key_groups ~ "Dominant/Key Group",
        TRUE ~ "Rare Group"
      ),
      Rank = dplyr::row_number()
    )
  print("=== MOB Group Classification ===")
  print(mob_classification)
  write.csv(mob_classification, paste0(output_prefix, "_MOB_group_classification.csv"), row.names = FALSE)
  # Key species identification by site
  mob_by_site <- plot_data %>%
    dplyr::group_by(Group, Taxon) %>%
    dplyr::summarise(Total_RelAbundance = sum(RelativeAbundance), .groups = "drop") %>%
    dplyr::arrange(Group, dplyr::desc(Total_RelAbundance))
  key_species_by_site <- mob_by_site %>%
    dplyr::group_by(Group) %>%
    dplyr::slice(1:3) %>%
    dplyr::mutate(Rank = dplyr::row_number()) %>%
    dplyr::ungroup()
  print("\n=== Key MOB Species by Site (Top 3 per Surrounding) ===")
  print(key_species_by_site)
  write.csv(key_species_by_site, paste0(output_prefix, "_MOB_key_species_by_site.csv"), row.names = FALSE)
  cat("\n", dataset_name, "analysis complete\n")
}
