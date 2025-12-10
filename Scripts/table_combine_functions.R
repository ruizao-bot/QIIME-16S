# table_combine_functions.R
# Functions for merging feature table and taxonomy, and parsing taxonomy ranks

combine_feature_taxonomy <- function(ft_path, taxonomy_path, output_dir) {
  feature_table_lines <- readLines(ft_path)
  header_idx <- which(grepl('^#OTU ID', feature_table_lines))[1]
  if (is.na(header_idx)) stop('Could not find header line beginning with "#OTU ID" in feature-table.tsv')
  header_raw <- sub('^#', '', feature_table_lines[header_idx])
  col_names <- strsplit(header_raw, '\t')[[1]]
  if (length(col_names) == 1L) col_names <- strsplit(header_raw, '\\s+')[[1]]
  feature_table <- read.table(ft_path, sep='\t', header=FALSE, row.names=1, skip=header_idx, check.names=FALSE, comment.char='')
  if (length(col_names) - 1L != ncol(feature_table)) stop(sprintf('Column name count mismatch: header has %d data columns, table has %d', length(col_names) - 1L, ncol(feature_table)))
  colnames(feature_table) <- col_names[-1]
  taxonomy <- read.table(taxonomy_path, sep='\t', header=TRUE, row.names=1)
  taxonomy_only <- data.frame(Taxon = taxonomy$Taxon, row.names = rownames(taxonomy))
  common_features <- intersect(rownames(feature_table), rownames(taxonomy_only))
  feature_table_matched <- feature_table[common_features, , drop=FALSE]
  taxonomy_matched <- taxonomy_only[common_features, , drop=FALSE]
  merged <- cbind(taxonomy_matched, feature_table_matched)
  write.table(merged, file.path(output_dir, 'taxonomy-abundance-table.tsv'), sep='\t', quote=FALSE, row.names=TRUE)
  parse_taxonomy <- function(tax_string) {
    ranks <- c(Domain=NA, Phylum=NA, Class=NA, Order=NA, Family=NA, Genus=NA, Species=NA)
    if (is.na(tax_string) || tax_string == '') {
      ranks[is.na(ranks)] <- 'Unassigned'
      return(ranks)
    }
    parts <- unlist(strsplit(tax_string, ';'))
    for (p in parts) {
      p <- trimws(p)
      if (p == '' ) next
      if (!grepl('^[dpcofgs]__', p)) next
      prefix <- substr(p, 1, 1)
      value <- sub('^[dpcofgs]__', '', p)
      if (value == '' || value == 'uncultured' || value == 'metagenome') value <- 'Unassigned'
      rankName <- switch(prefix, d='Domain', p='Phylum', c='Class', o='Order', f='Family', g='Genus', s='Species')
      ranks[rankName] <- ifelse(value == '', 'Unassigned', value)
    }
    ranks[is.na(ranks)] <- 'Unassigned'
    return(ranks)
  }
  rank_matrix <- t(vapply(taxonomy_matched$Taxon, parse_taxonomy, character(7)))
  rank_df <- as.data.frame(rank_matrix, stringsAsFactors = FALSE)
  rownames(rank_df) <- rownames(taxonomy_matched)
  merged_ranks <- cbind(rank_df, feature_table_matched)
  write.table(merged_ranks, file.path(output_dir, 'final-table-with-ranks.tsv'), sep='\t', quote=FALSE, row.names=TRUE)
  return(list(feature_table=feature_table_matched, taxonomy=taxonomy_matched, merged_ranks=merged_ranks))
}
