# ==============================================================================
# Configuration File for MOB Analysis
# ==============================================================================
# This file contains shared parameters and settings used across all analysis scripts
# ==============================================================================

# ------------------------------------------------------------------------------
# File Paths
# ------------------------------------------------------------------------------
BASE_DIR <- "/Users/jiayi/Desktop/metagenomic_pipeline/QIIME"
INPUT_FILE <- "Results/denoise_mode/final-table-with-ranks.tsv"
OUTPUT_DIR <- "Results/denoise_mode"

# ------------------------------------------------------------------------------
# Sample Filtering
# ------------------------------------------------------------------------------
EXCLUDE_SAMPLES <- c("A1", "A2")

# ------------------------------------------------------------------------------
# MOB Target Genera
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# MOB Color Palette
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# Analysis Parameters
# ------------------------------------------------------------------------------
TOP_N_FAMILIES <- 75  # Number of top families to display
TOP_N_DOMINANT <- 4   # Number of dominant MOB groups

