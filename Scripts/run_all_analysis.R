#!/usr/bin/env Rscript
# ==============================================================================
# Master Script to Run All MOB Analyses

# 直接设置工作目录到项目根目录
setwd("/Users/jiayi/Desktop/metagenomic_pipeline/QIIME")
# (All config variables are now inlined in each script)
# ==============================================================================
# This script runs all analysis steps in the correct order
# ==============================================================================

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("MOB ANALYSIS PIPELINE\n")
cat(rep("=", 80), "\n", sep = "")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Record start time
start_time <- Sys.time()


# ------------------------------------------------------------------------------
# Step 0: Combine Feature Table and Taxonomy
# ------------------------------------------------------------------------------
cat("\n[STEP 0/3] Combining Feature Table and Taxonomy...\n")
cat(rep("-", 80), "\n", sep = "")

tryCatch({
  source("Scripts/table_combine.R")
  cat("\n✓ Step 0 completed successfully\n")
}, error = function(e) {
  cat("\n✗ Step 0 failed with error:\n")
  cat(e$message, "\n")
  stop("Pipeline halted due to error in Step 0")
})

# ------------------------------------------------------------------------------
# Step 1: Family-Level Stacked Bar Plot
# ------------------------------------------------------------------------------
cat("\n[STEP 1/3] Running Family-Level Analysis...\n")
cat(rep("-", 80), "\n", sep = "")

tryCatch({
  source("Scripts/01_stacked_barplot.R")
  cat("\n✓ Step 1 completed successfully\n")
}, error = function(e) {
  cat("\n✗ Step 1 failed with error:\n")
  cat(e$message, "\n")
  stop("Pipeline halted due to error in Step 1")
})

# ------------------------------------------------------------------------------
# Step 2: MOB Analysis (Identification, Visualization, Statistics)
# ------------------------------------------------------------------------------
cat("\n[STEP 2/3] Running MOB Analysis...\n")
cat(rep("-", 80), "\n", sep = "")

tryCatch({
  source("Scripts/02_mob_analysis.R")
  cat("\n✓ Step 2 completed successfully\n")
}, error = function(e) {
  cat("\n✗ Step 2 failed with error:\n")
  cat(e$message, "\n")
  stop("Pipeline halted due to error in Step 2")
})

# ------------------------------------------------------------------------------
# Pipeline Summary
# ------------------------------------------------------------------------------
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "secs")

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("PIPELINE COMPLETED SUCCESSFULLY\n")
cat(rep("=", 80), "\n", sep = "")
cat("End time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total elapsed time:", round(elapsed_time, 2), "seconds\n")
cat("\nAll analysis steps completed successfully!\n")
cat("Check the 'denoise_mode' directory for outputs.\n\n")
