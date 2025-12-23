#!/usr/bin/env bash

# QIIME2 Shotgun Data Analysis Pipeline
# This script processes shotgun metagenomic data using QIIME2
# Ensure QIIME2 is installed and the environment is activated before running

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_DIR="${PROJECT_DIR}/Data"
RAW_DATA_DIR="${DATA_DIR}/raw_data/shotgun"
PROCESSED_DATA_DIR="${DATA_DIR}/processed_data"
REFERENCE_DB_DIR="${DATA_DIR}/reference_dbs"
METADATA_FILE="${DATA_DIR}/metadata/metadata.tsv"
MANIFEST_FILE="${RAW_DATA_DIR}/manifest.tsv"
LOG_FILE="${PROJECT_DIR}/Logs/shotgun_pipeline.log"
ENV_NAME="qiime2-moshpit"

# Create necessary directories
mkdir -p "${PROCESSED_DATA_DIR}"
mkdir -p "${PROJECT_DIR}/Logs"

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_FILE}"
}

# Error handling function
error_exit() {
    log "ERROR: $1"
    exit 1
}

# Check if output file exists and skip step
check_skip() {
    local output_file="$1"
    local step_name="$2"
    if [[ -f "${output_file}" ]]; then
        log "${step_name}: Output file already exists (${output_file}). Skipping..."
        return 0  # Skip
    fi
    return 1  # Don't skip
}

# Step 1: Import Data
step1_import_data() {
    log "Starting Step 1: Importing shotgun data"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/demux-paired-end.qza" "Step 1"; then
        return 0
    fi

    # Check if manifest file exists
    if [[ ! -f "${MANIFEST_FILE}" ]]; then
        error_exit "Manifest file not found at ${MANIFEST_FILE}. Please create it before running the pipeline."
    fi

    # Import data
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "${MANIFEST_FILE}" \
        --output-path "${PROCESSED_DATA_DIR}/demux-paired-end.qza" \
        --input-format PairedEndFastqManifestPhred33V2 || \
        error_exit "Failed to import data."

    log "Data imported successfully. Output: ${PROCESSED_DATA_DIR}/demux-paired-end.qza"
}

# Step 2: Quality Control (Updated to trim adapters)
step2_quality_control() {
    log "Starting Step 2: Quality control (Trimming adapters)"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/trimmed-seqs.qza" "Step 2"; then
        return 0
    fi

    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences "${PROCESSED_DATA_DIR}/demux-paired-end.qza" \
        --p-cores 4 \
        --o-trimmed-sequences "${PROCESSED_DATA_DIR}/trimmed-seqs.qza" || \
        error_exit "Failed to trim adapters."

    log "Adapters trimmed successfully. Output: ${PROCESSED_DATA_DIR}/trimmed-seqs.qza"
}

# Step 3: Remove Host DNA
step3_remove_host() {
    log "Starting Step 3: Removing host DNA"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/host-removed-seqs.qza" "Step 3"; then
        return 0
    fi

    # Check if host reference database exists
    HOST_DB="${REFERENCE_DB_DIR}/host_genome"
    if [[ ! -d "${HOST_DB}" ]]; then
        log "WARNING: Host genome database not found at ${HOST_DB}. Skipping host removal."
        log "To enable host removal, build a Bowtie2 index and place it in ${HOST_DB}"
        # Copy input to output if skipping
        cp "${PROCESSED_DATA_DIR}/trimmed-seqs.qza" "${PROCESSED_DATA_DIR}/host-removed-seqs.qza"
        return 0
    fi

    qiime quality-control filter-reads \
        --i-demultiplexed-sequences "${PROCESSED_DATA_DIR}/trimmed-seqs.qza" \
        --i-database "${HOST_DB}" \
        --p-n-threads 16 \
        --o-filtered-sequences "${PROCESSED_DATA_DIR}/host-removed-seqs.qza" || \
        error_exit "Failed to remove host DNA."

    log "Host DNA removed successfully. Output: ${PROCESSED_DATA_DIR}/host-removed-seqs.qza"
}

# Step 4: Assemble Contigs
step4_assemble_contigs() {
    log "Starting Step 4: Assembling contigs"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/assembled-contigs.qza" "Step 4"; then
        return 0
    fi

    qiime assembly assemble-megahit \
        --i-reads "${PROCESSED_DATA_DIR}/host-removed-seqs.qza" \
        --o-contigs "${PROCESSED_DATA_DIR}/assembled-contigs.qza" || \
        error_exit "Failed to assemble contigs."

    log "Contigs assembled successfully. Output: ${PROCESSED_DATA_DIR}/assembled-contigs.qza"
}

# Step 4b: Map reads to contigs (required for binning)
step4b_map_reads() {
    log "Starting Step 4b: Mapping reads to contigs"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/mapped-reads.qza" "Step 4b"; then
        return 0
    fi

    # Index contigs
    qiime assembly index-contigs \
        --i-contigs "${PROCESSED_DATA_DIR}/assembled-contigs.qza" \
        --o-index "${PROCESSED_DATA_DIR}/contigs-index.qza" || \
        error_exit "Failed to index contigs."

    # Map reads to contigs
    qiime assembly map-reads \
        --i-indexed-contigs "${PROCESSED_DATA_DIR}/contigs-index.qza" \
        --i-reads "${PROCESSED_DATA_DIR}/host-removed-seqs.qza" \
        --o-alignment-map "${PROCESSED_DATA_DIR}/mapped-reads.qza" || \
        error_exit "Failed to map reads to contigs."

    log "Reads mapped successfully. Output: ${PROCESSED_DATA_DIR}/mapped-reads.qza"
}

# Step 5: Bin Contigs
step5_bin_contigs() {
    log "Starting Step 5: Binning contigs"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/binned-contigs.qza" "Step 5"; then
        return 0
    fi

    qiime annotate bin-contigs-metabat \
        --i-contigs "${PROCESSED_DATA_DIR}/assembled-contigs.qza" \
        --i-maps "${PROCESSED_DATA_DIR}/mapped-reads.qza" \
        --o-mags "${PROCESSED_DATA_DIR}/binned-contigs.qza" || \
        error_exit "Failed to bin contigs."

    log "Contigs binned successfully. Output: ${PROCESSED_DATA_DIR}/binned-contigs.qza"
}

# Step 6: Evaluate MAGs
step6_evaluate_mags() {
    log "Starting Step 6: Evaluating MAGs"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/mags-evaluation.qzv" "Step 6"; then
        return 0
    fi

    qiime annotate evaluate-busco \
        --i-bins "${PROCESSED_DATA_DIR}/binned-contigs.qza" \
        --o-visualization "${PROCESSED_DATA_DIR}/mags-evaluation.qzv" || \
        error_exit "Failed to evaluate MAGs."

    log "MAGs evaluated successfully. Output: ${PROCESSED_DATA_DIR}/mags-evaluation.qzv"
}

# Step 7a: Predict genes from MAGs
step7a_predict_genes() {
    log "Starting Step 7a: Predicting genes"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/predicted-genes.qza" "Step 7a"; then
        return 0
    fi

    qiime annotate predict-genes-prodigal \
        --i-input "${PROCESSED_DATA_DIR}/binned-contigs.qza" \
        --o-gene-sequences "${PROCESSED_DATA_DIR}/predicted-genes.qza" || \
        error_exit "Failed to predict genes."

    log "Genes predicted successfully. Output: ${PROCESSED_DATA_DIR}/predicted-genes.qza"
}

# Step 7b: Search orthologs using DIAMOND
step7b_search_orthologs() {
    log "Starting Step 7b: Searching orthologs"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/orthologs.qza" "Step 7b"; then
        return 0
    fi

    qiime annotate search-orthologs-diamond \
        --i-gene-sequences "${PROCESSED_DATA_DIR}/predicted-genes.qza" \
        --o-ortholog-annotations "${PROCESSED_DATA_DIR}/orthologs.qza" || \
        error_exit "Failed to search orthologs."

    log "Orthologs searched successfully. Output: ${PROCESSED_DATA_DIR}/orthologs.qza"
}

# Step 7c: Annotate MAGs with eggNOG
step7c_annotate_mags() {
    log "Starting Step 7c: Annotating MAGs with eggNOG"
    
    # Check if already completed
    if check_skip "${PROCESSED_DATA_DIR}/mags-annotations.qza" "Step 7c"; then
        return 0
    fi

    qiime annotate map-eggnog \
        --i-ortholog-annotations "${PROCESSED_DATA_DIR}/orthologs.qza" \
        --o-eggnog-annotations "${PROCESSED_DATA_DIR}/mags-annotations.qza" || \
        error_exit "Failed to annotate MAGs."

    log "MAGs annotated successfully. Output: ${PROCESSED_DATA_DIR}/mags-annotations.qza"
}

# Main Pipeline Execution (Updated)
main() {
    log "Starting QIIME2 Shotgun Data Analysis Pipeline"

    step1_import_data
    step2_quality_control
    step3_remove_host
    step4_assemble_contigs
    step4b_map_reads
    step5_bin_contigs
    step6_evaluate_mags
    step7a_predict_genes
    step7b_search_orthologs
    step7c_annotate_mags

    log "Pipeline completed successfully."
}

main "$@"
