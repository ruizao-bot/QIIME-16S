
#!/usr/bin/env bash

# QIIME2 Pipeline with Checkpoint System
# This script processes demultiplexed paired-end FASTQ files using QIIME2 to generate an OTU table and representative sequences
# If unfamiliar with QIIME, see https://amplicon-docs.qiime2.org/en/latest/explanations/getting-started.html

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CHECKPOINT_DIR="${PROJECT_DIR}/Logs/checkpoints"
LOG_FILE="${PROJECT_DIR}/Logs/qiime2_pipeline.log"
ENV_NAME="qiime"
# Mode: 'denoise' (DADA2) or 'cluster' (OTU clustering from demux)
# Users can override with -m|--mode when running the script
MODE="denoise"

# Create checkpoint directory
mkdir -p "${CHECKPOINT_DIR}"

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_FILE}"
}

# Error handling function
error_exit() {
    log "ERROR: $1"
    exit 1
}

# Checkpoint functions
create_checkpoint() {
    local step_name="$1"
    echo "$(date '+%Y-%m-%d %H:%M:%S')" > "${CHECKPOINT_DIR}/${step_name}.checkpoint"
    log "Checkpoint created for step: ${step_name}"
}

check_checkpoint() {
    local step_name="$1"
    if [[ -f "${CHECKPOINT_DIR}/${step_name}.checkpoint" ]]; then
        log "Step ${step_name} already completed. Skipping..."
        return 0
    else
        return 1
    fi
}

# Function to pause and wait for user input
pause_script() {
    local step_name="$1"
    local message="${2:-Press Enter to continue to the next step, or Ctrl+C to exit}"
    
    log "Completed step: ${step_name}"
    echo ""
    echo "=================================================="
    echo "Step completed: ${step_name}"
    echo "Checkpoint saved at: ${CHECKPOINT_DIR}/${step_name}.checkpoint"
    echo "=================================================="
    echo "${message}"
    read -r
}

# Step 1: Environment Setup
step1_environment_setup() {
    if check_checkpoint "step1_environment_setup"; then
        return 0
    fi
    
    log "Starting Step 1: Environment Setup"
    
    # Check if conda is available
    if ! command -v conda &> /dev/null; then
        error_exit "Conda is not installed or not in PATH"
    fi
    
    # Initialize conda for this shell session
    source "$(conda info --base)/etc/profile.d/conda.sh"
    
    # Check if environment already exists
    if conda env list | grep -q "^${ENV_NAME} "; then
        log "QIIME2 environment '${ENV_NAME}' already exists"
        echo ""
        echo "=================================================="
        echo "Existing environment detected: ${ENV_NAME}"
        echo "=================================================="
        echo "1) Use existing environment"
        echo "2) Remove and recreate environment"
        echo "3) Exit and use a different environment name"
        echo ""
        read -p "Enter choice (1-3): " ENV_CHOICE
        
        case ${ENV_CHOICE} in
            1)
                log "Using existing QIIME2 environment: ${ENV_NAME}"
                ;;
            2)
                log "Removing existing environment and recreating..."
                conda env remove --name "${ENV_NAME}" -y || error_exit "Failed to remove existing environment"
                
                log "Creating new QIIME2 environment..."
                if ! conda env create \
                    --name "${ENV_NAME}" \
                    --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml; then
                    
                    log "Initial environment creation failed. Trying with flexible channel priority..."
                    conda config --set channel_priority flexible
                    
                    conda env create \
                        --name "${ENV_NAME}" \
                        --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml || \
                        error_exit "Failed to create QIIME2 environment"
                fi
                ;;
            3)
                log "Exiting. Please set ENV_NAME variable in the script to use a different environment name."
                exit 0
                ;;
            *)
                log "Invalid choice. Using existing environment."
                ;;
        esac
    else
        # Environment doesn't exist, create it
        log "Creating QIIME2 environment..."
        if ! conda env create \
            --name "${ENV_NAME}" \
            --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml; then
            
            log "Initial environment creation failed. Trying with flexible channel priority..."
            conda config --set channel_priority flexible
            
            conda env create \
                --name "${ENV_NAME}" \
                --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml || \
                error_exit "Failed to create QIIME2 environment"
        fi
    fi
    
    log "Activating QIIME2 environment..."
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"
    
    # Verify installation
    qiime --version || error_exit "QIIME2 installation verification failed"
    
    create_checkpoint "step1_environment_setup"
    pause_script "Environment Setup" "Environment is ready. Make sure your manifest.tsv file is properly formatted before continuing."
}
# Step 2: Import Data
step2_import_data() {
    if check_checkpoint "step2_import_data"; then
        return 0
    fi
    
    log "Starting Step 2: Import demultiplexed paired-end FASTQ files"
    
    # Ensure environment is activated
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"
    
    # Check if manifest file exists - look in Data/raw_data first, then current directory
    MANIFEST_PATH=""
    if [[ -f "Data/raw_data/manifest.tsv" ]]; then
        MANIFEST_PATH="Data/raw_data/manifest.tsv"
        log "Using manifest file: Data/raw_data/manifest.tsv"
    elif [[ -f "manifest.tsv" ]]; then
        MANIFEST_PATH="manifest.tsv"
        log "Using manifest file: manifest.tsv (in working directory)"
    else
        error_exit "manifest.tsv file not found. Please create this file in Data/raw_data/ or current directory according to QIIME2 specifications."
    fi
    
    log "Importing paired-end data from ${MANIFEST_PATH}..."
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "${MANIFEST_PATH}" \
        --output-path Data/processed_data/demux-paired-end.qza \
        --input-format PairedEndFastqManifestPhred33V2 || \
        error_exit "Failed to import data"
    
    # Verify output file was created
    if [[ ! -f "Data/processed_data/demux-paired-end.qza" ]]; then
        error_exit "Output file Data/processed_data/demux-paired-end.qza was not created"
    fi
    
    create_checkpoint "step2_import_data"
    pause_script "Data Import" "Data imported successfully. Review the Data/processed_data/demux-paired-end.qza file before proceeding."
}

# Step 3: Visualize Demux Data
step3_visualize_demux() {
    if check_checkpoint "step3_visualize_demux"; then
        return 0
    fi
    
    log "Starting Step 3: Visualize demux data for quality check"
    
    # Ensure environment is activated
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"
    
    # Check if input file exists - look in Data/processed_data first, then current directory
    DEMUX_PATH=""
    if [[ -f "Data/processed_data/demux-paired-end.qza" ]]; then
        DEMUX_PATH="Data/processed_data/demux-paired-end.qza"
    elif [[ -f "demux-paired-end.qza" ]]; then
        DEMUX_PATH="demux-paired-end.qza"
    else
        error_exit "Input file demux-paired-end.qza not found in Data/processed_data/ or current directory. Please run step 2 first."
    fi
    
    log "Creating demux visualization from ${DEMUX_PATH}..."
    qiime demux summarize \
        --i-data "${DEMUX_PATH}" \
        --o-visualization Data/processed_data/demux-paired-end.qzv || \
        error_exit "Failed to create demux visualization"
    
    # Verify output file was created
    if [[ ! -f "Data/processed_data/demux-paired-end.qzv" ]]; then
        error_exit "Output file Data/processed_data/demux-paired-end.qzv was not created"
    fi
    
    create_checkpoint "step3_visualize_demux"
    pause_script "Demux Visualization" "Visualization created at Data/processed_data/demux-paired-end.qzv. Please drag this file to https://view.qiime2.org to check quality plots and determine truncation lengths for the next step."
}

# Step 4: Remove Primers/Adapters (Optional)
step4_remove_primers() {
    if check_checkpoint "step4_remove_primers"; then
        return 0
    fi
    
    log "Starting Step 4: Remove primers/adapters (optional)"
    
    # Ensure environment is activated
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"
    
    echo ""
    echo "=================================================="
    echo "Primer/Adapter Removal"
    echo "=================================================="
    echo ""
    echo "Do you need to remove primers/adapters from your sequences?"
    echo ""
    echo "Common 16S primer sets:"
    echo "1) 515F/806R (V4 region) - Earth Microbiome Project"
    echo "   Forward: GTGCCAGCMGCCGCGGTAA"
    echo "   Reverse: GGACTACHVGGGTWTCTAAT"
    echo ""
    echo "2) Custom primers (you will enter them)"
    echo "3) Skip primer removal (already removed or not needed)"
    echo ""
    read -p "Enter choice (1-3): " PRIMER_CHOICE
    
    case ${PRIMER_CHOICE} in
        1)
            PRIMER_F="GTGCCAGCMGCCGCGGTAA"
            PRIMER_R="GGACTACHVGGGTWTCTAAT"
            log "Using 515F/806R primers (V4 region)"
            ;;
        2)
            echo ""
            read -p "Enter forward primer sequence (5' to 3'): " PRIMER_F
            read -p "Enter reverse primer sequence (5' to 3'): " PRIMER_R
            log "Using custom primers: F=${PRIMER_F}, R=${PRIMER_R}"
            ;;
        3)
            log "Skipping primer removal"
            create_checkpoint "step3b_remove_primers"
            return 0
            ;;
        *)
            log "Invalid choice. Skipping primer removal"
            create_checkpoint "step3b_remove_primers"
            return 0
            ;;
    esac
    
    # Check if input file exists
    DEMUX_PATH=""
    if [[ -f "Data/processed_data/demux-paired-end.qza" ]]; then
        DEMUX_PATH="Data/processed_data/demux-paired-end.qza"
    elif [[ -f "demux-paired-end.qza" ]]; then
        DEMUX_PATH="demux-paired-end.qza"
    else
        error_exit "Input file demux-paired-end.qza not found. Please run step 2 first."
    fi
    
    log "Removing primers using cutadapt..."
    log "Forward primer: ${PRIMER_F}"
    log "Reverse primer: ${PRIMER_R}"
    
    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences "${DEMUX_PATH}" \
        --p-front-f "${PRIMER_F}" \
        --p-front-r "${PRIMER_R}" \
        --p-discard-untrimmed \
        --p-no-indels \
        --o-trimmed-sequences Data/processed_data/demux-trimmed.qza \
        --verbose || error_exit "Primer removal failed"
    
    # Verify output file was created
    if [[ ! -f "Data/processed_data/demux-trimmed.qza" ]]; then
        error_exit "Output file Data/processed_data/demux-trimmed.qza was not created"
    fi
    
    # Create visualization of trimmed data
    log "Creating visualization of trimmed sequences..."
    qiime demux summarize \
        --i-data Data/processed_data/demux-trimmed.qza \
        --o-visualization Data/processed_data/demux-trimmed.qzv || \
        log "Warning: Failed to create trimmed data visualization"
    
    create_checkpoint "step4_remove_primers"
    
    echo ""
    echo "=================================================="
    echo "Primer Removal Complete!"
    echo "=================================================="
    echo ""
    echo "Generated files:"
    echo "- Data/processed_data/demux-trimmed.qza (primer-trimmed sequences)"
    echo "- Data/processed_data/demux-trimmed.qzv (visualization)"
    echo ""
    echo "IMPORTANT: Use demux-trimmed.qza for downstream analysis (Step 5)"
    echo ""
    
    pause_script "Primer Removal" "Primers removed. Review demux-trimmed.qzv to verify primer removal before proceeding."
}

# Step 5: DADA2 Denoising
step5_dada2_denoising() {
    if check_checkpoint "step5_dada2_denoising"; then
        return 0
    fi
    
    log "Starting Step 5: DADA2 denoising and feature table generation"
    
    # Ensure environment is activated
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"
    
    # Check if input file exists - prefer trimmed version if available
    DEMUX_PATH=""
    if [[ -f "Data/processed_data/demux-trimmed.qza" ]]; then
        DEMUX_PATH="Data/processed_data/demux-trimmed.qza"
        log "Using primer-trimmed sequences: ${DEMUX_PATH}"
    elif [[ -f "Data/processed_data/demux-paired-end.qza" ]]; then
        DEMUX_PATH="Data/processed_data/demux-paired-end.qza"
        log "Using original sequences (primers not removed): ${DEMUX_PATH}"
    elif [[ -f "demux-paired-end.qza" ]]; then
        DEMUX_PATH="demux-paired-end.qza"
    else
        error_exit "Input file not found. Please run step 2 (and optionally step 3b) first."
    fi
    
    # Get truncation lengths from user
    echo ""
    echo "Based on the quality plots from Data/processed_data/demux-paired-end.qzv:"
    read -p "Enter forward truncation length (default: 250): " TRUNC_LEN_F
    read -p "Enter reverse truncation length (default: 250): " TRUNC_LEN_R
    
    TRUNC_LEN_F=${TRUNC_LEN_F:-250}
    TRUNC_LEN_R=${TRUNC_LEN_R:-250}
    
    log "Using truncation lengths: Forward=${TRUNC_LEN_F}, Reverse=${TRUNC_LEN_R}"
    
    # Get number of threads to use
    echo ""
    echo "DADA2 can use multiple CPU cores to speed up processing."
    read -p "Enter number of threads (0 = all available, default: 0): " N_THREADS
    N_THREADS=${N_THREADS:-0}
    
    log "Using ${N_THREADS} threads (0 = all available cores)"
    
    # Create Results directory for denoise mode
    mkdir -p Results/denoise_mode
    
    log "Running DADA2 denoising (this may take a while)..."
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "${DEMUX_PATH}" \
        --p-trunc-len-f "${TRUNC_LEN_F}" \
        --p-trunc-len-r "${TRUNC_LEN_R}" \
        --p-n-threads "${N_THREADS}" \
        --o-table Results/denoise_mode/table.qza \
        --o-representative-sequences Results/denoise_mode/rep-seqs.qza \
        --o-denoising-stats Results/denoise_mode/denoising-stats.qza \
        --o-base-transition-stats Results/denoise_mode/base-transition-stats.qza || \
        error_exit "DADA2 denoising failed"
    
    # Verify output files were created
    for file in Results/denoise_mode/table.qza Results/denoise_mode/rep-seqs.qza Results/denoise_mode/denoising-stats.qza Results/denoise_mode/base-transition-stats.qza; do
        if [[ ! -f "${file}" ]]; then
            error_exit "Output file ${file} was not created"
        fi
    done
    
    # Create visualizations of DADA2 results
    log "Creating visualizations of DADA2 results..."
    
    log "Creating denoising stats visualization..."
    qiime metadata tabulate \
        --m-input-file Results/denoise_mode/denoising-stats.qza \
        --o-visualization Results/denoise_mode/denoising-stats.qzv || \
        log "Warning: Failed to create denoising stats visualization"
    
    log "Creating feature table summary..."
    qiime feature-table summarize \
        --i-table Results/denoise_mode/table.qza \
        --o-visualization Results/denoise_mode/table.qzv || \
        log "Warning: Failed to create feature table visualization"
    
    log "Creating representative sequences summary..."
    qiime feature-table tabulate-seqs \
        --i-data Results/denoise_mode/rep-seqs.qza \
        --o-visualization Results/denoise_mode/rep-seqs.qzv || \
        log "Warning: Failed to create rep-seqs visualization"
    
    create_checkpoint "step5_dada2_denoising"
    
    echo ""
    echo "=================================================="
    echo "DADA2 Denoising Complete!"
    echo "=================================================="
    echo ""
    echo "Generated files:"
    echo "1. Results/denoise_mode/table.qza - Feature table (ASV counts)"
    echo "2. Results/denoise_mode/rep-seqs.qza - Representative sequences"
    echo "3. Results/denoise_mode/denoising-stats.qza - Denoising statistics"
    echo ""
    echo "Visualizations (open at https://view.qiime2.org):"
    echo "1. Results/denoise_mode/denoising-stats.qzv - See how many reads passed filters"
    echo "2. Results/denoise_mode/table.qzv - Feature table summary (# features, sampling depth)"
    echo "3. Results/denoise_mode/rep-seqs.qzv - Browse representative sequences"
    echo ""
    
    pause_script "DADA2 Denoising" "DADA2 processing completed. Review the .qzv files at https://view.qiime2.org before proceeding."
}

# Step 5 Alternative: Build OTU table and rep-seqs directly from demux
step5_cluster_from_demux() {
    if check_checkpoint "step5_cluster_from_demux"; then
        return 0
    fi

    log "Starting Step 5 (cluster mode): Generate OTU table and representative sequences from demux"

    # Ensure environment is activated
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"

    # Check if input file exists - prefer trimmed version if available
    DEMUX_PATH=""
    if [[ -f "Data/processed_data/demux-trimmed.qza" ]]; then
        DEMUX_PATH="Data/processed_data/demux-trimmed.qza"
        log "Using primer-trimmed sequences: ${DEMUX_PATH}"
    elif [[ -f "Data/processed_data/demux-paired-end.qza" ]]; then
        DEMUX_PATH="Data/processed_data/demux-paired-end.qza"
        log "Using original sequences (primers not removed): ${DEMUX_PATH}"
    elif [[ -f "demux-paired-end.qza" ]]; then
        DEMUX_PATH="demux-paired-end.qza"
    else
        error_exit "Input file not found. Please run step 2 (and optionally step 3b) first."
    fi

    # Create Results directory for cluster mode
    mkdir -p Results/cluster_mode

    log "Merging paired reads using QIIME2 vsearch..."
    qiime vsearch merge-pairs \
        --i-demultiplexed-seqs "${DEMUX_PATH}" \
        --o-merged-sequences Results/cluster_mode/joined.qza \
        --o-unmerged-sequences Results/cluster_mode/unmerged.qza || error_exit "merge-pairs failed"

    log "Quality-filtering merged reads..."
    qiime quality-filter q-score \
        --i-demux Results/cluster_mode/joined.qza \
        --o-filtered-sequences Results/cluster_mode/filtered-seqs.qza \
        --o-filter-stats Results/cluster_mode/filtered-stats.qza || error_exit "quality-filter failed"

    log "Dereplicating sequences..."
    qiime vsearch dereplicate-sequences \
        --i-sequences Results/cluster_mode/filtered-seqs.qza \
        --o-dereplicated-table Results/cluster_mode/derep-table.qza \
        --o-dereplicated-sequences Results/cluster_mode/derep-seqs.qza || error_exit "dereplicate-sequences failed"
    
    log "Removing chimeras from dereplicated sequences..."
    qiime vsearch uchime-denovo \
        --i-table Results/cluster_mode/derep-table.qza \
        --i-sequences Results/cluster_mode/derep-seqs.qza \
        --o-chimeras Results/cluster_mode/chimeras.qza \
        --o-nonchimeras Results/cluster_mode/rep-seqs.qza \
        --o-stats Results/cluster_mode/chimera-stats.qza || error_exit "chimera removal failed"
    
    log "Filtering feature table to remove chimeric sequences..."
    qiime feature-table filter-features \
        --i-table Results/cluster_mode/derep-table.qza \
        --m-metadata-file Results/cluster_mode/rep-seqs.qza \
        --o-filtered-table Results/cluster_mode/table.qza || error_exit "table filtering failed"

    # At this point we have table.qza and rep-seqs.qza; offer clustering options
    echo ""
    echo "Choose clustering method to produce OTUs (from dereplicated sequences):"
    echo "1) de novo clustering (vsearch cluster-features-de-novo, 97%)"
    echo "2) closed-reference clustering (requires reference database)"
    echo "3) skip clustering (keep dereplicated sequences as features)"
    read -p "Enter choice (1-3): " CLUSTER_CHOICE

    case ${CLUSTER_CHOICE} in
        1)
            log "Running de novo clustering (97%)..."
            qiime vsearch cluster-features-de-novo \
                --i-table Results/cluster_mode/table.qza \
                --i-sequences Results/cluster_mode/rep-seqs.qza \
                --p-perc-identity 0.97 \
                --o-clustered-table Results/cluster_mode/table-dn-97.qza \
                --o-clustered-sequences Results/cluster_mode/rep-seqs-dn-97.qza || error_exit "de-novo clustering failed"
            # Promote outputs to standard names for downstream steps
            mv Results/cluster_mode/table-dn-97.qza Results/cluster_mode/table.qza
            mv Results/cluster_mode/rep-seqs-dn-97.qza Results/cluster_mode/rep-seqs.qza
            ;;
        2)
            if [[ ! -f "Data/reference_dbs/silva_97_otus.qza" ]] && [[ ! -f "silva_97_otus.qza" ]]; then
                log "Reference database not found in Data/reference_dbs/ or current directory. Cannot run closed-reference clustering."
            else
                REF_DB_PATH=""
                if [[ -f "Data/reference_dbs/silva_97_otus.qza" ]]; then
                    REF_DB_PATH="Data/reference_dbs/silva_97_otus.qza"
                else
                    REF_DB_PATH="silva_97_otus.qza"
                fi
                log "Running closed-reference clustering with ${REF_DB_PATH}..."
                qiime vsearch cluster-features-closed-reference \
                    --i-table Results/cluster_mode/table.qza \
                    --i-sequences Results/cluster_mode/rep-seqs.qza \
                    --i-reference-sequences "${REF_DB_PATH}" \
                    --p-perc-identity 0.97 \
                    --o-clustered-table Results/cluster_mode/table-cr-97.qza \
                    --o-clustered-sequences Results/cluster_mode/rep-seqs-cr-97.qza \
                    --o-unmatched-sequences Results/cluster_mode/unmatched.qza || error_exit "Closed-reference clustering failed"
                mv Results/cluster_mode/table-cr-97.qza Results/cluster_mode/table.qza
                mv Results/cluster_mode/rep-seqs-cr-97.qza Results/cluster_mode/rep-seqs.qza
            fi
            ;;
        3)
            log "Skipping additional clustering; using dereplicated features"
            ;;
        *)
            log "Invalid choice. Using dereplicated features"
            ;;
    esac

    # Verify outputs exist
    for file in Results/cluster_mode/table.qza Results/cluster_mode/rep-seqs.qza; do
        if [[ ! -f "${file}" ]]; then
            error_exit "Expected output ${file} was not created"
        fi
    done

    create_checkpoint "step5_cluster_from_demux"
    pause_script "Cluster-from-demux" "Clustering-from-demux completed. Outputs in Results/cluster_mode/: table.qza, rep-seqs.qza"
}

# Step 6: Taxonomic Classification
step6_taxonomic_classification() {
    if check_checkpoint "step6_taxonomic_classification"; then
        return 0
    fi
    
    log "Starting Step 6: Taxonomic Classification"
    
    # Ensure environment is activated
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"
    
    # Determine which mode we're in and set paths
    if [[ "${MODE}" == "denoise" ]]; then
        TABLE_PATH="Results/denoise_mode/table.qza"
        REP_SEQS_PATH="Results/denoise_mode/rep-seqs.qza"
        OUTPUT_DIR="Results/denoise_mode"
    elif [[ "${MODE}" == "cluster" ]]; then
        TABLE_PATH="Results/cluster_mode/table.qza"
        REP_SEQS_PATH="Results/cluster_mode/rep-seqs.qza"
        OUTPUT_DIR="Results/cluster_mode"
    else
        error_exit "Unknown mode: ${MODE}"
    fi
    
    # Check if required files exist
    if [[ ! -f "${TABLE_PATH}" ]] || [[ ! -f "${REP_SEQS_PATH}" ]]; then
        error_exit "Required files not found. Please run step 5 first."
    fi
    
    # Check for classifier
    CLASSIFIER_PATH=""
    if [[ -f "Data/reference_dbs/classifier.qza" ]]; then
        CLASSIFIER_PATH="Data/reference_dbs/classifier.qza"
    elif [[ -f "classifier.qza" ]]; then
        CLASSIFIER_PATH="classifier.qza"
    fi
    
    if [[ -z "${CLASSIFIER_PATH}" ]]; then
        echo ""
        echo "=================================================="
        echo "Taxonomic Classifier Not Found"
        echo "=================================================="
        echo ""
        echo "Choose an option:"
        echo "1) Use an existing classifier file (you provide the path)"
        echo "2) Build a custom classifier (requires reference sequences)"
        echo "3) Download a pre-trained classifier"
        echo "4) Skip taxonomic classification"
        echo ""
        read -p "Enter choice (1-4): " CLASSIFIER_CHOICE
        
        case ${CLASSIFIER_CHOICE} in
            1)
                read -p "Enter the full path to your classifier file: " USER_CLASSIFIER
                if [[ -f "${USER_CLASSIFIER}" ]]; then
                    CLASSIFIER_PATH="${USER_CLASSIFIER}"
                    log "Using classifier: ${CLASSIFIER_PATH}"
                else
                    error_exit "Classifier file not found: ${USER_CLASSIFIER}"
                fi
                ;;
            2)
                echo ""
                echo "To build a custom classifier, you need:"
                echo "1. Reference sequences (e.g., silva-138-99-seqs.qza)"
                echo "2. Taxonomy file (e.g., silva-138-99-tax.qza)"
                echo ""
                echo "You can download these from:"
                echo "https://docs.qiime2.org/2024.10/data-resources/"
                echo ""
                
                read -p "Enter path to reference sequences (.qza): " REF_SEQS
                read -p "Enter path to taxonomy file (.qza): " REF_TAX
                
                if [[ ! -f "${REF_SEQS}" ]]; then
                    error_exit "Reference sequences not found: ${REF_SEQS}"
                fi
                if [[ ! -f "${REF_TAX}" ]]; then
                    error_exit "Taxonomy file not found: ${REF_TAX}"
                fi
                
                echo ""
                echo "Do you want to extract a specific region using primers?"
                echo "This is recommended if you used specific primers (e.g., 515F/806R)"
                read -p "Extract primers? (y/n): " EXTRACT_PRIMERS
                
                if [[ "${EXTRACT_PRIMERS}" =~ ^[Yy]$ ]]; then
                    read -p "Enter forward primer sequence: " PRIMER_F
                    read -p "Enter reverse primer sequence: " PRIMER_R
                    
                    log "Extracting amplicon region from reference sequences..."
                    qiime feature-classifier extract-reads \
                        --i-sequences "${REF_SEQS}" \
                        --p-f-primer "${PRIMER_F}" \
                        --p-r-primer "${PRIMER_R}" \
                        --o-reads Data/reference_dbs/ref-seqs-extracted.qza || \
                        error_exit "Failed to extract reads"
                    
                    REF_SEQS="Data/reference_dbs/ref-seqs-extracted.qza"
                    log "Using extracted sequences: ${REF_SEQS}"
                fi
                
                log "Training classifier (this may take 30+ minutes)..."
                mkdir -p Data/reference_dbs
                qiime feature-classifier fit-classifier-naive-bayes \
                    --i-reference-reads "${REF_SEQS}" \
                    --i-reference-taxonomy "${REF_TAX}" \
                    --o-classifier Data/reference_dbs/custom-classifier.qza || \
                    error_exit "Classifier training failed"
                
                CLASSIFIER_PATH="Data/reference_dbs/custom-classifier.qza"
                log "Custom classifier created: ${CLASSIFIER_PATH}"
                ;;
            3)
                echo ""
                echo "To perform taxonomic classification, you need a pre-trained classifier."
                echo ""
                echo "You can:"
                echo "1. Download a pre-trained classifier from:"
                echo "   https://docs.qiime2.org/2024.10/data-resources/"
                echo ""
                echo "2. For 16S data, common classifiers:"
                echo "   - Silva 138 99% OTUs"
                echo "   - Greengenes2 2022.10"
                echo "   - GTDB"
                echo ""
                echo "3. Place the classifier in Data/reference_dbs/classifier.qza"
                echo ""
                echo "Example download command:"
                echo "  wget -O Data/reference_dbs/classifier.qza \\"
                echo "    https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"
                echo ""
                log "Skipping taxonomic classification. Please download a classifier and run step 6 again."
                create_checkpoint "step6_taxonomic_classification"
                return 0
                ;;
            4)
                log "Skipping taxonomic classification."
                create_checkpoint "step6_taxonomic_classification"
                return 0
                ;;
            *)
                log "Invalid choice. Skipping taxonomic classification."
                create_checkpoint "step6_taxonomic_classification"
                return 0
                ;;
        esac
    fi
    
    log "Using classifier: ${CLASSIFIER_PATH}"
    
    log "Running taxonomic classification (this may take a while)..."
    qiime feature-classifier classify-sklearn \
        --i-classifier "${CLASSIFIER_PATH}" \
        --i-reads "${REP_SEQS_PATH}" \
        --o-classification "${OUTPUT_DIR}/taxonomy.qza" || \
        error_exit "Taxonomic classification failed"
    
    log "Creating taxonomy visualization..."
    qiime metadata tabulate \
        --m-input-file "${OUTPUT_DIR}/taxonomy.qza" \
        --o-visualization "${OUTPUT_DIR}/taxonomy.qzv" || \
        error_exit "Failed to create taxonomy visualization"
    
    log "Creating taxonomic bar plots..."
    qiime taxa barplot \
        --i-table "${TABLE_PATH}" \
        --i-taxonomy "${OUTPUT_DIR}/taxonomy.qza" \
        --o-visualization "${OUTPUT_DIR}/taxa-bar-plots.qzv" || \
        error_exit "Failed to create taxonomic bar plots"
    
    # Optional: Export taxonomy table to TSV
    log "Exporting taxonomy table to TSV..."
    qiime tools export \
        --input-path "${OUTPUT_DIR}/taxonomy.qza" \
        --output-path "${OUTPUT_DIR}/exported-taxonomy" || \
        log "Warning: Failed to export taxonomy table"
    
    create_checkpoint "step6_taxonomic_classification"
    
    echo ""
    echo "=================================================="
    echo "Taxonomic Classification Complete!"
    echo "=================================================="
    echo ""
    echo "Generated files:"
    echo "1. ${OUTPUT_DIR}/taxonomy.qza"
    echo "   - Taxonomic assignments for each feature"
    echo ""
    echo "2. ${OUTPUT_DIR}/taxonomy.qzv"
    echo "   - Taxonomy table visualization"
    echo ""
    echo "3. ${OUTPUT_DIR}/taxa-bar-plots.qzv"
    echo "   - Interactive taxonomic bar plots (like your example image!)"
    echo "   - View at https://view.qiime2.org"
    echo ""
    if [[ -f "${OUTPUT_DIR}/exported-taxonomy/taxonomy.tsv" ]]; then
        echo "4. ${OUTPUT_DIR}/exported-taxonomy/taxonomy.tsv"
        echo "   - Taxonomy table in TSV format for further analysis"
        echo ""
    fi
    echo "The taxa-bar-plots.qzv file will show relative abundances"
    echo "of different taxa across your samples, just like your example!"
    echo ""
    
    pause_script "Taxonomic Classification" "Taxonomic classification complete. Open taxa-bar-plots.qzv to see your taxonomic composition!"
}

# Main execution function
main() {
    log "Starting QIIME2 Pipeline with Checkpoint System"
    log "Script directory: ${SCRIPT_DIR}"
    log "Project directory: ${PROJECT_DIR}"
    log "Checkpoint directory: ${CHECKPOINT_DIR}"
    log "Log file: ${LOG_FILE}"
    
    echo ""
    echo "=============================================="
    echo "QIIME2 Pipeline with Checkpoint System"
    echo "=============================================="
    echo ""
    
    # Ask user for environment name if not already set via command line
    if [[ -z "${ENV_NAME_SET:-}" ]]; then
        echo "Current conda environment name: ${ENV_NAME}"
        read -p "Enter your QIIME2 environment name (or press Enter to use '${ENV_NAME}'): " USER_ENV_NAME
        if [[ -n "${USER_ENV_NAME}" ]]; then
            ENV_NAME="${USER_ENV_NAME}"
            log "Using user-specified environment name: ${ENV_NAME}"
        else
            log "Using default environment name: ${ENV_NAME}"
        fi
        echo ""
    fi
    
    echo "This script will run the following steps:"
    echo "1. Environment Setup"
    echo "2. Import Data"
    echo "3. Visualize Demux Data"
    echo "4. Remove Primers/Adapters (optional)"
    echo "5. Denoising or Clustering (mode: ${MODE})"
    echo "6. Taxonomic Classification"
    echo ""
    echo "You can:"
    echo "- Run all steps: Press Enter"
    echo "- Run from specific step: Type step number (1-6)"
    echo "- Exit: Ctrl+C"
    echo ""
    if [[ -z "${START_STEP:-}" ]]; then
        read -p "Enter step number to start from (or press Enter for all steps): " START_STEP
        START_STEP=${START_STEP:-1}
    else
        log "Starting from step ${START_STEP} (provided via CLI)"
    fi

    # Validate start step
    if ! [[ "${START_STEP}" =~ ^[1-6]$ ]]; then
        error_exit "Invalid step number. Please enter 1-6."
    fi
    
    # Run steps based on start step
    if [[ ${START_STEP} -le 1 ]]; then
        step1_environment_setup
    fi
    
    if [[ ${START_STEP} -le 2 ]]; then
        step2_import_data
    fi
    
    if [[ ${START_STEP} -le 3 ]]; then
        step3_visualize_demux
    fi
    
    if [[ ${START_STEP} -le 4 ]]; then
        step4_remove_primers
    fi
    
    if [[ ${START_STEP} -le 5 ]]; then
        if [[ "${MODE}" == "denoise" ]]; then
            step5_dada2_denoising
        elif [[ "${MODE}" == "cluster" ]]; then
            step5_cluster_from_demux
        else
            error_exit "Unknown mode: ${MODE}. Supported: denoise, cluster"
        fi
    fi
    
    if [[ ${START_STEP} -le 6 ]]; then
        step6_taxonomic_classification
    fi
    
    log "Pipeline completed successfully!"
    echo ""
    echo "=============================================="
    echo "Pipeline completed successfully!"
    echo "=============================================="
    echo ""
    echo "Generated files:"
    echo "- Data/processed_data/demux-paired-end.qza (imported data)"
    echo "- Data/processed_data/demux-paired-end.qzv (quality visualization)"
    if [[ "${MODE}" == "denoise" ]]; then
        echo "- Results/denoise_mode/table.qza (feature table)"
        echo "- Results/denoise_mode/rep-seqs.qza (representative sequences)"
        echo "- Results/denoise_mode/denoising-stats.qza (DADA2 statistics)"
    elif [[ "${MODE}" == "cluster" ]]; then
        echo "- Results/cluster_mode/table.qza (feature table)"
        echo "- Results/cluster_mode/rep-seqs.qza (representative sequences)"
        echo "- Results/cluster_mode/joined.qza, filtered-seqs.qza (intermediate files)"
    fi
    if [[ -f "Results/${MODE}_mode/table-cr-97.qza" || -f "Results/${MODE}_mode/table-or-97.qza" ]]; then
        echo "- Results/${MODE}_mode/table-*-97.qza (additional clustered tables, if VSEARCH was used)"
    fi
    echo ""
    echo "Checkpoints saved in: ${CHECKPOINT_DIR}"
    echo "Log file: ${LOG_FILE}"
    echo ""
    echo "Next steps:"
    echo "1. Review the quality plots in Data/processed_data/demux-paired-end.qzv"
    echo "2. Consider running taxonomic classification"
    echo "3. Generate phylogenetic tree"
    echo "4. Perform diversity analysis"
    echo ""
}

# Cleanup function for graceful exit
cleanup() {
    log "Script interrupted. Checkpoints preserved for resuming."
    exit 130
}

# Set up signal handlers
trap cleanup SIGINT SIGTERM

# Usage function
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -h, --help         Show this help message"
    echo "  -c, --clean        Remove all checkpoints and start fresh"
    echo "  -s, --status       Show status of completed steps"
    echo "  -r, --remove       Remove a specific checkpoint"
    echo "                     Example: -r step5_dada2_denoising"
    echo "  -d, --delete-intermediate"
    echo "                     Delete intermediate outputs from step 4 and 5"
    echo "  -m, --mode         Pipeline mode: 'denoise' (DADA2) or 'cluster' (generate OTUs from demux). Default: denoise"
    echo "                     Example: -m cluster"
    echo "  -e, --env          Conda environment name to use. Default: qiime"
    echo "                     Example: -e qiime2-2024.5"
    echo ""
    echo "Examples:"
    echo "  $0                 Run the pipeline (resume from last checkpoint)"
    echo "  $0 --clean         Start the pipeline from scratch"
    echo "  $0 --status        Check which steps have been completed"
    echo "  $0 -r step5_dada2_denoising  Remove checkpoint for step 5"
    echo "  $0 -d              Delete intermediate outputs and re-run step 4/5"
    echo "  $0 -e myqiime      Use conda environment named 'myqiime'"
    echo "  $0 -m cluster -e qiime2-amplicon  Use cluster mode with custom environment"
    echo ""
    echo "Available checkpoints to remove:"
    echo "  - step1_environment_setup"
    echo "  - step2_import_data"
    echo "  - step3_visualize_demux"
    echo "  - step4_remove_primers"
    echo "  - step5_dada2_denoising"
    echo "  - step5_cluster_from_demux"
    echo "  - step6_taxonomic_classification"
}

# Status function
show_status() {
    echo "Pipeline Status:"
    echo "================"
    
    local steps=("step1_environment_setup" "step2_import_data" "step3_visualize_demux" "step4_remove_primers" "step5_dada2_denoising" "step6_taxonomic_classification")
    local step_names=("Step 1: Environment Setup" "Step 2: Import Data" "Step 3: Visualize Demux" "Step 4: Remove Primers" "Step 5: DADA2 Denoising/Clustering" "Step 6: Taxonomic Classification")
    
    for i in "${!steps[@]}"; do
        local step="${steps[$i]}"
        local name="${step_names[$i]}"
        if [[ -f "${CHECKPOINT_DIR}/${step}.checkpoint" ]]; then
            local timestamp=$(cat "${CHECKPOINT_DIR}/${step}.checkpoint")
            echo "✓ ${name} (completed: ${timestamp})"
        else
            echo "✗ ${name} (not completed)"
        fi
    done
    echo ""
}

# Clean function
clean_checkpoints() {
    echo "Removing all checkpoints..."
    rm -rf "${CHECKPOINT_DIR}"
    rm -f "${LOG_FILE}"
    echo "Checkpoints and log file removed. Pipeline will start from the beginning."
}

# Remove specific checkpoint function
remove_checkpoint() {
    local step_name="$1"
    local checkpoint_file="${CHECKPOINT_DIR}/${step_name}.checkpoint"
    
    if [[ -f "${checkpoint_file}" ]]; then
        rm "${checkpoint_file}"
        echo "Checkpoint removed: ${step_name}"
        log "Checkpoint manually removed: ${step_name}"
    else
        echo "Checkpoint not found: ${step_name}"
        echo "Available checkpoints:"
        if [[ -d "${CHECKPOINT_DIR}" ]]; then
            for checkpoint in "${CHECKPOINT_DIR}"/*.checkpoint; do
                if [[ -f "${checkpoint}" ]]; then
                    basename "${checkpoint}" .checkpoint
                fi
            done
        else
            echo "  (none)"
        fi
    fi
}

# Remove intermediate outputs function
remove_intermediate_outputs() {
    echo ""
    echo "=================================================="
    echo "Remove Intermediate Outputs"
    echo "=================================================="
    echo ""
    echo "This will delete intermediate files from:"
    echo ""
    echo "Step 4 (Remove Primers):"
    echo "  - demux-trimmed.qza"
    echo "  - demux-trimmed.qzv"
    echo ""
    echo "Step 5 (DADA2 Denoising):"
    echo "  - table.qza"
    echo "  - rep-seqs.qza"
    echo "  - denoising-stats.qza"
    echo "  - base-transition-stats.qza"
    echo "  - rep-seqs.qzv"
    echo ""
    echo "And reset checkpoints for these steps."
    echo ""
    read -p "Are you sure you want to delete these files? (y/n): " CONFIRM_DELETE
    
    if [[ "${CONFIRM_DELETE}" =~ ^[Yy]$ ]]; then
        echo "Deleting intermediate outputs..."
        
        # Remove Step 4 files
        if [[ -f "demux-trimmed.qza" ]]; then
            rm "demux-trimmed.qza"
            echo "✓ Deleted: demux-trimmed.qza"
        fi
        if [[ -f "demux-trimmed.qzv" ]]; then
            rm "demux-trimmed.qzv"
            echo "✓ Deleted: demux-trimmed.qzv"
        fi
        
        # Remove Step 5 files
        if [[ -f "table.qza" ]]; then
            rm "table.qza"
            echo "✓ Deleted: table.qza"
        fi
        if [[ -f "rep-seqs.qza" ]]; then
            rm "rep-seqs.qza"
            echo "✓ Deleted: rep-seqs.qza"
        fi
        if [[ -f "denoising-stats.qza" ]]; then
            rm "denoising-stats.qza"
            echo "✓ Deleted: denoising-stats.qza"
        fi
        if [[ -f "base-transition-stats.qza" ]]; then
            rm "base-transition-stats.qza"
            echo "✓ Deleted: base-transition-stats.qza"
        fi
        if [[ -f "rep-seqs.qzv" ]]; then
            rm "rep-seqs.qzv"
            echo "✓ Deleted: rep-seqs.qzv"
        fi
        
        # Remove checkpoints
        echo "Removing checkpoints..."
        remove_checkpoint "step4_remove_primers"
        remove_checkpoint "step5_dada2_denoising"
        
        echo ""
        echo "Intermediate outputs removed successfully."
        echo "You can now re-run these steps with different parameters."
        log "Intermediate outputs from step4 and step5 removed by user"
    else
        echo "Cancelled. No files were deleted."
    fi
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        -c|--clean)
            clean_checkpoints
            exit 0
            ;;
        -s|--status)
            show_status
            exit 0
            ;;
        -d|--delete-intermediate)
            remove_intermediate_outputs
            exit 0
            ;;
        -r|--remove)
            if [[ -n "${2:-}" ]]; then
                remove_checkpoint "$2"
                exit 0
            else
                echo "Error: --remove requires an argument (checkpoint name)"
                echo "Example: $0 --remove step4_dada2_denoising"
                exit 1
            fi
            ;;
        -m|--mode)
            if [[ -n "${2:-}" ]]; then
                MODE="$2"
                shift 2
            else
                echo "Error: --mode requires an argument (denoise|cluster)"
                exit 1
            fi
            ;;
        -e|--env)
            if [[ -n "${2:-}" ]]; then
                ENV_NAME="$2"
                ENV_NAME_SET="true"
                shift 2
            else
                echo "Error: --env requires an argument (environment name)"
                exit 1
            fi
            ;;
        # positional: start step (1-6)
        [1-6])
            START_STEP="$1"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

main

