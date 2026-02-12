#!/bin/bash

################################################################################
# Quick Metagenomic Search Pipeline - Methane Gene Detection
################################################################################
# 
# PURPOSE:
#   Fast detection and quantification of methane metabolism genes in shotgun
#   metagenomic data. Identifies methanotrophs (methane oxidizers) and 
#   methanogens (methane producers).
#
# WORKFLOW:
#   1. BBduk - Quality control and adapter trimming
#   2. SingleM - Taxonomic profiling using 16S marker genes
#   3. DIAMOND - Functional gene search (pmoA, mmoX, mcrA, etc.)
#   4. Quantification - Gene abundance and RPKM normalization
#
# REQUIREMENTS:
#   • Conda environment: quick_search
#     - BBMap (bbduk.sh)
#     - SingleM (with hmmer, orfm dependencies)
#     - DIAMOND
#   
#   • Databases (auto-downloaded or manual):
#     - SingleM metapackage: Data/reference_dbs/S5.4.0.GTDB_r226.metapackage_*.smpkg.zb
#     - DIAMOND database: Data/reference_dbs/DIAMOND/methane_master_db.dmnd
#
#   • Input data:
#     - Paired-end shotgun reads: Data/raw_data/shotgun/<SAMPLE>_R1*.fastq.gz
#                                                        <SAMPLE>_R2*.fastq.gz
#
# SETUP:
#   1. Create conda environment:
#      conda create -n quick_search -c bioconda -c conda-forge \
#          python=3.10 singlem bbmap diamond
#   
#   2. Download SingleM database (15-20GB):
#      conda activate quick_search
#      singlem data --output-directory Data/reference_dbs/
#   
#   3. Copy DIAMOND database to Data/reference_dbs/DIAMOND/
#
# USAGE EXAMPLES:
#   # Activate environment first
#   conda activate quick_search
#
#   # Single sample
#   bash Scripts/01_pipelines/quick_search.sh 53394
#
#   # Multiple samples
#   bash Scripts/01_pipelines/quick_search.sh 53394 53395 53396
#
#   # From sample list file
#   bash Scripts/01_pipelines/quick_search.sh --sample-list Config/samples.txt
#
#   # Auto-detect all samples
#   bash Scripts/01_pipelines/quick_search.sh --auto
#
#   # Custom threads (default: 4)
#   bash Scripts/01_pipelines/quick_search.sh --threads 16 53394
#
# HPC USAGE:
#   sbatch Scripts/02_hpc/submit_quick_search.sbatch 53394
#
# OUTPUT:
#   Data/processed_data/bbduk_cleaned/         - Quality-controlled reads
#   Data/processed_data/singlem_output/        - Taxonomic profiles
#   Data/functional_analysis/methane_genes/    - Gene search results
#   Logs/                                      - Sample-specific logs
#
# AUTHOR: Jiayi
# DATE: 2026-02-01
################################################################################

# ============================================================================
# Environment Check and Configuration
# ============================================================================

# Check if conda environment is activated
if [[ "$CONDA_DEFAULT_ENV" != "quick_search" ]] && [[ -z "$SLURM_JOB_ID" ]]; then
    echo "========================================================================"
    echo "WARNING: Conda environment 'quick_search' is not activated"
    echo "======================================================================="
    echo ""
    echo "Please activate the environment first:"
    echo "  conda activate quick_search"
    echo ""
    echo "If you don't have the environment yet, create it:"
    echo "  conda create -n quick_search -c bioconda -c conda-forge \\"
    echo "      python=3.10 singlem bbmap diamond"
    echo ""
    exit 1
fi

# Default configuration
THREADS=${THREADS:-4}              # Number of threads (can be set via environment variable)
MAX_PARALLEL_SAMPLES=${MAX_PARALLEL_SAMPLES:-1}  # Number of samples to process in parallel
CONDA_BASE="${CONDA_BASE:-/opt/anaconda3}"  # Conda installation path
RAW_DATA_SUBDIR="${RAW_DATA_SUBDIR:-shotgun}"  # Subdirectory for raw data (shotgun or empty for root)

# Auto-detect base directory if not set
if [ -z "$BASE_DIR" ]; then
    # Try to find the project root by looking for Scripts directory
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    BASE_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
fi

# Parse command line arguments
SAMPLE_LIST=()
SAMPLE_FILE=""
AUTO_DETECT=false

while [[ $# -gt 0 ]]; do
    case $1 in
        step0|--setup)
              echo "==== 环境与依赖包准备（可自定义环境名） ===="
              echo "1. 创建环境（任选名称）："
              echo "   conda create -n <your_env_name> -c bioconda -c conda-forge \\"
              echo "      python=3.10 singlem bbmap diamond"
             echo "2. Activate environment:"
              echo "   conda activate <your_env_name>"
             echo "3. Download SingleM reference database:"
              echo "   singlem data --output-directory Data/reference_dbs/"
             echo "4. Place your DIAMOND database at: Data/reference_dbs/DIAMOND/"
              echo ""
             echo "After these steps, you can run this script."
              exit 0
              ;;
        --sample-list)
            SAMPLE_FILE="$2"
            shift 2
            ;;
        --auto)
            AUTO_DETECT=true
            shift
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --parallel)
            MAX_PARALLEL_SAMPLES="$2"
            shift 2
            ;;
        --base-dir)
            BASE_DIR="$2"
            shift 2
            ;;
        --help|-h)
            cat << 'EOF'
========================================================================
  Quick Metagenomic Search - Methane Gene Detection Pipeline
========================================================================

USAGE:
  bash Scripts/01_pipelines/quick_search.sh [OPTIONS] [SAMPLE_IDs...]

OPTIONS:

    step0, --setup        Show recommended conda/package setup (custom env name allowed)
    --sample-list FILE    Read sample IDs from file (one per line)
    --auto                Auto-detect all samples in raw_data directory
    --threads N           Number of threads (default: 4)
    --base-dir PATH       Base directory path (default: auto-detect)
    --help, -h            Show this help message

EXAMPLES:
  # Single sample
  bash Scripts/01_pipelines/quick_search.sh 53394

  # Multiple samples
  bash Scripts/01_pipelines/quick_search.sh 53394 53395 53396

  # From sample list file
  bash Scripts/01_pipelines/quick_search.sh --sample-list Config/samples.txt

  # Auto-detect all samples
  bash Scripts/01_pipelines/quick_search.sh --auto

  # Custom thread count
  bash Scripts/01_pipelines/quick_search.sh --threads 16 53394

  # HPC submission
  sbatch Scripts/02_hpc/submit_quick_search.sbatch 53394

REQUIREMENTS:
  1. Conda environment: quick_search (must be activated)
  2. Input data: Data/raw_data/shotgun/<SAMPLE>_R1*.fastq.gz
  3. Databases:
     - SingleM: Data/reference_dbs/S5.4.0.GTDB_r226.metapackage_*.smpkg.zb
     - DIAMOND: Data/reference_dbs/DIAMOND/methane_master_db.dmnd

SETUP (first time):
  # Create environment
  conda create -n quick_search -c bioconda -c conda-forge \\
      python=3.10 singlem bbmap diamond

  # Activate environment
  conda activate quick_search

  # Download SingleM database
  singlem data --output-directory Data/reference_dbs/

OUTPUT:
  Data/processed_data/bbduk_cleaned/      - Quality-controlled reads
  Data/processed_data/singlem_output/     - Taxonomic profiles
  Data/functional_analysis/methane_genes/ - Gene detection results
  Logs/                                   - Execution logs

For more information, see: QUICK_SEARCH_USAGE.md
========================================================================
EOF
            exit 0
            ;;
        -*)
            echo "ERROR: Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
        *)
            SAMPLE_LIST+=("$1")
            shift
            ;;
    esac
done

# Determine samples to process
if [ "$AUTO_DETECT" = true ]; then
    echo "Auto-detecting samples from Data/raw_data/${RAW_DATA_SUBDIR}..."
    RAW_DATA_DIR="${BASE_DIR}/Data/raw_data/${RAW_DATA_SUBDIR}"
    if [ -d "$RAW_DATA_DIR" ]; then
        # Find all unique sample IDs by looking for _R1.fastq.gz files
        while IFS= read -r file; do
            SAMPLE_ID=$(basename "$file" | sed 's/_R1.fastq.gz//')
            SAMPLE_LIST+=("$SAMPLE_ID")
        done < <(find "$RAW_DATA_DIR" -name "*_R1.fastq.gz" -type f)
    fi
elif [ -n "$SAMPLE_FILE" ]; then
    echo "Reading samples from sample list file..."
    if [ ! -f "$SAMPLE_FILE" ]; then
        echo "ERROR: Sample file not found"
        echo "  Specified file: $(basename "$SAMPLE_FILE")"
        echo "  Please check:"
        echo "    1. File path is correct"
        echo "    2. File exists and is readable"
        exit 1
    fi
    while IFS= read -r line || [ -n "$line" ]; do
        # Skip empty lines and comments
        [[ -z "$line" || "$line" =~ ^#.*$ ]] && continue
        SAMPLE_LIST+=("$line")
    done < "$SAMPLE_FILE"
fi

# Check if any samples provided
if [ ${#SAMPLE_LIST[@]} -eq 0 ]; then
    echo "========================================================================"
    echo "ERROR: No samples provided"
    echo "======================================================================="
    echo ""
    echo "Usage: bash Scripts/01_pipelines/quick_search.sh [OPTIONS] [SAMPLE_IDs...]"
    echo ""
    echo "Examples:"
    echo "  bash Scripts/01_pipelines/quick_search.sh 53394"
    echo "  bash Scripts/01_pipelines/quick_search.sh --auto"
    echo "  bash Scripts/01_pipelines/quick_search.sh --sample-list Config/samples.txt"
    echo ""
    echo "Use --help for more information"
    exit 1
fi

echo "========================================================================"
echo "Quick Metagenomic Search Pipeline"
echo "========================================================================"
echo "Samples to process: ${#SAMPLE_LIST[@]}"
echo "Threads: ${THREADS}"
echo "Base directory: ${BASE_DIR}"
echo "Conda environment: ${CONDA_DEFAULT_ENV:-not detected}"
echo "======================================================================="
echo ""

# ============================================================================
# Prerequisites Check
# ============================================================================

echo "Checking prerequisites..."
echo ""

# Check required tools
MISSING_TOOLS=()

if ! command -v bbduk.sh &> /dev/null; then
    MISSING_TOOLS+=("bbduk.sh (BBMap)")
fi

if ! command -v singlem &> /dev/null; then
    MISSING_TOOLS+=("singlem")
fi

if ! command -v diamond &> /dev/null; then
    MISSING_TOOLS+=("diamond")
fi

if ! command -v orfm &> /dev/null; then
    MISSING_TOOLS+=("orfm (orfM)")
fi

if ! command -v hmmsearch &> /dev/null; then
    MISSING_TOOLS+=("hmmsearch (HMMER)")
fi


if ! command -v bc &> /dev/null; then
    MISSING_TOOLS+=("bc")
fi

if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
    echo "========================================================================"
    echo "ERROR: Missing required tools"
    echo "========================================================================"
    echo ""
    echo "The following tools are not found:"
    for tool in "${MISSING_TOOLS[@]}"; do
        echo "  - $tool"
    done
    echo ""
    echo "Please install them in the quick_search environment:"
    echo "  conda activate quick_search"
    echo "  conda install -c bioconda -c conda-forge bbmap singlem diamond orfm hmmer"
    echo "  # or use mamba for faster installs:"
    echo "  # mamba install -n quick_search -c bioconda -c conda-forge bbmap singlem diamond orfm hmmer"
    echo ""
    exit 1
fi

echo "  ✓ bbduk.sh found: $(which bbduk.sh)"
echo "  ✓ singlem found: $(singlem --version 2>&1 | head -1)"
echo "  ✓ diamond found: $(diamond --version 2>&1 | head -1)"
echo ""

# Set paths
RAW_DATA_DIR="${BASE_DIR}/Data/raw_data/${RAW_DATA_SUBDIR}"
BBDUK_DIR="${BASE_DIR}/Data/processed_data/bbduk_cleaned"
SINGLEM_DIR="${BASE_DIR}/Data/processed_data/singlem_output"
DIAMOND_DIR="${BASE_DIR}/Data/functional_analysis/methane_genes"
LOG_DIR="${BASE_DIR}/Logs"

# Create output directories
mkdir -p ${BBDUK_DIR}
mkdir -p ${SINGLEM_DIR}
mkdir -p ${DIAMOND_DIR}
mkdir -p ${LOG_DIR}

# Set database paths - auto-detect or use environment variables
if [ -z "$ADAPTERS" ]; then
    # Try to find BBMap adapters
    if [ -f "${CONDA_BASE}/envs/qiime2/opt/bbmap-39.10-0/resources/adapters.fa" ]; then
        ADAPTERS="${CONDA_BASE}/envs/qiime2/opt/bbmap-39.10-0/resources/adapters.fa"
    elif command -v bbduk.sh &> /dev/null; then
        BBMAP_PATH=$(dirname $(which bbduk.sh))
        ADAPTERS="${BBMAP_PATH}/../resources/adapters.fa"
    else
        echo "WARNING: Could not find BBMap adapters, will skip adapter trimming"
        ADAPTERS=""
    fi
fi

SINGLEM_METAPACKAGE="${SINGLEM_METAPACKAGE:-${BASE_DIR}/Data/reference_dbs/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb}"
METHANE_DB="${METHANE_DB:-${BASE_DIR}/Data/reference_dbs/DIAMOND/methane_master_db}"

# Check databases
echo "Checking databases..."
echo ""

DB_MISSING=false

# Check SingleM database
if [ -e "$SINGLEM_METAPACKAGE" ] || [ -d "$SINGLEM_METAPACKAGE" ]; then
    echo "  ✓ SingleM database found"
else
    echo "  ❌ SingleM database NOT found"
    echo "     Expected: $SINGLEM_METAPACKAGE"
    DB_MISSING=true
fi

# Check DIAMOND database
if [ -f "${METHANE_DB}.dmnd" ]; then
    echo "  ✓ DIAMOND database found: ${METHANE_DB}.dmnd"
elif [ -f "${METHANE_DB}" ]; then
    echo "  ✓ DIAMOND database found: ${METHANE_DB}"
else
    echo "  ❌ DIAMOND database NOT found"
    echo "     Expected: ${METHANE_DB}.dmnd"
    DB_MISSING=true
fi

if [ "$DB_MISSING" = true ]; then
    echo ""
    echo "========================================================================"
    echo "ERROR: Required databases are missing"
    echo "======================================================================="
    echo ""
    echo "To download SingleM database (15-20GB):"
    echo "  conda activate quick_search"
    echo "  singlem data --output-directory Data/reference_dbs/"
    echo ""
    echo "For DIAMOND database:"
    echo "  Copy or create methane_master_db.dmnd in Data/reference_dbs/DIAMOND/"
    echo ""
    exit 1
fi

echo ""
echo "All prerequisites checked - ready to start!"
echo ""

# ============================================================================
# Sample Processing
# ============================================================================

################################################################################
# Function: Process single sample
################################################################################
process_sample() {
    local SAMPLE_ID=$1
    local START_TIME=$(date +%s)
    local LOG_FILE
    LOG_FILE="${LOG_DIR}/${SAMPLE_ID}_quick_search.log"

    (
    echo "========================================================================"
    echo "Processing Sample: ${SAMPLE_ID}"
    echo "========================================================================"
    echo "Start time: $(date)"
    echo "======================================================================="
    
    # Input files - support both naming patterns
    # Pattern 1: SAMPLE_R1.fastq.gz (standard)
    # Pattern 2: SAMPLE_R1_LABEL.fastq.gz (with sample label)
    if [ -f "${RAW_DATA_DIR}/${SAMPLE_ID}_R1.fastq.gz" ]; then
        R1_RAW="${RAW_DATA_DIR}/${SAMPLE_ID}_R1.fastq.gz"
        R2_RAW="${RAW_DATA_DIR}/${SAMPLE_ID}_R2.fastq.gz"
    else
        # Try to find files with pattern SAMPLE_R1_*.fastq.gz
        R1_RAW=$(find "${RAW_DATA_DIR}" -name "${SAMPLE_ID}_R1_*.fastq.gz" -o -name "${SAMPLE_ID}_R1.fastq.gz" | head -1)
        R2_RAW=$(find "${RAW_DATA_DIR}" -name "${SAMPLE_ID}_R2_*.fastq.gz" -o -name "${SAMPLE_ID}_R2.fastq.gz" | head -1)
    fi

    # BBduk output files
    R1_CLEAN="${BBDUK_DIR}/${SAMPLE_ID}_R1_bbduk.fastq.gz"
    R2_CLEAN="${BBDUK_DIR}/${SAMPLE_ID}_R2_bbduk.fastq.gz"
    BBDUK_STATS="${BBDUK_DIR}/${SAMPLE_ID}_bbduk_stats.txt"

    # SingleM output files
    SINGLEM_OUTPUT="${SINGLEM_DIR}/singlem_output_${SAMPLE_ID}"
    SINGLEM_OTU="${SINGLEM_OUTPUT}/singlem_otu_table.csv"
    SINGLEM_PROFILE="${SINGLEM_OUTPUT}/singlem_profile.txt"

    # DIAMOND输出文件
    DIAMOND_OUTPUT="${DIAMOND_DIR}/${SAMPLE_ID}_combined_methane_hits.txt"
    DIAMOND_SUMMARY="${DIAMOND_DIR}/${SAMPLE_ID}_methane_summary.txt"
    
    # Log file
    :

################################################################################
# STEP 1: BBduk Quality Control
################################################################################
echo "=========================================="
echo "STEP 1: Quality Control with BBduk"
echo "Sample: ${SAMPLE_ID}"
echo "=========================================="

if [ ! -f "${R1_RAW}" ] || [ -z "${R1_RAW}" ]; then
    echo "ERROR: R1 raw data not found"
    echo "  Sample ID: ${SAMPLE_ID}"
    echo "  Location: Data/raw_data/${RAW_DATA_SUBDIR}/"
    echo "  Looking for: ${SAMPLE_ID}_R1*.fastq.gz"
    echo ""
    echo "  Available files in ${RAW_DATA_DIR}:"
    ls -la "${RAW_DATA_DIR}"/${SAMPLE_ID}*.fastq.gz 2>/dev/null || echo "    No matching files found"
    echo ""
    echo "  Please check:"
    echo "    1. Sample ID is correct"
    echo "    2. Files are in the raw_data/${RAW_DATA_SUBDIR} directory"
    echo "    3. File naming matches: <SAMPLE_ID>_R1*.fastq.gz"
    return 1
fi

if [ ! -f "${R2_RAW}" ] || [ -z "${R2_RAW}" ]; then
    echo "ERROR: R2 raw data not found"
    echo "  Sample ID: ${SAMPLE_ID}"
    echo "  Location: Data/raw_data/${RAW_DATA_SUBDIR}/"
    echo "  Looking for: ${SAMPLE_ID}_R2*.fastq.gz"
    echo ""
    echo "  Available files in ${RAW_DATA_DIR}:"
    ls -la "${RAW_DATA_DIR}"/${SAMPLE_ID}*.fastq.gz 2>/dev/null || echo "    No matching files found"
    echo ""
    echo "  Please check:"
    echo "    1. Sample ID is correct"
    echo "    2. Files are in the raw_data/${RAW_DATA_SUBDIR} directory"
    echo "    3. File naming matches: <SAMPLE_ID>_R2*.fastq.gz"
    return 1
fi

echo "Found input files:"
echo "  R1: $(basename ${R1_RAW})"
echo "  R2: $(basename ${R2_RAW})"
echo ""

echo "Running BBduk..."

# Wrapper to run BBduk safely and retry with constrained Java options if heap-init errors occur
run_bbduk() {
    local LOGFILE="${LOG_DIR}/${SAMPLE_ID}_bbduk.err"
    # First attempt: capture both stdout and stderr
    if bbduk.sh "$@" >"${LOGFILE}" 2>&1; then
        return 0
    fi

    # Check for Java heap-init errors in combined output
    if grep -q "Initial heap size set to a larger value than the maximum heap size" "${LOGFILE}" >/dev/null 2>&1; then
        echo "WARNING: bbduk Java heap error detected; retrying with constrained Java options (JAVA_TOOL_OPTIONS)..."
        # Retry using JAVA_TOOL_OPTIONS to override any internal JVM settings
        if JAVA_TOOL_OPTIONS="-Xmx4g -Xms512m" bbduk.sh "$@" >"${LOGFILE}" 2>&1; then
            return 0
        fi

        echo "WARNING: Retry with JAVA_TOOL_OPTIONS failed; attempting direct java invocation as last resort..."
        # Attempt to find BBMap installation dir
        BBMAP_BIN=$(which bbduk.sh 2>/dev/null || true)
        BBMAP_DIR=""
        if [ -n "${BBMAP_BIN}" ]; then
            BBMAP_DIR="$(cd "$(dirname "${BBMAP_BIN}")/.." && pwd 2>/dev/null || echo "")"
        fi
        # If BBMAP_DIR still empty, fall back to conda location
        if [ -z "${BBMAP_DIR}" ]; then
            BBMAP_DIR="/opt/miniconda3/envs/quick_search/opt/bbmap-*/current/"
        fi

        # Direct java invocation; append stderr to log
        java -ea -Xmx4g -Xms512m -cp "${BBMAP_DIR}" jgi.BBDuk "$@" >>"${LOGFILE}" 2>&1 || return 1
        return 0
    fi

    # No specific heap-init error found; leave log for inspection and return failure
    return 1
}

if [ -n "$ADAPTERS" ] && [ -f "$ADAPTERS" ]; then
    run_bbduk \
        in1=${R1_RAW} \
        in2=${R2_RAW} \
        out1=${R1_CLEAN} \
        out2=${R2_CLEAN} \
        ref=${ADAPTERS} \
        ktrim=r \
        k=23 \
        mink=11 \
        hdist=1 \
        tpe \
        tbo \
        qtrim=rl \
        trimq=20 \
        minlen=50 \
        stats=${BBDUK_STATS} \
        threads=${THREADS}
else
    echo "WARNING: Adapters not found, running BBduk without adapter trimming"
    run_bbduk \
        in1=${R1_RAW} \
        in2=${R2_RAW} \
        out1=${R1_CLEAN} \
        out2=${R2_CLEAN} \
        qtrim=rl \
        trimq=20 \
        minlen=50 \
        stats=${BBDUK_STATS} \
        threads=${THREADS}
fi

if [ $? -ne 0 ]; then
    echo "ERROR: BBduk failed for sample ${SAMPLE_ID}"
    return 1
fi

echo "BBduk completed!"
echo "Clean reads saved to: ${BBDUK_DIR}"
echo ""

################################################################################
# STEP 2: SingleM Taxonomic Profiling
################################################################################
echo "=========================================="
echo "STEP 2: Taxonomic Profiling with SingleM"
echo "=========================================="

# Activate singlem environment
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate quick_search 2>/dev/null || echo "WARNING: Could not activate quick_search environment"
fi

# Set environment variables
export SINGLEM_METAPACKAGE_PATH=${SINGLEM_METAPACKAGE}

# Set temporary directory for SingleM (avoid /tmp issues on HPC)
if [ -z "$TMPDIR" ]; then
    TMPDIR="${BASE_DIR}/Data/temp"
    mkdir -p "${TMPDIR}"
    export TMPDIR
    export TEMP="${TMPDIR}"
    export TMP="${TMPDIR}"
fi

# Create a unique temporary directory for this sample
SAMPLE_TMPDIR="${TMPDIR}/singlem_${SAMPLE_ID}_$$"
mkdir -p "${SAMPLE_TMPDIR}"
export TMPDIR="${SAMPLE_TMPDIR}"

# Set Python environment to keep temp files longer
export PYTHONDONTWRITEBYTECODE=1

echo "Using temporary directory: ${TMPDIR}"
# Use more threads for SingleM/DIAMOND to improve performance
# Only cap on low-memory systems
SINGLEM_THREADS=${SINGLEM_THREADS:-${THREADS}}
export OMP_NUM_THREADS=1
# Try to increase stack size but ignore permission warnings
ulimit -s 65536 2>/dev/null || true

echo "Running SingleM pipe..."
singlem pipe \
    --forward ${R1_CLEAN} \
    --reverse ${R2_CLEAN} \
    --otu-table ${SINGLEM_OTU} \
    --threads ${SINGLEM_THREADS} \
    --output-extras

if [ $? -ne 0 ]; then
    echo "ERROR: SingleM pipe failed for sample ${SAMPLE_ID}"

    # Diagnostic: run a small DIAMOND test on a subset of reads to see if DIAMOND itself fails
    echo "Running DIAMOND diagnostic (subset of reads, threads=1)..."
    mkdir -p "${DIAMOND_DIR}" 2>/dev/null || true
    DIAMOND_DIAG_ERR="${DIAMOND_DIR}/${SAMPLE_ID}_diamond_diag.err"
    set -o pipefail
    gunzip -c ${R1_CLEAN} ${R2_CLEAN} | head -n 400000 | diamond blastx \
        -d ${METHANE_DB} \
        -q - \
        -o /dev/null \
        --threads 1 \
        --tmpdir "${SAMPLE_TMPDIR}" \
        --sensitive 2> "${DIAMOND_DIAG_ERR}" || echo "DIAMOND diagnostic exit code: $? (see ${DIAMOND_DIAG_ERR})"
    set +o pipefail

    echo "See diagnostic file: ${DIAMOND_DIAG_ERR}"
    return 1
fi

echo "Running SingleM summarise..."
singlem summarise \
    --input-otu-tables ${SINGLEM_OTU} \
    --output-taxonomic-profile ${SINGLEM_PROFILE}

if [ $? -ne 0 ]; then
    echo "ERROR: SingleM summarise failed for sample ${SAMPLE_ID}"
    return 1
fi

echo "SingleM completed!"
echo "OTU table: ${SINGLEM_OTU}"
echo "Taxonomic profile: ${SINGLEM_PROFILE}"

# Clean up temporary directory for this sample
rm -rf "${SAMPLE_TMPDIR}" 2>/dev/null || true

echo ""

################################################################################
# STEP 3: DIAMOND Functional Gene Search
################################################################################
echo "=========================================="
echo "STEP 3: Methane Gene Search with DIAMOND"
echo "=========================================="

echo "Searching for methane metabolism genes..."
# Use a dedicated writable temp dir for DIAMOND
DIAMOND_TMPDIR="${TMPDIR}/diamond_${SAMPLE_ID}_$$"
mkdir -p "${DIAMOND_TMPDIR}" 2>/dev/null || true

gunzip -c ${R1_CLEAN} ${R2_CLEAN} | \
    DIAMOND_TMPDIR="${DIAMOND_TMPDIR}" diamond blastx \
    -d ${METHANE_DB} \
    -q - \
    -o ${DIAMOND_OUTPUT} \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
    --threads ${THREADS} \
    --sensitive \
    --query-cover 80 \
    --min-score 40 \
    --max-target-seqs 5 \
    --evalue 1e-5

# Clean DIAMOND temp dir
rm -rf "${DIAMOND_TMPDIR}" 2>/dev/null || true

if [ $? -ne 0 ]; then
    echo "ERROR: DIAMOND search failed for sample ${SAMPLE_ID}"
    return 1
fi

echo "DIAMOND search completed!"
echo "Results saved to: ${DIAMOND_OUTPUT}"
echo ""

################################################################################
# STEP 4: Generate Summary Report
################################################################################
echo "=========================================="
echo "STEP 4: Generating Summary Report"
echo "=========================================="

# Count results
TOTAL_HITS=$(wc -l < ${DIAMOND_OUTPUT})
MMO_HITS=$(grep -ic "methane monooxygenase" ${DIAMOND_OUTPUT})
PMO_HITS=$(grep -iEc "PmoA|pMMO" ${DIAMOND_OUTPUT})
TOTAL_METHANE=$((MMO_HITS + PMO_HITS))

# Generate report
cat > ${DIAMOND_SUMMARY} << EOF
===================================================================
METHANE METABOLISM GENE SEARCH RESULTS - Sample ${SAMPLE_ID}
===================================================================
Analysis Date: $(date +%Y-%m-%d)
Search Tool: DIAMOND blastx
Database: methane_master_db
Search Parameters:
  - Query coverage: >80%
  - Bit score: >40
  - E-value: <1e-5
  - Max target sequences: 5
  - Mode: sensitive

===================================================================
RESULTS SUMMARY
===================================================================

Combined Reads (R1 + R2 from bbduk_cleaned):
  Total alignments: ${TOTAL_HITS}
  
  Methane-specific genes:
    - Methane monooxygenase (MMO): ${MMO_HITS} hits
    - PmoA family proteins (pMMO): ${PMO_HITS} hits
    - Total methanotroph genes: ${TOTAL_METHANE} hits

===================================================================
BIOLOGICAL INTERPRETATION
===================================================================

EOF

if [ ${TOTAL_METHANE} -gt 0 ]; then
    cat >> ${DIAMOND_SUMMARY} << EOF
METHANOTROPHS DETECTED
  Sample ${SAMPLE_ID} contains aerobic methane-oxidizing bacteria.
  Key enzymes detected:
  - pMMO (particulate methane monooxygenase): Primary methane oxidation
  - MMO (soluble methane monooxygenase): Alternative pathway
  
  These organisms consume methane (CH4) as carbon and energy source.

===================================================================
ECOLOGICAL SIGNIFICANCE
===================================================================

The presence of methanotrophs suggests:
1. Oxic environment (aerobic methane oxidation)
2. Potential methane sink rather than source
3. Common in wetlands, rice paddies, landfill covers

===================================================================
EOF
else
    cat >> ${DIAMOND_SUMMARY} << EOF
NO SIGNIFICANT METHANE METABOLISM GENES DETECTED
  
  This sample shows low or no methanotroph/methanogen activity.

===================================================================
EOF
fi

echo "Summary report saved to: ${DIAMOND_SUMMARY}"
echo ""

################################################################################
# STEP 5: Gene Abundance Quantification and Normalization
################################################################################
echo "=========================================="
echo "STEP 5: Gene Abundance Quantification"
echo "=========================================="

# Output files for quantification
GENE_COUNTS="${DIAMOND_DIR}/${SAMPLE_ID}_gene_counts.txt"
GENE_RPKM="${DIAMOND_DIR}/${SAMPLE_ID}_gene_rpkm.txt"
FILTERED_HITS="${DIAMOND_DIR}/${SAMPLE_ID}_filtered_methane_hits.txt"

# Filter high-quality hits only (remove ribosomal proteins and low-quality hits)
echo "Filtering high-quality methane gene hits..."
grep -iE "methane|pmo|mmo" ${DIAMOND_OUTPUT} | \
    grep -viE "ribosom|ribokinase|ribonucleoside|ammonia" > ${FILTERED_HITS}

FILTERED_COUNT=$(wc -l < ${FILTERED_HITS})
echo "Filtered hits: ${FILTERED_COUNT} (removed non-methane genes)"

# Count reads per gene
echo "Counting reads per gene..."
awk '{print $7}' ${FILTERED_HITS} | sort | uniq -c | sort -rn > ${GENE_COUNTS}

echo "Top 10 most abundant genes:"
head -10 ${GENE_COUNTS}
echo ""

# Calculate normalization factors
# Get total mapped reads (from original DIAMOND output)
TOTAL_MAPPED_READS=$(awk '{print $1}' ${DIAMOND_OUTPUT} | sort -u | wc -l)

# Get total sequencing reads (from both R1 and R2) - optimized single pass
TOTAL_READS=$(( $(zcat ${R1_CLEAN} ${R2_CLEAN} | wc -l) / 4 ))
TOTAL_READS_MILLIONS=$(echo "scale=4; ${TOTAL_READS} / 1000000" | bc)

echo "Normalization factors:"
echo "  Total reads (R1+R2): ${TOTAL_READS}"
echo "  Total reads in millions: ${TOTAL_READS_MILLIONS}M"
echo "  Total mapped reads: ${TOTAL_MAPPED_READS}"
echo ""

# Calculate RPKM for each gene
# RPKM = (Read counts * 10^9) / (Gene length in bp * Total mapped reads)
# For protein-level alignment, we use alignment length as proxy for gene length
echo "Calculating RPKM values..."

awk -v total_reads="${TOTAL_READS_MILLIONS}" -v filtered_hits="${FILTERED_HITS}" '
BEGIN {
    print "Gene_Name\tRead_Count\tAvg_Length\tRPKM"
}
{
    gene_name = $0
    gsub(/^[[:space:]]+/, "", gene_name)
    split(gene_name, arr, " ")
    count = arr[1]
    gene = ""
    for(i=2; i<=length(arr); i++) gene = gene " " arr[i]
    gsub(/^[[:space:]]+/, "", gene)
    
    # Get average alignment length for this gene
    cmd = "grep -F \"" gene "\" " filtered_hits " | awk '"'"'{sum+=$4; n++} END {if(n>0) print sum/n; else print 1000}'"'"'"
    cmd | getline avg_length
    close(cmd)
    
    # Calculate RPKM: (count * 1000) / (avg_length * total_reads_millions)
    rpkm = (count * 1000) / (avg_length * total_reads)
    
    printf "%s\t%d\t%.1f\t%.4f\n", gene, count, avg_length, rpkm
}
' ${GENE_COUNTS} > ${GENE_RPKM}

echo "RPKM calculation completed!"
echo "Results saved to: ${GENE_RPKM}"
echo ""

echo "Top 10 genes by RPKM:"
head -11 ${GENE_RPKM} | tail -10
echo ""

################################################################################
# STEP 6: Display Results
################################################################################
echo "=========================================="
echo "ANALYSIS COMPLETE - Sample ${SAMPLE_ID}"
echo "=========================================="
echo ""
echo "Output Files:"
echo "  BBduk cleaned reads:"
echo "    - ${R1_CLEAN}"
echo "    - ${R2_CLEAN}"
echo ""
echo "  SingleM results:"
echo "    - ${SINGLEM_OTU}"
echo "    - ${SINGLEM_PROFILE}"
echo ""
echo "  DIAMOND results:"
echo "    - ${DIAMOND_OUTPUT}"
echo "    - ${DIAMOND_SUMMARY}"
echo ""
echo "  Gene abundance analysis:"
echo "    - ${FILTERED_HITS}"
echo "    - ${GENE_COUNTS}"
echo "    - ${GENE_RPKM}"
echo ""
echo "Quick Results:"
echo "  Total DIAMOND hits: ${TOTAL_HITS}"
echo "  Filtered methane gene hits: ${FILTERED_COUNT}"
echo "  Methane monooxygenase: ${MMO_HITS} hits"
echo "  PmoA family proteins: ${PMO_HITS} hits"
echo "  Total methanotroph genes: ${TOTAL_METHANE} hits"
echo "  Total sequencing reads: ${TOTAL_READS}"
echo ""
    local END_TIME=$(date +%s)
    local ELAPSED=$((END_TIME - START_TIME))
    echo "All steps completed successfully!"
    echo "Time elapsed: $((ELAPSED / 60)) minutes $((ELAPSED % 60)) seconds"
    echo "=========================================="
    echo ""
    ) > >(tee -a "$LOG_FILE") 2>&1
    return $?
}

################################################################################
# Main execution loop
################################################################################
MAIN_START_TIME=$(date +%s)
SUCCESS_COUNT=0
FAIL_COUNT=0
FAILED_SAMPLES=()

# Parallel processing function
if [ ${MAX_PARALLEL_SAMPLES} -gt 1 ]; then
    echo "Processing ${#SAMPLE_LIST[@]} samples with ${MAX_PARALLEL_SAMPLES} parallel jobs..."
    echo ""
    
    # Use GNU parallel or xargs for parallel processing
    if command -v parallel &> /dev/null; then
        export -f process_sample
        export BASE_DIR THREADS RAW_DATA_DIR BBDUK_DIR SINGLEM_DIR DIAMOND_DIR LOG_DIR
        export ADAPTERS SINGLEM_METAPACKAGE METHANE_DB SINGLEM_THREADS RAW_DATA_SUBDIR
        
        printf '%s\n' "${SAMPLE_LIST[@]}" | \
            parallel -j ${MAX_PARALLEL_SAMPLES} --line-buffer \
            'if process_sample {}; then echo "✓ {} completed"; else echo "✗ {} failed"; fi'
        
        # Count successes/failures from logs
        for SAMPLE_ID in "${SAMPLE_LIST[@]}"; do
            if grep -q "All steps completed successfully" "${LOG_DIR}/${SAMPLE_ID}_quick_search.log" 2>/dev/null; then
                SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
            else
                FAIL_COUNT=$((FAIL_COUNT + 1))
                FAILED_SAMPLES+=("$SAMPLE_ID")
            fi
        done
    else
        echo "WARNING: GNU parallel not found, falling back to sequential processing"
        echo "Install with: conda install -c conda-forge parallel"
        echo ""
        MAX_PARALLEL_SAMPLES=1
    fi
fi

# Sequential processing (default or fallback)
if [ ${MAX_PARALLEL_SAMPLES} -eq 1 ]; then
    for SAMPLE_ID in "${SAMPLE_LIST[@]}"; do
        echo ""
        echo "################################################################################"
        echo "# Processing sample ${SAMPLE_ID} ($(($SUCCESS_COUNT + $FAIL_COUNT + 1))/${#SAMPLE_LIST[@]})"
        echo "################################################################################"
        
        if process_sample "$SAMPLE_ID"; then
            SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
            echo "✓ Sample ${SAMPLE_ID} completed successfully"
        else
            FAIL_COUNT=$((FAIL_COUNT + 1))
            FAILED_SAMPLES+=("$SAMPLE_ID")
            echo "✗ Sample ${SAMPLE_ID} failed"
        fi
    done
fi

################################################################################
# Final Summary
################################################################################
MAIN_END_TIME=$(date +%s)
TOTAL_ELAPSED=$((MAIN_END_TIME - MAIN_START_TIME))

echo ""
echo "========================================================================"
echo "PIPELINE SUMMARY"
echo "========================================================================"
echo "Total samples processed: ${#SAMPLE_LIST[@]}"
echo "Successful: ${SUCCESS_COUNT}"
echo "Failed: ${FAIL_COUNT}"
echo "Total time: $((TOTAL_ELAPSED / 3600))h $((TOTAL_ELAPSED % 3600 / 60))m $((TOTAL_ELAPSED % 60))s"
echo ""

if [ ${FAIL_COUNT} -gt 0 ]; then
    echo "Failed samples:"
    for SAMPLE in "${FAILED_SAMPLES[@]}"; do
        echo "  - ${SAMPLE}"
    done
    echo ""
    echo "Check logs in: ${LOG_DIR}/"
    echo "========================================================================"
    exit 1
fi

echo "All samples completed successfully!"
echo ""
echo "Output locations:"
echo "  Quality-controlled reads: Data/processed_data/bbduk_cleaned/"
echo "  Taxonomic profiles:       Data/processed_data/singlem_output/"
echo "  Gene detection results:   Data/functional_analysis/methane_genes/"
echo "  Logs:                     Logs/"
echo "======================================================================="
