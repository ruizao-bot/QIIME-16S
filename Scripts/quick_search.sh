#!/bin/bash

################################################################################
# Quick Metagenomic Analysis Pipeline
# Automated workflow for shotgun metagenomic analysis
################################################################################
# Steps:
# 1. BBduk quality control
# 2. SingleM taxonomic profiling (OTU generation)
# 3. DIAMOND functional gene search (methane metabolism genes)
################################################################################

# Default configuration
THREADS=${THREADS:-4}              # Number of threads (can be set via environment variable)
CONDA_BASE="${CONDA_BASE:-/opt/anaconda3}"  # Conda installation path
RAW_DATA_SUBDIR="${RAW_DATA_SUBDIR:-shotgun}"  # Subdirectory for raw data (shotgun or empty for root)

# Auto-detect base directory if not set
if [ -z "$BASE_DIR" ]; then
    # Try to find the project root by looking for Scripts directory
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    BASE_DIR="$(dirname "$SCRIPT_DIR")"
fi

# Parse command line arguments
SAMPLE_LIST=()
SAMPLE_FILE=""
AUTO_DETECT=false

while [[ $# -gt 0 ]]; do
    case $1 in
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
        --base-dir)
            BASE_DIR="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: bash quick_search.sh [OPTIONS] [SAMPLE_IDs...]"
            echo ""
            echo "Options:"
            echo "  --sample-list FILE    Read sample IDs from file (one per line)"
            echo "  --auto                Auto-detect all samples in raw_data directory"
            echo "  --threads N           Number of threads (default: 4)"
            echo "  --base-dir PATH       Base directory path (default: current)"
            echo "  --help, -h            Show this help message"
            echo ""
            echo "Examples:"
            echo "  bash quick_search.sh 53394                      # Single sample"
            echo "  bash quick_search.sh 53394 53395 53396          # Multiple samples"
            echo "  bash quick_search.sh --sample-list samples.txt  # From file"
            echo "  bash quick_search.sh --auto                     # Auto-detect all"
            echo "  bash quick_search.sh --threads 16 53394         # With 16 threads"
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
    echo "ERROR: No samples provided"
    echo "Usage: bash quick_search.sh [OPTIONS] [SAMPLE_IDs...]"
    echo "Use --help for more information"
    exit 1
fi

echo "========================================="
echo "Processing ${#SAMPLE_LIST[@]} sample(s)"
echo "Threads: ${THREADS}"
echo "Base directory: ${BASE_DIR}"
echo "========================================="
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

################################################################################
# Function: Process single sample
################################################################################
process_sample() {
    local SAMPLE_ID=$1
    local START_TIME=$(date +%s)
    
    echo "========================================="
    echo "Processing Sample: ${SAMPLE_ID}"
    echo "Start time: $(date)"
    echo "Base directory: ${BASE_DIR}"
    echo "Raw data directory: ${RAW_DATA_DIR}"
    echo "========================================="
    
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
    LOG_FILE="${LOG_DIR}/${SAMPLE_ID}_quick_search.log"
    exec > >(tee -a "$LOG_FILE") 2>&1

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
if [ -n "$ADAPTERS" ] && [ -f "$ADAPTERS" ]; then
    bbduk.sh \
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
    bbduk.sh \
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
echo "Running SingleM pipe..."
singlem pipe \
    --forward ${R1_CLEAN} \
    --reverse ${R2_CLEAN} \
    --otu-table ${SINGLEM_OTU} \
    --threads ${THREADS} \
    --output-extras

if [ $? -ne 0 ]; then
    echo "ERROR: SingleM pipe failed for sample ${SAMPLE_ID}"
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
gunzip -c ${R1_CLEAN} ${R2_CLEAN} | \
    diamond blastx \
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

# Get total sequencing reads (from both R1 and R2)
TOTAL_R1_READS=$(zcat ${R1_CLEAN} | wc -l | awk '{print $1/4}')
TOTAL_R2_READS=$(zcat ${R2_CLEAN} | wc -l | awk '{print $1/4}')
TOTAL_READS=$(echo "${TOTAL_R1_READS} + ${TOTAL_R2_READS}" | bc)
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

awk -v total_reads="${TOTAL_READS_MILLIONS}" '
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
    cmd = "grep -F \"" gene "\" '${FILTERED_HITS}' | awk '"'"'{sum+=$4; n++} END {if(n>0) print sum/n; else print 1000}'"'"'"
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
}

################################################################################
# Main execution loop
################################################################################
MAIN_START_TIME=$(date +%s)
SUCCESS_COUNT=0
FAIL_COUNT=0
FAILED_SAMPLES=()

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

################################################################################
# Final Summary
################################################################################
MAIN_END_TIME=$(date +%s)
TOTAL_ELAPSED=$((MAIN_END_TIME - MAIN_START_TIME))

echo ""
echo "################################################################################"
echo "# PIPELINE SUMMARY"
echo "################################################################################"
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
    exit 1
fi

echo "All samples completed successfully!"
echo "################################################################################"
