#!/bin/bash -l
#SBATCH --job-name=quick_search        # Job name
#SBATCH --output=quick_search_%j.out   # Standard output log (%j = job ID)
#SBATCH --error=quick_search_%j.err    # Standard error log
#SBATCH --ntasks=1                     # Number of tasks (processes)
#SBATCH --cpus-per-task=8              # Number of CPU cores per task (reduced for stability)
#SBATCH --mem=64G                      # Memory per node
#SBATCH --time=24:00:00                # Time limit (24 hours)
#SBATCH --partition=large_336          # Partition name

################################################################################
# SLURM Job Script for Quick Metagenomic Search Pipeline
# 
# Usage on HPC:
#   Single sample:
#     sbatch submit_quick_search.sbatch 53394
#   
#   Multiple samples:
#     sbatch submit_quick_search.sbatch 53394 53395 53396
#   
#   From sample list file:
#     sbatch submit_quick_search.sbatch --sample-list samples.txt
#   
#   Auto-detect all samples:
#     sbatch submit_quick_search.sbatch --auto
#
################################################################################

# Print job information
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE"
echo "Start Time: $(date)"
echo "=========================================="
echo ""

# Set environment variables for HPC
export BASE_DIR="/home/jiayi/QIIME"
export THREADS=$SLURM_CPUS_PER_TASK
export CONDA_BASE="/home/jiayi/miniconda3"  # Adjust if your conda is elsewhere
export RAW_DATA_SUBDIR="shotgun"  # Specify shotgun subdirectory

# Set temporary directory to avoid /tmp space issues
export TMPDIR="${BASE_DIR}/Data/temp"
export TEMP="${TMPDIR}"
export TMP="${TMPDIR}"
mkdir -p "${TMPDIR}"

# Initialize Conda
if [ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
else
    eval "$(conda shell.bash hook)"
fi

# Activate the main environment (bbduk and diamond should be here)
conda activate quick_search

# Set database paths (adjust these for your HPC environment)
export SINGLEM_METAPACKAGE="${BASE_DIR}/Data/reference_dbs/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb"
export METHANE_DB="${BASE_DIR}/Data/reference_dbs/DIAMOND/methane_master_db"

# Find BBMap adapters
if [ -f "${CONDA_BASE}/envs/qiime2-moshpit/opt/bbmap*/resources/adapters.fa" ]; then
    export ADAPTERS=$(ls ${CONDA_BASE}/envs/qiime2-moshpit/opt/bbmap*/resources/adapters.fa 2>/dev/null | head -1)
fi

# Navigate to the scripts directory
cd ${BASE_DIR}/Scripts

# Run the quick search script with all arguments passed to sbatch
echo "Running quick_search.sh with parameters: $@"
echo ""

bash quick_search.sh "$@"

# Capture exit status
EXIT_STATUS=$?

# Deactivate Conda environment
conda deactivate

echo ""
echo "=========================================="
echo "Job completed with exit status: $EXIT_STATUS"
echo "End Time: $(date)"
echo "=========================================="

exit $EXIT_STATUS
