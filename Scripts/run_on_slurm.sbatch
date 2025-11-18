#!/bin/bash
#SBATCH --job-name=qiime_pipeline
#SBATCH --output=Logs/slurm-pipeline-%j.out
#SBATCH --error=Logs/slurm-pipeline-%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=large_336

# Example Slurm wrapper to run this pipeline non-interactively.
# Customize resource requests, partition, and paths for your cluster.

set -euo pipefail

# Determine submission/project directory. When submitted via sbatch, SLURM sets
# SLURM_SUBMIT_DIR to the directory where sbatch was invoked. Use that if present
# to avoid issues where $0 or relative paths differ on compute nodes.
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
cd "${SUBMIT_DIR}"

# Load modules if your cluster requires them (optional)
# module load anaconda

# Initialize conda and activate the environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qiime

# Export non-interactive variables here to control behavior.
# Examples (uncomment and edit as needed):
# export START_STEP=5
# export START_STEP=6
# export START_STEP=7
# NOTE: If you set START_STEP in the environment (uncomment above), the pipeline
# will start from that step instead of running all steps. Leave commented to run
# the full pipeline (steps 1-7) by default.
# export MODE=denoise
# export ENV_NAME=qiime
# export N_THREADS=8
# export TRUNC_LEN_F=240
# export TRUNC_LEN_R=200
# Use an existing classifier (recommended in batch):
# export CLASSIFIER_PATH=Data/reference_dbs/classifier.qza
# For decontamination step (Step 6):
# export DECONTAMINATION_CHOICE=1
# export NEGATIVE_CONTROLS=A1,A2

# If you want to build a classifier in batch, set BUILD_CLASSIFIER="silva" and
# allow RESCRIPt automatic install (set INSTALL_RESCRIPT=y) and provide PRIMER_F/PRIMER_R
# export BUILD_CLASSIFIER=silva
# export INSTALL_RESCRIPT=y
# export EXTRACT_PRIMERS=y
# export PRIMER_F=GTGCCAGCMGCCGCGGTAA
# export PRIMER_R=GGACTACHVGGGTWTCTAAT

# Run the pipeline in non-interactive (batch) mode. Provide a start step if desired.
# Examples:
# sbatch this_script -> runs default START_STEP (if set in env) or step 1
# sbatch this_script && will execute according to exported START_STEP

# By default run the full pipeline (all steps) in non-interactive batch mode.
# If you want to run a single start step, set START_STEP in the environment before
# submitting (uncomment and set above) or edit this file.
echo "Starting full pipeline (steps 1-7) in non-interactive batch mode"

# Prefer a specialized HPC wrapper if present; otherwise fall back to the main script
HPC_WRAPPER_PATH="${SUBMIT_DIR}/Scripts/main_hpc.sh"
MAIN_SCRIPT_PATH="${SUBMIT_DIR}/Scripts/main.sh"

if [[ -f "${HPC_WRAPPER_PATH}" ]]; then
	echo "Found ${HPC_WRAPPER_PATH} - using HPC wrapper"
	bash "${HPC_WRAPPER_PATH}" -b "$@"
elif [[ -f "${MAIN_SCRIPT_PATH}" ]]; then
	echo "Warning: ${HPC_WRAPPER_PATH} not found. Falling back to ${MAIN_SCRIPT_PATH} -b"
	bash "${MAIN_SCRIPT_PATH}" -b "$@"
else
	echo "ERROR: Neither ${HPC_WRAPPER_PATH} nor ${MAIN_SCRIPT_PATH} exist in ${SUBMIT_DIR}" >&2
	exit 1
fi
