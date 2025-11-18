#!/usr/bin/env bash

# HPC wrapper for running the QIIME pipeline on shared clusters (Slurm)
# - Creates a job-local workspace under $SCRATCH (or $TMPDIR) when available
# - Initializes conda robustly for non-interactive batch runs
# - Exports environment variables consumed by `Scripts/main.sh` and runs it in batch mode
# - Captures logs to the job-local directory and copies them back to the project Logs

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Determine scratch base: prefer $SCRATCH, then $TMPDIR, fall back to /tmp
SCRATCH_BASE="${SCRATCH:-${TMPDIR:-/tmp}}"
export SCRATCH="${SCRATCH_BASE}"
export TMPDIR="${TMPDIR:-${SCRATCH}/tmp}"

# Job id: use SLURM_JOB_ID if present; otherwise create a manual timestamp id
JOB_ID="${SLURM_JOB_ID:-manual-$(date +%s)}"

# Job-local working directories
JOB_DIR="${SCRATCH}/QIIME/${USER:-$(whoami)}/job-${JOB_ID}"
LOG_DIR="${JOB_DIR}/Logs"

# Small retrying mkdir to avoid NFS race conditions
mkdir_retry() {
    local d="$1"
    local max=5
    local i=0
    umask 002
    until mkdir -p "$d" 2>/dev/null; do
        i=$((i+1))
        if (( i >= max )); then
            echo "ERROR: failed to create directory $d" >&2
            return 1
        fi
        sleep $((RANDOM % 3 + 1))
    done
    chmod g+rwx "$d" 2>/dev/null || true
}

mkdir_retry "${JOB_DIR}"
mkdir_retry "${LOG_DIR}"

echo "[${JOB_ID}] Job workspace: ${JOB_DIR}"
echo "[${JOB_ID}] Log dir: ${LOG_DIR}"

# Ensure non-interactive mode for the pipeline wrapper
export NON_INTERACTIVE="true"

# Preserve any START_STEP provided by the user; otherwise default to 1
export START_STEP="${START_STEP:-1}"

# Make sure conda is initialised for this shell. On many HPC systems conda is available but
# activation requires eval "$(conda shell.bash hook)" or sourcing conda.sh from conda base.
if command -v conda >/dev/null 2>&1; then
    # Try the preferred shell hook first
    if ! eval "$(conda shell.bash hook)" 2>/dev/null; then
        # Fall back to sourcing the conda profile script
        CONDA_BASE="$(conda info --base 2>/dev/null || echo '')"
        if [[ -n "${CONDA_BASE}" && -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
            # shellcheck disable=SC1090
            source "${CONDA_BASE}/etc/profile.d/conda.sh"
        fi
    fi
else
    echo "Warning: conda not found in PATH. If your QIIME environment lives in a module, load it before running this wrapper." >&2
fi

# Export a location for the pipeline to write logs and checkpoints (main.sh will prefer SCRATCH and SLURM_JOB_ID)
export SCRATCH

# Run the existing main pipeline script in batch mode and capture its output to the job-local log
PIPELINE_LOG="${LOG_DIR}/qiime2_pipeline.log"
echo "[${JOB_ID}] Running pipeline: ${SCRIPT_DIR}/main.sh -b (logs -> ${PIPELINE_LOG})"

# Run pipeline and tee output to job-local log file
"${SCRIPT_DIR}/main.sh" -b "$@" 2>&1 | tee "${PIPELINE_LOG}"

# After completion, copy logs back to the project Logs directory so results are persistent
DEST_LOG_DIR="${PROJECT_DIR}/Logs/job-${JOB_ID}"
mkdir -p "${DEST_LOG_DIR}"
cp -a "${LOG_DIR}/." "${DEST_LOG_DIR}/" || echo "Warning: failed to copy logs back to project directory"

echo "[${JOB_ID}] Pipeline finished. Logs copied to: ${DEST_LOG_DIR}"

exit 0
