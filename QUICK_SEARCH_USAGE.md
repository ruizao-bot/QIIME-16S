# Quick Metagenomic Search Pipeline - Usage Guide

## Overview
Fast metagenomic analysis pipeline for methane metabolism gene detection with support for batch processing and HPC execution.

## Features
- ✅ Batch processing multiple samples
- ✅ HPC/SLURM support with resource management
- ✅ Flexible thread/CPU configuration
- ✅ Automated environment detection
- ✅ Comprehensive error handling and logging

---

## Local Usage (Mac/Linux)

### Single Sample
```bash
bash Scripts/quick_search.sh 53394
```

### Multiple Samples
```bash
bash Scripts/quick_search.sh 53394 53395 53396
```

### From Sample List File
```bash
bash Scripts/quick_search.sh --sample-list Config/samples.txt
```

### Auto-detect All Samples
```bash
bash Scripts/quick_search.sh --auto
```

### Custom Thread Count
```bash
bash Scripts/quick_search.sh --threads 8 53394
```

### Custom Base Directory
```bash
bash Scripts/quick_search.sh --base-dir /path/to/project 53394
```

---

## HPC Usage (SLURM)

### Single Sample
```bash
sbatch Scripts/submit_quick_search.sbatch 53394
```

### Multiple Samples
```bash
sbatch Scripts/submit_quick_search.sbatch 53394 53395 53396
```

### From Sample List File
```bash
sbatch Scripts/submit_quick_search.sbatch --sample-list Config/samples.txt
```

### Auto-detect All Samples
```bash
sbatch Scripts/submit_quick_search.sbatch --auto
```

### Custom Resources
Edit `submit_quick_search.sbatch` and modify:
```bash
#SBATCH --cpus-per-task=32     # More CPUs
#SBATCH --mem=128G             # More memory
#SBATCH --time=48:00:00        # Longer time
```

---

## Configuration

### Environment Variables
You can set these before running:

```bash
export THREADS=16                    # Number of threads
export BASE_DIR=/home/user/QIIME     # Project directory
export CONDA_BASE=/opt/conda         # Conda installation path
```

### HPC-Specific Settings
Edit `submit_quick_search.sbatch`:
- Job name, partition, resources
- Conda environment paths
- Database paths

---

## Sample List Format

Create `Config/samples.txt`:
```
# One sample ID per line
# Lines starting with # are comments
53394
53395
53396
```

---

## Output Structure

For each sample (e.g., `53394`):
```
Data/processed_data/bbduk_cleaned/
  ├── 53394_R1_bbduk.fastq.gz
  ├── 53394_R2_bbduk.fastq.gz
  └── 53394_bbduk_stats.txt

Data/processed_data/singlem_output/
  └── singlem_output_53394/
      ├── singlem_otu_table.csv
      └── singlem_profile.txt

Data/functional_analysis/methane_genes/
  ├── 53394_combined_methane_hits.txt
  ├── 53394_methane_summary.txt
  ├── 53394_filtered_methane_hits.txt
  ├── 53394_gene_counts.txt
  └── 53394_gene_rpkm.txt

Logs/
  └── 53394_quick_search.log
```

---

## Monitoring Jobs on HPC

### Check job status
```bash
squeue -u $USER
```

### View output in real-time
```bash
tail -f quick_search_<JOB_ID>.out
```

### Cancel a job
```bash
scancel <JOB_ID>
```

---

## Troubleshooting

### Environment activation fails
Check conda paths in `submit_quick_search.sbatch`:
```bash
export CONDA_BASE="/home/jiayi/miniconda3"
```

### Database not found
Set database paths explicitly:
```bash
export SINGLEM_METAPACKAGE="/path/to/metapackage"
export METHANE_DB="/path/to/diamond/db"
```

### BBMap adapters not found
The script will auto-detect or you can set:
```bash
export ADAPTERS="/path/to/adapters.fa"
```

### Low memory errors
Increase memory in SBATCH directives:
```bash
#SBATCH --mem=128G
```

---

## Pipeline Steps

1. **BBduk Quality Control** - Trim adapters and low-quality bases
2. **SingleM Taxonomic Profiling** - OTU/taxonomy classification  
3. **DIAMOND Functional Search** - Methane metabolism genes
4. **Summary Report** - Biological interpretation
5. **Gene Abundance** - Quantification and RPKM normalization
6. **Results Display** - Comprehensive output summary

---

## Help

```bash
bash Scripts/quick_search.sh --help
```
