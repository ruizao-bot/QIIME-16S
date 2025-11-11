# QIIME2 Metagenomic Pipeline

This repository contains a QIIME2 amplicon/metagenomic analysis pipeline with an opinionated directory layout and a checkpointing wrapper script at `Scripts/main.sh`.


## Directory structure (recommended)

```
QIIME/
├── README.md
├── Data/
│   ├── raw_data/             # Original FASTQ files and manifest.tsv
│   ├── processed_data/       # Imported QIIME2 artifacts (demux, trimmed .qza/.qzv)
│   ├── reference_dbs/        # Pre-trained classifiers and reference artifacts (.qza)
│   └── metadata/             # Optional: sample metadata TSV(s) for group comparison
├── Scripts/
│   └── main.sh               # Main pipeline (run this script)
├── Results/
│   ├── denoise_mode/         # DADA2 outputs (table.qza, rep-seqs.qza, .qzv)
│   └── cluster_mode/         # Clustering-mode outputs
├── Logs/
│   ├── checkpoints/          # checkpoint files written by the script
│   └── qiime2_pipeline.log   # pipeline run log
└── Config/
        ├── manifest_template.tsv
        └── metadata_template.tsv
```

## Notes about the layout

- `Data/raw_data/` — place your raw FASTQ and `manifest.tsv` here before running the pipeline. The manifest must follow QIIME2 format.
- `Data/processed_data/` — this is where the script writes imported/visualized demux artifacts (e.g. `demux-paired-end.qza`, `demux-paired-end.qzv`)
- `Data/reference_dbs/` — put any downloaded or trained classifier (e.g. `classifier.qza`) here; the script will look here first when classifying.
- `Data/metadata/` — optional location for sample metadata TSVs (`sample-id` column required). The current `main.sh` does not automatically inject metadata into every downstream command; you can supply metadata when running commands like `qiime taxa barplot`.
- `Results/denoise_mode/` — DADA2 outputs (feature table, representative sequences, denoising stats and visualizations).
- `Logs/checkpoints/` — the script writes `.checkpoint` files here to track completed steps; removing a checkpoint allows re-running that step.

## What the script does (high level)

- Provides a stepwise pipeline with checkpoints so long-running steps can be resumed.
- Supports two modes: `denoise` (DADA2) and `cluster` (OTU clustering from demux).
- Creates results under `Results/denoise_mode/` or `Results/cluster_mode/`.
- Writes logs to `Logs/qiime2_pipeline.log` and checkpoints to `Logs/checkpoints/`.

## Step-by-step explanation (Steps 1–6)

This pipeline is split into six interactive steps. Below is a short plain-language explanation of each step, what it needs, what it produces, and common notes/tips.

1) Environment Setup
    - Purpose: Create or activate the Conda environment with QIIME2 and required plugins.
    - Inputs: None (the script will ask for the Conda env name or use the default `qiime`).
    - Outputs: A ready-to-use Conda environment and a checkpoint file `Logs/checkpoints/step1_environment_setup.checkpoint`.
    - Notes: If the environment already exists the script offers to reuse it or recreate it from the official QIIME2 amplicon YAML. This step verifies `qiime` is callable before proceeding.

2) Import Data
    - Purpose: Import your demultiplexed paired-end FASTQ files into a QIIME2 artifact so downstream plugins can read them.
    - Inputs: `Data/raw_data/manifest.tsv` (or `manifest.tsv` in project root). The manifest must follow QIIME2 format and use Phred33V2 if that is your data.
    - Outputs: `Data/processed_data/demux-paired-end.qza` and checkpoint `step2_import_data`.
    - Notes: The script fails early if the manifest is not found or malformed. Verify your sample IDs match your metadata file if you plan to use metadata later.

3) Visualize Demux Data
    - Purpose: Produce an interactive summary of per-sample sequence quality to choose truncation values for trimming/denoising.
    - Inputs: `Data/processed_data/demux-paired-end.qza` (from step 2).
    - Outputs: `Data/processed_data/demux-paired-end.qzv` and checkpoint `step3_visualize_demux`.
    - Notes: Open the `.qzv` at https://view.qiime2.org and inspect per-base quality plots to choose sensible truncation lengths for DADA2 (step 5).

4) Remove Primers/Adapters (optional)
    - Purpose: Trim primers/adapters using cutadapt if your reads still contain primer sequences.
    - Inputs: The demux artifact from step 2/3 and user-specified primer sequences (or choose the common 515F/806R set in the prompt).
    - Outputs: `Data/processed_data/demux-trimmed.qza` (and `.qzv`) and checkpoint `step4_remove_primers` (skips if user chooses to skip).
    - Notes: Only run this if primers are present; removing primers before denoising improves downstream accuracy. The script offers common primers or allows custom entry.

5) Denoising or Clustering
    - Purpose: Produce the feature table and representative sequences using either DADA2 (default, higher-resolution ASVs) or a clustering workflow (OTUs).
    - Inputs: `demux-trimmed.qza` if primers were removed, otherwise `demux-paired-end.qza`. User-provided truncation lengths for DADA2.
    - Outputs (denoise mode): `Results/denoise_mode/table.qza`, `Results/denoise_mode/rep-seqs.qza`, denoising stats and visualizations, plus checkpoint `step5_dada2_denoising`.
      Outputs (cluster mode): `Results/cluster_mode/table.qza` and `rep-seqs.qza`, with optional clustered variants and checkpoint `step5_cluster_from_demux`.
    - Notes: DADA2 is recommended for modern 16S workflows; clustering is available for legacy comparisons or specific workflows. This is typically the slowest step; ensure you have enough CPU/threads and memory.

6) Taxonomic Classification
    - Purpose: Assign taxonomy to representative sequences using a classifier (pre-trained or custom-trained) and generate visualizations.
    - Inputs: `Results/<mode>_mode/rep-seqs.qza` and a classifier artifact (e.g., `Data/reference_dbs/classifier.qza`) or the option to build one using RESCRIPt or your own FASTA+taxonomy.
    - Outputs: `${OUTPUT_DIR}/taxonomy.qza`, `${OUTPUT_DIR}/taxonomy.qzv`, `${OUTPUT_DIR}/taxa-bar-plots.qzv`, exported TSVs, and checkpoint `step6_taxonomic_classification`.
    - Notes: You can supply a pre-trained classifier (fast), or build one (slower). If you build, the script supports RESCRIPt-based SILVA workflows and keeping interactive primer extraction so the training region matches your amplicon.

These explanations should help you decide which steps to run or re-run, and what inputs are required for each.
## Usage

## Prerequisites

Before using this pipeline, ensure you have:This pipeline automates microbiome data analysis using QIIME2. It identifies and visualizes bacteria and other microorganisms in your samples through a simple, step-by-step process.## 📚 What is This?This repository contains a QIIME2 amplicon/metagenomic analysis pipeline with an opinionated directory layout and a checkpointing wrapper script at `Scripts/main.sh`.

- Git installed on your computer

- Basic familiarity with Terminal/Command Line


### Clone the Repository

Open Terminal and navigate to where you want to store the pipeline, then clone:

git clone https://github.com/ruizao-bot/QIIME-16S.git

## Running guidance
Run the script from the project root (recommended) or let the script change to the project root if you modify it accordingly. Example invocations:

# Run interactively (will prompt for choices)
bash Scripts/main.sh

# Run in cluster mode
bash Scripts/main.sh -m cluster

# Show pipeline status (which steps completed)
bash Scripts/main.sh --status

# Remove a specific checkpoint so a step will run again
bash Scripts/main.sh --remove 5

# Delete intermediate outputs from primer removal and denoising (interactive confirmation)
bash Scripts/main.sh --delete-intermediate

# Remove all checkpoints and log file (start fresh)
bash Scripts/main.sh --clean
```

Options supported by the script include:
- `-h, --help` — show help
- `-c, --clean` — remove all checkpoints and logs
- `-s, --status` — show status of completed steps
- `-r, --remove <checkpoint>` — remove a single checkpoint file (e.g. `step5_dada2_denoising`)
- `-d, --delete-intermediate` — interactively delete intermediate files produced by primer removal and DADA2
- `-m, --mode <denoise|cluster>` — pipeline mode
- `-e, --env <envname>` — conda environment name to use

## Metadata usage

If you want group-aware visualizations (taxa bar plots, group-based diversity tests), create a metadata TSV and use it with QIIME2 commands such as `qiime taxa barplot`:

```bash
qiime taxa barplot \
    --i-table Results/denoise_mode/table.qza \
    --i-taxonomy Results/denoise_mode/taxonomy.qza \
    --m-metadata-file Data/metadata/sample-metadata.tsv \
    --o-visualization Results/denoise_mode/taxa-bar-plots.qzv
```
