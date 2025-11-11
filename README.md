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

## Usage

Run the script from the project root (recommended) or let the script change to the project root if you modify it accordingly. Example invocations:

```bash
# Run interactively (will prompt for choices)
bash Scripts/main.sh

# Run in cluster mode
bash Scripts/main.sh -m cluster

# Show pipeline status (which steps completed)
bash Scripts/main.sh --status

# Remove a specific checkpoint so a step will run again
bash Scripts/main.sh --remove step5_dada2_denoising

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

## Typical file flow

```
Data/raw_data/manifest.tsv -> Scripts/main.sh (import) -> Data/processed_data/demux-paired-end.qza
                                                         -> (visualize) -> Data/processed_data/demux-paired-end.qzv
                                                         -> (trim primers) -> Data/processed_data/demux-trimmed.qza
                                                         -> (denoise/cluster) -> Results/[mode]/table.qza, rep-seqs.qza, *.qzv
```

## Metadata usage

If you want group-aware visualizations (taxa bar plots, group-based diversity tests), create a metadata TSV and use it with QIIME2 commands such as `qiime taxa barplot`:

```bash
qiime taxa barplot \
    --i-table Results/denoise_mode/table.qza \
    --i-taxonomy Results/denoise_mode/taxonomy.qza \
    --m-metadata-file Data/metadata/sample-metadata.tsv \
    --o-visualization Results/denoise_mode/taxa-bar-plots.qzv
```

## Suggested next steps

- Ensure the `Data/` and `Results/` folders exist before running (or modify `Scripts/main.sh` to `cd` to project root and create them automatically).
- Add any download/training steps for reference classifiers under `Data/reference_dbs/`.
- Add your sample metadata (if needed) at `Data/metadata/sample-metadata.tsv` and use it with taxa barplots and diversity commands.

If you'd like, I can update `Scripts/main.sh` to automatically `cd` to the project root and create the common folders at startup — say the word and I'll apply that change and run a quick syntax check.