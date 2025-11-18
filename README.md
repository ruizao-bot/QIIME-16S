# QIIME2 Metagenomic Pipeline

This repository contains a QIIME2 amplicon/metagenomic analysis pipeline with an opinionated directory layout and a checkpointing wrapper script at `Scripts/main.sh`.


## Directory structure (recommended)

```
QIIME/
├── README.md
├── Data/
│   ├── raw_data/             # Original FASTQ files and manifest.tsv
│   ├── processed_data/       # Imported QIIME2 artifacts (demux, trimmed .qza/.qzv)
│   └── reference_dbs/        # Pre-trained classifiers and reference artifacts (.qza)
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
  - For decontamination (step 6), create a metadata file with columns: `sample-id` and `control_status` (values: 'control' or 'not_control')
  - Example decontamination metadata:
    ```
    sample-id	control_status
    Sample1	not_control
    Sample2	not_control
    NegControl1	control
    NegControl2	control
    ```
- `Results/denoise_mode/` — DADA2 outputs (feature table, representative sequences, denoising stats and visualizations).
- `Results/cluster_mode/` — OTU clustering outputs (if using cluster mode).
- `Logs/checkpoints/` — the script writes `.checkpoint` files here to track completed steps; removing a checkpoint allows re-running that step.

## What the script does (high level)

- Provides a stepwise pipeline with checkpoints so long-running steps can be resumed.
- Supports two modes: `denoise` (DADA2) and `cluster` (OTU clustering from demux).
- Includes optional decontamination step using prevalence-based decontam algorithm with negative controls.
- Creates results under `Results/denoise_mode/` or `Results/cluster_mode/`.
- Writes logs to `Logs/qiime2_pipeline.log` and checkpoints to `Logs/checkpoints/`.
- Automatically detects and uses decontaminated files for downstream analysis when available.

## Step-by-step explanation (Steps 1–7)

This pipeline is split into seven interactive steps. Below is a short plain-language explanation of each step, what it needs, what it produces, and common notes/tips.

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

6) Decontamination (optional)
    - Purpose: Identify and remove contaminant sequences using negative control samples and the decontam algorithm (prevalence-based).
    - Inputs: 
      - Feature table from step 5 (`Results/<mode>_mode/table.qza`)
      - Metadata file with `control_status` column indicating which samples are controls (values: 'control' or 'not_control')
      - The script looks for metadata in: `DECONTAM_METADATA` env variable, `Data/metadata/decontam-metadata.tsv`, or `Data/metadata/metadata.tsv`
    - Outputs:
      - `${OUTPUT_DIR}/decontam-scores.qza` — decontam scores for all features
      - `${OUTPUT_DIR}/decontam-scores.qzv` — visualization of score distribution
      - `${OUTPUT_DIR}/decontam_scores_export/stats.tsv` — exported scores (columns: `#OTU ID`, `prev`, `p`)
      - `${OUTPUT_DIR}/decontam_scores_export/contaminant-ids.txt` — list of identified contaminants
      - `${OUTPUT_DIR}/table-no-contam.qza` — feature table with contaminants removed
      - `${OUTPUT_DIR}/rep-seqs-no-contam.qza` — filtered representative sequences
      - `${OUTPUT_DIR}/table-clean.qza` (optional) — table with contaminants AND negative controls removed
      - Checkpoint: `step6_decontamination`
    - Notes:
      - Decontam score threshold (default: 0.5) determines which features are contaminants (higher p-value = more likely contaminant)
      - If no contaminants are identified, the original table is used for downstream analysis
      - The script prompts whether to remove negative control samples after decontamination (recommended)
      - View `decontam-scores.qzv` at https://view.qiime2.org to assess score distribution and adjust threshold if needed
      - Environment variables for non-interactive mode: `DECONTAM_THRESHOLD` (default: 0.5), `REMOVE_CONTROLS` (default: y)

7) Taxonomic Classification
    - Purpose: Assign taxonomy to representative sequences using a classifier (pre-trained or custom-trained) and generate visualizations.
    - Inputs: 
      - Representative sequences from step 5 or 6 (decontaminated if applicable)
      - Feature table from step 5 or 6 (decontaminated if applicable)
      - A classifier artifact (e.g., `Data/reference_dbs/classifier.qza`) or option to build one using RESCRIPt
    - Outputs: 
      - `${OUTPUT_DIR}/taxonomy.qza` — taxonomic assignments for each feature
      - `${OUTPUT_DIR}/taxonomy.qzv` — taxonomy table visualization
      - `${OUTPUT_DIR}/taxa-bar-plots.qzv` — interactive taxonomic bar plots
      - `${OUTPUT_DIR}/exported-taxonomy/taxonomy.tsv` — taxonomy in TSV format
      - Checkpoint: `step8_taxonomic_classification`
    - Notes: 
      - You can supply a pre-trained classifier (fast), or build one (slower)
      - If building, the script supports RESCRIPt-based SILVA workflows with optional primer extraction
      - The script automatically uses decontaminated files if they exist from step 6
      - Taxa bar plots show relative abundances across samples

These explanations should help you decide which steps to run or re-run, and what inputs are required for each.
## Usage

### Clone the Repository

Open Terminal and navigate to where you want to store the pipeline, then clone:
```
git clone https://github.com/ruizao-bot/QIIME-16S.git
```
### Running guidance

Run the script from the project root (recommended) or let the script change to the project root if you modify it accordingly. Example invocations:

#### Run interactively (will prompt for choices)
```bash
bash Scripts/main.sh
```
#### Run in cluster mode
```bash
bash Scripts/main.sh -m cluster
```
#### Run in non-interactive/batch mode (for HPC)
```bash
bash Scripts/main.sh --non-interactive
```
#### Show pipeline status (which steps completed)
```bash
bash Scripts/main.sh --status
```
#### Remove checkpoint for step 6 (decontamination)
```bash
bash Scripts/main.sh --remove 6
```
#### Remove checkpoint for step 7 (taxonomic classification)
```bash
bash Scripts/main.sh --remove 7
```
#### Delete intermediate outputs from primer removal and denoising (interactive confirmation)
```bash
bash Scripts/main.sh --delete-intermediate
```
#### Remove all checkpoints and log file (start fresh)
```bash
bash Scripts/main.sh --clean
```

Options supported by the script include:
- `-h, --help` — show help
- `-c, --clean` — remove all checkpoints and logs
- `-s, --status` — show status of completed steps
- `-r, --remove <step#>` — remove a single checkpoint by step number (1-7)
- `-d, --delete-intermediate` — interactively delete intermediate files produced by primer removal and DADA2
- `-m, --mode <denoise|cluster>` — pipeline mode (default: denoise)
- `-e, --env <envname>` — conda environment name to use (default: qiime)
- `-b, --batch, --non-interactive` — run in non-interactive mode (uses env vars or defaults for prompts)

### Non-interactive mode environment variables

When running with `--non-interactive`, the script uses these environment variables:

- `START_STEP` — step number to start from (1-7)
- `PRIMER_CHOICE` — primer removal choice (1=515F/806R, 2=custom, 3=skip)
- `TRUNC_LEN_F` — forward read truncation length for DADA2 (default: 250)
- `TRUNC_LEN_R` — reverse read truncation length for DADA2 (default: 250)
- `N_THREADS` — number of threads for DADA2 (default: 0 = all available)
- `DECONTAMINATION_CHOICE` — decontamination (1=yes, 2=no/skip)
- `DECONTAM_METADATA` — path to metadata file with control_status column
- `DECONTAM_THRESHOLD` — score threshold for identifying contaminants (default: 0.5)
- `REMOVE_CONTROLS` — remove negative control samples after decontam (y/n, default: y)
- `BUILD_CLASSIFIER` — classifier option (silva, download, skip, custom)
- `CLASSIFIER_CHOICE` — numeric choice for classifier (1-4)
- `INSTALL_RESCRIPT` — install RESCRIPt if needed (y/n)
- `EXTRACT_PRIMERS` — extract primer region when building classifier (y/n)
- `PRIMER_F` — forward primer sequence (if EXTRACT_PRIMERS=y)
- `PRIMER_R` — reverse primer sequence (if EXTRACT_PRIMERS=y)
