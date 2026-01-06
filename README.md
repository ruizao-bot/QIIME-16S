# QIIME2 Metagenomic Pipeline

This repository contains a QIIME2 amplicon/metagenomic analysis pipeline with an opinionated directory layout and a checkpointing wrapper script at `Scripts/main.sh`.


## Directory structure (recommended)

```
QIIME/
├── README.md
├── Data/
│   ├── raw_data/             
│   │   ├── shotgun/          # Shotgun metagenomic FASTQ files and manifest.tsv
│   │   └── manifest.tsv      # Amplicon FASTQ files and manifest.tsv
│   ├── processed_data/       # Imported QIIME2 artifacts (demux, trimmed .qza/.qzv)
│   ├── reference_dbs/        # Pre-trained classifiers and reference artifacts (.qza)
│   ├── fastqc_raw_reports/   # FastQC reports for raw data
│   └── fastqc_trimmed_reports/ # FastQC reports for trimmed data
├── Scripts/
│   ├── main.sh               # Amplicon analysis pipeline (16S/18S)
│   └── shotgun.sh            # Shotgun metagenomic analysis pipeline
├── Results/
│   ├── denoise_mode/         # DADA2 outputs (table.qza, rep-seqs.qza, .qzv)
│   └── cluster_mode/         # Clustering-mode outputs
├── Logs/
│   ├── checkpoints/          # checkpoint files written by the script
│   ├── qiime2_pipeline.log   # amplicon pipeline run log
│   └── shotgun_pipeline.log  # shotgun pipeline run log
└── Config/
        ├── manifest_template.tsv
        └── metadata_template.tsv
```

## Notes about the layout

- `Data/raw_data/` — place your raw FASTQ and `manifest.tsv` here before running the pipeline. The manifest must follow QIIME2 format.
  - For amplicon data: place files in the root of `raw_data/`
  - For shotgun metagenomic data: place files in `raw_data/shotgun/` subdirectory
- `Data/processed_data/` — this is where the script writes imported/visualized demux artifacts (e.g. `demux-paired-end.qza`, `demux-paired-end.qzv`)
- `Data/reference_dbs/` — put any downloaded or trained classifier (e.g. `classifier.qza`) here; the script will look here first when classifying.
  - For shotgun pipeline: place host genome Bowtie2 index in `reference_dbs/host_genome/` for host DNA removal
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

---

## Shotgun Metagenomic Analysis Pipeline

The `shotgun.sh` script provides a comprehensive pipeline for analyzing shotgun metagenomic sequencing data using QIIME2's moshpit plugin. This pipeline is designed for whole-genome sequencing (WGS) data and performs assembly, binning, and functional annotation of microbial genomes.

### Environment Setup

The shotgun pipeline requires the `qiime2-moshpit` conda environment and the `qc_env` environment for FastQC. Ensure both environments are installed before running the pipeline.

### Directory Structure for Shotgun Data

Place your shotgun metagenomic FASTQ files in:
```
Data/raw_data/shotgun/
```

Create a manifest file at:
```
Data/raw_data/shotgun/manifest.tsv
```

The manifest should follow QIIME2's PairedEndFastqManifestPhred33V2 format.

### Shotgun Pipeline Steps (Steps 0-7c)

The shotgun pipeline consists of the following steps:

**Step 0: FastQC on Raw Data**
- **Purpose**: Quality assessment of raw sequencing reads
- **Inputs**: Raw FASTQ.gz files in `Data/raw_data/shotgun/`
- **Outputs**: FastQC HTML reports in `Data/fastqc_raw_reports/`
- **Environment**: `qc_env`
- **Notes**: Helps identify quality issues before processing. Review reports to assess read quality and determine if preprocessing is needed.

**Step 1: Import Data**
- **Purpose**: Import paired-end shotgun metagenomic reads into QIIME2 format
- **Inputs**: Manifest file (`Data/raw_data/shotgun/manifest.tsv`) listing sample IDs and paths to forward/reverse FASTQ files
- **Outputs**: `Data/processed_data/demux-paired-end.qza`
- **Notes**: The manifest must follow QIIME2's PairedEndFastqManifestPhred33V2 format. Script will fail if manifest is not found or malformed.

**Step 2: Quality Control (Trim Primers/Adapters)**
- **Purpose**: Remove primer sequences from reads using cutadapt
- **Inputs**: `demux-paired-end.qza` from Step 1
- **Outputs**: `Data/processed_data/trimmed-seqs.qza`
- **Primers Used**: 
  - Forward: 799F (AACMGGATTAGATACCCKG)
  - Reverse: 1193R (ACGTCATCCCCACCTTCC)
- **Parameters**: `--p-discard-untrimmed` (reads without primers are discarded), `--p-cores 4`
- **Notes**: Uses cutadapt for primer trimming. Reads without both primers are discarded to ensure only target amplicons are retained.

**Step 2b: FastQC on Trimmed Data**
- **Purpose**: Quality assessment of trimmed reads to verify primer removal
- **Inputs**: `trimmed-seqs.qza` from Step 2
- **Outputs**: FastQC reports in `Data/fastqc_trimmed_reports/`
- **Process**: Exports trimmed sequences to FASTQ format, runs FastQC, then switches back to qiime2-moshpit environment
- **Notes**: Compare with raw FastQC reports to confirm primer removal and assess remaining read quality.

**Step 3: Remove Host DNA**
- **Purpose**: Filter out host (e.g., human, mouse) DNA sequences to focus on microbial content
- **Inputs**: `trimmed-seqs.qza` from Step 2, Host genome Bowtie2 index at `Data/reference_dbs/host_genome/`
- **Outputs**: `Data/processed_data/host-removed-seqs.qza`
- **Parameters**: `--p-n-threads 16`
- **Notes**: If host genome database is not found, this step is skipped and the input file is copied to output. To enable host removal, build a Bowtie2 index of the host genome and place it in the reference directory.

**Step 4: Assemble Contigs**
- **Purpose**: Assemble short reads into longer contiguous sequences (contigs) using MEGAHIT assembler
- **Inputs**: `host-removed-seqs.qza` from Step 3
- **Outputs**: `Data/processed_data/assembled-contigs.qza`
- **Assembler**: MEGAHIT (memory-efficient, fast assembler for metagenomes)
- **Notes**: This is often a time-consuming step depending on data size. MEGAHIT is optimized for metagenomic data and handles complex communities well.

**Step 4b: Map Reads to Contigs**
- **Purpose**: Align original reads back to assembled contigs to calculate coverage depth (required for binning)
- **Inputs**: `assembled-contigs.qza` from Step 4, `host-removed-seqs.qza` from Step 3
- **Outputs**: 
  - `Data/processed_data/contigs-index.qza` (Bowtie2 index of contigs)
  - `Data/processed_data/mapped-reads.qza` (alignment maps)
- **Process**: 
  1. Index contigs using `qiime assembly index-contigs`
  2. Map reads to contigs using `qiime assembly map-reads`
- **Notes**: Coverage information is essential for binning in Step 5. Higher coverage contigs are more reliable.

**Step 5: Bin Contigs (MAG Generation)**
- **Purpose**: Group contigs into Metagenome-Assembled Genomes (MAGs) based on sequence composition and coverage patterns using MetaBAT2
- **Inputs**: `assembled-contigs.qza` from Step 4, `mapped-reads.qza` from Step 4b
- **Outputs**: 
  - `Data/processed_data/binned-contigs.qza` (MAGs)
  - `Data/processed_data/contig-map.qza` (mapping of contigs to bins)
  - `Data/processed_data/unbinned-contigs.qza` (contigs that couldn't be binned)
- **Algorithm**: MetaBAT2 (uses tetranucleotide frequency and coverage depth)
- **Default Parameters**: minContig=2500 bp, minCV=1.0, minCVSum=1.0, maxP=95%, minS=60, maxEdges=200, minClsSize=200000
- **Notes**: 
  - MAG quality depends on contig length and coverage. Low-coverage samples may produce few or no MAGs.
  - If no MAGs are formed, check: (1) contig length distribution, (2) read mapping percentage, (3) sample complexity
  - Common causes of failure: insufficient sequencing depth, highly diverse communities, poor assembly quality

**Step 6: Evaluate MAGs**
- **Purpose**: Assess completeness and contamination of MAGs using BUSCO (Benchmarking Universal Single-Copy Orthologs)
- **Inputs**: `binned-contigs.qza` from Step 5
- **Outputs**: `Data/processed_data/mags-evaluation.qzv` (visualization)
- **Metrics**: 
  - Completeness: percentage of expected single-copy genes present
  - Contamination: percentage of duplicated single-copy genes
- **Notes**: High-quality MAGs typically have >90% completeness and <5% contamination. View the .qzv file at https://view.qiime2.org to assess MAG quality.

**Step 7a: Predict Genes**
- **Purpose**: Identify protein-coding genes in MAGs using Prodigal gene prediction tool
- **Inputs**: `binned-contigs.qza` from Step 5
- **Outputs**: `Data/processed_data/predicted-genes.qza` (gene sequences)
- **Tool**: Prodigal (specialized for microbial gene prediction)
- **Notes**: Predicts both gene locations and translates to protein sequences for downstream annotation.

**Step 7b: Search Orthologs**
- **Purpose**: Find homologous genes in reference databases using DIAMOND (fast protein sequence alignment)
- **Inputs**: `predicted-genes.qza` from Step 7a
- **Outputs**: `Data/processed_data/orthologs.qza` (ortholog annotations)
- **Tool**: DIAMOND (100-20,000x faster than BLAST)
- **Notes**: Identifies known orthologs to enable functional annotation. DIAMOND is preferred over BLAST for large metagenomic datasets.

**Step 7c: Annotate MAGs**
- **Purpose**: Assign functional annotations to genes using eggNOG database (clusters of orthologous groups)
- **Inputs**: `orthologs.qza` from Step 7b
- **Outputs**: `Data/processed_data/mags-annotations.qza` (functional annotations)
- **Database**: eggNOG (evolutionary genealogy of genes: Non-supervised Orthologous Groups)
- **Annotations Include**: COG categories, GO terms, KEGG pathways, enzyme predictions
- **Notes**: Final step provides comprehensive functional profiles of MAGs, enabling metabolic pathway analysis and functional comparisons.

### Running the Shotgun Pipeline

#### Basic Usage
```bash
bash Scripts/shotgun.sh
```

#### Test Mode (for testing with small datasets)
```bash
bash Scripts/shotgun.sh --test
```

#### Clear All Outputs and Restart Pipeline
To remove all output files and re-run from the beginning:
```bash
rm -f Data/processed_data/demux-paired-end.qza \
      Data/processed_data/trimmed-seqs.qza \
      Data/processed_data/host-removed-seqs.qza \
      Data/processed_data/assembled-contigs.qza \
      Data/processed_data/contigs-index.qza \
      Data/processed_data/mapped-reads.qza \
      Data/processed_data/binned-contigs.qza \
      Data/processed_data/contig-map.qza \
      Data/processed_data/unbinned-contigs.qza \
      Data/processed_data/mags-evaluation.qzv \
      Data/processed_data/predicted-genes.qza \
      Data/processed_data/orthologs.qza \
      Data/processed_data/mags-annotations.qza
```

### Pipeline Features

- **Checkpointing**: Each step checks if output already exists and skips if present
- **Error Handling**: Pipeline stops on errors with informative messages
- **Logging**: All operations logged to `Logs/shotgun_pipeline.log` with timestamps
- **Environment Management**: Automatically switches between conda environments as needed

### Troubleshooting

**No MAGs formed during binning:**
- Check contig length: MetaBAT2 requires contigs ≥2500 bp by default
- Verify read mapping percentage: Low percentage (<1%) indicates assembly issues
- Assess sequencing depth: Deep sequencing improves assembly and binning
- Review assembly quality: Poor assembly produces fragmented contigs that cannot be binned

**Export contigs for inspection:**
```bash
qiime tools export \
    --input-path Data/processed_data/assembled-contigs.qza \
    --output-path exported_contigs/
```

**Check contig statistics:**
```bash
# Install seqkit if not available
conda install -c bioconda seqkit

# Get contig length statistics
seqkit stats exported_contigs/*.fa
```

---

## Amplicon Analysis Pipeline (16S/18S rRNA)



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
