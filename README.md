# Metabarcoding of Environmental DNA from Soil Samples

A modular, HPC-ready Bash pipeline for identifying mammal species from environmental DNA (eDNA) in soil samples using Oxford Nanopore long-read sequencing and MEGABLAST taxonomic assignment.

---

## Project Overview

This pipeline converts raw Nanopore basecalled reads into taxonomic identifications by:

1. Quality-trimming reads with Cutadapt (16S mammalian primers)
2. Generating high-fidelity consensus sequences via NGSpeciesID clustering and MAFFT alignment
3. Deduplicating cross-barcode consensus sequences with CD-HIT-EST
4. Assigning taxonomy with BLASTn against a curated Mammalia-only NCBI nucleotide database
5. Visualizing species detections and IUCN range maps interactively

Designed to run on the Lonestar6 HPC system (TACC) using Apptainer containers, with a local WSL development path.

---

## Repository Structure

```text
Metabarcoding/
├── scripts/               # All pipeline and utility scripts
│   ├── pipeline.sh        # Master pipeline (all 4 phases)
│   ├── setup.sh           # Environment validation and directory init
│   ├── blastdb_init.sh    # SLURM job: download full NCBI nt database
│   ├── mammalia_nt_db.sh  # SLURM job: build Mammalia-only BLAST DB
│   ├── consensus.pl       # BioPerl MSA-to-consensus script
│   ├── preprocess_reads.sh
│   ├── run_consensus.sh
│   ├── blast_taxonomy.sh
│   ├── qc_nanoplot_multiqc.sh
│   └── ...                # Additional standalone phase scripts
├── iucnranges/            # R package for IUCN range data
│   └── R/iucnranges.R
├── misc/
│   └── primers.txt        # 6 mammalian primer sets
├── examples/
│   └── 2ER9WHC1013-Alignment-Descriptions.csv   # Sample BLAST output for VisualizationApp
├── Docker/                # Docker container definition
├── Apptainer/             # Apptainer/Singularity container definition
├── VisualizationApp.qmd   # Interactive Quarto/Shiny species range map app
├── requirements.txt       # Conda environment dependencies
└── eDNA/                  # Generated data (excluded from git — produced by pipeline)
```

---

## Dependencies

All bioinformatics tools are packaged in the Apptainer container (`Apptainer/edna-pipeline.def`) or the conda environment (`requirements.txt`).

### Core Tools

| Tool | Purpose |
| --- | --- |
| Cutadapt | Primer trimming and quality filtering |
| FastQC / NanoPlot / MultiQC | Read quality assessment |
| NGSpeciesID (v0.3.1) | Clustering nanopore reads into species groups |
| MAFFT | Multiple sequence alignment |
| `consensus.pl` | BioPerl consensus from MSA |
| CD-HIT-EST | Deduplication of consensus sequences |
| BLASTn / MEGABLAST | Taxonomic assignment |
| KronaTools (`ktImportText`) | Interactive taxonomy charts |

### R Packages (Visualization)

| Package | Purpose |
| --- | --- |
| shiny / quarto | Interactive app framework |
| rredlist | IUCN Red List API client |
| sf | Spatial data handling (shapefiles) |
| leaflet | Interactive maps |
| tidyverse | Data manipulation |

---

## Setup

### 1. Container (HPC — Apptainer)

```bash
# Build the Apptainer container from the definition file
apptainer build Apptainer/edna-pipeline.sif Apptainer/edna-pipeline.def
```

### 2. Conda (local / WSL)

```bash
conda create -n edna-pipeline --file requirements.txt
conda activate edna-pipeline
```

### 3. Initialize Working Environment

```bash
bash scripts/setup.sh <RAW_DIR> <DB_DIR> [THREADS] [SETTINGS_FILE]
```

This validates all dependencies, checks executable availability, creates output directories, and logs tool versions for reproducibility. Optionally generates SLURM job scripts if a settings file is provided.

---

## BLAST Database Setup

### Option A: Full NCBI nt database (HPC — SLURM)

Submit the SLURM job to download the complete nucleotide database:

```bash
sbatch scripts/blastdb_init.sh
```

This downloads `nt` via `update_blastdb.pl`, syncing to `$WORK/Metabarcoding/blastdb/nt/`.

### Option B: Mammalia-only database (recommended for efficiency)

Build a filtered database containing only Mammalia sequences (taxid 40674, excluding Hominidae):

```bash
sbatch scripts/mammalia_nt_db.sh
```

This downloads all `nt_euk` volumes, filters to Mammalia using `blastdbcmd -taxids 40674`, and builds a custom database as `nt_euk_mammalia.*`. This dramatically reduces search space and runtime.

---

## Running the Pipeline

```bash
bash scripts/pipeline.sh <RAW_DIR> <TRIMMED_DIR> <CONS_BASE> <DB_DIR> [THREADS] [RECLUSTER] [PRIMER_FILE]
```

### Arguments

| Argument | Description |
| --- | --- |
| `RAW_DIR` | Directory containing raw `.fastq` / `.fastq.gz` files |
| `TRIMMED_DIR` | Output directory for trimmed reads |
| `CONS_BASE` | Base output directory for consensus sequences |
| `DB_DIR` | Path to BLAST database directory, or remote DB name (e.g. `nt`) |
| `THREADS` | Number of CPU threads (default: 8) |
| `RECLUSTER` | `true`/`false` — run CD-HIT-EST reclustering at 100% identity (default: `true`) |
| `PRIMER_FILE` | Path to primer definitions file (default: built-in 16Smam primers) |

### Example

```bash
bash scripts/pipeline.sh \
  eDNA/2024 \
  eDNA/trimmed \
  eDNA/consensus \
  scripts/mammalia_nt_db \
  16 true misc/primers.txt
```

### Pipeline Phases

**Phase 1 — Preprocess**: Cutadapt trims adapter/primer sequences (min quality 12, length 10–1100 bp). Output: `TRIMMED_DIR/`.

**Phase 2 — QC & Merge**: FastQC, NanoPlot, and MultiQC generate quality reports. Reads are merged per barcode for downstream analysis. Output: `TRIMMED_DIR/fastqc_results_trimmed/`, `nanoplot_results_trimmed/`, `multiqc_report_trimmed/`, `merged_fastq/`.

**Phase 3 — Consensus**: NGSpeciesID clusters reads at 90% identity. Each cluster is aligned with MAFFT and reduced to a consensus with `consensus.pl`. Optional CD-HIT-EST reclustering at 100% removes cross-barcode duplicates. Primers are trimmed from final consensus FASTAs. Output: `CONS_BASE/per_barcode/`, `CONS_BASE/combined/`, `CONS_BASE/by_marker/`.

**Phase 4 — Taxonomy**: BLASTn assigns taxonomy at ≥97% identity and ≥90% query coverage. Results are generated per barcode and for the combined pool. Krona charts are built for interactive exploration. Output: per-barcode CSVs and Krona `.html` files.

---

## Primer Reference

Defined in `misc/primers.txt` (format: `name:forward:reverse`):

| Name | Target | Amplicon Size |
| --- | --- | --- |
| 16Smam | 16S rRNA (mammal) | 150–170 bp |
| MiMammal | 12S rRNA | ~180 bp |
| Mammal16S | 16S rRNA | variable |
| Taberlat | 16S rRNA | ~60 bp |
| Riaz12S-V5 | 12S V5 | ~10 bp |
| Franklin360bp | 16S rRNA | ~360 bp |

---

## Output Files

| File / Directory | Description |
| --- | --- |
| `eDNA/trimmed/` | Primer-trimmed FASTQs, QC reports |
| `eDNA/consensus/per_barcode/` | Per-sample consensus FASTAs, cluster metadata, Medaka polishing outputs |
| `eDNA/consensus/combined/final_consensus.fasta` | Deduplicated combined consensus sequences |
| `eDNA/taxonomy/*.csv` | BLAST hit tables (per sequence) |
| `eDNA/taxonomy/*.html` | Krona interactive taxonomy charts |

---

## Visualization App

`VisualizationApp.qmd` is a Quarto/Shiny web application that:

1. Accepts an uploaded BLAST Alignment Descriptions CSV (see `examples/` for a sample)
2. Filters hits by Query Cover > 90% and Per. Identity > 90%
3. Ranks species by detection frequency
4. Fetches IUCN Red List range shapefiles via the `iucnranges` R package
5. Renders an interactive Leaflet map with toggleable species ranges

### Running the App

```r
# In R or RStudio, with the IUCN API key set:
Sys.setenv(IUCN_REDLIST_KEY = "your_api_key_here")
quarto::quarto_serve("VisualizationApp.qmd")
```

An IUCN Red List API key is required. Register at <https://api.iucnredlist.org>.

A sample input CSV is provided in `examples/2ER9WHC1013-Alignment-Descriptions.csv`.

---

## iucnranges R Package

The `iucnranges/` directory contains a lightweight R package for IUCN spatial data:

- `rl_init(key)` — initialize API key
- `fetch_species_ranges(species, path)` — download range shapefiles
- `combine_ranges(sf_list)` — merge into single sf object
- `crop_to_study_area(ranges_sf, area_sf)` — clip to study region
- `plot_range_map(ranges_sf, area_sf)` — Leaflet map
- `get_top_species(summary_df, barcode, n)` — filter top N detections

Install locally with `devtools::install_local("iucnranges")`.

---

## FAIR Data Principles

- Pipeline scripts are fully parameterized and version-controlled
- Tool versions are logged by `setup.sh` at runtime
- Raw sequencing data is excluded from this repository (generated outputs in `eDNA/`)
- Container definitions (`Docker/`, `Apptainer/`) ensure computational reproducibility
- Example data provided in `examples/` for testing the visualization app
