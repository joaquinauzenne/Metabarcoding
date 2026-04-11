# Metabarcoding of Environmental DNA from Soil Samples

A modular, HPC-ready pipeline for identifying mammal species from environmental DNA (eDNA) in soil samples using Oxford Nanopore long-read sequencing and MEGABLAST taxonomic assignment.

---

## Project Overview

This pipeline converts raw Nanopore basecalled reads into taxonomic identifications and interactive visualizations by:

1. Quality-trimming reads with Cutadapt (16S/12S mammalian primers)
2. Generating consensus sequences via NGSpeciesID clustering + MAFFT alignment (OTU mode) or dereplicate + chimera removal (ASV mode)
3. Optionally running DADA2 as an experimental parallel ASV inference branch
4. Assigning taxonomy with BLASTn against a curated Mammalia-only NCBI nucleotide database
5. Applying Lowest Common Ancestor (LCA) post-processing for robust taxonomic resolution
6. Visualizing species detections and IUCN range maps interactively via the `edna-viz` Shiny app

Designed to run on the Lonestar6 HPC system (TACC) using Apptainer containers, with local WSL development.

---

## Repository Structure

```text
Metabarcoding/
├── scripts/
│   ├── pipeline.sh              # Master 4-phase pipeline (OTU or ASV mode)
│   ├── lca_assignment.R         # LCA post-processing for BLAST results
│   ├── pipeline_dada2.R         # DADA2 experimental ASV branch (Nanopore)
│   ├── asv_otu_comparison.R     # Concordance comparison across pipeline branches
│   ├── setup.sh                 # Environment validation and directory init
│   ├── blastdb_init.sh          # SLURM: download full NCBI nt database
│   ├── mammalia_nt_db.sh        # SLURM: build Mammalia-only BLAST DB
│   ├── consensus.pl             # BioPerl MSA-to-consensus
│   ├── blast_taxonomy.sh        # Standalone taxonomy phase
│   └── ...                      # Additional standalone phase scripts
├── edna-viz/                    # Golem Shiny app (replaces VisualizationApp.qmd)
│   ├── DESCRIPTION
│   ├── R/
│   │   ├── app_config.R         # Golem config + global options
│   │   ├── app_server.R         # Module wiring
│   │   ├── app_ui.R             # shinydashboard layout
│   │   ├── run_app.R            # Entry point: ednaviz::run_app()
│   │   ├── mod_upload.R         # File upload + format detection
│   │   ├── mod_filter.R         # Threshold sliders + detection bar chart
│   │   ├── mod_map.R            # Leaflet map (IUCN + GBIF + eDNA sites)
│   │   ├── fct_detect_format.R  # Auto-detect web CSV vs CLI TSV
│   │   ├── fct_parse_blast.R    # Dual-format BLAST parser → canonical schema
│   │   ├── fct_lca.R            # In-app LCA assignment
│   │   └── fct_gbif_query.R     # GBIF occurrence fetcher
│   ├── data-raw/
│   │   └── prep_iucn_ranges.R   # One-time IUCN shapefile preprocessing
│   ├── inst/app/app_data/       # Pre-processed IUCN range polygons (.rds)
│   ├── tests/testthat/          # Unit tests
│   └── Dockerfile               # Production container
├── misc/
│   └── primers.txt              # 6 mammalian primer sets
├── examples/
│   └── 2ER9WHC1013-Alignment-Descriptions.csv   # Sample BLAST output (demo input)
├── Docker/Dockerfile            # Pipeline container
├── Apptainer/edna-pipeline.def  # HPC container definition
├── VisualizationApp.qmd         # Legacy Quarto/Shiny app (superseded by edna-viz/)
├── requirements.txt             # Conda environment dependencies
└── eDNA/                        # Generated pipeline outputs (excluded from git)
```

---

## Dependencies

All bioinformatics tools are packaged in the Apptainer container (`Apptainer/edna-pipeline.def`) or the conda environment (`requirements.txt`).

### Core Bioinformatics Tools

| Tool | Purpose |
| --- | --- |
| Cutadapt | Primer trimming and quality filtering |
| FastQC / NanoPlot / MultiQC | Read quality assessment |
| NGSpeciesID (v0.3.1) | Clustering nanopore reads into species groups |
| vsearch | Dereplication and chimera removal (ASV mode) |
| MAFFT | Multiple sequence alignment (OTU mode) |
| `consensus.pl` | BioPerl consensus from MSA (OTU mode) |
| CD-HIT-EST | Deduplication of consensus sequences |
| BLASTn / MEGABLAST | Taxonomic assignment |
| KronaTools (`ktImportText`) | Interactive taxonomy charts |
| DADA2 (R/Bioconductor) | Experimental ASV inference |

### R Packages (Visualization — `edna-viz`)

| Package | Purpose |
| --- | --- |
| golem / shiny / shinydashboard | App framework |
| leaflet / leaflet.extras | Interactive maps |
| sf | Spatial data (IUCN shapefiles) |
| rgbif | GBIF occurrence data |
| plotly / DT | Charts and tables |
| dplyr / tidyr / readr | Data manipulation |

---

## Setup

### 1. Container (HPC — Apptainer)

```bash
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

---

## BLAST Database Setup

### Option A: Full NCBI nt database (HPC — SLURM)

```bash
sbatch scripts/blastdb_init.sh
```

### Option B: Mammalia-only database (recommended)

```bash
sbatch scripts/mammalia_nt_db.sh
```

Builds a filtered database (Mammalia taxid 40674, excluding Hominidae) — dramatically reduces search space and runtime.

---

## Running the Pipeline

```bash
bash scripts/pipeline.sh <RAW_DIR> <TRIMMED_DIR> <CONS_BASE> <DB_DIR> \
  [THREADS] [RECLUSTER] [PRIMER_FILE] [MODE]
```

### Arguments

| Argument | Description | Default |
| --- | --- | --- |
| `RAW_DIR` | Raw `.fastq` / `.fastq.gz` files | — |
| `TRIMMED_DIR` | Output for trimmed reads | — |
| `CONS_BASE` | Base dir for consensus sequences | — |
| `DB_DIR` | BLAST DB path or remote name (`nt`) | — |
| `THREADS` | CPU threads | 8 |
| `RECLUSTER` | CD-HIT-EST reclustering at 100% (OTU mode) | `true` |
| `PRIMER_FILE` | Primer definitions file | built-in 16Smam |
| `MODE` | `otu` or `asv` — pipeline mode | `otu` |

### Example (OTU mode — backward compatible)

```bash
bash scripts/pipeline.sh \
  eDNA/2024 eDNA/trimmed eDNA/consensus scripts/mammalia_nt_db \
  16 true misc/primers.txt otu
```

### Example (ASV mode)

```bash
bash scripts/pipeline.sh \
  eDNA/2024 eDNA/trimmed eDNA/consensus_asv scripts/mammalia_nt_db \
  16 true misc/primers.txt asv
```

### Pipeline Phases

**Phase 1 — Preprocess**: Cutadapt trims adapters (min quality 12, length 10–1100 bp). Identical for both modes.

**Phase 2 — QC & Merge**: FastQC, NanoPlot, MultiQC quality reports. Reads merged per barcode. Identical for both modes.

**Phase 3 — Consensus / ASV inference**:

- *OTU mode*: NGSpeciesID @90% → MAFFT → consensus.pl → CD-HIT-EST @100% → primer trim
- *ASV mode*: vsearch dereplicate → NGSpeciesID (`--abundance_ratio 0.01 --consensus_n_rounds 3`) → vsearch chimera removal → CD-HIT-EST @100% → primer trim

**Phase 4 — Taxonomy**: BLASTn megablast at ≥97% identity, ≥90% query coverage. Produces:

- `blast.tsv` — 15-column format for Krona charts
- `blast_results.tsv` — 10-column format for `edna-viz` Shiny app
- `taxonomy.tsv` — LCA-assigned taxonomy (`lca_assignment.R`)
- `krona.html` — interactive taxonomy chart

---

## DADA2 Experimental Branch

Runs DADA2 ASV inference in parallel on the same quality-filtered reads from Phase 1. **Experimental**: DADA2's error model was designed for Illumina; Nanopore compatibility depends on read quality.

```bash
Rscript scripts/pipeline_dada2.R eDNA/trimmed eDNA/dada2_output 16
```

Inspect `eDNA/dada2_output/error_rates.pdf` before trusting results. A clear decreasing trend with quality score indicates the model learned useful error structure.

---

## LCA Post-Processing

```bash
Rscript scripts/lca_assignment.R blast_results.tsv taxonomy.tsv
```

Reads 10-column BLAST TSV, applies top-hit LCA logic within a 2% bitscore margin, and assigns `lca_rank` (species | genus | family_or_higher) with `confidence` (high | medium | low).

---

## ASV vs OTU Comparison

```bash
Rscript scripts/asv_otu_comparison.R \
  otu_taxonomy.tsv asv_taxonomy.tsv dada2_taxonomy.tsv \
  comparison_output/
```

Produces `species_concordance.csv` with flags: `concordant`, `asv_only`, `dada2_only`, `otu_only`. Upload this file to the **Comparison** tab in `edna-viz`.

---

## edna-viz Shiny App

Replaces the legacy `VisualizationApp.qmd`. A production-quality Golem app with:

- Dual-format BLAST upload (web CSV or CLI TSV — auto-detected)
- LCA taxonomy assignment in-app
- Threshold sliders with flagging (not hard-dropping) of below-threshold hits
- Leaflet map: IUCN range polygons + GBIF occurrence points + eDNA sampling sites

### Run locally

```r
# Install the package
devtools::install_local("edna-viz")

# Run the app
ednaviz::run_app()
```

### Prepare IUCN spatial data (one-time)

```r
# 1. Download TERRESTRIAL_MAMMALS.shp from:
#    https://www.iucnredlist.org/resources/spatial-data-download
# 2. Update SHAPEFILE_PATH in the script, then run:
source("edna-viz/data-raw/prep_iucn_ranges.R")
```

### Run tests

```r
devtools::test("edna-viz")
```

### Docker deployment

```bash
cd edna-viz
docker build -t edna-viz .
docker run -p 3838:3838 edna-viz
# Open http://localhost:3838
```

---

## Primer Reference

Defined in `misc/primers.txt` (format: `name:forward:reverse`):

| Name | Target | Amplicon Size |
| --- | --- | --- |
| 16Smam | 16S rRNA (mammal) | 150–170 bp |
| MiMammal | 12S rRNA | ~180 bp |
| Mammal16S | 16S rRNA | variable |
| Taberlat | 16S rRNA | ~60 bp |
| Riaz12S-V5 | 12S V5 | ~107 bp |
| Franklin360bp | 16S rRNA | ~330 bp |

---

## Output Files

| File / Directory | Description |
| --- | --- |
| `eDNA/trimmed/` | Primer-trimmed FASTQs, QC reports |
| `eDNA/consensus/per_barcode/` | Per-sample consensus FASTAs |
| `eDNA/consensus/combined/final_consensus.fasta` | Deduplicated consensus sequences |
| `eDNA/consensus/combined/taxonomy/blast_results.tsv` | 10-column BLAST (edna-viz app input) |
| `eDNA/consensus/combined/taxonomy/taxonomy.tsv` | LCA-assigned taxonomy table |
| `eDNA/consensus/combined/taxonomy/krona.html` | Interactive Krona chart |
| `eDNA/dada2_output/` | DADA2 ASVs, abundance table, error plots |

---

## FAIR Data Principles

- Pipeline scripts are fully parameterized and version-controlled
- Tool versions are logged by `setup.sh` at runtime
- Raw sequencing data is excluded from this repository (generated outputs in `eDNA/`)
- Container definitions (`Docker/`, `Apptainer/`) ensure computational reproducibility
- Example data provided in `examples/` for testing the visualization app

*Note: This repository has been modified recently with the usage of agentic AI, specifically Claude Code, to ameliorate its structure. This repository remains in progress under my personal direction.*
