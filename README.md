# Unified eDNA Metabarcoding Pipeline (Nanopore)

A single, consolidated Bash script implementing a full end‑to‑end eDNA metabarcoding workflow for Oxford Nanopore data. Inspired by the PIMENTA pipeline and recent protocols (Van der Vorst et al. 2024; Erkenswick et al. 2024), this script is organized into four clear phases with built‑in checks to skip completed steps, ensuring reproducibility and ease of reruns.

---

## Table of Contents

1. [Requirements](#requirements)  
2. [Installation](#installation)  
3. [Usage](#usage)  
4. [Pipeline Phases](#pipeline-phases)  
    1. [Phase 1 – Preprocessing & Quality Filtering](#phase-1--preprocessing--quality-filtering)  
    2. [Phase 2 – Per‑Sample Clustering & Consensus (NGSpeciesID)](#phase-2--per-sample-clustering--consensus-ngspeciesid)  
    3. [Phase 3 – Cross‑Sample Deduplication & Primer Trimming](#phase-3--cross-sample-deduplication--primer-trimming)  
    4. [Phase 4 – Taxonomic Assignment (BLAST)](#phase-4--taxonomic-assignment-blast)  
5. [Directory Layout](#directory-layout)  
6. [References](#references)

---

## Requirements

- **Bash** (WSL or Linux shell)  
- [Cutadapt](https://cutadapt.readthedocs.io)  
- [NGSpeciesID](https://github.com/<NGSpeciesID>)  
- [CD-HIT-EST](http://weizhong-lab.ucsd.edu/cd-hit/)  
- [MAFFT](https://mafft.cbrc.jp/)  
- NCBI BLAST+ (`blastn`/`megablast`)  
- Python 3.x (for filtering scripts/consensus.pl)  

## Installation

Clone this repository and ensure all dependencies are in your `$PATH`:

```bash
git clone https://github.com/joaquinauzenne/Metabarcoding.git
cd Metabarcoding
# Install dependencies, e.g. via conda:
# conda install -c bioconda cutadapt ngspeciesid cd-hit mafft blast
# or; Docker files, .yml, and .tar.gz provided in main directory
```

## Usage

```bash
bash pipeline.sh \
  /path/to/raw_fastq_dir \
  /path/to/trimmed_dir \
  /path/to/consensus_base_dir \
  [THREADS]
```

- `<RAW_FASTQ_DIR>`: demultiplexed `.fastq` files  
- `<TRIMMED_DIR>`: output of Phase 1  
- `<CONS_BASE_DIR>`: root output for consensus & clustering  
- `[THREADS]`: (optional) number of CPU threads (default: 8)

Each phase checks for existing output files (e.g., `_trimmed.fastq`, `*_consensus.fasta`) and skips if present.

---

## Pipeline Phases

### Phase 1 – Preprocessing & Quality Filtering

- **Input**: demultiplexed FASTQ per sample (Guppy‐basecalled).  
- **Poly‑A/T Removal**: trim homopolymer tails (`cutadapt --poly-a --poly-t`).  
- **Quality & Length**: enforce mean Q‐score, min/max length (`cutadapt -q <Q> -m <min> -M <max>`).  
- **Output**: `<sample>_trimmed.fastq`.

> *Tip:* wrap with `if [ -f "${sample}_trimmed.fastq" ]; then continue; fi` to skip.

### Phase 2 – Per‑Sample Clustering & Consensus (NGSpeciesID)

- **Clustering**: group reads ≥90% similar into clusters.  
- **Consensus**: align within each cluster (SPOA) and call majority-rule consensus.  
- **Output**: `consensus_dir/<sample>_consensus.fasta`, one or more sequences per sample.

> NGSpeciesID replaces manual CD‑HIT + MAFFT + consensus steps for initial clustering.

### Phase 3 – Cross‑Sample Deduplication & Primer Trimming

1. **Final Clustering**: pool all per-sample FASTAs → `cd-hit-est -c 0.98`.  
2. **Consensus per Cluster**: MAFFT alignment + `consensus.pl` to collapse nearly identical.  
3. **Primer Trimming**: loop over `PRIMER_PAIRS` array with Cutadapt (`-g` front, `-a` adapter).  
4. **Output**: primer-trimmed FASTA per marker in `Consensus/trimmed_by_marker/`.

### Phase 4 – Taxonomic Assignment (BLAST)

- **megablast**: query each final consensus against NCBI nt (`-max_target_seqs 10 -evalue 1e-5`).  
- **Filter**: keep hits ≥90% identity & ≥90% coverage; exclude uncultured/environmental sequences.  
- **Rank Adjustment**: derive lowest common taxonomic rank when multiple hits conflict.  
- **Output**: `blast_results.tsv` and summary reports (e.g., Krona chart).

---

## Directory Layout

**REWORKING**

---

## References

---

*Designed for reproducible, modular eDNA workflows on Nanopore data.*

