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
7. [Project Summary](#project-summary)

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

## Project Summary

### Unified eDNA Metabarcoding Pipeline (Nanopore)

A single consolidated Bash script that integrates the user’s eDNA metabarcoding workflow with improvements inspired by the PIMENTA pipeline and recent protocols (e.g., Van der Vorst et al. 2024, Erkenswick et al. 2024). The pipeline is organized into clearly labeled phases and implements best practices (skipping completed steps, robust Bash syntax for WSL compatibility, etc.) for ease of use and reproducibility.

---

### Phase 1: Read Preprocessing and Quality Filtering

> In Phase 1, raw Nanopore reads (already basecalled and demultiplexed by Guppy) are preprocessed. This includes poly-A/T tail removal and stringent quality/length filtering, while deferring any primer trimming to a later stage.

**Input:**  
Demultiplexed FASTQ files per sample (Guppy has trimmed adapters during basecalling). Reads from all targeted loci remain combined at this stage.

**Poly-A/T Trimming:**  
Use Cutadapt to trim homopolymer tails (common in Nanopore reads) – e.g., remove runs of AAAA... or TTTT... at read ends (github.com). This step parallels PIMENTA’s initial adapter and poly-A trimming (PIMENTA uses Guppy + PRINSEQ for similar effect (github.com).

**Quality Filtering:**  
Apply minimum quality thresholds (e.g., mean Q-score cutoff) and length range filters. For example, require reads within an expected amplicon length range and discard low-quality reads (Cutadapt’s `-q` and `-m`/`-M` options can enforce this). In PIMENTA, PRINSEQ performed an analogous trimming of low-quality bases and enforced length bounds (github.com).

**Output:**  
Cleaned FASTQ files per sample (postfix e.g. `_trimmed.fastq`). These files contain full-length amplicon reads with adapters removed and homopolymers trimmed, but still include primer sequences at their ends (primer removal is postponed to Phase 3).

> **Note:** Deferring primer trimming until after consensus generation helps preserve maximum overlap during clustering/alignment. This ensures that even reads of different marker loci remain separated by inherent sequence dissimilarity rather than by manual partitioning, and maximizes the data retained for consensus building (github.com).

Each sample’s preprocessing can be wrapped in a check – for example, the script can skip this phase for a sample if its `_trimmed.fastq` already exists (`if [ -f "$sample_trimmed.fastq" ]` guard). This way, reruns don’t redo Phase 1 unnecessarily.

---

### Phase 2: Per-Sample Clustering & Consensus (NGSpeciesID)

> In Phase 2, each sample’s filtered reads are clustered by sequence similarity and collapsed into high-accuracy consensus sequences per barcode (marker) using NGSpeciesID. This tool implements reference-free clustering (e.g., via isONclust or a similar algorithm) followed by multiple sequence alignment and consensus calling, analogous to the first clustering step of PIMENTA (github.com).

- **Clustering per Sample:**  
  NGSpeciesID groups reads into clusters based on sequence similarity (e.g., reads ≥90% similar might form a cluster). This is similar to other long-read metabarcoding approaches which cluster reads at ~90% identity to delineate putative Molecular OTUs (bmcgenomics.biomedcentral.com). Each cluster ideally represents a unique DNA barcode (species or taxon) in that sample.

- **Consensus Calling:**  
  For each cluster, NGSpeciesID builds a multiple alignment (using an internal aligner like SPOA) and calls a majority-rule consensus sequence. This consensus corrects random ONT errors, yielding sequences with accuracy often >99% (bmcgenomics.biomedcentral.com). Recent studies show that consensus barcodes from Nanopore reads can be as accurate and indel-free as Illumina sequences when such clustering and consensus algorithms are applied (bmcgenomics.biomedcentral.com).

- **Multiple Consensus per Sample:**  
  NGSpeciesID may output multiple consensus sequences per sample, one per distinct cluster/taxon detected (bmcbiol.biomedcentral.com). This is crucial for eDNA/metabarcoding, as a single sample (environmental sample) can contain DNA from many species.

**Output:**  
For each sample and each primer-defined locus, a FASTA file of consensus sequences (e.g., `sampleX_consensus.fasta`). Each FASTA header can carry the sample ID and cluster ID (e.g., `>SampleX_cluster1 readCount=50`).

> NGSpeciesID automates clustering and consensus in one step (bmcbiol.biomedcentral.com), so it replaces the manual CD-HIT → MAFFT → consensus steps that some pipelines use for initial clustering. If an expected per-sample consensus FASTA already exists, the script will skip rerunning NGSpeciesID for that sample (preventing redundant computation).

---

### Phase 3: Cross-Sample Consensus Deduplication (Final Clustering & Primer Trimming)

> Phase 3 takes all consensus sequences from Phase 2 (pooled across samples, for all barcode loci) and performs a final clustering and deduplication, followed by primer trimming. This step collapses identical or nearly identical sequences across samples and ensures each unique sequence (presumably representing a unique taxon DNA barcode) is represented once. We use a two-step approach inspired by the PIMENTA “reclustering” phase (github.com).

1. **All-vs-All Clustering (CD-HIT):**  
   Combine all per-sample consensus FASTAs into one file and run `cd-hit-est` at a high identity threshold (e.g., 98–100%). This clusters redundant consensus sequences across samples (github.com).

2. **Consensus per Cluster:**  
   For each cluster produced by CD-HIT, align the member sequences with MAFFT and build a final consensus sequence using the user’s `consensus.pl` script (github.com). This yields one representative sequence per cluster.

3. **Primer Removal (Cutadapt):**  
   Trim primer sequences from these final consensus sequences. Define a Bash array `PRIMER_PAIRS` listing each locus and its forward/reverse primer sequences, then run Cutadapt (`-g/--front` and `-a/--adapter`) for each pair. Non-matching sequences go to an “untrimmed” file, ensuring each sequence is captured by its correct primer set.

**Output:**  
Final primer-trimmed consensus FASTA files (or one combined FASTA) containing unique high-quality barcode sequences ready for taxonomic assignment.

> Throughout Phase 3, the script uses clear subdirectories (e.g., `Clustering/All_consensus.fasta`, `Clustering/CDHIT_clusters.clstr`, `Consensus/final_consensus.fasta`, `Consensus/trimmed_by_marker/`) to organize intermediate files. Skipping logic can be applied at each sub-step.

---

### Phase 4: Taxonomic Assignment (BLAST + Post-processing)

> In Phase 4, each unique consensus sequence is subjected to taxonomic identification using a BLAST search against NCBI’s nucleotide database, followed by filtering of results to ensure reliable assignments.

- **BLAST (megablast) Search:**  
  Use `blastn` in megablast mode (`-max_target_seqs 10 -evalue 1e-5`) to query each consensus sequence against NCBI nt (github.com).

- **Identity and Coverage Filters:**  
  Parse the BLAST output to retain only hits with ≥90% sequence identity **AND** ≥90% query coverage (github.com).

- **Exclude Unwanted Taxa:**  
  Filter out hits corresponding to environmental or unplaced sequences using user-provided exclusion lists (github.com).

- **Consensus Taxonomy & Rank Adjustment:**  
  Determine the lowest common taxonomic rank among top hits. If all top hits are the same species, assign species; if they differ at species but share genus, assign genus; and so on.

**Output:**  
A tabular report (and/or annotated FASTA) listing each consensus sequence and its assigned taxonomy, plus BLAST metrics. The script can also produce summaries or visualizations (e.g., Krona charts).

> The script can incorporate a skip check for BLAST: if a previous run’s results file exists (e.g., `blast_results.tsv`), the BLAST search will be skipped or results merged, avoiding redundant queries.


*Designed for reproducible, modular eDNA workflows on Nanopore data.*

