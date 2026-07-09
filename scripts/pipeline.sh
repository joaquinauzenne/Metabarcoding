#!/usr/bin/env bash
set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Master eDNA pipeline (all phases in one script):
#   1) Preprocess        (Cutadapt trimming → *_trimmed.fastq.gz)
#   2) QC & merge        (FastQC, NanoPlot, MultiQC, per-barcode merge)
#   3) Consensus         (NGSpeciesID @90%, optional reclustering @100%, MAFFT + consensus.pl, primer trimming)
#   4) Taxonomy          (BLASTn → per-barcode + combined, Krona charts)
#
# Usage:
#   bash pipeline.sh \
#        <RAW_DIR> <TRIMMED_DIR> <CONS_BASE> <DB_DIR> \
#        [THREADS] [RECLUSTER] [PRIMER_FILE] [MODE] [MODULES]
#
#   RAW_DIR      Directory with raw .fastq(.gz) files
#   TRIMMED_DIR  Output dir for trimmed FASTQs
#   CONS_BASE    Base output dir for consensus step
#   DB_DIR       BLAST DB directory (local nt*)
#   THREADS      (optional) Number of threads [default: 8]
#   RECLUSTER    (optional) true|false — whether to recluster combined consensus [default: true]
#   PRIMER_FILE  (optional) Path to a file of primers (name:fwd:rev per line);
#                if omitted, uses built-in 16Smam primer
#   MODE         (optional) otu|asv — pipeline mode [default: otu]
#                  otu: NGSpeciesID @90% → MAFFT → consensus.pl → CD-HIT-EST @100% (existing behavior)
#                  asv: vsearch dereplicate → NGSpeciesID (high sensitivity) → chimera removal → CD-HIT-EST @100%
#   MODULES      (optional) Comma-separated phases to run: Phase1,Phase2,Phase3,Phase4
#                Default "all" — runs every phase. Use e.g. "Phase3,Phase4" to resume
#                from consensus without re-trimming. Mirrors PIMENTA's Runmodules param.
# ─────────────────────────────────────────────────────────────────────────────

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Load pipeline configuration (single source of truth) ─────────────────────
# All tunable parameters live in config/pipeline.yaml. Override the path with
# PIPELINE_CONFIG=/path/to/config.yaml if needed. Parsed into CFG_* shell vars.
PIPELINE_CONFIG="${PIPELINE_CONFIG:-$SCRIPTS_DIR/../config/pipeline.yaml}"
if [[ ! -f "$PIPELINE_CONFIG" ]]; then
  printf "[ERROR] Pipeline config not found: %s\n" "$PIPELINE_CONFIG" >&2
  printf "        Set PIPELINE_CONFIG or create config/pipeline.yaml\n" >&2
  exit 1
fi
# shellcheck source=lib/parse_yaml.sh
source "$SCRIPTS_DIR/lib/parse_yaml.sh"
eval "$(parse_yaml "$PIPELINE_CONFIG" CFG_)"
printf "[INFO] Loaded config: %s\n" "$PIPELINE_CONFIG"

# ─────────────────────────────────────────────────────────────────────────────
# Phase 1 — Adapter trimming and quality filtering (Cutadapt).
# Parameters match PIMENTA's Prinseq settings: Q>=12, length 10-1100 bp (§2.2.1).
# We use Cutadapt instead of Prinseq for Nanopore-specific poly-A/T removal and
# correct handling of gzipped per-barcode FASTQ files already demultiplexed by MinKNOW.
#
# Function: preprocess_reads
# Usage   : preprocess_reads <input_fastq_dir> <output_trimmed_dir>
# ─────────────────────────────────────────────────────────────────────────────
preprocess_reads() {
  local INPUT_DIR="${1%/}"
  local OUTPUT_DIR="${2%/}"
  local MIN_QUAL="$CFG_cutadapt_min_qual" MIN_LEN="$CFG_cutadapt_min_len" MAX_LEN="$CFG_cutadapt_max_len"

  mkdir -p "$OUTPUT_DIR"
  printf "[INFO] Phase 1: preprocessing reads -> %s\n" "$OUTPUT_DIR"
  printf "[INFO]   Quality cutoff: %s\n" "$MIN_QUAL"
  printf "[INFO]   Length range  : %s-%s bp\n" "$MIN_LEN" "$MAX_LEN"

  find "$INPUT_DIR" -type f \( -name '*.fastq' -o -name '*.fastq.gz' \) | \
  while IFS= read -r IN; do
    local REL="${IN#"$INPUT_DIR"/}"
    local base="${REL%.fastq.gz}"; base="${base%.fastq}"
    local OUT="$OUTPUT_DIR/${base}_trimmed.fastq.gz"
    mkdir -p "$(dirname "$OUT")"

    if [[ -s "$OUT" ]]; then
      printf "[SKIP] %s (already exists)\n" "$REL"
      continue
    fi

    printf "[INFO] Trimming & filtering: %s\n" "$REL"
    python -m cutadapt \
      --poly-a --trim-n --max-n "$CFG_cutadapt_max_n" \
      --minimum-length "$MIN_LEN" --maximum-length "$MAX_LEN" \
      --quality-cutoff "$MIN_QUAL" --cores 0 \
      -o "$OUT" "$IN"

    if [[ -s "$OUT" ]]; then
      printf "[OK]   Trimmed: %s_trimmed.fastq.gz\n" "$base"
    else
      printf "[ERROR] No output for %s\n" "$REL" >&2
      rm -f "$OUT"
    fi
  done

  printf "[OK]   Phase 1 complete.\n"
}

# ─────────────────────────────────────────────────────────────────────────────
# Phase 2 — QC visualization and per-barcode merging.
# FastQC + NanoPlot provide Nanopore-specific quality diagnostics (base quality,
# read length distribution, N50). MultiQC aggregates across all barcodes.
# PIMENTA skips this step; we retain it for QC reporting required by the ELN.
#
# Function: quality_check_merge
# Usage   : quality_check_merge <trimmed_dir> <threads>
# ─────────────────────────────────────────────────────────────────────────────
quality_check_merge() {
  local TRIMMED_DIR="${1%/}"
  local THREADS="${2:-8}"

  local FASTQC_OUT="$TRIMMED_DIR/fastqc_results_trimmed"
  local MERGED_DIR="$TRIMMED_DIR/merged_fastq"
  local NANOPLOT_OUT="$TRIMMED_DIR/nanoplot_results_trimmed"
  local MULTIQC_OUT="$TRIMMED_DIR/multiqc_report_trimmed"

  mkdir -p "$FASTQC_OUT" "$MERGED_DIR" "$NANOPLOT_OUT" "$MULTIQC_OUT"

  # FastQC
  if [[ -n "$(ls -A "$FASTQC_OUT")" ]]; then
    printf "[SKIP] FastQC outputs already present.\n"
  else
    printf "[INFO] Running FastQC on all *_trimmed.fastq(.gz)...\n"
    find "$TRIMMED_DIR" -type f \( -name "*_trimmed.fastq" -o -name "*_trimmed.fastq.gz" \) | \
    while IFS= read -r fq; do
      printf "[INFO] FastQC on %s\n" "$fq"
      fastqc --threads "$THREADS" --outdir "$FASTQC_OUT" "$fq" \
        && printf "[OK]   FastQC succeeded for %s\n" "$fq" \
        || printf "[WARN] FastQC FAILED for %s\n" "$fq" >&2
    done
  fi

  # Merge by barcode
  printf "[INFO] Merging trimmed FASTQs by barcode into %s...\n" "$MERGED_DIR"
  find "$TRIMMED_DIR" -type f \( -name "*_trimmed.fastq" -o -name "*_trimmed.fastq.gz" \) \
    -printf '%h\n' | sort -u | while IFS= read -r BDIR; do
    local SAMPLE_TAG="${BDIR#${TRIMMED_DIR}/}"
    SAMPLE_TAG="${SAMPLE_TAG//\//_}"
    local MERGED_FASTQ="$MERGED_DIR/${SAMPLE_TAG}_merged.fastq"

    printf "[INFO]   Merging %s -> %s\n" "$BDIR" "$MERGED_FASTQ"
    {
      find "$BDIR" -type f \( -name "*_trimmed.fastq" -o -name "*_trimmed.fastq.gz" \) \
        | sort | \
        while IFS= read -r in; do
          [[ "$in" == *.gz ]] && zcat "$in" || cat "$in"
        done
    } > "$MERGED_FASTQ"
  done

  # NanoPlot
  printf "[INFO] Running NanoPlot on merged FASTQs...\n"
  find "$MERGED_DIR" -type f -name "*_merged.fastq" | \
  while IFS= read -r mfq; do
    local SAMPLE="$(basename "$mfq" _merged.fastq)"
    printf "[INFO] NanoPlot on %s\n" "$mfq"
    NanoPlot --fastq "$mfq" \
             --outdir "$NANOPLOT_OUT/$SAMPLE" \
             --threads "$THREADS" \
             --N50 \
             --plots hex dot kde \
      && printf "[OK]   NanoPlot complete for %s\n" "$mfq" \
      || printf "[WARN] NanoPlot FAILED for %s\n" "$mfq" >&2
  done

  # MultiQC
  printf "[INFO] Aggregating FastQC with MultiQC...\n"
  multiqc "$FASTQC_OUT" --outdir "$MULTIQC_OUT" \
    && printf "[OK]   MultiQC report: %s/multiqc_report.html\n" "$MULTIQC_OUT" \
    || printf "[WARN] MultiQC failed\n" >&2

  printf "[OK]   QC & merge complete.\n"
  printf "       FastQC:    %s\n" "$FASTQC_OUT"
  printf "       Merged:    %s\n" "$MERGED_DIR"
  printf "       NanoPlot:  %s\n" "$NANOPLOT_OUT"
  printf "       MultiQC:   %s/multiqc_report.html\n" "$MULTIQC_OUT"
}

# ─────────────────────────────────────────────────────────────────────────────
# Phase 3 — Clustering and consensus generation.
#
# [OTU] Per-sample clustering and consensus (NGSpeciesID).
# NGSpeciesID clusters Nanopore reads at >=90% similarity then runs Medaka polishing,
# producing >99%-accurate consensus sequences from error-prone raw reads.
# This replaces PIMENTA's two-step PRINSEQ->CD-HIT-EST->MAFFT->consensus.pl chain
# with a purpose-built Nanopore tool (Paper §3.1: "clustering now done at 93%
# rather than 80% identity" -- we use NGSpeciesID's internal 90% threshold).
# Per-sample processing mirrors PIMENTA Phase 2 ("Analysis per sample").
#
# [OTU recluster] Cross-sample deduplication at 100% identity.
# Mirrors PIMENTA Phase 3 ("Analysis of all samples combined"): per-barcode consensus
# sequences are pooled and reclustered at 100% identity to merge identical OTUs
# sequenced in multiple barcodes. MAFFT L-INS-i + consensus.pl at 30% threshold
# (same as PIMENTA's RandSize=100, -t 30 settings; Paper §2.2.1).
#
# [ASV] Dereplicate -> chimera removal -> NGSpeciesID (high sensitivity).
# vsearch --derep_fulllength collapses exact-duplicate reads and discards singletons,
# reducing noise before clustering. uchime3_denovo removes chimeric sequences.
# NGSpeciesID runs with abundance_ratio=0.01 to resolve rare variants.
# Parallels the OTU->ASV methodological transition from PIMENTA-style OTU clustering
# to ASV-level resolution (Callahan et al. 2017 DADA2; PIMENTA §3.1 discussion).
#
# Function: run_consensus
# Usage   : run_consensus <merged_fastq_dir> <output_base_dir> <threads> <true|false> [primer_file] [mode]
# ─────────────────────────────────────────────────────────────────────────────
run_consensus() {
  local MERGED_DIR="${1%/}"
  local OUT_BASE="${2%/}"
  local THREADS="${3:-8}"
  local RECLUSTER_FLAG="${4:-true}"
  local PRIMER_FILE="${5:-}"
  local MODE="${6:-otu}"

  [[ "$RECLUSTER_FLAG" != "true" ]] && RECLUSTER=false || RECLUSTER=true

  # load primers
  declare -a PRIMERS
  if [[ -n "$PRIMER_FILE" && -f "$PRIMER_FILE" ]]; then
    mapfile -t PRIMERS < "$PRIMER_FILE"
    printf "[INFO] Loaded primers from %s\n" "$PRIMER_FILE"
  else
    PRIMERS=( "16Smam:CGGGTTGGGGTGACCTCGGA:GCTGTTATCCCTAGGGTAACT" )
    printf "[INFO] Using default primer: %s\n" "${PRIMERS[0]}"
  fi

  local PER_BARCODE_OUT="$OUT_BASE/per_barcode"
  local COMBINED_OUT="$OUT_BASE/combined"
  local ALL_CONS="$COMBINED_OUT/all_consensus.fasta"
  local FINAL_CONS_FASTA="$COMBINED_OUT/final_consensus.fasta"
  local MARKER_OUT="$OUT_BASE/by_marker"
  local TAXON_DIR="$OUT_BASE/../taxonomy"
  local CONSENSUS_PL="$SCRIPTS_DIR/consensus.pl"

  mkdir -p "$PER_BARCODE_OUT" "$COMBINED_OUT" "$MARKER_OUT" "$TAXON_DIR"

  # 3a) Per-barcode clustering & consensus
  if [[ "$MODE" == "asv" ]]; then
    printf "\n[INFO] Phase 3 [ASV]: per-sample dereplicate + NGSpeciesID (high sensitivity)\n"
    shopt -s nullglob
    mapfile -t FQS < <(printf '%s\n' "$MERGED_DIR"/*_merged.fastq)
    (( ${#FQS[@]} )) || { printf "[ERROR] No merged FASTQ in %s\n" "$MERGED_DIR" >&2; exit 1; }
    set -x
    for fq in "${FQS[@]}"; do
      local sample="$(basename "$fq" _merged.fastq)"
      local outdir="$PER_BARCODE_OUT/$sample"
      mkdir -p "$outdir"
      if compgen -G "${outdir}/medaka_cl_id_*/consensus.fasta" > /dev/null; then
        printf "[SKIP] [per-barcode ASV] %s\n" "$sample"
        continue
      fi

      # Dereplicate: collapse identical reads, track abundances, discard singletons
      printf "[INFO]   [ASV] Dereplicating: %s\n" "$sample"
      vsearch --derep_fulllength "$fq" \
              --output "$outdir/${sample}_derep.fasta" \
              --sizeout --minuniquesize "$CFG_vsearch_minuniquesize"

      # NGSpeciesID with higher sensitivity for ASV resolution
      NGSpeciesID --fastq "$fq" \
                  --q "$CFG_ngspeciesid_q" --ont --consensus --medaka \
                  --m "$CFG_ngspeciesid_m" --s "$CFG_ngspeciesid_s" --sample_size "$CFG_ngspeciesid_sample_size" \
                  --consensus_n_rounds "$CFG_ngspeciesid_consensus_n_rounds" \
                  --abundance_ratio "$CFG_ngspeciesid_abundance_ratio" \
                  --outfolder "$outdir"
    done
    set +x
  else
    printf "\n[INFO] Phase 3 [OTU]: per-sample clustering & consensus (NGSpeciesID @90%%)\n"
    shopt -s nullglob
    mapfile -t FQS < <(printf '%s\n' "$MERGED_DIR"/*_merged.fastq)
    (( ${#FQS[@]} )) || { printf "[ERROR] No merged FASTQ in %s\n" "$MERGED_DIR" >&2; exit 1; }
    set -x
    for fq in "${FQS[@]}"; do
      local sample="$(basename "$fq" _merged.fastq)"
      local outdir="$PER_BARCODE_OUT/$sample"
      mkdir -p "$outdir"
      if compgen -G "${outdir}/medaka_cl_id_*/consensus.fasta" > /dev/null; then
        printf "[SKIP] [per-barcode] %s\n" "$sample"
        continue
      fi
      NGSpeciesID --fastq "$fq" \
                  --q "$CFG_ngspeciesid_q" --ont --consensus --medaka \
                  --m "$CFG_ngspeciesid_m" --s "$CFG_ngspeciesid_s" --sample_size "$CFG_ngspeciesid_sample_size" \
                  --outfolder "$outdir"
    done
    set +x
  fi

  # 3b) Collapse + combined consensus / ASV deduplication
  printf "\n"
  if [[ "$MODE" == "asv" ]]; then
    printf "[INFO] Phase 3 [ASV]: chimera removal + CD-HIT-EST @100%% deduplication\n"
  else
    printf "[INFO] Phase 3 [OTU]: collapse @100%% + MSA & final consensus\n"
  fi

  if [[ -s "$FINAL_CONS_FASTA" ]]; then
    printf "[SKIP] Phase 3b (already done)\n"
  else
    # concatenate per-barcode consensus
    if [[ ! -s "$ALL_CONS" ]]; then
      printf "[INFO]   Building %s\n" "$ALL_CONS"
      : > "$ALL_CONS"
      while IFS= read -r -d '' f; do
        cat "$f" >> "$ALL_CONS"; printf "\n" >> "$ALL_CONS"
      done < <(find "$PER_BARCODE_OUT" -mindepth 2 -path '*/medaka_cl_id_*/consensus.fasta' -print0)
    fi

    if [[ "$MODE" == "asv" ]]; then
      # ASV: chimera removal → CD-HIT @100% → primer trim (no MAFFT/consensus.pl)
      local NOCHIM_FASTA="$COMBINED_OUT/all_consensus_nochim.fasta"
      local CHIMERA_FASTA="$COMBINED_OUT/chimeras.fasta"
      local CDHIT_OUT="$COMBINED_OUT/asv_representatives.fasta"

      printf "[INFO]   [ASV] Chimera removal -> %s\n" "$NOCHIM_FASTA"
      vsearch --uchime3_denovo "$ALL_CONS" \
              --nonchimeras "$NOCHIM_FASTA" \
              --chimeras "$CHIMERA_FASTA"

      printf "[INFO]   [ASV] CD-HIT-EST @100%% deduplication -> %s\n" "$CDHIT_OUT"
      cd-hit-est -i "$NOCHIM_FASTA" -o "$CDHIT_OUT" \
                 -c "$CFG_cdhit_identity" -aS "$CFG_cdhit_asv_coverage" -T "$THREADS" -M 0 -d 0

      # primer trimming on ASV representatives
      printf "[INFO]   [ASV] Primer trimming of %s\n" "$CDHIT_OUT"
      for e in "${PRIMERS[@]}"; do
        IFS=":" read -r name fwd rev <<<"$e"
        cutadapt -g "^${fwd}" -a "${rev}\$" -e "$CFG_primer_trim_error_rate" -O "$CFG_primer_trim_min_overlap" --fasta \
          -o "$MARKER_OUT/${name}_asvs.fasta" \
          "$CDHIT_OUT"
      done

      # symlink to expected final_consensus.fasta path for downstream compat
      cp "$CDHIT_OUT" "$FINAL_CONS_FASTA"

    elif $RECLUSTER; then
      # OTU: CD-HIT-EST @100% → MAFFT → consensus.pl → primer trim
      local CDHIT_OUT="$COMBINED_OUT/combined_cdhit100.fasta"
      if [[ ! -f "$CDHIT_OUT" ]]; then
        printf "[INFO]   CD-HIT-EST @100%% -> %s\n" "$CDHIT_OUT"
        cd-hit-est -i "$ALL_CONS" -o "$CDHIT_OUT" -c "$CFG_cdhit_identity" -T "$THREADS" -M 0 -d 0
      fi

      # per-cluster MSA & consensus
      printf "[INFO]   Per-cluster MSA & consensus -> %s\n" "$FINAL_CONS_FASTA"
      : > "$FINAL_CONS_FASTA"
      local CLSTR="$CDHIT_OUT.clstr"
      local cluster_id=-1 seq_ids=()

      process_cluster() {
        local cid cfasta afasta consf n
        cid=$(printf "%03d" "$1")
        cfasta="$COMBINED_OUT/cluster_${cid}.fasta"
        afasta="$COMBINED_OUT/cluster_${cid}_aligned.fasta"
        consf="$COMBINED_OUT/cluster_${cid}_consensus.fasta"
        : > "$cfasta"
        for id in "${seq_ids[@]}"; do
          awk -v id=">$id" '$0==id{p=1;print;next}p&&/^>/{exit}p{print}' "$ALL_CONS" >> "$cfasta"
        done
        n=$(grep -c '^>' "$cfasta")
        if (( n>1 )); then
          # Subsample to 100 reads before MAFFT to keep the MSA tractable on large clusters.
          # Mirrors PIMENTA's RandSize=100 setting (§2.2.1). Seed 42 ensures reproducibility.
          local msa_input="$cfasta"
          if (( n > 100 )); then
            local sub_fasta="${cfasta%.fasta}_sub100.fasta"
            seqtk sample -s42 "$cfasta" 100 > "$sub_fasta"
            msa_input="$sub_fasta"
          fi
          mafft --auto --thread "$THREADS" --localpair --maxiterate 1000 "$msa_input" > "$afasta"
        else
          cp "$cfasta" "$afasta"
        fi
        perl "$CONSENSUS_PL" -in "$afasta" -out "$consf" -t 30
        [[ -s "$consf" ]] && sed -e "1s/.*/>Cluster${cid}_consensus/" "$consf" >> "$FINAL_CONS_FASTA"
      }

      while IFS= read -r line; do
        if [[ $line =~ ^\>Cluster\ ([0-9]+) ]]; then
          (( ${#seq_ids[@]} )) && process_cluster "$cluster_id"
          cluster_id="${BASH_REMATCH[1]}" seq_ids=()
        elif [[ $line =~ ,\ >(consensus_cl_id_[^[:space:]]+) ]]; then
          seq_ids+=( "${BASH_REMATCH[1]}" )
        fi
      done < "$CLSTR"
      (( ${#seq_ids[@]} )) && process_cluster "$cluster_id"

      # primer trimming on final consensus
      printf "[INFO]   Primer-based trimming of %s\n" "$FINAL_CONS_FASTA"
      for e in "${PRIMERS[@]}"; do
        IFS=":" read -r name fwd rev <<<"$e"
        cutadapt -g "^${fwd}" -a "${rev}\$" -e "$CFG_primer_trim_error_rate" -O "$CFG_primer_trim_min_overlap" --fasta \
          -o "$MARKER_OUT/${name}_consensus.fasta" \
          "$FINAL_CONS_FASTA"
      done

    else
      # OTU: skip recluster — pool & trim per-barcode
      printf "[INFO]   Skipping recluster -- using Phase 2 consensus\n"
      cp "$ALL_CONS" "$FINAL_CONS_FASTA"

      printf "[INFO] Phase 3d: per-barcode pooling & primer trimming\n"
      for sample_dir in "$PER_BARCODE_OUT"/*; do
        local sample="$(basename "$sample_dir")"
        local pbdir="$sample_dir/per-barcode-consensus"
        mkdir -p "$pbdir"
        local pooled="$pbdir/${sample}_pooled.fasta"
        : > "$pooled"
        find "$sample_dir" -path '*/medaka_cl_id_*/consensus.fasta' -print0 | \
        while IFS= read -r -d '' f; do
          cat "$f" >> "$pooled"; printf "\n" >> "$pooled"
        done

        for e in "${PRIMERS[@]}"; do
          IFS=":" read -r name fwd rev <<<"$e"
          cutadapt -g "^${fwd}" -a "${rev}\$" -e "$CFG_primer_trim_error_rate" -O "$CFG_primer_trim_min_overlap" --fasta \
            -o "$pbdir/${sample}_${name}_trimmed.fasta" \
            "$pooled"
        done
      done
    fi
  fi

  printf "\n[OK]   Phase 3 complete.\n"
  printf "       Per-barcode dirs:     %s/*\n" "$PER_BARCODE_OUT"
  printf "       Combined consensus:   %s\n" "$FINAL_CONS_FASTA"
  if [[ "$MODE" == "asv" ]]; then
    printf "       ASV representatives:  %s/asv_representatives.fasta\n" "$COMBINED_OUT"
    printf "       Chimeras (excluded):  %s/chimeras.fasta\n" "$COMBINED_OUT"
  fi
  printf "       Primer outputs:       %s/*\n" "$MARKER_OUT"
  printf "       Taxonomy inputs:      %s/\n" "$TAXON_DIR"
}

# ─────────────────────────────────────────────────────────────────────────────
# Function: count_cluster_size
# Purpose : Sum total_supporting_reads across all clusters for each barcode
# Usage   : count_cluster_size <CONS_BASE> [OUT_FILE]
# Example : count_cluster_size "$CONS_BASE" "$CONS_BASE/cluster_totals.tsv"
# ─────────────────────────────────────────────────────────────────────────────
count_cluster_size() {
  local CONS_BASE="${1%/}"
  local OUT_FILE="${2:-$CONS_BASE/cluster_totals.tsv}"

  # write header
  printf 'Barcode\tTotalSupportReads\n' > "$OUT_FILE"

  # find every trimmed-FASTA under per_barcode/*/per-barcode-consensus
  find "$CONS_BASE"/per_barcode/*/per-barcode-consensus -type f -name '*_trimmed.fasta' | \
  while IFS= read -r fasta; do
    # derive sample/barcode name
    local sample_dir
    sample_dir="$(dirname "$fasta")"                 # .../per-barcode-consensus
    local barcode_dir
    barcode_dir="$(dirname "$sample_dir")"           # .../per_barcode/<barcode>
    local barcode
    barcode="$(basename "$barcode_dir")"             # e.g. eDNA_2024_Barcode01

    # extract and sum the numbers after '_total_supporting_reads_'
    local total
    total="$(grep '^>' "$fasta" \
      | sed 's/.*_total_supporting_reads_//' \
      | awk '{sum+=$1} END{print sum}')"

    printf '%s\t%s\n' "$barcode" "${total:-0}" >> "$OUT_FILE"
  done

  printf "[OK]   Cluster support totals written to %s\n" "$OUT_FILE"
}

# ─────────────────────────────────────────────────────────────────────────────
# Phase 4 — BLASTn taxonomic assignment and Krona/edna-viz reporting.
# Single 17-column megablast run, split by awk into:
#   (1) 15-col blast.tsv  -> ktImportText -> interactive Krona HTML per barcode/combined
#   (2) 10-col blast_results.tsv -> edna-viz Shiny app (Golem format)
# Mammalia-only DB (taxid 40674, excluding Hominidae 9604).
# Thresholds: pident>=97, qcov>=90 (stricter than PIMENTA's >=90%/>=90% -- appropriate
# for mammal 16S amplicons where 97% is the species-level boundary).
# LCA resolution follows PIMENTA §2.2.2: top-bitscore hits within 2% margin;
# downgrade to genus when species disagree, family_or_higher when genera disagree.
# MIN_PIDENT=97, MIN_QCOV=90, BITSCORE_MARGIN=0.98 passed to lca_assignment.R.
#
# Function: blast_taxonomy
# Usage   : blast_taxonomy <cons_base> <db_dir> [threads] [mode]
# ─────────────────────────────────────────────────────────────────────────────
blast_taxonomy() {
  local CONS_BASE="${1%/}"
  local DB_DIR="${2%/}"
  local THREADS="${3:-8}"
  local MODE="${4:-otu}"

  printf "[INFO] Phase 4: Taxonomic Assignment & Reporting\n"
  printf "       Consensus base: %s\n" "$CONS_BASE"
  printf "       BLAST DB dir:   %s\n" "$DB_DIR"
  printf "       Threads:        %s\n" "$THREADS"
  printf "       Mode:           %s\n" "$MODE"
  printf "\n"

  if [[ -d "$DB_DIR" ]]; then
    DB_OPT="-db $DB_DIR/$CFG_blast_db_name"
    THREAD_OPT="-num_threads $THREADS"
    FILTER_OPT="-taxids $CFG_blast_taxid -negative_taxids $CFG_blast_negative_taxid"
  else
    DB_OPT="-remote -db $DB_DIR"
    THREAD_OPT=""
    FILTER_OPT="-entrez_query '$CFG_blast_entrez_query'"
  fi

  printf "[INFO] Using BLAST options: %s %s\n\n" "$DB_OPT" "$THREAD_OPT"

  # 4a: per-barcode
  printf "[INFO] Phase 4a: Assigning taxonomy per barcode\n"
  for fasta in "$CONS_BASE"/per_barcode/*/per-barcode-consensus/*_trimmed.fasta; do
    local sample
    sample=$(basename "$(dirname "$(dirname "$fasta")")")
    local outdir="$CONS_BASE/per_barcode/$sample/taxonomy"
    mkdir -p "$outdir"

    printf "[INFO]   [%s] BLASTn megablast (17-col) -> %s/blast_all.tsv\n" "$sample" "$outdir"
    blastn -task megablast \
           $DB_OPT \
           $FILTER_OPT \
           -query "$fasta" \
           -perc_identity "$CFG_blast_perc_identity" \
           -qcov_hsp_perc "$CFG_blast_qcov_hsp_perc" \
           -max_target_seqs "$CFG_blast_max_target_seqs" \
           -evalue "$CFG_blast_evalue" \
           -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames stitle qcovs' \
           $THREAD_OPT \
           -out "$outdir/blast_all.tsv"

    # Derive Krona format (15-col: cols 1-14 + 16, skipping scomnames at col 15)
    awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$16}' \
        "$outdir/blast_all.tsv" > "$outdir/blast.tsv"

    # Derive Golem app format (10-col: qseqid sseqid pident qcovs evalue bitscore staxids sscinames scomnames stitle)
    awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,$17,$11,$12,$13,$14,$15,$16}' \
        "$outdir/blast_all.tsv" > "$outdir/blast_results.tsv"

    printf "[INFO]         Generating Krona chart -> %s/krona.html\n" "$outdir"
    ktImportText -o "$outdir/krona.html" "$outdir/blast.tsv"
  done

  # 4b: combined
  printf "\n[INFO] Phase 4b: Assigning taxonomy to pooled final consensus\n"
  mkdir -p "$CONS_BASE/combined/taxonomy"
  printf "[INFO]   BLASTn megablast (17-col) -> combined/taxonomy/blast_all.tsv\n"
  blastn -task megablast \
         $DB_OPT \
         $FILTER_OPT \
         -query "$CONS_BASE/combined/final_consensus.fasta" \
         -perc_identity "$CFG_blast_perc_identity" \
         -qcov_hsp_perc "$CFG_blast_qcov_hsp_perc" \
         -max_target_seqs "$CFG_blast_max_target_seqs" \
         -evalue "$CFG_blast_evalue" \
         -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames stitle qcovs' \
         $THREAD_OPT \
         -out "$CONS_BASE/combined/taxonomy/blast_all.tsv"

  # Derive Krona format (15-col: cols 1-14 + 16, skipping scomnames at col 15)
  awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$16}' \
      "$CONS_BASE/combined/taxonomy/blast_all.tsv" > "$CONS_BASE/combined/taxonomy/blast.tsv"

  # Derive Golem app format (10-col: qseqid sseqid pident qcovs evalue bitscore staxids sscinames scomnames stitle)
  awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,$17,$11,$12,$13,$14,$15,$16}' \
      "$CONS_BASE/combined/taxonomy/blast_all.tsv" > "$CONS_BASE/combined/taxonomy/blast_results.tsv"

  printf "[INFO]   LCA assignment -> combined/taxonomy/taxonomy.tsv\n"
  Rscript "$(dirname "${BASH_SOURCE[0]}")/lca_assignment.R" \
    "$CONS_BASE/combined/taxonomy/blast_results.tsv" \
    "$CONS_BASE/combined/taxonomy/taxonomy.tsv" \
    || printf "[WARN]   lca_assignment.R not found or failed -- skipping LCA step\n" >&2

  printf "[INFO]   Krona chart -> combined/taxonomy/krona.html\n"
  ktImportText -o "$CONS_BASE/combined/taxonomy/krona.html" \
               "$CONS_BASE/combined/taxonomy/blast.tsv"

  # 4c: summary guidance
  printf "\n[INFO] Phase 4c: Summaries & downstream\n"
  printf "       View combined taxonomy:    %s/combined/taxonomy/taxonomy.tsv\n" "$CONS_BASE"
  printf "       Load blast_results.tsv into the edna-viz Shiny app for interactive mapping\n"
  printf "       Merge per-barcode blast.tsv files into a master table\n"
  printf "       Filter hits by bitscore/e-value or percent identity\n"
  printf "       Use Krona charts for interactive taxonomy review\n"
  printf "       Convert blast.tsv to BIOM or taxonomy tables for R/Python plotting\n"
  printf "\n[OK]   Phase 4 complete.\n"
}

# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline entrypoint
# ─────────────────────────────────────────────────────────────────────────────
RAW_DIR="${1%/}"
TRIMMED_DIR="${2%/}"
CONS_BASE="${3%/}"
DB_DIR="${4%/}"
THREADS="${5:-8}"
RECLUSTER="${6:-true}"
PRIMER_FILE="${7:-}"
MODE="${8:-otu}"
MODULES="${9:-all}"

# Write run settings for reproducibility (mirrors PIMENTA's $MID.settings.txt).
mkdir -p "$CONS_BASE"
SETTINGS_FILE="$CONS_BASE/run_settings.txt"
{
  printf "eDNA Metabarcoding Pipeline -- Run Settings\n"
  printf "Date:        %s\n" "$(date)"
  printf "RAW_DIR:     %s\n" "$RAW_DIR"
  printf "TRIMMED_DIR: %s\n" "$TRIMMED_DIR"
  printf "CONS_BASE:   %s\n" "$CONS_BASE"
  printf "DB_DIR:      %s\n" "$DB_DIR"
  printf "THREADS:     %s\n" "$THREADS"
  printf "RECLUSTER:   %s\n" "$RECLUSTER"
  printf "MODE:        %s\n" "$MODE"
  printf "MODULES:     %s\n" "$MODULES"
  [[ -n "$PRIMER_FILE" ]] && printf "PRIMER_FILE: %s\n" "$PRIMER_FILE"
  printf "CONFIG:      %s\n" "$PIPELINE_CONFIG"
  printf "BLAST thresholds: pident>=%s, qcov>=%s, evalue<=%s, bitscore_margin=%s\n" \
    "$CFG_blast_perc_identity" "$CFG_blast_qcov_hsp_perc" "$CFG_blast_evalue" "$CFG_blast_bitscore_margin"
} > "$SETTINGS_FILE"
printf "[INFO] Settings written to %s\n" "$SETTINGS_FILE"

printf "[INFO] Starting pipeline:\n"
printf "       RAW_DIR     = %s\n" "$RAW_DIR"
printf "       TRIMMED_DIR = %s\n" "$TRIMMED_DIR"
printf "       CONS_BASE   = %s\n" "$CONS_BASE"
printf "       DB_DIR      = %s\n" "$DB_DIR"
printf "       THREADS     = %s\n" "$THREADS"
printf "       RECLUSTER   = %s\n" "$RECLUSTER"
printf "       MODE        = %s\n" "$MODE"
printf "       MODULES     = %s\n" "$MODULES"
[[ -n "$PRIMER_FILE" ]] && printf "       PRIMER_FILE = %s\n" "$PRIMER_FILE"
printf "\n"

# Phase 1
if [[ "$MODULES" == *"all"* ]] || [[ "$MODULES" == *"Phase1"* ]]; then
  if find "$TRIMMED_DIR" -type f -name "*_trimmed.fastq.gz" | grep -q .; then
    printf "[SKIP] Phase 1: trimmed reads already present in %s\n" "$TRIMMED_DIR"
  else
    preprocess_reads "$RAW_DIR" "$TRIMMED_DIR"
  fi
fi

# Phase 2
if [[ "$MODULES" == *"all"* ]] || [[ "$MODULES" == *"Phase2"* ]]; then
  if find "$TRIMMED_DIR/merged_fastq" -name '*_merged.fastq' 2>/dev/null | grep -q .; then
    printf "\n[SKIP] Phase 2: merged FASTQs already present in %s/merged_fastq\n" "$TRIMMED_DIR"
  else
    printf "\n[INFO] Phase 2: QC & merge trimmed FASTQs\n"
    quality_check_merge "$TRIMMED_DIR" "$THREADS"
  fi
fi

# Phase 3
if [[ "$MODULES" == *"all"* ]] || [[ "$MODULES" == *"Phase3"* ]]; then
  printf "\n[INFO] Phase 3: clustering & consensus (mode=%s, recluster=%s)\n" "$MODE" "$RECLUSTER"
  run_consensus "$TRIMMED_DIR/merged_fastq" "$CONS_BASE" "$THREADS" "$RECLUSTER" "$PRIMER_FILE" "$MODE"
fi

# Phase 4
if [[ "$MODULES" == *"all"* ]] || [[ "$MODULES" == *"Phase4"* ]]; then
  # Pre-flight guard: final_consensus.fasta must exist before BLAST can run.
  # If Phase 3 was skipped but the file is absent, the BLAST query will silently fail.
  if [[ ! -f "$CONS_BASE/combined/final_consensus.fasta" ]]; then
    printf "[ERROR] final_consensus.fasta missing -- Phase 3 may have failed or was skipped.\n" >&2
    printf "        Expected: %s/combined/final_consensus.fasta\n" "$CONS_BASE" >&2
    printf "        Re-run with MODULES=Phase3 or MODULES=all to generate it.\n" >&2
    exit 1
  fi
  printf "\n"
  blast_taxonomy "$CONS_BASE" "$DB_DIR" "$THREADS" "$MODE"
fi

printf "\n[OK]   Pipeline complete.\n"
exit 0
