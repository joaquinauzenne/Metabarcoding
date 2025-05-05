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
#        [THREADS] [RECLUSTER] [PRIMER_FILE]
#
#   RAW_DIR      Directory with raw .fastq(.gz) files
#   TRIMMED_DIR  Output dir for trimmed FASTQs
#   CONS_BASE    Base output dir for consensus step
#   DB_DIR       BLAST DB directory (local nt*)
#   THREADS      (optional) Number of threads [default: (1/2)npoc for local, all for HPC]
#   RECLUSTER    (optional) true|false — whether to recluster combined consensus [default: true]
#   PRIMER_FILE  (optional) Path to a file of primers (name:fwd:rev per line);
#                if omitted, uses built-in 16Smam primer
# ─────────────────────────────────────────────────────────────────────────────

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")")/scripts"

# ─────────────────────────────────────────────────────────────────────────────
# Function: preprocess_reads
# Purpose : Phase 1 of the eDNA pipeline — adapter trimming, poly-A/T removal,
#           quality- and length-filtering using Cutadapt.
# Usage   : preprocess_reads <input_fastq_dir> <output_trimmed_dir>
# ─────────────────────────────────────────────────────────────────────────────
preprocess_reads() {
  local INPUT_DIR="${1%/}"
  local OUTPUT_DIR="${2%/}"
  local MIN_QUAL=12 MIN_LEN=10 MAX_LEN=1100

  mkdir -p "$OUTPUT_DIR"
  echo "▶ Phase 1: preprocessing reads → $OUTPUT_DIR"
  echo "   Quality cutoff: $MIN_QUAL"
  echo "   Length range  : $MIN_LEN–$MAX_LEN bp"

  find "$INPUT_DIR" -type f \( -name '*.fastq' -o -name '*.fastq.gz' \) | \
  while IFS= read -r IN; do
    local REL="${IN#"$INPUT_DIR"/}"
    local base="${REL%.fastq.gz}"; base="${base%.fastq}"
    local OUT="$OUTPUT_DIR/${base}_trimmed.fastq.gz"
    mkdir -p "$(dirname "$OUT")"

    if [[ -s "$OUT" ]]; then
      echo "✔ SKIP: $REL (already exists)"
      continue
    fi

    echo "⟳ Trimming & filtering: $REL"
    python -m cutadapt \
      --poly-a --trim-n --max-n 0.05 \
      --minimum-length "$MIN_LEN" --maximum-length "$MAX_LEN" \
      --quality-cutoff "$MIN_QUAL" --cores 0 \
      -o "$OUT" "$IN"

    if [[ -s "$OUT" ]]; then
      echo "✔ Trimmed: ${base}_trimmed.fastq.gz"
    else
      echo "✖ ERROR: no output for $REL" >&2
      rm -f "$OUT"
    fi
  done

  echo "✔ Phase 1 complete."
}

# ─────────────────────────────────────────────────────────────────────────────
# Function: quality_check_merge
# Purpose : Phase 2 — FastQC, merge by barcode, NanoPlot, MultiQC
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
    echo "✅ FastQC outputs already present; skipping."
  else
    echo "🔍 Running FastQC on all *_trimmed.fastq(.gz)..."
    find "$TRIMMED_DIR" -type f \( -name "*_trimmed.fastq" -o -name "*_trimmed.fastq.gz" \) | \
    while IFS= read -r fq; do
      echo "▶ FastQC on $fq"
      fastqc --threads "$THREADS" --outdir "$FASTQC_OUT" "$fq" \
        && echo "✔ FastQC succeeded for $fq" \
        || echo "✖ FastQC FAILED for $fq" >&2
    done
  fi

  # Merge by barcode
  echo "🔗 Merging trimmed FASTQs by barcode into $MERGED_DIR..."
  find "$TRIMMED_DIR" -type f \( -name "*_trimmed.fastq" -o -name "*_trimmed.fastq.gz" \) \
    -printf '%h\n' | sort -u | while IFS= read -r BDIR; do
    local SAMPLE_TAG="${BDIR#${TRIMMED_DIR}/}"
    SAMPLE_TAG="${SAMPLE_TAG//\//_}"
    local MERGED_FASTQ="$MERGED_DIR/${SAMPLE_TAG}_merged.fastq"

    echo "  • Merging $BDIR → $MERGED_FASTQ"
    {
      find "$BDIR" -type f \( -name "*_trimmed.fastq" -o -name "*_trimmed.fastq.gz" \) \
        | sort | \
        while IFS= read -r in; do
          [[ "$in" == *.gz ]] && zcat "$in" || cat "$in"
        done
    } > "$MERGED_FASTQ"
  done

  # NanoPlot
  echo "📊 Running NanoPlot on merged FASTQs..."
  find "$MERGED_DIR" -type f -name "*_merged.fastq" | \
  while IFS= read -r mfq; do
    local SAMPLE="$(basename "$mfq" _merged.fastq)"
    echo "▶ NanoPlot on $mfq"
    NanoPlot --fastq "$mfq" \
             --outdir "$NANOPLOT_OUT/$SAMPLE" \
             --threads "$THREADS" \
             --N50 \
             --plots hex dot kde \
      && echo "✔ NanoPlot complete for $mfq" \
      || echo "✖ NanoPlot FAILED for $mfq" >&2
  done

  # MultiQC
  echo "📈 Aggregating FastQC with MultiQC..."
  multiqc "$FASTQC_OUT" --outdir "$MULTIQC_OUT" \
    && echo "✔ MultiQC report: $MULTIQC_OUT/multiqc_report.html" \
    || echo "✖ MultiQC failed" >&2

  echo "✅ QC & merge complete."
  echo "   FastQC:    $FASTQC_OUT"
  echo "   Merged:    $MERGED_DIR"
  echo "   NanoPlot:  $NANOPLOT_OUT"
  echo "   MultiQC:   $MULTIQC_OUT/multiqc_report.html"
}

# ─────────────────────────────────────────────────────────────────────────────
# Function: run_consensus
# Purpose : Phase 3 — NGSpeciesID clustering/consensus + optional reclustering
# Usage   : run_consensus <merged_fastq_dir> <output_base_dir> <threads> <true|false>
# ─────────────────────────────────────────────────────────────────────────────
run_consensus() {
  local MERGED_DIR="${1%/}"
  local OUT_BASE="${2%/}"
  local THREADS="${3:-8}"
  local RECLUSTER_FLAG="${4:-true}"
  local PRIMER_FILE="${5:-}"

  [[ "$RECLUSTER_FLAG" != "true" ]] && RECLUSTER=false || RECLUSTER=true

  # load primers
  declare -a PRIMERS
  if [[ -n "$PRIMER_FILE" && -f "$PRIMER_FILE" ]]; then
    mapfile -t PRIMERS < "$PRIMER_FILE"
    echo "▶ Loaded primers from $PRIMER_FILE"
  else
    PRIMERS=( "16Smam:CGGGTTGGGGTGACCTCGGA:GCTGTTATCCCTAGGGTAACT" )
    echo "▶ Using default primer: ${PRIMERS[0]}"
  fi

  local PER_BARCODE_OUT="$OUT_BASE/per_barcode"
  local COMBINED_OUT="$OUT_BASE/combined"
  local ALL_CONS="$COMBINED_OUT/all_consensus.fasta"
  local FINAL_CONS_FASTA="$COMBINED_OUT/final_consensus.fasta"
  local MARKER_OUT="$OUT_BASE/by_marker"
  local TAXON_DIR="$OUT_BASE/../taxonomy"
  local CONSENSUS_PL="$SCRIPTS_DIR/consensus.pl"

  mkdir -p "$PER_BARCODE_OUT" "$COMBINED_OUT" "$MARKER_OUT" "$TAXON_DIR"

  # 3a) Per-barcode NGSpeciesID
  echo
  echo "▶ Phase 3: per-sample clustering & consensus (NGSpeciesID @90%)"
  shopt -s nullglob
  mapfile -t FQS < <(printf '%s\n' "$MERGED_DIR"/*_merged.fastq)
  (( ${#FQS[@]} )) || { echo "[ERROR] No merged FASTQ in $MERGED_DIR" >&2; exit 1; }
  set -x
  for fq in "${FQS[@]}"; do
    local sample="$(basename "$fq" _merged.fastq)"
    local outdir="$PER_BARCODE_OUT/$sample"
    mkdir -p "$outdir"
    if compgen -G "${outdir}/medaka_cl_id_*/consensus.fasta" > /dev/null; then
      echo "  ✔ SKIP [per-barcode] $sample"
      continue
    fi
    NGSpeciesID --fastq "$fq" \
                --q 10.0 --ont --consensus --medaka \
                --m 550 --s 550 --sample_size 100 \
                --outfolder "$outdir"
  done
  set +x

  # 3b) Collapse + optional recluster + MAFFT + consensus.pl
  echo
  echo "▶ Phase 3: collapse @100% + MSA & final consensus"
  if [[ -s "$FINAL_CONS_FASTA" ]]; then
    echo "  ✔ SKIP Phase 3 (already done)"
  else
    # concatenate
    if [[ ! -s "$ALL_CONS" ]]; then
      echo "  ⟳ Building $ALL_CONS"
      : > "$ALL_CONS"
      while IFS= read -r -d '' f; do
        cat "$f" >> "$ALL_CONS"; echo "" >> "$ALL_CONS"
      done < <(find "$PER_BARCODE_OUT" -mindepth 2 -path '*/medaka_cl_id_*/consensus.fasta' -print0)
    fi

    if $RECLUSTER; then
      # CD-HIT-EST @100%
      local CDHIT_OUT="$COMBINED_OUT/combined_cdhit100.fasta"
      if [[ ! -f "$CDHIT_OUT" ]]; then
        echo "  ⟳ CD-HIT-EST @100% → $CDHIT_OUT"
        cd-hit-est -i "$ALL_CONS" -o "$CDHIT_OUT" -c 1.00 -T "$THREADS" -M 0 -d 0
      fi

      # per-cluster MSA & consensus
      echo "  ⟳ Per-cluster MSA & consensus → $FINAL_CONS_FASTA"
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
          mafft --auto --thread "$THREADS" --localpair --maxiterate 1000 "$cfasta" > "$afasta"
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
      echo "  ⟳ Primer-based trimming of $FINAL_CONS_FASTA"
      for e in "${PRIMERS[@]}"; do
        IFS=":" read -r name fwd rev <<<"$e"
        cutadapt -g "^${fwd}" -a "${rev}\$" -e 0.15 -O 25 --fasta \
          -o "$MARKER_OUT/${name}_consensus.fasta" \
          "$FINAL_CONS_FASTA"
      done

    else
      # skip recluster: pool & trim per-barcode
      echo "  ⟳ Skipping recluster — using Phase 2 consensus"
      cp "$ALL_CONS" "$FINAL_CONS_FASTA"

      echo "▶ Phase 3d: per-barcode pooling & primer trimming"
      for sample_dir in "$PER_BARCODE_OUT"/*; do
        local sample="$(basename "$sample_dir")"
        local pbdir="$sample_dir/per-barcode-consensus"
        mkdir -p "$pbdir"
        local pooled="$pbdir/${sample}_pooled.fasta"
        : > "$pooled"
        find "$sample_dir" -path '*/medaka_cl_id_*/consensus.fasta' -print0 | \
        while IFS= read -r -d '' f; do
          cat "$f" >> "$pooled"; echo "" >> "$pooled"
        done

        for e in "${PRIMERS[@]}"; do
          IFS=":" read -r name fwd rev <<<"$e"
          cutadapt -g "^${fwd}" -a "${rev}\$" -e 0.15 -O 25 --fasta \
            -o "$pbdir/${sample}_${name}_trimmed.fasta" \
            "$pooled"
        done
      done
    fi
  fi

  echo
  echo "✔ Phase 3 complete."
  echo "  • Per-barcode dirs:     $PER_BARCODE_OUT/*"
  echo "  • Combined consensus:   $FINAL_CONS_FASTA"
  echo "  • Primer outputs:       $MARKER_OUT/*"
  echo "  • Taxonomy inputs:      $TAXON_DIR/"
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

  echo "✔ Cluster support totals written to $OUT_FILE"
}


# ─────────────────────────────────────────────────────────────────────────────
# Function: blast_taxonomy
# (Phase 4) — BLASTn taxonomic assignment & reporting (per‐barcode + combined)
# Usage    : blast_taxonomy <cons_base> <db_dir> [threads]
# ─────────────────────────────────────────────────────────────────────────────
blast_taxonomy() {
  local CONS_BASE="${1%/}"
  local DB_DIR="${2%/}"
  local THREADS="${3:-8}"

  echo "▶ Phase 4: Taxonomic Assignment & Reporting"
  echo "  • Consensus base: $CONS_BASE"
  echo "  • BLAST DB dir:   $DB_DIR"
  echo "  • Threads:        $THREADS"
  echo

  if [[ -d "$DB_DIR" ]]; then
    DB_OPT="-db $DB_DIR/nt"
    THREAD_OPT="-num_threads $THREADS"
    FILTER_OPT="-taxids 40674 -negative_taxids 9604"
  else
    DB_OPT="-remote -db $DB_DIR"
    THREAD_OPT=""
    FILTER_OPT="-entrez_query 'Mammalia[Organism] NOT Hominidae[Organism]'"
  fi

  echo "↳ Using BLAST options: $DB_OPT $THREAD_OPT"
  echo

  # 4a: per-barcode
  echo "▶ Phase 4a: Assigning taxonomy per barcode"
  for fasta in "$CONS_BASE"/per_barcode/*/per-barcode-consensus/*_trimmed.fasta; do
    local sample
    sample=$(basename "$(dirname "$(dirname "$fasta")")")
    local outdir="$CONS_BASE/per_barcode/$sample/taxonomy"
    mkdir -p "$outdir"

    echo "  ⟳ [$sample] BLASTn megablast → $outdir/blast.tsv"
    blastn -task megablast \
           $DB_OPT \
           -query "$fasta" \
           -perc_identity 97 \
           -qcov_hsp_perc 90 \
           -max_target_seqs 100 \
           -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle' \
           $THREAD_OPT \
           -out "$outdir/blast.tsv"

    echo "      Generating Krona chart → $outdir/krona.html"
    ktImportText -o "$outdir/krona.html" "$outdir/blast.tsv"
  done

  # 4b: combined
  echo
  echo "▶ Phase 4b: Assigning taxonomy to pooled final consensus"
  mkdir -p "$CONS_BASE/combined/taxonomy"
  echo "  ⟳ BLASTn megablast → combined/taxonomy/blast.tsv"
  blastn -task megablast \
         $DB_OPT \
         $FILTER_OPT \
         -query "$CONS_BASE/combined/final_consensus.fasta" \
         -perc_identity 97 \
         -qcov_hsp_perc 90 \
         -max_target_seqs 100 \
         -outfmt '6 qseqid sseqid pident length evalue bitscore staxids sscinames stitle' \
         $THREAD_OPT \
         -out "$CONS_BASE/combined/taxonomy/blast.tsv"

  echo "  ⟳ Krona chart → combined/taxonomy/krona.html"
  ktImportText -o "$CONS_BASE/combined/taxonomy/krona.html" \
               "$CONS_BASE/combined/taxonomy/blast.tsv"

  # 4c: summary guidance
  echo
  echo "▶ Phase 4c: Summaries & downstream"
  echo "You can now:"
  echo "  • Merge per-barcode blast.tsv files into a master table"
  echo "  • Filter hits by bitscore/e-value or percent identity"
  echo "  • Use Krona charts for interactive taxonomy review"
  echo "  • Convert blast.tsv to BIOM or taxonomy tables for R/Python plotting"
  echo
  echo "✅ Phase 4 complete."
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

echo "▶ Starting pipeline:"
echo "    RAW_DIR     = $RAW_DIR"
echo "    TRIMMED_DIR = $TRIMMED_DIR"
echo "    CONS_BASE   = $CONS_BASE"
echo "    DB_DIR      = $DB_DIR"
echo "    THREADS     = $THREADS"
echo "    RECLUSTER   = $RECLUSTER"
[[ -n "$PRIMER_FILE" ]] && echo "    PRIMER_FILE = $PRIMER_FILE"
echo

# Phase 1
if find "$TRIMMED_DIR" -type f -name "*_trimmed.fastq.gz" | grep -q .; then
  echo "✔ SKIP Phase 1: trimmed reads already exist."
else
  preprocess_reads "$RAW_DIR" "$TRIMMED_DIR"
fi

# Phase 2
echo
echo "▶ Phase 2: QC & merge trimmed FASTQs"
quality_check_merge "$TRIMMED_DIR" "$THREADS"

# Phase 3
echo
echo "▶ Phase 3: clustering & consensus (recluster=${RECLUSTER})"
run_consensus "$TRIMMED_DIR/merged_fastq" "$CONS_BASE" "$THREADS" "$RECLUSTER" "$PRIMER_FILE"

# Phase 4
echo
blast_taxonomy "$CONS_BASE" "$DB_DIR" "$THREADS"

echo
echo "✅ Pipeline complete."
exit 0
