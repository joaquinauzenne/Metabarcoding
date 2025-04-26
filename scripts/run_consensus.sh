#!/usr/bin/env bash
set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Script : run_consensus.sh
# Purpose: Phase 2 & Phase 3 of the eDNA pipeline, with optional reclustering:
#   • Per-sample clustering & consensus (NGSpeciesID @90%)
#   • Cross-sample collapse @100% + MSA & consensus
#   • Primer-based trimming of final consensus (Phase 3b)
#   • OR, if no recluster: per-barcode pooling + primer trimming (Phase 3d)
#
# Usage:
#   bash run_consensus.sh \
#        /path/to/merged_fastq_dir \
#        /path/to/output_base_dir \
#        [threads] \
#        [--no-recluster]
# ─────────────────────────────────────────────────────────────────────────────

MERGED_DIR="${1%/}"
OUT_BASE="${2%/}"
THREADS="${3:-8}"
RECLUSTER=true
if [[ "${4:-}" == "--no-recluster" ]]; then
  RECLUSTER=false
fi

echo "[DEBUG] MERGED_DIR = $MERGED_DIR"
echo "[DEBUG] OUT_BASE   = $OUT_BASE"
echo "[DEBUG] THREADS    = $THREADS"
echo "[DEBUG] RECLUSTER  = $RECLUSTER"

# Verify input
if [[ ! -d "$MERGED_DIR" ]]; then
  echo "[ERROR] Directory not found: $MERGED_DIR" >&2
  exit 1
fi

# Paths
PER_BARCODE_OUT="$OUT_BASE/per_barcode"
COMBINED_OUT="$OUT_BASE/combined"
ALL_CONS="$COMBINED_OUT/all_consensus.fasta"
FINAL_CONS_FASTA="$COMBINED_OUT/final_consensus.fasta"
MARKER_OUT="$OUT_BASE/by_marker"
TAXON_DIR="$OUT_BASE/../taxonomy"
CONSENSUS_PL="$(dirname "$0")/consensus.pl"

mkdir -p \
  "$PER_BARCODE_OUT" \
  "$COMBINED_OUT" \
  "$MARKER_OUT" \
  "$TAXON_DIR"

# ─────────────────────────────────────────────────────────────────────────────
# Phase 2: per-sample clustering & Medaka polishing
# ─────────────────────────────────────────────────────────────────────────────
echo
echo "▶ Phase 2: per-sample clustering & consensus (NGSpeciesID @90%)"
shopt -s nullglob
mapfile -t FQS < <(printf '%s\n' "$MERGED_DIR"/*_merged.fastq)

if (( ${#FQS[@]} == 0 )); then
  echo "[ERROR] No merged FASTQ files found in $MERGED_DIR" >&2
  exit 1
fi
echo "[DEBUG] Found ${#FQS[@]} FASTQ files to process"

set -x
for fq in "${FQS[@]}"; do
  sample="$(basename "$fq" _merged.fastq)"
  outdir="$PER_BARCODE_OUT/$sample"
  mkdir -p "$outdir"

  # Skip if already polished
  if compgen -G "${outdir}/medaka_cl_id_*/consensus.fasta" > /dev/null; then
    echo "  ✔ SKIP [per-barcode] $sample"
    continue
  fi

  NGSpeciesID \
    --fastq "$fq" \
    --q 10.0 --ont \
    --consensus --medaka \
    --m 550 --s 550 \
    --sample_size 100 \
    --outfolder "$outdir"
done
set +x


# ─────────────────────────────────────────────────────────────────────────────
# Phase 3: concatenate + (optional) recluster + consensus
# ─────────────────────────────────────────────────────────────────────────────
echo
echo "▶ Phase 3: collapse @100% + MSA & final consensus"

if [[ -s "$FINAL_CONS_FASTA" ]]; then
  echo "  ✔ SKIP Phase 3 (already done)"
else
  # 3a) Concatenate all per-barcode medaka consensuses
  if [[ ! -s "$ALL_CONS" ]]; then
    echo "  ⟳ 3a) Building $ALL_CONS"
    : > "$ALL_CONS"
    while IFS= read -r -d '' f; do
      cat "$f" >> "$ALL_CONS"
      echo ""      >> "$ALL_CONS"
    done < <(
      find "$PER_BARCODE_OUT" -mindepth 2 -path '*/medaka_cl_id_*/consensus.fasta' -print0
    )
  fi

  if $RECLUSTER; then
    # 3b) recluster @100%
    CDHIT_OUT="$COMBINED_OUT/combined_cdhit100.fasta"
    CLSTR="$CDHIT_OUT.clstr"
    if [[ ! -f "$CDHIT_OUT" ]]; then
      echo "  ⟳ 3b) CD-HIT-EST @100% → $CDHIT_OUT"
      cd-hit-est -i "$ALL_CONS" -o "$CDHIT_OUT" -c 1.00 -T "$THREADS" -M 0 -d 0
    fi

    # 3c) per-cluster MSA & consensus
    echo "  ⟳ 3c) Building per-cluster consensus → $FINAL_CONS_FASTA"
    : > "$FINAL_CONS_FASTA"
    cluster_id=-1; seq_ids=()
    process_cluster() {
      local cid=$1 cfasta afasta consf n
      cid=$(printf "%03d" "$cid")
      cfasta="$COMBINED_OUT/cluster_${cid}.fasta"
      afasta="$COMBINED_OUT/cluster_${cid}_aligned.fasta"
      consf="$COMBINED_OUT/cluster_${cid}_consensus.fasta"
      : > "$cfasta"
      for id in "${seq_ids[@]}"; do
        awk -v id=">$id" '
          $0==id {p=1;print;next}
          p && /^>/ {exit}
          p{print}
        ' "$ALL_CONS" >> "$cfasta"
      done
      n=$(grep -c '^>' "$cfasta")
      if (( n > 1 )); then
        mafft --auto --thread "$THREADS" "$cfasta" > "$afasta"
      else
        cp "$cfasta" "$afasta"
      fi
      perl "$CONSENSUS_PL" -in "$afasta" -out "$consf" -t 30
      [[ -s "$consf" ]] && sed -e "1s/.*/>Cluster${cid}_consensus/" "$consf" >> "$FINAL_CONS_FASTA"
    }
    while IFS= read -r line; do
      if [[ $line =~ ^\>Cluster\ ([0-9]+) ]]; then
        (( ${#seq_ids[@]} )) && process_cluster "$cluster_id"
        cluster_id="${BASH_REMATCH[1]}"; seq_ids=()
      elif [[ $line =~ ,\ >(consensus_cl_id_[^[:space:]]+) ]]; then
        seq_ids+=( "${BASH_REMATCH[1]}" )
      fi
    done < "$CLSTR"
    (( ${#seq_ids[@]} )) && process_cluster "$cluster_id"

    # 3b) Primer removal from final consensus
    echo "  ⟳ 3b) Primer-based trimming of $FINAL_CONS_FASTA"
    declare -a PRIMERS=("16Smam:CGGGTTGGGGTGACCTCGGA:GCTGTTATCCCTAGGGTAACT")
    for e in "${PRIMERS[@]}"; do
      IFS=":" read -r name fwd rev <<<"$e"
      cutadapt -g "^${fwd}" -a "${rev}\$" -e 0.15 -O 25 --fasta \
        -o "$MARKER_OUT/${name}_consensus.fasta" \
        "$FINAL_CONS_FASTA"
    done

  else
    # No recluster: copy ALL_CONS → FINAL_CONS and do per-barcode trimming
    echo "  ⟳ Skipping recluster — using Phase 2 consensus"
    cp "$ALL_CONS" "$FINAL_CONS_FASTA"

    echo "▶ Phase 3d: per-barcode pooling & primer trimming"
    for sample_dir in "$PER_BARCODE_OUT"/*; do
      sample=$(basename "$sample_dir")
      pbdir="$sample_dir/per-barcode-consensus"
      mkdir -p "$pbdir"
      pooled="$pbdir/${sample}_pooled.fasta"

      # Pool medaka consensuses for this barcode
      : > "$pooled"
      while IFS= read -r -d '' f; do
        cat "$f" >> "$pooled"
        echo ""    >> "$pooled"
      done < <(
        find "$sample_dir" -path '*/medaka_cl_id_*/consensus.fasta' -print0
      )

      # Trim primers on each pooled barcode
      declare -a PRIMERS=("16Smam:CGGGTTGGGGTGACCTCGGA:GCTGTTATCCCTAGGGTAACT")
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
echo "✔ All done."
echo "  • Per-barcode dirs:     $PER_BARCODE_OUT/*"
echo "  • Combined consensus:   $FINAL_CONS_FASTA"
echo "  • Primer outputs:       $MARKER_OUT/*"
echo "  • Taxonomy inputs:      $TAXON_DIR/"
