#!/usr/bin/env bash
set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Script : blast_taxonomy.sh
# Purpose: Phase 4 of the eDNA pipeline: assign taxonomy to per-barcode and
#          pooled consensus sequences, then generate summary tables + Krona
#
# Usage:
#   bash phase4_taxonomy.sh <CONS_BASE> <DB_DIR> [threads]
#
#   CONS_BASE:  base dir of Phase 3 output (contains per_barcode/ and combined/)
#   DB_DIR:     directory containing formatted BLAST nt DB (nt.*) OR a
#               remote DB name (e.g. nt, nr, refseq_rna)
#   threads:    number of threads to use (default: 4)
# ─────────────────────────────────────────────────────────────────────────────

CONS_BASE="${1%/}"
DB_DIR="${2%/}"
THREADS="${3:-4}"

echo "▶ Phase 4: Taxonomic Assignment & Reporting"
echo "  • Consensus base: $CONS_BASE"
echo "  • BLAST DB dir:   $DB_DIR"
echo "  • Threads:        $THREADS"
echo

# 0) Select local vs. remote BLAST options
if [[ -d "$DB_DIR" ]]; then
  DB_OPT="-db $DB_DIR/nt"
  THREAD_OPT="-num_threads $THREADS"
else
  DB_OPT="-remote -db $DB_DIR"
  THREAD_OPT=""  # remote mode cannot use threads :contentReference[oaicite:1]{index=1}
fi
echo "↳ Using BLAST options: $DB_OPT $THREAD_OPT"
echo

# 1) Format local DB if needed
if [[ "$DB_OPT" =~ ^-db ]] && [[ ! -f "$DB_DIR/nt.00.nsq" ]]; then
  echo "↳ Formatting nt DB in $DB_DIR…"
  makeblastdb -in "$DB_DIR/nt.fa" \
              -dbtype nucl \
              -parse_seqids \
              -out "$DB_DIR/nt"
fi

# 2) Per-barcode taxonomic assignment
echo
echo "▶ Phase 4a: Assigning taxonomy per barcode"
for fasta in "$CONS_BASE"/per_barcode/*/per-barcode-consensus/*_trimmed.fasta; do
  sample=$(basename "$(dirname "$fasta")")
  outdir="$CONS_BASE/per_barcode/$sample/taxonomy"
  mkdir -p "$outdir"

  echo "  ⟳ [$sample] BLASTn megablast → $outdir/blast.tsv"
  blastn -task megablast \
         $DB_OPT \
         -query "$fasta" \
         -perc_identity 97 \
         -qcov_hsp_perc 80 \
         -max_target_seqs 10 \
         -outfmt '6 qseqid sseqid pident length evalue bitscore staxids sscinames stitle' \
         $THREAD_OPT \
         -out "$outdir/blast.tsv"

  echo "      Generating Krona chart → $outdir/krona.html"
  ktImportText -o "$outdir/krona.html" "$outdir/blast.tsv"
done

# 3) Pooled consensus taxonomic assignment
echo
echo "▶ Phase 4b: Assigning taxonomy to pooled final consensus"
pooled_fasta="$CONS_BASE/combined/final_consensus.fasta"
mkdir -p "$CONS_BASE/combined/taxonomy"

echo "  ⟳ BLASTn megablast → combined/taxonomy/blast.tsv"
blastn -task megablast \
       $DB_OPT \
       -query "$pooled_fasta" \
       -perc_identity 97 \
       -qcov_hsp_perc 80 \
       -max_target_seqs 10 \
       -outfmt '6 qseqid sseqid pident length evalue bitscore staxids sscinames stitle' \
       $THREAD_OPT \
       -out "$CONS_BASE/combined/taxonomy/blast.tsv"

echo "  ⟳ Krona chart → combined/taxonomy/krona.html"
ktImportText -o "$CONS_BASE/combined/taxonomy/krona.html" \
             "$CONS_BASE/combined/taxonomy/blast.tsv"

# 4) Summaries
echo
echo "▶ Phase 4c: Summaries & downstream"
echo "You can now:"
echo "  • Merge per-barcode blast.tsv files into a master table"
echo "  • Filter hits by bitscore/e-value or percent identity"
echo "  • Use Krona charts for interactive taxonomy review"
echo "  • Convert blast.tsv to BIOM or taxonomy tables for R/Python plotting"
echo
echo "✅ Phase 4 complete. Next steps: taxon tables, visualization in Python/R, Krona review."
