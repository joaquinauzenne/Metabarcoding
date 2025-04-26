#!/usr/bin/env bash
set -euo pipefail

# Root dir containing exactly these four subfolders
RAW_ROOT="/mnt/c/Users/joaau/Metagenomics/eDNA"

# Primer metadata: NAME:FWD:REV
declare -a PRIMER_PAIRS=(
  "16Smam:CGGGTTGGGGTGACCTCGGA:GCTGTTATCCCTAGGGTAACT"
)

# Initialize per-primer tallies
declare -A total_fwd total_rev
total_reads_global=0

# Header for per‐file report (optional; comment out if not needed)
printf "File\tPrimer\tReads_with_fwd\tReads_with_fwd(%%)\tReads_with_rev\tReads_with_rev(%%)\n"

# Loop over only the four allowed subdirs
find "$RAW_ROOT"/{2024,Qia25,TA12,TA25}/*/ -type f -name "*.fastq.gz" -print0 \
| while IFS= read -r -d '' fq; do
  # Count reads in this file (4 lines per read)
  reads_in_file=$(( $(zcat "$fq" | wc -l) / 4 ))
  total_reads_global=$(( total_reads_global + reads_in_file ))
  sample=$(basename "$fq")

  for meta in "${PRIMER_PAIRS[@]}"; do
    IFS=":" read -r name fwd rev <<<"$meta"

    # Count reads containing each primer (forward / reverse)
    # We look for the primer in the sequence lines only:
    reads_with_fwd=$(zgrep -A1 -m1 -c "$fwd" "$fq" || true)
    reads_with_rev=$(zgrep -A1 -m1 -c "$rev" "$fq" || true)

    # Accumulate into global tallies
    total_fwd[$name]=$(( ${total_fwd[$name]:-0} + reads_with_fwd ))
    total_rev[$name]=$(( ${total_rev[$name]:-0} + reads_with_rev ))

    # Per‐file output (optional; keep or remove)
    pct_fwd=$(awk -v a=$reads_with_fwd -v b=$reads_in_file 'BEGIN{printf "%.2f", b?100*a/b:0}')
    pct_rev=$(awk -v a=$reads_with_rev -v b=$reads_in_file 'BEGIN{printf "%.2f", b?100*a/b:0}')
    printf "%s\t%s\t%d\t%s%%\t%d\t%s%%\n" \
      "$sample" "$name" \
      "$reads_with_fwd" "$pct_fwd" \
      "$reads_with_rev" "$pct_rev"
  done
done

# After processing all files, print overall summary
echo
printf "Primer\tTotal_Reads_with_fwd\tTotal_Reads_with_rev\tPct_fwd_overall\tPct_rev_overall\n"
for meta in "${PRIMER_PAIRS[@]}"; do
  IFS=":" read -r name fwd rev <<<"$meta"
  tf=${total_fwd[$name]:-0}
  tr=${total_rev[$name]:-0}
  pctf=$(awk -v a=$tf -v b=$total_reads_global 'BEGIN{printf "%.2f", b?100*a/b:0}')
  pctr=$(awk -v a=$tr -v b=$total_reads_global 'BEGIN{printf "%.2f", b?100*a/b:0}')
  printf "%s\t%d\t%d\t%s%%\t%s%%\n" \
    "$name" "$tf" "$tr" \
    "$pctf" "$pctr"
done
