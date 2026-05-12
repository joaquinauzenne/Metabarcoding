#!/usr/bin/env bash
# set -eo pipefail before module load; -u deferred (TACC bash-completion references
# BASH_COMPLETION_DEBUG without a default, which aborts under set -u).
set -eo pipefail
IFS=$'\n\t'

# ─────────────────────────────────────────────────────────────────────────────
# Script : extract_mammalia_from_partial_nt_euk.sh
# Purpose: Phase A one-shot recovery — build mammalia_nt BLAST DB from the
#          partial nt_euk volumes stranded on disk after SLURM job 3138888
#          hit the /work quota at volume 186/254, then delete nt_euk to reclaim
#          ~500+ GB.
#
# Run inside an idev compute session (NOT on a login node):
#   idev -t 02:00:00 -p development
#   module load tacc-apptainer
#   bash scripts/extract_mammalia_from_partial_nt_euk.sh
#
# Coverage note: only the volumes already extracted on disk are used. Coverage
# is partial; some Mammalia sequences from volumes not on disk will be absent.
# See Phase A ELN for reactivation criteria for a full strategic rebuild.
#
# Exit codes:
#   0  — mammalia_nt DB built and nt_euk cleaned up
#   1  — preflight failure, extraction error, or verification failure
# ─────────────────────────────────────────────────────────────────────────────

# ─── Login-node guard ────────────────────────────────────────────────────────
HOSTNAME_VAL="$(hostname)"
if [[ "$HOSTNAME_VAL" == login* ]] && [[ -z "${SLURM_JOB_ID:-}" ]]; then
  printf "[ERROR] Refusing to run on login node '%s'.\n" "$HOSTNAME_VAL" >&2
  printf "        TACC forbids apptainer on login nodes. Start a compute session:\n" >&2
  printf "          idev -t 02:00:00 -p development\n" >&2
  printf "          module load tacc-apptainer\n" >&2
  printf "          bash scripts/extract_mammalia_from_partial_nt_euk.sh\n" >&2
  exit 1
fi

if ! command -v apptainer &>/dev/null; then
  printf "[ERROR] apptainer not found in PATH. Run: module load tacc-apptainer\n" >&2
  exit 1
fi

set -u

SIF="${SIF:-$WORK/Metabarcoding/Apptainer/edna-pipeline.sif}"
BLASTDB_DIR="${BLASTDB_DIR:-$WORK/Metabarcoding/blastdb}"
OUT_PREFIX="${BLASTDB_DIR}/mammalia_nt"
PARTIAL_NAL="${BLASTDB_DIR}/nt_euk_partial.nal"
MAMMALIA_FASTA="${BLASTDB_DIR}/nt_euk_mammalia.fasta"
MAMMALIA_TAXID=40674
FASTA_MIN_KB=1048576   # 1 GB lower bound on extracted FASTA size

if [[ ! -f "$SIF" ]]; then
  printf "[ERROR] Container not found: %s\n" "$SIF" >&2
  exit 1
fi

printf "╔══════════════════════════════════════════════════════════════╗\n"
printf "║  extract_mammalia_from_partial_nt_euk.sh — Phase A recovery ║\n"
printf "╚══════════════════════════════════════════════════════════════╝\n"
printf "[INFO] BLASTDB dir : %s\n" "$BLASTDB_DIR"
printf "[INFO] Container   : %s\n" "$SIF"

export BLASTDB="$BLASTDB_DIR"

# ═══════════════════════════════════════════════════════════════════════════════
# Step 1: Inventory — verify extracted volumes and taxonomy aux files
# ═══════════════════════════════════════════════════════════════════════════════

printf "\n[INFO] Step 1: Inventory\n"
cd "$BLASTDB_DIR"

mapfile -t EXTRACTED_VOLS < <(ls nt_euk.*.nhr 2>/dev/null | sed 's/\.nhr$//' | sort)
N_VOLS="${#EXTRACTED_VOLS[@]}"

if (( N_VOLS == 0 )); then
  printf "[ERROR] No extracted nt_euk volumes (.nhr) found in %s\n" "$BLASTDB_DIR" >&2
  printf "        Nothing to extract Mammalia from. Aborting.\n" >&2
  exit 1
fi

printf "[OK]   %d extracted volume(s): %s … %s\n" \
  "$N_VOLS" "${EXTRACTED_VOLS[0]}" "${EXTRACTED_VOLS[-1]}"

MISSING_AUX=()
for f in taxdb.btd taxdb.bti taxonomy4blast.sqlite3; do
  [[ -f "${BLASTDB_DIR}/${f}" ]] || MISSING_AUX+=("$f")
done
if (( ${#MISSING_AUX[@]} > 0 )); then
  printf "[ERROR] Missing taxonomy aux files required for -taxids %d:\n" "$MAMMALIA_TAXID" >&2
  for f in "${MISSING_AUX[@]}"; do printf "          - %s\n" "$f" >&2; done
  printf "        These were downloaded in job 3138888 but appear absent. Aborting.\n" >&2
  exit 1
fi
printf "[OK]   Taxonomy aux files present\n"

# Remove leftover .tar.gz for volumes that are already extracted
printf "[INFO] Removing leftover tarballs for extracted volumes...\n"
TARBALLS_REMOVED=0
for vol in "${EXTRACTED_VOLS[@]}"; do
  if [[ -f "${vol}.tar.gz" ]]; then
    rm -f "${vol}.tar.gz"
    TARBALLS_REMOVED=$((TARBALLS_REMOVED + 1))
  fi
done
# Remove any remaining unextracted/partial tarballs (vol 186 was truncated)
for tgz in nt_euk.*.tar.gz; do
  [[ -f "$tgz" ]] || continue
  printf "[INFO] Removing unextracted/partial tarball: %s\n" "$tgz"
  rm -f "$tgz"
  TARBALLS_REMOVED=$((TARBALLS_REMOVED + 1))
done
printf "[OK]   Removed %d tarball(s)\n" "$TARBALLS_REMOVED"

# Remove stale intermediate FASTA from prior partial run (misleads _nt_euk_present checks)
if [[ -f "$MAMMALIA_FASTA" ]]; then
  printf "[INFO] Removing stale intermediate FASTA from prior run\n"
  rm -f "$MAMMALIA_FASTA"
fi

# ═══════════════════════════════════════════════════════════════════════════════
# Step 2: Build partial alias file
# ═══════════════════════════════════════════════════════════════════════════════

printf "\n[INFO] Step 2: Write nt_euk_partial.nal\n"

{
  printf "TITLE Partial nt_euk (%d volumes — 2026-05-02 disk-quota recovery)\n" "$N_VOLS"
  printf "DBLIST"
  for vol in "${EXTRACTED_VOLS[@]}"; do printf " %s" "$vol"; done
  printf "\n"
} > "$PARTIAL_NAL"

printf "[OK]   %s written (%d volumes)\n" "$PARTIAL_NAL" "$N_VOLS"

# ═══════════════════════════════════════════════════════════════════════════════
# Step 3: Extract Mammalia FASTA
# ═══════════════════════════════════════════════════════════════════════════════

printf "\n[INFO] Step 3: Extract Mammalia (taxid %d) — per-volume loop (%d volumes); 30–90 min\n" \
  "$MAMMALIA_TAXID" "$N_VOLS"

# Per-volume loop avoids the BLAST4/BLAST5 alias blocker: the container's blastdbcmd
# (BLAST+ ≥2.12) looks for .ndb files when resolving a .nal alias, but these volumes
# were created in BLAST4 format (.nhr/.nin/.nsq only). Individual volumes are addressed
# directly and do not trigger the alias resolver.
: > "$MAMMALIA_FASTA"

VOLS_OK=0
VOLS_FAIL=0
for vol in "${EXTRACTED_VOLS[@]}"; do
  apptainer exec \
    --bind "$WORK":"$WORK" \
    --env "BLASTDB=$BLASTDB_DIR" \
    "$SIF" \
    blastdbcmd \
      -db "${BLASTDB_DIR}/${vol}" \
      -taxids "$MAMMALIA_TAXID" \
      -outfmt "%f" \
    >> "$MAMMALIA_FASTA" 2>/dev/null \
    && VOLS_OK=$((VOLS_OK + 1)) \
    || { printf "[WARN] blastdbcmd failed on %s — skipping\n" "$vol" >&2
         VOLS_FAIL=$((VOLS_FAIL + 1)); }
done

printf "[INFO] Volumes processed: %d OK, %d skipped\n" "$VOLS_OK" "$VOLS_FAIL"

if (( VOLS_FAIL == N_VOLS )); then
  printf "[ERROR] All %d volumes failed — blastdbcmd cannot read BLAST4 individual volumes.\n" "$N_VOLS" >&2
  printf "        Phase A is not viable. Proceed to Phase B (MIDORI2 / EDirect).\n" >&2
  exit 1
fi

NSEQ=$(grep -c '^>' "$MAMMALIA_FASTA" 2>/dev/null || echo 0)
FASTA_KB=$(du -sk "$MAMMALIA_FASTA" 2>/dev/null | cut -f1 || echo 0)
printf "[OK]   Extracted %s sequences (FASTA: %s KB)\n" "$NSEQ" "$FASTA_KB"

if (( NSEQ == 0 )); then
  printf "[ERROR] Extraction produced 0 sequences — all volumes may lack Mammalia entries\n" >&2
  exit 1
fi
if (( FASTA_KB < FASTA_MIN_KB )); then
  printf "[WARN] FASTA is %d KB, below the 1 GB lower bound (%d KB).\n" "$FASTA_KB" "$FASTA_MIN_KB" >&2
  printf "       Proceeding, but review blast hits carefully before P3 smoke test.\n" >&2
fi

# ═══════════════════════════════════════════════════════════════════════════════
# Step 4: Build mammalia_nt BLAST DB
# ═══════════════════════════════════════════════════════════════════════════════

printf "\n[INFO] Step 4: makeblastdb — mammalia_nt\n"

if [[ -f "${OUT_PREFIX}.nhr" ]] && [[ "${OUT_PREFIX}.nhr" -nt "$MAMMALIA_FASTA" ]]; then
  printf "[SKIP] mammalia_nt.nhr exists and is newer than the FASTA — skipping makeblastdb\n"
else
  apptainer exec \
    --bind "$WORK":"$WORK" \
    "$SIF" \
    makeblastdb \
      -in "$MAMMALIA_FASTA" \
      -dbtype nucl \
      -parse_seqids \
      -out "$OUT_PREFIX" \
      -title "Mammalia_nt (partial coverage; Hominidae excluded at query time)"
  printf "[OK]   makeblastdb complete\n"
fi

# ═══════════════════════════════════════════════════════════════════════════════
# Step 5: Verify DB
# ═══════════════════════════════════════════════════════════════════════════════

printf "\n[INFO] Step 5: Verify mammalia_nt\n"

MISSING_DB=0
for ext in nhr nin nsq; do
  if [[ -f "${OUT_PREFIX}.${ext}" ]]; then
    SIZE=$(du -sh "${OUT_PREFIX}.${ext}" | cut -f1)
    printf "[OK]   mammalia_nt.%s (%s)\n" "$ext" "$SIZE"
  else
    printf "[ERROR] mammalia_nt.%s not found\n" "$ext" >&2
    MISSING_DB=$((MISSING_DB + 1))
  fi
done

if (( MISSING_DB > 0 )); then
  printf "[ERROR] Verification failed — %d index file(s) missing. NOT cleaning up nt_euk.\n" "$MISSING_DB" >&2
  exit 1
fi

# ═══════════════════════════════════════════════════════════════════════════════
# Step 6: Cleanup nt_euk (only reached after verification passes)
# ═══════════════════════════════════════════════════════════════════════════════

printf "\n[INFO] Step 6: Deleting nt_euk volumes and intermediates\n"

FREED=0
for ext in nhr nin nnd nni nog nsq; do
  for f in nt_euk."$ext" nt_euk.*."$ext"; do
    [[ -f "$f" ]] || continue
    rm -f "$f"
    FREED=$((FREED + 1))
  done
done
for f in nt_euk.nal nt_euk_partial.nal nt_euk.files.list nt_euk.links.list; do
  [[ -f "$f" ]] && { rm -f "$f"; FREED=$((FREED + 1)); }
done
[[ -f "$MAMMALIA_FASTA" ]] && { rm -f "$MAMMALIA_FASTA"; FREED=$((FREED + 1)); }

printf "[OK]   Removed %d file(s)\n" "$FREED"
printf "[INFO] Kept: taxdb.btd, taxdb.bti, taxonomy4blast.sqlite3 (required by blastn -negative_taxids)\n"

printf "\n[INFO] Post-cleanup disk usage:\n"
du -sh "$BLASTDB_DIR" 2>/dev/null || true
df -h "$WORK" 2>/dev/null | tail -1 || true

printf "\n[OK]   Phase A complete.\n"
printf "       DB prefix : %s\n" "$OUT_PREFIX"
printf "       Sequences : %s\n" "$NSEQ"
printf "\n"
printf "Next — P3 smoke test:\n"
printf "  bash scripts/pipeline.sh \\\n"
printf "    eDNA/2024/barcode01 \\\n"
printf "    eDNA/2024_barcode01_trimmed_asv \\\n"
printf "    eDNA/2024_barcode01_consensus_asv \\\n"
printf "    \$WORK/Metabarcoding/blastdb/mammalia_nt \\\n"
printf "    16 false misc/primers.txt asv\n"
