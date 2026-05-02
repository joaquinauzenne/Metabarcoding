#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ─────────────────────────────────────────────────────────────────────────────
# Script : diag_edirect_tacc.sh
# Purpose: Phase B0 de-risk probe — verify that NCBI EDirect works on a TACC
#          compute node before committing to an EDirect-based Mammalia DB
#          rebuild strategy (Phase B2/B3).
#
#          EDirect uses HTTPS to eutils.ncbi.nlm.nih.gov (port 443), which is
#          a distinct channel from the FTP (port 21) that update_blastdb.pl
#          used and that TACC blocked. This script confirms the HTTPS path is
#          open.
#
# Run inside an idev compute session (NOT on a login node):
#   idev -t 00:15:00 -p development
#   bash scripts/diag_edirect_tacc.sh
#
# Exit codes:
#   0  — EDirect installed and esearch probe returned a positive result
#   1  — installation failed, network blocked, or probe returned 0/empty
# ─────────────────────────────────────────────────────────────────────────────

EDIRECT_DIR="${EDIRECT_DIR:-$HOME/edirect}"

# ─── Login-node guard ────────────────────────────────────────────────────────
HOSTNAME_VAL="$(hostname)"
if [[ "$HOSTNAME_VAL" == login* ]] && [[ -z "${SLURM_JOB_ID:-}" ]]; then
  printf "[ERROR] Refusing to run on login node '%s'.\n" "$HOSTNAME_VAL" >&2
  printf "        TACC policy: run network-bound tools on compute nodes only.\n" >&2
  printf "          idev -t 00:15:00 -p development\n" >&2
  printf "          bash scripts/diag_edirect_tacc.sh\n" >&2
  exit 1
fi

if ! command -v curl &>/dev/null; then
  printf "[ERROR] curl not found. EDirect install requires curl.\n" >&2
  exit 1
fi

# ─── Install EDirect if not already present ──────────────────────────────────
printf "[INFO] EDirect dir: %s\n" "$EDIRECT_DIR"

if [[ -x "${EDIRECT_DIR}/esearch" ]]; then
  printf "[OK]   EDirect already installed\n"
else
  printf "[INFO] Installing EDirect to %s...\n" "$EDIRECT_DIR"
  sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)" \
    2>&1 | grep -v "^$" || true
  if [[ ! -x "${EDIRECT_DIR}/esearch" ]]; then
    printf "[ERROR] EDirect install failed — esearch not found in %s\n" "$EDIRECT_DIR" >&2
    printf "        Check: curl https://ftp.ncbi.nlm.nih.gov to verify HTTPS reachability\n" >&2
    exit 1
  fi
  printf "[OK]   EDirect installed\n"
fi

export PATH="${EDIRECT_DIR}:${PATH}"

# ─── NCBI API key (optional, reduces rate-limit risk) ────────────────────────
if [[ -n "${NCBI_API_KEY:-}" ]]; then
  export NCBI_API_KEY
  printf "[INFO] NCBI_API_KEY set — using authenticated requests\n"
else
  printf "[INFO] NCBI_API_KEY not set — anonymous requests (may hit rate limits)\n"
fi

# ─── Probe: count mammalian complete genomes in GenBank ──────────────────────
printf "[INFO] Probing: esearch for Mammalia complete genomes in nucleotide DB...\n"

PROBE_RESULT=$(
  esearch -db nucleotide \
          -query "txid40674[Organism] AND complete genome[Title]" \
  2>&1
)

COUNT=$(printf '%s\n' "$PROBE_RESULT" \
  | xtract -pattern ENTREZ_DIRECT -element Count 2>/dev/null || echo "")

if [[ -z "$COUNT" ]]; then
  printf "[ERROR] esearch returned no Count — xtract parse failed or network timed out\n" >&2
  printf "        Raw output:\n" >&2
  printf "%s\n" "$PROBE_RESULT" | head -20 >&2
  exit 1
fi

printf "[INFO] Mammalia complete genome hits: %s\n" "$COUNT"

if (( COUNT > 0 )); then
  printf "[OK]   B0-EDirect PASS — HTTPS to eutils.ncbi.nlm.nih.gov is open\n"
  printf "       EDirect is usable on this compute node for Phase B2/B3.\n"
  exit 0
else
  printf "[FAIL] B0-EDirect FAIL — Count returned 0\n" >&2
  printf "       Either the query matched nothing (unexpected) or the result was malformed.\n" >&2
  printf "       Raw output:\n" >&2
  printf "%s\n" "$PROBE_RESULT" | head -20 >&2
  exit 1
fi
