#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ─────────────────────────────────────────────────────────────────────────────
# Script : check_asv_deps.sh
# Purpose: Verify that the ASV-mode runtime dependencies of `pipeline.sh`
#          (MODE=asv) and the standalone DADA2 path (`pipeline_dada2.R`)
#          are present inside the Apptainer container. Run before launching
#          ASV-mode jobs on Lonestar6 so missing tools surface here, not
#          mid-pipeline.
#
# Usage:
#   module load tacc-apptainer
#   bash scripts/check_asv_deps.sh
#
# Exit code:
#   0  — all required tools and R packages found
#   1  — one or more dependencies missing (counted in the summary)
# ─────────────────────────────────────────────────────────────────────────────

SIF="${SIF:-Apptainer/edna-pipeline.sif}"

if ! command -v apptainer &>/dev/null; then
  printf "[ERROR] apptainer not found in PATH. Run: module load tacc-apptainer\n" >&2
  exit 1
fi
if [[ ! -f "$SIF" ]]; then
  printf "[ERROR] Container image not found: %s\n" "$SIF" >&2
  exit 1
fi

printf "[INFO] Checking ASV runtime deps in %s\n" "$SIF"

MISSING=0

# Shell tools used by pipeline.sh ASV branch (vsearch, NGSpeciesID, medaka)
# and the OTU branch tools the ASV path also relies on (mafft, cd-hit-est, cutadapt).
printf "[INFO] Shell tools:\n"
for t in vsearch NGSpeciesID medaka mafft cd-hit-est cutadapt; do
  if apptainer exec "$SIF" bash -lc "command -v '$t' >/dev/null"; then
    printf "[OK]      %s\n" "$t"
  else
    printf "[MISSING] %s\n" "$t" >&2
    MISSING=$((MISSING + 1))
  fi
done

# R packages used by pipeline_dada2.R and lca_assignment.R.
printf "[INFO] R packages (Bioconductor):\n"
R_REPORT=$(apptainer exec "$SIF" Rscript -e '
suppressMessages({
  for (pkg in c("dada2", "ShortRead", "Biostrings")) {
    tag <- if (requireNamespace(pkg, quietly = TRUE)) "OK" else "MISSING"
    cat(sprintf("%s\t%s\n", tag, pkg))
  }
})
' 2>/dev/null)

while IFS=$'\t' read -r tag pkg; do
  [[ -z "$tag" ]] && continue
  if [[ "$tag" == "OK" ]]; then
    printf "[OK]      %s\n" "$pkg"
  else
    printf "[MISSING] %s\n" "$pkg" >&2
    MISSING=$((MISSING + 1))
  fi
done <<< "$R_REPORT"

printf "[INFO] Summary: %d missing dependency(ies)\n" "$MISSING"
if (( MISSING > 0 )); then
  printf "[FAIL] Resolve missing items before launching ASV-mode jobs.\n" >&2
  printf "       - Shell tools: extend Apptainer/edna-pipeline.def and rebuild the .sif.\n" >&2
  printf "       - R packages : add via BiocManager::install() inside the container,\n" >&2
  printf "                      or extend the .def with the install step and rebuild.\n" >&2
  exit 1
fi

printf "[OK]   ASV runtime dependencies satisfied.\n"