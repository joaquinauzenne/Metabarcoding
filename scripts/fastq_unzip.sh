#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------
# Script: fastq_unzip.sh
# Purpose: Recursively unzip all .fastq.gz under a given directory (e.g. eDNA)
#          into a parallel '_unzipped' directory tree.
#
# Usage:
#   bash fastq_unzip.sh /path/to/eDNA
#   → creates /path/to/eDNA/_unzipped/[same subdirs]/[files.fastq]
# ----------------------------------------------------------------------

# 1) Check args
if [[ "$#" -ne 1 ]]; then
    echo "Usage: $0 <eDNA_root_directory>"
    exit 1
fi

SOURCE_DIR="$1"
DEST_DIR="$SOURCE_DIR/_unzipped"

# 2) Create the destination base
mkdir -p "$DEST_DIR"

# 3) Find and unzip
find "$SOURCE_DIR" -type f -name "*.fastq.gz" | while IFS= read -r file; do
    # Strip SOURCE_DIR/ prefix to get the relative path
    REL_PATH="${file#${SOURCE_DIR}/}"
    
    # Compute where to write in DEST_DIR
    TARGET_DIR="$DEST_DIR/$(dirname "$REL_PATH")"
    mkdir -p "$TARGET_DIR"
    
    # Unzip into the mirrored location
    OUT_FILE="$TARGET_DIR/$(basename "${REL_PATH%.gz}")"
    echo "Extracting: $file → $OUT_FILE"
    gunzip -c "$file" > "$OUT_FILE"
done

echo "✅ All .fastq.gz files extracted under $DEST_DIR"
