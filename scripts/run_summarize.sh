#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
NANO_DIR="$PROJECT_ROOT/eDNA/trimmed/nanoplot_results_trimmed"

python3 "$SCRIPT_DIR/summarize_nanostats.py" "$NANO_DIR"