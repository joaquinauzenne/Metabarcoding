#!/usr/bin/env python3
import os
import sys
import glob
import re
import pandas as pd
import matplotlib.pyplot as plt

def parse_nanostats(path):
    """
    Parse the NanoStats.txt at `path` and extract key numeric metrics.
    """
    metrics = {}
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            # Match lines like "Mean read length:    125.0"
            m = re.match(r'^([\w \-]+):\s*([\d,]+(?:\.\d+)?)$', line)
            if m:
                key = m.group(1).strip()
                val = m.group(2).replace(',', '')
                try:
                    metrics[key] = float(val)
                except ValueError:
                    pass
    return metrics

def main(base_dir):
    # 1) Recursively find all NanoStats.txt files
    pattern = os.path.join(base_dir, '**', 'NanoStats.txt')
    paths = glob.glob(pattern, recursive=True)

    if not paths:
        print(f"No NanoStats.txt found under {base_dir}", file=sys.stderr)
        sys.exit(1)

    records = []
    for p in sorted(paths):
        sample = os.path.basename(os.path.dirname(p))
        stats = parse_nanostats(p)
        if stats:
            stats['Sample'] = sample
            records.append(stats)

    if not records:
        print(f"No metrics parsed from NanoStats.txt under {base_dir}", file=sys.stderr)
        sys.exit(1)

    # 2) Build DataFrame and aggregate
    df = pd.DataFrame(records).set_index('Sample')
    summary = df.describe().T.round(2)

    # 3) Plot as table
    fig, ax = plt.subplots(figsize=(10, len(summary)*0.4 + 1))
    ax.axis('off')
    tbl = ax.table(
        cellText=summary.values,
        rowLabels=summary.index,
        colLabels=summary.columns,
        cellLoc='center',
        loc='center'
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1, 1.5)
    plt.tight_layout()

    # 4) Save image
    out_path = os.path.join(base_dir, 'nanostats_summary.png')
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"Saved summary table image to {out_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <NanoStats_directory>", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1])
