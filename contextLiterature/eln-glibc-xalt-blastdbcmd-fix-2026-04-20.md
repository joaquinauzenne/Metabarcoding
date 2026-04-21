# 2026-04-20 — Run BLAST+ Inside Apptainer to Avoid TACC XALT/GLIBC Error

## Problem

After the `set -u` fix landed (commit `288143c`) and job 3095653 ran, the
`.err` log showed a new failure after module load:

```text
tacc-singularity module has been deprecated, loading tacc-apptainer instead.
[WARN] nt_euk index files are 420 days old (> 90 days)
/usr/bin/blastdbcmd: /lib/x86_64-linux-gnu/libc.so.6: version `GLIBC_2.25' not found (required by /opt/apps/xalt/xalt/lib64/libxalt_init.so)
```

The `.out` log confirms the script progressed past module load, validated
`nt_euk` was present (though 420 days old), and then died at the `blastdbcmd`
extraction step.

## Diagnosis

TACC injects XALT (an accounting/tracking wrapper) via `LD_PRELOAD` on every
job. When we load `ncbi-blast-plus` and invoke `blastdbcmd`, the binary that
gets picked up is `/usr/bin/blastdbcmd` (system path, not the module's private
bin). XALT's `libxalt_init.so` requires `GLIBC_2.25`, but the system `libc.so.6`
on this node is older and does not export that symbol version. Result:
dynamic-linker fatal error before `blastdbcmd` ever executes.

This is TACC-environment-specific and cannot be fixed by patching our script
in-place (we can't modify system libc, and stripping `LD_PRELOAD` is fragile
and TACC-unfriendly). The pragmatic fix is to stop relying on TACC's
`ncbi-blast-plus` module altogether and run BLAST+ inside our own Apptainer
container, which is already the pattern used by
[scripts/blastdb_init.sh](../scripts/blastdb_init.sh) and
[scripts/batch_job.sh](../scripts/batch_job.sh).

## Fix

In [scripts/mammalia_nt_db.sh](../scripts/mammalia_nt_db.sh):

1. **Drop** `module load biocontainers` and `module load ncbi-blast-plus`.
   Only `tacc-apptainer` is needed — the container provides BLAST+.
2. **Wrap** `blastdbcmd` and `makeblastdb` calls in `apptainer exec` with
   `--bind "$WORK":"$WORK"` and `--env "BLASTDB=$BLASTDB_DIR"` (so the
   in-container BLAST+ can resolve taxid DBs in `$WORK/Metabarcoding/blastdb`).
3. **Pre-flight** check that `$WORK/Metabarcoding/Apptainer/edna-pipeline.sif`
   exists before trying to use it, with a clear error if not.

Key change, `scripts/mammalia_nt_db.sh`:

```bash
module load tacc-apptainer

set -u

SIF="$WORK/Metabarcoding/Apptainer/edna-pipeline.sif"
if [[ ! -f "$SIF" ]]; then
  printf "[ERROR] edna-pipeline.sif not found at %s\n" "$SIF" >&2
  exit 1
fi

# ... later ...

apptainer exec \
  --bind "$WORK":"$WORK" \
  --env "BLASTDB=$BLASTDB_DIR" \
  "$SIF" \
  blastdbcmd -db "${DB_PREFIX}" -taxids 40674 -negative_taxids 9604 \
    -outfmt "%f" -out "${BLASTDB_DIR}/nt_euk_mammalia.fasta"

apptainer exec \
  --bind "$WORK":"$WORK" \
  "$SIF" \
  makeblastdb -in "${BLASTDB_DIR}/nt_euk_mammalia.fasta" \
    -dbtype nucl -parse_seqids \
    -out "${OUT_PREFIX}" -title "Mammalia_nt_excluding_Hominidae"
```

## Notes

- `batch_job.sh` and `blastdb_init.sh` already run BLAST+ through Apptainer,
  so they are not affected by this bug.
- The stale-index warning (`420 days old`) is informational only; the user
  explicitly chose option `[1]` (use existing) in the interactive prompt,
  and `blastdbcmd` was about to succeed against the existing indexes until
  XALT/GLIBC aborted it.

## Git

- Scope: `scripts/mammalia_nt_db.sh` only
- Commit pushed 2026-04-20

## Next Steps

1. `git pull` on Lonestar6
2. Resubmit: `bash scripts/mammalia_nt_db.sh` (account + choice `[1]`)
3. Confirm `.err` no longer contains `GLIBC_2.25` or any XALT error, and
   the job completes through the `makeblastdb` step producing
   `blastdb/mammalia_nt.*` index files
4. Proceed to F1 preflight (`setup.sh`) and F2 single-barcode test run
