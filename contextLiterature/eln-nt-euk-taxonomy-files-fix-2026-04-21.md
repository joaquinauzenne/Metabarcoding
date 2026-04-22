# ELN — `nt_euk.nal` + taxonomy aux files missing in `mammalia_nt_db.sh`

**Date**: 2026-04-21
**Script**: `scripts/mammalia_nt_db.sh`
**Host**: Lonestar6 (TACC), `login2.ls6`
**Job**: `3095766` (failed after the `-negative_taxids` fix landed)
**Related**: `eln-mammalia-db-negative-taxids-fix-2026-04-20.md`

## Observed errors

After commit `88e64f9` removed the invalid `-negative_taxids` flag, the job reached the `blastdbcmd -taxids 40674` step and failed with two stacked errors:

```
The -taxids command line option requires additional data files. Please see the section
'Taxonomic filtering for BLAST databases' in https://www.ncbi.nlm.nih.gov/books/NBK569839/
BLAST Database error: No alias or index file found for nucleotide database
[/work/.../blastdb/nt_euk] in search path [...]
```

## Diagnosis

Directory listing of `$WORK/Metabarcoding/blastdb/` showed **none** of:

- `nt_euk.nal` — the multi-volume alias (without it, `blastdbcmd -db nt_euk` can't open the DB even though volume files exist)
- `taxdb.btd`, `taxdb.bti` — scientific/common name lookup
- `taxonomy4blast.sqlite3` — lineage DB required for taxid subtree expansion (i.e. resolving "all descendants of Mammalia = 40674")

Two independent bugs:

1. **`_download_nt_euk` filter too narrow**: `awk '/nt_euk.*\.tar\.gz$/'` grabs only per-volume archives, never the standalone `nt_euk.nal` alias and never any taxonomy files.
2. **Verification too shallow**: `_nt_euk_present` only globs `.nhr` files. Presence of volume indices was misreported as "verified" even though the DB was unusable.

`blastdb_init.sh` does fetch `taxdb` but writes into `blastdb/nt/` (a sibling directory used by a different workflow). `mammalia_nt_db.sh` uses `blastdb/` directly, so those files weren't reachable via `BLASTDB`.

## Fix applied

`scripts/mammalia_nt_db.sh`:

- **`_download_nt_euk`** — broadened the awk filter to four classes of files:
  - `nt_euk.NN.tar.gz` (volume archives, as before)
  - `nt_euk.nal` (alias, standalone)
  - `taxdb.tar.gz` → extracts `taxdb.btd` + `taxdb.bti`
  - `taxonomy4blast.sqlite3.gz` → `gunzip` to `taxonomy4blast.sqlite3`
- **`_missing_aux_files`** (new helper) — returns the list of required auxiliary files not present in `$BLASTDB_DIR`. Guards against `set -u` by checking array length before `printf`.
- **Step 1 verification** — after the `nt_euk.*.nhr` check, also calls `_missing_aux_files`; if any are missing, prints each by name with actionable remediation (`DOWNLOAD_DB=true sbatch ...`) and exits.

Net effect: the script now either (a) verifies a complete, queryable DB install, or (b) fails up front with the exact list of missing pieces — never silently proceeds into a broken `blastdbcmd` call.

## Re-run command

```bash
cd $WORK/Metabarcoding
git pull
DOWNLOAD_DB=true sbatch scripts/mammalia_nt_db.sh
```

`DOWNLOAD_DB=true` will rsync the alias, both taxonomy archives, and any missing volumes in one pass. After decompression, the new verification step should report:

```
[OK]   nt_euk volumes, alias, and taxonomy files verified
```

before the `blastdbcmd` extraction starts.

## Design notes

- **Why not switch to `update_blastdb.pl`?** Would also solve this, but the script has an established rsync-based download pattern and the filter fix is a ~15-line change. Kept the existing mechanism; revisit if we hit further edge cases.
- **Why not share files with `blastdb_init.sh`'s `blastdb/nt/`?** That path serves full-nt workflows; Mammalia-specific work keeps its inputs in `blastdb/` top-level. Merging the two is a separate refactor.
- **Forward-looking**: `taxonomy4blast.sqlite3` was flagged in `.gitignore` as a planned file — now it's an actual runtime dependency of Phase 4's taxid filtering, so the gitignore entry is no longer aspirational.
