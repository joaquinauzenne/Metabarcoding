# ELN — `blastdbcmd -negative_taxids` fix in `mammalia_nt_db.sh`

**Date**: 2026-04-20
**Script**: `scripts/mammalia_nt_db.sh`
**Host**: Lonestar6 (TACC), node `c303-006.ls6`
**Related**: `eln-mammalia-db-module-fix-2026-04-13.md`, `eln-glibc-xalt-blastdbcmd-fix-2026-04-20.md`

## Observed error

After pulling the XALT/GLIBC fix (commit `64faffc`), `bash scripts/mammalia_nt_db.sh` on Lonestar6 aborted at the Mammalia-extraction step with:

```
Error: Unknown argument: "negative_taxids"
Error:  (CArgException::eInvalidArg) Unknown argument: "negative_taxids"
```

The USAGE printed by `blastdbcmd` listed only `-taxids` and `-taxidlist` — no negative variant.

## Diagnosis

`-negative_taxids` is a **search-time** flag supported by `blastn` (and other search programs), not by `blastdbcmd`. BLAST+ 2.17.0's `blastdbcmd` exposes only:

- `-taxids <comma-list>`
- `-taxidlist <file>`
- `-no_taxid_expansion`

There is no built-in way to subtract a taxon subtree at the DB-extraction stage. Doing it pre-extraction would require expanding the Hominidae subtree offline (taxdump/taxonkit/Entrez) — adding state, network, or container dependencies — for a filter that `blastn` already applies at query time.

## Decision

Delegate Hominidae exclusion entirely to search time. `scripts/pipeline.sh` (Phase 4) already sets:

```bash
FILTER_OPT="-taxids 40674 -negative_taxids 9604"
```

and passes it to `blastn`, where `-negative_taxids` is valid. The Mammalia DB therefore contains all of Mammalia; Hominidae sequences are filtered when queried. Hominidae is a tiny fraction of Mammalia in nt_euk, so the size cost is negligible, and the design is forward-flexible — adding future exclusions (e.g., *Bos taurus* contamination) requires only a flag tweak in `FILTER_OPT`, no DB rebuild.

Rejected alternatives: taxdump/taxonkit-based taxid expansion (container rebuild + ~200MB data), Entrez API (requires compute-node internet + API key), post-extraction FASTA filtering (same taxid-expansion problem elsewhere), Hominidae mask DB (over-engineered).

## Diff applied

`scripts/mammalia_nt_db.sh`:

- **Step 2 header + log line**: clarifies that Hominidae is excluded downstream by `blastn`, not here.
- **blastdbcmd args**: removed `-negative_taxids 9604`. Kept `-taxids 40674`.
- **makeblastdb `-title`**: `"Mammalia_nt_excluding_Hominidae"` → `"Mammalia_nt (Hominidae excluded at query time)"`.
- Added a 3-line comment above Step 2 documenting the design choice and pointing to `pipeline.sh`'s `FILTER_OPT`.

No SLURM, module, or apptainer-wrapper changes — those were stabilized by the 2026-04-13 and 2026-04-20 fixes and are working.

## Re-run command

```bash
cd $WORK/Metabarcoding
git pull
sbatch scripts/mammalia_nt_db.sh          # reuses existing nt_euk
# or, to refresh nt_euk first:
DOWNLOAD_DB=true sbatch scripts/mammalia_nt_db.sh
```

## Expected log progression

```
[INFO] Extracting Mammalia (taxid 40674); Hominidae (9604) excluded at search time via blastn -negative_taxids
[OK]   Extracted <N> sequences
[INFO] Building BLAST DB: mammalia_nt...
[OK]   Mammalia NT database ready: <BLASTDB_DIR>/mammalia_nt.*
```

## End-to-end verification

1. After DB build completes, run Phase 4 of `pipeline.sh` on one barcode.
2. Confirm the resulting BLAST TSV contains **zero** rows with taxid `9604`:
   ```bash
   awk -F'\t' '$7==9604' <phase4.tsv> | wc -l   # expect 0
   ```
3. Confirm expected Mammalia hits appear (e.g., *Alouatta seniculus* in the 2024 test barcodes).
4. Load the TSV into `edna-viz/` — the 10-column CLI schema is unchanged, so the parser accepts it without modification.
