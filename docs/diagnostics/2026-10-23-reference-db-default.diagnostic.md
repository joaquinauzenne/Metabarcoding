# Reference DB default — diagnostic record

- **Date:** 2026-10-23
- **Status:** Decided (default restored); MIDORI2 validation pending
- **Branch:** `claude/repo-cleanup`

## Finding

`config/pipeline.yaml` had `blast.db_name: mammalia_midori` since the 2026-07-09
checkpoint commit (`33e1d46`), and the README had been rewritten around
`build_midori_db.sh`. No chat, ELN, or diagnostic record documents a decision to
adopt MIDORI2. The only recorded analysis (2026-06-11) recommended an
EDirect/CRABS-built 16S database as primary with MIDORI2 as a cross-check, and a
2026-09-23 status check concluded MIDORI2 was "evaluated but never adopted".
The Phase B0 probes (`diag_midori_blast.sh`, `diag_edirect_tacc.sh`) have no
recorded pass/fail. The validated species calls to date came from NCBI web
BLAST; Phase 4 did not run locally until the 2026-09-07 `-taxids` fix, and the
2026-08-04 end-to-end check used the smoke-test fixture DB, not MIDORI2.

## Decision

- Default stays the nt_euk → Mammalia subset (`mammalia_nt`,
  `scripts/mammalia_nt_db.sh`). The public README documents this method.
- MIDORI2 (`scripts/build_midori_db.sh`) is a **candidate**, not the default.
- Leaning toward EDirect/CRABS over MIDORI2 because it keeps native NCBI
  taxids and matches the existing nt-based method; not yet decided.

## Gate to switch `db_name`

1. `scripts/diag_midori_blast.sh` (or `diag_edirect_tacc.sh`) exits 0, with
   the `staxids`/`sscinames` columns populated — record the run output here.
2. `scripts/smoke_test.sh` passes against the candidate DB.
3. Assignments on the validated sample set (`examples/` web-BLAST export)
   are compared with the candidate DB's calls and the discordances explained.
4. The 16Smam amplicon-length discrepancy (~135 bp nominal vs ~300–330 bp
   observed) is resolved before any amplicon-trimmed (CRABS) build.

## Not verified in this session

The MIDORI2 probe could not be run from the cloud session or the linked
device VM: reference-midori.info and ftp.ncbi.nlm.nih.gov are unreachable
from both, and BLAST+ is not installed in either. Run it in WSL or an
`idev` session.
