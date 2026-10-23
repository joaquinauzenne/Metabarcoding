# Dependabot alerts in requirements.txt — diagnostic record

- **Date:** 2026-10-23
- **Status:** Proposed in a pull request (not merged when written)

## Trigger

After the branch restructure (see `2026-10-23-branch-restructure.diagnostic.md`), GitHub reported 3 open Dependabot alerts on `main`: 1 high and 2 moderate. All of them are in the root `requirements.txt`. Dependabot scans pinned (`==`) pip entries. It does not scan conda `environment.yml`, R `DESCRIPTION` files or Dockerfile package installs.

## Findings

| Package | Pinned version | Issue | Fixed in |
| --- | --- | --- | --- |
| fonttools | 4.60.2 | `eval()` in `cffLib.parseBlendList()`: arbitrary code execution from a crafted TTX file (fonttools #4039, fixed by PR #4040) | 4.62.0 |
| fonttools | 4.60.2 | `fontTools.misc.etree` resolves XML external entities when lxml is older than 5.0, which allows a local file read (AIKIDO-2026-755272) | 4.64.0 |
| hpack | 4.0.0 | Unbounded variable-integer decoding gives O(n²) run-away computation on malformed input (Dependabot alert #3, affects versions ≤ 4.1.0) | 4.2.0 |
| pyarrow | 23.0.1 | Already patched for CVE-2026-25087 / GHSA-rgxp-2hwp-jwgg (fixed in 23.0.1) | — |

The two fonttools rows are inferred from public advisories. This session could not read the alert list itself; only the hpack alert was confirmed by the owner.

## Why requirements.txt is retired rather than kept

- **It is a `pip freeze` snapshot, not a spec.** 92 of its 109 lines point at conda build folders (`file:///home/conda/feedstock_root/...`, `file:///opt/conda/conda-bld/...`) that do not exist on any other machine.
- **It records Python 3.8,** through the `scipy` cp38 wheel line. `environment.yml` specifies Python 3.10, and `pyarrow==23.0.1` (added by the earlier security fix) does not support 3.8.
- **No script or container uses it.**
  - `Apptainer/edna-pipeline.def` and `Docker/Dockerfile` install their tools directly with mamba/micromamba.
  - `scripts/` does not reference it.
  - Only the README referenced it, and the README's command, `conda create --file requirements.txt`, cannot read a pip-freeze file.
- **`environment.yml` calls itself the canonical conda spec** and documents its own lock-file workflow (`environment.lock.txt` generated from the built container).

## Pre-commit assertion

- **Change:**
  1. Move `requirements.txt` to `misc/legacy/requirements.txt` and add a header that marks it legacy and non-installable.
  2. In that file, bump `fonttools` from 4.60.2 to 4.64.0 and `hpack` from 4.0.0 to 4.2.0.
  3. Point the README's tree, Dependencies line and conda setup step at `environment.yml`.
  4. Add this record and the branch-restructure record.
- **Canonical artifacts:** `requirements.txt` (moved to `misc/legacy/`), `README.md`, and `docs/diagnostics/`.
- **Edit mode:** `git mv` plus two surgical pin changes and a header comment; 4 surgical regions in the README; 2 new files.
- **Branch/phase scope:** Docs and packaging metadata only. No pipeline script, parameter, container definition or Golem module is changed, and no Phase 1–4 behaviour or OTU/ASV/DADA2 mode changes.
- **Tools and their limits:**
  - Dependabot raises alerts from version-pinned manifests. It cannot see unpinned conda specs or `file://` entries, so clearing these alerts says nothing about the container's actual package versions.
  - `mamba env create -f environment.yml` solves the environment fresh. Because the spec is unpinned, the versions it installs are not reproducible until `environment.lock.txt` is committed.
- **Execution environment:** No runtime change. The containers do not read either file.
- **Validation gate:**
  - After merge, Dependabot's alert list for `requirements.txt` should show 0 open alerts. They auto-close when the manifest path disappears; Dependabot still scans the moved file with the patched pins.
  - If any alert remains open, compare its package against `misc/legacy/requirements.txt`.
- **Cross-document consistency:**
  - Updated: the README sections (tree, Dependencies, Setup step 2).
  - Not changed: `environment.yml`, `Apptainer/edna-pipeline.def`, `Docker/Dockerfile`, and `edna-viz/` (R, which Dependabot does not scan).
- **Doubts:**
  1. The patched pins in the legacy file have not been install-tested. The file cannot be installed anyway, and recent `hpack` and `fonttools` releases may not support the Python 3.8 the snapshot came from.
  2. Dependabot's own hpack PR (#2) was merged into `main` while this change was in review. It made the same `hpack` 4.2.0 bump at the old path, so this change was rebased onto it, and the move to `misc/legacy/` now carries that bump.
  3. Whether the built container ships a vulnerable `fonttools` (pulled in through `matplotlib`/`NanoPlot`) is not checked here. It would show up in `environment.lock.txt` once one is committed.

## Defer (suggestions only, not part of this change)

- Commit `environment.lock.txt` from the built container, so pinned versions exist that can be audited.
