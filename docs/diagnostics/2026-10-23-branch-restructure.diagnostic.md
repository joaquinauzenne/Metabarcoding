# Branch restructure — diagnostic record

- **Date:** 2026-10-23
- **Status:** Executed and validated
- **Run by:** joaquinauzenne, from a local clone (Git Bash on Windows), using `restructure-branches.sh`

## Goal

Make the branch layout reflect pipeline versions:

- `feature/platform-consolidation` becomes `main`.
- The previous `main` is kept as version `v1.0`.
- `claude/project-plan-review-lpa4ln` is merged in and then deleted.

## Decisions

| Question | Decision |
| --- | --- |
| Where should the review branch be merged before deletion? | Into both `main` and `feature/platform-consolidation` |
| How is the old `main` preserved? | An annotated tag only (`v1.0`), with no archive branch |
| How does platform-consolidation become `main`? | A merge that keeps `main`'s history, not a force-reset |
| What should `main`'s files be afterwards? | Exactly `feature/platform-consolidation`'s files (a reverse-ours merge) |
| Where does `v1.0` point? | At `main` *before* the review merge, so it is the pure old pipeline |
| Do any remote clones or job scripts track `main`? | No |

## Pre-commit assertion (as confirmed)

- **Change:**
  1. Tag the old `main` as `v1.0`.
  2. Merge the review branch into both branches.
  3. Make `main` identical to `feature/platform-consolidation` with a reverse-ours merge that keeps history.
  4. Delete the review branch.
- **Canonical artifact:** The git repository (`origin` = `github.com/joaquinauzenne/Metabarcoding`). No project documents were edited.
- **Edit mode:** Branch operations only.
  1. Record the old `main` SHA and tag it `v1.0`.
  2. Run `git merge --no-ff` to merge the review branch into each of the two branches.
  3. Run `git merge -s ours main` on `feature/platform-consolidation`.
  4. Fast-forward `main` to it.
  5. Delete the review branch.
- **Branch/phase scope:** Affects only `main`, `feature/platform-consolidation`, the review branch, and the new tag. Other branches (`develop`, `feature/asv-pipeline`, `feature/edna-viz`, `claude/github-io-redesign-282f18`) and the `backup` remote were not touched. The pipeline's OTU/ASV/DADA2 modes are pipeline settings, not git branches, and were out of scope.
- **Tools and their limits:**
  - `git merge -s ours` records a parent commit without taking any of its changes. Content found only on the old `main` is therefore dropped from the tip and stays reachable only through `v1.0` and the history. This differs from `-X ours`.
  - `--ff-only` refuses to run if the branches have diverged.
  - `git tag -a` creates a pointer that does not move.
- **Validation gate:** The script stops before any deletion unless all four checks pass:
  - `main`'s tree equals `feature/platform-consolidation`'s tree.
  - `v1.0` points at the recorded old `main` SHA.
  - The old `main` is an ancestor of the new `main`.
  - The review branch is an ancestor of the new `main`.

## Pre-run findings

- **Branch name:** The branch is spelled `feature/platform-consolidation`, with a hyphen. The first draft of the script used an underscore and stopped harmlessly at its preflight check.
- **Review branch contents:** One commit, `63d0103`, which bumps `fonttools` and `pyarrow` in `requirements.txt` to fix security issues. It was based on the old `main` tip, `d1969fa`, so merging it also brought `d1969fa`'s README edit into `feature/platform-consolidation`. As a result, nothing from the old `main` was lost from the new `main`.
- **Trial merge:** Merging the review branch into `feature/platform-consolidation` produced no conflicts.
- **Existing tag:** `pre-consolidation-2026-09-07` already pointed at `d1969fa`. `v1.0` now names the same commit as a version.
- **Where it ran:** A cloud session could not push the tag (HTTP 403 policy block), so it was run from the local clone with the owner's own credentials.

## Result

| Ref | Commit | Meaning |
| --- | --- | --- |
| `v1.0` | `d1969fa` | Previous pipeline (pure old `main`) |
| — | `992f13c` | Merge of the review branch into the old `main` |
| — | `a2df0a6` | Merge of the review branch into `feature/platform-consolidation` |
| `main` | `2e3d041` | "Promote feature/platform-consolidation to main"; parents `a2df0a6` and `992f13c` |
| `feature/platform-consolidation` | `2e3d041` | Kept, and identical to `main` |
| `claude/project-plan-review-lpa4ln` | — | Deleted, both locally and on `origin` |

The validation gate passed all four checks.

## Follow-ups

- After promotion, GitHub reported 3 Dependabot alerts on `main` (1 high, 2 moderate). They are addressed in the same pull request as this record; see `2026-10-23-dependabot-requirements.diagnostic.md`.
- `feature/platform-consolidation` can be deleted whenever convenient. It is identical to `main`.
