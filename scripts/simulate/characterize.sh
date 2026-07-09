#!/usr/bin/env bash
set -eo pipefail
IFS=$'\n\t'

# ─────────────────────────────────────────────────────────────────────────────
# scripts/simulate/characterize.sh  — M1 step 1 of 4
#
# Purpose: make a NanoSim ONT error model available for the validation harness.
#          Per the M1 decision, we use a PRE-TRAINED (bundled) ONT model rather
#          than characterizing from a pilot run. NanoSim's pre-trained models are
#          NOT shipped inside the pip/conda package — they are downloaded
#          separately from the NanoSim project. This script locates an existing
#          model, or downloads one if you provide MODEL_URL, and records the model
#          prefix for simulate.sh.
#
# Capability / limitation (paired, per project convention):
#   DOES verify NanoSim is callable and pin a model prefix for reproducible sims.
#   Does NOT judge whether the chosen ONT model matches your flowcell/basecaller
#   chemistry — a pre-trained model is an approximation, not your run's true error
#   profile. Swap to a pilot-characterized model later for run-specific realism.
#
# Run inside an idev compute session (download needs network):
#   idev -t 00:20:00 -p development
#   module load tacc-apptainer
#   MODEL_URL="<nanosim pre-trained model tarball url>" bash scripts/simulate/characterize.sh
#
# Env vars:
#   SIF           Container image [default: Apptainer/edna-pipeline.sif]
#   MODELS_DIR    Where models live/download [default: $WORK/Metabarcoding/nanosim_models]
#   MODEL_PREFIX  Path prefix of an already-present model (skips download)
#   MODEL_URL     Tarball URL of a NanoSim pre-trained ONT model (used if no prefix)
#   NANOSIM_ENV   In-container env path [default: /opt/conda/envs/nanosim]
#
# Output: writes the resolved model prefix to $MODELS_DIR/MODEL_PREFIX.txt
# Exit:   0 model ready; 1 NanoSim missing, no model, or download failed.
# ─────────────────────────────────────────────────────────────────────────────

SIF="${SIF:-Apptainer/edna-pipeline.sif}"
MODELS_DIR="${MODELS_DIR:-${WORK:-$PWD}/Metabarcoding/nanosim_models}"
MODEL_PREFIX="${MODEL_PREFIX:-}"
MODEL_URL="${MODEL_URL:-}"
NANOSIM_ENV="${NANOSIM_ENV:-/opt/conda/envs/nanosim}"

# ─── Login-node guard (download/simulation belong on compute nodes) ──────────
HOSTNAME_VAL="$(hostname)"
if [[ "$HOSTNAME_VAL" == login* ]] && [[ -z "${SLURM_JOB_ID:-}" ]]; then
  printf "[ERROR] Refusing to run on login node '%s'. Start an idev session:\n" "$HOSTNAME_VAL" >&2
  printf "          idev -t 00:20:00 -p development\n          module load tacc-apptainer\n" >&2
  exit 1
fi

if ! command -v apptainer &>/dev/null; then
  printf "[ERROR] apptainer not found. Run: module load tacc-apptainer\n" >&2
  exit 1
fi
[[ -f "$SIF" ]] || { printf "[ERROR] Container not found: %s\n" "$SIF" >&2; exit 1; }

mkdir -p "$MODELS_DIR"

# ─── Verify NanoSim is callable in its env ───────────────────────────────────
printf "[INFO] Verifying NanoSim in %s\n" "$NANOSIM_ENV"
if ! apptainer exec "$SIF" "$NANOSIM_ENV/bin/simulator.py" --version >/dev/null 2>&1; then
  printf "[ERROR] simulator.py not runnable at %s/bin/. Rebuild the container from\n" "$NANOSIM_ENV" >&2
  printf "        Apptainer/edna-pipeline.def (adds the 'nanosim' env).\n" >&2
  exit 1
fi
printf "[OK]   NanoSim: %s\n" "$(apptainer exec "$SIF" "$NANOSIM_ENV/bin/simulator.py" --version 2>&1 | head -1)"

# ─── Resolve a model prefix ──────────────────────────────────────────────────
# A NanoSim "model prefix" is the common path stem shared by the model's files
# (e.g. .../human_NA12878_DNA_FAB49712_guppy/training  → prefix ".../training").
if [[ -n "$MODEL_PREFIX" ]]; then
  if compgen -G "${MODEL_PREFIX}*" >/dev/null; then
    printf "[OK]   Using provided MODEL_PREFIX: %s\n" "$MODEL_PREFIX"
  else
    printf "[ERROR] No files found for MODEL_PREFIX=%s*\n" "$MODEL_PREFIX" >&2
    exit 1
  fi
elif [[ -n "$MODEL_URL" ]]; then
  printf "[INFO] Downloading pre-trained model: %s\n" "$MODEL_URL"
  tarball="$MODELS_DIR/$(basename "$MODEL_URL")"
  curl -fSL --retry 3 --retry-delay 5 -o "$tarball" "$MODEL_URL"
  tar -xzf "$tarball" -C "$MODELS_DIR"
  # Heuristic: the model prefix is the stem of the '*_model_profile' file, if any,
  # else the extracted directory + '/training'. Confirm against the model's README.
  profile="$(find "$MODELS_DIR" -name '*_model_profile' | head -1 || true)"
  if [[ -n "$profile" ]]; then
    MODEL_PREFIX="${profile%_model_profile}"
  else
    extracted_dir="$(find "$MODELS_DIR" -maxdepth 1 -type d ! -path "$MODELS_DIR" | head -1 || true)"
    MODEL_PREFIX="${extracted_dir}/training"
  fi
  printf "[INFO] Inferred MODEL_PREFIX: %s (verify against the model's README)\n" "$MODEL_PREFIX"
else
  printf "[ERROR] No model available. Provide one of:\n" >&2
  printf "  MODEL_PREFIX=/path/to/model/stem   (already-downloaded model)\n" >&2
  printf "  MODEL_URL=<pre-trained model tarball url>\n" >&2
  printf "Pre-trained ONT models are listed in the NanoSim project (bcgsc/NanoSim)\n" >&2
  printf "under 'pre-trained_models'. Pick an ONT genomic model (e.g. a guppy DNA model).\n" >&2
  exit 1
fi

printf "%s\n" "$MODEL_PREFIX" > "$MODELS_DIR/MODEL_PREFIX.txt"
printf "[OK]   Model prefix recorded -> %s/MODEL_PREFIX.txt\n" "$MODELS_DIR"
printf "       simulate.sh will read it from there (or pass MODEL_PREFIX explicitly).\n"
