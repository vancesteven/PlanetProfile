#!/bin/bash
# Shared environment for PlanetProfile SLURM jobs. Sourced by every
# *.sbatch in this directory; not meant to be run directly.
#
# EDIT THESE FOR YOUR CLUSTER (or export them before sbatch to override --
# every value below is `${VAR:-default}` so the submit environment wins):
#
#   PP_REPO          repo checkout the jobs run from
#   PP_CONDA_ENV     conda/mamba env name (PPcl locally)
#   PP_CONDA_SH      path to conda.sh for `conda activate` in batch shells
#   PP_MODULES       space-separated `module load` args, empty to skip
#   PP_SCRATCH       fast node-local or parallel scratch for big artifacts
#
# The account/partition/QoS are NOT set here: they are per-cluster and per
# allocation, so they live in each .sbatch header (edit the #SBATCH lines)
# or come from your SLURM defaults.

set -euo pipefail

PP_REPO="${PP_REPO:-$SLURM_SUBMIT_DIR}"
PP_CONDA_ENV="${PP_CONDA_ENV:-PPcl}"
PP_CONDA_SH="${PP_CONDA_SH:-$HOME/mamba/etc/profile.d/conda.sh}"
PP_MODULES="${PP_MODULES:-}"

# Scratch: prefer SLURM's per-job node-local dir when the site provides it
# (fastest, auto-cleaned), else a per-job dir under a cluster scratch root.
if [[ -n "${SLURM_TMPDIR:-}" ]]; then
  PP_SCRATCH="${PP_SCRATCH:-$SLURM_TMPDIR/pp_$SLURM_JOB_ID}"
else
  PP_SCRATCH="${PP_SCRATCH:-${SCRATCH:-/tmp}/pp_${SLURM_JOB_ID:-$$}}"
fi
mkdir -p "$PP_SCRATCH"

if [[ -n "$PP_MODULES" ]]; then
  # shellcheck disable=SC2086
  module load $PP_MODULES
fi

if [[ -f "$PP_CONDA_SH" ]]; then
  # shellcheck disable=SC1090
  source "$PP_CONDA_SH"
  conda activate "$PP_CONDA_ENV"
else
  echo "WARNING: $PP_CONDA_SH not found; assuming '$PP_CONDA_ENV' is " \
       "already on PATH" >&2
fi

# --- hazards this project has actually been bitten by ------------------
# libomp double-load: PlanetProfile (numba) and torch each ship their own
# OpenMP runtime. Loading both in ONE process aborts. Keep dataset/cache
# builds and flow training in SEPARATE jobs, and set this so a job that
# does import both fails soft rather than dying.
export KMP_DUPLICATE_LIB_OK="${KMP_DUPLICATE_LIB_OK:-TRUE}"

# numba needs a writable cache dir; the default under a batch scheduler is
# often read-only or shared-and-racing between array tasks. Give every task
# its OWN dir -- two tasks compiling into one cache dir can corrupt it.
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-$PP_SCRATCH/numba_cache}"
mkdir -p "$NUMBA_CACHE_DIR"

# Thread pinning: these jobs get their parallelism from SLURM array tasks
# or an explicit --n-jobs, so per-process BLAS threading is oversubscription
# (measurably slower, and it fights the scheduler's cpu binding).
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

export PYTHONPATH="${PP_REPO}:${PYTHONPATH:-}"
export MPLBACKEND="${MPLBACKEND:-Agg}"   # no display on compute nodes

cd "$PP_REPO"

echo "=== PlanetProfile SLURM job ==============================="
echo "  job/array : ${SLURM_JOB_ID:-n/a} / ${SLURM_ARRAY_TASK_ID:-n/a}"
echo "  node      : $(hostname)"
echo "  repo      : $PP_REPO"
echo "  git       : $(git -C "$PP_REPO" rev-parse --short HEAD 2>/dev/null \
                      || echo 'not a git checkout')"
echo "  env       : $PP_CONDA_ENV ($(python -c 'import sys; print(sys.version.split()[0])'))"
echo "  scratch   : $PP_SCRATCH"
echo "  cpus/task : ${SLURM_CPUS_PER_TASK:-1}"
echo "==========================================================="
