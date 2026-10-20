#!/usr/bin/env bash
# Local (no-SLURM) launcher for the Enceladus (zb, w) production cache.
#
# Committed reconstruction of the uncommitted ~/scratch/enc_zbw_shards/run.sh
# used on 2026-10-20, which was never committed; this one lives in the repo.
# Same build as enceladus_zbw_shards.sbatch (11 shards, frozen axis on shard
# 0), run as local background processes capped at MAX_JOBS concurrent.
#
# Resumable at shard granularity: a shard whose .meta.json already exists is
# skipped, so after a reboot just re-run. An unfinished shard restarts from
# zero. Writes $SHARD_DIR/DONE when every shard has a .meta.json.
#
#   PP_SCRATCH=~/scratch nohup plans/cluster/enceladus_zbw_local.sh \
#       > ~/scratch/enc_zbw_shards.run.out 2>&1 &
#
# Then merge (see plans/STATUS.md for the gate checks):
#   HOME=$PP_SCRATCH/pphome python plans/scripts/enceladus_zbw_merge_shards.py \
#       --shard-dir $PP_SCRATCH/enc_zbw_shards \
#       --out $PP_SCRATCH/enceladus_zbw_production_cache.pkl
#
# Keep PP_SCRATCH OUT of /tmp (wiped on reboot) and out of the Dropbox-synced
# repo (large pickles mid-write).

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PP_SCRATCH="${PP_SCRATCH:-$HOME/scratch}"
SHARD_DIR="${SHARD_DIR:-$PP_SCRATCH/enc_zbw_shards}"
N_SHARDS="${N_SHARDS:-11}"
MAX_JOBS="${MAX_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"
PP_CONDA_ENV="${PP_CONDA_ENV:-PPcl}"
PP_CONDA_SH="${PP_CONDA_SH:-$HOME/mamba/etc/profile.d/conda.sh}"

# Resolve the REAL Perple_X/PlanetProfile cache before HOME is redirected
# (platformdirs.user_cache_dir: ~/Library/Caches on macOS, ~/.cache on Linux).
if [[ "$(uname)" == "Darwin" ]]; then
  REAL_CACHE="$HOME/Library/Caches/PlanetProfile"
  PPHOME_CACHE_PARENT="$PP_SCRATCH/pphome/Library/Caches"
else
  REAL_CACHE="${XDG_CACHE_HOME:-$HOME/.cache}/PlanetProfile"
  PPHOME_CACHE_PARENT="$PP_SCRATCH/pphome/.cache"
fi
if [[ ! -d "$REAL_CACHE" ]]; then
  echo "ERROR: PlanetProfile cache not found at $REAL_CACHE;" \
       "every Perple_X lookup would fail. Aborting." >&2
  exit 1
fi

if [[ -f "$PP_CONDA_SH" ]]; then
  # shellcheck disable=SC1090
  source "$PP_CONDA_SH"
  conda activate "$PP_CONDA_ENV"
fi

# HOME redirect: TidalPy's config under ~/Documents is EPERM for agent
# shells. The redirected HOME must still see the real cache.
mkdir -p "$SHARD_DIR" "$PPHOME_CACHE_PARENT"
ln -sfn "$REAL_CACHE" "$PPHOME_CACHE_PARENT/PlanetProfile"
export HOME="$PP_SCRATCH/pphome"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export KMP_DUPLICATE_LIB_OK=TRUE MPLBACKEND=Agg
export NUMBA_CACHE_DIR="$PP_SCRATCH/numba_cache"
mkdir -p "$NUMBA_CACHE_DIR"

cd "$REPO"
git rev-parse HEAD > "$SHARD_DIR/BUILD_COMMIT"
echo "repo $REPO @ $(cat "$SHARD_DIR/BUILD_COMMIT"); shards -> $SHARD_DIR;" \
     "max $MAX_JOBS concurrent; $(date)"
python plans/scripts/enceladus_zbw_shard_build.py --n-shards "$N_SHARDS" --plan \
  || { echo "ERROR: --plan failed" >&2; exit 1; }

for ((i = 0; i < N_SHARDS; i++)); do
  tag=$(printf 'shard_%03d_of_%03d' "$i" "$N_SHARDS")
  if [[ -f "$SHARD_DIR/$tag.meta.json" ]]; then
    echo "skip $tag (meta exists)"
    continue
  fi
  # poll, not `wait -n`: macOS /bin/bash is 3.2
  while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do sleep 30; done
  frozen=""
  [[ $i -eq 0 ]] && frozen="--with-frozen"
  echo "launch $tag $frozen $(date)"
  python plans/scripts/enceladus_zbw_shard_build.py \
    --shard-index "$i" --n-shards "$N_SHARDS" --out-dir "$SHARD_DIR" $frozen \
    > "$SHARD_DIR/$tag.log" 2>&1 &
done
wait

n_meta=$(ls "$SHARD_DIR"/shard_*_of_*.meta.json 2>/dev/null | wc -l)
if (( n_meta == N_SHARDS )); then
  date > "$SHARD_DIR/DONE"
  echo "all $N_SHARDS shards have meta; DONE $(date)"
else
  echo "only $n_meta/$N_SHARDS shards have meta; check $SHARD_DIR/*.log" >&2
  exit 1
fi
