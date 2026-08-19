# Running PlanetProfile inference workloads on a SLURM cluster

Batch scripts for the compute-heavy campaign steps. Everything here is
cluster-agnostic: no account, partition, or module name is hardcoded —
they come from env vars or the `#SBATCH` lines you edit once.

## First-time setup

1. **Fill in your allocation.** Each `.sbatch` has a commented block:
   ```
   ##SBATCH --account=CHANGEME
   ##SBATCH --partition=CHANGEME
   ```
   Uncomment (single `#SBATCH`) and set them, or rely on your site's
   defaults if it has them.

2. **Point `pp_env.sh` at your environment** (or export these before
   `sbatch` — the submit environment always wins):
   | var | default | meaning |
   |---|---|---|
   | `PP_REPO` | `$SLURM_SUBMIT_DIR` | repo checkout to run from |
   | `PP_CONDA_ENV` | `PPcl` | conda/mamba env name |
   | `PP_CONDA_SH` | `$HOME/mamba/etc/profile.d/conda.sh` | for `conda activate` in batch shells |
   | `PP_MODULES` | *(empty)* | `module load` args, e.g. `"python/3.11 gcc/12"` |
   | `PP_SCRATCH` | `$SLURM_TMPDIR` or `$SCRATCH` | fast scratch for big artifacts |

3. `mkdir -p slurm_logs` in the repo root (the scripts write there).

`pp_env.sh` also sets the things this project has actually been bitten by:
`KMP_DUPLICATE_LIB_OK` (numba+torch libomp double-load), a **per-task**
`NUMBA_CACHE_DIR` (array tasks sharing one cache dir can corrupt it),
BLAS thread pinning to 1 (parallelism comes from the array or an explicit
`--n-jobs`, so per-process threading is oversubscription), and
`MPLBACKEND=Agg`.

> **Keep PlanetProfile and torch in separate jobs.** They each ship an
> OpenMP runtime and loading both in one process aborts. Cache builds and
> dataset generation are one job; flow training is another.

## Jobs

### 1. Enceladus (zb, w) production cache — sharded

The big one: 88 zb × 40 w = 3520 ocean nodes + a 39-node frozen axis.
~10 h as a single serial process; ~1 h wall as an 11-task array.

```bash
# 0. sanity-check the decomposition — costs nothing, no allocation used
python plans/scripts/enceladus_zbw_shard_build.py --n-shards 11 --plan

# 1. build the shards, 2. merge when they all succeed
mkdir -p slurm_logs
JID=$(sbatch --parsable plans/cluster/enceladus_zbw_shards.sbatch)
SHARD_DIR=$SCRATCH/enc_zbw_shards_$JID \
  sbatch --dependency=afterok:$JID plans/cluster/enceladus_zbw_merge.sbatch
```

Pass the **same** `SHARD_DIR` to both. The shard job prints the directory
it chose; `afterok` means the merge never starts if any task failed.

**Why sharding does not weaken the r7 gates.** This is the load-bearing
claim — it rests on the gates' own structure, not on convenience:

- **C2 reachability (two-sided)** and the **zb_placement invariant** are
  enforced inside `build_zbw_grid_cache` with strictly **per-node**
  predicates (each node's own zb/w against the declared floor, each node's
  own residual against `zb_tol_km`; `cache_builder.py` ~2139-2165). No
  cross-node or cross-segment term exists, so running them over a
  partition of the zb axis enforces the *identical* predicate on the
  *identical* 3520 nodes. Distributed, not relaxed.
- **`ocean_moi_window`** is the one aggregate and it aggregates *exactly*:
  a max over nodes (merge = max of maxes) and a count (merge = sum).
- The merge **re-implements no gate**. It asserts every shard actually ran
  C2 (a shard with a `None` restriction table is rejected), that the shard
  set is an exact index-level partition of the canonical 88-node axis,
  that the reassembled axis is *bit-identical* to the config's, and that
  the contract fields (`ocean_comp`, `zb_tol_km`, restriction table, w
  axis) agree bit-for-bit across shards. Any failure → `MERGE REFUSED`,
  exit 1, no cache written.

**The tolerance trap** (why the shard script *requires* an explicit
`--zb-tol-km`): the builder's `zb_tol_km=None` default is "half the
smallest spacing in `zb_km_grid`", and it also sets the Tb root-find's
`solve_tol_km = 0.4 × zb_tol_km`. Left to default, a shard of the coarse
5-20@1.0 km segment would solve to 0.5 km and a shard of the fine
22-30@0.25 segment to 0.125 km — shards solved to *different precision*,
not assemblable into what any single-node run produces. The script
therefore demands the explicit canonical 0.125 and refuses anything else.

Output: the cache plus
`validation_reports/enceladus_isostasy/r7_production_build_report.json`
carrying the six §0.29 acceptance gates. **Report them to Machine A
whether they pass or fail** — gates are interpreted, never tuned.

### 2. Corrected-pipeline SBC (§0.28)

Single node, many cores. The driver already forks a process pool over
independent pairs, and keeping it in one process is deliberate: the
BH-FDR step needs every surviving pair at once.

```bash
# preregistered §5 precondition first (~12 min): fiducial reproduction
VERIFY_ONLY=1 sbatch plans/cluster/nh3_corrected_sbc.sbatch
# then the real run
sbatch plans/cluster/nh3_corrected_sbc.sbatch
```

Measured NH3 reference: 613 pairs ≈ 63 CPU-h → ~2 h at 32 cores (was ~14 h
at 8). **Scale `--mem` with `--cpus-per-task`** — each worker held ~1.4 GB,
so 32 workers need ~45 GB; under-requesting gets workers OOM-killed and
their pairs recorded as errors.

Point `DRIVER` at another composition's driver to reuse the harness.

### 3. Multi-seed reference sampler runs

One array task per (config, seed) — built for r7's post-build step (2
branch configs × 3 seeds = 6 tasks, the default `--array=0-5`), and it
also serves the B3 reference-wander protocol (1 config × ≥3 seeds).

```bash
CONFIGS="path/to/ocean_5param.json path/to/frozen_3param.json" \
SEEDS="8001 8002 8003" \
  sbatch plans/cluster/reference_runs.sbatch
```

`run_inference_cli` takes its seed from `config.random_state`, not a flag,
so each task writes a **scratch-only** seed-patched copy of the committed
config (`plans/scripts/cluster_seeded_config.py`). The committed config
stays the source of truth; every derived config records its parent and the
exact keys changed under `metadata.cluster_seed_patch`, and the config is
copied next to its result for provenance. `--set` refuses to *create* a
key, so a mistyped override fails loudly instead of becoming a silently
ignored setting.

Extra sampler overrides:
```bash
EXTRA_SET="--set sampler_settings.n_effective=8000" \
CONFIGS=... SEEDS=... sbatch plans/cluster/reference_runs.sbatch
```

> `--time=24:00:00` is extrapolated from this project's n_eff=2000 runs;
> n_eff=8000 with n_active=4096 is substantially heavier. Check the first
> task's rate early — **a sampler killed at the wall leaves a partial
> chain, which must not be pooled as if it had converged.**

## Wall-time and memory figures

Every number in these scripts is either measured on this project's own
runs or extrapolated from them, and says which in the script comments.
Re-measure on your hardware before trusting the array-wide requests:
node-local scratch vs a parallel filesystem alone moves the cache-build
per-node cost noticeably.

## Staging and Dropbox

The workstation checkout lives in Dropbox, so the jobs build large
artifacts in `$PP_SCRATCH` and `mv` them into place at the end. Writing a
multi-GB pickle straight into a synced folder invites a partially-synced
file. Harmless on a pure-cluster path, load-bearing on the shared one.
