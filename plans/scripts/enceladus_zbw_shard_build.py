#!/usr/bin/env python
"""Build ONE zb-shard of the Enceladus (zb, w) production cache.

Cluster counterpart to `enceladus_zbw_production_build.py` (which builds
the whole grid in one ~10h serial call). r7's build instructions
themselves say "ocean ~19.5 h serial (parallelize by segment)"; this is
that parallelization, as a SLURM array task.

WHY SHARDING PRESERVES THE r7 GATES (the load-bearing claim -- read this
before trusting a sharded cache):

- **C2 reachability (two-sided)** is enforced INSIDE
  `build_zbw_grid_cache`, and its predicate is strictly PER-NODE: for
  each (i_zb, i_w) it compares that node's own zb against
  `_reachability_zb_min_interp(w)` and that node's own built/None state
  (cache_builder.py ~2139-2165). There is no cross-node or cross-segment
  term. So running it on a partition of the zb axis enforces the
  IDENTICAL predicate on the IDENTICAL 3520 nodes -- it is distributed,
  not weakened. The merge step asserts every shard actually ran it (a
  shard whose `ocean_reachability_restriction_table` is None did not,
  and is rejected).
- **zb_placement invariant** is likewise per-node (each node's own
  |zb_achieved - zb_target| against zb_tol_km).
- **ocean_moi_window** is the one AGGREGATE, and it aggregates exactly:
  `max_abs_cmr2_deviation` is a max over nodes (merge = max of maxes),
  `n_nodes_measured` is a count (merge = sum). No approximation.
- **frozen branch** is NOT sharded (39 nodes, ~35 min, and its I-F1..I-F6
  invariants are per-node within one axis). Exactly ONE array task
  carries `--with-frozen`; the merge asserts exactly one shard did.

THE TOLERANCE TRAP (why --zb-tol-km is mandatory here, not optional):
`build_zbw_grid_cache`'s `zb_tol_km=None` default is "half the smallest
spacing in zb_km_grid", and it also sets the Tb root-find's own
`solve_tol_km = 0.4 * zb_tol_km`. A shard containing only the coarse
5-20@1.0 km segment would derive 0.5 km, a shard of the fine 22-30@0.25
segment 0.125 km -- so shards would be solved to DIFFERENT precision and
the assembled cache would not be the grid any single-node run produces.
This script therefore requires an explicit --zb-tol-km and asserts it
equals the canonical production value, so every shard is byte-comparable
to the whole-grid build.

Usage (one shard):
  python plans/scripts/enceladus_zbw_shard_build.py \
      --shard-index 0 --n-shards 11 --out-dir /scratch/enc_shards
  # add --with-frozen on exactly one shard (the SLURM wrapper uses 0)

Verify the sharding scheme without building anything:
  python plans/scripts/enceladus_zbw_shard_build.py --n-shards 11 --plan
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from enceladus_zbw_production_build import (  # noqa: E402
    CONFIG_PATH, ZB_TOL_KM, FROZEN_ZB_TOL_KM, FROZEN_CUNCERTAINTY,
    FROZEN_MASS_TOL, FROZEN_RHO_CLOSURE_TOL_KGM3, FROZEN_MAX_ITER,
    _load_config, _verify_axes_match_config, _zb_axis, _w_axis,
    _frozen_axis)

MIN_SHARD_NODES = 2  # build_zbw_grid_cache requires a >=2-point zb axis


def _shard_slices(n_zb, n_shards):
    """Contiguous, gapless, non-overlapping zb index blocks.

    Contiguity matters: zb-adjacent nodes share EOS/spline warm state in
    PlanetProfile's solver path, and a contiguous block also makes each
    shard's own residual profile inspectable as a zb sweep. np.array_split
    gives near-equal blocks (sizes differ by at most 1).
    """
    idx = np.arange(n_zb)
    blocks = np.array_split(idx, n_shards)
    small = [len(b) for b in blocks if len(b) < MIN_SHARD_NODES]
    if small:
        raise SystemExit(
            f"--n-shards {n_shards} yields shard(s) with < "
            f"{MIN_SHARD_NODES} zb nodes ({small}); build_zbw_grid_cache "
            f"requires a >=2-point zb axis. Use --n-shards <= "
            f"{n_zb // MIN_SHARD_NODES}.")
    return blocks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--shard-index', type=int)
    ap.add_argument('--n-shards', type=int, required=True)
    ap.add_argument('--out-dir', default=None)
    ap.add_argument('--with-frozen', action='store_true',
                    help='Also build the 39-node frozen axis in THIS shard. '
                         'Exactly one shard in the array must set it.')
    ap.add_argument('--zb-tol-km', type=float, default=ZB_TOL_KM,
                    help='Explicit placement invariant (default: the '
                         'canonical production 0.125). Must not be left to '
                         'the builder default -- see module docstring.')
    ap.add_argument('--plan', action='store_true',
                    help='Print the shard decomposition and exit.')
    args = ap.parse_args()

    cfg = _load_config()
    zb_all, w_all, frozen_all = _verify_axes_match_config(cfg)
    blocks = _shard_slices(len(zb_all), args.n_shards)

    if args.plan:
        print(f'{args.n_shards} shards over {len(zb_all)} zb nodes '
              f'x {len(w_all)} w nodes = {len(zb_all)*len(w_all)} ocean '
              f'nodes (+{len(frozen_all)} frozen on shard 0):')
        for i, b in enumerate(blocks):
            print(f'  shard {i:2d}: zb[{b[0]:2d}:{b[-1]+1:2d}] '
                  f'= {len(b):2d} nodes, {zb_all[b[0]]:.2f}-'
                  f'{zb_all[b[-1]]:.2f} km, '
                  f'{len(b)*len(w_all)} ocean nodes')
        total = sum(len(b) for b in blocks)
        assert total == len(zb_all), 'shard decomposition lost nodes'
        print(f'  total {total} zb nodes -- partition verified complete')
        return

    if args.shard_index is None:
        raise SystemExit('--shard-index is required unless --plan is given')
    if not (0 <= args.shard_index < args.n_shards):
        raise SystemExit(f'--shard-index {args.shard_index} out of range '
                         f'[0, {args.n_shards})')
    if abs(args.zb_tol_km - ZB_TOL_KM) > 1e-12:
        raise SystemExit(
            f'--zb-tol-km {args.zb_tol_km} != the canonical production '
            f'{ZB_TOL_KM}. Shards MUST share one explicit tolerance or '
            f'they are not assemblable into the whole-grid product (see '
            f'module docstring, THE TOLERANCE TRAP).')

    from PlanetProfile.Inference.cache_builder import build_zbw_grid_cache

    block = blocks[args.shard_index]
    zb_shard = zb_all[block]
    out_dir = Path(args.out_dir or (REPO / 'PlanetProfile/Test/mcmc_results/'
                                           'Enceladus/Cassini_isostasy/'
                                           'shards'))
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f'shard_{args.shard_index:03d}_of_{args.n_shards:03d}.pkl'

    print(f'shard {args.shard_index}/{args.n_shards}: '
          f'zb[{block[0]}:{block[-1]+1}] = {len(zb_shard)} nodes '
          f'({zb_shard[0]:.4f}-{zb_shard[-1]:.4f} km) x {len(w_all)} w '
          f'= {len(zb_shard)*len(w_all)} ocean nodes'
          + (f' + {len(frozen_all)} frozen' if args.with_frozen else ''))

    t0 = time.time()
    cache = build_zbw_grid_cache(
        planet_template_module='PlanetProfile.Default.Enceladus.PPEnceladus',
        zb_km_grid=zb_shard.tolist(),
        wOcean_ppt_grid=w_all.tolist(),
        output_path=str(out_path),
        ocean_overrides=dict(cfg['ocean_overrides']),
        extrap_ocean=False,
        tb_placeholder_K=272.0,
        zb_tol_km=args.zb_tol_km,
        frozen_zb_km_grid=(frozen_all.tolist() if args.with_frozen else None),
        frozen_Cuncertainty=FROZEN_CUNCERTAINTY,
        frozen_mass_tol=FROZEN_MASS_TOL,
        frozen_rho_closure_tol_kgm3=FROZEN_RHO_CLOSURE_TOL_KGM3,
        frozen_zb_tol_km=FROZEN_ZB_TOL_KM,
        frozen_moi_nonconditioning_window=None,
        frozen_max_iter=FROZEN_MAX_ITER,
        config=cfg,
        progress=True,
    )
    wall_s = time.time() - t0

    # Sidecar so the merge step can verify the partition without unpickling
    # every (large) shard cache first.
    meta = {
        'shard_index': args.shard_index,
        'n_shards': args.n_shards,
        'zb_index_lo': int(block[0]),
        'zb_index_hi_exclusive': int(block[-1] + 1),
        'zb_km_first': float(zb_shard[0]),
        'zb_km_last': float(zb_shard[-1]),
        'n_zb': int(len(zb_shard)),
        'n_w': int(len(w_all)),
        'zb_tol_km': float(args.zb_tol_km),
        'with_frozen': bool(args.with_frozen),
        'n_zb_placement_rejected': int(cache['n_zb_placement_rejected']),
        'n_ocean_none': int(sum(1 for s in cache['structures'] if s is None)),
        'reachability_check_ran': (
            cache.get('ocean_reachability_restriction_table') is not None),
        'schema_version': cache.get('schema_version'),
        'wall_s': wall_s,
        'cache_path': str(out_path),
    }
    with open(out_path.with_suffix('.meta.json'), 'w') as f:
        json.dump(meta, f, indent=2)
    print(f'\nshard complete in {wall_s/3600:.2f}h -> {out_path}')
    print(json.dumps(meta, indent=2))


if __name__ == '__main__':
    main()
