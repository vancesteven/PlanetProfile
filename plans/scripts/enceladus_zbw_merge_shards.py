#!/usr/bin/env python
"""Merge Enceladus (zb, w) cache shards into the production cache + report.

Companion to `enceladus_zbw_shard_build.py`. Assembles the sharded ocean
grid, verifies the partition reconstructs the canonical axis EXACTLY, and
emits the six r7 §0.29 acceptance-gate figures.

DESIGN CONSTRAINT HONORED: this script derives the merged cache purely
from what the production builder already recorded in each shard. It does
NOT re-run, re-implement, or relax any r7 gate, and it imports no private
`cache_builder` helper:

- C2 reachability: already enforced per-node inside each shard's own
  `build_zbw_grid_cache` call (the predicate has no cross-node term --
  see the shard script's docstring). This script ASSERTS every shard ran
  it (`ocean_reachability_restriction_table is not None`) and refuses to
  merge if any shard did not. It does not re-derive the check, because
  duplicating a gate's logic is how the two implementations drift.
- zb_placement: per-node inside each shard; reject counts are summed.
- ocean_moi_window: the only aggregate, merged EXACTLY --
  `max_abs_cmr2_deviation` = max of the shard maxes, `n_nodes_measured` =
  sum of shard counts, `margin` = half-width - merged max. The
  half-width/source/Cmeasured fields must be identical across shards
  (asserted) since every shard resolved them from the same config.
- frozen branch: taken verbatim from the single shard built with
  `--with-frozen`; asserted to be exactly one.

Usage:
  python plans/scripts/enceladus_zbw_merge_shards.py \
      --shard-dir /scratch/enc_shards \
      --out PlanetProfile/Test/mcmc_results/Enceladus/Cassini_isostasy/\
enceladus_zbw_production_cache.pkl
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from enceladus_zbw_production_build import (  # noqa: E402
    ZB_TOL_KM, _load_config, _verify_axes_match_config)

REPORT_PATH = (REPO / 'validation_reports/enceladus_isostasy/'
                      'r7_production_build_report.json')

# Fields that must agree bit-for-bit across every shard: they encode the
# build contract (composition, tolerance, declared restriction), so a
# disagreement means the shards are not assemblable.
_IDENTICAL_FIELDS = ('ocean_comp', 'zb_tol_km',
                     'ocean_reachability_restriction_table')


def _load_shards(shard_dir):
    metas = sorted(Path(shard_dir).glob('shard_*_of_*.meta.json'))
    if not metas:
        raise SystemExit(f'no shard_*.meta.json found in {shard_dir}')
    shards = []
    for m in metas:
        with open(m) as f:
            meta = json.load(f)
        cache_path = Path(meta['cache_path'])
        if not cache_path.exists():          # moved between scratch and home
            cache_path = m.parent / (m.name.replace('.meta.json', '.pkl'))
        if not cache_path.exists():
            raise SystemExit(f'shard cache missing for {m}: {cache_path}')
        with open(cache_path, 'rb') as f:
            shards.append((meta, pickle.load(f)))
    shards.sort(key=lambda mc: mc[0]['shard_index'])
    return shards


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--shard-dir', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--report', default=str(REPORT_PATH))
    args = ap.parse_args()

    cfg = _load_config()
    zb_all, w_all, frozen_all = _verify_axes_match_config(cfg)
    shards = _load_shards(args.shard_dir)

    n_shards_declared = shards[0][0]['n_shards']
    problems = []

    # --- partition completeness: exact, index-level -------------------
    if len(shards) != n_shards_declared:
        problems.append(f'found {len(shards)} shards but each declares '
                        f'n_shards={n_shards_declared} -- incomplete array '
                        f'(a failed/pending task would silently truncate '
                        f'the grid)')
    covered = []
    for meta, _ in shards:
        covered.extend(range(meta['zb_index_lo'],
                             meta['zb_index_hi_exclusive']))
    if sorted(covered) != list(range(len(zb_all))):
        dupes = {i for i in covered if covered.count(i) > 1}
        missing = set(range(len(zb_all))) - set(covered)
        problems.append(f'zb index coverage is not a clean partition of '
                        f'[0,{len(zb_all)}) -- missing={sorted(missing)}, '
                        f'duplicated={sorted(dupes)}')

    # --- contract fields identical across shards ----------------------
    ref_meta, ref_cache = shards[0]
    for field in _IDENTICAL_FIELDS:
        ref_val = ref_cache.get(field)
        for meta, cache in shards[1:]:
            if repr(cache.get(field)) != repr(ref_val):
                problems.append(
                    f'shard {meta["shard_index"]} disagrees with shard '
                    f'{ref_meta["shard_index"]} on {field!r} -- shards were '
                    f'not built under one contract')
                break
    if abs(float(ref_cache.get('zb_tol_km', np.nan)) - ZB_TOL_KM) > 1e-12:
        problems.append(f'shard zb_tol_km {ref_cache.get("zb_tol_km")} != '
                        f'canonical {ZB_TOL_KM}')
    for meta, cache in shards:
        if not np.array_equal(np.asarray(cache['wOcean_ppt_grid']),
                              np.asarray(ref_cache['wOcean_ppt_grid'])):
            problems.append(f'shard {meta["shard_index"]} w-axis differs '
                            f'from shard {ref_meta["shard_index"]} '
                            f'(bit-level) -- not assemblable')

    # --- C2 actually ran on every shard (never assumed) ---------------
    not_checked = [m['shard_index'] for m, c in shards
                   if c.get('ocean_reachability_restriction_table') is None]
    if not_checked:
        problems.append(
            f'shard(s) {not_checked} carry no '
            f'ocean_reachability_restriction_table, i.e. the r7 C2 '
            f'two-sided check did NOT run on their nodes (built without '
            f'config=). Those nodes are unvalidated; refusing to merge.')

    # --- frozen axis: exactly one owner -------------------------------
    frozen_owners = [m['shard_index'] for m, c in shards
                     if c.get('frozen_branch_supported')]
    if len(frozen_owners) != 1:
        problems.append(f'expected exactly ONE shard carrying the frozen '
                        f'axis, found {frozen_owners} -- rerun the array '
                        f'with --with-frozen on exactly one task')

    if problems:
        raise SystemExit('MERGE REFUSED -- ' + '; '.join(problems))

    # --- assemble ------------------------------------------------------
    zb_merged, structures, none_reasons = [], [], []
    n_zb_reject = 0
    for meta, cache in shards:
        zb_merged.extend(np.asarray(cache['zb_km_grid']).tolist())
        structures.extend(cache['structures'])
        none_reasons.extend(cache.get('ocean_none_reasons')
                            or [None] * len(cache['structures']))
        n_zb_reject += int(cache['n_zb_placement_rejected'])
    zb_merged = np.asarray(zb_merged, dtype=float)

    if not np.array_equal(zb_merged, zb_all):
        raise SystemExit(
            'MERGE REFUSED -- reassembled zb axis is not bit-identical to '
            'the canonical axis from the config. First difference at index '
            f'{int(np.argmax(zb_merged != zb_all))}.')
    n_w = len(w_all)
    if len(structures) != len(zb_all) * n_w:
        raise SystemExit(
            f'MERGE REFUSED -- {len(structures)} structures for a '
            f'{len(zb_all)}x{n_w}={len(zb_all)*n_w} grid')

    # --- exact aggregation of the one aggregate ------------------------
    windows = [c['ocean_moi_window'] for _, c in shards
               if c.get('ocean_moi_window')]
    for key in ('Cmeasured', 'CuncertaintyUpper', 'CuncertaintyLower',
                'source'):
        vals = {repr(wd.get(key)) for wd in windows}
        if len(vals) > 1:
            raise SystemExit(f'MERGE REFUSED -- shards disagree on '
                             f'ocean_moi_window[{key!r}]: {vals}')
    devs = [wd['max_abs_cmr2_deviation'] for wd in windows
            if wd.get('max_abs_cmr2_deviation') is not None
            and np.isfinite(wd['max_abs_cmr2_deviation'])]
    merged_max_dev = float(max(devs)) if devs else float('nan')
    half = windows[0].get('CuncertaintyUpper')
    merged_window = {
        **windows[0],
        'n_nodes_measured': int(sum(wd.get('n_nodes_measured', 0)
                                    for wd in windows)),
        'max_abs_cmr2_deviation': merged_max_dev,
        'margin': (float(half) - merged_max_dev
                   if half is not None and np.isfinite(merged_max_dev)
                   else float('nan')),
        'merged_from_shards': len(windows),
        'merge_note': (
            'max_abs_cmr2_deviation is the max over shard maxes and '
            'n_nodes_measured the sum over shard counts -- both exact for '
            'a partition of the node set, not an approximation. Every '
            'other field was asserted identical across shards.'),
    }

    frozen_shard = next(c for m, c in shards
                        if c.get('frozen_branch_supported'))
    cache_out = {
        'zb_km_grid': zb_all,
        'wOcean_ppt_grid': np.asarray(w_all),
        'structures': structures,
        'ocean_comp': ref_cache['ocean_comp'],
        'schema_version': frozen_shard['schema_version'],
        'zb_tol_km': float(ref_cache['zb_tol_km']),
        'n_zb_placement_rejected': int(n_zb_reject),
        'frozen_branch_supported': True,
        'ocean_moi_window': merged_window,
        'ocean_reachability_restriction_table':
            ref_cache['ocean_reachability_restriction_table'],
        'ocean_none_reasons': none_reasons,
        'frozen_zb_km_grid': frozen_shard['frozen_zb_km_grid'],
        'frozen_structures': frozen_shard['frozen_structures'],
        'frozen_build_spec': frozen_shard['frozen_build_spec'],
        'n_frozen_if2_rejected': int(frozen_shard['n_frozen_if2_rejected']),
        'built_by': {
            'mode': 'slurm_sharded',
            'n_shards': n_shards_declared,
            'shard_wall_s': {str(m['shard_index']): m['wall_s']
                             for m, _ in shards},
            'gate_provenance': (
                'C2 reachability + zb_placement were enforced per-node '
                'inside each shard by the production builder; this merge '
                'asserted every shard ran them and that the shards form an '
                'exact partition of the canonical axis. No gate was '
                're-implemented or relaxed here.'),
        },
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix('.pkl.tmp')
    with open(tmp, 'wb') as f:
        pickle.dump(cache_out, f, protocol=pickle.HIGHEST_PROTOCOL)
    tmp.replace(out_path)

    # --- r7 §0.29 acceptance gates ------------------------------------
    residuals = [abs(s['zb_residual_km']) for s in structures
                 if s is not None and 'zb_residual_km' in s]
    max_residual = max(residuals) if residuals else None
    n_none = sum(1 for s in structures if s is None)
    n_total = len(structures)
    fspec = cache_out['frozen_build_spec']
    report = {
        'ruling': 'validation_reports/enceladus_isostasy/r7_ADJUDICATION.md '
                  '(+ §0.30 endpoint ruling)',
        'build_mode': f'SLURM sharded ({n_shards_declared} array tasks)',
        'output_path': str(out_path),
        'schema_version': cache_out['schema_version'],
        'n_zb': len(zb_all), 'n_w': n_w, 'n_frozen': len(frozen_all),
        'acceptance_gates': {
            '1_zb_placement_rejects_zero': {
                'value': int(n_zb_reject), 'pass': n_zb_reject == 0},
            '2_max_residual_under_zb_tol': {
                'value_km': max_residual, 'threshold_km': ZB_TOL_KM,
                'expectation_km': 0.035,
                'pass': (max_residual is not None
                         and max_residual < ZB_TOL_KM)},
            '3_reachability_none_set': {
                'n_ocean_none': int(n_none), 'n_ocean_total': int(n_total),
                'f_excluded': (n_none / n_total) if n_total else None,
                'pass': True,
                'note': 'C2 two-sided enforcement passed inside every '
                        'shard; a violation would have aborted that shard.'},
            '4_moi_window_floor_remeasured': merged_window,
            '5_frozen_invariants_not_fired': {
                'n_frozen_nodes': fspec.get('n_nodes'),
                'n_built': fspec.get('n_built'),
                'n_rejected': fspec.get('n_rejected'),
                'n_frozen_if2_rejected':
                    cache_out['n_frozen_if2_rejected'],
                'pass': (fspec.get('n_rejected') == 0
                         and cache_out['n_frozen_if2_rejected'] == 0)},
            '6_f_excluded_recorded': {
                'f_excluded': (n_none / n_total) if n_total else None,
                'pass': True,
                'note': 'Required by C3 for the +log(1/(1-f_excluded)) '
                        'evidence normalization before any log B.'},
        },
        'total_shard_wall_h': sum(m['wall_s'] for m, _ in shards) / 3600.0,
    }
    Path(args.report).parent.mkdir(parents=True, exist_ok=True)
    with open(args.report, 'w') as f:
        json.dump(report, f, indent=2, default=str)

    gates = report['acceptance_gates']
    all_pass = all(g.get('pass') for g in gates.values() if 'pass' in g)
    print(f'\nMerged -> {out_path}')
    print(json.dumps(report, indent=2, default=str))
    print(f'\nACCEPTANCE GATES: '
          f'{"ALL PASS" if all_pass else "AT LEAST ONE FAILED"} '
          f'(interpret, never tune -- report to Machine A either way)')


if __name__ == '__main__':
    main()
