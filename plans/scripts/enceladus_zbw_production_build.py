#!/usr/bin/env python
"""Production Enceladus (zb, w) + frozen-branch structure cache build.

Authorized by validation_reports/enceladus_isostasy/r7_ADJUDICATION.md
("r7 — RATIFIED WITH CONDITIONS. BUILD AUTHORIZED after C1-C3.") and
plans/MACHINE-B-HANDOFF.md §0.29 ("the moment the C1-C3 commit lands, B
RUNS THE BUILD per the r7 instructions"). C1-C3 landed in commit
2f9c6774 (ocean_reachability_restriction table + two-sided builder
enforcement + evidence-protocol normalization record) and are covered
by 44 passing tests (enceladus_zb_placement_test.py,
enceladus_ocean_zbw_lookup_test.py, enceladus_moi_window_test.py) run
before this build was launched.

Per r7's literal build instructions (validation_reports/
enceladus_isostasy/r7_ADJUDICATION.md, "MACHINE B BUILD INSTRUCTIONS"):
- zb segmented axis: 5-20@1.0 / 20-22@0.5 / 22-30@0.25 / 30-42@0.5 /
  42-45@0.25 = 87 nodes (matches the config's own
  metadata.structure_cache_spec.zb_km_grid.n_total == 87, verified
  below before the build starts).
- w = 10**linspace(-1, 2, 40) (matches
  metadata.structure_cache_spec.log10_w_grid {lo:-1, hi:2, n:40}).
- ocean_overrides = {comp: Seawater, deltaT: 0.002} (config's own
  top-level ocean_overrides, unchanged).
- zb_tol_km=0.125 EXPLICIT (r7's stated invariant; the builder derives
  its root-find solve_tol_km = 0.4 * zb_tol_km = 0.05 internally).
- config=<the candidate config dict>, NO bulk_overrides argument — the
  builder resolves Cuncertainty=0.08 from the config's own
  metadata['bulk_overrides'] (cache_builder.bulk_overrides_from_config);
  passing bulk_overrides explicitly here would silently win over the
  config's declared MoI window per the documented precedence rule.
- frozen axis: arange(46.5, 65.8001, 0.5) = 39 nodes,
  frozen_Cuncertainty=0.015, mass_tol=1e-6, rho_closure_tol=12.0,
  frozen_zb_tol_km=0.25, moi_nonconditioning_window=None, max_iter=10.
- extrap_ocean=False, tb_placeholder_K=272.0. Schema v3.2-zbw-joint
  (automatic from build_zbw_grid_cache when frozen_zb_km_grid is
  supplied).

ACCEPTANCE (r7): 0 placement rejects, max residual < 0.125 km (expect
<=0.035); None set == C1's predicted reachability-restriction table
within one cell (any node outside it halts the build — enforced INSIDE
build_zbw_grid_cache itself per C2, not re-checked here); mass invariant
on every frozen node; frozen I-F1..I-F4 + I-F6 must not fire (also
enforced inside the builder — a violating node is a hard reject to
None, never a silent pass); f_excluded recorded on the returned cache.

Usage:
  python plans/scripts/enceladus_zbw_production_build.py --dry-run
  python plans/scripts/enceladus_zbw_production_build.py
"""
import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

CONFIG_PATH = str(REPO / 'PlanetProfile/Inference/configs/'
                        'enceladus_cassini_isostasy_7D.json')
OUT_DIR = REPO / 'PlanetProfile/Test/mcmc_results/Enceladus/Cassini_isostasy'
OUT_PATH = OUT_DIR / 'enceladus_zbw_production_cache.pkl'
REPORT_PATH = (REPO / 'validation_reports/enceladus_isostasy/'
                      'r7_production_build_report.json')

ZB_TOL_KM = 0.125
FROZEN_ZB_TOL_KM = 0.25
FROZEN_CUNCERTAINTY = 0.015
FROZEN_MASS_TOL = 1e-6
FROZEN_RHO_CLOSURE_TOL_KGM3 = 12.0
FROZEN_MAX_ITER = 10


def _zb_axis():
    segs = [(5.0, 1.0, 15), (20.0, 0.5, 4), (22.0, 0.25, 32),
            (30.0, 0.5, 24), (42.0, 0.25, 12)]
    zb = []
    for lo, step, n in segs:
        zb.extend(lo + step * i for i in range(n))
    return np.asarray(zb, dtype=float)


def _w_axis():
    return 10 ** np.linspace(-1.0, 2.0, 40)


def _frozen_axis():
    return np.arange(46.5, 65.8001, 0.5)


def _load_config():
    with open(CONFIG_PATH) as f:
        return json.load(f)


def _verify_axes_match_config(cfg):
    spec = cfg['metadata']['structure_cache_spec']
    zb = _zb_axis()
    w = _w_axis()
    frozen = _frozen_axis()
    n_zb_declared = spec['zb_km_grid']['n_total']
    n_w_declared = spec['log10_w_grid']['n']
    n_frozen_declared = spec['frozen_zb_km_grid']['n']
    problems = []
    if len(zb) != n_zb_declared:
        problems.append(f'zb axis len {len(zb)} != declared {n_zb_declared}')
    if len(w) != n_w_declared:
        problems.append(f'w axis len {len(w)} != declared {n_w_declared}')
    if len(frozen) != n_frozen_declared:
        problems.append(f'frozen axis len {len(frozen)} != declared '
                        f'{n_frozen_declared}')
    if abs(zb[0] - spec['zb_km_grid']['lo']) > 1e-9 or \
            abs(zb[-1] - spec['zb_km_grid']['hi']) > 1e-9:
        problems.append(f'zb axis endpoints {zb[0]}/{zb[-1]} != declared '
                        f"{spec['zb_km_grid']['lo']}/"
                        f"{spec['zb_km_grid']['hi']}")
    if problems:
        raise SystemExit('Axis construction does NOT match the config\'s '
                         'declared structure_cache_spec -- refusing to '
                         'build. ' + '; '.join(problems))
    print(f'Axes verified against config: zb={len(zb)} nodes '
          f'[{zb[0]}, {zb[-1]}] km, w={len(w)} nodes '
          f'[{w[0]:.3f}, {w[-1]:.3f}] ppt, frozen={len(frozen)} nodes '
          f'[{frozen[0]}, {frozen[-1]}] km.')
    return zb, w, frozen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true',
                    help='Verify axes + config wiring, print the exact '
                         'call that would be made, and exit without '
                         'building.')
    ap.add_argument('--out', default=str(OUT_PATH))
    args = ap.parse_args()

    from PlanetProfile.Inference.cache_builder import build_zbw_grid_cache

    cfg = _load_config()
    zb_arr, w_arr, frozen_arr = _verify_axes_match_config(cfg)
    ocean_overrides = dict(cfg['ocean_overrides'])

    call_kwargs = dict(
        planet_template_module='PlanetProfile.Default.Enceladus.PPEnceladus',
        zb_km_grid=zb_arr.tolist(),
        wOcean_ppt_grid=w_arr.tolist(),
        output_path=args.out,
        ocean_overrides=ocean_overrides,
        extrap_ocean=False,
        tb_placeholder_K=272.0,
        zb_tol_km=ZB_TOL_KM,
        frozen_zb_km_grid=frozen_arr.tolist(),
        frozen_Cuncertainty=FROZEN_CUNCERTAINTY,
        frozen_mass_tol=FROZEN_MASS_TOL,
        frozen_rho_closure_tol_kgm3=FROZEN_RHO_CLOSURE_TOL_KGM3,
        frozen_zb_tol_km=FROZEN_ZB_TOL_KM,
        frozen_moi_nonconditioning_window=None,
        frozen_max_iter=FROZEN_MAX_ITER,
        config=cfg,
        progress=True,
    )

    if args.dry_run:
        printable = {k: (v if not isinstance(v, list) else
                         f'<list len={len(v)}>')
                    for k, v in call_kwargs.items()}
        print(json.dumps(printable, indent=2, default=str))
        print('\nDRY RUN: no build performed.')
        return

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    cache = build_zbw_grid_cache(**call_kwargs)
    wall_s = time.time() - t0

    n_zb_reject = cache.get('n_zb_placement_rejected')
    residuals = [s['zb_residual_km'] for s in cache.get('structures', [])
                if s is not None and 'zb_residual_km' in s]
    max_residual = max((abs(r) for r in residuals), default=None)
    n_none = sum(1 for s in cache.get('structures', []) if s is None)
    n_total = len(cache.get('structures', []))

    report = {
        'ruling': 'validation_reports/enceladus_isostasy/r7_ADJUDICATION.md',
        'commit_c1_c3': '2f9c6774',
        'wall_s': wall_s,
        'output_path': args.out,
        'n_zb': len(zb_arr), 'n_w': len(w_arr), 'n_frozen': len(frozen_arr),
        'n_ocean_total': n_total,
        'n_ocean_none': n_none,
        'n_ocean_built': n_total - n_none,
        'n_zb_placement_rejected': n_zb_reject,
        'max_abs_zb_residual_km': max_residual,
        'acceptance_zb_reject_zero': (n_zb_reject == 0),
        'acceptance_residual_under_0125': (
            max_residual is not None and max_residual < 0.125),
        'ocean_moi_window': cache.get('ocean_moi_window'),
        'ocean_reachability_restriction_table_present':
            cache.get('ocean_reachability_restriction_table') is not None,
        'f_excluded': (n_none / n_total) if n_total else None,
        'schema_version': cache.get('schema_version'),
        'frozen_n_structures': len(cache.get('frozen_structures', []))
            if 'frozen_structures' in cache else None,
    }
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(REPORT_PATH, 'w') as f:
        json.dump(report, f, indent=2, default=str)
    print(f'\nWrote {REPORT_PATH}')
    print(json.dumps(report, indent=2, default=str))


if __name__ == '__main__':
    main()
