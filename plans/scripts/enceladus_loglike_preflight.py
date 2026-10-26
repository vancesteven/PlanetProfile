#!/usr/bin/env python
"""Preflight for the Enceladus ocean6D rho_rock posterior run, on the REAL
cache (reviewer 2026-10-25, before the production pilot):

  1. one log_likelihood_fn call at a central theta -> finite, no exception;
  2. N prior draws through log_likelihood_fn -> accepted / rejected tally,
     reject reasons from runner.last_mass_closure['reason'], ms per call.

Exit 1 if the central call is not finite or any call raises.

  python plans/scripts/enceladus_loglike_preflight.py --config CFG.json \
      [--n 1000] [--seed 0] [--json-out FILE]
"""
import argparse
import collections
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402
from PlanetProfile.Inference.mcmc_runner import MCMCRunner  # noqa: E402
from PlanetProfile.Inference.forward_models import (  # noqa: E402
    apply_parameters, UnservableSampleError)


def classify_reject(runner, th):
    """Why log_likelihood_fn rejected th: re-run its stages in order."""
    td = runner._expand_theta(th)
    try:
        apply_parameters(td, runner.structure_data)
    except UnservableSampleError:
        return 'unservable_cache_cell (None corner)'
    if 'zb_km' in td and runner._ocean_zbw_struct_for_theta(td) is None:
        return 'zbw_cell_has_None_corner'
    runner.last_mass_closure = None
    if runner._derive_gravity_isostatic(td) is None:
        mc = runner.last_mass_closure or {}
        return f"gravity_none: {mc.get('reason') or 'no closure reason'}"
    if runner._derive_libration_deg(td) is None:
        return 'libration_none'
    return 'other (gravity+libration finite)'


def draw_prior(param_space, rng):
    x = []
    for nm, p in param_space.items():
        lo, hi = p['bounds']
        if p['prior_type'] == 'truncated_gaussian':
            m, s = p['mean'], p['std']
            x.append(stats.truncnorm.rvs((lo - m) / s, (hi - m) / s, loc=m,
                                         scale=s, random_state=rng))
        elif p['prior_type'] == 'uniform':
            x.append(rng.uniform(lo, hi))
        else:
            raise SystemExit(f'{nm}: prior_type {p["prior_type"]!r} not handled')
    return np.array(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--n', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--json-out')
    a = ap.parse_args()

    cfg = InferenceConfig.from_json(a.config)
    runner = MCMCRunner(cfg)
    ps = cfg.param_space
    names = list(ps)
    print('params:', names)

    centre = {'zb_km': 25.0, 'log10_wOcean_ppt': 0.5, 'rho_ice_kgm3': 925.0,
              'rho_rock_kgm3': 2400.0, 'compensation_C2': 0.75,
              'libration_sys_frac': 0.0}
    th0 = np.array([centre.get(n, np.mean(ps[n]['bounds'])) for n in names])
    ok = True
    try:
        ll0 = float(runner.log_likelihood_fn(th0))
    except Exception:
        traceback.print_exc()
        ll0, ok = float('nan'), False
    print(f'central theta {dict(zip(names, th0))} -> loglike {ll0}')
    if not (np.isfinite(ll0) and ll0 > -1e29):
        ok = False

    rng = np.random.default_rng(a.seed)
    reasons = collections.Counter()
    n_acc = n_exc = 0
    best = (-np.inf, None)
    rej = []
    t0 = time.perf_counter()
    for _ in range(a.n):
        th = draw_prior(ps, rng)
        runner.last_mass_closure = None
        try:
            ll = float(runner.log_likelihood_fn(th))
        except Exception as e:
            n_exc += 1
            reasons[f'EXCEPTION {type(e).__name__}: {e}'[:120]] += 1
            continue
        if np.isfinite(ll) and ll > -1e29:
            n_acc += 1
            if ll > best[0]:
                best = (ll, th)
        else:
            reasons[classify_reject(runner, th)] += 1
            rej.append(dict(zip(names, map(float, th))))
    dt = (time.perf_counter() - t0) / max(a.n, 1) * 1e3
    if n_exc:
        ok = False
    out = {'config': a.config, 'cache': cfg.structure_cache_path,
           'central_loglike': ll0, 'n': a.n, 'accepted': n_acc,
           'rejected': a.n - n_acc - n_exc, 'exceptions': n_exc,
           'reject_reasons': dict(reasons), 'ms_per_call': dt,
           'best_prior_draw': ({'loglike': best[0],
                                'theta': dict(zip(names, map(float, best[1])))}
                               if best[1] is not None else None),
           'rejected_zb_logw': [[round(r['zb_km'], 3),
                                 round(r['log10_wOcean_ppt'], 4)] for r in rej
                                if 'zb_km' in r],
           'pass': ok}
    print(json.dumps(out, indent=1))
    if a.json_out:
        Path(a.json_out).write_text(json.dumps(out, indent=1))
    print('PREFLIGHT', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
