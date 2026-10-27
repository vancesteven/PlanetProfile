#!/usr/bin/env python
"""Enceladus branch comparison by an evidence UPPER BOUND on the frozen branch
(reviewer pre-run ruling 2026-10-27, GO-with-bound).

Because the frozen prior integrates to 1, Z_frozen <= max_theta L(theta), so
log Z_frozen <= max log L over the declared frozen prior box. This script
finds that maximum on the PRODUCTION cache with the PRODUCTION likelihood
(the ocean6D_rhorock config reduced to the frozen parameters -- identical
observables, sigmas, sigma_model_add, correlations, libration factor), by
  * a dense grid (rho_rock x rho_ice x sys_frac), and
  * a multistart bounded Nelder-Mead,
records per-channel residuals at the maximum, a LIBRATION-ONLY bound
(support-model independent: holds under rigid or Airy gravity), the
out-of-grid prior fraction, and feeds the bound + run (a)'s ocean log Z
through branch_evidence.branch_bayes_factor (exercising the observable-
identity check). Report is one-sided by MODEL DEPENDENCE, not underflow.

  python plans/scripts/enceladus_frozen_branch_bound.py --cache CACHE \
      --ocean-result validation_reports/enceladus_isostasy/posterior_a/rep1_result.pkl \
      --out-dir validation_reports/enceladus_isostasy/branch_bound
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402
from PlanetProfile.Inference.mcmc_runner import MCMCRunner  # noqa: E402
from PlanetProfile.Inference.branch_evidence import branch_bayes_factor  # noqa: E402

CFG = REPO / 'PlanetProfile/Inference/configs/enceladus_cassini_isostasy_ocean6D_rhorock.json'
FROZEN = ('rho_ice_kgm3', 'rho_rock_kgm3', 'libration_sys_frac')
OBS = ('C20', 'C22', 'C30', 'libration_deg')


def build(cache):
    raw = json.loads(CFG.read_text())
    raw['structure_cache_path'] = str(cache)
    for k in [k for k in raw['param_space'] if k not in FROZEN]:
        raw['param_space'].pop(k)
    run = MCMCRunner(InferenceConfig(**raw))
    return raw, run


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cache', required=True)
    ap.add_argument('--ocean-result', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--grid', default='161,21,13',
                    help='n_rho_rock,n_rho_ice,n_sys')
    ap.add_argument('--starts', type=int, default=12)
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    raw, run = build(a.cache)
    names = list(run.param_names)
    ps = raw['param_space']
    lo = np.array([ps[n]['bounds'][0] for n in names])
    hi = np.array([ps[n]['bounds'][1] for n in names])
    ll = run.log_likelihood_fn

    # --- dense grid ---
    ng = dict(zip(('rho_rock_kgm3', 'rho_ice_kgm3', 'libration_sys_frac'),
                  map(int, a.grid.split(','))))
    axes = [np.linspace(lo[i], hi[i], ng[n]) for i, n in enumerate(names)]
    G = np.array(np.meshgrid(*axes, indexing='ij')).reshape(len(names), -1).T
    L = np.array([float(ll(t)) for t in G])
    rej = L <= -1e29
    k = int(np.argmax(np.where(rej, -np.inf, L)))
    grid_best = {'loglike': float(L[k]), 'theta': dict(zip(names, map(float, G[k])))}

    # --- multistart bounded Nelder-Mead ---
    rng = np.random.default_rng(75)
    starts = []
    for _ in range(a.starts):
        x0 = lo + (hi - lo) * rng.uniform(0.05, 0.95, len(names))
        f = (lambda x: 1e6 if (np.any(x < lo) or np.any(x > hi)) else -float(ll(x)))
        res = minimize(f, x0, method='Nelder-Mead',
                       options={'xatol': 1e-4, 'fatol': 1e-6, 'maxiter': 2000})
        starts.append({'loglike': float(-res.fun),
                       'theta': dict(zip(names, map(float, res.x)))})
    opt_best = max(starts, key=lambda s: s['loglike'])
    best = max((grid_best, opt_best), key=lambda s: s['loglike'])

    # --- per-channel residuals at the maximum ---
    th = np.array([best['theta'][n] for n in names])
    td = run._expand_theta(th)
    run.last_mass_closure = None
    iso = run._derive_gravity_isostatic(td)
    lib = run._derive_libration_deg(td)
    md = raw['metadata']
    pred = {'C20': float(iso[(2, 0)]), 'C22': float(iso[(2, 2)]),
            'C30': float(iso[(3, 0)]), 'libration_deg': float(lib)}
    s_inf = {o: float(np.hypot(raw['observables'][o][1],
                               md.get('sigma_model_add', {}).get(o, 0.0)))
             for o in OBS}
    resid = {o: (pred[o] - raw['observables'][o][0]) / s_inf[o] for o in OBS}

    # --- libration-only bound (support-model independent) ---
    # max over the box of the libration prediction; the libration term alone
    # then bounds log L from above (the gravity terms are <= 0).
    lib_vals = []
    for t in G:
        try:
            v = run._derive_libration_deg(run._expand_theta(t))
        except Exception:
            v = None
        lib_vals.append(np.nan if v is None else float(v))
    lib_vals = np.array(lib_vals)
    lib_max = float(np.nanmax(lib_vals))
    d = raw['observables']['libration_deg'][0] - lib_max
    lib_only_bound = -0.5 * (d / s_inf['libration_deg']) ** 2 if d > 0 else 0.0

    # --- out-of-grid prior fraction (frozen zb beyond the cache axis) ---
    from PlanetProfile.Gravity.isostasy import frozen_zb_from_mass
    spec = run.structure_data.get('frozen_build_spec', {}) if isinstance(
        run.structure_data, dict) else {}
    fz = np.asarray(run.structure_data.get('frozen_zb_km_grid', []), float) \
        if isinstance(run.structure_data, dict) else np.array([])
    M = float(spec.get('Mtot_kg', 1.08022e20)) if spec else 1.08022e20
    R = float(spec.get('R_m', 252.1e3)) if spec else 252.1e3
    rr = np.random.default_rng(1)
    rs = rr.uniform(lo[names.index('rho_rock_kgm3')], hi[names.index('rho_rock_kgm3')], 20000)
    ri = rr.uniform(lo[names.index('rho_ice_kgm3')], hi[names.index('rho_ice_kgm3')], 20000)
    zbs = np.array([frozen_zb_from_mass(M, R, x, y) for x, y in zip(rs, ri)]) / 1e3
    out_of_grid = (float(np.mean((zbs < fz.min() - 1e-6) | (zbs > fz.max() + 1e-6)))
                   if fz.size else None)

    # --- ocean evidence + Bayes factor via the utility ---
    oc = pickle.load(open(a.ocean_result, 'rb'))
    oc_cfg = json.loads((Path(a.ocean_result).parent / 'rep1_config.json').read_text())
    obs_f = {k: tuple(v) for k, v in raw['observables'].items()}
    obs_o = {k: tuple(v) for k, v in oc_cfg['observables'].items()}
    bf = branch_bayes_factor(best['loglike'], float(oc.metadata['log_Z']), obs_f, obs_o,
                             log_Z_frozen_err=None,
                             log_Z_ocean_err=float(oc.metadata['log_Z_err']))
    bf_lib = branch_bayes_factor(lib_only_bound, float(oc.metadata['log_Z']), obs_f, obs_o,
                                 log_Z_frozen_err=None,
                                 log_Z_ocean_err=float(oc.metadata['log_Z_err']))

    rep = {
        'method': 'log Z_frozen <= max log L over the declared frozen prior (prior integrates to 1)',
        'config': str(CFG), 'cache': a.cache, 'frozen_params': names,
        'frozen_prior_bounds': dict(zip(names, map(list, zip(lo.tolist(), hi.tolist())))),
        'grid': ng, 'n_grid': int(len(G)), 'grid_rejects': int(rej.sum()),
        'grid_best': grid_best, 'multistart': starts, 'multistart_best': opt_best,
        'max_loglike': best,
        'prediction_at_max': pred, 'inflated_sigma': s_inf,
        'residual_sigma_at_max': resid,
        'libration_max_over_box_deg': lib_max,
        'libration_only_bound_logL': lib_only_bound,
        'frozen_out_of_grid_prior_fraction': out_of_grid,
        'ocean': {'log_Z': float(oc.metadata['log_Z']),
                  'log_Z_err': float(oc.metadata['log_Z_err']),
                  'result': a.ocean_result},
        'branch_bayes_factor_full_bound': bf,
        'branch_bayes_factor_libration_only_bound': bf_lib,
        'labels': ('log_Z_frozen entries are UPPER BOUNDS, not estimates; log B '
                   'values are upper bounds on log Z_frozen - log Z_ocean'),
    }
    (out / 'frozen_branch_bound.json').write_text(json.dumps(rep, indent=1, default=float))
    print(json.dumps({k: rep[k] for k in ('max_loglike', 'residual_sigma_at_max',
                                          'libration_max_over_box_deg',
                                          'libration_only_bound_logL',
                                          'frozen_out_of_grid_prior_fraction',
                                          'grid_rejects')}, indent=1, default=float))
    print('log B <= %.2f (full bound); <= %.2f (libration-only bound)' % (
        bf['log_bayes_factor'], bf_lib['log_bayes_factor']))
    print('utility one_sided flag:', bf['one_sided'], '| notes:', bf['notes'][:2])


if __name__ == '__main__':
    main()
