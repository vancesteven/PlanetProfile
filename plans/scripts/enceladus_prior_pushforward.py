#!/usr/bin/env python
"""Run (h): prior-only pushforward for the Enceladus ocean6D rho_rock campaign,
by DIRECT prior Monte Carlo (no sampler).

Reviewer Fisher-checkpoint condition 3 / posterior_a review:
  * prior pushforward of derived C/MR^2 -- "do not present C/MR^2 as a
    measurement until (h) shows the prior pushforward is much wider";
  * the R_s < R_o support edge (declared physical edge of the ocean branch);
  * the None-corner hole (zb <= ~13 km) must be TALLIED SEPARATELY from
    physical rejects, not mixed into the effective prior.

Per prior draw, through the SAME runner methods the likelihood and the
posterior derived-quantity loop use (no re-implemented physics):
  apply_parameters / _ocean_zbw_struct_for_theta  -> cache hole?
  _derive_gravity_isostatic + last_mass_closure   -> C20/C22/C30, C/MR^2,
                                                      R_s, ocean thickness,
                                                      physical reject reason
  _derive_libration_deg                           -> libration prediction
                                                      (incl. model budget)
Each draw is classified exactly once:
  'ok'                 all observables finite
  'hole'               cache cell unservable / has a None corner
  'physical:<reason>'  mass-closure / support reject (e.g. R_s_not_in_0_Ro)
  'libration_none'     gravity finite but libration None
  'exception:<type>'   anything raising (counted, never silently dropped)

Outputs a JSON with the tallies, the hole and physical-reject prior-volume
fractions, the effective joint prior's marginals vs the nominal prior, and
pushforward quantiles of every derived quantity; plus an .npz of all draws.

  python plans/scripts/enceladus_prior_pushforward.py --config CFG.json \
      [--n 20000] [--seed 0] [--out-dir DIR]
"""
import argparse
import collections
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402
from PlanetProfile.Inference.mcmc_runner import MCMCRunner  # noqa: E402
from PlanetProfile.Inference.forward_models import (  # noqa: E402
    apply_parameters, UnservableSampleError)
from enceladus_loglike_preflight import draw_prior  # noqa: E402

DERIVED = ('C20', 'C22', 'C30', 'libration_deg', 'cmr2_derived', 'R_s_km',
           'ocean_thickness_km', 'rho_s')
Q = (0.005, 0.025, 0.16, 0.5, 0.84, 0.975, 0.995)


def classify_and_derive(runner, th):
    """Return (category, derived dict). Mirrors log_likelihood_fn's stage
    order (cache -> isostatic gravity + mass closure -> libration)."""
    td = runner._expand_theta(th)
    try:
        apply_parameters(td, runner.structure_data)
    except UnservableSampleError:
        return 'hole', {}
    if 'zb_km' in td and runner._ocean_zbw_struct_for_theta(td) is None:
        return 'hole', {}
    runner.last_mass_closure = None
    iso = runner._derive_gravity_isostatic(td)
    mc = runner.last_mass_closure or {}
    if iso is None:
        return f"physical:{mc.get('reason') or 'gravity_none_no_reason'}", {}
    lib = runner._derive_libration_deg(td)
    d = {'C20': float(iso.get((2, 0), np.nan)),
         'C22': float(iso.get((2, 2), np.nan)),
         'C30': float(iso.get((3, 0), np.nan)),
         'libration_deg': float(lib) if lib is not None else np.nan,
         'cmr2_derived': float(mc.get('cmr2_derived', np.nan)),
         'R_s_km': float(mc.get('R_s_m', np.nan)) / 1e3,
         'ocean_thickness_km': float(mc.get('ocean_thickness_m', np.nan)) / 1e3,
         'rho_s': float(mc.get('rho_s', np.nan))}
    if lib is None:
        return 'libration_none', d
    return 'ok', d


def quantiles(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return None
    return {'n': int(v.size), 'q': dict(zip([str(q) for q in Q],
                                            np.quantile(v, Q).tolist())),
            'mean': float(v.mean()), 'sd': float(v.std(ddof=1))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--n', type=int, default=20000)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    cfg = InferenceConfig.from_json(a.config)
    runner = MCMCRunner(cfg)
    ps = cfg.param_space
    names = list(ps)
    assert names == list(runner.param_names), (names, runner.param_names)
    obs = {k: (float(v[0]), float(v[1])) for k, v in cfg.observables.items()}

    rng = np.random.default_rng(a.seed)
    TH = np.empty((a.n, len(names)))
    cat = []
    DV = {k: np.full(a.n, np.nan) for k in DERIVED}
    t0 = time.perf_counter()
    for i in range(a.n):
        th = draw_prior(ps, rng)
        TH[i] = th
        try:
            c, d = classify_and_derive(runner, th)
        except Exception as e:  # counted, never dropped silently
            c, d = f'exception:{type(e).__name__}', {}
        cat.append(c)
        for k, v in d.items():
            DV[k][i] = v
        if (i + 1) % 2000 == 0:
            el = time.perf_counter() - t0
            print(f'{i + 1}/{a.n}  {el:.0f}s  {collections.Counter(cat).most_common(4)}',
                  flush=True)
    wall = time.perf_counter() - t0
    cat = np.array(cat)
    tally = dict(collections.Counter(cat.tolist()))
    ok = cat == 'ok'
    hole = cat == 'hole'
    phys = np.char.startswith(cat.astype(str), 'physical:')
    n = a.n

    # Effective joint prior = nominal prior restricted to the declared
    # physical support (hole excluded from BOTH numerator and denominator:
    # it is a cache artefact, not physics).
    servable = ~hole
    rep = {
        'config': a.config, 'cache': cfg.structure_cache_path,
        'n': n, 'seed': a.seed, 'wall_s': wall, 'ms_per_draw': wall / n * 1e3,
        'tally': tally,
        'fractions': {
            'ok': float(ok.mean()),
            'hole_cache_artefact': float(hole.mean()),
            'physical_reject_of_all': float(phys.mean()),
            'physical_reject_of_servable': (float(phys[servable].mean())
                                            if servable.any() else None),
            'exceptions': float(np.char.startswith(cat.astype(str),
                                                   'exception:').mean()),
        },
        'hole_zb_range_km': ([float(TH[hole, 0].min()), float(TH[hole, 0].max())]
                             if hole.any() else None),
        'physical_reject_theta_ranges': {
            nm: [float(TH[phys, j].min()), float(TH[phys, j].max())]
            for j, nm in enumerate(names)} if phys.any() else None,
        'nominal_prior_marginals': {nm: quantiles(TH[:, j])
                                    for j, nm in enumerate(names)},
        'effective_prior_marginals_ok': {nm: quantiles(TH[ok, j])
                                         for j, nm in enumerate(names)},
        'pushforward_ok': {k: quantiles(DV[k][ok]) for k in DERIVED},
        'observed': obs,
    }
    # Pushforward width vs observation sigma, per observable.
    rep['pushforward_sd_over_sigma_obs'] = {
        k: (rep['pushforward_ok'][k]['sd'] / s) if rep['pushforward_ok'].get(k)
        else None for k, (_, s) in obs.items()}
    (out / 'prior_pushforward.json').write_text(json.dumps(rep, indent=1))
    np.savez_compressed(out / 'prior_pushforward_draws.npz', theta=TH,
                        param_names=np.array(names), category=cat,
                        **{k: DV[k] for k in DERIVED})
    print(json.dumps({k: rep[k] for k in ('tally', 'fractions', 'hole_zb_range_km',
                                          'pushforward_sd_over_sigma_obs')},
                     indent=1))
    print('cmr2_derived prior pushforward:', rep['pushforward_ok']['cmr2_derived'])
    print('wrote', out)


if __name__ == '__main__':
    main()
