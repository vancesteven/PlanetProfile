#!/usr/bin/env python
"""G0 data-integrity gate + merge for the Enceladus ocean-branch SBI training
set (DESIGN_PREREGISTRATION.md). Runs BEFORE training.

G0 (preregistered):
  * servable fraction of the training set in [0.887, 0.900]
  * empirical noise sd within 1% of sigma_infl, corr(C20,C22) in [0.46, 0.48]
    -- measured on (x - x0), i.e. the injected noise itself, pooled over all
    rows (x0 is the noise-free regeneration through theta_override)
  * 1,000 random rows' x0 match the likelihood path's own predictions
    (_derive_gravity_isostatic / _derive_libration_deg) to rtol 1e-10
  * shard k used noise_seed = 7575 + k (and data seed 75 + k)
Writes the merged training npz (theta, x, plus provenance) only if G0 passes.

  python plans/scripts/enceladus_sbi_g0.py --dir ~/scratch/enc_sbi --cache CACHE
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402
from PlanetProfile.Inference.mcmc_runner import MCMCRunner  # noqa: E402

CFG = REPO / 'PlanetProfile/Inference/configs/enceladus_cassini_isostasy_ocean6D_rhorock.json'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', required=True)
    ap.add_argument('--cache', required=True)
    ap.add_argument('--n-check', type=int, default=1000)
    a = ap.parse_args()
    D = Path(a.dir)
    shards = sorted(D.glob('train_shard_*.npz'))
    assert len(shards) == 11, f'expected 11 training shards, found {len(shards)}'

    TH, X, X0, metas = [], [], [], []
    for f in shards:
        d = np.load(f)
        m = json.loads(str(d['meta']))
        metas.append(m)
        TH.append(d['theta']); X.append(d['x']); X0.append(d['x0'])
    TH, X, X0 = map(np.concatenate, (TH, X, X0))
    res = {'n_kept': int(len(TH)),
           'n_request': int(sum(m['n_request'] for m in metas))}

    seeds_ok = all(m['data_seed'] == 75 + m['k'] and m['noise_seed'] == 7575 + m['k']
                   for m in metas) and sorted(m['k'] for m in metas) == list(range(11))
    hashes = {m['config_hash'] for m in metas}
    caches = {m['cache_sha256'] for m in metas}
    res['seeds_ok'] = bool(seeds_ok)
    res['single_config_hash'] = sorted(hashes)
    res['single_cache_sha'] = sorted(caches)
    servable = res['n_kept'] / res['n_request']
    res['servable_fraction'] = servable
    res['servable_ok'] = bool(0.887 <= servable <= 0.900)

    raw = json.loads(CFG.read_text())
    obs = list(raw['observables'].keys())
    md = raw['metadata']
    s_inf = np.array([np.hypot(raw['observables'][k][1], md['sigma_model_add'].get(k, 0.0))
                      for k in obs])
    E = X - X0
    sd = E.std(axis=0, ddof=1)
    rho = float(np.corrcoef(E[:, obs.index('C20')], E[:, obs.index('C22')])[0, 1])
    res['noise_sd_over_sigma_infl'] = dict(zip(obs, (sd / s_inf).tolist()))
    res['noise_corr_C20_C22'] = rho
    res['noise_ok'] = bool(np.all(np.abs(sd / s_inf - 1) <= 0.01) and 0.46 <= rho <= 0.48)

    raw['structure_cache_path'] = a.cache
    run = MCMCRunner(InferenceConfig(**raw))
    rng = np.random.default_rng(20261027)
    idx = rng.choice(len(TH), size=min(a.n_check, len(TH)), replace=False)
    worst = 0.0
    for i in idx:
        td = run._expand_theta(TH[i])
        run.last_mass_closure = None
        iso = run._derive_gravity_isostatic(td)
        lib = run._derive_libration_deg(td)
        ref = np.array([iso[(2, 0)], iso[(2, 2)], iso[(3, 0)], lib])
        worst = max(worst, float(np.max(np.abs(X0[i] - ref) / np.abs(ref))))
    res['likelihood_path_max_rel_diff'] = worst
    res['likelihood_path_ok'] = bool(worst <= 1e-10)
    res['x_finite'] = bool(np.isfinite(X).all() and np.isfinite(TH).all())

    res['G0_pass'] = bool(res['seeds_ok'] and len(hashes) == 1 and len(caches) == 1
                          and res['servable_ok'] and res['noise_ok']
                          and res['likelihood_path_ok'] and res['x_finite'])
    (D / 'G0_report.json').write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
    if res['G0_pass']:
        np.savez(D / 'enceladus_ocean6D_train_1m.npz', theta=TH, x=X,
                 param_names=np.array(metas[0]['param_names']),
                 obs_names=np.array(obs),
                 provenance=json.dumps({'shards': metas, 'G0': res}, default=str))
        print('merged ->', D / 'enceladus_ocean6D_train_1m.npz')
    else:
        print('G0 FAILED -- no merged training set written')
        sys.exit(1)


if __name__ == '__main__':
    main()
