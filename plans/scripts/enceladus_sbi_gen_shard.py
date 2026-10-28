#!/usr/bin/env python
"""One shard of the Enceladus ocean-branch SBI training set
(validation_reports/enceladus_isostasy/sbi/DESIGN_PREREGISTRATION.md).

Shard k: data seed BASE_DATA + k, noise seed BASE_NOISE + k (training:
75 / 7575; held-out G3/G5 pairs: 900 / 9900). Draws from the nominal prior
through SBIRunner.generate_training_set -> MCMCRunner.generate_sbi_dataset
(drop_nonfinite=True drops cache-hole draws: the declared EFFECTIVE prior),
obs_noise=True with the likelihood's inflated, correlated covariance.
Writes theta, x (noisy), and the NOISE-FREE x0 for the same theta (needed by
G0 and the G2 pushforward; regenerated through theta_override).

  python plans/scripts/enceladus_sbi_gen_shard.py --k K --n-request N \
      --out-dir DIR [--held-out]
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402
from PlanetProfile.Inference.mcmc_runner import MCMCRunner  # noqa: E402

CFG = REPO / 'PlanetProfile/Inference/configs/enceladus_cassini_isostasy_ocean6D_rhorock.json'


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--n-request', type=int, required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--cache', required=True)
    ap.add_argument('--held-out', action='store_true')
    ap.add_argument('--zb-max', type=float, default=None,
                    help='restrict the zb prior upper bound (G5 thin-shell set)')
    a = ap.parse_args()
    base_d, base_n = (900, 9900) if a.held_out else (75, 7575)
    if a.zb_max is not None:
        base_d, base_n = 950, 9950
    ds, ns = base_d + a.k, base_n + a.k

    raw = json.loads(CFG.read_text())
    raw['structure_cache_path'] = a.cache
    if a.zb_max is not None:
        raw['param_space']['zb_km']['bounds'] = [raw['param_space']['zb_km']['bounds'][0],
                                                 float(a.zb_max)]
    cfg = InferenceConfig(**raw)
    run = MCMCRunner(cfg)
    obs = list(cfg.observables.keys())
    t0 = time.time()
    theta, x = run.generate_sbi_dataset(
        n_samples=a.n_request, seed=ds, drop_nonfinite=True,
        imag_convention='abs', obs_noise=True, noise_seed=ns)
    stats = dict(getattr(run, 'last_sbi_dataset_stats', {}) or {})
    _, x0 = run.generate_sbi_dataset(theta_override=theta, drop_nonfinite=False,
                                     imag_convention='abs')
    dt = time.time() - t0
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tag = ('heldout' if a.held_out else 'train') if a.zb_max is None else 'thin'
    f = out / f'{tag}_shard_{a.k:02d}.npz'
    meta = {'k': a.k, 'data_seed': ds, 'noise_seed': ns, 'n_request': a.n_request,
            'n_kept': int(len(theta)), 'kept_fraction': len(theta) / a.n_request,
            'obs_names': obs, 'param_names': list(run.param_names),
            'config': str(CFG), 'config_hash': cfg.generate_hash(),
            'cache': a.cache, 'cache_sha256': sha256(a.cache),
            'zb_max': a.zb_max, 'wall_s': dt, 'generator_stats': stats}
    np.savez(f, theta=theta, x=x, x0=x0, meta=json.dumps(meta, default=str))
    print(json.dumps({k: meta[k] for k in ('k', 'data_seed', 'noise_seed', 'n_kept',
                                           'kept_fraction', 'wall_s')}))


if __name__ == '__main__':
    main()
