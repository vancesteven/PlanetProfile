#!/usr/bin/env python
"""Train the Enceladus ocean-branch nsf NPE artifact (DESIGN_PREREGISTRATION.md).

Deployed defaults (hidden 50, 5 transforms, 10 bins, z-score independent),
train seed 75. Refuses to run unless G0_report.json says G0_pass. MUST run in
a separate process from dataset generation / PlanetProfile forward models
(libomp hazard: torch vs numba). Capacity and N are FIXED by the
preregistration -- do not change them after any gate result.

  python plans/scripts/enceladus_sbi_train.py --dir ~/scratch/enc_sbi \
      [--n-sub N --name NAME]   # optional non-gating 300k diagnostic
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
from PlanetProfile.Inference.sbi_runner import SBIRunner  # noqa: E402

CFG = REPO / 'PlanetProfile/Inference/configs/enceladus_cassini_isostasy_ocean6D_rhorock.json'
TRAIN_SEED = 75


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', required=True)
    ap.add_argument('--cache', required=True)
    ap.add_argument('--n-sub', type=int, default=None)
    ap.add_argument('--name', default='enceladus_ocean6D_rhorock_posterior_1m.pt')
    a = ap.parse_args()
    D = Path(a.dir)
    g0 = json.loads((D / 'G0_report.json').read_text())
    if not g0.get('G0_pass'):
        raise SystemExit('G0 did not pass -- refusing to train')
    ds = D / 'enceladus_ocean6D_train_1m.npz'
    d = np.load(ds)
    theta = np.asarray(d['theta'], np.float64)
    x = np.asarray(d['x'], np.float64)
    if a.n_sub:
        theta, x = theta[:a.n_sub], x[:a.n_sub]

    raw = json.loads(CFG.read_text())
    raw['mode'] = 'sbi'
    raw['structure_cache_path'] = a.cache
    runner = SBIRunner(InferenceConfig.from_dict(raw))
    assert list(runner.obs_names) == [str(o) for o in d['obs_names']]
    assert list(runner.param_names) == [str(p) for p in d['param_names']]

    t0 = time.time()
    runner.train(theta, x, seed=TRAIN_SEED, density_estimator='nsf')
    runner._train_info['dataset'] = {'path': str(ds), 'sha256': sha256(ds),
                                     'n_train': int(len(theta)), 'G0': g0}
    out = D / a.name
    runner.save_artifact(str(out))
    dt = time.time() - t0
    man = {'artifact': str(out), 'artifact_sha256': sha256(out),
           'train_seed': TRAIN_SEED, 'density_estimator': 'nsf',
           'capacity': 'deployed defaults (hidden 50, 5 transforms, 10 bins)',
           'n_train': int(len(theta)), 'train_seconds': dt,
           'config': str(CFG), 'config_hash': runner.config.generate_hash(),
           'cache_sha256': g0['single_cache_sha'][0],
           'dataset_sha256': runner._train_info['dataset']['sha256'],
           'train_info': {k: v for k, v in runner._train_info.items()
                          if isinstance(v, (int, float, str, bool, list))}}
    (D / (a.name + '.manifest.json')).write_text(json.dumps(man, indent=1, default=str))
    print(json.dumps({k: man[k] for k in ('artifact', 'n_train', 'train_seconds',
                                          'artifact_sha256')}, indent=1))


if __name__ == '__main__':
    main()
