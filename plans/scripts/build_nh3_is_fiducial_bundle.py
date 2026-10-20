#!/usr/bin/env python
"""Build the shipped Titan NH3 fiducial importance-sampling (IS) bundle.

Re-runs the validated fiducial conditioning (config-central x_obs, N=20000,
seed = cfg.random_state + seed_offset, n_derived=None) exactly as
plans/scripts/is_correction_validate.py --comp nh3 does, computes the IS
correction, and ASSERTS the regenerated log-weights equal the committed
validated draws (validation_reports/titan_freegrav_nh3_1m/is_correction/
log_weights.npy). On mismatch it aborts without writing anything: the
shipped bundle must be the validated draws.

Writes PlanetProfile/Inference/sbi_artifacts/titan_freegrav_nh3_is_fiducial.npz
plus a JSON manifest sidecar (.json). Loader: PlanetProfile/Inference/
is_fiducial_bundle.py.

Usage: python plans/scripts/build_nh3_is_fiducial_bundle.py [--n-post 20000]
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

CONFIG = 'PlanetProfile/Inference/configs/test54_titan_nh3_freegrav.json'
REPORT_DIR = REPO / 'validation_reports/titan_freegrav_nh3_1m/is_correction'
OUT_NPZ = REPO / 'PlanetProfile/Inference/sbi_artifacts/' \
    'titan_freegrav_nh3_is_fiducial.npz'


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-post', type=int, default=20000)
    ap.add_argument('--seed-offset', type=int, default=0)
    ap.add_argument('--out', default=str(OUT_NPZ))
    args = ap.parse_args()

    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    from PlanetProfile.Inference.is_correction import (
        compute_is_correction, assert_byte_identity)

    cfg = InferenceConfig.from_json(str(REPO / CONFIG))
    cfg.mode = 'sbi'
    run_seed = int(cfg.random_state) + args.seed_offset
    t0 = time.time()
    runner = SBIRunner(cfg)
    base_hash = cfg.generate_hash()
    art = runner._find_matching_artifact(base_hash)
    if art is None:
        raise SystemExit(f'no artifact matches base hash {base_hash}')
    artifact = runner._load_artifact_dict(art)
    runner._posterior = artifact['posterior']
    runner._train_info = {k: artifact.get(k) for k in
                          ('density_estimator', 'seed', 'n_train_effective',
                           'rejection_stats', 'theta_norm', 'x_norm')}
    x_obs = {n: float(s[0]) for n, s in cfg.observables.items()}
    res = runner._condition_and_package(
        x_obs=x_obs, n_posterior_samples=args.n_post, seed=run_seed,
        n_reeval=int(cfg.sampler_settings.get('n_reeval', 500)),
        artifact_path=art, reused=True,
        rejection_stats=artifact.get('rejection_stats'),
        config_hash=base_hash, t0=t0, n_derived=None)
    identity = assert_byte_identity(runner, res)
    k2_bounds = cfg.sampler_settings.get('k2_support_bounds') or {}
    corr = compute_is_correction(res, k2_support_bounds=k2_bounds)

    # ---- hard assert: regenerated draws == committed validated draws ----
    ref_lw = np.load(REPORT_DIR / 'log_weights.npy')
    lw = corr.log_w
    if lw.shape != ref_lw.shape:
        raise SystemExit(f'ABORT: log_w shape {lw.shape} != {ref_lw.shape}')
    if not np.array_equal(np.isneginf(lw), np.isneginf(ref_lw)):
        raise SystemExit('ABORT: -inf mask differs from committed log_w')
    fin = np.isfinite(ref_lw)
    if not np.allclose(lw[fin], ref_lw[fin], rtol=1e-10, atol=0):
        d = np.abs(lw[fin] - ref_lw[fin])
        raise SystemExit(
            f'ABORT: log_w differs from committed; max abs diff {d.max():.3e}')
    max_abs = float(np.max(np.abs(lw[fin] - ref_lw[fin])))
    max_rel = float(np.max(np.abs(lw[fin] - ref_lw[fin])
                           / np.maximum(np.abs(ref_lw[fin]), 1e-300)))
    print(f'log_w match: OK  n_finite={int(fin.sum())} '
          f'n_-inf={int((~fin).sum())} max_abs={max_abs:.3e} '
          f'max_rel={max_rel:.3e}')

    # ---- assemble arrays ------------------------------------------------
    N = len(lw)
    arrays = {
        'samples': np.asarray(res.samples, float),
        'log_w': lw,
        'log_likelihoods': np.asarray(res.log_likelihoods, float),
        'flow_log_prob': np.asarray(res.metadata['flow_log_prob'], float),
        'param_names': np.asarray(list(res.param_names)),
    }
    carried = []
    for name in ('k2_results', 'h2_results', 'cmr2_results',
                 'D_ocean_results', 'D_iceIh_results', 'D_hsphere_results',
                 'D_iceHP_results', 'D_clath_results', 'c20_results',
                 'c22_results', 'cmr2_hydro_results'):
        v = getattr(res, name, None)
        if v is not None and np.size(v):
            v = np.asarray(v, float)
            if len(v) == N:
                arrays[name] = v
                carried.append(name)
    dph = getattr(res, 'D_icePhase_results', None)
    if dph:
        for k, v in dph.items():
            arrays[f'D_icePhase_{k}'] = np.asarray(v, float)
            carried.append(f'D_icePhase_{k}')
    dr = arrays.get('D_ocean_results')
    if dr is not None:
        arrays['has_ocean'] = np.nan_to_num(dr, nan=0.0) > 0.5
    hr = getattr(res, 'heating_results', None)
    heating_note = None
    if hr is not None:
        heating_note = (f'heating_results: list of {len(hr)} dicts '
                        f'(len != N={N}; not shipped)')
        if len(hr) == N:
            keys = sorted({k for d in hr for k in d})
            for k in keys:
                arrays[f'heating_{k}'] = np.array(
                    [d.get(k, np.nan) for d in hr], float)
            carried.append('heating_*')
            heating_note = None
    np.savez_compressed(args.out, **arrays)
    npz_sha = sha256_file(args.out)

    val = json.load(open(REPORT_DIR / 'is_validation_nh3.json'))
    import torch
    try:
        import sbi
        sbi_v = sbi.__version__
    except Exception:
        sbi_v = None
    obs_sigma = {n: float(s[1]) for n, s in cfg.observables.items()}
    diag = corr.summary()
    manifest = {
        'bundle': 'titan_freegrav_nh3_is_fiducial',
        'npz_file': Path(args.out).name,
        'npz_sha256': npz_sha,
        'artifact_sha256': sha256_file(art),
        'artifact_name': Path(art).name,
        'config_file': CONFIG,
        'config_hash': base_hash,
        'structure_cache_sha256': identity['structure_cache_sha256'],
        'x_obs': x_obs,
        'obs_sigma': obs_sigma,
        'observables': sorted(x_obs),
        'imag_convention': identity['imag_convention'],
        'param_names': list(res.param_names),
        'n_samples': N,
        'config_random_state': int(cfg.random_state),
        'seed_offset': args.seed_offset,
        'effective_seed': run_seed,
        'n_derived': None,
        'arrays': sorted(arrays),
        'derived_carried': carried,
        'heating_note': heating_note,
        'diagnostics': {
            'ess': diag['ess'], 'pareto_k': diag['pareto_k'],
            'w_max_frac': diag['w_max_frac'],
            'frac_rejected': diag['frac_rejected'],
            'n_rejected': diag['n_rejected'],
            'verdict': diag['verdict'],
            'ocean_prob_flow': diag['branch']['ocean']['prob_flow'],
            'ocean_prob_corrected': diag['branch']['ocean']['prob_corrected'],
        },
        'c5_3': val['gates']['c5_3_reverse_coverage'],
        'logw_match': {'n_finite': int(fin.sum()),
                       'n_neg_inf': int((~fin).sum()),
                       'max_abs_diff': max_abs, 'max_rel_diff': max_rel,
                       'rtol': 1e-10},
        'versions': {'torch': torch.__version__, 'sbi': sbi_v,
                     'numpy': np.__version__},
    }
    mpath = Path(args.out).with_suffix('.json')
    with open(mpath, 'w') as f:
        json.dump(manifest, f, indent=1, default=float)
    print(f'npz -> {args.out}  ({Path(args.out).stat().st_size/1e6:.2f} MB) '
          f'sha256 {npz_sha}')
    print(f'manifest -> {mpath}')
    print('arrays:', sorted(arrays))
    if heating_note:
        print(heating_note)


if __name__ == '__main__':
    main()
