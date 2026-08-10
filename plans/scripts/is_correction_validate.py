#!/usr/bin/env python
"""Track 1 validation driver: IS-corrected amortized posterior vs the
pooled reference MCMC (plans/active/tidal-sector-remedy-plan.md,
reviewer conditions C1-C16).

Per campaign: condition the deployed flow at the fiducial x with FULL
derived recompute (n_derived=None), form importance weights
(is_correction.compute_is_correction), run the preregistered diagnostics
and — where the pooled reference pkl is present — the C3 likelihood
consistency test, the weighted pushforward gates (median-to-median +
weighted KS/W1, C10), the C20/C22 no-regression gate (C11), and the C16
ocean-fraction comparison. The crosscheck gate (all 13 params) runs via
systematic resampling to L = int(ESS) draws fed to the ratified
validate_sbi crosscheck (preregistered: resampled-set size is ESS, never
N — review 6e).

Usage:
  python plans/scripts/is_correction_validate.py --comp nh3 \
      [--n-post 20000] [--reference /path/to/pooled_reference.pkl] \
      [--out validation_reports/<campaign>/is_correction]

Never tune a threshold after seeing a failure: stop and surface.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

CAMPAIGNS = {
    'nh3': {
        'config': 'PlanetProfile/Inference/configs/'
                  'test54_titan_nh3_freegrav.json',
        'report_dir': 'validation_reports/titan_freegrav_nh3_1m',
        'reference': 'validation_reports/titan_freegrav_nh3_1m/reference/'
                     'titan_freegrav_nh3_reference_result.pkl',
    },
    'mgso4': {
        'config': 'PlanetProfile/Inference/configs/'
                  'test54_titan_mgso4_freegrav.json',
        'report_dir': 'validation_reports/titan_freegrav_mgso4_1m',
        'reference': 'validation_reports/titan_freegrav_mgso4_1m/'
                     'reference/titan_freegrav_mgso4_reference_pooled.pkl',
    },
    'nacl': {
        'config': 'PlanetProfile/Inference/configs/'
                  'test54_titan_nacl_freegrav.json',
        'report_dir': 'validation_reports/titan_freegrav_nacl_1m',
        'reference': 'validation_reports/titan_freegrav_nacl_1m/'
                     'reference/titan_freegrav_nacl_reference_pooled.pkl',
    },
}

# Gate 2 (C10) preregistered tolerances (reviewer 2026-08-11):
MEDIAN_GAP_SIGMA_MAX = 0.5      # corrected-SBI vs MCMC pp, median-to-median
WKS_D_MAX = 0.15                # weighted two-sample KS on Re/Im k2
GRAV_NOREG_SIGMA = 0.1          # C11: C20/C22 pp median shift budget


def weighted_ks(v1, w1, v2, w2):
    """Two-sample KS distance between weighted samples."""
    v1, v2 = np.asarray(v1, float), np.asarray(v2, float)
    w1 = np.asarray(w1, float) / np.sum(w1)
    w2 = np.asarray(w2, float) / np.sum(w2)
    grid = np.sort(np.concatenate([v1, v2]))
    c1 = np.cumsum(w1[np.argsort(v1)])
    c2 = np.cumsum(w2[np.argsort(v2)])
    f1 = np.interp(grid, np.sort(v1), c1, left=0.0, right=1.0)
    f2 = np.interp(grid, np.sort(v2), c2, left=0.0, right=1.0)
    return float(np.max(np.abs(f1 - f2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--comp', required=True, choices=sorted(CAMPAIGNS))
    ap.add_argument('--n-post', type=int, default=20000)
    ap.add_argument('--reference', default=None)
    ap.add_argument('--out', default=None)
    ap.add_argument('--seed-offset', type=int, default=0,
                    help='C13 3-seed stability: offset added to the '
                         'config random_state for the conditioning draw')
    args = ap.parse_args()

    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    from PlanetProfile.Inference.is_correction import (
        compute_is_correction, assert_byte_identity,
        reference_likelihood_consistency, weighted_quantile,
        systematic_resample)

    spec = CAMPAIGNS[args.comp]
    # Config stays byte-identical to the training config (mode flip only,
    # matching the artifact hash): n_post and seed are ARGUMENTS of the
    # conditioning call, never config mutations — a mutated config would
    # trip the C2 hash assert, by design.
    cfg = InferenceConfig.from_json(str(REPO / spec['config']))
    cfg.mode = 'sbi'
    run_seed = int(cfg.random_state) + args.seed_offset

    out_dir = Path(args.out) if args.out else (
        REPO / spec['report_dir'] / 'is_correction')
    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    runner = SBIRunner(cfg)
    base_hash = cfg.generate_hash()
    art = runner._find_matching_artifact(base_hash)
    if art is None:
        raise SystemExit(f'no artifact matches base hash {base_hash}')
    artifact = runner._load_artifact_dict(art)
    runner._posterior = artifact['posterior']
    runner._train_info = {k: artifact.get(k) for k in
                          ('density_estimator', 'seed',
                           'n_train_effective', 'rejection_stats',
                           'theta_norm', 'x_norm')}
    x_obs = {n: s[0] for n, s in cfg.observables.items()}
    res = runner._condition_and_package(
        x_obs=x_obs, n_posterior_samples=args.n_post,
        seed=run_seed,
        n_reeval=int(cfg.sampler_settings.get('n_reeval', 500)),
        artifact_path=art, reused=True,
        rejection_stats=artifact.get('rejection_stats'),
        config_hash=base_hash, t0=t0, n_derived=None)
    wall_condition = time.time() - t0

    # --- C2 byte identity ------------------------------------------------
    identity = assert_byte_identity(runner, res)

    # --- weights + core diagnostics (C4-C7, C16, C6) ---------------------
    k2_bounds = cfg.sampler_settings.get('k2_support_bounds') or {}
    corr = compute_is_correction(res, k2_support_bounds=k2_bounds)

    report = {
        'campaign': args.comp,
        'config_hash': base_hash,
        'artifact': str(art),
        'n_post': args.n_post,
        'seed': int(cfg.random_state),
        'wall_condition_s': wall_condition,
        'byte_identity': identity,
        'is_diagnostics': corr.summary(),
        'gates': {},
    }

    k2 = np.asarray(res.k2_results, float)
    re_v, im_v = k2[:, 0], np.abs(k2[:, 1])
    w = corr.weights
    flow_w = np.full_like(w, 1.0 / len(w))

    def _pp(v, wt):
        return float(weighted_quantile(v, wt, 0.5)[0])

    pp = {
        'im_flow': _pp(im_v, flow_w), 'im_corrected': _pp(im_v, w),
        're_flow': _pp(re_v, flow_w), 're_corrected': _pp(re_v, w),
    }
    report['pushforward_medians'] = pp

    # --- reference-dependent gates --------------------------------------
    ref_path = Path(args.reference) if args.reference else (
        REPO / spec['reference'])
    if ref_path.exists():
        import pickle
        with open(ref_path, 'rb') as f:
            ref = pickle.load(f)
        ref_w = np.asarray(getattr(ref, 'weights', None) if
                           getattr(ref, 'weights', None) is not None
                           else np.ones(len(ref.samples)), float)
        ref_w = ref_w / ref_w.sum()

        # C3: likelihood recompute consistency
        report['c3_consistency'] = reference_likelihood_consistency(
            runner, np.asarray(ref.samples, float),
            np.asarray(ref.log_likelihoods, float))

        # reference pushforward (k2 recomputed columns must exist)
        ref_k2 = np.asarray(getattr(ref, 'k2_results', []), float)
        gates = report['gates']
        if ref_k2.size:
            rim = np.abs(ref_k2[:, 1]); rre = ref_k2[:, 0]
            sig_im = float(cfg.observables['Im_k2'][1])
            sig_re = float(cfg.observables['Re_k2'][1])
            gap_im = abs(_pp(im_v, w) - _pp(rim, ref_w)) / sig_im
            gap_re = abs(_pp(re_v, w) - _pp(rre, ref_w)) / sig_re
            ok_f = np.isfinite(im_v); rok = np.isfinite(rim)
            wks_im = weighted_ks(im_v[ok_f], w[ok_f], rim[rok], ref_w[rok])
            wks_re = weighted_ks(re_v[ok_f], w[ok_f], rre[rok], ref_w[rok])
            gates['pushforward'] = {
                'gap_im_sigma': gap_im, 'gap_re_sigma': gap_re,
                'wks_im': wks_im, 'wks_re': wks_re,
                'median_gap_max': MEDIAN_GAP_SIGMA_MAX,
                'wks_max': WKS_D_MAX,
                'pass': bool(gap_im <= MEDIAN_GAP_SIGMA_MAX
                             and gap_re <= MEDIAN_GAP_SIGMA_MAX
                             and wks_im <= WKS_D_MAX
                             and wks_re <= WKS_D_MAX),
            }
            # C11 gravity no-regression
            g = {}
            for j, name in ((2, 'C20'), (3, 'C22')):
                if ref_k2.shape[1] > j and k2.shape[1] > j:
                    sig = float(cfg.observables[name][1])
                    g[name] = abs(_pp(k2[:, j], w)
                                  - _pp(ref_k2[:, j], ref_w)) / sig
            if g:
                gates['gravity_noregression'] = {
                    **g, 'max_sigma': GRAV_NOREG_SIGMA,
                    'pass': bool(all(v <= GRAV_NOREG_SIGMA
                                     for v in g.values()))}

        # C16 ocean fraction vs reference
        ref_doc = np.asarray(getattr(ref, 'D_ocean_results', []), float)
        if ref_doc.size:
            ref_frac = float(np.sum(
                ref_w[np.nan_to_num(ref_doc, nan=0.0) > 0.5]))
            gates['ocean_fraction'] = {
                'reference': ref_frac,
                'corrected': corr.branch['ocean']['prob_corrected'],
                'note': 'compare against the reference 3-seed spread '
                        '(preregistered on B before this result is read)',
            }

        # crosscheck via systematic resampling to L = int(ESS)
        L = int(corr.ess)
        idx = systematic_resample(w, L, seed=run_seed)
        np.save(out_dir / 'resampled_indices.npy', idx)
        report['crosscheck_resample'] = {
            'L': L, 'note': 'feed samples[idx] to validate_sbi crosscheck '
            'as the corrected-SBI sample set; report ESS, never N'}
    else:
        report['reference'] = f'ABSENT at {ref_path} — reference-dependent '\
            'gates (C3, pushforward, C11, C16, crosscheck) must run on ' \
            'the machine holding the pooled reference.'

    np.save(out_dir / 'log_weights.npy', corr.log_w)
    out_path = out_dir / f'is_validation_{args.comp}.json'
    with open(out_path, 'w') as f:
        json.dump(report, f, indent=1, default=float)
    print(json.dumps({k: report[k] for k in
                      ('is_diagnostics', 'pushforward_medians', 'gates')},
                     indent=1, default=float))
    print(f'report -> {out_path}')


if __name__ == '__main__':
    main()
