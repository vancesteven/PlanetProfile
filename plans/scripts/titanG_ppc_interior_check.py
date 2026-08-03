"""Posterior-predictive coverage + prior-predictive interior check.

Two reviewer-recommended DIAGNOSTICS (NOT gates — nothing is tuned to these)
for the Titan JOINT no-ocean+ocean free-gravity SBI campaign, folded in after
the scientific-reviewer adjudicated the user's "k2 has an extraordinarily broad
range" observation (2026-08-03):

  (1) Posterior-predictive coverage at the fiducial.
      Draw N samples from the trained SBI posterior at x_obs, push them BACK
      through the identical theta->x forward loop the training set uses
      (MCMCRunner.generate_sbi_dataset(theta_override=...), so the gravity /
      k2 / support-guard physics is byte-identical — no re-implementation),
      and report the fraction of the noiseless posterior-predictive
      {C20,C22,Re_k2,Im_k2} landing within 1/2/3 sigma of x_obs, plus the
      posterior-predictive median vs observed. This converts the user's
      qualitative PRIOR-predictive worry into a quantitative POSTERIOR-
      predictive statement: after conditioning, do retained models reproduce
      the data?

  (2) Prior-predictive interior check.
      From the training .npz, report the percentile of each observed value
      within the prior-predictive marginal. x_obs interior (~5-95 pctile on
      all channels) => the flow is INTERPOLATING at the fiducial, not
      extrapolating to an out-of-envelope datum. This is the ~30-second guard
      that MgSO4 (1-194 ppt) / NaCl (1-300 ppt) — whose ocean EOS shifts the
      envelope — are still asked to interpolate.

Reusable across compositions via CLI flags. Interpreted by the
scientific-reviewer; never tuned to pass.

Run (after gen + train complete for the composition):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/titanG_ppc_interior_check.py \
      --artifact PlanetProfile/Inference/sbi_artifacts/titan_freegrav_nh3_posterior_1m.pt \
      --config   PlanetProfile/Inference/configs/test54_titan_nh3_freegrav.json \
      --dataset  /tmp/titanG_build/datasets/titanG_nh3_1m.npz \
      --output-dir validation_reports/titan_freegrav_nh3_1m/ppc \
      --n-post 4000 --seed 72
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np


def _pctile_of(value, samples):
    """Percentile rank of `value` within finite `samples` (0-100)."""
    s = np.asarray(samples, dtype=float)
    s = s[np.isfinite(s)]
    if s.size == 0:
        return float("nan")
    return float(100.0 * np.mean(s <= value))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--artifact", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--dataset", required=True,
                    help="training .npz (theta,x) for the prior-predictive check")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--n-post", type=int, default=4000,
                    help="posterior draws pushed through the forward model")
    ap.add_argument("--seed", type=int, default=72)
    ap.add_argument("--x-obs", default=None,
                    help="optional JSON dict override; default = config central values")
    args = ap.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.sbi_runner import SBIRunner

    cfg_obj = InferenceConfig.from_json(args.config)
    obs_names = list(cfg_obj.observables.keys())
    sigmas = {n: float(cfg_obj.observables[n][1]) for n in obs_names}

    # x_obs: config centrals unless overridden.
    if args.x_obs:
        x_obs = {k: float(v) for k, v in json.loads(args.x_obs).items()}
    else:
        x_obs = {n: float(cfg_obj.observables[n][0]) for n in obs_names}
    print(f"[ppc] obs ({len(obs_names)}): {obs_names}")
    print(f"[ppc] x_obs = {x_obs}")
    print(f"[ppc] sigma = {sigmas}")

    # -----------------------------------------------------------------
    # (2) Prior-predictive interior check (cheap; do first).
    # -----------------------------------------------------------------
    dat = np.load(args.dataset, allow_pickle=True)
    x_train = np.asarray(dat["x"], dtype=float)  # (n, K), obs order
    assert x_train.shape[1] == len(obs_names), \
        f"dataset x has {x_train.shape[1]} cols, config has {len(obs_names)} obs"
    interior = {}
    for j, name in enumerate(obs_names):
        col = x_train[:, j]
        pct = _pctile_of(x_obs[name], col)
        finite = col[np.isfinite(col)]
        interior[name] = {
            "obs": x_obs[name],
            "pctile_in_prior_predictive": pct,
            "interior_5_95": bool(5.0 <= pct <= 95.0),
            "train_median": float(np.median(finite)) if finite.size else float("nan"),
            "train_5_95": [float(np.percentile(finite, 5)),
                           float(np.percentile(finite, 95))] if finite.size else None,
            "n_finite": int(finite.size),
        }
        flag = "interior" if interior[name]["interior_5_95"] else "TAIL(!)"
        print(f"[ppc][prior] {name:8s} obs={x_obs[name]:+.4e}  "
              f"pctile={pct:5.1f}  [{flag}]")
    all_interior = all(v["interior_5_95"] for v in interior.values())
    print(f"[ppc][prior] x_obs interior on ALL channels: {all_interior}")

    # -----------------------------------------------------------------
    # (1) Posterior-predictive coverage at the fiducial.
    # -----------------------------------------------------------------
    # Sampling uses the artifact-loaded runner (flow only; deploy/GUI path,
    # reject_outside_prior=True). The forward model needs a CONFIG-built
    # runner (the artifact carries no config/cache), exactly as gen_dataset
    # builds it — this is the SAME runner class + config that produced the
    # training set, so the theta->x physics is identical.
    runner = SBIRunner.load_artifact(args.artifact)
    # Build the forward runner in sbi mode from the same config file
    # (gen_dataset does the identical mode override).
    _cfg_dict = json.load(open(args.config))
    _cfg_dict["mode"] = "sbi"
    fwd_runner = SBIRunner(InferenceConfig.from_dict(_cfg_dict))
    assert list(fwd_runner.param_names) == list(runner.param_names), (
        f"param order mismatch: artifact {list(runner.param_names)} vs "
        f"config {list(fwd_runner.param_names)}")
    assert list(fwd_runner.obs_names) == obs_names, "config obs order drift"

    print(f"[ppc] drawing {args.n_post} posterior samples at x_obs...")
    theta_post = runner.sample_posterior(
        x_obs, n_samples=args.n_post, seed=args.seed)  # reject_outside_prior=True (deploy path)
    print(f"[ppc] posterior theta: {theta_post.shape}")

    # Push posterior theta BACK through the identical theta->x forward loop.
    # obs_noise=False: we want the NOISELESS posterior-predictive (does the
    # model reproduce the data), not data + a fresh noise realization.
    mcmc = fwd_runner._get_mcmc_runner()
    theta_pp, x_pp = mcmc.generate_sbi_dataset(
        n_samples=len(theta_post),
        apply_support_guard=fwd_runner._support_guard_active(),
        imag_convention=fwd_runner.imag_convention,
        drop_nonfinite=True,
        obs_noise=False,
        theta_override=theta_post,
    )
    x_pp = np.asarray(x_pp, dtype=float)
    n_forward = len(x_pp)
    print(f"[ppc] posterior-predictive x: {x_pp.shape} "
          f"(dropped {len(theta_post) - n_forward} nonfinite/support)")

    # Per-channel coverage + median.
    ppc = {}
    for j, name in enumerate(obs_names):
        col = x_pp[:, j]
        col = col[np.isfinite(col)]
        if col.size == 0:
            ppc[name] = {"n_finite": 0}
            continue
        sig = sigmas[name]
        dev = np.abs(col - x_obs[name]) / sig
        ppc[name] = {
            "obs": x_obs[name],
            "sigma": sig,
            "pp_median": float(np.median(col)),
            "pp_5_95": [float(np.percentile(col, 5)),
                        float(np.percentile(col, 95))],
            "median_dev_sigma": float(np.median(dev)),
            "frac_within_1sig": float(np.mean(dev <= 1.0)),
            "frac_within_2sig": float(np.mean(dev <= 2.0)),
            "frac_within_3sig": float(np.mean(dev <= 3.0)),
            "n_finite": int(col.size),
        }
        print(f"[ppc][post ] {name:8s} obs={x_obs[name]:+.4e}  "
              f"pp_med={ppc[name]['pp_median']:+.4e}  "
              f"med_dev={ppc[name]['median_dev_sigma']:.2f}sig  "
              f"within1/2/3sig={ppc[name]['frac_within_1sig']:.2f}/"
              f"{ppc[name]['frac_within_2sig']:.2f}/{ppc[name]['frac_within_3sig']:.2f}")

    # Joint coverage (all channels simultaneously within k-sigma).
    finite_all = np.all(np.isfinite(x_pp), axis=1)
    xa = x_pp[finite_all]
    obs_vec = np.array([x_obs[n] for n in obs_names])
    sig_vec = np.array([sigmas[n] for n in obs_names])
    dev_all = np.abs(xa - obs_vec) / sig_vec
    joint = {
        "n_finite": int(len(xa)),
        "frac_all_within_1sig": float(np.mean(np.all(dev_all <= 1.0, axis=1))),
        "frac_all_within_2sig": float(np.mean(np.all(dev_all <= 2.0, axis=1))),
        "frac_all_within_3sig": float(np.mean(np.all(dev_all <= 3.0, axis=1))),
    }
    print(f"[ppc][post ] JOINT all-channels within 1/2/3sig="
          f"{joint['frac_all_within_1sig']:.3f}/{joint['frac_all_within_2sig']:.3f}/"
          f"{joint['frac_all_within_3sig']:.3f}")

    report = {
        "kind": "posterior_predictive_coverage_and_interior_check",
        "note": ("DIAGNOSTIC, not a gate — interpreted by scientific-reviewer, "
                 "never tuned to pass. Posterior-predictive is NOISELESS "
                 "(obs_noise=False); coverage vs config sigma."),
        "artifact": os.path.abspath(args.artifact),
        "config": os.path.abspath(args.config),
        "dataset": os.path.abspath(args.dataset),
        "obs_names": obs_names,
        "x_obs": x_obs,
        "sigma": sigmas,
        "n_post_requested": args.n_post,
        "n_post_forward_finite": n_forward,
        "seed": args.seed,
        "prior_predictive_interior": interior,
        "prior_predictive_all_interior_5_95": all_interior,
        "posterior_predictive": ppc,
        "posterior_predictive_joint": joint,
    }
    out = outdir / "ppc_interior_report.json"
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n[ppc] report -> {out}")


if __name__ == "__main__":
    main()
