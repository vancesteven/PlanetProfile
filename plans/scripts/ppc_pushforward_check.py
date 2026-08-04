"""Body/artifact-agnostic posterior-predictive check (PPC) with four-way table.

Generalizes ``titanG_ppc_interior_check.py`` to ANY SBI artifact + reference
MCMC, per the manager advisory (2026-08-04, MACHINE-B-HANDOFF §0.5/§0.6): the
Titan NH3 campaign proved SBC-PASS + per-param crosscheck can coexist with a
datum-local pushforward under-update the per-param gate is structurally blind
to. This script is the diagnostic that catches it, run across the deployed
Europa artifacts (v4, v1.1 first — the public HF app) and then v5/v6/v7.

For each observable channel it reports FOUR medians plus the preregistered flag:

  * prior-predictive median  — draws from the PRIOR pushed through the forward
    model (from the training .npz if present, else regenerated from the config
    prior through the identical theta->x loop; same physics, noiseless).
  * SBI posterior-predictive  — SBI posterior draws at x_obs pushed forward.
  * MCMC posterior-predictive  — reference-MCMC samples (importance-WEIGHTED)
    pushed through the SAME forward loop.
  * observed  — the datum (config central, or --x-obs override).

PREREGISTERED FLAG (fixed BEFORE running, never moved after seeing results):
  A channel is FLAGGED when
      | median_pp(SBI) - median_pp(MCMC) | > FLAG_K * sigma_obs
  with FLAG_K = 0.5. Sanity target (NH3, validated 2026-08-03/04): flags BOTH
  k2 channels (Re ~0.83 sigma, Im ~1.46 sigma) and NEITHER gravity channel
  (~0.03 sigma). Run this script on the NH3 artifact to confirm the statistic
  reproduces that pattern before trusting it on a new artifact.

All pushforwards reuse ``MCMCRunner.generate_sbi_dataset(theta_override=...)``
(the theta->x path validated exact against the reference, max|fwd-stored|=0) —
no forward physics is re-implemented here. NOISELESS (obs_noise=False): we ask
whether retained models reproduce the data, not data + a fresh noise draw.

DIAGNOSTIC, not a gate. Interpreted by the scientific-reviewer; never tuned.

Run (v4 example):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/ppc_pushforward_check.py \
      --artifact PlanetProfile/Inference/sbi_artifacts/europa_clipper_v4_geodesy_11D_posterior_1m.pt \
      --config   PlanetProfile/Inference/configs/europa_clipper_v4_geodesy_11D.json \
      --mcmc     PlanetProfile/Test/mcmc_results/Europa/Test53_geodesy_v4/europa_clipper_v4_reference_result.pkl \
      --output-dir validation_reports/europa_clipper_v4_1m/ppc \
      --n-post 4000 --n-prior 10000 --seed 49
"""
import argparse
import json
import os
import pickle
from pathlib import Path

import numpy as np

FLAG_K = 0.5  # PREREGISTERED: |SBI_pp - MCMC_pp| median gap > FLAG_K * sigma_obs => flag


def _pctile_of(value, samples):
    """Percentile rank of ``value`` within finite ``samples`` (0-100)."""
    s = np.asarray(samples, dtype=float)
    s = s[np.isfinite(s)]
    if s.size == 0:
        return float("nan")
    return float(100.0 * np.mean(s <= value))


def _weighted_median_1d(values, weights):
    """Importance-weighted median (midpoint rule, matches validate_sbi)."""
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    ok = np.isfinite(v) & np.isfinite(w) & (w >= 0)
    v, w = v[ok], w[ok]
    if v.size == 0 or w.sum() <= 0:
        return float("nan")
    order = np.argsort(v)
    cw = np.cumsum(w[order])
    cw = cw / cw[-1]
    return float(v[order][np.searchsorted(cw, 0.5)])


def _weighted_quantile(values, weights, q):
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    ok = np.isfinite(v) & np.isfinite(w) & (w >= 0)
    v, w = v[ok], w[ok]
    if v.size == 0 or w.sum() <= 0:
        return float("nan")
    order = np.argsort(v)
    cw = np.cumsum(w[order]) - 0.5 * w[order]
    cw = cw / w.sum()
    return float(np.interp(q, cw, v[order]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--artifact", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--mcmc", default=None,
                    help="reference MCMC InferenceResult .pkl (for the MCMC-pp column)")
    ap.add_argument("--dataset", default=None,
                    help="optional training .npz for the prior-predictive column; "
                         "if absent, prior-pred is regenerated from the config prior")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--n-post", type=int, default=4000,
                    help="SBI posterior draws pushed through the forward model")
    ap.add_argument("--n-prior", type=int, default=10000,
                    help="prior draws for the prior-predictive column when no --dataset")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--x-obs", default=None,
                    help="optional JSON dict override; default = config central values")
    ap.add_argument("--label", default=None, help="artifact label for the report header")
    args = ap.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.sbi_runner import SBIRunner

    cfg_obj = InferenceConfig.from_json(args.config)
    obs_names = list(cfg_obj.observables.keys())
    sigmas = {n: float(cfg_obj.observables[n][1]) for n in obs_names}

    if args.x_obs:
        x_obs = {k: float(v) for k, v in json.loads(args.x_obs).items()}
    else:
        x_obs = {n: float(cfg_obj.observables[n][0]) for n in obs_names}

    label = args.label or Path(args.artifact).stem
    print(f"[ppc] artifact: {label}")
    print(f"[ppc] obs ({len(obs_names)}): {obs_names}")

    # Sampling runner (deploy/GUI path, reject_outside_prior=True). Forward
    # runner is config-built in sbi mode (the artifact carries no config/cache),
    # exactly as gen_dataset builds it -> identical theta->x physics.
    runner = SBIRunner.load_artifact(args.artifact)
    _cfg_dict = json.load(open(args.config))
    _cfg_dict["mode"] = "sbi"
    fwd_runner = SBIRunner(InferenceConfig.from_dict(_cfg_dict))
    assert list(fwd_runner.obs_names) == obs_names, "config obs order drift"
    mcmc = fwd_runner._get_mcmc_runner()
    guard = fwd_runner._support_guard_active()
    conv = fwd_runner.imag_convention

    def push(theta):
        """Push theta (N,P) through the identical forward loop; drop_nonfinite=False
        preserves row order so external weights stay aligned. Returns (theta_used, x)."""
        th, x = mcmc.generate_sbi_dataset(
            n_samples=len(theta), apply_support_guard=guard,
            imag_convention=conv, drop_nonfinite=False,
            obs_noise=False, theta_override=np.asarray(theta, dtype=float))
        return np.asarray(th, dtype=float), np.asarray(x, dtype=float)

    # ---- prior-predictive column -----------------------------------------
    prior_x = None
    prior_src = None
    if args.dataset and os.path.exists(args.dataset):
        dat = np.load(args.dataset, allow_pickle=True)
        prior_x = np.asarray(dat["x"], dtype=float)
        prior_src = f"dataset:{args.dataset}"
    else:
        print(f"[ppc] no dataset -> regenerating prior-predictive from config prior "
              f"(n={args.n_prior}, seed={args.seed})")
        theta_prior = fwd_runner.prior.rvs(args.n_prior) if hasattr(fwd_runner, "prior") \
            else mcmc.prior.rvs(args.n_prior)
        _, prior_x = push(theta_prior)
        prior_src = f"regenerated_from_prior:n={args.n_prior}"
    assert prior_x.shape[1] == len(obs_names)

    # ---- SBI posterior-predictive ----------------------------------------
    print(f"[ppc] drawing {args.n_post} SBI posterior samples at x_obs...")
    theta_sbi = runner.sample_posterior(x_obs, n_samples=args.n_post, seed=args.seed)
    _, x_sbi = push(theta_sbi)

    # ---- MCMC posterior-predictive (weighted) ----------------------------
    mcmc_x = None
    mcmc_w = None
    if args.mcmc and os.path.exists(args.mcmc):
        res = pickle.load(open(args.mcmc, "rb"))
        m_theta = np.asarray(res.samples, dtype=float)
        m_w = np.asarray(res.weights, dtype=float) if res.weights is not None \
            else np.ones(len(m_theta))
        # Align param order artifact <- reference MCMC.
        ref_names = list(res.param_names)
        art_names = list(fwd_runner.param_names)
        if ref_names != art_names:
            idx = [ref_names.index(p) for p in art_names]
            m_theta = m_theta[:, idx]
            print(f"[ppc] reordered MCMC params {ref_names} -> {art_names}")
        m_theta_used, mcmc_x = push(m_theta)   # order preserved (drop_nonfinite=False)
        mcmc_w = m_w
        print(f"[ppc] MCMC-pp: pushed {len(mcmc_x)} weighted reference samples")
    else:
        print("[ppc] WARNING: no --mcmc reference -> MCMC-pp column unavailable, "
              "flag statistic cannot be computed")

    # ---- four-way table + preregistered flag -----------------------------
    channels = {}
    n_flagged = 0
    for j, name in enumerate(obs_names):
        sig = sigmas[name]
        pri = prior_x[:, j]
        pri_med = float(np.median(pri[np.isfinite(pri)])) if np.isfinite(pri).any() else float("nan")
        sbi = x_sbi[:, j]
        sbi_med = float(np.median(sbi[np.isfinite(sbi)])) if np.isfinite(sbi).any() else float("nan")
        rec = {
            "obs": x_obs[name],
            "sigma_obs": sig,
            "prior_pred_median": pri_med,
            "sbi_pp_median": sbi_med,
            "obs_pctile_in_prior_pred": _pctile_of(x_obs[name], pri),
        }
        if mcmc_x is not None:
            mc = mcmc_x[:, j]
            mc_med = _weighted_median_1d(mc, mcmc_w)
            rec["mcmc_pp_median"] = mc_med
            gap = abs(sbi_med - mc_med)
            rec["sbi_minus_mcmc_abs"] = float(gap)
            rec["sbi_minus_mcmc_sigma"] = float(gap / sig) if sig > 0 else float("nan")
            rec["flagged"] = bool(gap > FLAG_K * sig)
            if rec["flagged"]:
                n_flagged += 1
        channels[name] = rec
        flg = "FLAG" if rec.get("flagged") else "  ok"
        gaps = f"{rec.get('sbi_minus_mcmc_sigma', float('nan')):.2f}"
        print(f"[ppc][{flg}] {name:26s} obs={x_obs[name]:+.4e} sig={sig:.3e}  "
              f"prior={pri_med:+.3e} sbi={sbi_med:+.3e} "
              f"mcmc={rec.get('mcmc_pp_median', float('nan')):+.3e}  "
              f"|sbi-mcmc|={gaps}sig")

    all_interior = all(5.0 <= c["obs_pctile_in_prior_pred"] <= 95.0
                       for c in channels.values()
                       if np.isfinite(c["obs_pctile_in_prior_pred"]))

    report = {
        "kind": "ppc_pushforward_four_way",
        "note": ("DIAGNOSTIC, not a gate. Noiseless pushforward (obs_noise=False). "
                 f"Preregistered flag: |median_pp(SBI)-median_pp(MCMC)| > {FLAG_K}*sigma_obs. "
                 "FLAG_K fixed before running; never moved after seeing results."),
        "label": label,
        "artifact": os.path.abspath(args.artifact),
        "config": os.path.abspath(args.config),
        "mcmc_reference": os.path.abspath(args.mcmc) if args.mcmc else None,
        "prior_predictive_source": prior_src,
        "flag_k": FLAG_K,
        "obs_names": obs_names,
        "x_obs": x_obs,
        "sigma": sigmas,
        "n_post": args.n_post,
        "seed": args.seed,
        "prior_predictive_all_interior_5_95": all_interior,
        "n_channels_flagged": n_flagged,
        "flagged_channels": [n for n, c in channels.items() if c.get("flagged")],
        "channels": channels,
    }
    out = outdir / "ppc_pushforward_report.json"
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n[ppc] {n_flagged}/{len(obs_names)} channels FLAGGED "
          f"(|SBI-MCMC| > {FLAG_K}sig); x_obs interior on all: {all_interior}")
    print(f"[ppc] report -> {out}")


if __name__ == "__main__":
    main()
