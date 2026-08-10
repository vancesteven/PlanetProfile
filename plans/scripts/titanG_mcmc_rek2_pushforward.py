"""MCMC -> Re_k2 pushforward (closes reviewer required-validation #1).

Scientific-reviewer verdict on the Titan free-gravity MgSO4/NaCl §0.16 gates
(2026-08-10): the crosscheck validated SBI-vs-MCMC in parameter MARGINAL space,
but the Re_k2 under-prediction the split-status deployment rests on is a
pushforward-of-the-JOINT statistic. Two posteriors can match on all 13
marginals yet differ in correlation structure enough to shift a 1.3-1.5 sigma
pushforward. So before deploying Re_k2 as "informative-with-caveat", the
reference MCMC's Re_k2 must be demonstrated on data, not inferred.

This pushes the POOLED REPAIRED reference-MCMC posterior through the IDENTICAL
theta->x forward loop the training set + the SBI ppc use
(MCMCRunner.generate_sbi_dataset(theta_override=...), obs_noise=False), so the
gravity / k2 / support-guard physics is byte-identical to titanG_ppc_interior_
check.py — no re-implementation. The pooled samples carry pocoMC importance
weights; we resample to unweighted with those weights (systematic/stratified
resampling, seeded) before the pushforward.

Decision rule (reviewer):
  MCMC under-predicts Re_k2 within ~0.5 sigma of the SBI pushforward value
    => model-data tension is reproduced by an independent sampler => Re_k2 ACCEPT
       (deploy informative-with-caveat).
  MCMC pushforward centers on the observed 0.608
    => the SBI under-prediction is a flow offset => escalate Re_k2 to quarantine.

DIAGNOSTIC/validation, interpreted by the scientific-reviewer; never tuned.

Run (after gen + train + reference MCMC + repool complete):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
    python plans/scripts/titanG_mcmc_rek2_pushforward.py --comp MgSO4
"""
import argparse
import json
import os
import pickle
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
COMPS = {
    "MgSO4": {"cfg": "PlanetProfile/Inference/configs/test54_titan_mgso4_freegrav.json",
              "tag": "mgso4", "seed": 73},
    "NaCl":  {"cfg": "PlanetProfile/Inference/configs/test54_titan_nacl_freegrav.json",
              "tag": "nacl", "seed": 74},
}


def _systematic_resample(weights, n_out, rng):
    """Systematic (low-variance) resample of indices from importance weights."""
    w = np.asarray(weights, dtype=float)
    w = w / w.sum()
    positions = (rng.random() + np.arange(n_out)) / n_out
    cumsum = np.cumsum(w)
    cumsum[-1] = 1.0  # guard fp drift
    return np.searchsorted(cumsum, positions)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--comp", required=True, choices=list(COMPS))
    ap.add_argument("--n-push", type=int, default=8000,
                    help="unweighted resampled MCMC draws pushed through the forward model")
    ap.add_argument("--seed", type=int, default=None,
                    help="resample+forward seed; defaults to the campaign seed")
    args = ap.parse_args()

    spec = COMPS[args.comp]
    tag = spec["tag"]
    seed = args.seed if args.seed is not None else spec["seed"]
    cfg_path = ROOT / spec["cfg"]
    pooled_pkl = (ROOT / f"validation_reports/titan_freegrav_{tag}_1m/reference"
                        / f"titan_freegrav_{tag}_reference_pooled.pkl")
    outdir = ROOT / f"validation_reports/titan_freegrav_{tag}_1m/rek2_pushforward"
    outdir.mkdir(parents=True, exist_ok=True)
    ppc_report = (ROOT / f"validation_reports/titan_freegrav_{tag}_1m/ppc"
                        / "ppc_interior_report.json")
    assert pooled_pkl.exists(), f"missing pooled reference {pooled_pkl}"
    assert cfg_path.exists(), f"missing config {cfg_path}"

    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.sbi_runner import SBIRunner

    cfg_obj = InferenceConfig.from_json(str(cfg_path))
    obs_names = list(cfg_obj.observables.keys())
    sigmas = {n: float(cfg_obj.observables[n][1]) for n in obs_names}
    x_obs = {n: float(cfg_obj.observables[n][0]) for n in obs_names}
    assert obs_names == ["C20", "C22", "Re_k2", "Im_k2"], f"unexpected obs: {obs_names}"

    # Build the forward runner in sbi mode from the same config file
    # (gen_dataset + titanG_ppc_interior_check do the identical mode override).
    _cfg_dict = json.load(open(cfg_path))
    _cfg_dict["mode"] = "sbi"
    fwd_runner = SBIRunner(InferenceConfig.from_dict(_cfg_dict))
    fwd_param_names = list(fwd_runner.param_names)
    assert list(fwd_runner.obs_names) == obs_names, "config obs order drift"

    # Load pooled MCMC posterior (samples + importance weights).
    o = pickle.load(open(pooled_pkl, "rb"))
    mcmc_param_names = list(o.param_names)
    samples = np.asarray(o.samples, dtype=float)     # (N, 13)
    weights = np.asarray(o.weights, dtype=float)     # (N,)
    assert samples.shape[0] == weights.shape[0], \
        f"weights {weights.shape} != samples {samples.shape}"
    assert abs(weights.sum() - 1.0) < 1e-6, f"pooled weights sum={weights.sum()}"

    # Reorder MCMC columns into the forward runner's param order (defensive:
    # both are the joint 13-param set, but never assume identical ordering).
    assert sorted(mcmc_param_names) == sorted(fwd_param_names), (
        f"param-set mismatch MCMC {mcmc_param_names} vs fwd {fwd_param_names}")
    col_for = [mcmc_param_names.index(p) for p in fwd_param_names]
    samples_ord = samples[:, col_for]

    # Resample to unweighted with the importance weights, then push through the
    # identical forward loop (obs_noise=False => noiseless pushforward).
    rng = np.random.default_rng(int(seed))
    idx = _systematic_resample(weights, args.n_push, rng)
    theta_push = samples_ord[idx]
    print(f"[mcmc-push] {args.comp}: pooled N={samples.shape[0]} "
          f"-> resampled {theta_push.shape[0]} unweighted draws (seed={seed})")

    mcmc = fwd_runner._get_mcmc_runner()
    theta_fwd, x_fwd = mcmc.generate_sbi_dataset(
        n_samples=len(theta_push),
        apply_support_guard=fwd_runner._support_guard_active(),
        imag_convention=fwd_runner.imag_convention,
        drop_nonfinite=True,
        obs_noise=False,
        theta_override=theta_push,
    )
    x_fwd = np.asarray(x_fwd, dtype=float)
    n_forward = len(x_fwd)
    print(f"[mcmc-push] forward x: {x_fwd.shape} "
          f"(dropped {len(theta_push) - n_forward} nonfinite/support)")

    # Per-channel pushforward stats (Re_k2 is the target; report all four).
    per_channel = {}
    for j, name in enumerate(obs_names):
        col = x_fwd[:, j]
        col = col[np.isfinite(col)]
        if col.size == 0:
            per_channel[name] = {"n_finite": 0}
            continue
        sig = sigmas[name]
        signed = (np.median(col) - x_obs[name]) / sig  # + = over-predict
        dev = np.abs(col - x_obs[name]) / sig
        per_channel[name] = {
            "obs": x_obs[name], "sigma": sig,
            "pp_median": float(np.median(col)),
            "pp_5_95": [float(np.percentile(col, 5)), float(np.percentile(col, 95))],
            "signed_median_dev_sigma": float(signed),
            "median_abs_dev_sigma": float(np.median(dev)),
            "n_finite": int(col.size),
        }

    # Re_k2 decision arithmetic vs the SBI pushforward.
    rek2 = per_channel["Re_k2"]
    sbi_rek2 = None
    if ppc_report.exists():
        pj = json.load(open(ppc_report))
        sbi_rek2 = pj.get("posterior_predictive", {}).get("Re_k2", {})
    decision = {
        "mcmc_pp_median_Re_k2": rek2["pp_median"],
        "mcmc_signed_dev_sigma": rek2["signed_median_dev_sigma"],
        "mcmc_median_abs_dev_sigma": rek2["median_abs_dev_sigma"],
        "sbi_pp_median_Re_k2": sbi_rek2.get("pp_median") if sbi_rek2 else None,
        "sbi_median_dev_sigma": sbi_rek2.get("median_dev_sigma") if sbi_rek2 else None,
    }
    if sbi_rek2 and sbi_rek2.get("pp_median") is not None:
        # |MCMC median - SBI median| in units of obs sigma.
        d_mcmc_vs_sbi = abs(rek2["pp_median"] - sbi_rek2["pp_median"]) / sigmas["Re_k2"]
        decision["mcmc_vs_sbi_median_gap_sigma"] = float(d_mcmc_vs_sbi)
        decision["mcmc_within_0p5sigma_of_sbi"] = bool(d_mcmc_vs_sbi <= 0.5)
        # under-prediction preserved (both below obs)?
        decision["mcmc_under_predicts"] = bool(rek2["pp_median"] < x_obs["Re_k2"])
        if d_mcmc_vs_sbi <= 0.5 and rek2["pp_median"] < x_obs["Re_k2"]:
            decision["verdict_arithmetic"] = "REPRODUCES_TENSION (Re_k2 ACCEPT candidate)"
        elif abs(rek2["signed_median_dev_sigma"]) < 0.5:
            decision["verdict_arithmetic"] = "MCMC_CENTERS_ON_OBS (flow-offset => quarantine Re_k2)"
        else:
            decision["verdict_arithmetic"] = "INTERMEDIATE (scientific-reviewer adjudicates)"

    report = {
        "kind": "mcmc_rek2_pushforward",
        "note": ("Closes reviewer required-validation #1: pushes the pooled "
                 "repaired reference-MCMC posterior through the identical "
                 "theta->x forward loop (obs_noise=False) used by the SBI ppc "
                 "and training set. Importance-weighted pooled samples resampled "
                 "to unweighted (systematic, seeded) before pushforward. "
                 "DIAGNOSTIC; interpreted by scientific-reviewer, never tuned."),
        "comp": args.comp,
        "config": os.path.relpath(str(cfg_path), ROOT),
        "pooled_pkl": os.path.relpath(str(pooled_pkl), ROOT),
        "obs_names": obs_names, "x_obs": x_obs, "sigma": sigmas,
        "seed": int(seed),
        "n_pooled": int(samples.shape[0]),
        "n_push_requested": int(args.n_push),
        "n_push_forward_finite": int(n_forward),
        "per_channel": per_channel,
        "rek2_decision": decision,
    }
    out = outdir / "mcmc_rek2_pushforward_report.json"
    with open(out, "w") as f:
        json.dump(report, f, indent=2)

    print(f"\n[mcmc-push] Re_k2: obs={x_obs['Re_k2']:.4f} sigma={sigmas['Re_k2']:.4f}")
    print(f"[mcmc-push]   MCMC pp_median = {rek2['pp_median']:.4f} "
          f"(signed dev {rek2['signed_median_dev_sigma']:+.2f}sig, "
          f"|dev| {rek2['median_abs_dev_sigma']:.2f}sig)")
    if sbi_rek2:
        print(f"[mcmc-push]   SBI  pp_median = {sbi_rek2.get('pp_median'):.4f} "
              f"(|dev| {sbi_rek2.get('median_dev_sigma'):.2f}sig)")
        print(f"[mcmc-push]   MCMC-vs-SBI median gap = "
              f"{decision.get('mcmc_vs_sbi_median_gap_sigma'):.3f}sig "
              f"(within 0.5sig: {decision.get('mcmc_within_0p5sigma_of_sbi')})")
        print(f"[mcmc-push]   arithmetic verdict: {decision.get('verdict_arithmetic')}")
    print(f"[mcmc-push] report -> {out}")


if __name__ == "__main__":
    main()
