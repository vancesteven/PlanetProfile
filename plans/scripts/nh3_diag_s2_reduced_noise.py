"""NH3 under-update diagnosis — SEPARATOR S2: reduced-noise k2 pilot.

Manager decision (MACHINE-B-HANDOFF §0.8, 2026-08-04): run S1 + S2 before any
MgSO4/NaCl compute. This is S2 (diagnosis hypothesis #1, noise-augmentation
swamping, for the DOMINANT Im_k2 miss — the no-ocean control cleared #1 for
Re_k2 only, Im_k2 untested).

Question it separates. Training x pairs get Gaussian noise at sigma_obs; if that
noise swamps the k2 signal spread the flow rationally down-weights k2 and the
pushforward under-updates. S2 retrains on the SAME data with the k2 noise
REDUCED (sigma/4) and, as a second arm, ZERO k2 noise. Then:
  - Im_k2 gap shrinks materially -> the abs-fold/additive-noise convention
                                    contributes to the miss;
  - gap unchanged                -> #1 is closed for Im_k2 too (mechanism is the
                                    ocean-admitting apparatus, not the noise).

NO NEW FORWARD SIMS. The saved dataset stores the NOISED x, but the noise was
added as a single one-shot draw AFTER the |Im| fold, NOT re-folded:
    x_saved = x_clean + default_rng(noise_seed).normal(0,1,(N,4)) * sigmas
(mcmc_runner.generate_sbi_dataset ~L2579-2597; stats.obs_noise.refold_im=False).
So x_clean is recovered EXACTLY by reproducing the identical draw and subtracting
(verified: reconstructed |Im_k2| has min 0.0 and 0.0 frac<0 — the perfect
abs-fold signature; the saved x has ~22% negative Im_k2). We then re-noise the k2
channels at the reduced sigma and keep C20/C22 noise at the original level (the
gravity channels are not under study and assimilate cleanly).

Run:
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/nh3_diag_s2_reduced_noise.py
"""
import argparse
import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.sbi_runner import SBIRunner

ROOT = Path(__file__).resolve().parents[2]
CFG = ROOT / "PlanetProfile/Inference/configs/test54_titan_nh3_freegrav.json"
DATASET = Path("/tmp/titanG_build/datasets/titanG_nh3_1m.npz")
OUTDIR = ROOT / "validation_reports/nh3_diagnosis/s2_reduced_noise"
TMP_ART = Path("/tmp/nh3_diag/s2")
TRAIN_SEED = 72
K2_CHANNELS = ("Re_k2", "Im_k2")


def _recover_clean(x, stats, obs_names):
    """x_clean = x_saved - reproduced_noise. Requires diagonal gaussian noise
    (the NH3 dataset), refold_im=False."""
    on = stats.get("obs_noise")
    if on is None:
        raise SystemExit("[s2] dataset has no obs_noise metadata — cannot recover clean signal")
    if on.get("correlations"):
        raise SystemExit("[s2] correlated noise not supported by this recovery path")
    assert on.get("refold_im") is False, "recovery assumes noise NOT re-folded"
    sigmas = np.array([float(on["sigma"][n]) for n in obs_names])
    noise = np.random.default_rng(on["noise_seed"]).normal(0.0, 1.0, size=x.shape) * sigmas
    x_clean = x - noise
    # Validate the abs-fold signature on any |Im| channel.
    for j, n in enumerate(obs_names):
        if n.startswith("Im"):
            imk = x_clean[:, j]
            assert imk.min() >= -1e-9, f"{n} clean min {imk.min()} < 0 — recovery failed"
    return x_clean, sigmas


def _renoise(x_clean, sigmas, obs_names, k2_factor, seed):
    """Re-noise all channels at their original sigma, EXCEPT the k2 channels
    which are scaled by k2_factor (0.0 = zero-noise k2 arm)."""
    eff = sigmas.copy()
    for j, n in enumerate(obs_names):
        if n in K2_CHANNELS:
            eff[j] = sigmas[j] * k2_factor
    noise = np.random.default_rng(seed).normal(0.0, 1.0, size=x_clean.shape) * eff
    return x_clean + noise, eff


def _train_and_ppc(arm_label, x_arm, theta, runner_cfg, stats, args):
    tmp = TMP_ART / arm_label
    tmp.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / arm_label
    out.mkdir(parents=True, exist_ok=True)

    runner = SBIRunner(InferenceConfig.from_dict(runner_cfg))
    print(f"[s2:{arm_label}] training: theta={theta.shape} x={x_arm.shape} "
          f"seed={TRAIN_SEED} {args.density_estimator} max_epochs={args.max_epochs}")
    t0 = time.time()
    runner.train(theta, x_arm, seed=TRAIN_SEED,
                 density_estimator=args.density_estimator,
                 max_num_epochs=args.max_epochs)
    runner._train_info["rejection_stats"] = stats
    art = tmp / f"nh3_s2_{arm_label}_pilot.pt"
    runner.save_artifact(str(art))
    _ts = getattr(runner, "_last_train_summary", {}) or {}
    _ep = _ts.get("epochs_trained") or _ts.get("epochs")
    if isinstance(_ep, list):
        _ep = _ep[-1] if _ep else None
    print(f"[s2:{arm_label}] trained {(time.time()-t0)/60:.1f} min -> {art}; "
          f"epochs_trained={_ep} (early-stop fired if < {args.max_epochs})")

    # PPC uses the ARM's x so the prior-predictive matches what was trained.
    npz = tmp / f"nh3_s2_{arm_label}_dataset.npz"
    np.savez_compressed(npz, theta=theta, x=x_arm,
                        stats=json.dumps(stats, default=str))
    cmd = [
        "python", str(ROOT / "plans/scripts/titanG_ppc_interior_check.py"),
        "--artifact", str(art), "--config", str(CFG),
        "--dataset", str(npz), "--output-dir", str(out),
        "--n-post", str(args.n_post), "--seed", str(TRAIN_SEED),
    ]
    env = dict(os.environ, PYTHONPATH=str(ROOT),
               NUMBA_CACHE_DIR="/tmp/pp_numba_cache", KMP_DUPLICATE_LIB_OK="TRUE")
    r = subprocess.run(cmd, env=env, cwd=str(ROOT))
    print(f"[s2:{arm_label}] PPC exit {r.returncode}; report under {out}")
    return {"artifact": str(art), "ppc_report": str(out / "ppc_interior_report.json"),
            "epochs_trained": _ep}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k2-factor", type=float, default=0.25,
                    help="reduced-noise k2 sigma multiplier (sigma/4 default)")
    ap.add_argument("--zero-arm", action="store_true", default=True,
                    help="also train a zero-k2-noise arm")
    ap.add_argument("--no-zero-arm", dest="zero_arm", action="store_false")
    ap.add_argument("--max-epochs", type=int, default=60)
    ap.add_argument("--density-estimator", default="nsf")
    ap.add_argument("--n-post", type=int, default=4000)
    ap.add_argument("--renoise-seed", type=int, default=7272)
    args = ap.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    TMP_ART.mkdir(parents=True, exist_ok=True)
    if not DATASET.exists():
        raise SystemExit(f"[s2] missing dataset {DATASET}")

    with np.load(DATASET, allow_pickle=True) as d:
        theta = np.asarray(d["theta"], np.float64)
        x = np.asarray(d["x"], np.float64)
        stats = json.loads(str(d["stats"].item()))

    cfg = json.load(open(CFG)); cfg["mode"] = "sbi"
    obs_names = list(InferenceConfig.from_dict(cfg).observables.keys())
    x_clean, sigmas = _recover_clean(x, stats, obs_names)
    print(f"[s2] recovered clean signal; k2 channels {K2_CHANNELS}; "
          f"reduced factor={args.k2_factor}")

    arms = {}
    # Arm A: k2 noise reduced by factor (default sigma/4).
    x_red, eff_red = _renoise(x_clean, sigmas, obs_names, args.k2_factor, args.renoise_seed)
    arms[f"k2_sigma_over_{int(round(1/args.k2_factor))}"] = x_red
    # Arm B: zero k2 noise (optional).
    if args.zero_arm:
        x_zero, _ = _renoise(x_clean, sigmas, obs_names, 0.0, args.renoise_seed + 1)
        arms["k2_zero_noise"] = x_zero

    results = {}
    for label, x_arm in arms.items():
        results[label] = _train_and_ppc(label, x_arm, theta, cfg, stats, args)

    manifest = {
        "kind": "nh3_diag_s2_reduced_noise_pilot",
        "config": str(CFG), "dataset": str(DATASET),
        "orig_sigma": {n: float(s) for n, s in zip(obs_names, sigmas)},
        "k2_channels": list(K2_CHANNELS), "k2_factor": args.k2_factor,
        "arms": results, "train_seed": TRAIN_SEED,
        "density_estimator": args.density_estimator, "max_epochs": args.max_epochs,
        "clean_recovery": ("x_clean = x_saved - default_rng(noise_seed).normal*sigmas; "
                           "verified abs-fold signature (Im min>=0)."),
        "preregistered_reading": (
            "Im_k2 pushforward gap shrinks materially vs the CAPPED full-joint "
            "anchor (nh3_diag_capped_anchor.py, cap-vs-cap) -> noise/abs-fold "
            "convention contributes; unchanged -> #1 closed for Im_k2."),
        "anchor_note": ("compare against the capped full-joint anchor at the SAME "
                        "max_epochs/n_post/seed, NOT the historical converged "
                        "0.042; a 'gap unchanged' reading is only trustworthy if "
                        "the anchor also under-updates at this cap."),
    }
    with open(OUTDIR / "s2_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, default=str)
    print(f"[s2] manifest -> {OUTDIR/'s2_manifest.json'}")


if __name__ == "__main__":
    main()
