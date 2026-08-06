"""NH3 under-update diagnosis — #4 ARCHITECTURE PILOT (capacity / embedding).

Machine A §0.10 item 2 (PRODUCTION AUTHORIZATION, 2026-08-06). Runs in parallel
with the v5/v6 priority-(a) gates; reuses the EXISTING 1M NH3 dataset (NO new
forward sims). PILOT ONLY — the production-architecture decision is
manager + reviewer + user, made WITH this pilot's four-way table in hand. Do NOT
fold any arm into an MgSO4/NaCl production build from this script.

Question
--------
The NH3 diagnosis eliminated observation-vector width, x-normalization, the
noise convention, the joint no-ocean mixture, the salinity axis + Tb-w
degeneracy (#1), training-set size (#1 matched-N control), and reference-MCMC
resolution (#2/B3). The mechanism is concentration-FAILURE (#3): the deployed
flow's tidal posterior sits near the prior. The last candidate is flow
CAPACITY / EMBEDDING for this problem class (13 params, strong Tb-w degeneracy,
weakly identified tidal sector). #4 tests whether a larger flow and/or an
x-embedding closes the |Im_k2| under-update.

Frozen design — scientific-reviewer PASS WITH CONCERNS (2026-08-06)
------------------------------------------------------------------
Arms (all on the FULL 689,845-row ocean-admitting dataset; the deployed flow's
footing — comparability to the DEPLOYED flow is what matters here, NOT to the
S1/#1 ocean-only pilots):
  D0 : deployed-arch nsf h50/t5/b10, NO embedding  (converged control /
       no-regression check — MUST reproduce ~0.042; the harness extension must
       not perturb the deployed architecture)
  A  : capacity-up nsf h128/t10/b16, no embedding
  B  : deployed nsf h50/t5/b10 + embedding MLP (4->64->64)
  C  : capacity-up h128/t10/b16 + embedding MLP (4->64->64)

CRITICAL fixes folded in:
- The deployed flow trained UNCAPPED (titanG_nh3_train_all.py: no max_num_epochs
  -> early-stopping convergence). So the 0.043 CAPPED anchor is the WRONG
  reference for these uncapped arms; the converged deployed flow (~0.042) is.
  D0 provides that reference IN THIS SCRIPT at matched seed/harness, removing
  harness/seed drift. Read A/B/C against D0, not against 0.043.
- All arms UNCAPPED with the DEPLOYED early-stopping criterion
  (stop_after_epochs=20, validation_fraction=0.1) held fixed; only a
  non-binding safety ceiling max_num_epochs=500 (FLAG any arm that hits it —
  its convergence claim is then void).
- z_score_theta = z_score_x = 'independent' held fixed across ALL arms so
  capacity/embedding is cleanly isolated (the dim-3 Im_k2 z-scoring outlier
  persists in every arm; addressing it is a SEPARATE diagnostic, out of scope).
- >=3 seeds per arm (72, 172, 272). A single-seed PASS is only a screen; any arm
  meeting the target must be confirmed at >=3 seeds BEFORE it enters the
  manager+reviewer+user production decision (target shift 0.042->0.086 ~ +1.26
  sigma_obs is within plausible NPE seed variance).

Pre-registered JOINT success criterion (no single metric alone is a PASS):
  (a) SBI-pp |Im_k2| median >= 0.0862  (= 0.1037 - 0.5*0.035; primary), AND
  (b) concentration_ratio < 1 (pp narrower than prior; #3 found 1.215 > 1), AND
  (c) frac(pushforward Im_k2 >= obs) materially increased from the prior
      baseline toward ~0.5.
  OVERSHOOT flag: SBI-pp median > 0.121 (= 0.1037 + 0.5*0.035) signals
  over-update / mode-collapse, NOT success.
Per arm/seed also report: epochs_trained (+ hit-ceiling flag), best validation
log-prob, train-vs-validation gap, and the pushforward nonfinite/support-drop
fraction (a median shift driven by more mass in support-rejected regions is an
artifact, not assimilation).

Reference anchors (provenance; NOT recomputed): capped anchor 0.043 (historical,
wrong reference for uncapped arms — kept for the table only), converged deployed
~0.042 (regenerated here as D0), matched MCMC-pp ceiling 0.1037 (#2), obs 0.135,
sigma_obs 0.035.

libomp: cache/dataset read + torch train in ONE process (numba forward model NOT
imported here); each PPC step is a SEPARATE subprocess (imports the forward
model). No new forward sims.

Run:
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/nh3_diag_4_architecture_pilot.py
Single-seed screen first (then confirm passers at >=3 seeds):
  ... nh3_diag_4_architecture_pilot.py --seeds 72
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
OUTDIR = ROOT / "validation_reports/nh3_diagnosis/f4_architecture"
TMP_ART = Path("/tmp/nh3_diag/f4")

# Deployed provenance (verified from the artifact): nsf, seed 72,
# config_hash e596574d1e81567c, trained UNCAPPED.
DEPLOYED_SEED = 72
SEEDS_DEFAULT = [72, 172, 272]
MAX_EPOCHS_CEILING = 500          # non-binding safety ceiling; FLAG if hit
STOP_AFTER_EPOCHS = 20            # deployed early-stopping criterion (fixed)
VALIDATION_FRACTION = 0.1         # deployed (fixed)

# Comparison anchors (provenance; NOT recomputed here).
CAPPED_ANCHOR_IMK2 = 0.043        # historical 60-epoch-capped full-joint anchor
DEPLOYED_CONVERGED_IMK2 = 0.042   # deployed uncapped flow (D0 must reproduce ~this)
MATCHED_MCMC_PP_CEILING = 0.1037  # #2 matched-resolution (n_eff=2000) ceiling
OBS_IMK2 = 0.135
SIGMA_IMK2 = 0.035
TARGET_IMK2 = MATCHED_MCMC_PP_CEILING - 0.5 * SIGMA_IMK2   # 0.0862 (primary)
OVERSHOOT_IMK2 = MATCHED_MCMC_PP_CEILING + 0.5 * SIGMA_IMK2  # 0.121

# Frozen arms. embedding_net None -> Identity (deployed); "mlp:64-64" -> FCEmbedding.
ARMS = [
    {"tag": "D0", "hidden_features": 50, "num_transforms": 5, "num_bins": 10,
     "embedding_net": None,
     "role": "converged deployed-arch control / no-regression check"},
    {"tag": "A", "hidden_features": 128, "num_transforms": 10, "num_bins": 16,
     "embedding_net": None, "role": "capacity-up"},
    {"tag": "B", "hidden_features": 50, "num_transforms": 5, "num_bins": 10,
     "embedding_net": "mlp:64-64", "role": "embedding-only"},
    {"tag": "C", "hidden_features": 128, "num_transforms": 10, "num_bins": 16,
     "embedding_net": "mlp:64-64", "role": "capacity + embedding"},
]


def _train_one(arm, seed, theta, x, stats, runner_cfg, args):
    """Train one (arm, seed) flow uncapped with the deployed stopping rule,
    save the artifact, and run PPC in a separate subprocess. Fresh runner per
    call so arms/seeds never share fit state."""
    tag = f"{arm['tag']}_seed{seed}"
    runner = SBIRunner(InferenceConfig.from_dict(runner_cfg))
    print(f"[f4:{tag}] training: theta={theta.shape} x={x.shape} seed={seed} "
          f"nsf h{arm['hidden_features']}/t{arm['num_transforms']}/"
          f"b{arm['num_bins']} embed={arm['embedding_net']}")
    t1 = time.time()
    runner.train(
        theta, x, seed=seed, density_estimator="nsf",
        hidden_features=arm["hidden_features"],
        num_transforms=arm["num_transforms"],
        num_bins=arm["num_bins"],
        embedding_net=arm["embedding_net"],
        stop_after_epochs=STOP_AFTER_EPOCHS,
        validation_fraction=VALIDATION_FRACTION,
        max_num_epochs=args.max_epochs,
    )
    runner._train_info["rejection_stats"] = stats
    runner._train_info["f4_arm"] = arm
    art = TMP_ART / f"nh3_f4_{tag}.pt"
    runner.save_artifact(str(art))

    # Convergence + overfit diagnostics from the sbi training summary.
    ts = getattr(runner, "_last_train_summary", {}) or {}
    ep = ts.get("epochs_trained") or ts.get("epochs")
    if isinstance(ep, list):
        ep = ep[-1] if ep else None
    hit_ceiling = (ep is not None) and (ep >= args.max_epochs)
    # best validation log-prob and train-vs-val gap (schema varies; be defensive)
    val_curve = ts.get("validation_log_probs") or ts.get("validation_log_prob")
    train_curve = ts.get("training_log_probs") or ts.get("training_log_prob")
    best_val = ts.get("best_validation_log_prob")
    if best_val is None and isinstance(val_curve, list) and val_curve:
        best_val = max(val_curve)
    train_val_gap = None
    if (isinstance(train_curve, list) and isinstance(val_curve, list)
            and train_curve and val_curve):
        n = min(len(train_curve), len(val_curve))
        train_val_gap = float(train_curve[n - 1] - val_curve[n - 1])
    print(f"[f4:{tag}] trained in {(time.time()-t1)/60:.1f} min -> {art}; "
          f"epochs={ep} hit_ceiling={hit_ceiling} best_val={best_val} "
          f"train_val_gap={train_val_gap}")
    if hit_ceiling:
        print(f"[f4:{tag}] WARNING: hit max_num_epochs ceiling -> "
              f"convergence claim VOID for this arm/seed")

    out = OUTDIR / tag
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, PYTHONPATH=str(ROOT),
               NUMBA_CACHE_DIR="/tmp/pp_numba_cache", KMP_DUPLICATE_LIB_OK="TRUE")

    # (1) Pushforward-shape diagnostic — the SOURCE of the three joint-criterion
    #     inputs (pp-median, concentration ratio vs prior, frac-mass>=obs and its
    #     prior baseline) plus the nonfinite drop-bias fraction. Same script as #3.
    print(f"[f4:{tag}] launching pushforward-shape subprocess ...")
    pf_cmd = [
        "python", str(ROOT / "plans/scripts/nh3_diag_pushforward_plot.py"),
        "--artifact", str(art), "--config", str(CFG), "--dataset", str(DATASET),
        "--n-post", str(args.n_post), "--seed", str(seed),
        "--outdir", str(out / "pushforward"), "--label", tag,
    ]
    rpf = subprocess.run(pf_cmd, env=env, cwd=str(ROOT))
    print(f"[f4:{tag}] pushforward exit {rpf.returncode}")

    # (2) Full-interior PPC (all channels) for the standard interior report.
    print(f"[f4:{tag}] launching interior PPC subprocess ...")
    ppc_cmd = [
        "python", str(ROOT / "plans/scripts/titanG_ppc_interior_check.py"),
        "--artifact", str(art), "--config", str(CFG),
        "--dataset", str(DATASET), "--output-dir", str(out),
        "--n-post", str(args.n_post), "--seed", str(seed),
    ]
    rppc = subprocess.run(ppc_cmd, env=env, cwd=str(ROOT))
    print(f"[f4:{tag}] interior PPC exit {rppc.returncode}; reports under {out}")

    # Read verdict statistics from the pushforward report (correct key names).
    pp_imk2 = concentration = frac_ge_obs = prior_frac_ge_obs = drop_frac = None
    pf_rep = out / "pushforward" / "pushforward_shape_report.json"
    if pf_rep.exists():
        try:
            pr = json.load(open(pf_rep))
            pp_imk2 = pr.get("sbi_posterior_predictive_spread", {}).get("median")
            concentration = pr.get("sbi_concentration_ratio_vs_prior")
            frac_ge_obs = pr.get("sbi_frac_mass_ge_obs")
            prior_frac_ge_obs = pr.get("prior_frac_mass_ge_obs")
            drop_frac = pr.get("nonfinite_drop_bias_check", {}).get("drop_frac")
        except Exception as e:
            print(f"[f4:{tag}] WARNING: could not parse pushforward report: {e}")

    return {
        "tag": tag, "arm": arm["tag"], "seed": seed,
        "artifact": str(art), "output_dir": str(out),
        "epochs_trained": ep, "hit_ceiling": bool(hit_ceiling),
        "best_validation_log_prob": best_val, "train_val_gap": train_val_gap,
        "architecture_tag": runner._train_info.get("density_estimator"),
        "config_hash": runner.config.generate_hash(),
        "pp_imk2_median": pp_imk2, "concentration_ratio": concentration,
        "frac_ge_obs": frac_ge_obs, "prior_frac_ge_obs": prior_frac_ge_obs,
        "pushforward_drop_fraction": drop_frac,
    }


def _verdict(r):
    """Pre-registered JOINT success verdict for one (arm, seed) run."""
    m = r.get("pp_imk2_median")
    c = r.get("concentration_ratio")
    if m is None:
        return {"status": "PPC_MISSING"}
    overshoot = m > OVERSHOOT_IMK2
    passes = (m >= TARGET_IMK2) and (c is not None and c < 1.0) and (not overshoot)
    return {
        "meets_primary_median": bool(m >= TARGET_IMK2),
        "concentrates": bool(c is not None and c < 1.0),
        "overshoot": bool(overshoot),
        "joint_pass": bool(passes),
        "median_dev_sigma_from_ceiling": (m - MATCHED_MCMC_PP_CEILING) / SIGMA_IMK2,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=SEEDS_DEFAULT,
                    help="training seeds per arm (>=3 for a decision; 1 = screen)")
    ap.add_argument("--arms", nargs="+", default=[a["tag"] for a in ARMS],
                    help="subset of arm tags to run (D0 A B C)")
    ap.add_argument("--max-epochs", type=int, default=MAX_EPOCHS_CEILING,
                    help="non-binding safety ceiling; arms flagged if hit")
    ap.add_argument("--n-post", type=int, default=4000)
    args = ap.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    TMP_ART.mkdir(parents=True, exist_ok=True)
    if not DATASET.exists():
        raise SystemExit(f"[f4] missing dataset {DATASET}")

    with np.load(DATASET, allow_pickle=True) as d:
        theta = np.asarray(d["theta"], np.float64)
        x = np.asarray(d["x"], np.float64)
        stats = json.loads(str(d["stats"].item())) if "stats" in d.files else {}
    print(f"[f4] dataset: theta={theta.shape} x={x.shape} (FULL ocean-admitting)")

    cfg = json.load(open(CFG)); cfg["mode"] = "sbi"
    probe = SBIRunner(InferenceConfig.from_dict(cfg))
    obs_names = list(probe.obs_names)
    assert theta.shape[1] == len(probe.param_names), "theta col mismatch"

    arms = [a for a in ARMS if a["tag"] in set(args.arms)]
    screen = len(args.seeds) < 3

    runs = []
    for arm in arms:
        for seed in args.seeds:
            runs.append(_train_one(arm, seed, theta, x, stats, cfg, args))

    # Attach per-run verdicts and aggregate median-of-seeds per arm.
    for r in runs:
        r["verdict"] = _verdict(r)
    per_arm = {}
    for arm in arms:
        rs = [r for r in runs if r["arm"] == arm["tag"]
              and r.get("pp_imk2_median") is not None]
        meds = [r["pp_imk2_median"] for r in rs]
        per_arm[arm["tag"]] = {
            "n_seeds": len(rs),
            "pp_imk2_median_of_seeds": float(np.median(meds)) if meds else None,
            "pp_imk2_min": float(np.min(meds)) if meds else None,
            "pp_imk2_max": float(np.max(meds)) if meds else None,
            "pp_imk2_spread": float(np.max(meds) - np.min(meds)) if len(meds) > 1 else None,
            "any_seed_joint_pass": any(r["verdict"].get("joint_pass") for r in rs),
            "all_seeds_joint_pass": bool(rs) and all(r["verdict"].get("joint_pass") for r in rs),
        }

    manifest = {
        "kind": "nh3_diag_f4_architecture_pilot",
        "config": str(CFG), "dataset": str(DATASET),
        "note": ("PILOT ONLY. Production-architecture decision is "
                 "manager+reviewer+user WITH this table. Do NOT fold any arm "
                 "into an MgSO4/NaCl production build from this script."),
        "footing": "FULL ocean-admitting 1M NH3 dataset (deployed flow's footing)",
        "seeds": args.seeds, "is_single_seed_screen": screen,
        "stop_after_epochs": STOP_AFTER_EPOCHS,
        "validation_fraction": VALIDATION_FRACTION,
        "max_epochs_ceiling": args.max_epochs,
        "z_score_theta": "independent", "z_score_x": "independent",
        "obs_names": obs_names,
        "arms": arms,
        "comparison_anchors": {
            "capped_anchor_imk2_WRONG_ref_for_uncapped": CAPPED_ANCHOR_IMK2,
            "deployed_converged_imk2_D0_target": DEPLOYED_CONVERGED_IMK2,
            "matched_mcmc_pp_ceiling": MATCHED_MCMC_PP_CEILING,
            "obs_imk2": OBS_IMK2, "sigma_imk2": SIGMA_IMK2,
            "primary_target_imk2": TARGET_IMK2,
            "overshoot_flag_imk2": OVERSHOOT_IMK2,
        },
        "success_criterion": (
            "JOINT: pp-median |Im_k2| >= 0.0862 AND concentration_ratio < 1 AND "
            "frac_ge_obs materially increased vs prior; overshoot > 0.121 flags "
            "over-update. D0 must reproduce ~0.042 (no-regression check on the "
            "harness extension). Read A/B/C against D0, not against 0.043. A "
            "single-seed PASS is a SCREEN only; confirm at >=3 seeds before the "
            "production decision."),
        "reviewer": ("PASS WITH CONCERNS 2026-08-06: additive train() capacity/"
                     "embedding plumbing (CRITICAL); D0 converged control + read "
                     "vs D0 not 0.043 (MAJOR); >=3 seeds before decision (MAJOR); "
                     "joint criterion not median-alone (MAJOR); z-scoring fixed "
                     "'independent' across arms (MODERATE); epochs/best-val/"
                     "train-val-gap/drop-fraction per arm (MODERATE)."),
        "runs": runs,
        "per_arm_summary": per_arm,
    }
    with open(OUTDIR / "f4_pilot_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, default=str)
    print(f"[f4] manifest -> {OUTDIR/'f4_pilot_manifest.json'}")
    print(f"[f4] arms={[a['tag'] for a in arms]} seeds={args.seeds} "
          f"screen={screen}")
    for tag, s in per_arm.items():
        print(f"[f4] {tag}: pp_imk2 med-of-seeds={s['pp_imk2_median_of_seeds']} "
              f"(spread {s['pp_imk2_spread']}) any_pass={s['any_seed_joint_pass']}")


if __name__ == "__main__":
    main()
