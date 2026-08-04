"""NH3 diagnosis — capped full-JOINT reference anchor (reviewer-required).

The scientific-reviewer flagged that comparing capped pilots (max_num_epochs=60)
against the full flow's converged SBI-pp |Im_k2|=0.042 is asymmetric: a
"gap unchanged" pilot reading could be under-training, not the mechanism. The fix
is a head-to-head anchor — retrain the FULL joint dataset (unfiltered, original
noise) at the SAME epoch cap, same architecture/seed, and PPC it identically. The
pilots then compare against THIS number (cap-vs-cap), not the historical 0.042.

If this capped full-joint anchor already under-updates like the deployed flow
(|Im_k2| ~0.04), the cap is not the cause and the pilots' "unchanged" readings
are trustworthy. If the capped anchor assimilates BETTER than 0.042, the cap
matters and the pilots need convergence (raise --max-epochs).

NO NEW FORWARD SIMS — reuses /tmp/titanG_build/datasets/titanG_nh3_1m.npz as-is.

Run:
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/nh3_diag_capped_anchor.py
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
OUTDIR = ROOT / "validation_reports/nh3_diagnosis/capped_full_joint_anchor"
TMP_ART = Path("/tmp/nh3_diag/anchor")
TRAIN_SEED = 72


def _epochs_from(runner):
    """Best-effort (epochs_trained, best_val_logprob) from the train summary."""
    ts = getattr(runner, "_last_train_summary", {}) or {}
    epochs = ts.get("epochs_trained") or ts.get("epochs")
    if isinstance(epochs, list):
        epochs = epochs[-1] if epochs else None
    bvl = ts.get("best_validation_log_prob") or ts.get("best_validation_log_probs")
    if isinstance(bvl, list):
        bvl = bvl[-1] if bvl else None
    return epochs, bvl, ts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-epochs", type=int, default=60)
    ap.add_argument("--density-estimator", default="nsf")
    ap.add_argument("--n-post", type=int, default=4000)
    args = ap.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    TMP_ART.mkdir(parents=True, exist_ok=True)
    if not DATASET.exists():
        raise SystemExit(f"[anchor] missing dataset {DATASET}")

    with np.load(DATASET, allow_pickle=True) as d:
        theta = np.asarray(d["theta"], np.float64)
        x = np.asarray(d["x"], np.float64)
        stats = json.loads(str(d["stats"].item()))

    cfg = json.load(open(CFG)); cfg["mode"] = "sbi"
    runner = SBIRunner(InferenceConfig.from_dict(cfg))
    print(f"[anchor] capped full-joint: theta={theta.shape} x={x.shape} "
          f"seed={TRAIN_SEED} {args.density_estimator} max_epochs={args.max_epochs}")
    t0 = time.time()
    runner.train(theta, x, seed=TRAIN_SEED,
                 density_estimator=args.density_estimator,
                 max_num_epochs=args.max_epochs)
    runner._train_info["rejection_stats"] = stats
    epochs, bvl, ts = _epochs_from(runner)
    art = TMP_ART / "nh3_capped_full_joint_anchor.pt"
    runner.save_artifact(str(art))
    print(f"[anchor] trained {(time.time()-t0)/60:.1f} min; epochs_trained={epochs} "
          f"best_val_logprob={bvl} (early-stop fired if < {args.max_epochs})")

    manifest = {
        "kind": "nh3_diag_capped_full_joint_anchor",
        "artifact": str(art), "config": str(CFG),
        "max_epochs": args.max_epochs, "epochs_trained": epochs,
        "best_validation_log_prob": bvl, "train_summary_keys": list(ts.keys()),
        "train_seed": TRAIN_SEED, "density_estimator": args.density_estimator,
        "n_train": int(len(theta)), "config_hash": runner.config.generate_hash(),
        "purpose": ("head-to-head cap-vs-cap anchor for S1/S2 pilots; if this "
                    "under-updates like the deployed flow (|Im_k2|~0.04) the cap "
                    "is not the cause and pilot 'unchanged' readings are trusted."),
    }
    with open(OUTDIR / "anchor_train_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, default=str)

    print("[anchor] launching PPC subprocess ...")
    cmd = [
        "python", str(ROOT / "plans/scripts/titanG_ppc_interior_check.py"),
        "--artifact", str(art), "--config", str(CFG),
        "--dataset", str(DATASET), "--output-dir", str(OUTDIR),
        "--n-post", str(args.n_post), "--seed", str(TRAIN_SEED),
    ]
    env = dict(os.environ, PYTHONPATH=str(ROOT),
               NUMBA_CACHE_DIR="/tmp/pp_numba_cache", KMP_DUPLICATE_LIB_OK="TRUE")
    r = subprocess.run(cmd, env=env, cwd=str(ROOT))
    print(f"[anchor] PPC exit {r.returncode}; report + manifest under {OUTDIR}")


if __name__ == "__main__":
    main()
