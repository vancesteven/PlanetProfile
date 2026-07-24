"""Diagnostic retrain of the Titan no-ocean nsf flow: capture the NLL curve
and test the reviewer's under-training hypothesis (crosscheck a2d7f121 step 3).

The cleaned 1M manifold (k2 guard) makes SBC pass but the crosscheck still
shows a UNIFORM ~1.17x over-dispersion on all 12 params (worst on dC20_nh/
dC22_nh). Reviewer step 3: "if over-dispersion persists after cleaning, the
cause is under-training, not data -- increase epochs / check the loss curve."

This trains once with a LONGER early-stop patience (stop_after_epochs=40 vs
sbi default 20) and DUMPS sbi's training/validation log-prob curve so we can
see whether validation NLL was still improving at the default stop (=> genuine
under-training) or had plateaued (=> the reference MCMC is the suspect). Writes
a SEPARATE diagnostic artifact; does NOT overwrite the production .pt until the
re-gate confirms improvement.

Run (separate process; libomp):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/titanG_train_diag.py
"""
import json, os, time, shutil
import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.sbi_runner import SBIRunner

DSDIR = "/tmp/titanG_build/datasets"
ARTDIR = "/tmp/titanG_build/artifacts"
CFG = "PlanetProfile/Inference/configs/titan_freegrav_noocean.json"
ARTNAME = "titan_freegrav_noocean_posterior_1m_diag.pt"
TRAIN_SEED = 71
STOP_AFTER_EPOCHS = 40   # longer patience than sbi default (20)


def main():
    os.makedirs(ARTDIR, exist_ok=True)
    npz = os.path.join(DSDIR, "titanG_noocean_1m.npz")
    with np.load(npz, allow_pickle=True) as d:
        theta = np.asarray(d["theta"], np.float64)
        x = np.asarray(d["x"], np.float64)

    cfg = json.load(open(CFG)); cfg["mode"] = "sbi"
    runner = SBIRunner(InferenceConfig.from_dict(cfg))
    print(f"[diag] theta={theta.shape} x={x.shape} seed={TRAIN_SEED} nsf "
          f"stop_after_epochs={STOP_AFTER_EPOCHS}")
    t0 = time.time()
    runner.train(theta, x, seed=TRAIN_SEED, density_estimator="nsf",
                 stop_after_epochs=STOP_AFTER_EPOCHS)
    dt = time.time() - t0

    summary = getattr(runner, "_last_train_summary", {}) or {}
    # sbi summary keys vary; print whatever curve-like entries exist.
    print(f"[diag] trained in {dt/60:.1f} min; summary keys: {list(summary.keys())}")
    for key in ("training_log_probs", "validation_log_probs",
                "best_validation_log_prob", "epochs_trained", "epochs"):
        if key in summary:
            v = summary[key]
            if isinstance(v, list):
                n = len(v)
                tail = ", ".join(f"{e:.4f}" for e in v[-8:])
                print(f"[diag]   {key}: n={n}, last8=[{tail}]")
                if key == "validation_log_probs" and n >= 22:
                    # Was val log-prob still improving at the sbi default stop
                    # (epoch = best + 20)? Compare best epoch to n.
                    best_ep = int(np.argmax(v))
                    print(f"[diag]   best val-logprob at epoch {best_ep}/{n-1}; "
                          f"default-stop would trigger at ~{best_ep+20}; "
                          f"{'STILL IMPROVING past default stop' if best_ep+20 < n-1 else 'plateaued by default stop'}")
            else:
                print(f"[diag]   {key}: {v}")

    out = os.path.join(ARTDIR, ARTNAME)
    runner.save_artifact(out)
    dropbox = os.path.join("PlanetProfile/Inference/sbi_artifacts", ARTNAME)
    shutil.copy2(out, dropbox)

    with open(os.path.join(ARTDIR, "titanG_train_diag_summary.json"), "w") as f:
        json.dump({"stop_after_epochs": STOP_AFTER_EPOCHS,
                   "train_seconds": dt, "n_train": int(len(theta)),
                   "artifact": dropbox, "train_summary": summary}, f, indent=2)
    print(f"[diag] diagnostic artifact -> {dropbox}")


if __name__ == "__main__":
    main()
