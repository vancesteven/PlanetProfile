"""NH3 under-update diagnosis — SEPARATOR S1: ocean-only-with-salinity pilot.

Manager decision (MACHINE-B-HANDOFF §0.8, 2026-08-04): run S1 + S2 before any
MgSO4/NaCl compute. This is S1.

Question it separates. The full NH3 JOINT flow under-updates k2 (SBI-pp |Im_k2|
median 0.042 vs MCMC-pp 0.100, obs 0.135). The diagnosis localized the mechanism
to the "ocean-admitting apparatus" — the joint no-ocean+ocean MIXTURE *and/or*
the co-varying salinity axis — but could not separate the two (they are bundled
in the no-ocean control). S1 breaks the mixture out: train a pilot flow on the
OCEAN-BRANCH rows only, still carrying the log10_wOcean_ppt salinity axis. Then:
  - clean k2 assimilation  -> the joint MIXTURE (bimodal conditional) is the driver;
  - persistent k2 miss     -> the salinity degeneracy (or another shared element)
                              is implicated even without the frozen branch.

NO NEW FORWARD SIMS. Reuses the existing 1M NH3 dataset
(/tmp/titanG_build/datasets/titanG_nh3_1m.npz).

has_ocean recovery (dataset carries NO per-row tag). For each row we look up its
(Tb_K = theta[:,7], log10_wOcean_ppt = theta[:,8]) in the SAME NH3 2D cache the
forward model used, via the SAME bilinear operator (grid_interp_2d), and read the
DOMINANT-WEIGHT (nearest) node's has_ocean tag. A row is "ocean" iff its nearest
grid node has_ocean=True. This mirrors how the forward model selected structure
(nearest/blend over the 12x12 grid) — so the filter matches the trained support.

libomp hazard: cache read (numpy/pickle) + torch training happen in ONE process
here but PlanetProfile's numba forward model is NOT imported (we only read the
cache dict and train torch). Env still pins KMP_DUPLICATE_LIB_OK=TRUE. The PPC
step is a SEPARATE subprocess (it imports the forward model).

Run:
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/nh3_diag_s1_ocean_only.py
"""
import argparse
import json
import os
import pickle
import subprocess
import time
from pathlib import Path

import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.sbi_runner import SBIRunner
from PlanetProfile.Inference import grid_interp_2d as g2

ROOT = Path(__file__).resolve().parents[2]
CFG = ROOT / "PlanetProfile/Inference/configs/test54_titan_nh3_freegrav.json"
DATASET = Path("/tmp/titanG_build/datasets/titanG_nh3_1m.npz")
CACHE = (ROOT / "PlanetProfile/Test/mcmc_results/Titan/Test54_nh3_ocean/"
         "titan_nh3_joint_structure_grid_2d.pkl")
OUTDIR = ROOT / "validation_reports/nh3_diagnosis/s1_ocean_only"
TMP_ART = Path("/tmp/nh3_diag/s1")
TRAIN_SEED = 72


def _has_ocean_grid(cache):
    """(n_Tb*n_w,) bool array: does each flat grid node have an ocean?"""
    S = cache["structures"]
    ho = np.zeros(len(S), dtype=bool)
    for i, node in enumerate(S):
        ho[i] = bool(isinstance(node, dict) and node.get("has_ocean"))
    return ho


def _row_is_ocean(theta, param_names, cache):
    """Per-row ocean-branch membership from the nearest grid node.

    Uses grid_interp_2d.bilinear_weights (the forward model's operator) and takes
    the dominant-weight corner's has_ocean tag — the node the nearest/blend path
    would key on."""
    Tb_grid = np.asarray(cache["Tb_K_grid"], float)
    w_grid = np.asarray(cache["wOcean_ppt_grid"], float)
    ho = _has_ocean_grid(cache)
    i_tb = param_names.index("Tb_K")
    i_lw = param_names.index("log10_wOcean_ppt")
    out = np.zeros(len(theta), dtype=bool)
    for r in range(len(theta)):
        Tb = float(theta[r, i_tb])
        w_ppt = float(10.0 ** theta[r, i_lw])
        corners, weights = g2.bilinear_weights(Tb_grid, w_grid, Tb, w_ppt)
        dom = corners[int(np.argmax(weights))]
        out[r] = bool(ho[dom])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-epochs", type=int, default=60,
                    help="cap NPE epochs (pilot; full run is uncapped)")
    ap.add_argument("--density-estimator", default="nsf")
    ap.add_argument("--n-post", type=int, default=4000)
    args = ap.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    TMP_ART.mkdir(parents=True, exist_ok=True)
    if not DATASET.exists():
        raise SystemExit(f"[s1] missing dataset {DATASET}")

    with np.load(DATASET, allow_pickle=True) as d:
        theta = np.asarray(d["theta"], np.float64)
        x = np.asarray(d["x"], np.float64)
        stats = json.loads(str(d["stats"].item())) if "stats" in d.files else {}

    cfg = json.load(open(CFG)); cfg["mode"] = "sbi"
    runner = SBIRunner(InferenceConfig.from_dict(cfg))
    param_names = list(runner.param_names)
    obs_names = list(runner.obs_names)
    assert theta.shape[1] == len(param_names), "theta col mismatch"

    cache = pickle.load(open(CACHE, "rb"))
    print(f"[s1] recovering has_ocean for {len(theta):,} rows via nearest grid node ...")
    t0 = time.time()
    is_ocean = _row_is_ocean(theta, param_names, cache)
    n_ocean = int(is_ocean.sum())
    print(f"[s1] ocean rows = {n_ocean:,} / {len(theta):,} "
          f"({100.0*n_ocean/len(theta):.1f}%); recover {time.time()-t0:.1f}s")

    theta_o = theta[is_ocean]
    x_o = x[is_ocean]

    print(f"[s1] training ocean-only pilot: theta={theta_o.shape} x={x_o.shape} "
          f"seed={TRAIN_SEED} {args.density_estimator} max_epochs={args.max_epochs}")
    t1 = time.time()
    runner.train(theta_o, x_o, seed=TRAIN_SEED,
                 density_estimator=args.density_estimator,
                 max_num_epochs=args.max_epochs)
    runner._train_info["rejection_stats"] = stats
    runner._train_info["s1_filter"] = {
        "kind": "ocean_branch_only", "n_ocean": n_ocean,
        "n_total_dataset": int(len(theta)),
        "has_ocean_recovery": "nearest-grid-node via grid_interp_2d.bilinear_weights",
        "cache": str(CACHE),
    }
    art = TMP_ART / "nh3_s1_ocean_only_pilot.pt"
    runner.save_artifact(str(art))
    _ts = getattr(runner, "_last_train_summary", {}) or {}
    _ep = _ts.get("epochs_trained") or _ts.get("epochs")
    if isinstance(_ep, list):
        _ep = _ep[-1] if _ep else None
    print(f"[s1] trained in {(time.time()-t1)/60:.1f} min -> {art}; "
          f"epochs_trained={_ep} (early-stop fired if < {args.max_epochs})")

    # Persist the filtered dataset so the PPC subprocess uses the pilot's
    # actual prior-predictive (ocean-branch x), not the joint dataset.
    filt_npz = TMP_ART / "nh3_s1_ocean_only_dataset.npz"
    np.savez_compressed(filt_npz, theta=theta_o, x=x_o,
                        stats=json.dumps(stats, default=str))

    manifest = {
        "kind": "nh3_diag_s1_ocean_only_pilot",
        "artifact": str(art), "config": str(CFG), "cache": str(CACHE),
        "n_ocean": n_ocean, "n_total": int(len(theta)),
        "train_seed": TRAIN_SEED, "density_estimator": args.density_estimator,
        "max_epochs": args.max_epochs, "epochs_trained": _ep,
        "obs_names": obs_names,
        "config_hash": runner.config.generate_hash(),
        "anchor_note": ("compare SBI-pp |Im_k2| against the CAPPED full-joint "
                        "anchor (nh3_diag_capped_anchor.py), NOT the historical "
                        "converged 0.042 — cap-vs-cap head-to-head."),
        "preregistered_reading": (
            "clean k2 assimilation -> joint MIXTURE is the driver; persistent "
            "k2 miss -> salinity degeneracy (or shared element) implicated."),
    }
    with open(OUTDIR / "s1_train_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, default=str)

    # PPC in a separate process (imports the forward model).
    print("[s1] launching PPC subprocess at Titan datum ...")
    cmd = [
        "python", str(ROOT / "plans/scripts/titanG_ppc_interior_check.py"),
        "--artifact", str(art), "--config", str(CFG),
        "--dataset", str(filt_npz), "--output-dir", str(OUTDIR),
        "--n-post", str(args.n_post), "--seed", str(TRAIN_SEED),
    ]
    env = dict(os.environ, PYTHONPATH=str(ROOT),
               NUMBA_CACHE_DIR="/tmp/pp_numba_cache", KMP_DUPLICATE_LIB_OK="TRUE")
    r = subprocess.run(cmd, env=env, cwd=str(ROOT))
    print(f"[s1] PPC exit code {r.returncode}; report under {OUTDIR}")
    print(f"[s1] manifest -> {OUTDIR/'s1_train_manifest.json'}")


if __name__ == "__main__":
    main()
