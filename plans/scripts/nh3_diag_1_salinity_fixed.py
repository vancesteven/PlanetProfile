"""NH3 under-update diagnosis — FOLLOW-UP #1: salinity-fixed ocean-only pilot
(+ a sample-size-matched, salinity-VARYING control).

Manager AUTHORIZED (MACHINE-B-HANDOFF §0.9, 2026-08-06). #2's outcome satisfied
the preregistered branch (SBI-pp-vs-matched-MCMC-pp gap SURVIVED at +1.76 sigma_obs
>> 0.5), so #1 runs. It is a PILOT retrain on the EXISTING 1M NH3 dataset with NO
artifact-design change and NO new forward sims — no further sign-off needed to RUN
it (reviewer + user sign-off attaches to REMEDY selection AFTER #1).

What it separates
-----------------
S1 (ocean-only, salinity STILL varying) left the strong Tb<->w salinity
degeneracy (corr ~ -0.986) intact. #1 removes it: fix salinity at the reference
posterior's median and train an ocean-only pilot on just that salinity slice.
  - Im_k2 gap COLLAPSES  -> the salinity axis / its degeneracy is the driver;
                            remedy discussion (manager + reviewer + user) starts
                            from degeneracy-aware options.
  - Im_k2 gap PERSISTS   -> salinity is NOT the driver; remaining candidates are
                            flow capacity / embedding (#4); elimination restarts
                            from the widened residual.
Either way: STOP after #1 and surface. Remedy selection is not Machine B's call.

SAMPLE-SIZE CONFOUND CONTROL (scientific-reviewer, PASS WITH CONCERNS 2026-08-05)
--------------------------------------------------------------------------------
The fixed-salinity band keeps only ~9% of the dataset (~55-62k rows) vs the
capped anchor's ~690k and S1's ~643k — an ~11x reduction. Under-resourced
training generically biases a flow toward UNDER-concentration (staying near the
prior), which is the SAME |Im_k2|~0.04 signature we are diagnosing. So reading
#1 against the 690k anchor alone confounds "salinity axis removed" with "11x less
data": a PERSIST could be a genuine capacity/embedding limit OR mere data
starvation, and the preregistered reading maps PERSIST -> capacity/embedding (#4).
To isolate the salinity axis at matched N we ALSO train a CONTROL: ocean-only,
salinity STILL VARYING, subsampled to the SAME n_kept, same seed/arch/cap.
Then:
  #1(banded, N)  vs  control(varying, N)   -> isolates the salinity axis at matched N
  control(varying, N)  vs  S1(varying, 643k) -> isolates the pure size effect
Do NOT read #1 against the 690k anchor alone; read the banded-vs-control delta.

Design (manager-pinned)
-----------------------
- Salinity FIXED at the reference posterior weighted median
  log10_wOcean_ppt = 1.1007 (= 12.61 ppt). Implemented as a symmetric +-0.084-dex
  band in log10_w about the median (half the 0.1677 log10 w-grid spacing), i.e. a
  geometric-mean-centred ~+-21% salinity slice [10.4, 15.3] ppt about 12.6 ppt.
  NOTE: this band is NOT the nearest-node Voronoi cell (the nearest w-node,
  14.93 ppt, has cell [1.090, 1.258]); the band is centred on the true median to
  avoid a +18% salinity bias, and the forward model blends BILINEARLY across
  nodes (grid_interp_2d.bilinear_weights), not by nearest-node snapping. The band
  is a symmetric physical salinity slice about the median, nothing more.
- Ocean-only rows (reuse the validated S1 nearest-node has_ocean recovery).
- SAME architecture / seed family / 60-epoch cap as the S-pilots + capped anchor.
  Cap-vs-cap comparison; but see the sample-size control above — the primary
  read is banded-vs-control, not banded-vs-690k-anchor.
- Verdict statistic: SBI-pp |Im_k2| median for #1 AND the control vs the
  matched-resolution MCMC-pp ceiling 0.1037 (#2) and the capped anchor (0.043),
  plus the four-way table and the flag statistic in sigma_obs (0.035).

NO NEW FORWARD SIMS. Reuses /tmp/titanG_build/datasets/titanG_nh3_1m.npz.
libomp: cache read + torch train in ONE process (numba forward model NOT imported
here); each PPC step is a SEPARATE subprocess (it imports the forward model).

Run:
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/nh3_diag_1_salinity_fixed.py
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
OUTDIR = ROOT / "validation_reports/nh3_diagnosis/f1_salinity_fixed"
TMP_ART = Path("/tmp/nh3_diag/f1")
TRAIN_SEED = 72

# Manager-pinned comparison anchors (provenance; NOT recomputed here).
REF_LOG10_W_MEDIAN = 1.1007         # reference posterior weighted median (12.61 ppt)
HALF_CELL_LOG10_W = 0.0839          # half of the 0.1677 log10 w-grid spacing
CAPPED_ANCHOR_IMK2 = 0.043          # capped full-joint anchor SBI-pp |Im_k2| median
MATCHED_MCMC_PP_CEILING = 0.1037    # #2 matched-resolution (n_eff=2000) ceiling
OBS_IMK2 = 0.135
SIGMA_IMK2 = 0.035


def _has_ocean_grid(cache):
    S = cache["structures"]
    return np.array([bool(isinstance(n, dict) and n.get("has_ocean")) for n in S],
                    dtype=bool)


def _row_is_ocean(theta, param_names, cache):
    """Per-row ocean-branch membership from the nearest grid node (S1 logic)."""
    Tb_grid = np.asarray(cache["Tb_K_grid"], float)
    w_grid = np.asarray(cache["wOcean_ppt_grid"], float)
    ho = _has_ocean_grid(cache)
    i_tb = param_names.index("Tb_K")
    i_lw = param_names.index("log10_wOcean_ppt")
    out = np.zeros(len(theta), dtype=bool)
    for r in range(len(theta)):
        corners, weights = g2.bilinear_weights(
            Tb_grid, w_grid, float(theta[r, i_tb]), float(10.0 ** theta[r, i_lw]))
        out[r] = bool(ho[corners[int(np.argmax(weights))]])
    return out


def _train_and_ppc(tag, theta_sub, x_sub, stats, runner_cfg, args, extra_info):
    """Train one pilot flow on (theta_sub, x_sub), save it, and run PPC in a
    separate subprocess. Returns a dict of run metadata (incl. epochs_trained).
    A FRESH runner is built per call so the two pilots never share fit state."""
    runner = SBIRunner(InferenceConfig.from_dict(runner_cfg))
    i_lw = list(runner.param_names).index("log10_wOcean_ppt")
    kept_lw = theta_sub[:, i_lw]
    print(f"[f1:{tag}] training: theta={theta_sub.shape} x={x_sub.shape} "
          f"seed={TRAIN_SEED} {args.density_estimator} max_epochs={args.max_epochs}")
    print(f"[f1:{tag}] log10_w 5/50/95 = {np.percentile(kept_lw, [5, 50, 95])} "
          f"(ppt {10.0 ** np.percentile(kept_lw, [5, 50, 95])})")
    t1 = time.time()
    runner.train(theta_sub, x_sub, seed=TRAIN_SEED,
                 density_estimator=args.density_estimator,
                 max_num_epochs=args.max_epochs)
    runner._train_info["rejection_stats"] = stats
    runner._train_info[f"f1_filter_{tag}"] = extra_info
    art = TMP_ART / f"nh3_f1_{tag}_pilot.pt"
    runner.save_artifact(str(art))
    _ts = getattr(runner, "_last_train_summary", {}) or {}
    _ep = _ts.get("epochs_trained") or _ts.get("epochs")
    if isinstance(_ep, list):
        _ep = _ep[-1] if _ep else None
    early = (_ep is not None) and (_ep < args.max_epochs - 5)
    print(f"[f1:{tag}] trained in {(time.time()-t1)/60:.1f} min -> {art}; "
          f"epochs_trained={_ep} (early-stop parity flag: {early})")

    filt_npz = TMP_ART / f"nh3_f1_{tag}_dataset.npz"
    np.savez_compressed(filt_npz, theta=theta_sub, x=x_sub,
                        stats=json.dumps(stats, default=str))

    ppc_out = OUTDIR / tag
    ppc_out.mkdir(parents=True, exist_ok=True)
    print(f"[f1:{tag}] launching PPC subprocess at Titan datum ...")
    cmd = [
        "python", str(ROOT / "plans/scripts/titanG_ppc_interior_check.py"),
        "--artifact", str(art), "--config", str(CFG),
        "--dataset", str(filt_npz), "--output-dir", str(ppc_out),
        "--n-post", str(args.n_post), "--seed", str(TRAIN_SEED),
    ]
    env = dict(os.environ, PYTHONPATH=str(ROOT),
               NUMBA_CACHE_DIR="/tmp/pp_numba_cache", KMP_DUPLICATE_LIB_OK="TRUE")
    r = subprocess.run(cmd, env=env, cwd=str(ROOT))
    print(f"[f1:{tag}] PPC exit code {r.returncode}; report under {ppc_out}")
    return {
        "tag": tag, "artifact": str(art), "ppc_output_dir": str(ppc_out),
        "n_train": int(len(theta_sub)), "epochs_trained": _ep,
        "early_stop_parity_flag": bool(early),
        "config_hash": runner.config.generate_hash(),
        "extra": extra_info,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-epochs", type=int, default=60,
                    help="cap NPE epochs (pilot; matches S-pilots + capped anchor)")
    ap.add_argument("--density-estimator", default="nsf")
    ap.add_argument("--n-post", type=int, default=4000)
    ap.add_argument("--half-cell", type=float, default=HALF_CELL_LOG10_W,
                    help="half-width in log10_w of the fixed-salinity band")
    ap.add_argument("--skip-control", action="store_true",
                    help="(debug) skip the matched-N varying control")
    args = ap.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    TMP_ART.mkdir(parents=True, exist_ok=True)
    if not DATASET.exists():
        raise SystemExit(f"[f1] missing dataset {DATASET}")

    with np.load(DATASET, allow_pickle=True) as d:
        theta = np.asarray(d["theta"], np.float64)
        x = np.asarray(d["x"], np.float64)
        stats = json.loads(str(d["stats"].item())) if "stats" in d.files else {}

    cfg = json.load(open(CFG)); cfg["mode"] = "sbi"
    probe = SBIRunner(InferenceConfig.from_dict(cfg))
    param_names = list(probe.param_names)
    obs_names = list(probe.obs_names)
    assert theta.shape[1] == len(param_names), "theta col mismatch"

    cache = pickle.load(open(CACHE, "rb"))
    i_lw = param_names.index("log10_wOcean_ppt")
    lw = theta[:, i_lw]

    print(f"[f1] recovering has_ocean for {len(theta):,} rows via nearest grid node ...")
    t0 = time.time()
    is_ocean = _row_is_ocean(theta, param_names, cache)
    n_ocean = int(is_ocean.sum())
    # Salinity fixed: keep ocean rows within a symmetric half-grid-cell band in
    # log10_w of the reference posterior median (a ~+-21% salinity slice about
    # 12.6 ppt). This is a physical salinity slice, NOT a nearest-node cell.
    in_band = np.abs(lw - REF_LOG10_W_MEDIAN) < args.half_cell
    keep = is_ocean & in_band
    n_keep = int(keep.sum())
    print(f"[f1] ocean={n_ocean:,}  in-band={int(in_band.sum()):,}  "
          f"ocean&fixed-salinity={n_keep:,} / {len(theta):,} "
          f"({100.0*n_keep/len(theta):.2f}%); recover {time.time()-t0:.1f}s")
    if n_keep < 5000:
        print(f"[f1] WARNING: only {n_keep} rows — pilot statistics may be thin")

    runs = {}

    # --- Pilot #1: salinity-fixed (banded) ocean-only ---
    theta_b = theta[keep]
    x_b = x[keep]
    runs["banded"] = _train_and_ppc(
        "banded", theta_b, x_b, stats, cfg, args,
        extra_info={
            "kind": "ocean_branch_only_AND_salinity_fixed",
            "ref_log10_w_median": REF_LOG10_W_MEDIAN,
            "ref_w_ppt_median": float(10.0 ** REF_LOG10_W_MEDIAN),
            "half_cell_log10_w": args.half_cell,
            "n_kept": n_keep, "n_total_dataset": int(len(theta)),
            "has_ocean_recovery": "nearest-grid-node via grid_interp_2d.bilinear_weights",
            "band_note": ("symmetric +-half-cell log10_w slice about the median; "
                          "NOT a nearest-node Voronoi cell; forward model blends "
                          "bilinearly across nodes."),
            "cache": str(CACHE),
        })

    # --- Control: matched-N, salinity-VARYING ocean-only (reviewer MAJOR) ---
    if not args.skip_control:
        ocean_idx = np.flatnonzero(is_ocean)
        rng = np.random.default_rng(TRAIN_SEED)
        pick = rng.choice(ocean_idx, size=n_keep, replace=False)
        pick.sort()
        theta_c = theta[pick]
        x_c = x[pick]
        runs["control_varying"] = _train_and_ppc(
            "control_varying", theta_c, x_c, stats, cfg, args,
            extra_info={
                "kind": "ocean_branch_only_salinity_VARYING_matched_N",
                "n_kept": int(n_keep), "n_ocean_pool": n_ocean,
                "subsample_seed": TRAIN_SEED,
                "purpose": ("matched-N control isolating the salinity axis from "
                            "the ~11x training-set-size reduction; read "
                            "banded-vs-control at matched N, NOT banded-vs-690k."),
                "cache": str(CACHE),
            })
    else:
        print("[f1] --skip-control set: matched-N varying control NOT trained")

    manifest = {
        "kind": "nh3_diag_f1_salinity_fixed_pilot",
        "config": str(CFG), "cache": str(CACHE),
        "ref_log10_w_median": REF_LOG10_W_MEDIAN,
        "ref_w_ppt_median": float(10.0 ** REF_LOG10_W_MEDIAN),
        "half_cell_log10_w": args.half_cell,
        "n_ocean": n_ocean, "n_kept_banded": n_keep, "n_total": int(len(theta)),
        "train_seed": TRAIN_SEED, "density_estimator": args.density_estimator,
        "max_epochs": args.max_epochs,
        "obs_names": obs_names,
        "runs": runs,
        "comparison_anchors": {
            "capped_anchor_imk2": CAPPED_ANCHOR_IMK2,
            "matched_mcmc_pp_ceiling": MATCHED_MCMC_PP_CEILING,
            "obs_imk2": OBS_IMK2, "sigma_imk2": SIGMA_IMK2,
        },
        "anchor_note": ("PRIMARY read = banded-vs-control at matched N (isolates "
                        "the salinity axis from the ~11x size reduction). SECONDARY "
                        "= each pilot's SBI-pp |Im_k2| vs the capped anchor (0.043) "
                        "and matched MCMC-pp ceiling (0.1037); cap-vs-cap."),
        "preregistered_reading": (
            "banded-vs-control Im_k2 gap COLLAPSES -> salinity axis/degeneracy is "
            "the driver; gap PERSISTS (banded ~= control, both under-updated) -> "
            "capacity/embedding (#4). If the CONTROL alone under-updates vs S1's "
            "643k, that is the pure size effect and must be subtracted before "
            "reading the salinity axis. STOP after #1 and surface; remedy "
            "selection is manager + reviewer + user."),
        "reviewer": ("PASS WITH CONCERNS (2026-08-05): matched-N varying control "
                     "added (MAJOR); band rationale corrected to a symmetric "
                     "salinity slice, not a Voronoi cell (MODERATE); epochs_trained "
                     "+ early-stop parity flag recorded per pilot (MINOR)."),
    }
    with open(OUTDIR / "f1_train_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, default=str)
    print(f"[f1] manifest -> {OUTDIR/'f1_train_manifest.json'}")
    print(f"[f1] pilots: {list(runs.keys())}; PPC reports under {OUTDIR}/<tag>/")


if __name__ == "__main__":
    main()
