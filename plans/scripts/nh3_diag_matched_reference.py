"""NH3 under-update diagnosis — FOLLOW-UP #2: matched-resolution reference MCMC.

Manager judgement (MACHINE-B-HANDOFF §0.9, 2026-08-05) reordered the reviewer's
follow-ups so #2 runs after #3. #2 is the DECISIVE, cheap step: the entire
"flow under-updates by 0.042-vs-0.100" magnitude is measured against an NH3
reference MCMC run at n_effective=500 — exactly the resolution class B3 just
discredited (the 1.06 km v5/v7 reference wander was an n_eff=500 artifact that
shrank ~5.6x at n_eff=2000). Before chasing the gap with an expensive retrain
(#1), re-measure the TARGET: re-run the NH3 JOINT reference MCMC at n_eff=2000
across >=2 seeds and re-measure the weighted |Im_k2| MCMC-pp "ceiling".

Methodology (identical to B3, reviewer-ratified)
-----------------------------------------------
pocoMC 1.2.6 is a preconditioned SMC sampler with NO n_live; its ensemble knob
is n_effective (default 512) with n_active (default 256) the active count. The
existing NH3 reference ran n_effective=500 (n_active pocoMC-default 256). Raising
n_effective ALONE changes the sampler REGIME (the n_effective/n_active ratio
drives flow-train cadence + annealing). To quadruple the ensemble while keeping
the originals' regime (ratio ~2), raise n_active in step: n_active=1024 with
n_effective=2000. Both overridable. The TRACKED config is NEVER mutated (its
config_hash is referenced by the deployed flow + existing reference); the knobs
are set on the runner instance, exactly as b3_reference_wander does.

The ceiling
-----------
The target statistic is the WEIGHTED |Im_k2| median from each run's k2_results
(col1 is SIGNED Im; |Im_k2| = abs(col1); weights are the importance weights).
This is the "0.100 ceiling" the flow is judged against. We report per-seed
weighted 5/50/95, the between-seed scatter of the median, and — for the release
decision — the gap between the deployed-flow SBI-pp median (0.0423, from #3) and
the matched-resolution MCMC-pp median, expressed in units of sigma_obs (0.035).

DIAGNOSTIC STUDY. It does NOT replace the committed n_eff=500 reference pkl
(that stays as the ratification-time artifact); it re-measures the ceiling so the
manager's preregistered MgSO4/NaCl release criterion can be evaluated. Per-seed
pickles + manifests to /tmp (Dropbox sync-race avoidance); the aggregate report
is copied into validation_reports/.

Run (separate process; env pinned per repo convention):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/nh3_diag_matched_reference.py \
    --seeds 72 172 --n-effective 2000 --n-active 1024
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.mcmc_runner import MCMCRunner

ROOT = Path(__file__).resolve().parents[2]
CFG = ROOT / "PlanetProfile/Inference/configs/test54_titan_nh3_freegrav.json"
OUT_DIR = ROOT / "validation_reports/nh3_diagnosis/matched_reference"
TMP_DIR = Path("/tmp/nh3_diag/matched_reference")

# From follow-up #3 (deployed 1M flow pushforward): the SBI-pp |Im_k2| median the
# ceiling is compared against. Recorded here as provenance; NOT recomputed.
DEPLOYED_SBI_PP_IMK2_MEDIAN = 0.0423
LEGACY_NEFF500_CEILING = 0.0999   # weighted |Im_k2| median of the n_eff=500 ref
OBS_IMK2 = 0.135
SIGMA_IMK2 = 0.035
FLAG_THRESHOLD_SIGMA = 0.5        # manager release criterion (< 0.5 sigma_obs -> PROCEED)


def _sha256(path):
    if not path or not os.path.exists(str(path)):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_sha():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip()
    except Exception:
        return None


def _env_versions():
    vers = {}
    for mod in ("scipy", "pocomc", "numpy", "torch"):
        try:
            vers[mod] = __import__(mod).__version__
        except Exception:
            vers[mod] = None
    return vers


def _weighted_quantiles(x, w, qs):
    i = np.argsort(x)
    x = x[i]; w = w[i]
    c = np.cumsum(w)
    if c[-1] <= 0:
        return np.full(len(qs), np.nan)
    c = c / c[-1]
    return np.interp(qs, c, x)


def _run_one(cfg, seed, n_effective, n_active, tmp_dir):
    runner = MCMCRunner(cfg)
    runner.random_state = int(seed)
    runner.n_effective = int(n_effective)     # override config n_effective
    if n_active is not None:
        runner.n_active = int(n_active)        # raise in step (regime-preserving)
    jsonl = tmp_dir / f"prog_seed{seed}.jsonl"
    t0 = time.time()
    print(f"[ref2] seed={seed} n_effective={n_effective} n_active={n_active} "
          f"starting ...", flush=True)
    res = runner.run(progress_jsonl_path=str(jsonl))
    dt = time.time() - t0
    md = getattr(res, "metadata", {}) or {}
    kk = np.asarray(res.k2_results, dtype=float)   # (N,2) [Re, Im signed]
    ww = np.asarray(res.weights, dtype=float) if res.weights is not None \
        else np.ones(kk.shape[0])
    fin = np.all(np.isfinite(kk), axis=1) & np.isfinite(ww)
    im = np.abs(kk[fin, 1]); w = ww[fin]
    q = _weighted_quantiles(im, w, [0.05, 0.5, 0.95])
    print(f"[ref2] seed={seed} done {dt/60:.1f} min: n_k2={im.size} "
          f"log_Z={md.get('log_Z')} |Im_k2| wtd 5/50/95={q}", flush=True)
    res.save(str(tmp_dir / f"nh3_reference_seed{seed}_neff{n_effective}.pkl"))
    return {
        "seed": int(seed), "n_effective": int(n_effective),
        "n_active": (int(n_active) if n_active is not None else None),
        "elapsed_s": dt, "n_samples": int(np.asarray(res.samples).shape[0]),
        "n_k2_finite": int(im.size),
        "log_Z": md.get("log_Z"), "log_Z_err": md.get("log_Z_err"),
        "imk2_weighted_p5": float(q[0]),
        "imk2_weighted_median": float(q[1]),
        "imk2_weighted_p95": float(q[2]),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(CFG))
    ap.add_argument("--seeds", type=int, nargs="+", default=[72, 172])
    ap.add_argument("--n-effective", type=int, default=2000)
    ap.add_argument("--n-active", type=int, default=1024,
                    help="pocoMC active-particle count raised in step with "
                         "n_effective to preserve the sampler regime (ratio ~2)")
    ap.add_argument("--out", default=str(OUT_DIR))
    args = ap.parse_args()

    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    cfg = InferenceConfig.from_json(args.config)
    cache_path = cfg.structure_cache_path
    prov = {
        "config": os.path.relpath(args.config, ROOT),
        "config_hash": cfg.generate_hash(),
        "structure_cache_path": str(cache_path),
        "structure_cache_sha256": _sha256(cache_path),
        "n_effective": int(args.n_effective),
        "n_active": int(args.n_active),
        "config_n_effective_original": cfg.sampler_settings.get("n_effective"),
        "obs_names": list(cfg.observables),
        "param_names": list(cfg.param_space.keys()),
        "env_versions": _env_versions(),
        "git_sha": _git_sha(),
    }
    print(f"[ref2] config_hash={prov['config_hash']} "
          f"(original n_eff={prov['config_n_effective_original']}); "
          f"cache sha256={prov['structure_cache_sha256']}")
    print(f"[ref2] observables={prov['obs_names']}")

    runs = []
    for sd in args.seeds:
        runs.append(_run_one(cfg, sd, args.n_effective, args.n_active, TMP_DIR))

    meds = np.array([r["imk2_weighted_median"] for r in runs])
    pooled_median = float(np.mean(meds))
    between_seed_std = float(np.std(meds, ddof=1)) if meds.size > 1 else None
    between_seed_range = float(meds.max() - meds.min()) if meds.size else None

    # Release-criterion arithmetic (manager §0.9): the reported target is
    # SBI-pp vs matched-resolution MCMC-pp, expressed in sigma_obs.
    gap = pooled_median - DEPLOYED_SBI_PP_IMK2_MEDIAN
    gap_in_sigma = gap / SIGMA_IMK2
    ceiling_move_vs_legacy = pooled_median - LEGACY_NEFF500_CEILING

    report = {
        "kind": "nh3_diag_matched_resolution_reference_followup2",
        "note": ("DIAGNOSTIC, not a gate. Re-measures the weighted |Im_k2| "
                 "MCMC-pp ceiling at n_eff=2000 (B3-matched resolution) so the "
                 "flow under-update is judged SBI-pp vs matched MCMC-pp, not vs "
                 "the n_eff=500 ceiling. Does NOT replace the committed reference."),
        "provenance": prov,
        "seeds": list(args.seeds),
        "runs": runs,
        "matched_resolution_ceiling": {
            "pooled_weighted_imk2_median": pooled_median,
            "between_seed_std": between_seed_std,
            "between_seed_range": between_seed_range,
            "per_seed_median": {int(r["seed"]): r["imk2_weighted_median"]
                                for r in runs},
        },
        "legacy_neff500_ceiling": LEGACY_NEFF500_CEILING,
        "ceiling_move_vs_legacy": ceiling_move_vs_legacy,
        "deployed_sbi_pp_imk2_median": DEPLOYED_SBI_PP_IMK2_MEDIAN,
        "obs_imk2": OBS_IMK2, "sigma_imk2": SIGMA_IMK2,
        "sbi_vs_matched_mcmc_gap": gap,
        "sbi_vs_matched_mcmc_gap_in_sigma_obs": gap_in_sigma,
        "flag_threshold_sigma": FLAG_THRESHOLD_SIGMA,
        "release_reading": (
            "manager §0.9: if |SBI-pp - matched MCMC-pp| < 0.5*sigma_obs the "
            "under-update is substantially a reference artifact -> MgSO4/NaCl "
            "PROCEED (standard gates + pushforward gate) and NH3 split "
            "ratification returns to the manager; if the gap survives -> run #1 "
            "(salinity-fixed retrain) with reviewer+user sign-off before any "
            "MgSO4/NaCl compute. This script reports the number; the "
            "scientific-reviewer + manager adjudicate."),
    }
    with open(out_dir / "matched_reference_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\n[ref2] report -> {out_dir/'matched_reference_report.json'}")
    print(f"[ref2] matched-resolution ceiling (weighted |Im_k2| median) = "
          f"{pooled_median:.4f}  (legacy n_eff=500 = {LEGACY_NEFF500_CEILING}; "
          f"move {ceiling_move_vs_legacy:+.4f})")
    print(f"[ref2] between-seed: std={between_seed_std} range={between_seed_range}")
    print(f"[ref2] SBI-pp {DEPLOYED_SBI_PP_IMK2_MEDIAN} vs matched MCMC-pp "
          f"{pooled_median:.4f}  gap={gap:+.4f} = {gap_in_sigma:+.2f} sigma_obs "
          f"(flag threshold {FLAG_THRESHOLD_SIGMA})")

    # Copy per-seed pickles into the report dir for reader confirmation.
    for sd in args.seeds:
        src = TMP_DIR / f"nh3_reference_seed{sd}_neff{args.n_effective}.pkl"
        if src.exists():
            shutil.copy2(src, out_dir / src.name)
        prg = TMP_DIR / f"prog_seed{sd}.jsonl"
        if prg.exists():
            shutil.copy2(prg, out_dir / prg.name)


if __name__ == "__main__":
    main()
