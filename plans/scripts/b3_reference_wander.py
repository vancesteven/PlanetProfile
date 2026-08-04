"""B3 reference-wander study — multi-seed v5 + v7 reference MCMC (pocoMC).

Adjudication context (plans/active/europa-v5v6v7-gate-adjudication.md, item #3):
the v5 and v7 reference posteriors differ by 1.06 km in D_iceIh_km mean (~14x
the quoted MC error), a near-rigid translation (sigma ratio 0.975) with the
WRONG sign for a support effect — the signature of nested-sampling correlated
log-volume error, NOT of the added |Ae| support. Which reference wandered is
OPEN. This driver measures the EMPIRICAL between-seed reference-noise floor by
running each config at >=3 fresh seeds, so we can judge whether the observed
1.06 km wander is bracketed by seed-to-seed scatter (§0.7 step 1; feeds the v5
shape-excess re-evaluation in step 3).

Resolution mapping (methodology decision — reviewer-ratified before launch)
--------------------------------------------------------------------------
pocoMC 1.2.6 is a PRECONDITIONED SMC sampler. It has NO ``n_live`` argument:
its ensemble-size knob is ``n_effective`` (default 512) with ``n_active``
(default 256) the active-particle count. The originals ran ``n_effective=500``
(n_active left at pocoMC default 256), which the adjudication doc itself labels
"n_live=500" (line 58). The manager's B3 spec "n_live >= 2000" therefore maps to
``n_effective=2000`` — 4x the originals' ensemble; n_effective is the knob that
reduces the correlated annealing/evidence error blamed for the wander.

Reviewer (scientific-reviewer, PASS-WITH-CONCERNS): raising n_effective ALONE
while n_active stays pinned at 256 changes the sampler REGIME (the
n_effective/n_active ratio drives flow train cadence + annealing), not a clean
4x resolution scale. To keep the originals' regime (ratio ~2, train_frequency=1)
while quadrupling the ensemble, raise n_active in step: ``n_active=1024`` with
``n_effective=2000``. Both are overridable (--n-effective / --n-active).

Determinism / provenance
------------------------
- Each run pins ``runner.random_state`` to the CLI seed (independent of the
  config's random_state, which is deliberately ignored here so seeds sweep).
- The forward model reads the pinned structure cache + Ae sidecar named in each
  config; their sha256 are recorded per config for reader confirmation.
- Per-seed pickles + manifests go to /tmp (repo memory: avoid Dropbox sync
  races); the aggregate scatter report is copied into validation_reports/.

DIAGNOSTIC STUDY — reports whatever the scatter gives. The committed primary
references (seed-51 v5, seed-71 v7) are NOT replaced by this script; success is
a measurement (does seed scatter bracket 1.06 km?), never a reference selected
to pass a gate.

Run (separate process; env pinned per repo convention):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/b3_reference_wander.py \
    --which both --seeds 101 202 303 --n-effective 2000
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

from PlanetProfile.Inference.inference_core import InferenceConfig, InferenceResult
from PlanetProfile.Inference.mcmc_runner import MCMCRunner

ROOT = Path(__file__).resolve().parents[2]

CONFIGS = {
    "v5": {
        "config": ROOT / "PlanetProfile/Inference/configs/europa_clipper_v5_geodesy_11D.json",
        "primary": (ROOT / "PlanetProfile/Test/mcmc_results/Europa/Test52_seawater_v5/"
                    "europa_clipper_v5_reference_result.pkl"),
        "primary_seed": 51,
    },
    "v7": {
        "config": ROOT / "PlanetProfile/Inference/configs/europa_clipper_v7_openae_11D.json",
        "primary": (ROOT / "PlanetProfile/Test/mcmc_results/Europa/Test54_seawater_v7/"
                    "europa_clipper_v7_reference_result.pkl"),
        "primary_seed": 71,
    },
}
# Parameters the adjudication tracks as the wander/degeneracy signal.
KEY_PARAMS = ["D_iceIh_km", "log10_wOcean_ppt", "Tb_K"]
OBSERVED_WANDER_KM = 1.06  # v5-vs-v7 D_iceIh mean gap the scatter must bracket


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
    """Record the versions the B3 spec's 'pinned environment' clause needs
    (scipy straddled 1.16.3->1.17.1 across the originals)."""
    vers = {}
    for mod in ("scipy", "pocomc", "numpy", "torch"):
        try:
            vers[mod] = __import__(mod).__version__
        except Exception:
            vers[mod] = None
    return vers


def _run_one(cfg, seed, n_effective, n_active, tmp_dir):
    runner = MCMCRunner(cfg)
    runner.random_state = int(seed)
    runner.n_effective = int(n_effective)  # override config's n_effective
    if n_active is not None:
        runner.n_active = int(n_active)     # raise in step (regime-preserving)
    jsonl = tmp_dir / f"prog_seed{seed}.jsonl"
    t0 = time.time()
    print(f"[b3]   seed={seed} n_effective={n_effective} n_active={n_active} "
          f"starting ...", flush=True)
    res = runner.run(progress_jsonl_path=str(jsonl))
    dt = time.time() - t0
    s = np.asarray(res.samples, dtype=np.float64)
    md = getattr(res, "metadata", {}) or {}
    print(f"[b3]   seed={seed} done {dt/60:.1f} min samples={s.shape} "
          f"log_Z={md.get('log_Z')}", flush=True)
    res.save(str(tmp_dir / f"reference_seed{seed}.pkl"))
    return s, md, dt, list(res.param_names)


def _study(which, seeds, n_effective, n_active, include_primary):
    spec = CONFIGS[which]
    cfg = InferenceConfig.from_json(str(spec["config"]))
    param_names = list(cfg.param_space.keys())
    tmp_dir = Path(f"/tmp/b3_build/{which}_reference_wander")
    tmp_dir.mkdir(parents=True, exist_ok=True)

    cache_path = cfg.structure_cache_path
    sidecar_path = str(cache_path) + ".ae_sidecar.pkl"
    prov = {
        "config": os.path.relpath(str(spec["config"]), ROOT),
        "config_hash": cfg.generate_hash(),
        "structure_cache_path": str(cache_path),
        "structure_cache_sha256": _sha256(cache_path),
        "ae_sidecar_sha256": _sha256(sidecar_path),
        "n_effective": int(n_effective),
        # n_active plumbed through MCMCRunner; None -> pocoMC default 256.
        "n_active": (int(n_active) if n_active is not None else "pocomc_default_256"),
        "n_total_termination_ess": 4096,   # hardcoded in mcmc_runner.run()
        "n_evidence_default": 4096,         # pocoMC evidence-IS sample count
        "env_versions": _env_versions(),
        "param_names": param_names,
    }
    print(f"[b3] === {which} === config_hash={prov['config_hash']}")
    print(f"[b3] cache sha256={prov['structure_cache_sha256']}")

    runs = []
    # Optionally seed the study with the committed primary reference (as-is; it
    # was run at n_effective=500, so it is labelled distinctly, NOT pooled into
    # the equal-resolution scatter statistic).
    if include_primary and os.path.exists(str(spec["primary"])):
        prim = InferenceResult.load(str(spec["primary"]))
        runs.append({
            "label": f"seed{spec['primary_seed']}_primary_neff500",
            "seed": spec["primary_seed"], "n_effective_run": 500,
            "samples": np.asarray(prim.samples, dtype=np.float64),
            "metadata": getattr(prim, "metadata", {}) or {},
            "equal_res": False,
        })

    for sd in seeds:
        s, md, dt, pnames = _run_one(cfg, sd, n_effective, n_active, tmp_dir)
        runs.append({
            "label": f"seed{sd}_neff{n_effective}", "seed": int(sd),
            "n_effective_run": int(n_effective), "samples": s,
            "metadata": md, "elapsed_s": dt, "equal_res": True,
            "param_names": pnames,
        })

    # Per-parameter mean + std per run; scatter computed over EQUAL-RESOLUTION
    # runs only (the fresh n_effective seeds), primary reported alongside.
    per_param = {}
    equal = [r for r in runs if r["equal_res"]]
    for name in param_names:
        j = param_names.index(name)
        means = {r["label"]: float(np.mean(r["samples"][:, j])) for r in runs}
        stds = {r["label"]: float(np.std(r["samples"][:, j])) for r in runs}
        eq_means = np.array([np.mean(r["samples"][:, j]) for r in equal])
        # Per-seed means over equal-resolution runs, keyed by seed, so the
        # cross-config aggregation in main() can pair v5 vs v7 at matched seeds
        # and matched resolution (the number step 3 needs).
        eq_means_by_seed = {int(r["seed"]): float(np.mean(r["samples"][:, j]))
                            for r in equal}
        entry = {
            "means": means, "stds": stds,
            "equal_res_means_by_seed": eq_means_by_seed,
            "between_seed_mean_range": (float(eq_means.max() - eq_means.min())
                                        if eq_means.size else None),
            "between_seed_mean_std": (float(eq_means.std(ddof=1))
                                      if eq_means.size > 1 else None),
        }
        per_param[name] = entry

    report = {
        "kind": "b3_reference_wander",
        "which": which,
        "provenance": prov,
        "git_sha": _git_sha(),
        "observed_v5_v7_D_iceIh_wander_km": OBSERVED_WANDER_KM,
        "seeds": list(seeds),
        "runs": [{"label": r["label"], "seed": r["seed"],
                  "n_effective_run": r["n_effective_run"],
                  "n_samples": int(r["samples"].shape[0]),
                  "log_Z": (r["metadata"] or {}).get("log_Z"),
                  "log_Z_err": (r["metadata"] or {}).get("log_Z_err"),
                  "elapsed_s": r.get("elapsed_s"),
                  "equal_res": r["equal_res"]} for r in runs],
        "per_param": per_param,
        "key_params": KEY_PARAMS,
    }
    return report


def _paired_cross_config_gap(reports, out_dir):
    """Matched-resolution v5-vs-v7 gap (reviewer-required for step 3).

    The 1.06 km anchor was a v5(n_eff=500)-vs-v7(n_eff=500) gap. Comparing it
    to a within-config floor measured at n_eff=2000 is apples-to-oranges. This
    computes the v5-v7 D_iceIh (and log10_w) mean gap at MATCHED seeds and
    MATCHED resolution, and reports it against each config's within-config
    between-seed floor at the same resolution — the number that decides whether
    the wander is a low-resolution artifact (gap shrinks toward floor) or a real
    config/support effect (gap persists while floor collapses)."""
    if "v5" not in reports or "v7" not in reports:
        return None
    v5, v7 = reports["v5"], reports["v7"]
    gap = {"kind": "b3_matched_resolution_v5_v7_gap",
           "legacy_neff500_gap_km": OBSERVED_WANDER_KM,
           "note": ("legacy gap was v5(n_eff=500) vs v7(n_eff=500); the matched "
                    "gaps below are at the fresh resolution and ARE the step-3 "
                    "comparison — the legacy number is NOT directly comparable."),
           "per_param": {}}
    for name in ("D_iceIh_km", "log10_wOcean_ppt"):
        # In-memory reports: keys are ints (not JSON-stringified).
        e5 = v5["per_param"].get(name, {}).get("equal_res_means_by_seed", {})
        e7 = v7["per_param"].get(name, {}).get("equal_res_means_by_seed", {})
        shared = sorted(set(e5) & set(e7))
        per_seed = {int(s): float(e5[s]) - float(e7[s]) for s in shared}
        gaps = np.array(list(per_seed.values())) if per_seed else np.array([])
        gap["per_param"][name] = {
            "per_seed_v5_minus_v7": per_seed,
            "mean_gap": float(gaps.mean()) if gaps.size else None,
            "gap_std": float(gaps.std(ddof=1)) if gaps.size > 1 else None,
            "v5_within_config_floor_std": v5["per_param"].get(name, {}).get("between_seed_mean_std"),
            "v7_within_config_floor_std": v7["per_param"].get(name, {}).get("between_seed_mean_std"),
        }
    out_path = out_dir / "matched_resolution_v5_v7_gap.json"
    with open(out_path, "w") as f:
        json.dump(gap, f, indent=2, default=str)
    print(f"\n[b3] matched-resolution v5-v7 gap -> {out_path}")
    for name, e in gap["per_param"].items():
        print(f"[b3]   {name:20s} matched mean gap={e['mean_gap']} "
              f"(v5 floor std={e['v5_within_config_floor_std']}, "
              f"v7 floor std={e['v7_within_config_floor_std']}); "
              f"legacy n_eff=500 gap {OBSERVED_WANDER_KM} km NOT comparable")
    return gap


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", choices=["v5", "v7", "both"], default="both")
    ap.add_argument("--seeds", type=int, nargs="+", default=[101, 202, 303])
    ap.add_argument("--n-effective", type=int, default=2000)
    ap.add_argument("--n-active", type=int, default=1024,
                    help="pocoMC active-particle count raised in step with "
                         "n_effective to preserve the sampler regime (ratio ~2)")
    ap.add_argument("--no-primary", action="store_true",
                    help="exclude the committed n_eff=500 primary from the report")
    ap.add_argument("--out", default=str(ROOT / "validation_reports/b3_reference_wander"))
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    whichs = ["v5", "v7"] if args.which == "both" else [args.which]

    reports = {}
    for which in whichs:
        rep = _study(which, args.seeds, args.n_effective, args.n_active,
                     include_primary=not args.no_primary)
        reports[which] = rep
        out_path = out_dir / f"{which}_reference_wander.json"
        with open(out_path, "w") as f:
            json.dump(rep, f, indent=2, default=str)
        print(f"\n[b3] {which} report -> {out_path}")
        for name in KEY_PARAMS:
            e = rep["per_param"].get(name)
            if not e:
                continue
            print(f"[b3]   {name:20s} within-config between-seed mean range="
                  f"{e['between_seed_mean_range']}  std={e['between_seed_mean_std']}")

    # Matched-resolution cross-config gap (the step-3 comparison).
    if len(reports) == 2:
        _paired_cross_config_gap(reports, out_dir)


if __name__ == "__main__":
    main()
