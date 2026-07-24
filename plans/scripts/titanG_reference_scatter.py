"""Reference-side diagnostic for the Titan no-ocean crosscheck residual
(scientific-reviewer a2d7f121 required-validation item 2).

The re-gated flow is uniformly ~1.17x over-dispersed vs a SINGLE pocoMC
reference (n_eff=500). A single reference carries its own few-percent
sigma scale scatter. This runs N INDEPENDENT pocoMC references at different
seeds and reports the per-parameter sigma_MCMC spread, so we can judge whether
the 1.17x flow/reference ratio sits within reference-to-reference scatter plus
the known mild NSF over-dispersion.

DIAGNOSTIC ONLY — run once, report, DO NOT iterate or select a reference to
pass the gate (that would be gate-tuning-by-reference-selection). The primary
crosscheck reference remains the seed-71 pickle.

Run (separate process):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/titanG_reference_scatter.py \
    --seeds 171 271
"""
import argparse, json, time
from pathlib import Path
import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig, InferenceResult
from PlanetProfile.Inference.mcmc_runner import MCMCRunner

ROOT = Path(__file__).resolve().parents[2]
CFG = ROOT / "PlanetProfile/Inference/configs/titan_freegrav_noocean.json"
PRIMARY = (ROOT / "validation_reports/titan_freegrav_noocean_1m/reference/"
           "titan_freegrav_noocean_reference_result.pkl")
OUT = ROOT / "validation_reports/titan_freegrav_noocean_1m/reference_scatter"
TMP = Path("/tmp/titanG_build/reference_scatter")


def _std_by_param(samples):
    return np.asarray(samples, dtype=np.float64).std(axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=[171, 271])
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)

    cfg = InferenceConfig.from_json(str(CFG))
    param_names = list(cfg.param_space.keys())

    # Primary (seed-71) reference stds as the baseline.
    primary = InferenceResult.load(str(PRIMARY))
    runs = [("seed71_primary", 71, np.asarray(primary.samples, np.float64))]

    for sd in args.seeds:
        runner = MCMCRunner(cfg)
        runner.random_state = sd
        t0 = time.time()
        print(f"[scatter] running reference seed={sd} ...", flush=True)
        res = runner.run(progress_jsonl_path=str(TMP / f"prog_{sd}.jsonl"))
        dt = time.time() - t0
        s = np.asarray(res.samples, np.float64)
        print(f"[scatter]   seed={sd} done {dt/60:.1f} min, samples={s.shape}")
        res.save(str(TMP / f"reference_seed{sd}.pkl"))
        runs.append((f"seed{sd}", sd, s))

    # Per-parameter sigma across the references.
    stds = {name: [] for name in param_names}
    order = list(runs[0][2].shape[1:]) and param_names
    for label, sd, s in runs:
        sp = _std_by_param(s)
        for j, name in enumerate(param_names):
            stds[name].append(float(sp[j]))

    report = {"kind": "titanG_reference_sigma_scatter",
              "runs": [{"label": l, "seed": sd, "n": int(s.shape[0])}
                       for l, sd, s in runs],
              "param_names": param_names, "per_param": {}}
    print(f"\n{'param':14s} {'sig_min':>10s} {'sig_max':>10s} {'spread%':>8s}")
    for name in param_names:
        arr = np.array(stds[name])
        spread = (arr.max() - arr.min()) / arr.mean() * 100.0
        report["per_param"][name] = {
            "sigmas": stds[name], "min": float(arr.min()),
            "max": float(arr.max()), "mean": float(arr.mean()),
            "spread_pct": float(spread)}
        print(f"{name:14s} {arr.min():10.4g} {arr.max():10.4g} {spread:8.2f}")

    with open(OUT / "reference_sigma_scatter.json", "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n[scatter] report -> {OUT/'reference_sigma_scatter.json'}")
    print("[scatter] Interpretation: compare the flow/reference ratio (~1.17x, "
          "i.e. +17%) against these reference-to-reference spread% values plus "
          "the known mild NSF over-dispersion. DIAGNOSTIC ONLY — not iterated.")


if __name__ == "__main__":
    main()
