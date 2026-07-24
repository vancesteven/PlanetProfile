"""Limits non-monotonicity attribution (scientific-reviewer a2d7f121 item 3).

The SBI limits gate fails: log10_eta_Ih median vs Im_k2 rises then falls
(low-Im_k2 end 0.05->0.15 INCREASES, which is backwards for a
dissipation->viscosity mapping). Question: is this a FLOW artifact or genuine
forward-model PHYSICS (a degeneracy where multiple eta_Ih give the same low
Im_k2 at fixed Re_k2)?

Test: run the SAME sweep on the REFERENCE MCMC — condition the pocoMC posterior
on each Im_k2 value (all other observables fixed at the SBI sweep's fixed_obs)
and record the log10_eta_Ih posterior median. If the MCMC is ALSO non-monotonic
at the low-Im_k2 end, it is physics (and the gate's monotonicity premise is
wrong for Titan). If the MCMC is monotone-decreasing and only the flow is not,
it is a flow limitation on the Im_k2 channel.

This is a targeted 6-point MCMC sweep (~1 h). DIAGNOSTIC — reports the
attribution; does not change the gate.

Run (separate process):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/titanG_limits_attribution.py
"""
import copy, json, time
from pathlib import Path
import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.mcmc_runner import MCMCRunner

ROOT = Path(__file__).resolve().parents[2]
CFG = ROOT / "PlanetProfile/Inference/configs/titan_freegrav_noocean.json"
OUT = ROOT / "validation_reports/titan_freegrav_noocean_1m/limits_attribution"
TMP = Path("/tmp/titanG_build/limits_attribution")

SWEEP_OBS = "Im_k2"
SWEEP_VALUES = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
FIXED_OBS = {"C20": -3.3511e-05, "C22": 1.0107e-05, "Re_k2": 0.608}
MONOTONE_PARAM = "log10_eta_Ih"
SEED = 71
# SBI medians for side-by-side (from limits_report.json this session).
SBI_MEDIANS = [12.196, 12.546, 12.607, 12.452, 11.529, 10.961]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
    base = InferenceConfig.from_json(str(CFG))
    pidx = list(base.param_space.keys()).index(MONOTONE_PARAM)

    medians = []
    grid = []
    for val in SWEEP_VALUES:
        # Rebuild a config whose observable centrals match the sweep point:
        # sweep Im_k2, pin the others to FIXED_OBS, keep each observable's sigma.
        cfg = copy.deepcopy(base)
        obs = dict(cfg.observables)
        for k, v in FIXED_OBS.items():
            obs[k] = [float(v), float(obs[k][1])]
        obs[SWEEP_OBS] = [float(val), float(obs[SWEEP_OBS][1])]
        cfg.observables = obs

        runner = MCMCRunner(cfg)
        runner.random_state = SEED
        t0 = time.time()
        print(f"[attr] MCMC at {SWEEP_OBS}={val} ...", flush=True)
        res = runner.run(progress_jsonl_path=str(TMP / f"prog_{val}.jsonl"))
        dt = time.time() - t0
        s = np.asarray(res.samples, np.float64)
        med = float(np.median(s[:, pidx]))
        medians.append(med)
        grid.append({"sweep_value": val, f"{MONOTONE_PARAM}_median": med,
                     "n_samples": int(s.shape[0])})
        print(f"[attr]   {SWEEP_OBS}={val}: {MONOTONE_PARAM}_median={med:.3f} "
              f"({dt/60:.1f} min, n={s.shape[0]})")

    diffs = [medians[i + 1] - medians[i] for i in range(len(medians) - 1)]
    n_increasing = sum(1 for d in diffs if d > 0.01)
    low_end_increasing = medians[2] > medians[0]  # 0.05 -> 0.15 rises?
    mcmc_monotone_dec = all(d <= 0.01 for d in diffs)

    report = {
        "kind": "titanG_limits_attribution",
        "sweep_obs": SWEEP_OBS, "sweep_values": SWEEP_VALUES,
        "fixed_obs": FIXED_OBS, "monotone_param": MONOTONE_PARAM,
        "seed": SEED, "grid": grid,
        "mcmc_medians": medians, "mcmc_diffs": diffs,
        "sbi_medians": SBI_MEDIANS,
        "mcmc_low_end_increasing_0p05_to_0p15": bool(low_end_increasing),
        "mcmc_monotone_decreasing": bool(mcmc_monotone_dec),
        "attribution": (
            "PHYSICS (reference also non-monotone at low Im_k2; gate "
            "monotonicity premise wrong for Titan)" if low_end_increasing
            else "FLOW (reference monotone-decreasing; flow non-monotonicity "
                 "is an Im_k2-channel flow limitation)"),
    }
    with open(OUT / "limits_attribution.json", "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n[attr] MCMC medians: {[round(m,3) for m in medians]}")
    print(f"[attr] SBI  medians: {SBI_MEDIANS}")
    print(f"[attr] low-end (0.05->0.15) increasing in MCMC? {low_end_increasing}")
    print(f"[attr] ATTRIBUTION: {report['attribution']}")
    print(f"[attr] report -> {OUT/'limits_attribution.json'}")


if __name__ == "__main__":
    main()
