"""Run the v7 open-|Ae| SBI gates (single arm).

Gates via PlanetProfile/Inference/validate_sbi.py:
  - sbc        : simulation-based calibration (fresh held-out pairs via --config)
  - limits     : limiting-behavior monotonicity + prior containment (|Im k2| sweep)
  - crosscheck : SBI vs the v7 reference MCMC pkl on the same obs
                 (ratification-blocking). MUST use the open-bounds reference.

Gates are interpreted, never tuned to pass. Same pass criteria as v5/v6
(multiplicity-aware, SHAPE_EXCESS residual clause).

Known v7-specific risks to watch (spec):
  (a) posterior mass piling at the |Ae| saturation edge for tight synthetic
      conditioning (flow edge behavior near the physical support edge);
  (b) the near-frozen corner is ~0 ocean — verify D_ocean -> 0 nodes are
      handled (they are built and legitimate for the open config).

Reports -> validation_reports/europa_clipper_v7_openae_1m/.
Run:
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/v7_run_gates.py
"""
import json, os, subprocess, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "PlanetProfile/Inference/sbi_artifacts"
CFG = ROOT / "PlanetProfile/Inference/configs"
REPORTS = ROOT / "validation_reports"
REF_MCMC = (ROOT / "PlanetProfile/Test/mcmc_results/Europa/Test54_seawater_v7/"
            "europa_clipper_v7_reference_result.pkl")
VALIDATE = ["python", "-m", "PlanetProfile.Inference.validate_sbi"]
SEED = "71"

ARTNAME = "europa_clipper_v7_openae_11D_posterior_1m.pt"
CFGNAME = "europa_clipper_v7_openae_11D.json"

_IM_K2_ALIASES = ("Im_k2", "im_k2", "abs_Im_k2")


def _fixed_obs_json(cfgpath):
    obs = json.load(open(cfgpath))["observables"]
    swept = next((a for a in _IM_K2_ALIASES if a in obs), None)
    fixed = {k: float(v[0]) for k, v in obs.items() if k != swept}
    return json.dumps(fixed)


def run(cmd, logpath):
    print(f"\n[gate] $ {' '.join(str(c) for c in cmd)}", flush=True)
    with open(logpath, "w") as f:
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=ROOT)
    print(f"[gate]   exit={p.returncode} log={logpath}", flush=True)
    return p.returncode


def main():
    art = ART / ARTNAME
    cfg = CFG / CFGNAME
    if not art.exists():
        raise SystemExit(f"[gate] artifact missing: {art} — train first")
    outdir = REPORTS / "europa_clipper_v7_openae_1m"
    outdir.mkdir(parents=True, exist_ok=True)
    rc = {}

    rc["sbc"] = run(VALIDATE + [
        "sbc", "--artifact", str(art), "--config", str(cfg),
        "--n-sbc", "300", "--num-posterior-samples", "1000",
        "--seed", SEED, "--output-dir", str(outdir / "sbc")], outdir / "sbc.log")

    rc["limits"] = run(VALIDATE + [
        "limits", "--artifact", str(art),
        "--fixed-obs", _fixed_obs_json(cfg),
        "--seed", SEED, "--output-dir", str(outdir / "limits")],
        outdir / "limits.log")

    if REF_MCMC.exists():
        rc["crosscheck"] = run(VALIDATE + [
            "crosscheck", "--artifact", str(art), "--mcmc", str(REF_MCMC),
            "--seed", SEED, "--output-dir", str(outdir / "crosscheck")],
            outdir / "crosscheck.log")
    else:
        rc["crosscheck"] = f"MISSING reference MCMC: {REF_MCMC} (run v7_reference_mcmc.py)"

    with open(REPORTS / "v7_gate_summary.json", "w") as f:
        json.dump({"seed": int(SEED), "artifact": str(art),
                   "config": str(cfg), "gates": rc,
                   "ref_mcmc": str(REF_MCMC)}, f, indent=2)
    print(f"\n[gate] summary -> {REPORTS/'v7_gate_summary.json'}")
    print(json.dumps({"gates": rc}, indent=2))


if __name__ == "__main__":
    main()
