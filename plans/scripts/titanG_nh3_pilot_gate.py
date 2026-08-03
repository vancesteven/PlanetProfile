"""Reviewer item-4 PILOT gate for the Titan free-gravity NH3 JOINT artifact.

Binding directive (scientific-reviewer afc74fe6, 2026-08-02):
  "Run SBC + SBI-vs-MCMC crosscheck on a pilot (e.g. 50-100k) *after* the
   [support-guard] fix, before committing the 1M — Phase A precedent says the
   box is what makes these pass; verify empirically rather than by analogy."

This trains a flow on a PILOT SUBSET of the already-generated 1M NH3 dataset,
then runs the two ratification-relevant gates (SBC + crosscheck vs the fresh
NH3 joint reference MCMC). If both look healthy, the full 1M flow train
(titanG_nh3_train_all.py) is justified; if not, we stop before wasting the
long train.

Determinism: pilot subset is the FIRST n_pilot rows of the shuffled-at-gen
dataset (gen already applied seed=72); train seed = 72. The pilot artifact is
a DISTINCT filename so it never shadows the production 1M artifact.

libomp hazard: torch-only process (no PlanetProfile forward-model import here;
the dataset is read from npz). Safe to run standalone.

Run (after gen + reference MCMC complete):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/titanG_nh3_pilot_gate.py
"""
import argparse, json, os, subprocess, time
from pathlib import Path

import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.sbi_runner import SBIRunner

ROOT = Path(__file__).resolve().parents[2]
DSDIR = Path("/tmp/titanG_build/datasets")
NPZ = DSDIR / "titanG_nh3_1m.npz"
PILOT_DIR = Path("/tmp/titanG_build/pilot")
CFG = ROOT / "PlanetProfile/Inference/configs/test54_titan_nh3_freegrav.json"
REF_MCMC = (ROOT / "validation_reports/titan_freegrav_nh3_1m/reference/"
            "titan_freegrav_nh3_reference_result.pkl")
REPORTS = ROOT / "validation_reports/titan_freegrav_nh3_1m/pilot"
PILOT_ART = PILOT_DIR / "titan_freegrav_nh3_posterior_pilot.pt"
VALIDATE = ["python", "-m", "PlanetProfile.Inference.validate_sbi"]
TRAIN_SEED = 72
SEED = "72"

_IM_K2_ALIASES = ("Im_k2", "im_k2", "abs_Im_k2")


def _fixed_obs_json(cfgpath):
    obs = json.load(open(cfgpath))["observables"]
    swept = next((a for a in _IM_K2_ALIASES if a in obs), None)
    fixed = {k: float(v[0]) for k, v in obs.items() if k != swept}
    return json.dumps(fixed)


def run(cmd, logpath):
    print(f"\n[pilot-gate] $ {' '.join(str(c) for c in cmd)}", flush=True)
    with open(logpath, "w") as f:
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=ROOT)
    print(f"[pilot-gate]   exit={p.returncode} log={logpath}", flush=True)
    return p.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-pilot", type=int, default=100000)
    args = ap.parse_args()

    PILOT_DIR.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    if not NPZ.exists():
        raise SystemExit(f"[pilot-gate] missing dataset {NPZ}")

    with np.load(NPZ, allow_pickle=True) as d:
        theta = np.asarray(d["theta"], np.float64)
        x = np.asarray(d["x"], np.float64)
    n = min(args.n_pilot, len(theta))
    theta_p, x_p = theta[:n], x[:n]

    cfg = json.load(open(CFG)); cfg["mode"] = "sbi"
    runner = SBIRunner(InferenceConfig.from_dict(cfg))
    assert x_p.shape[1] == len(runner.obs_names), \
        f"x cols {x_p.shape[1]} != obs {len(runner.obs_names)}"

    print(f"[pilot-gate] PILOT train theta={theta_p.shape} x={x_p.shape} "
          f"(of {len(theta)} available) seed={TRAIN_SEED} nsf")
    t0 = time.time()
    runner.train(theta_p, x_p, seed=TRAIN_SEED, density_estimator="nsf")
    runner.save_artifact(str(PILOT_ART))
    print(f"[pilot-gate] pilot flow trained in {(time.time()-t0)/60:.1f} min "
          f"-> {PILOT_ART}")

    rc = {"n_pilot": n}
    rc["sbc"] = run(VALIDATE + [
        "sbc", "--artifact", str(PILOT_ART), "--config", str(CFG),
        "--n-sbc", "300", "--num-posterior-samples", "1000",
        "--seed", SEED, "--output-dir", str(REPORTS / "sbc")],
        REPORTS / "sbc.log")

    if REF_MCMC.exists():
        rc["crosscheck"] = run(VALIDATE + [
            "crosscheck", "--artifact", str(PILOT_ART), "--mcmc", str(REF_MCMC),
            "--seed", SEED, "--output-dir", str(REPORTS / "crosscheck")],
            REPORTS / "crosscheck.log")
    else:
        rc["crosscheck"] = f"MISSING reference MCMC: {REF_MCMC}"

    with open(REPORTS / "titanG_nh3_pilot_gate_summary.json", "w") as f:
        json.dump({"seed": int(SEED), "n_pilot": n, "artifact": str(PILOT_ART),
                   "config": str(CFG), "ref_mcmc": str(REF_MCMC),
                   "gates": rc}, f, indent=2)
    print(f"\n[pilot-gate] summary -> {REPORTS/'titanG_nh3_pilot_gate_summary.json'}")
    print(json.dumps(rc, indent=2))


if __name__ == "__main__":
    main()
