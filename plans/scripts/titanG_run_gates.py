"""Run the ratified SBI gates for the Titan free-gravity (no-ocean) artifact.

SINGLE-ARM (no ablation comparison — Titan has no induction/h2 siblings).
Gates via PlanetProfile/Inference/validate_sbi.py:
  - sbc        : simulation-based calibration (fresh held-out pairs via --config)
  - limits     : limiting-behavior monotonicity + prior containment
                 (|Im k2| sweep; Titan keeps the k2 channel)
  - crosscheck : SBI vs the fresh Titan free-gravity reference MCMC pkl on the
                 same {C20,C22,Re_k2,Im_k2} obs (ratification-blocking)

Reports -> validation_reports/titan_freegrav_noocean_1m/.
Gates are interpreted by the scientific-reviewer, NEVER tuned to pass.

Run (after gen + train + reference MCMC complete):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/titanG_run_gates.py
"""
import json, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "PlanetProfile/Inference/sbi_artifacts/titan_freegrav_noocean_posterior_1m.pt"
CFG = ROOT / "PlanetProfile/Inference/configs/titan_freegrav_noocean.json"
REPORTS = ROOT / "validation_reports/titan_freegrav_noocean_1m"
# Reference MCMC lives OUTSIDE Test/ (repo rule: no Test/ writes without
# permission) -- under the campaign's own validation_reports tree.
REF_MCMC = (ROOT / "validation_reports/titan_freegrav_noocean_1m/reference/"
            "titan_freegrav_noocean_reference_result.pkl")
VALIDATE = ["python", "-m", "PlanetProfile.Inference.validate_sbi"]
SEED = "71"

# The |Im k2| channel the limits sweep walks; every OTHER observable must be
# pinned via --fixed-obs. Aliases mirror validate_sbi._IM_K2_ALIASES.
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
    REPORTS.mkdir(parents=True, exist_ok=True)
    if not ART.exists():
        raise SystemExit(f"[gate] artifact missing: {ART} — run titanG_train_all.py first")
    rc = {}

    rc["sbc"] = run(VALIDATE + [
        "sbc", "--artifact", str(ART), "--config", str(CFG),
        "--n-sbc", "300", "--num-posterior-samples", "1000",
        "--seed", SEED, "--output-dir", str(REPORTS / "sbc")], REPORTS / "sbc.log")

    rc["limits"] = run(VALIDATE + [
        "limits", "--artifact", str(ART),
        "--fixed-obs", _fixed_obs_json(CFG),
        "--seed", SEED, "--output-dir", str(REPORTS / "limits")],
        REPORTS / "limits.log")

    if REF_MCMC.exists():
        rc["crosscheck"] = run(VALIDATE + [
            "crosscheck", "--artifact", str(ART), "--mcmc", str(REF_MCMC),
            "--seed", SEED, "--output-dir", str(REPORTS / "crosscheck")],
            REPORTS / "crosscheck.log")
    else:
        rc["crosscheck"] = f"MISSING reference MCMC: {REF_MCMC} — run titanG_reference_mcmc.py"

    with open(REPORTS / "titanG_gate_summary.json", "w") as f:
        json.dump({"seed": int(SEED), "artifact": str(ART), "config": str(CFG),
                   "ref_mcmc": str(REF_MCMC), "gates": rc}, f, indent=2)
    print(f"\n[gate] summary -> {REPORTS/'titanG_gate_summary.json'}")
    print(json.dumps(rc, indent=2))


if __name__ == "__main__":
    main()
