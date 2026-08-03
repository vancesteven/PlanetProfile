"""Run the ratified SBI gates for the Titan free-gravity NH3 JOINT artifact.

SINGLE-ARM (no ablation comparison — Titan has no induction/h2 siblings).
Gates via PlanetProfile/Inference/validate_sbi.py:
  - sbc        : simulation-based calibration (fresh held-out pairs via --config)
  - limits     : limiting-behavior monotonicity + prior containment
                 (|Im k2| sweep; Titan keeps the k2 channel)
  - crosscheck : SBI vs the fresh Titan NH3 joint reference MCMC pkl on the
                 same {C20,C22,Re_k2,Im_k2} obs (ratification-blocking)

Reports -> validation_reports/titan_freegrav_nh3_1m/.
Gates are interpreted by the scientific-reviewer, NEVER tuned to pass.

Run (after gen + train + reference MCMC complete):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/titanG_nh3_run_gates.py
"""
import json, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "PlanetProfile/Inference/sbi_artifacts/titan_freegrav_nh3_posterior_1m.pt"
CFG = ROOT / "PlanetProfile/Inference/configs/test54_titan_nh3_freegrav.json"
REPORTS = ROOT / "validation_reports/titan_freegrav_nh3_1m"
REF_MCMC = (ROOT / "validation_reports/titan_freegrav_nh3_1m/reference/"
            "titan_freegrav_nh3_reference_result.pkl")
# Training dataset for the posterior-predictive + prior-predictive-interior
# diagnostic (built by titanG_nh3_gen_dataset.py in /tmp). Optional: the PPC
# step is skipped with a note if it is absent.
DATASET = Path("/tmp/titanG_build/datasets/titanG_nh3_1m.npz")
PPC = [str(ROOT / "plans/scripts/titanG_ppc_interior_check.py")]
VALIDATE = ["python", "-m", "PlanetProfile.Inference.validate_sbi"]
SEED = "72"

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
        raise SystemExit(f"[gate] artifact missing: {ART} — run titanG_nh3_train_all.py first")
    rc = {}

    # SBC uses >=400 pairs per scientific-reviewer pilot adjudication (2026-08-03,
    # validation item 3): n=212 had modest power to detect the ~1.2x pilot
    # over-dispersion; the production SBC widens power and watches log10_wOcean_ppt
    # (the only param near the raw-KS threshold at pilot, p=0.058).
    rc["sbc"] = run(VALIDATE + [
        "sbc", "--artifact", str(ART), "--config", str(CFG),
        "--n-sbc", "500", "--num-posterior-samples", "1000",
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
        rc["crosscheck"] = f"MISSING reference MCMC: {REF_MCMC} — run titanG_nh3_reference_mcmc.py"

    # Posterior-predictive coverage + prior-predictive interior check.
    # DIAGNOSTIC (not a gate — never tuned to pass). Answers the user's
    # "k2 has an extraordinarily broad range" observation quantitatively:
    # does the conditioned posterior reproduce the data, and is x_obs
    # interior to the prior-predictive envelope? Skipped (noted) if the
    # training .npz is absent (it lives in /tmp, may be cleared).
    if DATASET.exists():
        rc["ppc"] = run(["python"] + PPC + [
            "--artifact", str(ART), "--config", str(CFG),
            "--dataset", str(DATASET),
            "--output-dir", str(REPORTS / "ppc"),
            "--n-post", "4000", "--seed", SEED], REPORTS / "ppc.log")
    else:
        rc["ppc"] = f"SKIPPED (diagnostic): training dataset absent at {DATASET}"

    with open(REPORTS / "titanG_nh3_gate_summary.json", "w") as f:
        json.dump({"seed": int(SEED), "artifact": str(ART), "config": str(CFG),
                   "ref_mcmc": str(REF_MCMC), "gates": rc}, f, indent=2)
    print(f"\n[gate] summary -> {REPORTS/'titanG_nh3_gate_summary.json'}")
    print(json.dumps(rc, indent=2))


if __name__ == "__main__":
    main()
