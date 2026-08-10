"""Run the ratified SBI gates for the Titan free-gravity MgSO4 / NaCl JOINT artifacts — parameterized.

Phase B of the Titan free-gravity campaign; §0.16 post-training gate sequence
steps 2-4 (step 1 = pushforward + tidal-quarantine, run separately via
titanG_ppc_interior_check.py). Parameterized analogue of titanG_run_gates.py
(no-ocean) for the MgSO4/NaCl joint configs.

Gates via PlanetProfile/Inference/validate_sbi.py:
  - sbc        : simulation-based calibration (fresh held-out pairs via --config).
                 §0.16 preregisters n_sbc=1500 (NOT the no-ocean driver's 300).
  - limits     : limiting-behavior monotonicity + prior containment (|Im k2| sweep;
                 the tidal channel is QUARANTINED split-status — limits on the tidal
                 sweep are advisory, gravity/Re_k2 containment is the binding read).
  - crosscheck : SBI vs the POOLED B3 multi-seed reference MCMC on the same
                 {C20,C22,Re_k2,Im_k2} obs.

Reports -> validation_reports/titan_freegrav_<comp>_1m/. Gates are interpreted by
the scientific-reviewer, NEVER tuned to pass.

Run (after gen + train + reference MCMC + pushforward complete):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
    python plans/scripts/titanG_ocean_run_gates.py --comp NaCl
"""
import argparse, json, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
COMPS = {
    "MgSO4": {"cfg": "PlanetProfile/Inference/configs/test54_titan_mgso4_freegrav.json",
              "tag": "mgso4", "seed": "73"},
    "NaCl":  {"cfg": "PlanetProfile/Inference/configs/test54_titan_nacl_freegrav.json",
              "tag": "nacl", "seed": "74"},
}
VALIDATE = ["python", "-m", "PlanetProfile.Inference.validate_sbi"]
N_SBC = "1500"  # §0.16 preregistered
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--comp", required=True, choices=list(COMPS))
    ap.add_argument("--gates", nargs="+", default=["sbc", "limits", "crosscheck"],
                    choices=["sbc", "limits", "crosscheck"])
    args = ap.parse_args()
    spec = COMPS[args.comp]
    tag = spec["tag"]
    ART = ROOT / f"PlanetProfile/Inference/sbi_artifacts/titan_freegrav_{tag}_posterior_1m.pt"
    CFG = ROOT / spec["cfg"]
    REPORTS = ROOT / f"validation_reports/titan_freegrav_{tag}_1m"
    REF_MCMC = REPORTS / "reference" / f"titan_freegrav_{tag}_reference_pooled.pkl"
    SEED = spec["seed"]
    REPORTS.mkdir(parents=True, exist_ok=True)
    assert ART.exists(), f"missing artifact {ART}"
    assert CFG.exists(), f"missing config {CFG}"

    rc = {}
    if "sbc" in args.gates:
        rc["sbc"] = run(VALIDATE + [
            "sbc", "--artifact", str(ART), "--config", str(CFG),
            "--n-sbc", N_SBC, "--num-posterior-samples", "1000",
            "--seed", SEED, "--output-dir", str(REPORTS / "sbc")], REPORTS / "sbc.log")
    if "limits" in args.gates:
        rc["limits"] = run(VALIDATE + [
            "limits", "--artifact", str(ART),
            "--fixed-obs", _fixed_obs_json(CFG),
            "--seed", SEED, "--output-dir", str(REPORTS / "limits")],
            REPORTS / "limits.log")
    if "crosscheck" in args.gates:
        if REF_MCMC.exists():
            rc["crosscheck"] = run(VALIDATE + [
                "crosscheck", "--artifact", str(ART), "--mcmc", str(REF_MCMC),
                "--seed", SEED, "--output-dir", str(REPORTS / "crosscheck")],
                REPORTS / "crosscheck.log")
        else:
            rc["crosscheck"] = f"MISSING reference MCMC: {REF_MCMC}"

    print(f"\n[gate] {args.comp} exit codes: {rc}")
    manifest = {"kind": f"titanG_{tag}_gate_runner", "comp": args.comp,
                "artifact": str(ART), "config": spec["cfg"],
                "reference_mcmc": str(REF_MCMC), "n_sbc": N_SBC,
                "seed": SEED, "gate_exit_codes": {k: v for k, v in rc.items()},
                "tidal_split_status": ("tidal k2 sector QUARANTINED (pushforward "
                    "step-1 verdict); limits on the |Im k2| sweep are advisory, "
                    "gravity+Re_k2 containment is the binding read.")}
    with open(REPORTS / "gate_run_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[gate] manifest -> {REPORTS / 'gate_run_manifest.json'}")


if __name__ == "__main__":
    main()
