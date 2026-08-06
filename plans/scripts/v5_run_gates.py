"""Run the ratified v5 SBI gates for all three arms + the ablation comparison.

Per-arm gates via PlanetProfile/Inference/validate_sbi.py:
  - sbc     : simulation-based calibration (fresh held-out pairs via --config)
  - limits  : limiting-behavior monotonicity + prior containment
  - crosscheck (BASELINE ONLY): SBI vs the reference MCMC pkl on the same obs
              (ratification-blocking). Ablations have no reference MCMC.

Also writes the ablation comparison (D_iceIh median + corr(D,w) per arm at the
fiducial conditioning) — the headline deliverable: which observable set drives
the thick-ice pull and the D<->w degeneracy.

Reports -> validation_reports/europa_clipper_v5_{baseline,noinduction,nok2}_1m/.
Invoked by the caller with the PPcl env + PYTHONPATH/NUMBA/KMP env vars.
"""
import json, os, subprocess, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "PlanetProfile/Inference/sbi_artifacts"
CFG = ROOT / "PlanetProfile/Inference/configs"
REPORTS = ROOT / "validation_reports"
# §0.7 step-2: crosscheck targets the FRESH pooled n_eff~2000 reference (B3
# seeds 101/202/303), NOT the OLD n_eff=500 pkl (whose 1.06 km v5-v7 wander was
# a resolution artifact). Build it with plans/scripts/pool_v5_reference_neff2000.py.
POOLED_REF = Path("/tmp/b3_build/europa_clipper_v5_reference_pooled_neff2000.pkl")
OLD_REF = (ROOT / "PlanetProfile/Test/mcmc_results/Europa/Test52_seawater_v5/"
           "europa_clipper_v5_reference_result.pkl")
REF_MCMC = POOLED_REF
VALIDATE = ["python", "-m", "PlanetProfile.Inference.validate_sbi"]
SEED = "51"
# B2: request enough SBC pairs that >=500 survive the ~0.64 support rejection
# (500/0.36 ~ 1389; use 1500 for margin). The old n=300 kept only 108 pairs.
N_SBC = "1500"
# §0.7 step-3 (reviewer-ratified 2026-08-06): D_iceIh_km reference-wander floor
# = 2x the 0.18 km B3 empirical floor. Injected into the crosscheck shape-excess
# mean-shift budget (relax-only). D_iceIh ONLY; NOT applied to v6 (no measured
# wander) or other params. Passed to `crosscheck --empirical-floor`.
EMPIRICAL_FLOOR = json.dumps({"D_iceIh_km": 0.36})


def _validate_sbi_sha():
    """HEAD SHA touching validate_sbi.py, for B1's per-run provenance record."""
    try:
        r = subprocess.run(
            ["git", "log", "-1", "--format=%H", "--",
             "PlanetProfile/Inference/validate_sbi.py"],
            cwd=ROOT, capture_output=True, text=True)
        return r.stdout.strip() or None
    except Exception:
        return None

# The |Im k2| channel the limits sweep walks; every OTHER observable must be
# pinned via --fixed-obs or validate_sbi raises (fixed values required for
# non-swept observables). Aliases mirror validate_sbi._IM_K2_ALIASES.
_IM_K2_ALIASES = ("Im_k2", "im_k2", "abs_Im_k2")


def _fixed_obs_json(cfgpath):
    """Central values for every observable except the swept |Im k2| channel,
    as a JSON string for `limits --fixed-obs`. Without this the limits gate
    errors on all non-swept channels (reviewer finding, 2026-07-21)."""
    obs = json.load(open(cfgpath))["observables"]
    swept = next((a for a in _IM_K2_ALIASES if a in obs), None)
    fixed = {k: float(v[0]) for k, v in obs.items() if k != swept}
    return json.dumps(fixed)

ARMS = {
    "baseline":    ("europa_clipper_v5_geodesy_11D_posterior_1m.pt",
                    "europa_clipper_v5_geodesy_11D.json"),
    "noinduction": ("europa_clipper_v5_noinduction_7obs_posterior_1m.pt",
                    "europa_clipper_v5_noinduction_7obs.json"),
    "nok2":        ("europa_clipper_v5_nok2_17obs_posterior_1m.pt",
                    "europa_clipper_v5_nok2_17obs.json"),
}


def run(cmd, logpath):
    print(f"\n[gate] $ {' '.join(str(c) for c in cmd)}", flush=True)
    with open(logpath, "w") as f:
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=ROOT)
    print(f"[gate]   exit={p.returncode} log={logpath}", flush=True)
    return p.returncode


def _ablation_comparison(seed=51, n_draw=20000):
    """Headline deliverable: per arm, condition the flow on that arm's fiducial
    observable central values and report the D_iceIh posterior median + the
    (D_iceIh, log10_w) correlation — the quantity the ablations exist to
    attribute to observable groups. Not an automated pass/fail gate; a reported
    diagnostic (reviewer finding, 2026-07-21)."""
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    out = {}
    for tag, (artname, cfgname) in ARMS.items():
        art = ART / artname
        cfg = CFG / cfgname
        if not art.exists():
            out[tag] = {"error": f"artifact missing: {art}"}
            continue
        obs = json.load(open(cfg))["observables"]
        x_obs = {k: float(v[0]) for k, v in obs.items()}
        runner = SBIRunner.load_artifact(str(art))
        names = list(runner.param_names)
        if "D_iceIh_km" not in names or "log10_wOcean_ppt" not in names:
            out[tag] = {"error": f"expected params absent in {names}"}
            continue
        s = runner.sample_posterior(x_obs, n_samples=n_draw, seed=seed)
        di = s[:, names.index("D_iceIh_km")]
        lw = s[:, names.index("log10_wOcean_ppt")]
        out[tag] = {
            "n_draw": int(s.shape[0]),
            "D_iceIh_median_km": float(np.median(di)),
            "D_iceIh_q16_km": float(np.percentile(di, 16)),
            "D_iceIh_q84_km": float(np.percentile(di, 84)),
            "corr_D_logw": float(np.corrcoef(di, lw)[0, 1]),
            "logw_median": float(np.median(lw)),
        }
        print(f"[ablation] {tag}: D_iceIh med={out[tag]['D_iceIh_median_km']:.2f} km "
              f"corr(D,logw)={out[tag]['corr_D_logw']:.3f}")
    return out


def main():
    if not REF_MCMC.exists():
        raise SystemExit(
            f"[gate] pooled fresh reference missing: {REF_MCMC}\n"
            f"       build it first: python plans/scripts/pool_v5_reference_neff2000.py")
    summary = {}
    for tag, (artname, cfgname) in ARMS.items():
        art = ART / artname
        cfg = CFG / cfgname
        outdir = REPORTS / f"europa_clipper_v5_{tag}_1m"
        outdir.mkdir(parents=True, exist_ok=True)
        rc = {}

        # SBC: generate fresh held-out pairs from the config forward model.
        rc["sbc"] = run(VALIDATE + [
            "sbc", "--artifact", str(art), "--config", str(cfg),
            "--n-sbc", N_SBC, "--num-posterior-samples", "1000",
            "--seed", SEED, "--output-dir", str(outdir / "sbc")], outdir / "sbc.log")

        # limits: |Im k2| sweep only makes sense where Im_k2 is an observable
        # (baseline + nok2 keep it; noinduction keeps it too — all three have k2
        #  EXCEPT nok2 drops k2. So limits' default Im-k2 sweep applies to
        #  baseline + noinduction; nok2 has no k2 channel to sweep.)
        if tag != "nok2":
            rc["limits"] = run(VALIDATE + [
                "limits", "--artifact", str(art),
                "--fixed-obs", _fixed_obs_json(cfg),
                "--seed", SEED, "--output-dir", str(outdir / "limits")],
                outdir / "limits.log")
        else:
            rc["limits"] = "skipped (no k2 channel to sweep)"

        # crosscheck: baseline only (only arm with a reference MCMC).
        if tag == "baseline":
            rc["crosscheck"] = run(VALIDATE + [
                "crosscheck", "--artifact", str(art), "--mcmc", str(REF_MCMC),
                "--empirical-floor", EMPIRICAL_FLOOR,
                "--seed", SEED, "--output-dir", str(outdir / "crosscheck")],
                outdir / "crosscheck.log")
        else:
            rc["crosscheck"] = "n/a (no reference MCMC for this ablation)"

        summary[tag] = rc

    print("\n[gate] ablation comparison (D_iceIh + D<->w degeneracy per arm) ...")
    ablation = _ablation_comparison(seed=int(SEED))

    with open(REPORTS / "v5_gate_summary.json", "w") as f:
        json.dump({"seed": int(SEED), "arms": summary,
                   "ablation_comparison": ablation,
                   "ref_mcmc": str(REF_MCMC),
                   "ref_mcmc_kind": "fresh pooled n_eff~2000 (B3 seeds 101/202/303)",
                   "n_sbc_requested": int(N_SBC),
                   "crosscheck_empirical_floor": json.loads(EMPIRICAL_FLOOR),
                   "derived_params_sbc_na": {
                       "Tb_K": ("N/A (derived from sampled D_iceIh_km via "
                                "per-salinity PCHIP inversion; not a sampled "
                                "parameter, no SBC rank-uniformity test "
                                "applicable) — reviewer Item-2 2026-08-06")},
                   "validate_sbi_sha": _validate_sbi_sha()}, f, indent=2)
    print(f"\n[gate] summary -> {REPORTS/'v5_gate_summary.json'}")
    print(json.dumps({"arms": summary, "ablation_comparison": ablation}, indent=2))


if __name__ == "__main__":
    main()
