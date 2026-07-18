#!/usr/bin/env python
"""Clipper v2 same-version sampling reproduction check (Machine B, 2026-07-18).

Last item of the v2 ratification package (handoff 2026-07-18 queue). Machine A
started and killed this because SBIRunner(config) triggers the 3-frequency Ae
grid precompute. This version avoids that: load_artifact() is the pure-sampling
path (no config, no cache, no Ae grid). The conditioning x_obs and reference
sbi_mean/sbi_std are already committed in crosscheck_report.json.

Two checks:
  1. Same-version reproduction: draw 4000 samples at config-central x_obs
     (seed 42, reject_outside_prior=True = the public/GUI contract path) and
     compare per-param mean/sigma against the committed flow stats. Same torch
     2.8 / sbi 0.26.1 / seed -> expect near-exact.
     Gate: |dmean| < 0.1*sigma_ref ; sigma_new/sigma_ref in [0.9, 1.1].
  2. Signed-Im semantics probe. Flipping a channel's sign conditions the flow
     on off-manifold data, so (per validate_sbi's off-manifold convention,
     lines 903/1012) these draws use reject_outside_prior=False — otherwise the
     posterior legitimately spreads outside the prior box and rejection
     sampling stalls at 0% acceptance (this is what killed the first attempt).
       - flip a *signed* Bind imag channel -> posterior MUST move;
       - condition on +v vs -v of the *abs*-folded Im_k2 -> _x_obs_vector's
         abs() maps both to |v|, so same seed -> byte-identical draws.

The spline resample guard in sample_posterior absorbs the intermittent nflows
discriminant assert (guard active on the reject=True path only).
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
warnings.filterwarnings("ignore")

ARTIFACT = "PlanetProfile/Inference/sbi_artifacts/europa_seawater_andrade_clipper_v2.pt"
REPORT = ("PlanetProfile/Inference/sbi_artifacts/validation_reports/"
          "europa_clipper_v2_1m/crosscheck_report.json")

N_SAMPLES = 4000
SEED = 42
DMEAN_TOL_FRAC = 0.1
SIGMA_RATIO_LO = 0.9
SIGMA_RATIO_HI = 1.1


def main():
    from PlanetProfile.Inference.sbi_runner import SBIRunner

    ref = json.loads(Path(REPORT).read_text())
    x_obs = dict(ref["x_obs"])
    per_param_ref = {p["param"]: p for p in ref["results"]["per_parameter"]}

    runner = SBIRunner.load_artifact(ARTIFACT)
    param_names = list(runner.param_names)
    print(f"loaded {Path(ARTIFACT).name}: params={param_names}")
    print(f"conditioning x_obs has {len(x_obs)} observables\n")

    # --- Check 1: same-version reproduction (public reject=True path) ------
    samples = runner.sample_posterior(x_obs, n_samples=N_SAMPLES, seed=SEED,
                                      reject_outside_prior=True)
    assert samples.shape == (N_SAMPLES, len(param_names)), samples.shape
    assert np.isfinite(samples).all(), "non-finite draws"

    print("=== Check 1: same-version reproduction (seed 42, reject=True) ===")
    print(f"{'param':<16}{'new_mean':>14}{'ref_mean':>14}{'dmean':>12}"
          f"{'tol':>12}{'sig_ratio':>11}  pass")
    repro = []
    all_pass = True
    for i, p in enumerate(param_names):
        col = samples[:, i]
        new_mean = float(col.mean())
        new_std = float(col.std(ddof=1))
        r = per_param_ref[p]
        ref_mean = r["sbi_mean"]
        ref_std = r["sbi_std"]
        dmean = new_mean - ref_mean
        tol = DMEAN_TOL_FRAC * ref_std
        sig_ratio = new_std / ref_std if ref_std else float("nan")
        mean_ok = abs(dmean) < tol
        sig_ok = SIGMA_RATIO_LO <= sig_ratio <= SIGMA_RATIO_HI
        ok = mean_ok and sig_ok
        all_pass = all_pass and ok
        repro.append({
            "param": p, "new_mean": new_mean, "ref_mean": ref_mean,
            "new_std": new_std, "ref_std": ref_std, "dmean": dmean,
            "dmean_tol": tol, "sigma_ratio": sig_ratio,
            "mean_pass": mean_ok, "sigma_pass": sig_ok, "pass": ok,
        })
        print(f"{p:<16}{new_mean:>14.5g}{ref_mean:>14.5g}{dmean:>12.4g}"
              f"{tol:>12.4g}{sig_ratio:>11.4f}  {'OK' if ok else 'FAIL'}")
    print(f"\nCheck 1 verdict: {'PASS' if all_pass else 'FAIL'}\n")

    # --- Check 2: signed-Im semantics probe (off-manifold => reject=False) --
    signed_ch = "Bind_synodic_x_imag"
    abs_ch = "Im_k2"
    assert runner.channel_conventions.get(signed_ch) == "signed", signed_ch
    assert runner.channel_conventions.get(abs_ch) == "abs", abs_ch

    def draw(xd, seed=SEED):
        return runner.sample_posterior(xd, n_samples=N_SAMPLES, seed=seed,
                                       reject_outside_prior=False)

    base = draw(x_obs)
    base_mean = base.mean(axis=0)
    base_std = base.std(axis=0, ddof=1)

    # (a) flip a signed Bind imag channel -> expect movement
    x_signed = dict(x_obs); x_signed[signed_ch] = -x_signed[signed_ch]
    s_signed = draw(x_signed)
    shift_signed = np.abs(s_signed.mean(axis=0) - base_mean) / base_std

    # (b) abs-fold: condition on +v vs -v of Im_k2 -> abs() => identical draws.
    #     (x_obs Im_k2 is 0.0 here, so a bare sign flip is a vacuous no-op;
    #      use a nonzero synthetic |Im_k2| and compare the two signs.)
    v = 0.08
    xp = dict(x_obs); xp[abs_ch] = +v
    xm = dict(x_obs); xm[abs_ch] = -v
    sp = draw(xp)
    sm = draw(xm)
    shift_abs = np.abs(sp.mean(axis=0) - sm.mean(axis=0)) / base_std
    abs_probe_desc = (f"Im_k2 = +{v} vs -{v} (abs-fold -> identical |Im_k2|); "
                      f"same seed + reject=False -> byte-identical draws expected")
    identical = bool(np.array_equal(sp, sm))

    # (c) Im_h2 abs-fold probe. Im_h2 is marked 'abs' in channel_conventions
    #     and folded to |Im_h2| during training, but _x_obs_vector only folds
    #     the Im_k2 aliases (sbi_runner.py:674) -> a signed Im_h2 conditions
    #     the flow off-manifold. This probe (+v vs -v) exposes that gap: if the
    #     draws are NOT byte-identical the fold is missing for Im_h2.
    h2_ch = "Im_h2"
    h2_conv = runner.channel_conventions.get(h2_ch)
    xhp = dict(x_obs); xhp[h2_ch] = +v
    xhm = dict(x_obs); xhm[h2_ch] = -v
    shp = draw(xhp)
    shm = draw(xhm)
    h2_identical = bool(np.array_equal(shp, shm))
    h2_shift = float((np.abs(shp.mean(axis=0) - shm.mean(axis=0)) / base_std).max())
    h2_fold_ok = (h2_conv == "abs") and h2_identical

    print("=== Check 2: signed-Im semantics probe (reject_outside_prior=False) ===")
    print(f"signed channel flipped: {signed_ch} "
          f"({x_obs[signed_ch]:+.3f} -> {-x_obs[signed_ch]:+.3f})")
    print(f"abs probe: {abs_probe_desc}\n")
    print(f"{'param':<16}{'signed_shift':>14}{'abs_shift':>14}  (units: sigma)")
    for i, p in enumerate(param_names):
        print(f"{p:<16}{shift_signed[i]:>14.4f}{shift_abs[i]:>14.3e}")

    max_signed = float(shift_signed.max())
    max_abs = float(shift_abs.max())
    signed_moves = max_signed > 0.1
    abs_static = identical and max_abs < 1e-9
    probe_pass = signed_moves and abs_static
    print(f"\nmax signed-flip shift: {max_signed:.4f} sigma "
          f"(expect > 0.1: {'OK' if signed_moves else 'FAIL'})")
    print(f"Im_k2 abs-fold +v/-v draws byte-identical: {identical} "
          f"(max shift {max_abs:.2e} sigma; expect exactly 0: "
          f"{'OK' if abs_static else 'FAIL'})")
    print(f"Im_h2 abs-fold +v/-v draws byte-identical: {h2_identical} "
          f"(conv={h2_conv}, max shift {h2_shift:.2e} sigma; "
          f"{'OK' if h2_fold_ok else 'GAP: _x_obs_vector folds Im_k2 only'})")
    print(f"\nCheck 2 verdict: {'PASS' if probe_pass else 'FAIL'}\n")

    report = {
        "check": "v2_same_version_sampling_reproduction",
        "machine": "B",
        "date": "2026-07-18",
        "artifact": Path(ARTIFACT).name,
        "reference_report": "crosscheck_report.json",
        "n_samples": N_SAMPLES,
        "seed": SEED,
        "sampling_path": "SBIRunner.load_artifact (no config, no Ae precompute)",
        "spline_guard_active_on_reject_true_path": True,
        "reproduction": {
            "gate": "|dmean| < 0.1*sigma_ref AND sigma_ratio in [0.9,1.1]",
            "reject_outside_prior": True,
            "per_parameter": repro,
            "verdict": "PASS" if all_pass else "FAIL",
        },
        "signed_im_probe": {
            "reject_outside_prior": False,
            "note": ("flips drive x off-manifold; reject=False matches "
                     "validate_sbi off-manifold convention (lines 903/1012). "
                     "reject=True here legitimately stalls at 0% acceptance."),
            "signed_channel": signed_ch,
            "abs_channel": abs_ch,
            "abs_probe_desc": abs_probe_desc,
            "abs_draws_byte_identical": identical,
            "im_h2_fold": {
                "channel_convention": h2_conv,
                "draws_byte_identical": h2_identical,
                "max_shift_sigma": h2_shift,
                "fold_applied": h2_fold_ok,
                "note": ("Im_h2 is marked 'abs' and folded in training but "
                         "_x_obs_vector (sbi_runner.py:674) folds only Im_k2 "
                         "aliases. Immaterial here (Im_h2=0 at x_obs) but a "
                         "latent train/sample mismatch on the public path if a "
                         "signed Im_h2 is ever passed."),
            },
            "signed_shift_sigma": {param_names[i]: float(shift_signed[i])
                                   for i in range(len(param_names))},
            "abs_shift_sigma": {param_names[i]: float(shift_abs[i])
                                for i in range(len(param_names))},
            "max_signed_shift_sigma": max_signed,
            "max_abs_shift_sigma": max_abs,
            "verdict": "PASS" if probe_pass else "FAIL",
        },
        "overall_verdict": "PASS" if (all_pass and probe_pass) else "FAIL",
    }
    out = Path("/tmp/v2_sampling_check_report.json")
    out.write_text(json.dumps(report, indent=2))
    print(f"wrote {out}")
    print(f"\nOVERALL: {report['overall_verdict']}")
    return 0 if report["overall_verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
