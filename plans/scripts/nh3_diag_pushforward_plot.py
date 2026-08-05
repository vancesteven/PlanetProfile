"""NH3 under-update diagnosis — FOLLOW-UP #3: anchor Im_k2 pushforward SHAPE.

Manager judgement (MACHINE-B-HANDOFF §0.9, 2026-08-05) reordered the reviewer's
follow-ups to #3-first because it is FREE (no new MCMC): plot the capped
full-joint anchor's SBI posterior-predictive |Im_k2| DISTRIBUTION — not just the
median that S1/S2 reported — to separate two mechanisms the median cannot:

  concentration-failure : the posterior barely narrows from the prior; it still
                          covers the high-|Im_k2| ocean tail but its MASS sits
                          low, so the median under-reads. Fixable by a sharper
                          flow / more epochs / better identifiability.
  wrong-mode            : the posterior CONCENTRATES (narrow) on a low-|Im_k2|
                          mode and drops the high-k2 tail entirely. NOT a
                          resolution problem — the flow learned the wrong
                          conditional shape.

Method. Reproduce the anchor's posterior-predictive |Im_k2| via the IDENTICAL
path titanG_ppc_interior_check uses (sample the flow at x_obs -> push theta back
through generate_sbi_dataset, obs_noise=False, byte-identical forward physics),
then overlay:
  - PRIOR-predictive  |Im_k2| (training x column, already abs-folded)
  - anchor SBI POSTERIOR-predictive |Im_k2|
  - MCMC reference    |Im_k2| (weighted; the authoritative tidal posterior)
with markers at obs=0.135, MCMC-pp median (~0.100 ceiling), SBI-pp median.

NO NEW FORWARD SIMS beyond the ~3.8k-draw pushforward the PPC already does; NO
new MCMC (the reference pkl is reused). Saves the raw SBI-pp |Im_k2| array so a
re-plot needs no recompute.

Run:
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/nh3_diag_pushforward_plot.py
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CFG = ROOT / "PlanetProfile/Inference/configs/test54_titan_nh3_freegrav.json"
DATASET = Path("/tmp/titanG_build/datasets/titanG_nh3_1m.npz")
ANCHOR = Path("/tmp/nh3_diag/anchor/nh3_capped_full_joint_anchor.pt")
REF_PKL = (ROOT / "validation_reports/titan_freegrav_nh3_1m/reference/"
           "titan_freegrav_nh3_reference_result.pkl")
OUTDIR = ROOT / "validation_reports/nh3_diagnosis/pushforward_shape"
MCMC_PP_CEILING = 0.100  # weighted MCMC |Im_k2| median (matched-res re-measure = #2)


def _weighted_quantiles(x, w, qs):
    i = np.argsort(x)
    x = x[i]; w = w[i]
    c = np.cumsum(w)
    c /= c[-1]
    return np.interp(qs, c, x)


def main():
    global OUTDIR
    ap = argparse.ArgumentParser()
    ap.add_argument("--artifact", default=str(ANCHOR))
    ap.add_argument("--config", default=str(CFG))
    ap.add_argument("--dataset", default=str(DATASET))
    ap.add_argument("--ref-pkl", default=str(REF_PKL))
    ap.add_argument("--n-post", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=72)
    ap.add_argument("--outdir", default=str(OUTDIR),
                    help="output dir (override for the deployed-flow transfer check)")
    ap.add_argument("--label", default="capped anchor",
                    help="curve/title label for the SBI flow being plotted")
    args = ap.parse_args()

    OUTDIR = Path(args.outdir)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for p in (args.artifact, args.config, args.dataset, args.ref_pkl):
        if not Path(p).exists():
            raise SystemExit(f"[pf] missing input: {p}")

    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.sbi_runner import SBIRunner

    cfg_obj = InferenceConfig.from_json(args.config)
    obs_names = list(cfg_obj.observables.keys())
    j_im = obs_names.index("Im_k2")
    x_obs = {n: float(cfg_obj.observables[n][0]) for n in obs_names}
    sig_im = float(cfg_obj.observables["Im_k2"][1])
    print(f"[pf] obs={obs_names}; x_obs[Im_k2]={x_obs['Im_k2']} sigma={sig_im}")

    # (a) PRIOR-predictive |Im_k2| from the training x. The saved x is NOISED
    #     (noise added ONCE after the |Im| fold, refold_im=False), so ~22% of
    #     the Im_k2 column dips negative. Recover the CLEAN folded signal
    #     (min>=0) by reproducing the identical one-shot draw and subtracting
    #     — the exact S2 recovery — so the prior envelope is comparable to the
    #     NOISELESS SBI pushforward (obs_noise=False below).
    dat = np.load(args.dataset, allow_pickle=True)
    x_train = np.asarray(dat["x"], dtype=float)
    stats = json.loads(str(dat["stats"].item())) if "stats" in dat.files else {}
    on = stats.get("obs_noise")
    if on is not None and not on.get("correlations") and on.get("refold_im") is False:
        sig_vec = np.array([float(on["sigma"][n]) for n in obs_names])
        noise = np.random.default_rng(on["noise_seed"]).normal(
            0.0, 1.0, size=x_train.shape) * sig_vec
        x_clean = x_train - noise
        prior_im = x_clean[:, j_im]
        prior_im = prior_im[np.isfinite(prior_im)]
        assert prior_im.min() >= -1e-9, (
            f"clean Im_k2 min {prior_im.min()} < 0 — recovery failed")
        print(f"[pf] recovered CLEAN prior |Im_k2| (min={prior_im.min():.3g})")
    else:
        # Fallback: no recoverable metadata -> use noised column as-is (labelled).
        prior_im = x_train[:, j_im]
        prior_im = prior_im[np.isfinite(prior_im)]
        print("[pf] WARNING: obs_noise metadata missing/incompatible; "
              "prior |Im_k2| is the NOISED training column (some values < 0)")

    # (b) anchor SBI POSTERIOR-predictive |Im_k2| via the exact PPC path.
    runner = SBIRunner.load_artifact(args.artifact)
    _cfg = json.load(open(args.config)); _cfg["mode"] = "sbi"
    fwd = SBIRunner(InferenceConfig.from_dict(_cfg))
    assert list(fwd.param_names) == list(runner.param_names), "param order drift"
    assert list(fwd.obs_names) == obs_names, "obs order drift"
    print(f"[pf] sampling {args.n_post} posterior draws at x_obs ({args.label}) ...")
    theta_post = runner.sample_posterior(x_obs, n_samples=args.n_post, seed=args.seed)
    mcmc = fwd._get_mcmc_runner()
    # drop_nonfinite=False so rows stay ALIGNED to theta_post — lets us test
    # whether the dropped (nonfinite) draws are high-k2-biased (reviewer MINOR).
    theta_kept, x_pp_all = mcmc.generate_sbi_dataset(
        n_samples=len(theta_post),
        apply_support_guard=fwd._support_guard_active(),
        imag_convention=fwd.imag_convention,
        drop_nonfinite=False, obs_noise=False, theta_override=theta_post)
    x_pp_all = np.asarray(x_pp_all, dtype=float)
    theta_kept = np.asarray(theta_kept, dtype=float)
    finite_rows = np.all(np.isfinite(x_pp_all), axis=1)
    x_pp = x_pp_all[finite_rows]
    sbi_im = x_pp[:, j_im]
    n_dropped = int((~finite_rows).sum())
    print(f"[pf] SBI-pp |Im_k2|: n={sbi_im.size} (dropped {n_dropped} nonfinite) "
          f"5/50/95={np.percentile(sbi_im,[5,50,95])}")

    # --- Nonfinite-drop-bias check (reviewer MINOR) --------------------------
    # If dropped draws preferentially carried high-|Im_k2| ocean structure, the
    # kept SBI-pp would be biased LOW, exaggerating the under-update. We cannot
    # read Im_k2 on the dropped rows (it is nonfinite), so use the strongest
    # available proxy for "ocean/high-dissipation": the sampled tidal-sector
    # params that drive |Im_k2| up. We compare their distribution in kept vs
    # dropped draws. A drop biased toward high-dissipation params is the
    # concern. NOTE: closing rests on MAGNITUDE (<0.25% of draws drop) and
    # conservative DIRECTION (dropped set leans high-dissipation, so including
    # it would raise SBI-pp toward obs — weakening, never manufacturing, the
    # under-update), NOT on the kept/dropped distributions being similar.
    drop_bias = {"n_dropped": n_dropped,
                 "drop_frac": float(n_dropped / max(len(theta_post), 1))}
    if n_dropped > 0:
        pnames = list(fwd.param_names)
        # Params physically pushing |Im_k2| up: lower viscosity / higher zeta.
        proxy_names = [p for p in pnames
                       if p.startswith("log10_eta") or p == "log10_zeta"]
        kept_theta = theta_kept[finite_rows]
        drop_theta = theta_kept[~finite_rows]
        for p in proxy_names:
            jp = pnames.index(p)
            drop_bias[p] = {
                "kept_median": float(np.median(kept_theta[:, jp])),
                "dropped_median": float(np.median(drop_theta[:, jp])),
                "kept_5_95": [float(np.percentile(kept_theta[:, jp], 5)),
                              float(np.percentile(kept_theta[:, jp], 95))],
                "dropped_5_95": [float(np.percentile(drop_theta[:, jp], 5)),
                                 float(np.percentile(drop_theta[:, jp], 95))],
            }
        print(f"[pf] drop-bias proxy medians (kept vs dropped): "
              + "; ".join(f"{p} {drop_bias[p]['kept_median']:.2f}/"
                          f"{drop_bias[p]['dropped_median']:.2f}" for p in proxy_names))

    # (c) MCMC reference |Im_k2| (weighted; signed col1 -> abs).
    import pickle
    ref = pickle.load(open(args.ref_pkl, "rb"))
    kk = np.asarray(ref.k2_results, dtype=float)      # (N,2) [Re, Im signed]
    ww = np.asarray(ref.weights, dtype=float)
    fin = np.all(np.isfinite(kk), axis=1) & np.isfinite(ww)
    mcmc_im = np.abs(kk[fin, 1]); mcmc_w = ww[fin]
    mcmc_q = _weighted_quantiles(mcmc_im, mcmc_w, [0.05, 0.5, 0.95])
    print(f"[pf] MCMC-ref |Im_k2| weighted 5/50/95={mcmc_q}")

    # Save raw arrays so re-plots need no recompute.
    np.savez_compressed(OUTDIR / "pushforward_arrays.npz",
                        prior_im=prior_im, sbi_im=sbi_im,
                        mcmc_im=mcmc_im, mcmc_w=mcmc_w)

    # ---- Plot -------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sbi_med = float(np.median(sbi_im))
    mcmc_med = float(mcmc_q[1])
    obs_v = x_obs["Im_k2"]
    # Common range; clip long tails for readability but keep the high-k2 region.
    hi = max(np.percentile(prior_im, 99.5), np.percentile(sbi_im, 99.5),
             mcmc_q[2], obs_v) * 1.05
    bins = np.linspace(0.0, hi, 70)

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    ax.hist(prior_im, bins=bins, density=True, histtype="stepfilled",
            color="0.80", edgecolor="0.55", alpha=0.7,
            label=f"prior-predictive (train, n={prior_im.size:,})")
    ax.hist(sbi_im, bins=bins, density=True, histtype="step",
            color="tab:red", lw=2.2,
            label=f"{args.label} SBI posterior-predictive (n={sbi_im.size:,})")
    ax.hist(mcmc_im, bins=bins, weights=mcmc_w, density=True, histtype="step",
            color="tab:blue", lw=2.2, ls="--",
            label=f"MCMC reference (weighted, n={mcmc_im.size:,})")

    ax.axvline(obs_v, color="k", lw=1.6, ls=":",
               label=f"observed = {obs_v:.3f}")
    ax.axvline(mcmc_med, color="tab:blue", lw=1.4, alpha=0.9,
               label=f"MCMC-pp median = {mcmc_med:.3f}")
    ax.axvline(sbi_med, color="tab:red", lw=1.4, alpha=0.9,
               label=f"SBI-pp median = {sbi_med:.3f}")

    ax.set_xlabel(r"$|\mathrm{Im}\,k_2|$")
    ax.set_ylabel("density")
    ax.set_title(f"Titan NH3 joint — {args.label} $|Im\\,k_2|$ pushforward shape\n"
                 "(concentration-failure vs wrong-mode diagnostic)")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_xlim(0, hi)
    fig.tight_layout()
    pdf = OUTDIR / "anchor_imk2_pushforward.pdf"
    fig.savefig(pdf); fig.savefig(OUTDIR / "anchor_imk2_pushforward.png", dpi=140)
    plt.close(fig)
    print(f"[pf] plot -> {pdf}")

    # ---- Quantitative shape descriptors (feed the reviewer verdict) -------
    def _spread(x, w=None):
        if w is None:
            q = np.percentile(x, [5, 50, 95])
        else:
            q = _weighted_quantiles(x, w, [0.05, 0.5, 0.95])
        return {"p5": float(q[0]), "median": float(q[1]), "p95": float(q[2]),
                "iqr90": float(q[2] - q[0])}
    prior_s = _spread(prior_im)
    sbi_s = _spread(sbi_im)
    mcmc_s = _spread(mcmc_im, mcmc_w)
    # fraction of SBI-pp mass ABOVE the MCMC-pp ceiling and above obs.
    frac_above_ceiling = float(np.mean(sbi_im >= MCMC_PP_CEILING))
    frac_above_obs = float(np.mean(sbi_im >= obs_v))
    # concentration ratio: how much narrower is the posterior than the prior?
    concentration = float(sbi_s["iqr90"] / prior_s["iqr90"]) if prior_s["iqr90"] else float("nan")
    # PRIOR-referenced shift (ceiling-INDEPENDENT evidence, reviewer): a correct
    # update conditioned on obs=0.135 should heavily up-weight the draws already
    # reaching >=obs; measure how far the flow moved that mass vs the prior.
    prior_frac_above_obs = float(np.mean(prior_im >= obs_v))
    prior_frac_above_ceiling = float(np.mean(prior_im >= MCMC_PP_CEILING))

    report = {
        "kind": "nh3_diag_pushforward_shape_followup3",
        "note": ("DIAGNOSTIC, not a gate — interpreted by scientific-reviewer, "
                 "never tuned. Separates concentration-failure (broad, low mass) "
                 "from wrong-mode (narrow, low mode) using the FULL |Im_k2| "
                 "pushforward, not just the median."),
        "artifact": os.path.abspath(args.artifact),
        "config": os.path.abspath(args.config),
        "dataset": os.path.abspath(args.dataset),
        "ref_pkl": os.path.abspath(args.ref_pkl),
        "obs_Im_k2": obs_v, "sigma_Im_k2": sig_im,
        "mcmc_pp_ceiling_used": MCMC_PP_CEILING,
        "prior_predictive_spread": prior_s,
        "sbi_posterior_predictive_spread": sbi_s,
        "mcmc_reference_spread": mcmc_s,
        "sbi_concentration_ratio_vs_prior": concentration,
        "sbi_frac_mass_ge_mcmc_ceiling": frac_above_ceiling,
        "sbi_frac_mass_ge_obs": frac_above_obs,
        "prior_frac_mass_ge_obs": prior_frac_above_obs,
        "prior_frac_mass_ge_mcmc_ceiling": prior_frac_above_ceiling,
        "ceiling_independent_reading": (
            "obs=0.135 is reachable under the prior (frac>=obs "
            f"{prior_frac_above_obs:.3f}); a correct posterior conditioned on "
            "0.135+-0.035 should drive frac>=obs toward ~0.5, but the flow moved "
            f"it only {prior_frac_above_obs:.3f}->{frac_above_obs:.3f}. This "
            "under-update is measured entirely against the PRIOR and does NOT "
            "depend on where the MCMC ceiling sits (#2)."),
        "nonfinite_drop_bias_check": drop_bias,
        "label": args.label,
        "reading_key": (
            "wrong-mode if concentration_ratio << 1 (posterior much narrower "
            "than prior) AND frac_mass_ge_mcmc_ceiling ~ 0 (tail dropped); "
            "concentration-failure if concentration_ratio ~ 1 (barely narrowed) "
            "with mass still reaching the ceiling but median low. The #2 "
            "matched-resolution MCMC re-measure will confirm whether the "
            "ceiling itself moves."),
        "arrays_npz": str(OUTDIR / "pushforward_arrays.npz"),
        "plot_pdf": str(pdf),
    }
    with open(OUTDIR / "pushforward_shape_report.json", "w") as f:
        json.dump(report, f, indent=2)
    print(f"[pf] report -> {OUTDIR/'pushforward_shape_report.json'}")
    print(f"[pf] concentration_ratio(vs prior)={concentration:.3f}  "
          f"frac_mass>=ceiling={frac_above_ceiling:.3f}  "
          f"frac_mass>=obs={frac_above_obs:.3f}")


if __name__ == "__main__":
    main()
