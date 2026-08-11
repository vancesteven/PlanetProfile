#!/usr/bin/env python
"""Track 2a — conditional mutual information I(theta; Im_k2 | C20,C22,Re_k2)
in nats, per composition. Estimator RATIFIED by scientific-reviewer
2026-08-10 (APPROVE-WITH-CHANGES, agent a872477da30620eb0);
design + ruling in plans/active/track2a-infogain-design.md.

Decision quantity (Track-1 remedy plan Track 2, sequencing item 1): the
incremental information Im_k2 carries about the interior theta beyond the
other three observables. If the CMI's combined-band UPPER bound < 0.1 nat,
CANCEL Track 2b (transform retrain) and 2c (Re/Im-only ablation) for that
composition. The gate is the NH3 CMI; NaCl/MgSO4 are free context.

Ratified estimand + reduction (EXACT for this dataset):
  x = g(theta) + eps, eps ~ N(0, diag sigma^2), eps _|_ theta, and EVERY
  rejection (k2-support box on the NOISELESS forward, no-ocean/density-
  inversion guards, drop_nonfinite) is a deterministic function of theta
  alone -> x_S|theta = N(g_S(theta), diag sigma_S^2) survives rejection, so
  H(x_S|theta) collapses to the analytic noise entropy and rejection only
  redefines the marginal to the EFFECTIVE PRIOR p'(theta) prop
  p(theta)*1[theta in A]. Hence
    CMI = H(x4) - H(x3) - 0.5*ln(2*pi*e*sigma_Im^2)
        = H(Im_k2 | C20,C22,Re_k2) - 0.5*ln(2*pi*e*sigma_Im^2)
  a <=4-dim differential-entropy problem on the kept-sample marginal; theta
  is integrated out exactly by the noise identity.

Reviewer-required estimator details (folded in):
  - Kozachenko-Leonenko k-NN differential entropy; SAME k for H(x4),H(x3)
    (the difference does NOT cancel k-NN bias since d differs); sweep
    k in {3,5,10,20} as the systematic (bias) band.
  - Rank-Gaussian transform the THREE CONDITIONING channels ONLY (invariant
    for H(Im|rest); tames Re_k2's heavy right tail). Im_k2 gets a LINEAR
    z-score only (its scale constant cancels analytically in the CMI); the
    subtracted noise entropy uses sigma_Im in the SAME z-scored units.
  - Structural check CMI >= 0 (Var(Im|rest) >= sigma_Im^2 always); a
    negative estimate flags k-NN bias.
  - CI by m-out-of-n subsampling WITHOUT replacement (20 x 200k), NOT
    with-replacement bootstrap (duplicate rows -> zero neighbor distances
    -> downward entropy bias). Combined band = variance CI (+/-) k-sweep
    systematic. CANCEL only if combined UPPER bound < 0.1 nat.
  - Gaussian moment-matched CMI reported as a BRACKET, not a bound.

Run (spare cores; reads ONLY training (theta,x) pairs, touches NO
corrected/IS-reweighted quantity -> does not interact with the Track-1
manager gate):
  mamba run -n PPcl python plans/scripts/track2a_infogain.py \
      --comp nacl mgso4 [--n-sub 200000 --n-draws 20]
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import digamma, gammaln
from scipy.stats import norm, rankdata
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[2]
DATA = {
    "nacl": "/tmp/titanG_build/datasets/titanG_nacl_1m.npz",
    "mgso4": "/tmp/titanG_build/datasets/titanG_mgso4_1m.npz",
    "nh3": "/tmp/titanG_build/datasets/titanG_nh3_1m.npz",
}
OBS_ORDER = ["C20", "C22", "Re_k2", "Im_k2"]
COND_IDX = [0, 1, 2]   # C20, C22, Re_k2
IM_IDX = 3
KS = [3, 5, 10, 20]
CANCEL_THRESHOLD_NAT = 0.10


def _ln_unit_ball_vol(d):
    """ln volume of the unit Euclidean d-ball."""
    return 0.5 * d * np.log(np.pi) - gammaln(d / 2.0 + 1.0)


def kl_entropy(X, ks):
    """Kozachenko-Leonenko differential entropy (nats) for each k in ks.

    Ĥ = psi(N) - psi(k) + ln V_d + (d/N) sum_i ln r_{i,k}, r = Euclidean
    distance to the k-th nearest neighbour (self excluded). One neighbour
    query (max(ks)+1) serves every k.
    """
    X = np.ascontiguousarray(X, dtype=np.float64)
    n, d = X.shape
    kmax = max(ks)
    nn = NearestNeighbors(n_neighbors=kmax + 1, algorithm="auto")
    nn.fit(X)
    dist, _ = nn.kneighbors(X)          # column 0 is self (distance 0)
    ln_vd = _ln_unit_ball_vol(d)
    out = {}
    for k in ks:
        r = dist[:, k]                  # k-th neighbour, self excluded
        # guard against exact ties at zero distance (duplicate rows)
        r = np.where(r > 0, r, np.finfo(float).tiny)
        out[k] = float(digamma(n) - digamma(k) + ln_vd
                       + d * np.mean(np.log(r)))
    return out


def rank_gaussian(col):
    """Monotone rank->Gaussian map (invariant reparam of a conditioner)."""
    r = rankdata(col, method="average")
    u = r / (len(col) + 1.0)
    return norm.ppf(u)


def transform_block(x):
    """Rank-Gaussian the 3 conditioners; linear z-score Im by its std.

    Returns (cond_rg [n,3], im_z [n], s_im) where s_im is the Im std used
    for the z-score (needed to place sigma_Im in the same units)."""
    cond = np.column_stack([rank_gaussian(x[:, j]) for j in COND_IDX])
    im = x[:, IM_IDX]
    s_im = float(im.std(ddof=0))
    im_z = (im - im.mean()) / s_im
    return cond, im_z, s_im


def cmi_from_block(x, sigma_im, ks):
    """CMI(nats) per k from a sample block, plus the two entropies."""
    cond_rg, im_z, s_im = transform_block(x)
    x4 = np.column_stack([cond_rg, im_z])
    x3 = cond_rg
    h4 = kl_entropy(x4, ks)
    h3 = kl_entropy(x3, ks)
    # noise entropy in the z-scored Im units: sigma_Im/s_im
    sig_z = sigma_im / s_im
    noise_ent = 0.5 * np.log(2.0 * np.pi * np.e * sig_z ** 2)
    cmi = {k: h4[k] - h3[k] - noise_ent for k in ks}
    return cmi, h4, h3, noise_ent, s_im


def gaussian_bracket(x, sigma_im):
    """Moment-matched Gaussian CMI = 0.5 ln(Var(Im|rest)/sigma_Im^2)."""
    cond = x[:, COND_IDX]
    im = x[:, IM_IDX]
    A = np.column_stack([np.ones(len(im)), cond])
    beta, *_ = np.linalg.lstsq(A, im, rcond=None)
    resid = im - A @ beta
    var_cond = float(resid.var(ddof=cond.shape[1] + 1))
    return 0.5 * np.log(var_cond / sigma_im ** 2), var_cond


def run_comp(comp, n_sub, n_draws, seed0):
    path = DATA[comp]
    if not Path(path).exists():
        return {"comp": comp, "status": "MISSING", "path": path}
    d = np.load(path, allow_pickle=True)
    x = np.asarray(d["x"], float)
    stats = json.loads(str(d["stats"]))
    sigma_im = float(stats["obs_noise"]["sigma"]["Im_k2"])
    n = x.shape[0]

    # --- point estimate on the full kept-sample marginal ---------------
    cmi_pt, h4_pt, h3_pt, noise_ent, s_im = cmi_from_block(x, sigma_im, KS)
    # raw (no rank transform; linear z-score all channels) cross-check that
    # the conditioner reparam does not move the CMI (reviewer validation 2)
    xz = (x - x.mean(0)) / x.std(0, ddof=0)
    sig_z_raw = sigma_im / float(x[:, IM_IDX].std(ddof=0))
    noise_raw = 0.5 * np.log(2.0 * np.pi * np.e * sig_z_raw ** 2)
    h4_raw = kl_entropy(xz, KS)
    h3_raw = kl_entropy(xz[:, COND_IDX], KS)
    cmi_raw = {k: h4_raw[k] - h3_raw[k] - noise_raw for k in KS}

    gcmi, var_cond = gaussian_bracket(x, sigma_im)

    # --- m-out-of-n subsampling WITHOUT replacement for the CI ---------
    rng = np.random.default_rng(seed0)
    m = min(n_sub, n)
    sub_cmi = {k: [] for k in KS}
    for _ in range(n_draws):
        idx = rng.choice(n, size=m, replace=False)
        c, _, _, _, _ = cmi_from_block(x[idx], sigma_im, KS)
        for k in KS:
            sub_cmi[k].append(c[k])

    # per-k variance CI (2-sigma) from subsamples; combined band folds the
    # k-sweep systematic spread of the point estimate
    per_k = {}
    for k in KS:
        arr = np.asarray(sub_cmi[k], float)
        per_k[str(k)] = {
            "point_full": cmi_pt[k],
            "point_full_raw": cmi_raw[k],
            "H_x4": h4_pt[k], "H_x3": h3_pt[k],
            "sub_mean": float(arr.mean()),
            "sub_std": float(arr.std(ddof=1)),
            "sub_ci95": [float(arr.mean() - 2 * arr.std(ddof=1)),
                         float(arr.mean() + 2 * arr.std(ddof=1))],
        }
    pts = np.asarray([cmi_pt[k] for k in KS], float)
    ksweep_lo, ksweep_hi = float(pts.min()), float(pts.max())
    # combined upper bound: widest per-k upper CI, widened by the k-sweep
    # systematic half-range on top of the k-sweep max
    max_sub_std = max(per_k[str(k)]["sub_std"] for k in KS)
    ksweep_half = 0.5 * (ksweep_hi - ksweep_lo)
    combined_upper = ksweep_hi + 2 * max_sub_std + ksweep_half
    combined_lower = ksweep_lo - 2 * max_sub_std - ksweep_half

    cmi_min = float(pts.min())
    return {
        "comp": comp,
        "status": "OK",
        "n_kept": int(n),
        "sigma_Im": sigma_im,
        "im_std": s_im,
        "noise_entropy_zunits": float(noise_ent),
        "ks": KS,
        "cmi_point_per_k": {str(k): cmi_pt[k] for k in KS},
        "cmi_point_raw_per_k": {str(k): cmi_raw[k] for k in KS},
        "per_k": per_k,
        "k_sweep_range_nat": [ksweep_lo, ksweep_hi],
        "gaussian_bracket_nat": float(gcmi),
        "var_im_given_rest": float(var_cond),
        "structural_cmi_nonneg": bool(cmi_min >= -1e-6),
        "cmi_min_over_k_nat": cmi_min,
        "combined_band_nat": [combined_lower, combined_upper],
        "cancel_threshold_nat": CANCEL_THRESHOLD_NAT,
        "cancel_2b2c": bool(combined_upper < CANCEL_THRESHOLD_NAT),
        "subsample": {"m": int(m), "n_draws": int(n_draws)},
        "note": ("CMI = I(theta; Im_k2 | C20,C22,Re_k2), effective-prior "
                 "(support-guarded) marginal. Reads training (theta,x) "
                 "only; no corrected/IS quantity. CANCEL 2b/2c iff combined "
                 "band upper < 0.1 nat. Gaussian value is a BRACKET, not a "
                 "bound. Estimand is an UPPER bound on tidal-parameter-"
                 "specific recoverable info (>0.1 nat is necessary-not-"
                 "sufficient for 2b/2c to help)."),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--comp", nargs="+", default=["nacl", "mgso4"],
                    choices=sorted(DATA))
    ap.add_argument("--n-sub", type=int, default=200000)
    ap.add_argument("--n-draws", type=int, default=20)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    results = {}
    for i, comp in enumerate(args.comp):
        print(f"[track2a] {comp}: estimating CMI ...", flush=True)
        r = run_comp(comp, args.n_sub, args.n_draws, args.seed + 1000 * i)
        results[comp] = r
        if r["status"] != "OK":
            print(f"[track2a] {comp}: {r['status']} ({r.get('path')})")
            continue
        band = r["combined_band_nat"]
        print(f"[track2a] {comp}: CMI k-sweep {r['k_sweep_range_nat']} nat "
              f"| combined band [{band[0]:.3f},{band[1]:.3f}] "
              f"| gaussian bracket {r['gaussian_bracket_nat']:.3f} "
              f"| CMI>=0 {r['structural_cmi_nonneg']} "
              f"| CANCEL 2b/2c {r['cancel_2b2c']}")

    out_dir = ROOT / "validation_reports/track2a_infogain"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "track2a_cmi.json"
    payload = {
        "kind": "track2a_conditional_mutual_information",
        "estimand": "I(theta; Im_k2 | C20,C22,Re_k2) [nats], effective-prior",
        "ratified_by": "scientific-reviewer 2026-08-10 (a872477da30620eb0)",
        "cancel_rule": "CANCEL 2b/2c iff combined band upper < 0.1 nat; "
                       "gate composition is NH3 (salts are context)",
        "results": results,
    }
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"[track2a] -> {out_path}")


if __name__ == "__main__":
    main()
