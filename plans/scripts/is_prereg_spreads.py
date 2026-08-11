#!/usr/bin/env python
"""C9/C16 PREREGISTRATION spreads for the Track-1 IS correction
(plans/active/tidal-sector-remedy-plan.md, conditions C9 + C16).

Reviewer rule (§0.17 task 2 / §0.18 Phase-1 item 1): BEFORE reading ANY
corrected result, compute and commit each reference MCMC's 3-seed
- ocean-fraction spread (C16 acceptance band: the corrected has_ocean mass
  must match the reference within THIS spread), and
- per-parameter weighted-median spread across seeds (C9 unidentified-
  nuisance adjudication band).

These are computed from the REFERENCE MCMC seed pkls ONLY — no corrected /
IS-reweighted quantity is read here, by construction; that is what makes
this a preregistration and not a post-hoc tolerance.

Ocean fraction per seed = sum of importance weights where D_ocean > 0.5 km
(same >0.5 km liquid-layer threshold the is_correction C16 branch uses).
Between-seed spread reported as [min, max], mean, sample-std, and half-range
(the band half-width the corrected value must fall within).

Compositions:
  nh3   : matched-resolution n_eff=2000 seeds 72/172/272 (this campaign's
          B3-class reference; the single-seed n_eff=500 production pkl is the
          ratification artifact but NOT the 3-seed spread source).
  mgso4 : B3 seeds 73/173/273.
  nacl  : B3 seeds 74/174/274.

Run:
  mamba run -n PPcl python plans/scripts/is_prereg_spreads.py
"""
import json
import pickle
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
D_OCEAN_LIQUID_KM = 0.5  # matches is_correction.compute_is_correction C16

# Reference seed pkls per composition. NH3 uses the matched-resolution
# n_eff=2000 diagnosis references (same joint 2D cache, sha 3d837c...) so it
# is B3-class like the salts; the salts use their committed B3 seed pkls.
REFS = {
    "nh3": [
        ROOT / "validation_reports/nh3_diagnosis/matched_reference"
             / "nh3_reference_seed72_neff2000.pkl",
        ROOT / "validation_reports/nh3_diagnosis/matched_reference"
             / "nh3_reference_seed172_neff2000.pkl",
        ROOT / "validation_reports/nh3_diagnosis/matched_reference"
             / "nh3_reference_seed272_neff2000.pkl",
    ],
    "mgso4": [
        ROOT / f"validation_reports/titan_freegrav_mgso4_1m/reference"
             / f"titan_freegrav_mgso4_reference_seed{s}.pkl"
        for s in (73, 173, 273)
    ],
    "nacl": [
        ROOT / f"validation_reports/titan_freegrav_nacl_1m/reference"
             / f"titan_freegrav_nacl_reference_seed{s}.pkl"
        for s in (74, 174, 274)
    ],
}

# C9 "unidentified nuisances": the dissipation/rheology block is the
# poorly-identified sector (NaCl crosscheck FAILed on log10_eta_V; §0.16).
# We report the spread for ALL 13 params but flag these as the C9 set.
C9_NUISANCE_PARAMS = [
    "log10_eta_Ih", "log10_eta_III", "log10_eta_V", "log10_eta_VI",
    "log10_eta_sil", "log10_zeta", "alpha",
]


def _wquantile(v, w, q):
    """Weighted quantile (linear interp on the weighted CDF)."""
    v = np.asarray(v, float)
    w = np.asarray(w, float)
    ok = np.isfinite(v) & np.isfinite(w) & (w > 0)
    v, w = v[ok], w[ok]
    order = np.argsort(v)
    v, w = v[order], w[order]
    cw = np.cumsum(w)
    cw = cw / cw[-1]
    # midpoint convention
    p = (cw - 0.5 * w / w.sum())
    return float(np.interp(q, p, v))


def _load(p):
    with open(p, "rb") as f:
        return pickle.load(f)


def _spread(vals):
    a = np.asarray(vals, float)
    return {
        "per_seed": [float(x) for x in a],
        "min": float(a.min()),
        "max": float(a.max()),
        "mean": float(a.mean()),
        "std": float(a.std(ddof=1)) if a.size > 1 else 0.0,
        "half_range": float(0.5 * (a.max() - a.min())),
    }


def main():
    out = {
        "kind": "is_correction_preregistration_spreads",
        "note": ("C9/C16 preregistration. Computed from REFERENCE MCMC seed "
                 "pkls ONLY (no corrected/IS-reweighted quantity read). The "
                 "C16 ocean-fraction acceptance band and the C9 nuisance "
                 "median band are FROZEN by this commit BEFORE any corrected "
                 "result is read. D_ocean liquid threshold "
                 f"{D_OCEAN_LIQUID_KM} km."),
        "d_ocean_liquid_km": D_OCEAN_LIQUID_KM,
        "c9_nuisance_params": C9_NUISANCE_PARAMS,
        "compositions": {},
    }
    for comp, paths in REFS.items():
        missing = [str(p) for p in paths if not p.exists()]
        if missing:
            out["compositions"][comp] = {"status": "INCOMPLETE",
                                         "missing": missing}
            print(f"[prereg] {comp}: MISSING {missing}")
            continue
        objs = [_load(p) for p in paths]
        pnames = list(objs[0].param_names)
        seeds_ocean, seeds_medians = [], {p: [] for p in pnames}
        seed_ids = []
        for o, path in zip(objs, paths):
            w = np.asarray(o.weights, float)
            w = w / w.sum()
            doc = np.nan_to_num(np.asarray(o.D_ocean_results, float), nan=0.0)
            seeds_ocean.append(float(w[doc > D_OCEAN_LIQUID_KM].sum()))
            S = np.asarray(o.samples, float)
            for j, p in enumerate(pnames):
                seeds_medians[p].append(_wquantile(S[:, j], w, 0.5))
            # recover seed id from filename
            seed_ids.append(path.stem)
        entry = {
            "status": "OK",
            "seed_pkls": [p.name for p in paths],
            "n_samples_per_seed": [int(o.samples.shape[0]) for o in objs],
            "ocean_fraction": _spread(seeds_ocean),
            "param_median_spread": {
                p: _spread(seeds_medians[p]) for p in pnames},
        }
        out["compositions"][comp] = entry
        of = entry["ocean_fraction"]
        print(f"[prereg] {comp}: ocean_frac per-seed {of['per_seed']} "
              f"=> band [{of['min']:.4f}, {of['max']:.4f}] "
              f"half-range {of['half_range']:.4f}")

    out_dir = ROOT / "validation_reports/is_correction_prereg"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "is_prereg_spreads.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"[prereg] -> {out_path}")


if __name__ == "__main__":
    main()
