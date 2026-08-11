"""Pool the three matched-resolution NH3 reference seeds (72/172/272,
n_eff=2000) into an NH3 pooled reference, using the RATIFIED salt pooling
scheme (plans/scripts/titanG_repool_reference.py): concatenate samples AND
weights, renormalizing each seed to 1/n_seeds for equal per-seed posterior
mass (pocoMC weights are importance weights, not uniform).

Why the matched seeds (not the Aug-2 single-seed n_eff=500 pkl): the
committed §0.18 P1.1 preregistration C16 ocean-fraction band and C9
nuisance-median bands (commit 817c5919) are computed from exactly these
three seeds. The Track-1 driver's C16 gate compares the corrected ocean
mass against BOTH the reference point value AND the 3-seed spread; those
two must describe the SAME reference, which forces the pooled matched-seed
reference. Mixing resolutions would compare a point value to a spread that
does not correspond to it.

Writes validation_reports/titan_freegrav_nh3_1m/reference/
titan_freegrav_nh3_reference_pooled.pkl (salt naming convention). Does NOT
overwrite the existing single-seed titan_freegrav_nh3_reference_result.pkl.
Pass the pooled path to the driver via --reference.

Run:
  mamba run -n PPcl env PYTHONPATH=. python plans/scripts/pool_nh3_matched_reference.py
"""
import pickle
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SEEDS = [72, 172, 272]
SRC = ROOT / "validation_reports/nh3_diagnosis/matched_reference"
OUT = (ROOT / "validation_reports/titan_freegrav_nh3_1m/reference"
       / "titan_freegrav_nh3_reference_pooled.pkl")


def main():
    per_seed = []
    last = None
    for sd in SEEDS:
        p = SRC / f"nh3_reference_seed{sd}_neff2000.pkl"
        o = pickle.load(open(p, "rb"))
        s = np.asarray(o.samples)
        w = np.asarray(o.weights, dtype=float)
        assert w.shape[0] == s.shape[0], \
            f"seed{sd}: weights {w.shape} != samples {s.shape}"
        per_seed.append((s, w, sd))
        last = o
    n_seeds = len(per_seed)
    pooled = np.concatenate([s for s, _, _ in per_seed], axis=0)
    pooled_w = np.concatenate(
        [w / w.sum() / n_seeds for _, w, _ in per_seed], axis=0)
    assert pooled_w.shape[0] == pooled.shape[0]
    assert abs(pooled_w.sum() - 1.0) < 1e-9, f"wsum={pooled_w.sum()}"

    # concatenate the auxiliary per-sample arrays the driver's gates read,
    # in the SAME seed order as samples/weights
    for attr in ("log_likelihoods", "D_ocean_results"):
        arrs = [np.asarray(getattr(o, attr))
                for o in [pickle.load(
                    open(SRC / f"nh3_reference_seed{sd}_neff2000.pkl", "rb"))
                    for sd in SEEDS]]
        setattr(last, attr, np.concatenate(arrs, axis=0))
    k2s = [np.asarray(pickle.load(
        open(SRC / f"nh3_reference_seed{sd}_neff2000.pkl", "rb")).k2_results)
        for sd in SEEDS]
    last.k2_results = np.concatenate(k2s, axis=0)

    last.samples = pooled
    last.weights = pooled_w
    OUT.parent.mkdir(parents=True, exist_ok=True)
    last.save(str(OUT))

    chk = pickle.load(open(OUT, "rb"))
    cs = np.asarray(chk.samples)
    cw = np.asarray(chk.weights)
    doc = np.nan_to_num(np.asarray(chk.D_ocean_results, float), nan=0.0)
    ocean_frac = float(cw[doc > 0.5].sum())
    print(f"[pool_nh3] seeds={SEEDS} pooled samples={cs.shape} "
          f"weights={cw.shape} wsum={cw.sum():.6f} "
          f"ll={np.asarray(chk.log_likelihoods).shape} "
          f"k2={np.asarray(chk.k2_results).shape} "
          f"{'OK' if cs.shape[0] == cw.shape[0] else 'MISMATCH'}")
    print(f"[pool_nh3] pooled ocean fraction (D_ocean>0.5km) = {ocean_frac:.4f} "
          f"(prereg band [0.9101,0.9222])")
    print(f"[pool_nh3] -> {OUT}")


if __name__ == "__main__":
    main()
