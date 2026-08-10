"""Repair the Titan free-gravity MgSO4/NaCl POOLED reference pkl (no MCMC recompute).

Root cause (2026-08-10): titanG_ocean_reference_mcmc.py pooled the multi-seed
reference by setting ``res.samples = concat(all seed samples)`` but never updating
``res.weights`` — so the pooled pkl carried only the LAST seed's weight array
(len ~6365) against all-seed samples (len ~19325). validate_sbi crosscheck then
raised ``weights length != n_samples``. The per-seed pkls are intact (weights
match samples, each sums to 1); pocoMC weights are importance weights (NOT
uniform), so pooling must concatenate BOTH samples and weights, renormalizing
each seed to 1/n_seeds for equal per-seed posterior mass.

This rebuilds the pooled pkl from the committed per-seed pkls using the corrected
scheme (mirrors the fixed titanG_ocean_reference_mcmc.py pooling block). No Test/
write; reads/writes only validation_reports/.../reference/.

Run:
  mamba run -n PPcl env PYTHONPATH=. python plans/scripts/titanG_repool_reference.py
"""
import pickle
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
COMPS = {"mgso4": [73, 173, 273], "nacl": [74, 174, 274]}


def repool(tag, seeds):
    d = ROOT / f"validation_reports/titan_freegrav_{tag}_1m/reference"
    per_seed = []
    last = None
    for sd in seeds:
        o = pickle.load(open(d / f"titan_freegrav_{tag}_reference_seed{sd}.pkl", "rb"))
        s = np.asarray(o.samples)
        w = np.asarray(o.weights, dtype=float)
        assert w.shape[0] == s.shape[0], f"{tag} seed{sd}: weights {w.shape} != samples {s.shape}"
        per_seed.append((s, w))
        last = o
    n_seeds = len(per_seed)
    pooled = np.concatenate([s for s, _ in per_seed], axis=0)
    pooled_w = np.concatenate([w / w.sum() / n_seeds for _, w in per_seed], axis=0)
    assert pooled_w.shape[0] == pooled.shape[0]
    assert abs(pooled_w.sum() - 1.0) < 1e-9, f"{tag} pooled weights sum={pooled_w.sum()}"
    last.samples = pooled
    last.weights = pooled_w
    out = d / f"titan_freegrav_{tag}_reference_pooled.pkl"
    last.save(str(out))
    # verify round-trip
    chk = pickle.load(open(out, "rb"))
    cs = np.asarray(chk.samples); cw = np.asarray(chk.weights)
    print(f"[repool] {tag}: n_seeds={n_seeds} pooled samples={cs.shape} weights={cw.shape} "
          f"wsum={cw.sum():.6f}  {'OK' if cs.shape[0]==cw.shape[0] else 'MISMATCH'}")


if __name__ == "__main__":
    for tag, seeds in COMPS.items():
        repool(tag, seeds)
