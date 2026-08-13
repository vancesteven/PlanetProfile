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
    """Pool ALL per-sample arrays (samples, weights, AND every derived array)
    in seed order.

    2026-08-12 completion: the 2026-08-10 repair concatenated only ``samples``
    and ``weights`` and left every derived per-sample array (``log_likelihoods``,
    ``k2_results``, ``D_ocean_results``, ``c20/c22_results``, ``h2_results``,
    ``cmr2_results``, the ``D_*`` thickness arrays, ``cmr2_hydro_results``) at
    the LAST seed's length. Any consumer indexing a derived array against the
    pooled weights (the IS-validation C3/pushforward/C16 path) then crashed with
    an IndexError. This pools every ndarray attribute whose leading dim equals
    that seed's sample count, in the SAME seed order as samples, so derived[i]
    stays aligned to sample[i]. Only ``weights`` gets the 1/n_seeds renorm; the
    derived arrays are plain-concatenated (they are per-sample quantities, not
    posterior mass). No MCMC recompute; reads/writes only validation_reports/.
    """
    d = ROOT / f"validation_reports/titan_freegrav_{tag}_1m/reference"
    objs = []
    for sd in seeds:
        o = pickle.load(open(d / f"titan_freegrav_{tag}_reference_seed{sd}.pkl", "rb"))
        s = np.asarray(o.samples)
        w = np.asarray(o.weights, dtype=float)
        assert w.shape[0] == s.shape[0], f"{tag} seed{sd}: weights {w.shape} != samples {s.shape}"
        objs.append(o)
    n_seeds = len(objs)
    last = objs[-1]
    n_per = [np.asarray(o.samples).shape[0] for o in objs]

    # Discover every per-sample array attribute (leading dim == that seed's
    # n_samples on ALL seeds) other than weights, which is handled specially.
    per_sample_attrs = []
    for k, v in vars(last).items():
        if k == "weights":
            continue
        if not isinstance(v, np.ndarray) or v.ndim < 1:
            continue
        if all(isinstance(getattr(o, k, None), np.ndarray)
               and np.asarray(getattr(o, k)).shape[0] == n
               for o, n in zip(objs, n_per)):
            per_sample_attrs.append(k)

    pooled_w = np.concatenate(
        [np.asarray(o.weights, float) / np.asarray(o.weights, float).sum() / n_seeds
         for o in objs], axis=0)
    n_total = pooled_w.shape[0]
    assert n_total == sum(n_per)
    assert abs(pooled_w.sum() - 1.0) < 1e-9, f"{tag} pooled weights sum={pooled_w.sum()}"

    for k in per_sample_attrs:
        pooled_arr = np.concatenate([np.asarray(getattr(o, k)) for o in objs], axis=0)
        assert pooled_arr.shape[0] == n_total, f"{tag} {k}: {pooled_arr.shape[0]} != {n_total}"
        setattr(last, k, pooled_arr)
    last.weights = pooled_w

    out = d / f"titan_freegrav_{tag}_reference_pooled.pkl"
    last.save(str(out))
    # verify round-trip: every per-sample array must match the pooled length.
    chk = pickle.load(open(out, "rb"))
    lens = {"weights": np.asarray(chk.weights).shape[0]}
    for k in per_sample_attrs:
        lens[k] = np.asarray(getattr(chk, k)).shape[0]
    ok = all(v == n_total for v in lens.values())
    print(f"[repool] {tag}: n_seeds={n_seeds} n_total={n_total} wsum={np.asarray(chk.weights).sum():.6f} "
          f"pooled_attrs={sorted(per_sample_attrs)}  {'OK' if ok else 'MISMATCH ' + str(lens)}")


if __name__ == "__main__":
    for tag, seeds in COMPS.items():
        repool(tag, seeds)
