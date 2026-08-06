"""Pool the three FRESH B3 v5 reference-MCMC seeds into one n_eff~2000 reference
for the v5 crosscheck gate (Machine A §0.7 step-2 / recon Step A).

WHY: the deployed v5 crosscheck (v5_run_gates.py REF_MCMC) still targets the OLD
n_eff=500 europa_clipper_v5_reference_result.pkl. §0.7 step-2 forbids that stale
target and requires the FRESH pooled n_eff=2000 references produced by B3. The
1.06 km v5-v7 wander was an n_eff=500 resolution artifact; at matched n_eff=2000
the gap collapsed to -0.19+-0.22 km (B3). This pools seeds 101/202/303.

POOLING (concatenate, do NOT average): concatenate each seed's `samples`,
concatenate `weights` and renormalize to sum=1 (each seed already sums to 1, so
the pool weights all three seeds equally), concatenate `log_likelihoods`, reuse
seed-101's config/param_names/param_labels, construct one InferenceResult.
Concatenation preserves each seed's full weighted marginal; equal-weight-per-seed
follows from the renormalization.

Output: /tmp/b3_build/europa_clipper_v5_reference_pooled_neff2000.pkl
(a build artifact in /tmp per repo convention; NOT written into Test/).

Run:
  mamba run -n PPcl env PYTHONPATH=. KMP_DUPLICATE_LIB_OK=TRUE \
    python plans/scripts/pool_v5_reference_neff2000.py
"""
import numpy as np
from pathlib import Path

from PlanetProfile.Inference.inference_core import InferenceResult

SRC = Path("/tmp/b3_build/v5_reference_wander")
SEEDS = [101, 202, 303]
OUT = Path("/tmp/b3_build/europa_clipper_v5_reference_pooled_neff2000.pkl")


def _ess(w):
    w = np.asarray(w, float)
    return float(w.sum() ** 2 / np.sum(w ** 2))


def main():
    pkls = [SRC / f"reference_seed{s}.pkl" for s in SEEDS]
    for p in pkls:
        if not p.exists():
            raise SystemExit(f"[pool] missing fresh B3 reference {p}")

    results = [InferenceResult.load(str(p)) for p in pkls]
    base = results[0]
    pnames = list(base.param_names)
    for r, s in zip(results, SEEDS):
        assert list(r.param_names) == pnames, f"param_names mismatch seed {s}"
        print(f"[pool] seed {s}: samples={r.samples.shape} "
              f"ESS={_ess(r.weights):.0f} w_sum={np.sum(r.weights):.6f}")

    samples = np.concatenate([r.samples for r in results], axis=0)
    # renormalize each seed's weights to 1 (defensive) then pool with equal
    # per-seed mass, then renormalize the pool to sum=1.
    w_parts = [np.asarray(r.weights, float) / np.sum(r.weights) for r in results]
    weights = np.concatenate(w_parts, axis=0)
    weights = weights / weights.sum()
    loglik = np.concatenate([np.asarray(r.log_likelihoods, float)
                             for r in results], axis=0)

    pooled = InferenceResult(
        samples=samples, weights=weights, log_likelihoods=loglik,
        param_names=pnames, param_labels=list(base.param_labels),
        config=base.config,
    )
    print(f"[pool] pooled: samples={pooled.samples.shape} "
          f"ESS={_ess(pooled.weights):.0f} config_hash={pooled.config.generate_hash()}")

    # sanity: pooled D_iceIh weighted mean/std vs per-seed
    if "D_iceIh_km" in pnames:
        j = pnames.index("D_iceIh_km")
        wm = np.average(samples[:, j], weights=weights)
        wstd = np.sqrt(np.average((samples[:, j] - wm) ** 2, weights=weights))
        per = [np.average(r.samples[:, j], weights=r.weights) for r in results]
        print(f"[pool] D_iceIh_km pooled wmean={wm:.3f} wstd={wstd:.3f}; "
              f"per-seed wmeans={[f'{v:.3f}' for v in per]} "
              f"between-seed std={np.std(per):.3f}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    pooled.save(str(OUT))
    print(f"[pool] saved -> {OUT}")


if __name__ == "__main__":
    main()
