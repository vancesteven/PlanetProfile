# Corrected-pipeline SBC — preregistration (before any result is read)

Spec source: `plans/MACHINE-B-HANDOFF.md` §0.17 item 6 (canonical wording,
C14/C15) and §0.28 (manager ruling, 2026-08-16) answering
`corrected_sbc_SPEC_QUESTIONS.md`. This document fixes every remaining
design choice before the driver runs, per this campaign's "never tune
after seeing a failure" discipline — same posture as
`c12_sweep_PREREGISTRATION.md`.

## 1. Ruling inputs (§0.28, binding, not re-derived here)

1. N = 20,000 (deployed/fiducial N — the estimator under validation is the
   deployed one). n_sbc = 500 SURVIVING pairs, NH3 only. Budget corrected
   to ~52 CPU-h; run in gaps.
2. theta_0, x pairs drawn via `SBIRunner.generate_training_set` (effective
   prior, all training support cuts) — CONFIRMED as C15, exactly the C12
   sweep's prior-predictive construction.
3. Low-ESS/Pareto-k policy: pairs failing Pareto-k<=0.7 OR the ESS floor
   are RECORDED, EXCLUDED from the rank statistics, and their fraction is
   REPORTED with the excluded pairs' x z-distance distribution (same
   z-metric as `c12_sweep_report.json`, distance from the artifact's
   `x_norm` mean in std units). Expected ~18% (C12 precedent). If excluded
   fraction > 25%: verdict = INCONCLUSIVE-BY-COVERAGE (not FAIL), routed
   to the C12 scope ruling rather than treated as a calibration failure.
4. Deploy-time N is fixed at 20,000 (this ruling). The GUI's n_derived=1000
   default cannot produce IS weights — out of scope here, queued as a
   Codex guard task by the manager, not a B action item.

## 2. Over-request target

n_sbc target is 500 SURVIVING pairs (after the ~30% support-guard
attrition AND after the ~18% expected Pareto-k/ESS exclusion). Solving
0.70 * 0.82 * n_request >= 500 => n_request ~= 872. Request **900** pairs
from `generate_training_set` up front (round number, small margin over
872) rather than iterating in batches — a single over-request avoids a
second seed-management decision. If fewer than 500 pairs survive both
filters from 900, the shortfall is reported explicitly (never silently
topped up with a second draw at a different seed, which would break the
single-preregistered-draw discipline).

## 3. Rank statistic (C14, weighted-rank "option i")

Per pair i (theta_0_i, x_i): condition the deployed flow on x_i with full
recompute (`n_derived=None`), form IS weights via
`compute_is_correction` exactly as the fiducial/C12 driver does, then the
weighted rank for parameter j is:

    r_ij = sum_k w_ik * 1[ posterior_sample_ijk < theta_0_ij ]

over the posterior draws k=1..N for pair i, using the SAME normalized IS
weights w_ik (sums to 1 over k) already computed for that pair — i.e. the
weighted empirical CDF of the corrected posterior evaluated at the true
theta_0, per parameter. This is continuous on [0,1] (not an integer count
out of N), so `sbi.diagnostics.check_sbc` (which assumes integer ranks on
[0,N] and its own internal bootstrap) is NOT used. Per-parameter
uniformity test: `scipy.stats.kstest(r_j, 'uniform')` across the n_sbc
surviving pairs, exactly as recorded in `corrected_sbc_SPEC_QUESTIONS.md`
"what is settled."

BH-FDR multiplicity correction across the 13 sampled parameters
(`alpha, log10_zeta, log10_eta_{Ih,III,V,VI,sil}, Tb_K, log10_wOcean_ppt,
R_core_km, rho_core_kgm3, dC20_nh, dC22_nh`) reuses
`validate_sbi.py`'s existing BH block (lines ~534-542) refactored into a
shared helper — NOT duplicated. Verdict uses `bh_all_pass`; raw per-param
`ks_pass` is retained informationally, matching the existing SBC gate's
own precedent for the reason recorded there (a single param grazing 0.05
is an expected false positive at K=13 under perfect calibration, not a
finding).

## 4. Per-pair recipe (reuses the C12 per-point-rebuild fix)

For each surviving (theta_0_i, x_i):
1. Locate + load the artifact against the **unmutated fiducial** config
   hash (same as the C12 driver) — artifact identity must not depend on
   the swept x.
2. Build a FRESH `InferenceConfig` with `observables` centrals set to
   x_i (sigmas unchanged), fresh `SBIRunner`, inject the loaded
   `posterior`/`_train_info` — this is the exact fix that closed the
   MCMCRunner likelihood-closure trap identified before the C12 sweep
   ran (the closure is bound to `config.observables` at construction and
   does not track a swept `x_obs` argument).
3. `runner._condition_and_package(x_obs=x_i, n_posterior_samples=20000,
   n_derived=None, ...)` — full recompute, required per Finding 1 of the
   spec-questions doc (no partial-derived path exists).
4. `compute_is_correction(result, k2_support_bounds=...)` for weights,
   Pareto-k, ESS. Record `pareto_pass = (pareto_k<=0.7) and not nan`;
   `ess_pass = (ess>=1000)` (absolute floor per the reliability set in
   §0.17, not ESS/N).
5. Compute r_ij per §3 using `result.samples` (theta draws, columns in
   `runner.param_names` order) and the pair's weights.
6. `c2_hash_expected_mismatch` recorded explicitly (expected True
   off-fiducial), never silently skipped — same as C12.

## 5. Verification before the real run

Before the 900-pair draw, the per-pair point-builder is run on ~20
fixed-seed pairs and checked for:
- Exact reproduction of the committed fiducial Pareto-k/ESS
  (`is_validation_nh3.json`) when x_i is forced to the fiducial x and the
  same seed as that report — same check as the C12 sweep's
  `verify_at_fiducial()`, reused verbatim (it already exists as
  importable logic in `nh3_c12_amortized_sweep.py`).
- Non-degenerate rank values (not all 0 or all 1) on a small synthetic
  check where theta_0 is drawn from the flow's own posterior at x_i
  (should give ranks ~ Uniform under a correct, non-buggy implementation
  even before any real calibration claim is made) — this is a
  software-correctness smoke test on the r_ij formula itself, not a
  calibration result, and is reported separately from the real SBC.

If either check fails, the driver does not proceed to the 900-pair draw.

## 6. Seeds

- Pair draw (theta_0, x): `generate_training_set(n_simulations=900,
  seed=530170, obs_noise=True, noise_seed=530171)` — dedicated seeds, not
  reused from C12's `PRIOR_SEED=90210` (a distinct preregistered draw).
- Per-pair posterior sampling seed: `531000 + i` (i = pair index in the
  order returned by `generate_training_set`), deterministic and
  reproducible.

## 7. Deliverables

- Per-pair record: theta_0, x, pareto_k, ess, pareto_pass, ess_pass,
  survived (bool), z-distance from x_norm mean, wall_s.
- Excluded-pair fraction + z-distance distribution of the excluded set.
- Per-parameter: `ks_pval`, `ks_pass`, `ks_pval_bh_adj`, `ks_pass_bh`,
  n_sbc actually used.
- Verdict: PASS (bh_all_pass and excluded_frac<=25%) / FAIL
  (bh fails on >=1 param and excluded_frac<=25%) /
  INCONCLUSIVE-BY-COVERAGE (excluded_frac>25%, regardless of the KS
  numbers — coverage governs before calibration is even asked).
- This document, committed alongside the report, so the design record
  predates the result (matching C12's own discipline).

Not self-adjudicated: PASS/FAIL/INCONCLUSIVE is reported per the rule
above; the scientific-reviewer interprets what the result means for the
quarantine-lift path. Gates are interpreted, never tuned to pass.
