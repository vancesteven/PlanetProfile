# Corrected-pipeline SBC (§0.17 item 6 / §0.27 standing work) — spec questions before compute

Per §0.27, NH3 corrected-SBC is my only standing work item while the
Enceladus build trigger is rescinded (r6 NOT RATIFIED). Before building a
driver and spending real compute, sent the underspecified spec to an
independent Opus pass — same discipline as the C12 sweep, which caught a
real design trap before it burned hours. This one found the budget itself
is inconsistent with the stated method, plus a factual correction to my
own assumption about the NH3 config. Routing to Machine A for a ruling
rather than picking defaults myself, per direction to conserve compute and
only run heavy tasks as directed.

## Task as stated (`plans/MACHINE-B-HANDOFF.md` line 451-454)

"SBC of the corrected pipeline (C14/C15): rank definition = weighted ranks
r = sum_i w_i 1[theta_i < theta_0] vs Uniform(0,1); theta_0 from the
EFFECTIVE prior (all training support cuts); n>=500 pairs at the
deploy-time N; BH-FDR across 13 params. Budget ~14 CPU-h/comp."

## Finding 1 — the budget is wrong as stated, not just tight

A full `_condition_and_package(..., n_derived=None)` recompute is
structurally required per SBC pair (compute_is_correction hard-fails on
any NaN-padded log_likelihood, and n_derived=None is the only way to avoid
padding). Cost is empirically linear and consistent across four existing
runs (fiducial, two C13 seeds, R2 100k) at ~8.8 ms/work-unit
(work_units = 2N + n_reeval), plus ~19s fixed overhead/pair (runner
construction, cache unpickle, artifact load).

At the fiducial N=20,000: **~52 CPU-h for 500 pairs — 3.7x the stated
~14 CPU-h budget.** The budget is only satisfiable at N<=4,400. My own C12
sweep data (149 points, N=5000, 4.62 CPU-h) independently confirms this
arithmetic to two digits (500 pairs at N=5000 extrapolates to 15.5 CPU-h).

A leaner likelihood-only recompute path exists (SBC needs only
log_likelihoods + flow_log_prob, not k2/heating/derived quantities;
skipping those roughly halves cost) but is not reachable through existing
knobs — n_derived is the only lever and any value below N triggers the
same hard-fail. Building that path would be new code on the
Track-1-correction-critical conditioning function, not a preregistration
choice.

**Decision needed:** pin N, and rule whether to (a) accept the true cost
at the chosen N. (b) raise the budget, or (c) authorize building the
lean likelihood-only path first.

## Finding 2 — correction to my own prior assumption: NH3 DOES have active support cuts

I had assumed the NH3 joint (no-ocean+ocean) design meant no support-guard
cuts applied, since the campaign deliberately drops
`phase_stability.enforce='no_ocean_Ih'` so frozen nodes are retained
rather than rejected. That assumption is WRONG: the config independently
sets `apply_support_guard: true`, and `_support_guard_active()` ORs the
two conditions rather than requiring both. The committed uncorrected SBC
report already shows this cost: 500 requested pairs -> 353 kept (29.4%
rejected) from the k2-support box, density-inversion guard, and
cache-servability checks.

Consequence: "theta_0 from the effective prior (all training support
cuts)" is NOT moot for NH3 as I'd assumed — theta_0 must be drawn via
`generate_training_set` (which applies the guard automatically), never a
raw prior box draw, and pairs must be over-requested (~750-800 for 500
survivors) to hit n>=500 after the guard's ~30% attrition.

One structural asymmetry flagged, not yet acted on: the k2-support box is
enforced on the training/generation side but does NOT appear in the
likelihood itself (`_make_flexible_log_likelihood` has no k2-box
rejection) — so the module docstring's claim that every guard is folded
into the likelihood sentinel is not exactly true for this one box. Small
expected effect (guarded theta draws land in the likelihood's own
near-zero tail anyway) but worth a reviewer note since SBC is precisely
the test that would notice a support/likelihood mismatch if one exists.

## Finding 3 — low-ESS pairs need a preregistered policy before results exist

SBC's rank resolution tracks effective sample size, not N. From the C12
sweep's own ESS distribution at N=5000: median ESS 1550, but 28% below
500, 12% below 50, 5th percentile ESS=12. A quick simulation under perfect
calibration shows the weighted-rank KS test's false-positive rate stays
near-nominal (4.5-6.5%) down to ESS~50, but inflates to ~24% at ESS=12.
~18% of pairs will fail the C12 Pareto-k>0.7 gate outright (that finding is
already pushed and under review, see c12_sweep_FINDINGS.md).

**Decision needed:** whether to exclude failing/low-ESS pairs (introduces
a scope-bound conditioning echo of §0.25's disposition), keep them
(conservative — biases toward FAIL, does not risk false PASS), or report
both. This must be fixed in writing BEFORE any pair is run, given how
directly the C12 outcome bears on it.

## Finding 4 — "deploy-time N" is not a well-defined quantity in the code

Candidates found: GUI default n_post=10000 (public-mode default 5000);
the fiducial validation's N=20000; GUI max 50000. The GUI's own
`n_derived` default is 1000, which CANNOT produce IS weights at all
(same C4 hard-fail) — so the literal current deploy configuration is one
where the correction is not computable, independent of any SBC question.
Worth flagging as its own finding regardless of the SBC ruling.

## What is settled, ready to build the moment a ruling lands

- Per-pair recipe: `generate_training_set` for (theta_0, x) pairs
  (over-requested for the support-guard attrition), then per pair rebuild
  a fresh config+runner with observables set to that pair's x (the same
  per-point-runner-rebuild fix already verified in
  `plans/scripts/nh3_c12_amortized_sweep.py`), full `n_derived=None`
  recompute, `compute_is_correction` for weights.
- Rank statistic: `scipy.stats.kstest(r, 'uniform')` per parameter, r
  computed via `scipy.stats.kstest`-compatible continuous weighted rank
  (NOT sbi's `check_sbc`, which assumes integer ranks on [0,N] and does
  not accept a continuous weighted statistic without abuse).
- BH-FDR: reuse `validate_sbi.py`'s existing block verbatim (the only
  BH-FDR implementation in the repo); refactor into a shared helper rather
  than duplicate.
- Pre-check before the real run: verify the driver reproduces
  `is_validation_nh3.json`'s committed fiducial Pareto-k/ESS exactly at
  n_pairs~20, same pattern as the C12 sweep's verify-only step.

Not self-adjudicated, no compute spent. Routing per instruction to
conserve tokens and run heavy tasks only as directed.
