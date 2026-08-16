# C12 amortized sweep — result and interpretation request

Raw report: `c12_sweep_report.json`. Design fixed in advance:
`c12_sweep_PREREGISTRATION.md`. Not self-adjudicated — routed to the
scientific-reviewer.

## Headline

**POOLED: 29/149 failed (19.5%) -> FAIL** against the preregistered
threshold (<=5%). 141/200 requested prior-predictive draws survived the
support guard (the driver reports this; it is expected sampler behavior,
not a defect).

- prior-predictive-only: 26/141 failed (18.4%)
- endpoints-only: 3/8 failed (37.5%), including one hard C4 REJECT
  (75.6% sentinel-rejected at the C22 hi endpoint — recorded as a failure,
  not silently excluded, per the preregistration's NaN/error handling rule)

This is a NEGATIVE result for amortization breadth as specified. It does
NOT bear on the fiducial-point validations (C3/C5.3/C10/C11/C13/C16),
which already passed and are unaffected by this sweep.

## The failures are NOT randomly scattered — they concentrate at high z-distance from the training center

Computed each point's Euclidean z-distance in the artifact's x_norm space
(the same normalization the flow trains on):

| | mean z-distance | median z-distance | n |
|---|---|---|---|
| prior-predictive PASSED | 1.59 | 1.61 | 115 |
| prior-predictive FAILED | 3.48 | 2.89 | 26 |

The 10 worst Pareto-k values (up to 16.27) all sit at z-distance 3.3-7.5,
overwhelmingly driven by extreme Re_k2 and/or Im_k2 (individual
per-channel z as high as +7.2). These are legitimate prior-predictive
draws (theta ~ effective prior, forward-modeled, noised) — not constructed
adversarially — but they land in the tails of the JOINT x-distribution the
flow was trained on, where amortized posterior mass and true likelihood
mass can diverge sharply.

**This is exactly the mechanism the ratified `x_obs_limits` guard exists
for.** The deployed GUI hard limit is `Im_k2 in (0.0, 0.30)`
(`PlanetProfileApp/pages/Inference.py`), and the C12 endpoint construction
independently derived a near-identical value (Im_k2 hi = 0.3008) from the
±2σ_train box — see the preregistration. So the deployed conditioning
domain is ALREADY narrower than the full prior-predictive support this
sweep probes.

## Endpoint failures, itemized

- **C22 hi**: HARD REJECT (75.6% sentinel-rejected, C4). This point is
  ~46.7σ_obs from the fiducial C22 (per the earlier reviewer measurement
  in `b2prime_ADJUDICATED_drho_weighting.json`'s sibling analysis of this
  artifact's x_norm) — an extreme corner even by the ±2σ_train construction.
- **Re_k2 lo** (0.204): Pareto-k 0.948, fails by a small margin (threshold
  0.7).
- **Im_k2 hi** (0.301): Pareto-k 1.113 — just outside the deployed 0.30
  GUI limit and fails.
- The other 5 endpoints (C20 lo/hi, C22 lo, Re_k2 hi, Im_k2 lo) all PASS
  cleanly (Pareto-k 0.36-0.62).

So endpoint failures concentrate on the Re_k2/Im_k2 axes specifically —
consistent with the k2 channel being the one already under a separate
support-guard/quarantine discussion in this campaign, not a new finding
about C20/C22.

## ESS/N and amortization pricing (reported per the task spec, not gated)

- ESS/N distribution across all 149 points: min 0.00052, median 0.287,
  max 0.563.
- N_required (at the 1000-ESS floor) distribution: min 1776, median 3487,
  max 1,915,498 (driven by the same tail points that failed Pareto-k).
- At the fiducial, ESS/N was 0.0635 (N_required ~15,739) — the SWEPT
  points have a HIGHER median ESS/N (0.287) than the fiducial itself,
  because most swept points are drawn from the effective prior (closer to
  the bulk of the flow's own training distribution) rather than the
  specific, somewhat off-center fiducial observation. The amortization
  price is therefore data-point-dependent and the fiducial is not
  representative of the median cost — it is worse than the median but far
  better than the worst-case tail.

## 0.9326 sanity line (task spec item)

Across the 148 points with recoverable branch data, the ocean-branch
flow-vs-corrected probability moves: mean |flow - corrected| = 0.0115,
MAX |flow - corrected| = 0.127 (12.7 percentage points). This is a much
larger swing than the fiducial's own 0.93255 -> 0.92925 (0.33 pp), already
cleared as real reweighting under §0.18 Phase-1 item 2. The fiducial's
near-fixed-point behavior does NOT generalize across x — off-fiducial the
branch correction can be substantial, which is expected (the branch mass
depends on where x sits relative to the ocean/no-ocean boundary in
x-space) and is not itself evidence of a defect, but it means the fiducial
plumbing check's stability does not certify branch-correction stability
elsewhere.

## What this does and does not mean — reviewer questions

1. **Does the C12 FAIL block the NH3 quarantine-lift path**, or does it
   only bound the domain over which amortized IS correction can be
   trusted (i.e. informative about scope, not a blocker on the deployed
   fiducial-adjacent conditioning that GUI users actually exercise)?
2. **Is the ±2σ_train endpoint construction the right domain to gate on**,
   given it already reproduces the deployed Im_k2 hard limit almost
   exactly — should the C12 pass/fail be restricted to the IN-SUPPORT
   region the GUI already permits (`x_obs_limits`), rather than the full
   ±2σ box which admits some combinations (e.g. Re_k2 hi=0.854, Im_k2
   hi=0.301, simultaneously) the deployed app would never actually let a
   user request?
3. **Should the 200-prior-predictive-draw denominator be reported at
   n=141** (what the support guard actually delivered) or does the
   guard's 30% attrition itself need investigation before the sweep is
   considered representative?
4. Given the failures concentrate in the Re_k2/Im_k2 tail exactly where
   this campaign's separate k2-quarantine discussion already lives, is
   this FAIL evidence FOR extending that quarantine's scope, or is it a
   separate, amortization-specific finding that should be tracked on its
   own?

Not tuned, not re-run with a different threshold after seeing this
result — reported as designed. Gates are interpreted, never tuned to
pass.
