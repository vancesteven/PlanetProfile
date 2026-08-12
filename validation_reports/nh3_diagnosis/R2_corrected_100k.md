# §0.20 R2 — N=100k corrected conditioning at the fiducial (tighter corrected SE)

**Status: complete — scientific-reviewer verdict PASS-WITH-CONCERNS, proceed to
R3 (2026-08-11).** The reviewer confirmed the "persists" determination is
*decisive*: the corrected fraction moved 0.92925 → 0.93299, i.e. AWAY from the
reference at higher N, which *falsifies* (not merely fails to support) the
finite-N corrected-side-bias hypothesis — any O(1/ESS) SNIS bias would relax
toward the reference, and the flow-corrected fixed point (~0.933) lies above it.
R2 consumes no gate; C16 status unchanged. Three carry-forward concerns folded
in below (comparator asymmetry; between-seed corrected SE for the R3 gate;
delta-method SE cross-check) — none alter routing. Second step in the C16
resolution plan (Machine A §0.20). R1 (PASS-WITH-CONCERNS,
2026-08-11) established that neither side's *legitimate* weighting treatment
explains the residual and that the corrected side is stable. R2 tightens the
**corrected-side** standard error at the fiducial by scaling the deployed-flow
conditioning from N=20k (R1) to **N=100k**, and asks whether the fiducial
residual shrinks (→ 1/ESS-ceiling / finite-N-bias explanation) or persists (→
tension is not corrected-side bias; proceed to R3). Diagnostic only — consumes
no gate, does NOT re-ratify C16. R3 (n_eff≥8000 reference recompute) remains
decisive. The scientific-reviewer, not this session, adjudicates the reading.

Script (repo copy): `plans/scripts/r2_corrected_100k.py`. Output:
`validation_reports/nh3_diagnosis/r2_corrected_100k.json`. Reuses the exact R1
SNIS conditioning machinery (`_condition_and_package` → `compute_is_correction`,
seed 72, fiducial x, ocean mask `nan_to_num(D_ocean,0)>0.5`), scaled to 100k,
plus a 2000-rep bootstrap SE on the weighted corrected fraction.

## Result

| quantity | R1 (N=20k) | **R2 (N=100k)** |
|---|---|---|
| corrected ocean fraction | 0.92925 | **0.93299** |
| bootstrap SE | — | **0.00197** |
| residual vs reference A (0.91725) | +0.01200 | **+0.01573** |
| residual / bootstrap SE | — | **7.98** |
| Pareto-k (tail) | 0.191 | 0.289 (< 0.5 gate; PSIS off) |
| Kish ESS | — | 5976 |
| n_finite / n_rejected | 19774 / 226 | 98837 / 1163 |

Bootstrap 95% CI on the corrected fraction: **[0.92902, 0.93674]** — the lower
bound (0.92902) sits **+0.0118 above** reference A (0.91725).

## Preregistered reading (for the reviewer — NOT self-adjudicated)

Manager's R2 spec (§0.20): *"residual PERSISTS at tighter SE → confirmation the
tension is not finite-N corrected-side bias, proceed to R3; if it SHRINKS
materially → the 1/ESS-ceiling argument fails, the SNIS/corrected side reopens,
R3 pauses."*

Observed:

1. **The residual did not shrink.** At 5× the draws it moved +0.01200 → +0.01573
   (slightly *up*, not toward zero). The corrected fraction rose 0.92925 →
   0.93299 — consistent with the C13 seed values (0.9293/0.9310/0.9363) the
   reviewer cited, not with a finite-N bias that would relax toward the
   reference at higher N.
2. **The residual is now 8.0 bootstrap SE above zero**, and the entire 95% CI
   ([0.92902, 0.93674]) lies above reference A. The corrected side is a
   precisely-resolved fixed point, not a noisy estimate that happens to sit
   high.
3. **The tail did not turn heavy.** Pareto-k rose only 0.191 → 0.289, still well
   below the 0.5 clean gate — PSIS would still not fire. The corrected fraction
   is not being propped up by a few extreme weights (consistent with R1's
   top-0.1% ablation, which removed the heaviest mass and still left the
   residual standing).

**Advisory recommendation (reviewer's call, not this session's):** R2 lands on
the preregistered "residual persists / tightens" branch. It does not reopen the
corrected/SNIS side and does not re-ratify C16. It supports proceeding to R3
(the decisive reference recompute at n_eff≥8000, ≥3 seeds), where the
outstanding suspect — the reference-side seed scatter (R1: span 0.0121, std
0.0063 at n_eff=2000) — is resolved and the |corrected − reference| ≤ 2×
combined-SE re-ratification gate is evaluated. Per the R1 carry-forward, R3
should report the pooled reference SE (≈0.0063/√3 ≈ 0.0036) under treatment-A
pooling so the gate compares like with like.

**Scope caveat (carried from R1):** this residual is the **fiducial**
single-seed residual (corrected − reference A). The committed **pooled** 3-seed
residual (+0.0149) is R3's remit; R2 only tightens the corrected-side SE at the
fiducial.

## Reviewer corrections folded in (2026-08-11, PASS-WITH-CONCERNS → proceed to R3)

- **Comparator asymmetry — the +0.01573 headline is NOT "the tension grew past
  0.0149."** It compares seed-72 (single-seed, high-side) corrected against the
  pooled-3-seed (low-side) reference REF_A. Using R1's *seed-72* reference
  (0.92223) instead, the like-with-like seed-72 residual is **+0.01076** —
  *below* the committed pooled +0.0149. The fiducial residual exceeding the
  pooled value is a comparator artifact (single-seed-high corrected vs
  pooled-low reference on both ends), not evidence the tension grew. This is the
  answer to §0.20 R2-question 3: the report's scoping is correct and the
  fiducial-exceeds-pooled observation changes nothing.
- **The single-seed bootstrap SE (0.00197) understates corrected-side
  uncertainty.** It captures only within-conditioning Monte-Carlo variance at
  one flow seed. The honest corrected-side SE includes between-seed scatter
  (C13: sample std 0.00365, ≈2× the bootstrap SE). The residual is still 4.3 SE
  even against the between-seed std, and all three C13 corrected values sit above
  the reference — so routing is unchanged — but **R3 must form the corrected-side
  SE from between-seed pooling (≥3 flow seeds)**, not the single-seed bootstrap.
- **Bootstrap SE < ESS-based delta-method proxy** (0.00197 vs
  √(p(1−p)/ESS)=0.00323). The nonparametric bootstrap is the consistent estimator
  for a self-normalized IS fraction (the ESS-binomial heuristic is conservative
  when indicator and weights correlate); the smaller value is defensible.
  Routing is robust either way (8.0 SE bootstrap / 4.9 SE ESS-based). R3 should
  report both as a cross-check.

## R3 gate specification (reviewer required-validation, carried forward)

1. Compute BOTH sides at n_eff≥8000, ≥3 seeds; evaluate
   |corrected_pooled − reference_pooled| ≤ 2× combined SE with **like-with-like
   treatment-A pooling on both sides**.
2. combined SE = √(SE_ref² + SE_corr²), BOTH between-seed pooled SEs. Carry-forward
   estimates: SE_ref ≈ 0.0063/√3 ≈ 0.0036 (R1); SE_corr ≈ 0.00365/√3 ≈ 0.0021
   (C13). Illustrative 2× combined ≈ **0.0083** — the committed pooled residual
   0.0149 exceeds this, so unless R3's higher-n_eff reference recompute moves the
   reference materially upward, **the gate fails and the tension is real.** Do
   NOT substitute the single-seed bootstrap SE for the corrected side.
3. Report the delta-method SNIS SE alongside the bootstrap SE as a cross-check.
