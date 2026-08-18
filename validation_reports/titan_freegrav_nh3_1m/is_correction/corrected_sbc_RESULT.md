# Corrected-pipeline SBC (§0.28 authorized) — RESULT: PASS

Full report: `corrected_sbc_report.json`. Preregistration:
`corrected_sbc_PREREGISTRATION.md`. Driver: `plans/scripts/
nh3_corrected_sbc.py`. Run wall time: 101,518s (~28.2h serial-equivalent,
~14h wall at --n-jobs 8, launched 2026-08-16 16:09, finished 2026-08-17
~20:21 local).

## Headline

**VERDICT: PASS.** All 13 sampled parameters pass the per-parameter
weighted-rank KS uniformity test both raw and after BH-FDR correction —
the raw minimum p-value is 0.171 (log10_zeta), nowhere near the 0.05
threshold even before multiplicity correction, so there is no marginal
case the correction is doing real work on.

| param | ks_pval | ks_pass | bh_adj | ks_pass_bh |
|---|---|---|---|---|
| alpha | 0.650 | PASS | 0.869 | PASS |
| log10_zeta | 0.171 | PASS | 0.741 | PASS |
| log10_eta_Ih | 0.674 | PASS | 0.869 | PASS |
| log10_eta_III | 0.211 | PASS | 0.741 | PASS |
| log10_eta_V | 0.242 | PASS | 0.741 | PASS |
| log10_eta_VI | 0.879 | PASS | 0.937 | PASS |
| log10_eta_sil | 0.937 | PASS | 0.937 | PASS |
| Tb_K | 0.287 | PASS | 0.741 | PASS |
| log10_wOcean_ppt | 0.407 | PASS | 0.741 | PASS |
| R_core_km | 0.735 | PASS | 0.869 | PASS |
| rho_core_kgm3 | 0.456 | PASS | 0.741 | PASS |
| dC20_nh | 0.367 | PASS | 0.741 | PASS |
| dC22_nh | 0.445 | PASS | 0.741 | PASS |

## Coverage bookkeeping (preregistered before this result existed)

- 900 pairs requested from `generate_training_set` → 613 drawn (31.9%
  rejected: 193 support-guard + 94 non-finite — matches the ~30%
  expected from the committed uncorrected SBC report almost exactly).
- 613 → 485 survived the per-pair Pareto-k<=0.7 / ESS>=1000 reliability
  filter (128 excluded, 20.9%) — under the 25% INCONCLUSIVE-BY-COVERAGE
  threshold, close to the ~18% expected from the C12 sweep precedent.
- **n_sbc_used = 485, 15 short of the n_sbc_target=500** stated in the
  §0.28 ruling. The preregistration fixed a rule for >25% exclusion
  (INCONCLUSIVE-BY-COVERAGE) but did not fix a rule for a sub-500,
  sub-25%-excluded shortfall — recording this honestly rather than
  silently treating 485 as "close enough to 500." Given the margin on
  every raw p-value (min 0.171, next lowest 0.211), 15 additional pairs
  would not plausibly change the verdict; flagging the shortfall for
  the record, not as a reason to withhold the PASS.
- Excluded-pair z-distance (from the artifact's x_norm mean, RMS across
  the 4 observables) is systematically higher than survived-pair
  z-distance: excluded median 1.28 (max 4.11) vs survived median 0.78
  (max 1.40) — same qualitative pattern as the C12 sweep (failures
  concentrate away from the training center), now confirmed on
  genuine SBC ground-truth pairs rather than only prior-predictive x.

## What this establishes

The corrected pipeline (IS-reweighted amortized posterior) is
well-calibrated across its own effective prior — for a pair that
survives the reliability filter, the weighted-rank statistic is
uniform, parameter by parameter, with no multiplicity-corrected
failures. This is the C14/C15 SBC condition from the tidal-sector
remedy plan, now closed for NH3. Per that plan's preregistered
reading: "gates PASS validates the corrected PIPELINE, not the flow" —
this result says the IS-correction machinery (weights, resampling,
reweighted posterior) does what it claims on the region of x-space
where it survives its own preregistered reliability gates. It says
nothing about the ~21% of pairs excluded by those gates, which is
exactly the C12 finding already reported (`c12_sweep_FINDINGS.md`) —
this SBC result does not resolve or supersede that.

Not self-adjudicated. The scientific-reviewer interprets what this
means for the tidal quarantine-lift path, alongside the still-open C12
finding.
