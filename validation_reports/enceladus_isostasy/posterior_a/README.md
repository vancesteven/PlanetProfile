# Run (a) full posterior — 2026-10-26 (Mac, repo 060ef826)

Config: enceladus_cassini_isostasy_ocean6D_rhorock.json, seed 75, pocoMC n_effective 8000 / n_active 4096 /
n_posterior_samples 20000, production cache ~/scratch/enceladus_zbw_production_cache.pkl.
Launcher: `plans/cluster/enceladus_posterior_a.sh`. Summary + figure: `plans/scripts/enceladus_posterior_a_report.py`.
Result pkls and the 88 MB log are local only (untracked).

## Plumbing and sampling (verified)
- Preflight PASS: 884/1000 prior draws accepted, 0 exceptions, 5.7 ms/call. **All 116 rejects are
  `zbw_cell_has_None_corner`, and all have zb 5.1–13.0 km** (any w). This is the support gap left by the
  121 failed builds. The prior is effectively truncated at thin shells; the posterior (zb ≈ 25–28 km) is far from it.
- rep1 and rep2 (same seed): samples, log-likelihoods and weights **bit-identical**. 180,224 calls and 1511 s each.
- 21,271 weighted samples, Kish ESS 16,213, log Z = −7.096 ± 0.009. pocoMC gives importance-weighted
  independent samples, so there is no chain autocorrelation to report.
- tidal_solve_skipped = True; isostatic_derived finite on 100% of draws; max |mass_residual_after| = 1.2e-15.

## Posterior (rep1; 16/50/84, weighted)
| quantity | 16 / 50 / 84 | sd | Fisher C σ (best fit) |
|---|---|---|---|
| zb_km | 25.87 / 26.76 / 27.74 | 0.94 | 0.79 |
| log10 w (ppt) | prior-like | 0.87 | 0.87 (prior) |
| rho_ice | 918.2 / 924.9 / 931.8 | 5.78 | 5.77 (prior) |
| rho_rock | 2373 / 2411 / 2453 | 40.2 | 37 |
| C2 | 0.718 / 0.814 / 0.906 | 0.091 | 0.086 |
| libration_sys_frac | −0.0039 / 0 / 0.0039 | 0.0039 | prior |
| C/MR² (derived) | 0.3306 / 0.3324 / 0.3340 | — | 0.0016 |
| R_s (km) | 190.4 / 192.4 / 194.3 | | |
| ocean thickness (km) | 30.8 / 32.9 / 35.2 | | |
| R_s − PP core radius (km) | −4.9 / −2.9 / −1.1 | | |

Max log-likelihood −0.52 at zb 26.47, log10 w 0.28, rho_ice 920.6, rho_rock 2402, C2 0.81.
Correlations: rho_s–C2 −0.877; zb–log10 w 0.20; zb–C2 0.12; zb–sys 0.11; rho_ice–rho_s 0.14.

## zb jointly with salinity (`zb_w_joint.png`)
| log10 w bin | mass | zb 16/50/84 | corr(rho_s,C2) |
|---|---|---|---|
| [−1, 0.5) | 0.498 | 25.76 / 26.61 / 27.56 | −0.894 |
| [0.5, 1.5) | 0.331 | 25.87 / 26.74 / 27.68 | −0.892 |
| [1.5, 2] | 0.171 | 26.33 / 27.26 / 28.24 | −0.882 |

The zb median is flat at ≈26.6 km for log10 w < 1. It rises to 27.45 km in [1.75, 2). The total shift
(≈0.85 km) is under 1σ.

## Fisher comparison
The w-sweep Fisher (variant C) predicted that corr(rho_s,C2) weakens at high salinity
(−0.887/−0.872/−0.752/−0.388 at log10 w −0.5/1.0/1.5/1.9). The posterior shows no weakening within w bins.
Conditioning the same Fisher matrices on w (dropping the w row/column) gives −0.887/−0.888/−0.888/−0.887
and σ(zb|w) 0.79/0.78/0.77/0.74 km. Those match the posterior's within-bin values: corr −0.88 to −0.90, and
σ(zb) ≈ 0.91, a ratio of ~1.17. The w-marginal weakening came from the local steep dzb/dw slope being
extrapolated across the whole w box. It is not a property of the posterior.

## Status
Sampling, reproducibility and closure: `verified` (artifacts above).
Scientific acceptance: `implemented, unverified` as a paper result.

Reviewer (2026-10-27): **PASS WITH CONCERNS.** All comparison criteria are met (σ ratios vs Fisher C:
zb 1.19, rho_rock 1.07, C2 1.05, C/MR² 1.05), and the Fisher-correlation resolution above is confirmed.
- MAJOR: same-seed bit-identity shows determinism, not Monte-Carlo error. A different-seed replicate is
  needed before `verified`: `SEED=2027 REPS=1 plans/cluster/enceladus_posterior_a.sh`.
- The None-corner hole (zb ≤ 13 km) has negligible effect on (a). Confirm with
  `THIN_CHECK=1 plans/cluster/enceladus_posterior_a.sh`; the best loglike at zb 13–15 km should be far below −50.
  log Z is over the full prior with L = 0 in the hole (not renormalised by the 0.884 acceptance).
  [2026-10-27, run (h)] Correct as stated: this IS the nominal-prior evidence and is not renormalised. 0.884 is the 1,000-draw preflight (± 0.010); the 20k-draw servable fraction is 0.893. See ../prior_h/README.md.
- The hole WILL distort the prior-only run (h) and SBI prior draws. Fill the corners, or tally the
  hole separately from physical rejects.
- Wording: "zb = 26.8 +1.0/−0.9 km (16–84%), marginalised over an unconstrained ocean salinity,
  log10 w ~ U[−1, 2] (w in ppt). Conditional on w, the median is 26.6 km for w ≲ 10 ppt and rises to
  about 27.3–27.5 km for w ≳ 30 ppt; the salinity-prior dependence (≲0.9 km) is below the 1σ width.
  Salinity and ice density are not constrained by these data." Qualify it with the ocean branch,
  isostatic_hm2019, uniform rock ρ_s ~ U[2200, 2600], and libration 0.092° (Park erratum; source not yet in papers/).
  Do not present C/MR² as a measurement until (h) shows the prior pushforward is much wider.
- Next, in order: seed replicate + thin check; (g) w uniform in ppt (expect zb +0.5–0.8 km);
  (h) via direct prior MC; injections at log10 w ≈ 1.7 and rho_ice 918/933; rho_ice prior variants; Fisher at the MAP.

## Different-seed replicate and thin check (2026-10-27): both PASS, reviewer MAJOR closed

**Different-seed replicate (seed 2027).** Produced with `SEED=2027 REPS=1 plans/cluster/enceladus_posterior_a.sh` (1472 s). The result is in `../posterior_a_seed2027/`.
- Kish ESS: 16,061 vs 16,213 for seed 75.
- log Z: −7.1106 ± 0.0065 vs −7.0964 ± 0.0088. The difference is 0.014, within 1.3 combined σ.

Per-parameter comparison, seed 75 vs seed 2027:

| param | 16/50/84%, seed 75 | 16/50/84%, seed 2027 | Δmean/SE | Δmean/sd | sd ratio | weighted KS |
|---|---|---|---|---|---|---|
| zb_km | 25.87 / 26.76 / 27.74 | 25.87 / 26.78 / 27.76 | +1.25 | +0.014 | 1.009 | 0.0085 |
| log10_wOcean_ppt | −0.52 / 0.51 / 1.53 | −0.50 / 0.52 / 1.55 | +1.35 | +0.015 | 0.996 | 0.0100 |
| rho_ice_kgm3 | 918.2 / 924.9 / 931.8 | 918.2 / 924.9 / 931.8 | −0.07 | −0.001 | 0.997 | 0.0084 |
| rho_rock_kgm3 | 2372.9 / 2410.7 / 2453.2 | 2373.5 / 2410.7 / 2453.3 | +0.13 | +0.001 | 1.001 | 0.0091 |
| compensation_C2 | 0.718 / 0.814 / 0.906 | 0.718 / 0.815 / 0.908 | +0.76 | +0.009 | 1.009 | 0.0085 |
| libration_sys_frac | −0.0039 / 0 / 0.0039 | −0.0039 / 0 / 0.0040 | +0.34 | +0.004 | 1.009 | 0.0064 |

- Every mean shift is ≤ 1.35 Monte-Carlo SE and ≤ 0.015 posterior sd.
- Every sd ratio is within 1% of 1.
- Every weighted KS distance is below the 1% two-sample critical value, about 0.018 at the combined ESS.
- The zb headline shifts by 0.015 km (seed 2027 − seed 75), far below its ~0.9 km width.

So the Monte-Carlo error is negligible against the posterior widths.

**Thin check.** `THIN_CHECK=1` drew 5000 prior draws at zb 13–15 km. All were accepted (0 rejects, 0 exceptions). The best loglike was **−318.2**, far below the reviewer's −50 bar and the central −3.67. The None-corner hole at zb ≤ 13 km sits where the likelihood is negligible, so it does not affect (a).

Status of the run (a) posterior:
- Sampling, reproducibility and Monte-Carlo error: **`verified`**. The evidence is the table above and `thin_zb13_15.json`.
- Scientific acceptance as a paper result: still `implemented, unverified`, pending (g), (h) and the remaining reviewer items. The hole still matters for (h) and for SBI prior draws.
