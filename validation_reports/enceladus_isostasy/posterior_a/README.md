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
  log10 w ~ U[−1, 2] (w in ppt). Conditional on w, the median is 26.6 km for w ≲ 3 ppt and rises to
  about 27.3–27.5 km for w ≳ 30 ppt; the salinity-prior dependence (≲0.9 km) is below the 1σ width.
  Salinity and ice density are not constrained by these data." [2026-10-27: threshold corrected 10 → 3 ppt, matching the bin edge log10 w = 0.5; superseded by the adjudicated headline in ../injections/README.md.] Qualify it with the ocean branch,
  isostatic_hm2019, uniform rock ρ_s ~ U[2200, 2600], and libration 0.092° (Park et al. 2024 erratum, filed verbatim in papers/park2024global_ERRATUM.md, bcead860).
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
- Scientific acceptance: **scientifically ACCEPTED (ocean branch), Opus final-acceptance reviewer 2026-10-27.** The adjudicated headline and caveats are in `../injections/README.md`. The frozen-branch evidence comparison is a separate later step.

---

# Run (a′) re-run under the repaired libration — 2026-11-11 (Mac, repo b736860b)

Same config file, seed (75), sampler settings, cache and launcher as run (a),
but the libration channel now runs `formulation='vh2008_repaired'` × the
per-node elastic-factor sidecar F(zb,w)·(1 + c·(ρ_rock − ρ_ref)); the frozen
×1.008 factor is retired (commit 2a75a910; validation in
`../kint_switchover/README.md`). The launcher auto-built the production-cache
sidecar before sampling. THE PKLS IN THIS DIRECTORY ARE NOW THE RE-RUN'S
(the launcher reuses the directory); run (a)'s numbers survive in the tables
above and in `../posterior_a_seed2027/`.

## Plumbing and sampling
- Preflight PASS (same support hole as (a): thin-shell None corners only).
- rep1 and rep2 (same seed 75): **bit-identical posteriors** (same-seed
  determinism PASS; distinct pkl md5s are timing metadata only).
- 21,165 weighted samples, Kish ESS 16,130, log Z = −7.1127 ± 0.0062
  (run (a): −7.0964 ± 0.0088; Δ = −0.016, negligible against the
  branch-comparison margin log B ≤ −280).
- isostatic_derived finite on 100% of draws; max |mass_residual_after| = 1.2e-15.
- The posterior zb window (25.6–27.5 km, 16–84%) sits inside the F-table's
  best-validated window (zb 24–30: F 0.978–0.984, separability ≤ 0.11%),
  far from the thin-shell (zb < 13–14 km) degradation.

## Posterior (rep1; 16/50/84, weighted) vs run (a)
| quantity | re-run 16 / 50 / 84 | sd | run (a) 50% | shift |
|---|---|---|---|---|
| zb_km | 25.57 / 26.47 / 27.47 | 0.95 | 26.76 | **−0.29 km** |
| log10 w (ppt) | −0.52 / 0.52 / 1.55 (prior-like) | 0.88 | 0.51 | +0.01 |
| rho_ice | 918.2 / 924.8 / 931.7 | 5.75 | 924.9 | −0.1 |
| rho_rock | 2373 / 2411 / 2455 | 40.7 | 2411 | 0 |
| C2 | 0.712 / 0.811 / 0.904 | 0.092 | 0.814 | −0.003 |
| libration_sys_frac | −0.0040 / −0.0001 / 0.0039 | 0.0039 | 0 | ≈0 (prior-like) |
| C/MR² (derived) | 0.3306 / 0.3324 / 0.3341 | 0.0017 | 0.3324 | 0 |
| R_s (km) | 190.3 / 192.4 / 194.3 | 2.0 | 192.4 | 0 |
| ocean thickness (km) | 31.1 / 33.3 / 35.5 | 2.2 | 32.9 | +0.4 |

Max log-likelihood −0.518 (χ² ≈ 1.04 for 4 data) at zb 26.22, log10 w 0.02,
rho_ice 915.8, rho_rock 2394, C2 0.818, sys −0.0053.
Correlations: rho_s–C2 −0.879; zb–log10 w 0.24; zb–C2 0.14; zb–sys 0.12;
rho_ice–rho_s 0.15 — all within 0.04 of run (a).

zb by log10 w bin (vs run (a) medians 26.61 / 26.74 / 27.26):
| bin | mass | zb 16/50/84 | shift |
|---|---|---|---|
| [−1, 0.5) | 0.492 | 25.44 / 26.31 / 27.25 | −0.30 |
| [0.5, 1.5) | 0.334 | 25.58 / 26.44 / 27.41 | −0.30 |
| [1.5, 2] | 0.174 | 26.12 / 27.01 / 28.02 | −0.25 |

The shift is uniform across salinity bins; the joint zb–w structure is
unchanged.

## Reading (as adjudicated 2026-11-12)
The −0.29 km zb shift matches the switch-over prediction (−0.2..−0.3 km from
the +0.8–0.9% frozen-channel bias at C2 ≈ 0.81, 0.28 km per 1%). The clean
forecast comparison is the reviewer's weighted regression of zb on
(C2, ρ_s, w, ρ_ice, sys), (a′) against (a): the intercept at C2 = 0.81,
ρ_s = 2411 shifts **−0.278 km** (forecast −0.2..−0.3) and the C2 slope
changes by **−0.46 km per unit C2** (forecast −0.34..−0.38), so both the
magnitude and the C2 structure of the repair are confirmed. (Marginal C2-bin
shifts look flat because of the ρ_s–C2 confound, corr −0.88; do not quote
them against the forecast.) The w-bin shifts are −0.30 km for w ≲ 30 ppt and
−0.25 km above (consistent with F's w-dependence, ~2σ of per-bin MC noise).
C2, ρ_rock, C/MR², R_s move ≤ 0.05 sd; sd ratios 1.00–1.01; correlations
within 0.03: the repair moved the thickness scale, not the inference
geometry. (sys_frac stays prior-like, but that is uninformative here —
sys_frac and zb are degenerate on the single libration datum, so any
multiplicative change goes into zb regardless; its ±0.4% widens zb only
~0.11 km in quadrature.)

## Binding review (2026-11-12): GO-WITH-CONDITIONS
Every submitted number independently reproduced from the pickles (quantiles,
sds, ESS 16,130.4, log Z −7.112727 ± 0.006207, max-logL point, correlations,
salinity bins, C/MR², bit-identity of rep1/rep2, sidecar md5 vs F-table
provenance, preflight support unchanged at 884/116 None-corner). Run (a′) is
VALID as the headline posterior under the repaired model. Δ log Z −0.016 is
seed-level noise; the branch bound log B ≤ −280 survives (the frozen bound,
computed with 1.008 on the frozen branch, becomes ~1.8 nat tighter under the
switched config's factor 1.0 — see ../branch_bound/README.md note).

**Adjudicated headline (supersedes ../injections/README.md wording; effective
once condition C1 below passes):**

> Under the ocean-branch model (uniform-density rock interior,
> ρ_s ~ U[2200, 2600] kg m⁻³; Airy-isostatic Hemingway & Mittal (2019)
> gravity with Tajeddine et al. (2017) shape; Δρ-consistent Eq.-12 libration
> with the repaired Van Hoolst et al. (2008) rigid formulation and a per-node
> elastic-shell factor), conditioning on Cassini C20, C22 and C30 and the
> Park et al. (2024, erratum) physical libration of 0.092° ± 0.003° gives an
> ice-shell thickness zb = 26.5 +1.0/−0.9 km (median, 16–84%) for a
> log-uniform ocean salinity prior (log10 w ~ U[−1, 2], w in ppt), and
> 26.9 +1.0/−0.9 km if salinity is uniform in w over 0.1–100 ppt. zb is
> constrained jointly with salinity, which the data do not constrain:
> conditional on w, zb is 26.3 km for w ≲ 3 ppt and 27.0 km for w ≳ 30 ppt.
> The prior dependence (≤ 0.5 km) is below the 1σ width. Rock density and
> gravitational compensation are degenerate (corr = −0.88):
> ρ_s = 2411 +44/−38 kg m⁻³ and C2 = 0.81 ± 0.09 individually, and the ρ_s
> 95% interval [2344, 2502] kg m⁻³ lies well above the prior floor. The
> derived C/MR² is 0.3324 ± 0.0017; under the uniform-rock model this
> constraint is carried by the rock density (corr −0.98, ~96% of its
> variance) and is a model-dependent derived quantity, not an independent
> moment-of-inertia measurement. Ice density is not constrained: widening its
> prior from U[915, 935] to U[880, 940] kg m⁻³ moves zb by about +0.3 km
> (0.3σ).

> **Change note.** This supersedes the provisional zb = 26.8 +1.0/−0.9 km
> (October 2026). That value used the shipped VH2008 rigid libration with a
> fixed ×1.008 correction, which over-predicted the libration by ~0.9% at the
> posterior C2 and biased zb thick by 0.3 km (0.3σ). The posterior widths,
> the ρ_s–C2 degeneracy, C/MR² and the salinity structure are unchanged.

Supporting numbers (reviewer-computed, exact prior reweighting of (a′),
method validated on (a): reweighted-(a)-seed2027 27.155 vs direct g1 27.169):
w-uniform 25.98 / 26.92 / 27.95 (ESS 4,741); g2 (log10 w ≤ 1.3) 26.33;
ρ_ice half-boxes 26.57 / 26.37 (±0.10); ρ_ice derivative unchanged
(−0.32 vs −0.31 km per 10 kg/m³). Optional: a direct g1 run (reweighting
ESS is 4.7k).

All seven caveats of ../injections/README.md carry over, plus a new one:

> 8. The libration elastic factor assumes the cache's elastic shell (F from
> ALMA elastic y1; void if the shell rheology changes). A softer or
> viscoelastic shell would make the inferred zb thinner, by up to ~0.7 km for
> Maxwell and ~2 km for Andrade/soft-base
> (../libration_rheology_bound/README.md; that bound is not yet
> reviewer-accepted). This one-sided model-form term is not in
> libration_sys_frac.

October acceptance artifacts, adjudicated carry-over: injections (noise-free
+ noisy) CARRY OVER (they test sampler fidelity and identifiability geometry,
unchanged within 1%/0.03); (g) variants REPLACED by the exact reweighting
above; (h) C/MR² pushforward carries over exactly (mass closure has no
libration), but its libration pushforward row is frozen-convention — labelled
in ../prior_h/README.md; Fisher-MAP carries over (widths within 1%); branch
comparison carries over as the inequality (note added).

## Different-seed replicate (seed 2027, 2026-11-11): C1 PASS

Run on the Mac with `SEED=2027 REPS=1 plans/cluster/enceladus_posterior_a.sh`
(repo 5806d86b; config confirmed formulation vh2008_repaired,
multiplicative_factor 1.0, sidecar declared). The frozen-convention baseline
was preserved beforehand in `../posterior_a_seed2027_frozenconv/`.
Kish ESS 15,979; log Z −7.1132 ± 0.0068. zb 25.545 / 26.478 / 27.462.

Against (a′) seed 75, per the pre-stated acceptance criteria:

| param | Δmean/SE (≤3) | Δmean/sd (≤0.03) | sd ratio (±2%) | wKS (<0.018) |
|---|---|---|---|---|
| zb_km | −0.62 | −0.007 | 1.000 | 0.0089 |
| log10_wOcean_ppt | +1.80 | +0.020 | 0.986 | 0.0172 |
| rho_ice_kgm3 | +0.74 | +0.008 | 1.005 | 0.0113 |
| rho_rock_kgm3 | +0.15 | +0.002 | 0.998 | 0.0095 |
| compensation_C2 | −0.01 | −0.000 | 0.998 | 0.0098 |
| libration_sys_frac | +1.68 | +0.019 | 1.008 | 0.0112 |

- |Δ median zb| = 0.011 km (criterion ≤ 0.05).
- |Δ log Z| = 0.0004 = 0.05 combined σ (criterion ≤ 3).
- corr(ρ_s, C2): −0.8789 vs −0.8776, diff 0.0013 (criterion ≤ 0.01).
- The one near-threshold value is wKS(log10_w) = 0.0172 against the ~0.018
  critical value; log10_w is the prior-dominated parameter, so its KS tests
  the two runs' prior sampling noise. All criteria pass.

**Every criterion PASSES. C1 is met; the adjudicated headline is in effect.**

## Status
- Same-seed reproducibility, closure, support: **`verified`** (artifacts above).
- Monte-Carlo error under the repaired model: **`verified`** (seed-2027
  replicate table above; every pre-stated criterion met).
- Revised headline: **IN EFFECT (2026-11-12).** zb = 26.5 +1.0/−0.9 km
  (adjudicated wording above) replaces the provisional 26.8 +1.0/−0.9 km.
