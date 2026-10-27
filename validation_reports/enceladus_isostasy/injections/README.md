# Injection recovery and the ρ_ice U[880, 940] sensitivity run (Fisher-checkpoint conditions 4 and 5)

- **Launched:** 2026-10-27, commit e4df7f95, on the Mac production cache. Same sampler settings as run (a) (seed 75, n_eff 8000, one repeat). Each run took 1667–1882 s, four running concurrently.
- **Injection configs:** made by `plans/scripts/enceladus_injection_config.py`, noise-free. Each truth scores log L = 0 under its own data.
- **Report:** `plans/scripts/enceladus_injection_report.py` writes `injection_report.json`.
- **Raw outputs:** `../posterior_{inj_w18,inj_ri918,inj_ri933,ri880}/`.

**Recovery criterion, stated in STATUS before the results:** each truth lies inside the 16–84 % interval for zb, ρ_s and C2, or the reason is understood. w and ρ_ice are prior-dominated, so their posteriors revert toward the prior rather than recovering the truth.

## Injections: all three PASS

Common truth: zb 27 km, ρ_s 2410, C2 0.81, sys_frac 0.
- **inj_w18:** log10 w 1.8, ρ_ice 925. This is the high-salinity truth where the cache hole and the w-dependence interact (run (h) review).
- **inj_ri918 / inj_ri933:** log10 w 0.5, with ρ_ice at the ends of its prior.

Each cell gives the posterior 16 / 50 / 84 %, the truth's percentile, and **IN** or **out** of the 16–84 % interval.

| param | inj_w18 | inj_ri918 | inj_ri933 |
|---|---|---|---|
| **zb (truth 27)** | 25.64 / 26.52 / 27.49 · pct 0.69 · **IN** | 26.16 / 27.10 / 28.11 · 0.46 · **IN** | 26.40 / 27.33 / 28.39 · 0.36 · **IN** |
| **ρ_s (2410)** | 2369 / 2408 / 2452 · 0.52 · **IN** | 2385 / 2424 / 2467 · 0.37 · **IN** | 2371 / 2408 / 2452 · 0.52 · **IN** |
| **C2 (0.81)** | 0.682 / 0.779 / 0.876 · 0.62 · **IN** | 0.708 / 0.805 / 0.899 · 0.52 · **IN** | 0.708 / 0.806 / 0.900 · 0.52 · **IN** |
| log10 w | truth 1.8: −0.51 / 0.54 / 1.54 · 0.93 · out (prior-dominated) | truth 0.5: · 0.50 · IN | truth 0.5: · 0.49 · IN |
| ρ_ice | truth 925: · 0.50 · IN | truth 918: 918.2 / 924.8 / 931.7 · 0.15 · out (prior-dominated) | truth 933: · 0.90 · out (prior-dominated) |
| sys_frac | · 0.50 · IN | · 0.50 · IN | · 0.49 · IN |
| C/MR² | truth 0.3329 · 0.58 · IN | truth 0.3320 · 0.55 · IN | truth 0.3326 · 0.53 · IN |

The best log L in each run is −0.004, −0.001 and −0.002, so the sampler finds the χ² = 0 truth.

**Why the out-of-band cells are understood.** In each case the nuisance's posterior reverts to its prior, as expected, and zb follows it through the known correlations. Conditioning on the truth's own nuisance value brings zb back to the truth:
- **inj_w18:** zb in the bin log10 w ∈ [1.7, 1.9) is 27.07 [26.18, 28.03]. The marginal is pulled to 26.52 because the log-uniform w prior puts 83 % of the mass at log10 w < 1.5, where zb is about 0.6 km thinner. This is the same mechanism run (g) quantified.
- **inj_ri918 / inj_ri933:** zb with |ρ_ice − truth| < 2 is 27.25 and 27.23 against a truth of 27. The marginals shift ∓0.2 km through dzb/dρ_ice ≈ −0.19 km per 10 kg/m³.

So the joint (zb, nuisance) posterior is unbiased. The marginal zb carries a prior-dependent shift bounded by the nuisance prior width, consistent with run (g).

## ρ_ice ~ U[880, 940] (sensitivity only)

This box is outside the declared EOS-nuisance semantics. It was narrowed to [915, 935] on 2026-08-13, and NEW_MODERATE_4, the question of how to justify the ρ_ice bound, is still OPEN.

| quantity | (a) U[915, 935] | U[880, 940] | shift / sd(a) | sd ratio |
|---|---|---|---|---|
| zb | 26.76 | **27.06** | +0.31 | 1.07 |
| ρ_ice | 924.9 | 909.4 (16–84 %: 889–930) | −2.65 | 3.02 |
| ρ_s | 2411 | 2397 | −0.34 | 1.04 |
| C2 | 0.814 | 0.810 | −0.05 | 1.03 |
| log10 w | 0.506 | 0.516 | 0.00 | 1.00 |
| C/MR² | 0.33238 | 0.33242 | +0.00 | 1.00 |
| log Z | −7.096 | −7.094 | | |

- **ρ_ice is unconstrained by the data.** Its posterior is the wider prior.
- **Effect on zb.** zb moves +0.3 km, widening 7 %, through dzb/dρ_ice ≈ −0.19 km per 10 kg/m³. In ρ_ice bins zb runs from 27.46 (880–900) to 26.54 (935–940).
- **Everything else is unaffected,** including C/MR² and the evidence.

So even an ice-density box three times wider moves the zb headline by less than a third of its 1σ.

This run is a nuisance-width sensitivity, not a test of a porosity interpretation of ρ_ice. A porous shell would also have to feed the hydrosphere mass balance, and ri880 changes ρ_ice only in the isostasy and libration paths (the declared EOS-nuisance semantics).

## Ice-density prior variants (condition 5), complete

| ρ_ice prior | zb median [16–84 %] | method |
|---|---|---|
| U[915, 935] (run a) | 26.76 [25.87, 27.74] | direct |
| U[915, 925] | 26.84 [25.94, 27.80] | exact reweighting of (a), ESS 8.1k, both seeds agree |
| U[925, 935] | 26.69 [25.80, 27.66] | exact reweighting of (a), ESS 8.1k |
| U[880, 940] | 27.06 [26.10, 28.09] | direct |

## Status

- Injection recovery (condition 4): **`verified` PASS** against the pre-stated criterion.
- ρ_ice prior sensitivity (condition 5): **`verified`**, and negligible at ≤ 0.3 km.
- **Ocean-branch headline (runs a/g/h, injections, ρ_ice variants, Fisher-MAP): scientifically ACCEPTED, reviewer 2026-10-27** (see the final review section below). The frozen-branch evidence comparison is a separate later step.

## Final acceptance review (Opus, 2026-10-27): ACCEPT WITH EDITS. The ocean branch is scientifically ACCEPTED (stands after the noisy-injection re-adjudication).

The reviewer independently reproduced every number in this README, the run (a)/(g)/(h) records, the half-box reweighting, the Fisher ratios and the g1 prior implementation.

**Applied edits:**
- The libration provenance text in both Enceladus configs now points to `papers/park2024global_ERRATUM.md`. The config hash is unchanged (metadata only).
- Added the "not a porosity test" sentence (ri880 section above).
- Corrected the run (a) wording threshold from 10 to 3 ppt.

**Calibration (reviewer re-adjudication 2026-10-27).** Noise-free injections test bias only.

A posterior-local importance-reweighted SBC was attempted and **REJECTED as a method**:
- Without subtracting L(x_obs), it is tautological. That was the reviewer's original script, and its KS ≤ 0.026 carries no information.
- With the subtraction (`plans/scripts/enceladus_posterior_calibration.py` → `calibration_report.json`), its truth distribution, posterior (a), does not match the inference prior.
- Posterior (a) is also an inadequate proposal for replicate data about √2σ away. A linear-Gaussian control fails worse, with ESS down to 13.
- Its output (rank KS 0.055–0.081, 68 % coverage 0.71) is therefore **not evidence either way**. It is kept, labelled as a rejected method, so the record shows why the approach was dropped.

**Sampler fidelity rests instead on:**
- the different-seed replicate (≤ 0.015 sd);
- Fisher at the MAP (widths within 0–8 %);
- agreement between direct and reweighted g1/g2 (≤ 0.02 km);
- ri880 restricted to [915, 935] reproducing run (a) within 0.01 km (reviewer check);
- noise-free recovery with max log L ≈ 0.

**Ruling: ocean-branch acceptance does not wait on calibration.** Global SBC is deferred to the SBI stage.

As a robustness RECORD (not a gate), three real noisy-injection pocoMC runs are under way: inj_noise1–3, truth zb 27, log10 w 0.5, ρ_ice 925, ρ_s 2410, C2 0.81, noise seeds 1–3, with truth log L −1.23, −3.30 and −10.10.
- Gate for the record: zb, ρ_s and C2 truths inside the 95 % interval in all three; −½χ²(truth) ≤ max log L; ESS > 10k.
- The reviewer reopens acceptance if a gated truth falls outside its 95 % interval.

**Noisy-injection result (2026-10-27, `noisy_injections.json`): the reopen condition tripped. Reviewer RULING: expected exact-posterior behaviour on extreme draws, not a sampler defect, so acceptance STANDS.**

| run | noise / σ_infl (C20, C22, C30, lib) | whitened χ²₄ (p) | ESS | max log L vs truth log L | zb truth pct | ρ_s truth pct | C2 truth pct |
|---|---|---|---|---|---|---|---|
| noise1 | +0.34, +0.82, +0.33, −1.30 | 2.46 (0.65) | 16,641 | −0.24 ≥ −1.23 ✓ | 0.069 IN | 0.607 IN | 0.341 IN |
| noise2 | +0.19, −0.52, −0.41, **−2.43** | 6.60 (0.16) | 12,767 | −0.14 ≥ −3.30 ✓ | **0.005 OUT** | 0.369 IN | 0.342 IN |
| noise3 | **+2.02, −2.54**, +0.41, −0.57 | **20.19 (0.0005)** | 16,030 | −0.01 ≥ −10.10 ✓ | 0.270 IN | **0.003 OUT** | 0.742 IN |

**Reviewer-adjudicated reading:**
- **Spec deviation.** The noise was drawn from the raw σ without correlation: < 1.1 % in σ, plus the missing C20/C22 correlation. `enceladus_injection_config.py` is now corrected to draw from the inflated, correlated covariance (empirical check: sd within 1 %, correlation 0.468 against a target of 0.47). No rerun is required.
- **The excursions.** noise2 (libration −2.43σ) and noise3 (χ²₄ = 20.2, p = 5e-4) put a gated truth at the 0.005 and 0.003 percentiles.
- **Where the sampler can be checked, it agrees.** Importance-reweighting run (a) to the same data reproduces the direct posterior where it has coverage. In noise1 the gated percentiles agree within 0.01, and the quantiles within 0.15 km and 1 kg m⁻³.
- **Where it can't, the direct runs are the better posteriors.** IS ESS is 1660 and 466 in noise2 and noise3. The direct runs sit 0.5 and 0.9 sd further toward the new data, with wider upper tails. noise3's max log L is −0.71 under IS against −0.013 direct. The IS estimate is therefore a lower bound on the shift and cannot adjudicate those two.
- **The excursions are driven by the data draws.** Under both methods the out-of-band truths sit at percentile ≤ 0.02.
- **Scope.** Three runs test sampler robustness, not coverage. A global SBC is deferred to SBI.
- **Linear prediction.** Percentiles more extreme than the linear prediction are not a concern. The linear prediction and IS share run (a)'s local geometry, and the extra extremity goes toward the new data with wider tails, the opposite of overconfidence.
- **Optional, non-gating: different-seed pocoMC replicate of noise3, DONE** (`posterior_inj_noise3s2027_seed2027`, same injected data, sampler seed 2027). It confirms the extreme-tail posterior directly:
  - ρ_s truth percentile 0.0032 against 0.0029; zb 0.272 against 0.270; C2 0.736 against 0.742;
  - every parameter's mean agrees to ≤ 0.02 posterior sd (≤ 1.8 MC-SE), and the 95 % bounds agree to ≤ 1 kg m⁻³ and 0.02 km;
  - ESS 16,011 against 16,030; max log L −0.012 against −0.013; log Z −6.386 against −6.388.

**Required caveats on the headline:**
1. uniform-rock assumption;
2. 2200 ρ_s floor, open but not binding;
3. ρ_ice EOS-nuisance semantics (NEW_MODERATE_4 open; a caveat, not a blocker);
4. ocean branch only;
5. H&M isostasy and shape σ_model are a model choice;
6. libration 0.092 ± 0.003, with the σ cited to Park 2024's original uncertainty statement because the erratum text omits it;
7. the cache hole, negligible for the posterior.

**Adjudicated headline wording:**

> Under the ocean-branch model (uniform-density rock interior, ρ_s ~ U[2200, 2600] kg m⁻³; Airy-isostatic Hemingway & Mittal (2019) gravity with Tajeddine et al. (2017) shape; Δρ-consistent Eq.-12 libration), conditioning on Cassini C20, C22 and C30 and the Park et al. (2024, erratum) physical libration of 0.092° ± 0.003° gives an ice-shell thickness zb = 26.8 +1.0/−0.9 km (median, 16–84 %) for a log-uniform ocean salinity prior (log10 w ~ U[−1, 2], w in ppt), and 27.2 ± 1.0 km if salinity is uniform in w over 0.1–100 ppt. zb is constrained jointly with salinity, which the data do not constrain: conditional on w, zb is 26.6 km for w ≲ 3 ppt and 27.3 km for w ≳ 30 ppt. The prior dependence (≤ 0.4 km) is below the 1σ width. Rock density and gravitational compensation are degenerate (corr = −0.88): ρ_s = 2411 +42/−38 kg m⁻³ and C2 = 0.81 ± 0.09 individually, and the ρ_s 95 % interval [2345, 2501] kg m⁻³ lies well above the prior floor. The derived C/MR² is 0.3324 ± 0.0017; under the uniform-rock model this constraint is carried by the rock density (corr −0.98, ~96 % of its variance) and is a model-dependent derived quantity, not an independent moment-of-inertia measurement (see ../prior_h/README.md). Ice density is not constrained: widening its prior from U[915, 935] to U[880, 940] kg m⁻³ moves zb by +0.3 km (0.3σ).
