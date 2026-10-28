# Enceladus ocean-branch SBI: design and preregistration (DRAFT, review pending)

Status: **PREREGISTERED (2026-10-27).** The scientific-reviewer (Opus 5.5) returned APPROVE WITH EDITS, and every edit is applied below. The joint-gate k and floor recipe are fixed from reference-only calibration (`joint_gate_calibration.json`, `plans/scripts/enceladus_sbi_joint_gate_calibration.py`), and the required permutation-mutation sanity check PASSES. This file is committed BEFORE any dataset is generated or any flow trained. Changing N or capacity after any gate result is forbidden.

## Scope

An amortized NPE posterior for the ACCEPTED ocean-branch model: config `enceladus_cassini_isostasy_ocean6D_rhorock.json`, the Mac production cache, and the 6 sampled parameters.
- Observable vector x = (C20, C22, C30, libration_deg), in the config's order.
- The reference posterior is run (a) (`posterior_a/`). Seeds 75 and 2027 are both accepted.

## Declared training prior: the EFFECTIVE prior (run (h) ruling)

- **Definition.** The nominal prior box, restricted to the servable (zb, w) region: cells whose four corners are all built.
- **Implementation.** `MCMCRunner.generate_sbi_dataset(drop_nonfinite=True)` drops draws whose cell is a hole. Holes are counted and discarded, not labelled as rejects.
- **Pilot check.** On the production cache, 364 of 400 draws were kept (91%, consistent with run (h)'s 0.893). The first 40 kept rows equal the likelihood path's `_derive_gravity_isostatic` / `_derive_libration_deg` to rtol 1e-12 (0 mismatches).
- **Recorded with the artifact:**
  - the servable fraction;
  - the hole (zb, w) region;
  - the statement that the flow's support is the effective prior.
- **Inside the hole.** A user-supplied x is conditioned on the effective prior, which is out of support inside the hole. This is declared in the GUI scope note.

## Dataset

- **Size.** N = 1,000,000 kept rows, over-requested about 1.12× to absorb hole drops.
- **Noise.** obs_noise=True: noise is added at simulation time from the likelihood's inflated, correlated covariance (σ ⊕ sigma_model_add, ρ(C20, C22) = 0.47), matching the NPE convention used for Titan.
  - **Verified and FIXED (2026-10-27).** `generate_sbi_dataset(obs_noise=True)` applied the correlation but used the RAW σ, not σ ⊕ sigma_model_add, so it did not match the likelihood.
    - Fixed in `mcmc_runner.py`. Configs without sigma_model_add, which is every non-Enceladus campaign, are byte-identical.
    - Test `tests/sbi_noise_sigma_model_test.py` drives the real generator through theta_override. It passes, and it fails on the old code (mutation-checked).
    - Suites: sbi/enceladus/noise/isostasy/libration, 251 passed, 1 xfail.
- **Seeds.** Shard k (k = 0…10) uses data seed 75 + k and noise seed 7575 + k. The G3 and G5 held-out pairs use data seeds 900 + k and noise seeds 9900 + k, disjoint from training.
- **Generation.** Single-threaded per process; 11 parallel shards with disjoint seed streams (seed 75 + k), merged in order; about 10 minutes wall.
- **Provenance.** Config hash, cache sha256, and per-shard hole counts.

## Training

- **Estimator.** nsf (`SBIRunner.train(density_estimator='nsf')`) at the deployed defaults: hidden 50, 5 transforms, 10 bins, z-score independent.
- **Seed.** Train seed 75.
- **Process isolation.** Training runs in a separate process from generation (the libomp hazard).
- **No tuning against gate results.** If a gate fails, the failure is recorded and escalated, never re-run with changed capacity to pass.

## Gates (preregistered; thresholds fixed here before any result)

- **G0. Data integrity (before training).**
  - The servable fraction of the 1M set is in [0.887, 0.900].
  - The empirical noise sd is within 1% of σ_infl, and corr(C20, C22) is in [0.46, 0.48]. Measure these on theta_override repeats of the same θ, as in `tests/sbi_noise_sigma_model_test.py`.
  - 1,000 random rows match the likelihood path's noise-free predictions to rtol 1e-10 (noise-free regeneration through theta_override).
  - Shard k used noise_seed = 7575 + k.

- **G1. Crosscheck at x_obs.**
  - Compare the flow posterior against the pooled run (a) reference (seeds 75 and 2027, ESS 32,273) on all 6 parameters. Use N_post = 40,000, so that n_ks = floor(ESS_ref) = 32,273 and the stock floor is built at the n actually compared.
  - Use `_run_crosscheck(..., matched_n_floor=False)` with the stock ratified gates, and record n_ks.
  - No C9 or other fallback applies. w, ρ_ice and sys_frac are prior-dominated (posterior/prior sd 1.001, 0.999 and the prior sd), so a correct flow reproduces their prior, and a failure is a defect to escalate.
- **G2. Derived pushforward at x_obs.**
  - Flow-sample C/MR² (through the runner's mass closure) and the 4 predicted observables, against run (a)'s.
  - Pass if the median-to-median gap is ≤ 0.25 posterior sd for zb, ρ_s, C2 and C/MR², AND the weighted KS on C/MR² and on each observable pushforward is ≤ 0.10.
  - Both sides are pushed forward noise-free through the same `theta_override` path. Run (a) has no stored C20/C22 columns, so these are recomputed, and the gate fails loudly if any column is absent.
  - Additional pass conditions:
    - (i) each observable's pushforward median gap is ≤ 0.1 σ_infl;
    - (ii) **joint gate (first class):** |corr_flow(ρ_s, C2) − corr_ref(ρ_s, C2)| ≤ 0.05, with corr_ref = −0.8764, AND the energy distance between flow and reference in reference-whitened (zb, ρ_s, C2) space is ≤ k × p99, with k = 2 and p99 = 0.006162. The p99 is the reference self-distance floor: 100 independent reference-vs-reference resamples at matched n (40,000 vs 32,272), energy distance on m = 2,000 subsamples per side, RNG 20261027. Gate value: 0.01232.
      - Null false-fail rate: 0/100 on both statistics.
      - **Mutation sanity (required, PASSES):** a reference sample with ρ_s and C2 independently permuted in physical space fails 20/20 (energy distance 0.28–0.31, Δcorr 0.86). A 20 % partial permutation fails 20/20 jointly (Δcorr ≥ 0.17; energy distance alone catches 75 %).
    - (iii) |P_flow(C2 > 0.95) − P_ref(C2 > 0.95)| ≤ 0.02, with a reference value of 0.062. Report the sbi out-of-box leakage fraction.
- **G3. Global SBC.**
  - Use n = 2,000 held-out (θ₀, x₀) pairs from the effective prior, with seeds disjoint from training (data 900 + k, noise 9900 + k).
  - Rank θ₀ among 1,000 flow draws per pair.
  - Pass if, per parameter, the BH-FDR-adjusted KS p is ≥ 0.05 over the 6 parameters, AND the 68 % central coverage is in [0.65, 0.71], which is ±2.9 binomial σ.
  - Also pass the joint rank statistic: the rank of log q(θ₀|x₀) among log q(θᵢ|x₀), i.e. HPD coverage, must be KS-uniform at p ≥ 0.05.
  - Also report, without gating, the same statistics on the zb₀ < 15 km stratum.
  - SBC is prior-averaged and does not certify local calibration at x_obs. That is the job of G1, G2 and G4.
- **G4. Noisy-injection recovery.**
  - Condition the flow on the 3 noisy-injection data sets plus the 3 noise-free injection data sets that already have pocoMC posteriors.
  - Pass if, for zb, ρ_s and C2 on each of the 6 data sets, the flow median is within 0.3 pocoMC-sd of the pocoMC median AND the sd ratio is in [0.7, 1.4].
  - The noise3 reference is the pooled pocoMC seeds 75 + 2027.
  - pocoMC's L = 0 in the hole makes its target identical to the effective-prior posterior for every x. No adjustment is made.
  - This is a direct flow-versus-sampler comparison on off-nominal data, and it removes the IS-coverage problem.
- **G5. Hole leakage.**
  - Measure on the G3 pairs with zb₀ < 15 km and on a dedicated 500-pair thin-shell set (the effective prior restricted to zb < 15 km).
  - The flow mass in hole cells, averaged over pairs, must be ≤ 1 %, and the 95th percentile per pair must be ≤ 5 %.
  - At x_obs the leakage is reported only. The gate is vacuous there because the reference's minimum zb is 23.6 km.

## Query-time guard (required for GUI and amortized use)

GUI and amortized queries draw 10k flow samples. If more than 1 % of them fall at zb < 15 km, or in a hole cell, the result is labelled "effective-prior posterior; servable cells only — the cache hole is excluded and the physical posterior there is unavailable to either SBI or pocoMC". No reachability input is added to the flow, because reachability is a property of θ, not of x.

**Outcome.**
- All of G1–G5 pass: deploy as split-status-free "SBI, validated against pocoMC at x_obs and globally calibrated (SBC)".
- Any FAIL: record it, do not deploy as validated, and escalate. The pocoMC run (a) stays authoritative.
- Changing N or capacity after any gate result is forbidden.
- Optional, non-gating diagnostic: train a 300k flow and report its G1 and G2 statistics against the 1M flow, to check for undertraining. It is never a substitute for the 1M artifact.

## Reviewer rulings on the open questions (2026-10-27)

1. **N = 1M is kept.** About 0.9 % of rows land in the posterior's identified volume (posterior/prior width 0.082 for zb, 0.35 for ρ_s, 0.32 for C2), which is about 9,000 rows at 1M against about 2,700 at 300k. NH3 failed by under-updating in an informative region.
2. **Thresholds.**
   - G2: 0.25 sd and KS 0.10 are consistent with each other; the σ_obs and joint gates are added.
   - G3: n = 2,000 and the window [0.65, 0.71].
   - G4: an sd-ratio condition is added.
   - G5: redefined as above.
3. **No C9 fallback.** The ρ_s–C2 ridge and the box edge are handled by the joint gate and the edge-mass check, not by G1.
4. **No reachability input.** The query-time guard above is required instead.

**Target equivalence.** pocoMC run (a) samples π_nom · L · 1_servable / Z, which equals π_eff · L / Z_eff for every x. So the flow's effective-prior target and pocoMC's target are the same posterior, and they differ only in the normalising constant (0.113 nat). Neither method can see the physical posterior inside the hole. That limitation is the scope of the query guard, not a gate.
