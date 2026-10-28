# Enceladus ocean-branch SBI: design and preregistration (DRAFT, review pending)

Status: **DRAFT**. Nothing in this file has been run. It becomes PREREGISTERED only after scientific-reviewer approval and a commit, and before any training or gate statistic is computed.

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
- **Seeds.** data = 75, noise = 7575.
- **Generation.** Single-threaded per process; 11 parallel shards with disjoint seed streams (seed 75 + k), merged in order; about 10 minutes wall.
- **Provenance.** Config hash, cache sha256, and per-shard hole counts.

## Training

- **Estimator.** nsf (`SBIRunner.train(density_estimator='nsf')`) at the deployed defaults: hidden 50, 5 transforms, 10 bins, z-score independent.
- **Seed.** Train seed 75.
- **Process isolation.** Training runs in a separate process from generation (the libomp hazard).
- **No tuning against gate results.** If a gate fails, the failure is recorded and escalated, never re-run with changed capacity to pass.

## Gates (preregistered; thresholds fixed here before any result)

- **G1. Crosscheck at x_obs.**
  - Compare the flow posterior (N_post = 20,000) against the pooled run (a) reference (seeds 75 and 2027) on all 6 parameters.
  - Use `validate_sbi._run_crosscheck(..., matched_n_floor=False)` with the stock ratified gates. The flow draws are many, so n_ks = ESS_ref and the stock floor is correct, unlike the NH3 IS case.
  - Use median materiality for log10_w, and the ratified sd-scaled rule for the others.
- **G2. Derived pushforward at x_obs.**
  - Flow-sample C/MR² (through the runner's mass closure) and the 4 predicted observables, against run (a)'s.
  - Pass if the median-to-median gap is ≤ 0.25 posterior sd for zb, ρ_s, C2 and C/MR², AND the weighted KS on C/MR² and on each observable pushforward is ≤ 0.10.
- **G3. Global SBC.**
  - This is the calibration test deferred from the headline review.
  - Use n = 500 (θ₀, x₀) pairs from the EFFECTIVE prior (servable draws, noise added).
  - Compute the rank of θ₀ among 1,000 flow draws per pair, per parameter.
  - Pass if the KS test against uniform has BH-FDR-adjusted p ≥ 0.05 for all 6 parameters, and the 68% coverage is in [0.63, 0.73] for each.
- **G4. Noisy-injection recovery.**
  - Condition the flow on the 3 noisy-injection data sets plus the 3 noise-free injection data sets that already have pocoMC posteriors.
  - Pass if, for zb, ρ_s and C2, the flow's median is within 0.3 pocoMC-sd of the pocoMC median on each.
  - This is a direct flow-versus-sampler comparison on off-nominal data, and it removes the IS-coverage problem.
- **G5. Support.** No flow mass inside the hole region beyond 1%, measured with flow draws at x_obs.

**Outcome.**
- All of G1–G5 pass: deploy as split-status-free "SBI, validated against pocoMC at x_obs and globally calibrated (SBC)".
- Any FAIL: record it, do not deploy as validated, and escalate. The pocoMC run (a) stays authoritative.

## Open questions for the reviewer

1. Is N = 1M justified? Run (a) has 6 parameters, 2 identified directions and 4 prior-dominated ones. Would 300k be adequate, and is the 1M default for Titan transferable?
2. Are the G2, G3 and G4 thresholds right? G3 uses 500 pairs at ~0.05 s each for flow sampling, which is cheap. G4 uses 0.3 sd.
3. The ρ_s–C2 ridge (corr −0.88) and the prior-box edges: should G1 use the C9-style rule for the prior-dominated w and ρ_ice?
4. Anything missing, for example a reachability (hole) indicator as an extra flow input or as a declared out-of-support guard at query time.
