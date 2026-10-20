# R1 / R2 preregistration: NH3 corrected-posterior crosscheck (G1), C10 re-run, C11 gravity no-regression

Status: **PREREGISTERED (2026-10-20).** This file and r1_r2_prereg_bands.json
are committed BEFORE any R1/R2 statistic is computed.

Reviewed by the scientific-reviewer (Opus 5.5), verdict APPROVE WITH EDITS.
All five edits are applied below.
- **Critical edit.** The stock shape floor is mis-calibrated when the
  corrected set is much smaller than the reference ESS. Reference self-draws
  of size L fail it on 4–8 of 8 replicates for every parameter. The floor is
  therefore evaluated at matched n.
- **Major edit.** The C9 band now includes the corrected side's own
  median sampling error.

Every calibration number here comes from reference-only computation.

Governing rulings: MACHINE-B-HANDOFF §0.32, RATIFICATION.md (2026-10-20 section).
Runs on Machine A, which holds the n_eff=8000 seed pkls (user directive 2026-10-20: "we can run everything on this machine").

## Inputs (fixed)

- **Corrected set.** Use the shipped fiducial bundle
  `PlanetProfile/Inference/sbi_artifacts/titan_freegrav_nh3_is_fiducial.npz`
  (sha256 `0818b4e5…`). It holds N=20000 draws at seed 72, and their
  log_w is bit-exact with the committed `log_weights.npy`, which
  `fc5f29f5` verified. The corrected sample set is
  `samples[resampled_indices]`, using the committed
  `resampled_indices.npy` (L=1270, systematic resampling at seed 72). It
  is used unweighted, exactly as validate_sbi expects. No new draws are
  made.
- **Reference.** The n_eff=8000 pocoMC seeds 72/172/272:
  `validation_reports/nh3_diagnosis/matched_reference_neff8000/nh3_reference_seed{72,172,272}_neff8000.pkl`.
  Pool them under treatment A, the convention the §0.21 C16
  re-ratification used: concatenate the three seeds, renormalize each
  seed's weights to w/Σw/3, and make the pooled weights sum to 1.
  Assert before use that config_hash `1611b65f…` is identical across all
  three seeds (`config.generate_hash()`; already confirmed at prereg
  time). The seed pkls do not carry the structure-cache sha, so take it
  from `matched_reference_neff8000/matched_reference_report.json`
  (`3d837cf8…`) and assert it equals the bundle manifest's cache sha.
- **Seeds (all stochastic steps):**
  - crosscheck RNG seed 72;
  - C9 SE_corr bootstrap RandomState(72), B=1000;
  - null-fail calibration RandomState(1072), 10 replicates;
  - R2 plumbing subsamples seed 72 (200 reference draws and 200 corrected
    draws, uniform without replacement).
- **Run-time provenance asserts (blocking).** The driver must check, before
  computing anything:
  - the bundle npz sha256 equals the manifest's `npz_sha256` (`0818b4e5…`);
  - the bundle `log_w` is bit-exact with the committed `log_weights.npy`
    (identical −inf mask, finite entries equal);
  - `systematic_resample(normalized w, 1270, seed=72)` reproduces the
    committed `resampled_indices.npy` exactly.

  The bundle stands in for regenerating the conditioning (§0.32). These
  asserts are what make that substitution valid.

## R1a: G1 crosscheck (13 parameters)

Run `validate_sbi._run_crosscheck(param_names, ref_samples, ref_weights,
corrected_samples, seed=72, matched_n_floor=True)` with the ratified gates:

| Gate | Criterion |
|---|---|
| mean | max(0.25σ, MC error) |
| σ-ratio | within ratified bounds |
| shape-D | excess within floor factor × p99 self-D |
| median materiality | ratified rule |

- `empirical_floor=None`: no wander floor is measured for this campaign.
- **Shape floor evaluated at matched n:**
  `_self_d_floor(..., ess=min(ESS_ref, L))`. Factor 1.5, p99 and 200
  bootstraps are unchanged. This is the opt-in `matched_n_floor=True`;
  the stock default is untouched for every other campaign.
  - Rationale: the stock floor assumes the KS compares at n_ks = ESS_ref.
    Here n_ks = L = 1270, far below ESS_ref (≈39,100, value in the bands
    JSON).
  - Reference-only null (bands JSON `null_fail_table`): the stock floor
    fails self-draws, and the matched floor does not.
- No other threshold is changed.
- Report a per-parameter verdict for all 13 parameters by name.
- The 13 parameters split into:
  - the 7 C9 nuisances below;
  - the **6 non-C9 parameters**: Tb_K, log10_wOcean_ppt, R_core_km,
    rho_core_kgm3, dC20_nh, dC22_nh. The last two are the gravity
    non-hydrostatic offsets, not structure, but are gated identically.

**C9 decision rule for unidentified nuisances.** It applies ONLY to
`c9_nuisance_params` = {log10_eta_Ih, log10_eta_III, log10_eta_V,
log10_eta_VI, log10_eta_sil, log10_zeta, alpha}, and ONLY when the
stock crosscheck FAILs for that parameter.
- s_seed is the sample std (ddof=1) of the three n_eff=8000 per-seed
  weighted medians.
- m_ref is the pooled weighted median.
- m_corr is the median of the corrected set.
- SE_corr is the std of the median of L=1270 weight-resampled
  pooled-reference draws (B=1000, RandomState(72)). It is the corrected
  side's own median sampling error at this L.
- The parameter is **PASS-by-C9** iff BOTH of these hold:
  - the stock σ-ratio gate passes. An unidentified nuisance with the
    wrong WIDTH is a real defect; C9 excuses only the mean and shape
    gates.
  - |m_corr − m_ref| ≤ tol = 2·sqrt(s_seed² + SE_corr²).
- Otherwise the parameter FAILs.
- s_seed, SE_corr, tol, m_ref and the null pass rate per parameter are
  computed from the reference pkls ONLY. They are committed in
  `r1_r2_prereg_bands.json` in the same commit as this file.
- The n_eff=2000 bands (13646f18) are SUPERSEDED because of the
  resolution mismatch: methods §3 forbids mixing resolutions. They are
  reported for information only.

**Lift mapping (per the 2026-10-20 ruling).**
- Re k₂, |Im k₂|: not crosscheck parameters. They depend on R1b and R2
  below.
- ζ, the five η's and alpha: each comes out of quarantine at the
  fiducial only on a stock PASS or a PASS-by-C9.
  - alpha (the rheology exponent) is on the frozen C9 list and is treated
    exactly like ζ and the η's, including quarantine by name on FAIL.
  - A PASS-by-C9 parameter is displayed as "unidentified nuisance:
    agrees with the reference within reference seed scatter" and never
    as constrained.
  - A FAIL parameter stays quarantined by name.
- Any FAIL among the 6 non-C9 parameters means the lift does NOT take
  effect: re-escalate.

## R1b: C10 re-run against the n_eff=8000 pool

- Compute the weighted two-sample KS between the bundle's
  importance-weighted Re k₂ and |Im k₂| pushforward (full N, weights =
  normalized exp(log_w)) and the pooled reference pushforward, as in
  the original C10 driver.
- Gate: KS ≤ 0.15 for each.
- Also report the median-to-median gaps in σ_obs (σ_Re = 0.048,
  σ_Im = 0.035). These stay held to the G2 ≤ 0.5 σ_obs gate.

## R2: C11 gravity no-regression

- **Reference side.** The n_eff=8000 seed pkls carry per-sample
  `c20_results` and `c22_results`, 100% finite; this was checked at
  prereg time. The "no C20/C22 columns" limitation applied only to the
  n_eff=2000 pooled pkl. Use the stored arrays with the treatment-A
  pooled weights.
- **Plumbing cross-check (blocking).** Recompute C20 and C22 through the
  config's `gravity_forward_model` on the structure cache for 200
  reference samples (seed 72, uniform without replacement). Require a
  max relative difference ≤ 1e-9 against the stored values. A mismatch
  means the two sides are not on the same forward path: FAIL R2.
- **Corrected-side plumbing check (blocking).** Do the same for 200
  corrected bundle draws: recompute at HEAD against the bundle's stored
  `c20_results`/`c22_results`, max relative difference ≤ 1e-9. It
  compares the corrected side with itself, so it does not break
  preregistration.
  - Both sides use `_derive_gravity_pair` (the SBI path calls it through
    `_get_mcmc_runner()`), but at different SHAs: the reference at
    0fe68b39, the bundle later.
- **Driver note.** `is_correction_validate.py` lines ~239–251 read
  C20/C22 from `k2[:, 2:4]`. Those columns don't exist, so that path
  silently skipped C11. The R2 driver must read
  `c20_results`/`c22_results` directly.
- Weighted medians are taken over the finite-`log_w` subset. All
  non-finite c20 values coincide with non-finite `log_w`; this was
  checked at review.
- **Corrected side.** Use the bundle's `c20_results` and `c22_results`
  with the IS weights. Assert both exist and are finite on the
  non-rejected draws. **If either side lacks C20/C22, FAIL loudly; never
  skip.**
- **Gate.** For C20 and C22 each: |weighted median_corr − weighted
  median_ref| ≤ 0.1 σ_obs, where σ_obs(C20) = 4.71e-7 and
  σ_obs(C22) = 6.2e-8 (config sigmas; assert they match the config).
  Also report the weighted KS for information only, not gated.
- **Scope.** The corrected-median SE is about 0.035 σ_obs, so the gate
  sits about 2.8 SE out and the null passes about 99%. Because
  dC20_nh/dC22_nh absorb structure, R2 is a gravity no-regression check
  ONLY. It is not a test of the structure sector.

## Outcome rule

| Outcome | Consequence |
|---|---|
| R1a structure sector PASS, R1b PASS, R2 PASS | The fiducial-scope lift takes effect per the 2026-10-20 ruling. ζ and the η's are listed by their R1a/C9 verdicts. |
| Any FAIL | Recorded, no lift, re-escalated to the scientific-reviewer. No threshold or rule is revisited after results are read. |

Outputs: `r1_r2_report.json` and `r1_r2_RESULT.md` in this directory.
