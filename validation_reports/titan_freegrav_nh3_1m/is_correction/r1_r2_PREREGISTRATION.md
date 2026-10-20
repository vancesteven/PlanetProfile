# R1 / R2 preregistration: NH3 corrected-posterior crosscheck (G1), C10 re-run, C11 gravity no-regression

Status: **DRAFT, REVIEW PENDING. NOT YET PREREGISTERED.** Two scientific-reviewer
passes (Opus) on 2026-10-20 died on network errors before returning a verdict.
Don't run R1/R2 until a reviewer approves this file (with any edits applied) and
the approved version is committed with status PREREGISTERED. No R1/R2 statistic
has been computed. The C9 bands in r1_r2_prereg_bands.json come from the
reference seed pkls only.

Open questions for the reviewer (manager's, not yet answered):
- Is 2·s_seed from 3 seeds too tight? For example, eta_VI has s_seed 0.014 at
  n_eff=8000 against a half-range of 0.116 at n_eff=2000. Options are a
  t-quantile, a floor, or the originally frozen n_eff=2000 bands (commit 13646f18).
- Is L=1270 enough power for the stock crosscheck gates?
- Does the corrected side's c20/c22 use the same forward path as the reference's
  stored arrays?
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
- **Seeds.** Crosscheck RNG seed 72. Nothing else is stochastic.

## R1a: G1 crosscheck (13 parameters)

Run `validate_sbi._run_crosscheck(param_names, ref_samples, ref_weights,
corrected_samples, seed=72)` with the ratified stock gates:

| Gate | Criterion |
|---|---|
| mean | max(0.25σ, MC error) |
| σ-ratio | within ratified bounds |
| shape-D | excess within floor factor × p99 self-D |
| median materiality | ratified rule |

- `empirical_floor=None`: no wander floor is measured for this campaign.
- No threshold is changed.
- Report a per-parameter verdict, with the 6 rheology parameters and ζ
  given by name.

**C9 decision rule for unidentified nuisances.** It applies ONLY to
`c9_nuisance_params` = {log10_eta_Ih, log10_eta_III, log10_eta_V,
log10_eta_VI, log10_eta_sil, log10_zeta, alpha}, and ONLY when the
stock crosscheck FAILs for that parameter.
- s_seed is the sample std (ddof=1) of the three n_eff=8000 per-seed
  weighted medians.
- m_ref is the pooled weighted median.
- m_corr is the median of the corrected set.
- The parameter is **PASS-by-C9** iff |m_corr − m_ref| ≤ 2·s_seed.
  - This is the same 2σ-of-reference-seed-scatter scale as the C16
    bound.
  - Rationale: a nuisance the reference itself cannot pin down shows
    seed-to-seed median wander. Disagreement inside twice that wander
    cannot be attributed to the corrected pipeline.
- Otherwise the parameter FAILs.
- The s_seed values are computed from the reference pkls only and
  committed in `r1_r2_prereg_bands.json` (n_eff=8000 seeds; values in that file) in the same commit as this
  file. No corrected quantity is read to produce them.

**Lift mapping (per the 2026-10-20 ruling).**
- Re k₂, |Im k₂|: not crosscheck parameters. They depend on R1b and R2
  below.
- ζ and the η's (layer viscosities; `alpha` is the rheology exponent):
  each comes out of quarantine at the fiducial only on a stock PASS or a
  PASS-by-C9.
  - A PASS-by-C9 parameter is displayed as "unidentified nuisance:
    agrees with the reference within reference seed scatter" and never
    as constrained.
  - A FAIL parameter stays quarantined by name.
- Any FAIL among the 6 structure-sector parameters (non-C9) means the
  lift does NOT take effect: re-escalate.

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
- **Corrected side.** Use the bundle's `c20_results` and `c22_results`
  with the IS weights. Assert both exist and are finite on the
  non-rejected draws. **If either side lacks C20/C22, FAIL loudly; never
  skip.**
- **Gate.** For C20 and C22 each: |weighted median_corr − weighted
  median_ref| ≤ 0.1 σ_obs, where σ_obs(C20) = 4.71e-7 and
  σ_obs(C22) = 6.2e-8 (config sigmas; assert they match the config).
  Also report the weighted KS for information only, not gated.

## Outcome rule

| Outcome | Consequence |
|---|---|
| R1a structure sector PASS, R1b PASS, R2 PASS | The fiducial-scope lift takes effect per the 2026-10-20 ruling. ζ and the η's are listed by their R1a/C9 verdicts. |
| Any FAIL | Recorded, no lift, re-escalated to the scientific-reviewer. No threshold or rule is revisited after results are read. |

Outputs: `r1_r2_report.json` and `r1_r2_RESULT.md` in this directory.
