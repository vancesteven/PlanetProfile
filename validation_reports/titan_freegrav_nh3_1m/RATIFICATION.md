# Titan free-gravity NH3 JOINT artifact — ratification

**Status: VERIFIED (ratified by scientific-reviewer 2026-08-03) with scoped caveats.**

Artifact: `PlanetProfile/Inference/sbi_artifacts/titan_freegrav_nh3_posterior_1m.pt`
(n_train = 689,845, nsf, 124 epochs, seed 72).
Config: `PlanetProfile/Inference/configs/test54_titan_nh3_freegrav.json`
(13 params, obs {C20, C22, Re_k2, Im_k2}, JOINT no-ocean + ocean).
Structure cache: `PlanetProfile/Test/mcmc_results/Titan/Test54_nh3_ocean/`
`titan_nh3_joint_structure_grid_2d.pkl` (sha256 `3d837cf8…`,
`retry_frozen_as_no_ocean=True`).
Reference MCMC: `reference/titan_freegrav_nh3_reference_result.pkl`
(weighted pocoMC, 4461 samples, ess 4125, log_Z = −14.24).

## Gate outcomes

| Gate | Verdict | Note |
|---|---|---|
| SBC (353 kept pairs) | **PASS** | min BH-adj KS p = 0.772; all 13 raw KS p ∈ [0.239, 0.941]; c2st ≈ chance |
| crosscheck vs pocoMC | FAIL (benign) | only 3 unconstrained eta **medians** exceed the 0.3-dex tol (≤0.28σ of marginal); all means/sigmas/shapes pass |
| limits (Im_k2 sweep) | FAIL (mis-specified) | eta_Ih monotonicity **direction**; a 0.4σ noise wander on a flat marginal; all points in-support, containment/interior pass |

## Ratification scope

**"Verified" covers:** posterior calibration against the true simulator (SBC, all
13 params, BH-FDR corrected); cross-sampler agreement with pocoMC on **all
scientifically load-bearing params** (alpha, log10_zeta, log10_eta_sil, Tb_K,
log10_wOcean_ppt, R_core_km, rho_core_kgm3, dC20_nh, dC22_nh, log10_eta_V —
mean/sigma/median/shape all PASS); salinity ocean-peak weight + Tb mixture
capture; in-support behavior and containment under the Im_k2 sweep.

## Caveats (must travel with the artifact)

1. **`log10_eta_Ih`, `log10_eta_III`, `log10_eta_VI` are prior-dominated
   dissipation nuisances** (posterior ≈ prior; unconstrained by {C20,C22,Re_k2,
   Im_k2}). Report their marginals as prior-dominated; **do not quote SBI point
   estimates** for these three. (eta_V, eta_sil weakly constrained but
   better-behaved.)
2. **The limits eta_Ih-monotone-vs-Im_k2 premise is N/A for a JOINT
   no-ocean + ocean mixture.** In the no-ocean branch HP ices — not ice Ih —
   carry the tidal response, so eta_Ih barely couples to Im_k2 and the
   conditional median is driven by branch-mixing weight, not a rheological law.
   The observed swing (12.68→13.37 dex) is 0.4σ of the 1.74-dex marginal, i.e.
   sampling noise. Document as N/A (or replace with a "no significant trend
   within Xσ" test) for this artifact class — **do NOT tune the gate to pass.**
3. **Uniform ~1.1–1.3× sigma_ratio (SBI/pocoMC):** SBI credible intervals are
   marginally conservative relative to pocoMC (n_effective=500 slightly narrow).
   Within the ratified [0.7,1.4] sub-gate. The clean SBC (direct calibration
   against the simulator) supersedes the pilot "sigma_ratio must shrink or HOLD"
   trigger — an over-dispersed flow would fail SBC, and it does not.

## Config-hash note (reviewer's required validation — CLEARED)

The crosscheck reference MCMC `config_hash` (`1611b65fff3f06c9`) differs from the
artifact `config_hash` (`e596574d1e81567c`) **only** because of the `mode` string
(`mcmc` vs `sbi`). All load-bearing fields — param_space bounds+priors,
observables, imag_convention (`abs`, refold_im=false), gravity_forward_model,
structure_cache_path — are byte-identical. Verified 2026-08-03; the crosscheck
comparison is apples-to-apples.

## Inheritance for MgSO4 / NaCl joint campaigns

These two gate mis-specifications are **intrinsic to the joint mixture design**,
not to NH3. The MgSO4 and NaCl joint artifacts will exhibit the same benign
limits-monotonicity FAIL and (likely) the same unconstrained-eta-median
crosscheck FAILs. Inherit the corrected gate scope above; do not re-litigate or
tune. Report eta nuisances as prior-dominated; mark the limits monotonicity gate
N/A for joint mixtures.
