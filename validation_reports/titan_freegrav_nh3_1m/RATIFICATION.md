# Titan free-gravity NH3 JOINT artifact — ratification

**Status: SPLIT (amended by scientific-reviewer 2026-08-03 after the
posterior-predictive check).**
- **VERIFIED for the gravity/structure sector** — C20, C22, R_core, rho_core,
  salinity (log10_wOcean_ppt), Tb, and the null non-hydrostatic offsets
  (dC20_nh, dC22_nh). These load-bearing inferences hold.
- **NOT VERIFIED for the tidal (k2) / dissipation sector** — Re_k2, Im_k2,
  log10_zeta, and the eta viscosities. The SBI flow under-updates this sector at
  the Titan datum (see caveat 1); the MCMC reference is authoritative for any
  k2- or dissipation-derived quantity.

The earlier "VERIFIED with benign eta-nuisance caveats" framing was too narrow
and is superseded — see the PPC finding below.

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
| crosscheck vs pocoMC | FAIL | 3 eta medians exceed the 0.3-dex tol on the marginals; **BUT** the per-param gate is blind to the joint tidal-sector under-update surfaced by the PPC (see below) — do not read this FAIL as "benign nuisances only" |
| limits (Im_k2 sweep) | FAIL (mis-specified) | eta_Ih monotonicity **direction**; a 0.4σ noise wander on a flat marginal; all points in-support, containment/interior pass |

## Ratification scope

**"Verified" covers:** posterior calibration against the true simulator (SBC, all
13 params, BH-FDR corrected); cross-sampler agreement with pocoMC on **all
scientifically load-bearing params** (alpha, log10_zeta, log10_eta_sil, Tb_K,
log10_wOcean_ppt, R_core_km, rho_core_kgm3, dC20_nh, dC22_nh, log10_eta_V —
mean/sigma/median/shape all PASS); salinity ocean-peak weight + Tb mixture
capture; in-support behavior and containment under the Im_k2 sweep.

## Posterior-predictive finding (2026-08-03) — tidal sector under-update

A posterior-predictive check (`ppc/ppc_interior_report.json`, built via the
`theta_override` pushforward validated exact against the MCMC reference) shows
the SBI flow **under-updates the entire tidal (k2) sector at the Titan datum**,
not merely the eta "nuisances":

| channel | observed | SBI posterior-pred. median | MCMC posterior-pred. median |
|---|---|---|---|
| Re_k2 | 0.608 ± 0.048 | 0.541 (1.74σ low) | 0.581 (0.6σ low) |
| Im_k2 | 0.135 ± 0.035 | 0.042 (2.9σ low) | 0.093 (1.2σ low) |

Both SBI k2 pushforwards sit essentially at the **prior-predictive median**
(Re_k2 0.526, Im_k2 0.037) — i.e. the flow behaves as if the two k2
observations carry almost no information — while the MCMC posterior-predictive
updates toward the data. Decomposition: of the 2.9σ Im_k2 gap, ~1.2σ is real
data–model tension shared with the MCMC (benign; physics+prior struggle to reach
the observed dissipation — the datum sits at the 86th prior-predictive
percentile), and ~1.7σ is the SBI-vs-MCMC flow deficiency. Mechanism: the three
dissipation-controlling params all shift toward *less* dissipation
(log10_zeta +0.28 dex, eta_Ih +0.31, eta_III +0.21), each passing the per-param
crosscheck by 90–92% of tolerance, compounding multiplicatively in the nonlinear
k2 pushforward. The per-param marginal crosscheck **structurally cannot see this**
(it never pushes parameters forward to an observable).

**SBC-vs-PPC reconciliation (no paradox):** SBC is a *prior-averaged* rank test,
dominated by draws near the prior-predictive median where "sit near the prior"
is correct. It has few draws in the informative upper-Im_k2 tail where the Titan
datum lives, so a localized under-update there washes out. The flow is globally
calibrated (SBC PASS) yet **locally miscalibrated at the Titan datum
specifically** — both are true simultaneously.

## Caveats (must travel with the artifact)

1. **The flow under-updates the tidal-dissipation sector.** SBI
   posterior-predictive Re_k2/Im_k2 sit at the prior-predictive median and are
   ~1.7σ / ~2.9σ below the observed values, vs ~0.6σ / ~1.2σ for the reference
   MCMC. **The MCMC reference is authoritative for any k2- or dissipation-derived
   quantity; do NOT report SBI posteriors for Re_k2, Im_k2, log10_zeta, or the
   eta viscosities at the Titan datum.** The constraint the flow misses lives in
   the joint/correlation structure and the pushforward, not in any single
   marginal — so this is a flow deficiency, NOT the previously-claimed "benign
   prior-dominated nuisance." (eta_V, eta_sil weakly constrained.)
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

The limits-monotonicity N/A and the config-hash mode-only diff are intrinsic to
the joint mixture design; inherit those unchanged. **But do NOT inherit the old
"eta issue is benign and prior-dominated" framing** — it propagates the
mischaracterization the PPC corrected. Instead, MgSO4 and NaCl must:

1. **Run the PPC pushforward** (`titanG_ppc_interior_check.py`, now wired into
   the gate runner) and report, per channel, SBI-PP vs MCMC-PP median + a
   distributional distance (KS / 1-Wasserstein) gated at k× a matched-n
   bootstrap self-distance floor of the **MCMC-PP** (anchor on the MCMC-PP
   dispersion, NOT the measurement sigma). Also re-run the prior-predictive
   interior check: their ocean EOS shifts the envelope, so confirm the datum is
   not pushed into the tail (extrapolation).
2. **Promote the pushforward-observable crosscheck to a first-class gate** for
   the joint campaign (the per-param marginal gate provably cannot see a
   joint/pushforward miss). Sanity requirement: the new gate MUST flag the known
   NH3 k2-sector under-update when run on this artifact — a gate that does not
   catch it is mis-specified.
3. **Diagnose the flow under-update before prescribing more compute.** Likely
   cause is the sharp nonlinearity / near-degeneracy of the dissipation sector,
   not global undertraining (SBC passes, structural params clean). Do NOT assume
   "more epochs" fixes it: rule out MC noise first, then consider training
   weighted toward the informative Im_k2 tail, then an embedding net / alternate
   density estimator. Interpret; never tune to pass.

The gravity/structure sector is expected to remain cleanly verifiable for both
compositions; the split-status discipline (verified gravity/structure, k2 sector
deferred to the MCMC reference) travels with every joint artifact.

---

## Tidal-sector quarantine-lift adjudication (2026-10-20, scientific-reviewer + manager): NOT LIFTED YET, fiducial-scope lift pre-authorized

The verdict is binding. It is recorded here because this file is what
the GUI warning cites.

**Blocking gates (never run):**
- **G1**: crosscheck of the IS-CORRECTED posterior against the
  reference on all 13 parameters. The only crosscheck on disk
  (`crosscheck/`, 2026-08-03) is of the UNCORRECTED flow;
  `is_validation_nh3.json` `crosscheck_resample` holds only the
  resample indices (L=1270).
- **C11**: C20/C22 gravity no-regression. The reference lacks C20/C22
  columns.

Both must PASS. On PASS, a fiducial-scope lift takes effect without
further adjudication. On either FAIL, re-escalate.

**Passing at the fiducial (seeds 72/+1000/+2000, N=20000):**
- G2 pushforward medians: Im gap ≤0.19 σ_obs, Re gap ≤0.10 σ_obs.
- C2, C3 (max rel diff 0.0), C4.
- C5/C7 reliability set: Pareto-k ≤0.26, ESS ≥1229, w_max ≤0.0062.
- C5.3 reverse coverage.
- C6.
- C10: weighted KS ≤0.074, run against the n_eff=2000 reference; re-run
  against n_eff=8000 with G1.
- C13: 3-seed spreads ≤0.074 σ_obs.
- C14/C15.
- C16: 1.20σ against the n_eff=8000 reference.
- G4 corrected SBC: PASS on all 13 parameters, 485 pairs, 20.9%
  excluded.

**Dispositions:**
- **Plan gate 3** (ESS/N ≥0.1):
  - Fiducial half: FAIL-ADJUDICATED (0.0635), superseded by the C7
    reliability set, which passes.
  - Sweep half: replaced by C12, which FAILED (29/149; 17/114 = 14.9%
    inside the GUI Im_k2 box).
  - Discharged ONLY as a scope restriction: corrected at the default
    datum, uncorrected and quarantined everywhere else (architecture
    C8(i)).
- **C8(iii)** (x-keyed correction for user-modified x): NOT authorized.
  It needs its own preregistered gate set, including a substitute for
  C5.3 and a Re_k2 `x_obs_limits`.
- **Scope of the lift, once it takes effect:**
  - Re k₂ and |Im k₂| pushforwards: lifted.
  - ζ and η: lifted only for the parameters that pass G1 under the C9
    rule.
  - Ocean fraction: shows the corrected value (≈0.929).
  - The model–data k₂ tension caveat stays in every case.
- **GUI conditions B1–B9** (hashed precomputed bundle, exact-match
  activation, load-time diagnostic recompute, n_derived guard, weights
  in every panel, ESS shown rather than N, manifest-sourced numbers,
  screenshot verification, C8(iii) disabled) and the approved
  warning texts are in MACHINE-B-HANDOFF §0.32.

### R1/R2 outcome (2026-10-20): ALL PASS, the fiducial-scope lift TAKES EFFECT

- G1: all 13 parameters pass the stock gates against the n_eff=8000 pool,
  so the C9 fallback was never needed.
- C10 re-run: Re KS 0.032, |Im| KS 0.048; median gaps 0.055 / 0.098 σ_obs.
- C11: C20 shift 0.007 σ_obs, C22 shift 0.0975 σ_obs against a 0.1 gate.
  The C22 margin is thin and is recorded rather than reinterpreted.

Full record: `is_correction/r1_r2_RESULT.md` and `r1_r2_report.json`.

Lifted at the default datum only: Re k₂, |Im k₂|, ζ, the η's and alpha.
The k₂ model–data tension caveat stays. The GUI flag stays OFF until the
B8 visual verification is done.
