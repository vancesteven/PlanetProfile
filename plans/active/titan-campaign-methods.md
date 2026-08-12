# Titan free-gravity SBI campaign — methods

**Scope.** A paper-facing methods section for the Titan free-gravity simulation-based
inference (SBI) campaign: the forward model and observable design, the amortized
posterior and its importance-sampling (IS) correction against reference MCMC, the
C1–C16 validation gates, and the Track-2a information-gain decision. Numbers,
priors, and gate outcomes are drawn from the deployed configs and validation
reports; each subsection cites the source artifact so the section can be audited
against the run record.

**Status honesty.** This document records methodology and *current* validation
state. One gate (C16, ocean-fraction) is OPEN — it was reopened after the 3-seed
stability check and is under a manager-gate STOP pending a larger-N corrected run
and a higher-resolution reference. It is documented below as an unresolved,
escalated item, not as a passed gate. Do not cite the corrected NH3 ocean fraction
as ratified.

---

## 1. Interior model and observable design

Each Titan artifact conditions an amortized posterior on the observable vector

    x = {C20, C22, Re k2, Im k2},

with no CMR² (moment of inertia), no magnetic-induction, and no h2 channels. The
rationale is body-specific (configs `test54_titan_{nh3,nacl,mgso4}_freegrav.json`,
`titan_freegrav_noocean.json`):

- **CMR² is dropped to avoid a double count.** The published MoI (0.343, Petricca
  et al. 2025) is derived from the gravity figure within the hydrostatic framework
  (the Radau–Darwin relation maps MoI to J2 = −C20, and C22 is tied to J2 through
  the hydrostatic 10/3 ratio). Conditioning on both the C20/C22 figure and the
  figure-derived MoI applies the same hydrostatic relation twice and biases the
  posterior toward hydrostaticity. MoI is
  retained only as a candidate *inference-time* reweight (a future, independent
  scalar constraint), never as a likelihood term in this campaign.
- **Induction and h2 are dropped for lack of a clean signal.** Saturn's field is
  spin-aligned and Titan's ionosphere screens the inducing signal, so there is no
  clean Cassini induced-dipole channel; there is no measured Titan h2.
- **k2 and C22 do not double-count.** The static figure (C22) and the periodic
  tidal response (k2) are physically distinct and information-additive; both are
  kept.

**Gravity forward model.** `gravity_forward_model = "clairaut_hydrostatic"`. The
hydrostatic figure is set by the Clairaut Love number and the rotational parameter
q_r = ω²R³/(GM); for Titan (Torb = 15.945 d, R = 2574.73 km, M = 1.3452×10²³ kg)
q_r = 3.956×10⁻⁵ and k_f (Clairaut) = 1.0138. Non-hydrostatic offsets enter as two
sampled nuisances (dC20_nh, dC22_nh; §2) so the gravity channel is interior-agnostic.

**Non-hydrostaticity is a null signal for Titan.** Petricca's ratio −C20/C22 =
33.511/10.107 = 3.316 ≈ the classical synchronous value 10/3 = 3.333. The offset
nuisances are sampled for agnosticism but are expected to sit near zero (contrast
Europa SOL-A, ratio 3.157). The forward map uses the classical ratio
`gravity_j2_over_c22 = 3.3333…` and `gravity_ref_radius_m = 2.575×10⁶ m`. The
Radau–Darwin vs. Clairaut k_f gap (−0.76%: 1.0216 vs. 1.0138) is the RD
one-parameter approximation error and is tracked as a known frame effect, not
absorbed into the nuisances.

**Gravity provenance.** C20 = −J2 = −3.3511×10⁻⁵ (J2 = 33.511 ± 0.471 ×10⁻⁶);
C22 = 1.0107×10⁻⁵ (10.107 ± 0.062 ×10⁻⁶); R_ref = 2575 km; MoI = 0.343 — Petricca
et al. (2025), *Nature* 648:558.

**Compositions.** Ocean composition cannot be a sampled axis within one flow (the
structure cache bakes in one EOS), so the ocean-bearing campaign is three separate
artifacts — NH3, NaCl, MgSO4 — each with its own 2D (Tb × log10 w) structure cache
and its own reference MCMC, plus a no-ocean artifact. The NH3 config additionally
admits frozen (no-ocean) interiors into the *same* posterior as ocean interiors
(joint no-ocean + ocean; frozen nodes build as real no-ocean structures rather than
being silently rejected), so the NH3 posterior is not implicitly conditioned on
"an ocean exists."

---

## 2. Sampled parameters and priors

Thirteen parameters are sampled (all uniform priors; log10 priors are Jeffreys in
the underlying quantity). Eleven are shared across all compositions; Tb_K and the
salinity ceiling are composition-specific.

| Parameter | Prior bounds | Meaning |
|---|---|---|
| alpha | [0.15, 0.45] | Andrade rheology exponent |
| log10_zeta | [−3.0, 2.0] | Andrade pseudo-period prefactor |
| log10_eta_Ih | [10.0, 16.0] | ice Ih viscosity (Pa·s) |
| log10_eta_III | [10.0, 16.0] | ice III viscosity |
| log10_eta_V | [10.0, 16.0] | ice V viscosity |
| log10_eta_VI | [10.0, 16.0] | ice VI viscosity |
| log10_eta_sil | [18.0, 22.0] | silicate viscosity |
| Tb_K | composition-specific | ocean-bottom / ice-base temperature (K) |
| log10_wOcean_ppt | composition-specific | ocean salinity (g/kg solution) |
| R_core_km | [0.0, 2000.0] | metallic-core radius |
| rho_core_kgm3 | [2591.0, 3600.0] | core density |
| dC20_nh | [−2.0×10⁻⁵, 2.0×10⁻⁵] | non-hydrostatic C20 offset |
| dC22_nh | [−5.0×10⁻⁶, 5.0×10⁻⁶] | non-hydrostatic C22 offset |

Composition-specific ranges:

| Composition | Tb_K | log10_wOcean_ppt (w ceiling) | Seeds (data/noise/train) |
|---|---|---|---|
| NH3 | [249.0, 263.0] | [0.0, 1.8451] (70 ppt, 7 wt%) | 72 / 7272 / 72 |
| NaCl | [233.0, 272.0] | [0.0, 2.4624] (290 ppt, SeaFreeze wMax) | 74 / 7474 / 74 |
| MgSO4 | [248.0, 258.0] | [0.0, 2.2878] (194 ppt, 2-molal) | 73 / 7373 / 73 |

Silicate density rho_sil_kgm3 is derived by mass conservation (bounds [1800, 3500]
kg/m³, `reject_if_outside_bounds`, `density_inversion_guard`) — not sampled.

**Observation noise** (identical across compositions, diagonal Gaussian, θ-independent):
σ_C20 = 4.71×10⁻⁷, σ_C22 = 6.2×10⁻⁸, σ_Re_k2 = 0.048, σ_Im_k2 = 0.035.
Imaginary convention is `abs` (the absolute fold is applied to the noiseless forward
Im k2, then independent Gaussian noise is added).

**Physics-anchored support guard** (not tuned to pass): k2_support_bounds
Re_k2 ∈ [−0.1, 1.5], |Im_k2| ∈ [0, 1.0], applied to the *noiseless* forward output
before observation noise. Re_k2 ≤ 3/2 is the fluid-limit secular Love ceiling;
|Im_k2| ≤ 1 because dissipation cannot exceed the elastic response. A guard rejection
is counted as a support rejection so the SBI training support equals the
reference-MCMC effective support. The rejection is a deterministic function of θ
alone (independent of the noise realization), which is what makes the
additive-Gaussian-noise reduction in §5 exact on the guarded prior.

---

## 3. Amortized posterior and reference MCMC

The deployed artifact is a normalizing-flow posterior q(θ | x) trained on ~1M
(θ, x) pairs generated from the priors of §2 with the guard of §2 applied. The
reference against which the amortized posterior is validated is a pocoMC run on the
*same* observable set and likelihood.

**NH3 reference pooling.** The NH3 reference is the pool of three matched-resolution
seeds (72, 172, 272), each n_eff = 2000. Samples are concatenated across seeds;
weights are each renormalized to w/Σw/n_seeds so every seed carries equal posterior
mass (pocoMC weights are importance weights, not uniform), and the pooled weights
sum to 1. The matched seeds are chosen because the preregistered C16 ocean-fraction
band and the C9 nuisance-median bands are computed from exactly these three seeds;
mixing resolutions would compare a point value against a spread that does not
correspond to it. Output:
`validation_reports/titan_freegrav_nh3_1m/reference/titan_freegrav_nh3_reference_pooled.pkl`.

**Known limitation (C11 unevaluated).** The pooled NH3 reference stores
`k2_results` as (N, 2) = [Re_k2, Im_k2] only; it carries no C20/C22 pushforward
columns. The C11 gravity-no-regression gate reads columns for C20/C22 and therefore
cannot execute against this reference — no `gravity_noregression` verdict appears in
the NH3 report. Full C1–C16 ratification requires recomputing the per-sample C20/C22
for the pooled reference and running C11.

---

## 4. Importance-sampling correction

The amortized flow is a fast approximation; the IS correction reweights flow draws
by the exact likelihood to recover the target posterior for a specific observation
and to quantify the flow's approximation error.

**Weights.** For flow draws θ_i ~ q(· | x_obs),

    log w_i = log p(x_obs | θ_i) − log q(θ_i | x_obs).

The uniform-prior density cancels (all priors are box-uniform). log p is the
reference-MCMC likelihood closure evaluated on the deploy path; log q is the flow
log-density with `norm_posterior=False` (the truncation constant cancels). Support
guards enter as a −1×10³⁰ sentinel that is masked before any weight arithmetic.

**Reliability set** (preregistered; never tuned after seeing results):

- Pareto-k of the top weights (Zhang & Stephens 2009): ≤ 0.5 clean; 0.5–0.7 →
  PSIS smoothing; > 0.7 FAIL.
- Absolute ESS ≥ 1000 at the N actually run.
- Max single-draw weight fraction ≤ 0.01.
- Rejected fraction < 0.5.
- Per-branch (ocean / no-ocean): probability ≥ 0.05 and ESS ≥ 50 to gate that
  branch.
- ESS/N is REPORTED-ONLY (used to size N via N_required = 1000 / (ESS/N)); it was
  struck as a gate on 2026-08-11 — it is a cost statistic, not a validity criterion.

**The verdict** is `clean`, `smoothed` (PSIS applied), or `fail` (any fail reason).

---

## 5. Validation gates (C3, C5.3, C10, C11, C13, C16)

The gates fall into two groups: reliability of the weights themselves (§4) and
agreement of the corrected posterior with the reference.

- **C3 — likelihood self-consistency (byte-exactness).** Recompute the
  log-likelihood for 200 reference draws through the deploy-path runner; require max
  relative difference < 10⁻⁹ and zero sentinel disagreements. This proves the flow
  path and the reference path evaluate the *same* likelihood.
- **C5.3 — reverse-direction coverage (BLOCKING).** Evaluate the flow log-density at
  the reference-MCMC draws; require < 1% of reference mass to fall below the flow's
  0.1st-percentile self-log-density. Pareto-k is blind to regions where the target
  greatly exceeds the flow but the flow density is so low no draw landed there; C5.3
  is not.
- **C10 — tidal pushforward.** Weighted two-sample KS on the Re/Im k2 pushforward
  (KS distance ≤ 0.15) plus median-to-median agreement (≤ 0.5 σ_obs). The
  full-distribution KS is used because a width defect is not caught by a median-only
  test.
- **C11 — gravity no-regression.** C20/C22 pushforward-median shift ≤ 0.1 σ_obs.
  (Currently unevaluated for NH3; see §3.)
- **C13 — 3-seed stability.** The corrected medians and branch fractions must be
  stable across three conditioning seeds to within 0.1 σ_obs (and, for ocean
  fraction, within the reference 3-seed spread).
- **C16 — ocean fraction.** The corrected has-ocean mass must agree with the
  reference within a properly combined standard error (§6).

---

## 6. NH3 results and the C16 reopen

NH3, n_post = 20000, seed 72, config_hash `e596574d1e81567c`
(`validation_reports/titan_freegrav_nh3_1m/is_correction/is_validation_nh3.json`).

**Weight reliability (clean).** ESS = 1271 (ESS/N = 0.0635), Pareto-k = 0.191
(no PSIS smoothing), max weight fraction = 0.0062, rejected fraction = 0.011.
Verdict `clean`.

**Tidal deliverable (the campaign's main scientific result).** The IS correction
moves the Im k2 pushforward median from 0.0449 (flow) to 0.1089 (corrected), and
Re k2 from 0.541 to 0.594. The corrected Im k2 recovers the reference-MCMC tidal
ceiling. At the fiducial seed, C10 passes: the weighted KS distances are 0.061 (Im)
/ 0.041 (Re) and the corrected-vs-reference median gaps are 0.147 σ_obs (Im) /
0.069 σ_obs (Re), both well within the 0.5 σ_obs median tolerance. C13 confirms the
tidal medians are seed-stable: across the three conditioning seeds the
corrected-vs-reference Im/Re median gaps span 0.130–0.188 / 0.027–0.101 σ_obs — a
seed-to-seed spread of ≤ 0.074 σ_obs, within the 0.1 σ_obs C13 stability threshold.
C3 (max rel-diff 0.0, 200 draws) and C5.3 (0.094% reference mass below the flow tail
vs. 1% cap) also pass. **The tidal-sector deliverable is sound.**

**C16 is OPEN (reopened → manager-gate STOP).** Corrected ocean fraction = 0.9293
vs. reference 0.9173, a residual of +0.0120 at the fiducial seed. Under the SNIS
delta-method branch-fraction standard error (SE_corr = 0.00412, SE_ref = 0.00365,
combined 0.00550) the residual is z = 2.18 — above the preregistered reopen
trigger. The pass rule is |corrected − reference| ≤ 2·√(SE_corr² + SE_ref²), not the
raw 3-seed min/max band. The 3-seed check (ocean fraction 0.9293 / 0.9310 / 0.9363;
residuals +0.012 / +0.014 / +0.019) shows the residual is stable, all-positive, and
*not shrinking*. The between-seed std (0.00367) matches the SNIS delta SE and
falsifies the cruder ESS-based SE (0.00767, ~2× too large).

The decisive argument that this is **structural, not finite-N**: the +0.0149 mean
residual is ~19× the O(1/ESS) finite-N bias scale (1/ESS = 7.9×10⁻⁴ at ESS ≈ 1270
in the clean-Pareto regime, k = 0.19; the SNIS bias is O(1/ESS) up to a
problem-dependent constant, so this is a scale rather than a hard ceiling, but a
factor of ~19 is robust to that constant). More draws would sharpen the tension,
not close it.
This is a ~1.5% branch-mass inconsistency between two same-target (C3 byte-identical)
estimators, and the escalation diagnosis points at the reference side (n_eff = 2000
is modest). It is not a catastrophic failure mode — the no-ocean branch is preserved
and C5.3 passes.

Per the manager gate ("NH3 fail on C3/C5.3/C10/C16 ⇒ STOP all further Track-1
compute, escalate"), all NH3 Track-1 corrected compute is stopped and the item is
escalated. The resolution — none of which is self-authorized — is:

1. A larger-N corrected conditioning run (N ≈ 100k–200k, ESS ≈ 6k–13k) at the same
   seed(s), to test whether the residual stays ~+0.015 (structural) or drops with
   ESS (finite-N). This is Track-1 compute under the STOP and must be authorized.
2. A higher-n_eff reference recomputation (n_eff ≳ 8000) with re-verified pooled
   weighting, reporting the reference ocean fraction and its SE.
3. Re-ratify C16 only on agreement of the two same-target estimators within a
   properly combined SE, with the +0.015 source identified.

---

## 7. Europa v5 positive control

To confirm the IS machinery does not distort a well-calibrated posterior, the same
correction was run on the deployed Europa v5 geodesy artifact (11 params, 21
observables), where the flow is already well-calibrated and the correction should be
a near-no-op (`validation_reports/europa_clipper_v5_baseline_1m/is_correction/is_validation_v5.json`,
seed 43).

Result: ESS = 13487 (ESS/N = 0.674), Pareto-k = 0.256, max weight fraction = 0.0013,
verdict `clean`. Pushforward medians are effectively unchanged by the correction
(Im k2 0.002637 → 0.002635; Re k2 0.2525 → 0.2497); ocean fraction 1.0 → 1.0. C3,
C5.3, and C10 all pass. The correction is a clean no-op on a good posterior — which
corroborates that the NH3 C16 tension (ESS/N 0.06, a real ~1.5% branch-mass shift)
is a property of the NH3 reference/target pair, not an artifact of the correction
machinery.

---

## 8. Track-2a — information-gain per observable channel

Before committing costly downstream work (transform-retrain and Re/Im-only
ablations), we quantified how much interior information the Im k2 channel carries
beyond {C20, C22, Re_k2}. The estimand is the conditional mutual information

    CMI = I(θ ; Im_k2 | C20, C22, Re_k2)   [nats, effective-prior],

with a cancel threshold of 0.1 nat (cancel downstream work only if the upper bound
of the combined uncertainty band is below threshold — i.e., require positive
evidence *against* information before cancelling).

**Estimator.** Because the observation noise is additive, diagonal-Gaussian, and
θ-independent, the conditional law x | θ has constant (analytic) entropy, so for any
channel subset S, I(θ; x_S) = H(x_S) − Σ_{c∈S} ½ ln(2πe σ_c²). By the chain rule the
CMI collapses to a low-dimensional differential-entropy difference on the
noised-observable marginal:

    CMI = [H(x₄) − H(x₃)] − ½ ln(2πe σ_Im²) = H(Im_k2 | C20,C22,Re_k2) − ½ ln(2πe σ_Im²).

The 13-dim θ never enters the estimator — it is integrated out exactly by the
noise-model identity (the guard rejections are θ-only, so the reduction holds on the
guarded prior). H is estimated with the Kozachenko–Leonenko k-NN estimator over
k ∈ {3, 5, 10, 20}, with the three conditioners rank-Gaussianized (an invertible
reparameterization that tames k-NN boundary bias without changing the estimand) and
σ_Im = 0.035 subtracted exactly in the same units. Uncertainty is from m-out-of-n
subsampling (20 draws of 200k, without replacement).

**Result.** NH3 gates the decision; the salts are context. The NH3 row is from the
ratified `track2a_cmi.json`; the NaCl/MgSO4 rows are from
`track2a-infogain-design.md` (the ratified JSON contains only NH3).

| Composition | n_kept | CMI band (nat) | Gaussian bracket | CMI ≥ 0 | Cancel 2b/2c |
|---|---|---|---|---|---|
| **NH3 (gate)** | 689,845 | [0.599, 0.615] | 1.211 | yes | **No** |
| NaCl | 631,214 | [0.616, 0.640] | 1.120 | yes | No |
| MgSO4 | 691,075 | [0.380, 0.402] | 0.737 | yes | No |

**Decision.** NH3 CMI ≈ 0.605 nat is ~6× the 0.1-nat threshold, so the
transform-retrain (2b) and Re/Im-only ablation (2c) are **not cancelled**: Im k2
carries substantial incremental information about the interior beyond the other
three channels. The estimand is an upper bound on the *any-θ-dimension* recoverable
information, so CMI > 0.1 forbids cancellation but is necessary-not-sufficient for
2b/2c to help. All reviewer cross-checks pass (CMI ≥ 0; raw vs. rank-transformed
agree to < 0.01 nat for k = 3, 5, 10; k-sweep range < 0.004 nat; subsample std
≈ 0.002 nat; the k-NN CMI sits below the moment-matched Gaussian bracket, as
expected). Track-2a reads only training (θ, x) pairs and no corrected/IS-reweighted
quantity, so it is independent of the C16 manager-gate STOP.

---

## Source artifacts

- Configs: `PlanetProfile/Inference/configs/test54_titan_{nh3,nacl,mgso4}_freegrav.json`,
  `titan_freegrav_noocean.json`.
- IS correction + gates: `PlanetProfile/Inference/is_correction.py`,
  `plans/scripts/is_correction_validate.py`.
- NH3 report: `validation_reports/titan_freegrav_nh3_1m/is_correction/is_validation_nh3.json`.
- Europa v5 control: `validation_reports/europa_clipper_v5_baseline_1m/is_correction/is_validation_v5.json`.
- Track-2a: `plans/active/track2a-infogain-design.md`,
  `validation_reports/track2a_infogain/track2a_cmi.json`.
- Reference pooling: `plans/scripts/pool_nh3_matched_reference.py`.
- Campaign master plan: `plans/fluffy-snacking-fountain.md`;
  remedy plan: `plans/active/tidal-sector-remedy-plan.md`.
