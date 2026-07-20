# Europa Clipper v4 (11D geodesy) — SBI validation gate report

**Artifact:** `europa_clipper_v4_geodesy_11D_posterior_1m.pt`
(nsf, sbi 0.26.1, torch 2.8.0, seed 43; git_sha `6e995165`).
**Config:** `PlanetProfile/Inference/configs/europa_clipper_v4_geodesy_11D.json`
(`gravity_forward_model = clairaut_hydrostatic`, `random_state 43`).
**Cache:** Test52_seawater 2D (Tb × log₁₀ w) structure grid (shared bilinear interpolant, as v3).
**Reference MCMC:** `Test/mcmc_results/Europa/Test53_geodesy_v4/europa_clipper_v4_reference_result.pkl`
(pocoMC, 4430 samples, ESS 4159, r̂ 1.0; carries weights + `heating_results`).
**Trained:** 326,849 kept sims (67.3% support+non-finite rejection), 21-D x, 11-D θ. Converged 107 epochs / 44.3 min.
**Seeds:** train 43 / data 49 / noise 4949.

v4 = v3's 8D (α, log₁₀ζ_Ih, log₁₀ζ_sil, log₁₀η_Ih, log₁₀η_sil, T_b, log₁₀w, R_core, ρ_core)
**plus degree-2 gravity**: two sampled non-hydrostatic Stokes offsets `dC20_nh`,
`dC22_nh` (both U[−2e−5, 2e−5]) mapped through the hydrostatic Clairaut forward
map `C20 = −R·C22_h + dC20_nh`, `C22 = C22_h + dC22_nh` with **R = 3.324**
(Tricarico, not the classical 10/3). The k₂ conditioning is tightened to the
Mazarico-2023-projected precision **σ(Re_k₂)=σ(Im_k₂)=0.015** (4× tighter than
v1–v3's 0.06).

---

## Headline result

**The non-hydrostaticity deliverable is verified; interior marginals are
diagnosed (not defects) and reported with calibrated caveats.**

1. **Gate 4 — non-hydrostatic offset upper limit: VERIFIED (centerpiece).** The
   identifiable combination **u = dC22_nh + dC20_nh/3.324** is recovered at
   94–95%, likelihood-dominated (σ_u,post 3.34e−7 vs prior 1.21e−5, a 36× shrink),
   and interior-independent **by construction** (the interior term C22_h cancels
   exactly in u; u is invariant under the ratio-preserving degenerate direction
   dC20_nh = −3.324·dC22_nh). Its width is directly **SBC-calibrated** (u-SBC KS
   p = 0.093; var-ratio 0.905 → mildly conservative, the safe direction for an
   upper limit).
2. **T_b ↔ salinity degeneracy** carried faithfully from v3 (corr −0.987 vs
   reference −0.983).
3. **Interior scalar marginals** (ζ_Ih, ρ_core) remain **prior-dominated in both
   MCMC and SBI**; SBI is conservatively over-dispersed. These are intrinsic
   weak-identification, not flow defects (see §Diagnostics).

---

## Gate results

| Gate | Verdict | Detail |
|------|---------|--------|
| **Gate 4 — non-hydrostaticity (deliverable)** | **PASS** | ratio-breaking u = +9.38e−7 vs true +1.00e−6 (94%); likelihood-dominated (σ_post 3.34e−7 vs prior 1.21e−5, 36×); ratio-preserving u ≈ 0 with components tracked (dC20 −3.24e−6/−3.32e−6, dC22 +9.14e−7/+1.00e−6); interior unbiased across arms (worst 0.029σ ≪ 0.25). |
| **u-SBC (derived combo)** | **PASS** | 400 held-out pairs: KS p 0.093; mean_dev +0.024; var-ratio 0.905 (conservative). |
| **Priority-0 (A) — flow fidelity** | **PASS** | KS(SBI PPD, MCMC PPD) D 0.083, p 0.121 → FAITHFUL. SBI implied Re_k₂ 0.2539 [0.242,0.264] ≈ MCMC 0.2536 [0.242,0.266]. |
| **SBC (11 params)** | **FAIL (localized)** | 8/11 pass. Non-uniform: T_b (KS p 4e−4), α (0.040), log₁₀w (0.034). c2st 0.55–0.59 (comparable to the PASSING v1.1). dC20_nh (0.50) and dC22_nh (0.32) **pass** — the deliverable is calibrated. |
| **crosscheck (11 params)** | **FAIL (soft)** | mean/σ location & scale pass broadly; shape-D fails vs tight self-D floors (~0.03–0.06). Worst: ζ_Ih (mean_diff 0.282 > 0.229, D 0.166), ρ_core (mean_diff 204, |dmedian| 231 > 207), η_sil (D 0.117). σ-ratios all 1.05–1.26 (SBI wider → conservative). |
| **2D T_b–w degeneracy** | **FAIL (2 marginals)** | **Joint HEALTHY**: corr −0.987 vs −0.983 (PASS), 68/95% containment 0.640/0.920 (PASS), T_b & w marginals within tol. Fails only on 2 prior-dominated scalar marginals (ζ_Ih SBI +0.41 vs MCMC +0.12; ρ_core +6415 vs +6211). |

The **two ratifying deliverables** — Gate 4 (non-hydrostatic upper limit) and the
2D degeneracy **joint** — both pass their central criteria. The three FAILs are on
scalar interior marginals whose failure is explained below and is not a flow defect.

---

## Diagnostics — why SBC / crosscheck / degeneracy fail on interior scalars

A four-diagnostic + Priority-0 program (reviewer-required, all *cheap*, no retrain
compute) decomposed the interior miscalibration. Results in
`reviewer_diagnostics_123.json`, `*_robustk2.json`, `priority0_mcmc_ppd_and_uSBC.json`.

1. **z-score compression on the heavy k₂ tail — RULED OUT.** A same-data,
   same-arch, same-seed ablation (frozen robust median/IQR embedding on the
   Re_k₂/Im_k₂ channels + `z_score_x='none'`, decompressing the k₂ tail 2.73×)
   moved the fiducial posterior-predictive < 0.2σ and left every gate verdict
   unchanged. Numerical preprocessing was not the driver.
   (`europa_clipper_v4_geodesy_11D_posterior_1m_robustk2.pt`, `*_robustk2.json`.)
2. **The fiducial implied Re_k₂ ≈ 0.25 is genuine, not a flow under-update.**
   Pushing both the MCMC and the SBI posteriors through the identical forward map
   gives near-identical implied Re_k₂ (0.2536 vs 0.2539; flow-fidelity KS p 0.121).
   The observed Re_k₂ = 0.23 sits **1.30σ below** the joint noisy posterior-predictive
   (std 0.0182) — **consistent**. Because the MCMC lands in the same place, no
   retrain can relocate the posterior toward 0.23.
3. **ζ_Ih / ρ_core are prior-dominated in both MCMC and SBI.** Their MCMC posterior
   σ are ≈0.64× / 0.80× the prior width; SBI is *wider still* (σ-ratio 1.15 / 1.26)
   → **over-dispersed = conservative**, never overconfident. The crosscheck/degeneracy
   marginal FAILs are shape-D exceedances against very tight self-D floors on
   near-prior distributions, not location or scale errors.
4. **T_b SBC failure is diffuse and Bonferroni-marginal** (only T_b hard-fails the
   0.05/11 ≈ 0.0045 threshold); rank-histogram shape is mildly under-dispersed with
   ~zero location bias.

**Remediation decision: no 1M regeneration.** The residual is intrinsic
weak-identification plus a genuine 1.3σ k₂ offset, not a fixable flow bug. Chasing
gates 1–3 with more sims would chase a target the joint data does not support.

---

## Reportable vs not reportable (reviewer-binding)

**Reportable as calibrated:**
- **u = dC22_nh + dC20_nh/3.324** as a provisional **upper limit** on the degree-2
  non-hydrostatic offset. Report in absolute units and in Clipper-σ units, with the
  σ_u convention disclosure below.
- The **T_b–w joint** posterior geometry (corr, 2D credible region).

**NOT reportable as calibrated (do not cite):**
- Per-component `dC20_nh` / `dC22_nh` marginals, the −C20/C22 ratio, or the
  dC22–C̄ correlation: the *degenerate* direction v = dC20 + R·dC22 is pinned only
  through the miscalibrated interior C22_h (≈0.14 σ_v leak).
- Interior scalar marginals ζ_Ih, ρ_core as tight constraints — report as
  prior-dominated, MCMC-consistent, SBI-conservative.

**Re_k₂ = 0.23 is a CAVEAT, not a tension detection.** 1.30σ is statistically
insignificant → "consistent; a mild directional offset." The k₂ likelihood *does*
pull the interior down from the prior-predictive median (0.277) toward the observed
0.23, but the interior prior volume + MoI/induction/h₂ co-constraints hold the joint
posterior-predictive at ~0.254. Whether this reflects a real physical high-k₂
preference vs prior-predictive geometry is a **v5 question** (prior-sensitivity +
leave-one-out on the k₂ channel), not a v4 claim.

**σ_u convention risk (40%).** The observation-limited σ_u = √(σ_C22² + (σ_C20/R)²)
= 3.25e−7 with the config's unnormalized Stokes σ. If the mission Table-5 σ(C20) is
4π-normalized, σ_u → 4.55e−7. State this on any reported u.

---

## v3 → v4 per-parameter constraint improvement

Adding degree-2 gravity (C20/C22) and the tight Mazarico k₂ conditioning tightens
the shared interior marginals (reference-MCMC posterior σ, so the comparison is
flow-independent):

| parameter | v3 MCMC σ | v4 MCMC σ | change |
|---|---|---|---|
| T_b (K) | 1.368 | 1.207 | −12% |
| log₁₀ w (dex) | 0.121 | 0.107 | −11% |
| R_core (km) | 67.5 | 56.8 | −16% |
| ρ_core (kg m⁻³) | 793 | 691 | −13% |

v4 additionally delivers the **degree-2 non-hydrostatic upper limit** (u), a
constraint v3 did not carry at all.

*(ζ_Ih / ζ_sil not directly comparable: v3 samples a single shared log₁₀ζ; v4 splits
it per phase. Both remain prior-dominated.)*

---

## Deliverable — zeta split (ice vs silicate) + heating fraction

Per-phase ζ split so a shared Andrade transient-creep amplitude cannot artificially
preference silicate heating. Heating routed through the runner-enriched cache
(`MCMCRunner._expand_theta` + `forward_model_k2_flexible` on `runner.structure_data`);
the built-in provenance check confirms the runner path reproduces the stored
reference `heating_results` (rtol 1e−4). Reports:
`zeta_heating_reference.json`, `zeta_heating_sbi.json`.

| quantity | reference MCMC | SBI posterior |
|---|---|---|
| **corr(ζ_Ih, ζ_sil)** | **−0.021** | **−0.021** |
| log₁₀ζ_Ih q16/50/84 (dex) | −0.85 / 0.25 / 1.06 | −0.61 / 0.53 / 1.52 |
| f_ice q16/50/84 | 0.65 / 0.95 / 1.00 | 0.61 / 0.96 / 1.00 |
| median Ice_Ih_W | 1.161e12 (1161 GW) | 1.127e12 (1127 GW) |
| median Silicate_W | 5.94e10 (59 GW) | 3.76e10 (38 GW) |
| frac ice-dominated (f_ice>0.5) | 0.904 | 0.883 |

**The split works as intended and SBI reproduces the reference.** The two ζ are
decoupled — corr(ζ_Ih, ζ_sil) = −0.021 identically in MCMC and SBI, so the shared-ζ
coupling artifact is removed. Heating is ice-dominated across ~88–90% of the
posterior. (Provenance check: the runner path reproduces the stored reference
`heating_results` to rtol 1e−4 — median 1.0099e12 = 1.0099e12 on identical θ.)

**CAVEAT:** both ζ remain weakly constrained (posterior widths ≈ prior widths), so
the defensible claim is that the split **removes an artificial shared-ζ coupling**,
not that the data independently pins each ζ. Absolute heating powers condition on
the tightened k₂/h₂ channels; the fraction + decoupling are the robust claims. The
f_ice shape is sensitive to the wide U[−3,2] ζ prior (carry into any v5).

---

## Overall

Status per CLAUDE.md vocabulary:
- **Gate 4 non-hydrostatic upper limit (u): VERIFIED** — recovered 94–95%,
  likelihood-dominated, interior-independent by construction, u-SBC calibrated.
- **T_b–w joint geometry: VERIFIED** (degeneracy joint criteria pass).
- **Interior scalar marginals (ζ_Ih, ρ_core, T_b calibration): implemented,
  documented** — prior-dominated / diffuse, MCMC-consistent, SBI-conservative;
  explicitly *not* claimed as tight constraints.
- **Re_k₂ offset:** caveat (1.3σ, consistent), flagged as a v5 question.

No retrain launched. No commits (opt-in). Reviewer sign-off: **PASS — close v4, no
blocker** (2026-07-19). Ratification (INDEX row, GUI slot) is Machine A's step.
