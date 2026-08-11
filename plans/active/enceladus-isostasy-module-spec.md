# Enceladus isostatic-gravity module — implementation spec (B11/B12/B14/B15/B16)

User-confirmed fit path (2026-08-12). Implements the Hemingway & Mittal
2019 (Icarus 332; papers/hemingway2019enceladus.pdf) forward model inside
the 1D campaign, per the reviewer's shape-channel ruling (freeze doc).
Target: `PlanetProfile/Gravity/isostasy.py` + wiring into
`PlanetProfile/Inference/gravity_obs.py` / `mcmc_runner`.

## Physics (H&M equation map)

1. **Hydrostatic figure per interface** — EXISTS: `Librations.compute_
   eccentricities` (full iterative multi-layer Tricarico 2014; B16
   structural verification done — it returns per-layer polar/equatorial
   eccentricities (ep, eq), NOT first-order theory, which H&M show
   underestimates H22_hyd by ~4% → ~30% error on H22_nh). Convert
   (ep_i, eq_i) at each interface radius R_i to unnormalized H20_hyd,
   H22_hyd (triaxial semi-axes → degree-2 shape coefficients; write the
   conversion once, unit-test on a known ellipsoid).
2. **Non-hydrostatic surface topography** (shape INPUT, not observable):
   H_nh_t(l,m) = H_obs(l,m) − H_hyd(l,m) at l=2; H30_obs is wholly
   non-hydrostatic. Source: Tajeddine et al. 2017 per H&M Table 1
   (unnormalized, R_ref = 252.22 km): H20 = −3510 ± 3.9 m,
   H22 = 857 ± 1.3 m, H30 = 420 ± 4.5 m. Nimmo 2011 as the
   preregistered ablation (−3846 ± 178.9, 917 ± 19.4, 384 ± 4.8 m at
   252.1 km).
3. **Airy root (equal-pressures isostasy, H&M Eq. 12 / Hemingway &
   Matsuyama 2017)**: basal relief
   H_nh_b(l,m) = C2 * H_nh_t(l,m) · (ρ_ice/Δρ) · (g_t/g_b) ·
   (R_t/R_b)^l-ish factor per Eq. 12 exactly — TRANSCRIBE Eq. 12 from
   the PDF during implementation, do not trust this sketch; the
   equal-pressures form differs from Lambeck equal-masses and from
   Iess/McKinnon/Čadek conventions (H&M §2.4 explicitly). C2 ∈ [0,1] is
   the sampled compensation fraction (C2=1 pure Airy).
4. **Gravity from topography** (H&M Eq. 4 layer sum): each interface's
   topography contributes δC_lm ∝ 4π ρ_contrast H_lm R_i^(l+2) /
   [M R_ref^l (2l+1)] with the **finite-amplitude correction**
   (Wieczorek & Phillips 1998; H&M Eqs. 5-7 define the adjusted H±
   coefficients) — REQUIRED (~5% on J2 = 8σ if omitted). Contributions:
   surface H_nh_t (ρ_ice against vacuum) + basal H_nh_b (Δρ = ρ_ocean −
   ρ_ice, opposite sign). Total predicted C_lm = C_lm_hyd (existing
   Tricarico/Radau path per config gravity model) + C_lm_nh(isostatic).
5. **Σ_model (H&M Eq. 22)**: added variance on predicted C20/C22/C30
   and libration from shape uncertainty: Tajeddine σ_model = 5.3e-6 /
   1.7e-6 / 4.4e-6 and 0.00025 deg (Nimmo: 244.8e-6 / 26.5e-6 / 4.7e-6
   and 0.004 deg). Implement as config-declared additive variances on
   those observables' likelihood terms.
6. **Reference radii (B14)**: Park gravity at R_ref = 256.6 km;
   Tajeddine shape at 252.22 km; PP body radius 252.1 km (Archinal).
   One explicit, unit-tested conversion (l=2 factor (256.6/252.22)^2 =
   1.035 → 3.6σ on C22 if mishandled; l=3 factor 1.053).

## Branch handling (pending reviewer consolidation — user ruled no-ocean STAYS)

Ocean nodes: full Airy coupling with sampled C2. Frozen nodes: NO
ice/ocean interface → rigid (uncompensated) support: C_lm_nh = surface
term only (equivalent to C2 = 0 with no basal root). Expected physics:
frozen branch predicts ~12σ-discrepant C20_nh AND ~10x-small libration
— doubly discriminated by real predictions, honoring the user's
keep-no-ocean ruling. FINAL FORM awaits the reviewer consolidation
ruling; implement the ocean branch first (it is invariant).

## Acceptance gate (B13 — the whole design rests on this)

With Iess 2014 gravity + Thomas 2016 libration (0.120 ± 0.014 deg) +
Tajeddine 2017 shape + equal-pressures Airy + finite-amplitude
correction + ρ_shell = 925, ρ_ocean = 1020 kg/m³ (H&M's stated example
values), the module must reproduce H&M's preferred solution: shell
19-24 km, ocean 35-39 km, core 192-195 km, ρ_core 2340-2410 kg/m³
(their §4/abstract). Implement as a committed test that grid-scans
shell/ocean thickness at H&M's 1 km steps and checks the misfit minimum
lands in their box. FALLBACK on failure: shape → display-only, restore
dC20_nh/dC22_nh free boxes (reviewer-preregistered).

## Sequencing

1. (this spec) → 2. conversion helpers + unit tests (ellipsoid → H_lm;
   radius rescaling) → 3. ocean-branch module + finite-amplitude →
   4. B13 gate vs H&M → 5. branch handling per reviewer consolidation →
   6. wire into gravity_obs/mcmc_runner behind
   gravity_forward_model = 'clairaut_hydrostatic_plus_airy_isostasy' →
   7. B1'/B2'/B5 scans on the assembled model → 8. config freeze.
Transcribe Eqs. 4-7, 12, 19-22 from the PDF verbatim at step 3 — the
sketches above are structural, not authoritative.
