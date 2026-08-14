# Prior + observational constraints per mission–body 1D model

Machine A, 2026-07-19 (user-directed: "prepare prior and posterior
constraints for the models I named based on the provided references").
Companion to plans/mission-body-1d-roadmap.md. Status: constraint spec —
configs authored from these tables follow the roadmap discipline
(mission-appropriate observables only, hypothetical channels labeled,
full gates + reference MCMC + user ratification per artifact).

Conventions: all gravity coefficients UNNORMALIZED (the 3GM simulation
papers state their sigmas unnormalized; Galileo/Cassini literature values
are unnormalized). "Hypothetical channel" = no measurement exists; the
channel is retained at theory width for GUI exploration, labeled as such
in config metadata + scope note (Galileo v1.1 pattern, user-ratified).
Zeta split (user 2026-07-19) applies to every model below:
`log10_zeta_ice` (param_groups → log10_zeta_Ih + log10_zeta_HP where HP
ices exist) + `log10_zeta_sil`, each U[-3, 2].

## Shared prior blocks (from the ratified Europa/Callisto configs)

| Parameter | Prior | Notes |
|---|---|---|
| alpha | U[0.15, 0.45] | Andrade exponent, all bodies |
| log10_zeta_ice, log10_zeta_sil | U[-3, 2] | split per user; grouped Ih+HP for ice |
| log10_eta_Ih | U[12, 17] | body configs may narrow (existing: Europa [12,17]) |
| log10_eta_HP | U[12, 18] | HP-ice bodies only |
| log10_eta_sil | U[18, 22] | |
| Tb_K | cache grid span | per-body box, set by the structure grid |
| R_core_km | U[0, 0.5 R_body] | no-core limit admitted |
| rho_core_kgm3 | U[5000, 8000] | metallic Fe/Fe-S; see Callisto caveat |
| dC20_nh, dC22_nh | body-specific (below) | only with gravity_forward_model |

## 1. Galileo–Ganymede (upgrade of ganymede_pureh2o_andrade_8D; cache EXISTS)

Data (real):
- CMR2 = **0.3115 ± 0.0028** (Anderson et al. 1996; HYDROSTATIC-derived:
  J2/C22 = 10/3 imposed, Radau–Darwin — same caveat class as Europa GC21.
  Gomez Casajus et al. 2022 (Juno G34 + reanalyzed Galileo) find degree-2
  compatible with hydrostatic at 1σ but assign a LARGER realistic MoI
  error once non-hydrostatic degree-2 is admitted — config knob
  `moi_prior` with the Anderson solution default, GC2022 variant noted.)
Hypothetical channels (labeled; no Galileo k2/h2 measurement):
- Re_k2 [0.45, 0.10], Im_k2 [0.008, 0.10] — ocean-consistent theory
  central (ocean-bearing Ganymede k2 ~0.3–0.6; no-ocean < 0.1, Van
  Hoolst et al. 2024 discussion); Im central = |k2|·sin(1°) parity with
  the Europa convention. WIDE sigmas: exploration channels.
- Re_h2 [1.3, 0.2], Im_h2 [0.0, 0.2].
Gaps (documented, not modeled): Galileo magnetometer induction evidence
(Kivelson et al. 2002) — NO Be1xyz_Ganymede excitation table in-tree, so
no Bind/Ae channels until one is added (future work; roadmap).
Params (9): alpha, log10_zeta_ice(group Ih+HP), log10_zeta_sil,
log10_eta_Ih, log10_eta_HP, log10_eta_sil, Tb_K, R_core_km,
rho_core_kgm3.

## 2. Galileo–Callisto (revision of callisto_nacl_andrade_8D; caches EXIST)

Data (real):
- CMR2 = **0.3549 ± 0.0042** (Anderson et al. 2001; hydrostatic
  assumption ESPECIALLY suspect for Callisto — slow rotator, small
  q_r, and the non-hydrostatic literature (McKinnon; Gao & Stevenson)
  allows a substantially larger true MoI, up to an undifferentiated
  body. Metadata caveat REQUIRED; optional inflated-sigma variant
  [0.3549, 0.010] pre-registered for a sensitivity check.)
- Galileo degree-2 (Anderson 2001, unnormalized): J2 = 32.7 ± 0.8 e-6,
  C22 = 10.2 ± 0.3 e-6 — usable as C20/C22 observables via
  gravity_forward_model='clairaut_hydrostatic' + nuisances INSTEAD of
  the RD-derived CMR2 (v4 pattern; preferred variant. dC20_nh/dC22_nh
  priors U[-5e-6, 5e-6]: brackets the hydrostatic-vs-measured tension
  at Callisto's scale with margin).
- Induction: Be1xyz_Callisto EXISTS; existing config's Ae channels +
  synodic support cut stay (Zimmer et al. 2000 near-conductor
  amplitude — keep the existing channel sigmas; they are
  Callisto-derived, not leaked from another body).
Hypothetical: k2/h2 wide (no measurement): Re_k2 [0.3, 0.1] (ocean) /
Im_k2 [0.005, 0.1]; h2 as Ganymede.
Params (10): as Ganymede + (compositional axes stay fixed per config —
NaCl and MgSO4 variants exist; salinity sampling is a non-goal here).

## 3. Cassini–Enceladus (NEW; template PPEnceladus.py exists; cache = Machine B build)

**PRIMARY SOURCE UPDATED (user 2026-07-19): Park et al. 2024 (JGR
Planets 129, e2023JE008054; papers/park2024global.pdf) and its erratum
supersede Iess et al. 2014 gravity + Thomas et al. 2016 libration.**
Their Case 2 (quadrupole +
J3) is the paper-recommended field; unnormalized, reference radius
**256.6 km**, formal 1σ:

Data (real, Park et al. 2024 Table 8 Case 2):
- **C20 = −5477.45 ± 36.99 e-6** (J2 = 5477.45e-6)
- **C22 = 1517.90 ± 14.70 e-6**
- **corr(C20,C22) = 0.47**
- Measured J2/C22 = 3.61 — non-hydrostaticity STRONGER than Iess's
  3.51. With Enceladus's body-specific hydrostatic ratio (~3.22–3.25
  at q_r ≈ 6.3e-3), hydrostatic J2 from the measured C22 is
  ~4890–4930e-6 → **ΔJ2_nh ≈ +5.5–5.9e-4**. Nuisance priors
  ASYMMETRIC per review, widened for the Park field:
  **dC20_nh U[-1.5e-3, 1.5e-3]** (~2.6× the detected J2 excess),
  **dC22_nh U[-1e-4, 1e-4]** (C22 near-hydrostatic; a wide common
  prior would destroy the C22 channel). Pre-register
  prior-sensitivity on both.
- **Forced libration amplitude = 0.092° ± 0.009° (3σ → 1σ ≈ 0.003°)**
  (ERRATUM value, 2026-07-19 — originally published as 0.091°; Eq. 6
  amplitude 0.092119°) — REVISED DOWN from Thomas et al. 2016's
  0.120° ± 0.014° (2σ).
  With the merged Librations.py forward model (main → genai,
  2026-07-19) this is now a WIRABLE observable channel — the
  ocean/shell-decoupling constraint enters the likelihood instead of
  being a documented gap. New observable 'libration_deg' + runner
  channel (Machine A test runs).
- Park interior context (their inference, for fiducial/crosscheck
  only, NOT observables; ERRATUM values): mean shell 25–29 km
  (preferred 27), ocean 26–30 km (preferred 28), core radius ~197 km,
  core density 2290–2350 kg/m3 (preferred 2320), MoI 0.337
  [0.335, 0.338], total conductive heat loss 20–30 GW.
- **ERRATUM note (authoritative version of record):** the erratum
  corrects the shape/topography model and everything derived from it
  (libration 0.092°, interior ranges above, shape harmonics in
  Table 8's shape-derived columns) but does NOT touch the measured
  Case 2 gravity field — the C20/C22 observables above stand. It also
  adds the coefficient rescaling C_nm(new) = (R_old/R_new)^n
  C_nm(old), confirming this pipeline's degree-2 (R/R_ref)^2
  reference-radius convention.
- Out of 1D scope (documented): S22 = −275.31 ± 10.87 e-6 (significant
  lateral structure — needs a 3D/lateral model); C30 = 177.82 ± 33.42
  e-6 (degree-3 zonal; our forward model is degree-2 — J3 modeling is
  future work alongside higher-degree Clipper terms).
- **Nuisance-dominated honesty statement (review-binding):** the
  structural signal in C22_h (≈ k_f q_r/4 ~ 1.6e-3, k_f swing a few
  e-4) is comparable to the nuisance box, so C20/C22 alone constrain
  Enceladus's interior WEAKLY — the ocean/shell evidence is the
  implemented libration channel. Say so in the config scope note; do
  not present the gravity channel as a CMR2 replacement.
- CMR2 ≈ 0.335 (Iess 2014 interpretation) is DERIVED (hydrostatic+
  corrections) — do NOT double-count; use only as a sanity cross-check,
  not an observable, when C20/C22 are conditioned directly.
Remaining gap:
- k2 is unmeasured and is **omitted**, not represented by a hypothetical
  channel. The rejection is user-ratified in
  `plans/active/enceladus-config-freeze.md`: the proposed
  Re_k2 [0.015, 0.02], Im_k2 [0.005, 0.02] channel carried essentially
  no information and risked fake-measurement optics.
Params (9): as Ganymede minus HP ices (no III/V/VI at Enceladus
pressures → no log10_eta_HP; zeta_ice = Ih only), plus dC20_nh/dC22_nh.
Implementation prerequisites (below): body-dependent hydrostatic ratio +
reference radius — Enceladus q_r ≈ 6.3e-3 (fast rotator) makes the
Europa-hardcoded 3.324 WRONG here.

## 4. JUICE–Callisto (same caches as Galileo–Callisto; new sigmas)

Projected data (Cappuccio et al. 2022, PSJ 3GM simulation, 21 flybys,
unnormalized; Van Hoolst et al. 2024):
- Re_k2 [ocean-theory central 0.3, **0.059**] (3GM σ_k2, HAA on).
- Im_k2: NOT constrained unless Q < ~10 (Van Hoolst) — keep as
  hypothetical-width channel [0.005, 0.1], labeled.
- C20 [fiducial_h, **1.3e-7**], C22 [fiducial_h, **1.9e-8**]
  (σ_J2/σ_C22 from the 3GM sim; gravity_forward_model + nuisances
  U[-5e-6, 5e-6]; fiducials computed from the production
  _derive_gravity_pair at the body fiducial, v4 pattern).
Gap: obliquity (15 arcmin projected) — unsupported observable class.
Params (10): as Galileo–Callisto.

## 5. JUICE–Ganymede (flagship; orbital phase GCO500)

Projected data (Cappuccio et al. 2020 PSS sim + Van Hoolst et al. 2024):
- Re_k2 [0.45, **1e-4**] AND Im_k2 [0.008, **1e-4**] — the orbital
  phase constrains BOTH parts at 1e-4 (Van Hoolst): Im k2 becomes a
  real shell-viscosity/thermal-state measurement, not a hypothetical.
- Re_h2 [1.3, **0.04**], Im_h2 [0.0, 0.04] ("a few percent" accuracy,
  Van Hoolst; PULL the exact GALA number from the special-issue GALA
  paper before freezing — flagged).
- C20/C22: GCO500 continuous tracking — sigmas from Cappuccio et al.
  2020 (NOT fetched this session; flagged TO PULL before config
  freeze; expect ≲1e-8 class) via gravity_forward_model + nuisances
  U[-5e-6, 5e-6].
Gaps: obliquity 0.2 arcsec, libration 0.4–0.8 arcsec — unsupported
observable class (same roadmap item as Enceladus libration).
Params (9): as Galileo–Ganymede + dC20_nh/dC22_nh (11 total).
NOTE the conditioning-central caveat: at σ(k2) = 1e-4 the flow's
validated domain will be razor-thin around the training central —
grid-walk anchors in k2 are MANDATORY at this precision, and the
central should be refreshed to the actual measured value when JUICE
delivers (retrain).

## 6. JUICE–Europa (gravity-only complement to Clipper v4)

Projected data (Cappuccio et al. 2022, 2 flybys ~400 km, unnormalized):
- C20 [fiducial_h, **3.8e-6**], C22 [fiducial_h, **5.1e-7**]
  (hydrostatic-ratio test at σ(J2/C22) = 0.016). WORSE than Clipper v4
  (σ(C22) 2e-7) — value is the independent cross-check, not new
  constraint. Low priority; share the v3 2D cache + v4 config as base.
- k2/h2/rotation: 2 flybys insufficient (Van Hoolst) — no channels
  beyond the v4 set.
Params: = v4's 11.

## Implementation prerequisites (before ANY of these configs freeze)

1. **Body-dependent hydrostatic machinery in gravity_obs.py** — both
   currently Europa-hardcoded:
   - `J2_OVER_C22 = 3.324` is Europa's Tricarico-corrected value. The
     correction grows with q_r (ordering verified: Enceladus 6.3e-3 ≫
     Europa 5.0e-4 > Ganymede 1.9e-4 > Callisto 3.7e-5): Enceladus
     uses the frozen ratio **3.25**; Ganymede/Callisto sit nearer
     10/3 (~3.332). Make it a config key (`gravity_j2_over_c22`,
     threaded through hydrostatic_c20_c22 alongside the existing
     R_ref_m parameter) with the body value computed/cited at
     config-authoring time (Tricarico 2014 formalism).
     **Review-binding: correct ratio required for ALL C20/C22 configs,
     not just Enceladus** — at JUICE sub-1e-7 sigmas, using Europa's
     3.324 for Callisto/Ganymede (true ~3.332) injects an ~8e-8 C20
     systematic > σ(C22); the nuisances would absorb it and bias the
     non-hydrostaticity posterior.
   - `R_REF_GC21_M = 1565 km` is the Europa/GC21 reference radius.
     Config key `gravity_ref_radius_m`; per-body literature reference
     radii (Callisto 2410.3 km Anderson 2001; Enceladus **256.6 km,
     Park 2024 Table 8 Case 2** — verify each remaining body's paper
     reference radius when transcribing).
   - **Correlated degree-2 conditioning (review-binding):** the
     pipeline currently treats observables as independent Gaussians;
     the Enceladus (and Callisto) J2–C22 published solutions are
     correlated — the Park Enceladus detection is RATIO-driven
     (J2/C22 = 3.61, corr = 0.47), so the off-diagonal is material. Add
     a per-pair correlation option to the likelihood (2×2 covariance
     for the C20/C22 block); the Enceladus value is frozen and the
     Callisto covariance remains to be obtained before its freeze. SBI
     side: generate the noise draw from the same 2×2 covariance.
2. **Libration CLOSED for Enceladus; obliquity remains open.** The Van
   Hoolst libration formalism is implemented in `Librations.py`, wired
   through `mcmc_runner` and the SBI observable vector, and guarded by
   `tests/librations_test.py`. JUICE rotation-state obliquity remains an
   unsupported observable class.
3. **Be1xyz_Ganymede excitation table** — required for any future
   Ganymede induction channel (Kivelson 2002 evidence; JUICE J-MAG).
   MoonMag-side work.
4. **Enceladus structure cache** — PPEnceladus template exists, no
   cache; Machine B build (small body, cheap grid).
5. Flagged numbers to pull before the respective config freezes:
   Cappuccio et al. 2020 Ganymede GCO C20/C22 sigmas; GALA h2 sigma
   (special-issue paper); each paper's reference radius.

## References

- Anderson et al. 1996, Nature 384, 541 (Ganymede gravity, MoI 0.3115)
- Anderson et al. 2001, Icarus 153, 157 (Callisto gravity/MoI)
- Park et al. 2024, JGR Planets 129, e2023JE008054 + erratum
  (authoritative Enceladus gravity and libration)
- Iess et al. 2014, Science 344, 78 (superseded Enceladus gravity)
- Thomas et al. 2016, Icarus 264, 37 (superseded Enceladus libration)
- Gomez Casajus et al. 2022, GRL 49 (Ganymede gravity after Juno)
- Cappuccio et al. 2020, PSS (Ganymede 3GM simulation)
- Cappuccio et al. 2022, PSJ 3, 208 (Callisto + Europa 3GM simulation)
- Van Hoolst et al. 2024, SSRv 220, 54 (JUICE geophysics overview)
- SSRv JUICE special issue: link.springer.com/collections/dficcaejhd

## Review trail

Scientific review 2026-07-19 (opus): PASS WITH CONCERNS; all quoted
literature values independently verified (incl. Enceladus q_r = 6.26e-3
and the q_r ordering). Four binding changes INCORPORATED above:
(1) Enceladus ΔJ2_nh resized with the body-specific ratio (~4.0-4.5e-4,
not 2.7e-4); (2) asymmetric Enceladus nuisances + nuisance-dominated
honesty statement; (3) correlated J2-C22 conditioning (2x2 covariance)
as an implementation prerequisite; (4) body-specific hydrostatic ratio
required for ALL C20/C22 configs (8e-8-class systematic at JUICE
sigmas). The later 2026-08-12 freeze adjudication rejects the
Enceladus hypothetical k2 channel and supersedes that optional flag.
Confirm JUICE-Europa ratio-test sigma and Anderson-2001 covariance from
the papers before freeze.
