# Approved Plan: Test52 — Differentiated No-Ocean Titan (k2 + CMR2, Petricca et al. 2025)

Approved by user 2026-07-08. Design by opus agent (grounded in PPTest50.py, Petricca 2025 PDF,
live compute_cmr2 evaluations). Companion to plans/monte-carlo-sbi-implementation-plan.md.

## Approved decisions

1. **CMR2 discretization-offset ANCHOR CORRECTION**: compute_cmr2 on the reduced cache
   under-reports CMR2 by a stable +0.00095 vs PP-native CMR2mean (pure hydrosphere
   radial-discretization error; core-independent, verified constant across grid and across
   (R_core, rho_core)). Fix: store per-grid-point `CMR2_offset = CMR2_pp_native − recon_nocore`
   at cache build time; add it in `_derive_cmr2_via_mass_conservation` before return.
   Regression test pins offset = 0.00095 ± 1e-4 for the Test52 grid. This CHANGES likelihood
   numerics for derived-rho configs — approved explicitly. Callisto impact: +0.23σ shift
   (σ=0.0042); disclose when Callisto runs are next regenerated.
2. **rho_core prior: rock differentiation [2591, 3600] kg/m³.** Lower bound = Petricca 2025
   bulk rock density. Upper 3600 (anhydrous silicate) NOT paper-sourced — user-blessed
   2026-07-08; cite as "user-adopted anhydrous silicate bound" until a citation is supplied.
3. **Density-inversion guard ON**: reject samples with rho_sil > rho_core (gravitationally
   unstable); count as support rejection like other guards. New logic (Callisto pattern has
   none); applies only when derived_params mass_conservation active AND guard enabled in config.
4. **Cache bodyname 'Test50' accepted** (PPTest50 template). Revisit when induction observables
   added (see memory: inference-roadmap-induction).

## Observables (Petricca et al. 2025, Nature 648:558, papers/petricca2025titan.pdf)

Re_k2 = 0.608 ± 0.048; Im_k2 = 0.135 ± 0.035; CMR2 = 0.343 ± 0.001 (Radau–Darwin).

## Config: configs/test52_titan_noocean_andrade_10D.json

10D = Test50 8D (alpha, log10_zeta, log10_eta_{Ih,III,V,VI,sil}, Tb_K — bounds verbatim from
test50_titan_noocean_andrade_8D.json) + R_core_km [0, 2000] + rho_core_kgm3 [2591, 3600];
derived_params rho_sil_kgm3 mass_conservation bounds [1800, 3500] reject_if_outside_bounds;
phase_stability no_ocean_Ih margin 0.1; planet_template_module PlanetProfile.Test.PPTest50;
random_state 42; structure_cache_path
PlanetProfile/Test/mcmc_results/Titan/Test52_andrade_noocean_diff/titan_diff_noocean_structure_grid.pkl.
Full JSON draft in opus design report (this session); implementer copies it verbatim, then adds
`density_inversion_guard: true` under derived_params (new key per decision 3).

## Physical sanity anchors (verify after build; offset-corrected)

no-core CMR2 = 0.34299 (0.0σ); R_core=600 km rock = −1.1σ; 900 km rock = −3.4σ;
600 km metallic = −4.1σ; 900 km metallic = −28.6σ. R_oceanbot ≈ 2081 km, M_hydro ≈ 3.81e22 kg,
Mtot = 1.3452e23 kg, R = 2574.73 km (PPTest50.py L22-23, Iess et al. 2010).
Petricca rock (R=2026 km, ρ=2591) reproduced by no-core reconstruction (2081 km, 2554 kg/m³).

## Registry additions (parameter_registry.py)

ParameterDefs R_core_km, rho_core_kgm3 (sketch in design report); preset
'andrade_titan_noocean_diff_10D' with observables BYTE-MATCHING the config JSON
(Im_k2 spelling — the documented landmine).

## Phasing

| Phase | Task | Agent | Gate |
|---|---|---|---|
| 1 | Config JSON + registry entries + grid build (build_phase_c1_cache, PPTest50 template, 9 Tb pts) + build-checklist verification | sonnet | checklist incl. offset stability |
| 2 | Offset storage in cache_builder + anchor correction & inversion guard in _derive_cmr2_via_mass_conservation + regression tests | sonnet | AFTER Callisto regeneration agent completes (same file); likelihood diff reviewed |
| 3 | Smoke MCMC (n_effective≈50, seed 42): rejection<20%, CMR2 posterior-predictive centered 0.343, small-core posterior | sonnet | physically sane |
| 4 | Numerics review + posterior vs Petricca + manuscript notes | opus | sign-off |

## Do NOT touch

structure_derivation.py, forward_models.py, any Test50 file/config/PPTest50.py, other body
configs, k2 likelihood numerics, sbi_runner.py, validate_sbi.py, Perple_X, manuscripts.

## Build checklist (Phase 1 gate)

1. phases unique {1,3,5,6,30,50}; hydrosphere contiguous; extract_hydrosphere_layers non-empty;
   R_oceanbot ≈ 2081 km; M_hydro ≈ 3.81e22 kg.
2. R_body_m ≈ 2.57473e6, Mtot_kg ≈ 1.3452e23 every grid point.
3. Cached no-core CMR2 ≈ 0.343 (PP native); recon-no-core ≈ 0.342; offset 0.00095 ± <1e-4
   across grid — RECORD per point.
4. CMR2 spread across (R_core, rho_core) reproduces sanity table (not inert).
