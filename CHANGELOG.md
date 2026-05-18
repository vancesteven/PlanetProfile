# Changelog

## [Unreleased] – genai branch

### Bug Fixes (narrow scope)
- **Corrected conduction profile added for clathrate depth calculation.** A new function `ConductiveTemperatureCorrect` in `PlanetProfile/Thermodynamics/ThermalProfiles/ThermalProfiles.py` implements the strict Turcotte & Schubert (2002) §4.9 eq. 4.40 form: `c1 = qTop·rTop²/k − Htot·rTop³/(3k)` with planar-limit ΔT = qTop·ΔR/k (Fourier's law). `GetPbConduct` now uses this function; the resulting clathrate layer thickness matches the user-specified `Bulk.clathMaxThick_m` (previously PP produced a clathrate 2× deeper than requested).
- **Silicate boundary-condition fix for solid-sphere bodies (T_center-anchored).** `SilRecursionSolid` and `SilRecursionPorous` in `PlanetProfile/Thermodynamics/Geophysical.py` now detect the solid-sphere case (no Fe core, no `CONSTANT_INNER_DENSITY`) and enforce the only physically admissible boundary condition with T finite at r=0: at each downward step they override the inherited `qTop` with the consistency value `qTop = Htot·rTop/3`, which makes `c1 = 0` by construction in `ConductiveTemperatureCorrect`. The closed-form profile is then `T(r) = T_top + Htot·(rTop² − r²)/(6k)`. Heat-flux consistency at the ice/silicate interface is a derived quantity rather than a prescribed one. Shell bodies (`Fe_CORE=True` or `CONSTANT_INNER_DENSITY=True`) continue to use the legacy `ConductiveTemperature` call to preserve existing inner-core behavior. This resolves the long-standing silicate BC issue that the legacy halved-c1 had masked: the legacy form attenuated a real divergence by a factor of 2, "passing" tests numerically while propagating an inconsistent `qTop` from the ice shell down through the silicate. Validated against Tests 1–36 of the BuildTest suite (no behavior changes for shell bodies; finite T(r=0) for solid-sphere bodies). Test 37 still fails with the pre-existing clathrate-underplate stub error documented at `LayerPropagators.py:76`, unrelated to this fix.

### Renames / Clarifications
- `Constants.triplePointT_K` → `Constants.TtripleIh_III_L_K` (with deprecation alias). Value unchanged at 251.165 K. The old name was misleading because the constant is the Ice Ih–III–liquid triple point, not the water Ih–liquid–vapour point at 273.16 K. Deprecation alias preserved in `defineStructs.py` for external scripts.

### Kalousova HP-ice clarifications
- **Ice VI under `KALOUSOVA_CONVECTION` now sets a melting-curve T profile** in the ice VI layer (`Convection.py::IceVIConvectSolid`). Previously, convection parameters (`eLid`, `Dconv`, `ΔTBL`, `Ra*`) were computed but the layer's `T_K` was left untouched and a runtime warning was emitted. The temperature is now uniformly set to `TconvVI_K` (which `ConvectionKalousova2018` already places on the melting curve), and layer properties are re-evaluated. The full three-segment construction (lid + adiabatic interior + lower TBL), as implemented for Ice III and V, remains a follow-up — it requires `Planet.Steps.nVbottom`, `nIceVILitho`, layer allocation in `IceLayers()`, and `IceVIConductSolid/Porous`.
- **`Planet.meltFractionVI` placeholder bumped from 0.01 to 0.5** when a temperate layer is detected, reflecting that an actively melting two-phase HP-ice layer is much more likely to sit near 50% melt than at the percolation threshold. `Constants.phiPercolation` (= 0.05) is intentionally *unchanged* — Kalousova & Sotin (2018) Eq. 10 outflow velocity / mass flux equations in `LayerPropagators.py` are conditioned on this value. Top-level `meltFractionIII` and `meltFractionV` placeholders remain at 0.01 pending the same two-phase solver. The in-ocean path (`LayerPropagators.py`) continues to set melt fraction = `phiPercolation` per Kalousova's outflow-balance assumption.
- **CLAUDE.md HP-ice section rewritten** to reflect the actual implementation status: temperate-layer detection done for III/V/VI; profile propagation done for III/V, simplified for VI; melt fraction is a placeholder (not a solver output). The previous wording overstated the implementation.

### Inference / MCMC — Phase C1 (cache v2.1, B-layered blend, Callisto)

- **fix(layer_propagators): HPIceConvection early-return when
  KALOUSOVA=False.** `HPIceConvection` in
  `PlanetProfile/Thermodynamics/LayerPropagators.py` previously
  dispatched HP ices (III/V/VI) to `ConvectionDeschampsSotin2001`
  when `KALOUSOVA_CONVECTION=False`, but D-S 2001 has a hardcoded
  `Tupper_K=274` melt-curve bracket
  (`ThermalProfiles.py:311–315`) that fails for HP ices at deep
  ocean pressures. The early-return now sets no-convection
  diagnostic defaults (eLid=zb, Dconv=0, Ra=0, RaCrit=∞,
  meltFraction=0) and leaves T on the melt curve as set by
  `SelfConsistentOceanLayer`. Pre-fix Titan / Ganymede caches are
  not corrupted — D-S 2001's phase-mismatch escape hatch
  (`ThermalProfiles.py:296–305`) returned a conductive-profile
  fallback, so T_K, phases, ρ, μ, K, η (the arrays MCMC
  consumes) were unaffected; only diagnostic fields (eLid, Ra,
  meltFraction) were silently wrong. Verified by Callisto MgSO₄
  rebuild: 6/9 Tb points previously failed with D-S 2001
  ValueError at low Tb; post-fix all 9 succeed (21 min vs 74 min
  that crashed).

- **feat(cache_builder): v2.1 transition-aware Tb grid.**
  `PlanetProfile/Inference/cache_builder.py` extended with
  `_bisect_transition` and a post-grid refinement phase. After the
  regular Tb grid is built, `build_tb_grid_cache` walks adjacent
  grid pairs and detects layer-set changes (different
  `region_phases` or `n_layers`). Each transition is bisected to
  ε_T = 0.01 K; the converged anchor pair flanks the transition by
  < ε_T. Cache schema bumped to `v2.1` with a `transitions`
  metadata list `[{Tb_lo, Tb_hi, regions_lo, regions_hi}, ...]`.
  Backward-compatible — old v2.0 caches load unchanged. Rationale:
  a Tb where a new layer first stabilises (HP ice III appearance,
  ocean appearance, Ih conv/cond split) is a physical
  discontinuity; bisecting it bounds the unblendable interval
  narrowly enough that nearest-neighbour fallback at MCMC time
  introduces no measurable bias.

- **feat(forward_models): B-layered structure blend.**
  `PlanetProfile/Inference/forward_models.py::_blend_b_layered`
  replaces the legacy element-wise blend. For each layer index k,
  boundary radii are blended linearly in `w`; per-cell continuous
  fields (ρ, K, μ, η, T, P, bulk_visc) are resampled onto s0's
  intra-layer normalised grid via `np.interp`; phases / discrete
  metadata are copied from s0 (identical to s1 by precondition);
  per-Tb scalars (CMR², D_iceIh_km, …) are linear-blended.
  `apply_bottom_temperature` dispatches: matching `region_phases` →
  B-layered blend; mismatched (across-transition bracket) →
  nearest-neighbour. Transition windows are < 0.01 K wide so the
  discontinuity is below any meaningful posterior precision. Why:
  naive element-wise blending across a moving boundary produces
  unphysical "mush" cells where ρ does not match phase identity.
  Includes a precondition check asserting every per-cell field has
  length equal to `len(r_m)`, surfacing cache-builder-padding bugs
  surgically by field name rather than from inside `np.interp`.
  12 unit tests in `tests/test_layered_blend.py` cover boundary
  linearity, variable cell counts, no-mush invariant, transition
  dispatch, scalar/meta handling, and w=0/w=1 endpoints.

- **fix(cache_builder): P_MPa included in MIN_POINTS padding pass.**
  `build_single_structure`'s thin-layer padding pass (lines 256–316)
  was extending r, ρ, K, μ, η, phases, bulk_visc, T_K via
  `np.interp` but not P_MPa. Caches with any padded layer therefore
  had `len(P_MPa) < len(r_m)` (typically by 4, the number of
  single-cell layers padded). The new B-layered blend accesses
  P_MPa per cell, surfacing the discrepancy; fixed by extending the
  padding pass to interpolate P_MPa identically to T_K. Verified
  zero field-length mismatches in rebuilt Callisto / Europa /
  Ganymede caches.

- **feat(cache_builder): ocean_overrides for composition switching.**
  New optional dict in v2.1 BodyConfig JSON: `ocean_overrides`.
  Applied to `Planet.Ocean.<key>` after deepcopy of the body's
  default template, before running PP. Enables composition switching
  and concentration sweeps without spawning N template Python files.
  Motivation: Callisto's MgSO₄ EOS only extends to P=800 MPa;
  Callisto's deep ocean frequently triggers extrapolation
  regeneration on every PP run, making MgSO₄ Callisto builds
  prohibitively slow. SeaFreeze's NaCl EOS extends to P=5000 MPa
  cleanly. New config `callisto_nacl_andrade_8D.json` switches to
  NaCl 100 ppt; the MgSO₄ config is retained as deprecated for
  record. Future SeaFreeze MgSO₄ release with wider P range can
  revisit.

- **chore(docs): inference toolkit methodology README.**
  New `PlanetProfile/Inference/README.md` (723 lines) — canonical
  methodology reference covering the pocoMC algorithm (preconditioned
  normalising flows + SMC tempering, n_effective semantics, ESS
  termination), v2.1 cache schema, B-layered blend mechanics,
  ocean_overrides usage, and end-to-end Phase C1 workflow
  (config → cache → smoke → production MCMC).

- **chore(caches): Phase C1 v2.1 caches built for Callisto, Europa,
  Ganymede.** Three structure caches with the v2.1 schema,
  transition refinement, and P_MPa fix:
  `callisto_nacl_structure_grid.pkl` (17 points, 1 transition, ~40 min
  build), `europa_seawater_structure_grid.pkl` (16 points, 1
  transition, 53 s build), and
  `ganymede_pureh2o_structure_grid.pkl` (23 points, 2 transitions,
  55 s build). Smoke MCMC against the Callisto cache completes
  cleanly: ESS=4096, acceptance=0.61, posterior appropriately
  prior-dominated for 1 observable in 8D. Note for future work:
  Gao & Stevenson 2013 argues slow rotators may not be hydrostatic,
  so Callisto's CMR² could be up to 10% lower than the Anderson 2001
  value. Inference variants at CMR² 5% and 10% below nominal (same
  σ=0.0042) should be built; the structure cache is independent of
  the observable and need not be rebuilt.

### Inference / MCMC — Phase B refactor
- **Test 50 refactored to thin wrapper** around `MCMCRunner(InferenceConfig.from_json(...))`. The 901-line monolithic script is now 528 lines (structure-building code unchanged; MCMC/forward/likelihood code replaced entirely by toolkit delegation). The 5D variant (`explore_test50_5D.py`) is now a 20-line shim that delegates to Test50 with `--config test50_titan_noocean_andrade_5D.json`.
- **`InferenceConfig` extended with `param_groups` and `fixed_params`** (`inference_core.py`). `param_groups` maps a sampled scalar (e.g., `log10_eta_HP`) to a list of member parameters that all receive the same value at runtime — enabling the 5D HP-locked variant without code changes. `fixed_params` injects constants (e.g., `Tb_K=250.965`) into every forward model call without entering the prior.
- **Body-agnostic no-ocean safeguard** added to `MCMCRunner._make_flexible_log_likelihood`. Enable with `sampler_settings.phase_stability = {"enforce": "no_ocean_Ih", "margin_K": 0.1}` in any config JSON. The guard rejects samples where any Ice Ih cell's interpolated T exceeds `Tm_Ih_lin(P) − margin_K`, where `Tm_Ih_lin = 273.16 − 0.068·P_MPa`. Independent of body name — usable whenever the no-ocean null hypothesis needs to be enforced.
- **Per-phase viscosity hooks**: `log10_eta_III`, `log10_eta_V`, `log10_eta_VI` now map to individual HP ice phases in `apply_viscosity_params`. The lumped `log10_eta_HP` still works; per-phase keys override it.
- **Single-zeta Andrade hook**: `log10_zeta` (no phase suffix) applies uniformly to all solid phases, matching Test50/Test48's original parameterisation. Per-phase `log10_zeta_Ih/HP/sil` overrides still take precedence.
- **Test50-format grid cache** (`{'Tb_K_grid': ..., 'structures': [...]}`) now supported natively by `MCMCRunner` and the `Tb_K` parameter hook. Linear interpolation is bit-for-bit identical to Test50's original `_interp_structure`.
- **Arrhenius reference temperature** is now dynamic: the `Tb_K_sampled` value stored by the `Tb_K` hook is used as `reference_temp_K` in `apply_arrhenius_viscosity`, matching Test50's `η(T) = η_Ih · exp(E/R · (1/T − 1/Tb))` semantics.
- **Two JSON configs** in `PlanetProfile/Inference/configs/`: `test50_titan_noocean_andrade_8D.json` (8D independent HP-ice η, Tb sampled) and `test50_titan_noocean_andrade_5D.json` (5D HP-locked, Tb fixed). Both include `arrhenius_params` (E_Ih=60 kJ/mol) and `phase_stability.no_ocean_Ih`.
- **Test 50 upgraded to 8D MCMC with sampled basal temperature** (`Test50_mcmc_andrade_noocean_yao2014.py`). Adds `T_b ∈ [249.0, 250.965] K` as an 8th parameter; structure cache pre-built on a 9-point Tb grid; forward model linearly interpolates structures between bracketing grid points per sample. All ice-phase viscosity priors broadened to `[10, 16]` (log₁₀, Pa·s) to admit Petricca-style low-η regimes. Forward model adds a no-ocean safeguard rejecting any Ih cell whose interpolated T crosses the linearized melt curve. PPTest50 default clathrate cap reduced from 10 km to 2 km; PPTest50 default `Tb_K = TtripleIh_III_L_K − 0.2 K`.

## [3.1.0] – 2026-01-13
**Author:** @Chang-Scott

### Optimizations
- Improved memory usage in the EOS framework.
- Implemented EOS pre-loading to support large-scale explorations.
- Updated computational pathways to reduce typical PlanetProfile runs to ~1–5 seconds.

### Plotting Extensibility
- Overhauled storage of large-scale exploration results, enabling:
  - Exploreograms to utilize inductogram plots.
  - Inductograms to utilize exploreogram results.
- Exploration data is now saved and reloaded as `.pkl` files instead of `.mat`.
  - `.mat` files are still generated for users who wish to post-process results in MATLAB.

### Monte Carlo
- Added an initial framework to enable Monte Carlo sampling of model properties.
- See individual commits for implementation details.
- *Note:* Monte Carlo functionality is a work in progress and has not yet been extensively tested.

### Non-Self-Consistent Modeling
- Established a framework to support layer modeling with prescribed constant properties rather than EOS-derived values.
  - For example, users can now directly specify ocean density, thermal expansivity, and related parameters without enforcing compositional self-consistency.
- *Note:* This functionality is a work in progress and has not yet been extensively tested.

### MgSO₄–SeaFreeze Coupling
- Implemented dynamic coupling between ice polymorph chemical potentials and MgSO₄ thermodynamic data to generate phase grids on the fly.
- This replaces the previous approach that relied on precomputed MgSO₄ phase lookup tables, which were:
  - Low resolution.
  - Memory intensive.

### Bug Fixes
- Numerous bug fixes addressing edge-case planetary configurations and plotting-related errors.

## [3.0.0] - 2024-05-01
**Author:** @Chang-Scott

This release implements a broad set of changes developed over the past year, with several major new scientific and modeling capabilities.

### Major Additions
- Added the ability to model an **NaCl ocean** using an in-development NaCl(aq) Equation of State from **SeaFreeze**.
- Introduced support for **speciated ocean chemistries** (*CustomSolutions*) using the **Frezchem** and **SUPCRT16** databases, adapted through the chemical modeling package **Reaktoro**.
- Enabled calculation of **chemical (metabolic) reaction affinities** up to **1 GPa** using the **SUPCRT16-organics** database via **Reaktoro**.
- Integrated **PyALMA3** (the Python implementation of **ALMA3**) to compute **tidal Love numbers**.

### Expanded Ocean Chemistry Support
- These updates extend ocean world modeling beyond the previously supported:
  - Seawater
  - MgSO₄
  - Pure H₂O  
- PlanetProfile can now explore a substantially broader and more realistic geochemical parameter space.

### Documentation and Usage
- Guidance on using the **CustomSolution** capability is available upon request.
- A dedicated **PlanetProfile tutorial webpage** is currently in development.

### Bug Fixes
- Fixed numerous incidental bugs across modeling and analysis workflows.