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