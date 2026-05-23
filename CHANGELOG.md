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

### Inference / MCMC — k₂/h₂ observables and three-body validation

Six bug fixes and three framework additions, all confined to
`PlanetProfile/Inference/`. End-to-end smoke MCMCs (`n_eff=100`)
against the v2.1 caches now produce **real Bayesian posteriors for
Europa, Ganymede, and Callisto** (4166 / 4232 / 4313 samples, 100 %
non-rejected, all observables finite), with `R_core_km` and
`log10_eta_Ih` showing genuine constraint by the data. The framework
generalizes the existing Titan deployment to ocean worlds with
ice + ocean + (HP ice) + silicate + Fe core layouts.

#### Bug fixes

- **`forward_models.py` — per-layer rheology mapping.** The Andrade
  branch in `forward_model_k2_flexible` previously assigned `Andrade`
  to every layer except `('0', 'Clath')` and `Elastic` to those two.
  PP's own `Params.Gravity.rheology_models` (defined in
  `Gravity/defaultConfigGravity.py:26`) maps `'Fe' → 'elastic'`,
  `'0' → 'newton'`, `'Clath' → 'elastic'`, ices and silicate to
  `'andrade'`. Treating the Fe core as Andrade with finite viscosity
  produced a near-singular complex shear modulus that TidalPy's
  `radial_solver` could not traverse, surfacing as
  `RadialSolver.ShootingMethod:: Integration problem at layer 0;
  solution 0: Required step size is less than spacing between
  numbers` and 100 % NaN k₂ for Europa during smoke MCMC. The fix
  routes Fe + Clath through `Elastic`, ocean (`'0'`) through
  `Newton`, ices/silicate through `Andrade` — matching PP's
  reference path. Confirmed against PP's direct
  `_run_tidalpy_backend` on default Europa: `k₂ = 0.2657 − 0.0008j`,
  `h₂ = 1.1896 − 0.0029j`, matching Mazarico+ 2023 within bounds.

- **`forward_models.py` — eager TidalPy import.** `forward_model_k2_flexible`
  used to lazy-import `TidalPy.rheology` and `TidalPy.RadialSolver`
  inside the function body. In sandboxed runtime environments where
  the default numba cache path is unwritable, the lazy import
  triggered `numba.core.caching.RuntimeError: cannot cache function
  'off': no locator available for file ...partial_melt/melting_models.py`
  on the first sample of the MCMC. Imports are now at module top so
  numba cache state is established once before any sampling. Combined
  with `run_inference_cli.py` setting `NUMBA_CACHE_DIR` (below) this
  is fully resolved.

- **`mcmc_runner.py:513` — pocoMC `posterior()` unpack swap.**
  Previous code: `samples, log_likes, log_post, _ = sampler.posterior()`.
  pocoMC's actual return is `(samples, weights, logl, logp)`, so
  `log_likes` was storing **importance weights** (≈ 1/N after pocoMC's
  trim) and `log_post` was storing the actual log-likelihoods. The
  `samples` array was always correct (the MCMC drew from the right
  posterior) — only the `r.log_likelihoods` field on `InferenceResult`
  was mislabeled. **Existing pre-fix Titan-era result pickles are
  scientifically unaffected**: corner plots, posterior medians, IQRs,
  and predicted-observable summaries all derive from `r.samples` or
  from the post-MCMC recomputation pass, neither of which touched the
  buggy unpack. Anything that read `r.log_likelihoods` directly
  (logl histograms, importance-weighted reweighting, BIC/DIC) would
  have been wrong.

- **`mcmc_runner.py` — v2.1 list-format cache in CMR² dispatcher.**
  The CMR² mass-conservation re-derivation block only handled the
  legacy dict-format cache (top-level `'grid_cache'` + `'grid_Tb_values'`)
  for picking the per-sample structure. The v2.1 list format produced
  by `cache_builder.build_phase_c1_cache` (top-level `'Tb_K_grid'` +
  `'structures'` + `'transitions'`) silently fell through to
  `struct_for_cmr2 = structure_data`, where the top-level dict has
  no `Mtot_kg` / `R_body_m` keys. `np.isfinite(np.nan)` → False → every
  sample rejected at `-1e30`. Posteriors against v2.1 caches were
  the prior, "shape-decorated by random rejection." The dispatcher
  now handles both layouts. Pre-Phase-C1 analyses (Test48, Test50,
  Titan) are unaffected because they used the dict format that the
  legacy branch already handled correctly.

- **`structure_derivation.py` — skip zero-thickness hydrosphere shells
  in `extract_hydrosphere_layers`.** PP caches legitimately duplicate
  radii at layer interfaces (TidalPy's `radial_solver` requires the
  boundary radius to appear twice — once for the layer below, once
  for the layer above). Those duplicates surface in the per-cell
  hydrosphere extraction as cells with `r_in == r_out`, which
  `compute_cmr2` then rejects with `Bad layer geometry`. A
  zero-thickness shell contributes zero mass and zero inertia, so
  dropping it is correct for the mass-conservation derivation.
  Surfaced as Callisto's smoke MCMC at 0/4096 non-rejected even
  after the v2.1-dispatcher fix, traced via direct
  `MCMCRunner.log_likelihood_fn` invocation.

- **`inference_core.py` — `ocean_overrides` field on `InferenceConfig`.**
  Added `ocean_overrides: Dict[str, Any] = field(default_factory=dict)`.
  This is a cache-build-only knob (consumed by
  `cache_builder.build_phase_c1_cache` to swap `Planet.Ocean.comp` /
  `Planet.Ocean.wOcean_ppt`, e.g. NaCl 100 ppt for Callisto). The
  MCMC runner ignores it. Required so a single body config file can
  drive both stages without `InferenceConfig.from_json` rejecting
  the unknown key.

- **`run_inference_cli.py` — set `NUMBA_CACHE_DIR` before TidalPy
  import.** Defaults to `${TMPDIR}/pp_numba_cache` if the env var is
  unset. Avoids the numba "no locator available" failure under
  sandboxed runtimes (see eager-import fix above).

- **`cache_builder.py` — `Texc_hr` label-zip via `ExcitationsList`.**
  `Planet.Magnetic.Texc_hr` is a numpy array of periods (canonical
  labels live in the body-keyed lookup table in
  `MagneticInduction.Moments.ExcitationsList`). The cache_builder
  was treating it as a dict (`dict(mag.Texc_hr).items()`), raising
  TypeError silently caught by the surrounding `except` clause and
  storing `Texc_hr = None`. Now matches each cached period to its
  closest labelled period via `ExcitationsList()` and builds the
  `{label: period_hr}` dict the induction forward model consumes.
  PureH2O bodies (e.g. the default-config Ganymede) skip
  `MagneticInduction` entirely (PP `Main.py:352` short-circuits when
  `Ocean.comp == 'PureH2O'`); for those, `Texc_hr` legitimately
  remains `None` — documented as expected, not a bug.

#### Framework additions

- **`probe_tb_prior.py` (new module)** — generalizable diagnostic
  for surfacing integrator-failing prior regions. Given an
  `InferenceConfig`, evaluates `forward_model_k2_flexible` at each
  Tb_K grid point with `n_fiducials` theta draws (default 5: prior
  median + four random draws), reports per-Tb success rate, and
  recommends a tightened `[Tb_min, Tb_max]` bound covering only
  grid points whose success rate meets the threshold. Use
  pre-flight on any new body cache before launching production MCMC
  to avoid silent posterior truncation. Headless CLI:
  `python -m PlanetProfile.Inference.probe_tb_prior --config
  <cfg.json>`.

- **k₂/h₂ observables (Re/Im k₂, Re/Im h₂) wired through the
  forward model and dispatcher.** `forward_model_k2_flexible` returns
  a 5-tuple `(Re_k2, Im_k2, Re_h2, Im_h2, perPhase_W)` from a single
  TidalPy `radial_solver` call — h₂ is read from `result.h` on the
  same RadialSolver result, not a separate solve. Dispatcher in
  `mcmc_runner.py` handles `Re_k2`/`Im_k2`/`abs_Im_k2`/`k2`/`Re_h2`/
  `Im_h2`/`abs_Im_h2`/`h2` keys. Default observable values for the
  three production configs use Mazarico+ 2023 Europa values
  (`k₂ ≈ 0.20 ± 0.05`, `h₂ ≈ +1.12 ± 0.10`, `Im k₂ = Im h₂ = 0 ±
  0.05`/`0.1`); these are **placeholders** for Ganymede and
  Callisto (Mazarico+ 2023 publishes only Europa values).
  Body-specific updates pending Clipper / JUICE measurements and
  any other relevant Ganymede / Callisto literature.

- **Magnetic induction observables wired (forward path only).** Cache builder
  extracts `sigma_Sm` (per-cell), `rSigChange_m` and
  `sigmaLayers_Sm` (compact-layer representation), and `Texc_hr`
  (labelled period dict) from `Planet.Magnetic`. Dispatcher in
  `mcmc_runner.py` accepts `Ae_<label>_real`/`Ae_<label>_imag`
  (Re/Im default) or `BiAmp_<label>`/`BiPhase_<label>_deg`
  (amp/phase legacy) observables; both pull from the same
  `forward_model_induction` call. **No production config currently
  enables induction observables** — the wiring is in place but
  awaiting body-specific Ae values.

- **J₂ / C₂₂ deferred.** Initial wiring assumed PP exposed J₂
  and C₂₂ as forward-model outputs; in fact `Planet.Bulk.J₂` and
  `Planet.Bulk.C₂₂` are observation *inputs* set in body config /
  test files, not predicted quantities. The cache_builder
  silently NaN'd these via `getattr(Planet, "J2", None)`. Dropped
  from all three production configs pending a generalized
  gravity-coefficient forward model (J\_n, C\_nm, S\_nm to
  arbitrary degree/order) for bodies with high-resolution gravity
  (Mars MRO 165×165, Moon GRAIL ~900×900, Enceladus, …). Tracked
  in user memory as future work.

- **Callisto `rho_sil_kgm3` lower bound 2200 → 1800.** Callisto's
  silicate is porous (φ\_vac ≈ 0.91 from PP MoI matching), giving
  cached `rho_sil ≈ 2025 kg/m³`. The Europa-tuned bound `[2200,
  3500]` rejected every MCMC sample at the mass-conservation gate,
  producing a degenerate prior-shaped posterior. Widened lower
  bound to 1800 for Callisto only; Europa and Ganymede continue
  to use `[2200, 3500]`. Documented inline in
  `callisto_nacl_andrade_8D.json::core_bounds_rationale`.

#### Validation summary

| Body | samples | non-rejected | Re(k₂) median (target ≈ 0.20 †) | CMR² median (target) | R_core_km IQR / prior bounds |
|---|---|---|---|---|---|
| Europa | 4166 | 100 % | 0.26 | 0.346 / 0.355 | [489, 587] / [0, 780] |
| Ganymede | 4232 | 100 % | 0.46 | 0.312 / 0.312 | [454, 759] / [0, 1316] |
| Callisto | 4313 | 100 % | 0.50 | 0.341 / 0.355 | [685, 776] / [0, 1205] |

† Re(k₂) target uses the Mazarico+ 2023 Europa value as a
placeholder for Ganymede and Callisto. The 5–7σ tension on those
two reflects priors-vs-data, not a code defect; revisit at
production `n_eff` after replacing placeholders with body-specific
Re/Im k₂ priors.

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