# Session State — 2026-05-23

## Branch: `genai`

## What Was Done (2026-05-21 to 2026-05-23 session — k₂/h₂ observables and three-body validation)

Took the existing toolbox from "Titan-only deployment with v2.1
caches" to **three-body validation against the Europa, Ganymede, and
Callisto v2.1 caches**. Six bug fixes in `PlanetProfile/Inference/`,
plus a generalizable Tb-prior probe and three production configs.

### Three-body validation outcome

End-to-end smoke MCMC at `n_eff=100` per body:

- **Europa** (Seawater, 7D): 4166 samples, 100% non-rejected, log_likelihoods [-5.4, -0.4]. R_core_km IQR [489, 587] / prior bounds [0, 780] — strong constraint.
- **Ganymede** (PureH2O, 8D, HP ices): 4232 samples, 100% non-rejected, log_likelihoods [-21.0, -13.2]. Tb_K IQR [252.7, 254.4] — strong constraint.
- **Callisto** (NaCl 100 ppt, 8D, HP ices): 4313 samples, 100% non-rejected, log_likelihoods [-26.9, -18.1]. R_core_km IQR [685, 776] / prior bounds [0, 1205] — strong constraint.

The high Re(k₂) values on Ganymede / Callisto (~0.5) vs the prior
target ~0.20 reflect a 5-7σ model-data tension driven by the
**placeholder priors**: Mazarico+ 2023 publishes k₂ values for Europa
only, and the Ganymede/Callisto configs were initialized with the
Europa values pending body-specific replacements. This is a
configuration matter to revisit at production `n_eff`, not a code
defect.

### Six bugs surfaced and fixed (all in `PlanetProfile/Inference/`)

1. **Per-layer rheology dispatch in `forward_models.py`** — was
   feeding Andrade rheology into the Fe core (and Elastic into the
   ocean instead of Newton), causing TidalPy `radial_solver` to fail
   at "layer 0" with step-size-floor errors on Europa-shaped
   structures. Now mirrors PP's
   `Params.Gravity.rheology_models` defaults (Fe→Elastic,
   ocean→Newton, ices/silicate→Andrade). Cross-validated against
   PP's `_run_tidalpy_backend` direct call: matches `k₂ = 0.27,
   h₂ = 1.19` for default Europa.

2. **Eager TidalPy import in `forward_models.py` + `NUMBA_CACHE_DIR`
   setup in `run_inference_cli.py`** — fixes a numba "no locator
   available" failure under sandboxed runtime environments. Imports
   are now resolved at module load time, before any sampling, with
   a writable cache directory established.

3. **pocoMC `posterior()` unpack swap in `mcmc_runner.py`** —
   pocoMC returns `(samples, weights, logl, logp)`; the runner was
   reading `(samples, logl, logp, _)`, so the
   `InferenceResult.log_likelihoods` field was storing importance
   weights and the actual log-likelihoods were dropped. **The
   `samples` array was always correct** — pre-fix Titan-era
   posteriors are unaffected for any analysis that uses `r.samples`
   or the post-MCMC recomputation pass. Anything that read
   `r.log_likelihoods` directly (logl histograms, BIC/DIC,
   importance-weighted reweighting) was wrong.

4. **v2.1 list-format cache support in CMR² dispatcher
   (`mcmc_runner.py`)** — the dispatcher only handled the legacy
   dict-format cache (`'grid_cache'` + `'grid_Tb_values'`). The v2.1
   list format produced by `cache_builder.build_phase_c1_cache`
   (`'Tb_K_grid'` + `'structures'` + `'transitions'`) silently fell
   through to the top-level dict, which has no `Mtot_kg` /
   `R_body_m` keys → NaN → -1e30 rejection of every sample.
   Posteriors against v2.1 caches were degenerate (prior-shaped).
   Pre-Phase-C1 caches in dict format were unaffected.

5. **Zero-thickness hydrosphere shells in `structure_derivation.py`**
   — PP caches duplicate radii at layer interfaces by design (TidalPy
   requires the boundary to appear twice). Those duplicates surfaced
   as `r_in == r_out` shells in `extract_hydrosphere_layers`, which
   `compute_cmr2` rejected with "Bad layer geometry". The
   extraction now skips zero-thickness shells (zero mass, zero
   inertia — correct for mass conservation). Surfaced as Callisto's
   posterior remaining 0/4096 non-rejected even after fix #4.

6. **`ocean_overrides` field in `inference_core.py` +
   `Texc_hr` label-zip in `cache_builder.py`** — schema-coverage
   fixes for the unified body-config / cache + Texc_hr extraction
   from `Planet.Magnetic` (was treating numpy-array-of-periods as a
   dict, raising TypeError silently caught and storing `None`).

### Framework additions

- **`probe_tb_prior.py` (new module)** — generalizable diagnostic
  for any future body cache. Surfaces integrator-failing prior
  regions before launching production MCMC.
- **k₂/h₂ wiring**: Re/Im k₂ + Re/Im h₂ as MCMC observables. Single
  TidalPy `radial_solver` call returns both k₂ and h₂.
- **Magnetic induction wiring** (forward path only): induction observables
  (`Ae_<label>_real`/`_imag` or `BiAmp_<label>`/`BiPhase_<label>_deg`).
  No config enables them yet — awaiting body-specific Ae values.
- **Three production configs updated** (J₂/C₂₂ dropped pending
  generalized gravity-coefficient forward model; Callisto
  `rho_sil_kgm3` lower bound widened 2200 → 1800 for porous
  interior; metadata descriptions updated inline).

### Deferred (tracked in user memory)

- **Generalized gravity-coefficient forward model** — produce J\_n,
  C\_nm, S\_nm to arbitrary degree/order from interior structure.
  Targets bodies with high-resolution gravity: Mars MRO 165×165,
  Moon GRAIL ~900×900, Enceladus, etc. PP currently treats J₂ /
  C₂₂ as observation inputs, not predictions.

- **Callisto CMR² Gao & Stevenson 2013 nonhydrostatic-rotator
  variants** — MCMC variants at 5 % and 10 % lower CMR² priors,
  reflecting nonhydrostatic correction. Cache is independent of
  the observable and need not be rebuilt.

- **Ganymede / Callisto body-specific Re/Im k₂ priors** —
  placeholder values from Mazarico+ 2023 (Europa) currently used.
  Replace with body-specific literature values when available;
  Clipper / JUICE will narrow these substantially.

- **Ice VI full lid+adiabat+TBL profile propagation** under
  `KALOUSOVA_CONVECTION` (Ice III and V already done; Ice VI
  currently uses uniform-T placeholder). Pending `Steps.nVbottom`,
  `nIceVILitho`, layer allocation in `IceLayers()`, and
  `IceVIConductSolid/Porous`.

### Validation hygiene

- **Source-tree edits this session**: 7 files in
  `PlanetProfile/Inference/` plus 3 production configs. Diagnostic
  scripts (`/tmp/claude-503/`) are not staged.
- **No tests broken** — the existing layered-blend / CMR² /
  per-phase-zeta / no-ocean-safeguard test suites continue to
  pass.

## What Was Done (2026-05-18 session)

### HPIceConvection early-return fix — layer_propagators

`HPIceConvection` in `LayerPropagators.py` previously fell through to
`ConvectionDeschampsSotin2001` for HP ices (III/V/VI) when
`KALOUSOVA_CONVECTION=False`. D-S 2001 has a hardcoded `Tupper_K=274`
bracket in `ThermalProfiles.py:311–315` that is appropriate for Ice Ih
but fails at the pressures HP ices inhabit, raising a ValueError at low
Tb. The fix inserts an early return that writes no-convection diagnostic
defaults (eLid=zb, Dconv=0, Ra=0, RaCrit=∞, meltFraction=0) and leaves
T_K on the melt curve as `SelfConsistentOceanLayer` set it — the correct
original PP behaviour for HP ices. Pre-existing Titan and Ganymede caches
were not corrupted because D-S 2001's phase-mismatch escape hatch
(`ThermalProfiles.py:296–305`) returned a conductive-profile fallback
without touching T_K, ρ, μ, K, or η (the fields MCMC consumes); only
diagnostic fields were silently wrong. The regression was visible on
Callisto: 6/9 Tb grid points crashed; the rebuilt cache hit all 9 in
21 min compared to the 74 min partial run that aborted.

### v2.1 transition-aware Tb grid — cache_builder

`cache_builder.py` gained `_bisect_transition` and a post-grid refinement
phase. After the regular Tb sweep, the builder walks adjacent grid-point
pairs and detects layer-set changes (different `region_phases` or
`n_layers`). Each detected transition is bisected to ε_T = 0.01 K, and
the converged anchor pair flanks the discontinuity by less than that
tolerance. The cache schema is bumped to `v2.1` with a `transitions`
metadata list `[{Tb_lo, Tb_hi, regions_lo, regions_hi}, ...]`. Old v2.0
caches load unchanged. Physical motivation: a Tb at which a new layer
first stabilises — HP ice III appearing, ocean surface forming, Ih
switching from convective to conductive — is a genuine discontinuity. The
0.01 K window means nearest-neighbour fallback at MCMC time is below any
posterior precision we could resolve.

### B-layered structure blend — forward_models

`forward_models.py::_blend_b_layered` replaces the legacy element-wise
blend that produced unphysical mush cells (ρ not matching phase identity)
when blending across a moving layer boundary. The new implementation
blends boundary radii linearly in `w`, then resamples each continuous
per-cell field (ρ, K, μ, η, T, P, bulk_visc) onto s0's intra-layer
normalised grid via `np.interp`; phase and discrete metadata are copied
from s0 (guaranteed identical to s1 across same-structure brackets). Per-
Tb scalars (CMR², D_iceIh_km, …) are linearly blended. Dispatch logic in
`apply_bottom_temperature`: matching `region_phases` → B-layered blend;
mismatched (cross-transition bracket) → nearest-neighbour. A precondition
check asserts every per-cell field has length equal to `len(r_m)` before
blending, surfacing padding bugs by field name rather than inside
`np.interp`. 12 unit tests in `tests/test_layered_blend.py` cover the
key invariants (boundary linearity, variable cell counts, no-mush, w=0/1
endpoints, transition dispatch, scalar/meta handling).

### P_MPa padding bug — cache_builder

The thin-layer padding pass in `build_single_structure` (lines 256–316)
extended all per-cell arrays via `np.interp` except `P_MPa`, which was
kept at its original length. Caches with any padded layer therefore had
`len(P_MPa) < len(r_m)` by the number of single-cell layers padded
(typically 4). The bug was latent — Test50's runner never blended P_MPa
cell-by-cell. The B-layered blend exercises P_MPa per cell, so the mismatch
surfaced immediately. Fixed by adding `P_MPa` to the padding interpolation,
identical to T_K. Zero field-length mismatches confirmed across rebuilt
Callisto, Europa, and Ganymede caches.

### ocean_overrides — cache_builder

A new optional `ocean_overrides` dict in the v2.1 BodyConfig JSON is
applied to `Planet.Ocean.<key>` after deepcopy of the body template and
before running PP. This enables composition switching and concentration
sweeps without maintaining N separate Python template files. The immediate
use case is Callisto: PP's MgSO₄ EOS tops out at P=800 MPa; Callisto's
deep ocean exceeds 200 MPa and triggers extrapolation regeneration on every
call, making MgSO₄ builds prohibitively slow. `callisto_nacl_andrade_8D.json`
uses `ocean_overrides` to switch to SeaFreeze's NaCl EOS (valid to 5 GPa).
The MgSO₄ config is kept as a deprecated reference; a future SeaFreeze
MgSO₄ extension can revisit.

### Inference toolkit README — docs

`PlanetProfile/Inference/README.md` (723 lines) written as the canonical
methodology reference for the Phase C1 workflow. Covers pocoMC algorithm
internals (preconditioned normalising flows, SMC tempering, n_effective
semantics, ESS termination), v2.1 cache schema and transition semantics,
B-layered blend derivation, ocean_overrides usage, and the full Phase C1
end-to-end sequence (config JSON → build cache → smoke test → production
MCMC).

### Phase C1 v2.1 caches — Callisto, Europa, Ganymede

Three production structure caches built with the v2.1 schema, transition
refinement, and P_MPa fix:

- `callisto_nacl_structure_grid.pkl`: 17 points, 1 transition (ocean
  appearance), ~40 min build time with NaCl 100 ppt.
- `europa_seawater_structure_grid.pkl`: 16 points, 1 transition, 53 s.
- `ganymede_pureh2o_structure_grid.pkl`: 23 points, 2 transitions, 55 s.

Smoke MCMC against the Callisto NaCl cache completed cleanly:
ESS=4096, acceptance=0.61, posterior appropriately prior-dominated for
1 observable in 8D (CMR² alone). Open future work: Gao & Stevenson 2013
argues that slow rotators may not be in hydrostatic equilibrium, placing
Callisto's CMR² up to 10% below the Anderson 2001 value. Inference runs at
CMR² nominal−5% and nominal−10% (same σ=0.0042) should be built; the
structure cache is observable-independent and need not be rebuilt.

---

## What Was Done (2026-05-13 session)

### Clathrate-thickness bug — narrow-scope fix

Diagnosed and corrected a factor-of-2 error in PP's spherical conduction integration constant `c1` (relative to Turcotte & Schubert 2002 §4.9 eq. 4.40). The bug caused `GetPbConduct` to integrate twice as deep as needed before reaching the target temperature, so user-specified `Bulk.clathMaxThick_m = 2 km` produced a 4 km clathrate layer.

- **Initial broad fix** (then reverted): corrected `c1` in `ConductiveTemperature` and `ConductiveTemperatureActual` directly. Confirmed clathrate fix at structural level (2 km) but full test suite (`python -m PlanetProfile.BuildTest`) regressed at Test 15: silicate T propagation diverges as `1/r` near `r=0` for solid silicate bodies (no inner core or `CONSTANT_INNER_DENSITY`), driving SilicateLayers mass-balance failure. The legacy halved-`c1` had been masking this for years.
- **Final narrow fix**: legacy `ConductiveTemperature` / `ConductiveTemperatureActual` retained with their halved `c1` (now clearly documented as a known discrepancy + pointer to the correct version). Added new `ConductiveTemperatureCorrect` with strict T&S 4.40 form. Pointed `GetPbConduct` (the only caller that needs the right answer for clathrate depth) to the new function. Test suite re-ran clean.
- **Silicate BC issue resolved (2026-05-14).** `SilRecursionSolid` and `SilRecursionPorous` now detect the solid-sphere case (no Fe core, no `CONSTANT_INNER_DENSITY`) and enforce the only physically admissible BC for T finite at r=0: at each downward step the inherited `qTop` is overridden with the consistency value `qTop = Htot·rTop/3`, making `c1 = 0` in `ConductiveTemperatureCorrect`. The closed-form profile is `T(r) = T_top + Htot·(rTop² − r²)/(6k)`. Shell bodies (Fe core or `CONSTANT_INNER_DENSITY`) retain the legacy `ConductiveTemperature` call. FIXME markers removed. Validated against BuildTest Tests 1–36 (all pass; Test 37 still fails with the pre-existing clathrate-underplate stub error, unrelated).

### Test50 MCMC upgraded 7D → 8D

`Test50_mcmc_andrade_noocean_yao2014.py`:
- Added `Tb_K` as 8th sampled parameter, uniform prior `[249.0, 250.965]` K. Physically marginalizes over solute-driven Ih-III-L triple-point depression up to ~15 ppt NaCl-equivalent.
- Structure cache pre-built on a 9-point Tb grid (Option 2 — grid + interpolate). `forward_model` linearly interpolates between bracketing grid points per sample, capturing eLid / TBL / PbI drift with Tb.
- All ice-phase viscosity priors broadened to `log10 η ∈ [10, 16]` to admit Petricca low-η regime across every ice phase.
- No-ocean safeguard: forward model rejects any Ih cell whose interpolated T crosses the linearized Ih-L melt curve `Tm(P) = 273.16 − 0.068·P_MPa` (with 0.1 K margin). Keeps the no-ocean assumption internally consistent across the sampled Tb band.

`PPTest50.py`:
- Default clathrate cap reduced from 10 km to 2 km.
- Default `Tb_K = TtripleIh_III_L_K − 0.2 K` (ε scaled to grid resolution; comment explains derivation `2·|dTm/dP|·ΔP_cell`).

### Constant rename

`Constants.triplePointT_K` (251.165 K) → `Constants.TtripleIh_III_L_K`. The original name was misleading; this is the Ice Ih-III-L triple point, not the water Ih-L-vapour point at 273.16 K. Deprecation alias preserved in `defineStructs.py`. 13 callers updated across PP source files and test configs (including PPTest21, which is normally protected by the `<28` rule but the user explicitly authorized inclusion).

### New diagnostic scripts (committable)

`PlanetProfile/Test/scripts/`:
- `probe_heat_flux.py` — top-down q profile through clathrate + Ih layers from a built grid cache.
- `probe_sil_moduli.py` — verifies silicate ρ/μ/K aren't contaminated by divergent T (used during the broad-fix audit to confirm Test50 MCMC was safe before we narrowed scope).

## What's Pending

### Immediate (this session)
1. **Full BuildTest re-run** — running in background with the narrow fix; awaiting completion to confirm no regression (Tests 1–14 pre-confirmed; Test 15 was the regression point and should now pass).
2. **Commit + push** to `genai` after BuildTest passes.

### Scheduled follow-up
3. **Silicate BC rework** (Task #11): COMPLETED 2026-05-14 — Option A (T_center-anchored, c1=0) implemented for solid-sphere bodies in `SilRecursionSolid` / `SilRecursionPorous`. Shell bodies preserved with legacy halved-c1.
4. **Launch Test50 MCMC** with the new 8D + grid setup. Compare posterior against Test48 (ocean-bearing) to test Petricca's no-ocean hypothesis.

### Toolkit continuation (carried from previous session)
5. Port Test48 `forward_model` to `mcmc_common.apply_arrhenius_ih` / `split_silicate_core` / `compute_per_phase_heating`. Currently only `make_plots` uses the toolkit.
6. Port Test50 plotting to `mcmc_plots`. Current Test50 has inline plots from initial port.
7. BodyConfig abstraction for Europa, Ganymede, Callisto.
8. PlanetProfileApp integration (Phase B–D): body selector, live progress, multi-run comparison.

### Kalousova HP-ice follow-ups (NEW 2026-05-14)
9. **Ice VI full thermal-profile propagation** — port the lid + adiabatic interior + lower TBL three-segment construction from `IceVConvectSolid` (lines ~764–839) to `IceVIConvectSolid`. Currently Ice VI uses the simplified uniform-T-along-melting-curve placeholder. Required scaffolding: `Planet.Steps.nVbottom`, `Planet.Steps.nIceVILitho`, layer allocation in `IceLayers()`, `IceVIConductSolid/Porous` in `IceConduction.py`.
10. **Ice VI porous variant** — `IceVIConvectPorous` is a stub that currently calls `IceVIConvectSolid` with a warning; needs a real porous implementation once #9 lands.
11. **Two-phase melt-fraction solver** — replace the placeholder `meltFractionIII/V/VI` values (0.01 for III/V, 0.5 for VI in top-level dispatchers; `phiPercolation` for in-ocean path) with a proper two-phase steady-state solver per Kalousova & Sotin 2018. `Constants.phiPercolation` (= 0.05) must remain unchanged — Kalousova's outflow equations are conditioned on it.

## Key Files Touched This Session

| File | Status |
|---|---|
| `PlanetProfile/Thermodynamics/ThermalProfiles/ThermalProfiles.py` | Added `ConductiveTemperatureCorrect`; legacy functions documented + retained; `GetPbConduct` rewired |
| `PlanetProfile/Thermodynamics/Geophysical.py` | Silicate BC fix: solid-sphere case now uses c1=0 form via `ConductiveTemperatureCorrect` with overridden qTop; shell case unchanged |
| `PlanetProfile/Utilities/defineStructs.py` | `triplePointT_K` → `TtripleIh_III_L_K` (with deprecation alias) |
| `PlanetProfile/Test/PPTest50.py` | 2 km clathrate, Tb default `TtripleIh_III_L_K − 0.2 K` |
| `PlanetProfile/Test/Test50_mcmc_andrade_noocean_yao2014.py` | 8D MCMC, grid + interp, no-ocean safeguard, broadened priors |
| `PlanetProfile/Test/PPTest{21,35,37,38,39,41,41NoClath,43,46_allice}.py` | `triplePointT_K` rename |
| `PlanetProfile/Test/Test42_mcmc_maxwell_ocean.py`, `PPTitanPureWaterNoOcean.py` | `triplePointT_K` rename (1 ref each) |
| `PlanetProfile/Test/scripts/probe_heat_flux.py`, `probe_sil_moduli.py` | New diagnostic scripts |
| `CHANGELOG.md` | New `[Unreleased]` entry |
| `MCMC_INFERENCE_GUIDE.md` | Test 50 section rewritten; staleness fixes |

## Cache Files

In `PlanetProfile/Test/mcmc_results/`:
- `titan_allice_yao2014_structure_grid.pkl` — current (post-narrow-fix; 2 km clathrate, 9-point Tb grid)
- `titan_allice_yao2014_structure_grid_q20mW.pkl` — pre-fix snapshot (4 km clathrate, kept for reference)

## Environment Notes

- Mamba env: `PPcl`
- Python: `/Users/svance/mamba/envs/PPcl/bin/python`
- All ad-hoc scripts in `PlanetProfile/Test/scripts/` (per project convention established this session). Mamba env directory is off-limits.
