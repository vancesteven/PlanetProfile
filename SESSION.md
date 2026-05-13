# Session State — 2026-05-13

## Branch: `genai`

## What Was Done (this session)

### Clathrate-thickness bug — narrow-scope fix

Diagnosed and corrected a factor-of-2 error in PP's spherical conduction integration constant `c1` (relative to Turcotte & Schubert 2002 §4.9 eq. 4.40). The bug caused `GetPbConduct` to integrate twice as deep as needed before reaching the target temperature, so user-specified `Bulk.clathMaxThick_m = 2 km` produced a 4 km clathrate layer.

- **Initial broad fix** (then reverted): corrected `c1` in `ConductiveTemperature` and `ConductiveTemperatureActual` directly. Confirmed clathrate fix at structural level (2 km) but full test suite (`python -m PlanetProfile.BuildTest`) regressed at Test 15: silicate T propagation diverges as `1/r` near `r=0` for solid silicate bodies (no inner core or `CONSTANT_INNER_DENSITY`), driving SilicateLayers mass-balance failure. The legacy halved-`c1` had been masking this for years.
- **Final narrow fix**: legacy `ConductiveTemperature` / `ConductiveTemperatureActual` retained with their halved `c1` (now clearly documented as a known discrepancy + pointer to the correct version). Added new `ConductiveTemperatureCorrect` with strict T&S 4.40 form. Pointed `GetPbConduct` (the only caller that needs the right answer for clathrate depth) to the new function. Test suite re-ran clean.
- **Silicate BC issue deferred** as a scheduled follow-up. FIXME markers placed in `Geophysical.py::SilRecursionSolid` and `SilRecursionPorous`.

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
3. **Silicate BC rework** (Task #11): rework `PropagateConductionProfilesSolid` / `...Porous` in `Geophysical.py` to use `T_center`-anchored inward propagation with closed-form `T(r) = T_top + Htot·(rTop² − r²)/(6k)`, or impose `q_top = Htot·rTop/3` consistency at the ice/silicate interface. Until done, silicate T profiles are physically wrong but stable (legacy halved-c1).
4. **Launch Test50 MCMC** with the new 8D + grid setup. Compare posterior against Test48 (ocean-bearing) to test Petricca's no-ocean hypothesis.

### Toolkit continuation (carried from previous session)
5. Port Test48 `forward_model` to `mcmc_common.apply_arrhenius_ih` / `split_silicate_core` / `compute_per_phase_heating`. Currently only `make_plots` uses the toolkit.
6. Port Test50 plotting to `mcmc_plots`. Current Test50 has inline plots from initial port.
7. BodyConfig abstraction for Europa, Ganymede, Callisto.
8. PlanetProfileApp integration (Phase B–D): body selector, live progress, multi-run comparison.

## Key Files Touched This Session

| File | Status |
|---|---|
| `PlanetProfile/Thermodynamics/ThermalProfiles/ThermalProfiles.py` | Added `ConductiveTemperatureCorrect`; legacy functions documented + retained; `GetPbConduct` rewired |
| `PlanetProfile/Thermodynamics/Geophysical.py` | FIXME markers at silicate conduction loops |
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
