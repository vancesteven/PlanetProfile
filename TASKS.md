# PlanetProfile Inference Framework - Task Status

**Last Updated:** 2026-05-02  
**Current Phase:** Maxwell Tb_K Grid Cache Complete ✓ — Next: induction likelihood

---

## Phase Status Overview

| Phase | Status | Notes |
|-------|--------|-------|
| **Phase 1** | ✅ Complete | GUI skeleton, cache utilities, validation |
| **Phase 4** | ✅ Complete | Core inference engine, forward models, structure cache |
| **Phase 5** | ✅ Complete | MCMCRunner (54% acceptance, 0.7 min), per-phase zeta |
| **Phase 2** | ✅ Complete | Parameter configuration UI (extensible, per-body) |
| **Phase 3** | ✅ Complete | Corner plot, k2 scatter, heating distribution, tab rename |
| **Grid Cache** | ✅ Complete | Test42 Tb_K grid: 121/121 pts in 255–270 K prior range |
| **Phase 7** | ⏸️ Pending | GUI integration (subprocess, progress, results) |
| **Phase 8** | ⏸️ Pending | Export handlers, replotting CLI |
| **Phase 6** | ⏸️ Deferred | SBI runner (neural posterior estimation) |

---

## Completed This Session: Maxwell Grid Cache + TidalPy Thin-Layer Fix

### 1. `build_structure_grid()` in `structure_cache.py` ✓

Builds a `dict[float → structure_data]` cache keyed by `Tb_K`, covering the full
prior range for Test42.  Key features:
- Incremental save after each point — a killed process loses no progress
- `force_rebuild=False` skips already-cached points on re-run
- `rheology` parameter selects the `configParams.Gravity.rheology_models` dict

### 2. TidalPy thin-layer fix ✓

**Root cause:** With `POROUS_ROCK=True`, a porosity phase transition occupying a single
radial point in the full PP model produces a 1-point layer in `Planet.Reduced`.  TidalPy
requires ≥ 5 slices per layer; `Gravity.py`'s existing padding only handles `n_pts ≥ 2`.

**Fix — two-part:**

a) `_expand_thin_reduced_layers(Planet)` in `structure_cache.py`: after `PlanetProfile()`
   returns, expand any layer with fewer than 5 points to exactly 5 by linear interpolation
   of radii and constant-fill of all other arrays (ρ, phase, seismic, η).

b) Inner try/except around the `PlanetProfile()` call catches the TidalPy
   "Layer X has 2 slices when at least 5 are required" error that fires from
   PlanetProfile's *internal* `SetupGravity` call.  Because `Planet` is mutated in-place,
   `Planet.Reduced` is fully populated at the point of failure; we expand and retry
   `SetupGravity` directly.

**Result:** 50 additional Tb_K grid points now succeed that previously failed.

### 3. `PPTest42.py` physical parameter updates ✓

| Parameter | Before | After | Reason |
|-----------|--------|-------|--------|
| Tb_K MCMC prior | [255, 275] K | [255, 270] K | > 273.15 K is unphysical (above 1-bar ice Ih melting point) |
| `clathMaxThick_m` | 10 km | 5 km | 10 km clathrate cap exceeds thin ice shell at high Tb_K, causing PlanetProfile to abort |
| `POROUS_ROCK` | False | True | Required for correct silicate MoI matching |
| `CuncertaintyLower/Upper` | ±0.01 | ±0.05 | ±0.01 window hit zero-init C_kgm2 entries, breaking iValid lookup |

### 4. Grid cache status ✓

- Total cached: **157/189** grid points (full 251–275 K sweep)
- MCMC prior range coverage: **121/121** points in 255–270 K, spacing ~0.13 K
- Remaining gaps are physical failures outside the prior:
  - Tb_K = 251.0 K: no valid ice-ocean phase transition at low pressure
  - Tb_K ≥ 271.2 K: clathrate transition pressure falls below ice Ih melting curve
    (these sit outside the revised 255–270 K prior — not a problem for MCMC)

---

## Immediate Next Task: Wire Induction into MCMC Likelihood

**Status:** Infrastructure present, chi² terms not yet connected.

`forward_model_induction()` exists in `forward_models.py` and can compute complex k2
and induction responses from a cached structure + rheology sample.  What's missing:

1. **In `forward_model_k2_flexible` / `compute_log_likelihood`:** add Im(k2) as a
   second observable alongside Re(k2).  Test42 constraints (Petricca et al. 2025):
   - Re(k2) = 0.608 ± 0.048
   - |Im(k2)| = 0.135 ± 0.035

2. **Grid lookup for Tb_K:** `forward_model_k2_flexible` currently assumes a single
   cached structure.  For Test42 with `Tb_K` as a free parameter, it must look up the
   nearest grid point:
   ```python
   tb_idx = np.argmin(np.abs(np.array(sorted(grid_cache.keys())) - Tb_K))
   structure = grid_cache[sorted(grid_cache.keys())[tb_idx]]
   ```

3. **MCMCRunner:** when `'Tb_K'` is in `param_names`, pass `grid_cache` (dict) instead
   of `structure_data` (single dict) to `forward_model_k2_flexible`.

4. **Log-likelihood function signature update:** accept optional `obs_Im_k2` and
   `obs_Im_k2_sigma` in addition to existing `obs_Re_k2` / `obs_Re_k2_sigma`.

**Files to edit:** `forward_models.py`, `mcmc_runner.py`, `inference_core.py`

---

## Pending: Phase 7 — GUI Integration

**Prerequisites:** Induction likelihood + grid-cache dispatch working end-to-end.

- Run control buttons (Run / Abort / Force Re-run)
- Subprocess execution via `run_inference_cli.py`
- Real-time progress display (poll JSON file written by MCMCRunner)
- Results display (tabbed: Corner | Trace | Diagnostics | k2 | Heating)
- Error handling with traceback display

---

## Pending: Phase 8 — Export & Replotting

**Prerequisites:** Phase 7 complete.

- Export handlers (PKL, HDF5, CSV, PNG)
- Standalone replotting CLI (`replot_inference.py`)
- Download buttons in GUI

---

## Deferred: Phase 6 — SBI Runner

MCMC sufficient for v1.  SBI adds ~350 lines + torch/sbi dependencies; defer until
MCMC workflow is end-to-end in the GUI.

---

## Known Issues

| Issue | Impact | Status |
|-------|--------|--------|
| `PPTest*NoClath.py` files incomplete | Blocks Test41 no-clath cache gen | Identified, fix not yet applied |
| Checkpoint resumption untested | Long runs only | Implementation exists |
| No multi-chain R-hat | Diagnostic accuracy | pocoMC single-chain; R-hat set to 1.0 |

---

## File Inventory

### Completed ✅

| File | Lines | Role |
|------|-------|------|
| `PlanetProfile/Inference/inference_core.py` | ~400 | Orchestration, config loading |
| `PlanetProfile/Inference/forward_models.py` | ~500 | k2/induction forward models |
| `PlanetProfile/Inference/structure_cache.py` | ~730 | Single + grid cache, thin-layer fix |
| `PlanetProfile/Inference/mcmc_runner.py` | ~380 | pocoMC wrapper, per-phase zeta |
| `PlanetProfile/Inference/parameter_registry.py` | ~350 | Prior definitions, param metadata |
| `PlanetProfile/Inference/run_inference_cli.py` | ~250 | CLI entry point |
| `PlanetProfile/Inference/prepare_structure_variants.py` | ~200 | Cache generation helpers |
| `PlanetProfileApp/pages/Inference.py` | ~950 | GUI page (all phases) |
| `PlanetProfile/Test/PPTest42.py` | ~145 | Maxwell ocean-bearing Titan config |

### Pending ⏸️

| File | Est. lines | Phase |
|------|-----------|-------|
| `PlanetProfile/Inference/export_results.py` | ~200 | 8 |
| `PlanetProfile/Inference/replot_inference.py` | ~300 | 8 |
| `PlanetProfileApp/pages/Inference.py` (extend) | ~500 | 7 |

---

## Session State

**Branch:** genai  
**Cache location:** `titan_cache/titan_maxwell_grid_cache.pkl`  
**Grid coverage:** 121 pts, 255.09–270.91 K, spacing ~0.13 K — ready for MCMC
