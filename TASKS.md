# PlanetProfile Inference Framework - Task Status

**Last Updated:** 2026-05-02  
**Current Phase:** PPTest45 Hybrid Hydrosphere MCMC — 6D run complete, plots verified

---

## Phase Status Overview

| Phase | Status | Notes |
|-------|--------|-------|
| **Phase 1** | ✅ Complete | GUI skeleton, cache utilities, validation |
| **Phase 4** | ✅ Complete | Core inference engine, forward models, structure cache |
| **Phase 5** | ✅ Complete | MCMCRunner (54% acceptance, 0.7 min), per-phase zeta |
| **Phase 2** | ✅ Complete | Parameter configuration UI (extensible, per-body) |
| **Phase 3** | ✅ Complete | Corner plot, k2 scatter, heating distribution, tab rename |
| **Grid Cache (Test42)** | ✅ Complete | 121/121 pts in 255–270 K prior range |
| **PPTest45 Hybrid Cache** | ✅ Complete | 72 pts (8×9), 255–270 K × 300–700 km |
| **PPTest45 MCMC harness** | ✅ Complete | 6D run done; Mtot≈1.345e23 kg, ρ_sil≈2570 kg/m³ posterior |
| **PPTest46 (Europa)** | ⏸️ Pending | Naming fixed; cache + MCMC script not yet written |
| **Phase 7** | ⏸️ Pending | GUI integration |
| **Phase 8** | ⏸️ Pending | Export handlers, replotting CLI |
| **Phase 6** | ⏸️ Deferred | SBI runner |

---

## Current Priority: PPTest45 Complete — Evaluate Results / Next Steps

## PPTest45 Hybrid Hydrosphere Non-Self-Consistent MCMC (DONE)

### Motivation

PPTest42 (`MCMCRunner` + Maxwell grid) treats Tb_K as free but still uses PlanetProfile's
MoI closure to place the silicate boundary — the structure is self-consistent.  PPTest45
relaxes this: the hydrosphere PT is self-consistent in Tb_K, but total hydrosphere
thickness D_hydro_km is a free MCMC parameter.  `Mtot_kg` and `CMR2mean` become model
*outputs* compared against observations, not closure targets.

### PPTest45 harness: `Test45_mcmc_maxwell_hybrid_hydro.py`

**Location:** `PlanetProfile/Test/Test45_mcmc_maxwell_hybrid_hydro.py`

**Parameter space (5D):**
| Parameter | Prior | Range |
|-----------|-------|-------|
| `log10_eta_Ih` | uniform | [12, 16] |
| `log10_eta_HP` | uniform | [10, 18] |
| `log10_eta_sil` | uniform | [12, 22] |
| `Tb_K` | uniform | [255, 270] K |
| `D_hydro_km` | uniform | [300, 700] km |

**Observational constraints:**
- Re(k2) = 0.608 ± 0.048
- |Im(k2)| = 0.135 ± 0.035
- CMR2 = 0.343 ± 0.001 (Petricca et al. 2025)
- Mtot_kg = 1.3452e23 ± 2×10²⁰ kg (Jacobson 2006)

**Grid:** 8 × 9 = 72 points (Tb_K × D_hydro_km).  Auto-built by harness if missing.

**To run:**
```bash
conda activate ./venvPP
python PlanetProfile/Test/Test45_mcmc_maxwell_hybrid_hydro.py
# To force grid rebuild:
python PlanetProfile/Test/Test45_mcmc_maxwell_hybrid_hydro.py --rebuild-grid
```

**Output:** `PlanetProfile/Test/mcmc_results/hybrid_hydro_mcmc.pkl` + 4 plots.

### Proof-of-concept grid (2×2, from SESSION_STATE.md)
- Tb_K=251.2, D_hydro_km=400: Mtot_kg=1.835e23, CMR2=0.4544
- Tb_K=251.2, D_hydro_km=500: Mtot_kg=1.727e23, CMR2=0.4178
- D_hydro_km variation changes both Mtot and CMR2 at fixed Tb_K ✓

---

## Completed This Session

### PPTest46 naming fixes
- `PPTest46.py` docstring: "PPTest45" → "PPTest46"
- `PPTest46.py` line 46: `PlanetStruct('Test45')` → `PlanetStruct('Test46')`
- `parameter_registry.py`: `andrade_europa` and `andrade_europa_nocore` presets → `PPTest46`
- `prepare_structure_variants.py`: added `elif 'PPTest46'` branch for Europa

### Prior session: Maxwell Grid Cache + TidalPy Thin-Layer Fix (Test42)

#### `build_structure_grid()` in `structure_cache.py` ✓
Builds `dict[float → structure_data]` keyed by Tb_K, incremental save per point.

#### TidalPy thin-layer fix ✓
`_expand_thin_reduced_layers(Planet)` expands any <5-point layer to 5; inner try/except
retries `SetupGravity` after expansion.  50 additional grid points now succeed.

#### `PPTest42.py` physical updates ✓
| Parameter | Before | After | Reason |
|-----------|--------|-------|--------|
| Tb_K prior | [255, 275] K | [255, 270] K | > 273.15 K unphysical |
| `clathMaxThick_m` | 10 km | 5 km | prevents overflow at thin-ice end |

#### Grid cache status ✓
- **121/121** points in 255–270 K MCMC prior range (spacing ~0.13 K)

---

## Pending: PPTest46 Europa Andrade MCMC

**Status:** Naming fixed; structure cache and harness not yet written.

Structure is a single point (no HP ice, fixed Tb_K ≈ 271 K), so cache generation is
fast.  Script modelled on `Test41_mcmc_andrade_no_ocean.py`.

To generate cache:
```bash
python -m PlanetProfile.Inference.prepare_structure_variants \
    --test-module PlanetProfile.Test.PPTest46 \
    --output-dir cache/ --rheology andrade
```

---

## Pending: Phase 7 — GUI Integration

**Prerequisites:** PPTest45 MCMC confirmed working end-to-end.

- Run control buttons (Run / Abort / Force Re-run)
- Subprocess execution via `run_inference_cli.py`
- Real-time progress display
- Results display (Corner | Trace | Diagnostics | k2 | Heating)

---

## Pending: Phase 8 — Export & Replotting

**Prerequisites:** Phase 7 complete.

---

## Deferred: Phase 6 — SBI Runner

---

## Known Issues

| Issue | Impact | Status |
|-------|--------|--------|
| `PPTest*NoClath.py` files incomplete | Blocks Test41 no-clath cache gen | Identified |
| Checkpoint resumption untested | Long runs only | Implementation exists |
| No multi-chain R-hat | Diagnostic accuracy | pocoMC single-chain; R-hat=1.0 |

---

## File Inventory

### Completed ✅

| File | Role |
|------|------|
| `PlanetProfile/Inference/inference_core.py` | Orchestration, config loading |
| `PlanetProfile/Inference/forward_models.py` | k2/induction forward models + hook system |
| `PlanetProfile/Inference/structure_cache.py` | Single + grid cache, thin-layer fix |
| `PlanetProfile/Inference/hybrid_structure_cache.py` | Tb_K × D_hydro_km hybrid grid |
| `PlanetProfile/Inference/mcmc_runner.py` | pocoMC wrapper, grid cache dispatch |
| `PlanetProfile/Inference/parameter_registry.py` | Prior definitions, param metadata |
| `PlanetProfile/Inference/prepare_structure_variants.py` | Cache generation CLI |
| `PlanetProfile/Test/PPTest42.py` | Maxwell ocean Titan config |
| `PlanetProfile/Test/PPTest45.py` | Hybrid hydrosphere Titan (deep-copy of Test42) |
| `PlanetProfile/Test/PPTest46.py` | Europa Andrade MCMC base (Variant B, Fe-FeS) |
| `PlanetProfile/Test/Test42_mcmc_maxwell_ocean.py` | Maxwell ocean MCMC harness |
| `PlanetProfile/Test/Test45_mcmc_maxwell_hybrid_hydro.py` | Hybrid hydrosphere MCMC harness |
| `PlanetProfileApp/pages/Inference.py` | GUI page (all phases) |

### Pending ⏸️

| File | Phase |
|------|-------|
| `Test46_mcmc_andrade_europa.py` | Europa Andrade inference |
| `PlanetProfile/Inference/export_results.py` | Phase 8 |
| `PlanetProfileApp/pages/Inference.py` (extend) | Phase 7 |

---

## Session State

**Branch:** genai  
**Test42 cache:** `titan_cache/titan_maxwell_grid_cache.pkl` — 121 pts, 255–270 K  
**Test45 cache:** `PlanetProfile/Test/mcmc_results/titan_maxwell_hybrid_hydro_grid_cache.pkl` — to be built on first run
