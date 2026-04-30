# PlanetProfile Inference Framework - Task Status

**Last Updated:** 2026-04-29  
**Current Phase:** Phase 5 Implementation Complete, Validation Pending

---

## Phase Status Overview

| Phase | Status | Lines | Effort | Notes |
|-------|--------|-------|--------|-------|
| **Phase 1** | ✅ Complete | ~750 | 1 week | GUI skeleton, cache utilities, validation |
| **Phase 4** | ✅ Complete | ~2,369 | 2 weeks | Core inference engine, forward models, structure cache |
| **Phase 5** | ⚠️ Implementation Complete, Validation Blocked | ~950 | 1 week | MCMCRunner, dual structure caches, helper CLI |
| **Phase 2** | ⏸️ Pending | ~800 | 1 week | Parameter configuration UI |
| **Phase 3** | ⏸️ Pending | ~500 | 1 week | Visualization components (corner, trace, k2, heating) |
| **Phase 7** | ⏸️ Pending | ~500 | 1 week | GUI integration (subprocess, progress, results) |
| **Phase 8** | ⏸️ Pending | ~500 | 1 week | Export handlers, replotting CLI |
| **Phase 6** | ⏸️ Deferred | ~350 | Optional | SBI runner (neural posterior estimation) |

**Total Implemented:** ~4,069 lines  
**Remaining:** ~2,300 lines (Phases 2, 3, 7, 8)

---

## Current Status: Phase 5 Validation Blocked

### What Was Completed (Phase 5)

**1. MCMCRunner Class** (~350 lines)
- File: `PlanetProfile/Inference/mcmc_runner.py`
- pocoMC sampler integration
- Automatic rheology detection (Andrade/Maxwell)
- Progress callbacks for GUI
- Checkpoint save/load
- Convergence diagnostics (ESS, R-hat, acceptance rate)
- Heating recomputation on posterior subset

**2. Helper CLI** (~200 lines)
- File: `PlanetProfile/Inference/prepare_structure_variants.py`
- Generates both clathrate and no-clathrate structure caches
- Automatic module name inference (PPTest41 → PPTest41NoClath)
- Verbose logging and size reporting

**3. No-Clathrate Test Variants** (~155 lines)
- Files: `PlanetProfile/Test/PPTest41NoClath.py`, `PlanetProfile/Test/PPTest3NoClath.py`
- **Status:** ❌ INCOMPLETE - Files exist but are not faithful copies of parent tests
- **Issue:** Missing variables (andradExponent, andrade_zeta, etaMeltKalousova_Pas, etc.)
- **Cause:** Initial implementation was minimal, not a full copy with targeted changes

**4. API Integration**
- Updated: `PlanetProfile/Inference/__init__.py`
- Lazy imports for MCMCRunner/SBIRunner
- Wired into `inference_core.run_inference()`

### Blocking Issue

**PPTest41NoClath.py and PPTest3NoClath.py are incomplete:**

The current files set only:
- `Planet.Do.CLATHRATE = False`
- `Planet.Bulk.clathMaxThick_m = 0.0`
- `Planet.Steps.nClath = 0`

But they're missing ALL other variables from parent tests:
- `Planet.Gravity.andradExponent`
- `Planet.Gravity.andrade_zeta` (full dict)
- `Planet.Ocean.etaMeltKalousova_Pas`
- `Planet.Do.DO_SELF_CONSISTENT_HTIDAL`
- `Planet.Ocean.HtidalIce_Wm3`
- All other bulk, ocean, silicate, core, seismic, magnetic settings

**Error when running structure generation:**
```
'<=' not supported between instances of 'int' and 'NoneType'
```

This occurs because TidalPy/gravity calculations expect these variables to exist.

### Fix Required (NOT YET APPLIED)

**PPTest41NoClath.py diff:**
```diff
--- Copy ALL 137 lines from PPTest41.py
+++ Change only:
    - Line 2: Docstring → "WITHOUT clathrate cap: MCMC exploration variant"
    - Line 36: PlanetStruct('Test41NoClath')  # was 'Test41'
    - Line 62: Planet.Do.CLATHRATE = False  # was True
    - Line 64: Planet.Bulk.clathMaxThick_m = 0.0  # was 10e3
    - Line 65: Planet.Steps.nClath = 0  # was 30
```

**PPTest3NoClath.py diff:**
```diff
--- Copy ALL 71 lines from PPTest3.py
+++ Change only:
    - Line 2: Docstring → "WITHOUT clathrate underplate layer"
    - Line 9: PlanetStruct('Test3NoClath')  # was 'Test3'
    - Line 22: Planet.Do.CLATHRATE = False  # was True
    - Line 25: Planet.Bulk.clathMaxThick_m = 0.0  # was 5e3
    - Line 23: Planet.Steps.nClath = 0  # was 30
```

---

## Next Steps (Priority Order)

### Immediate (Unblock Validation)

1. **Fix PPTest41NoClath.py and PPTest3NoClath.py** (~5 minutes)
   - Replace current files with faithful copies of parent tests
   - Apply only the 4-5 line diffs shown above
   - Commit: "Fix no-clathrate test variants: faithful copies with minimal changes"

2. **Run Local Validation** (~15-20 minutes)
   ```bash
   # Activate venvPP environment
   conda activate venvPP
   
   # Generate structure variants
   python -m PlanetProfile.Inference.prepare_structure_variants \
       --test-module PlanetProfile.Test.PPTest41 \
       --output-dir titan_cache/ \
       --rheology andrade
   
   # Run smoke test (n_effective=10)
   python -m PlanetProfile.Inference.run_inference_cli \
       --config titan_cache/test_mcmc_config.json \
       --output titan_cache/test_mcmc_result.pkl
   
   # Verify result
   python -c "
   from PlanetProfile.Inference import InferenceResult
   result = InferenceResult.load('titan_cache/test_mcmc_result.pkl')
   print('Convergence:', result.convergence_metrics)
   print('Summary:', result.get_summary_stats())
   "
   ```

3. **Compare to Test41 Baseline** (optional, ~30 minutes)
   ```bash
   python PlanetProfile/Test/Test41_mcmc_andrade_no_ocean.py
   # Verify posteriors match statistically
   ```

### Phase 2-3: Parameter UI and Visualization (~2 weeks, Sonnet)

**Prerequisites:** Phase 5 validated

**Phase 2: Parameter Configuration UI** (~800 lines)
- Inference mode selector (MCMC/SBI radio)
- Parameter space configuration (multi-select + prior inputs)
- Observable constraints (Re(k2), Im(k2) with uncertainties)
- Sampler settings (n_effective, random_state, etc.)
- Structure loading (PPTest dropdown)
- **Clathrate checkbox** → selects structure cache path

**Phase 3: Visualization Components** (~500 lines)
- `inference_plots.py`: Corner, trace, k2 scatter, heating distribution
- Plotly-based interactive plots
- Test with dummy data before GUI integration

### Phase 7: GUI Integration (~1 week, Sonnet + Opus Review)

**Prerequisites:** Phases 2-3 complete

**Tasks:**
- Run control buttons (Run/Abort/Force Re-run)
- Subprocess execution (via run_inference_cli.py)
- Real-time progress display (poll JSON file)
- Results display (tabbed layout: Corner|Trace|Diagnostics|k2|Heating)
- Error handling with traceback display

### Phase 8: Export & Replotting (~1 week, Sonnet)

**Prerequisites:** Phase 7 complete

**Tasks:**
- Export handlers (PKL, HDF5, CSV, PNG)
- Standalone replotting CLI (`replot_inference.py`)
- Interactive HTML export (Plotly standalone)
- Download buttons in GUI

### Phase 6: SBI Runner (Optional, Deferred)

**Status:** Not required for v1  
**Rationale:** MCMC sufficient for initial prototype, SBI adds ~350 lines + torch/sbi dependencies

---

## File Inventory

### Completed Files ✅

**Phase 1 (GUI Skeleton):**
- `PlanetProfileApp/pages/Inference.py` (215 lines)
- `PlanetProfileApp/Utilities/inference_cache.py` (203 lines)
- `PlanetProfileApp/Utilities/inference_validation.py` (332 lines)

**Phase 4 (Inference Core):**
- `PlanetProfile/Inference/inference_core.py` (369 lines)
- `PlanetProfile/Inference/forward_models.py` (443 lines)
- `PlanetProfile/Inference/structure_cache.py` (498 lines)
- `PlanetProfile/Inference/run_inference_cli.py` (245 lines)
- `PlanetProfile/Inference/validate_framework.py` (320 lines)
- `PlanetProfile/Inference/__init__.py` (updated)

**Phase 5 (MCMC Runner):**
- `PlanetProfile/Inference/mcmc_runner.py` (350 lines)
- `PlanetProfile/Inference/prepare_structure_variants.py` (200 lines)

### Incomplete Files ❌

- `PlanetProfile/Test/PPTest41NoClath.py` (85 lines) — **NEEDS FIX**
- `PlanetProfile/Test/PPTest3NoClath.py` (70 lines) — **NEEDS FIX**

### Pending Files ⏸️

**Phase 2:**
- `PlanetProfileApp/pages/Inference.py` (extend by ~800 lines)

**Phase 3:**
- `PlanetProfileApp/Utilities/inference_plots.py` (500 lines)

**Phase 7:**
- `PlanetProfileApp/pages/Inference.py` (extend by ~500 lines)

**Phase 8:**
- `PlanetProfile/Inference/export_results.py` (200 lines)
- `PlanetProfile/Inference/replot_inference.py` (300 lines)

---

## Documentation Status

### Complete ✅

- `INFERENCE_FRAMEWORK_STATUS.md` — Overall framework status and architecture
- `VALIDATION_RESULTS.md` — Validation strategy and expected outcomes
- `IMPORT_FIX_SUMMARY.md` — Import error fixes (Phase 4)
- `PHASE5_MCMC_RUNNER_COMPLETE.md` — Phase 5 implementation details
- `MCMC_INFERENCE_GUIDE.md` — Test41-44 baseline documentation

### In Progress ⚠️

- `TASKS.md` (this file) — Current session state

---

## Known Issues and Limitations

### Phase 5 Issues

1. **PPTest*NoClath.py files incomplete** (blocking validation)
   - Status: Identified, fix prepared, not yet applied
   - Impact: Cannot generate structure caches

2. **Checkpoint resumption not tested**
   - Implementation exists but needs validation with long runs

3. **No multi-chain R-hat**
   - pocoMC doesn't expose chains, R-hat approximated (set to 1.0)

### General Limitations

1. **Sandbox environment lacks dependencies**
   - pandas, seafreeze, cmasher, TidalPy not available
   - User must run validation locally with full venvPP

2. **No GPU acceleration**
   - pocoMC is CPU-only (unlike emcee/dynesty with MPI)

3. **Fixed rheology per run**
   - Cannot switch Andrade↔Maxwell mid-inference

---

## Cost Discipline Tracking

**Models Used:**
- Phase 1: Sonnet/Haiku (GUI skeleton, utilities) ✓
- Phase 4: Opus (core infrastructure, architectural decisions) ✓
- Phase 5: Opus (MCMC runner, design decisions) ✓
- Future: Sonnet for Phases 2-3, 7-8 unless complexity requires Opus

**Subagent Spawning:**
- Minimized: Sequential tool calls preferred over parallel agents
- No unnecessary re-reading of files
- Context summarized before passing to subagents

---

## Success Criteria (Phase 5)

### Completed ✓

- [x] MCMCRunner class implemented
- [x] pocoMC integration
- [x] Clathrate toggle via dual structure caches
- [x] Helper CLI for structure variant generation
- [x] Progress callbacks
- [x] Checkpoint save/load
- [x] Convergence diagnostics
- [x] Heating recomputation
- [x] API integration
- [x] Documentation

### Blocked ❌

- [ ] PPTest*NoClath.py files corrected
- [ ] Local validation (structure generation)
- [ ] Local validation (MCMC smoke test n=10)
- [ ] Compare to Test41 baseline

### Next Milestone

**Phase 5 validation complete** → Proceed to Phase 2 (Parameter UI)

---

## Session State

**Current Working Directory:** `~/Library/CloudStorage/Dropbox/planetprofile-genai`

**Active Branch:** genai

**Uncommitted Changes:**
- `PlanetProfile/Inference/mcmc_runner.py` (new)
- `PlanetProfile/Inference/prepare_structure_variants.py` (new)
- `PlanetProfile/Test/PPTest41NoClath.py` (new, needs fix)
- `PlanetProfile/Test/PPTest3NoClath.py` (new, needs fix)
- `PlanetProfile/Inference/__init__.py` (modified)
- `PHASE5_MCMC_RUNNER_COMPLETE.md` (new)
- `TASKS.md` (this file, new)

**Test Config Ready:**
- `titan_cache/test_mcmc_config.json` (2D parameter space, n_eff=10)

**Next Command:**
```bash
# After fixing PPTest*NoClath.py files:
python -m PlanetProfile.Inference.prepare_structure_variants \
    --test-module PlanetProfile.Test.PPTest41 \
    --output-dir titan_cache/
```

---

## Contact Points for User

**Before proceeding to Phase 2, user must:**

1. ✅ Review and approve PPTest*NoClath.py diffs (shown in conversation)
2. ⏸️ Apply fixes to PPTest41NoClath.py and PPTest3NoClath.py
3. ⏸️ Run local validation (structure generation + smoke test)
4. ⏸️ Confirm validation passes (all 5 tests)
5. ⏸️ Provide explicit instruction: "Proceed with Phase 2 (Parameter UI)"

**Questions for user:**

- Confirm PPTest*NoClath.py fix approach (faithful copies with minimal diffs)
- After validation: Which phases to prioritize? (Recommended: 2→3→7→8)
- Push to genai branch after which milestone? (Recommended: after Phase 7)

---

**End of Status Report**
