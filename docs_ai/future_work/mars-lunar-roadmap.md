# Plan: Lunar 1D Modeling and MCMC Refinement (Meta-Optimization with Era)

## Objective
Develop 1D interior structure models for the Lunar nearside and farside, update Perple_X integration, and leverage the `era` (Evaluation and Research Assistant) tool to meta-optimize the MCMC inference process. This plan is designed to be fully portable between Gemini and Claude.

## Context for Successors (Claude/Gemini)
This project, **PlanetProfile**, uses self-consistent thermodynamics to model planetary interiors. We are currently shifting focus from 3D Mars modeling to **Lunar 1D models** (Nearside vs. Farside). The goal is to use MCMC (specifically **pocoMC**) to constrain Lunar parameters while maintaining a balance between physical self-consistency and data-driven flexibility.

## Phase 1: Research & Perple_X Modernization
- [ ] **Data Extraction:** Extract crustal thickness, MoI, mass, and seismic constraints from user-provided Lunar papers.
- [ ] **Perple_X Refactor:** 
    - [ ] Update `PlanetProfile/Thermodynamics/InnerEOS.py` to handle the latest Perple_X output formats.
    - [ ] Generate new `.tab` look-up tables for Bulk Silicate Moon (BSM) compositions.
- [ ] **Baseline Models:** Create `Luna/PPLunaNearside.py` and `Luna/PPLunaFarside.py` using standard literature values as a starting point.

## Phase 2: MCMC Framework & Era Meta-Optimization
Instead of manually tuning MCMC priors and likelihoods, we will use the `era` tool's **Flat UCB Tree Search (FUTS)** algorithm to find the optimal modeling configuration.

### Era Integration Logic:
1. **Problem Definition:** "Find the Lunar model configuration (Nearside/Farside) that best fits seismic and gravity data while maintaining physical plausibility."
2. **Solution Space:** A "Solution" is a tuple of `(PPBody.py config, MCMC JSON config)`.
3. **`execute_fn` (Evaluation):** 
    - Run a short MCMC run (low `n_effective`) using the proposed solution.
    - Calculate a **Score** based on:
        - Goodness-of-fit ($\chi^2$ for MoI, mass, seismic).
        - Convergence diagnostics (Effective Sample Size).
        - Physicality penalties (e.g., density inversions where not expected).
4. **`generate_fn` (Iteration):** 
    - Use the LLM to analyze the results of the previous `execute_fn`.
    - Propose *modifications* to the priors or model structure (e.g., "The core size is hitting the upper bound; shift the log-uniform prior range higher").

## Phase 3: GUI & Educational Generalization
- [ ] **MCMC Tab Refinement:**
    - [ ] Integrate the "Structure Wedge" plots from `mcmc_plots.py` into the Streamlit dashboard.
    - [ ] Add an "Inference Mode" toggle: **Self-Consistent** (locked physics) vs. **Empirical** (flexible rheology).
- [ ] **Modular Registries:** Refactor `parameter_registry.py` and `forward_models.py` into a more extensible library format to support rapid integration of new bodies (like the Moon).

## Phase 4: Execution & Verification
- [ ] **Regression:** Ensure no regressions in existing body tests (`BuildTest.py`).
- [ ] **Lunar Benchmark:** Execute `PPTest_Luna_MCMC.py` to verify the `era`-optimized configuration.
- [ ] **Audit Trail:** Maintain a detailed log of numerical changes and scientific impact in the final session report.

## Verification Checklist
- [ ] `InnerEOS.py` supports new Perple_X version.
- [ ] Lunar Nearside/Farside models achieve stable MCMC convergence.
- [ ] GUI correctly visualizes sampled Lunar structures.
- [ ] Era-optimized priors show improved ESS vs. manual baselines.
