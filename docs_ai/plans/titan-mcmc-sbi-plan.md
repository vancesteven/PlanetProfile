# Plan: Titan MCMC & Tidal Heat Partitioning (SBI-Ready)

## Objective
Refine the MCMC framework using Titan as the primary test case to investigate tidal heat partitioning (Silicate, HP Ices, Ice Ih). Prepare the pipeline for **Simulation-Based Inference (SBI)** pre-training and enhance the GUI for educational, body-specific prior configuration.

## Context for Successors
We are focusing on **Titan** (with and without an ocean). The scientific goal is to track where tidal energy is dissipated. This MCMC work is the foundation for an "Amortized Inference" (SBI) tab in the GUI. The **era** tool will be used to optimize the MCMC parameter space (priors/likelihoods) for these specific Titan scenarios. **Mars** modeling is the next major body focus, but is on hold until this MCMC infrastructure is solid.

## Phase 1: Titan MCMC & Heat Partitioning
- [x] **Tidal Heat Tracking:**
    - [x] Ensure `PlanetProfile/Inference/forward_models.py` correctly reports dissipated power in three distinct reservoirs: **Silicate Core**, **HP Ice Layer**, and **Ice Ih Shell**.
    - [x] Update the MCMC likelihood to (optionally) include constraints on total tidal heating if literature values are available.
- [ ] **Ocean vs. No-Ocean Comparison:**
    - [ ] Refine the configurations for `Test50` (No-Ocean Titan) and an Ocean-bearing counterpart.
    - [ ] Use `era` to meta-optimize the prior bounds for both cases to ensure stable convergence and physically plausible heat partitioning.

## Phase 2: SBI Pre-training Pipeline
- [x] **Sampling for SBI:**
    - [x] Modify the MCMC runner to export (theta, observable) pairs in a format suitable for training normalizing flows (e.g., via the `sbi` Python package).
    - [x] Ensure the prior distributions used in MCMC are consistent with the ranges needed for robust SBI training.

## Phase 3: GUI & Educational Enhancements
- [x] **MCMC Tab Updates:**
    - [x] **Manual Bounds:** Allow users to override the `era`-optimized or default priors via the UI. (Already supported, verified)
    - [x] **Body-Specific Hints:** Implement a "Prior Recommendations" system that provides clues based on the selected body (e.g., "For Titan, Ice Ih viscosity is typically expected between $10^{13}$ and $10^{15}$ Pa·s").
    - [x] **Heat Visualization:** Add a pie chart or bar graph to the MCMC results showing the percentage of tidal heat partitioned into the three reservoirs.
- [x] **Inference (SBI) Tab Scaffolding:** Create the initial page for the SBI tab, which will eventually load pre-trained models for instant posterior estimation.

## Phase 4: Era Meta-Optimization (Titan Focus)

We will use `era`'s **Flat UCB Tree Search (FUTS)** to automate the discovery of optimal MCMC configurations. This avoids the "manual trial-and-error" phase of setting priors and likelihood weights.

- [x] **Infrastructure Implementation:** Created `planetprofile_era_titan.py` which implements the `execute_fn` and `generate_fn` wrappers for PlanetProfile.
- [x] **Scoring Logic:** Defined a multi-objective score based on ESS, data fit, and physical heat partitioning.
- [ ] **LLM Wiring:** Wire the `generate_fn` to a live LLM (Gemini/Claude) for automated prior refinement.

## Verification & Testing
- [ ] **Titan Regression:** Run `Test50` and verify the heat partitioning outputs match previous benchmarks.
- [ ] **Prior Sensitivity Test:** Verify that manual prior overrides in the GUI are correctly propagated to the `pocoMC` sampler.
- [ ] **SBI Export Test:** Verify that a short run correctly generates a training-ready dataset of parameters and observables.

## Post-MCMC Roadmap
- [ ] **Mars 1D Models:** Revisit Mars modeling (composition-based literature) once the MCMC/SBI pipeline is verified on Titan.
