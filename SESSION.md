# Project State & Session Handoff

**Authoritative reference** for the `genai` branch — supersedes the prior `MCMC_ARCHITECTURE.md` (now deleted). When this file conflicts with any other tracked doc, this file wins.

- **Branch:** `genai`
- **Env:** `mamba activate PPcl` → `/Users/svance/mamba/envs/PPcl/bin/python`
- **Repo:** `/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai`
- **Last updated:** 2026-07-07

---

## Part I — Architecture (canonical)

### 1. Sampler: Preconditioned Monte Carlo (pocoMC), not vanilla MCMC

The framework uses **pocoMC** (Karamanis et al. 2022) — Sequential Monte Carlo with normalizing-flow preconditioning. Inside the codebase "MCMC" is shorthand for "PMC sampler driven by pocoMC".

| | MCMC (emcee, MH) | PMC (pocoMC) |
|---|---|---|
| Strategy | Random walk | Tempered importance sampling + walker re-sampling |
| Posterior shape | Best: unimodal, near-Gaussian | Handles multimodal, non-Gaussian, curved |
| Convergence | Manual autocorrelation | Sequential β annealing → automatic stop |
| Log-evidence | Not provided | Returned as `sampler.logz` |
| Test50 5D wall-time at n_eff=500 | minutes–hours | ~1 minute |

### 2. Three-layer forward model

```
Layer 1 — STRUCTURE CONSTRUCTION (PlanetProfile)
  Build Planet from PPBody template + sampled overrides.
  Cache as N-D Tb grid; interpolate per sample.
  Outputs: r, ρ, K, μ, η_base, T, P, phases, layer indices

Layer 2 — FORWARD PHYSICS (overlays; NOT self-consistent)
  Apply sampled rheology overrides on cached structure.
  Plug-in convection / heating / observables (TidalPy, MoonMag).
  Outputs: predicted k₂, Ae, J₂, ...

Layer 3 — OBSERVATION LIKELIHOOD
  Gaussian χ² over observed quantities. log L = Σᵢ log Lᵢ.
  Each constraint is a pluggable module.
```

**The "no required self-consistency" principle.** Rheology overrides are *not* required to match the structure that produced them. PP's self-consistent η is a first guess; the MCMC overrides at the per-phase level. This is the only way to test hypotheses (e.g. Petricca+2025's low-η preference) outside what self-consistent PP models allow. Toggles like `force_self_consistent_eta_Ih`, `lock_HP_to_Ih`, `freeze_structure`, `enforce_no_ocean` lock specific axes when desired.

### 3. Plug-points and registries

| Plug-point | Examples | Interface |
|---|---|---|
| Convection | Yao 2014, Deschamps & Sotin 2001, Kalousová 2018 | `Planet → Planet` |
| Rheology | Andrade, Maxwell, Burgers (future) | `(structure, θ) → complex μ profile per layer` |
| Tidal heating | volumetric, depth-weighted, per-phase | `(μ profile, ω) → q(r)` |
| Parameter registry | `Inference/parameter_registry.py` | `ParameterDef`: id, prior, bounds, units, rebuild flag |
| Constraint registry | (proposed; currently inline per body config) | `ConstraintDef`: id, units, required outputs, obs+σ |

Adding a constraint is *additive* — observations don't know about each other, rheology doesn't know which observations are enabled. Communication is only through structured outputs from the forward layer.

### 4. GUI surfaces: proper vs amortized inference

The Streamlit app exposes **two distinct inference workflows**:

- **MCMC tab (proper inference, slow, exact):** real pocoMC run on local machine. Body selector, parameter table from `parameter_registry`, constraint panel, run button with live progress, posterior viewer (corner, k₂ scatter, heating). Time per run: ~1 min to hours.
- **SBI sub-tab (amortized, fast, approximate):** simulation-based inference with a pre-trained normalizing flow (`sbi` package). Offline: 10⁴–10⁶ (θ, obs) pairs train a conditional flow. Online: user enters obs values, gets posterior in milliseconds. Ties to a specific body + parameter set + structure assumptions baked into the training set.

Per [[project-gui-architecture]] (2026-06-25 confirmed), the layout is body-organized:
```
[Body: Europa | Ganymede | Callisto | Titan]
  └── [Run Config] [Browse Runs] [Posterior] [Model Explorer] [SBI]
```
SBI is a sub-tab inside MCMC, not a top-level page. Both workflows share the parameter registry, constraint registry, and posterior visualization toolkit.

### 5. Roadmap with current status

| Phase | Description | Status (2026-06-29) |
|---|---|---|
| A | Minimal MCMC tab — reproduce Test50 5D from GUI | **Delivered** (`MCMC.py`, configs, runner) |
| B | Refactor Test50 onto toolkit so CLI == GUI | **Not started** (low priority — Test50 still runs as bespoke script) |
| C | Multi-body + multi-constraint: Europa, Ganymede, Callisto, with k₂ + induction | **Mostly delivered** — bodies wired, T25 induction observables committed; production posteriors exist for Europa + Callisto |
| D | Amortized inference tab + 3D scaffolding (Mars), RunRegistry, StructureBackend/ForwardChannel | **In progress** — `mcmc_run_loader` merged (T07a, commit `46c0723d`); SBI runner missing (see Part II §3) |

Future hooks already designed: `OneDStructureBackend` / `ThreeDStructureBackend`, `ForwardChannel` capability declarations (refuse Mars seismic with 1D backend), `RunRegistry` (flat dir of pkls + JSON index, Streamlit-Cloud friendly).

### 6. Key files at a glance

| File | Role |
|---|---|
| `PlanetProfile/Inference/mcmc_runner.py` | `MCMCRunner` class; `generate_sbi_dataset()` stub; JSONL streaming |
| `PlanetProfile/Inference/inference_core.py` | `InferenceConfig`, `InferenceResult`; **dead import `from .sbi_runner import SBIRunner`** at lines 363, 370 |
| `PlanetProfile/Inference/forward_models.py` | Structure cache lookup; rheology overlay; TidalPy k₂/h₂; MoonMag Ae |
| `PlanetProfile/Inference/cache_builder.py` | v2.1 Tb-grid cache; transition refinement; B-layered blend |
| `PlanetProfile/Inference/parameter_registry.py` | `PARAMETER_REGISTRY`, `PARAMETER_PRESETS`, prior definitions |
| `PlanetProfile/Inference/configs/*.json` | 10 body configs (Titan 2, Europa 2, Ganymede 1, Callisto 4) |
| `PlanetProfile/Inference/diagnostics/` | `run_audit.py`, `evidence_comparison.py`, `progress_tail.py` |
| `PlanetProfile/Test/mcmc_results/{body}/{run_name}/*.pkl` | All MCMC result pkls live here |
| `PlanetProfileApp/pages/Inference.py` | MCMC inference tab (1699 lines — ripe for refactor into body-organized layout) |
| `PlanetProfileApp/pages/AmortizedInference.py` | SBI scaffold (75 lines, button `disabled=True`) |
| `PlanetProfileApp/pages/Exploreogram.py` | Parametric sweep + Plotly/matplotlib viz (2215 lines) |
| `PlanetProfileApp/pages/CompareRuns.py` | Multi-run comparison (not yet linked from Inference tab) |
| `PlanetProfileApp/utils/mcmc_run_loader.py` | Browse-runs glob + run-record builder |

### 7. Critical bug history (do not re-introduce)

- **pocoMC `posterior()` returns `(samples, weights, logl, logp)`** — not `(samples, logl, logp, weights)`. Fixed in Phase C1 bug #3. The `samples` array was always correct; pre-fix log_likelihoods were importance weights.
- **v2.1 cache list format** — new caches use `'Tb_K_grid'` + `'structures'` + `'transitions'`. Old v2.0 used `'grid_cache'` + `'grid_Tb_values'`. `mcmc_runner.py` CMR² dispatcher handles both. Any new cache-reading code must too.
- **NumPy 2.x migration done** (commit `a9f24584`) — `np.complex_` → `np.complex128` etc. across MoonMag and core modules. Do not reintroduce deprecated aliases.
- **matplotlib 3.9+ removed `matplotlib.cm.get_cmap`** — use `matplotlib.colormaps[name]` (or `mpl.colormaps[name]`). Fixed 2026-07-07 in `PlanetProfile/Plotting/ExplorationPlots.py` and `PlanetProfileApp/Utilities/exploreogram_plotly.py`. Env has matplotlib 3.11.0; do not reintroduce `get_cmap`.
- **numpy ABI split in PPcl env (2026-07-07)** — conda had clobbered numpy 2.4.6 files with numpy 1.26.4 while scipy/matplotlib/etc. were built against numpy 2.x, producing `numpy._core.multiarray failed to import` on any scientific import (surfaced as "Error loading planet module" in the GUI). Fix: `pip install --force-reinstall --no-deps numpy==2.4.6`. Leftover `numpy-1.26.4.dist-info/` can be deleted. If it recurs, look for a `mamba/conda install` that pins numpy <2.

---

## Part II — Current state of work (2026-06-29)

### 1. Just shipped (verified)

All four Exploreogram induction bugs are verified in the running app:

| Bug | Status | Commit |
|---|---|---|
| 1 — matplotlib component routing (Bx/By/Bz → Bi1{x,y,z}_nT zName) | `verified` (PDF `Europa_explore_..._MultiSubplot_Bi1x_nT.pdf`) | `b7da7807` |
| 2 — salinity↔conductivity axis sort | `verified` | `8a631cc4` |
| 3 — Plotly secondary axis (one derived) | `verified` | `8a631cc4` |
| 4 — Plotly component selector | `verified` | `8a631cc4` |

Verifier scenarios live at `.claude/skills/verifier-streamlit/scenarios/bug{1,2,3,4}_*.py`. CLAUDE.md verification-discipline section was restored in `b7da7807`.

Bug 1 patch detail (`PlanetProfileApp/pages/Exploreogram.py`):
- `~L1438`: added `_COMP_TO_ZNAME = {'Bx':'Bi1x_nT','By':'Bi1y_nT','Bz':'Bi1z_nT','Total':'Bi1Tot_nT'}`
- `~L1439`: `_pdf_append = results['zName']` (dropped misleading `f"{zName}_{_ic_mpl}"` suffix)
- `~L1475`: `real_imaginary` branch when `zName in ('Amp_nT','phase_deg')`: route via `_COMP_TO_ZNAME.get(_ic_mpl,'Bi1Tot_nT')`
- `~L1573`: added bare `Bi1x_nT/Bi1y_nT/Bi1z_nT` to `_will_plot_ri`
- `~L1417`: banner narrowed to `amplitude_phase` mode only

Backend supports this with no changes (`EssentialHelpers.py:848`, `defineStructs.py:2315-2318`).

### 2. T25 induction MCMC state

Production runs exist for Europa + Callisto at `PlanetProfile/Test/mcmc_results/{Europa,Callisto}/T25_production_*/`:

| Run | n_samples | log_likelihoods | log_evidence |
|---|---|---|---|
| Europa seawater 7D | 4374 × 7 | [-6.4, -0.4] (healthy) | **None** — predates log_Z metadata fix |
| Callisto NaCl 8D | 4451 × 8 | [-24.2, -17.1] (healthy) | **None** — same |
| Europa T25 smoke | 4096 × 7 | -1e30 everywhere | degenerate; ignore |

Δln Z Bayes-factor table (nominal vs −5% vs −10% CMR² Callisto) is blocked on filling in `log_evidence` — see Part III §1.D.

### 3. Open workstreams (priority order for Clipper deadline ~2026-07-13)

**A. Exploreogram student handoff** — Exploreogram is verified and ready for a student to take over for usability polish + new science. Need: `docs/EXPLOREOGRAM_STUDENT_GUIDE.md` capturing file map, session-state keys, verifier scenarios as regression suite, four-bug commit history, suggested "good first issues".

**B. Body-organized MCMC tab refactor** — `Inference.py` is currently flat at 1699 lines. Target structure per Part I §4. Split into `PlanetProfileApp/utils/mcmc_body_view.py` + `PlanetProfileApp/inference_tabs/{run_config,browse_runs,posterior,model_explorer,sbi}.py` + thin `Inference.py` wiring the body selector.

**C. Implement `sbi_runner.py`** — blocks all SBI work. Dead import at `Inference/inference_core.py:363,370`. Minimum viable: `SBIRunner.train(dataset_path, output_path)`, `SBIRunner.posterior(obs, artifact_path)`, `SBIRunner.run(progress_callback)` matching MCMCRunner interface. Use SNPE-C via `sbi` package. Then run `generate_sbi_dataset()` (stub in `mcmc_runner.py` ~L722) for Titan no-ocean 5D, target 50K–100K (θ, x) pairs. Train + serialize to `sbi_artifacts/titan_andrade_noocean_posterior.pt`.

**D. T25 log_Z backfill** — re-run Europa and Callisto T25 production configs at original n_eff with current MCMCRunner (post-metadata-fix) to populate `log_evidence`. Closes the Δln Z Bayes-factor table.

### 4. Important but not critical-path

- **Live JSONL progress display in Inference tab** — runner already writes JSONL every 2s; GUI just needs an `st.empty()` polling loop wired to `diagnostics/progress_tail.py`.
- **Multi-run comparison handoff** — `CompareRuns.py` exists; needs a button in Inference tab.
- **Body-specific k₂ priors for Ganymede / Callisto** (currently Europa placeholders from Mazarico+2023). Requires literature review (Anderson, Moore, Kamata, Hussmann) and approval gate.

### 5. Deferred

- Ice VI full-TBL Kalousová-style convection — out of scope per [[project-goals-confirmed]]; stub stays
- Gravity-coefficient forward model (J_n, C_nm, S_nm) — large scope, separate phase
- Ganymede induction (pure-H₂O config skips MagneticInduction by design)
- Phase B Test50 CLI refactor — low priority once SBI is moving
- Mars 3D StructureBackend / ForwardChannel — designed, not built

---

## Part III — Pickup guidance for incoming agents

### 1. Economical agent fan-out

| Workstream | Lead | Pattern |
|---|---|---|
| A. Exploreogram student handoff | Sonnet | One agent writes the guide; self-contained, no review needed |
| B. GUI body-organized refactor | Opus plans → Sonnet implements → Opus reviews | 1 plan turn + 2 implementation turns + verifier-streamlit regression check |
| C. SBI runner | Opus drafts API → Sonnet implements → scientific-reviewer + verifier | 1 plan + 2-3 implementation + 1 review |
| D. T25 log_Z backfill | Sonnet launches + summarizes | 1 launch turn + wall-time + 1 summary turn |

A, C, D are independent and can run in parallel on day 1. B and C combine downstream (SBI sub-tab needs SBI runner).

### 2. Devin (Sonnet 4.6 in IDE)

Devin maintains `DEVIN.html` as their working tracker. **Until 2026-06-29 their copy claimed Exploreogram bugs were unverified and cited a stale `mcmc_results/` path at repo root.** They have been asked to update to reflect the current state above. If you see them about to layer fixes on those, redirect.

Devin's strengths here: mechanical refactors (workstream B), documentation generation (workstream A), test/verifier scaffolding. They should NOT change scientific assumptions, touch `PlanetProfile/Test/` files, or push without instruction. Plan nontrivial work first (≥50 lines or multi-file).

### 3. Hard constraints (carry-forward, do not relax)

Per `CLAUDE.md` and confirmed memory:

- `python -m py_compile` is **not** verification. UI changes must be observed in the running app (screenshot, PDF, or `/verify` skill).
- Status vocabulary: `verified` | `implemented, unverified` | `not implemented`. Never `done`/`fixed`/`complete`/`syntax clean`.
- Do not layer fixes on `implemented, unverified` work.
- Do not alter `PlanetProfile/Test/` files without permission.
- Do not change fundamental scientific assumptions.
- Sonnet/Haiku for implementation; Opus for planning and review.
- Small commits; push only on explicit instruction.
- Preserve migration compatibility with Vellard/Niesyt/Chang branches (small, well-scoped commits).

### 4. Verifier infrastructure

```
Env: ~/.venvs/pp-verify/bin/python  (Python 3.14, Playwright + Chromium)
Skill: .claude/skills/verifier-streamlit/
Scenarios: .claude/skills/verifier-streamlit/scenarios/bug{1,2,3,4}_*.py
Output: /tmp/claude-503/verify/<UTC-stamp>/<scenario>/
Run all: ~/.venvs/pp-verify/bin/python .claude/skills/verifier-streamlit/run.py
Streamlit: localhost:8501; user keeps a long-running session.
Launch:    /Users/svance/mamba/envs/PPcl/bin/streamlit run PlanetProfileApp/PlanetProfileApp.py --server.port 8501
```

### 5. Pickup checklist

1. Read this SESSION.md fully.
2. Skim `DEVIN.html` for current Devin state; cross-check against Part II.
3. Verify with `git log --oneline -10` that the Exploreogram bug commits (`8a631cc4`, `b7da7807`) are present.
4. Confirm `PlanetProfile/Inference/sbi_runner.py` still missing — if Devin built it, skip workstream C step 1.
5. Pick workstream(s) per user priority for the day; fan out per §1.

---

## References (carried over from prior architecture doc)

- Karamanis et al. 2022, *Accelerating astronomical and cosmological inference with preconditioned Monte Carlo*, MNRAS 516, 1644.
- Petricca et al. 2025, *Constraints on Titan's interior from Cassini gravity and rotation*, JGR Planets.
- Yao et al. 2014, *Convective parameterization of stagnant lid in Ice Ih layers*.
- Kalousová & Sotin 2018, *Melting in high-pressure ice layers of large ocean worlds*, GRL 45, 8096.
- sbi package: https://www.mackelab.org/sbi/
- pocoMC: https://github.com/minaskar/pocomc

Related tracked docs (user-facing guides, linked from `README.md`):
- `MCMC_INFERENCE_GUIDE.md` — tutorial walkthrough of running Test48 / Test50
- `CONVECTION_GUIDE.md`, `ARRHENIUS_VISCOSITY_GUIDE.md`, `TIDAL_LOVE_NUMBERS_GUIDE.md` — physics module guides
- `CHANGELOG.md` — chronological change record
- `PP_ARCHITECTURE_TESTING.md` — repository architecture + test organization

Claude-only guidance (untracked, in `.claude/`):
- `CLAUDE.md` (at repo root — auto-loaded by harness)
- `.claude/CLAUDE.prior.md` — historical full version of CLAUDE.md (reference only)
