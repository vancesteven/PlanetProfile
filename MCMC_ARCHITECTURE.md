# PlanetProfile MCMC Architecture

This document records the design rationale for PlanetProfile's Bayesian-inference framework: what the forward model is, why we use **preconditioned** Monte Carlo rather than vanilla MCMC, what the toolkit's plug-points are, and how the GUI surfaces (proper inference + amortized inference) fit together.

It is intended as a contributor / reviewer reference, not a tutorial. For a tutorial-style walkthrough see `MCMC_INFERENCE_GUIDE.md`.

---

## 1. Sampler: Preconditioned Monte Carlo (pocoMC), not MCMC

The framework uses **pocoMC** (Karamanis et al. 2022, https://github.com/minaskar/pocomc), which implements **Preconditioned Monte Carlo (PMC)** — a Sequential Monte Carlo (SMC) sampler with normalizing-flow preconditioning.

This is **not** vanilla MCMC. The distinction matters:

| | MCMC (e.g. emcee, Metropolis-Hastings) | PMC (pocoMC) |
|---|---|---|
| Sampling strategy | Random walk in posterior | Tempered importance sampling, walker re-sampling |
| Posterior shape | Best for unimodal, near-Gaussian | Handles multimodal, non-Gaussian, curved |
| Burn-in / convergence | Manual; needs autocorrelation diagnostics | Sequential β annealing → automatic stop |
| Evidence (log Z) | Not provided; needs nested sampling for it | Returned directly |
| Wall-clock efficiency | Scales poorly with curvature / dimension | Strong in moderate D (~5–30) thanks to preconditioner |
| What's "preconditioned" | n/a | A normalizing flow learns the posterior shape and reshapes proposals to it |

For our problem class (5–10 D, weakly informative priors, observable depends non-linearly on viscosity), PMC is the right tool. Test50's full 8D run finished in **~1 minute** at N_eff=500 with the preconditioner doing the heavy lifting. A naive emcee run with comparable acceptance would be 10–100× slower.

When you read "MCMC" in the codebase or Test files, it's shorthand for "PMC sampler driven by pocoMC". The distinction is:
- The **outer protocol** is sequential importance sampling with annealing (β: 0 → 1).
- The **inner refresh** is short MCMC chains over the preconditioned proposal — those are MCMC-style moves, but they live inside the PMC scaffold.

**Reference:** Karamanis, M., Beutler, F., Peacock, J. A., Nabergoj, D., & Seljak, U. (2022). *Accelerating astronomical and cosmological inference with preconditioned Monte Carlo.* MNRAS 516, 1644. arXiv:2207.05652.

---

## 2. Forward modeling: three orthogonal layers

```
┌─────────────────────────────────────────────────────────────┐
│  Layer 1: STRUCTURE CONSTRUCTION (PlanetProfile)            │
│  • Build a Planet from PPBody.py + sampled overrides        │
│  • Some parameters force structure rebuild                  │
│    (Tb, ocean salinity, clathrate thickness, core size)     │
│  • Cache as N-D grid; interpolate per sample                │
│  Outputs: r, ρ, K, μ, η_base, T, P, phases, layer indices   │
├─────────────────────────────────────────────────────────────┤
│  Layer 2: FORWARD PHYSICS (overlays — NOT self-consistent)  │
│  • Apply sampled rheology overrides on cached structure:    │
│    per-phase η, α, ζ, μ, etc.                                │
│  • Plug in alternative convection / heating models          │
│  • Compute observables: k₂, induction, gravity, ...          │
│  Outputs: predicted observables                              │
├─────────────────────────────────────────────────────────────┤
│  Layer 3: OBSERVATION LIKELIHOOD                            │
│  • Gaussian χ² over observed quantities                      │
│  • Each constraint type is a pluggable module                │
│  Output: log L                                               │
└─────────────────────────────────────────────────────────────┘
```

### Why these are independent

- **Structure ≠ rheology.** PP's self-consistent η is a *first guess*. The MCMC overrides it at the per-phase level. The forward layer treats the cached structure as fixed geometry/density and sweeps rheology freely.
- **Rheology ≠ observation.** The same η + α + ζ samples generate k₂ *and* dissipation *and* tidal heating *and* (eventually) gravity / seismic predictions. Adding a new observable doesn't change how rheology is evaluated.
- **Constraints stack additively.** log L = Σᵢ log Lᵢ over enabled constraint modules. Adding Mars J₂/J₃/seismic = adding modules, not rewriting the inference.

### The "no required self-consistency" principle

A key design decision: **rheology overrides are not required to be consistent with the structure that produced them.** If PP would compute a particular η_Ih(P,T) from Arrhenius + grain-size assumptions, the MCMC can ignore that and sample η_Ih freely from the prior. This is the only way to test hypotheses like Petricca et al. (2025), where the data preference for low η lies outside what self-consistent ice rheology models naturally allow.

The framework supports a switch (`force_self_consistent_eta_Ih` etc.) to lock specific parameters to PP-derived values when you *do* want self-consistency. By default, sampling overrides PP.

---

## 3. Plug-points and registries

Three plug-points, each with a registered name and a uniform interface:

| Plug-point | Examples in tree | Interface |
|---|---|---|
| Convection model | Yao 2014 spherical, Deschamps & Sotin 2001, Kalousova 2018 two-phase | `Planet → Planet` (pre-structure-extraction) |
| Rheology | Andrade, Maxwell, Burgers (future) | `(structure, theta) → complex μ profile per layer` |
| Tidal-heating distribution | volumetric uniform, depth-weighted, per-phase | `(complex μ profile, ω) → q(r)` |

Adding a new model = registering it under a name; the GUI then offers it as a dropdown without code changes there.

### Parameter registry

Each MCMC parameter is described by a `ParameterDef` (`PlanetProfile/Inference/parameter_registry.py`):

```
id, label, latex_label, description
category            ∈ {rheology, structure, ocean, magnetic, ...}
prior_type          ∈ {uniform, log-uniform, normal}
default_bounds, default_mean, default_std
units
requires_structure_rebuild   # → which axis of the structure grid?
rheology_constraint          # only valid if rheology=='andrade' etc.
group                        # NEW (proposed): bundle params (e.g. all η_HP share a slider)
```

### Constraint registry (proposed addition)

Each observable is described by a `ConstraintDef`:

```
id                  e.g. 'k2_real', 'cmr2', 'gravity_J2', 'induction_amp', 'seismic_Vs'
label, units
category            ∈ {tidal, gravity, magnetic, seismic}
required_outputs    forward-model channels needed (e.g. ['k2'])
default_obs         {value, sigma}      # body-specific override
forward_module      Python qualified name; runs once per sample
```

This generalizes what is currently a body-specific OBS dict in each test script.

### Self-consistency switches

Top-level toggles exposed by the framework:

- `force_self_consistent_eta_Ih` — η_Ih is *not* a free parameter; recomputed from PP at each sample (forces structure rebuild).
- `lock_HP_to_Ih` — collapse η_III, η_V, η_VI into one parameter (what Test50 5D does).
- `freeze_structure` — sample only rheology, never rebuild structure (what most current Test4× scripts do).
- `enforce_no_ocean` — reject samples whose interpolated T crosses Tm_Ih(P) (Test50's no-ocean safeguard).

Each is a checkbox in the GUI Self-Consistency panel.

---

## 4. Worked example: Test50 in 8D and 5D

`Test50_mcmc_andrade_noocean_yao2014.py` is the reference no-ocean Titan inference. It is the toolkit-ready example.

**8D:** `[α, log₁₀ζ, log₁₀η_Ih, log₁₀η_III, log₁₀η_V, log₁₀η_VI, log₁₀η_sil, T_b]`

**5D (`scripts/explore_test50_5D.py`):** `[α, log₁₀ζ, log₁₀η_Ih, log₁₀η_HP, log₁₀η_sil]` with T_b fixed at the upper grid edge (PPTest50 default = `TtripleIh_III_L_K − 0.2 K`) and HP ices locked together.

Both use the same:
- structure cache (`titan_allice_yao2014_structure_grid.pkl`, 9-point Tb grid)
- forward_model (with HP ices independently overridable)
- OBS dict from Petricca et al. 2025
- Andrade rheology
- pocoMC sampler

**Posterior comparison:**

| | 8D median [16, 84%] | 5D median [16, 84%] |
|---|---|---|
| α | 0.306 [0.21, 0.40] | 0.298 [0.20, 0.39] |
| log₁₀ ζ | −1.74 [−2.50, +0.37] | −1.67 [−2.58, +0.64] |
| log₁₀ η_Ih | 12.72 [10.96, 14.78] | 12.70 [10.83, 14.77] |
| log₁₀ η_HP (combined) | — | 11.74 [11.15, 14.17] |
| log₁₀ η_III | 12.54 [10.73, 14.68] | (collapsed) |
| log₁₀ η_V | 11.79 [10.55, 14.29] | (collapsed) |
| log₁₀ η_VI | 12.50 [10.79, 14.75] | (collapsed) |
| log₁₀ η_sil | 19.68 [18.56, 21.07] | 19.61 [18.49, 21.05] |
| T_b (K) | 250.0 [249.4, 250.6] | (fixed at 250.965) |
| Re(k₂) | 0.555 [0.49, 0.62] | 0.560 [0.49, 0.62] |
| ‖Im(k₂)‖ | 0.095 [0.07, 0.13] | 0.088 [0.06, 0.12] |

**Reading the diff:**

- η_HP in 5D (10¹¹·⁷⁴) sits *below* the marginalized HP medians in 8D (~10¹²·⁵). When III/V/VI are locked together, the sampler can pick a single dissipation-optimal value rather than averaging over three weakly-coupled axes; that single value drifts down toward Petricca's low-η preference. This is informative: HP-ice dissipation is the dominant Im(k₂) contribution, and the data want it lower than Ih.
- T_b posterior in 8D was uniform-flat across [249.0, 250.59] — the 8th parameter added no information because structural changes across that 2 K range are <0.5% (eLid/TBL/D_hsphere). Removing it saved 1 dimension at no cost.
- α and η_Ih are essentially unchanged. The shell-Ih dissipation is well-constrained by Re(k₂); collapsing HP doesn't disturb it.
- η_sil ~10¹⁹·⁶ Pa·s in both. **This is independent of the silicate thermal-profile bug** because PPTest50 sets `CONSTANT_INNER_DENSITY = True` — silicate ρ/K/μ are constants, not derived from the (incorrect) T(r). And η_sil is a free MCMC parameter, overriding any PP value. The posterior reflects the data's preference given clean elastic moduli.

**Conclusion: 5D is the cleaner physical model when no T_b sensitivity exists** in the data, and locked HP gives a sharper read on HP dissipation than three independent axes that the data can't separate.

---

## 5. Constraint extensibility — what it costs to add a new observable

The architecture is designed so that adding constraints is *additive*, not *intrusive*.

**Adding tidal Im(k₂) at a different frequency** (multi-frequency Cassini fits): one new ConstraintDef pointing at the existing radial solver with a different ω. No rheology code touched.

**Adding induction (Europa, Ganymede)**: forward channel `forward_models.forward_model_induction` already exists. ConstraintDef `induction_amp_phase` consumes it. No structure-cache changes.

**Adding gravity J₂, J₃, ..., C_{nm}, S_{nm} (Mars 3D)**: each is a ConstraintDef whose forward routine integrates ρ(r,θ,φ) against a spherical-harmonic kernel. *Zero changes* to rheology code path. Requires the structure-construction layer to deliver 3D ρ; see §6.

**Adding seismic Vs(r) at station depths**: ConstraintDef with `required_outputs=['Vs_profile']`. Forward physics adds a Vs computation from existing μ and ρ in the cache. Each station/depth measurement is one χ² term.

**Adding seismic Vs(r,θ,φ) for InSight-style point measurements (Mars 3D)**: same as above but consumes the 3D structure backend.

The principle: **observations don't know about each other**, and **rheology doesn't know which observations are enabled**. They communicate only through structured outputs from the forward layer.

---

## 6. Future-proofing for 3D models (Mars and beyond)

Three interfaces are needed for the 3D direction:

**`StructureBackend`**:
- `OneDStructureBackend` — current Planet object; returns r, ρ, μ, η, T, P arrays.
- `ThreeDStructureBackend` (future) — returns the same plus a lateral basis (spherical-harmonic coefficients of ρ, μ, etc.).

**`ForwardChannel`**:
- Each physical predictor (k₂, J_n, induction, Vs(r,θ,φ)) declares its required `StructureBackend` capability.
- Mars seismic + high-order gravity demand 3D backend; refusing to enable them with a 1D backend is a UI-level guard.

**`RunRegistry`**:
- Every MCMC run gets a UUID, a config snapshot, and a result pkl.
- Posterior-comparison page lets users tag/compare/diff runs across sessions.
- Streamlit-cloud-friendly: the registry is a flat directory of pkls + a JSON index.

Two `InferenceConfig` fields make saved result pkls self-describing for post-hoc analysis:
- `arrhenius_params` — Arrhenius viscosity parameters (e.g. `{"E_Ih_J_per_mol": 60000}`). Stored on the config so `reanalyze_k2_from_pickle` auto-loads them without the caller re-supplying them. Set this as a **top-level JSON key** in new body configs (preferred over `sampler_settings["arrhenius_params"]`).
- `planet_template_module` — importable PPTest module (e.g. `"PlanetProfile.Test.PPTest48"`) used by `plot_structure_wedge_pp` to re-run PlanetProfile at the posterior median and call `PlotWedge`. Must be set in the config JSON for wedge plots to work without passing the module explicitly.

These are not built yet. Phase D (below) introduces them.

### Posterior-point wedge plot

`mcmc_plots.plot_structure_wedge_pp(result, grid_cache, output_path, ...)`
re-runs PlanetProfile's full forward model at a chosen posterior point
(`use='median' | 'best_fit' | 'sample'`) and invokes the canonical
`PlotWedge` to produce a manuscript-quality figure of the inferred interior.
This bypasses the cached-structure interpolation that is used during
likelihood evaluation, so the wedge reflects exactly what PP itself would
draw for the inferred parameters.

A subtle but important pitfall: PP defaults `FigMisc.figFormat='pdf'`, so a
file written to `*.png` can otherwise be a PDF byte stream in a `.png`
wrapper. The wrapper temporarily overrides `FigMisc.figFormat` to match the
actual extension on `output_path` (and restores it in a `finally` block).
It also resets `Params.SKIP_PLOTS=False` (cleared by the silent forward
run) and patches `tight_layout`/`savefig` to handle zero-width Wedge
patches without crashing. See `PlanetProfile/Inference/README.md` §"Wedge
plot (PP canonical)" for full details.

Typical call:

```python
plot_structure_wedge_pp(result, grid_cache,
    f"mcmc_results/{body}/{run_label}/wedge_pp.png", use='median')
```

---

## 7. GUI surfaces: proper vs amortized inference

The Streamlit app exposes **two distinct inference workflows**, not one.

### MCMC tab — **proper inference** (slow, exact)

What the doc-plan elsewhere calls Phases A–D. This tab drives a real pocoMC run on the user's local machine (or a backend worker).

- Body selector → BodyConfig
- Parameter table from `parameter_registry`, with sample/fix toggles, prior bounds, group-locking
- Constraint panel from `constraint_registry`, with checkboxes and editable obs values
- Self-consistency panel
- Run button → live progress (sampler iteration count, β value, ESS)
- Posterior viewer (corner plot, k₂ scatter, heating)

Time per run: ~1 min (Test50 5D) to ~hours (high-D, multi-constraint, structure-rebuilding).

This tab is for the **science user** doing a real fit.

### Inference tab — **amortized inference** (fast, approximate)

NEW concept (per user direction): a separate tab using **simulation-based inference (SBI)** with a normalizing-flow density estimator (e.g. via the `sbi` package).

The flow:

1. **Offline training (developer / CI):** sample 10⁴–10⁶ (θ, observable) pairs by running the forward model across the prior. Train a conditional normalizing flow `q(θ | observable)`.
2. **Online inference (Streamlit Cloud):** at runtime, the user enters obs values; the network returns a posterior in **milliseconds**.

Why this matters:
- Streamlit Cloud has no compute budget for full pocoMC runs.
- Casual users (mission scientists exploring "what if k₂=0.55±0.04?") get instant feedback.
- Exploration is interactive — drag an obs slider, watch the posterior reshape live.

Trade-offs:
- Approximate (NN error, training-set coverage limits).
- Tied to a specific body + parameter set + structure assumptions baked into the training set.
- Not a substitute for the MCMC tab when accuracy matters or the parameter space changes.

The Inference tab consumes **pre-trained models** distributed with the app (or downloaded from a model registry). The user flow is:

- Select a pre-trained model (e.g. "Titan no-ocean / Andrade / 5D / Petricca-2025 priors").
- Enter or load observed values.
- See posterior corner plot + k₂ scatter in <1 s.
- Optionally export the posterior; optionally launch a confirmation MCMC run on the same config (handoff to MCMC tab).

Both tabs share the same parameter registry, constraint registry, and posterior visualization toolkit. They differ only in the sampler.

---

## 8. Roadmap (Phases A–D)

### Phase A — minimal MCMC tab (Titan all-ice only)

Goal: reproduce Test50 5D from the GUI. Validate the toolkit ↔ GUI binding.

New: `PlanetProfileApp/pages/MCMC.py`, `PlanetProfileApp/Utilities/inference_config.py`.
Consume (no changes): `parameter_registry.py`, `forward_models.py`, `mcmc_runner.py`.

### Phase B — refactor Test50 onto the toolkit

Goal: prove GUI run path = CLI run path. Both go through `MCMCRunner.run(config)`.

- Move Test50's inline `forward_model` body-specific bits into `BodyConfig` for Titan no-ocean.
- Replace bespoke prior + likelihood builders with the toolkit's runner.
- Test50 becomes a thin script: `python Test50_mcmc.py` ≡ `MCMCRunner.run(load_config('test50.yaml'))`.
- 5D and 8D variants in `scripts/` differ only in config (no code duplication).

### Phase C — multi-body + multi-constraint (MCMC tab)

Goal: Europa, Ganymede, Callisto. Combine k₂ + induction + (where available) gravity.

- BodyConfigs for each body.
- ConstraintDefs for induction, gravity (J₂, C₂₂).
- Body selector becomes meaningful.
- Posterior viewer page: load .pkl, overlay multiple runs, compare across bodies.

### Phase D — Inference tab (amortized) + Mars 3D scaffolding

Goal: the fast inference tab (ships as separate workflow) + extensibility hooks.

- `StructureBackend` interface.
- `ForwardChannel` interface with capability declarations.
- `RunRegistry` for cross-session run management.
- Inference-tab page consuming pre-trained `sbi` posteriors.
- CI pipeline / Makefile target for retraining a pre-trained model when the underlying body config changes.

---

## 9. Wedge-plot sanitization in `plot_structure_wedge_pp`

`mcmc_plots.py::plot_structure_wedge_pp` re-runs PlanetProfile internally via
`_run_pp_with_overrides` (CALC_NEW, NO_SAVEFILE, SKIP_PLOTS=True) to render a
posterior-median wedge diagram. Because `SKIP_PLOTS=True` bypasses `GetLayerMeans`
(called in `Main.PlanetProfile` only when plotting), the scalar layer-summary
attributes that `PlotWedge` needs are not set. The sanitization block (lines
~1584–1750) repairs them in-place before calling `PlotWedge`.

### Attributes set by `GetLayerMeans` but missing after re-run

| Attribute | PP meaning | How sanitized |
|---|---|---|
| `dzIceI_km`, `zIceI_m` | Ice Ih total thickness / depth to top | Reconstructed from `(r_m, phase==1)` |
| `dzClath_km`, `zClath_km` | Clathrate thickness / depth to top | Reconstructed from `(r_m, phase==30)` |
| `dzIceII_km`, `zIceII_m` | Analogues for HP ices II/III/V/VI | Same reconstruction |
| `dzSilPorous_km`, `dzFeS_km`, etc. | Modes PP may not have exercised | NaN → 0 (PlotWedge guards with `> 0`) |
| `Dconv_m`, `deltaTBL_m`, `eLid_m` | Convective zone / conductive lid | Finite when PP ran; NaN → reconstruct from `dzIceI - eLid` |
| `D_km`, `zb_km` | Ocean thickness / depth | NaN → 0 for no-ocean models |
| `Core.Rmean_m`, `Sil.Rmean_m` | Mean radii for silicate/core boundaries | Reconstructed from phase masks |

### Grid convention

PP uses a cell-centred grid: `r_m` has N+1 boundary radii (descending from surface
to centre); `phase` has N cell values. Reconstruction accepts both `n_r == n_phase`
and `n_r == n_phase + 1`, using `r_m[last+1]` as the cell bottom when the N+1
boundary array is available.

### What `PlotWedge` expects vs what PP sets for `clathType='top'`

`eLid_m` in PP's thermal model is the **ice Ih conductive lid**, not the clathrate
layer thickness. `PlotWedge` (line 1157) uses it as the clathrate-lid width for
`clathType='top'` — this is intentional PP convention, not a bug. For strongly
convective parameter sets (e.g. Titan Test46 all-ice), `Dconv_m + deltaTBL_m ≈
dzIceI_km * 1e3`, leaving essentially no conductive ice Ih band; `iceIcond` pixels
will be absent from the rendered image. This is physically correct.

The clathrate lid (`clathCond` color) renders with width `eLid_m / rMax_km`. For
Titan (R≈2575 km) with a ~3–5 km clathrate lid this is <0.2% of the wedge radius —
sub-pixel at standard figure size. For Europa (R≈1561 km) the same lid is ~0.3% and
is marginally visible. No minimum-width floor is applied; the rendering is an
accurate representation of the physical layer proportions.

---

## 10. References

- Karamanis et al. 2022, *Accelerating astronomical and cosmological inference with preconditioned Monte Carlo*, MNRAS 516, 1644.
- Petricca et al. 2025, *Constraints on Titan's interior from Cassini gravity and rotation*, JGR Planets.
- Yao et al. 2014, *Convective parameterization of stagnant lid in Ice Ih layers*.
- Kalousová & Sotin 2018, *Melting in high-pressure ice layers of large ocean worlds*, GRL 45, 8096.
- Tejfel et al. (sbi package): https://www.mackelab.org/sbi/
- pocoMC: https://github.com/minaskar/pocomc

## 11. See also

- `MCMC_INFERENCE_GUIDE.md` — tutorial-style walkthrough of running Test48 / Test50.
- `KALOUSOVA_IMPLEMENTATION_GUIDE.md` — HP-ice convection physics.
- `SESSION.md` — current session state.
- `CHANGELOG.md` — what changed when.
