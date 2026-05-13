# MCMC Bayesian Inference for Titan's Rheology

**Purpose:** Constrain Titan's interior rheology parameters using observed tidal Love numbers  
**References:** Petricca et al. (2025) *Nature*; Durante et al. (2019)

## Overview

PlanetProfile integrates with **TidalPy** to perform Markov Chain Monte Carlo (MCMC) Bayesian inference of Titan's rheological parameters. Tests 41-46 explore Andrade and Maxwell rheologies with fixed-structure and hybrid-hydro (variable structure) models, comparing ocean-present and ocean-free scenarios.

### Observational Constraints (Cassini)
- **Re(k₂):** 0.608 ± 0.048
- **|Im(k₂)|:** 0.135 ± 0.035

### Key Physical Insight (Petricca et al. 2025)
A liquid ocean tidally decouples HP ice below it (ocean transmits no shear stress). With ocean present, only Ice Ih dissipates — requiring unphysically low viscosity. HP ice (III, V, VI) has greater volume and CAN account for tidal heating, but only without an ocean layer.

## Test Suite

### Test 41: Andrade No-Ocean (`Test41_mcmc_andrade_no_ocean.py`)
**Scenario:** Titan without subsurface ocean (450 km total hydrosphere)  
**Rheology:** Andrade (power-law frequency-dependent)  
**Parameter space (5D):**
- α (Andrade exponent): [0.2, 0.4]
- log₁₀(ζ) (Andrade compliance): [-2, 2]
- log₁₀(η_Ih) (Ice Ih viscosity, Pa·s): [12, 16]
- log₁₀(η_HP) (HP ice viscosity, Pa·s): [10, 18]
- log₁₀(η_sil) (Silicate viscosity, Pa·s): [18, 22]

**Outputs:**
- `andrade_no_ocean_corner.png` — Posterior parameter distributions (KDE corner plot)
- `andrade_no_ocean_k2_scatter.png` — Re(k₂) vs Im(k₂) with error ellipse, colored by silicate heating fraction
- `andrade_no_ocean_heating.png` — Per-phase tidal heating distribution across posterior
- `andrade_no_ocean_mcmc.pkl` — Full chain + heating arrays for re-analysis

### Test 42: Maxwell Ocean (`Test42_mcmc_maxwell_ocean.py`)
**Scenario:** Titan with 28 km subsurface ocean (Tᵦ = 255 K)  
**Rheology:** Maxwell (viscous relaxation, elastic at short timescales)  
**Parameter space (5D):**
- log₁₀(μ_Ih) (Ice Ih shear modulus, Pa): [9, 10]
- log₁₀(η_Ih) (Ice Ih viscosity, Pa·s): [12, 16]
- log₁₀(η_HP) (HP ice viscosity, Pa·s): [10, 18]
- log₁₀(η_sil) (Silicate viscosity, Pa·s): [18, 22]
- Tᵦ (Ocean-ice boundary temperature, K): [250, 260]

**Outputs:**
- `maxwell_ocean_corner.png`
- `maxwell_ocean_k2_scatter.png`
- `maxwell_ocean_heating.png`
- `maxwell_ocean_structure_profile.png` — Interior structure vs Tᵦ
- `maxwell_ocean_Tb_structure.png` — Layer thickness evolution with boundary temperature
- `maxwell_ocean_thickness_vs_heating.png` — Thickness-heating correlations
- `maxwell_ocean_mcmc.pkl`

### Test 43: Andrade Arrhenius No-Ocean (`Test43_mcmc_andrade_arrhenius_no_ocean.py`)
**Scenario:** No ocean, 4 km clathrate cap  
**Rheology:** Andrade + Arrhenius temperature-dependent viscosity for HP ice  
**Parameter space (5D):** Same as Test41

**Key difference:** HP ice viscosity η_HP(T) = η₀ exp[Q/R(1/T - 1/Tₘ)] with activation energy Q derived from laboratory measurements (Kalousova et al. 2018).

**Outputs:** Same as Test41 with `_arrhenius` suffix

### Test 44: Maxwell Arrhenius Ocean (`Test44_mcmc_maxwell_arrhenius_ocean.py`)
**Scenario:** Ocean present (28 km)  
**Rheology:** Maxwell + Arrhenius for HP ice  
**Parameter space (5D):** Same as Test42

**Outputs:** Same as Test42 with `_arrhenius` suffix

### Test 46: Andrade Hybrid-Hydro (`Test46_mcmc_andrade_hybrid_hydro.py`)
**Scenario:** Variable hydrosphere structure (Tb and D_hydro sampled)  
**Rheology:** Andrade with per-phase viscosity (Ih, III, V, VI, sil independent)  
**Parameter space (10D):**
- alpha (Andrade exponent): [0.2, 0.4]
- log10(zeta) (compliance ratio): [-3, 2]
- log10(eta_Ih): [12, 16]
- log10(eta_III): [10, 16]
- log10(eta_V): [10, 16]
- Tb_K (ice-ocean boundary temp): [252, 270]
- D_hydro_km (total hydrosphere): [50, 800]
- log10(eta_VI): [10, 16]
- log10(eta_sil): [18, 22]
- f_core (core mass fraction): [0.15, 0.35]

**Grid cache:** Pre-computed PlanetProfile structures on (Tb_K x D_hydro) grid (25 km spacing, 6789 points). Per-phase HP ice breakdown (III, V, VI) stored individually.

**Outputs:**
- `hybrid_hydro_andrade_corner.png` — 10D posterior
- `hybrid_hydro_andrade_k2_scatter.png` — k2 constraint match
- `hybrid_hydro_andrade_layers_vs_docean.png` — All posterior models sorted by D_ocean (stackplot + per-phase heating vs f_sil)
- `hybrid_hydro_andrade_structure_profile.png` — Wedge diagram at posterior median
- `hybrid_hydro_andrade_mcmc.pkl` — Full chain

### Test 48: Andrade Yao2014 Spherical Convection (`Test48_mcmc_andrade_yao2014.py`)
**Scenario:** Variable hydrosphere structure with Yao et al. (2014) spherical shell convection for Ice Ih  
**Rheology:** Andrade with per-phase viscosity  
**Convection model:** Yao et al. (2014) — 3D spherical shell stagnant-lid scaling (replaces DS2001 Cartesian)  
**Parameter space (10D):** Same as Test46

**Key differences from Test46:**
- Grid cache built with `SPHERICAL_CONVECTION=True` → thicker stagnant lids, lower basal heat flux
- Adds surface heat flux constraint: q_surface = 10 ± 5 mW/m² (Nimmo & Bills 2010)
- Inline `yao_heat_flux_mWm2()` couples sampled η_Ih to convective physics via analytical scaling law
- Uses Ra_fullDT for onset check, Ra_m (viscous-temperature) for heat flux scaling
- Penalizes samples where η_Ih produces non-physical heat flux (too high or no convection onset)

**Physical motivation:** DS2001 uses Cartesian scaling laws that underestimate heat flux for thick shells. Yao 2014 corrects with f-dependent 3D spherical geometry — most impactful for Titan (f≈0.96) where the Ice Ih shell fills most of the spherical gap.

**Grid cache:** `titan_yao2014_hybrid_hydro_grid_cache.pkl` (separate from DS2001 cache)

**Structure config:** `PPTest48.py` (deepcopy of PPTest46 with `SPHERICAL_CONVECTION=True`)

**Expected results vs Test46:**
- Thicker lids: ~33-45 km (Yao) vs ~6-10 km (DS2001) for same (Tb, D_hydro)
- Lower convective heat flux (~7% reduction)
- η_Ih posterior bounded below ~10¹³ Pa·s (Yao onset) and above ~10¹⁶ Pa·s (no convection)
- Im(k₂) posterior slightly shifted due to different lid geometry

**Outputs:**
- `yao2014_andrade_corner.png` — 10D posterior
- `yao2014_andrade_k2_scatter.png` — k2 constraint match
- `yao2014_andrade_mcmc.pkl` — Full chain

### Test 46 All-Ice Variant (`Test46_mcmc_allice.py`)
**Scenario:** No ocean, fixed structure (D_hsphere=493.7 km)  
**Rheology:** Andrade with per-phase viscosity (7D: alpha, zeta, eta_Ih, eta_III, eta_V, eta_VI, eta_sil)  
**Convection:** DS2001 for Ice Ih, Kalousova for HP, Arrhenius viscosity.  
**Key result:** Petricca-compatible model where HP ice can dissipate without ocean decoupling.

**Outputs:** `allice_andrade_corner.png`, `allice_andrade_k2_heating.png`

### Test 50: Andrade No-Ocean Titan with Yao 2014 (`Test50_mcmc_andrade_noocean_yao2014.py`)
**Scenario:** All-ice (no-ocean) Titan with Yao 2014 spherical-shell Ih convection.
**Rheology:** Andrade with per-phase viscosity + basal temperature (8D).
**Convection:** Yao for Ice Ih (replaces DS2001), Kalousova for HP, Arrhenius.
**Clathrate cap:** 2 km insulating layer at surface (Fourier-matched to ~40 mW/m² interior heat flux when clathrate is the full stagnant lid; Yao self-consistently places a ~15 km eLid encompassing the clathrate plus an Ih conductive sub-layer).

**Motivation:** Tests Petricca's central claim that no-ocean structures are required for HP ice dissipation to drive the observed |Im(k₂)|. The Test48 hybrid chain couldn't sample no-ocean structures fairly (0/4434 samples at D_hydro < 200 km due to CMR2 + ρ_sil-floor degeneracy) — this test samples the no-ocean branch cleanly via a fixed structure, following Petricca's methodology.

**Structure config:** `PPTest50.py` — `SPHERICAL_CONVECTION=True` added on top of `PPTest46_allice.py`; clathrate cap set to 2 km.  Tb default = `TtripleIh_III_L_K − 0.2 K` (ε sized to grid resolution; see PPTest50 comment).

**Structure handling (Option 2 — grid + interpolate):** cache-build produces 9 structures over Tb ∈ [249.0, 250.965] K.  Forward model linearly interpolates arrays between bracketing grid points per MCMC sample, capturing eLid/TBL/PbI drift.  Consistency check at build time asserts identical `changeIndices` / `layer_types` / `region_phases` across grid.

**Parameter priors (8D):**
- α ∈ [0.15, 0.45], log₁₀ζ ∈ [−3, 2]
- log₁₀η_{Ih, III, V, VI} ∈ [10, 16] — broadened to include Petricca low-η regime across all ice phases
- log₁₀η_sil ∈ [18, 22]
- T_b ∈ [249.0, 250.965] K — narrow band below the Ih–III–L triple point; physically equivalent to marginalizing over realistic solute-driven triple-point depressions (e.g. NaCl up to ~15 ppt).

**No-ocean safeguard:** forward model rejects any Ih cell whose (interpolated) T meets or exceeds the linearized Ih–L melt curve Tm(P) = 273.16 − 0.068·P_MPa (with 0.1 K margin).  Keeps the no-ocean assumption internally consistent across the sampled Tb band.

**Expected results:** Comparison against Test48 (ocean-bearing) posterior will reveal whether HP ices become dissipation-competitive when ocean decoupling is removed.

**Outputs:** `allice_yao2014_andrade_*.png` (same plot suite as Test48 via the mcmc_plots toolkit).

## Running the Tests

### Prerequisites
```bash
pip install TidalPy pocoMC corner seaborn
```

### Execution
Each test is standalone:
```bash
python PlanetProfile/Test/Test41_mcmc_andrade_no_ocean.py   # ~2-4 hours
python PlanetProfile/Test/Test42_mcmc_maxwell_ocean.py      # ~3-5 hours
python PlanetProfile/Test/Test43_mcmc_andrade_arrhenius_no_ocean.py
python PlanetProfile/Test/Test44_mcmc_maxwell_arrhenius_ocean.py
python PlanetProfile/Test/Test48_mcmc_andrade_yao2014.py    # ~3-5 hours (requires Yao2014 grid)
python PlanetProfile/Test/Test50_mcmc_andrade_noocean_yao2014.py   # ~30 min – 1 hr (9-point Tb grid build + MCMC)
```

### Re-plotting Saved Results
```bash
python PlanetProfile/Test/replot_mcmc.py
```
Regenerates all figures from saved `.pkl` files without re-running MCMC.

## Key Findings

### Ocean vs No-Ocean Scenarios
- **Ocean case (Test42):** Re(k₂) posterior centered at ~0.58-0.60, consistent with observations. Heating dominated by HP ice dissipation (70-80% Ice VI).
- **No-ocean case (Test41):** Re(k₂) posterior ~0.60-0.62, matching observations. Heating concentrated in Ice Ih (60-70% due to thin stagnant lid from 4 km clathrate cap).

### Arrhenius Temperature Dependence
- **With Arrhenius (Tests 43-44):** HP ice viscosity reduces by 1-2 orders of magnitude near melting point (~272 K for Ice VI). Posteriors tighten around η_HP ~ 10¹²-10¹³ Pa·s.
- **Without Arrhenius (Tests 41-42):** Constant viscosity assumption yields broader posteriors, η_HP ~ 10¹¹-10¹⁴ Pa·s.

### Rheology Comparison
- **Andrade:** Better fit to imaginary k₂ (Im(k₂) captures frequency-dependent dissipation). α ~ 0.3 (median), consistent with laboratory measurements.
- **Maxwell:** Simpler model, fewer parameters. Adequate for Re(k₂) but underpredicts Im(k₂) by 10-20%.

### Silicate Viscosity Constraint
All scenarios converge on log₁₀(η_sil) ~ 19-20 Pa·s, consistent with partially hydrated, warm silicate mantle. Values < 10¹⁸ Pa·s imply unrealistic partial melting; values > 10²¹ Pa·s prevent sufficient tidal heating.

## Outputs and Interpretation

### Corner Plots
- **1D marginals:** Median ± 1σ credible intervals for each parameter
- **2D contours:** 1σ and 2σ confidence regions (68% and 95%)
- **Color scheme:** Consistent phase colors across all plots

### k₂ Scatter Plots
- **Color scale:** Silicate heating fraction (0 = all ice, 1 = all silicate)
- **Error ellipse:** Observational constraints from Petricca et al. 2025
- **Interpretation:** Points within ellipse match observations; clustering shows degeneracies

### Heating Distribution Plots
- **Stacked bar chart:** Per-phase tidal heating (Fe, Sil, VI, V, III, Ocean, Ih, Clath)
- **Median + percentiles:** 16th/50th/84th percentiles across posterior
- **Key insight:** Ocean case → HP ice dissipation; no-ocean case → Ice Ih dissipation

### Structure Profile Plots (Ocean cases only)
- **Left panel:** Layer thickness vs Tᵦ
- **Right panel:** Heating rate vs Tᵦ
- **Correlation analysis:** Thicker Ice VI → higher total dissipation → higher Im(k₂)

## File Structure
```
PlanetProfile/Test/
├── Test41_mcmc_andrade_no_ocean.py          # Andrade no-ocean MCMC (5D)
├── Test42_mcmc_maxwell_ocean.py             # Maxwell ocean MCMC (5D)
├── Test43_mcmc_andrade_arrhenius_no_ocean.py
├── Test44_mcmc_maxwell_arrhenius_ocean.py
├── Test46_mcmc_andrade_hybrid_hydro.py      # Andrade hybrid-hydro (10D, DS2001)
├── Test46_mcmc_allice.py                    # All-ice variant (no ocean, 7D)
├── Test48_mcmc_andrade_yao2014.py           # Andrade hybrid-hydro (10D, Yao2014) — uses mcmc_plots toolkit
├── Test50_mcmc_andrade_noocean_yao2014.py   # Andrade no-ocean Yao+Kalousova (8D: +Tb)
├── PPTest41.py ... PPTest44.py              # Structural configs
├── PPTest46_allice.py                       # All-ice structural config
├── PPTest48.py                              # Titan Yao2014 config (SPHERICAL_CONVECTION)
├── PPTest50.py                              # Titan no-ocean Yao+Kalousova + 2 km clathrate
├── Test40_maxwell_sweep.py                  # Parameter sweep (precursor)
├── replot_mcmc.py                           # Regenerate figures from .pkl
└── mcmc_results/                            # Output directory
    ├── titan_maxwell_hybrid_hydro_grid_cache.pkl  # DS2001 grid cache
    ├── titan_yao2014_hybrid_hydro_grid_cache.pkl  # Yao2014 grid cache
    ├── hybrid_hydro_andrade_mcmc.pkl
    ├── hybrid_hydro_andrade_*.png
    ├── yao2014_andrade_mcmc.pkl
    ├── yao2014_andrade_*.png
    ├── allice_andrade_mcmc_results.pkl
    ├── allice_andrade_*.png
    └── [andrade_no_ocean_*, maxwell_ocean_*, etc.]

PlanetProfile/Inference/
├── hybrid_structure_cache.py                # Grid cache builder (DS2001 or Yao2014 via convection_model kwarg)
├── structure_cache.py                       # Fixed-structure cache for Test46_allice
├── inference_core.py                        # Shared MCMC/SBI dispatch logic
├── mcmc_common.py                           # Body-agnostic MCMC helpers (NEW)
└── mcmc_plots.py                            # Body-agnostic plot library (NEW)
```

## Inference Toolkit (`mcmc_common.py`, `mcmc_plots.py`)

Extracted from Test48 to allow reuse across bodies (Titan, Europa, Ganymede, ...) with consistent forward-model physics and plotting.

### `mcmc_common.py` — body-agnostic helpers

| Function | Purpose |
|---|---|
| `apply_arrhenius_ih(eta_mod, phases, ci, n_layers, T_K_profile, Tb_K, E_act_J_mol=60e3, R_gas=8.314462, ih_phase_id=1)` | Scale Ice Ih viscosity by Arrhenius law where sampled η is the *basal* value at Tb (Yao convention). |
| `split_silicate_core(...)` | Insert inner-core boundary at f_core · r_sil_top, apply two-layer density + elastic moduli, extend layering. |
| `build_andrade_shear_bulk(region_phases, alpha, log10_zeta)` | Per-layer Andrade shear + Elastic bulk (ocean/clathrate → Elastic). |
| `build_maxwell_shear_bulk(region_phases)` | Maxwell variant. |
| `compute_per_phase_heating(result, data)` | Integrate TidalPy volumetric dissipation by phase. Returns dict {phase name → W}. |
| `run_pocomc_sampler(prior, log_like_fn, n_effective=500, random_state=42)` | pocoMC wrapper with timing + logging. |
| `evaluate_posterior(samples, forward_fn, n_eval=500, random_state=42)` | Re-evaluate forward model on posterior subset for heating breakdowns. |
| `gaussian_chi2_terms(observed, predicted)` | Sum ((pred - μ)/σ)² over matched keys. |

### `mcmc_plots.py` — 9 body-agnostic plot functions

| Function | Figure |
|---|---|
| `plot_corner(samples, labels, title, output_path, quantiles=(0.16,0.5,0.84))` | Posterior corner with quantile titles. |
| `plot_k2_scatter_by(k2_results, color_values, colorbar_label, obs_re, obs_im, ...)` | k₂ Re/Im scatter, arbitrary per-sample colouring + Petricca 1σ/2σ ellipses. |
| `plot_ice_comparison(k2_results, heating_results, d_ocean_eval, obs_re, obs_im, ...)` | Two-panel Ih-vs-HP dominance diagnostics. |
| `plot_heating_vs_parameters(eval_samples, heating_results, param_labels, extra_xvals, ..., cumulative_bar=True, eval_d_ocean=...)` | 2×5 heating-vs-parameter scatter + optional full-width cumulative stacked-fraction bar (sorted by (f_sil, D_ocean)). |
| `plot_mass_cmr2_diagnostics(d_hydro_values, mtot_results, cmr2_results, obs_*)` | Mass and CMR² vs D_hydro with observation bands. |
| `plot_cmr2_surface(tb_vals, d_vals, grid_cache, output_path)` | pcolormesh of CMR² across (Tb, D_hydro). |
| `plot_tb_structure(tb_vals, d_vals, grid_cache, samples, output_path)` | Ice Ih thickness vs Tb + posterior Tb histogram. |
| `plot_layers_vs_docean(samples, eval_idx, grid_cache, ..., R_body_km, body_name='Titan', param_indices=None, equil_heating_GW=None)` | 3-panel: D_ocean density / cumulative-thickness stackplot / per-phase heating, sorted by (f_sil ASC, D_ocean ASC). |
| `plot_structure_wedge(samples, eval_idx, grid_cache, ..., R_body_km, body_name='Titan', param_indices=None)` | Wedge diagram of posterior median interior with 5/95-percentile arcs. |

`param_indices` defaults to Test48 layout `{'Tb': 5, 'D_hydro': 6, 'f_core': 9}` — override when adapting to a body with a different parameter ordering.

### Usage pattern (new body)

```python
from PlanetProfile.Inference import mcmc_common as mc, mcmc_plots as mp

# 1. Forward model uses helpers:
mc.apply_arrhenius_ih(eta_mod, phases, ci, n_layers, data['T_K'], Tb_K)
mc.split_silicate_core(r_m, rho_mod, K_Pa_mod, mu_Pa_mod, phases, ci, ...)
shear, bulk = mc.build_andrade_shear_bulk(region_phases, alpha, log10_zeta)
# ...TidalPy call...
heating = mc.compute_per_phase_heating(result, data)

# 2. MCMC run:
samples, log_likes, sampler = mc.run_pocomc_sampler(prior, log_like_fn)
eval_idx, eval_results = mc.evaluate_posterior(samples, forward_fn)

# 3. Plots — body-agnostic:
mp.plot_corner(samples, PARAM_LABELS, title, out_path)
mp.plot_layers_vs_docean(samples, eval_idx, grid_cache, tb_vals, d_vals,
                         heating_results, out_path,
                         R_body_km=EUROPA_R_M / 1e3, body_name='Europa',
                         param_indices={'Tb': 5, 'D_hydro': 6, 'f_core': 9})
```

## Technical Details

### MCMC Sampler
- **Algorithm:** pocoMC (Preconditioned Monte Carlo)
- **Target:** n_eff = 500 effective samples
- **Convergence:** Automatically stops when ESS > n_eff
- **Warm-up:** Dynamically determined (typically 1000-2000 iterations)

### Forward Model
1. **Build structure:** Run PlanetProfile once with PPTest configuration
2. **Cache layers:** Store P, T, ρ, r profiles (fixed across MCMC)
3. **TidalPy solver:** Radial solver with per-layer rheology
4. **Heating calculation:** `calc_radial_volumetric_tidal_heating_from_rs_solution()`
5. **Log-likelihood:** Gaussian with observational covariance

### Posterior Re-evaluation
After MCMC completes, 500 samples are re-evaluated with full structure + heating:
- Enables phase-resolved heating plots
- Allows structure-parameter correlation analysis
- Preserves computational efficiency (cached structures)

## Physics Notes

### Tidal Decoupling (Petricca et al. 2025)
- Ocean layer transmits no shear stress → HP ice below is tidally decoupled
- With ocean: only Ice Ih dissipates → requires eta_Ih ~ 10^13 (possibly unphysical heating rates)
- Without ocean: HP ice volume dominates dissipation → steady-state heating achievable
- Two-phase convection models justify low ice viscosity (porous medium transport)
- Thermal equilibrium constraint: Ice Ih alone can't sustain >10 TW without melting

### TidalPy Heating vs Simple Formula
- Simple: E_dot = (21/2)|Im(k2)| n^5 R^5 e^2 / G
- TidalPy volumetric integration gives ~1.74x the simple formula (additional forcing modes)
- Per-phase breakdown requires full radial solver evaluation

## References

- Petricca et al. (2025) *Nature* — Cassini gravity constraints, ocean-free Titan paradigm
- Yao, C., Deschamps, F., Lowman, J. P., Sanchez-Valle, C., & Tackley, P. J. (2014) *JGR Planets* — Stagnant-lid convection in spherical shell geometry
- Kalousova & Sotin (2018) *GRL* — HP ice convection with partial melt prediction
- Nimmo & Bills (2010) *Icarus* — Titan tidal heating and heat flux estimates
- Renaud & Henning (2018) *ApJ* — TidalPy framework for tidal heating calculations
- Andrade (1910) — Power-law frequency-dependent rheology
- Maxwell (1867) — Viscoelastic relaxation model

## Key Findings — Test48 Path B (Yao + Arrhenius + loose structure)

Five progressive MCMC runs (2026-05-05 to 2026-05-07) converged on the B-path result:

| Run | Change | Re(k₂) | \|Im(k₂)\| | χ to obs | f_Ih (median) | 1σ hit |
|-----|--------|--------|-----------|----------|---------------|--------|
| pre-fix (bug) | uniform η_Ih everywhere | 0.641 | 0.109 | 1.00σ | 100% | 26% |
| A1 | Arrhenius Ice Ih (basal ref) | 0.650 | 0.096 | 1.42σ | 100% | 19% |
| A1 + B2 | HP prior [10,18] → [12,16] | 0.666 | 0.086 | 1.86σ | 100% | 11% |
| A1 + B2 + no-qs | q_surface constraint removed | 0.660 | 0.090 | 1.68σ | 100% | 11% |
| **B (loose struct)** | **CMR2 err 0.001→0.005, ρ_sil floor 2000→1800** | **0.633** | **0.104** | **1.01σ** | **100%** | **28%** |

**Interpretation:** Under physically correct Yao + Arrhenius + Andrade, Titan tidal heating is dominated by ~100 km of warm basal Ice Ih at the Maxwell dissipation peak (η ≈ 10¹⁵ Pa·s, μ ≈ 3.5 GPa, ωτ ≈ 1 at the tidal frequency). This result is robust to prior choice — a D4 forward-model sweep (2025 rheology combinations at fixed structure) found zero HP-dominated solutions. The observed |Im(k₂)| = 0.135 sits on the edge of the reachable (Re, Im) curve; best samples reach 0.06σ but the posterior median is broader.

**Active ingredient:** CMR2 relaxation (0.001 → 0.005) was the key to improving fit — it opened the D_hydro × (ρ_sil, ρ_core, f_core) degeneracy that had pinned the chain at D_hydro ≈ 460 km.

## Next Steps

### Pending science
- **Test50 MCMC** — 8D MCMC script ready (rheology + Tb_K); launch pending. Comparison against Test48 posterior will reveal whether HP ices become dissipation-competitive when ocean decoupling is removed.
- **Test49: Yao + 4 km clathrate + ocean** — perturbation check of Path B result (plan documented; grid rebuild required).
- **Test51: Europa MCMC** — apply toolkit to Europa with variable (ρ_sil, ρ_core, f_core) and Galileo/Juno k₂ constraints.

### Toolkit continuation
- **Port Test48 forward_model** to `mcmc_common.apply_arrhenius_ih`, `split_silicate_core`, `compute_per_phase_heating` — currently only `make_plots` uses the toolkit.
- **Port Test50 plotting** to `mcmc_plots` — current Test50 has inline plots from initial port.
- **BodyConfig abstraction** — centralize OBS, R/M, n, e constants per body so Test51+ are thin wrappers.
- **PlanetProfileApp integration** — live progress, multi-run comparison, body selector.

### Advanced inference options
- **SBI (Simulation-Based Inference):** `sbi_runner.py` for continuous parameter sampling. Infrastructure in place (`inference_core.py` has dispatch).
- **Petricca constraint comparison:** `--petricca` flag for Re(k₂)=0.133 constraint set (Durante et al. 2019 phase-lag interpretation).
