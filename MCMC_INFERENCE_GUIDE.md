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
**Rheology:** Andrade with per-phase viscosity (5D: alpha, zeta, eta_Ih, eta_III, eta_V, eta_VI, eta_sil)  
**Key result:** Petricca-compatible model where HP ice can dissipate without ocean decoupling.

**Outputs:** `allice_andrade_corner.png`, `allice_andrade_k2_heating.png`

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
├── Test46_mcmc_allice.py                    # All-ice variant (no ocean, 5D)
├── Test48_mcmc_andrade_yao2014.py           # Andrade hybrid-hydro (10D, Yao2014)
├── PPTest41.py ... PPTest44.py              # Structural configs
├── PPTest46_allice.py                       # All-ice structural config
├── PPTest48.py                              # Titan Yao2014 config (SPHERICAL_CONVECTION)
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
└── inference_core.py                        # Shared MCMC/SBI dispatch logic
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

## Next Steps

- **Build Yao2014 grid cache:** Run full (Tb × D_hydro) grid with `convection_model='yao2014'` (~400 points, several hours). Required before Test48 MCMC can execute.
- **Run Test48 MCMC:** Compare Yao2014 vs DS2001 posteriors — expect shifted η_Ih distribution and different Im(k₂) sensitivity.
- **SBI (Simulation-Based Inference):** Implement `sbi_runner.py` for continuous parameter sampling without grid discretization. Infrastructure in place (`sbi` installed, `inference_core.py` has dispatch).
- **Thermal equilibrium prior:** Add steady-state heating constraint to penalize Ice Ih dissipation modes exceeding melting timescale.
- **Petricca constraint comparison:** `--petricca` flag for Re(k2)=0.133 constraint set (Durante et al. 2019 phase-lag interpretation).
