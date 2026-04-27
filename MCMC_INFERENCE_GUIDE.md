# MCMC Bayesian Inference for Titan's Rheology

**Status:** Complete (April 2026)  
**Purpose:** Constrain Titan's interior rheology parameters using observed tidal Love numbers
**Reference:** Petricca et al. (2025) *Nature*

## Overview

PlanetProfile integrates with **TidalPy** to perform Markov Chain Monte Carlo (MCMC) Bayesian inference of Titan's rheological parameters. Tests 41-44 explore Andrade and Maxwell rheologies with and without Arrhenius temperature-dependent viscosity, comparing ocean-present and ocean-free scenarios.

### Observational Constraints (Petricca et al. 2025)
- **Re(k₂):** 0.608 ± 0.048
- **|Im(k₂)|:** 0.135 ± 0.035

These constraints from Cassini gravity field measurements challenge the traditional ocean-world paradigm and suggest HP ice layers with reduced viscosities.

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
├── Test41_mcmc_andrade_no_ocean.py          # Andrade no-ocean MCMC
├── Test42_mcmc_maxwell_ocean.py             # Maxwell ocean MCMC
├── Test43_mcmc_andrade_arrhenius_no_ocean.py
├── Test44_mcmc_maxwell_arrhenius_ocean.py
├── PPTest41.py                              # Structural config (no-ocean)
├── PPTest42.py                              # Structural config (ocean)
├── PPTest43.py                              # Structural config (Arrhenius no-ocean)
├── PPTest44.py                              # Structural config (Arrhenius ocean)
├── Test40_maxwell_sweep.py                  # Parameter sweep (precursor to MCMC)
├── replot_mcmc.py                           # Regenerate all figures from .pkl
└── mcmc_results/                            # Output directory
    ├── andrade_no_ocean_corner.png
    ├── andrade_no_ocean_k2_scatter.png
    ├── andrade_no_ocean_heating.png
    ├── andrade_no_ocean_mcmc.pkl
    ├── maxwell_ocean_corner.png
    ├── maxwell_ocean_k2_scatter.png
    ├── maxwell_ocean_heating.png
    ├── maxwell_ocean_structure_profile.png
    ├── maxwell_ocean_Tb_structure.png
    ├── maxwell_ocean_thickness_vs_heating.png
    ├── maxwell_ocean_mcmc.pkl
    └── [andrade_arrhenius_no_ocean_*, maxwell_arrhenius_ocean_*]
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

## Manuscript Integration

Results documented in:
- `TITAN_RESULTS_SECTION_FINAL.md` — Results section for HP ice convection paper
- `~/Dropbox/planetprofile-genai-manuscript/PlanetProfile-genai/main.tex` — Draft manuscript

Key manuscript points:
- MCMC posteriors validate Petricca et al. 2025 no-ocean model
- Fluid-bearing HP ice viscosities (10¹²-10¹³ Pa·s) required for both scenarios
- Clathrate cap reduces Ice Ih stagnant lid by 93%, enabling efficient tidal dissipation

## References

- Petricca et al. (2025) *Nature* — Cassini gravity constraints, ocean-free Titan paradigm
- Kalousova & Sotin (2018) *GRL* — HP ice convection with partial melt prediction
- Renaud & Henning (2018) *ApJ* — TidalPy framework for tidal heating calculations
- Andrade (1910) — Power-law frequency-dependent rheology
- Maxwell (1867) — Viscoelastic relaxation model

## Future Work

- **3D tidal heating:** Extend to lateral heterogeneities (TidalPy 3D solver)
- **Joint inversion:** Combine k₂, moment of inertia, and magnetic induction
- **Time-dependent evolution:** Couple with thermal-orbital evolution models
- **Exoplanet application:** Scale to super-Earths and mini-Neptunes with thick ice shells
