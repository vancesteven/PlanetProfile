# 3D Ice Thickness Modes in PlanetProfile

This document describes the three modes for specifying ice shell thickness in 3D lateral structure calculations.

## Overview

PlanetProfile supports three modes for determining ice thickness distribution across the surface grid:

| Mode | Flag/Setting | Use Case | Physics |
|------|-------------|----------|---------|
| **Equilibrium** | `DO_EQUILIBRIUM_ICE=True` | Scientific studies | Self-consistent from heat balance |
| **Prescribed** | `dIce_Cpq_km` set | Exploratory, observational | User-specified pattern |
| **Uniform** | Neither set | Testing, 1D-like | Constant thickness |

## Mode 1: Equilibrium Ice Thickness (RECOMMENDED)

### When to Use
- **Scientific publications**: When you need physically self-consistent ice thickness
- **Thermal-tidal coupling studies**: Ice thickness affects tidal heating, which affects temperature, which affects viscosity, which affects tidal heating
- **Habitability studies**: Self-consistent ice-ocean interface topography
- **Predictive modeling**: No observational constraints on ice thickness available

### How It Works
Solves the steady-state heat balance equation at each grid pixel:

```
k * (Tb - Tsurf) / d_ice = q_basal + H_tidal(pixel) * d_ice
```

where:
- `k`: Ice thermal conductivity (W/m/K)
- `Tb`, `Tsurf`: Basal and surface temperatures (K)
- `d_ice`: Ice shell thickness (m) — **this is what we solve for**
- `q_basal`: Heat flux from silicate mantle (W/m²)
- `H_tidal(pixel)`: Volumetric tidal dissipation rate (W/m³), which depends on local ice viscosity

The solver iterates:
1. Start with uniform thickness from reference 1D model
2. Run lateral column models → compute thermal profiles
3. Compute 3D tidal heating → H_tidal depends on viscosity(T)
4. Solve heat balance → new thickness distribution
5. Repeat steps 2-4 until thickness converges

### Physics Background
Based on Tobie et al. (2003, JGR doi:10.1029/2003JE002099):
- Conductive heat transport through ice shell
- Volumetric tidal dissipation in ice (Maxwell or Andrade rheology)
- Basal heat flux from silicate mantle (radiogenic + tidal)
- Self-consistent coupling between thickness, temperature, and heating

**Expected behavior (Europa example)**:
- Thicker ice at sub/anti-Jovian points (low tidal dissipation)
- Thinner ice at poles and mid-latitudes (high dissipation)
- ~5 km peak-to-peak variation for ~20-25 km mean shell
- Nearly uniform surface heat flux (~35-40 mW/m²)

### Configuration
```python
# Enable 3D lateral structure
Planet.Lateral.DO_3D = True
Planet.Lateral.gridType = 'healpix'
Planet.Lateral.nSide = 4  # 192 pixels

# Enable 3D tidal heating (REQUIRED for equilibrium mode)
Planet.Lateral.DO_TIDAL_3D = True

# Enable equilibrium ice thickness solver
Planet.Lateral.DO_EQUILIBRIUM_ICE = True
Planet.Lateral.equilibriumTol_m = 100.0     # Convergence tolerance (m)
Planet.Lateral.equilibriumMaxIter = 10      # Max iterations
Planet.Lateral.kThermIce_WmK = 2.3          # Ice thermal conductivity

# Optional: Override basal heat flux (otherwise computed from Sil properties)
# Planet.Lateral.qBasal_Wm2 = 0.007  # W/m² (7 mW/m² for Europa)

# Do NOT set dIce_Cpq_km — equilibrium solver computes thickness

# Required orbital parameters
Planet.Bulk.eccentricity = 0.0094
Planet.Bulk.meanMotion_radps = 2 * np.pi / (3.551181 * 86400)  # Europa
```

### Outputs
After the run, `Planet.Lateral` contains:
- `dIce_m`: Equilibrium ice thickness at each pixel (m)
- `qSurf_Wm2`: Surface heat flux at each pixel (W/m²)
- `equilibriumIterations`: Number of iterations to converge
- `equilibriumResidual_m`: Final max residual (m)

## Mode 2: Prescribed Ice Thickness (Spherical Harmonics)

### When to Use
- **Exploratory studies**: Testing sensitivity to ice thickness patterns
- **Observational constraints**: Matching radar or gravity-derived thickness
- **Comparison studies**: Comparing prescribed vs equilibrium patterns
- **Quick prototyping**: Faster than equilibrium mode (no iteration)

### How It Works
Ice thickness is specified as a sum of spherical harmonics:

```
d_ice(θ, φ) = Σ_p Σ_q [ C_pq * Y_pq^cos(θ, φ) + S_pq * Y_pq^sin(θ, φ) ]
```

where:
- `p`: Degree (0 to pMax)
- `q`: Order (0 to p)
- `C_pq`, `S_pq`: Cosine and sine coefficients (km)
- `Y_pq`: Real spherical harmonics

**Common patterns**:
- `C_00`: Mean thickness (monopole)
- `C_20`: Pole-equator variation (oblateness)
- `C_22, S_22`: 4-fold longitudinal pattern (sub/anti-Jovian + leading/trailing)

### Configuration
```python
# Enable 3D lateral structure
Planet.Lateral.DO_3D = True
Planet.Lateral.gridType = 'healpix'
Planet.Lateral.nSide = 4

# Enable 3D tidal heating (optional but recommended)
Planet.Lateral.DO_TIDAL_3D = True

# Disable equilibrium solver
Planet.Lateral.DO_EQUILIBRIUM_ICE = False

# Prescribe ice thickness via spherical harmonics
Planet.Lateral.dIce_pMax = 2  # Maximum degree

Planet.Lateral.dIce_Cpq_km = np.array([
    [29.0,  0.0,  0.0],   # C_00 = 29 km mean
    [ 0.0,  0.0,  0.0],   # p=1: zero (synchronous rotator)
    [-3.5,  0.0, -1.5],   # C_20 = -3.5 km (polar thinning)
                          # C_22 = -1.5 km (sub/anti-Jovian thinning)
])
Planet.Lateral.dIce_Spq_km = np.array([
    [0.0,  0.0,  0.0],
    [0.0,  0.0,  0.0],
    [0.0,  0.0, -0.7],    # S_22 = -0.7 km (phase offset)
])
```

### Physical Interpretation of Coefficients

#### Europa Example (Tobie et al. 2003, 2005)
```python
# Mean thickness (C_00)
C_00 = 29.0  # km — Bolton et al. (2025) Juno MWR constraint

# Polar thinning (C_20 < 0)
C_20 = -3.5  # km — Thinner at poles due to higher tidal dissipation
             # Physical origin: Combination of tidal strain pattern
             # (Ojakangas & Stevenson 1989) and ice viscosity structure

# Sub/anti-Jovian thinning (C_22, S_22)
C_22 = -1.5  # km — 4-fold longitudinal pattern
S_22 = -0.7  # km — Phase offset from librational forcing
             # Physical origin: Eccentricity tide pattern
             # (Roberts & Nimmo 2008, Soderlund et al. 2023)
```

### Outputs
After the run, `Planet.Lateral` contains:
- `dIce_m`: Prescribed ice thickness at each pixel (m)
- `HtidalIce_Wm3`: Tidal heating given the prescribed thickness (W/m³)

## Mode 3: Uniform Ice Thickness (Fallback)

### When to Use
- **Testing**: Debugging 3D workflow with simplest case
- **1D-like runs**: Want 3D grid but no lateral variation
- **Baseline comparison**: Reference for comparing against equilibrium or prescribed

### How It Works
All pixels get the same thickness from the reference 1D model:

```python
d_ice(θ, φ) = Planet.zb_km * 1000  # meters
```

This is automatically selected if neither equilibrium mode nor prescribed coefficients are set.

### Configuration
```python
# Enable 3D lateral structure
Planet.Lateral.DO_3D = True
Planet.Lateral.gridType = 'healpix'
Planet.Lateral.nSide = 4

# Do NOT set DO_EQUILIBRIUM_ICE or dIce_Cpq_km
# Uniform mode is the default fallback
```

## Mode Selection Priority

The modes are checked in this order (first match wins):

1. **Equilibrium**: If `DO_EQUILIBRIUM_ICE=True`
2. **Prescribed SH**: If `dIce_Cpq_km is not None`
3. **Prescribed function**: If `dIce_func is not None`
4. **Uniform**: Otherwise (default fallback)

**Warning**: If you set both `DO_EQUILIBRIUM_ICE=True` and `dIce_Cpq_km`, the prescribed coefficients are **ignored** and a warning is logged.

## Validation and Error Checking

### Equilibrium Mode Requirements
The following parameters are **required** for equilibrium mode:
- `Planet.Lateral.DO_TIDAL_3D = True` (must compute tidal heating)
- `Planet.Bulk.eccentricity` (orbital eccentricity)
- `Planet.Bulk.meanMotion_radps` (mean motion = 2π / orbital period)
- `Planet.Bulk.Tsurf_K` (surface temperature)
- `Planet.Bulk.Tb_K` (basal temperature, used as initial guess)

If any are missing, `InitLateralStructure()` raises a `ValueError` with a helpful message.

## Examples

### Example 1: Equilibrium Mode (Recommended)
```python
from PlanetProfile import GetConfig
from PlanetProfile.Main import PlanetProfile
from PlanetProfile.Utilities.defineStructs import PlanetStruct
import numpy as np

# Create Europa configuration
Planet = PlanetStruct('Europa')
Planet.Bulk.R_m = 1560.8e3
Planet.Bulk.M_kg = 4.800e22
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Tb_K = 268.0
Planet.Bulk.eccentricity = 0.0094
Planet.Bulk.meanMotion_radps = 2 * np.pi / (3.551181 * 86400)

# Configure 3D with equilibrium ice thickness
Planet.Lateral.DO_3D = True
Planet.Lateral.DO_TIDAL_3D = True
Planet.Lateral.DO_EQUILIBRIUM_ICE = True
Planet.Lateral.nSide = 4

# Run PlanetProfile
Params = GetConfig.Params
Params.CALC_NEW = True
Planet, Params = PlanetProfile(Planet, Params)

# Results
print(f"Mean ice thickness: {np.mean(Planet.Lateral.dIce_m)/1e3:.2f} km")
print(f"Converged in {Planet.Lateral.equilibriumIterations} iterations")
```

### Example 2: Prescribed Mode
```python
# Same setup as above, but with prescribed thickness
Planet.Lateral.DO_EQUILIBRIUM_ICE = False  # Turn off equilibrium

Planet.Lateral.dIce_pMax = 2
Planet.Lateral.dIce_Cpq_km = np.array([
    [29.0,  0.0,  0.0],
    [ 0.0,  0.0,  0.0],
    [-3.5,  0.0, -1.5],
])
Planet.Lateral.dIce_Spq_km = np.array([
    [0.0,  0.0,  0.0],
    [0.0,  0.0,  0.0],
    [0.0,  0.0, -0.7],
])

# Run as before
Planet, Params = PlanetProfile(Planet, Params)
```

### Example 3: Switching Between Modes
```python
# Run both modes for comparison
def run_both_modes(base_config):
    # Equilibrium mode
    Planet_eq = deepcopy(base_config)
    Planet_eq.Lateral.DO_EQUILIBRIUM_ICE = True
    # Remove any prescribed coefficients
    Planet_eq.Lateral.dIce_Cpq_km = None
    Planet_eq, _ = PlanetProfile(Planet_eq, Params)
    
    # Prescribed mode
    Planet_presc = deepcopy(base_config)
    Planet_presc.Lateral.DO_EQUILIBRIUM_ICE = False
    Planet_presc.Lateral.dIce_Cpq_km = np.array([[29, 0, 0], [0, 0, 0], [-3.5, 0, -1.5]])
    Planet_presc.Lateral.dIce_Spq_km = np.array([[0, 0, 0], [0, 0, 0], [0, 0, -0.7]])
    Planet_presc, _ = PlanetProfile(Planet_presc, Params)
    
    return Planet_eq, Planet_presc
```

## Performance Considerations

| Mode | Relative Speed | Notes |
|------|---------------|-------|
| Uniform | 1× (baseline) | No ice thickness variation |
| Prescribed | 1× | Same as uniform (single column pass) |
| Equilibrium | 5-10× | Iterative (typically 5-10 iterations) |

**Optimization tips**:
- Use lower `nSide` for testing (e.g., nSide=4 → 192 pixels)
- Increase `equilibriumTol_m` for faster convergence (e.g., 200 m instead of 100 m)
- Reduce `equilibriumMaxIter` if you're okay with approximate solutions

## References

1. **Tobie et al. (2003)**: "Tidal dissipation within large icy satellites: Applications to Europa and Titan", *JGR* 108(E11), 5124, doi:10.1029/2003JE002099
   - Equilibrium ice thickness physics

2. **Tobie et al. (2005)**: "Solid tidal friction above a liquid water reservoir as the origin of the south polar hotspot on Enceladus", *Icarus* 196, 642-652
   - Geographic tidal heating patterns

3. **Ojakangas & Stevenson (1989)**: "Thermal state of an ice shell on Europa", *Icarus* 81, 220-241
   - Tidal strain pattern formalism

4. **Roberts & Nimmo (2008)**: "Tidal heating and the long-term stability of a subsurface ocean on Enceladus", *Icarus* 194, 675-689
   - Ice thickness-tidal heating coupling

5. **Beuthe (2013)**: "Crustal control of dissipative ocean tides in Enceladus and other icy moons", *Icarus* 223, 308-329
   - Thin-shell tidal dissipation theory

6. **Bolton et al. (2025)**: "Juno microwave radiometer observations of Europa", *Science*
   - Europa ice thickness constraints

## Frequently Asked Questions

### Q: Which mode should I use for my study?
**A**: Use **equilibrium mode** for scientific publications where you need physically self-consistent ice thickness. Use **prescribed mode** for exploratory studies or when matching observational constraints.

### Q: Why does equilibrium mode take so long?
**A**: Equilibrium mode iterates to self-consistency: thickness → thermal profile → viscosity → tidal heating → thickness. Each iteration requires running all lateral columns. Typical convergence takes 5-10 iterations.

### Q: Can I start with a prescribed pattern and then refine to equilibrium?
**A**: The equilibrium solver always starts from uniform thickness (from the reference 1D model). If you want to bias the initial guess, adjust `Planet.Bulk.Tb_K` to get the desired mean thickness in the reference model.

### Q: What if equilibrium mode doesn't converge?
**A**: Try:
1. Increase `equilibriumTol_m` (e.g., 200 m instead of 100 m)
2. Increase `equilibriumMaxIter` (e.g., 20 instead of 10)
3. Check that `kThermIce_WmK` and `qBasal_Wm2` are physically reasonable
4. Verify orbital parameters are correct

### Q: How do I know if my equilibrium solution is physical?
**A**: Check:
- Ice thickness is positive everywhere: `np.all(Planet.Lateral.dIce_m > 0)`
- Surface heat flux is realistic: `30-50 mW/m²` for Europa
- Variation matches expected patterns: thicker at sub/anti-Jovian, thinner at poles
- Converged residual is small: `< 100 m`

### Q: Can I use equilibrium mode for bodies other than Europa?
**A**: Yes! Equilibrium mode is body-agnostic. Just provide:
- Orbital parameters (`eccentricity`, `meanMotion_radps`)
- Thermal parameters (`Tsurf_K`, `Tb_K`, `kThermIce_WmK`)
- Silicate heating (Qrad_Wkg, Htidal_Wm3) or override `qBasal_Wm2`

Examples: Enceladus, Titan, Ganymede, Callisto, Triton, Pluto's moons

---

**Last updated**: 2026-06-18  
**Author**: Emma Vellard  
**Contact**: emma.vellard@outlook.fr
