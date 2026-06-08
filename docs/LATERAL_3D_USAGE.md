# 3D Lateral Structure Usage Guide

PlanetProfile v3.2+ supports spatially-varying (3D) interior structure calculations for reproducing geographic tidal heating distributions, asymmetric magnetic induction, and spatially-resolved habitability assessments.

## Quick Start

### Minimal Example

```python
# In your PPBody.py configuration file:

# Enable 3D lateral structure
Planet.Lateral.DO_3D = True
Planet.Lateral.gridType = 'healpix'  # Equal-area grid
Planet.Lateral.nSide = 4  # 192 pixels (low resolution for testing)

# Run PlanetProfile as usual
# python PlanetProfileCLI.py YourBody
```

This creates:
- 3D interior structure with 192 independent columns
- Geographic maps saved to `YourBody/figures/`
- 3D results saved to `YourBody/YourBody_lateral3D.pkl`

### Full Example (Europa with Geographic Tidal Heating)

```python
# PPEuropa3D.py

from PlanetProfile import *

Planet = PlanetStruct('Europa')

# Bulk parameters (as usual)
Planet.Bulk.R_m = 1560.8e3
Planet.Bulk.M_kg = 4.7998e22
Planet.Bulk.Tsurf_K = 110.0
Planet.Bulk.Tb_K = 269.8  # 1D reference ice-ocean boundary temperature
Planet.Bulk.Cmeasured = 0.346
Planet.Bulk.Cuncertainty = 0.005

# Ocean composition
Planet.Ocean.comp = 'MgSO4'
Planet.Ocean.wOcean_ppt = 10.0

# Enable 3D lateral structure
Planet.Lateral.DO_3D = True
Planet.Lateral.gridType = 'healpix'
Planet.Lateral.nSide = 8  # 768 pixels (medium resolution)

# Enable 3D tidal heating
Planet.Lateral.DO_TIDAL_3D = True
Planet.Bulk.eccentricity = 0.0094  # Europa's eccentricity
Planet.Bulk.meanMotion_radps = 2.05e-5  # 3.55 day orbital period

# Optionally: add ice thickness variation via spherical harmonics
# Example: Y_20 pattern (pole-equator variation)
Planet.Lateral.dIce_pMax = 2
Planet.Lateral.dIce_Cpq_km = np.array([
    [30.0, 0.0, 0.0],  # C_00 = mean thickness 30 km
    [0.0, 0.0, 0.0],   # p=1 (no degree-1)
    [5.0, 0.0, 0.0]    # C_20 = 5 km pole-equator variation
])
Planet.Lateral.dIce_Spq_km = np.zeros((3, 3))  # All sine coefficients zero

# Run PlanetProfile
PlanetProfile(Planet)
```

### Output Files

When `DO_3D=True`, PlanetProfile generates additional outputs:

**Data files:**
- `YourBody/YourBody_lateral3D.pkl` — All 3D results (ice thickness, heating, temperatures, etc.)
- Existing 1D files unchanged (`.txt`, `.csv`, `.pkl`)

**Figures (if `Params.PLOT_LATERAL=True`, default):**
- `YourBody_IceThickness.pdf` — Geographic ice shell thickness map
- `YourBody_TidalHeating.pdf` — Surface heat flux from tidal dissipation
- `YourBody_BasalTemperature.pdf` — Ice-ocean boundary temperature deviation
- `YourBody_OceanConductivity.pdf` — Mean ocean electrical conductivity
- `YourBody_ClathrateFraction.pdf` — Clathrate volume fraction (if `DO_CLATH_LATERAL=True`)
- `YourBody_EffectiveConductivity.pdf` — Effective thermal conductivity
- `YourBody_LateralSummary.pdf` — Multi-panel figure with 5 key fields
- `YourBody_TidalHeatingByLayer.pdf` — By-layer heating (cf. Tobie et al. 2005 Fig 10)
- `YourBody_TidalHeatingVsDepth.pdf` — H(r) profiles at 4 representative locations

## Configuration Options

### Grid Selection

**HEALPix (recommended for most cases):**
```python
Planet.Lateral.gridType = 'healpix'
Planet.Lateral.nSide = 8  # nPix = 12 * nSide^2
```
- Equal-area pixels (uniform spatial resolution)
- No pole singularities
- Efficient spherical harmonic transforms
- Requires `healpy` (install with `conda install -c conda-forge healpy`)
- Common values: nSide=1 (12 pix), 4 (192 pix), 8 (768 pix), 16 (3072 pix)

**Lat-lon (no dependencies, simpler):**
```python
Planet.Lateral.gridType = 'latlon'
Planet.Lateral.nLat = 37
Planet.Lateral.nLon = 72  # nPix = nLat * nLon = 2664
```
- Regular grid in latitude/longitude
- Non-equal-area (equator pixels ~24× larger than pole pixels)
- No additional dependencies
- Better for rectangular output formats

### Ice Thickness Specification

**Mode 1: Uniform (default)**
```python
# Uses 1D reference ice thickness from Planet.Bulk.Tb_K
Planet.Lateral.DO_3D = True
# No additional configuration needed
```

**Mode 2: Spherical Harmonics**
```python
Planet.Lateral.dIce_pMax = 4  # Use up to degree 4
Planet.Lateral.dIce_Cpq_km = np.array([...])  # (pMax+1, pMax+1) cosine coeffs
Planet.Lateral.dIce_Spq_km = np.array([...])  # (pMax+1, pMax+1) sine coeffs
```
- Cpq[0,0] = C_00 = mean ice thickness in km
- Higher degrees add spatial variation
- 4π normalization: Y_pq = P_pq(cos θ) × [C_pq cos(qφ) + S_pq sin(qφ)]

**Mode 3: Custom Function**
```python
def ice_thickness_func(theta_rad):
    """Ice thickness in meters as function of colatitude"""
    return 25e3 + 10e3 * np.cos(theta_rad)  # Thicker at poles

Planet.Lateral.dIce_func = ice_thickness_func
```

### Tidal Heating Options

**Disable tidal heating (fastest):**
```python
Planet.Lateral.DO_TIDAL_3D = False  # Default
```

**Enable 3D tidal heating:**
```python
Planet.Lateral.DO_TIDAL_3D = True
Planet.Bulk.eccentricity = 0.0094
Planet.Bulk.meanMotion_radps = 2 * np.pi / (3.55 * 86400)  # orbital period
```
- Computes H(r,θ,φ) from Ojakangas & Stevenson (1989) formulation
- Geographic pattern depends on synchronous rotation + eccentricity
- Requires orbital parameters (eccentricity, mean motion)

**Self-consistent thermal-tidal coupling (requires TidalPy):**
```python
Planet.Do.DO_SELF_CONSISTENT_HTIDAL = True
Params.Gravity.backend = 'tidalpy'
Planet.Lateral.DO_TIDAL_3D = True
```
- Iteratively converges heating with temperature-dependent viscosity
- Per-phase heating rates from viscoelastic theory
- Slower but more physically realistic

### Clathrate Variation

```python
Planet.Lateral.DO_CLATH_LATERAL = True
Planet.Lateral.fClath_pMax = 2
Planet.Lateral.fClath_Cpq = np.array([
    [0.15, 0.0, 0.0],  # C_00 = 15% mean clathrate fraction
    [0.0, 0.0, 0.0],
    [0.05, 0.0, 0.0]   # C_20 = 5% more at poles
])
Planet.Lateral.fClath_Spq = np.zeros((3, 3))
```
- Volume fraction of clathrate in mixed ice (0-1)
- Affects thermal conductivity: k_eff = f_clath * k_clath + (1 - f_clath) * k_ice
- k_clath ≈ 0.5 W/m/K (lower than ice Ih ~3 W/m/K at 200K)

### Mass Conservation

```python
Planet.Lateral.DO_MASS_CONSERVE = True  # Default
```
- Adjusts ocean floor radius to match total body mass
- Preserves ice thickness spatial pattern
- Typical residuals <0.01% after adjustment
- Disable only for diagnostic purposes

### Plotting Control

```python
Params.PLOT_LATERAL = True  # Default, create all geographic maps
Params.PLOT_LATERAL = False  # Skip plotting (faster, data still saved)
```

## Computational Requirements

### Time

| Grid Resolution | nPixels | Time (no tidal) | Time (with 3D tidal) |
|-----------------|---------|-----------------|----------------------|
| nSide=1 (test)  | 12      | ~10 sec         | ~30 sec              |
| nSide=4 (low)   | 192     | ~2 min          | ~5 min               |
| nSide=8 (medium)| 768     | ~8 min          | ~20 min              |
| nSide=16 (high) | 3072    | ~30 min         | ~90 min              |

*Times for Europa-like body on 8-core CPU with `Params.DO_PARALLEL=True`. Without parallelization, multiply by ~6×.*

### Memory

| Grid Resolution | RAM Usage |
|-----------------|-----------|
| nSide=1         | ~100 MB   |
| nSide=4         | ~300 MB   |
| nSide=8         | ~800 MB   |
| nSide=16        | ~2 GB     |

Memory scales linearly with nPix because each column stores full profile + EOS data.

### Parallelization

```python
Params.DO_PARALLEL = True  # Default, use all CPU cores
```
- Columns are embarrassingly parallel (no communication)
- Speedup ≈ min(nCores, nPixels)
- Disable for debugging: `Params.DO_PARALLEL = False`

## Use Cases

### 1. Tobie et al. (2005) Figure 10 Reproduction

**Goal:** Reproduce geographic tidal heating distribution in Titan's ice shell showing ocean decoupling effect.

```python
Planet = PlanetStruct('Titan')
Planet.Bulk.Tb_K = 255.0
Planet.Ocean.comp = 'MgSO4'
Planet.Ocean.wOcean_ppt = 100.0  # Concentrated ocean

Planet.Lateral.DO_3D = True
Planet.Lateral.gridType = 'healpix'
Planet.Lateral.nSide = 8  # 768 pixels

Planet.Lateral.DO_TIDAL_3D = True
Planet.Bulk.eccentricity = 0.0288  # Titan's eccentricity
Planet.Bulk.meanMotion_radps = 2 * np.pi / (15.945 * 86400)  # 15.945 day period

# Use Andrade rheology (more realistic than Maxwell)
Params.Gravity.backend = 'tidalpy'
Planet.Gravity.QS_Andrade = 100
Planet.Gravity.andrade_zeta = 1.0
```

**Expected output:**
- Geographic maps showing 67× heating ratio between ice above vs below ocean
- Peak heating at sub/anti-Saturn points and poles
- Validates ocean decoupling physics

### 2. Asymmetric Magnetic Induction (MoonMag Integration)

**Goal:** Compute induced magnetic field from non-uniform ice-ocean boundary.

```python
Planet = PlanetStruct('Europa')
# ... standard Europa configuration ...

Planet.Lateral.DO_3D = True
Planet.Lateral.nSide = 4  # 192 pixels sufficient for induction

# Add ice thickness variation
Planet.Lateral.dIce_pMax = 4
Planet.Lateral.dIce_Cpq_km = load_ice_thickness_sh_from_data()  # Your SH coefficients

# Run PlanetProfile (automatically populates Planet.Magnetic.asymShape_m)
PlanetProfile(Planet)

# Feed to MoonMag
from MoonMag import MoonMag
asymMag = MoonMag(Planet, Params)  # Uses asymShape_m automatically
```

**Expected output:**
- Asymmetric induced dipole Bx, By, Bz components
- Spatial variations ~10-50 nT for Europa
- Can validate against Galileo magnetometry

### 3. Spatial Habitability Mapping

**Goal:** Identify regions with highest astrobiological potential.

```python
Planet = PlanetStruct('Enceladus')
# ... standard Enceladus configuration ...

Planet.Lateral.DO_3D = True
Planet.Lateral.DO_TIDAL_3D = True  # Tidal heating → warm spots
Planet.Lateral.nSide = 8

# After running, analyze results:
lateral_data = Planet.Lateral
high_heat_flux_regions = lateral_data.qBase_Wm2 > np.percentile(lateral_data.qBase_Wm2, 90)
thin_ice_regions = lateral_data.dIce_m < np.percentile(lateral_data.dIce_m, 25)

# High-priority landing sites: thin ice + high basal heat flux
priority_sites = high_heat_flux_regions & thin_ice_regions
```

**Expected output:**
- Maps showing ocean thickness, basal heat flux, ice permeability
- Target regions for future missions

## Loading and Replotting Results

```python
from PlanetProfile.Lateral.LateralIO import ReloadLateralResults
from PlanetProfile.Plotting.LateralPlots import GenerateLateralPlots

# Reload from file (much faster than recomputing)
Planet = PlanetStruct('Europa')
Planet = ReloadLateralResults(Planet, 'Europa/Europa_lateral3D.pkl')

# Regenerate plots with new settings
Params = ConfigParams()
Params.FigMisc.figFormat = 'png'  # Change to PNG instead of PDF
GenerateLateralPlots(Planet, Params)
```

## Troubleshooting

### ImportError: No module named 'healpy'

**Solution:**
```bash
# Option 1: Install healpy
conda install -c conda-forge healpy

# Option 2: Use lat-lon grid instead
Planet.Lateral.gridType = 'latlon'
Planet.Lateral.nLat = 37
Planet.Lateral.nLon = 72
```

### PyALMA negative Love numbers (k₂ < 0)

**Symptom:** Unphysical negative k₂ for complex mantle EOS + high core FeS content.

**Solution:**
```python
# Option 1: Use TidalPy backend instead
Params.Gravity.backend = 'tidalpy'

# Option 2: Simplify parameters
Params.Sil.mantleEOS = 'CV3hy1wt_678_1.tab'  # Simpler EOS
Params.Core.xFeS = 0.55  # Lower FeS fraction
```

### Memory issues with large grids

**Solution:**
```python
# Reduce grid resolution
Planet.Lateral.nSide = 4  # Instead of 16

# Or: Disable parallel processing (uses less memory)
Params.DO_PARALLEL = False

# Or: Clear EOS between columns (not yet implemented, TODO)
```

### 3D mean doesn't match 1D reference

**Expected:** Small differences (<5%) due to:
- Mass conservation adjustment
- Numerical precision
- Ice thickness pattern

**Concerning (>10% difference):** Check:
- Ice thickness SH coefficients are reasonable
- No extreme spatial variations
- Mass conservation converged (`massResidual_frac < 1e-4`)

## Advanced: Custom Analysis

Access raw 3D data for custom analysis:

```python
# After running with DO_3D=True
nPix = Planet.Lateral.nPix
theta = Planet.Lateral.theta_rad  # Colatitude [rad]
phi = Planet.Lateral.phi_rad      # Longitude [rad]

# Ice thickness field
d_ice = Planet.Lateral.dIce_m  # (nPix,) [m]

# Basal temperature field
Tb = Planet.Lateral.Tb_K  # (nPix,) [K]

# Surface heat flux
q_surf = Planet.Lateral.qSurf_Wm2  # (nPix,) [W/m²]

# Mean ocean conductivity
sigma = Planet.Lateral.sigma_mean_Sm  # (nPix,) [S/m]

# Tidal heating (if DO_TIDAL_3D=True)
H_ice = Planet.Lateral.HtidalIce_Wm3  # (nPix,) [W/m³]

# Example: Compute spherical mean
from PlanetProfile.Lateral.SpatialGrid import IntegrateOverSphere
d_ice_mean_m = IntegrateOverSphere(d_ice, Planet.Lateral.pixArea_sr)

# Example: Convert to lat-lon for plotting
import numpy as np
lat_deg = 90 - np.degrees(theta)  # [-90, 90]
lon_deg = np.degrees(phi)  # [0, 360]

import matplotlib.pyplot as plt
plt.scatter(lon_deg, lat_deg, c=d_ice/1e3, cmap='viridis')
plt.xlabel('Longitude [°E]')
plt.ylabel('Latitude [°N]')
plt.colorbar(label='Ice thickness [km]')
plt.show()
```

## References

- **Ojakangas & Stevenson (1989):** Tidal heating with eccentricity forcing, *Icarus* 81:220-241, DOI: 10.1016/0019-1035(89)90052-3
- **Tobie et al. (2005):** Geographic tidal dissipation in Titan, *Icarus* 177:534-549, DOI: 10.1016/j.icarus.2005.04.006
- **Górski et al. (2005):** HEALPix equal-area pixelization, *ApJ* 622:759-771, DOI: 10.1086/427976
- **Styczinski et al. (2021):** MoonMag asymmetric induction method, *Icarus* 376:114840, DOI: 10.1016/j.icarus.2021.114840
- **Hobbs (1974):** Ice Ih thermal conductivity, *Ice Physics*, Oxford University Press

## Support

For questions or issues:
- GitHub: https://github.com/vancesteven/PlanetProfile/issues
- Email: steven.d.vance@jpl.nasa.gov
- Documentation: https://vancesteven.github.io/PlanetProfile
