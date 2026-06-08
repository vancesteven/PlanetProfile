# 3D Lateral Structure Examples

This directory contains working examples demonstrating PlanetProfile's 3D lateral structure capability.

## Quick Start

### Example 1: Simple Europa-like 3D Model (Fast)

```bash
# Minimal test (12 pixels, ~30 seconds)
python run_3d_example.py

# Medium resolution (48 pixels, ~2 minutes)
python run_3d_example.py --grid-size 2

# Higher resolution (192 pixels, ~10 minutes)
python run_3d_example.py --grid-size 4

# Fast mode without tidal heating (~10 seconds)
python run_3d_example.py --grid-size 1 --no-tidal
```

**Output:**
- Geographic maps in `Europa_3D_Example/figures/`
- 3D data in `Europa_3D_Example/Europa_3D_Example_lateral3D.pkl`
- Summary statistics printed to terminal

### Example 2: Titan Model (Tobie et al. 2005 Figure 10 Target)

```bash
# Default (768 pixels, ~20 minutes)
python run_3d_titan_tobie2005.py

# Lower resolution for testing (192 pixels, ~5 minutes)
python run_3d_titan_tobie2005.py --grid-size 4

# High resolution like Tobie 2005 (3072 pixels, ~2 hours)
python run_3d_titan_tobie2005.py --grid-size 16
```

**Output:**
- Geographic tidal heating maps showing ocean decoupling
- By-layer heating maps (cf. Tobie 2005 Figure 10 format)
- Comparison with Tobie 2005 published values

## What Each Example Demonstrates

### `run_3d_example.py` — Basic 3D Capability

**Features:**
- Europa-like configuration (Seawater ocean)
- Small grid for fast testing (12-192 pixels)
- 3D tidal heating H(r,θ,φ)
- Ice thickness variation via spherical harmonics (Y₂₀ pole-equator pattern)
- Mass conservation enforcement
- Geographic visualization (8 plot types)

**Physics:**
- Synchronous rotation + eccentricity → geographic tidal pattern
- Thicker ice at poles (5 km variation on 25 km mean)
- Ojakangas & Stevenson (1989) tidal heating formulation
- Maxwell rheology

**Use Case:**
- Verify 3D installation works
- Understand 3D workflow
- Test parameter variations quickly

### `run_3d_titan_tobie2005.py` — Scientific Validation

**Features:**
- Titan Model 2 configuration (Ice I above ocean)
- Medium to high resolution (768-3072 pixels)
- Ocean decoupling physics demonstration
- Comparison with Tobie et al. (2005) published results

**Physics:**
- Concentrated MgSO₄ ocean (100 ppt)
- Ice I above ocean → high tidal dissipation
- Ocean mechanically decouples surface ice from interior
- Target: ~3822 × 10⁻⁹ W/m³ at base of ice I (67× higher than Model 1)

**Use Case:**
- Validate 3D implementation against literature
- Reproduce Tobie et al. (2005) Figure 10
- Demonstrate ocean decoupling effect
- Publication-quality figures

## Command-Line Options

### `run_3d_example.py`

```
python run_3d_example.py [OPTIONS]

Options:
  --grid-size NSIDE    HEALPix nSide (1, 2, 4, 8, 16)
                       nPix = 12*nSide² 
                       Default: 2 (48 pixels)
  
  --no-tidal           Disable 3D tidal heating (faster)
  
  --no-plot            Skip plotting (data still saved)
  
  --help               Show help message
```

### `run_3d_titan_tobie2005.py`

```
python run_3d_titan_tobie2005.py [OPTIONS]

Options:
  --grid-size NSIDE    HEALPix nSide (4, 8, 16, 20)
                       Default: 8 (768 pixels)
  
  --model MODEL        Tobie 2005 model (1, 2, 3)
                       Default: 2 (Ice I above ocean)
                       [Models 1 & 3 not yet implemented]
  
  --help               Show help message
```

## Grid Size Guide

| nSide | nPixels | Time (no tidal) | Time (with tidal) | Use Case |
|-------|---------|-----------------|-------------------|----------|
| 1     | 12      | ~10 sec         | ~30 sec           | Quick test |
| 2     | 48      | ~30 sec         | ~2 min            | Testing |
| 4     | 192     | ~2 min          | ~10 min           | Development |
| 8     | 768     | ~8 min          | ~30 min           | Science |
| 16    | 3072    | ~30 min         | ~2 hours          | Publication |
| 20    | 4800    | ~1 hour         | ~4 hours          | High-res |

*Times approximate for Europa-like body on 8-core CPU with parallel processing.*

## Output Files

Both examples generate:

### Data Files
```
BodyName/
├── BodyName_Profile.txt         # 1D profile (full precision)
├── BodyName_Profile.csv         # 1D profile (spreadsheet format)
├── BodyName.pkl                 # 1D full state
└── BodyName_lateral3D.pkl       # 3D lateral data (NEW)
```

### Figures
```
BodyName/figures/
├── BodyName_IceThickness.pdf
├── BodyName_TidalHeating.pdf
├── BodyName_BasalTemperature.pdf
├── BodyName_OceanConductivity.pdf
├── BodyName_EffectiveConductivity.pdf
├── BodyName_LateralSummary.pdf          # Multi-panel overview
├── BodyName_TidalHeatingByLayer.pdf     # Tobie 2005 Fig 10 format
└── BodyName_TidalHeatingVsDepth.pdf     # H(r) profiles
```

## Interpreting Results

### Terminal Output

Both scripts print:
- Grid configuration (pixels, type)
- Ice thickness statistics (mean, min, max, range)
- Basal temperature statistics
- Surface heat flux statistics
- Tidal heating statistics (if enabled)
- Mass conservation residual (<0.01% is good)
- Execution time

### Geographic Maps

**Ice Thickness:** Shows spatial variation from SH coefficients or uniform
- Viridis colormap (blue = thin, yellow = thick)
- Contours at 8 levels
- Units: km

**Tidal Heating:** Shows surface heat flux from integrated H(r,θ,φ)
- Magma colormap (purple = low, yellow = high)
- Peaks at sub/anti-planetary points and poles
- Units: W/m²

**Basal Temperature:** Shows deviation from mean Tb_K
- Seismic diverging colormap (blue = cold, red = warm)
- Correlation with ice thickness (thicker → higher Pb → higher Tb)
- Units: K

**Ocean Conductivity:** Shows mean σ for each column
- Cividis colormap (sequential)
- Anticorrelation with ice thickness (thicker ice → colder ocean → lower σ)
- Units: S/m

**By-Layer Heating:** Following Tobie et al. (2005) Figure 10
- 2×2 panels: (top ice I, bottom ice I), (top HP ice, bottom HP ice)
- RdBu_r diverging colormap centered at zero
- Shows % deviation from mean heating
- Direct comparison with Tobie 2005 Figure 10

### Validation Checks

**Mass Conservation:**
- Residual should be <0.01% (0.0001)
- Larger residuals may indicate:
  - Mass conservation disabled
  - Convergence issues
  - Extreme spatial variations

**Physical Ranges:**
- Ice thickness: 1-100 km (Europa), 50-200 km (Titan)
- Basal temperature: 250-275 K (ice I melting range)
- Surface heat flux: 1-100 mW/m² (typical)
- Tidal heating: 10⁻⁹ to 10⁻⁵ W/m³ (depends on viscosity)

**3D vs 1D Consistency:**
- 3D spherical mean should approximate 1D reference
- Small differences (<5%) expected due to:
  - Mass conservation adjustment
  - Numerical precision
  - Spatial pattern effects

## Customization

### Modify Ice Thickness Pattern

Edit the spherical harmonic coefficients:

```python
# Example: Y₂₀ pole-equator + Y₂₂ sectoral pattern
Planet.Lateral.dIce_pMax = 2
Planet.Lateral.dIce_Cpq_km = np.array([
    [25.0, 0.0, 0.0],  # C₀₀ = mean 25 km
    [0.0, 0.0, 0.0],   # p=1
    [5.0, 0.0, 3.0]    # C₂₀ = 5 km (pole-eq), C₂₂ = 3 km (sectoral)
])
Planet.Lateral.dIce_Spq_km = np.zeros((3, 3))
```

### Add Clathrate Variation

```python
Planet.Lateral.DO_CLATH_LATERAL = True
Planet.Lateral.fClath_pMax = 2
Planet.Lateral.fClath_Cpq = np.array([
    [0.15, 0.0, 0.0],  # 15% mean clathrate
    [0.0, 0.0, 0.0],
    [0.05, 0.0, 0.0]   # 5% more at poles
])
Planet.Lateral.fClath_Spq = np.zeros((3, 3))
```

### Use TidalPy Backend

For self-consistent thermal-tidal coupling:

```python
Params.Gravity.backend = 'tidalpy'
Planet.Do.DO_SELF_CONSISTENT_HTIDAL = True
```

Requires: `pip install PlanetProfile[tidal]`

### Change Grid Type

Use lat-lon instead of HEALPix:

```python
Planet.Lateral.gridType = 'latlon'
Planet.Lateral.nLat = 37
Planet.Lateral.nLon = 72  # 2664 pixels
```

No healpy dependency required.

## Troubleshooting

### ImportError: No module named 'healpy'

**Solution:**
```bash
conda install -c conda-forge healpy
# or
pip install healpy
# or change gridType to 'latlon'
```

### PyALMA negative Love numbers

**Symptom:** k₂ < 0 (unphysical)

**Solution:**
```python
# Use TidalPy backend
Params.Gravity.backend = 'tidalpy'
# or simplify EOS
Planet.Sil.mantleEOS = 'CV3hy1wt_678_1.tab'
Planet.Core.xFeS = 0.55
```

### Memory issues

**Solution:**
- Reduce grid size: `--grid-size 2` instead of `--grid-size 8`
- Disable parallel: `Params.DO_PARALLEL = False` (slower but less memory)
- Close other applications

### 3D heating too low/high

**Possible causes:**
- Viscosity too high/low → check Planet.Sil.etaConv_Pas
- Eccentricity incorrect → check Planet.Bulk.eccentricity
- Rheology mismatch → Maxwell vs Andrade
- Ocean composition → affects Tb_K → affects viscosity

**Debug:**
1. Check 1D reference heating first
2. Verify orbital parameters (eccentricity, mean motion)
3. Compare with Tobie 2005 expected values
4. Check viscosity temperature dependence

## Next Steps

After running examples:

1. **Reload and replot** without recomputing:
   ```python
   from PlanetProfile.Lateral.LateralIO import ReloadLateralResults
   Planet = ReloadLateralResults(Planet, 'BodyName/BodyName_lateral3D.pkl')
   ```

2. **Custom analysis** with raw data:
   ```python
   theta = Planet.Lateral.theta_rad
   phi = Planet.Lateral.phi_rad
   d_ice = Planet.Lateral.dIce_m
   H_ice = Planet.Lateral.HtidalIce_Wm3
   # Your analysis here
   ```

3. **Feed to MoonMag** for asymmetric induction:
   ```python
   from MoonMag import MoonMag
   asymMag = MoonMag(Planet, Params)  # Uses Planet.Magnetic.asymShape_m automatically
   ```

4. **Create your own configuration** following `PPTest1_3D.py` template

## References

- **Ojakangas & Stevenson (1989):** Tidal heating formulation, *Icarus* 81:220-241
- **Tobie et al. (2005):** Geographic tidal dissipation in Titan, *Icarus* 177:534-549
- **Górski et al. (2005):** HEALPix equal-area pixelization, *ApJ* 622:759-771
- **Styczinski et al. (2021):** MoonMag asymmetric induction, *Icarus* 376:114840

## Support

- **Documentation:** `docs/LATERAL_3D_USAGE.md` (comprehensive guide)
- **Variable reference:** `VARIABLE_REFERENCES.md` (all Lateral fields)
- **GitHub Issues:** https://github.com/vancesteven/PlanetProfile/issues
- **Email:** steven.d.vance@jpl.nasa.gov
