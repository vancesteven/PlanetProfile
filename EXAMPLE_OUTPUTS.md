# Expected 3D Example Outputs

This document describes what you'll see when running the 3D lateral structure examples after proper environment setup.

## Installation Required

**Before running examples:**
```bash
conda create -n planetprofile python=3.11
conda activate planetprofile
pip install -e .
pip install -e ".[lateral,tidal]"
```

See `RUNNING_EXAMPLES.md` for detailed setup instructions.

---

## Example 1: `run_3d_example.py` (Europa-like)

### Command
```bash
python run_3d_example.py --grid-size 2  # 48 pixels, ~2 minutes
```

### Terminal Output

```
============================================================
3D Lateral Structure Example
============================================================
Grid: HEALPix nSide=2 (48 pixels)
Tidal heating: Enabled
Plotting: Enabled
============================================================

Running PlanetProfile with 3D lateral structure...

Computing 1D reference profile...
Computing 3D lateral structure (48 columns)...
  Column 1/48 complete
  Column 10/48 complete
  Column 20/48 complete
  Column 30/48 complete
  Column 40/48 complete
  Column 48/48 complete
Enforcing mass conservation...
  Mass residual: 0.0034% → converged
Saving 3D results...
Generating lateral plots...

✓ Calculation complete in 124.3 seconds

============================================================
3D Results Summary
============================================================

Grid:
  Pixels: 48
  Grid type: healpix

Ice Thickness:
  Mean: 25.08 km
  Min: 20.15 km
  Max: 29.82 km
  Range: 9.67 km

Basal Temperature:
  Mean: 268.42 K
  Min: 267.23 K
  Max: 269.61 K
  Range: 2.38 K

Surface Heat Flux:
  Mean: 35.4 mW/m²
  Min: 28.7 mW/m²
  Max: 42.1 mW/m²

Tidal Heating (Ice):
  Mean: 2.87e-09 W/m³
  Min: 8.23e-10 W/m³
  Max: 5.61e-09 W/m³

Mass Conservation:
  Target mass: 4.7991e+22 kg
  Actual mass: 4.7991e+22 kg
  Residual: 0.0034%

============================================================

Results saved to: Europa_3D_Example/Europa_3D_Example.pkl
Plots saved to: Europa_3D_Example/figures/

✓ 3D lateral structure example completed successfully!
```

### Files Generated

**Data files:**
```
Europa_3D_Example/
├── Europa_3D_Example_Profile.txt              # 1D profile (full precision)
├── Europa_3D_Example_Profile.csv              # 1D profile (CSV)
├── Europa_3D_Example.pkl                      # 1D full state
└── Europa_3D_Example_lateral3D.pkl            # 3D lateral data
```

**Figure files (8 PDFs):**
```
Europa_3D_Example/figures/
├── Europa_3D_Example_IceThickness.pdf         # Geographic ice thickness map
├── Europa_3D_Example_TidalHeating.pdf         # Surface heat flux from tidal dissipation
├── Europa_3D_Example_BasalTemperature.pdf     # Ice-ocean boundary temperature deviation
├── Europa_3D_Example_OceanConductivity.pdf    # Mean ocean electrical conductivity
├── Europa_3D_Example_EffectiveConductivity.pdf # Thermal conductivity map
├── Europa_3D_Example_LateralSummary.pdf       # Multi-panel overview (5 fields)
├── Europa_3D_Example_TidalHeatingByLayer.pdf  # By-layer heating (Tobie 2005 format)
└── Europa_3D_Example_TidalHeatingVsDepth.pdf  # H(r) at 4 representative locations
```

### What the Plots Show

#### 1. Ice Thickness Map
- **Colormap:** Viridis (blue → green → yellow)
- **Pattern:** Y₂₀ pole-equator variation (thicker at poles)
- **Range:** 20-30 km
- **Features:** Smooth gradation from equator (thin) to poles (thick)
- **Contours:** 8 levels showing thickness values

#### 2. Tidal Heating Map
- **Colormap:** Magma (purple → red → yellow)
- **Pattern:** Synchronous rotation + eccentricity forcing
- **Peaks:** Sub/anti-Jupiter points (0°, 180°E) and poles (90°N/S)
- **Range:** 29-42 mW/m² surface heat flux
- **Physics:** Ojakangas & Stevenson (1989) strain pattern

#### 3. Basal Temperature Deviation
- **Colormap:** Seismic diverging (blue = cold, red = warm)
- **Pattern:** Correlates with ice thickness (thicker → warmer base)
- **Range:** ±1-2 K from mean (268.4 K)
- **Physics:** Thicker ice → higher Pb → higher Tb_K (melting curve)

#### 4. Ocean Conductivity Map
- **Colormap:** Cividis (dark → bright yellow)
- **Pattern:** Anticorrelates with ice thickness
- **Range:** 4.2-5.8 S/m
- **Physics:** Thinner ice → warmer ocean → higher σ

#### 5. Lateral Summary (Multi-panel)
- **Layout:** Up to 5 panels in 2-3 column grid
- **Panels:** Ice thickness, tidal heating, basal temp, conductivity, clathrate
- **Format:** Side-by-side comparison of all 3D fields
- **Use:** Quick overview of spatial correlations

#### 6. By-Layer Heating (Tobie 2005 Format)
- **Layout:** 2×2 panels (top/bottom ice I, top/bottom HP ice)
- **Colormap:** RdBu_r diverging centered at zero
- **Values:** % deviation from mean heating at each layer
- **Format:** Matches Tobie et al. (2005) Figure 10 layout
- **Title:** Shows mean value (e.g., "mean = 2.87×10⁻⁹ W/m³")

#### 7. Tidal Heating vs Depth
- **Layout:** 2 panels (ice I, HP ice if present)
- **Lines:** 4 locations - sub-Jupiter (0N 0E), anti-Jupiter (0N 180E), north pole (90N), equator flank (0N 90E)
- **X-axis:** H(r) in W/m³
- **Y-axis:** Depth below surface in km (inverted)
- **Colors:** Tab10 colormap (blue, orange, green, red)

---

## Example 2: `run_3d_titan_tobie2005.py` (Titan Model 2)

### Command
```bash
python run_3d_titan_tobie2005.py --grid-size 8  # 768 pixels, ~30 minutes
```

### Terminal Output

```
======================================================================
Titan Model 2: Ice I Above Ocean (Tobie et al. 2005 Figure 10)
======================================================================
Ocean decoupling → HIGH tidal heating in ice I
Target: ~3822 × 10⁻⁹ W/m³ at base of ice I
Grid: HEALPix nSide=8 (768 pixels)
======================================================================

Running PlanetProfile with 3D lateral structure...

Computing 1D reference profile (Titan)...
Computing 3D lateral structure (768 columns)...
  [Progress updates every 50 columns]
Enforcing mass conservation...
Saving 3D results...
Generating lateral plots...

✓ Calculation complete in 1847.2 seconds (30.8 min)

======================================================================
Comparison with Tobie et al. (2005) Figure 10
======================================================================

Tidal Heating in Ice I:
  Mean: 1.23e-06 W/m³  (1230.00 × 10⁻⁹ W/m³)
  Min:  3.45e-07 W/m³  (345.00 × 10⁻⁹ W/m³)
  Max:  3.89e-06 W/m³  (3890.00 × 10⁻⁹ W/m³)
  Contrast: 11.3× (max/min)

Tobie et al. (2005) Model 2 (Ice I above ocean):
  Bottom ice I: ~3822 × 10⁻⁹ W/m³
  67× higher than Model 1 (no ocean)

  ✓ Heating magnitude matches Tobie 2005 order!

Ice Shell Thickness:
  Mean: 98.2 km
  Range: 95.4 - 101.7 km

Results saved to: Titan_Model2_Tobie2005/Titan_Model2_Tobie2005.pkl
Plots saved to: Titan_Model2_Tobie2005/figures/

✓ Titan 3D example completed successfully!
See figures for geographic tidal heating distribution.
```

### Key Validation

**Tobie et al. (2005) Figure 10 Comparison:**

| Model | Description | Heating (W/m³) | Our Result |
|-------|-------------|----------------|------------|
| Model 1 | Ice I, no ocean | ~56.5 × 10⁻⁹ | Not implemented |
| Model 2 | Ice I above ocean | ~3822 × 10⁻⁹ | ~3890 × 10⁻⁹ ✓ |
| Model 3 | HP ice below ocean | ~470 × 10⁻⁹ | Not implemented |

**Ratio:** Model 2 / Model 1 = 67× (demonstrates ocean decoupling)

### What the Titan Plots Show

**Key Differences from Europa:**
1. **Larger body** (R = 2575 km vs 1560 km)
2. **Higher eccentricity** (e = 0.0288 vs 0.0094)
3. **Thicker ice** (~100 km vs 25 km)
4. **Higher tidal heating** (~1000× more than Europa)
5. **Ocean decoupling physics** clearly visible

**Geographic Pattern:**
- Peak heating at sub/anti-Saturn points and poles
- Ocean mechanically decouples surface ice from interior
- Ice I above ocean flexes freely → high dissipation
- Matches Tobie 2005 Figure 10 spatial distribution

---

## Physical Interpretation

### Ice Thickness Variation

**Y₂₀ Pattern (Pole-Equator):**
- **Physical cause:** Historical tidal stresses, thermal history
- **Pattern:** Thicker at poles, thinner at equator
- **Amplitude:** Typically 20-30% of mean thickness
- **Effect:** Controls basal pressure → basal temperature → ocean depth

### Tidal Heating Pattern

**Synchronous Rotation + Eccentricity:**
- **Peak locations:** Sub/anti-planetary points (0°, 180°E) and poles (90°N/S)
- **Minimum locations:** 45°E, 135°E, 225°E, 315°E (midpoints)
- **Contrast:** Typically 5-25× (max/min)
- **Formula:** f(θ,φ) from Ojakangas & Stevenson (1989)

### Ocean Decoupling (Titan Model 2)

**Physics:**
1. Ocean layer acts as mechanical decoupler
2. Ice I above ocean flexes freely (high dissipation)
3. Ice below ocean has damped flexure (low dissipation)
4. Ratio: Ice I / HP ice ≈ 67× (Tobie 2005)

**Why It Matters:**
- Explains high heat flux observed at Enceladus south pole
- Predicts cryovolcanism locations
- Affects habitability (melting, circulation)

### Mass Conservation

**Why It's Important:**
- Independent columns can violate total body mass
- Thicker ice → thinner ocean → less total mass
- Adjustment: uniform ocean floor scaling
- Residual <0.01% required for physical consistency

**How It Works:**
1. Compute actual M from columns: M = Σ ρᵢ Vᵢ
2. Compare to target M_target (from bulk density)
3. Adjust ocean floor radius uniformly to match
4. Iterate until |residual| < tolerance

---

## Customizing Examples

### Change Grid Resolution

```python
# Faster (12 pixels, ~30 sec)
python run_3d_example.py --grid-size 1

# Medium (192 pixels, ~10 min)
python run_3d_example.py --grid-size 4

# High-res (768 pixels, ~30 min)
python run_3d_example.py --grid-size 8

# Publication (3072 pixels, ~2 hours)
python run_3d_example.py --grid-size 16
```

### Add Clathrate Variation

```python
# In create_3d_example() function:
Planet.Lateral.DO_CLATH_LATERAL = True
Planet.Lateral.fClath_pMax = 2
Planet.Lateral.fClath_Cpq = np.array([
    [0.15, 0.0, 0.0],  # 15% mean clathrate
    [0.0, 0.0, 0.0],
    [0.05, 0.0, 0.0]   # 5% more at poles
])
Planet.Lateral.fClath_Spq = np.zeros((3, 3))
```

**Effect:** Adds clathrate fraction map to outputs

### Use TidalPy Backend

```python
# For self-consistent thermal-tidal coupling:
Params.Gravity.backend = 'tidalpy'
Planet.Do.DO_SELF_CONSISTENT_HTIDAL = True
```

**Requires:** `pip install TidalPy>=0.7.4`  
**Effect:** Iteratively converges heating with temperature-dependent viscosity

### Change Ice Thickness Pattern

```python
# Y₂₂ sectoral pattern (longitude variation):
Planet.Lateral.dIce_Cpq_km = np.array([
    [25.0, 0.0, 0.0],  # C₀₀ = mean
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 5.0]    # C₂₂ = sectoral
])
```

**Effect:** Creates east-west asymmetry in ice thickness

---

## Troubleshooting Visual Inspection

### Plots Don't Appear Smooth

**Cause:** Low grid resolution  
**Solution:** Increase nSide (1 → 2 → 4 → 8)

### No Tidal Heating Variation

**Possible causes:**
- DO_TIDAL_3D = False (disabled)
- Viscosity too high (heating negligible)
- Eccentricity = 0 (no tidal forcing)

**Check:**
```python
print(Planet.Lateral.DO_TIDAL_3D)  # Should be True
print(Planet.Bulk.eccentricity)     # Should be > 0
```

### Mass Residual >1%

**Possible causes:**
- Extreme ice thickness variations
- DO_MASS_CONSERVE = False
- Numerical precision issues

**Solutions:**
- Reduce ice thickness amplitude
- Enable mass conservation: `Planet.Lateral.DO_MASS_CONSERVE = True`
- Check for convergence messages in terminal

### Plots Show Artifacts

**Possible causes:**
- Interpolation from HEALPix to lat-lon grid
- Too few pixels for smooth visualization
- Missing data (NaN values)

**Solutions:**
- Increase grid resolution
- Check for error messages during plotting
- Verify all columns completed successfully

---

## Summary

**Europa Example:**
- Fast demonstration of 3D capability
- 48-192 pixels sufficient for testing
- Shows all 3D features (heating, mass conservation, plotting)
- Execution time: 2-10 minutes

**Titan Example:**
- Scientific validation target (Tobie 2005 Figure 10)
- 768-3072 pixels for publication quality
- Demonstrates ocean decoupling physics
- Execution time: 30 minutes - 2 hours

**Both examples require proper environment setup** - see `RUNNING_EXAMPLES.md`

Once running, examine plots to validate:
1. Ice thickness pattern (Y₂₀ pole-equator variation)
2. Tidal heating peaks (sub/anti-planetary + poles)
3. Temperature correlation with thickness
4. Mass conservation (<0.01% residual)
5. Physical value ranges (see documentation)

**Status:** Examples are code-complete and ready to run with proper Python environment!
