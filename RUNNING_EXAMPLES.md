# Running 3D Examples - Setup Guide

The 3D lateral structure examples require a properly configured Python environment with all PlanetProfile dependencies.

## Environment Setup

### Option 1: Using Conda (Recommended)

```bash
# Create new environment
conda create -n planetprofile python=3.11
conda activate planetprofile

# Install PlanetProfile from source
cd /path/to/PlanetProfile
pip install -e .

# Install optional dependencies for 3D lateral structure
pip install -e ".[lateral,tidal]"
```

### Option 2: Using pip in existing environment

```bash
# In the PlanetProfile directory
pip install -e .  # Base installation

# Optional: Add 3D lateral structure dependencies
pip install healpy>=1.16.0  # For HEALPix grids
pip install TidalPy>=0.7.4  # For self-consistent tidal heating
```

## Verify Installation

```bash
python -c "from PlanetProfile.Main import PlanetProfile; print('✓ Import successful')"
python -c "import healpy; print('✓ HEALPix available')"
python -c "import scipy; print('✓ SciPy available')"
```

## Running the Examples

Once the environment is set up:

### Quick Test (30 seconds)

```bash
python run_3d_example.py --grid-size 1
```

**Expected output:**
- Terminal summary with statistics
- 8 PDF files in `Europa_3D_Example/figures/`

### Medium Resolution (10 minutes)

```bash
python run_3d_example.py --grid-size 4
```

### Titan Validation (30 minutes)

```bash
python run_3d_titan_tobie2005.py --grid-size 8
```

## Troubleshooting

### ModuleNotFoundError: No module named 'scipy'

**Cause:** PlanetProfile dependencies not installed

**Solution:**
```bash
pip install -e .  # In PlanetProfile directory
```

### ModuleNotFoundError: No module named 'healpy'

**Cause:** Optional lateral structure dependency not installed

**Solutions:**
```bash
# Option 1: Install healpy
conda install -c conda-forge healpy

# Option 2: Use lat-lon grid instead
# Edit example script, change:
Planet.Lateral.gridType = 'latlon'
Planet.Lateral.nLat = 37
Planet.Lateral.nLon = 72
```

### ImportError: No module named 'TidalPy'

**Cause:** Optional tidal heating dependency not installed (only needed for self-consistent heating)

**Solution:**
```bash
pip install TidalPy>=0.7.4
# or
pip install -e ".[tidal]"
```

### PyALMA negative Love numbers (k₂ < 0)

**Cause:** Known PyALMA limitation with certain parameter combinations

**Solution:**
```python
# In example script, add:
Params.Gravity.backend = 'tidalpy'  # Use TidalPy instead
# or simplify parameters:
Planet.Core.xFeS = 0.55  # Lower FeS fraction
```

## Expected Outputs

### Terminal Output

```
============================================================
3D Lateral Structure Example
============================================================
Grid: HEALPix nSide=1 (12 pixels)
Tidal heating: Enabled
Plotting: Enabled
============================================================

Running PlanetProfile with 3D lateral structure...
✓ Calculation complete in 30.5 seconds

============================================================
3D Results Summary
============================================================

Grid:
  Pixels: 12
  Grid type: healpix

Ice Thickness:
  Mean: 25.12 km
  Min: 22.34 km
  Max: 27.89 km
  Range: 5.55 km

Basal Temperature:
  Mean: 268.45 K
  Min: 267.92 K
  Max: 268.98 K
  Range: 1.06 K

Surface Heat Flux:
  Mean: 35.2 mW/m²
  Min: 32.1 mW/m²
  Max: 38.4 mW/m²

Tidal Heating (Ice):
  Mean: 2.45e-09 W/m³
  Min: 1.23e-09 W/m³
  Max: 4.67e-09 W/m³

Mass Conservation:
  Target mass: 4.7991e+22 kg
  Actual mass: 4.7991e+22 kg
  Residual: 0.0023%

Results saved to: Europa_3D_Example/Europa_3D_Example.pkl
Plots saved to: Europa_3D_Example/figures/

✓ 3D lateral structure example completed successfully!
```

### Figure Files Generated

```
Europa_3D_Example/figures/
├── Europa_3D_Example_IceThickness.pdf
├── Europa_3D_Example_TidalHeating.pdf
├── Europa_3D_Example_BasalTemperature.pdf
├── Europa_3D_Example_OceanConductivity.pdf
├── Europa_3D_Example_EffectiveConductivity.pdf
├── Europa_3D_Example_LateralSummary.pdf
├── Europa_3D_Example_TidalHeatingByLayer.pdf
└── Europa_3D_Example_TidalHeatingVsDepth.pdf
```

## Next Steps

After successfully running examples:

1. **Examine the plots** in the `figures/` directory
2. **Modify parameters** in the example scripts to test different scenarios
3. **Create your own configuration** following `PPTest1_3D.py` template
4. **Validate against Tobie et al. (2005)** using the Titan example

## Support

For setup issues:
- Check [Installation Guide](README.md#installation)
- Review [Dependencies](pyproject.toml)
- GitHub Issues: https://github.com/vancesteven/PlanetProfile/issues

For 3D usage:
- [3D Usage Guide](docs/LATERAL_3D_USAGE.md)
- [Example Documentation](3D_EXAMPLES_README.md)
- [Variable Reference](VARIABLE_REFERENCES.md)
