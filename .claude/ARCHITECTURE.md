# PlanetProfile Architecture Reference

This document provides detailed architecture and module structure for PlanetProfile. For the mission-critical workflow, see [../CLAUDE.md](../CLAUDE.md).

---

## Project Overview

PlanetProfile is a scientific software framework for constructing 1D interior structure models of planetary bodies (especially ocean worlds like Europa, Enceladus, Titan, etc.). It uses self-consistent thermodynamics for fluid, rock, and mineral phases, calculating properties like sound speeds, attenuation, and electrical conductivities.

**This is peer-reviewed scientific software used in planetary science research.**

---

## Repository Layout

```
PlanetProfile/              # Python package
  Main.py                   # Top-level entry: run(), PlanetProfile(), InductOgram(), ExploreOgram()
  __init__.py               # Package initialization, cache paths (_CACHE, _PERPLEXCACHE, _REFPROFILESCACHE)
  
  Utilities/
    defineStructs.py        # Core data structures: PlanetStruct, Constants, EOSlist, Timing
    SetupInit.py            # Initialization and setup
    ResultsIO.py            # Save/load results as .pkl, .mat, .txt
    PPversion.py            # Version info
    
  Thermodynamics/
    LayerPropagators.py     # Layer propagation: IceLayers(), OceanLayers(), InnerLayers()
    HydroEOS.py             # Ocean equation of state handling (GetOceanEOS)
    OceanProps.py           # Ocean chemistry, pH, speciation (LiquidOceanPropsCalcs)
    Seismic.py              # Seismic velocity calculations
    Electrical.py           # Electrical conductivity
    Viscosity.py            # Viscosity calculations
    Reaktoro/CustomSolution/  # Arbitrary ocean compositions via geochemical databases
    EOSdata/Perple_X/       # Perple_X EOS data files for silicates/core (~164 MB, downloaded during install)
    
  MagneticInduction/
    MagneticInduction.py    # Forward modeling with MoonMag
    Moments.py              # Magnetic moment calculations
    
  Gravity/
    Gravity.py              # Tidal Love numbers using PyALMA3
    
  Plotting/
    ProfilePlots.py         # Interior structure profiles
    MagPlots.py             # Magnetic induction plots
    ExplorationPlots.py     # Multi-model exploration visualizations
    
  Test/                     # Test modules: PPTest1.py through PPTest28.py
    TestBayes.py            # Bayesian inference tests
    
  Default/                  # Default body configurations
    <Body>/
      PP<Body>.py           # Default configuration for each body
      
  CustomSolution/           # Custom ocean solution support
  MonteCarlo/               # Monte Carlo sampling (work in progress)
  Inversion/                # Inversion routines
  TrajecAnalysis/           # Trajectory analysis tools
  SPICE/                    # SPICE kernel support
  install.py                # Installation helper
  BuildTest.py              # Full test suite runner (alternative to Testing.py)

<Body>/                     # Per-body directories (e.g., Europa/, Titan/)
  PP<Body>.py               # User configuration (if different from default)
  figures/                  # Generated plots
  inductionData/            # Magnetic induction results
  seismicData/              # Seismic model outputs
  *.pkl, *.mat, *.txt       # Results files

PlanetProfileCLI.py         # CLI entry point for developers
Testing.py                  # Full test suite runner (primary)
pyproject.toml              # Build config, dependencies, version
MANIFEST.in                 # Source distribution file inclusions
CHANGELOG.md                # Release history
VARIABLE_REFERENCES.md      # Canonical Planet/Params field names and units
REFERENCES.md               # Scientific papers underpinning PlanetProfile
```

---

## Data Flow

### 1. Entry Point
User runs:
- `python PlanetProfileCLI.py Europa` (developers)
- `python -m PlanetProfile.Main Europa` (pip users)

### 2. Configuration Loading
`Main.run()` loads:
- `Params` from `GetConfig.py` (global settings)
- `Planet` structure from body configuration file

### 3. Main Orchestration
`PlanetProfile()` function in `Main.py`:
1. **Initialize**: `SetupInit()` prepares structures
2. **Layer propagation**: 
   - `IceLayers()` — surface to ice-ocean boundary
   - `OceanLayers()` — ocean (if present)
   - `InnerLayers()` — silicate mantle and core
3. **Property calculations**:
   - Seismic properties via `SeismicCalcs()`
   - Electrical conductivity via `ElecConduct()`
   - Ocean chemistry via `LiquidOceanPropsCalcs()`
4. **Optional forward models**:
   - Gravity/Love numbers via `GravityParameters()`
   - Magnetic induction via `MagneticInduction()`
5. **Results output**:
   - Save: `WriteResults()` → `.pkl`, `.mat`, `.txt`
   - Plot: `GeneratePlots()` → `<body>/figures/`

### 4. Caching and Reloading
- `CALC_NEW` = False: Reload from cached `.pkl` files
- `CALC_NEW` = True: Recalculate everything
- Specific flags: `CALC_NEW_REF`, `CALC_NEW_INDUCT`, `CALC_NEW_GRAVITY`

---

## Body Configuration Files

Configuration files define planetary parameters. Structure:

### `Planet.Bulk.*` — Bulk Properties
- `R_m` — Mean radius (meters)
- `M_kg` — Mass (kilograms)
- `Tsurf_K` — Surface temperature (Kelvin)
- `Psurf_MPa` — Surface pressure (megapascals)
- `Cmeasured` — Measured moment of inertia factor (C/MR²)
- `Cuncertainty` — Uncertainty in C/MR²
- `Tb_K` — Ocean bottom temperature (Kelvin)
- `J2`, `C22` — Gravity field coefficients

### `Planet.Ocean.*` — Ocean Settings
- `comp` — Composition: `'Seawater'`, `'PureH2O'`, `'NaCl'`, `'MgSO4'`, `'CustomSolution'`
- `wOcean_ppt` — Salinity (parts per thousand)
- `deltaP` — Pressure step (MPa)
- `deltaT` — Temperature step (K)
- `PHydroMax_MPa` — Maximum hydrosphere pressure
- `THydroMax_K` — Maximum hydrosphere temperature

### `Planet.Sil.*` — Silicate Mantle
- `mantleEOS` — Perple_X EOS filename (e.g., `'CM_hydrous_differentiated_*.tab'`)
- `rhoSilWithCore_kgm3` — Silicate density at reference conditions
- `Qrad_Wkg` — Radiogenic heating rate (W/kg)
- `Htidal_Wm3` — Tidal heating rate (W/m³)

### `Planet.Core.*` — Core Properties
- `coreEOS` — Core EOS filename (e.g., `'Fe-S_3D_EOS.mat'`)
- `rhoFe_kgm3`, `rhoFeS_kgm3` — Core material densities
- `xFeS`, `xFeCore` — Iron-sulfur composition fractions
- `QScore` — Core heat flux

### `Planet.Steps.*` — Discretization
- `nIceI` — Number of Ice I steps
- `nSilMax` — Maximum silicate layer steps
- `nCore` — Core layer steps

---

## EOS (Equation of State) System

### Caching Mechanism
- EOSs stored in `EOSlist.loaded` dictionary (defined in `defineStructs.py`)
- Key prevents recomputation of expensive EOS evaluations
- **Critical**: Clear between tests to prevent memory bloat

### Ocean EOS Hierarchy
1. **SeaFreeze** — `'PureH2O'`, `'NaCl'` (ice polymorphs + liquid)
2. **GSW** — `'Seawater'` (standard seawater via TEOS-10)
3. **Custom MgSO4** — `'MgSO4'` (magnesium sulfate solutions)
4. **Reaktoro** — `'CustomSolution'` (arbitrary compositions via Frezchem/SUPCRT databases)

### Silicate/Core EOS
- **Perple_X** outputs: Pre-computed phase diagrams and properties
- Location: `PlanetProfile/Thermodynamics/EOSdata/Perple_X/`
- Downloaded during installation (~164 MB)
- Examples: `CM_hydrous_differentiated_*.tab`, `Fe-S_3D_EOS.mat`

### Reference Profiles
- Cached in user directory via `platformdirs.user_cache_dir('PlanetProfile', 'PlanetProfile')`
- Location stored in `_REFPROFILESCACHE` (see `__init__.py`)

---

## Exploration Modes

### InductOgram
- Explore magnetic induction responses across parameter space
- Function: `InductOgram()` in `Main.py`
- Output: Magnetic field amplitudes, phases at different frequencies
- Uses: `MagneticInduction()` repeatedly with varying parameters

### ExploreOgram
- 2D exploration across two model properties (e.g., ice thickness vs. salinity)
- Function: `ExploreOgram()` in `Main.py`
- Output: Contour plots of interior properties
- Results saved as `.pkl` and `.mat` files

### Monte Carlo (Work in Progress)
- Statistical sampling of model parameter space
- Located in `PlanetProfile/MonteCarlo/`
- **Not fully tested** — use with caution

---

## Parallel Processing

- Uses Python `multiprocessing` module
- Platform-specific context:
  - Windows: `'spawn'`
  - Linux/macOS: `'fork'`
- Gated by `Params.DO_PARALLEL` flag
- Set `DO_PARALLEL = False` to disable (useful for debugging)

---

## Testing Structure

### Test Suites
- **`Testing.py`** (primary) — Run from repo root
- **`PlanetProfile/BuildTest.py`** (alternative) — Run as module

### Test Files
- Location: `PlanetProfile/Test/`
- Naming: `PPTest1.py` through `PPTest28.py`
- Special tests:
  - `PPTestBayes.py` — Bayesian inference
  - `PPTest*InductOgram.py` — Magnetic induction explorations
  - `PPTest*Explore.py` — ExploreOgram tests

### Test Pattern
Each test defines a `Planet` object with specific configuration, then runs `PlanetProfile(Planet, Params)`.

### Adding Tests
When porting a feature:
1. Create `PPTest<N+1>.py` in `PlanetProfile/Test/`
2. Define `Planet` with feature-specific configuration
3. Import and run in `BuildTest.py` or `Testing.py`

---

## Key Dependencies

From `pyproject.toml`:

### Required (Thermodynamics/Geophysics)
- `numpy>=1.26.4,<2.0.0`
- `scipy>=1.16.3,<1.17.0`
- `SeaFreeze>=1.0.0` — Ice/ocean thermodynamics
- `gsw>=3.6.20` — Seawater (TEOS-10)
- `MoonMag>=1.7.5` — Magnetic induction
- `pyalma3>=1.0.1` — Tidal Love numbers

### Optional/Conditional
- `reaktoro` — Custom ocean compositions (not in pyproject.toml, installed separately)
- `obspy` — Seismic wave propagation (optional)
- `spiceypy>=8.0.0` — SPICE kernels for trajectory analysis

### Visualization
- `matplotlib>=3.10.6`
- `cmasher>=1.9.2` — Colormaps

### Utilities
- `mpmath>=1.3.0` — Arbitrary precision math
- `hdf5storage>=0.2.2` — MATLAB file I/O
- `platformdirs>=3.0.0` — Cross-platform cache directories

---

## Common Patterns and Idioms

### Logging
```python
import logging
log = logging.getLogger('PlanetProfile')

log.debug("Detailed diagnostic")
log.info("General information")
log.warning("Something unexpected")
log.error("Error occurred")
```

### Phase Indexing
```python
from PlanetProfile.Utilities.Indexing import GetPhaseIndices

indsLiq, indsI, indsII, indsIII, indsV, indsVI, indsSil, indsFe = GetPhaseIndices(Planet.phase)
```

### EOS Access
```python
# Ocean EOS
from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS

Planet.Ocean.EOS = GetOceanEOS(comp, wOcean_ppt, P_MPa, T_K, ...)

# Use EOS
rho = Planet.Ocean.EOS.fn_rho_kgm3(P_MPa, T_K)
Cp = Planet.Ocean.EOS.fn_Cp_JkgK(P_MPa, T_K)
```

### Timing
```python
from PlanetProfile.Utilities.defineStructs import Timing

Timing.setStartingTime(time.time())
# ... do work ...
Timing.printFunctionTimeDifference('MyFunction()', time.time())
```

---

## File Format Conventions

### Output Files
- **`.pkl`** — Python pickle (primary format for reloading)
- **`.mat`** — MATLAB format (for MATLAB interoperability)
- **`.txt`** — Human-readable text (summary data)

### Figures
- Default: `.png` at 300 DPI
- Can also save as `.pdf`, `.eps`, `.fig` (controlled by `Params.GRIDS_ONLY`, `Params.vectorFormats`)

---

## Documentation Build

### Prerequisites
```bash
pip install sphinx sphinxcontrib-matlabdomain sphinxcontrib.apidoc sphinx-rtd-theme myst-parser
```

### Build Process
```bash
cd docs
rm -rf stubs/      # Clear old API docs
make clean         # Clear build artifacts
make html          # Generate HTML docs
```

### View Results
Open `docs/_build/html/index.html` in browser.

### Sphinx Configuration
- `docs/conf.py` — Main config
- `docs/index.rst` — Documentation homepage
- Auto-generates API docs from NumPy-style docstrings

---

## Version and Release Management

### Version Numbers
- `pyproject.toml` line 7: `version = "X.X.X"`
- `PlanetProfile/Utilities/PPverNum.txt` — MATLAB compatibility

### Building Package
```bash
rm -rf dist/*
python -m build
```

Creates:
- `dist/PlanetProfile-X.X.X.tar.gz` (source distribution)
- `dist/PlanetProfile-X.X.X-py3-none-any.whl` (wheel)

### PyPI Upload
```bash
python -m twine upload dist/* --verbose
```

Requires PyPI API token for authentication.

---

## Scientific References

Key papers are listed in `REFERENCES.md`:
- Vance et al. (2018) — Geophysical investigations of habitability (JGR Planets)
- Styczinski et al. (2023) — PlanetProfile software description (Earth and Space Science)

When porting features that implement models from papers:
1. Add paper to `REFERENCES.md` if not present
2. Reference in commit message
3. Add inline comment with citation

---

## Contacts and Resources

- **Documentation:** https://vancesteven.github.io/PlanetProfile
- **Main Repository:** https://github.com/vancesteven/PlanetProfile
- **Mirror:** https://github.com/NASA-Planetary-Science/PlanetProfile
- **PyPI Package:** https://pypi.org/project/PlanetProfile/
- **Slack:** Development discussions (invitation on request)
