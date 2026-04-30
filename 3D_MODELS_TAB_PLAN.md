# Plan: 3D Models Tab for PlanetProfileApp

## Overview

Add a "3D Models" page to PlanetProfileApp below the Exploreogram tab. This page will let users construct spatially resolved (lat/lon) interior models by running PlanetProfile column-by-column across a spherical grid, with laterally varying ice thickness and tidal heating. The page leverages the existing `PlanetProfile/Lateral/` module.

## Prerequisites

A completed 1D PlanetProfile run is required as the reference model. The 3D tab takes the reference Planet object from session state and extends it to a spatial grid.

## Page Location

- **File:** `PlanetProfileApp/pages/ThreeDModels.py`
- **Navigation registration** in `PlanetProfileApp.py`:
  ```python
  three_d_models = st.Page("pages/ThreeDModels.py", title="3D Models", icon="🧊")
  ```
  Added after `exploreogram` in the navigation list.

## Page Layout

### Section 1: Grid Configuration

| Parameter | Widget | Default | Notes |
|-----------|--------|---------|-------|
| Grid type | `st.selectbox` | `"latlon"` | `"latlon"` or `"healpix"` |
| nLat | `st.number_input` | 36 | Latitude resolution (latlon only) |
| nLon | `st.number_input` | 72 | Longitude resolution (latlon only) |
| nSide | `st.number_input` | 8 | HEALPix resolution (healpix only; nPix = 12 * nSide^2) |

### Section 2: Ice Thickness Variation

Two modes:

**A. Spherical harmonic coefficients** (default for tidal bodies):
- Degree-2 ice thickness variation from tidal equilibrium
- User specifies amplitude `dIce_km` (peak-to-peak variation)
- Internally converted to C20, C22 coefficients via Ojakangas & Stevenson (1989) pattern
- Shows preview of ice thickness map before running

**B. Uniform** (no lateral variation):
- Uses reference model ice thickness everywhere
- Still computes spatially varying tidal heating

### Section 3: Tidal Heating Configuration

| Parameter | Widget | Default | Notes |
|-----------|--------|---------|-------|
| Heating pattern | `st.selectbox` | `"Ojakangas & Stevenson (1989)"` | Degree-2 spatial pattern |
| Rheology | `st.selectbox` | `"Maxwell"` | `"Maxwell"` or `"Andrade"` |
| Mean Htidal (W/m^3) | `st.number_input` | from reference Planet | Mean volumetric heating |
| Andrade alpha | `st.slider` | 0.2 | Only shown for Andrade |

### Section 4: Run Controls

- **Run button**: Launches `RunLateralColumns()` with progress bar
- **Parallel toggle**: Enable/disable multiprocessing
- **Estimated runtime** based on grid size and reference model time
- Progress: "Column X/N (lat, lon)" with Streamlit progress bar

### Section 5: Results Display

Three sub-tabs within results:

**A. Surface Maps** (Plotly interactive heatmaps):
- Ice shell thickness (km)
- Basal temperature (K)
- Surface heat flux (mW/m^2)
- Tidal heating distribution (W/m^3)
- Ocean thickness (km)
- Convective layer thickness (km)

**B. Cross Sections**:
- Equatorial slice (lat=0): radial profile of T, rho, phase vs longitude
- Polar slice (lon=0): radial profile vs latitude
- Uses Plotly heatmap with phase boundaries overlaid

**C. Statistics**:
- Table of min/mean/max for each field
- Histogram of ice thickness distribution
- Comparison to reference 1D model

### Section 6: Export

- Save results as PKL (full columnPlanets array)
- Save maps as interactive HTML (Plotly)
- Save maps as PNG/PDF (matplotlib fallback)
- CSV export of per-column summary table (lat, lon, zIce, Tb, qSurf, ...)

## Backend Integration

### Existing modules to use:

| Module | Function | Purpose |
|--------|----------|---------|
| `Lateral.LateralStructure` | `InitLateralStructure()` | Set up spatial grid and ice thickness field |
| `Lateral.LateralStructure` | `RunLateralColumns()` | Per-column PlanetProfile execution |
| `Lateral.TidalHeating3D` | `TidalStrainPattern()` | Compute spatial heating pattern |
| `Lateral.TidalHeating3D` | `_MaxwellDissipation()` | Maxwell rheology dissipation |
| `Lateral.TidalHeating3D` | `_AndradeDissipation()` | Andrade rheology dissipation |
| `Lateral.SpatialGrid` | `InitGrid()` | Create lat/lon or HEALPix grid |
| `Lateral.SpatialGrid` | `SHtoGrid()` | Evaluate SH on grid |
| `Plotting.LateralPlots` | `_PlotSurfaceMap()` | Generate lat/lon maps |

### New code needed:

1. **`ThreeDModels.py`** (~500-700 lines): The Streamlit page itself
2. **Plotly wrappers** for interactive maps: `_PlotlyHeatmap()` utility within the page or in a shared plotting module, following the Exploreogram pattern for Plotly vs matplotlib switching
3. **Cross-section extraction**: Helper to slice columnPlanets array along latitude or longitude and assemble radial profiles into 2D arrays

### Session state variables:

```python
st.session_state['3d_grid_type']          # 'latlon' or 'healpix'
st.session_state['3d_nLat']               # int
st.session_state['3d_nLon']               # int
st.session_state['3d_nSide']              # int
st.session_state['3d_dIce_km']            # float, ice thickness amplitude
st.session_state['3d_rheology']           # 'Maxwell' or 'Andrade'
st.session_state['3d_running']            # bool
st.session_state['3d_results']            # dict with columnPlanets, grid, maps
st.session_state['3d_error']              # str or None
st.session_state['3d_error_traceback']    # str or None
```

## Implementation Phases

### Phase 1: Core page + grid setup + tidal heating preview
- Page skeleton with navigation integration
- Grid configuration UI
- Tidal heating pattern preview (just the spatial pattern, no PlanetProfile runs)
- Ice thickness SH preview map

### Phase 2: Column execution + surface maps
- Wire `RunLateralColumns()` with progress tracking
- Surface map display (Plotly interactive)
- Statistics table
- PKL/CSV export

### Phase 3: Cross sections + advanced visualization
- Equatorial and polar cross-section views
- Phase boundary overlays
- Comparison with 1D reference

### Phase 4: Polish
- Session save/load integration
- Presets for common bodies (Europa, Enceladus, Titan)
- Runtime estimation
- Help text integration

### Phase 5: Magnetic Induction on 3D Models

#### Overview

After computing a spatially resolved interior model, run MoonMag magnetic induction to predict induced magnetic field signatures from the laterally varying conductivity structure. This extends the existing `PlanetProfile/MagneticInduction/` module to work with 3D conductivity distributions, not just 1D radial profiles.

#### Conductivity Structure from 3D Model

Each column in the 3D model produces a 1D radial conductivity profile (σ(r) at that lat/lon). The collection of profiles across the grid defines a 3D conductivity distribution σ(r, θ, φ). To use MoonMag:

1. **Expand conductivity to spherical harmonic coefficients**: At each radial shell, the laterally varying σ values on the grid are decomposed into SH coefficients (Clm, Slm) up to a user-specified maximum degree.
2. **Higher-order SH support**: The current Lateral module focuses on degree-2 (C20, C22) for ice thickness. Magnetic induction supports arbitrary SH degree, limited only by computational cost and available constraints.
3. **User-supplied coefficients**: Users can supply arbitrary SH coefficients for the conductivity structure rather than deriving them from column models.

#### UI Elements

| Parameter | Widget | Default | Notes |
|-----------|--------|---------|-------|
| Run induction | `st.checkbox` | False | Enable magnetic induction computation |
| Max SH degree | `st.number_input` | 4 | Maximum spherical harmonic degree for conductivity expansion |
| Excitation frequencies | `st.multiselect` | body-dependent | Orbital, synodic, etc. from existing MoonMag config |
| Coefficient source | `st.radio` | `"From 3D model"` | `"From 3D model"`, `"Manual entry"`, `"From file"` |
| Manual coefficients | `st.data_editor` | — | Shown when "Manual entry" selected; editable table of (l, m, Clm, Slm) |
| Coefficient file upload | `st.file_uploader` | — | Shown when "From file" selected; CSV/JSON with SH coefficients |

#### Results Display

**D. Magnetic Induction** (new sub-tab within results):
- Induced field maps at each excitation frequency (Bx, By, Bz components)
- Comparison of 3D-derived vs 1D-derived induced fields
- SH spectrum of induced response (power per degree)
- Amplitude and phase maps of induced dipole, quadrupole, octupole

#### Backend Integration

| Module | Function | Purpose |
|--------|----------|---------|
| `MagneticInduction.MagneticInduction` | `MagneticInduction()` | Core induction calculation (1D per column) |
| `MagneticInduction.Moments` | `Excitations()` | Compute excitation moments for the body |
| `MagneticInduction.MoonMag/` | vendored MoonMag library | Harmonic field propagation |
| New: `Lateral.ConductivitySH` | `ConductivityToSH()` | Decompose 3D σ(r,θ,φ) into SH coefficients per radial shell |
| New: `Lateral.ConductivitySH` | `SHFromFile()` | Load user-supplied SH coefficients |

#### Body-Specific Constraints

Available higher-order constraints by body:
- **Europa**: Degree-2 from ice shell thickness variation; higher degrees from crater/geology constraints
- **Ganymede**: Limited to degree-2 from tidal equilibrium; ionospheric conductivity from Hubble UV
- **Enceladus**: Strong degree-2 + 3 from south polar thermal anomaly
- **Titan**: Degree-2 from Cassini gravity; potential higher-order from topography

#### Implementation Notes

- MoonMag currently operates on 1D radial conductivity profiles. For true 3D induction, either:
  (a) Run MoonMag per-column and combine (approximate, fast), or
  (b) Extend MoonMag to accept SH conductivity coefficients directly (exact, requires MoonMag modification)
- Start with approach (a) for Phase 5; approach (b) is a stretch goal
- SH decomposition of grid data uses `scipy.special.sph_harm` or `pyshtools` if available

## Key Design Decisions

1. **Plotly-first visualization**: Interactive maps are critical for spatial data. Matplotlib fallback for publication-quality static exports.
2. **Reference model required**: The 3D tab reuses the 1D reference from session state (from the "Run PlanetProfile" page). Users must run a 1D model first.
3. **Progressive disclosure**: Start with simple defaults (latlon 36x72, tidal equilibrium ice, Maxwell), show advanced options in expanders.
4. **Computation cost**: At 36x72 = 2592 columns, each taking ~1-5s, total runtime is 30-120 min. The page should support background execution and partial results display. Consider a coarser default (18x36 = 648 columns, ~10-30 min).
