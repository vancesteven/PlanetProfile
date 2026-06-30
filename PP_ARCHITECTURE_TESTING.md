# PP_ARCHITECTURE_TESTING.md

Repository architecture, test organization, and benchmark guidance for
PlanetProfile.

Operational rules (autonomy, planning, scientific correctness, commit
discipline) live in `CLAUDE.md`. This file is the technical companion:
how the code is laid out, how to run it, and how it is tested.

---

## Project overview

PlanetProfile constructs 1-D interior structure models for ocean worlds
and rocky dwarf planets from planetary properties, using self-consistent
thermodynamics for fluid, rock, and mineral phases. Outputs include
sound speeds, attenuation, electrical conductivities, tidal Love
numbers, and magnetic induction responses.

**Capabilities:**

- Self-consistent ocean modeling for pure water, NaCl, seawater, MgSO4,
  or custom compositions via Reaktoro.
- Silicate and core EOS via PerpleX lookups.
- High-pressure ice convection with partial-melt prediction
  (Kalousova et al. 2018) for Ice III, V, VI.
- Forward modeling of tidal Love numbers (PyALMA3 or TidalPy) and
  magnetic induction (MoonMag, vendored).
- MCMC Bayesian inference for rheology constraints (TidalPy-backed,
  Andrade/Maxwell rheologies, Arrhenius viscosity).
- Large-scale parameter sweeps via ExploreOgrams.
- Output formats: `.txt`, `.pkl`, `.mat`.

**Supported bodies:** Europa, Ganymede, Callisto, Titan, Enceladus,
Ariel, Miranda, Titania, Oberon, Mimas, Tethys, Dione, Rhea, Iapetus,
Luna, Io, Pluto, Triton.

---

## Development environment

- **Python:** 3.10–3.12 (3.11 recommended).
- **Required env:** `mamba activate PPcl` — never use system Python.
  See `CLAUDE.md` for the operational rule.

### Installation

```bash
conda install "numpy>=2.0" scipy matplotlib mpmath pandas
conda install -c conda-forge gsw obspy spiceypy cmasher reaktoro
pip install SeaFreeze hdf5storage PyALMA3 TidalPy
python -m PlanetProfile.install PPinstall
```

---

## Running PlanetProfile

### CLI

```bash
python PlanetProfileCLI.py Europa                      # by body name
python PlanetProfileCLI.py path/to/PPBody.py           # by config file
python PlanetProfileCLI.py Europa PPEuropa1.py PPEuropa2.py
python PlanetProfileCLI.py Europa reload PPEuropa.txt  # reload saved
python -m PlanetProfile.Main Europa exploreogram       # 2-D sweep
```

### As a module

```bash
python -m PlanetProfile.Main Europa
```

```python
from PlanetProfile.Main import RunPPfile
RunPPfile('Europa', 'PPEuropa.py')
```

### GUI (PlanetProfileApp)

Streamlit-based interface. Install:

```bash
pip install streamlit pdf2image Pillow pandas
conda install poppler            # macOS
# Windows: install Poppler for Windows and add /bin to PATH
```

Run:

```bash
streamlit run PlanetProfileApp/PlanetProfileApp.py
```

**Page flow:** About → Main Settings → Bulk Planetary Settings →
Ocean Settings → Core and Silicate Settings → Layer Step Settings →
Run PlanetProfile → Outputs → Exploreogram.

**Architecture:**

- `PlanetProfileApp.py` — entry point with page navigation via
  `st.navigation()`.
- `pages/*.py` — individual page implementations.
- `Utilities/`:
  - `session_manager.py` — session save/load (`SessionManager` class,
    sessions stored under `sessions/`).
  - `app_helpers.py` — UI helpers and runtime estimation.
  - `get_planet.py` — planet object management across pages.
  - `planet_sidebar.py` — sidebar status display.
  - `CustomPlanetGenerator.py` — custom body creation.
  - `presets.py` — preset configuration loader.
  - `help_system.py` — in-app help.

State is held in `st.session_state`; planet configuration is stored
there and written to config files when running simulations.
`SessionManager` only persists JSON-serializable values; planet objects
need special handling.

---

## Testing

### Full suite

```bash
python -m PlanetProfile.BuildTest
```

Runs every body test in `PlanetProfile/Test/PPTest*.py` excluding
InductOgram, Bayes, and Explore tests. Run before any commit.

### Conventions

- Numbered `PPTest*.py` files (`PPTest1.py … PPTest28.py`) are part of
  the regression suite. Per `CLAUDE.md`, do **not** modify these files
  numbered ≤28 without explicit approval.
- Higher-numbered `PPTest*` files (>40) host MCMC / inference
  development tests; bug fixes in those are autonomous, but new logic
  still requires planning per `CLAUDE.md`.
- Add new bodies / scientific cases by creating a new `PPTest#.py`
  following existing patterns.

### GUI tests

```bash
python PlanetProfileApp/test_utilities.py
streamlit run PlanetProfileApp/PlanetProfileApp.py
```

Focus on session-state persistence, parameter validation, config-file
generation, and integration with `PlanetProfile.Main`.

### Benchmark / regression discipline

Per `CLAUDE.md`, passing tests is necessary but not sufficient for
PlanetProfile. When numerical results change, identify whether the
cause is a physical/modeling change, a bug fix, a default change, a
tolerance change, a data/benchmark change, or an unresolved issue.

Where applicable, check:

- unit and dimensional consistency,
- expected monotonicity or limiting behavior,
- benchmark agreement,
- reproducibility with fixed seeds,
- likelihood / prior behavior for inference,
- physically plausible outputs (e.g., Love-number signs and ranges,
  density monotonicity, energy balance),
- manuscript-reportable limitations.

For stochastic methods, fix random seeds where practical, report
sampling settings, and separate stochastic variability from
code-induced changes. Do not loosen tolerances without numerical
justification.

---

## Architecture

### Two interfaces: CLI vs GUI

1. **CLI (`PlanetProfileCLI.py`)** — direct Python execution for batch
   processing, scripting, and programmatic access. Users edit
   configuration files manually.
2. **GUI (`PlanetProfileApp/`)** — Streamlit interface that builds
   config files from web forms and calls the same core functions.

Both ultimately invoke `PlanetProfile.Main.PlanetProfile()`. The GUI is
a wrapper that simplifies configuration and provides visualization.

### Configuration system

Hierarchical override:

1. **Default configs (lowest priority):** `PlanetProfile/Default/Body/PPBody.py`.
2. **Default general settings:** `PlanetProfile/defaultConfig.py` and
   `PlanetProfile/*/defaultConfig*.py`.
3. **User configs (highest priority):** `Body/PPBody.py` in the user's
   working directory.

Defaults under `PlanetProfile/Default/` should only be changed for
development. User experimentation belongs in working-directory configs
copied by `python -m PlanetProfile.install PPinstall`.

### Main entry point

`PlanetProfile/Main.py` exposes:

- `PlanetProfile()` — main function for single-body models.
- `InductOgram()` — 2-D exploration of magnetic induction.
- `ExploreOgram()` — 2-D exploration across parameter space.
- `RunPPfile()` — convenience wrapper for running body configs.

### Module structure

- **`Thermodynamics/`** — core physics.
  - `LayerPropagators.py` — `IceLayers`, `OceanLayers`, `InnerLayers`
    propagate models through planetary layers.
  - `HydroEOS.py` — ocean EOS calculations.
  - `InnerEOS.py` — silicate / core EOS.
  - `Electrical.py` — electrical conductivity (and `PhaseConv` —
    canonical phase-code → label mapping).
  - `Seismic.py` — seismic properties.
  - `OceanProps.py` — liquid ocean property calculations.
  - `Reaktoro/` — custom solution thermodynamics via geochemical
    modeling.
  - `ThermalProfiles/` — thermal-profile and convection logic
    (including the Kalousova HP-ice model — see below).
  - `RefProfiles/` — cached reference profiles (ASCII).

- **`Plotting/`** — visualization.
  - `ProfilePlots.py` — interior structure plots.
  - `MagPlots.py` — magnetic induction plots.
  - `ExplorationPlots.py` — parameter exploration plots.

- **`MagneticInduction/`** — magnetic field calculations using MoonMag.
  - `MagneticInduction.py` — core induction calculations.
  - `Moments.py` — excitation moment calculations.

- **`Gravity/`** — gravitational field.
  - `Gravity.py` — Love numbers and gravity parameters.

- **`Inference/`** — MCMC and Bayesian inference toolkit.
  - `inference_core.py` — `InferenceConfig` and runner glue.
  - `cache_builder.py` — pre-builds per-Tb structure caches consumed by
    the forward models.
  - `forward_models.py` — k₂ / h₂ and induction forward models that
    operate against cached structures.
  - `mcmc_plots.py` — corner plots, layer-vs-D_ocean diagnostics,
    structure wedges.
  - `probe_tb_prior.py` — diagnostic probe over the cache to find
    Tb_K ranges where the integrator succeeds.
  - `find_tb_bounds.py` — bisection probe of any body's PP template to
    find Tb_K satisfying geological D_iceIh / D_ocean targets.
  - `configs/*.json` — per-body inference configurations (priors,
    fixed parameters, observables, derived parameters).

- **`Utilities/`** — helpers and data structures.
  - `defineStructs.py` — core data structures (`Constants`,
    `PlanetStruct`, etc.).
  - `SetupInit.py` — initialization and setup functions.
  - `ResultsIO.py` — I/O for model results.

- **`Default/`** — default configuration files for each body.

### Data flow (single-body run)

1. Body configuration file (e.g., `PPEuropa.py`) defines planetary
   parameters and settings.
2. `Main.py` imports the config and calls `PlanetProfile()`.
3. `SetupInit()` initializes the `Planet` structure with constants and
   parameters.
4. Layer propagators (`IceLayers`, `OceanLayers`, `InnerLayers`)
   compute properties from the surface inward.
5. Secondary calculations: electrical conductivity, seismic
   properties, magnetic induction.
6. Results saved via `WriteResults()`; plots generated via
   `GeneratePlots()`.

---

## High-pressure ice convection (Kalousova et al. 2018)

Implementation of the Kalousova et al. 2018 HP-ice convection model
with partial-melt prediction for Ice III, V, and VI.

### Quick start

```python
Planet.Do.KALOUSOVA_CONVECTION = True   # default: False
```

### Status

- **Temperate-layer detection (Ra* > Ra*c)** for Ice III, V, VI:
  implemented and identical across phases. `ConvectionKalousova2018`
  returns `Tconv_K, etaMelt_Pas, eLid_m, Dconv_m, deltaTBL_m, Qbot_W,
  Ra*, RaCrit` for all three.
- **Thermal-profile propagation** (lid + adiabatic interior + lower
  TBL): fully implemented for Ice **III** and **V**. Ice **VI**
  currently uses a uniform-T placeholder (T = `TconvVI_K`) along the
  melting curve, pending the lid+adiabat+TBL port. Tracked as a
  follow-up; full propagation requires `Planet.Steps.nVbottom`,
  `nIceVILitho`, layer allocation in `IceLayers()`, and
  `IceVIConductSolid/Porous`.
- **Triple-point boundaries:** interior temperatures approach
  Ice III–V–liquid (≈254 K) and V–VI–liquid (≈272 K) triple points.
- **Per-phase configuration:** enable / disable independently for
  Ice III, V, VI.
- **Melt fraction:** `Planet.meltFraction{III,V,VI}` is reported but is
  **not** the output of a two-phase solver. Top-level dispatchers
  (`Convection.py::Ice*ConvectSolid`) use a fixed placeholder when a
  temperate layer is detected: 0.01 for Ice III/V, 0.5 for Ice VI
  (placeholder, pending the full solver). The in-ocean path
  (`LayerPropagators.py::IceShellHydroSolid` and porous variant) sets
  melt fraction = `Constants.phiPercolation` (= 0.05) and uses that
  value inside Kalousova & Sotin 2018 Eq. 10 for outflow velocity and
  mass flux. **`Constants.phiPercolation` should not be edited** —
  Kalousova's outflow equations are conditioned on it.

### Configuration

```python
# Disable per phase
Planet.Do.NO_ICE_CONVECTION_III = True
Planet.Do.NO_ICE_CONVECTION_V = True
Planet.Do.NO_ICE_CONVECTION_VI = True

# Or disable all ice convection
Planet.Do.NO_ICE_CONVECTION = True
```

### Outputs

```python
if Planet.DO_HP_MELT:
    Planet.meltFractionIII, Planet.meltFractionV, Planet.meltFractionVI
Planet.eLidIII_m, Planet.eLidV_m, Planet.eLidVI_m       # lid thicknesses (m)
Planet.TconvIII_K, Planet.TconvV_K, Planet.TconvVI_K    # convection temps
Planet.RaConvectIII, Planet.RaConvectV, Planet.RaConvectVI
```

### Key files

- `PlanetProfile/Thermodynamics/ThermalProfiles/ThermalProfiles.py` —
  `ConvectionKalousova2018()`.
- `PlanetProfile/Thermodynamics/ThermalProfiles/Convection.py` —
  dispatch logic for Ice III / V / VI.
- `PlanetProfile/Utilities/defineStructs.py` — config flags and
  struct fields.

### Model comparison

| Feature              | Deschamps & Sotin (2001) | Kalousova et al. (2018)        |
|----------------------|--------------------------|--------------------------------|
| Applicable phases    | Ice I, III, V            | HP ice only (III, V, VI)       |
| Interior temperature | Adiabatic                | Melting curve (solidus)        |
| Melt formation       | Not modeled              | Predicted via Ra* > Ra*c       |
| Stagnant lid         | Top conductive layer     | Top temperate layer            |

**Reference:** Kalousová, K., & Sotin, C. (2018). Melting in
high-pressure ice layers of large ocean worlds — Implications for
volatiles transport. *Geophysical Research Letters*, 45(16),
8096–8103. https://doi.org/10.1029/2018GL078889

---

## Key patterns and conventions

### Streamlit session state (GUI)

```python
from Utilities.get_planet import get_planet
planet = get_planet()                          # st.session_state.get("Planet")

st.session_state['parameter_name'] = value
value = st.session_state.get('parameter_name', default_value)
```

- `SessionManager` handles save / load (only JSON-serializable values;
  planet objects need special handling).
- Each page is a separate Python file in `pages/`; navigation is
  configured in `PlanetProfileApp.py` via `st.navigation()`.
- Use `show_planet_status()` for sidebar display and
  `create_progress_indicator()` for step progress.

### Parallel processing

Python's `multiprocessing` is used for parallel computations. Disable
for cross-platform debugging:

```python
# In UserConfigs/configPP.py or body config
Params.DO_PARALLEL = False
```

- Windows uses `spawn`; macOS / Linux use `fork`.

### Equation of state (EOS) data

- Perple_X EOS files (~164 MB) are auto-downloaded during installation.
- Located in `PlanetProfile/Thermodynamics/EOSdata/Perple_X/`.
- Generated with Perple_X v6.7.9.

### `CALC_NEW` flags

Models use `CALC_NEW` flags to control recalculation vs cached reload.
Recalculate all parameters whenever:

- PlanetProfile is updated.
- Input parameters change that may affect layer thicknesses.
- Intermediate variables may be affected.

### Output files

Models generate multiple outputs in the body directory:

- `.txt` — ASCII results for reloading.
- `.pkl` — Python pickle with `Planet` objects.
- `.mat` — MATLAB-compatible.
- Figures under `Body/figures/`.

---

## Important files and locations

### Version

- `pyproject.toml` — package version (line 7).
- `PlanetProfile/Utilities/PPverNum.txt` — version number file.
- `PlanetProfile/Utilities/PPversion.py` — version retrieval functions.

### Configuration files

- `configPP.py` — general PlanetProfile settings.
- `configPPinduct.py` — magnetic induction settings.
- `configPPplots.py` — plotting settings.
- `configPPgravity.py` — gravity calculation settings.
- `configPPmodel.py` — model-specific settings.
- `configPPcustomsolution.py` — Reaktoro custom solution settings.

### Body-specific files

Each body directory (e.g., `Europa/`) contains:

- `PPBody.py` — main configuration file.
- `PPBodyInductOgram.py` — magnetic induction exploration (optional).
- `figures/` — output plots.
- `inductionData/`, `seismicData/` — data files (when applicable).

---

## Documentation

- Main documentation: https://vancesteven.github.io/PlanetProfile
- Build docs:
  ```bash
  cd docs && rm -rf stubs/ && make clean && make html
  ```
  Requires:
  ```bash
  pip install sphinx sphinxcontrib-matlabdomain sphinxcontrib.apidoc \
              sphinx-rtd-theme myst-parser
  ```

---

## Common issues

### SeaFreeze < v1.0.0

Manually remove the old install before reinstalling:

```bash
python -m site                       # find site-packages directory
# Delete: seafreeze.py, SeaFreeze_Gibbs.mat, SeaFreeze* directories
pip install SeaFreeze
```

### Perple_X data not found

```bash
python -m PlanetProfile.install
```

### PlanetProfileApp issues

**Poppler not found (PDF conversion fails):**

- macOS: `conda install poppler`
- Windows: download from
  https://github.com/oschwartz10612/poppler-windows and add `bin/` to
  PATH.
- Linux: `apt-get install poppler-utils` (or distro equivalent).

**Session state not persisting:**

- Check that `PlanetProfileApp/sessions/` exists.
- Verify values are JSON-serializable.
- Use `SessionManager.save_session()` / `load_session()`.

**Planet not loading in pages:**

- Ensure a planet is selected on the Main Settings page.
- Check `st.session_state.get("Planet")` is not None.
- Use `get_planet()` from `Utilities/get_planet.py`.

**Streamlit errors on page navigation:**

- Each page must call `st.set_page_config()` if used (must be the
  first Streamlit command).
- Session-state keys must be consistent across pages.
- Avoid `st.stop()` except when truly needed.

---

## Repository structure

```
PlanetProfile/
├── Main.py                      # Main entry point
├── GetConfig.py                 # Configuration loading
├── BuildTest.py                 # Test suite runner
├── Default/                     # Default body configurations
│   ├── Europa/PPEuropa.py
│   └── [other bodies]/
├── Thermodynamics/              # Physics calculations
│   ├── LayerPropagators.py
│   ├── HydroEOS.py
│   ├── InnerEOS.py
│   ├── Electrical.py
│   ├── ThermalProfiles/
│   ├── RefProfiles/
│   └── Reaktoro/
├── Plotting/                    # Visualization
├── MagneticInduction/           # Magnetic field modeling
├── Gravity/                     # Gravitational calculations
├── Inference/                   # MCMC and Bayesian inference
│   ├── inference_core.py
│   ├── cache_builder.py
│   ├── forward_models.py
│   ├── mcmc_plots.py
│   ├── probe_tb_prior.py
│   ├── find_tb_bounds.py
│   └── configs/*.json
├── Utilities/                   # Helper functions
└── Test/                        # Test configurations
    ├── PPTest1.py … PPTest28.py
    └── mcmc_results/            # Inference validation artifacts

PlanetProfileApp/                # Streamlit GUI application
├── PlanetProfileApp.py          # GUI entry point
├── pages/                       # Multi-page app
├── Utilities/                   # GUI-specific utilities
├── sessions/                    # Saved user sessions
└── config*.py                   # Config files for GUI runs

[Working directory for each body]/
├── PPBody.py                    # User configuration (overrides defaults)
├── configPP.py                  # User general settings
├── figures/                     # Output plots
└── [model outputs].txt/pkl/mat
```
