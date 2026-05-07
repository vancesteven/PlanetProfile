# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PlanetProfile is a software framework for constructing 1D interior structure models for ocean worlds and rocky dwarf planets based on planetary properties. It uses self-consistent thermodynamics for fluid, rock, and mineral phases, and computes outputs including sound speeds, attenuation, and electrical conductivities.

**Key capabilities:**
- Self-consistent ocean world modeling with various ocean compositions (pure water, NaCl, seawater, MgSO4, or custom via Reaktoro)
- Interior modeling for silicates and cores using PerpleX equation of states
- **High-pressure ice convection with partial melt prediction** (Kalousova et al. 2018) for Ice III, V, VI — see [KALOUSOVA_IMPLEMENTATION_GUIDE.md](KALOUSOVA_IMPLEMENTATION_GUIDE.md)
- Forward modeling of tidal Love numbers (PyALMA3 or TidalPy) and magnetic induction responses (MoonMag, included in the repository)
- **MCMC Bayesian inference** for rheology constraints using TidalPy (Tests 41-44, Andrade/Maxwell rheologies, Arrhenius viscosity) — see [MCMC_INFERENCE_GUIDE.md](MCMC_INFERENCE_GUIDE.md)
- Large-scale explorations (ExploreOgrams, MCMC inference, parameter sweeps)
- Multiple output formats (.txt, .pkl, .mat)

**Supported bodies:** Europa, Ganymede, Callisto, Titan, Enceladus, Ariel, Miranda, Titania, Oberon, Mimas, Tethys, Dione, Rhea, Iapetus, Luna, Io, Pluto, Triton

## Development Environment

### Python Version
- **Required:** Python 3.8-3.11
- **Recommended for users:** Python 3.11
- **Required for developers:** Python 3.11
- Python 3.14 support is planned but not yet available

### Installation Commands

For development work:
```bash
# Install prerequisites (use conda/mamba environment)
conda install numpy=1.26.4 scipy matplotlib mpmath pandas
conda install -c conda-forge gsw obspy spiceypy cmasher reaktoro
pip install SeaFreeze hdf5storage PyALMA3 TidalPy

# Clone repo and install
python -m PlanetProfile.install PPinstall
```

## Running PlanetProfile

### GUI Application (PlanetProfileApp)

PlanetProfileApp is a Streamlit-based graphical interface for PlanetProfile that simplifies model configuration and execution.

**Installation:**
```bash
pip install streamlit pdf2image Pillow pandas
conda install poppler  # macOS
# Windows: Download Poppler from https://github.com/oschwartz10612/poppler-windows
#          and add its /bin folder to PATH
```

**Running the GUI:**
```bash
# From PlanetProfileApp directory
streamlit run PlanetProfileApp.py
```

**GUI Structure:**
The app uses a multi-page navigation system:
1. **About** - Introduction and documentation
2. **PlanetProfile Main Settings** - Select planet or create custom body
3. **Bulk Planetary Settings** - Mass, radius, and bulk properties
4. **Ocean Settings** - Configure ocean composition (predefined or custom salts)
5. **Core and Silicate Settings** - Core and mantle configuration
6. **Layer Step Settings** - Set simulation resolution/granularity
7. **Run PlanetProfile** - Execute simulation with progress tracking
8. **PlanetProfile Outputs** - View figures and data tables
9. **Exploreogram** - 2D parameter space exploration (NEW)

**Key GUI Features:**
- Session management: Save/load configurations via `SessionManager` class
- Progress indicators for 6-step workflow
- Real-time simulation monitoring
- Interactive parameter adjustment
- Figure viewing with PDF/image conversion
- Preset configurations via `Utilities/presets.py`
- Help system integrated throughout (`Utilities/help_system.py`)
- **Exploreogram**: Interactive 2D parameter space exploration with heatmap visualization
  - 15+ exploration parameters (ocean salinity, temperature, core composition, etc.)
  - Configurable grid resolution (2×2 to 50×50)
  - Multiple output variables (ice thickness, temperature, magnetic field, etc.)
  - Parallel processing for efficiency
  - Save results as HTML (interactive) or PKL (data)

**GUI Architecture:**
- `PlanetProfileApp.py`: Main entry point with page navigation
- `pages/*.py`: Individual page implementations
- `Utilities/`:
  - `session_manager.py`: Session persistence
  - `app_helpers.py`: UI helpers and runtime estimation
  - `get_planet.py`: Planet object management across pages
  - `planet_sidebar.py`: Sidebar status display
  - `CustomPlanetGenerator.py`: Custom body creation
  - `presets.py`: Preset configuration loader

**Important:** The GUI uses Streamlit's `st.session_state` to maintain state across pages. All planet configuration is stored in session state and written to config files when running simulations.

### Command Line Interface
```bash
# Run with body name
python PlanetProfileCLI.py Europa

# Run with specific configuration file
python PlanetProfileCLI.py path/to/PPBody.py

# Run with multiple files
python PlanetProfileCLI.py Europa PPEuropa1.py PPEuropa2.py

# Reload from saved results
python PlanetProfileCLI.py Europa reload PPEuropa.txt

# Run exploreogram (2D parameter space exploration)
# Requires PPBodyExplore.py configuration file
python -m PlanetProfile.Main Europa exploreogram
```

### As a Python Module
```bash
# Run as module
python -m PlanetProfile.Main Europa

# Or in Python code
from PlanetProfile.Main import RunPPfile
RunPPfile('Europa', 'PPEuropa.py')
```

## Testing

### Run Full Test Suite
```bash
# From main PlanetProfile directory
python -m PlanetProfile.BuildTest
```

This runs all test bodies in `PlanetProfile/Test/PPTest*.py` (excluding InductOgram, Bayes, and Explore tests).

**Before any commits:** Run the full test suite to ensure no functionality is broken.

### Adding New Tests
When adding major new functionality, create a new `PPTest#.py` file in `PlanetProfile/Test/` directory following the existing test patterns.

### GUI Testing
For PlanetProfileApp development:
```bash
# Run basic utility tests
python PlanetProfileApp/test_utilities.py

# Run Streamlit app locally for testing
streamlit run PlanetProfileApp/PlanetProfileApp.py
```

The GUI testing focuses on:
- Session state persistence
- Parameter validation
- Configuration file generation
- Integration with PlanetProfile.Main

## Architecture

### Two Interfaces: CLI vs GUI

PlanetProfile has two main interfaces:

1. **Command Line (PlanetProfileCLI.py)**: Direct Python execution for batch processing, scripting, and programmatic access. Users manually edit configuration files.

2. **GUI (PlanetProfileApp/)**: Streamlit web interface for interactive configuration. The GUI:
   - Builds configuration files from user input through web forms
   - Calls the same `PlanetProfile.Main` functions as CLI
   - Provides visual feedback and output display
   - Saves/loads sessions for reproducibility

**Critical:** Both interfaces ultimately call the same core `PlanetProfile.Main.PlanetProfile()` function. The GUI is a wrapper that simplifies configuration file creation and provides a visual interface.

### Configuration System (Critical Pattern)

PlanetProfile uses a **hierarchical configuration override system**:

1. **Default configs** (lowest priority): `PlanetProfile/Default/Body/PPBody.py`
2. **Default general settings**: `PlanetProfile/defaultConfig.py`, `PlanetProfile/*/defaultConfig*.py`
3. **User configs** (highest priority): `Body/PPBody.py` in working directory

**Important:** Default files in `PlanetProfile/Default/` should only be changed for development. User experimentation should happen in the working directory configs copied by `python -m PlanetProfile.install PPinstall`.

### Main Entry Point

`PlanetProfile/Main.py` is the primary entry point containing:
- `PlanetProfile()`: Main function for single body models
- `InductOgram()`: 2D exploration of magnetic induction
- `ExploreOgram()`: 2D exploration across parameter space
- `RunPPfile()`: Convenience wrapper for running body configs

### Module Structure

- **`Thermodynamics/`**: Core physics calculations
  - `LayerPropagators.py`: IceLayers, OceanLayers, InnerLayers - propagates models through planetary layers
  - `HydroEOS.py`: Ocean equation of state calculations
  - `InnerEOS.py`: Silicate/core equation of state
  - `Electrical.py`: Electrical conductivity calculations
  - `Seismic.py`: Seismic properties
  - `OceanProps.py`: Liquid ocean property calculations
  - `Reaktoro/`: Custom solution thermodynamics via geochemical modeling

- **`Plotting/`**: Visualization modules
  - `ProfilePlots.py`: Interior structure plots
  - `MagPlots.py`: Magnetic induction plots
  - `ExplorationPlots.py`: Parameter exploration plots

- **`MagneticInduction/`**: Magnetic field calculations using MoonMag
  - `MagneticInduction.py`: Core induction calculations
  - `Moments.py`: Excitation moment calculations

- **`Gravity/`**: Gravitational field calculations
  - `Gravity.py`: Love numbers and gravity parameters

- **`Utilities/`**: Helper functions and data structures
  - `defineStructs.py`: Core data structures (Constants, PlanetStruct, etc.)
  - `SetupInit.py`: Initialization and setup functions
  - `ResultsIO.py`: Input/output for model results

- **`Default/`**: Default configuration files for each body
  - Each body has `PPBody.py` configuration file

### Data Flow

1. Body configuration file (e.g., `PPEuropa.py`) defines planetary parameters and settings
2. `Main.py` imports config and calls `PlanetProfile()`
3. `SetupInit()` initializes Planet structure with constants and parameters
4. Layer propagators (`IceLayers`, `OceanLayers`, `InnerLayers`) calculate properties from surface to center
5. Secondary calculations: electrical conductivity, seismic properties, magnetic induction
6. Results saved via `WriteResults()` and plots generated via `GeneratePlots()`

## High-Pressure Ice Convection (Kalousova et al. 2018)

PlanetProfile includes an implementation of the Kalousova et al. (2018) HP ice convection model with partial melt prediction for Ice III, V, and VI layers.

### Quick Start

Enable in configuration files:
```python
# In PPBody.py or config file
Planet.Do.KALOUSOVA_CONVECTION = True  # Enable Kalousova model (default: False)
```

### Key Features

- **Two-phase convection parameterization**: Predicts partial melting in HP ice layers when convection is vigorous
- **Temperate layer detection**: Uses modified Rayleigh number criterion (Ra* > Ra*c) to identify melt-bearing regions
- **Triple-point temperature boundaries**: Interior temperatures approach Ice III-V-liquid (254 K) and V-VI-liquid (272 K) triple points
- **Per-phase configuration**: Can enable/disable convection independently for Ice III, V, VI

### Configuration Options

```python
# Enable/disable per phase
Planet.Do.NO_ICE_CONVECTION_III = True   # Suppress Ice III convection
Planet.Do.NO_ICE_CONVECTION_V = True     # Suppress Ice V convection
Planet.Do.NO_ICE_CONVECTION_VI = True    # Suppress Ice VI convection

# Or disable all ice convection
Planet.Do.NO_ICE_CONVECTION = True       # All phases
```

### Output Variables

After running with Kalousova convection:
```python
# Check for partial melt
if Planet.DO_HP_MELT:
    print(f"Ice III melt fraction: {Planet.meltFractionIII}")
    print(f"Ice V melt fraction: {Planet.meltFractionV}")
    print(f"Ice VI melt fraction: {Planet.meltFractionVI}")
    
# Temperate layer thicknesses
Planet.eLidIII_m   # Ice III temperate layer thickness (m)
Planet.eLidV_m     # Ice V temperate layer thickness (m)
Planet.eLidVI_m    # Ice VI temperate layer thickness (m)

# Convection parameters
Planet.TconvIII_K, Planet.TconvV_K, Planet.TconvVI_K  # Interior temperatures
Planet.RaConvectIII, Planet.RaConvectV, Planet.RaConvectVI  # Rayleigh numbers
```

### Implementation Details

**Full technical documentation:** [KALOUSOVA_IMPLEMENTATION_GUIDE.md](KALOUSOVA_IMPLEMENTATION_GUIDE.md)

**Key files:**
- `PlanetProfile/Thermodynamics/ThermalProfiles/ThermalProfiles.py:70-214` — `ConvectionKalousova2018()` function
- `PlanetProfile/Thermodynamics/ThermalProfiles/Convection.py` — Dispatch logic for Ice III (line 463), Ice V (line 714), Ice VI (line 1173)
- `PlanetProfile/Utilities/defineStructs.py` — Config flags and struct fields

**Model comparison:**

| Feature | Deschamps & Sotin (2001) | Kalousova et al. (2018) |
|---------|--------------------------|-------------------------|
| Applicable phases | Ice I, III, V | HP ice only (III, V, VI) |
| Interior temperature | Adiabatic | Melting curve (solidus) |
| Melt formation | Not modeled | Predicted via Ra* > Ra*c |
| Stagnant lid | Top conductive layer | Top temperate layer |

**Reference:** Kalousová, K., & Sotin, C. (2018). Melting in high-pressure ice layers of large ocean worlds—Implications for volatiles transport. *Geophysical Research Letters*, 45(16), 8096-8103. https://doi.org/10.1029/2018GL078889

## Key Patterns and Conventions

### Streamlit Session State (GUI)

PlanetProfileApp uses Streamlit's `st.session_state` to maintain user configuration across pages:

**Critical pattern:**
```python
# In every page, retrieve planet object
from Utilities.get_planet import get_planet
planet = get_planet()  # Gets st.session_state.get("Planet")

# Store configuration values
st.session_state['parameter_name'] = value

# Access across pages
value = st.session_state.get('parameter_name', default_value)
```

**Session persistence:**
- `SessionManager` class (`Utilities/session_manager.py`) handles saving/loading
- Sessions stored in `PlanetProfileApp/sessions/` directory
- Only JSON-serializable values are saved
- Planet objects and complex structures require special handling

**Page navigation:**
- Each page is a separate Python file in `pages/`
- Navigation configured in `PlanetProfileApp.py` using `st.navigation()`
- Use `show_planet_status()` to display current planet in sidebar
- Progress indicators via `create_progress_indicator()` helper

### Parallel Processing
PlanetProfile uses Python's `multiprocessing` module for parallel computations. Disable if encountering cross-platform issues:
```python
# In UserConfigs/configPP.py or body config
Params.DO_PARALLEL = False
```

Platform-specific multiprocessing:
- Windows: uses 'spawn'
- macOS/Linux: uses 'fork'

### Equation of State (EOS) Data
- Perple_X EOS files (~164 MB) are automatically downloaded during installation
- Located in: `PlanetProfile/Thermodynamics/EOSdata/Perple_X/`
- Generated with Perple_X v6.7.9

### CALC_NEW Flags
Models use `CALC_NEW` flags to control whether to recalculate or reload cached results. It's recommended to recalculate all parameters whenever:
- PlanetProfile is updated
- Input parameters change that may affect layer thicknesses
- Intermediate variables may be affected

### Output Files
Models generate multiple output files in the body directory:
- `.txt`: ASCII results for reloading
- `.pkl`: Python pickle files with Planet objects
- `.mat`: MATLAB-compatible files
- Figures in `Body/figures/`

## Important Files and Locations

### Version Information
- `pyproject.toml`: Package version (line 7)
- `PlanetProfile/Utilities/PPverNum.txt`: Version number file
- `PlanetProfile/Utilities/PPversion.py`: Version retrieval functions

### Configuration Files
- `configPP.py`: General PlanetProfile settings
- `configPPinduct.py`: Magnetic induction settings
- `configPPplots.py`: Plotting settings
- `configPPgravity.py`: Gravity calculation settings
- `configPPmodel.py`: Model-specific settings
- `configPPcustomsolution.py`: Reaktoro custom solution settings

### Body-Specific Files
Each body directory (e.g., `Europa/`) contains:
- `PPBody.py`: Main configuration file
- `PPBodyInductOgram.py`: Magnetic induction exploration (optional)
- `figures/`: Output plots
- `inductionData/`: Magnetic field data (if applicable)
- `seismicData/`: Seismic data (if applicable)

## Documentation

- Main documentation: https://vancesteven.github.io/PlanetProfile
- Build docs: `cd docs && rm -rf stubs/ && make clean && make html`
- Requires: `pip install sphinx sphinxcontrib-matlabdomain sphinxcontrib.apidoc sphinx-rtd-theme myst-parser`

## Common Issues

### SeaFreeze Version
If using SeaFreeze < v1.0.0, manually remove old installation:
```bash
# Find site-packages directory
python -m site
# Delete: seafreeze.py, SeaFreeze_Gibbs.mat, SeaFreeze* directories
# Then: pip install SeaFreeze
```

### Perple_X Data Not Found
Run installation script to download EOS data:
```bash
python -m PlanetProfile.install
```

### PlanetProfileApp Issues

**Poppler not found (PDF conversion fails):**
- macOS: `conda install poppler`
- Windows: Download from https://github.com/oschwartz10612/poppler-windows and add bin/ to PATH
- Linux: `apt-get install poppler-utils` or equivalent

**Session state not persisting:**
- Check that `PlanetProfileApp/sessions/` directory exists
- Verify values are JSON-serializable (no custom objects without conversion)
- Use `SessionManager.save_session()` and `SessionManager.load_session()`

**Planet not loading in pages:**
- Ensure planet selected in Main Settings page
- Check `st.session_state.get("Planet")` is not None
- Use `get_planet()` helper function from `Utilities/get_planet.py`

**Streamlit errors on page navigation:**
- Each page must call `st.set_page_config()` if used (first Streamlit command)
- Session state keys must be consistent across pages
- Avoid using `st.stop()` except when truly needed to halt execution

## Repository Structure

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
│   └── Reaktoro/
├── Plotting/                    # Visualization
├── MagneticInduction/           # Magnetic field modeling
├── Gravity/                     # Gravitational calculations
├── Utilities/                   # Helper functions
└── Test/                        # Test configurations
    ├── PPTest1.py ... PPTest28.py

PlanetProfileApp/                # Streamlit GUI application
├── PlanetProfileApp.py          # GUI entry point
├── pages/                       # Multi-page app structure
│   ├── About.py
│   ├── PlanetProfileMainSettings.py
│   ├── BulkPlanetarySettings.py
│   ├── OceanSettings.py
│   ├── CoreSettings.py
│   ├── LayerStepSettings.py
│   ├── RunPlanetProfile.py
│   ├── PlanetProfileOutputs.py
│   └── CompareRuns.py
├── Utilities/                   # GUI-specific utilities
│   ├── session_manager.py       # Session save/load
│   ├── app_helpers.py           # UI components
│   ├── help_system.py           # Help documentation
│   ├── presets.py               # Preset configurations
│   ├── get_planet.py            # Planet state management
│   └── CustomPlanetGenerator.py
├── sessions/                    # Saved user sessions
└── config*.py                   # Config files for GUI runs

[Working directories for each body]/
├── PPBody.py                    # User configuration (overrides defaults)
├── configPP.py                  # User general settings
├── figures/                     # Output plots
└── [model outputs].txt/pkl/mat
```

## Cost discipline
- Default model: sonnet for all agents and subagents
- Never switch to opus without explicit per-task instruction from user
- Minimize subagent spawning — prefer sequential tool calls 
  over parallel agents unless tasks are truly independent
- Summarize context before passing to subagents — never 
  pass full conversation history
- Avoid re-reading files already in context
- Do not spin up a subagent for tasks completable in <5 tool calls
- Prefer grep/glob over full file reads when locating symbols
- If requirements are ambiguous, ask one clarifying question 
  before executing — do not assume and proceed

## graphify

This project has a graphify knowledge graph at graphify-out/.

graphify-out/GRAPH_REPORT.md is a navigation index only (500 community stubs) — do NOT read it.
Use `graphify query "topic"` instead: returns relevant nodes with file+line locations for ~30 tokens.

Rules:
- Before using Grep or reading multiple files to locate a symbol or module, run
  `graphify query "topic"` first — it replaces broad file searches at a fraction of the cost
- For targeted bug fixes within a known single module, read that file directly
- After modifying code files in this session, run `graphify update .` to keep the graph current (AST-only, no API cost)
- Do NOT read graphify-out/GRAPH_REPORT.md (it is a 42K-token Obsidian index with no useful content)


## Constraints
- Do NOT implement code unless explicitly instructed — summarize and propose first
- Perple_X directories: large lookup tables — DO NOT scan, grep, or modify
- Do not modify PPTest files numbered <28 unless explicitly instructed
- Do NOT re-run graphify without user confirmation
- Do not modify .claudeignore or .graphifyignore

## Directory boundaries
- All work in ~/Dropbox/planetprofile-genai/
- Manuscript: ~/Dropbox/manuscripts/ — read only unless explicitly told to commit/push

## Session resume protocol
Always run at every session start (no exceptions):
1. Run: `graphify query "inference mcmc forward model parameter space"` to orient to current module structure
2. Run: `cd ~/llm-wiki && python3 -m llmwiki query "PlanetProfile sbi PlanetProfileApp 3D viscosity convection" && cd -`
3. Read README.md

## Environment
- Always use the mamba environment: PPcl
- Activate with: mamba activate PPcl
- Never use system python

## Cost discipline
- Minimize subagent spawning — prefer sequential tool calls
  over parallel agents unless tasks are truly independent
- Summarize context before passing to subagents — never
  pass full conversation history
- Avoid re-reading files already in context
- Do not spin up a subagent for tasks completable in <5 tool calls
- If requirements are ambiguous, ask one clarifying question
  before executing — do not assume and proceed

## Model usage
- Default: Sonnet for all tasks and subagents
- Haiku for: file reads, grep/search, simple pattern-following
  edits, git status, summarizing known content
- Sonnet for: new logic, multi-file architecture, debugging
- If a task warrants Opus, say so and wait for confirmation 
  before proceeding
- Opus only when explicitly instructed per-task
- For any bug fix, read the relevant file with Haiku first —
  never guess at API names

## Autonomy levels
- Inference module (PlanetProfile/Inference/*): fix bugs
  and implement approved designs autonomously
- PPTest files >40: autonomous fixes only, no new logic
- Core PP modules: always summarize before changing
- Git operations: always ask before commit/push
