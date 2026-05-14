# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PlanetProfile is a software framework for constructing 1D interior structure models for ocean worlds and rocky dwarf planets based on planetary properties. It uses self-consistent thermodynamics for fluid, rock, and mineral phases, and computes outputs including sound speeds, attenuation, and electrical conductivities.

**Key capabilities:**
- Self-consistent ocean world modeling with various ocean compositions (pure water, NaCl, seawater, MgSO4, or custom via Reaktoro)
- Interior modeling for silicates and cores using PerpleX equation of states
- **High-pressure ice convection with partial melt prediction** (Kalousova et al. 2018) for Ice III, V, VI — see [KALOUSOVA_IMPLEMENTATION_GUIDE.md](KALOUSOVA_IMPLEMENTATION_GUIDE.md)
- Forward modeling of tidal Love numbers (PyALMA3 or TidalPy) and magnetic induction responses (MoonMag, included in the repository)
- **MCMC Bayesian inference** for rheology constraints using TidalPy (Tests 41–44, Andrade/Maxwell rheologies, Arrhenius viscosity) — see [MCMC_INFERENCE_GUIDE.md](MCMC_INFERENCE_GUIDE.md)
- Large-scale explorations (ExploreOgrams, MCMC inference, parameter sweeps)
- Multiple output formats (.txt, .pkl, .mat)

**Supported bodies:** Europa, Ganymede, Callisto, Titan, Enceladus, Ariel, Miranda, Titania, Oberon, Mimas, Tethys, Dione, Rhea, Iapetus, Luna, Io, Pluto, Triton

## Development Environment

### Python Version
- **Required:** Python 3.8–3.11
- **Recommended (users and developers):** Python 3.11
- Python 3.14 support is planned but not yet available

### Installation Commands

```bash
# Install prerequisites (use conda/mamba environment)
conda install numpy=1.26.4 scipy matplotlib mpmath pandas
conda install -c conda-forge gsw obspy spiceypy cmasher reaktoro
pip install SeaFreeze hdf5storage PyALMA3 TidalPy

# Clone repo and install
python -m PlanetProfile.install PPinstall
```

## Running PlanetProfile

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

# Run exploreogram (2D parameter space exploration; requires PPBodyExplore.py)
python -m PlanetProfile.Main Europa exploreogram
```

### As a Python Module
```bash
python -m PlanetProfile.Main Europa
```
```python
from PlanetProfile.Main import RunPPfile
RunPPfile('Europa', 'PPEuropa.py')
```

### GUI Application (PlanetProfileApp)

PlanetProfileApp is a Streamlit-based graphical interface that simplifies model configuration and execution.

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

**Page flow:** About → Main Settings → Bulk Planetary Settings → Ocean Settings → Core and Silicate Settings → Layer Step Settings → Run PlanetProfile → Outputs → Exploreogram.

**Architecture:**
- `PlanetProfileApp.py`: entry point with page navigation via `st.navigation()`
- `pages/*.py`: individual page implementations
- `Utilities/`:
  - `session_manager.py`: session save/load (`SessionManager` class, stored in `sessions/`)
  - `app_helpers.py`: UI helpers and runtime estimation
  - `get_planet.py`: planet object management across pages
  - `planet_sidebar.py`: sidebar status display
  - `CustomPlanetGenerator.py`: custom body creation
  - `presets.py`: preset configuration loader
  - `help_system.py`: in-app help

**Features:**
- Session save/load (JSON-serializable values only; planet objects need special handling)
- Progress indicators across a 6-step workflow
- Interactive parameter adjustment with real-time monitoring
- Figure viewing with PDF/image conversion
- **Exploreogram**: 15+ parameters, configurable grid (2×2 to 50×50), multiple output variables, parallel processing, HTML/PKL output

**Important:** The GUI uses `st.session_state` to maintain state across pages. All planet configuration is stored in session state and written to config files when running simulations.

## Testing

### Run Full Test Suite
```bash
# From main PlanetProfile directory
python -m PlanetProfile.BuildTest
```

Runs all test bodies in `PlanetProfile/Test/PPTest*.py` (excluding InductOgram, Bayes, and Explore tests). **Run before any commit.**

### Adding New Tests
When adding major new functionality, create a new `PPTest#.py` in `PlanetProfile/Test/` following existing patterns.

### GUI Testing
```bash
python PlanetProfileApp/test_utilities.py
streamlit run PlanetProfileApp/PlanetProfileApp.py
```
Focus: session state persistence, parameter validation, config file generation, integration with `PlanetProfile.Main`.

## Architecture

### Two Interfaces: CLI vs GUI

1. **CLI (`PlanetProfileCLI.py`)**: direct Python execution for batch processing, scripting, programmatic access. Users edit configuration files manually.
2. **GUI (`PlanetProfileApp/`)**: Streamlit interface that builds config files from web forms and calls the same core functions.

**Critical:** Both ultimately call `PlanetProfile.Main.PlanetProfile()`. The GUI is a wrapper that simplifies configuration and provides visualization.

### Configuration System (Critical Pattern)

PlanetProfile uses a **hierarchical configuration override system**:

1. **Default configs** (lowest priority): `PlanetProfile/Default/Body/PPBody.py`
2. **Default general settings**: `PlanetProfile/defaultConfig.py`, `PlanetProfile/*/defaultConfig*.py`
3. **User configs** (highest priority): `Body/PPBody.py` in working directory

**Important:** Default files in `PlanetProfile/Default/` should only be changed for development. User experimentation belongs in the working directory configs copied by `python -m PlanetProfile.install PPinstall`.

### Main Entry Point

`PlanetProfile/Main.py` contains:
- `PlanetProfile()`: main function for single body models
- `InductOgram()`: 2D exploration of magnetic induction
- `ExploreOgram()`: 2D exploration across parameter space
- `RunPPfile()`: convenience wrapper for running body configs

### Module Structure

- **`Thermodynamics/`**: core physics calculations
  - `LayerPropagators.py`: `IceLayers`, `OceanLayers`, `InnerLayers` — propagate models through planetary layers
  - `HydroEOS.py`: ocean EOS calculations
  - `InnerEOS.py`: silicate/core EOS
  - `Electrical.py`: electrical conductivity
  - `Seismic.py`: seismic properties
  - `OceanProps.py`: liquid ocean property calculations
  - `Reaktoro/`: custom solution thermodynamics via geochemical modeling

- **`Plotting/`**: visualization
  - `ProfilePlots.py`: interior structure plots
  - `MagPlots.py`: magnetic induction plots
  - `ExplorationPlots.py`: parameter exploration plots

- **`MagneticInduction/`**: magnetic field calculations using MoonMag
  - `MagneticInduction.py`: core induction calculations
  - `Moments.py`: excitation moment calculations

- **`Gravity/`**: gravitational field calculations
  - `Gravity.py`: Love numbers and gravity parameters

- **`Utilities/`**: helper functions and data structures
  - `defineStructs.py`: core data structures (`Constants`, `PlanetStruct`, etc.)
  - `SetupInit.py`: initialization and setup functions
  - `ResultsIO.py`: I/O for model results

- **`Default/`**: default configuration files for each body

### Data Flow

1. Body configuration file (e.g., `PPEuropa.py`) defines planetary parameters and settings
2. `Main.py` imports config and calls `PlanetProfile()`
3. `SetupInit()` initializes `Planet` structure with constants and parameters
4. Layer propagators (`IceLayers`, `OceanLayers`, `InnerLayers`) calculate properties from surface to center
5. Secondary calculations: electrical conductivity, seismic properties, magnetic induction
6. Results saved via `WriteResults()`; plots generated via `GeneratePlots()`

## High-Pressure Ice Convection (Kalousova et al. 2018)

Implementation of the Kalousova et al. (2018) HP ice convection model with partial melt prediction for Ice III, V, and VI layers.

### Quick Start

```python
# In PPBody.py or config file
Planet.Do.KALOUSOVA_CONVECTION = True  # Enable Kalousova model (default: False)
```

### Key Features (current implementation status)

- **Temperate-layer detection** (Ra\* > Ra\*c) for Ice III, V, VI: implemented and identical across phases. `ConvectionKalousova2018` returns `Tconv_K, etaMelt_Pas, eLid_m, Dconv_m, deltaTBL_m, Qbot_W, Ra*, RaCrit` for all three.
- **Thermal-profile propagation** (lid + adiabatic interior + lower TBL): fully implemented for Ice **III** and **V**. Ice **VI** currently uses a simplified placeholder — uniform T along the melting curve (T = `TconvVI_K`) across the layer — pending the lid+adiabat+TBL port. The full propagation requires `Planet.Steps.nVbottom`, `nIceVILitho`, layer allocation in `IceLayers()`, and `IceVIConductSolid/Porous`. Tracked as a follow-up.
- **Triple-point temperature boundaries**: interior temperatures approach Ice III–V–liquid (254 K) and V–VI–liquid (272 K) triple points.
- **Per-phase configuration**: enable/disable independently for Ice III, V, VI.
- **Melt fraction**: `Planet.meltFraction{III,V,VI}` is reported but **not** the output of a two-phase solver. Top-level dispatchers (`Convection.py::Ice*ConvectSolid`) use a fixed placeholder when a temperate layer is detected: 0.01 for Ice III/V, 0.5 for Ice VI (placeholder pending the full solver). The in-ocean path (`LayerPropagators.py::IceShellHydroSolid` and porous variant) sets melt fraction = `Constants.phiPercolation` (= 0.05) and uses that value inside Kalousova & Sotin 2018 Eq. 10 to compute outflow velocity and mass flux. **`Constants.phiPercolation` should not be edited** — Kalousova's outflow equations are conditioned on it.

### Configuration Options

```python
# Disable per phase
Planet.Do.NO_ICE_CONVECTION_III = True
Planet.Do.NO_ICE_CONVECTION_V = True
Planet.Do.NO_ICE_CONVECTION_VI = True

# Or disable all ice convection
Planet.Do.NO_ICE_CONVECTION = True
```

### Output Variables

```python
if Planet.DO_HP_MELT:
    print(f"Ice III melt fraction: {Planet.meltFractionIII}")
    print(f"Ice V melt fraction: {Planet.meltFractionV}")
    print(f"Ice VI melt fraction: {Planet.meltFractionVI}")

# Temperate layer thicknesses (m)
Planet.eLidIII_m, Planet.eLidV_m, Planet.eLidVI_m

# Convection parameters
Planet.TconvIII_K, Planet.TconvV_K, Planet.TconvVI_K
Planet.RaConvectIII, Planet.RaConvectV, Planet.RaConvectVI
```

### Implementation Details

**Full technical documentation:** [KALOUSOVA_IMPLEMENTATION_GUIDE.md](KALOUSOVA_IMPLEMENTATION_GUIDE.md)

**Key files:**
- `PlanetProfile/Thermodynamics/ThermalProfiles/ThermalProfiles.py:70-214` — `ConvectionKalousova2018()`
- `PlanetProfile/Thermodynamics/ThermalProfiles/Convection.py` — Dispatch logic: Ice III (L463), Ice V (L714), Ice VI (L1173)
- `PlanetProfile/Utilities/defineStructs.py` — Config flags and struct fields

**Model comparison:**

| Feature | Deschamps & Sotin (2001) | Kalousova et al. (2018) |
|---------|--------------------------|-------------------------|
| Applicable phases | Ice I, III, V | HP ice only (III, V, VI) |
| Interior temperature | Adiabatic | Melting curve (solidus) |
| Melt formation | Not modeled | Predicted via Ra\* > Ra\*c |
| Stagnant lid | Top conductive layer | Top temperate layer |

**Reference:** Kalousová, K., & Sotin, C. (2018). Melting in high-pressure ice layers of large ocean worlds — Implications for volatiles transport. *Geophysical Research Letters*, 45(16), 8096–8103. https://doi.org/10.1029/2018GL078889

## Key Patterns and Conventions

### Streamlit Session State (GUI)

```python
# In every page, retrieve planet object
from Utilities.get_planet import get_planet
planet = get_planet()  # Gets st.session_state.get("Planet")

# Store configuration values
st.session_state['parameter_name'] = value

# Access across pages
value = st.session_state.get('parameter_name', default_value)
```

- `SessionManager` handles save/load (only JSON-serializable values; planet objects need special handling)
- Each page is a separate Python file in `pages/`; navigation configured in `PlanetProfileApp.py` via `st.navigation()`
- Use `show_planet_status()` for sidebar display; `create_progress_indicator()` for step progress

### Parallel Processing

Python's `multiprocessing` module is used for parallel computations. Disable for cross-platform debugging:

```python
# In UserConfigs/configPP.py or body config
Params.DO_PARALLEL = False
```
- Windows: `spawn`
- macOS/Linux: `fork`

### Equation of State (EOS) Data
- Perple_X EOS files (~164 MB) auto-downloaded during installation
- Located in: `PlanetProfile/Thermodynamics/EOSdata/Perple_X/`
- Generated with Perple_X v6.7.9

### CALC_NEW Flags
Models use `CALC_NEW` flags to control recalculation vs cached reload. Recalculate all parameters whenever:
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
- `pyproject.toml`: package version (line 7)
- `PlanetProfile/Utilities/PPverNum.txt`: version number file
- `PlanetProfile/Utilities/PPversion.py`: version retrieval functions

### Configuration Files
- `configPP.py`: general PlanetProfile settings
- `configPPinduct.py`: magnetic induction settings
- `configPPplots.py`: plotting settings
- `configPPgravity.py`: gravity calculation settings
- `configPPmodel.py`: model-specific settings
- `configPPcustomsolution.py`: Reaktoro custom solution settings

### Body-Specific Files
Each body directory (e.g., `Europa/`) contains:
- `PPBody.py`: main configuration file
- `PPBodyInductOgram.py`: magnetic induction exploration (optional)
- `figures/`: output plots
- `inductionData/`, `seismicData/`: data files (if applicable)

## Documentation

- Main documentation: https://vancesteven.github.io/PlanetProfile
- Build docs: `cd docs && rm -rf stubs/ && make clean && make html`
- Requires: `pip install sphinx sphinxcontrib-matlabdomain sphinxcontrib.apidoc sphinx-rtd-theme myst-parser`

## Common Issues

### SeaFreeze Version
If using SeaFreeze < v1.0.0, manually remove old installation:
```bash
python -m site  # find site-packages directory
# Delete: seafreeze.py, SeaFreeze_Gibbs.mat, SeaFreeze* directories
pip install SeaFreeze
```

### Perple_X Data Not Found
```bash
python -m PlanetProfile.install
```

### PlanetProfileApp Issues

**Poppler not found (PDF conversion fails):**
- macOS: `conda install poppler`
- Windows: download from https://github.com/oschwartz10612/poppler-windows and add `bin/` to PATH
- Linux: `apt-get install poppler-utils` or equivalent

**Session state not persisting:**
- Check `PlanetProfileApp/sessions/` exists
- Verify values are JSON-serializable
- Use `SessionManager.save_session()` / `load_session()`

**Planet not loading in pages:**
- Ensure planet selected in Main Settings page
- Check `st.session_state.get("Planet")` is not None
- Use `get_planet()` helper from `Utilities/get_planet.py`

**Streamlit errors on page navigation:**
- Each page must call `st.set_page_config()` if used (first Streamlit command)
- Session state keys must be consistent across pages
- Avoid `st.stop()` except when truly needed

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
├── Utilities/                   # GUI-specific utilities
├── sessions/                    # Saved user sessions
└── config*.py                   # Config files for GUI runs

[Working directories for each body]/
├── PPBody.py                    # User configuration (overrides defaults)
├── configPP.py                  # User general settings
├── figures/                     # Output plots
└── [model outputs].txt/pkl/mat
```

---

# Working with Claude (operational rules)

## Environment
- mamba environment: `PPcl` — activate with `mamba activate PPcl`
- Never use system python
- Working directory: `~/Dropbox/planetprofile-genai/`
- Manuscripts: `~/Dropbox/manuscripts/` — read-only unless explicitly told to commit/push

## Session start protocol
At every session start, in order:
1. `graphify query "inference mcmc forward model parameter space"` — orient to current module structure
2. `cd ~/llm-wiki && python3 -m llmwiki query "PlanetProfile sbi PlanetProfileApp 3D viscosity convection" && cd -`
3. Read `README.md`

## graphify usage
Knowledge graph lives at `graphify-out/`. AST-only build, zero API cost. Measured average: **~7k tokens/query, ~15× reduction vs naive file reads** (`graphify benchmark`, 2026-05-14; per-query range 6×–270× depending on specificity).

- **Use `graphify query "topic"`** for architectural questions ("how does X flow", "what connects A to B") and symbol location across modules
- **Use grep/ripgrep directly** for specific symbol lookups in a known single module — often faster than the graph at that granularity
- After modifying code, run `graphify update .` to refresh (AST-only, free)
- `GRAPH_REPORT.md` is currently ~25 KB and readable for orientation — but if it ever exceeds ~100 KB, `.graphifyignore` is failing and the graph needs rebuilding
- Do **NOT** modify `.graphifyignore` without user confirmation — clustering quality depends on it

## Model usage
- **Default**: Opus 4.7 for primary reasoning (current actual usage)
- **Route down to Haiku 4.5 for**: file reads, grep/search, git status, classifying test failures, summarizing known content, simple pattern-following edits
- **Route down to local Gemma 4 (Ollama, via MCP wrapper)** for: log summarization, manuscript Q&A, batch text processing — when wired up
- **There is no further escalation** — Opus is the ceiling. If a task feels too hard for Opus, the problem is the prompt, not the model
- For any bug fix, read the relevant file with Haiku first — never guess at API names from memory

## Cost / context discipline
On Opus, every wasted token costs ~5× what it would on Sonnet. The discipline below matters more, not less.

- Summarize context before passing to subagents — never pass full conversation history
- Avoid re-reading files already in context
- Do not spawn a subagent for tasks completable in <5 tool calls
- Prefer sequential tool calls over parallel agents unless tasks are truly independent
- Prefer `graphify query` (~7k tokens) over Grep over full file reads, in that order — except where grep is the right tool (specific symbol in known module)
- **Biggest trap**: conflating "I need to understand this file" with "I need to read it myself." A Haiku summary is usually sufficient input for Opus reasoning — delegate the read before opening the file yourself

## Hard constraints
- Do **NOT** scan, grep, or modify Perple_X directories (large lookup tables)
- Do **NOT** modify `PPTest*` files numbered ≤28 unless explicitly instructed
- Do **NOT** modify `.claudeignore` or `.graphifyignore`
- Do **NOT** commit or push without explicit instruction

## Autonomy and checkpoints

The model: Claude runs autonomously **between two human checkpoints** — planning and review.

### Pre-execution checkpoint (planning)
For any non-trivial change (>1 file or >50 LOC, or any new logic), produce a plan **before writing code**:
- Files to touch, files NOT to touch
- Approach and key design decisions
- Tests that will exercise the change
- Success criteria

Wait for approval before implementing. For trivial bug fixes within a known single module, skip the plan and proceed.

### Autonomous execution
- Implement, run tests, self-debug
- **Re-plan and confirm** if scope expands materially beyond what was approved (e.g., new file types, new dependencies, touching modules outside the approved set)
- Stop and report if blocked on ambiguity — do not guess

### Post-execution checkpoint (review)
Before yielding control, produce a concise summary:
- Files changed (one line of intent per file)
- Test results (pass/fail counts, any new failures)
- Any deviations from the approved plan
- Open questions or follow-ups

Then **STOP**. Do not commit or push.

### Autonomy by module
- `PlanetProfile/Inference/*` — implement approved designs and fix bugs autonomously
- `PPTest*` files numbered >40 — autonomous fixes only, no new logic
- Core PP modules (`Thermodynamics/`, `Main.py`, `Plotting/`, `MagneticInduction/`, `Gravity/`) — always plan first
- `PlanetProfileApp/` — plan first for new pages or session-state changes; autonomous for cosmetic / help-text edits
- Git operations — always ask before commit/push, no exceptions

### Multi-agent failure modes to watch for
- **Context collapse**: subagents losing original task framing → keep summaries explicit and bounded
- **Specification gaming**: tests passing without real correctness → verify outputs against physical expectations (energy balance, monotonicity, known limits), not just exit codes
- **Goal misgeneralization**: agent over-extending scope → trigger re-planning when scope grows

If any of these fire, stop and surface to user before continuing.
