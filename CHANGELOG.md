# Changelog

## [Unreleased]

### Enhancements
- **Added 3D lateral structure module with optional healpy dependency.** PlanetProfile now supports spatially-varying (3D) interior structure calculations via the new `PlanetProfile.Lateral` module. Geographic grids can be initialized with HEALPix (equal-area, requires healpy) or lat-lon (simple, no dependencies). Spherical harmonic transforms use 4π-normalization consistent with MoonMag. Enables calculation of H(r,θ,φ) for reproducing geographic tidal heating distributions (Tobie et al. 2005 Figure 10) and asymmetric magnetic induction from spatially-varying ice-ocean boundary topography. **Optional dependency:** healpy>=1.16.0 (install with `pip install PlanetProfile[lateral]` or `conda install -c conda-forge healpy`). Graceful fallback to lat-lon grids if healpy unavailable. New classes: `LateralSubstruct` in `defineStructs.py`, modules: `SpatialGrid.py`, `TidalHeating3D.py`, `LateralStructure.py`, `MassConservation.py`, `ClathrateLateral.py`, `LateralIO.py`, tests: `PPTestLateralPhase1.py`, `PPTestLateralPhase2.py`, `PPTestLateralPhase3.py`, `PPTestLateralPhase4.py`, `PPTestLateralPhase5.py`, `PPTestLateralPhase6.py`, `PPTestLateralPhase7.py`. Literature: Ojakangas & Stevenson (1989), Tobie et al. (2005), Górski et al. (2005), Hobbs (1974), Styczinski et al. (2021). (Phase 1-7 of lateral heating port: +3055 lines, 40 tests)
- **Added TidalPy backend for self-consistent tidal heating calculations.** PlanetProfile now supports TidalPy as an alternative gravity backend alongside PyALMA3, enabling calculation of **Love numbers (k₂, h₂, l₂) AND tidal heating distribution H(r)** from viscoelastic theory. The TidalPy backend computes per-phase volumetric heating rates (ice I, III, V, VI, ocean, silicate) using proper radial integration ∫H(r)·4πr²dr, and supports multiple rheologies (Andrade, Maxwell, Elastic, Newton) with per-phase configuration. Self-consistent thermal-tidal coupling enabled via `Planet.Do.DO_SELF_CONSISTENT_HTIDAL = True`. Backend selection via `Params.Gravity.backend = 'tidalpy'` or `'pyalma'` (default). **Optional dependency:** TidalPy>=0.7.4 (install with `pip install PlanetProfile[tidal]`). Graceful fallback to PyALMA3 if TidalPy not installed. New attributes: `Planet.Bulk.eccentricity`, `Planet.Bulk.meanMotion_radps`, `Planet.Ocean.HtidalIce_Wm3`, `Planet.Gravity.andrade_zeta`, `Planet.Gravity.tidalpy_Htidal_perPhase_W`, `Planet.Gravity.tidalpy_Htidal_perPhase_Wm3`. Literature: Tobie et al. (2003, 2005), Renaud & Henning (2018). (+694 lines in `PlanetProfile/Gravity/Gravity.py`)
- **Added CSV export of profile data.** PlanetProfile now automatically exports profile data to CSV format alongside the full-precision `.txt` files. The CSV file contains 23 columns (P, T, r, phase, ρ, Cp, α, g, φ, σ, k, VP, VS, QS, KS, GS, Ppore, ρmatrix, ρpore, M, V, Htidal, η) with 8 significant figures for easy import into spreadsheet programs (Excel, Google Sheets, etc.). New module: `PlanetProfile/Utilities/HumanReadableOutput.py` with `WriteProfileCSV()` function and optional `WriteProfileSpreadsheet()` for Excel export (requires openpyxl). Original `.txt` files remain unchanged. **New dependency:** pandas>=2.0.0 added to support CSV/Excel export.
- **Renamed water triple point constant for clarity.** `Constants.triplePointT_K` is now `Constants.TtripleIh_III_L_K` (ice Ih-III-liquid triple point at 251.165 K, not the Ih-liquid-vapor point at 273.16 K). The old name remains as a deprecated alias. Added triple point constants for HP ice phase transitions: `TtripleIII_V_L_K` (254 K), `PtripleIII_V_L_MPa` (350 MPa), `TtripleV_VI_L_K` (272 K), `PtripleV_VI_L_MPa` (632 MPa) from Kalousová et al. (2018) and Journaux et al. (2020). Updated `PPTest21.py` to use the new constant name.

### Bug Fixes
- **Fixed missing pMax initialization in inductogram tests.** Added `Planet.Magnetic.pMax = 0` to `PPTest7InductOgram.py`, `PPTest11InductOgram.py`, and `PPTest7Explore.py` to explicitly set spherically symmetric induction (monopole only) for inductogram calculations. Without this initialization, inductogram tests would crash.
- **Replaced deprecated NumPy types for forward compatibility.** Replaced `np.int_` with `np.intp` (18 instances) and `np.float_` with `np.float64` (1 instance) across 7 files. These types were deprecated in NumPy 2.0 and will be removed in NumPy 2.1+. Also fixed typo in `PTPlots.py` line 646 that checked `ymin` twice instead of both `ymin` and `ymax`. Ensures compatibility with future NumPy versions.
- **Corrected conduction profile for clathrate depth calculation.** Added `ConductiveTemperatureCorrect()` function in `PlanetProfile/Thermodynamics/ThermalProfiles/ThermalProfiles.py` implementing the strict Turcotte & Schubert (2002) §4.9 eq. 4.40 form: `c1 = qTop·rTop²/k − Htot·rTop³/(3k)` with planar-limit ΔT = qTop·ΔR/k (Fourier's law). `GetPbConduct()` now uses this function; the resulting clathrate layer thickness matches the user-specified `Bulk.clathMaxThick_m` (previously PP produced a clathrate 2× deeper than requested).
- **Fixed silicate boundary condition for solid-sphere bodies.** `SilRecursionSolid()` and `SilRecursionPorous()` in `PlanetProfile/Thermodynamics/Geophysical.py` now enforce T finite at r=0 for solid-sphere silicate bodies (no Fe core, no CONSTANT_INNER_DENSITY) by overriding `qTop` with `Htot·rTop/3` (the c1=0 value consistent with the T&S 4.40 closed form). Shell bodies (Fe_CORE=True or CONSTANT_INNER_DENSITY=True) retain the legacy halved-c1 behavior to preserve existing test-suite behavior. This resolves the long-standing 1/r divergence issue at the body center for solid-sphere cases.

## [3.1.0] – 2026-01-13
**Author:** @Chang-Scott

### Optimizations
- Improved memory usage in the EOS framework.
- Implemented EOS pre-loading to support large-scale explorations.
- Updated computational pathways to reduce typical PlanetProfile runs to ~1–5 seconds.

### Plotting Extensibility
- Overhauled storage of large-scale exploration results, enabling:
  - Exploreograms to utilize inductogram plots.
  - Inductograms to utilize exploreogram results.
- Exploration data is now saved and reloaded as `.pkl` files instead of `.mat`.
  - `.mat` files are still generated for users who wish to post-process results in MATLAB.

### Monte Carlo
- Added an initial framework to enable Monte Carlo sampling of model properties.
- See individual commits for implementation details.
- *Note:* Monte Carlo functionality is a work in progress and has not yet been extensively tested.

### Non-Self-Consistent Modeling
- Established a framework to support layer modeling with prescribed constant properties rather than EOS-derived values.
  - For example, users can now directly specify ocean density, thermal expansivity, and related parameters without enforcing compositional self-consistency.
- *Note:* This functionality is a work in progress and has not yet been extensively tested.

### MgSO₄–SeaFreeze Coupling
- Implemented dynamic coupling between ice polymorph chemical potentials and MgSO₄ thermodynamic data to generate phase grids on the fly.
- This replaces the previous approach that relied on precomputed MgSO₄ phase lookup tables, which were:
  - Low resolution.
  - Memory intensive.

### Bug Fixes
- Numerous bug fixes addressing edge-case planetary configurations and plotting-related errors.

## [3.0.0] - 2024-05-01
**Author:** @Chang-Scott

This release implements a broad set of changes developed over the past year, with several major new scientific and modeling capabilities.

### Major Additions
- Added the ability to model an **NaCl ocean** using an in-development NaCl(aq) Equation of State from **SeaFreeze**.
- Introduced support for **speciated ocean chemistries** (*CustomSolutions*) using the **Frezchem** and **SUPCRT16** databases, adapted through the chemical modeling package **Reaktoro**.
- Enabled calculation of **chemical (metabolic) reaction affinities** up to **1 GPa** using the **SUPCRT16-organics** database via **Reaktoro**.
- Integrated **PyALMA3** (the Python implementation of **ALMA3**) to compute **tidal Love numbers**.

### Expanded Ocean Chemistry Support
- These updates extend ocean world modeling beyond the previously supported:
  - Seawater
  - MgSO₄
  - Pure H₂O  
- PlanetProfile can now explore a substantially broader and more realistic geochemical parameter space.

### Documentation and Usage
- Guidance on using the **CustomSolution** capability is available upon request.
- A dedicated **PlanetProfile tutorial webpage** is currently in development.

### Bug Fixes
- Fixed numerous incidental bugs across modeling and analysis workflows.