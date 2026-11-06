# Contributions to PlanetProfile — Emma Vellard (JPL, May–November 2026)

This note documents the contributions on three PlanetProfile development branches during my internship in Steve Vance's group: where each change lives, what was ported from `genai`, what was added or hardened, and what is needed before upstream integration. It is written for the PlanetProfile developers.

Status vocabulary follows the repository rule: **verified** (observed under documented steps, artifact cited),
**implemented, unverified** (code written, targeted behaviour not observed for this note), **not implemented**.

Repository state inspected on 4 November and rechecked on 5 November 2026: `origin/main` at `2ae66345` (2 July 2026), `origin/genai` at `60d348a5` (27 October 2026), local lunar branch at `8e3efe07` after the seven November additions. This note incorporates the factual and targeted numerical checks of 5 November; it does not certify every historical campaign or all-body regression.

**How to read the two notes.** Start with [`LUNAR_PROJECT_PLANETPROFILE_USAGE.md`](LUNAR_PROJECT_PLANETPROFILE_USAGE.md):
it says which parts of PlanetProfile the lunar results rest on and how each was checked. This note is the reference for merging the branches and includes work the lunar production study did not use. The cited local November verification scripts and JSON results are not committed in the lunar repository; these two notes are also untracked in PlanetProfile at the time of this check. They need to accompany an upstream review (usage note, access paragraph and Section 6).

### Decisions requested from Steve

1. ~~Push `lunar-structure-development`~~ **Done on 5 November 2026** (`origin/lunar-structure-development`, head `5d4b052c`). The lunar profile now ships with `PPLuna.py` instead of a path on my machine, and `PPLuna` runs on a fresh clone with an empty cache (Section 2.8).
2. **Review the bar-conversion fix in its own pull request.** It changes the conversion applied to raw Perple_X tables and can affect predictions for every body using them. For the lunar reference, it moves the predicted moment of inertia from outside to inside the adopted Williams model envelope (Section 2.1), which is not an independent observational confidence interval. Across the PlanetProfile test bodies, it changes I/MR² for 31 of the 35 that run, by up to 7.0e-4, without changing any validity (Section 2.6).
3. **Replace the shipped core table `Fe-S_3D_EOS.mat`.** Its Fe95S5 slice is broken and its pressure axis is biased. The lunar branch now has a builder for a traceable table (Section 2.7).
4. Decide the fate of **`integrate-genai-to-main` and `genai-clean-port`** (Section 6.1).
5. **Confirm the attribution of ports, new features and hardening in Sections 3 and 4.** The comparison includes the original `genai` code written by Steve; absence of a function name alone does not establish a new capability.
6. **Agree how to name or document the lunar use of `Bulk.Cmeasured`** so that its mean-moment meaning is explicit (Section 2.5), before considering a wider API rename.
7. **Review the seven additions of 5 November** (Section 2.7): which of them should go to `genai` or `main`.

---

## 1. Summary

| Branch                        | Period              | My commits | Pushed              | Base                               | Content                                                                                                                                     |
| ----------------------------- | ------------------- | ---------- | ------------------- | ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `integrate-genai-to-main`     | 8–21 May 2026       | 25         | yes                 | `main` at `d8e08dd0` (4 Feb 2026)  | High-pressure ice rheology: optional Arrhenius viscosity, diagnostic-only Kalousová convection; new Ice VI diagnostic framework |
| `genai-clean-port`            | 26 May–30 June 2026 | 48         | yes                 | `main` at `d8e08dd0`               | Port of the 3-D lateral structure, the TidalPy backend and several fixes from `genai`, with new modes, hardening and tests                  |
| `lunar-structure-development` | 24 Aug–5 Nov 2026   | 17         | yes (5 Nov 2026)    | `genai` at `493205ac` (4 Aug 2026) | Real lunar body: explicit crust and Moho, EOS-domain policy, prescribed silicate temperature, exact bar conversion, regression tests; on 5 Nov, seven additions taken from the lunar study (Section 2.7) |

Counts are commits authored under Emma's Git identity in each stated base-to-tip range, including documentation and the HP-ice branch's merge commit. None of these branch-specific commits is on the inspected `main` revision.

Three points need a decision from the team (Section 6):

1. **`lunar-structure-development` was pushed on 5 November 2026** (new branch on GitHub). It holds the real lunar configuration; before that, it existed only on my machine.
2. **The conversion fix has cross-body scope:** the bar-to-MPa/GPa constants are corrected to the exact SI value (Section 2.1). The cache correction also applies automatically (Section 2.7); neither is behind a feature flag.
3. **Both pushed branches conflict with today's `main`** (9 and 5 conflicting files, Section 6.2).

---

## 2. Lunar structure (`lunar-structure-development`, pushed 5 November 2026)

These fifteen commits sit on top of `genai` at `493205ac`: eight from August (Sections 2.1–2.5) and seven added on 5 November from what the lunar study had built outside PlanetProfile (Section 2.7). Two later commits make `PPLuna` run on a fresh clone (Section 2.8). The framework defaults for explicit crust, prescribed temperature, predictive MoI and density correction preserve the legacy path unless selected; the Luna configuration explicitly enables its lunar options. The table builders and 3-D evaluator are called explicitly. The bar conversion and cached-EOS fixes are automatic corrections, without feature flags.

### 2.1 Exact bar-to-MPa/GPa conversion (`c70c1bf7`): cross-body scope

- `Constants.bar2GPa` and `Constants.bar2MPa` were `1.01325e-4` and `1.01325e-1`, the standard-atmosphere factor, instead of the exact `1e-4` and `1e-1`. The inspected `main` and `genai` revisions still carry the old values (`defineStructs.py`). The [NIST conversion table](https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b9) gives 1 bar = 100000 Pa and 1 standard atmosphere = 101325 Pa.
- Effects of the old values:
  - Raw Perple_X tables, when `PerplexEOSStruct.__init__` loads them: the pressure axis and the moduli `KS_GPa` and `GS_GPa` converted from bar were scaled by +1.325 %. Already converted arrays in the shipped core MAT require separate repair or replacement.
  - PREM pressures, read in the same constructor, were scaled by the same factor.
  - The output-unit multiplier `PmultHydro`, used when `PHYDRO_IN_bar` is set (`SetUnits`), was also affected.
- `MPa2seaPressure` in `SwProps.py` used `bar2MPa` as the atmospheric-pressure offset. It now uses `Constants.P0` (1 atm) explicitly, so the seawater result does not change; only the two quantities are no longer conflated.
- **Size of the effect** on one lunar structure, measured on 5 November 2026:
  - set-up: Warren Preferred mantle, temperature field BU16, Garcia et al. (2019) M3 crust, `nSilMax = 1200`;
  - two closures, identical except for the two constants;
  - result:

| Constants                      | Core radius | Predicted I/MR² | Inside adopted GPM model envelope 0.39311–0.39338? |
| ------------------------------ | ----------- | --------------- | ----------------------- |
| exact (`0.1`, `1e-4`)          | 357.56 km   | 0.393165        | yes                     |
| old (`0.101325`, `1.01325e-4`) | 361.91 km   | 0.393080        | no                      |

The shift of the moment of inertia, 8.6e-5, is about a third of the adopted envelope's width and changes this reference's classification under that criterion. The endpoints span the five GPM models of [Williams et al. (2014), Table 8](https://doi.org/10.1002/2013JE004559); they are not an independent observational confidence interval for the whole-Moon moment. The classifications are conditional on that adopted comparison (usage note, Section 3.1). Script and numbers: `scripts/planetprofile/bar_constant_effect.py` and `results/planetprofile/bar_constant_effect/bar_constant_effect.json` in the lunar repository. Both closures were reproduced in the 5 November audit.

- **Not behind a flag:** every raw-table load using these constants receives the corrected conversion; the resulting change in model predictions depends on the body and inputs. I suggest a separate pull request so that the cross-body effects can be reviewed and tested independently.

### 2.2 Explicit structural crust and Moho handling (`d6b6351d`, tests `787013e2`)

- Optional explicit crust for waterless, non-porous models with an EOS mantle and an iron core; off by default.
- `BuildMohoAwareSilicateGrid` inserts exactly one node at the Moho radius and keeps every legacy boundary bit for bit.
- `GetSilicateMaterialMask` assigns each silicate shell to crust or mantle and raises if a shell straddles the Moho.
- Crust density and thermal conductivity are prescribed constants. The mantle EOS is queried only for mantle shells.
- Structural material IDs (`Constants.material*`) use their own namespace, separate from phase IDs.
- The original crust addition is structure-only: it contributes to mass, pressure, gravity and moment of inertia. By default its seismic properties are set to NaN and flagged `CRUST_SEISMIC_PROPERTIES_UNAVAILABLE`; Section 2.7 adds optional prescribed crust seismic values.
- Tests: `tests/lunar_crust_boundary_test.py` (no Perple_X table or full run needed).

### 2.3 EOS-domain policy (`b8d37102`, `f26debc2`, `e495e3cb`)

- With `EXTRAP=False`, the inner EOS silently replaced out-of-domain (P, T) by the nearest table edge.
- New `Params.EOS_DOMAIN_POLICY`:
  - `'legacy_clamp'` (default) keeps the historical, silent behaviour;
  - `'warn'` gives the same numbers plus a warning with the EOS name, the requested and tabulated ranges, the number of affected nodes and the largest exceedance;
  - `'error'` raises `EOSDomainError` before any out-of-domain evaluation.
- `Planet.Do.EOS_DOMAIN_POLICY` overrides the run-wide value per body (`f26debc2`).
- The policy is not part of `EOSlabel` and is re-applied when a cached EOS is reused.
- Strictness at the physical level (`e495e3cb`): candidate profiles in the core-radius search run down to r = 0 as scratch storage, and those scratch nodes exceed the table ceiling.
  - `EOSCandidatePropagation()` downgrades `'error'` to the clamp during candidate propagation only.
  - `ValidateCandidateEOSDomains()` then applies the declared policy to each candidate's active extent (`[0, nSilFinal)`) and to the core nodes.
  - Crust nodes are excluded from the mantle check and counted separately.
- `EOSIntrospection()` suspends enforcement only around the sentinel probe `fn_phiEOS_frac(0, 0)` in `GetphiCalc`.
- Luna declares `EOS_DOMAIN_POLICY = 'error'`. The historical record for `e495e3cb` reports that enforcing it on physical nodes reproduces the earlier solution at the same resolution. The current Warren integration checks exercised the strict configuration, but the historical before/after strictness comparison was not re-run in this audit.
- Tests: `tests/eos_domain_guard_test.py`, `tests/eos_physical_node_guard_test.py`.

### 2.4 Prescribed silicate temperature (`f26debc2`, `e495e3cb`)

- `Planet.Do.PRESCRIBED_SILICATE_T` with `Planet.Sil.T_of_z_K`, a callable T(depth), so that a body can take its silicate temperature from an external field. Generic, not tied to a body or a dataset.
- The temperature assignment changes in `InitSil` (layer 0) and in the `SilRecursionSolid` loop. The algorithms for pressure, density, gravity, layer mass, moment of inertia and the core-radius search are retained, with the EOS evaluated at the prescribed temperature; their numerical predictions can therefore change. The core keeps its own adiabat.
- A model built this way is hydrostatically consistent given the prescribed field, but not thermally consistent.
- New module `PlanetProfile/Thermodynamics/PrescribedT.py`: linear interpolation only; queries outside the table raise.
- Deep continuation below the data, needed because candidate profiles run to r = 0:
  - `deepContinuation='isothermal'` (opt-in) holds the deepest value and counts every use;
  - `'linear'` (`e495e3cb`) was added to test sensitivity;
  - the default is strict (`None`).
- Tests: `tests/prescribed_silicate_t_test.py`.

### 2.5 A real lunar body (`13015d6b`, guards `17d6ddaa`, resolution `e495e3cb`)

- `PPLuna.py` was Io-derived on the inspected `main`, `genai` and two port branches, retaining the "PPIo" header and Io's bulk radius and mass. It was byte-identical to `PPIo.py` on the two port branches, but not on `main`/`genai`, where `PPIo.py` additionally sets two orbital properties.
- It now defines the Moon, each value labelled as observation, assumption, provisional input or derived value:
  - R = 1737.15 km (LOLA, via Williams et al. 2014);
  - M = 7.3458e22 kg (GM/G);
  - `Cmeasured` = 0.3931 ± 0.0002, the whole-Moon mean moment of Konopliv et al. (1998);
  - a crust of Garcia et al. (2019) model M3, structure-only, with one bulk density of 2751.4 kg m⁻³ that reproduces the layered M3 crust mass;
  - the mantle EOS is a Perple_X table built for the Warren (2005) Preferred composition with the Stixrude and Lithgow-Bertelloni (2024) database;
  - the core is marked `PROVISIONAL_CORE_CONFIGURATION_PENDING_FINAL_SELECTION`; in the lunar production drivers its radius is selected for mass closure on the radial grid and the moment of inertia is a prediction. The legacy `Main.PlanetProfile()` route also applies a moment-of-inertia selection unless the new predictive mode is selected (Section 2.7).
- Resolution: the silicate recursion is first order in the radial step. Refining `nSilMax` from 300 to 600, 1200 and 2400 changes the core radius by 2.895, 1.448 and 0.724 km, so `Steps.nSilMax` is raised from 300 to 1200. At 1200, one radial interval (1.448 km) is the floor on resolvable core-radius differences.
- `tests/luna_config_test.py` checks the sourced values and fails if any Io value or token comes back.
- **Change committed on 5 November (`9ccaac07`)** in `PPLuna.py` (comments only, no code line changed):
  - `Cmeasured` holds the mean moment I = (A+B+C)/3, not the polar moment C;
  - the quoted ±0.0002 is "five times the formal error", not 1σ;
  - the adopted whole-Moon comparator is the model envelope 0.39311–0.39338 of Williams et al. (2014, Table 8), with the model-dependence explained in the usage note, Section 3.1.
  
  Renaming `Cmeasured` is left as a TODO: it appears about 135 times across the code, the app and the tests.

### 2.6 Verification of this branch

- **verified** (5 November 2026, working tree of `lunar-structure-development`, with the comment-only `PPLuna.py` change): the 12 test files in `tests/` listed below → **122 passed** (environment `planetprofile`):  `lunar_crust_boundary`, `eos_domain_guard`, `eos_physical_node_guard`, `luna_config`, `prescribed_silicate_t`,  `perplex_eos_cache`, `fe_s_traceable_core`, `layered_mantle_table`, `angular_mean_density`, `predict_moi`,  `silicate_3d` and `crust_seismic` (`*_test.py`).
- **verified for the Warren BU16/M3 reference:** its closure on the November revision reproduces the recorded production core radius (357.563375 km) and I/MR² (0.39316529266712535). This is a regression check on that case, not a guarantee for every body or every cached-EOS call.
- **Bar conversion across the test bodies: verified** (5 November 2026).
  - Method: every body of `PlanetProfile/Test/PPTest*.py` (58) was run twice through `Main.PlanetProfile()`, with
    the exact constants and with the old ones, using the `BuildTest` settings without plots, files or induction.
  - Results, for the 35 bodies that run:
    - none changes validity;
    - I/MR² changes for 31 of them, by up to 7.0e-4 (PPTest42, 45, 46);
    - the core radius changes for 14, by up to 36.6 km (PPTest9).
  - The other 23 bodies fail identically in both modes and at the branch point `493205ac` of `genai`, so their
    failures are pre-existing:
    - 12 need porosity or a core for the moment-of-inertia match;
    - 5 need reaktoro;
    - 3 call `GetTfreeze` with an argument it does not take;
    - 1 needs numpy ≥ 2;
    - 1 has `Tb_K` too high;
    - 1 fails with a `NoneType` error.
  - Numbers, scripts and logs: `results/planetprofile/branch_verification_2026-11-05/` in the lunar repository.
- The full `BuildTest` (inductograms and exploreograms included) was not run.

### 2.7 Additions of 5 November 2026, from the lunar study

The new modelling options are opt-in; the builders and 3-D evaluator are invoked explicitly. The cache fix is automatic. Each addition has targeted tests, with the scope of the comparison stated below. The local integration scripts are in the lunar repository, `scripts/planetprofile/`; the four scripts for density correction, predictive MoI, 3-D fields and crust seismic values were re-run in the audit. They and their JSON outputs are not yet committed.

| Commit | Addition | Check against the study and verification scope |
| --- | --- | --- |
| `1b227ca0` | `PerplexEOSStruct.__init__` keeps the loaded state when its table is already cached. A second direct construction used to return an object without domain limits or property functions | **verified:** 4 tests pass; the commit records two failures on the previous code, which was not re-tested in this audit |
| `93f7f1cb` | `Thermodynamics/Perple_X/stackFeSTraceable.py` builds an optional traceable Fe–S core table from the five shipped werami tables with exact unit conversions; the distributed MAT is left unchanged | **verified:** every physical array equals the study product; semantic SHA-256 `452147a82062916219dc8b83e53dcb53114ccfe8f60030b8afb4a9db72834445` under the documented array-shape convention |
| `6cf4cf32` | `Thermodynamics/LayeredMantle.BuildLayeredPerplexTable`: a mantle of two compositions, as one table split at an interface pressure | **verified** for synthetic splits and byte equality when reusing the final L500 table as its upper input. Independent reconstruction from the original U500 + M0 inputs is **implemented, unverified**; U500 source builds were unavailable in the expected directory |
| `bfed82be` | `Thermodynamics/AngularMeanDensity`: structure built on the angular-mean density of a 3-D mantle (`Planet.Sil.densityCorrection`, offset calculator) | **verified:** identical core radius, I/MR² and profiles to the external wrapper for a synthetic density correction on Warren BU16 (`verify_native_density_correction.py`) |
| `6a08e9f5` | `Planet.Do.PREDICT_MOI`: for a waterless iron-core body, C/MR² is reported as a prediction and compared with `Bulk.CMR2comparison`, without invalidating the model | **verified** on Warren BU16 through `Main.PlanetProfile()`: identical I/MR² and core radius to the study's shell sum (`verify_predict_moi.py`) |
| `7a7e6a3c` | `Thermodynamics/Silicate3D.EvaluateSilicate3D`: rho, VP, VS, KS, GS of the mantle on a 3-D temperature field around a radial structure | **verified** on `W_LOW_HET` and `ALT5_B1`: radial references identical; maximum relative 3-D property difference below 5.8e−7 against the stored float32 cubes (`verify_silicate_3d.py`) |
| `8e3efe07` | `Seismic.AssignCrustSeismic` with optional `Crust.VP_kms`, `VS_kms`, `QS`: prescribed crust velocities and moduli instead of NaN | **verified** on Warren BU16 through `Main.PlanetProfile()`: prescribed VP/VS set, mantle seismic values and structure unchanged (`verify_crust_seismic.py`) |

Not added: the study's hybrid interpolation (bilinear within one mineral assemblage, nearest across a phase
boundary), left for later. The export of a PlanetProfile structure to LOV3D already exists on the LOV3D side
(`planetstruct_to_interior_model` in `LOV3D_multi/pylov3d/compat.py`, branch `python-conversion`).

### 2.8 Running `PPLuna` on a fresh clone (5 November 2026)

Before these commits, `PPLuna` ran only on my machine: `PerplexEOSStruct` reads tables only from the user
Perple_X cache, and neither lunar table was in the repository or downloaded by the installer.

| Commit | Change | Verification |
| --- | --- | --- |
| `c505de49` | Ships `moon_w05_preferred_stx24_unified.tab` (12,150,425 bytes, sha256 `eb254ac5d454cd59f5ca1c29220214eb8d036cc21b44251b6db558c9df2d8348`) in `PlanetProfile/Thermodynamics/EOStables/Perple_X/`. `InnerEOS.EnsurePerplexTableInCache` copies a shipped table to the cache on first use, and builds `Fe-S_3D_EOS_5slice_traceable_v1.mat` with `stackFeSTraceable.py`, checked against its semantic SHA-256. A table already in the cache is never replaced | **verified:** `python PlanetProfile/Main.py Luna` on a clean checkout with an empty cache (temporary `HOME`) copies the mantle table, builds the core table, and gives a profile identical byte for byte to the run with my existing cache (C/MR² = 0.3931, R_core = 0.21 R). `tests/perplex_table_provisioning_test.py`: 4 tests pass |
| `5d4b052c` | `EOSDomainExceedance` now counts the P × T grid when a P axis and a T axis of different lengths are passed (`grid=True`), instead of raising a broadcast error. `PlotPvTPerpleX` leaves blank the plot nodes outside a table under the `'warn'` and `'error'` policies; `'legacy_clamp'` is unchanged | **verified:** the run used to stop after saving the profile, because the core plot grid reaches 475 K, below the 500 K edge of the core table; it now completes, and the `InnerPTprops` figure was inspected. A grid test was added to `tests/eos_domain_guard_test.py` |

Full `tests/` folder, before and after these two commits: the same 23 failures and 13 errors, which already exist at
`7d569fc2`; 5 more tests pass (161 → 166). Only the Warren Preferred mantle table is shipped: the other lunar
tables (Section 2.2 of the usage note) are still provisioned by the lunar study.

---

## 3. 3-D lateral structure and TidalPy (`genai-clean-port`, pushed)

Goal of the branch: bring the 3-D lateral structure of `genai` (Steve, `272499b0` and `fc9a5c93`, 21–22 April 2026) back to `main` as small, reviewable commits with tests. My commits prefixed `[port]` carry code from `genai`, but some also add capability or harden existing behaviour. The tables distinguish those roles using the inspected source behaviour, not only function names. In particular, column parallelism was already present in April; the later refactor is credited as hardening.

### 3.1 Ported from `genai`

| Feature | Commit(s) | Notes |
| --- | --- | --- |
| Lateral module: `SpatialGrid`, `TidalHeating3D`, `LateralStructure`, `MassConservation`, `ClathrateLateral`, `LateralIO`, `LateralPlots` | Phases 1–10 (`b9ef6346` … `dfc3d8bd`) | Split into ten phases; HEALPix grid (optional `healpy`) or lat-lon fallback; 4π-normalised spherical harmonics as in MoonMag; `Main.py` integration behind `Planet.Lateral.DO_3D` (default off); documentation in `docs/LATERAL_3D_USAGE.md`. The original `genai` module already had `_RunColumnsParallel` with Unix `fork` / Windows `spawn`; its later rework is listed in Section 3.2 |
| TidalPy backend | `0efa263a` | `Params.Gravity.backend = 'tidalpy'` or `'pyalma'` (default); per-phase heating; Andrade, Maxwell, elastic and Newton rheologies. Explicitly selecting TidalPy when it is unavailable raises `ImportError`; there is no automatic switch to PyALMA3 on that path |
| CSV export of profiles | `d57bc365` (from `73ac449f`) | `HumanReadableOutput.WriteProfileCSV()`, 23 columns; adds `pandas>=2.0.0`; `pMax = 0` added to three inductogram tests |
| Clathrate depth fix | `34ce46ea` | `ConductiveTemperatureCorrect()`, Turcotte & Schubert (2002) eq. 4.40; the old factor 2 in c1 put clathrate layers at twice the requested depth |
| Silicate boundary condition for solid spheres | `beb53112` (from `80853f4f`) | c1 = 0 so T stays finite at r = 0 when there is no iron core; shell bodies keep the legacy call |
| Triple-point constant renamed | `d09b2641` | `triplePointT_K` → `TtripleIh_III_L_K`, old name kept as alias; HP-ice triple points added |

### 3.2 Additions and hardening on `genai-clean-port`

| Feature | Commit(s) | Notes |
| --- | --- | --- |
| Equilibrium ice-thickness mode | `14b27215` | `CalcEquilibriumIceThickness()` and `DO_EQUILIBRIUM_ICE`; mode hierarchy equilibrium → prescribed spherical harmonics → function → uniform; `_ValidateEquilibriumMode()`; documented in `docs/3D_ICE_THICKNESS_MODES.md` |
| Obliquity tides | `61ac88b1`, `61c01f8c` | `TidalStrainPattern()` includes obliquity forcing; the obliquity parameter is connected to the heating calculation |
| Column orchestration and parallelism hardening | `61ac88b1` | Refactors the existing parallel-column capability with `_PrepareColumn()` and `_RunSingleColumn()` and fixes NaNs from a missing `Params.nModels`. Parallel execution and the Unix/Windows process contexts already existed in `272499b0`; they are not new capabilities of this commit |
| NPZ I/O | `61ac88b1` | `np.savez_compressed` replaces pickle; old `.pkl` files still load with a warning |
| Tobie et al. (2005) regression test | `61ac88b1` | `PPTestTobie2005Regression.py`: Titan heating against 3e-8 W m⁻³, 50 % tolerance |
| Heating-profile transfer and export | `ee1fd747`, `2ebd3c67` | `TransferTidalPyHeatingToProfile()` maps the TidalPy heating onto `Planet.Htidal_Wm3` for profile/CSV output. In the main workflow this follows the thermal calculation and does not itself re-solve or converge it; a converged thermal feedback calculation needs a separate iteration |
| Numerical guards at the poles | `bde1c067` | Epsilon in `TidalStrainPattern`, `_MaxwellDissipation`, `_AndradeDissipation` |
| Failed-column repair | `dc24fc28` | Failed 3-D columns are treated as invalid data, not zero heating; isolated failures repaired by nearest-neighbour interpolation; QA metadata saved; equilibrium solves guarded against non-finite heating |
| Love-number normalisation for 3-D heating | `25b98950` | With TidalPy support |
| NumPy 2 and LaTeX compatibility | `679d362f`, `f2cc373e`, `21b7bb40` | `np.int_` → `np.intp`, `np.float_` → `np.float64`, `np.trapz` → `np.trapezoid`; a `ymin`/`ymax` typo in `PTPlots.py` |
| Tests | Phases 1–9, `PPTest3DIceThicknessModes.py`, `test_lateral_failure_repair.py` | The CHANGELOG reports 59 tests for the lateral port |

### 3.3 Verification of this branch

- **verified in part** (5 November 2026, test environment below): of the 72 tests of the 14 test files of the port,
  **65 pass and 7 fail**. The 7 failures are defects of the tests themselves; none points to the module code:
  - `PPTestLateralPhase4` imports `_RunColumnsSerial`, removed when the column dispatch was refactored (`61ac88b1`);
  - three `PPTestLateralPhase7` tests still expect a `.pkl` file, but the module now writes NPZ (`61ac88b1`);
  - `PPTestTobie2005Regression` builds a mock planet without `T_K`, which `ComputeTidalHeating3D` needs;
  - two `PPTestEuropaStructure` tests call `GetConfig`, a module, as a function.
- **Test fixes, local commit `97b229cd`** on `genai-clean-port` (pushed 5 November 2026): the four test files are brought up to
  date with the module. Result: **71 of 72 tests pass**.
  - The Europa tests now run the default Europa body and pass: ice shell 13.4 km, H2O layer 101.6 km.
  - The Tobie 2005 test now runs and **fails on physics**: the mean ice heating of the simplified path is
    6.2e-10 W m⁻³, against the 3e-8 W m⁻³ reference (50 % tolerance), about 50 times too low. The inputs were not
    tuned; this is an open point for Steve.
  - Seen during these runs: the TidalPy heating transfer rewrites the profile CSV even with `Params.NO_SAVEFILE`.
- Test environment (temporary; the `planetprofile` environment itself is unchanged):
  - a detached worktree of the branch;
  - `healpy` 1.16.6, with a `scipy.integrate.trapz` shim because the environment has numpy 1.26 and scipy 1.16;
  - the MoonMag copy vendored on the lunar branch, with version metadata 1.7.5;
  - an import-only `reaktoro` stub that raises on any use;
  - PlanetProfile version metadata;
  - a temporary Perple_X cache linked to the repository tables.
  
  Logs and harness files: `results/planetprofile/branch_verification_2026-11-05/` in the lunar repository.
- Source inspection and an isolated execution of the missing-TidalPy branch confirm the explicit `ImportError` described above. That check does not validate the TidalPy solver or a converged thermal coupling. The transfer call order and output behaviour were inspected in `Gravity.py` and `Main.py`.
- `beb53112` records that its tests could not be run at the time (spiceypy installation failed).
- The branch also contains working files at the repository root that should not go to `main`: session handoffs, status and summary files, example scripts, `CLAUDE.md` and `.claude/`.

---

## 4. High-pressure ice rheology (`integrate-genai-to-main`, pushed)

Goal: bring the HP-ice rheology of `genai` (`272499b0`) to `main` without changing default results.

### 4.1 Ported from `genai`

- **Arrhenius viscosity:** eta(T) = eta_melt exp[(Eact/R)(1/T − 1/Tmelt)]. Here it is opt-in and per phase, with
  the flags `ARRHENIUS_VISCOSITY`, `_Ih`, `_III`, `_V` and `_VI`, all `False` by default (`145ddec0`, `20eff971`; `docs/rheology/arrhenius_viscosity.md`).
- **Kalousová and Sotin (2018) scalings:** kept under the `genai` flag name `KALOUSOVA_CONVECTION`, but
  **diagnostic-only** here. They do not change T, phase, mass, gravity, heat flux or viscosity profiles (`1b6416a4`).

### 4.2 Added on `integrate-genai-to-main`

These flags and records are absent from `genai`:
- **Model selector:** `HP_ICE_CONVECTION_MODEL` with the values `"none"` (default), `"DS2001_diagnostic"`,
  `"Kalousova2018_diagnostic"` and `"Kalousova2018_production_experimental"`. The legacy flags still work as shims  (`3a4d5113`).
- **Production gate:** `ALLOW_EXPERIMENTAL_HP_KALOUSOVA_PRODUCTION`. It enables a dry run that records candidate statuses and residuals without applying any profile update (`a3fac645`).
- **Phase-local state:** a `Planet.HPIceDiagnostics` record for Ice III, V and VI, with geometry, heat-through-stack bookkeeping, residuals and validity status. `_SetHPIceDiagnosticFields()` is the single writer (`1aee2889`, `af53bf59`, `c5af50bf`).
- **Ice VI candidate checks:** melt-curve checks, acceptance criteria, residual evaluator, margins and risk,
  rollback scaffolding, and a thermal reconstruction on the original boundary (`3a2ba639` … `f4e46daf`).
- **Documentation and tests:** `docs/rheology/hp_ice_convection.md` and `PlanetProfile/Test/test_rheology_features.py`.

### 4.3 Verification of this branch

- **verified** (5 November 2026, same temporary test environment as Section 3.3, with `sqlalchemy` added):
  `PlanetProfile/Test/test_rheology_features.py` → **151 passed**, no failure.
- By design, nothing here activates production feedback. Defaults leave existing results unchanged; this was the contract the branch was written to, but it was not re-checked for this note.

---

## 5. Smaller changes

- CI and documentation workflow: artifact actions v3 → v4 and a newer Python for the docs build (`d805ef51`, `ece745e4`, `99edb509` on `integrate-genai-to-main`; `eb3e3200` on `genai-clean-port`).
- Missing `SPEC_FILE` default (`8b313541`, `integrate-genai-to-main`).
- `GetTfreeze` handles array returns from `fn_phase` (`f9f84f81`, `genai-clean-port`).

---

## 6. Status and what is needed to reach `main`

### 6.1 Coordination questions in `docs_ai/DEVELOPMENT_PATHWAY.md` (on `genai`)

- **(a) Is `PlanetProfileApp/` on `genai-clean-port`?** No. It is on neither `genai-clean-port` nor
  `integrate-genai-to-main` (checked with `git ls-tree`), so there is no duplicate app directory to reconcile. Integration with Mara's port can still touch shared dependencies and needs its own merge check.
- **How the two pushed branches relate.** Both start from `main` at `d8e08dd0`. The HP-ice branch ports work from `272499b0`; the lateral/TidalPy branch ports work from `272499b0` and `fc9a5c93`. Their main feature sets are distinct, although some files overlap:
  - they share the base history, but no branch-specific commits beyond `d8e08dd0`;
  - the high-pressure ice rheology (Arrhenius, Kalousová, selector, Ice VI diagnostics) exists only on `integrate-genai-to-main`;
  - the lateral module, TidalPy and the other fixes exist only on `genai-clean-port`;
  - they touch eight common files and conflict with each other in three (`sphinx.yml`, `.gitignore`, `defineStructs.py`).
  
  Dropping either branch would therefore drop its feature set.
- **(b) Is `integrate-genai-to-main` abandoned or still needed?**
- **(c) Does the HP-ice work there supersede the corresponding `genai` commits?** Factually, it keeps the `genai` physics but makes it opt-in and diagnostic-only, and adds the selector, the phase-local state and the Ice VI checks.

### 6.2 Merge check (simulated with `git merge-tree --write-tree`, rechecked 5 November 2026)

| Merge | Conflicting files |
| --- | --- |
| `genai-clean-port` into `origin/main` | 9: `Gravity/Gravity.py`, `Main.py`, `Clathrates/ClathrateProps.py`, `Geophysical.py`, `HydroEOS.py`, `IronCore.py`, `LayerPropagators.py`, `ThermalProfiles/Convection.py`, `Utilities/defineStructs.py` |
| `integrate-genai-to-main` into `origin/main` | 5: `HydroEOS.py`, `LayerPropagators.py`, `ThermalProfiles/IceConduction.py`, `Utilities/SetupInit.py`, `Utilities/defineStructs.py` |
| `lunar-structure-development` into `origin/genai` | 1: `Utilities/defineStructs.py` (`genai` has moved 312 commits since the branch point; still one file after the additions of Section 2.7, checked on 5 November) |

### 6.3 Proposed next steps (proposals, for the team to decide)

The decisions themselves are listed at the top of this note. Once they are taken, the practical steps are:

1. ~~Commit the comment-only `PPLuna.py` change, then push `lunar-structure-development`~~ Done on 5 November (`9ccaac07`, then `d6a34a24` for the shipped profile, and `c505de49`, `5d4b052c` so that `PPLuna` runs on a fresh clone).
2. **Recover the 14 targeted bar-unit tests of August** (`pp_bar_unit_test.py`, not in the current tree). The
   before/after comparison on the test bodies is done (Section 2.6).
3. **Rebase the lunar commits on current `genai`**: one conflict, in `defineStructs.py`.
4. **Clean `genai-clean-port` before a pull request to `main`:** remove the session and handoff files, fix the 7
   out-of-date tests (fixed in `97b229cd`, pushed) and look into the Tobie 2005 result
   (Section 3.3), then resolve the 9 conflicts.
5. **Merge the two pushed branches one after the other.** Their pairwise merge currently conflicts in 3 files; recompute conflicts after the first integration rather than treating the separate counts as additive.
6. **Include the review evidence:** version these two notes and the intended local verification scripts/results, or provide them with the review. Recover the U500 source builds and perform an independent U500 + M0 layered-table reconstruction before claiming that source-to-product reproduction is verified. The usage note, Section 6, separates that input preparation from core-table construction and records the remaining prerequisites.

---

## 7. Where to look

- Lunar branch: `git log 493205ac..lunar-structure-development`; tests in `tests/` (`*_test.py`); module docstrings of `AngularMeanDensity.py`, `LayeredMantle.py`, `Silicate3D.py` and `Thermodynamics/Perple_X/stackFeSTraceable.py`.
- Lateral port: `docs/LATERAL_3D_USAGE.md`, `docs/3D_ICE_THICKNESS_MODES.md`, `CHANGELOG.md` [Unreleased] on `genai-clean-port`.
- HP ice: `docs/rheology/arrhenius_viscosity.md`, `docs/rheology/hp_ice_convection.md` on `integrate-genai-to-main`.
- How the lunar structures were used: the repository `lunar-bsm-interior-models` (not part of PlanetProfile).
