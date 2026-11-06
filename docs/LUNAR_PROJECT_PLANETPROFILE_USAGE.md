# How the lunar interior project used PlanetProfile — Emma Vellard (JPL, 2026)

This note describes what the lunar project (repository `lunar-bsm-interior-models`) actually used in PlanetProfile. It covers the version, the functions and options, the inputs we supplied, what we computed outside PlanetProfile, and how each piece was validated. It also lists the problems we found in PlanetProfile along the way.

A companion note, [`CONTRIBUTIONS_EMMA_VELLARD_2026.md`](CONTRIBUTIONS_EMMA_VELLARD_2026.md), documents the contributions on the three development branches. The work on `genai-clean-port` and `integrate-genai-to-main` played **no part** in the lunar production results.

Status vocabulary: **verified** (observed, artifact cited), **implemented, unverified**, **not implemented**.

**Reading order.** Read this note first: it is what the lunar results rest on. The contributions note is the reference for merging the branches, and it opens with the list of decisions requested from the team.

**Access to the cited files (checked 5 November 2026).** Paths starting with `lunar:` refer to the private project repository [`EmmaVellard/lunar-bsm-interior-models`](https://github.com/EmmaVellard/lunar-bsm-interior-models) or its local working tree. The historical audit documents are tracked, but the October audit snapshot `results/audits/park2025_k3m_20261020/initial_git_state.json`, the November scripts `scripts/planetprofile/verify_*.py` and `bar_constant_effect.py`, their result JSON files, and the internship report appendix are **local, uncommitted artifacts**. These two PlanetProfile notes are also untracked at the time of this check. Sharing the committed repositories alone therefore does not supply every cited artifact. Each row of Section 4 states its result; the local evidence must accompany a review that needs to reproduce it.

---

## 1. Version used

| Work | PlanetProfile version |
| --- | --- |
| Production, September–October 2026: structures rebuilt on each published temperature field, selection under the adopted admissibility criteria, three-dimensional property fields, the 19 retained interiors and their tidal backgrounds | branch `lunar-structure-development` at `e495e3cb`, PlanetProfile 3.1.0; the 20 October snapshot records an uncommitted, comment-only `PPLuna.py` change |
| First diagnostics, August 2026: first lunar structure, audit of the thermal formulation, bar-unit audit, Fe–S core audit | recorded HEAD `493205ac`, the `genai` commit on which the lunar branch was started; some diagnostic changes were still uncommitted, so the individual run records are needed as well as the SHA |

The commit is recorded in the project, for example in:
- `lunar:results/audits/park2025_k3m_20261020/initial_git_state.json`;
- `lunar:docs/STAGE_D_STRUCTURAL_ADMISSIBILITY.md`;
- `lunar:docs/STAGE_C_THERMAL_UNCERTAINTY_AND_COMPOSITION_INTERACTION.md`.

The production drivers use external wrappers rather than editing PlanetProfile source files (Section 3.2). The recorded working tree was not strictly clean: `initial_git_state.json` lists the comment-only `PPLuna.py` difference. No executable change appears in that recorded diff.

Two points to know:
- **The branch was pushed on 5 November 2026** (`origin/lunar-structure-development`). The production results were computed at `e495e3cb`, which is part of it.
- **Two later commits touch `PPLuna.py` without changing any result:** `9ccaac07` (comments on the meaning and uncertainty of `Cmeasured`) and `d6a34a24` (the temperature profile now ships in `PlanetProfile/Default/Luna/` instead of being read from a local path).

---

## 2. What PlanetProfile computes for us

The project uses PlanetProfile for three things:
1. loading the Perple_X tables;
2. building a spherically symmetric, hydrostatic, mass-closed crust–mantle–core structure;
3. evaluating density and seismic properties at given (P, T).

### 2.1 Body configuration (`PlanetProfile/Default/Luna/PPLuna.py`)

| Quantity | Value | Source or role |
| --- | --- | --- |
| `Bulk.R_m` | 1737.15 km | LOLA mean radius, via Williams et al. (2014) |
| `Bulk.M_kg` | 7.3458e22 kg | GM / G |
| `Bulk.Cmeasured`, `Cuncertainty` | 0.3931, 0.0002 | Konopliv et al. (1998), whole-Moon mean moment. Not used for selection in production (Section 3.1) |
| `Bulk.Tsurf_K`, `Psurf_MPa` | 250 K, 0 | Assumptions; the mantle temperature is prescribed (below) |
| `Do.Fe_CORE`, `Do.POROUS_ROCK` | `True`, `False` | Waterless body with a metallic core, no porosity model |
| `Do.PRESCRIBED_SILICATE_T` | `True` | Mantle temperature from an external field |
| `Do.EOS_DOMAIN_POLICY` | `'error'` | No out-of-domain EOS value at a physical node |
| `Steps.nSilMax`, `Steps.nCore` | 1200, 10 | Silicate resolution 1200 selected by the convergence test in Section 4; core resolution 10 retained from the configuration, without a separate core-grid convergence test cited here |
| `Sil.Qrad_Wkg` | 0 | The temperature is prescribed, so no conductive heating is needed |
| `Crust.*` | enabled, structure-only; thickness and bulk density set per crust | Explicit crust (Section 2.3) |
| `Core.coreEOS`, `Core.wFe_ppt` | `Fe-S_3D_EOS_5slice_traceable_v1.mat`, 950 | Fe–5 wt% S core (Section 2.4) |

### 2.2 Mantle equation of state

- **Loader:** `InnerEOS.PerplexEOSStruct` reads each Perple_X table, which supplies density, K_S, G_S, V_P and V_S.
- **Tables:** built in the project (planetary-eos-lab, Stixrude and Lithgow-Bertelloni 2024 database, "stx24") and copied into PlanetProfile's cache `~/Library/Caches/PlanetProfile/Perple_X/`. The N2 installer `lunar:scripts/dps_2026_followup/fu_n2_a_install_tables.py` handles the FeO-rich tables and assembles the layered tables; it assumes the earlier M0/M1 tables are already available (Section 6). There are three families:
	  - homogeneous residual mantles: `moon_{hau05_lowFe, w05_preferred, war05_alt5}_M1<crust>_stx24_unified.tab`;
	  - layered mantles (`..._L500<crust>_...`): stitched in the project from an upper-mantle residual table above a fixed pressure and the bulk-composition table below it;
	  - FeO-rich mantles: `moon_{Khan06S, ONeill91S}_M1<crust>_stx24_unified.tab`.
- **Interpolation:** `EOSinterpMethod='nearest'`, the PlanetProfile default (`Params.lookupInterpMethod`), and `EOSdomainPolicy='error'`.

### 2.3 Explicit crust

- The crust is structure-only (`Crust.enabled = True`, `structure_only = True`). It contributes mass, pressure, gravity and moment of inertia, but no seismic property.
- Its thickness and bulk density are set for each crust tested: Garcia et al. (2019) model M3, 40 km, 2751.4 kg m⁻³, the GRAIL crust of Wieczorek et al. (2013), 34 km, 2550 kg m⁻³; and the documented alternatives.
- The Moho is a single structural node, and no shell mixes crust and mantle.

### 2.4 Core equation of state

- The core table shipped with PlanetProfile, `Fe-S_3D_EOS.mat`, failed our audit (Section 5), so it was not used.
- The project built `Fe-S_3D_EOS_5slice_traceable_v1.mat` instead (`lunar:scripts/build_fe_s_5slice_core_eos.py`). It stacks, without resampling, the five werami tables shipped in the PlanetProfile repository (`80Fe20S` … `100Fe0S_core.tab`).
- The Fe95S5 slice (`wFe_ppt = 950`) is used: 5 wt% S, inside the 0–6 wt% range of the literature.
- Product record: `lunar:docs/FE_S_CORE_EOS_PRODUCT.md`.

### 2.5 Prescribed mantle temperature

- `PrescribedT.LoadPrescribedTProfileCSV()` reads T(depth): the area-weighted mean of each three-dimensional field of Santangelo et al. (2025).
- `deepContinuation='isothermal'` holds the deepest temperature below the data. This is needed because candidate profiles are propagated down to r = 0 (Section 3.1).

### 2.6 Radial closure

- **Calls:** the production drivers (`lunar:scripts/stage_d/d_lib.py`, `lunar:scripts/stage_e9a5/e9a5_lib.py`, `lunar:scripts/stage_e9a6/e9a6_lib.py`) call PlanetProfile's layer functions directly, in this order:
	  1. `SetupInit`;
	  2. `GetPlanetInnerEOS` (mantle and core);
	  3. `LayerPropagators.GetHtidalFunc` and `GetphiCalc`;
	  4. `Silicates.SilicateLayers` and `IronCore.IronCoreLayers`, inside the `InnerEOS.EOSCandidatePropagation()` context.
- **Outputs:** the core radius selected for lunar mass closure on the radial grid, and the profiles of P, T, r, ρ, g and mass of the active nodes. Closure is numerical, not exact: the Warren BU16 reference has a mass residual of −0.006743 %, with a 1.448 km silicate radial interval.
- **Run options:** `CALC_NEW = True`, `SKIP_PLOTS`, `NO_SAVEFILE`, `SKIP_INDUCTION`, `SKIP_GRAVITY`,
  `CALC_CONDUCT = False`.

### 2.7 Seismic properties of the three-dimensional fields

- `PerplexEOSStruct.fn_Seismic(P, T)` returns V_P, V_S, K_S and G_S at every point of the 3-D temperature field.
- The pressure is that of the closed radial structure, so these 3-D fields are a P/T decomposition around it, not a 3-D hydrostatic solution.
- They feed the seismic comparisons and the lateral shear-modulus coefficients passed to the tidal solver.

### 2.8 Not used

The following parts of PlanetProfile did not enter the lunar production results described here:
- the gravity and Love-number modules (PyALMA3, TidalPy), which were skipped (the tides were computed with LOV3D);
- magnetic induction and electrical conductivity;
- the porosity model and all hydrosphere code;
- PlanetProfile's own selection on the moment of inertia (Section 3.1);
- the `Lateral` module, which computes the hydrosphere only and keeps silicates and core at the reference;
- the HP-ice rheology;
- plotting and output files.

The first diagnostics of August and the first comparison of three compositions (`lunar:scripts/run_r1_*.py`,
`run_r2b_path_case.py`, `scripts/planetprofile/three_bsm_run.py`) ran the full `Main.PlanetProfile()`, but no
production result depends on them.

---

## 3. What the project computes outside PlanetProfile

The lunar results were produced with the external versions below, at `e495e3cb`. Since 5 November 2026, corresponding implementations also exist on the lunar branch. Predictive MoI, density correction and prescribed crust seismic properties are opt-in; the 3-D evaluator and table builders are called explicitly. The cache correction is **automatic**, with no off-by-default flag. The scope of each comparison, including the limited layered-table check, is stated in the contributions note, Section 2.7.

| Done outside PlanetProfile for the results | Now in PlanetProfile (lunar branch) |
| --- | --- |
| Moment of inertia as a shell sum, never fitted (3.1) | `Planet.Do.PREDICT_MOI`, `Planet.Bulk.CMR2comparison` |
| Angular-mean density wrapper (3.2) | `Planet.Sil.densityCorrection`, `Thermodynamics/AngularMeanDensity.py` |
| 3-D mantle properties with `fn_Seismic` (2.7) | `Thermodynamics/Silicate3D.EvaluateSilicate3D` |
| Layered mantle tables stitched by a script (2.2) | `Thermodynamics/LayeredMantle.BuildLayeredPerplexTable` |
| Traceable Fe–S core table built in the project (2.4) | `Thermodynamics/Perple_X/stackFeSTraceable.py` |
| Crust moduli from published velocities, in the tidal chain | `Planet.Crust.VP_kms`, `VS_kms`, `QS` |
| One memoised EOS object per table (5, item on cached tables) | fixed in `PerplexEOSStruct.__init__` |

The export of a PlanetProfile structure to LOV3D is on the LOV3D side: `planetstruct_to_interior_model` in
`LOV3D_multi/pylov3d/compat.py` (branch `python-conversion`).

### 3.1 Moment of inertia

- The moment of inertia is computed by the project from the closed structure, as a shell sum
  I = Σ (8π/15) ρ (r_top⁵ − r_bot⁵) over all active silicate shells (crust included) and core shells. It is a **prediction**. The project compares it with the whole-Moon model envelope 0.39311–0.39338 in [Williams et al. (2014), Table 8](https://doi.org/10.1002/2013JE004559), and separately compares core radius with Weber et al. (2011).
- **Meaning of the Williams comparator:** those endpoints span five GPM interior models that satisfy the adopted density, solid-Moon moment and Love number. They are not the bounds of an independent observational confidence interval for the whole-Moon moment. Williams et al. report the solid-Moon mean moment as 0.393112 ± 0.000012 at R = 1737.151 km; converting that quantity to a whole-Moon moment requires a fluid-core model. The project's classifications are conditional on its adopted model envelope and other declared criteria. Here, "admissible" means passing those criteria, not an independent validation of every model ingredient.
- PlanetProfile's selection by `Cmeasured` (`FindInnerWithMoIAndEOS`) is not used, so the moment of inertia is never fitted.
- For a waterless body the first valid candidate (`indsSilValid[0]`) is the structure; its active extent is
  `[0, nSilFinal)`.

### 3.2 Density of a three-dimensional mantle

- The later angular-mean reclosures (Stages E9a.5/E9a.6 and their follow-ups) build structures on the angular mean of the 3-D density. Earlier closures such as Stage D used the density at the prescribed mean temperature. The two differ where density is not linear in temperature.
- This is done with an external wrapper, `AngularMeanDensityEOS` (`lunar:scripts/stage_e9a5/e9a5_lib.py`). It overrides only `fn_rho_kgm3`, adds a prescribed correction Δρ(P), and is clamped to zero above the Moho. Every other accessor of the PlanetProfile EOS is delegated unchanged.
- With Δρ = 0, it reproduces the unwrapped closure.

### 3.3 Everything downstream

- The lateral shear-modulus fields, their spherical-harmonic coefficients and the degree-3 tidal response are computed outside PlanetProfile. The tidal response uses pyLOV3D on a radial background reduced from the PlanetProfile structure.
- The seismic comparisons with the Apollo models are also outside PlanetProfile.

---

## 4. Validation of what we used

| What | Result | Status and artifact |
| --- | --- | --- |
| Lunar regression tests: crust and Moho, EOS-domain guard, physical-node guard, Luna configuration, prescribed T | 83 passed | **verified** 4 Nov 2026, `pytest tests/{lunar_crust_boundary,eos_domain_guard,eos_physical_node_guard,luna_config,prescribed_silicate_t}_test.py`, environment `planetprofile` |
| Bar-unit correction (`c70c1bf7`) | The historical audit reports 14 targeted tests passed (10 Aug 2026), and 32 discovery errors classified there as environment or incompatible-cache failures. This is a historical record, not a current full-suite pass | **verified** at the time according to `lunar:docs/PLANETPROFILE_BAR_UNIT_AUDIT.md`. The test file `pp_bar_unit_test.py` cited there is not in today's tree: it must be recovered before upstream review |
| Bar-unit correction across the PlanetProfile test bodies | 35 of 58 bodies run with the exact and the old constants: no change of validity; I/MR² changes by up to 7.0e-4 and the core radius by up to 36.6 km. The 23 others fail identically in both modes and at `493205ac` (pre-existing) | **verified** 5 Nov 2026, `lunar:results/planetprofile/branch_verification_2026-11-05/` |
| Core EOS product | Reproduces every source cell after unit conversion; `PerplexEOSStruct` loads it exactly | **verified** 10 Aug 2026, `lunar:docs/FE_S_CORE_EOS_PRODUCT.md` (SHA-256 recorded) |
| Radial resolution | First order in the step: `nSilMax` 300 → 600 → 1200 → 2400 changes R_core by 2.895, 1.448, 0.724 km; criteria registered before the runs met at 1200 | **verified**, commit `e495e3cb` and `lunar:docs/LUNAR_PLANETPROFILE_WARREN_NUMERICAL_CLOSURE.md` |
| EOS domain at physical nodes | `'error'` holds on every active mantle and core node of the Warren reference | **verified**, `lunar:docs/LUNAR_PLANETPROFILE_WARREN_NUMERICAL_CLOSURE.md` §2–4 |
| Mantle table resolution | Unified Warren table at 300 × 300 against an independent 450 × 450 build: radial, mineralogical, shallow-boundary and mid-mantle convergence pass | **verified**, `lunar:docs/UNIFIED_STX24_RESOLUTION_CONVERGENCE.md` |
| Core-model sensitivity | The ordering of the compositions in predicted moment of inertia is unchanged over the tested core range | **verified**, `lunar:docs/LUNAR_PLANETPROFILE_CORE_SENSITIVITY.md` |
| Angular-mean density wrapper | With Δρ = 0 the closure is reproduced exactly; pressure feedback is 3.7–6.3 % of the frozen-pressure first-order shift on four structures | **verified**, `lunar:docs/STAGE_E9A5_ANGULAR_MEAN_DENSITY_RECLOSURE.md` §8; also summarised in report appendix A.4 (`lunar:docs/supaero_report/appendices/appendix_a.tex`) |
| Rebuilt structures | The 12 layered and FeO-rich structures rebuilt for the seismic and tidal analyses reproduce the moment of inertia, the core radius and the shell description of the earlier runs exactly | **verified**, report appendix A.4 |

**Audit recheck, 5 November 2026:** the 12 test files listed in Section 6 passed again (**122 passed**). The four November integration scripts reproduced the reported checks: predictive MoI and core radius on Warren BU16; native versus external density correction for a synthetic offset on that structure; prescribed crust seismic values without changing the mantle or structure; and 3-D property comparisons on `W_LOW_HET` and `ALT5_B1` (maximum relative difference below 5.8e−7 against the stored cubes). The exact/old bar-constant experiment was also reproduced. These are targeted checks on the stated cases. The historical resolution campaigns, the 12-structure rebuild campaign and `BuildTest` were not re-run in this audit.

The lunar production checks do not validate the unused components in Section 2.8 or the lateral/TidalPy and HP-ice branches as a whole. Those branches retain the verification limits stated in the contributions note.

---

## 5. Problems found in PlanetProfile (for upstream review)

1. **Bar conversion constants:** `Constants.bar2MPa` and `bar2GPa` use the atmosphere factor (1.01325). This displaces raw Perple_X table pressure axes and scales K_S and G_S converted from bar by +1.325 %. The problem is still on the inspected `main` and `genai` revisions. It is fixed on the lunar branch (`c70c1bf7`), where the seawater module now uses `Constants.P0` for the atmospheric offset. The [NIST conversion table](https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b9) distinguishes 1 bar = 100000 Pa from 1 standard atmosphere = 101325 Pa.

   Measured effect on one lunar structure (Warren Preferred mantle, field BU16, Garcia M3 crust, 5 November 2026):
   with the old constants, the core radius is 361.91 km instead of 357.56 km, and the predicted I/MR² is 0.393080
   instead of 0.393165. The old-constant result falls outside the project's adopted Williams model envelope 0.39311–0.39338 (Section 3.1).

2. **Shipped `Fe-S_3D_EOS.mat`**, statuses in `lunar:docs/FE_S_CORE_EOS_AUDIT.md` §4:
   - its pressure axis carries the same 1.01325 factor (confirmed);
   - K_S and G_S appear scaled by the same factor in 7 of 8 slices (supported numerically);
   - the Fe95S5 slice is broken: its moduli are 10⁵–10⁶ times too small;
   - the six slices Fe65 through Fe90 have local negative or outlier values; Fe100 has no such anomalies in the reported checks, although its conversion bias remains;
   - the original werami source set for this eight-slice MAT was not found in the checkout. The five tables used by the traceable alternative are a separate source set with different grid coverage.
   
   The bar fix does not repair this file, because its arrays are pre-scaled. The lunar branch now has a builder that
   rebuilds a traceable table from the werami tables of the repository (`stackFeSTraceable.py`, `93f7f1cb`).
3. **Silent clamping:** with `EXTRAP=False`, an out-of-domain (P, T) returned the table-edge value without warning. The lunar branch addresses it with `EOS_DOMAIN_POLICY` (`'legacy_clamp'` default, `'warn'`, `'error'`).
4. **Cached tables:** at production revision `e495e3cb`, `PerplexEOSStruct` short-circuits `__init__` for a table already in its cache and returns an object without its domain limits. The project works around this with one memoised object per table (`lunar:scripts/stage_e9a5/e9a5_lib.py`, `eos()`). Fixed automatically on the lunar branch by `1b227ca0`, without a feature flag.
5. **`PPLuna.py` was Io-derived** on the inspected `main`, `genai` and two port branches, including its header, radius and mass. It was not byte-identical to the contemporaneous `PPIo.py` on `main`/`genai`, which additionally set two orbital properties. The real lunar configuration is on the lunar branch.
6. **`Bulk.Cmeasured` is misnamed for the Moon:** it stores the mean moment (A+B+C)/3, not the polar moment C. Its ±0.0002 is "five times the formal error" (Konopliv et al. 1998), not 1σ.
7. **Core radius is quantised by the radial step:** the core radius is the grid node where the mass match falls, so it is quantised to one radial interval (1.448 km at `nSilMax = 1200`). Differences smaller than one step are not resolved.
8. **Nearest-neighbour lookup:** the default `'nearest'` lookup produces step artefacts where a phase boundary crosses a profile (`lunar:docs/STAGE_E3_M1_RESIDUAL_MANTLE.md`).

---

## 6. Reproducing

1. **Version and environment.** Use an isolated checkout at `e495e3cb` for the production code, retaining the recorded comment-only diff as provenance. Use `8e3efe07` or later (branch head `d6a34a24` on 5 November) for the November additions and the checks below. Record the SHA and working-tree diff, rather than only the moving branch name. The checked environment is `planetprofile`; the lunar study's scripts contain absolute paths to its own checkout that must be adapted on another machine; `PPLuna.py` no longer does since `d6a34a24`, which ships the temperature profile in `PlanetProfile/Default/Luna/`. The Warren reference reproduces its recorded production result on the November revision; this is not an all-body regression guarantee.

2. **Mantle tables and thermal inputs.** Provision the M0/M1 tables, the selected Santangelo inputs and their radial temperature reductions before running a closure. For the N2 additions, `lunar:scripts/dps_2026_followup/fu_n2_a_install_tables.py` copies the four Khan06S/ONeill91S M1 tables and assembles six L500 tables from **U500 upper-mantle residual tables plus M0 lower-mantle tables**. It assumes the earlier tables and thermal inputs are already present. It is not a general installer: its `tmp/n2_perplex/outputs` source directory was absent during the audit. Those source builds must be recovered or regenerated before that route is reproducible. The layered-table unit test reuses the final L500 table as its upper input; an independent rebuild from U500 + M0 remains **implemented, unverified** in this audit.

3. **Core table, separately.** The N2 installer does **not** install the core MAT. On the November PlanetProfile revision, run from the PlanetProfile repository root:

   ```bash
   python Thermodynamics/Perple_X/stackFeSTraceable.py
   ```

   This builds `Fe-S_3D_EOS_5slice_traceable_v1.mat` from the five shipped werami tables into the Perple_X cache. Alternatively, copy `lunar:results/derived_eos/core/Fe-S_3D_EOS_5slice_traceable_v1.mat` there. The rebuilt arrays match the study product exactly; their semantic SHA-256, using the documented array-shape convention, is `452147a82062916219dc8b83e53dcb53114ccfe8f60030b8afb4a9db72834445`. The same input product can be used with the production code revision.

4. **Closure and numerical reference.** `close_corrected()` imported by `lunar:scripts/stage_e9a6/e9a6_lib.py` performs a closure given the case's EOS, radial temperature, crust and supplied density correction. The correction must be derived by the corresponding angular-mean workflow; omitting it runs the uncorrected closure. The simple Warren BU16/M3 reference from `d_lib.close('Warren', 16)` gives R_core = 357.563375 km and I/MR² = 0.39316529266712535, with mass residual −0.006742574 %.

5. **Targeted tests.** From the PlanetProfile root, on the November revision, the following command gave **122 passed** on 5 November 2026:

   ```bash
   python -m pytest -q \
     tests/{lunar_crust_boundary,eos_domain_guard,eos_physical_node_guard,luna_config,prescribed_silicate_t,perplex_eos_cache,fe_s_traceable_core,layered_mantle_table,angular_mean_density,predict_moi,silicate_3d,crust_seismic}_test.py
   ```

   The brace expansion requires a shell such as bash or zsh. Some tests require local EOS tables and skip if they are absent; a checkout without those inputs must not report a skipped check as a reproduced result.

6. **Integration comparisons.** From the lunar repository root, using the November PlanetProfile revision, run the local scripts `scripts/planetprofile/verify_predict_moi.py`, `verify_native_density_correction.py`, `verify_crust_seismic.py` and `verify_silicate_3d.py`. The last script also needs the stored Stage-E7 cubes. Run `scripts/planetprofile/bar_constant_effect.py exact` and then `scripts/planetprofile/bar_constant_effect.py atm` in separate Python processes for the bar experiment. These scripts and their JSON results need to accompany the review because they are not committed (access note above).

`BuildTest`, an independent layered-table source rebuild and recovery of the historical 14 bar tests remain outstanding before claiming broader upstream validation.
