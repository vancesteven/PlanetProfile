# SeaFreeze direct evaluation vs PlanetProfile lookup tables (2026-11-02)

User request: "seafreeze can now compute individual properties quickly. Test
doing that instead of using lookup tables, including using whichphase
instead of the lookup table."

Status: measurements `verified` (scripts below). No default changed.

## What PlanetProfile does today

- **Properties.** `HydroEOS.OceanEOSStruct` / `IceEOSStruct` call
  `seafreeze()` ONCE on a (P, T) grid and wrap each property in a
  `RectBivariateSpline`. Every later evaluation is a spline call.
- **Phase.** With `Ocean.phaseType = 'lookup'` (the default), `whichphase`
  runs once on a grid and the result is wrapped in a nearest-neighbour
  `Nearest2DInterpolator`. With `phaseType = 'calc'`, `HydroEOS.SFphase`
  calls `whichphase` directly at every query. That option already existed.
- **Composition.** The Seawater composition (the Europa and Enceladus
  defaults) uses GSW and its own phase finder (`SwPhase`), so `phaseType`
  only affects its melting-curve EOS. The direct-phase change matters
  mainly for PureH2O, NaCl and MgSO4 oceans and for the ice phases.

## Measurements (SeaFreeze 1.1.3, this Mac)

Script: `plans/scripts/seafreeze_direct_benchmark.py`.

| Call | n = 1,000 | n = 10,000 |
|---|---|---|
| `seafreeze(PT, 'Ih')`, all properties | 333 ms | 3.1 s |
| `getProp(PT, 'Ih', path, 'rho')` | 61 ms | 0.56 s |
| `getProp(PT, 'Ih', path, 'rho','Cp','alpha')` | 146 ms | 1.4 s |
| `seafreeze(PT, 'water1')`, all properties | 1.07 s | 10.9 s |
| `getProp(PT, 'water1', path, 'rho')` | 162 ms | 1.66 s |
| `whichphase(PT)` (scatter) | 477 ms | 4.6 s (~0.46 ms/pt) |
| PP spline lookup, `rho` (`RectBivariateSpline`) | — | **1.4 ms** |

- **Per-property speedup is real:** requesting only rho is 5–7× faster than
  all properties, and three properties 2.3×. Values are identical
  (`rho` matches to 1e-9).
- **Direct evaluation is still ~1,000× slower per point than PP's spline
  lookup** (1.66 s vs 1.4 ms for 10k water rho). The spline reproduces the
  direct value to 2.3e-11 relative on a 1 MPa × 0.2 K grid, i.e. the
  lookup costs no accuracy for properties. Building the table takes 16 ms.

Phase finder (script `plans/scripts/seafreeze_phase_lookup_vs_direct.py`):

| Lookup grid (dP, dT) | Build | Disagreement with direct `whichphase`, 20k random points |
|---|---|---|
| 1.0 MPa, 0.5 K | 0.1 s | 0.21 % |
| 0.5 MPa, 0.25 K | 0.2 s | 0.07 % |
| 0.25 MPa, 0.1 K | 0.9 s | 0.03 % |

- Every disagreement lies within half a grid step of a phase boundary. The
  Ih upper boundary at 250 / 260 / 268 K is shifted by −0.19 / −0.41 /
  +0.02 MPa at dP = 1 MPa.

Full PP default runs, `phaseType` lookup → calc (script
`plans/scripts/seafreeze_phase_mode_runs.py`, output `phase_mode_runs.json`):

| Body (ocean) | lookup | calc | Change |
|---|---|---|---|
| Europa (Seawater) | 3.0 s | 1.0 s | identical outputs |
| Enceladus (Seawater) | 5.2 s | 4.4 s | identical outputs |
| Ganymede (PureH2O) | 2.7 s | **38.3 s** | zb +0.013 km, D_ocean +2.2 km (0.9 %), CMR2 −6e-5, k2 +0.0002 |

- **Defect found and fixed:** `phaseType = 'calc'` crashed on default
  Ganymede (`TypeError: only 0-dimensional arrays…`). `SFphase` returned a
  1-element array for a scalar query, and `GetPfreeze`'s brentq multiplies
  two results. It now returns a scalar for scalar input, matching the lookup
  finder. Test: `tests/seafreeze_direct_phase_test.py` (3 passed).
- The Seawater runs are faster under 'calc' only because they skip building
  the melt-EOS phase grid; their phase finder is GSW either way.

## Recommendation

1. **Keep spline lookup for PROPERTIES.** It is exact to 1e-11 and ~1,000×
   faster per point. SeaFreeze's per-property calls are the right tool for
   BUILDING the tables (5–7× faster when only some properties are needed),
   not for replacing them. Follow-up: switch the table builds in
   `HydroEOS` to `getProp(..., 'rho', 'Cp', 'alpha', ...)` with the
   property list each EOS actually uses. That is a pure speedup with
   identical values, worth a test that pins equality.
2. **Phase: offer direct `whichphase` as an accuracy option, not the
   default.**
   - It removes the half-grid-step boundary error (≤ 0.5 MPa at the default
     resolution).
   - Costs: 14× slower for the PureH2O Ganymede run (38 s vs 2.7 s).
   - Effect: it moved the Ganymede ocean thickness by 2.2 km, a real effect
     of boundary resolution on HP-ice bodies.
   - A cheaper middle path is a finer lookup grid: 0.25 MPa costs 0.9 s to
     build and cuts disagreement to 0.03 %. Alternatively, a hybrid that
     uses the table except near boundaries and calls `whichphase` only
     there.
3. **Inference caches** are built once and reused, so 'calc' (or a finer
   grid) is affordable there and would remove a small systematic. Decide
   per campaign; do not change cached artifacts retroactively.
