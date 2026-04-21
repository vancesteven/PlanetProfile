# Session Handoff

## 2026-04-20 (Session 4 of Arrhenius/Tidal work)

### Branch
`genai` -- Arrhenius HP ice viscosity, self-consistent tidal heating, and related bug fixes

### Current Objective
All planned Arrhenius/tidal tasks are COMPLETE, including Ice V etaMelt correction.

### Progress This Session (Session 4)

**Completed:**
1. Lowered Ice V `etaMelt_Pas` from 5e14 to 2.8e14 Pa·s in defineStructs.py, citing Durham & Stern (2001) flow law analysis. This pushes Ice V past the critical Rayleigh number threshold.
2. Re-ran PPTest35: Ice V now convects (Dconv=134.1 km, eLid=17.8 km). All three HP ice phases are supercritical.

### Progress Sessions 1-3

**Completed:**
1. Wired Arrhenius viscosity through `GetOceanHPIceEOS` in LayerPropagators.py -- both MgSO4 and general code paths now pass `**_arrheniusKwargs[phaseID]` to `GetIceEOS()` for phases III, V, VI when `Planet.Do.ARRHENIUS_VISCOSITY = True`.
2. Fixed ConvectionPlots.py tidal formulas: eps0 = (3/2)*e*n^2*R/g and Maxwell D = omega^2*eta*mu^2/(mu^2+omega^2*eta^2).
3. Implemented self-consistent tidal heating switch in Gravity.py: when `DO_SELF_CONSISTENT_HTIDAL=True`, overrides `HtidalIce_Wm3` with k2-implied value.
4. Updated PPTest35 with `ARRHENIUS_VISCOSITY = True` and `etaRock_Pas = [1e22, 1e22]` (Petricca 2025).
5. Fixed subcritical convection eLid bug: eLid = zb_m (not 0) in subcritical case.
6. Fixed spurious "temperate layer" message for subcritical ice layers (added `Dconv_m > 0` guard).
7. Fixed waterless body `CalcMoIConstantRho` crash (empty dChydro guard).
8. Fixed nan Pmid crash in `ConvectionDeschampsSotin2001` (early return when Pmid becomes nan).
9. Fixed whole-shell clathrate `surfIceEOS['Ih']` KeyError (fall back to any loaded EOS).
10. Ran full BuildTest: 78+ profiles, 0 errors.

### PPTest35 Results (Titan no-ocean, Arrhenius HP ice)
**After Ice V etaMelt correction (2.8e14 Pa·s):**
- Ice III: RaCrit=7.86e7, eLid=7.6 km, Dconv=74.3 km -- SUPERCRITICAL
- Ice V: RaCrit=1.12e8, eLid=17.8 km, Dconv=134.1 km -- SUPERCRITICAL (was subcritical at 5e14)
- Ice VI: RaCrit=7.03e8, eLid=21.3 km, Dconv=513.4 km -- SUPERCRITICAL

### Decision: Tmelt_K for Arrhenius
Used fixed reference melting temperatures: Ice III=253K, Ice V=264K, Ice VI=290K.

### Remaining Opportunities
- Verify Arrhenius viscosity is actually temperature-dependent in HP ice layers (check pkl output)
- The surface heat flux (5.8 mW/m^2) is below Petricca's ~40 mW/m^2 -- may need higher tidal heating
- Full self-consistent iteration (currently single-pass): re-run with updated HtidalIce_Wm3
- Test 28 (the last standard test) was not reached during BuildTest timeout

### Files Modified (All Sessions)

| File | Changes |
|------|---------|
| `PlanetProfile/Utilities/defineStructs.py` | Ice V etaMelt_Pas: 5e14 → 2.8e14 (Durham & Stern 2001) |
| `PlanetProfile/Thermodynamics/LayerPropagators.py` | +Arrhenius kwargs block in GetOceanHPIceEOS, +Dconv>0 guard for melt, +nan/thin layer guard, +CalcMoI empty dChydro fix, +surfIceEOS fallback |
| `PlanetProfile/Plotting/ConvectionPlots.py` | Fixed eps0 and Maxwell D formulas |
| `PlanetProfile/Gravity/Gravity.py` | Added DO_SELF_CONSISTENT_HTIDAL override branch |
| `PlanetProfile/Thermodynamics/ThermalProfiles/ThermalProfiles.py` | Fixed subcritical eLid=zb_m, +nan Pmid early return |
| `PlanetProfile/Test/PPTest35.py` | Added ARRHENIUS_VISCOSITY=True, etaRock_Pas=[1e22,1e22] |
