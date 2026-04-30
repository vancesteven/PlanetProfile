# HP Ice Convection Plan: Kalousova et al. (2018)

## Overview

Add Kalousova et al. (2018) parameterization for high-pressure (HP) ice convection as an **optional alternative** to the existing Deschamps & Sotin (2001) model. The new model is controlled by `Planet.Do.KALOUSOVA_CONVECTION = False` (default off), preserving all existing behavior.

Additionally, per-phase convection toggles have been added:
- `Planet.Do.NO_ICE_CONVECTION_Ih` — suppress ice Ih convection only
- `Planet.Do.NO_ICE_CONVECTION_III` — suppress ice III convection only
- `Planet.Do.NO_ICE_CONVECTION_V` — suppress ice V convection only
- `Planet.Do.NO_ICE_CONVECTION_VI` — suppress ice VI convection only
- `Planet.Do.NO_ICE_CONVECTION` — global suppression (existing, unchanged)

These are already implemented in `defineStructs.py` and `LayerPropagators.py`.

---

## Physical Basis

### Current Model: Deschamps & Sotin (2001)

The existing model (`ConvectionDeschampsSotin2001` in `ThermalProfiles.py:67-224`) treats each HP ice layer independently with:
- Stagnant lid convection with Arrhenius temperature-dependent viscosity
- Rayleigh number: Ra = αCpρ²gΔTz³/(ηk)
- Critical Ra from Solomatov (1995): Ra_crit = 20.9(E_act·ΔT/(RT_conv²))⁴
- If Ra > Ra_crit: convecting layer develops between stagnant lid (eLid) and lower TBL (deltaTBL)
- Returns: T_conv, η_conv, eLid, D_conv, deltaTBL, Q_bot, Ra, Ra_crit

### Proposed Model: Kalousova et al. (2018)

Kalousova et al. parameterize HP ice convection as a single mixed layer to determine conditions for partial melt formation. Key differences from Deschamps & Sotin:

1. **Two-phase convection insight**: When HP ices convect vigorously, temperatures at layer interfaces approach triple-point temperatures:
   - Ice III-V-liquid triple point: ~350 MPa, ~254 K (Journaux et al. 2020)
   - Ice V-VI-liquid triple point: ~632 MPa, ~272 K (Journaux et al. 2020)

2. **Melt formation criterion**: If the convective adiabat intersects the solidus, partial melt forms within the HP ice layer. The melt fraction depends on the Rayleigh number and layer thickness.

3. **Effective viscosity**: Uses composite rheology (diffusion + dislocation creep) rather than pure Arrhenius dependence. Parameters per ice phase:
   - Ice III: E_diff ~ 98 kJ/mol, E_disl ~ 190 kJ/mol
   - Ice V: E_diff ~ 136 kJ/mol, E_disl ~ 195 kJ/mol
   - Ice VI: E_diff ~ 110 kJ/mol, E_disl ~ 150 kJ/mol

4. **Temperature scaling**: Nu-Ra relationship modified for temperature-dependent viscosity with internal heating, following Deschamps & Sotin but with updated scaling coefficients.

---

## Implementation Design

### New Function: `ConvectionKalousova2018`

Location: `PlanetProfile/Thermodynamics/ThermalProfiles/ThermalProfiles.py`

**Signature** (mirrors Deschamps & Sotin for drop-in replacement):
```python
def ConvectionKalousova2018(Ttop_K, rTop_m, kTop_WmK, Tb_K, zb_m, 
                             gtop_ms2, Pmid_MPa, oceanEOS, iceEOS, 
                             phaseBot, EQUIL_Q, Eact_kJmol):
    """HP ice convection with melt-formation criterion (Kalousova et al. 2018).
    
    Returns: Tconv_K, etaConv_Pas, eLid_m, Dconv_m, deltaTBL_m, 
             Qbot_W, Ra, RaCrit
             
    Additional attributes set on return dict (or via Planet):
        meltFraction: Volume fraction of melt in convecting region
        DO_HP_MELT: Boolean, whether melt forms
    """
```

**Same 8 return values** as `ConvectionDeschampsSotin2001` for backward compatibility. Additional melt-related quantities stored on Planet struct.

### New Struct Fields

In `defineStructs.py`, add to the relevant ice/HP section:
```python
# HP ice melt (Kalousova 2018)
Planet.Ocean.meltFractionIII = 0.0  # Melt fraction in ice III layer
Planet.Ocean.meltFractionV = 0.0    # Melt fraction in ice V layer
Planet.Ocean.meltFractionVI = 0.0   # Melt fraction in ice VI layer
Planet.Ocean.DO_HP_MELT = False     # Whether HP ice partial melting occurs
```

### Dispatch Logic

In `Convection.py`, the existing functions `IceIIIConvectSolid`, `IceVConvectSolid`, `IceVIConvectSolid` call `ConvectionDeschampsSotin2001`. The dispatch changes to:

```python
# In IceIIIConvectSolid:
if Planet.Do.KALOUSOVA_CONVECTION:
    Planet.Tconv_K, Planet.etaConvIII_Pas, ... = \
        ConvectionKalousova2018(Ttop, rTop, kTop, TbIII, zbIII, ...)
else:
    Planet.Tconv_K, Planet.etaConvIII_Pas, ... = \
        ConvectionDeschampsSotin2001(Ttop, rTop, kTop, TbIII, zbIII, ...)
```

Same pattern for ice V and VI. **Ice Ih is NOT affected** — Kalousova applies only to HP ices.

### Triple-Point Temperature Boundaries

When `KALOUSOVA_CONVECTION = True`, the hot-side boundary temperature for each HP ice layer is set to the relevant triple-point temperature rather than the value propagated from the ocean:

| Layer | Hot boundary (Deschamps) | Hot boundary (Kalousova) |
|-------|--------------------------|--------------------------|
| Ice III | TbIII_K (from ocean EOS) | min(TbIII_K, 254 K) |
| Ice V | TbV_K (from ocean EOS) | min(TbV_K, 272 K) |
| Ice VI | TbVI_K (from ocean EOS) | TbVI_K (no change — deepest layer) |

Triple-point values should be stored as constants:
```python
# In Constants class (defineStructs.py)
Constants.TtripleIII_V_L_K = 254.0   # Ice III-V-liquid triple point
Constants.PtripleIII_V_L_MPa = 350.0
Constants.TtripleV_VI_L_K = 272.0    # Ice V-VI-liquid triple point
Constants.PtripleV_VI_L_MPa = 632.0
```

---

## Files to Modify

| File | Change |
|------|--------|
| `Utilities/defineStructs.py` | `Do.KALOUSOVA_CONVECTION` flag (DONE), melt fraction fields, triple-point constants |
| `ThermalProfiles/ThermalProfiles.py` | New `ConvectionKalousova2018` function |
| `ThermalProfiles/Convection.py` | Dispatch to Kalousova in IceIII/V/VI ConvectSolid/Porous |
| `Thermodynamics/LayerPropagators.py` | Per-phase convection toggles (DONE) |

---

## Implementation Steps

1. **Add triple-point constants** to `Constants` in `defineStructs.py`
2. **Add melt fraction fields** to `PlanetStruct` 
3. **Implement `ConvectionKalousova2018`** in `ThermalProfiles.py`:
   - Read Kalousova et al. 2018 for exact parameterization coefficients
   - Implement composite rheology (diffusion + dislocation creep)
   - Implement melt formation criterion
   - Return same 8 values as Deschamps & Sotin
4. **Add dispatch logic** in `Convection.py` IceIII/V/VI functions
5. **Add test** — PPTest with `KALOUSOVA_CONVECTION = True` for a Ganymede-like body (thick HP ice layers)
6. **Verify** — Run with flag off to confirm no change to existing results; run with flag on for Ganymede

---

## Key Differences Summary

| Feature | Deschamps & Sotin (2001) | Kalousova et al. (2018) |
|---------|--------------------------|-------------------------|
| Rheology | Arrhenius (single mechanism) | Composite (diffusion + dislocation) |
| Boundary T | From ocean EOS | Triple-point temperatures |
| Melt formation | Not modeled | Predicted from convective state |
| Applicable phases | All ice phases | HP ices only (III, V, VI) |
| Default | Yes (existing) | No (opt-in via config flag) |

---

## Status

- [x] `Planet.Do.KALOUSOVA_CONVECTION` config flag added
- [x] Per-phase convection toggles (`NO_ICE_CONVECTION_Ih/III/V/VI`) added
- [x] `LayerPropagators.py` updated to use per-phase flags
- [ ] Triple-point constants in `defineStructs.py`
- [ ] Melt fraction struct fields
- [ ] `ConvectionKalousova2018` function
- [ ] Dispatch logic in `Convection.py`
- [ ] Test case (Ganymede-like)
- [ ] Verification
