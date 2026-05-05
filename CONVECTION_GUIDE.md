# Ice Shell Convection Models in PlanetProfile

## 1. Overview

PlanetProfile offers three convection parameterizations for modeling heat transport in ice shells. The choice of model affects stagnant lid thickness, interior temperature profile, and (for high-pressure ice) whether partial melting is predicted.

| Model | Reference | Applicable Phases | Enable Flag |
|---|---|---|---|
| Deschamps & Sotin (2001) | DS2001 | Ice Ih, III, V | Default (no flag needed) |
| Yao et al. (2014) | Yao2014 | Ice Ih only | `Planet.Do.SPHERICAL_CONVECTION = True` |
| Kalousova & Sotin (2018) | Kalousova | Ice III, V, VI | `Planet.Do.KALOUSOVA_CONVECTION = True` |

**Important note on Ice VI:** Deschamps & Sotin is not implemented for Ice VI. If `KALOUSOVA_CONVECTION` is not enabled, an Ice VI layer is treated as purely conductive.

---

## 2. Deschamps & Sotin (2001) — Default Model

### Description

DS2001 is the default convection parameterization and applies Cartesian stagnant-lid scaling laws to ice layers. It divides the ice shell into three thermal regions:

- **Stagnant (conductive) lid** at the top, where viscosity contrasts suppress convection
- **Convecting sublayer** in the interior, with near-adiabatic temperature profile
- **Bottom thermal boundary layer** adjacent to the ocean or silicate interface

### Applicable Phases

Ice Ih, Ice III, Ice V (Ice VI is not supported; see note above).

### Configuration

No special flags are required. DS2001 is used by default for all ice phases unless overridden.

```python
# In PPBody.py — DS2001 is active unless you enable an alternative
# No flag required; the following are the defaults:
# Planet.Do.SPHERICAL_CONVECTION = False  (default)
# Planet.Do.KALOUSOVA_CONVECTION = False  (default)
```

### Key outputs

```python
Planet.eLid_m        # Stagnant lid thickness (m)
Planet.Dconv_m       # Convecting layer thickness (m)
Planet.deltaTBL_m    # Bottom thermal boundary layer thickness (m)
Planet.Tconv_K       # Interior convective temperature (K)
Planet.RaConvect     # Rayleigh number
Planet.RaCrit        # Critical Rayleigh number
```

### Reference

Deschamps, F., & Sotin, C. (2001). Thermal convection in the outer shell of large icy satellites. *Journal of Geophysical Research*, 106(E3), 5107–5121. https://doi.org/10.1029/2000JE001253

---

## 3. Yao et al. (2014) — Spherical Shell Scaling

### Description

The Yao2014 model replaces the default Cartesian scaling laws with 3D spherical shell scaling for Ice Ih. It introduces a curvature parameter f = R_base / R_top that accounts for the ratio of inner to outer shell radius. When f is significantly less than 1 (thick shells), spherical geometry substantially changes the predicted lid thickness and surface heat flux relative to DS2001.

Yao2014 is most impactful for bodies with thick ice shells:

- **Pluto** — f ≈ 0.75 (very thick shell, largest correction)
- **Titan** — f ≈ 0.96 (moderately thick shell)
- **Europa** — f ≈ 0.998 (thin shell, Cartesian and spherical agree closely)

For thin shells, Yao2014 and DS2001 produce nearly identical results. For thick shells, Yao2014 predicts **thicker stagnant lids** and **lower surface heat flux** than DS2001.

### Applicable Phases

Ice Ih only. DS2001 continues to be used for Ice III and V even when `SPHERICAL_CONVECTION` is enabled.

### Configuration

```python
# In PPBody.py
Planet.Do.SPHERICAL_CONVECTION = True
```

Enabling this flag automatically activates Arrhenius (temperature-dependent) viscosity for Ice Ih via the Frank-Kamenetskii approximation. You do not need to set `ARRHENIUS_VISCOSITY_Ih` separately.

You can also enable Arrhenius viscosity independently per phase without enabling spherical convection:

```python
Planet.Do.ARRHENIUS_VISCOSITY_Ih = True   # Arrhenius for Ice Ih (DS2001 geometry)
Planet.Do.ARRHENIUS_VISCOSITY_III = True  # Arrhenius for Ice III
Planet.Do.ARRHENIUS_VISCOSITY_V = True    # Arrhenius for Ice V
Planet.Do.ARRHENIUS_VISCOSITY_VI = True   # Arrhenius for Ice VI
```

### Key outputs

Same output variables as DS2001 (see Section 2). The curvature parameter f is logged at debug level when the model runs.

### Reference

Yao, C., Deschamps, F., Lowman, J. P., Sanchez-Valle, C., & Tackley, P. J. (2014). Stagnant lid convection in bottom-heated thin spherical shells: Influence of curvature and implications for dwarf planets and icy moons. *Journal of Geophysical Research: Planets*, 119(8), 1895–1913. https://doi.org/10.1002/2014JE004653

---

## 4. Kalousova & Sotin (2018) — HP Ice with Partial Melt

### Description

The Kalousova model is designed specifically for high-pressure ice phases (III, V, VI) in large ocean worlds. It implements a two-phase convection parameterization in which vigorous convection drives interior temperatures toward the local melting curve (solidus) rather than an adiabat. When convection is sufficiently strong (modified Rayleigh number Ra* exceeds a critical value Ra*c), a **temperate layer** forms at the top of the HP ice layer containing a mixture of ice and liquid water.

Key differences from DS2001:

- Interior temperature follows the ice melting curve, not an adiabat
- Temperate (partially molten) layer can form if Ra* > Ra*c
- Scaling laws depend on heat flux from the silicate layer
- Boundary temperatures are limited by ice-water triple-point temperatures (Ice III-V-liquid at 254 K; Ice V-VI-liquid at 272 K)

### Applicable Phases

Ice III, Ice V, Ice VI.

### Configuration

```python
# In PPBody.py
Planet.Do.KALOUSOVA_CONVECTION = True
```

To use a custom reference viscosity for melt-bearing ice:

```python
# Optional: override reference viscosity per phase (Pa*s)
Planet.Ocean.etaMeltKalousova_Pas = {3: 5e12, 5: 5e12, 6: 5e12}
```

### Key outputs

```python
# Global melt flag
Planet.DO_HP_MELT          # True if any HP ice layer has a temperate zone

# Per-phase temperate layer thickness (m); 0 if no melt
Planet.eLidIII_m
Planet.eLidV_m
Planet.eLidVI_m

# Per-phase melt fractions (nominal 0.01 when temperate layer detected)
Planet.meltFractionIII
Planet.meltFractionV
Planet.meltFractionVI

# Per-phase convection diagnostics
Planet.TconvIII_K, Planet.TconvV_K, Planet.TconvVI_K   # Interior temperatures (K)
Planet.RaConvectIII, Planet.RaConvectV, Planet.RaConvectVI  # Modified Rayleigh numbers
Planet.RaCritIII, Planet.RaCritV, Planet.RaCritVI          # Critical Rayleigh numbers
```

Checking for partial melt in post-processing:

```python
if Planet.DO_HP_MELT:
    print(f"Ice III temperate layer: {Planet.eLidIII_m / 1e3:.1f} km")
    print(f"Ice V temperate layer:   {Planet.eLidV_m / 1e3:.1f} km")
    print(f"Ice VI temperate layer:  {Planet.eLidVI_m / 1e3:.1f} km")
```

### Note on Ice VI

Deschamps & Sotin is not implemented for Ice VI. If `KALOUSOVA_CONVECTION = False` (the default), an Ice VI layer is treated as purely conductive. To model Ice VI convection, set `KALOUSOVA_CONVECTION = True`.

### Reference

Kalousová, K., & Sotin, C. (2018). Melting in high-pressure ice layers of large ocean worlds—Implications for volatiles transport. *Geophysical Research Letters*, 45(16), 8096–8103. https://doi.org/10.1029/2018GL078889

---

## 5. Suppressing Convection

Any ice phase can be forced into a purely conductive profile by setting the appropriate flag. This is useful for sensitivity studies or when convection onset is uncertain.

```python
# Suppress convection in all ice phases
Planet.Do.NO_ICE_CONVECTION = True

# Suppress convection in specific HP ice phases only
Planet.Do.NO_ICE_CONVECTION_III = True   # Ice III treated as conductive
Planet.Do.NO_ICE_CONVECTION_V = True     # Ice V treated as conductive
Planet.Do.NO_ICE_CONVECTION_VI = True    # Ice VI treated as conductive

# Suppress convection in Ice Ih only
Planet.Do.NO_ICE_CONVECTION_Ih = True
```

When a phase-specific flag is set, the global `NO_ICE_CONVECTION` flag is overridden for that phase. The layer receives a conductive (linear or parameterized) temperature profile from top to bottom boundary temperature.

---

## 6. Model Comparison

| Feature | DS2001 | Yao2014 | Kalousova |
|---|---|---|---|
| Applicable phases | Ih, III, V | Ih only | III, V, VI |
| Geometry | Cartesian | 3D spherical shell | Cartesian |
| Interior T profile | Adiabatic | Adiabatic | Melting curve (solidus) |
| Melt prediction | Not modeled | Not modeled | Yes, via Ra* > Ra*c criterion |
| Scaling basis | Cartesian stagnant-lid laws | Spherical curvature-corrected laws | Heat-flux-dependent two-phase laws |
| Key advantage | Broad phase coverage, fast | Accurate for thick shells (Pluto, Titan) | HP ice melt tracking for large ocean worlds |
| Arrhenius viscosity | Optional (per-phase flags) | Automatic for Ih | Implicit in reference viscosity |

---

## 7. References

Deschamps, F., & Sotin, C. (2001). Thermal convection in the outer shell of large icy satellites. *Journal of Geophysical Research*, 106(E3), 5107–5121. https://doi.org/10.1029/2000JE001253

Yao, C., Deschamps, F., Lowman, J. P., Sanchez-Valle, C., & Tackley, P. J. (2014). Stagnant lid convection in bottom-heated thin spherical shells: Influence of curvature and implications for dwarf planets and icy moons. *Journal of Geophysical Research: Planets*, 119(8), 1895–1913. https://doi.org/10.1002/2014JE004653

Kalousová, K., & Sotin, C. (2018). Melting in high-pressure ice layers of large ocean worlds—Implications for volatiles transport. *Geophysical Research Letters*, 45(16), 8096–8103. https://doi.org/10.1029/2018GL078889
