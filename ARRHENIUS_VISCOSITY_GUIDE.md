# Arrhenius Viscosity Guide

## Overview

PlanetProfile supports temperature-dependent (Arrhenius) viscosity for ice layers. Rather than
using a constant reference viscosity, the Arrhenius model scales viscosity exponentially with
temperature:

    η(T) = η_melt · exp[ Eact/(R·T) − Eact/(R·T_melt) ]

where `η_melt` is the viscosity at the melting temperature, `Eact` is the activation energy,
`R` is the gas constant, and `T_melt` is a reference melting temperature for the phase.

This model can be enabled independently of which convection model is active. It affects both
conductive and convective layer calculations: conductive layers use temperature-dependent
viscosity in their thermal profiles, and convective layers use it when computing Rayleigh
numbers and stagnant lid thicknesses.

---

## Configuration

Arrhenius viscosity is controlled by per-phase flags in `Planet.Do`. All default to `False`.

| Flag | Phase |
|------|-------|
| `Planet.Do.ARRHENIUS_VISCOSITY_Ih` | Ice Ih |
| `Planet.Do.ARRHENIUS_VISCOSITY_III` | Ice III |
| `Planet.Do.ARRHENIUS_VISCOSITY_V` | Ice V |
| `Planet.Do.ARRHENIUS_VISCOSITY_VI` | Ice VI |
| `Planet.Do.ARRHENIUS_VISCOSITY_sil` | Silicates (defined, not yet implemented) |
| `Planet.Do.ARRHENIUS_VISCOSITY` | Legacy: enables all ice phases at once |

The legacy flag `ARRHENIUS_VISCOSITY` is OR'd with each per-phase flag, so setting it `True`
enables Arrhenius for Ice Ih, III, V, and VI simultaneously. Prefer the per-phase flags for
new configurations.

**Example PPBody.py snippet — enable only for HP ice phases:**

```python
# Enable Arrhenius viscosity for high-pressure ice phases only
Planet.Do.ARRHENIUS_VISCOSITY_III = True
Planet.Do.ARRHENIUS_VISCOSITY_V   = True
Planet.Do.ARRHENIUS_VISCOSITY_VI  = True

# Leave Ice Ih at constant reference viscosity
Planet.Do.ARRHENIUS_VISCOSITY_Ih  = False
```

**Example — enable for all ice phases (legacy shorthand):**

```python
Planet.Do.ARRHENIUS_VISCOSITY = True
```

---

## Physics

The activation energies `Eact` (in kJ/mol) are drawn from laboratory measurements and are
stored in `Constants.Eact_kJmol`. Default values used by PlanetProfile:

| Phase | Eact (kJ/mol) | Source |
|-------|--------------|--------|
| Ice Ih | 59.4 | Lab creep experiments |
| Ice III | 127 | Lab creep experiments |
| Ice V | 136 | Lab creep experiments |
| Ice VI | 110 | Lab creep experiments |
| Clathrate | 90 | Durham et al. (2003) |

Reference melting temperatures used internally for HP ice phases:

| Phase | T_melt reference (K) |
|-------|---------------------|
| Ice III | ~253 K (at ~300 MPa) |
| Ice V | ~264 K (at ~500 MPa) |
| Ice VI | ~290 K (at ~1100 MPa) |

These reference temperatures are approximate mid-range values. The Arrhenius formula is
primarily sensitive to the ratio `(1/T − 1/T_melt)`, so the effect is strongest in cold,
thick conductive lids far below the melting point.

You can override the default activation energies per phase in the body config:

```python
# Override activation energy for Ice Ih (in kJ/mol)
Planet.Ocean.Eact_kJmol['Ih'] = 65.0
```

---

## Interaction with Convection

When a convecting layer is present, Arrhenius viscosity affects two key quantities:

- **Rayleigh number**: A larger viscosity contrast between the cold lid and the warm
  interior (larger `Eact`) increases the effective Rayleigh number, making convection onset
  easier.
- **Stagnant lid thickness**: Stronger temperature dependence produces a thicker stagnant
  (conductive) lid at the top of the convecting layer.

**Note on SPHERICAL_CONVECTION (Yao 2014):** Enabling `Planet.Do.SPHERICAL_CONVECTION`
automatically activates Arrhenius viscosity for Ice Ih. You do not need to set
`ARRHENIUS_VISCOSITY_Ih` separately when using the spherical convection model.

**Note on Kalousova convection:** When `Planet.Do.KALOUSOVA_CONVECTION = True` and
partial melt is predicted in HP ice, the effective reference viscosity `η_melt` is
automatically reduced to melt-bearing values. Arrhenius scaling is then applied on top of
this lower reference, giving a physically consistent treatment of partially molten HP layers.

---

## Interaction with Conduction

Even in purely conductive ice layers (convection off or below the Rayleigh threshold),
enabling Arrhenius viscosity causes the layer thermal profile to use temperature-dependent
viscosity rather than a constant value. This is why the flag is separate from the convection
flags: it controls viscosity representation across the full layer stack, not only in actively
convecting regions.

---

## Usage Examples

**Enable Arrhenius for all ice (simplest approach):**
```python
Planet.Do.ARRHENIUS_VISCOSITY = True
```

**Enable only for HP ice (common for large ocean worlds with thick HP ice shells):**
```python
Planet.Do.ARRHENIUS_VISCOSITY_III = True
Planet.Do.ARRHENIUS_VISCOSITY_V   = True
Planet.Do.ARRHENIUS_VISCOSITY_VI  = True
```

**Enable with Kalousova convection (recommended combination for HP ice):**
```python
Planet.Do.KALOUSOVA_CONVECTION    = True
Planet.Do.ARRHENIUS_VISCOSITY_III = True
Planet.Do.ARRHENIUS_VISCOSITY_V   = True
Planet.Do.ARRHENIUS_VISCOSITY_VI  = True
```

**Enable with spherical convection for Ice Ih (Arrhenius is set automatically):**
```python
Planet.Do.SPHERICAL_CONVECTION = True
# ARRHENIUS_VISCOSITY_Ih is forced True internally; no need to set it explicitly
```
