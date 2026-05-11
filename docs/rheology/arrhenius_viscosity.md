# Arrhenius Ice Viscosity

PlanetProfile can optionally use a temperature-dependent Arrhenius viscosity
law for supported ice EOS objects:

```text
eta(T) = eta_melt * exp[(Eact / R) * (1 / T - 1 / Tmelt)]
```

where `eta_melt` is the reference viscosity at the phase melting temperature,
`Eact` is activation energy in J/mol, `R` is the gas constant, and `Tmelt` is the
reference melting temperature.

## Flags

All flags default to `False`, so standard behavior is preserved unless a body
config enables the feature.

```python
Planet.Do.ARRHENIUS_VISCOSITY = False      # Legacy/global ice switch
Planet.Do.ARRHENIUS_VISCOSITY_Ih = False   # Ice Ih
Planet.Do.ARRHENIUS_VISCOSITY_III = False  # Ice III
Planet.Do.ARRHENIUS_VISCOSITY_V = False    # Ice V
Planet.Do.ARRHENIUS_VISCOSITY_VI = False   # Ice VI
```

The legacy `ARRHENIUS_VISCOSITY` flag enables every supported ice phase. For new
work, prefer the per-phase flags so Ice Ih and high-pressure ice behavior can be
reviewed independently.

## Parameters

Default activation energies come from `Constants.Eact_kJmol`. A PPBody config can
override one phase without changing the constants:

```python
Planet.Ocean.Eact_kJmol['Ih'] = 65.0  # kJ/mol
```

EOS/profile reference viscosities come from `Constants.etaMelt_Pas`. Reference
melting temperatures are currently approximate defaults for HP ices, with
computed `Bulk.TbIII_K` and `Bulk.TbV_K` used when available.

Arrhenius viscosity is kept separate from convection diagnostic viscosity. That
separation is intentional: changing diagnostic rheology should not silently
change EOS/profile `eta_Pas`, and Arrhenius profile flags should not silently
drive convection diagnostic `etaConv*` values.

## Example PPBody Config

Enable Arrhenius viscosity only for Ice Ih:

```python
Planet.Do.ARRHENIUS_VISCOSITY_Ih = True
Planet.Do.ARRHENIUS_VISCOSITY_III = False
Planet.Do.ARRHENIUS_VISCOSITY_V = False
Planet.Do.ARRHENIUS_VISCOSITY_VI = False
```

Enable Arrhenius viscosity only for high-pressure ices:

```python
Planet.Do.ARRHENIUS_VISCOSITY_Ih = False
Planet.Do.ARRHENIUS_VISCOSITY_III = True
Planet.Do.ARRHENIUS_VISCOSITY_V = True
Planet.Do.ARRHENIUS_VISCOSITY_VI = True
```

## Current Limitations

Arrhenius viscosity is an EOS/profile viscosity law. It changes saved `eta_Pas`
for enabled phases, but it does not overwrite the convection reference viscosity
inside the Arrhenius EOS object.
