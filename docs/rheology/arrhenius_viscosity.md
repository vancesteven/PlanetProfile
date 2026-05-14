# Arrhenius Ice Viscosity

PlanetProfile can optionally use a temperature-dependent Arrhenius viscosity
law for supported ice EOS objects:

```text
eta(T) = eta_melt * exp[(Eact / R) * (1 / T - 1 / Tmelt)]
```

where `eta_melt` is the reference viscosity at the phase melting temperature,
`Eact` is activation energy in J/mol, `R` is the gas constant, and `Tmelt` is the
reference melting temperature.

## Validation Status

Arrhenius viscosity has passed integration-safety validation for the current
extraction: default behavior remains off by default, enabled cases leave
non-viscosity propagated arrays unchanged, and per-phase flags are phase-local.
This is not yet a full validation of self-consistent Arrhenius convection.

The current safety contract is:

- enabling Ice Ih Arrhenius changes `eta_Pas` only in Ice Ih layers;
- enabling high-pressure ice Arrhenius changes `eta_Pas` only in enabled HP ice
  phases;
- the legacy global flag changes all supported ice phases;
- HP ice diagnostics and Kalousova diagnostics remain disabled/empty in
  Arrhenius-only validation cases.

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

Per-phase flags are expected to be phase-local. A disabled phase should retain
the default viscosity profile even if another phase has Arrhenius enabled.

## Parameters

Default activation energies come from `Constants.Eact_kJmol`. A PPBody config can
override one phase without changing the constants:

```python
Planet.Ocean.Eact_kJmol['Ih'] = 65.0  # kJ/mol
```

EOS/profile reference viscosities come from `Constants.etaMelt_Pas`. The
`Tmelt_K` values used by the Arrhenius helper are fallback reference
temperatures, not pressure-dependent melting curves:

```text
Ice Ih:  Constants.T0
Ice III: 253 K
Ice V:   264 K
Ice VI:  290 K
```

Body-specific `Bulk.TbIII_K` and `Bulk.TbV_K` values override the Ice III and
Ice V fallbacks when available. Ice VI currently remains a fixed fallback
reference pending fuller Ice VI thermal/profile support.

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

Ice Ih frequently reaches the current clipping ceiling
`etaMax_Pas = 1e21 Pa s` in validation cases. That ceiling is a scientific
sensitivity item and should be reviewed before interpreting cold-Ice Ih
viscosity magnitudes.

Ice VI is wired for high-pressure ice EOS construction, but the current surface
thermal-profile code does not have a separate production Ice VI lithosphere
conduction/convection path. Its Arrhenius `Tmelt_K` value is therefore a fixed
fallback reference, not a pressure-dependent Ice VI melting curve.

PlanetProfile now resets mutable ice-viscosity and run-local EOS cache state at
independent `PlanetProfile()` run boundaries. This prevents one body run from
contaminating the next while preserving legacy within-run viscosity handoff
behavior.
