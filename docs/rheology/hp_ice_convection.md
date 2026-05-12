# High-Pressure Ice Convection Diagnostics

PlanetProfile can optionally report high-pressure ice convection diagnostics for
in-ocean Ice III, Ice V, and Ice VI layers. The diagnostics are disabled by
default:

```python
Planet.Do.HP_ICE_CONVECTION_DIAGNOSTICS = False
```

When enabled, the model computes diagnostic-only quantities after the ocean and
HP ice layers have already been propagated. These calculations do not update the
propagated thermal profile, phase structure, heat flux, mass, gravity, or layer
thickness solution.

## Diagnostic Method

The current diagnostic path uses the existing Deschamps and Sotin (2001)-style
PlanetProfile convection scaling. Diagnostics are evaluated phase by phase in
bottom-to-top HP ice order:

```text
Ice VI -> Ice V -> Ice III
```

Ice III and Ice V use the existing DS2001 helper with a phase-fixed EOS wrapper
so the diagnostic convective temperature cannot alter the saved phase profile.
Ice VI uses a phase-local DS2001 fallback because the production DS2001 melt
lookup assumes a low-temperature pure-water search range that is not valid for
these Ice VI diagnostic blocks.

Ice VI values should be interpreted as phase-local inspection values only. They
are not a validated production Ice VI convection model.

## Reported Fields

The diagnostics populate top-level fields where available:

```text
TconvIII_K, TconvV_K, TconvVI_K
etaConvIII_Pas, etaConvV_Pas, etaConvVI_Pas
etaMeltIII_Pas, etaMeltV_Pas, etaMeltVI_Pas
eLidIII_m, eLidV_m, eLidVI_m
DconvIII_m, DconvV_m, DconvVI_m
deltaTBLIII_m, deltaTBLV_m, deltaTBLVI_m
RaConvectIII, RaConvectV, RaConvectVI
RaCritIII, RaCritV, RaCritVI
```

The same values are also recorded in:

```python
Planet.HPIceDiagnostics
```

`_SetHPIceDiagnosticFields()` is the single writer for both the top-level fields
and `Planet.HPIceDiagnostics`, which avoids the two storage paths drifting apart.

## Excluded From This Diagnostic Block

This feature does not include Kalousova two-phase HP ice convection, melt
transport, mass-flux or outflow diagnostics, latent-heat transport, self-consistent
tidal heating, Yao 2014 Ice Ih spherical convection, or Ice VI production
propagation stubs. Those require separate feature reviews.
