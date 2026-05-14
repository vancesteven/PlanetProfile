# High-Pressure Ice Convection Diagnostics

PlanetProfile includes opt-in diagnostic calculations for in-ocean high-pressure
ice convection. These diagnostics inspect Ice III, Ice V, and Ice VI behavior
without changing propagated profiles or scalar model outputs.

These are diagnostic calculators only. They do not modify propagated thermal,
phase, mass, gravity, heat-flux, conductivity, or viscosity profiles, including
`T_K`, `phase`, `MLayer_kg`, `g_ms2`, `qSurf_Wm2`, `qCon_Wm2`,
`kTherm_WmK`, or `eta_Pas`.

## Flags

```python
Planet.Do.HP_ICE_CONVECTION_DIAGNOSTICS = False
Planet.Do.KALOUSOVA_CONVECTION = False
```

`HP_ICE_CONVECTION_DIAGNOSTICS=True` enables diagnostic reporting. When
`KALOUSOVA_CONVECTION=False`, diagnostics use a Deschamps and Sotin (2001)-style
path where possible, with a phase-local Ice VI fallback that is explicitly
diagnostic-only.

In this extraction, `KALOUSOVA_CONVECTION=True` enables Kalousova HP ice
convection diagnostics only. It does not update the propagated HP ice thermal
profile, heat flux, mass, gravity, or phase structure. The flag name is kept for
traceability with `origin/genai`, where the same name is used. In this branch,
the flag selects diagnostic-only calculations and should not be interpreted as a
production Kalousova convection model.

## Kalousova Rheology Parameters

Kalousova diagnostics use an explicit parameter set:

```python
Planet.Ocean.etaMeltKalousova_Pas = {
    'III': 5e12,
    'V': 2.8e14,
    'VI': 5e14,
}
Planet.Ocean.phiPercolationKalousova_frac = 0.05
```

The Ice V value follows the `genai` Kalousova reference value of
`2.8e14 Pa s`, but it is local to Kalousova diagnostics. The global
Ice V EOS/profile viscosity remains unchanged unless a user or future reviewed
feature explicitly changes it.

Kalousova does not consume saved profile `eta_Pas`. Arrhenius HP ice flags can
change profile `eta_Pas`, but they do not silently change Kalousova `etaConvIII`,
`etaConvV`, or `etaConvVI`.

The Kalousova scalings currently implemented for diagnostics are:

```text
RaCrit = 19.965e3 * qs_mWm2**3.690
Ht_km = (0.145e-3 * qs_mWm2 + 0.015) * etaMelt_Pas**0.21
delta_prime = 2.746 * Ra**-0.271
```

The heat flux `qs_mWm2` is converted from W/m2 to mW/m2 before applying the
Kalousova and Sotin scaling, and `Ht_km` is converted back to meters for the
stored diagnostic fields.

`phiPercolationKalousova_frac = 0.05` is a configured diagnostic marker assigned
when the diagnostic Kalousova scaling is supercritical. It is not an evolved
porosity field and it is not computed melt transport.

`DO_HP_MELT` and `meltFractionIII`, `meltFractionV`, and `meltFractionVI` are
diagnostic markers only. They do not trigger melt transport, latent-heat
feedback, layer remeshing, or profile mutation in this extraction.

## Diagnostic Outputs

The diagnostic fields include:

```text
TconvIII_K, TconvV_K, TconvVI_K
etaConvIII_Pas, etaConvV_Pas, etaConvVI_Pas
eLidIII_m, eLidV_m, eLidVI_m
DconvIII_m, DconvV_m, DconvVI_m
deltaTBLIII_m, deltaTBLV_m, deltaTBLVI_m
RaConvectIII, RaConvectV, RaConvectVI
RaCritIII, RaCritV, RaCritVI
meltFractionIII, meltFractionV, meltFractionVI
```

Absent phases are reported as absent/NaN rather than creating new layers.

For the Kalousova diagnostic path, the existing `TconvIII_K`, `TconvV_K`, and
`TconvVI_K` fields store the phase-top reference temperature used in the scaling.
They are kept under the existing field names for compatibility with the
PlanetProfile convection diagnostics, but they should not be interpreted as a
propagated DS2001-style convective temperature.

The helper `_SetHPIceDiagnosticFields()` is the single writer for both the
top-level HP diagnostic fields (`TconvV_K`, `etaConvVI_Pas`,
`meltFractionVI`, etc.) and the structured `Planet.HPIceDiagnostics` record.
Future changes should keep those two representations synchronized through that
helper.

## Future API Note

A future cleaner API could replace the two HP ice diagnostic booleans with a
single model selector, for example:

```python
Planet.Do.HP_ICE_CONVECTION_MODEL = "none" | "DS2001_diagnostic" | "Kalousova2018_diagnostic"
```

That refactor is intentionally deferred so this extraction remains close to the
`genai` flag names and minimal review scope.

## Limitations

These calculations are diagnostics. They do not implement production Ice VI
underplate propagation, HP ice mass-transport feedback, self-consistent tidal
heating, or unrelated broader `genai` features.

Ice VI diagnostics remain science-review pending. The current implementation
reports finite inspection values for existing Ice VI blocks, but it should not be
treated as a complete Ice VI convection model.

The Kalousova and Sotin mass-transport constants used in broader `genai`
experiments, such as latent heat `306e3 J/kg` and melt density `1270 kg/m3`,
are intentionally excluded here because this extraction does not calculate melt
transport.

Kalousova and Sotin's published two-phase work focuses primarily on Ice VI. The
Ice III and Ice V diagnostic calculations follow `genai` and should be treated
as extrapolative inspection values until reviewed. They are useful for
traceability and sensitivity inspection, but they should not be used as
validated production predictions without a separate science review.
