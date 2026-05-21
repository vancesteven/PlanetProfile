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
Planet.Do.HP_ICE_CONVECTION_MODEL = "none"
Planet.Do.ALLOW_EXPERIMENTAL_HP_KALOUSOVA_PRODUCTION = False
Planet.Do.HP_ICE_CONVECTION_DIAGNOSTICS = False
Planet.Do.KALOUSOVA_CONVECTION = False
```

## HP_ICE_CONVECTION_MODEL selector

`HP_ICE_CONVECTION_MODEL` is the preferred future-facing selector for HP ice
convection handling. Supported values are:

```text
"none"
"DS2001_diagnostic"
"Kalousova2018_diagnostic"
"Kalousova2018_production_experimental"
```

The default is `"none"`, so existing default behavior is unchanged.

The existing `HP_ICE_CONVECTION_DIAGNOSTICS` and `KALOUSOVA_CONVECTION` flags are
still supported for backward compatibility. When
`HP_ICE_CONVECTION_MODEL == "none"`, those flags act as shims:

```text
HP_ICE_CONVECTION_DIAGNOSTICS=False and KALOUSOVA_CONVECTION=False -> "none"
HP_ICE_CONVECTION_DIAGNOSTICS=True  and KALOUSOVA_CONVECTION=False -> "DS2001_diagnostic"
KALOUSOVA_CONVECTION=True -> "Kalousova2018_diagnostic"
```

If `HP_ICE_CONVECTION_MODEL` is explicitly set to a value other than `"none"`,
the selector wins over the legacy flags.

The only implemented selector values are diagnostic-only. The reserved
`"Kalousova2018_production_experimental"` value requires:

```python
Planet.Do.ALLOW_EXPERIMENTAL_HP_KALOUSOVA_PRODUCTION = True
```

With that gate enabled, the selector enters an experimental dry-run path. The
dry run records candidate production statuses and residuals, but it does not
accept or apply profile updates. No active production heat-stack coupling,
latent-heat feedback, transported melt, or profile feedback is activated by
this selector scaffold.

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

## Phase-local diagnostic state

Each HP ice phase diagnostic is also represented as a phase-local state record
inside `Planet.HPIceDiagnostics`. The record mirrors the legacy top-level fields
while adding bookkeeping fields such as:

```text
phaseName, phaseID, present
iTop, iBot
rTop_m, rBot_m, zTop_m, zBot_m, thickness_m
Ttop_K, Tbot_K, Ptop_MPa, Pbot_MPa, Pmid_MPa
qBot_Wm2, Qbot_W
Q_in_W, Q_out_W, q_in_Wm2, q_out_Wm2
Q_internal_W, Q_latent_W
energyResidual_W, energyResidual_frac, heatFluxResidual_Wm2
heatBookkeepingStatus
productionCandidate, productionMode, updateAccepted
candidateStatus, candidateReason
massResidual_kg, massResidual_frac
phaseBoundaryResidual_K
Tconv_K, etaConv_Pas, etaMelt_Pas
eLid_m, Dconv_m, deltaTBL_m
RaConvect, RaCrit
meltFraction, DO_HP_MELT
validityStatus, skipReason
```

This is a scaffold for future reviewable HP ice work. It does not change current
diagnostic results and it does not activate any production feedback. The helper
`_SetHPIceDiagnosticFields()` remains the synchronization point between the
legacy top-level fields and the structured `Planet.HPIceDiagnostics` record.

The `validityStatus` and `skipReason` fields are status checks only. They can
identify absent phases, too-thin layers, zero thermal contrast, nonfinite
diagnostic values, invalid geometry, subcritical regimes, and valid computed
diagnostics. These statuses do not alter model execution in this extraction.

## Diagnostic heat bookkeeping

The phase-local state includes diagnostic heat-through-stack bookkeeping fields
for future HP ice model review:

```text
Q_in_W:
  Total heat power entering the bottom of an HP ice phase block.

Q_out_W:
  Total heat power leaving the top of that phase block.

q_in_Wm2:
  Area-normalized heat flux at the bottom boundary.

q_out_Wm2:
  Area-normalized heat flux at the top boundary.
```

`Q_internal_W` and `Q_latent_W` are placeholders for future source/sink
accounting and default to `0.0`. They do not activate internal heating or
latent-heat physics.

The diagnostic residual convention is:

```text
energyResidual_W = Q_out_W - Q_in_W - Q_internal_W + Q_latent_W
energyResidual_frac = energyResidual_W / max(abs(Q_in_W), abs(Q_out_W), 1.0)
heatFluxResidual_Wm2 = q_out_Wm2 - q_in_Wm2
```

When no source or sink term is active, the diagnostic bookkeeping preserves heat
throughput by default: `Q_out_W = Q_in_W` when an incoming heat value is
available. Zero-temperature-contrast diagnostics do not force finite incoming
heat to zero. Absent or skipped phases leave heat bookkeeping unavailable rather
than inventing heat values.

These fields are inspection records only. They do not feed back into the
propagated thermal profile, phase structure, mass, gravity, viscosity, or heat
flux outputs.

## Experimental production dry run

`HP_ICE_CONVECTION_MODEL = "Kalousova2018_production_experimental"` activates a
dry-run-only scaffold when
`ALLOW_EXPERIMENTAL_HP_KALOUSOVA_PRODUCTION = True`. The dry run reuses the
Kalousova diagnostic calculations and records whether a future production update
would be considered, rejected, or left as diagnostic-only.

This path is fail-closed:

```text
updateAccepted = False unless all Ice VI acceptance criteria are explicitly met
```

Even when a synthetic candidate state can satisfy those criteria,
`updateAccepted` is recorded only in the phase-local candidate state. It does
not modify propagated profiles, heat fluxes, viscosity, mass, gravity, or phase
boundaries.

Ice VI is the only phase that can be evaluated for future production
acceptance. A supercritical, finite, positive-contrast Ice VI diagnostic is not
enough by itself: real-body acceptance remains blocked until pressure-dependent
Ice VI melting and phase-boundary checks are available. In normal body runs
without those checks, Ice VI is rejected with a machine-readable status such as:

```text
candidateStatus = "missing_melt_curve_rejected"
updateAccepted = False
```

When a tightly controlled synthetic state supplies all required acceptance
inputs, the candidate state may record:

```text
candidateStatus = "accepted_candidate_state_only"
acceptanceCriteriaPassed = True
updateAccepted = True
```

This accepted state is still candidate-only and does not alter propagated model
outputs.

Rejected Ice VI cases use machine-readable statuses such as
`"subcritical_rejected"`, `"zero_contrast_rejected"`, `"too_thin_rejected"`,
`"invalid_geometry_rejected"`, `"invalid_viscosity_rejected"`,
`"missing_melt_curve_rejected"`, `"phase_boundary_rejected"`,
`"energy_residual_rejected"`, `"mass_residual_rejected"`,
`"boundary_layer_exceeds_layer_rejected"`, and `"nonfinite_rejected"`.

Ice III and Ice V remain diagnostic and extrapolative in this mode. When they
are evaluated, their phase states are marked:

```text
candidateStatus = "diagnostic_only_extrapolative"
```

The dry run also records zero/no-op mass residuals for evaluated phase states.
It does not add latent heat, melt density, mass transport, or active Ice VI
profile propagation. The purpose is to make future production decisions
auditable before any accepted profile update is introduced.

## Ice VI production acceptance criteria

The experimental production path evaluates candidate-state acceptance criteria
for Ice VI only. The criteria are intentionally stricter than the diagnostic
Kalousova formula checks. Before an Ice VI candidate can record
`updateAccepted=True` in its phase-local state, all of the following must pass:

```text
phaseID == 6
Ice VI is present
finite geometry
positive layer thickness
layer thickness >= 1 km
positive temperature contrast
finite positive material properties
finite positive viscosity
RaConvect > RaCrit
boundary layers close within the layer
energy residual is within tolerance
mass residual is exactly zero
Q_latent_W == 0.0
Tmelt_top_K, Tmelt_mid_K, and Tmelt_bot_K are finite
phaseBoundaryResidual_K is finite and <= 0.1 K
```

The current real-body diagnostic path does not yet provide pressure-dependent
`TmeltVI(P)` or a phase-boundary residual, so Ganymede/Titan-style body runs
remain rejected even when their Ice VI Kalousova diagnostics are otherwise
supercritical.

The acceptance thresholds are:

```text
abs(energyResidual_frac) <= 1e-10
abs(energyResidual_W) <= max(1e-10 * Qscale_W, 1e-6 W)
massResidual_kg == 0.0
massResidual_frac == 0.0
abs(phaseBoundaryResidual_K) <= 0.1 K
abs(layerClosureResidual_m) <= max(1e-8 * thickness_m, 1e-3 m)
```

Ice III and Ice V remain diagnostic-only extrapolations in this mode. The
acceptance criteria do not add latent heat, melt density, mass transport, or
active Ice VI profile propagation. Active profile mutation remains future work
and requires separate validation.

## Ice VI melt-curve candidate checks

Future Ice VI production acceptance requires an EOS-backed pressure-dependent
melting curve. The candidate dry-run now records the following phase-local
fields when the experimental production selector is enabled:

```text
Tmelt_top_K
Tmelt_mid_K
Tmelt_bot_K
TmeltSource
meltCurveStatus
phaseBoundaryResidual_K
phaseBoundaryStatus
eosBoundaryContext
eosBoundaryStatus
eosBoundaryReason
candidateBoundarySource
finalProfileCoverageStatus
```

The candidate check uses `GetTfreeze()` with the existing ocean EOS phase
convention when that lookup is available and finite. The diagnostic fixed
Ice VI fallback temperature is not allowed for production acceptance. If the
lookup is unavailable, outside the EOS domain, or nonfinite, Ice VI remains
rejected with a machine-readable status such as:

```text
missing_melt_curve_rejected
melt_curve_nonfinite_rejected
outside_eos_domain_rejected
```

An `outside_eos_domain_rejected` status can refer specifically to the
provisional in-run production-candidate bounds. For example, a candidate Ice VI
bottom boundary may slightly exceed the declared EOS pressure domain even when
the finalized propagated Ice VI nodes are later found to be inside the EOS
domain. The phase-local diagnostic state records this distinction with
machine-readable fields such as:

```text
candidateBoundarySource = in_run_candidate_bounds
eosBoundaryContext = in_run_candidate
eosBoundaryStatus = in_run_candidate_boundary_outside_eos_domain
finalProfileCoverageStatus = final_profile_nodes_in_domain
```

This distinction is diagnostic only. It does not make the candidate acceptable,
and it does not change the fail-closed behavior of the experimental production
selector. Silent clamping to the EOS boundary is intentionally rejected. A
post-hoc finalized-profile candidate check may be useful as a future diagnostic,
but it is not active production modeling and does not feed back into the
propagated profile.

## Post-hoc Ice VI production candidate evaluation

Post-hoc candidate evaluation is a finalized-profile diagnostic path for Ice VI.
It is separate from the provisional in-run dry-run candidate checks. The
post-hoc path reads finalized Ice VI nodes from the propagated profile, selects
top/mid/bottom Ice VI points, and evaluates `GetTfreeze()` through the existing
EOS only for those finalized points.

Results are recorded under the phase-local diagnostics, for example:

```text
Planet.HPIceDiagnostics["VI"]["posthocProductionCandidate"]
```

The post-hoc record may report `posthoc_candidate_passed` and
`posthocUpdateAccepted = True` as candidate-state metadata when finalized Ice VI
nodes are in-domain and the residual checks pass. This is not an accepted model
update. It does not modify temperature, phase, pressure, heat flux, viscosity,
mass, gravity, or layer geometry. Active profile mutation remains future work.

The post-hoc path rejects silent EOS-domain clamping. If finalized Ice VI points
fall outside the declared EOS domain, the record remains fail-closed with a
status such as `posthoc_outside_eos_domain_rejected`. Ice III and Ice V remain
diagnostic-only extrapolations.

The phase-boundary residual is a candidate-state stability check: it records
whether the current Ice VI candidate temperatures stay within the validated
Ice VI side of the returned melting curve under the current convention. Missing
or failed checks use statuses such as:

```text
phase_boundary_ok
phase_boundary_rejected
phase_boundary_unavailable_rejected
```

These checks only populate phase-local candidate state. They do not modify the
propagated profile, heat fluxes, viscosity, mass, gravity, or phase structure.
Ganymede/Titan-style real-body acceptance remains blocked unless the
pressure-dependent melt curve and phase-boundary checks are available and pass.
Ice III and Ice V remain diagnostic-only extrapolations.

## Post-hoc margin and risk diagnostics

The post-hoc Ice VI evaluator also records margin and risk diagnostics before
any future active production update is considered. These fields are stored only
in the phase-local candidate record:

```text
posthocAllNodesInEOSDomain
posthocAllNodesIceVI
posthocEOSPressureMargin_MPa
posthocEOSTemperatureMargin_K
posthocMinPhaseBoundaryMargin_K
posthocTemperatureContrastStatus
posthocRayleighRegimeStatus
posthocThicknessStatus
posthocLayerClosureStatus
posthocViscosityStatus
posthocSensitivityRiskStatus
posthocRiskReasons
```

Unlike the initial top/mid/bottom candidate check, these diagnostics inspect all
finalized Ice VI nodes for EOS-domain coverage and Ice VI phase identity. They
also record margins to the EOS pressure/temperature bounds and the minimum
`Tfreeze(P) - T` phase-boundary margin where the EOS lookup is available.

The `posthocSensitivityRiskStatus` field summarizes whether the candidate is
nominal, close to an EOS or phase boundary, near critical Rayleigh conditions,
or high-risk/rejected. Machine-readable reasons include:

```text
near_eos_boundary
near_phase_boundary
near_critical
invalid_contrast
subcritical
too_thin
boundary_layer_exceeds_layer
invalid_viscosity
```

High-risk synthetic states such as zero or negative temperature contrast,
subcritical Rayleigh number, too-thin Ice VI geometry, boundary layers that do
not fit inside the layer, or nonfinite viscosity are rejected at the
candidate-state level. Near-boundary states can still be recorded as candidate
metadata, but they are explicitly marked for review before any future active
profile mutation.

A nominal post-hoc Ice VI candidate pass is still not active production. It does
not modify propagated temperatures, pressures, phases, heat fluxes, viscosity,
mass, gravity, or layer geometry. Ice III and Ice V remain diagnostic-only
extrapolations.

## Active Ice VI candidate profile copy

The active Ice VI scaffold can now create an isolated candidate profile copy from
finalized post-hoc Ice VI nodes. The copy is stored only in the phase-local
candidate metadata record:

```text
Planet.HPIceDiagnostics["VI"]["activeProductionCandidate"]
```

The candidate record includes independent copies of the finalized Ice VI
pressure, temperature, phase, radius/depth, and viscosity arrays where those
fields are available. It also records the source node range, node count, heat
bookkeeping summaries, rollback placeholders, and protected-field checks.

This is still not active production. The candidate copy is not written back into
the propagated profile, and the metadata keeps:

```text
candidateAppliedToProfile = False
candidateAccepted = False
```

Creating or editing the candidate copy does not modify propagated temperatures,
pressures, phases, heat fluxes, viscosity, mass, gravity, or layer geometry. Ice
VI is the only phase with an active-production candidate copy; Ice III and Ice V
remain diagnostic-only extrapolations. This scaffold prepares later conservation
residual and rollback checks before any active profile update is considered.

## Active Ice VI candidate residual evaluator

The active Ice VI scaffold can also evaluate conservation and consistency
residuals on the isolated candidate copy. The evaluator reads the candidate copy
and phase-local diagnostics, then writes residual metadata back into:

```text
Planet.HPIceDiagnostics["VI"]["activeProductionCandidate"]
```

The recorded metadata includes:

```text
candidateEnergyResidual_W
candidateEnergyResidual_frac
candidateHeatPowerResidual_W
candidateHeatFluxResidual_Wm2
candidateExpected_q_in_Wm2
candidateExpected_q_out_Wm2
candidateSphericalFluxScalingStatus
candidateMassResidual_kg
candidateMassResidual_frac
candidatePhaseBoundaryResidual_K
candidateLayerClosureResidual_m
candidateEOSPressureMargin_MPa
candidateEOSTemperatureMargin_K
candidateMinPhaseBoundaryMargin_K
candidateRaOverRaCrit
candidateResidualStatus
candidateResidualReasons
candidateResidualsPassed
```

Passing residual checks is metadata only. It does not accept or apply the
candidate copy to the propagated profile, and the scaffold keeps:

```text
candidateAppliedToProfile = False
candidateAccepted = False
```

The hard heat-conservation check is based on total heat power through the Ice VI
candidate (`Q_in_W` and `Q_out_W`). The raw area-normalized flux difference
`q_out_Wm2 - q_in_Wm2` is still recorded as `candidateHeatFluxResidual_Wm2`,
but spherical shells need not have equal area-normalized flux at their bottom
and top boundaries when total power is conserved. When radius data are
available, the evaluator records the expected spherical area-scaled fluxes and
sets `candidateSphericalFluxScalingStatus`; a direct `q` difference is not a
blocker if the finalized candidate's top flux is explained by `Q/(4*pi*r^2)`
geometry. The incoming bottom flux may come from a deeper in-run or provisional
HP boundary, so a mismatch against the finalized candidate-copy bottom radius is
recorded as metadata rather than treated as heat-power loss.

The evaluator rejects missing required inputs, nonconserved total heat power,
unexplained heat-flux geometry mismatches, nonzero mass residuals,
phase-boundary residual failures, layer-closure failures, nonpositive EOS or
phase-boundary margins, subcritical Rayleigh state, invalid temperature
contrast, and invalid viscosity. It does not compute a new thermal solution and
does not add latent heat, melt-density, mass-transport, or active Ice III/V
production behavior.

## Active Ice VI rollback and discard metadata

The active Ice VI chain records rollback/discard metadata after candidate-copy
residual evaluation. This policy is metadata-only. It does not restore a changed
Planet profile, because the active candidate chain has not modified the
propagated profile.

Rollback metadata is stored in the same phase-local record:

```text
Planet.HPIceDiagnostics["VI"]["activeProductionCandidate"]
```

The policy records:

```text
candidateDiscarded
candidateDiscardReason
rollbackStatus
rollbackReasons
rollbackPolicyVersion
rollbackRequired
rollbackApplied
rollbackReason
```

If residuals pass, the candidate remains available for later metadata-only
inspection and the policy records `rollback_not_required`. If residuals fail, the
candidate is marked discarded and cannot be used for a future active step without
an explicit future reset/rebuild helper. Re-running the rollback policy on an
already-discarded candidate preserves the original residual status and rollback
reason so failed candidates are terminal and idempotent. If protected-field
checks indicate that propagated fields changed,
the policy records `rollback_protected_fields_changed` and does not mark rollback
as applied, because that case requires investigation rather than metadata cleanup.

In all paths, the scaffold keeps:

```text
candidateAppliedToProfile = False
candidateAccepted = False
```

Thus rollback/discard still does not apply an Ice VI thermal update, does not
change phase structure, and does not alter heat flux, viscosity, mass, gravity,
or layer geometry.

## Active Ice VI thermal-update candidate reconstruction

The active Ice VI scaffold now has a no-op thermal-update reconstruction step
for the isolated candidate copy. This step creates a thermal-update namespace
inside:

```text
Planet.HPIceDiagnostics["VI"]["activeProductionCandidate"]["thermalUpdate"]
```

The reconstruction requires an existing Ice VI candidate copy with residuals
already passed. It rejects missing, discarded, or residual-failed candidates and
does not attempt to rebuild them implicitly. For a clean candidate it writes
copy-only metadata such as:

```text
candidateThermalUpdateStrategy = "no_op_reconstruction"
candidateThermalUpdateStatus
candidateThermalUpdateReasons
candidateThermalUpdateAppliedToCopy
candidateThermalUpdateAppliedToPlanet
candidateThermalUpdateAccepted
candidateUpdatedT_K
candidateUpdatedQ_W
candidateUpdatedq_Wm2
candidateUpdatedPhaseArray
candidateThermalHeatPowerResidual_W
candidateThermalRiskStatus
candidateThermalRiskReasons
```

`candidateUpdatedT_K` is an independent copy of `candidateT_K`; it is not a
view of `Planet.T_K` or the candidate source array. The no-op reconstruction
does not compute a new thermal solution. It preserves Ice VI phases, carries
total heat power through the copied nodes, and reconstructs area-normalized
candidate flux metadata with:

```text
q = Q / (4*pi*r**2)
```

This is still candidate-state only:

```text
candidateThermalUpdateAppliedToPlanet = False
candidateThermalUpdateAccepted = False
candidateAppliedToProfile = False
candidateAccepted = False
```

No propagated temperatures, pressures, phases, heat fluxes, viscosity, mass,
gravity, or layer geometry are modified. The purpose of this no-op step is to
validate the thermal-update namespace, independent-array behavior, Q/q metadata,
and future residual/rollback hooks before any science-bearing Ice VI thermal
proposal is added.

The same namespace also supports an explicitly requested
`linear_conservative_reconstruction` strategy. This is still a copied-array
candidate calculation, not a propagated model update. It requires copied Ice VI
radius, depth, pressure, temperature, phase, and thermal-conductivity metadata,
then reconstructs a trial copied temperature profile with the existing local
conductive sign convention:

```text
dT/dz = q/k
q(r) = Q / (4*pi*r**2)
```

Here `z` increases downward through the copied Ice VI nodes and `Q` is the
through-going total heat power already accepted by the candidate residual
checks. The strategy rejects missing or nonfinite conductivity, invalid copied
geometry, nonfinite updated temperatures, EOS-domain failures, non-Ice-VI
phase results, or updated temperatures that cross the pressure-dependent
freezing curve. Passing the linear reconstruction only means:

```text
candidateThermalUpdateAppliedToCopy = True
candidateThermalUpdateAppliedToPlanet = False
candidateThermalUpdateAccepted = False
```

The linear reconstruction is intended to test copied-array plumbing,
conservation metadata, and fail-closed phase-boundary checks before a future
reviewed Ice VI thermal proposal is considered.

An explicitly requested `original_boundary_reconstruction` strategy provides a
more bounded copied-profile diagnostic. It preserves the finalized copied Ice VI
top and bottom temperatures, interpolates interior copied temperatures against
`candidateZ_m`, and records Q/q metadata with the same spherical area scaling.
It is not a conduction solve and does not claim to redistribute heat. It is a
diagnostic reconstruction that asks whether a bounded copied profile tied to
the finalized Ice VI boundary temperatures remains compatible with the EOS and
pressure-dependent freezing curve.

The original-boundary strategy rejects missing or nonmonotonic copied depth
metadata instead of silently switching coordinates. A successful reconstruction
still means only:

```text
candidateThermalUpdateAppliedToCopy = True
candidateThermalUpdateAppliedToPlanet = False
candidateThermalUpdateAccepted = False
candidateAppliedToProfile = False
candidateAccepted = False
```

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
