---
source_file: "PlanetProfile/Thermodynamics/HydroEOS.py"
type: "rationale"
community: "Clathrate & Custom EOS"
location: "L981"
tags:
  - graphify/rationale
  - graphify/INFERRED
  - community/Clathrate_&_Custom_EOS
---

# Determine if we need to load a new EOS, or if we can reuse one that's already be

## Connections
- [[CheckIfEOSLoaded()]] - `rationale_for` [EXTRACTED]
- [[ClathSeismic]] - `uses` [INFERRED]
- [[EOSLookupTableLoader]] - `uses` [INFERRED]
- [[EOSwrapper]] - `uses` [INFERRED]
- [[GetphiCalc]] - `uses` [INFERRED]
- [[MgSO4Conduct]] - `uses` [INFERRED]
- [[MgSO4PhaseLookup]] - `uses` [INFERRED]
- [[MgSO4PhaseMargules]] - `uses` [INFERRED]
- [[MgSO4Seismic]] - `uses` [INFERRED]
- [[Nearest2DInterpolator]] - `uses` [INFERRED]
- [[ReturnConstantSpecies]] - `uses` [INFERRED]
- [[ReturnZeros]] - `uses` [INFERRED]
- [[RktConduct]] - `uses` [INFERRED]
- [[RktHydroSpecies]] - `uses` [INFERRED]
- [[RktPhaseLookup]] - `uses` [INFERRED]
- [[RktPhaseOnDemand]] - `uses` [INFERRED]
- [[RktSeismic]] - `uses` [INFERRED]
- [[SwConduct]] - `uses` [INFERRED]
- [[SwPhase]] - `uses` [INFERRED]
- [[SwSeismic]] - `uses` [INFERRED]

#graphify/rationale #graphify/INFERRED #community/Clathrate_&_Custom_EOS