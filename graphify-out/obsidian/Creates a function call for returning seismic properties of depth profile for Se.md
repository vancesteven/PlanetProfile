---
source_file: "PlanetProfile/Thermodynamics/HydroEOS.py"
type: "rationale"
community: "Clathrate & Custom EOS"
location: "L1177"
tags:
  - graphify/rationale
  - graphify/INFERRED
  - community/Clathrate_&_Custom_EOS
---

# Creates a function call for returning seismic properties of depth profile for Se

## Connections
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
- [[SFSeismic]] - `rationale_for` [EXTRACTED]
- [[SwConduct]] - `uses` [INFERRED]
- [[SwPhase]] - `uses` [INFERRED]
- [[SwSeismic]] - `uses` [INFERRED]

#graphify/rationale #graphify/INFERRED #community/Clathrate_&_Custom_EOS