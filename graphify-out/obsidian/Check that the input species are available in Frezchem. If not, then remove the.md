---
source_file: "PlanetProfile/Thermodynamics/Reaktoro/reaktoroProps.py"
type: "rationale"
community: "Clathrate & Custom EOS"
location: "L188"
tags:
  - graphify/rationale
  - graphify/INFERRED
  - community/Clathrate_&_Custom_EOS
---

# Check that the input species are available in Frezchem. If not, then remove the

## Connections
- [[Nearest2DInterpolator]] - `uses` [INFERRED]
- [[ReturnConstantPTw]] - `uses` [INFERRED]
- [[checkSpeciesCompatibleWithFrezchem()]] - `rationale_for` [EXTRACTED]

#graphify/rationale #graphify/INFERRED #community/Clathrate_&_Custom_EOS