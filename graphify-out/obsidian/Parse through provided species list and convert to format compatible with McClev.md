---
source_file: "PlanetProfile/Thermodynamics/Reaktoro/reaktoroProps.py"
type: "rationale"
community: "Clathrate & Custom EOS"
location: "L1062"
tags:
  - graphify/rationale
  - graphify/INFERRED
  - community/Clathrate_&_Custom_EOS
---

# Parse through provided species list and convert to format compatible with McClev

## Connections
- [[.McClevskyIonParser()]] - `rationale_for` [EXTRACTED]
- [[Nearest2DInterpolator]] - `uses` [INFERRED]
- [[ReturnConstantPTw]] - `uses` [INFERRED]

#graphify/rationale #graphify/INFERRED #community/Clathrate_&_Custom_EOS