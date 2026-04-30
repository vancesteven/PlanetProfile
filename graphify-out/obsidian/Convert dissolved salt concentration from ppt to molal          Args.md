---
source_file: "PlanetProfile/Thermodynamics/MgSO4/MgSO4Props.py"
type: "rationale"
community: "Clathrate & Custom EOS"
location: "L33"
tags:
  - graphify/rationale
  - graphify/INFERRED
  - community/Clathrate_&_Custom_EOS
---

# Convert dissolved salt concentration from ppt to molal          Args:

## Connections
- [[Nearest2DInterpolator]] - `uses` [INFERRED]
- [[Ppt2molal()]] - `rationale_for` [EXTRACTED]
- [[ReturnConstantPTw]] - `uses` [INFERRED]

#graphify/rationale #graphify/INFERRED #community/Clathrate_&_Custom_EOS