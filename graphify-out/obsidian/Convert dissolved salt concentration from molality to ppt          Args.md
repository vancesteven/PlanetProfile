---
source_file: "PlanetProfile/Thermodynamics/MgSO4/MgSO4Props.py"
type: "rationale"
community: "Clathrate & Custom EOS"
location: "L18"
tags:
  - graphify/rationale
  - graphify/INFERRED
  - community/Clathrate_&_Custom_EOS
---

# Convert dissolved salt concentration from molality to ppt          Args:

## Connections
- [[Molal2ppt()]] - `rationale_for` [EXTRACTED]
- [[Nearest2DInterpolator]] - `uses` [INFERRED]
- [[ReturnConstantPTw]] - `uses` [INFERRED]

#graphify/rationale #graphify/INFERRED #community/Clathrate_&_Custom_EOS