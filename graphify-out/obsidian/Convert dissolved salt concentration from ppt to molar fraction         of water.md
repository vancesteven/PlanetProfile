---
source_file: "PlanetProfile/Thermodynamics/MgSO4/MgSO4Props.py"
type: "rationale"
community: "Clathrate & Custom EOS"
location: "L49"
tags:
  - graphify/rationale
  - graphify/INFERRED
  - community/Clathrate_&_Custom_EOS
---

# Convert dissolved salt concentration from ppt to molar fraction         of water

## Connections
- [[Massppt2molFrac()]] - `rationale_for` [EXTRACTED]
- [[Nearest2DInterpolator]] - `uses` [INFERRED]
- [[ReturnConstantPTw]] - `uses` [INFERRED]

#graphify/rationale #graphify/INFERRED #community/Clathrate_&_Custom_EOS