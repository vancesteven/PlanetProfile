---
source_file: "PlanetProfile/Thermodynamics/Reaktoro/reaktoroProps.py"
type: "rationale"
community: "Clathrate & Custom EOS"
location: "L24"
tags:
  - graphify/rationale
  - graphify/INFERRED
  - community/Clathrate_&_Custom_EOS
---

# Convert the ocean species string from g/kg of H2O to mol/kg of H2O

## Connections
- [[MolalConverter()]] - `rationale_for` [EXTRACTED]
- [[Nearest2DInterpolator]] - `uses` [INFERRED]
- [[ReturnConstantPTw]] - `uses` [INFERRED]

#graphify/rationale #graphify/INFERRED #community/Clathrate_&_Custom_EOS