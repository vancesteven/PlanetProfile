---
source_file: "PlanetProfile/Thermodynamics/InnerEOS.py"
type: "rationale"
community: "Constant EOS Framework"
location: "L43"
tags:
  - graphify/rationale
  - graphify/INFERRED
  - community/Constant_EOS_Framework
---

# Loads in Perple_X table and creates interpolators that can be called to

## Connections
- [[ConstantEOSStruct]] - `uses` [INFERRED]
- [[EOSwrapper]] - `uses` [INFERRED]
- [[PerplexEOSStruct]] - `rationale_for` [EXTRACTED]
- [[ReturnZeros]] - `uses` [INFERRED]

#graphify/rationale #graphify/INFERRED #community/Constant_EOS_Framework