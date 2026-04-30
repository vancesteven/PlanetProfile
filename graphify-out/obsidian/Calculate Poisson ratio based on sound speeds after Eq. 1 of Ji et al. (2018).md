---
source_file: "PlanetProfile/Thermodynamics/InnerEOS.py"
type: "rationale"
community: "Constant EOS Framework"
location: "L492"
tags:
  - graphify/rationale
  - graphify/INFERRED
  - community/Constant_EOS_Framework
---

# Calculate Poisson ratio based on sound speeds after Eq. 1 of Ji et al. (2018):

## Connections
- [[ConstantEOSStruct]] - `uses` [INFERRED]
- [[EOSwrapper]] - `uses` [INFERRED]
- [[ReturnZeros]] - `uses` [INFERRED]
- [[nuPoisson()]] - `rationale_for` [EXTRACTED]

#graphify/rationale #graphify/INFERRED #community/Constant_EOS_Framework