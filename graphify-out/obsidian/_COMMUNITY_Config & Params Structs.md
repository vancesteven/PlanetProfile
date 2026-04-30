---
type: community
cohesion: 0.12
members: 24
---

# Config & Params Structs

**Cohesion:** 0.12 - loosely connected
**Members:** 24 nodes

## Members
- [[.__init__()_76]] - code - PlanetProfile/Utilities/defineStructs.py
- [[.__init__()_71]] - code - PlanetProfile/Utilities/defineStructs.py
- [[ExploreParamsStruct]] - code - PlanetProfile/Utilities/defineStructs.py
- [[General runtime configuration parameters. Overridden by any settings contained w]] - rationale - UserConfigs/configPP.py
- [[General runtime configuration parameters. Overridden by any settings contained w_1]] - rationale - PlanetProfile/defaultConfig.py
- [[ParamsStruct]] - code - PlanetProfile/Utilities/defineStructs.py
- [[Quick test to verify exploreogram imports work correctly]] - rationale - test_exploreogram_imports.py
- [[Test that Ganymede profile can be loaded without errors.]] - rationale - test_ganymede_fix.py
- [[Test that Main module can be imported.]] - rationale - test_ganymede_fix.py
- [[Test that all required Params attributes exist.]] - rationale - test_ganymede_fix.py
- [[Test that all required imports work]] - rationale - test_exploreogram_imports.py
- [[configAssign()]] - code - UserConfigs/configPP.py
- [[configAssign()_1]] - code - PlanetProfile/defaultConfig.py
- [[configPP.py_1]] - code - output/exploreograms/UserConfigs/configPP.py
- [[configPP.py]] - code - PlanetProfileApp/configPP.py
- [[configPP.py_2]] - code - UserConfigs/configPP.py
- [[defaultConfig.py]] - code - PlanetProfile/defaultConfig.py
- [[main()]] - code - test_ganymede_fix.py
- [[test_exploreogram_imports.py]] - code - test_exploreogram_imports.py
- [[test_ganymede_fix.py]] - code - test_ganymede_fix.py
- [[test_ganymede_load()]] - code - test_ganymede_fix.py
- [[test_imports()]] - code - test_exploreogram_imports.py
- [[test_main_import()]] - code - test_ganymede_fix.py
- [[test_params_attributes()]] - code - test_ganymede_fix.py

## Live Query (requires Dataview plugin)

```dataview
TABLE source_file, type FROM #community/Config_&_Params_Structs
SORT file.name ASC
```

## Connections to other communities
- 2 edges to [[_COMMUNITY_Gravity & Induction Config]]
- 1 edge to [[_COMMUNITY_Parameter Import & Runs]]
- 1 edge to [[_COMMUNITY_Clathrate Properties]]
- 1 edge to [[_COMMUNITY_Test Suite & Body Configs]]

## Top bridge nodes
- [[.__init__()_71]] - degree 3, connects to 2 communities
- [[configAssign()]] - degree 8, connects to 1 community
- [[ExploreParamsStruct]] - degree 6, connects to 1 community
- [[ParamsStruct]] - degree 6, connects to 1 community