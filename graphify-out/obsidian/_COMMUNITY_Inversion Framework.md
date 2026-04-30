---
type: community
cohesion: 0.09
members: 42
---

# Inversion Framework

**Cohesion:** 0.09 - loosely connected
**Members:** 42 nodes

## Members
- [[.__init__()_78]] - code - PlanetProfile/Utilities/defineStructs.py
- [[.assignRealPlanetModel()]] - code - PlanetProfile/Utilities/defineStructs.py
- [[.setSpaceCraft()]] - code - PlanetProfile/Utilities/defineStructs.py
- [[Assign inversion parameters for PlanetProfile runs.]] - rationale - UserConfigs/configPPinversion.py
- [[Assign inversion parameters for PlanetProfile runs._1]] - rationale - PlanetProfile/Inversion/defaultConfigInversion.py
- [[CalcGridWithinUncertainty()]] - code - PlanetProfile/Inversion/Inversion.py
- [[Calculate the grid of models within the uncertainty of the best-fit model.     F]] - rationale - PlanetProfile/Inversion/Inversion.py
- [[Compute interpolation-aware uncertainty regions for combined contours.]] - rationale - PlanetProfile/Inversion/Inversion.py
- [[Configuration settings specific to inversion calculations]] - rationale - UserConfigs/configPPinversion.py
- [[Configuration settings specific to inversion calculations_1]] - rationale - PlanetProfile/Inversion/defaultConfigInversion.py
- [[Create a higher resolution interpolated grid for smoother contours.          Arg]] - rationale - PlanetProfile/Inversion/Inversion.py
- [[Create uncertainty plots for multiple best-fit planets in a single multiplot fig]] - rationale - PlanetProfile/Inversion/Inversion.py
- [[Extract and process data for contour plotting, handling complex and 3D data.]] - rationale - PlanetProfile/Inversion/Inversion.py
- [[Fit the best-fit model within the uncertainty of the best-fit model.]] - rationale - PlanetProfile/Inversion/Inversion.py
- [[FitWithinUncertainty()]] - code - PlanetProfile/Inversion/Inversion.py
- [[Helper function to check if grid values are within uncertainty bounds.]] - rationale - PlanetProfile/Inversion/Inversion.py
- [[Interpolate z data onto higher resolution grid using linear interpolation.]] - rationale - PlanetProfile/Inversion/Inversion.py
- [[Inversion.py]] - code - PlanetProfile/Inversion/Inversion.py
- [[InversionParamsStruct]] - code - PlanetProfile/Utilities/defineStructs.py
- [[Invert for best-fit interior structureconstraints based on a set of input param]] - rationale - PlanetProfile/Inversion/Inversion.py
- [[InvertBestPlanet()]] - code - PlanetProfile/Inversion/Inversion.py
- [[InvertBestPlanetList()]] - code - PlanetProfile/Inversion/Inversion.py
- [[InvertBestPlanetMultiplot()]] - code - PlanetProfile/Inversion/Inversion.py
- [[Plot a contour from a boolean region mask.          Args         ax Matplotlib]] - rationale - PlanetProfile/Inversion/Inversion.py
- [[Plot a single uncertainty contour with proper handling of complex3D data.     U]] - rationale - PlanetProfile/Inversion/Inversion.py
- [[Plot the uncertainty of the best-fit model.          Creates a plot showing unce]] - rationale - PlanetProfile/Inversion/Inversion.py
- [[PlotUncertainty()]] - code - PlanetProfile/Inversion/Inversion.py
- [[PlotUncertaintySubplot()]] - code - PlanetProfile/Inversion/Inversion.py
- [[Set which observations to invert]] - rationale - PlanetProfile/Utilities/defineStructs.py
- [[Subplot version of PlotUncertainty - simply calls PlotUncertainty with ax parame]] - rationale - PlanetProfile/Inversion/Inversion.py
- [[_check_within_bounds()]] - code - PlanetProfile/Inversion/Inversion.py
- [[_compute_interpolated_uncertainty_regions()]] - code - PlanetProfile/Inversion/Inversion.py
- [[_create_interpolated_grid()]] - code - PlanetProfile/Inversion/Inversion.py
- [[_extract_contour_data()]] - code - PlanetProfile/Inversion/Inversion.py
- [[_interpolate_z_data()]] - code - PlanetProfile/Inversion/Inversion.py
- [[_plot_boolean_contour()]] - code - PlanetProfile/Inversion/Inversion.py
- [[_plot_uncertainty_contour()]] - code - PlanetProfile/Inversion/Inversion.py
- [[configPPinversion.py]] - code - output/exploreograms/UserConfigs/configPPinversion.py
- [[configPPinversion.py_1]] - code - UserConfigs/configPPinversion.py
- [[defaultConfigInversion.py]] - code - PlanetProfile/Inversion/defaultConfigInversion.py
- [[inversionAssign()]] - code - UserConfigs/configPPinversion.py
- [[inversionAssign()_1]] - code - PlanetProfile/Inversion/defaultConfigInversion.py

## Live Query (requires Dataview plugin)

```dataview
TABLE source_file, type FROM #community/Inversion_Framework
SORT file.name ASC
```

## Connections to other communities
- 3 edges to [[_COMMUNITY_Plotting & Visualization]]
- 2 edges to [[_COMMUNITY_Test Suite & Body Configs]]
- 1 edge to [[_COMMUNITY_MgSO4 EOS API]]
- 1 edge to [[_COMMUNITY_Gravity & Induction Config]]

## Top bridge nodes
- [[InversionParamsStruct]] - degree 23, connects to 1 community
- [[PlotUncertainty()]] - degree 10, connects to 1 community
- [[_compute_interpolated_uncertainty_regions()]] - degree 7, connects to 1 community
- [[InvertBestPlanet()]] - degree 6, connects to 1 community
- [[_interpolate_z_data()]] - degree 5, connects to 1 community