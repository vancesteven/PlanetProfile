---
type: community
cohesion: 0.04
members: 71
---

# Lateral Variation & Gravity

**Cohesion:** 0.04 - loosely connected
**Members:** 71 nodes

## Members
- [[3D laterally-varying interior structure driver for PlanetProfile.  Orchestrates]] - rationale - PlanetProfile/Lateral/LateralStructure.py
- [[3D tidal heating for laterally-varying ice shell structure.  Computes spatially]] - rationale - PlanetProfile/Lateral/TidalHeating3D.py
- [[Adjust ocean floor radius uniformly to enforce mass conservation.          The p]] - rationale - PlanetProfile/Lateral/MassConservation.py
- [[AdjustForMassConservation()]] - code - PlanetProfile/Lateral/MassConservation.py
- [[Area-weighted integration of a field over the sphere.          Args]] - rationale - PlanetProfile/Lateral/SpatialGrid.py
- [[Arrhenius viscosity for ice.          eta = eta0  exp(QR  (1T - 1T_ref))]] - rationale - PlanetProfile/Lateral/TidalHeating3D.py
- [[Calculate induced gravity responses for the body and prints them to disk.]] - rationale - PlanetProfile/Gravity/Gravity.py
- [[CheckMassConservation()]] - code - PlanetProfile/Lateral/MassConservation.py
- [[ClathrateLateral.py]] - code - PlanetProfile/Lateral/ClathrateLateral.py
- [[Compute 3D tidal heating from local viscoelastic properties.          Supports M]] - rationale - PlanetProfile/Lateral/TidalHeating3D.py
- [[Compute 4pi-normalized associated Legendre function P_nm(cos(theta)).]] - rationale - PlanetProfile/Lateral/SpatialGrid.py
- [[Compute actual total mass from 3D column models and compare to target.]] - rationale - PlanetProfile/Lateral/MassConservation.py
- [[Compute effective thermal conductivity for mixed ice-clathrate.          Uses si]] - rationale - PlanetProfile/Lateral/ClathrateLateral.py
- [[Compute the degree-2 tidal strain heating pattern f(theta, phi).          For ec]] - rationale - PlanetProfile/Lateral/TidalHeating3D.py
- [[Compute volumetric dissipation factor for Andrade rheology.          The Andrade]] - rationale - PlanetProfile/Lateral/TidalHeating3D.py
- [[Compute volumetric dissipation factor for Maxwell rheology.          Uses the co]] - rationale - PlanetProfile/Lateral/TidalHeating3D.py
- [[ComputeEffectiveConductivity()]] - code - PlanetProfile/Lateral/ClathrateLateral.py
- [[ComputeTidalHeating3D()]] - code - PlanetProfile/Lateral/TidalHeating3D.py
- [[ConvergeTidalFeedback()]] - code - PlanetProfile/Lateral/TidalHeating3D.py
- [[Convert 3D ice-ocean boundary topography to SH coefficients         for MoonMag]] - rationale - PlanetProfile/Lateral/LateralStructure.py
- [[Decompose a grid field into spherical harmonic coefficients.          Uses numer]] - rationale - PlanetProfile/Lateral/SpatialGrid.py
- [[Estimate self-consistent ice thickness from lateral clathrate variation.]] - rationale - PlanetProfile/Lateral/ClathrateLateral.py
- [[EstimateIceThicknessFromClathrate()]] - code - PlanetProfile/Lateral/ClathrateLateral.py
- [[Evaluate spherical harmonic coefficients on a grid.          Uses 4pi-normalized]] - rationale - PlanetProfile/Lateral/SpatialGrid.py
- [[Extract summary fields from column models into Lateral arrays.]] - rationale - PlanetProfile/Lateral/LateralStructure.py
- [[Full 3D lateral structure pipeline.          1. Initialize grid and ice thicknes]] - rationale - PlanetProfile/Lateral/LateralStructure.py
- [[Gravity.py]] - code - PlanetProfile/Gravity/Gravity.py
- [[GravityParameters()]] - code - PlanetProfile/Gravity/Gravity.py
- [[GridToSH()]] - code - PlanetProfile/Lateral/SpatialGrid.py
- [[InitClathrateLateral()]] - code - PlanetProfile/Lateral/ClathrateLateral.py
- [[InitGrid()]] - code - PlanetProfile/Lateral/SpatialGrid.py
- [[InitLateralStructure()]] - code - PlanetProfile/Lateral/LateralStructure.py
- [[Initialize 3D lateral structure from configuration.          Reads ice thickness]] - rationale - PlanetProfile/Lateral/LateralStructure.py
- [[Initialize lateral clathrate variation from configuration.          Reads clathr]] - rationale - PlanetProfile/Lateral/ClathrateLateral.py
- [[Initialize spatial grid on Lateral substruct.          Sets theta_rad, phi_rad,]] - rationale - PlanetProfile/Lateral/SpatialGrid.py
- [[IntegrateOverSphere()]] - code - PlanetProfile/Lateral/SpatialGrid.py
- [[Iterate tidal heating - thermal structure to convergence.          Tidal heati]] - rationale - PlanetProfile/Lateral/TidalHeating3D.py
- [[Lateral clathrate variation for 3D interior models.  Specifies spatially varying]] - rationale - PlanetProfile/Lateral/ClathrateLateral.py
- [[LateralIO.py]] - code - PlanetProfile/Lateral/LateralIO.py
- [[LateralStructure.py]] - code - PlanetProfile/Lateral/LateralStructure.py
- [[Mass conservation enforcement for 3D laterally-varying interior models.  Compute]] - rationale - PlanetProfile/Lateral/MassConservation.py
- [[MassConservation.py]] - code - PlanetProfile/Lateral/MassConservation.py
- [[Reconfigure layer boundaries and gravity model into a format         usable by g]] - rationale - PlanetProfile/Gravity/Gravity.py
- [[Reload 3D lateral structure results from a pickle file.          Args]] - rationale - PlanetProfile/Lateral/LateralIO.py
- [[Reload gravity paramters from disk]] - rationale - PlanetProfile/Gravity/Gravity.py
- [[ReloadGravityParameters()]] - code - PlanetProfile/Gravity/Gravity.py
- [[ReloadLateralResults()]] - code - PlanetProfile/Lateral/LateralIO.py
- [[Run column models in parallel using multiprocessing.]] - rationale - PlanetProfile/Lateral/LateralStructure.py
- [[Run column models in serial.]] - rationale - PlanetProfile/Lateral/LateralStructure.py
- [[Run hydrosphere calculations for each grid column.          Each column gets the]] - rationale - PlanetProfile/Lateral/LateralStructure.py
- [[RunLateral3D()]] - code - PlanetProfile/Lateral/LateralStructure.py
- [[RunLateralColumns()]] - code - PlanetProfile/Lateral/LateralStructure.py
- [[SHtoGrid()]] - code - PlanetProfile/Lateral/SpatialGrid.py
- [[Save 3D lateral structure results to a pickle file.          Output file is plac]] - rationale - PlanetProfile/Lateral/LateralIO.py
- [[Save and reload 3D lateral structure fields.  Stores lateral results as pickle f]] - rationale - PlanetProfile/Lateral/LateralIO.py
- [[SaveLateralResults()]] - code - PlanetProfile/Lateral/LateralIO.py
- [[SetupGravity()]] - code - PlanetProfile/Gravity/Gravity.py
- [[Spatial grid management and spherical harmonic transforms for 3D lateral structu]] - rationale - PlanetProfile/Lateral/SpatialGrid.py
- [[SpatialGrid.py]] - code - PlanetProfile/Lateral/SpatialGrid.py
- [[TidalHeating3D.py]] - code - PlanetProfile/Lateral/TidalHeating3D.py
- [[TidalStrainPattern()]] - code - PlanetProfile/Lateral/TidalHeating3D.py
- [[UpdateAsymShapeFrom3D()]] - code - PlanetProfile/Lateral/LateralStructure.py
- [[Write the gravity parameters to disk]] - rationale - PlanetProfile/Gravity/Gravity.py
- [[WriteGravityParameters()]] - code - PlanetProfile/Gravity/Gravity.py
- [[_AndradeDissipation()]] - code - PlanetProfile/Lateral/TidalHeating3D.py
- [[_ArrheniusViscosity()]] - code - PlanetProfile/Lateral/TidalHeating3D.py
- [[_ExtractColumnSummaries()]] - code - PlanetProfile/Lateral/LateralStructure.py
- [[_MaxwellDissipation()]] - code - PlanetProfile/Lateral/TidalHeating3D.py
- [[_RunColumnsParallel()]] - code - PlanetProfile/Lateral/LateralStructure.py
- [[_RunColumnsSerial()]] - code - PlanetProfile/Lateral/LateralStructure.py
- [[_assoc_legendre_4pi()]] - code - PlanetProfile/Lateral/SpatialGrid.py

## Live Query (requires Dataview plugin)

```dataview
TABLE source_file, type FROM #community/Lateral_Variation_&_Gravity
SORT file.name ASC
```

## Connections to other communities
- 8 edges to [[_COMMUNITY_Test Suite & Body Configs]]
- 3 edges to [[_COMMUNITY_Ice Convection Models]]
- 1 edge to [[_COMMUNITY_Magnetic Induction Core]]
- 1 edge to [[_COMMUNITY_Asymmetry Functions]]
- 1 edge to [[_COMMUNITY_Community 22]]

## Top bridge nodes
- [[RunLateral3D()]] - degree 16, connects to 3 communities
- [[GravityParameters()]] - degree 10, connects to 2 communities
- [[ComputeTidalHeating3D()]] - degree 7, connects to 1 community
- [[InitLateralStructure()]] - degree 6, connects to 1 community
- [[InitClathrateLateral()]] - degree 5, connects to 1 community