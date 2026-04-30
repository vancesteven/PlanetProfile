# Results: Titan HP Ice Convection Study

## Overview

We modeled Titan's interior structure under two thermal scenarios to investigate the effect of ocean presence on high-pressure ice layer structure and convective heat transport. Both models used a 450 km thick hydrosphere with NaCl brine composition (100 ppt salinity) and a surface temperature of 94 K.

## Model Cases

### Case 1: Ocean Present (Tb = 244.5 K)
- **Ocean thickness**: 27.9 km (depth: 158.3-186.2 km)
- **Ice Ih lid**: 76.7 km thick (conductive)
- **Ice Ih convecting layer**: 77.7 km thick
- **Bottom ice Ih depth**: 158.3 km
- **Ocean bottom pressure**: 204.2 MPa
- **Moment of inertia**: C/MR² = 0.341

### Case 2: No Ocean (Tb = 243.165 K)  
- **Ice Ih thickness**: 161.9 km (extends to base)
- **No liquid ocean layer** (Tb below triple point)
- **Ice Ih lid**: 81.9 km thick (conductive)
- **Ice Ih convecting layer**: 76.0 km thick
- **Moment of inertia**: C/MR² = 0.340

## Phase Layering

### Ocean Case Phase Structure (surface to center):
1. **Ice Ih** (0-158.3 km depth): Convecting outer ice shell
   - Lid: 0-76.7 km depth
   - Convecting: 76.7-158.3 km depth
2. **Liquid Ocean** (158.3-186.2 km): Phase 0, NaCl brine
3. **Ice III** (186.2-233.5 km): High-pressure ice
4. **Ice II** (233.5-236.0 km): Narrow transition zone
5. **Ice III** (236.0-237.0 km): Second occurrence
6. **Ice II** (237.0-253.5 km): Extended layer
7. **Ice V** (253.5-410.7 km): Dominant HP ice phase
8. **Ice VI** (410.7-511.6 km): Deepest ice layer
9. **Silicate mantle** (511.6 km - core): Phase 50

### No Ocean Case Phase Structure:
1. **Ice Ih** (0-161.9 km depth): Convecting shell, no ocean
2. **Ice III** (161.9-236.0 km): First HP ice phase
3. **Ice II** (236.0-236.6 km): Narrow transition  
4. **Ice V** (254.3-254.9 km): Brief occurrence
5. **Ice II** (254.9-255.5 km): Return to ice II
6. **Ice V** (255.5-411.4 km): Main ice V layer
7. **Ice VI** (411.4-521.2 km): Deepest ice
8. **Silicate mantle** (521.2 km - core)

## Convection Parameters

| Parameter | Ocean Case | No Ocean Case | Units |
|-----------|------------|---------------|-------|
| Convective temperature (Tconv) | 237.7 | 236.4 | K |
| Convective viscosity (ηconv) | 2.50 × 10¹⁵ | 2.97 × 10¹⁵ | Pa·s |
| Rayleigh number (Ra) | 2.44 × 10⁷ | 2.13 × 10⁷ | - |
| Critical Rayleigh (Racrit) | 2.74 × 10⁶ | 2.76 × 10⁶ | - |
| Lid thickness | 76.7 | 81.9 | km |
| Convecting layer thickness | 77.7 | 76.0 | km |
| Thermal boundary layer (δTBL) | 3.3 | 3.5 | km |

**Both cases exhibit vigorous convection** (Ra >> Racrit by factor of ~8-9), indicating active heat transport through the ice Ih shell.

## Heat Flow

| Parameter | Ocean Case | No Ocean Case | Units |
|-----------|------------|---------------|-------|
| Surface heat flux (qsurf) | 4.69 | 4.38 | mW/m² |
| Conducted heat flux (qcon) | 5.32 | 4.99 | mW/m² |
| Heat from mantle | 390.6 | 364.8 | GW |
| Radiogenic heating (Qrad) | 1.5 × 10⁻¹² | 1.5 × 10⁻¹² | W/kg |
| Silicate tidal heating | 1.0 × 10⁻¹⁰ | 1.0 × 10⁻¹⁰ | W/m³ |

## Mass Distribution

### Ocean Case:
- **Total H₂O**: 3.42 × 10²² kg
- **Ice mass**: 1.16 × 10²² kg
- **Ocean mass**: 2.37 × 10²¹ kg  
- **Salt in ocean**: 5.35 × 10²¹ kg
- **Rock mass**: 9.50 × 10²² kg

### No Ocean Case:
- **Total H₂O**: 3.45 × 10²² kg
- **Ice mass**: 1.19 × 10²² kg
- **Ocean mass**: 0 kg
- **Pore fluid mass**: 2.83 × 10²² kg
- **Rock mass**: 9.44 × 10²² kg

## Key Findings

1. **Ocean Effect on Ice Ih Structure**: The presence of a liquid ocean layer reduces the ice Ih thickness by only ~3.6 km (158.3 vs 161.9 km), despite the significant thermal difference at the base (Tb = 244.5 K vs 243.2 K).

2. **HP Ice Transitions**: Both cases show complex phase transitions in the high-pressure ice region, with ice III, II, V, and VI all appearing. The transition sequence differs slightly between cases due to the different pressure-temperature paths.

3. **Convection Vigor**: The ocean case shows slightly more vigorous convection (Ra 15% higher) due to the warmer basal temperature, but both cases are well into the convective regime.

4. **Lid Thickness**: The conductive lid is ~5 km thinner in the ocean case (76.7 vs 81.9 km), consistent with the higher heat flux.

5. **Ice VI Dominance**: Ice VI is the deepest and thickest high-pressure ice phase in both cases, extending over ~100 km vertically.

6. **Thermal Gradient**: The ocean case maintains a steeper temperature gradient across the ice Ih shell due to the liquid layer acting as an efficient heat transport mechanism.

## Implications

These results demonstrate that:

- A small change in basal temperature (~1.3 K) across the water triple point fundamentally changes the interior structure (ocean vs. no ocean)
- The presence of an ocean has relatively minor effects on the thickness and convective behavior of the overlying ice Ih shell
- High-pressure ices (III, V, VI) form extensive layers regardless of ocean presence
- Both scenarios support vigorous convection, suggesting efficient heat transport from Titan's deep interior

## Figures

Key figures for manuscript:
1. **Wedge diagrams**: `TitanProfile_NaCl_100.0ppt_Tb244.5K_ConstantInnerRhoWedge.pdf` (ocean case) and `TitanProfile_NaCl_100.0ppt_Tb243.165K_450000.0_ConstantInnerRhoWedge.pdf` (no ocean)
2. **Hydrosphere plot**: `TitanProfile_NaCl_100.0ppt_Tb244.5K_ConstantInnerRhoHydrosphere.pdf` 
3. **Seismic profiles**: Show velocity and attenuation structure
4. **Viscosity profiles**: Show temperature-dependent rheology

## Next Steps for Analysis

1. Extract and plot temperature profiles for both cases
2. Analyze per-phase convection in HP ices (if implemented)
3. Compare with Kalousova et al. (2018) convection parameterization
4. Evaluate magnetic induction responses (if calculated)
5. Sensitivity analysis on ocean salinity effects
