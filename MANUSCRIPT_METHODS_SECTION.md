# Methods Section for HP Ice Convection Manuscript

## 2. Methods

### 2.1 PlanetProfile Framework

We use PlanetProfile v3.1.0+, a software framework for constructing 1D interior structure models for ocean worlds and rocky dwarf planets (Vance et al. 2018). PlanetProfile employs self-consistent thermodynamics for fluid, rock, and mineral phases, propagating models from surface to center through layered shells. The code calculates physical properties at each depth including temperature, pressure, density, thermal conductivity, seismic velocities, electrical conductivity, and gravitational acceleration.

For ice phases, PlanetProfile uses the SeaFreeze thermodynamic database (Journaux et al. 2020), which provides equations of state for ice Ih, II, III, V, and VI based on Gibbs energy parameterizations fitted to laboratory measurements and theoretical calculations. For ocean compositions, the framework supports pure water, NaCl, MgSO4, and NH₃ solutions using self-consistent EOS implementations. Silicate and core properties are calculated using PerpleX equation of state tables (Connolly 2009).

Layer propagation proceeds sequentially: surface ice layers (Ih → III → V → VI) are calculated first using ocean-ice interface temperature as a boundary condition, followed by ocean properties (if present), then HP ice underplating, silicate mantle, and metallic core. At each layer boundary, continuity of temperature and heat flux is enforced. Convection is parameterized using scaling laws rather than solving full Navier-Stokes equations, allowing rapid exploration of parameter space while maintaining physical fidelity.

### 2.2 High-Pressure Ice Convection: Kalousová & Sotin (2018)

We implement the two-phase convection parameterization of Kalousová & Sotin (2018) for HP ice layers (III, V, VI). This model accounts for partial melt formation at ice-silicate and ice-ice phase boundaries when convection becomes sufficiently vigorous. The key insight is that strong convection drives temperatures toward the solidus (melting curve), creating "temperate" layers where ice and liquid coexist in thermodynamic equilibrium.

The model is governed by four relationships derived from boundary layer theory and validated against numerical convection simulations:

**Modified Rayleigh number (Ra\*):**
The far-field Rayleigh number characterizes convective vigor:

$$\text{Ra}^* = \frac{\alpha \rho^2 g c_p \Delta T_c H_c^3}{k \mu_0}$$

where α is thermal expansivity [K⁻¹], ρ is density [kg m⁻³], g is gravitational acceleration [m s⁻²], cₚ is specific heat [J kg⁻¹ K⁻¹], ΔTc = Tb - Tinterior is temperature contrast [K], Hc is layer thickness [m], k is thermal conductivity [W m⁻¹ K⁻¹], and μ₀ is reference viscosity at the melting temperature [Pa s]. This formulation differs from the classical Rayleigh number by using the reference viscosity at the melting point rather than the bottom boundary temperature, reflecting the expectation that vigorous two-phase convection maintains near-solidus conditions throughout the layer.

**Critical Rayleigh number (Ra\*c):**
Temperate layer formation occurs when Ra\* exceeds a critical threshold that depends on heat flux from the silicate layer:

$$\text{Ra}_c^* = 1.9965 \times 10^4 \times (q_s[\text{mW m}^{-2}])^{3.690}$$

where qs is the heat flux from the silicate mantle in mW m⁻². This strong power-law dependence (exponent ~ 3.7) reflects the nonlinear feedback between heat transport efficiency and melt formation: higher silicate heat flux increases convective vigor, which drives temperatures closer to the melting point, which reduces viscosity, which further enhances convection.

**Temperate layer thickness (Ht):**
When Ra\* > Ra\*c, a temperate (partially molten) layer forms at the ocean-ice or ice-ice interface with thickness:

$$H_t[\text{km}] = (1.45 \times 10^{-4} \times q_s[\text{mW m}^{-2}] + 0.015) \times \mu_0[\text{Pa s}]^{0.21}$$

This relationship shows that temperate layer thickness increases linearly with silicate heat flux and scales as a weak power of reference viscosity (exponent 0.21). For typical Titan conditions (qs ~ 5–10 mW m⁻², μ₀ ~ 10¹¹–10¹² Pa s), Ht ranges from a few kilometers to tens of kilometers.

**Thermal boundary layer thickness (δ):**
The lower thermal boundary layer, which connects the convecting interior to the hot silicate boundary, scales as:

$$\delta^* = 2.746 \times (\text{Ra}^*)^{-0.271}$$
$$\delta = \delta^* \times H_c$$

This represents a ~27% decrease in boundary layer thickness per order-of-magnitude increase in Rayleigh number, consistent with turbulent thermal convection scaling theory.

**Implementation in PlanetProfile:**
The Kalousová model is implemented in `ConvectionKalousova2018()` (ThermalProfiles.py:70-214) and dispatched to HP ice layers via configuration flag `Planet.Do.KALOUSOVA_CONVECTION = True`. The function signature maintains compatibility with the existing Deschamps & Sotin (2001) convection model, returning eight values: convecting interior temperature, reference viscosity, temperate layer thickness, convecting layer thickness, thermal boundary layer thickness, heat flux, Ra\*, and Ra\*c.

Physical properties (ρ, cₚ, α, k) are evaluated at mid-layer pressure and melting temperature using the SeaFreeze database. The convecting interior temperature follows the melting curve T = Tmelt(P) at each depth rather than an adiabatic profile, reflecting the assumption that two-phase convection maintains ice at the solidus throughout the convecting region.

When Ra\* > Ra\*c (supercritical regime), the layer structure consists of:
1. **Temperate layer** (thickness Ht): Partial melt at ocean-ice interface
2. **Convecting interior** (thickness Hc - Ht - δ): Solid ice at melting curve
3. **Thermal boundary layer** (thickness δ): Conductive gradient from Tmelt to Tb

When Ra\* < Ra\*c (subcritical regime), no temperate layer forms and the entire layer convects as single-phase solid ice.

### 2.3 Per-Phase Convection Controls

PlanetProfile allows independent control of convection in each ice phase via configuration flags:
- `Planet.Do.NO_ICE_CONVECTION_Ih = True`: Suppress ice Ih convection
- `Planet.Do.NO_ICE_CONVECTION_III = True`: Suppress ice III convection
- `Planet.Do.NO_ICE_CONVECTION_V = True`: Suppress ice V convection
- `Planet.Do.NO_ICE_CONVECTION_VI = True`: Suppress ice VI convection

This granular control enables investigation of scenarios where specific phases are conductive (e.g., stagnant lid in ice Ih) while others convect (e.g., mobile HP ice). When convection is suppressed, the layer temperature profile is calculated conductively using Fourier's law with depth-dependent thermal conductivity from SeaFreeze.

For ice Ih, the Deschamps & Sotin (2001) stagnant-lid convection model is always used regardless of the Kalousová flag, as two-phase convection is not expected in ice Ih (melting point is at lower pressure). For HP ice layers (III, V, VI), the choice between Deschamps & Sotin and Kalousová models is controlled by `Planet.Do.KALOUSOVA_CONVECTION`.

Triple-point temperatures for HP ice phase boundaries are taken from Journaux et al. (2020):
- Ice III-V-liquid: T = 254 K, P = 350 MPa
- Ice V-VI-liquid: T = 272 K, P = 632 MPa

When the Kalousová model predicts temperate layer formation (Ra\* > Ra\*c), the flag `Planet.DO_HP_MELT = True` is set and melt fractions are stored per phase in `Planet.meltFractionIII/V/VI`. In the current implementation, a nominal melt fraction of 0.01 (1%) is assigned; future work will implement full two-phase flow modeling to calculate melt fractions from enthalpy balance.

### 2.4 Rheological Parameters

Reference viscosities at the melting point (μ₀) are phase-specific and stored in `Constants.etaMelt_Pas`:

| Phase | μ₀ [Pa s] | Source |
|-------|-----------|--------|
| Ice Ih | 10¹⁴ | Goldsby & Kohlstedt (2001) |
| Ice III | 10¹³ | Durham et al. (1997) |
| Ice V | 10¹³ | Durham et al. (1997) |
| Ice VI | 10¹³ | Durham et al. (1997) |

These values represent crystalline ice without fluid inclusions. For models exploring fluid-bearing HP ice (motivated by Petricca et al. 2025), we reduce μ₀ by factors of 100–1000 to achieve effective viscosities of 10¹⁰–10¹² Pa s. This reduction is consistent with laboratory measurements showing that small amounts of liquid (~0.01–0.1% by volume) can reduce ice viscosity by 2–4 orders of magnitude (Durham et al. 1996).

Activation energies for temperature-dependent viscosity (used in ice Ih convection with Deschamps & Sotin model) are phase-specific:
- Ice Ih: Eₐ = 60 kJ mol⁻¹ (solid-state diffusion creep)
- Ice III/V/VI: Eₐ = 90–110 kJ mol⁻¹ (estimated from high-pressure rheology experiments)

The current implementation uses a reference viscosity approach where μ₀ encodes both the intrinsic ice rheology and the effect of fluid inclusions. Future work will implement explicit composite rheology combining diffusion and dislocation creep mechanisms (see Kalousová et al. 2018 supplementary material).

### 2.5 Titan Model Configuration

We model Titan's interior with two end-member scenarios:

**Scenario 1: Ocean present**
- Total hydrosphere thickness: 450 km
- Ocean composition: NaCl 100 ppt
- Ocean layer thickness: ~28 km (dynamically determined from phase boundaries)
- Ice layers: Ih (surface) → Ocean → III → V → VI → Silicates
- Boundary temperature (ocean-ice III): Tb = 244.5 K (from phase diagram)

**Scenario 2: No ocean**
- Total hydrosphere thickness: 450 km
- No liquid water layer
- Ice layers: Ih (surface) → III → V → VI → Silicates
- Boundary temperature (ice Ih-III): Tb = 243.2 K (from phase diagram)

Both scenarios use:
- Surface temperature: Ts = 94 K (Cassini observations)
- Surface pressure: Ps = 0.15 MPa
- Moment of inertia constraint: C/MR² = 0.33–0.341 (Iess et al. 2012)
- Clathrate hydrate in upper ice Ih layer: 5% by volume (enables ice Ih convection, Choukroun & Sotin 2012)

The silicate core and rocky mantle are modeled using PerpleX equation of state tables for hydrated silicates at appropriate pressure-temperature conditions. Core radius and density are adjusted to satisfy the measured C/MR² constraint.

HP ice viscosities are varied parametrically in the range 10¹⁰–10¹⁴ Pa s to explore the fluid-bearing vs crystalline end-members and assess consistency with Petricca et al. (2025) tidal response constraints.
