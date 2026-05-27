# Scientific References

This file lists the key papers and textbooks underpinning PlanetProfile.
When porting a feature that implements a specific model, add its reference
here if not already present. Follow the existing format: number, full
citation, DOI or identifier, and a short note on what it is used for.

---

## Core PlanetProfile papers

**[9]** Styczinski, M. J., Vance, S. D., & Melwani Daswani, M. (2023).
*PlanetProfile: Self-consistent interior structure modeling for ocean worlds
and rocky dwarf planets in Python.*
Earth and Space Science, 10, e2022EA002748.
DOI: [10.1029/2022EA002748](https://doi.org/10.1029/2022EA002748)
> Core paper for PlanetProfile. Describes the Python framework, EOS-based
> radial profiles, parameter-space searches, magnetic induction coupling,
> seismic properties, and default models for ocean worlds.

---

## Ice rheology and grain-size dynamics

**[3]** Durham, W. B., & Stern, L. A. (2001).
*Rheological properties of water ice: Applications to satellites of the outer
planets.*
Annual Review of Earth and Planetary Sciences, 29, 295–330.
DOI: [10.1146/annurev.earth.29.1.295](https://doi.org/10.1146/annurev.earth.29.1.295)
> Foundational review for water-ice rheology, high-pressure ice phases,
> grain-size-sensitive creep, brittle-to-ductile transition, and icy
> satellite applications. Baseline reference for flow laws used in icy-world
> modeling.

**[7]** Kubo, T., Durham, W. B., Stern, L. A., & Kirby, S. H. (2006).
*Grain size-sensitive creep in ice II.*
Science, 311, 1267–1269.
DOI: [10.1126/science.1121296](https://doi.org/10.1126/science.1121296)
> Key experimental rheology paper for high-pressure ice II. Shows that
> reducing grain size from ~40 µm to ~6 µm changes the stress exponent from
> ~5 to ~2.5, making grain-size-sensitive creep relevant at low stresses in
> icy moons.

**[1]** Barr, A. C., & McKinnon, W. B. (2007).
*Convection in ice I shells and mantles with self-consistent grain size.*
Journal of Geophysical Research: Planets, 112, E02012.
DOI: [10.1029/2006JE002781](https://doi.org/10.1029/2006JE002781)
> Important for ice I shell convection and grain-size-dependent viscosity.
> Key result: grain size may evolve to 30–80 mm by dynamic
> recrystallization, while impurities may keep grains near 1–5 mm and allow
> thinner shells to convect.

---

## High-pressure ice layers and convection

**[6]** Kalousová, K., Sotin, C., Choblet, G., Tobie, G., & Grasset, O.
(2018).
*Two-phase convection in Ganymede's high-pressure ice layer: Implications for
its geological evolution.*
Icarus, 299, 133–147.
DOI: [10.1016/j.icarus.2017.07.018](https://doi.org/10.1016/j.icarus.2017.07.018)
> Central paper for HP ice convection in Ganymede. Models a two-phase
> ice–water mixture in 2D Cartesian geometry; studies melt production at the
> silicate/HP-ice interface and upward migration toward the ocean.

**[5, 4]** Kalousová, K., & Sotin, C. (2018).
*Melting in high-pressure ice layers of large ocean worlds: Implications for
volatiles transport.*
Geophysical Research Letters, 45.
DOI: [10.1029/2018GL078889](https://doi.org/10.1029/2018GL078889)
> Key for high-pressure ice layers in Ganymede and Titan. Derives scaling
> laws for volatile and water transport through HP ice; supporting information
> (no separate DOI) provides governing equations, numerical implementation,
> thermal boundary layer derivations, and simulation tables. Includes triple
> point temperatures and pressures for ice III-V-liquid and V-VI-liquid
> transitions used in PlanetProfile constants.

**[10]** Journaux, B., Brown, J. M., Pakhomova, A., Collings, I. E.,
Petitgirard, S., Espinoza, P., Ballaran, T. B., Vance, S. D., Ott, J.,
Cova, F., Garbarino, G., & Hanfland, M. (2020).
*Holistic approach for studying planetary hydrospheres: Gibbs representation of
ices and aqueous electrolytes – Application to Earth, Europa, and Ceres.*
Journal of Geophysical Research: Planets, 125, e2019JE006176.
DOI: [10.1029/2019JE006176](https://doi.org/10.1029/2019JE006176)
> Comprehensive thermodynamic framework for water ice phases and aqueous
> solutions. Provides Gibbs energy parameterizations, phase diagrams including
> triple points (ice Ih-III-liquid at 251.165 K, III-V-liquid at 254 K,
> V-VI-liquid at 272 K), and applications to ocean world interiors.

---

## Magnetic induction

**[11]** Vance, S. D., Styczinski, M. J., Bills, B. G., Cochrane, C. J.,
Soderlund, K. M., Gómez-Pérez, N., & Paty, C. (2021).
*Magnetic induction responses of Jupiter's ocean moons including effects from
adiabatic convection.*
Journal of Geophysical Research: Planets, 126, e2020JE006418.
DOI: [10.1029/2020JE006418](https://doi.org/10.1029/2020JE006418)
> Important for Europa, Ganymede, and Callisto magnetic sounding. Shows that
> depth-dependent ocean conductivity and adiabatic temperature profiles can
> significantly change induction responses compared with uniform-conductivity
> oceans.

**[2]** Cochrane, C. J., Vance, S. D., Castillo-Rogez, J. C., Styczinski,
M. J., & Liuzzo, L. (2025).
*Stronger evidence of a subsurface ocean within Callisto from a
multifrequency investigation of its induced magnetic field.*
AGU Advances, 6, e2024AV001237.
DOI: [10.1029/2024AV001237](https://doi.org/10.1029/2024AV001237)
> Relevant for magnetic induction and Callisto ocean constraints. Jointly
> considers multiple Galileo flybys and argues that Callisto's induction
> signal is better explained by a conductive ocean plus ionosphere than by
> the ionosphere alone.

---

## Europa Clipper science

**[8]** Petricca, F., Cochrane, C. J., Cascioli, G., Mazarico, E.,
Chang, S. D., Vance, S. D., Nimmo, F., & Castillo-Rogez, J. (2026).
*Characterization of Europa's interior through synthesis of Europa Clipper
data.*
The Planetary Science Journal, 7, 86.
DOI: [10.3847/PSJ/ae5225](https://doi.org/10.3847/PSJ/ae5225)
> Very relevant for Europa Clipper interpretation. Combines static gravity,
> magnetic induction, tides, rotation/orientation, and compositional
> constraints to recover Europa's ice shell thickness, ocean thickness, and
> salinity.

**[12]** Vance, S. D., Craft, K. L., Shock, E., et al. (2023).
*Investigating Europa's habitability with the Europa Clipper.*
Space Science Reviews, 219, 81.
DOI: [10.1007/s11214-023-01025-2](https://doi.org/10.1007/s11214-023-01025-2)
> Broad mission-context paper for Europa Clipper habitability science. Useful
> for connecting interior structure, composition, geology, ocean salinity,
> exchange processes, and habitability assessment.

---

## Tidal Heating and Viscoelastic Theory

**[11]** Tobie, G., Mocquet, A., & Sotin, C. (2003).
*Tidal dissipation within large icy satellites: Applications to Europa and Titan.*
Icarus, 177(2), 534–549. DOI: [10.1016/j.icarus.2005.04.006](https://doi.org/10.1016/j.icarus.2005.04.006)
> Foundational paper for viscoelastic tidal theory in icy satellites. Derives Love
> numbers for Maxwell and Andrade rheologies, and volumetric tidal heating H(r)
> from Im(k₂). Used by TidalPy backend for self-consistent tidal-thermal coupling.

**[12]** Tobie, G., Grasset, O., Lunine, J. I., Mocquet, A., & Sotin, C. (2005).
*Titan's internal structure inferred from a coupled thermal-orbital model.*
Icarus, 175(2), 496-502. DOI: [10.1016/j.icarus.2004.12.007](https://doi.org/10.1016/j.icarus.2004.12.007)
> Application of viscoelastic tidal theory to Titan. Demonstrates self-consistent
> coupling between interior structure, tidal heating, and orbital evolution. 
> Benchmark calculations for TidalPy backend validation.

**[13]** Renaud, J. P., & Henning, W. G. (2018).
*Increased tidal dissipation using advanced rheological models: Implications for Io and tidally active exoplanets.*
The Astrophysical Journal, 857(2), 98. DOI: [10.3847/1538-4357/aab784](https://doi.org/10.3847/1538-4357/aab784)
> Demonstrates superiority of Andrade rheology over Maxwell for ice deformation.
> Shows frequency-dependent tidal quality factor Q(ω) matches laboratory data.
> Justifies Andrade model use in TidalPy backend.

**[14]** Roberts, J. H., & Nimmo, F. (2008).
*Tidal heating and the long-term stability of a subsurface ocean on Enceladus.*
Icarus, 194(2), 675-689. DOI: [10.1016/j.icarus.2007.11.010](https://doi.org/10.1016/j.icarus.2007.11.010)
> Calculation of tidal heating distribution H(r) in Enceladus using viscoelastic
> propagator-matrix method. Validation reference for TidalPy per-phase heating rates.

---

## Geophysics textbook

**[15]** Turcotte, D. L., & Schubert, G. (2014).
*Geodynamics.* 3rd ed. Cambridge University Press.
ISBN: 978-1-107-00653-9 (hardback), 978-0-521-18623-0 (paperback).
> Key textbook reference for heat transfer, gravity, stress/strain, fluid
> mechanics, thermal convection, porous flow, and geophysical modeling
> fundamentals. No DOI; cite by ISBN.

---

## How to add a new reference

1. Choose the appropriate thematic section above, or create a new one.
2. Use the same format: bold number tag, full author list, italicised title,
   journal, volume/page, DOI as a hyperlink.
3. Write a one- to three-sentence note explaining what the paper is used for
   in PlanetProfile.
4. Add a `[port]`-flagged commit that updates this file alongside the feature
   that cites it.
