# Radionuclide-inventory heating + post-hoc Perple_X mineralogy

User direction 2026-07-25: Perple_X mineralogy and real radionuclide
estimates can be POST-HOC for the current runs; a self-consistent
Perple_X fit to a given (density-rescaled) model requires invoking
porosity. PlanetProfile should optionally compute radiogenic heating
from a radionuclide inventory, adapting McDonough, Šrámek & Wipperfurth
(2020, G-cubed 21, doi:10.1029/2019GC008865;
papers/mcdonough2020radiogenic.pdf; supplement
https://data.4tu.nl/articles/_/12681371/2) and scaling to likely
radionuclide fractions in other bodies.

## Done (2026-07-25, this commit)

- `PlanetProfile/Thermodynamics/RadiogenicHeating.py`: LLR systems
  (238U, 235U, 232Th, 40K, 87Rb, 147Sm) from McDonough2020 Table 3;
  `specific_heating_Wkg(inventory, t_beforePresent_Ga)`; presets
  BSE (McDonough2020) and CI (McDonough & Sun 1995); uniform
  `scale_inventory` for leaching/depletion.
- Opt-in wiring: `Planet.Sil.radionuclideInventory` (+
  `Sil.inventoryTime_Ga`) → SetupInit derives `Sil.Qrad_Wkg`.
  Default None = existing constant behavior, no science change.
- `tests/radiogenic_heating_test.py` (`verified`): reproduces the
  paper's Eq-1 coefficients (<1%) and 19.9 TW BSE; paleo monotonicity
  + early 40K dominance.

## Next: post-hoc analysis of existing inference runs (Machine A)

1. Per-body inventory choices (science decisions, need scientific-review
   sign-off): CI×leach-factor for Europa/Titan/Callisto rock; document
   the chosen fraction + rationale per body in the run configs.
2. GUI Heating tab: optional selector [body-config constant Qrad_Wkg |
   inventory (CI, scaled)] so posteriors can be re-expressed under the
   inventory assumption post hoc — display only, no likelihood change.
3. Compare inventory-based radiogenic power against the constant-Qrad
   values currently baked into body configs; flag discrepancies.

## Next: Perple_X mineralogy refit (bigger; blocked on porosity)

- Goal: for a posterior draw with sampled/derived rho_sil, find the
  Perple_X assemblage + porosity phi that matches the draw's silicate
  density profile — restoring a mineralogy (and hence a
  radionuclide-partitioning basis) for density-rescaled interiors.
- Requires: Do.POROUS_ROCK handling in the refit path (inference caches
  are built porosity-off), choice of mantle EOS table family to search,
  and a fit criterion (rho profile? rho + VP/VS?).
- The geotherm tab already plots T(P) with provenance caption; the
  refit would add Perple_X phase boundaries along that path (roadmap
  noted in the tab caption).

## SLR note

McDonough2020 also treats short-lived radionuclides (26Al, 60Fe) —
relevant only for accretion-era thermal history, not present-day
structure; out of scope unless a thermal-evolution use case appears.
