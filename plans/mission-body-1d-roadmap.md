# Mission–body 1D artifact roadmap

User direction 2026-07-19: current GUI models are 1D (spherically
symmetric radial interior structure) and are labeled so in the GUI.
Plan ONE 1D amortized model per mission–body pair. Parameter counts
(7D/8D/10D) refer to sampled-parameter dimensionality, NOT structure
dimensionality — the GUI model-assumptions text disambiguates.

## Current (labeled 1D in GUI)

| Mission–body | Artifact | Status |
|---|---|---|
| Cassini–Titan | titan_andrade_noocean_posterior.pt (Test50 8D) | deployed |
| Galileo–Europa | europa_seawater_andrade_posterior_1m.pt | RETIRING → Galileo v1.1 (honest observables; see v4 plan) |
| Clipper–Europa | europa_seawater_andrade_clipper_v2.pt | RETIRING → v4 geodesy (Mazarico projections) |

## Planned (one per mission–body; priority per user)

| Mission–body | Key measured/projected constraints | Notes |
|---|---|---|
| Galileo–Ganymede | CMR2 (Anderson 1996 hydrostatic), no k2 | Test/ has ganymede_pureh2o structure grid (C2 era) — cache exists |
| Galileo–Callisto | CMR2 (Anderson 2001; possibly non-hydrostatic!), no k2 | Test/ has callisto nacl/mgso4 grids; non-hydrostaticity literature strong for Callisto — v4-style C20/C22 nuisance treatment especially apt |
| Cassini–Enceladus | libration (Thomas 2016), gravity (Iess 2014), no k2 measured | no PP cache yet; small-body regime |
| JUICE–Callisto | k2 ± 0.06 (21 flybys); Im k2 unconstrained unless Q < ~10; obliquity 15 arcmin | Van Hoolst et al. 2024 |
| JUICE–Ganymede | k2 ± 1e-4 (REAL AND IMAGINARY, orbital phase GCO500); h2 to a few %; obliquity 0.2 arcsec; libration 0.4–0.8 arcsec | the flagship: k2+h2 joint inversion → ice-shell thickness ~few-km class; Im k2 → shell viscosity/thermal state |
| JUICE–Europa | 2 close flybys — insufficient for h2/rotation; modest k2 at best | complements Clipper v4, not a substitute |

Projected-accuracy source: Van Hoolst, T., Tobie, G., Vallat, C., et
al. (2024), "Geophysical Characterization of the Interiors of
Ganymede, Callisto and Europa by ESA's JUpiter ICy moons Explorer",
SSRv 220, 54, doi:10.1007/s11214-024-01085-y. The SSRv JUICE special
issue (https://link.springer.com/collections/dficcaejhd) has
per-instrument sensitivity papers (3GM, GALA, JANUS) — consult when
freezing each config's sigmas, per-mission normalization/convention
gates as in the v4 plan (learned the hard way: Titan sigmas leaked
into Europa configs).

Per-artifact discipline (unchanged): mission-appropriate observables
ONLY (no cross-body sigma reuse); hypothetical channels labeled as
such; config metadata cites the sensitivity source; full gate suite +
reference MCMC; INDEX row; user ratification before slot deploy.

## Future work: 3D models + mission visualization

For eventual 3D (laterally varying) models, missionwidget.com is the
candidate visualization layer — both for 3D interior rendering and for
spacecraft trajectory/measurement-geometry context (which flyby
measures what, where). Same dependency class as the Clipper-trajectory
simulation needed for higher-degree gravity and directional Bind
inversion (v4 plan future work). Not scheduled.
