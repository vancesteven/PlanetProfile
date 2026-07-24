# v5 ice-thickness coverage validation

Cache: `PlanetProfile/Test/mcmc_results/Europa/Test52_seawater_v3/europa_seawater_structure_grid_v3_2d.pkl`
Grid: n_Tb=93  (Tb 259.500–271.000 K, step ~0.1250 K)
      n_w=16   (w 0.100–100.000 ppt, log-spaced)
Prior: truncated Gaussian mu=29.0, sigma=10.0, trunc [5.0,60.0] km; nominal band [19.0,39.0] km; target thin ice 10.0 km.

## 1. Monotonicity of D_iceIh(Tb) per salinity column

| w (ppt) | valid nodes | strictly decreasing? | max +delta (km) violation |
|---|---|---|---|
| 0.100 | 80 | YES | 0.0000 |
| 0.158 | 80 | YES | 0.0000 |
| 0.251 | 80 | YES | 0.0000 |
| 0.398 | 80 | YES | 0.0000 |
| 0.631 | 80 | YES | 0.0000 |
| 1.000 | 80 | YES | 0.0000 |
| 1.585 | 81 | YES | 0.0000 |
| 2.512 | 81 | YES | 0.0000 |
| 3.981 | 82 | YES | 0.0000 |
| 6.310 | 83 | YES | 0.0000 |
| 10.000 | 84 | YES | 0.0000 |
| 15.849 | 87 | YES | 0.0000 |
| 25.119 | 91 | YES | 0.0000 |
| 39.811 | 92 | YES | 0.0000 |
| 63.096 | 81 | YES | 0.0000 |
| 100.000 | 61 | YES | 0.0000 |

**All columns strictly decreasing: YES** (inversion D->Tb well-posed per column: confirmed).

## 2. Achievable ice-thickness band [D_min, D_max] per salinity

| w (ppt) | D_min (km) | D_max (km) | reaches 10 km? | covers [19,39]? | trunc [5,60] gap |
|---|---|---|---|---|---|
| 0.100 | 22.50 | 101.05 | NO | NO | thin<22.5 |
| 0.158 | 22.47 | 101.03 | NO | NO | thin<22.5 |
| 0.251 | 22.41 | 101.00 | NO | NO | thin<22.4 |
| 0.398 | 22.33 | 100.94 | NO | NO | thin<22.3 |
| 0.631 | 22.20 | 100.85 | NO | NO | thin<22.2 |
| 1.000 | 22.00 | 100.71 | NO | NO | thin<22.0 |
| 1.585 | 21.69 | 101.32 | NO | NO | thin<21.7 |
| 2.512 | 21.19 | 100.97 | NO | NO | thin<21.2 |
| 3.981 | 20.41 | 101.25 | NO | NO | thin<20.4 |
| 6.310 | 19.16 | 101.22 | NO | NO | thin<19.2 |
| 10.000 | 17.18 | 100.66 | NO | yes | thin<17.2 |
| 15.849 | 13.98 | 100.94 | NO | yes | thin<14.0 |
| 25.119 | 8.70 | 100.68 | yes | yes | thin<8.7 |
| 39.811 | 1.11 | 96.41 | yes | yes | ok |
| 63.096 | 0.78 | 85.99 | yes | yes | ok |
| 100.000 | 0.46 | 65.89 | yes | yes | ok |

## 3. Warm-edge behaviour & Tb needed for 10 km ice (informs v5 Tb bounds)

Linear fit of D vs Tb over the warmest ~4 valid nodes per column; extrapolate to D=10 km. Tb_10 above the current warm edge (271.0 K) => the v5 cache must extend the warm edge to reach 10 km ice at that salinity.

| w (ppt) | warm-edge D (km) | slope (km/K) | Tb for 10 km (K) | beyond 271 K? |
|---|---|---|---|---|
| 0.100 | 22.50 | -9.70 | 272.29 | YES |
| 0.158 | 22.47 | -9.71 | 272.28 | YES |
| 0.251 | 22.41 | -9.71 | 272.28 | YES |
| 0.398 | 22.33 | -9.71 | 272.27 | YES |
| 0.631 | 22.20 | -9.72 | 272.26 | YES |
| 1.000 | 22.00 | -9.73 | 272.23 | YES |
| 1.585 | 21.69 | -9.74 | 272.20 | YES |
| 2.512 | 21.19 | -9.77 | 272.15 | YES |
| 3.981 | 20.41 | -9.81 | 272.06 | YES |
| 6.310 | 19.16 | -9.87 | 271.93 | YES |
| 10.000 | 17.18 | -9.97 | 271.72 | YES |
| 15.849 | 13.98 | -10.14 | 271.39 | YES |
| 25.119 | 8.70 | -10.42 | 270.88 | no |
| 39.811 | 1.11 | -10.83 | 270.05 | no |
| 63.096 | 0.78 | -10.73 | 268.64 | no |
| 100.000 | 0.46 | -10.40 | 266.08 | no |

**Max Tb needed for 10 km ice across all w: ~272.29 K** (current warm edge 271.0 K). v5 warm-edge target must reach this per-column, capped at the pressure-corrected Ih melting temperature Tmelt(P_base, w).

## 4. Recommended v5 Tb grid bounds (draft)

- Cold edge: keep 259.500 K unless [19,39] km unreachable at high w (check the 'covers [19,39]?' column above; extend colder if the thick end is missing).
- Warm edge: extend to ~272.29 K so 10 km ice is reachable at low salinity. Nodes past each column's Tmelt(P_base, w) will build to None and be excluded by the inversion's edge rejection (per-column cap, not a ragged grid).
- Keep n_w=16 log-spaced; increase n_Tb to hold ~0.125 K spacing over the wider range: n_Tb ~= 104.

NOTE: 10 km ice at low salinity sits ~0.8 K below the 273.15 K bulk Ih melting point and ON the pressure-corrected melt curve; per reviewer, the D prior is truncated at ~5 km to stay off the near-melt regime — do not chase thinner ice there.
