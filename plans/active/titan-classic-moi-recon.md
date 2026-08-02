# Titan classic-MoI reachability reconnaissance

Status: **verified** reconnaissance (2026-08-02). This is evidence for the
Machine A model manager; it does not adjudicate or change a scientific
assumption. No code or `PlanetProfile/Test/` file was changed. The targeted
runs took less than one minute in aggregate, well below the 30-minute limit.

## Bottom line

The shipped self-consistent Perple_X families do not reach
`C/MR^2 = 0.341` under their ordinary Titan configurations. At the default
MgSO4-100-ppt, `Tb=255 K` hydrosphere, the best core-bearing table is
`CV3hy1wt_678_1.tab`, with a resolution-converged upper edge of about
0.3395--0.3398. No-core porous families are lower (at most 0.3212).

Two constructions can cross 0.341, but neither closes the classic-model gap
within the current defensible bounds:

1. The reduced, uniform-density Test52 no-ocean stack reaches 0.341 for many
   core densities, but requires a mean silicate density near 2525 kg/m3. The
   shipped Perple_X tables give volume-weighted grain densities of
   3148--3602 kg/m3 on the same approximate geotherm. Reducing them to 2525
   would require pervasive mean porosity of about 0.29--0.41. The configured
   Han2014 law with even `phiRockMax_frac=0.9` and
   `Pclosure_MPa=2000` retains only 0.011 volume-mean porosity on that path.
2. A warm, dense MgSO4 hydrosphere plus CV3 rock and an Fe-S core crosses the
   target numerically. The closest resolved model was 0.3412935. It selected
   the `PHydroMax_MPa=1800` boundary, however, where the MgSO4 property table
   has already ended at 800 MPa and values are nearest-clamped. It also used a
   roughly 947 km hydrosphere and a 3663 kg/m3 mean mantle, above the current
   differentiated-Titan `rho_sil <= 3500 kg/m3` support. Restricting
   `PHydroMax_MPa` to 600--1200 MPa produced no mass-compatible mantle in this
   family.

Accordingly, this reconnaissance found **no physically defensible current
configuration that closes the gap**. The numerical boundary case and the
reduced-stack solutions should not be promoted as classic-MoI solutions
without manager-led review of the hydrosphere EOS domain, silicate-density
bounds, and porosity/compaction assumptions.

## Self-consistent Perple_X screen

All rows use the default Titan MgSO4-100-ppt ocean at `Tb=255 K`, target
`C/MR^2=0.341 +/- 1e-6`, `nIceI=80`, `nSilMax=80`, and `nCore=12`. “Core” uses
`Fe-S_3D_EOS.mat`, `wFe_ppt=700`; changing `xFeS` from 0.0405 to 0.55 did not
change the reported four-decimal range because this route selects the 3-D EOS
with `wFe_ppt`, not `xFeS`.

| Silicate table | No core, porous range (`phiTop=0.1--0.9`) | Fe-S core range | Maximum |
|---|---:|---:|---:|
| `Comet_67P_CG_v7_excluding_fluid_properties.tab` | 0.3177--0.3212 | 0.3057--0.3353 | 0.3353 |
| `CI_undifferentiated_hhph_DEW17_nofluid_nomelt_685.tab` | 0.3153--0.3176 | 0.3058--0.3357 | 0.3357 |
| `CM_undifferentiated_hhph_DEW17_nofluid_nomelt_685.tab` | 0.3166--0.3183 | 0.3052--0.3348 | 0.3348 |
| `CV_undifferentiated_v4_687_DEW17_nofluid_nomelt_v2.tab` | 0.3129--0.3148 | 0.3092--0.3389 | 0.3389 |
| `CV3hy1wt_678_1.tab` | 0.3126--0.3146 | 0.3100--0.3396 | **0.3396** |

The CV3 maximum was stable to the tested numerical refinements:

| `nIceI, nSilMax, nCore, Ocean.deltaP` | C/MR2 range |
|---|---:|
| `80, 80, 12, 8 MPa` | 0.3100--0.3396 |
| `200, 200, 30, 8 MPa` | 0.3101--0.3398 |
| `200, 300, 50, 4 MPa` | 0.3096--0.3395 |

## Hydrosphere dependence

The following screen used the strongest table families from the first pass.
The listed maxima are the output of the exact classic
`FindInnerWithMoIAndEOS` path.

| Interior family | MgSO4 100 ppt, 255 K | Pure H2O ocean, 270 K | NaCl 100 ppt, 244.5 K | Pure-H2O no-ocean, 250 K |
|---|---:|---:|---:|---:|
| Comet, porous (`phiTop=0.9`) | 0.3212 | 0.3188 | 0.3187 | 0.3174 |
| Comet + Fe-S core | 0.3353 | 0.3270 | 0.3322 | 0.3309 |
| CV + Fe-S core | 0.3389 | 0.3297 | 0.3353 | 0.3339 |
| CV3 + Fe-S core | **0.3396** | 0.3302 | 0.3363 | 0.3344 |

For CV3 + Fe-S, increasing the MgSO4 loading and `Tb` raises the upper edge:

| MgSO4 (ppt) | Tb (K) | C/MR2 range or outcome |
|---:|---:|---|
| 50 | 260 | 0.3076--0.3370 |
| 100 | 255 | 0.3100--0.3397 |
| 100 | 260 | 0.3118--0.3412 (coarse); 0.3115--0.3413 (refined) |
| 150 | 250 | 0.3094--0.3394 |
| 150 | 260 | 0.3166--0.3443 |
| 200 | 250--260 | no mass-compatible silicate mantle |

With the normal `0.001` MoI tolerance, the 100-ppt, 260-K refined run selected:

| C/MR2 | Rsil (km) | Rcore (km) | mean rho_sil (kg/m3) | mean rho_core (kg/m3) | Mmodel/MTitan | seafloor P (MPa) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.3412935 | 1627.8 | 667.4 | 3662.9 | 5065.1 | 0.999984 | **1800.0 (upper boundary)** |

This is a numerical reachability example, not a defensible result: the run
warns that the MgSO4 lookup ends at 800 MPa, and its mantle density and
hydrosphere thickness sit outside the current Titan inference support.

## Direct core-radius screen

The pure-math helpers in
`PlanetProfile/Inference/structure_derivation.py` were applied to the middle
Test52 cache node (`Tb=249.9825 K`) plus its committed CMR2 discretization
offset. Silicate density was derived by mass conservation and restricted to
1800--3500 kg/m3 with `rho_sil <= rho_core`. A 1-km core-radius grid gives:

| rho_core (kg/m3) | Rcore nearest 0.341 (km) | derived rho_sil (kg/m3) | C/MR2 |
|---:|---:|---:|---:|
| 3000 | 849 | 2521.3 | 0.341000 |
| 3600 | 630 | 2524.0 | 0.340998 |
| 5150 | 461 | 2525.3 | 0.340997 |
| 6500 | 400 | 2525.6 | 0.340994 |
| 8000 | 358 | 2526.0 | 0.341004 |

At 2591 kg/m3, even a 2000-km core remains slightly high at 0.341343. The
other density families cross the target cleanly. This demonstrates that the
gap is not a moment-integration limitation; it is the incompatibility between
the low effective mantle density required by the reduced stack and the
self-consistent EOS/porosity profiles.

On the same approximate Test52 geotherm, the table-loader
`GetInnerEOS(...).fn_rho_kgm3` gives these volume-weighted grain densities:

| Table | volume-weighted grain rho (kg/m3) | uniform mean porosity needed for rho=2525 (pore rho=1000) | Han2014 `Pclosure` needed with `phiTop=0.9` on this path |
|---|---:|---:|---:|
| Comet 67P | 3148 | 0.290 | 11.2 GPa |
| CI | 3353 | 0.352 | 13.8 GPa |
| CM | 3308 | 0.339 | 13.2 GPa |
| CV | 3581 | 0.409 | 16.7 GPa |
| CV3 | 3602 | 0.414 | 17.0 GPa |

The configured `Pclosure=2 GPa` produces only 0.011 volume-mean porosity;
closure pressures of 11--17 GPa would keep substantial porosity through a
mantle whose central pressure is only about 4.8 GPa. That is why the
parametrically reachable reduced solutions do not survive the classic EOS
route.

## Reproduction

Run from the repository root in the Machine A environment:

```bash
MPLCONFIGDIR=/tmp/ppgenai-mpl \
NUMBA_CACHE_DIR=/tmp/ppgenai-numba \
KMP_DUPLICATE_LIB_OK=TRUE \
PYTHONPATH=. \
/opt/miniconda3/envs/PP/bin/python -
```

The full-profile harness imported:

```python
from copy import deepcopy
from PlanetProfile.Default.Titan.PPTitan import Planet as template
from PlanetProfile.GetConfig import Params as base
from PlanetProfile.Main import PlanetProfile

p, q = deepcopy(template), deepcopy(base)
q.CALC_NEW = True
q.NO_SAVEFILE = q.SKIP_PLOTS = True
q.SKIP_INDUCTION = q.SKIP_GRAVITY = True
q.CALC_CONDUCT = q.CALC_SEISMIC = q.CALC_VISCOSITY = False
q.PRINT_COMPLETION = q.PRELOAD_EOS = False
p.Bulk.Cmeasured = 0.341
p.Bulk.Cuncertainty = 1e-6
p.Bulk.CuncertaintyLower = p.Bulk.CuncertaintyUpper = 1e-6
p.Do.NONHYDROSTATIC = False
p.Do.ConstantProps['Inner'] = False
```

For each table row, the harness set `p.Sil.mantleEOS`, porosity/core flags,
the hydrosphere values printed above, and the stated step counts, then called
`PlanetProfile(p, q)`. When the target was not selected, the documented range
is the `Min`/`Max` emitted by `FindInnerWithMoIAndEOS`. The warm MgSO4 refined
run used `nIceI=200`, `nSilMax=300`, `nCore=50`, and `Ocean.deltaP=4 MPa`.

The direct screen loaded
`PlanetProfile/Test/mcmc_results/Titan/Test52_andrade_noocean_diff/titan_diff_noocean_structure_grid.pkl`
and its `_offsets.json` sidecar, then called
`extract_hydrosphere_layers`, `derive_silicate_density`, and `compute_cmr2`
from `PlanetProfile.Inference.structure_derivation` over
`R_core_km = np.arange(0, 2001, 1)` and the six core-density values
`[2591, 3000, 3600, 5150, 6500, 8000] kg/m3`.

