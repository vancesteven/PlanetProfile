# Test52 Production MCMC — Summary

**Status:** implemented, unverified (pending scientific-reviewer sign-off).

## Config

`PlanetProfile/Inference/configs/test52_titan_noocean_andrade_10D.json` (n_effective=500, seed=42, unchanged).

## Convergence

| Metric | Value |
|---|---|
| n_samples | 4321 |
| ESS | 4180 |
| Acceptance | 0.4404 |
| R̂ | 1.000 |
| Wall time | ~110 s (sampler 100.6 s per metadata) |
| Structure cache | `titan_diff_noocean_structure_grid.pkl` (9 Tb points, rebuilt locally 2026-07-08; offset sidecar reproduced to <1e-6 vs committed sidecar) |

## Posterior medians ± std (all 10 params)

| Parameter | median | std | Prior bounds |
|---|---|---|---|
| alpha | 0.2905 | 0.0821 | [0.15, 0.45] |
| log10_zeta | −1.883 | 1.258 | [−3.0, 2.0] |
| log10_eta_Ih | 12.67 | 1.66 | [10, 16] |
| log10_eta_III | 12.67 | 1.65 | [10, 16] |
| log10_eta_V | 11.52 | 1.58 | [10, 16] |
| log10_eta_VI | 12.43 | 1.72 | [10, 16] |
| log10_eta_sil | 19.59 | 1.11 | [18, 22] |
| Tb_K | 250.0 | 0.563 | [249.0, 250.965] |
| R_core_km | 388 | 293 | [0, 2000] |
| rho_core_kgm3 | 2965 | 288 | [2591, 3600] |

## Posterior-predictive check vs Petricca et al. 2025 (Nature 648:558)

| Observable | Posterior mean | σ_post | Petricca μ ± σ | χ (post_mean − obs)/σ_obs |
|---|---|---|---|---|
| CMR2 | 0.342599 | 0.000519 | 0.343 ± 0.001 | −0.40σ |
| Re_k2 | 0.5705 | 0.0498 | 0.608 ± 0.048 | −0.78σ |
| \|Im_k2\| | 0.0957 | 0.0271 | 0.135 ± 0.035 | −1.12σ |

## Differentiation

- **R_core < 600 km posterior mass: 76.14%** (smoke run: 76.5% — consistent).
- Petricca weak-differentiation regime reproduced.

## Plots (this directory)

- `test52_production_corner.png` — 10D corner (977 KB)
- `test52_production_cmr2_ppc.png` — CMR2 PPC vs Petricca Gaussian (62 KB)
- `test52_production_k2_ppc.png` — Re_k2 + |Im_k2| PPC (241 KB)

## Provenance

- Result pickle: `test52_production_result.pkl` (0.6 MB, 4321 samples)
- Log: `/tmp/test52_production_mcmc.log`
- Progress JSON: `test52_production_progress.json`
- Plot script: `/tmp/test52_production_plots.py` (KDE corner style from `PlanetProfile/Test/replot_mcmc.py`)
