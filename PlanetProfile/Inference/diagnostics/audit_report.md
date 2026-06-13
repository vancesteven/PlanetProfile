# MCMC Run Audit Report

Generated: 2026-06-12 19:51 UTC  
Results dir: `/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai/PlanetProfile/Test/mcmc_results`  
Total pkls scanned: 9  
Unreadable: 0  

## Flag legend

| Flag | Meaning |
|---|---|
| ⚠ ESS<500 | Effective sample size below convergence threshold |
| ⚠ log-Z err>0.5 | Log-evidence MC error too large for Bayes-factor comparison |
| ⚠ seed missing | Random seed not recorded; run not reproducible |
| ⚠ pre-bugfix | git SHA absent or predates ca1b600 (2026-05-23 pocoMC unpack-swap fix); samples ok, weights/histograms may be wrong |
| ℹ smoke | n_samples ≤ 200; exploratory only, not for Bayes-factor comparison |
| ok | No issues detected |

## Runs

| Run name | Body | SHA | Seed | n_samples | ESS | log-Z | log-Z err | Walltime (h) | Sampler | Flags |
|---|---|---|---|---:|---:|---:|---:|---:|---|---|
| Test51_seawater | Europa | — | 42 | 4152 | 4152 | — | — | 0.024 | pocoMC | ℹ smoke ⚠ pre-bugfix |
| Test51_seawater | Europa | — | 42 | 4227 | 4227 | — | — | 0.017 | pocoMC | ⚠ pre-bugfix |
| Test46_allice | Titan | — | — | 4864 | — | — | — | — | legacy-dict | ⚠ seed missing ⚠ pre-bugfix |
| Test46_allice | Titan | — | — | 5120 | — | — | — | — | legacy-dict | ⚠ seed missing ⚠ pre-bugfix |
| Test46_allice | Titan | — | — | 5632 | — | — | — | — | legacy-dict | ⚠ seed missing ⚠ pre-bugfix |
| Test46_allice | Titan | — | — | 4864 | — | — | — | — | legacy-dict | ⚠ seed missing ⚠ pre-bugfix |
| Test46_allice | Titan | — | — | 5376 | — | — | — | — | legacy-dict | ⚠ seed missing ⚠ pre-bugfix |
| Test49_clathrate2km | Titan | — | — | 4530 | — | — | — | — | legacy-dict | ⚠ seed missing ⚠ pre-bugfix |
| Test49_clathrate4km | Titan | — | — | 4421 | — | — | — | — | legacy-dict | ⚠ seed missing ⚠ pre-bugfix |
