# MCMC Evidence Comparison Report

Generated: 2026-06-12 20:03 UTC  
Results dir: `/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai/PlanetProfile/Test/mcmc_results`  
Pairs evaluated: 4  
Pairs with log_Z available: 0  
Pairs skipped (log_Z unavailable): 4  

## Kass & Raftery (1995) classification

| |ln B₁₂| | Grade |
|---|---|
| < 1 | not worth mentioning |
| 1 – 3 | positive evidence |
| 3 – 5 | strong evidence |
| ≥ 5 | very strong evidence |

ln B₁₂ = log Z₁ − log Z₂ > 0 favours Model 1.  
σ(ln B) = √(σ_Z1² + σ_Z2²).  
Pairs with σ(ln B) > 1.0 are flagged as unresolvable.  

## Pairwise comparisons

| Model 1 | Model 2 | ln B₁₂ | σ(ln B) | K&R grade | Flag |
|---|---|---:|---:|---|---|
| Test43_andrade_arrhenius_no_ocean | Test44_maxwell_arrhenius_ocean | — | — | — | log_Z not available (pre-fix run — rerun to populate): Test43_andrade_arrhenius_no_ocean, Test44_maxwell_arrhenius_ocean |
| Test45_maxwell_hybrid_hydro | Test46_allice | — | — | — | log_Z not available (pre-fix run — rerun to populate): Test45_maxwell_hybrid_hydro, Test46_allice |
| Test43_andrade_arrhenius_no_ocean | Test46_allice | — | — | — | log_Z not available (pre-fix run — rerun to populate): Test43_andrade_arrhenius_no_ocean, Test46_allice |
| Test44_maxwell_arrhenius_ocean | Test45_maxwell_hybrid_hydro | — | — | — | log_Z not available (pre-fix run — rerun to populate): Test44_maxwell_arrhenius_ocean, Test45_maxwell_hybrid_hydro |

## Skipped pairs — detail

- **Andrade no-ocean vs Maxwell ocean** (Test43_andrade_arrhenius_no_ocean vs Test44_maxwell_arrhenius_ocean):  
  log_Z not available (pre-fix run — rerun to populate): Test43_andrade_arrhenius_no_ocean, Test44_maxwell_arrhenius_ocean
  pkl 1: `/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai/PlanetProfile/Test/mcmc_results/Titan/Test43_andrade_arrhenius_no_ocean/andrade_arrhenius_no_ocean_mcmc.pkl`
  pkl 2: `/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai/PlanetProfile/Test/mcmc_results/Titan/Test44_maxwell_arrhenius_ocean/maxwell_arrhenius_ocean_mcmc.pkl`

- **Maxwell hybrid-hydro vs all-ice** (Test45_maxwell_hybrid_hydro vs Test46_allice):  
  log_Z not available (pre-fix run — rerun to populate): Test45_maxwell_hybrid_hydro, Test46_allice
  pkl 1: `/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai/PlanetProfile/Test/mcmc_results/Titan/Test45_maxwell_hybrid_hydro/hybrid_hydro_mcmc.pkl`
  pkl 2: `/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai/PlanetProfile/Test/mcmc_results/Titan/Test46_allice/allice_yao2014_andrade_5D_results.pkl`

- **Andrade no-ocean vs all-ice** (Test43_andrade_arrhenius_no_ocean vs Test46_allice):  
  log_Z not available (pre-fix run — rerun to populate): Test43_andrade_arrhenius_no_ocean, Test46_allice
  pkl 1: `/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai/PlanetProfile/Test/mcmc_results/Titan/Test43_andrade_arrhenius_no_ocean/andrade_arrhenius_no_ocean_mcmc.pkl`
  pkl 2: `/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai/PlanetProfile/Test/mcmc_results/Titan/Test46_allice/allice_yao2014_andrade_5D_results.pkl`

- **Maxwell ocean vs Maxwell hybrid-hydro** (Test44_maxwell_arrhenius_ocean vs Test45_maxwell_hybrid_hydro):  
  log_Z not available (pre-fix run — rerun to populate): Test44_maxwell_arrhenius_ocean, Test45_maxwell_hybrid_hydro
  pkl 1: `/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai/PlanetProfile/Test/mcmc_results/Titan/Test44_maxwell_arrhenius_ocean/maxwell_arrhenius_ocean_mcmc.pkl`
  pkl 2: `/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai/PlanetProfile/Test/mcmc_results/Titan/Test45_maxwell_hybrid_hydro/hybrid_hydro_mcmc.pkl`

---

> **Note:** All comparisons require production-quality runs (n_eff ≥ 500, ESS ≥ 500).  
> Smoke runs (n_eff = 100) are listed for pipeline validation only and **must not** be used for Bayes-factor conclusions.  
> Runs flagged ⚠ pre-bugfix (before commit ca1b600, 2026-05-23) have unreliable importance weights.  
