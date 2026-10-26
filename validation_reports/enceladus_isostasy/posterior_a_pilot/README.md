# Run (a) timing pilot — 2026-10-25 (Mac, repo 6c929d54)

Config: enceladus_cassini_isostasy_ocean6D_rhorock.json, seed 75, n_effective 500 / n_active 256 /
n_posterior_samples 2000 (pilot override), production cache ~/scratch/enceladus_zbw_production_cache.pkl.
Pilot only — NOT an inference result (n_eff 500).

## Plumbing (verified)
- Preflight PASS: central loglike −3.67; 1000 prior draws, 884 accepted, 116 rejected, 0 exceptions; 6.0 ms/call.
  The rejects were logged as "other (not mass closure)"; the preflight now classifies them
  (None-corner cache cell / gravity / libration) and records their (zb, log10 w). Re-measure at the full run.
- Both repeats completed (291 s, 309 s; 34,560 likelihood calls each; ESS 4158, 4339 posterior samples).
- **Same-seed reproducibility: samples and log-likelihoods bit-identical between rep1 and rep2.**
- tidal_solve_skipped = True; isostatic_derived finite on 100% of draws; mass_residual_after = 0.

## Pilot posterior (rep1) vs Fisher variant C at best fit (reviewer criterion: within ~30%)
| quantity | pilot 68% | pilot σ | Fisher C σ |
|---|---|---|---|
| zb_km | 26.73 (+1.05/−0.95) | ~1.0 | 0.79 |
| rho_rock | 2413 (+48/−40) | 43 | 37 |
| C2 | 0.81 (+0.09/−0.11) | 0.095 | 0.086 |
| C/MR² derived | 0.3323 (+0.0017/−0.0020) | ~0.0019 | 0.0016 |
| rho_ice | near-uniform | 5.6 | 5.77 (prior) |
| log10 w | near-uniform | 0.85 | 0.866 (prior) |

corr(rho_s, C2) = −0.875 (Fisher −0.886); corr(zb, log10 w) = 0.19.

zb by log10 w bin (mass fraction; zb 16/50/84): [−1,0.5) 0.49: 25.7/26.6/27.6; [0.5,1.5) 0.35: 25.8/26.7/27.7;
[1.5,2) 0.16: 26.4/27.3/28.4 — shift 0.66 km across bins, < 1σ.

Other derived (16/50/84): R_s 190.1/192.3/194.3 km; ocean thickness 30.8/33.0/35.4 km.
