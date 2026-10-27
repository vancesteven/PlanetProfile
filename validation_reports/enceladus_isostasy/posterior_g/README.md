# Run (g): salinity-prior sensitivity (Fisher-checkpoint condition 2)

Launched 2026-10-27 at commit 1a53dfb6, on the same Mac production cache and the same config as run (a). Only the salinity prior changes:
- g1: w ~ U[0.1, 100] ppt, via the new prior_type `uniform-in-linear`.
- g2: log10 w ~ U[−1, 1.3].

Each run used seed 75 with n_eff 8000 and one repeat. g1 took 1896 s and g2 1648 s. Raw outputs are in `../posterior_g1/` and `../posterior_g2/`.

## Results (16 / 50 / 84 %)

| quantity | (a) log10 w ~ U[−1, 2] | g1: w ~ U ppt | g2: log10 w ≤ 1.3 |
|---|---|---|---|
| **zb (km)** | 25.87 / **26.76** / 27.74 | 26.22 / **27.17** / 28.18 | 25.78 / **26.67** / 27.60 |
| log10 w | −0.52 / 0.51 / 1.53 | 1.21 / 1.70 / 1.93 | −0.63 / 0.15 / 0.93 |
| ρ_ice | 918.2 / 924.9 / 931.8 | 918.2 / 924.9 / 931.8 | 918.2 / 924.8 / 931.7 |
| ρ_rock | 2373 / 2411 / 2453 | 2379 / 2417 / 2461 | 2371 / 2408 / 2451 |
| C2 | 0.718 / 0.814 / 0.906 | 0.730 / 0.828 / 0.918 | 0.714 / 0.811 / 0.903 |
| C/MR² | 0.3306 / 0.3324 / 0.3340 | 0.3307 / 0.3325 / 0.3340 | 0.3306 / 0.3324 / 0.3340 |
| Kish ESS | 16,213 | 16,303 | 16,024 |
| log Z | −7.096 ± 0.009 | −7.095 ± 0.006 | −7.119 ± 0.009 |

**zb within log10 w bins** (median, 16–84 %):

| bin | (a) | g1 | g2 |
|---|---|---|---|
| [−1, 0.5) | 26.61 [25.76, 27.56] | 26.56 [25.69, 27.55] | 26.63 [25.75, 27.56] |
| [0.5, 1.5) | 26.74 [25.87, 27.68] | 26.78 [25.92, 27.75] | 26.73 [25.83, 27.67] |
| [1.5, 2] | 27.26 [26.33, 28.24] | 27.34 [26.43, 28.33] | — |

## Reading

- **Size of the shift.** zb moves by **+0.41 km** under g1 and −0.10 km under g2. Both are well below the ~0.95 km 1σ half-width.
  - The g1 shift is somewhat smaller than the reviewer's expected +0.5–0.8 km.
  - The cause is that "w uniform in ppt" puts 69% of the posterior mass at w ≳ 30 ppt (log10 w ≥ 1.5).
- **Same likelihood, different weighting.** The w-conditioned zb is the same across all three priors: bin medians agree to within ≤ 0.08 km.
  - The marginal zb shift is therefore entirely the prior re-weighting the salinity bins.
  - This confirms the reporting rule: zb is constrained jointly with ocean salinity.
- **Insensitive to the w prior.** ρ_ice, ρ_rock, C2 and C/MR² barely move (C/MR² within 0.0001), and the evidence is unchanged within errors.
- **Check against prediction.** The prediction was made by importance-reweighting run (a) and recorded in STATUS before these runs. It matches both direct runs:
  - g1: zb 27.17 predicted vs 27.17 actual.
  - g2: 26.65 predicted vs 26.67 actual.

  So run (a)'s posterior already contained g1 and g2, and the reweighting shortcut is validated for later prior variants.

## Suggested wording, to update the run (a) wording

> zb = 26.8 +1.0/−0.9 km (16–84 %) under a log-uniform salinity prior (log10 w ~ U[−1, 2], w in ppt). It is 27.2 ± 1.0 km if salinity is instead uniform in ppt over 0.1–100 ppt. Conditional on salinity, zb is 26.6 km for w ≲ 3 ppt and 27.3 km for w ≳ 30 ppt. The prior dependence (≤ 0.4 km) is below the 1σ width, and salinity is not constrained by these data.

## Status

- Sampling for g1 and g2: `verified`. Both completed with ESS > 16k.
- Single seed only. Run (a)'s two-seed check showed MC error ≪ 0.02 sd, and the reweighting prediction matches each run to ≤ 0.02 km, which independently corroborates both.
- Scientific acceptance: still pending (h) and the remaining reviewer items.
