# Enceladus ocean-branch Fisher identifiability

Config: `/Users/svance/src/planetprofile-genai/PlanetProfile/Inference/configs/enceladus_cassini_isostasy_ocean6D_rhorock.json`  
Cache: `/Users/svance/scratch/enceladus_zbw_production_cache.pkl`  
Closure: `mass_rho_rock`; covariance sigmas (inflated): {'C20': 3.7367768196669176e-05, 'C22': 1.4797972834141844e-05, 'C30': 3.370840251332003e-05, 'libration_deg': 0.003010398644698074}, rho(C20,C22)=0.47  
Steps multipliers: [1.0, 2.0]; zb/log10w base step = largest node spacing (1 km, 0.07692 dex).  
Self-test: no

**Linearised (Fisher/Cramer-Rao) widths at a point; not posterior widths. Headline variant = B (Gaussian prior on libration_sys_frac only).** Four observables cannot identify six parameters: the data Fisher has rank <= 4 (variant A) and <= 5 with the sys_frac prior (B), so A/B are SINGULAR with all six free; their sigmas are inf along the null directions listed and correlations on those directions are the limiting +/-1. Finite widths in the free case come only from the prior box (variant C proxy). The rho_s-fixed columns are well posed in B.

## Summary (variant B_gauss_prior_sys (headline))

| point | rank free/fixed | corr(rho_s,C2) | flag | sig(C2) free | sig(C2) rho_s fixed | sig(zb) free | sig(zb) rho_s fixed | sig(rho_s) | sig(C/MR2) free | sig(C/MR2) rho_s fixed |
|---|---|---|---|---|---|---|---|---|---|---|
| maxL_a | 5/6, 5/5 | n/a (limit, singular F) | ok | inf | 1.221 | inf | 26.75 | inf | inf | 0.01123 |
| postmean_a | 5/6, 5/5 | n/a (limit, singular F) | ok | inf | 1.217 | inf | 24.68 | inf | inf | 0.01153 |

## Summary (variant C_B_plus_uniform_width_proxy)

| point | rank free/fixed | corr(rho_s,C2) | flag | sig(C2) free | sig(C2) rho_s fixed | sig(zb) free | sig(zb) rho_s fixed | sig(rho_s) | sig(C/MR2) free | sig(C/MR2) rho_s fixed |
|---|---|---|---|---|---|---|---|---|---|---|
| maxL_a | 6/6, 5/5 | -0.8859 | ok | 0.08957 | 0.04156 | 0.9036 | 0.9031 | 37.15 | 0.001611 | 0.0002288 |
| postmean_a | 6/6, 5/5 | -0.8844 | ok | 0.08922 | 0.04165 | 0.9277 | 0.9272 | 37.97 | 0.001618 | 0.0002367 |


## Point `maxL_a`

x = zb=26.4733, log10w=0.27578, rho_ice=920.573, rho_s=2402.21, C2=0.81224, sys=-0.00741

pred = C20=-0.00545673, C22=0.00152625, C30=0.000204309, libration_deg=0.0918694, cmr2_derived=0.332495

FD stencil: zb:central, log10w:central, rho_ice:central, rho_s:central, C2:central, sys:central

Step stability (relative difference of covariance-whitened Jacobian columns, step x1.0 vs x2.0): zb=0.00604, log10w=0.0235, rho_ice=0.00034, rho_s=0.000154, C2=2.12e-05, sys=0

||whitened J col||: zb=1.15, log10w=0.0923, rho_ice=0.067, rho_s=0.062, C2=25.7, sys=30.7

### A_likelihood_only

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=inf

sigma (rho_s fixed): zb=inf, log10w=inf, rho_ice=inf, C2=inf, sys=inf

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed inf)

rank 4/6, nullity 2; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): [{"zb": -0.466, "log10w": 0.471, "rho_ice": 0.014, "rho_s": 0.239, "C2": 0.034, "sys": -0.709}, {"zb": 0.157, "log10w": 0.212, "rho_ice": -0.731, "rho_s": -0.597, "C2": 0.1, "sys": -0.172}]

null directions (raw parameter units, unit-normalised): [{"zb": -0.0633, "log10w": 0.796, "rho_ice": 0.0316, "rho_s": 0.601, "C2": 0.000205, "sys": -0.0036}, {"zb": 0.00926, "log10w": 0.156, "rho_ice": -0.741, "rho_s": -0.653, "C2": 0.000264, "sys": -0.000381}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a      n/a      n/a      n/a      n/a
     n/a    1.000      n/a      n/a      n/a      n/a
     n/a      n/a    1.000      n/a      n/a      n/a
     n/a      n/a      n/a    1.000      n/a      n/a
     n/a      n/a      n/a      n/a    1.000      n/a
     n/a      n/a      n/a      n/a      n/a    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): -2.183e-16, 1.249e-16, 0.000295, 0.1153, 2.325, 3.56; cond = 3.56e+300

Weakest eigenvector: zb=-0.466, log10w=+0.471, rho_ice=+0.014, rho_s=+0.239, C2=+0.034, sys=-0.709

### B_gauss_prior_sys (headline)

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=0.004

sigma (rho_s fixed): zb=26.75, log10w=502, rho_ice=237.5, C2=1.221, sys=0.004

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed 0.01123)

rank 5/6, nullity 1; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): [{"zb": 0.262, "log10w": 0.095, "rho_ice": -0.714, "rho_s": -0.636, "C2": 0.089, "sys": -0.0}]

null directions (raw parameter units, unit-normalised): [{"zb": 0.0154, "log10w": 0.0695, "rho_ice": -0.718, "rho_s": -0.692, "C2": 0.000234, "sys": -8.88e-21}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a   -1.000   -1.000      n/a    0.000
     n/a    1.000      n/a      n/a      n/a    0.000
  -1.000      n/a    1.000    1.000      n/a    0.000
  -1.000      n/a    1.000    1.000      n/a    0.000
     n/a      n/a      n/a      n/a    1.000    0.000
   0.000    0.000    0.000    0.000    0.000    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): -1.103e-16, 0.0002767, 0.1146, 0.9481, 1.38, 3.557; cond = 3.56e+300

Weakest eigenvector: zb=+0.262, log10w=+0.095, rho_ice=-0.714, rho_s=-0.636, C2=+0.089, sys=-0.000

### C_B_plus_uniform_width_proxy

sigma (marginal): zb=0.9036, log10w=0.866, rho_ice=5.767, rho_s=37.15, C2=0.08957, sys=0.004

sigma (rho_s fixed): zb=0.9031, log10w=0.8659, rho_ice=5.71, C2=0.04156, sys=0.004

corr(rho_s,C2) = -0.8859; sigma(C/MR2) = 0.001611 (rho_s fixed 0.0002288)

rank 6/6, nullity 0; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): []

null directions (raw parameter units, unit-normalised): []

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    0.043   -0.108   -0.036    0.120    0.121
   0.043    1.000   -0.000    0.012    0.022   -0.000
  -0.108   -0.000    1.000    0.140    0.009    0.000
  -0.036    0.012    0.140    1.000   -0.886   -0.003
   0.120    0.022    0.009   -0.886    1.000    0.018
   0.121   -0.000    0.000   -0.003    0.018    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 0.09812, 0.7859, 0.8304, 0.9903, 1.093, 2.202; cond = 22.4

Weakest eigenvector: zb=+0.030, log10w=+0.002, rho_ice=-0.026, rho_s=-0.708, C2=+0.705, sys=+0.004


## Point `postmean_a`

x = zb=26.8016, log10w=0.50495, rho_ice=924.978, rho_s=2413.55, C2=0.81101, sys=0

pred = C20=-0.00545265, C22=0.00152379, C30=0.000206368, libration_deg=0.0912293, cmr2_derived=0.332158

FD stencil: zb:central, log10w:central, rho_ice:central, rho_s:central, C2:central, sys:central

Step stability (relative difference of covariance-whitened Jacobian columns, step x1.0 vs x2.0): zb=0.00518, log10w=0.0324, rho_ice=0.000403, rho_s=0.000152, C2=2.34e-05, sys=0

||whitened J col||: zb=1.12, log10w=0.15, rho_ice=0.0675, rho_s=0.0607, C2=25.8, sys=30.3

### A_likelihood_only

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=inf

sigma (rho_s fixed): zb=inf, log10w=inf, rho_ice=inf, C2=inf, sys=inf

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed inf)

rank 4/6, nullity 2; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): [{"zb": 0.279, "log10w": 0.035, "rho_ice": -0.692, "rho_s": -0.656, "C2": 0.087, "sys": 0.065}, {"zb": -0.415, "log10w": 0.498, "rho_ice": -0.218, "rho_s": 0.017, "C2": 0.067, "sys": -0.726}]

null directions (raw parameter units, unit-normalised): [{"zb": 0.0167, "log10w": 0.0158, "rho_ice": -0.689, "rho_s": -0.725, "C2": 0.000225, "sys": 0.000145}, {"zb": -0.0793, "log10w": 0.712, "rho_ice": -0.695, "rho_s": 0.0616, "C2": 0.000561, "sys": -0.00515}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a      n/a      n/a      n/a      n/a
     n/a    1.000      n/a      n/a      n/a      n/a
     n/a      n/a    1.000      n/a      n/a      n/a
     n/a      n/a      n/a    1.000      n/a      n/a
     n/a      n/a      n/a      n/a    1.000      n/a
     n/a      n/a      n/a      n/a      n/a    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 1.832e-16, 3.421e-16, 0.0003166, 0.1187, 2.291, 3.59; cond = 1.96e+16

Weakest eigenvector: zb=+0.279, log10w=+0.035, rho_ice=-0.692, rho_s=-0.656, C2=+0.087, sys=+0.065

### B_gauss_prior_sys (headline)

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=0.004

sigma (rho_s fixed): zb=24.68, log10w=291.6, rho_ice=204.8, C2=1.217, sys=0.004

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed 0.01153)

rank 5/6, nullity 1; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): [{"zb": 0.24, "log10w": 0.08, "rho_ice": -0.709, "rho_s": -0.651, "C2": 0.092, "sys": -0.0}]

null directions (raw parameter units, unit-normalised): [{"zb": 0.0142, "log10w": 0.0355, "rho_ice": -0.699, "rho_s": -0.714, "C2": 0.000238, "sys": -2.63e-20}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a   -1.000   -1.000      n/a    0.000
     n/a    1.000      n/a      n/a      n/a    0.000
  -1.000      n/a    1.000    1.000      n/a    0.000
  -1.000      n/a    1.000    1.000      n/a    0.000
     n/a      n/a      n/a      n/a    1.000    0.000
   0.000    0.000    0.000    0.000    0.000    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 1.825e-16, 0.0002995, 0.1181, 0.9462, 1.348, 3.588; cond = 1.97e+16

Weakest eigenvector: zb=+0.240, log10w=+0.080, rho_ice=-0.709, rho_s=-0.651, C2=+0.092, sys=-0.000

### C_B_plus_uniform_width_proxy

sigma (marginal): zb=0.9277, log10w=0.8659, rho_ice=5.766, rho_s=37.97, C2=0.08922, sys=0.004

sigma (rho_s fixed): zb=0.9272, log10w=0.8658, rho_ice=5.707, C2=0.04165, sys=0.004

corr(rho_s,C2) = -0.8844; sigma(C/MR2) = 0.001618 (rho_s fixed 0.0002367)

rank 6/6, nullity 0; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): []

null directions (raw parameter units, unit-normalised): []

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    0.067   -0.102   -0.033    0.125    0.119
   0.067    1.000   -0.000    0.018    0.037   -0.000
  -0.102   -0.000    1.000    0.143    0.008    0.000
  -0.033    0.018    0.143    1.000   -0.884   -0.003
   0.125    0.037    0.008   -0.884    1.000    0.019
   0.119   -0.000    0.000   -0.003    0.019    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 0.09797, 0.788, 0.8071, 0.9809, 1.104, 2.222; cond = 22.7

Weakest eigenvector: zb=+0.030, log10w=+0.004, rho_ice=-0.027, rho_s=-0.708, C2=+0.706, sys=+0.004

