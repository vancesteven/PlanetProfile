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
| lw-0.5 | 5/6, 5/5 | n/a (limit, singular F) | ok | inf | 1.144 | inf | 28.01 | inf | inf | 0.01186 |
| lw1.0 | 5/6, 5/5 | n/a (limit, singular F) | ok | inf | 1.196 | inf | 35.33 | inf | inf | 0.01215 |
| lw1.5 | 5/6, 5/5 | n/a (limit, singular F) | ok | inf | 1.288 | inf | 55.06 | inf | inf | 0.01245 |
| lw1.9 | 5/6, 5/5 | n/a (limit, singular F) | ok | inf | 1.405 | inf | 101.5 | inf | inf | 0.0121 |

## Summary (variant C_B_plus_uniform_width_proxy)

| point | rank free/fixed | corr(rho_s,C2) | flag | sig(C2) free | sig(C2) rho_s fixed | sig(zb) free | sig(zb) rho_s fixed | sig(rho_s) | sig(C/MR2) free | sig(C/MR2) rho_s fixed |
|---|---|---|---|---|---|---|---|---|---|---|
| lw-0.5 | 6/6, 5/5 | -0.8872 | ok | 0.08627 | 0.0398 | 0.7867 | 0.7863 | 37.41 | 0.001619 | 0.0002173 |
| lw1.0 | 6/6, 5/5 | -0.872 | ok | 0.08726 | 0.04271 | 0.8089 | 0.8088 | 37.76 | 0.001621 | 0.0002984 |
| lw1.5 | 6/6, 5/5 | -0.7523 | ok | 0.09084 | 0.05985 | 1.035 | 1.023 | 39.67 | 0.001629 | 0.0006553 |
| lw1.9 | 6/6, 5/5 | -0.3883 | ok | 0.09873 | 0.09099 | 1.76 | 1.514 | 48.44 | 0.001636 | 0.001197 |


## Point `lw-0.5`

x = zb=25, log10w=-0.5, rho_ice=925, rho_s=2400, C2=0.75, sys=0

pred = C20=-0.00549857, C22=0.00153182, C30=0.000213029, libration_deg=0.0978242, cmr2_derived=0.332839

FD stencil: zb:central, log10w:central, rho_ice:central, rho_s:central, C2:central, sys:central

Step stability (relative difference of covariance-whitened Jacobian columns, step x1.0 vs x2.0): zb=0.00639, log10w=0.0376, rho_ice=0.000409, rho_s=0.000154, C2=2.45e-05, sys=7.08e-14

||whitened J col||: zb=1.31, log10w=0.0149, rho_ice=0.0665, rho_s=0.0613, C2=26.5, sys=32.5

### A_likelihood_only

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=inf

sigma (rho_s fixed): zb=inf, log10w=inf, rho_ice=inf, C2=inf, sys=inf

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed inf)

rank 4/6, nullity 2; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): [{"zb": -0.084, "log10w": -0.27, "rho_ice": 0.73, "rho_s": 0.54, "C2": -0.065, "sys": 0.301}, {"zb": -0.497, "log10w": 0.404, "rho_ice": 0.133, "rho_s": 0.326, "C2": 0.016, "sys": -0.682}]

null directions (raw parameter units, unit-normalised): [{"zb": -0.0028, "log10w": -0.79, "rho_ice": 0.478, "rho_s": 0.384, "C2": -0.000106, "sys": 0.000403}, {"zb": -0.0137, "log10w": 0.979, "rho_ice": 0.0722, "rho_s": 0.191, "C2": 2.17e-05, "sys": -0.000756}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a      n/a      n/a      n/a      n/a
     n/a    1.000      n/a      n/a      n/a      n/a
     n/a      n/a    1.000      n/a      n/a      n/a
     n/a      n/a      n/a    1.000      n/a      n/a
     n/a      n/a      n/a      n/a    1.000      n/a
     n/a      n/a      n/a      n/a      n/a    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): -5.78e-16, 5.434e-16, 0.0002386, 0.11, 2.392, 3.498; cond = 3.5e+300

Weakest eigenvector: zb=-0.084, log10w=-0.270, rho_ice=+0.730, rho_s=+0.540, C2=-0.065, sys=+0.301

### B_gauss_prior_sys (headline)

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=0.004

sigma (rho_s fixed): zb=28.01, log10w=3474, rho_ice=302.3, C2=1.144, sys=0.004

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed 0.01186)

rank 5/6, nullity 1; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): [{"zb": 0.278, "log10w": 0.084, "rho_ice": -0.722, "rho_s": -0.626, "C2": 0.053, "sys": 0.0}]

null directions (raw parameter units, unit-normalised): [{"zb": 0.0133, "log10w": 0.354, "rho_ice": -0.681, "rho_s": -0.641, "C2": 0.000125, "sys": 7.76e-21}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a   -1.000   -1.000      n/a    0.000
     n/a    1.000      n/a      n/a      n/a    0.000
  -1.000      n/a    1.000    1.000      n/a    0.000
  -1.000      n/a    1.000    1.000      n/a    0.000
     n/a      n/a      n/a      n/a    1.000    0.000
   0.000    0.000    0.000    0.000    0.000    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 1.617e-16, 0.000228, 0.1097, 0.9486, 1.451, 3.491; cond = 2.16e+16

Weakest eigenvector: zb=+0.278, log10w=+0.084, rho_ice=-0.722, rho_s=-0.626, C2=+0.053, sys=+0.000

### C_B_plus_uniform_width_proxy

sigma (marginal): zb=0.7867, log10w=0.866, rho_ice=5.767, rho_s=37.41, C2=0.08627, sys=0.004

sigma (rho_s fixed): zb=0.7863, log10w=0.866, rho_ice=5.713, C2=0.0398, sys=0.004

corr(rho_s,C2) = -0.8872; sigma(C/MR2) = 0.001619 (rho_s fixed 0.0002173)

rank 6/6, nullity 0; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): []

null directions (raw parameter units, unit-normalised): []

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    0.007   -0.117   -0.034    0.092    0.127
   0.007    1.000   -0.000    0.002    0.003   -0.000
  -0.117   -0.000    1.000    0.136    0.013    0.000
  -0.034    0.002    0.136    1.000   -0.887   -0.003
   0.092    0.003    0.013   -0.887    1.000    0.016
   0.127   -0.000    0.000   -0.003    0.016    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 0.09926, 0.7776, 0.8624, 0.9997, 1.103, 2.158; cond = 21.7

Weakest eigenvector: zb=+0.024, log10w=+0.000, rho_ice=-0.025, rho_s=-0.709, C2=+0.704, sys=+0.004


## Point `lw1.0`

x = zb=25, log10w=1, rho_ice=925, rho_s=2400, C2=0.75, sys=0

pred = C20=-0.00550291, C22=0.00153263, C30=0.000213843, libration_deg=0.0982333, cmr2_derived=0.332944

FD stencil: zb:central, log10w:central, rho_ice:central, rho_s:central, C2:central, sys:central

Step stability (relative difference of covariance-whitened Jacobian columns, step x1.0 vs x2.0): zb=0.00641, log10w=0.00503, rho_ice=0.00027, rho_s=0.000156, C2=2.01e-05, sys=7.08e-14

||whitened J col||: zb=1.32, log10w=0.531, rho_ice=0.0677, rho_s=0.061, C2=26.3, sys=32.6

### A_likelihood_only

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=inf

sigma (rho_s fixed): zb=inf, log10w=inf, rho_ice=inf, C2=inf, sys=inf

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed inf)

rank 4/6, nullity 2; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): [{"zb": -0.181, "log10w": -0.326, "rho_ice": 0.74, "rho_s": 0.494, "C2": -0.069, "sys": 0.254}, {"zb": 0.45, "log10w": -0.456, "rho_ice": -0.1, "rho_s": -0.339, "C2": -0.022, "sys": 0.682}]

null directions (raw parameter units, unit-normalised): [{"zb": -0.0101, "log10w": -0.045, "rho_ice": 0.803, "rho_s": 0.594, "C2": -0.000193, "sys": 0.000571}, {"zb": 0.0587, "log10w": -0.147, "rho_ice": -0.255, "rho_s": -0.954, "C2": -0.000143, "sys": 0.00359}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a      n/a      n/a      n/a      n/a
     n/a    1.000      n/a      n/a      n/a      n/a
     n/a      n/a    1.000      n/a      n/a      n/a
     n/a      n/a      n/a    1.000      n/a      n/a
     n/a      n/a      n/a      n/a    1.000      n/a
     n/a      n/a      n/a      n/a      n/a    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): -1.638e-16, 2.009e-16, 0.0001968, 0.1048, 2.418, 3.477; cond = 3.48e+300

Weakest eigenvector: zb=-0.181, log10w=-0.326, rho_ice=+0.740, rho_s=+0.494, C2=-0.069, sys=+0.254

### B_gauss_prior_sys (headline)

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=0.004

sigma (rho_s fixed): zb=35.33, log10w=116.6, rho_ice=418.7, C2=1.196, sys=0.004

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed 0.01215)

rank 5/6, nullity 1; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): [{"zb": 0.327, "log10w": 0.146, "rho_ice": -0.729, "rho_s": -0.581, "C2": 0.057, "sys": -0.0}]

null directions (raw parameter units, unit-normalised): [{"zb": 0.0173, "log10w": 0.0192, "rho_ice": -0.749, "rho_s": -0.662, "C2": 0.000151, "sys": -2.09e-20}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    1.000   -1.000   -1.000      n/a    0.000
   1.000    1.000   -1.000   -1.000      n/a    0.000
  -1.000   -1.000    1.000    1.000      n/a    0.000
  -1.000   -1.000    1.000    1.000      n/a    0.000
     n/a      n/a      n/a      n/a    1.000    0.000
   0.000    0.000    0.000    0.000    0.000    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): -6.445e-17, 0.0001784, 0.104, 0.9501, 1.477, 3.469; cond = 3.47e+300

Weakest eigenvector: zb=+0.327, log10w=+0.146, rho_ice=-0.729, rho_s=-0.581, C2=+0.057, sys=-0.000

### C_B_plus_uniform_width_proxy

sigma (marginal): zb=0.8089, log10w=0.8651, rho_ice=5.767, rho_s=37.76, C2=0.08726, sys=0.004

sigma (rho_s fixed): zb=0.8088, log10w=0.8624, rho_ice=5.716, C2=0.04271, sys=0.004

corr(rho_s,C2) = -0.872; sigma(C/MR2) = 0.001621 (rho_s fixed 0.0002984)

rank 6/6, nullity 0; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): []

null directions (raw parameter units, unit-normalised): []

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    0.258   -0.120   -0.016    0.112    0.123
   0.258    1.000   -0.001    0.079    0.102   -0.000
  -0.120   -0.001    1.000    0.133    0.018    0.000
  -0.016    0.079    0.133    1.000   -0.872   -0.004
   0.112    0.102    0.018   -0.872    1.000    0.016
   0.123   -0.000    0.000   -0.004    0.016    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 0.09899, 0.6152, 0.772, 0.957, 1.267, 2.29; cond = 23.1

Weakest eigenvector: zb=+0.025, log10w=+0.005, rho_ice=-0.023, rho_s=-0.708, C2=+0.706, sys=+0.004


## Point `lw1.5`

x = zb=25, log10w=1.5, rho_ice=925, rho_s=2400, C2=0.75, sys=0

pred = C20=-0.00551219, C22=0.00153458, C30=0.000215202, libration_deg=0.0992163, cmr2_derived=0.333185

FD stencil: zb:central, log10w:central, rho_ice:central, rho_s:central, C2:central, sys:central

Step stability (relative difference of covariance-whitened Jacobian columns, step x1.0 vs x2.0): zb=0.00665, log10w=0.0178, rho_ice=0.00012, rho_s=0.000158, C2=1.37e-05, sys=0

||whitened J col||: zb=1.33, log10w=1.71, rho_ice=0.0696, rho_s=0.0602, C2=26.1, sys=33

### A_likelihood_only

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=inf

sigma (rho_s fixed): zb=inf, log10w=inf, rho_ice=inf, C2=inf, sys=inf

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed inf)

rank 4/6, nullity 2; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): [{"zb": 0.091, "log10w": 0.639, "rho_ice": -0.561, "rho_s": -0.112, "C2": 0.07, "sys": -0.501}, {"zb": -0.477, "log10w": 0.222, "rho_ice": 0.468, "rho_s": 0.548, "C2": -0.017, "sys": -0.452}]

null directions (raw parameter units, unit-normalised): [{"zb": 0.00829, "log10w": 0.0451, "rho_ice": -0.973, "rho_s": -0.225, "C2": 0.000323, "sys": -0.00184}, {"zb": -0.0316, "log10w": 0.0114, "rho_ice": 0.594, "rho_s": 0.804, "C2": -5.7e-05, "sys": -0.00121}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a      n/a      n/a      n/a      n/a
     n/a    1.000      n/a      n/a      n/a      n/a
     n/a      n/a    1.000      n/a      n/a      n/a
     n/a      n/a      n/a    1.000      n/a      n/a
     n/a      n/a      n/a      n/a    1.000      n/a
     n/a      n/a      n/a      n/a      n/a    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): -6.374e-17, 2.543e-16, 0.0001518, 0.09987, 2.471, 3.429; cond = 3.43e+300

Weakest eigenvector: zb=+0.091, log10w=+0.639, rho_ice=-0.561, rho_s=-0.112, C2=+0.070, sys=-0.501

### B_gauss_prior_sys (headline)

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=0.004

sigma (rho_s fixed): zb=55.06, log10w=52.03, rho_ice=713.5, C2=1.288, sys=0.004

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed 0.01245)

rank 5/6, nullity 1; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): [{"zb": -0.416, "log10w": -0.264, "rho_ice": 0.723, "rho_s": 0.482, "C2": -0.059, "sys": 0.0}]

null directions (raw parameter units, unit-normalised): [{"zb": -0.0238, "log10w": -0.0117, "rho_ice": 0.792, "rho_s": 0.61, "C2": -0.000173, "sys": 2.95e-20}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    1.000   -1.000   -1.000      n/a    0.000
   1.000    1.000   -1.000   -1.000      n/a    0.000
  -1.000   -1.000    1.000    1.000      n/a    0.000
  -1.000   -1.000    1.000    1.000      n/a    0.000
     n/a      n/a      n/a      n/a    1.000    0.000
   0.000    0.000    0.000    0.000    0.000    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 2.472e-16, 0.0001178, 0.09778, 0.9528, 1.534, 3.416; cond = 1.38e+16

Weakest eigenvector: zb=-0.416, log10w=-0.264, rho_ice=+0.723, rho_s=+0.482, C2=-0.059, sys=+0.000

### C_B_plus_uniform_width_proxy

sigma (marginal): zb=1.035, log10w=0.8573, rho_ice=5.767, rho_s=39.67, C2=0.09084, sys=0.004

sigma (rho_s fixed): zb=1.023, log10w=0.8242, rho_ice=5.722, C2=0.05985, sys=0.004

corr(rho_s,C2) = -0.7523; sigma(C/MR2) = 0.001629 (rho_s fixed 0.0006553)

rank 6/6, nullity 0; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): []

null directions (raw parameter units, unit-normalised): []

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    0.669   -0.102    0.151    0.235    0.095
   0.669    1.000   -0.004    0.275    0.261   -0.001
  -0.102   -0.004    1.000    0.125    0.023    0.000
   0.151    0.275    0.125    1.000   -0.752   -0.005
   0.235    0.261    0.023   -0.752    1.000    0.015
   0.095   -0.001    0.000   -0.005    0.015    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 0.09874, 0.193, 0.7577, 0.9552, 1.399, 2.596; cond = 26.3

Weakest eigenvector: zb=-0.023, log10w=+0.003, rho_ice=+0.021, rho_s=+0.710, C2=-0.703, sys=-0.004


## Point `lw1.9`

x = zb=25, log10w=1.9, rho_ice=925, rho_s=2400, C2=0.75, sys=0

pred = C20=-0.00553164, C22=0.00153936, C30=0.000216955, libration_deg=0.101543, cmr2_derived=0.333741

FD stencil: zb:central, log10w:backward, rho_ice:central, rho_s:central, C2:central, sys:central

Step stability (relative difference of covariance-whitened Jacobian columns, step x1.0 vs x2.0): zb=0.0068, log10w=0.139, rho_ice=2.96e-05, rho_s=0.000165, C2=7.14e-06, sys=6.85e-14

||whitened J col||: zb=1.38, log10w=4.16, rho_ice=0.0723, rho_s=0.0586, C2=25.7, sys=33.7

### A_likelihood_only

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=inf

sigma (rho_s fixed): zb=inf, log10w=inf, rho_ice=inf, C2=inf, sys=inf

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed inf)

rank 4/6, nullity 2; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): [{"zb": -0.056, "log10w": 0.722, "rho_ice": -0.135, "rho_s": 0.36, "C2": 0.031, "sys": -0.571}, {"zb": -0.544, "log10w": -0.243, "rho_ice": 0.663, "rho_s": 0.428, "C2": -0.042, "sys": -0.143}]

null directions (raw parameter units, unit-normalised): [{"zb": -0.00634, "log10w": 0.027, "rho_ice": -0.292, "rho_s": 0.956, "C2": 0.00019, "sys": -0.00264}, {"zb": -0.0336, "log10w": -0.00497, "rho_ice": 0.781, "rho_s": 0.623, "C2": -0.00014, "sys": -0.000362}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a      n/a      n/a      n/a      n/a
     n/a    1.000      n/a      n/a      n/a      n/a
     n/a      n/a    1.000      n/a      n/a      n/a
     n/a      n/a      n/a    1.000      n/a      n/a
     n/a      n/a      n/a      n/a    1.000      n/a
     n/a      n/a      n/a      n/a      n/a    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 2.711e-16, 9.003e-16, 0.0001284, 0.09861, 2.509, 3.392; cond = 1.25e+16

Weakest eigenvector: zb=-0.056, log10w=+0.722, rho_ice=-0.135, rho_s=+0.360, C2=+0.031, sys=-0.571

### B_gauss_prior_sys (headline)

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=0.004

sigma (rho_s fixed): zb=101.5, log10w=37.93, rho_ice=1425, C2=1.405, sys=0.004

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed 0.0121)

rank 5/6, nullity 1; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): [{"zb": -0.514, "log10w": -0.411, "rho_ice": 0.676, "rho_s": 0.328, "C2": -0.049, "sys": -0.0}]

null directions (raw parameter units, unit-normalised): [{"zb": -0.0342, "log10w": -0.00905, "rho_ice": 0.857, "rho_s": 0.514, "C2": -0.000173, "sys": -1.22e-19}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    1.000   -1.000   -1.000      n/a    0.000
   1.000    1.000   -1.000   -1.000      n/a    0.000
  -1.000   -1.000    1.000    1.000      n/a    0.000
  -1.000   -1.000    1.000    1.000      n/a    0.000
     n/a      n/a      n/a      n/a    1.000    0.000
   0.000    0.000    0.000    0.000    0.000    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 2.638e-16, 7.637e-05, 0.09459, 0.9535, 1.584, 3.368; cond = 1.28e+16

Weakest eigenvector: zb=-0.514, log10w=-0.411, rho_ice=+0.676, rho_s=+0.328, C2=-0.049, sys=-0.000

### C_B_plus_uniform_width_proxy

sigma (marginal): zb=1.76, log10w=0.8187, rho_ice=5.767, rho_s=48.44, C2=0.09873, sys=0.004

sigma (rho_s fixed): zb=1.514, log10w=0.6618, rho_ice=5.739, C2=0.09099, sys=0.004

corr(rho_s,C2) = -0.3883; sigma(C/MR2) = 0.001636 (rho_s fixed 0.001197)

rank 6/6, nullity 0; null directions (unit-diagonal-scaled Fisher coordinates, i.e. raw component x sqrt(F_kk)): []

null directions (raw parameter units, unit-normalised): []

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    0.907   -0.073    0.510    0.426    0.054
   0.907    1.000   -0.012    0.589    0.437   -0.001
  -0.073   -0.012    1.000    0.099    0.025    0.000
   0.510    0.589    0.099    1.000   -0.388   -0.006
   0.426    0.437    0.025   -0.388    1.000    0.013
   0.054   -0.001    0.000   -0.006    0.013    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 0.04555, 0.1001, 0.7363, 0.9555, 1.443, 2.719; cond = 59.7

Weakest eigenvector: zb=-0.487, log10w=-0.745, rho_ice=-0.001, rho_s=-0.409, C2=-0.200, sys=-0.003

