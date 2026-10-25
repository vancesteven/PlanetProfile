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
| rho_s2300 | 5/6, 5/5 | n/a (limit, singular F) | ok | inf | 0.8935 | inf | 24.31 | inf | inf | 0.007266 |
| best_fit_rho_s2400 | 5/6, 5/5 | n/a (limit, singular F) | ok | inf | 1.16 | inf | 29.39 | inf | inf | 0.01196 |
| rho_s2500 | 5/6, 5/5 | -1 (limit, singular F) | ok | inf | 1.599 | inf | 37.76 | inf | inf | 0.01991 |
| zb20 | 5/6, 5/5 | n/a (limit, singular F) | ok | inf | 1.001 | inf | 44.67 | inf | inf | 0.01272 |
| zb30 | 5/6, 5/5 | n/a (limit, singular F) | ok | inf | 1.294 | inf | 31.41 | inf | inf | 0.01091 |

## Summary (variant C_B_plus_uniform_width_proxy)

| point | rank free/fixed | corr(rho_s,C2) | flag | sig(C2) free | sig(C2) rho_s fixed | sig(zb) free | sig(zb) rho_s fixed | sig(rho_s) | sig(C/MR2) free | sig(C/MR2) rho_s fixed |
|---|---|---|---|---|---|---|---|---|---|---|
| rho_s2300 | 6/6, 5/5 | -0.8782 | ok | 0.08855 | 0.04236 | 0.7784 | 0.7774 | 31.05 | 0.001518 | 0.000203 |
| best_fit_rho_s2400 | 6/6, 5/5 | -0.8857 | ok | 0.08649 | 0.04016 | 0.7881 | 0.7877 | 37.49 | 0.001619 | 0.0002274 |
| rho_s2500 | 6/6, 5/5 | -0.8894 | ok | 0.08421 | 0.03849 | 0.7969 | 0.7968 | 44.12 | 0.001699 | 0.0002501 |
| zb20 | 6/6, 5/5 | -0.8888 | ok | 0.07849 | 0.03596 | 0.4773 | 0.4769 | 36.94 | 0.00158 | 0.0002088 |
| zb30 | 6/6, 5/5 | -0.8792 | ok | 0.09573 | 0.04561 | 1.195 | 1.195 | 37.93 | 0.001654 | 0.0002451 |


## Point `rho_s2300`

x = zb=25, log10w=0.5, rho_ice=925, rho_s=2300, C2=0.75, sys=0

pred = C20=-0.0056437, C22=0.00157506, C30=0.000213832, libration_deg=0.0989193, cmr2_derived=0.337508

FD stencil: zb:central, log10w:central, rho_ice:central, rho_s:central, C2:central, sys:central

Step stability (relative difference of covariance-whitened Jacobian columns, step x1.0 vs x2.0): zb=0.007, log10w=0.00709, rho_ice=0.000347, rho_s=0.00018, C2=2.15e-05, sys=7.01e-14

||whitened J col||: zb=1.33, log10w=0.155, rho_ice=0.0627, rho_s=0.071, C2=24.8, sys=32.9

### A_likelihood_only

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=inf

sigma (rho_s fixed): zb=inf, log10w=inf, rho_ice=inf, C2=inf, sys=inf

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed inf)

rank 4/6, nullity 2; null directions (components in raw-parameter units, scaled by sqrt(F_kk)): [{"zb": -0.537, "log10w": 0.338, "rho_ice": 0.118, "rho_s": 0.264, "C2": -0.021, "sys": -0.717}, {"zb": -0.072, "log10w": -0.165, "rho_ice": 0.742, "rho_s": 0.568, "C2": 0.012, "sys": 0.307}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a      n/a      n/a      n/a      n/a
     n/a    1.000      n/a      n/a      n/a      n/a
     n/a      n/a    1.000      n/a      n/a      n/a
     n/a      n/a      n/a    1.000      n/a      n/a
     n/a      n/a      n/a      n/a    1.000      n/a
     n/a      n/a      n/a      n/a      n/a    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 1.309e-18, 1.423e-16, 0.0003458, 0.1189, 2.551, 3.33; cond = 2.54e+18

Weakest eigenvector: zb=-0.537, log10w=+0.338, rho_ice=+0.118, rho_s=+0.264, C2=-0.021, sys=-0.717

### B_gauss_prior_sys (headline)

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=0.004

sigma (rho_s fixed): zb=24.31, log10w=256, rho_ice=221.9, C2=0.8935, sys=0.004

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed 0.007266)

rank 5/6, nullity 1; null directions (components in raw-parameter units, scaled by sqrt(F_kk)): [{"zb": 0.278, "log10w": 0.019, "rho_ice": -0.728, "rho_s": -0.626, "C2": -0.003, "sys": 0.0}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a   -1.000   -1.000      n/a    0.000
     n/a    1.000      n/a      n/a      n/a    0.000
  -1.000      n/a    1.000    1.000      n/a    0.000
  -1.000      n/a    1.000    1.000      n/a    0.000
     n/a      n/a      n/a      n/a    1.000    0.000
   0.000    0.000    0.000    0.000    0.000    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 2.711e-16, 0.0003454, 0.1187, 0.9562, 1.603, 3.322; cond = 1.23e+16

Weakest eigenvector: zb=+0.278, log10w=+0.019, rho_ice=-0.728, rho_s=-0.626, C2=-0.003, sys=+0.000

### C_B_plus_uniform_width_proxy

sigma (marginal): zb=0.7784, log10w=0.8659, rho_ice=5.769, rho_s=31.05, C2=0.08855, sys=0.004

sigma (rho_s fixed): zb=0.7774, log10w=0.8658, rho_ice=5.719, C2=0.04236, sys=0.004

corr(rho_s,C2) = -0.8782; sigma(C/MR2) = 0.001518 (rho_s fixed 0.000203)

rank 6/6, nullity 0; null directions (components in raw-parameter units, scaled by sqrt(F_kk)): []

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    0.089   -0.126   -0.050    0.109    0.128
   0.089    1.000   -0.000    0.020    0.032   -0.000
  -0.126   -0.000    1.000    0.132    0.012    0.000
  -0.050    0.020    0.132    1.000   -0.878   -0.003
   0.109    0.032    0.012   -0.878    1.000    0.016
   0.128   -0.000    0.000   -0.003    0.016    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 0.108, 0.7784, 0.8447, 0.9803, 1.155, 2.134; cond = 19.8

Weakest eigenvector: zb=+0.032, log10w=+0.003, rho_ice=-0.025, rho_s=-0.708, C2=+0.705, sys=+0.004


## Point `best_fit_rho_s2400`

x = zb=25, log10w=0.5, rho_ice=925, rho_s=2400, C2=0.75, sys=0

pred = C20=-0.00549986, C22=0.00153205, C30=0.000213286, libration_deg=0.0979418, cmr2_derived=0.33287

FD stencil: zb:central, log10w:central, rho_ice:central, rho_s:central, C2:central, sys:central

Step stability (relative difference of covariance-whitened Jacobian columns, step x1.0 vs x2.0): zb=0.00693, log10w=0.00653, rho_ice=0.00036, rho_s=0.000155, C2=2.31e-05, sys=7.08e-14

||whitened J col||: zb=1.31, log10w=0.169, rho_ice=0.0669, rho_s=0.0612, C2=26.4, sys=32.5

### A_likelihood_only

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=inf

sigma (rho_s fixed): zb=inf, log10w=inf, rho_ice=inf, C2=inf, sys=inf

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed inf)

rank 4/6, nullity 2; null directions (components in raw-parameter units, scaled by sqrt(F_kk)): [{"zb": -0.373, "log10w": 0.503, "rho_ice": -0.253, "rho_s": 0.028, "C2": 0.048, "sys": -0.735}, {"zb": -0.335, "log10w": -0.046, "rho_ice": 0.701, "rho_s": 0.62, "C2": -0.049, "sys": -0.083}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a      n/a      n/a      n/a      n/a
     n/a    1.000      n/a      n/a      n/a      n/a
     n/a      n/a    1.000      n/a      n/a      n/a
     n/a      n/a      n/a    1.000      n/a      n/a
     n/a      n/a      n/a      n/a    1.000      n/a
     n/a      n/a      n/a      n/a      n/a    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): -9.931e-17, 2.054e-16, 0.000228, 0.1083, 2.386, 3.505; cond = 3.5e+300

Weakest eigenvector: zb=-0.373, log10w=+0.503, rho_ice=-0.253, rho_s=+0.028, C2=+0.048, sys=-0.735

### B_gauss_prior_sys (headline)

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=0.004

sigma (rho_s fixed): zb=29.39, log10w=320.1, rho_ice=335.4, C2=1.16, sys=0.004

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed 0.01196)

rank 5/6, nullity 1; null directions (components in raw-parameter units, scaled by sqrt(F_kk)): [{"zb": 0.291, "log10w": 0.102, "rho_ice": -0.725, "rho_s": -0.613, "C2": 0.054, "sys": 0.0}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    1.000   -1.000   -1.000      n/a    0.000
   1.000    1.000   -1.000   -1.000      n/a    0.000
  -1.000   -1.000    1.000    1.000      n/a    0.000
  -1.000   -1.000    1.000    1.000      n/a    0.000
     n/a      n/a      n/a      n/a    1.000    0.000
   0.000    0.000    0.000    0.000    0.000    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 1.463e-16, 0.0002151, 0.1078, 0.9481, 1.446, 3.498; cond = 2.39e+16

Weakest eigenvector: zb=+0.291, log10w=+0.102, rho_ice=-0.725, rho_s=-0.613, C2=+0.054, sys=+0.000

### C_B_plus_uniform_width_proxy

sigma (marginal): zb=0.7881, log10w=0.8659, rho_ice=5.767, rho_s=37.49, C2=0.08649, sys=0.004

sigma (rho_s fixed): zb=0.7877, log10w=0.8657, rho_ice=5.714, C2=0.04016, sys=0.004

corr(rho_s,C2) = -0.8857; sigma(C/MR2) = 0.001619 (rho_s fixed 0.0002274)

rank 6/6, nullity 0; null directions (components in raw-parameter units, scaled by sqrt(F_kk)): []

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    0.082   -0.119   -0.033    0.094    0.127
   0.082    1.000   -0.000    0.024    0.035   -0.000
  -0.119   -0.000    1.000    0.135    0.015    0.000
  -0.033    0.024    0.135    1.000   -0.886   -0.003
   0.094    0.035    0.015   -0.886    1.000    0.016
   0.127   -0.000    0.000   -0.003    0.016    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 0.09917, 0.774, 0.8288, 0.9773, 1.145, 2.176; cond = 21.9

Weakest eigenvector: zb=+0.024, log10w=+0.002, rho_ice=-0.024, rho_s=-0.709, C2=+0.705, sys=+0.004


## Point `rho_s2500`

x = zb=25, log10w=0.5, rho_ice=925, rho_s=2500, C2=0.75, sys=0

pred = C20=-0.0053753, C22=0.00149485, C30=0.00021281, libration_deg=0.0970984, cmr2_derived=0.328755

FD stencil: zb:central, log10w:central, rho_ice:central, rho_s:central, C2:central, sys:central

Step stability (relative difference of covariance-whitened Jacobian columns, step x1.0 vs x2.0): zb=0.00686, log10w=0.00603, rho_ice=0.000371, rho_s=0.000135, C2=2.45e-05, sys=7.16e-14

||whitened J col||: zb=1.29, log10w=0.186, rho_ice=0.0706, rho_s=0.0532, C2=27.9, sys=32.3

### A_likelihood_only

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=inf

sigma (rho_s fixed): zb=inf, log10w=inf, rho_ice=inf, C2=inf, sys=inf

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed inf)

rank 4/6, nullity 2; null directions (components in raw-parameter units, scaled by sqrt(F_kk)): [{"zb": -0.408, "log10w": 0.535, "rho_ice": 0.069, "rho_s": 0.354, "C2": 0.081, "sys": -0.641}, {"zb": -0.218, "log10w": -0.365, "rho_ice": 0.723, "rho_s": 0.495, "C2": -0.154, "sys": 0.165}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a      n/a      n/a      n/a      n/a
     n/a    1.000      n/a      n/a      n/a      n/a
     n/a      n/a    1.000      n/a      n/a      n/a
     n/a      n/a      n/a    1.000      n/a      n/a
     n/a      n/a      n/a      n/a    1.000      n/a
     n/a      n/a      n/a      n/a      n/a    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): -3.699e-17, 3.527e-16, 0.0001634, 0.1005, 2.273, 3.627; cond = 3.63e+300

Weakest eigenvector: zb=-0.408, log10w=+0.535, rho_ice=+0.069, rho_s=+0.354, C2=+0.081, sys=-0.641

### B_gauss_prior_sys (headline)

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=0.004

sigma (rho_s fixed): zb=37.76, log10w=425, rho_ice=511.9, C2=1.599, sys=0.004

corr(rho_s,C2) = -1; sigma(C/MR2) = inf (rho_s fixed 0.01991)

rank 5/6, nullity 1; null directions (components in raw-parameter units, scaled by sqrt(F_kk)): [{"zb": 0.313, "log10w": 0.22, "rho_ice": -0.717, "rho_s": -0.568, "C2": 0.129, "sys": -0.0}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    1.000   -1.000   -1.000    1.000    0.000
   1.000    1.000   -1.000   -1.000    1.000    0.000
  -1.000   -1.000    1.000    1.000   -1.000    0.000
  -1.000   -1.000    1.000    1.000   -1.000    0.000
   1.000    1.000   -1.000   -1.000    1.000    0.000
   0.000    0.000    0.000    0.000    0.000    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): -1.968e-16, 0.0001285, 0.09966, 0.9385, 1.339, 3.622; cond = 3.62e+300

Weakest eigenvector: zb=+0.313, log10w=+0.220, rho_ice=-0.717, rho_s=-0.568, C2=+0.129, sys=-0.000

### C_B_plus_uniform_width_proxy

sigma (marginal): zb=0.7969, log10w=0.8659, rho_ice=5.764, rho_s=44.12, C2=0.08421, sys=0.004

sigma (rho_s fixed): zb=0.7968, log10w=0.8656, rho_ice=5.71, C2=0.03849, sys=0.004

corr(rho_s,C2) = -0.8894; sigma(C/MR2) = 0.001699 (rho_s fixed 0.0002501)

rank 6/6, nullity 0; null directions (components in raw-parameter units, scaled by sqrt(F_kk)): []

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    0.076   -0.112   -0.018    0.081    0.126
   0.076    1.000   -0.001    0.027    0.039   -0.000
  -0.112   -0.001    1.000    0.137    0.020    0.000
  -0.018    0.027    0.137    1.000   -0.889   -0.003
   0.081    0.039    0.020   -0.889    1.000    0.016
   0.126   -0.000    0.000   -0.003    0.016    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 0.09394, 0.7655, 0.8184, 0.9735, 1.136, 2.212; cond = 23.5

Weakest eigenvector: zb=+0.018, log10w=+0.002, rho_ice=-0.023, rho_s=-0.709, C2=+0.704, sys=+0.004


## Point `zb20`

x = zb=20, log10w=0.5, rho_ice=925, rho_s=2400, C2=0.75, sys=0

pred = C20=-0.00545503, C22=0.00153134, C30=0.000191408, libration_deg=0.122951, cmr2_derived=0.333264

FD stencil: zb:central, log10w:central, rho_ice:central, rho_s:central, C2:central, sys:central

Step stability (relative difference of covariance-whitened Jacobian columns, step x1.0 vs x2.0): zb=0.0112, log10w=0.08, rho_ice=0.000451, rho_s=0.000154, C2=2.26e-05, sys=0

||whitened J col||: zb=2.22, log10w=0.337, rho_ice=0.0628, rho_s=0.062, C2=29.1, sys=40.8

### A_likelihood_only

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=inf

sigma (rho_s fixed): zb=inf, log10w=inf, rho_ice=inf, C2=inf, sys=inf

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed inf)

rank 4/6, nullity 2; null directions (components in raw-parameter units, scaled by sqrt(F_kk)): [{"zb": -0.324, "log10w": 0.056, "rho_ice": 0.698, "rho_s": 0.635, "C2": -0.013, "sys": -0.029}, {"zb": 0.165, "log10w": -0.589, "rho_ice": 0.204, "rho_s": -0.055, "C2": -0.023, "sys": 0.762}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a      n/a      n/a      n/a      n/a
     n/a    1.000      n/a      n/a      n/a      n/a
     n/a      n/a    1.000      n/a      n/a      n/a
     n/a      n/a      n/a    1.000      n/a      n/a
     n/a      n/a      n/a      n/a    1.000      n/a
     n/a      n/a      n/a      n/a      n/a    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): -1.037e-16, 1.431e-16, 6.111e-05, 0.1074, 2.893, 3; cond = 3e+300

Weakest eigenvector: zb=-0.324, log10w=+0.056, rho_ice=+0.698, rho_s=+0.635, C2=-0.013, sys=-0.029

### B_gauss_prior_sys (headline)

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=0.004

sigma (rho_s fixed): zb=44.67, log10w=289.1, rho_ice=546.9, C2=1.001, sys=0.004

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed 0.01272)

rank 5/6, nullity 1; null directions (components in raw-parameter units, scaled by sqrt(F_kk)): [{"zb": -0.318, "log10w": 0.034, "rho_ice": 0.706, "rho_s": 0.632, "C2": -0.014, "sys": -0.0}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a   -1.000   -1.000      n/a    0.000
     n/a    1.000      n/a      n/a      n/a    0.000
  -1.000      n/a    1.000    1.000      n/a    0.000
  -1.000      n/a    1.000    1.000      n/a    0.000
     n/a      n/a      n/a      n/a    1.000    0.000
   0.000    0.000    0.000    0.000    0.000    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 2.357e-16, 5.391e-05, 0.1073, 0.9492, 1.994, 2.949; cond = 1.25e+16

Weakest eigenvector: zb=-0.318, log10w=+0.034, rho_ice=+0.706, rho_s=+0.632, C2=-0.014, sys=-0.000

### C_B_plus_uniform_width_proxy

sigma (marginal): zb=0.4773, log10w=0.8659, rho_ice=5.768, rho_s=36.94, C2=0.07849, sys=0.004

sigma (rho_s fixed): zb=0.4769, log10w=0.8653, rho_ice=5.723, C2=0.03596, sys=0.004

corr(rho_s,C2) = -0.8888; sigma(C/MR2) = 0.00158 (rho_s fixed 0.0002088)

rank 6/6, nullity 0; null directions (components in raw-parameter units, scaled by sqrt(F_kk)): []

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    0.232   -0.167   -0.041    0.060    0.154
   0.232    1.000   -0.000    0.038    0.037   -0.000
  -0.167   -0.000    1.000    0.125    0.008    0.000
  -0.041    0.038    0.125    1.000   -0.889   -0.003
   0.060    0.037    0.008   -0.889    1.000    0.013
   0.154   -0.000    0.000   -0.003    0.013    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 0.09961, 0.6851, 0.8179, 0.96, 1.342, 2.095; cond = 21

Weakest eigenvector: zb=+0.020, log10w=+0.000, rho_ice=-0.023, rho_s=-0.710, C2=+0.704, sys=+0.003


## Point `zb30`

x = zb=30, log10w=0.5, rho_ice=925, rho_s=2400, C2=0.75, sys=0

pred = C20=-0.00554396, C22=0.0015333, C30=0.000233362, libration_deg=0.0825126, cmr2_derived=0.33255

FD stencil: zb:central, log10w:central, rho_ice:central, rho_s:central, C2:central, sys:central

Step stability (relative difference of covariance-whitened Jacobian columns, step x1.0 vs x2.0): zb=0.00416, log10w=0.0143, rho_ice=0.000294, rho_s=0.000155, C2=2.34e-05, sys=8.41e-14

||whitened J col||: zb=0.879, log10w=0.122, rho_ice=0.0716, rho_s=0.0605, C2=24, sys=27.4

### A_likelihood_only

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=inf

sigma (rho_s fixed): zb=inf, log10w=inf, rho_ice=inf, C2=inf, sys=inf

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed inf)

rank 4/6, nullity 2; null directions (components in raw-parameter units, scaled by sqrt(F_kk)): [{"zb": -0.494, "log10w": 0.371, "rho_ice": 0.258, "rho_s": 0.39, "C2": 0.003, "sys": -0.632}, {"zb": 0.023, "log10w": 0.396, "rho_ice": -0.693, "rho_s": -0.472, "C2": 0.108, "sys": -0.359}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000      n/a      n/a      n/a      n/a      n/a
     n/a    1.000      n/a      n/a      n/a      n/a
     n/a      n/a    1.000      n/a      n/a      n/a
     n/a      n/a      n/a    1.000      n/a      n/a
     n/a      n/a      n/a      n/a    1.000      n/a
     n/a      n/a      n/a      n/a      n/a    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): -5.013e-16, 2.105e-16, 0.0003631, 0.1131, 2.296, 3.591; cond = 3.59e+300

Weakest eigenvector: zb=-0.494, log10w=+0.371, rho_ice=+0.258, rho_s=+0.390, C2=+0.003, sys=-0.632

### B_gauss_prior_sys (headline)

sigma (marginal): zb=inf, log10w=inf, rho_ice=inf, rho_s=inf, C2=inf, sys=0.004

sigma (rho_s fixed): zb=31.41, log10w=337.7, rho_ice=163.4, C2=1.294, sys=0.004

corr(rho_s,C2) = n/a; sigma(C/MR2) = inf (rho_s fixed 0.01091)

rank 5/6, nullity 1; null directions (components in raw-parameter units, scaled by sqrt(F_kk)): [{"zb": 0.263, "log10w": 0.161, "rho_ice": -0.73, "rho_s": -0.603, "C2": 0.093, "sys": -0.0}]

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    1.000   -1.000   -1.000      n/a    0.000
   1.000    1.000   -1.000   -1.000      n/a    0.000
  -1.000   -1.000    1.000    1.000      n/a    0.000
  -1.000   -1.000    1.000    1.000      n/a    0.000
     n/a      n/a      n/a      n/a    1.000    0.000
   0.000    0.000    0.000    0.000    0.000    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 1.679e-16, 0.0003425, 0.1121, 0.9553, 1.344, 3.588; cond = 2.14e+16

Weakest eigenvector: zb=+0.263, log10w=+0.161, rho_ice=-0.730, rho_s=-0.603, C2=+0.093, sys=-0.000

### C_B_plus_uniform_width_proxy

sigma (marginal): zb=1.195, log10w=0.8659, rho_ice=5.766, rho_s=37.93, C2=0.09573, sys=0.004

sigma (rho_s fixed): zb=1.195, log10w=0.8659, rho_ice=5.707, C2=0.04561, sys=0.004

corr(rho_s,C2) = -0.8792; sigma(C/MR2) = 0.001654 (rho_s fixed 0.0002451)

rank 6/6, nullity 0; null directions (components in raw-parameter units, scaled by sqrt(F_kk)): []

Posterior correlation matrix (order: zb, log10w, rho_ice, rho_s, C2, sys):

```
   1.000    0.060   -0.088   -0.023    0.139    0.107
   0.060    1.000   -0.000    0.013    0.032   -0.000
  -0.088   -0.000    1.000    0.142    0.022    0.000
  -0.023    0.013    0.142    1.000   -0.879   -0.003
   0.139    0.032    0.022   -0.879    1.000    0.019
   0.107   -0.000    0.000   -0.003    0.019    1.000
```

Eigenvalues of unit-diagonal-scaled Fisher (ascending): 0.09899, 0.7694, 0.7962, 0.9866, 1.075, 2.274; cond = 23

Weakest eigenvector: zb=+0.032, log10w=+0.004, rho_ice=-0.024, rho_s=-0.706, C2=+0.707, sys=+0.004

