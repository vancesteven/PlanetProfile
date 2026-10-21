# C3 mass-deficit study

Status vocabulary: results below are `implemented, unverified` until the user runs them on the Mac cache.

- cache: `/Users/svance/scratch/enceladus_zbw_production_cache.pkl`
- refine: Ocean.deltaP /[2.0, 4.0] (gate x4); nSilMax x2 diagnostic
- sigma_obs (raw): {'C20': 3.699e-05, 'C22': 1.47e-05, 'C30': 3.342e-05, 'libration_deg': 0.003}
- sigma (model-inflated): {'C20': 3.7367768196669176e-05, 'C22': 1.4797972834141844e-05, 'C30': 3.370840251332003e-05, 'libration_deg': 0.003010398644698074}
- gate: shift extrapolated to zero deficit <= 0.1 sigma_obs (max over compensation_C2 in [0.0, 0.5, 1.0])

## Map / roughness

- finite nodes 3399, masked 121, range [-0.005079502934320468, -0.0003722761790641016]
- jagged_flag: **True** (max slope-sign-change fraction > 0.4 (smooth ~0, white noise ~0.67). A quantised-step deficit is sawtooth by construction, so judge this against the observable shifts, not in isolation.)
- along_zb: abs dev median 0.0010224647417304555, p90 0.002338574465520583, max 0.004240885178729281 (0.900930717644664 of range); rel median 0.33542125218276914, max 11.028178323920036; slope sign changes 2187/3319
- along_log10w: abs dev median 4.009096790103016e-05, p90 0.0017370364378379258, max 0.004119403669418431 (0.8751232697295627 of range); rel median 0.01511981254478893, max 10.121867258336156; slope sign changes 1034/3159

## Nodes (gate: extrapolated-to-zero-deficit shift at x4, worst over C2, in sigma_obs)

| node | role | zb | w | def cached | repro | x2 | x4 | bound | C20 | C22 | C30 | libration_deg | PASS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| (70,18) | worst | 39.500 | 2.424 | -0.00508 | -0.00508 | -0.00161 | -0.00158 | 0.00441 | +1.7822 | -1.2930 | -0.0614 | -0.0988 | FAIL |
| (77,10) | worst | 42.500 | 0.588 | -0.00508 | -0.00508 | -0.00164 | -0.00045 | 0.00442 | +1.3765 | -0.9885 | -0.0610 | -0.0707 | FAIL |
| (80,1) | worst | 43.250 | 0.119 | -0.00507 | -0.00507 | -0.00164 | -0.00045 | 0.00443 | +1.3728 | -0.9854 | -0.0609 | -0.0685 | FAIL |
| (80,2) | worst | 43.250 | 0.143 | -0.00507 | -0.00507 | -0.00164 | -0.00045 | 0.00443 | +1.3717 | -0.9847 | -0.0609 | -0.0684 | FAIL |
| (80,3) | worst | 43.250 | 0.170 | -0.00506 | -0.00506 | -0.00163 | -0.00044 | 0.00443 | +1.3704 | -0.9838 | -0.0608 | -0.0684 | FAIL |
| (61,15) | control_median | 35.000 | 1.425 | -0.00315 | -0.00315 | -0.00046 | -0.00043 | 0.00445 | +0.1370 | -0.0794 | -0.0373 | -0.0195 | FAIL |
| (66,34) | control_median | 37.500 | 41.246 | -0.00315 | -0.00315 | -0.00090 | -0.00087 | 0.00419 | +0.1386 | -0.0818 | -0.0373 | -0.0199 | FAIL |

## Findings

- (70,18): |deficit| 5.080e-03 exceeds predicted one-layer bound 4.408e-03 (1% slack): deficit is NOT explained by quantisation alone.
- (70,18): FLOOR: refined (x4) |residual| 1.578e-03 did not fall well below the repro 5.080e-03 (or stays >= 3e-4): a second, non-quantisation mass error source exists.
- (77,10): |deficit| 5.078e-03 exceeds predicted one-layer bound 4.424e-03 (1% slack): deficit is NOT explained by quantisation alone.
- (77,10): FLOOR: refined (x4) |residual| 4.461e-04 did not fall well below the repro 5.078e-03 (or stays >= 3e-4): a second, non-quantisation mass error source exists.
- (77,10): LINEARITY: extrapolations from factors 2 and 4 disagree by >20% (see linearity_f2_over_f4); extrapolation unreliable.
- (80,1): |deficit| 5.072e-03 exceeds predicted one-layer bound 4.429e-03 (1% slack): deficit is NOT explained by quantisation alone.
- (80,1): FLOOR: refined (x4) |residual| 4.504e-04 did not fall well below the repro 5.072e-03 (or stays >= 3e-4): a second, non-quantisation mass error source exists.
- (80,1): LINEARITY: extrapolations from factors 2 and 4 disagree by >20% (see linearity_f2_over_f4); extrapolation unreliable.
- (80,2): |deficit| 5.068e-03 exceeds predicted one-layer bound 4.428e-03 (1% slack): deficit is NOT explained by quantisation alone.
- (80,2): FLOOR: refined (x4) |residual| 4.465e-04 did not fall well below the repro 5.068e-03 (or stays >= 3e-4): a second, non-quantisation mass error source exists.
- (80,2): LINEARITY: extrapolations from factors 2 and 4 disagree by >20% (see linearity_f2_over_f4); extrapolation unreliable.
- (80,3): |deficit| 5.064e-03 exceeds predicted one-layer bound 4.428e-03 (1% slack): deficit is NOT explained by quantisation alone.
- (80,3): FLOOR: refined (x4) |residual| 4.419e-04 did not fall well below the repro 5.064e-03 (or stays >= 3e-4): a second, non-quantisation mass error source exists.
- (80,3): LINEARITY: extrapolations from factors 2 and 4 disagree by >20% (see linearity_f2_over_f4); extrapolation unreliable.
- (61,15): FLOOR: refined (x4) |residual| 4.262e-04 did not fall well below the repro 3.153e-03 (or stays >= 3e-4): a second, non-quantisation mass error source exists.
- (66,34): FLOOR: refined (x4) |residual| 8.722e-04 did not fall well below the repro 3.152e-03 (or stays >= 3e-4): a second, non-quantisation mass error source exists.

## Diagnostics

- nSilMax x2 on top of deltaP/4, node (70,18) (not gated): {'knob': {'Ocean.deltaP': (0.1, 0.025), 'Steps.nSilMax': (120, 240)}, 'mass_residual': -0.0006152190577932704, 'dsigma_vs_deltaP4_alone': {'C20': 0.010297221283316704, 'C22': 0.01572745246003876, 'C30': 0.011338371179765518, 'libration_deg': 0.003493686741587454}, 'note': 'DIAGNOSTIC, not gated'}

**Overall (worst-deficit nodes only): FAIL / incomplete**
