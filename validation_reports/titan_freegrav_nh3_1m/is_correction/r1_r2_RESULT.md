# R1 / R2 result: NH3 corrected-posterior crosscheck, C10 re-run, C11: ALL PASS. Fiducial-scope lift TAKES EFFECT.

- **Run:** 2026-10-20 on Machine A, driver `plans/scripts/nh3_r1r2_gates.py` at 7f52c758.
- **Preregistration:** `r1_r2_PREREGISTRATION.md` and `r1_r2_prereg_bands.json` (634ea0cd), committed before this run.
- **Report:** `r1_r2_report.json`, which holds every input hash and every number.
- **Rules:** applied exactly as preregistered. No threshold, seed or rule was revisited after results were read.

## Provenance (all blocking asserts passed)

- Bundle sha256 is `0818b4e5…`.
- `log_w` is bit-exact against `log_weights.npy`: 19774 finite, 226 −inf, identical mask.
- `systematic_resample(w, 1270, seed=72)` reproduces `resampled_indices.npy` exactly.
- The reference config_hash `1611b65f…` is identical across the three seeds and equals the current config.
- The cache sha is `3d837cf8…` on both sides.
- The config sigmas match.
- The pooled reference ESS is 39097.
- The C9 bands re-derive to better than 1e-9.
- R2 plumbing: 200 reference draws and 200 corrected draws were each recomputed at HEAD against their own stored C20/C22. The max relative diff is 0.0 on both sides.

## R1a: G1 crosscheck, 13 parameters

The corrected set is L=1270 against the n_eff=8000 pooled reference, with a matched-n floor (floor_ess = 1270).

**All 13 parameters pass the stock gates**, so the C9 fallback was never invoked.

The table gives each value against its tolerance. Columns:
- **|Δmean|**: absolute difference of the means.
- **σ-ratio**: ratio of the posterior standard deviations.
- **D-excess**: KS shape excess.
- **|Δmedian|**: absolute difference of the medians.

| Parameter | |Δmean| / tol | σ-ratio | D-excess / tol | |Δmedian| / tol |
|---|---|---|---|---|
| alpha | 0.0022 / 0.021 | 1.042 | 0.034 / 0.094 | 0.0006 / 0.025 |
| log10_zeta | 0.034 / 0.36 | 1.030 | 0.016 / 0.094 | 0.032 / 0.3 |
| log10_eta_Ih | 0.046 / 0.42 | 1.023 | 0.034 / 0.097 | 0.060 / 0.3 |
| log10_eta_III | 0.019 / 0.42 | 1.035 | 0.018 / 0.092 | 0.032 / 0.3 |
| log10_eta_V | 0.070 / 0.43 | 1.016 | 0.032 / 0.099 | 0.136 / 0.3 |
| log10_eta_VI | 0.043 / 0.42 | 1.008 | 0.029 / 0.094 | 0.035 / 0.3 |
| log10_eta_sil | 0.003 / 0.28 | 1.025 | 0.049 / 0.099 | 0.102 / 0.3 |
| Tb_K | 0.28 / 0.97 | 1.016 | 0.031 / 0.097 | 0.46 / 1.17 |
| log10_wOcean_ppt | 0.010 / 0.13 | 1.050 | 0.019 / 0.090 | 0.013 / 0.3 |
| R_core_km | 12.8 / 125 | 1.018 | 0.024 / 0.090 | 15.5 / 150 |
| rho_core_kgm3 | 2.3 / 57.6 | 1.005 | 0.032 / 0.099 | 0.36 / 69.1 |
| dC20_nh | 5.2e-8 / 2.3e-7 | 1.070 | 0.023 / 0.094 | 1.6e-8 / 2.8e-7 |
| dC22_nh | 2.4e-8 / 6.1e-8 | 1.095 | 0.004 / 0.102 | 9.2e-9 / 7.3e-8 |

**Contrast with the uncorrected flow.** The uncorrected crosscheck of 2026-08-03 FAILED the three η medians. All three pass after correction.

## R1b: C10 against the n_eff=8000 pool, PASS

| Quantity | Corrected median | Reference median | Gap | Weighted KS |
|---|---|---|---|---|
| Re k₂ | 0.5937 | 0.5911 | 0.055 σ_obs | 0.032 |
| \|Im k₂\| | 0.1089 | 0.1055 | 0.098 σ_obs | 0.048 |

The gates are KS ≤ 0.15 and gap ≤ 0.5 σ_obs.

## R2: C11 gravity no-regression, PASS

| Quantity | Corrected median | Reference median | Shift | Gate |
|---|---|---|---|---|
| C20 | −3.350e-5 | −3.351e-5 | 0.0071 σ_obs | ≤ 0.1 |
| C22 | 1.011e-5 | 1.011e-5 | **0.0975 σ_obs** | ≤ 0.1 |

**Thin margin, recorded rather than reinterpreted.** The C22 shift passes at 0.0975 against the 0.1 gate. The preregistration put the null pass rate at about 99%, so a value this close to the line is uncommon but expected occasionally.
- Under the outcome rule it is a PASS, and the rule is not revisited after results are read.
- **Scope, as preregistered:** R2 is a gravity no-regression check only. dC20_nh and dC22_nh absorb structure, so R2 says nothing about the structure sector. R1a carries that.
- The dC22_nh crosscheck row itself is comfortably inside every gate.

## Outcome

**LIFT_TAKES_EFFECT.** This follows the binding 2026-10-20 adjudication (RATIFICATION.md): on an R1/R2 PASS, the fiducial-scope lift takes effect without further adjudication.

**Lifted at the DEFAULT Cassini datum only (C8(i)):**
- Re k₂ and |Im k₂| pushforwards.
- ζ and all five η's and alpha. All pass the stock gates, so none needs the "unidentified nuisance" qualifier and none is quarantined by name.

**Unchanged:**
- Any modified datum falls back to the uncorrected, quarantined tidal sector. C12's FAIL bounds the scope, and C8(iii) is not authorized.
- The model–data k₂ tension caveat stays: corrected |Im k₂| ≈ 0.109 against observed 0.135 ± 0.035, shared with the MCMC.

**Deployment status: lift ratified, GUI NOT yet switched on.**
- `is_fiducial_bundle.FEATURE_ENABLED` stays False until B8 visual verification is done:
  - at the default datum, the corrected warning shows, |Im k₂| ≈ 0.109 and ocean fraction ≈ 0.929;
  - after a one-field change, the fallback warning shows.
- Per CLAUDE.md, the GUI change stays `implemented, unverified` until then.
