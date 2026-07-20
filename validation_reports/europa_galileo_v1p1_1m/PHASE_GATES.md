# Europa Galileo v1.1 (8D) — SBI validation gates

**Artifact:** `europa_galileo_v1p1_8D_posterior_1m.pt` (nsf, sbi 0.26.1, seed 44)
**Config:** `PlanetProfile/Inference/configs/europa_galileo_v1p1_8D.json`
**Cache:** Test51_seawater 1D (Tb) structure grid
**Trained:** 831,903 sims (16.8% support-cut rejection), 5D x (CMR2, Re_k2, Im_k2, Re_h2, Im_h2), 8D theta. Converged 51 epochs / 45.1 min.
**Reference MCMC:** `galileo_v1p1_reference_result.pkl` (pocoMC, 4267 samples, r_hat 1.0, ESS 4179.7).
**Seeds:** train 44 / data 50 / noise 5050.

Honest-data framing (config metadata): the ONLY Galileo-era measurement is CMR2 [0.3547, 0.0024] (GC21 MoI) + the synodic |Ae|>0.7 support cut. k2/h2 are LABELED hypothetical-conditioning channels (Re_k2 0.23±0.05, Im_k2 0.004±0.05, Re_h2 1.2±0.1, Im_h2 0±0.1), NOT measurements.

Gate suite: **SBC + crosscheck** only (no degeneracy/2D gate — v1.1 has no salinity parameter). Plus the user-directed zeta-joint + heating-fraction deliverable.

---

## Gate 1 — SBC (simulation-based calibration): **PASS**

1000 SBC pairs, 1000 posterior samples each, α=0.05. **min KS p = 0.06742 ≥ 0.05.** All 8 params pass; c2st rank-accuracy 0.56–0.59 (≈0.5 ideal → well-calibrated).

| param | KS p | c2st |
|---|---|---|
| alpha | 0.956 | 0.563 |
| log10_zeta_Ih | 0.452 | 0.585 |
| log10_zeta_sil | 0.093 | 0.564 |
| log10_eta_Ih | 0.500 | 0.576 |
| log10_eta_sil | 0.067 | 0.568 |
| Tb_K | 0.406 | 0.592 |
| R_core_km | 0.323 | 0.566 |
| rho_core_kgm3 | 0.169 | 0.560 |

Report: `reports/sbc/sbc_report.json` + rank plot.

---

## Gate 2 — Cross-check (SBI-vs-MCMC at x_obs = config centrals): **PASS**

8/8 params pass — **no soft fails** (tighter than v3, which soft-failed rho_core). ESS_MCMC 4179.7.

| param | pass | Δmean/σ | σ_SBI/σ_MCMC | shape |
|---|---|---|---|---|
| alpha | ✓ | +0.086 | 1.034 | ✓ |
| log10_zeta_Ih | ✓ | +0.003 | 1.019 | ✓ |
| log10_zeta_sil | ✓ | +0.047 | 1.015 | ✓ |
| log10_eta_Ih | ✓ | +0.036 | 1.008 | ✓ |
| log10_eta_sil | ✓ | +0.016 | 1.001 | ✓ |
| Tb_K | ✓ | +0.020 | 1.031 | ✓ |
| R_core_km | ✓ | +0.008 | 1.062 | ✓ |
| rho_core_kgm3 | ✓ | +0.013 | 1.011 | ✓ |

Report: `reports/crosscheck/crosscheck_report.json` + overlay plot.

---

## Deliverable — zeta split (ice vs silicate) + heating fraction

User-directed: split `log10_zeta` → `log10_zeta_Ih` + `log10_zeta_sil` so a shared Andrade transient-creep amplitude cannot systematically preference silicate heating. **Report the joint posterior (corr) + pre-registered ice-vs-silicate heating-fraction comparison.**

| quantity | reference MCMC | SBI posterior |
|---|---|---|
| **corr(zeta_Ih, zeta_sil)** | **−0.021** | **−0.040** |
| log10_zeta_Ih q16/50/84 | −2.05 / −0.44 / +1.17 | −2.13 / −0.37 / +1.18 |
| log10_zeta_sil q16/50/84 | −1.40 / +0.03 / +1.34 | −1.38 / +0.10 / +1.39 |
| f_ice q16/50/84 | 0.22 / 0.85 / 0.99 | 0.25 / 0.87 / 0.99 |
| median Ice_Ih_W | 455.6 GW | 451.8 GW |
| median Silicate_W | 91.0 GW | 69.9 GW |
| frac ice-dominated (f_ice>0.5) | 0.708 | 0.736 |

**Result: the split works as intended, and SBI reproduces the reference on both the fraction AND the absolute powers.** The two zetas are decoupled (corr ≈ −0.03, consistent between MCMC and SBI). Heating is ice-dominated across the majority of the posterior (~74% SBI, 71% MCMC) with a wide tail toward silicate-dominated.

**Provenance note (reviewer item 1, resolved):** heating MUST be recomputed through the runner-enriched cache (`MCMCRunner._expand_theta` + `forward_model_k2_flexible` on `runner.structure_data`). A first draft of the SBI report loaded the structure cache with a bare `pickle.load`, yielding a wrong forward map (137 GW / f_ice 0.70) and a spurious 4.5σ gap vs the reference. Routing through a real `MCMCRunner` (which enriches the cache via `_load_grid_cache`) is validated to reproduce the stored reference `heating_results` exactly (455.57 GW = 455.57 GW on identical thetas), and the corrected SBI posterior then agrees with the reference. Script `/tmp/galileo_v1p1/sbi_zeta_heating_v2.py`; validation `/tmp/galileo_v1p1/heating_provenance_check.py`.

**CAVEAT:** absolute heating powers still condition on the HYPOTHETICAL k2/h2 channels (no Galileo tidal measurement exists). Additionally (reviewer): both zetas remain weakly constrained (posterior widths ≈ prior widths), so the split's defensible claim is that it removes an *artificial* shared-zeta coupling — not that the data independently pins each zeta; and the f_ice distribution shape is sensitive to the wide U[−3,2] zeta prior. **Carry a zeta-prior-sensitivity check and a derived-quantity (heating-fraction) crosscheck gate into Campaign 2** (reviewer required-validation items 2 and 3).

Reports: `reports/zeta_heating_reference.json`, `reports/zeta_heating_sbi.json`.

---

## Overall: **PASS** (SBC + crosscheck + deliverable)

All gates in the v1.1 suite pass. Status per CLAUDE.md vocabulary: the artifact is **verified** against its reference MCMC under the documented reproduction (SBC + crosscheck reports cited above). Ratification (INDEX row, GUI slot) is Machine A's step.
