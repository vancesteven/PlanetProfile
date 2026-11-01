# Enceladus k₂ campaign — plan (DRAFT, not approved)

Status: `not implemented`. Drafted 2026-11-01 for user approval and
scientific-reviewer pre-registration review. No compute has been run.

## Why

- **User direction (2026-10-31):** "we want k2 for Enceladus because future
  missions will measure it". This reverses the 2026-08-12 ruling that
  omitted k₂ (freeze doc point 1), and with it blocker B6, which dropped the
  rheology parameters because no conditioned channel depended on them.
- **Basis:** Bagheri et al. 2025 (PSJ 6:245, `papers/bagheri2025exploring.pdf`).
  - Predicted |k₂| is 0.024 ± 0.010 (core-dissipative) and 0.032 ± 0.013
    (shell-dissipative), 3σ (Table 3, p. 13).
  - Future tracking recovers k₂ to < 0.002 with X-band DSN and < 0.0002 with
    Ka-band (p. 17).
  - k₂ alone constrains the total dissipation, not its radial distribution.
    The amplitude and phase of h₂/l₂ discriminate between core and shell
    dissipation (p. 2, p. 9).
- **Physics coupling found 2026-10-31**
  (`validation_reports/enceladus_isostasy/libration_rheology_bound/`):
  - Shell rheology moves the libration by up to ~6 % (Andrade, or a soft
    base under an elastic lid). That is up to ~2 km in zb, comparable to the
    current posterior 1σ.
  - k₂ is strongly rheology-sensitive: |k₂| 0.012 → 0.075 over the plausible
    range at 26.7 km.
  - A k₂ channel therefore constrains the term that is now a systematic on
    the published shell thickness.

## Deliverable

An amortized model, one GUI slot, that conditions on the current Cassini
observables (C₂₀, C₂₂, C₃₀, libration) plus **k₂ (Re, |Im|) as a
future-mission observable**.

- **Conditioning value:** x_obs for k₂ is free (the user enters a
  hypothetical measurement). The default is the posterior-predictive median
  at the Cassini datum, not a fabricated "observed" value.
- **σ(k₂) presets:** two, matching Bagheri p. 17: 0.002 (X-band) and 0.0002
  (Ka-band). Each preset is a separately trained artifact, or one artifact
  trained across a σ range (decision D3).
- **Scope:** this is a mission-design / forecasting model. It is NOT a
  re-analysis of Cassini data: Cassini measured no k₂. The current
  scope-bound ocean-branch model stays the Cassini-only slot. GUI wording:
  "k₂ is a hypothetical future measurement".

## Parameter space

The base is Bagheri Table 2 (p. 9), adapted to PP's self-consistent
structure. Sampled log-uniformly where Bagheri does.

| Parameter | Prior | Basis |
|---|---|---|
| zb_km | U[5, 45] (as now) | current campaign |
| log10_wOcean_ppt | U[−1, 2] (as now) | current campaign |
| rho_ice_kgm3, rho_rock_kgm3, compensation_C2 | as now | current campaign |
| log10_eta_Ih (shell effective viscosity, Pa·s) | U[10, 25] | Bagheri η_shell 10¹⁰–10²⁵ |
| log10_mu_Ih (shell shear modulus) | log U[2, 4] GPa | Bagheri μ_shell 2–4 GPa |
| log10_eta_sil (core viscosity, Pa·s) | U[9, 22] | Bagheri η_core 10⁹–10²² |
| core shear modulus | log U[1, 50] GPa | Bagheri μ_core (new parameter; see D2) |
| Andrade α (shell and core) | fixed 0.2, β ≈ μ^(α−1)/η^α | Bagheri p. 21 (α = 0.2, β = 1; Castillo-Rogez 2011) |
| libration_sys_frac | re-derived (see Physics, step 2) | replaces the frozen rigid-branch budget |

Notes:
- **Fixed by PP's self-consistency, not sampled:**
  - Ocean density and thickness: sampled by Bagheri, but set here by
    (zb, w) and mass closure.
  - Core radius: set by ρ_rock and mass closure.
  - Bulk moduli: Bagheri fixes them (shell 10 GPa, core 100 GPa); PP's EOS
    supplies them.
- **Heterogeneous shell:** Bagheri treats a laterally uniform shell, as here.
  Their heterogeneous-shell results are out of scope.
- **Heating sanity band:** a declared display check against Bagheri's
  E_total 15–40 GW (p. 10), not a conditioned observable, unless the user
  rules otherwise (D4).

## Physics and solver changes (before any cache build)

1. **ALMA viscosity units.**
   - Fixed in the working tree 2026-11-01 (`Gravity.AlmaViscosity`;
     `tests/alma_viscosity_units_test.py`). The fix is pending binding review.
   - The inference forward model uses TidalPy
     (`forward_models.forward_model_k2_flexible`, η in Pa·s, SI frequency)
     and is NOT affected. The campaign uses that path.
2. **Libration on the deformable branch.**
   - With rheology sampled, libration must use the same rheology as k₂: the
     Van Hoolst 2013 deformable branch, with y₁ from the SAME TidalPy solve.
     The rigid branch plus a frozen correction is not valid once rheology
     varies (freeze ruling: the +0.8 % correction is "VOID if shell rheology
     changes").
   - **Prerequisite:** the reviewed rigid-branch K_int repair (missing 8π/15
     plus the pressure-figure linearization and elastic terms, which must land
     together; pinned by `tests/librations_test.py`, strict xfail). Until it
     lands, the deformable branch's absolute normalization is unverified.
   - **Gate P1:** the deformable branch reproduces Van Hoolst et al. 2013 /
     2016 published values for an elastic shell to a stated tolerance, AND
     the strict-xfail physical-direction test flips to pass (elasticity
     lowers libration).
3. **TidalPy y₁ for libration.** TidalPy returns the radial functions. Wire
   y₁ through `forward_model_k2_flexible` and the Van Hoolst Eq. 24
   partial-Love-number sum, as `test_partial_love_numbers_reproduce_alma_k2`
   does for ALMA.
   - **Gate P2:** the TidalPy y₁ partial sum reproduces TidalPy k₂ to 1e-6.
4. **Andrade β convention.** Bagheri uses β ≈ μ^(α−1)/η^α. PP maps Andrade
   ζ between pyALMA and TidalPy (`Gravity._pyalma_zeta_to_tidalpy`).
   - **Gate P3:** match the convention to Bagheri's at one parameter point
     (k₂ within 1 %).
5. **Cross-check against Bagheri.**
   - **Gate P4:** at Bagheri's shell-dissipative and core-dissipative
     posterior parameter points (Tables 5/6), PP's |k₂| and k₂ phase land
     inside Bagheri's Table 3 ranges.

## Structure cache

- **Today:** the (zb, w) production cache has structure only. Rheology is
  applied per draw by the TidalPy forward model (`apply_parameters` hooks:
  `log10_eta_Ih`, `log10_mu_Ih`, `log10_eta_sil`, `alpha`/`zeta`), so a
  rheology axis needs NO cache rebuild.
- **Needed:**
  - Core shear modulus is not a registered parameter
    (`parameter_registry.py`). Add a `log10_mu_sil` hook, with a test.
  - Per-draw cost: one TidalPy solve (~1–3 ms at these node counts, measured
    2026-11-01) plus the libration. Training simulation cost is dominated by
    the solve: ~10⁶ draws ≈ 1 h on 10 workers.

## Inference campaign (pre-registered before running)

1. **Reference posterior, pocoMC,** at the Cassini datum PLUS one
   hypothetical k₂ value per σ preset (X-band, Ka-band), n_eff 8000, 2 seeds.
   This is the authoritative result.
2. **Prior pushforward:** Bagheri-range check. Does the prior predictive
   |k₂| cover 0.01–0.05, and the heating 15–40 GW?
3. **Injection-recovery** at Bagheri's two regimes (shell- vs
   core-dissipative). Can a k₂ measurement at each σ separate them? This is
   the science deliverable: Bagheri says k₂ alone cannot (p. 2). Test that
   claim in PP.
4. **SBI training:** 1M simulations, nsf, effective-prior support as now.
5. **Gates G1–G5,** same battery as the ocean-branch artifact, plus:
   - G6, k₂-sector crosscheck, flow vs pocoMC: Re k₂, |Im k₂|, η_Ih, μ_Ih,
     η_sil within the stock tolerances;
   - G7, libration–rheology joint: corr(zb, μ_Ih) and corr(zb, η_Ih) within
     0.05 of pocoMC.
6. **Scope guard** on the same pattern as the current artifact, if needed.

## Open decisions for the user

- **D1.** Confirm the deliverable is a forecasting model (hypothetical k₂)
  alongside the Cassini-only slot. Alternatively, should it replace that slot?
- **D2.** Add core shear modulus as a sampled parameter (Bagheri samples it),
  or fix it?
- **D3.** σ(k₂): two artifacts (X-band, Ka-band), or one artifact trained over
  a σ range?
- **D4.** Heating: display check only, or a conditioned observable
  (15–40 GW)?
- **D5.** Add h₂/l₂ as optional future observables now? Bagheri: they are the
  discriminating channels.
- **D6.** Andrade (Bagheri) or Maxwell? Proposed: Andrade α = 0.2, with a
  Maxwell ablation.

## Order of work

1. Binding review of the ALMA fix (in progress). Commit it.
2. K_int repair, as its own reviewed change (it moves the Cassini-only
   libration channel too; that has to be re-evaluated and recorded).
3. Gates P1–P4.
4. User decisions D1–D6, then the reviewer pre-registers.
5. pocoMC references, pushforward and injections.
6. SBI training, gates, GUI slot and deploy.
