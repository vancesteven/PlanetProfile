# Enceladus k₂ campaign — plan (DRAFT, not approved)

Status: `not implemented`. Drafted 2026-11-01; **decisions D1–D6 and plan
corrections C1–C8 APPROVED by the user 2026-11-02** (expert input: Opus
scientific-reviewer 2026-11-02, see "Approved decisions" below). The
pre-registration review is still required before any compute.

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

## Approved decisions (user, 2026-11-02)

- **D1. Alongside, not replacing.** The k2 forecast model sits next to the
  Cassini-only slot, with a "k2 not measured" (very large σ) mode. That mode
  gives the rheology-marginalized Cassini posterior; record its difference
  from pocoMC run (a) as the size of the frozen-rheology understatement of
  the published zb. A material difference sends the Cassini slot back for
  re-adjudication.
- **D2. Sample core rigidity:** log10_mu_sil U[log10 1, log10 50] GPa
  (Bagheri Table 2), as a nuisance. Reviewer toy model: ΔRe k2 is
  0.0009–0.0021 over that range, ~1σ at X-band and 5–10σ at Ka-band. It is
  expected to be prior-dominated (Bagheri p. 13), so exclude it from
  per-parameter SBC pass/fail verdicts.
- **D3. One artifact, σ(k2) as a conditioning input,** sampled
  log-uniformly over a range that spans Ka-band to "not measured". Gated
  separately at the two presets, with SBC stratified by σ. If the Ka-band
  stratum fails, scope-bound the artifact to X-band and use pocoMC for Ka.
- **D4. Heating is display plus a pocoMC ablation only;** never conditioned
  together with Im k2. E = (21/2)(−Im k2)(nR)⁵e²/G ≈ 1.49e12 W × (−Im k2)
  (checked 2026-11-02), so 15–40 GW ≡ −Im k2 0.010–0.027. A pocoMC
  "steady-state" ablation imposes that band, because without it the shell-η
  prior favours the elastic plateau (Im k2 ≈ 0).
- **D5. Compute and store Re/Im h2 and l2 in every simulation row;** use
  them in pocoMC injections. Conditioning the SBI on them is a follow-on
  artifact (lateral heterogeneity, Bagheri §5.4, is out of scope). Needs
  wiring: forward_model_k2_flexible returns no l2 today.
- **D6. Andrade, alpha sampled U[0.15, 0.35],** log10_zeta = 0 (Bagheri's
  β = μ^(α−1)/η^α; the reviewer showed TidalPy Andrade(α, 1) matches it to
  machine precision). Maxwell is a pocoMC ablation only.

## Approved plan corrections (user, 2026-11-02)

- **C1. Rheology-selection bug (CRITICAL, verified).** The log10_mu_Ih hook
  (`forward_models.py:246`) unconditionally sets rheology_type = 'maxwell',
  so with Andrade the result depends on hook order. Fix it with a test
  before any cache or training work.
- **C2. 1σ, not 2σ.** Bagheri's expected uncertainties are 2σ (Table 3 note,
  p. 13, verified). Use 1σ(k2) = 0.001 (X-band) and 0.0001 (Ka-band); lag
  2° (0.2°).
- **C3. Shell-η prior U[12, 20].** Bagheri's discriminating analysis uses
  η_shell > 1e12 (p. 10); k2 is flat above ~1e17. Core η stays U[9, 22]
  (Bagheri uses > 1e11 for discrimination; record why if widened).
- **C4. Default hypothetical x_obs:** Bagheri Table 3 favoured values as
  two labelled presets: core-dissipative |k2| 0.024 at 39°, and
  shell-dissipative 0.032 at 29°. Not a prior-dependent predictive median.
- **C5. Observable −Im k2 (signed), not |Im k2|.** A sign-violating
  simulation is flagged, not folded. Document the (|k2|, lag) ↔
  (Re, −Im) mapping.
- **C6. Deformable libration uses only Re(y1).** Validate it, or restrict
  it, where the lag is large (Andrade η ≲ 1e14, lag 23–47°).
- **C7. Gate P1b:** after the K_int repair, the deformable branch in the
  elastic limit (η → ∞) reproduces the frozen Cassini libration. Gate P3 is
  already met (TidalPy Andrade at ζ = 1). P4 uses Bagheri's Fig. 2
  fixed-parameter case rather than the broad Table 3 ranges.
- **C8. Assumption change, recorded:** sampling μ_Ih U[2, 4] GPa overrides
  PP's EOS ice rigidity (3.5–4.7 GPa at the as-built nodes). That also
  affects the Cassini channels via the deformable libration.

The parameter table above is superseded where it conflicts: shell η
U[12, 20] (C3), alpha sampled (D6), core rigidity sampled (D2), σ(k2) as a
conditioning input (D3).

## Order of work

1. Binding review of the ALMA fix (in progress). Commit it.
2. K_int repair, as its own reviewed change (it moves the Cassini-only
   libration channel too; that has to be re-evaluated and recorded).
3. Gates P1–P4.
4. ~~User decisions D1–D6~~ (approved 2026-11-02); C1 fix first, then
   the reviewer pre-registers.
5. pocoMC references, pushforward and injections.
6. SBI training, gates, GUI slot and deploy.
