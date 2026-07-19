# Europa Clipper v3 (8D + ocean salinity) — Phase-3 gate report

Machine B, 2026-07-19. Artifact `europa_seawater_v3_clipper_8D_posterior_1m.pt`
(nsf, 336,537 kept training sims, seeds train 42 / data 48 / noise 4848).
Reference MCMC `Test52_seawater_v3/europa_seawater_v3_reference_result.pkl`
(pocoMC, 4419 samples). Both use the shared bilinear (Tb, log10 w) interpolant,
so the SBI training support == reference-MCMC support by construction.

## Headline science result
**Strong Tb↔salinity degeneracy, faithfully reproduced by the flow.**
- corr(Tb, log10_w): MCMC −0.986, SBI −0.988 (the iso-|Ae_synodic| ridge).
- w weakly identified near seawater: MCMC 38.7 ppt [28.2, 48.6],
  SBI 39.4 ppt [28.2, 51.3]; Tb 264.3 K.

## Gate results

| Gate | Verdict | Detail |
|------|---------|--------|
| **2D DEGENERACY (central deliverable)** | **PASS** | corr |Δ|=0.002 ≤ 0.15; 68% containment 0.650, 95% 0.943 (±0.1); all 8 marginals within tol (d/thr ≤ 0.45, σ-ratios 1.03–1.10) |
| **SBC (8 params)** | **PASS** | 497 held-out pairs, min KS p = 0.118 ≥ 0.05 |
| **crosscheck** | **FAIL (soft, 1/8)** | 7/8 params full PASS; `rho_core_kgm3` shape-only fail: KS D 0.0565 > tol 0.0507 (mean/median/σ all PASS: mean_diff 68.6 ≪ tol 198, σ-ratio 1.09) |
| **limits (Im_k2 → eta_Ih monotone)** | **FAIL (soft, premise)** | medians flat [13.13, 13.12, 13.20, 13.19] over |Im_k2| 0.15–0.30; containment 0.985 |

## Interpretation of the two soft failures (NOT gate-tuned)

Scientific-reviewer verdict (2026-07-19): **PASS WITH CONCERNS — RATIFIABLE.**
Both soft failures **NON-BLOCKING** (documented caveats). Refinements below are
the reviewer's corrections to the mechanism.

**crosscheck `rho_core_kgm3` shape → NON-BLOCKING.** `rho_core` is the
least-identified parameter — posterior σ 793 ≈ 26% of its [5000,8000] prior
width, i.e. nearly prior-dominated (CMR² + mass conservation weakly constrain
core density given a sampled R_core). The KS *p-value is informational only* in
this gate; the fail is the D-statistic (0.0565) exceeding the Monte-Carlo self-D
floor (0.0338) by ~11% — a *small but genuine* shape offset beyond MC noise, NOT
p-value hypersensitivity. It is non-blocking because: (i) mean (Δ=69 ≪ 198 tol),
median, and σ-ratio (1.09) all pass — location and scale agree; (ii) SBC passes
rho_core cleanly (KS p 0.284, c2st 0.543 ≈ 0.5) — the flow's *calibration* over
the prior-predictive is sound, so the offset is local to this one fiducial;
(iii) the mismatch is **conservative** — SBI is *wider* (σ 866 vs 793), so
credible intervals do not over-state the core-density constraint; (iv) rho_core
is not the v3 deliverable, and the degeneracy gate passed it at d/thr 0.45.
(Minor systematic: SBI is over-dispersed vs MCMC across all 8 params, σ-ratio
1.03–1.10, largest at rho_core — expected, safe NSF behavior.)

**limits monotonicity → NON-BLOCKING (inapplicable premise).** The gate's premise
(eta_Ih decreases as |Im k2| rises) is a Titan-regime, ocean-free relationship;
in ocean-bearing v3 the induction channels + |Ae_synodic| support cut dominate
and eta_Ih is decoupled from Im_k2 at the fiducial (medians flat to 0.083 dex ≈
0.05σ against the 1.67-dex eta_Ih marginal — statistically indistinguishable from
flat at n=2000). The containment 0.985 sub-fail is a *second symptom of the same
issue*: the sweep conditions the flow at Im_k2 = 0.15–0.30 while the real
observation is Im_k2 = 0.0, so the flow is evaluated off its conditioning in a
sparse-training region where ~1.5% NSF boundary leakage is expected — it does NOT
reflect containment at the real x_obs (the degeneracy gate's frac68/frac95 at the
true fiducial are within tol). This gate is best relabeled "not run / inapplicable
for ocean-bearing bodies," not "soft fail."

A v3-appropriate limiting test was probed (sweep the dominant orbital channel
`Bind_orbital_x_imag`, expect log10_w to rise): w barely moved (36.5→36.0 ppt,
−0.005 dex) — confirming that a SINGLE w-carrying channel varied in isolation
carries almost no salinity information (the identification is weak and *joint*,
per the −0.986 Tb↔w degeneracy). Per the reviewer, the correct future
falsification test is to sweep the synodic induction AMPLITUDE and verify ocean
**conductance** (σ×thickness, the quantity the induction likelihood actually
constrains) rises monotonically — NOT to sweep w-carrying channels expecting a
monotone w. Adding that test is future work; its absence does not block this
artifact because SBC (global calibration) and the degeneracy gate (fiducial
joint) already exercise the induction channels. See limits_induction.json.

The two ratifying gates for the v3 salinity deliverable — the 2D degeneracy gate
and SBC — both PASS. The two parameters carrying the v3 degeneracy
(log10_wOcean_ppt and Tb) pass every crosscheck criterion cleanly.

## Pre-training gates (see plans/europa-clipper-v3-phase1-validation.md)
C4 transition-cell jump 0.000σ PASS; edge Tb-discretization 0.0851 K PASS;
edge w-interpolation 0.285 K seawater / 0.099 K worst-case (reviewer
GO-WITH-CONDITIONS Option A, documented as config
`salinity_w_interp_caveat_2026_07_18`).

## Product fixes made during Phase 3 (committed)
- `UnservableSampleError` reject for unbuilt tilted-band corners in the MCMC
  likelihood + SBI dataset generator (commit 60f22bb2).
- `_structure_R_body_km` + run() post-processing made 2D-cache-safe (aab41fd7).

## Dual-space + prior sensitivity
Reported above in both log10 w and linear w. The w-posterior is prior-shape
sensitive (log-uniform Jeffreys prior; weak identification), as pre-registered;
the near-seawater mode is driven by the orbital Bind channels against the |Ae|
support cut, and carries the documented 0.285 K support-edge w-interpolation
bias. Do not over-interpret the w mode to better than ~0.3 K / ~0.05 dex in Tb–w.
