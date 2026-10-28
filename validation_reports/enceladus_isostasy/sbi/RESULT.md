# Enceladus ocean-branch SBI artifact: gate results (2026-10-27) — **FAIL (G3, G5)**

**Artifact**
- `~/scratch/enc_sbi/enceladus_ocean6D_rhorock_posterior_1m.pt`, sha256 `43255b3c…`.
- nsf at deployed defaults, train seed 75, 1,000,953 rows, converged at 91 epochs in 93 min.
- Config hash `8aeb033e` in sbi mode. The same file in mcmc mode hashes to `afcb9901`; the difference is only the mode key.
- Dataset sha `8b07744f…`. G0 PASS (`G0_report.json`).

**Gates and verdict**
- Gates: `plans/scripts/enceladus_sbi_gates.py` (53f333fa, committed before training finished). Full numbers in `gates_report.json`.
- Preregistered outcome applied: **record it, do not deploy as validated, escalate. pocoMC run (a) stays authoritative.** Nothing has been retrained or re-tuned.

| gate | verdict | key numbers |
|---|---|---|
| G0 data integrity | PASS | servable 0.8930; noise sd/σ_infl 0.999–1.001; corr 0.472; likelihood path exact |
| G1 crosscheck at x_obs (6 params, n_ks 32,272) | **PASS** | every parameter far inside tolerance (D-excess ≤ 0.009 vs ≈ 0.019; σ-ratio 0.995–1.036) |
| G2 pushforward + joint at x_obs | **PASS** | medians ≤ 0.047 sd; KS ≤ 0.035; observable medians ≤ 0.056 σ_infl; **joint Δcorr 0.0003** (gate 0.05), **energy 0.0047** (gate 0.0123); P(C2 > 0.95) 0.066 vs 0.062; out-of-box leakage 0.3 % |
| G3 global SBC (2,000 held-out pairs) | **FAIL** | per-parameter coverage68 0.673–0.694, all in [0.65, 0.71]; joint HPD-rank KS p = 0.82; **log10_w rank KS BH-p = 0.014 (< 0.05)**; the other 5 parameters have BH-p 0.12–0.98 |
| G4 flow vs pocoMC on 6 off-nominal data sets | **PASS** | all 18 cells: median gap ≤ 0.066 pocoMC-sd, sd ratio 1.015–1.048 (incl. noise2/noise3 extremes) |
| G5 hole leakage (thin-shell pairs) | **FAIL** | per-pair flow mass in hole cells: mean **7.9 %** (gate 1 %), p95 **20 %** (gate 5 %), over n = 1,028 pairs; at x_obs: 0.0 % |
| query guard at x_obs | no label needed | 0.0 % of 10k draws at zb < 15 km or in a hole cell |

The thin-shell G3 stratum (zb₀ < 15 km, n = 353, not gated) shows coverage 0.66–0.73, HPD p = 0.054, and C2 BH-p 0.023.

## Manager's initial reading (superseded by the adjudication below)

My hypothesis was that the log10_w G3 defect was the hole seen through the prior. **The reviewer REJECTED this hypothesis.**

## Reviewer adjudication (Opus, 2026-10-27): SCOPE-BOUND deploy; G3 and G5 remain recorded FAILs

- **G3 (log10_w, BH-p 0.014) is not evidence of a flow defect.**
  - log10_w is prior-dominated, so its SBC rank equals the effective-prior CDF of the held-out truth (corr 0.984).
  - The KS of F_eff(w₀) alone, with no flow involved, gives p = 0.0026, against 0.0023 for the ranks. The held-out seed-900 sample is itself non-uniform in w (raw-prior KS p = 0.011).
  - The defect is absent in the thin-shell stratum (p 0.32) and in leaking pairs (p 0.77). It persists for zb₀ ≥ 25 km (p 0.011) and after hole truncation (p 0.0052).
  - **The hole hypothesis is rejected.**
  - The manager reproduced the key numbers independently: KS of F_train(w₀) p = 0.0026; held-out vs training 2-sample p = 0.0026; mass below log10_w −0.5 is 0.214 held-out vs 0.184 training. Among 300 random training blocks of 2,000, 5.3 % have p < 0.05 and 0.3 % have p ≤ 0.0026, so this is a rare but legitimate null draw.
- **G5 is a support defect of the smooth NSF across the jagged (zb, w) cache hole.**
  - Leakage tracks posterior zb width, not w.
  - In the thin-shell stratum the flow is also mis-shaped: zb coverage68 0.734 with a ∩-shaped rank histogram (over-dispersed), and C2 rank p 0.004. Truncation does not fix this.
- **The artifact is NOT validated globally and NOT valid for thin-shell (zb < 15 km) data.**
- It reproduces pocoMC at x_obs (G1, G2) and on 6 off-nominal injections (G4).
- **pocoMC run (a) remains authoritative for all quoted results.** Nothing was retrained or re-thresholded.

**Findings recorded for future preregistrations:**
- **MAJOR.** The per-parameter SBC KS has no power against flow error for prior-dominated parameters: it reduces to a uniformity test of the held-out θ₀ sample. New preregistrations must report a null-control KS of F_eff(θ₀) on the same pairs and gate conditional on it, or use stratified or inverse-CDF held-out sampling.
- **MAJOR.** Thin-shell queries are invalid. Handled by the refusal guard.
- **Forbidden under 08659c67:** re-running G3 on another held-out seed (data-shopping), and changing capacity or N.
- **(c) Retraining** is permitted only under a NEW preregistration that changes the estimator class or support handling and carries the G3 null-control fix.

## Deploy scope (GUI and model card text, verbatim)

> **Enceladus ocean-branch SBI (scope-bound).** Validated against pocoMC at the Cassini datum and on six off-nominal injection data sets. It is not globally validated: preregistered gates G3 (log10_w rank uniformity, diagnosed as held-out sampling fluctuation) and G5 (thin-shell cache-hole leakage) FAILED.
> - Queries whose posterior puts more than 1 % mass at zb < 15 km or in a cache-hole cell are **refused**. Use pocoMC for those.
> - Below that threshold, draws falling in hole cells are dropped and the fraction is shown.
> - Results at other x are "effective-prior posterior; servable cells only; not individually validated".
> - Cite pocoMC run (a), not this artifact, for published values.

**Conditions:**
1. pocoMC is authoritative.
2. The guard is a REFUSAL tier at the 1 % threshold.
3. Drop-and-disclose below the threshold is support enforcement, not a gate fix: G5 stays FAIL and is never re-evaluated on truncated draws.
4. "Validated" applies only to x_obs and the 6 G4 sets.

**Required before the GUI is enabled:**
- **Refusal-guard test: `verified`.** `PlanetProfile/Inference/enceladus_sbi_guard.py` (`guarded_posterior`, with `SCOPE_TEXT` carrying the verbatim scope text) is tested by `tests/enceladus_sbi_guard_test.py`, 3 passed:
  - x_obs is accepted and validated, with 0 % dropped and 0 % thin;
  - the thinnest thin-shard x is refused (more than 1 %, "use pocoMC");
  - an off-datum x is accepted with the not-validated label.
- **Scope text on the model card and in the GUI: NOT YET.** The artifact is not wired into the GUI. Until it is, the status is `implemented, unverified` for GUI use.
