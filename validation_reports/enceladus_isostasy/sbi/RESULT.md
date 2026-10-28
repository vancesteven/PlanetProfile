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

## Reading (manager; NOT self-adjudicated, routed to the reviewer)

- **At and around the Cassini datum the artifact reproduces pocoMC.** That covers G1, G2 (including the joint ρ_s–C2 ridge) and G4. G4 includes the two extreme noise draws that tripped the headline review.
- **G5 is a support defect on the thin-shell edge.** The flow's density is smooth over the (zb, w) box. It cannot represent the hole, a jagged, grid-quantised region at zb 5–13 km, so for thin-shell data about 8 % of its mass lands in unbuildable cells. At x_obs this is 0 %.
- **G3 fails on log10_w alone.** The coverage is fine, and the joint HPD rank is uniform. log10_w's prior-effective distribution is tilted by the hole, whose fraction depends strongly on w. A flow that leaks into the hole would mis-rank w specifically. **Hypothesis:** the w rank defect is the same hole effect seen through the prior. This is untested.
- **Consequence for use.** As preregistered, the artifact is NOT deployable as "validated". The pocoMC posterior stays the result.
