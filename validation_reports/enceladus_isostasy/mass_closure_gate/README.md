# Mass-closure jitter gate — 2026-10-22

Run: `plans/cluster/enceladus_mass_closure_gate.sh 6` on the Mac (commit ad03a68b). 90/90 builds. Overall **FAIL**.
Details are in `jitter_gate.md` and `jitter_gate.json`.

| | σ_model / σ_obs (max(std, \|bias\|)) | limit |
|---|---|---|
| C20 | 1.27 | 0.5 |
| C22 | 0.95 | 0.5 |
| libration | 0.107 | 0.1 |

## Scientific-reviewer ruling (2026-10-22)

**Verdict.** PASS WITH CONCERNS on the closure. BLOCK on using C20/C22 as they stand.

**What the jitter is.** After closure, the leftover jitter is quantisation of R_sil:
- With ρ_sil closed, C_hyd depends only on R_sil.
- R_sil candidates sit on the hydrosphere layer grid, about 0.79 km per step at deltaP/1.
- C20 lies on a lattice of steps: 0.577σ at /1, 0.29σ at /2 and 0.144σ at /4.
- C20 and C22 are a single mode (correlation −1.000, ratio −1.34, which matches the hydrostatic 10/3 prediction).
- Libration follows the same mode, and C30 does not move.

**The finest build is not converged.**
- φ choice alone moves it by 0–0.58σ at deltaP/4.
- nPoros 8 is a nested subset of 16, while 12 is not nested.
- So max|prod − finest| overstates the error. The production error across nodes is about zero mean, with RMS ≈ 0.57σ (C20) and 0.43σ (C22).

**Closure ON vs OFF.** ON has the smaller spread at 8 of 10 nodes. At (25,30), OFF benefited from an accidental cancellation with the mass deficit.

**Rulings.**
- (a) **Keep closure ON.**
- (b) **Recommended remedy: Option 4.** A runner-level closure of the reduced stack in two constraints, solved in closed form for (R_s, ρ_s) so that M = Bulk.M and C/MR² = Cmeasured:
  - a = 3M/4π − ρ_o R_o³ − ρ_i (R³ − R_o³)
  - b = (15/8π) C M R² − ρ_o R_o⁵ − ρ_i (R⁵ − R_o⁵)
  - R_s² = b/a, and ρ_s = ρ_o + a/R_s³
  - It needs no rebuild and no PP core edit. **It requires user approval**, because it makes MoI an exact equality in the gravity and libration forward model and overrides R_sil by about 1–2 km.
  - Option 2 (a finer rebuild) is rejected: roughly 4–7× the r7 cost with convergence unproven.
  - Option 3 is unnecessary.
  - Fallback if the user declines: a rank-1 σ_model along the measured mode with s = RMS (C20 ≈ 2.1e-5), added to the full covariance Σ_obs(ρ) + s² u uᵀ, not as diagonal inflation.
- (c) Libration at 0.107σ is an acceptable marginal fail (RMS 0.044σ).

**Validation required for Option 4.** Re-run the gate on the same nodes and builds. It passes if:
- the C20/C22 spread is ≤ 0.1σ;
- the libration spread is ≤ 0.05σ;
- the frozen residual stays at about 1e-15;
- every node is feasible (a, b > 0; R_s < R_o; ρ_s inside the cache's rock-density range);
- the report gives R_s − R_sil(PP) and stack MoI vs PP CMR2 per node;
- the mean offset from the deltaP/4 ON values is ≤ 0.3σ;
- the posterior-median check is run;
- the local-quadratic check is re-run.
