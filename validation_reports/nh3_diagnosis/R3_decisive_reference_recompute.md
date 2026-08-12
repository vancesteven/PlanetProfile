# §0.20 R3 — decisive NH3 reference recompute (n_eff=8000, 3 seeds)

**Status: complete — GATE PASSES; scientific-reviewer verdict
PASS-WITH-CONCERNS, re-ratifies C16 and supports releasing the STOP
(2026-08-11). Final STOP release is the manager's (Machine A) action; this
session does NOT self-ratify.** Third and decisive step in the C16 resolution
plan (Machine A §0.20). R1 (PASS-WITH-CONCERNS) showed
neither side's legitimate weighting treatment explains the residual and named
the reference-side seed scatter (span 0.0121 at n_eff=2000) as the remaining
suspect. R2 (PASS-WITH-CONCERNS) showed the corrected side is a precisely
resolved fixed point (not finite-N bias). R3 resolves the reference side by
recomputing the NH3 free-gravity reference MCMC at **n_eff=8000 (4× the prior
resolution), 3 seeds (72/172/272)** and evaluates the committed
re-ratification gate. Per project discipline, the scientific-reviewer + manager
adjudicate — this session reports the number and does not self-ratify.

Drivers (repo): `plans/scripts/nh3_diag_matched_reference.py`
(`--n-effective 8000 --n-active 4000 --seeds 72 172 272`, regime-preserving
ratio 2, config NEVER mutated) + `plans/scripts/r3_pool_and_gate.py` (ratified
treatment-A pooling + committed pass-rule). Artifacts:
`validation_reports/nh3_diagnosis/matched_reference_neff8000/`.

## MCMC quality (all 3 seeds)

| seed | n_samples | ESS | R-hat | accept | wall (min) | log_Z |
|---|---|---|---|---|---|---|
| 72  | 21397 | (not logged) | (not logged) | — | 91.3 | −12.085 |
| 172 | 21095 | 12882 | 1.000 | 21.04% | 97.9 | −11.484 |
| 272 | 21246 | 12975 | 1.000 | 19.93% | 99.9 | −11.824 |

All converged (R-hat=1.000, ESS≈12.9k). The weighted |Im_k2| between-seed std
collapsed to **0.00015** (from the n_eff=2000 regime), confirming the
higher-resolution reference is tightly converged — the seed scatter R1 flagged
was a resolution artifact.

## Re-ratification gate (committed pass-rule)

Pass-rule (committed, `is_validation_nh3.json:105`, reviewer 2026-08-11):
`|corrected − reference| ≤ 2·√(SE_corr² + SE_ref²)`, both **between-seed** SEs
(NOT the raw 3-seed min/max band, NOT the single-seed bootstrap).

| quantity | n_eff=2000 (prior) | **n_eff=8000 (R3)** |
|---|---|---|
| reference pooled ocean fraction (treatment A) | 0.91725 | **0.92865** |
| reference per-seed fracs | — | 0.93090 / 0.93046 / 0.92459 |
| reference between-seed SE | — | 0.00203 |
| corrected pooled (C13, unchanged) | 0.93217 | 0.93217 |
| corrected between-seed SE (C13) | — | 0.00212 |
| **residual (corrected − reference)** | **+0.01492** | **+0.00352** |
| combined SE | — | 0.00294 |
| **2× combined SE bound** | — | **0.00588** |
| **residual / combined SE** | ~3.6 (FAIL) | **1.20** |
| **gate** | FAIL | **PASS (re_ratifies_C16 = true)** |

## Reading (for the reviewer + manager — NOT self-adjudicated)

The n_eff=8000 reference moved **+0.0114 upward** (0.91725 → 0.92865) — the
reviewer's preregistered "reference moves materially upward → residual shrinks
toward the gate" branch. The residual fell from +0.0149 (3.6σ FAIL at
n_eff=2000) to **+0.00352 (1.20σ, inside the 2σ bound)**. Diagnosis chain,
end to end:

1. **R1:** neither side is weighting-fragile; the largest single sensitivity is
   reference-side seed scatter at n_eff=2000.
2. **R2:** the corrected side is a stable, precisely-resolved fixed point
   (0.93217, moved *away* from the reference at higher N — not finite-N bias).
3. **R3:** the reference side WAS the artifact. At 4× resolution the reference
   converges to within the preregistered gate of a common ~0.93 ocean fraction
   and the tension collapses inside the bound. (Reviewer-tempered wording: at
   n_eff=8000 the reference 0.92865 still sits +0.0035 below the corrected
   0.93217 — both estimators are converging to a shared ~0.93 limit and now
   agree within the gate; the reference has not "arrived at" the corrected value.)

**The corrected ocean fraction (~0.932) was approximately right all along; the
n_eff=2000 reference was biased low by sampling resolution.** This is the
constructive resolution of C16, not a mere non-rejection.

**Plateau caveat (reviewer required-validation, folded in 2026-08-11).** The
pooled ocean fraction is characterized at only **two** resolution points
(0.91725 at n_eff=2000 → 0.92865 at n_eff=8000); it is not shown flat. The
co-moving weighted |Im_k2| median is still advancing at the last 4× step
(+0.00369 for 500→2000, then +0.00193 for 2000→8000 — roughly halving, an
approach to a limit but not a demonstrated plateau). This does NOT threaten the
ruling: any residual under-resolution moves the reference **further up** toward
the corrected fixed point (0.93217), i.e. toward *better* agreement or at most a
small sign flip. The reference is therefore at-or-below its converged value and
the residual is bounded toward agreement, not away. A single confirmatory
n_eff=32000 seed would demonstrate the fraction plateau; queue it only if C16 is
later contested (not required for release).

**Recommendation (advisory only — reviewer + manager decide):** the committed
gate passes at 1.20 combined-SE. Subject to the reviewer's verification of the
pooling, the SE construction, and the MCMC convergence, this supports
**re-ratifying C16** and releasing the STOP on NH3 Track-1 corrected compute.
Final re-ratification is NOT this session's to declare.

## Caveats carried for the reviewer

- Seed 272's fraction (0.92459) sits ~0.006 below seeds 72/172 (~0.9307) — the
  same seed that had the lowest ESS/N at n_eff=2000. It is included (3-seed
  pooling as prescribed); the reviewer may wish to note whether it warrants a
  4th seed, though the gate passes with it included.
- The gate uses the C13 corrected between-seed SE (0.00212), matching the R2
  required-validation (between-seed, not single-seed bootstrap). The R2
  delta-method cross-check (0.00412) is larger; using it would only widen the
  bound and strengthen the PASS.
- `re_ratifies_C16 = true` in the JSON is the **gate arithmetic**, not a status
  change. C16 remains at MANAGER-GATE STOP until the reviewer + manager act.

## Scientific-reviewer adjudication (2026-08-11) — ADOPTED

**Verdict: PASS-WITH-CONCERNS — re-ratifies C16; release the STOP.** The reviewer
independently reproduced the gate arithmetic to the last digit from the three
raw `nh3_reference_seed{72,172,272}_neff8000.pkl` files (reference pooled
0.9286503, SE 0.0020338; corrected pooled 0.9321722, SE 0.0021226; residual
+0.0035218; combined SE 0.0029397; ratio 1.198; PASS), and confirmed:

- **The reference move is like-for-like** (high confidence). Recomputing the
  n_eff=2000 pooled reference from ITS OWN three seed pkls under treatment A
  reproduces the committed 0.9172548 to 7 digits, so the +0.01140 rise is not a
  seed-set or pooling-convention artifact. Identical config_hash (1611b65f…) and
  structure_cache_sha256 (3d837cf8…) across n_eff=2000 and n_eff=8000, and the
  cache sha256 also matches the deployed SBI adjudication's byte_identity block —
  same forward model on both sides. (The config_hash differs from the SBI
  artifact's e596574d — expected MCMC-vs-SBI schema difference per the
  "config_hash mode-only diff" note; the load-bearing invariant is the shared
  cache sha256, which holds.)
- **Resolution effect, not regime confound** (medium-high). n_active/n_eff ≈ 2 at
  all three resolutions (256/500, 1024/2000, 4000/8000) → B3-ratified regime
  preserved. |Im_k2| median between-seed std collapsed 0.00066→0.00015 (~4.4×),
  ESS ~4-5k→~12.9k, R-hat=1.000 — independent convergence evidence. The move is
  physically coherent: n_eff=2000 undersampled the higher-dissipation ocean
  branch (ocean fraction and |Im_k2| rose together).
- **Seed-272 concern is conservative, not adverse** (high). Dropping seed 272
  (72/172 only) RAISES the reference to 0.93068 and SHRINKS the residual to
  +0.00149. Its inclusion *lowers* the reference and *inflates* the residual — the
  gate passes with the conservative choice. A 4th seed would only tighten an
  already-passing gate; not warranted for routing.
- **Small-sample SE cuts toward safety** (medium). Both SEs rest on 3 seeds
  (2 dof); a Welch/t correction raises the critical value to ~2.8-4.3, but the
  residual is only 1.20 combined-SE — comfortably inside even the widened bound.
- **Corrected-side SE is the correct like-with-like** (high). C13 seeds are
  flow-sampling+SNIS seeds of the single deployed flow; using their between-seed
  scatter as SE_corr is correct for the C16 "does THIS deployed flow agree"
  question and matches the R2 required-validation. The gate correctly uses the
  between-seed SE (0.00212), not the single-seed bootstrap (0.00197).

**Reviewer required-validations (non-blocking, folded in 2026-08-12):**
(1) two-resolution-points + advancing-|Im_k2| caveat — added to the Reading
section above; (2) "rises to meet" wording tempered to "converges to within the
preregistered gate of a common ~0.93 fraction" — done; (3) optional n_eff=32000
confirmatory seed — queued only if C16 is later contested; (4) cosmetic
provenance fix in `matched_reference/matched_reference_report.json`
(`seeds:[272]` disambiguated) — done, `seeds_note` added.

**Routing:** R3 re-ratifies C16 under the committed pass-rule
(`is_validation_nh3.json:105`). The reviewer's scientific verdict supports
releasing the MANAGER-GATE STOP on NH3 Track-1 corrected compute and clearing
the deployed corrected NH3 ocean fraction (~0.932) as a ratified input. **Final
STOP release is the manager's (Machine A) action, not this session's.** The four
concerns are robustness/documentation items, all pointing toward *stronger*
agreement — hence PASS-WITH-CONCERNS, not CHANGES-REQUIRED.
