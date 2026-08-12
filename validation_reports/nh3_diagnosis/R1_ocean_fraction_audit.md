# §0.20 R1 — Ocean-fraction weighting audit (ZERO / minimal compute)

**Status: complete — scientific-reviewer verdict PASS-WITH-CONCERNS, proceed to
R2 (2026-08-11).** This is the first, cheapest discriminator in the C16
resolution plan (Machine A §0.20, 2026-08-12). It asks whether the committed
corrected-vs-reference NH3 ocean-fraction residual (**+0.0149**, pooled 3-seed)
can be an artifact of *weighting treatment* on either side. Nothing here
re-ratifies C16; R3 (n_eff≥8000 reference recompute) remains the decisive step.
Per project discipline the scientific-reviewer, not this session, adjudicated
the reading.

**Reviewer ruling (verified against committed artifacts + the `0acff866`
diff):** methodology faithful (byte-identical `log_weights`,
`max_abs_logw_diff=0.0`); neither side's *legitimate* weighting treatment
explains the residual; the −0.0048 top-0.1% ablation is an adversarial deletion
of 20 valid draws (0.33× residual), NOT a candidate estimator, and does NOT
reopen the corrected side; R1 consumes no gate and **C16 status is unchanged**
(remains MANAGER-GATE STOP/REOPEN). Two non-blocking corrections were required
and are folded in above: (1) treatment B demoted as ill-defined; (2) SNIS audit
scoped to the fiducial residual, not pooled. **Cleared to proceed to R2.**

Scripts (repo-tracked copies): `/tmp/r1_reference_audit.py` (zero-compute,
reference seed pkls only) and `/tmp/r1_snis_audit.py` (one flow conditioning at
the fiducial). JSON outputs `/tmp/r1_reference_audit.json`,
`/tmp/r1_snis_audit.json`.

## Reference side (zero compute — committed seed pkls 72/172/272, n_eff=2000)

Per-seed ocean fraction, weighted (pocoMC importance weights) vs unweighted:

| seed | n | Kish ESS | ESS/N | frac (weighted) | frac (unweighted) |
|---|---|---|---|---|---|
| 72  | 7187 | 5218 | 0.726 | 0.92223 | 0.91652 |
| 172 | 7187 | 5214 | 0.725 | 0.91939 | 0.91540 |
| 272 | 6230 | 4142 | 0.665 | 0.91014 | 0.90433 |

Pooled ocean fraction under two legitimate weighting treatments:

| treatment | pooled frac |
|---|---|
| **A — correct pooling** (concat samples+weights, each seed → 1/n_seeds) [committed reference] | **0.91725** |
| C — unweighted pooled (uniform) | 0.91244 |

- **A vs C differ by only 0.0048**, i.e. **0.32× the +0.0149 residual**
  → `weighting_fragile = false` (threshold 0.5×residual = 0.00745). The pocoMC
  importance weights are near-uniform on the ocean/no-ocean split (ESS/N ≈
  0.67–0.73); the reference fraction is NOT a weighting artifact.
- **On treatment B (pre-repair "emulation") — reviewer correction (2026-08-11):**
  the `0acff866` bug did NOT produce an ocean fraction. It left `res.weights` as
  the last seed's array (length 6230) against the re-pooled samples (length
  20604), which *raised a length-mismatch assertion* — no statistic ever came
  out. The audit's treatment B (tiling seed 272's normalized weights across all
  three seeds' draws) is an arbitrary construct with no physical meaning; it
  lands near C only because the weights are near-uniform, so tiling washes out.
  **Treatment B is ill-defined and NOT load-bearing.** The
  weighting-robustness conclusion rests on **A vs C alone**; dropping B in fact
  *shrinks* the treatment swing (0.0056 → 0.0048), strengthening
  `weighting_fragile = false`. The earlier "repair moved the reference the
  correct direction / pre-repair ≈+0.0176" inference is withdrawn.
- **Caveat that motivates R3:** the three per-seed weighted fractions span
  **0.0121** (std 0.0063, ddof1) — comparable in size to the +0.0149 residual
  being adjudicated. At n_eff=2000 the reference comparator is under-resolved at
  exactly the level of the effect. This is the reference-side seed scatter R3
  (n_eff≥8000, ≥3 seeds) is designed to shrink.

## SNIS / corrected side (one flow conditioning, seed 72, fiducial x)

**Scope (reviewer note, 2026-08-11):** this single-seed audit bounds the
weighting sensitivity of the **fiducial** corrected residual (+0.0120 vs
reference A), NOT the **pooled** 3-seed residual (+0.0149) that the C16
escalation named. The pooled/multi-seed corrected behavior is R2's remit; the
"moves materially" reasoning below is against the fiducial residual it actually
probes.

Determinism check: re-conditioning reproduced the committed
`log_weights.npy` **byte-for-byte** (`max_abs_logw_diff = 0.0`, 226 rejected
draws match). The reconstructed per-draw ocean mask is therefore validated —
the `raw` treatment reproduces the committed corrected fiducial 0.9292539
exactly. Pareto-k = 0.191 (well below the 0.5 clean gate; PSIS would not fire).

| treatment | corrected frac | residual vs ref A |
|---|---|---|
| raw (committed) | 0.92925 | +0.01200 |
| PSIS-forced (tail smoothing forced on despite k<0.5) | 0.92935 | +0.01210 |
| top-0.1% ablated (20 of ~19774 finite weights zeroed) | 0.92449 | +0.00724 |

- raw → PSIS-forced swing = **+0.0001** (near-no-op; the tail is not heavy).
- top-0.1% ablation moves the corrected fraction by **−0.0048** (≈0.33× the
  residual). Removing the 20 heaviest draws shrinks the residual to +0.0072 but
  does not eliminate it. This is a *deliberately aggressive* stress (zeroing the
  single heaviest bin of ocean-favoring mass); it does not, on its own, meet a
  "moves materially" bar that would reopen the corrected side, but the
  reviewer should rule on that threshold.

## Preregistered reading (for the reviewer — NOT self-adjudicated)

Manager's R1 spec: *"if the +0.015 residual moves materially under any of these,
the corrected side reopens."* Observed:

1. **Reference side is weighting-robust** (swing 0.38× residual; correct vs
   unweighted differ 0.0048; the repair helped, not hurt).
2. **Corrected side is weighting-robust to raw/PSIS** (0.0001 swing) and only
   partially sensitive to the most aggressive tail ablation (0.33× residual,
   residual persists).
3. **Neither side's weighting treatment explains the full +0.0149.** The largest
   single sensitivity is the reference-side *seed scatter* (0.0121), pointing to
   reference-side sampling resolution — exactly R3's target — rather than a
   correction-machinery defect.

**Recommendation to the reviewer (advisory only):** R1 does not reopen the
corrected/SNIS side and does not re-ratify C16. It supports proceeding to R2
(one N=100k corrected run — tighter corrected SE) and R3 (decisive reference
recompute) per §0.20. Final routing is the reviewer's.

## Reviewer routing (2026-08-11) — ADOPTED

- **Verdict: PASS-WITH-CONCERNS. Proceed to R2.** Nothing needs fixing before
  R2; the two documentation corrections above are recorded, not blocking.
- **Corroboration cited by the reviewer:** existing C13 evidence (corrected
  0.9293/0.9310/0.9363, stable and not shrinking; residual ≈19× the 1/ESS bias
  ceiling) confirms the tension is NOT finite-N corrected-side bias. The largest
  single sensitivity in R1 is the reference-side **seed scatter (0.0121, std
  0.0063)** at n_eff=2000 — exactly R3's target.
- **Carried to R3:** report the **pooled reference SE** (≈0.0063/√3 ≈ 0.0036) so
  the |corrected − reference| ≤ 2× combined-SE re-ratification gate is evaluated
  against the same treatment-A pooling used here.
