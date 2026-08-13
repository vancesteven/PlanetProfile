# Resume note — Salt IS validations (Machine B)

_Written 2026-08-12 before a planned shutdown._

## What is done (durable on disk)

- **Pooled reference repair COMPLETE for MgSO4 + NaCl.** The 2026-08-10 repool
  had left the 12 derived per-sample arrays (log_likelihoods, k2_results,
  D_ocean_results, c20/c22, h2, cmr2, D_* thicknesses, cmr2_hydro) at the last
  seed's length; only samples+weights were pooled. `titanG_repool_reference.py`
  now pools ALL 13 per-sample arrays in seed order (verified aligned: MgSO4
  19325, NaCl 18514, wsum=1.0). See `../reference/repool_repair_note.json`
  (`repool_completion_2026_08_12`). No MCMC recompute; no science-assumption
  change. §0.16 crosscheck verdicts are NOT affected (that path reads only
  samples/weights).

## MgSO4 result — COMPLETE (2026-08-12, is_validation_mgso4.json)

Ran clean against the repaired reference. EVIDENCE, not adjudicated.
- C3 byte-identity: max_rel_diff 0.0, PASS (confirms the reference repair).
- IS diagnostics: clean, Pareto-k 0.326, ESS 1592, ESS/N 0.080, no PSIS.
- Reverse coverage (C5.3): PASS. Pushforward: PASS (gap_im 0.33s, gap_re 0.21s;
  WKS_im 0.135, WKS_re 0.075).
- **C16 ocean fraction: reference 0.4901 · corrected 0.5913 · residual +0.1012.**
  SAME-SIGN as NH3 (+0.015), larger. Per the preregistered rule this corroborates
  the reference-side sampling-resolution systematic (§0.20 R3), NOT an
  NH3-specific proposal cause. ESS/N 0.080 > NH3's 0.06. -> reviewer, not self.

## NaCl result — COMPLETE (2026-08-12, is_validation_nacl.json)

Ran clean against the repaired reference (18514 aligned). EVIDENCE, not adjudicated.
- C3 byte-identity: max_rel_diff 0.0, PASS.
- IS diagnostics: clean, Pareto-k 0.261, ESS 1516, ESS/N 0.076, no PSIS.
- Reverse coverage (C5.3): PASS. Pushforward: PASS (gap_im 0.34s, gap_re 0.17s;
  WKS_im 0.130, WKS_re 0.067).
- **C16 ocean fraction: reference 0.5318 · corrected 0.6498 · residual +0.1180.**
  SAME-SIGN as NH3 (+0.015) and MgSO4 (+0.1012). -> reviewer, not self.

## How to read the result (manager's preregistered rule — EVIDENCE, not a gate)

Compare each salt's C16 ocean-fraction residual SIGN against NH3's +0.015:
- SAME-SIGN (salts also corrected > reference) -> corroborates the R3
  reference-side sampling-resolution systematic (§0.20 resolution).
- near-zero / mixed-sign -> points at an NH3-specific proposal cause
  (NH3 ESS/N 0.06 vs salts' higher).
Record as EVIDENCE. NaCl crosscheck is expected to FAIL on log10_eta_V only
(poorly-identified nuisance, non-blocking, seen in §0.16).

Route scientific adjudication of the salt evidence to the scientific-reviewer;
never self-adjudicate. Do NOT commit/push without explicit user instruction
(Machine A owns pushes).
