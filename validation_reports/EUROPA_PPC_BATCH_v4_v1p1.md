# Europa deployed-artifact PPC batch — v4 + v1.1 (2026-08-04)

**Purpose.** Machine A manager advisory (MACHINE-B-HANDOFF §0.5/§0.6): the Titan
NH3 campaign proved SBC-PASS + per-param crosscheck can coexist with a
datum-local pushforward under-update the per-param gate is structurally blind
to. Run the same posterior-predictive check (PPC) on the DEPLOYED Europa
artifacts (v4 + v1.1 first — they are on the public HF app) to decide whether
the public slots need a sector warning before the held HF redeploy.

**Method.** `plans/scripts/ppc_pushforward_check.py` (body/artifact-agnostic
generalization of `titanG_ppc_interior_check.py`). Per observable channel it
reports a four-way median: prior-predictive / SBI posterior-predictive / MCMC
posterior-predictive / observed. All pushforwards reuse the exact
`MCMCRunner.generate_sbi_dataset(theta_override=...)` theta->x loop (validated
exact, max|fwd-stored|=0); NOISELESS (obs_noise=False). MCMC-pp median is
importance-WEIGHTED by the reference weights.

**Preregistered flag statistic (fixed before running, never moved):**
a channel is FLAGGED when `|median_pp(SBI) - median_pp(MCMC)| > 0.5 * sigma_obs`.

**Sanity validation of the statistic (the LICENSE for the whole exercise).** The
IDENTICAL generalized script (`ppc_pushforward_check.py`, same `|SBI−MCMC|`
statistic, n_post=4000, seed 72) was run on the NH3 artifact + its reference MCMC
and the report is PERSISTED at
`validation_reports/titan_freegrav_nh3_1m/ppc/ppc_pushforward_report.json`
(2/4 flagged). It reproduces the known pattern: **FLAGS both k2 channels (Re
0.94σ, Im 1.64σ) and NEITHER gravity channel (C20 0.04σ, C22 0.07σ).** The
statistic demonstrably fires on the known-bad tidal case and stays silent on the
well-assimilated gravity case — so "0 flagged" on an *informative* channel is a
meaningful null. (Earlier drafts cited an n_post=2000 sanity run that was not
persisted; this archived n_post=4000 run supersedes it and matches the Europa-run
statistic exactly.)

---

## Result: neither deployed artifact exhibits the NH3-style under-update

| Artifact | channels | flagged | max |SBI−MCMC| | verdict |
|---|---|---|---|---|
| Clipper–Europa v4 geodesy 11D | 21 | **0** | 0.09σ (Bind_synodic_x_imag) | no sector warning needed |
| Galileo–Europa v1.1 8D | 5 | **0** | 0.08σ (CMR2) | no sector warning needed |

Reports:
- `validation_reports/europa_clipper_v4_1m/ppc/ppc_pushforward_report.json`
- `validation_reports/europa_galileo_v1p1_1m/ppc/ppc_pushforward_report.json`

**IMPORTANT — read "0 flagged" correctly (scientific-reviewer, 2026-08-04).**
The flag `|median_pp(SBI) − median_pp(MCMC)| > 0.5·σ` can only fire on a channel
that the MCMC reference *itself* materially updates away from the prior-predictive
(i.e. an *informative* channel) AND on which the SBI flow then fails to follow.
On a non-informative channel — one where the datum carries little information and
the posterior barely moves off the prior — the flag is **structurally
insensitive**: "0 flagged" there means "no signal to test," NOT "flow verified."
So the batch's clean result is carried by the channels that ARE informative, not
by the count of zeros. See the per-channel breakdown below.

### v4 (21 channels) — the decisive channel is Re_k2
The one informative tidal channel in the deployed v4 set behaves correctly:

| channel | σ | prior-pp | MCMC-pp | SBI-pp | MCMC update | \|SBI−MCMC\| |
|---|---|---|---|---|---|---|
| **Re_k2** | 0.015 | 0.2741 | 0.2527 | 0.2534 | **1.43σ** (datum at 0.3 pctile) | **0.04σ** |
| Im_k2 | 0.015 | 0.0022 | 0.0029 | 0.0030 | 0.04σ | 0.01σ |
| Re_h2 | 0.10 | 1.224 | 1.136 | 1.137 | ~0.13σ | 0.01σ |
| Im_h2 | 0.10 | 0.0065 | 0.0110 | 0.0113 | ~0.03σ | 0.00σ |

**Re_k2 is genuine positive evidence.** The datum sits at the 0.3 prior-predictive
percentile — an *extreme* informative tail, more extreme than NH3's 86th — the
MCMC reference updates 1.43σ off the prior, and the SBI flow tracks that update to
0.04σ. On the one deployed tidal channel where the NH3 pathology *could* show up,
it demonstrably does not, on a harder datum than NH3's. Im_k2 / Re_h2 / Im_h2 are
barely updated by the MCMC (≤0.13σ), so their 0-flagged status is expected, not
informative. C20/C22 <0.1σ; all 14 induction channels ≤0.09σ.

### v1.1 (5 channels) — all tidal channels non-informative by construction
CMR2 0.08σ, Re_k2 0.01σ, Im_k2 0.00σ, Re_h2 0.01σ, Im_h2 0.00σ. The MCMC reference
updates every v1.1 k2/h2 channel by ≤0.13σ off the prior-predictive (e.g. Re_k2
prior 0.2742 → MCMC 0.2695), because v1.1's k2/h2 are **explicitly-labeled
hypothetical-conditioning channels** (σ 0.05/0.10; only CMR2 + the synodic |Ae|>0.7
support cut are honest Galileo data). A flow frozen at the prior would still score
≤0.09σ < 0.5σ on these and would NOT flag — so 0-flagged on v1.1's tidal sector is
the *expected* result of a non-informative datum, and carries no flow-fidelity
assurance. This is fine for the deployment question: **the NH3 pathology is by
definition a failure to assimilate informative data, and cannot arise where the
datum is non-informative at the deployed default.**

---

## Why NH3 differs — empirical, not a "tail" rule

An earlier draft of this report attributed NH3's k2 miss to "a tight measurement
in the informative upper tail (86th percentile)." **That causal story is retracted:
it is contradicted by v4's own data** — v4 Re_k2 sits at the 0.3 percentile (a
*more* extreme tail) yet the flow assimilates it correctly (0.04σ). Tail-location
alone does not produce the under-update.

The honest, empirical statement is narrower:
- **v4** assimilates its one informative tidal datum (Re_k2, extreme-tail,
  MCMC-updated 1.43σ) correctly — SBI tracks MCMC to 0.04σ.
- **v1.1**'s tidal channels are non-informative by construction (hypothetical
  conditioning), so the NH3 pathology cannot arise at the deployed default.

The distinction from NH3 is most likely training-pair density / flow behavior at
the specific NH3 joint datum (a no-ocean+ocean *mixture* posterior), NOT a general
"tail → under-update" law. The under-update diagnosis (queued, pilot-scale) is
where that mechanism gets pinned down; this report does not lean on it.

### Model-adequacy caveat: out-of-envelope tidal defaults (separate from the flag)
Four channels place the deployed *default* datum outside the 5–95 prior-predictive
band. This is a **model-adequacy** note about the deployed defaults, distinct from
the flow-agreement (tidal-warning) question:
- **Im_h2 = 0 exactly** (v4 & v1.1, pctile 0.0): the forward model never produces
  Im_h2 ≤ 0 — tidal dissipation is positive, so Im_h2 > 0 always. Conditioning on
  Im_h2 = 0 conditions on a physically-unreachable boundary; with σ = 0.1 the
  Gaussian likelihood puts ~half its mass on unphysical negative Im_h2 and the
  channel is effectively one-sided/uninformative. Both samplers agree at the
  model's achievable floor (~0.005–0.011) — expected.
- **Re_k2 = 0.23** (v4 pctile 0.3 / v1.1 0.8): below the lower prior-predictive
  tail; the posterior is pulled to the model's achievable edge (0.253 / 0.269) and
  cannot fully reach 0.23 — a mild model-adequacy tension.

On all four channels the SBI–MCMC agreement is ≤0.013σ: the flow tracks the
reference even at the envelope edge, so this triggers **no** NH3-style tidal-sector
warning. It is recorded here so the deployed-default tidal values are known to lie
partly beyond what the interior model can produce.

---

## Verdict for the held HF redeploy

**The public v4 and v1.1 slots do NOT need a tidal-sector warning.** The honest
basis (scientific-reviewer PASS WITH CONCERNS, 2026-08-04) is:
1. **v4 assimilates its one informative tidal datum correctly** — Re_k2 at the 0.3
   percentile, MCMC-updated 1.43σ, SBI tracked to 0.04σ; the NH3 under-update does
   not appear on the deployed tidal channel where it could.
2. **v1.1's tidal channels are explicitly-labeled non-informative hypotheticals**,
   where the NH3 pathology cannot arise at the deployed default.

This is NOT "the flow is verified clean across the whole tidal sector" — on the
non-informative channels the flag is insensitive and 0-flagged is expected. The
deployment decision is defensible on (1)+(2). The v4/v1.1 PPC is the gate the
manager set on the redeploy; both pass on that basis.

**Untested regime (GUI caveat, reviewer optional item).** The PPC tests only the
deployed *default* conditioning. It gives no assurance about flow behavior if a GUI
user dials an *off-default, informative* k2 (a tight σ at an extreme theory value)
— the exact regime where NH3 failed. A follow-up PPC at one or two off-default
informative k2 conditionings is recommended before relying on the flow for
user-supplied extreme tidal hypotheticals; it does not gate the current redeploy
(deployed defaults only).

## Provenance
- v4 reference MCMC: committed `Test53_geodesy_v4/europa_clipper_v4_reference_result.pkl`
  (4430 samples, ESS 4159/94%).
- v1.1 reference MCMC: was `/tmp`-only and lost; REGENERATED this session from
  `europa_galileo_v1p1_8D.json` seed 50 via `titanG_reference_mcmc.py`
  (4267 samples, ESS ~98%, 1.5 min) →
  `validation_reports/europa_galileo_v1p1_1m/reference/`. Reproducibility gap
  closed by COMMITTING THE REGENERATION RECIPE — `titanG_reference_manifest.json`
  (seed 50, config path, sample count) + `titanG_reference_progress.jsonl` — which
  regenerate the pkl deterministically. The pkl itself stays gitignored per repo
  convention (`*.pkl`; same as the NH3 reference), but is no longer only-in-`/tmp`.
  (An earlier draft cited "log_Z −3.02"; the pickle exposes `log_evidence = None`,
  so that number is unconfirmed and immaterial to the median-based PPC — dropped.)
- Training `.npz` for both were `/tmp`-only and lost; prior-predictive column
  regenerated from the config prior through the identical forward loop (noted in
  each report's `prior_predictive_source`).
- NH3 license (statistic sanity): PERSISTED
  `validation_reports/titan_freegrav_nh3_1m/ppc/ppc_pushforward_report.json`
  (generalized script, n_post=4000 seed 72; flags Re_k2 0.94σ + Im_k2 1.64σ,
  neither gravity) — closes the reviewer's MODERATE-3 (statistic must have
  demonstrably fired on the known-bad case, from an archived artifact).
- Both PPC reports interpreted by the scientific-reviewer (PASS WITH CONCERNS,
  2026-08-04) before STATUS update. Required corrections MAJOR-1/2 (report carried
  by v4 Re_k2, not "wide σ"; "tail → under-update" causal claim retracted),
  MODERATE-3 (NH3 license persisted), and MODERATE-4 (out-of-envelope model-adequacy
  caveat) are all folded into this report. Optional item (off-default informative-k2
  PPC) noted as a GUI follow-up; does not gate the deployed-default redeploy.
