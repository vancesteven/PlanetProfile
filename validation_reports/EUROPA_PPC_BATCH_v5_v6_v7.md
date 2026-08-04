# Europa v5/v6/v7 baseline PPC batch (2026-08-04)

**Purpose.** Second half of the Machine A manager directive (MACHINE-B-HANDOFF
§0.5 B5 EXTENDED): after the deployed v4 + v1.1 artifacts came back clean, run
the identical posterior-predictive check (PPC) on the remaining Europa campaigns
— v5 (thick-ice geodesy), v6 (free-gravity), v7 (open-|Ae|). These are
**delivered but NOT ratified** (GUI-gated pending the v5/v6/v7 gate adjudication);
the PPC is a flow-fidelity diagnostic, NOT a ratification action — those verdicts
stay with Machine A + the adjudication doc.

**Method.** `plans/scripts/ppc_pushforward_check.py` (the same body/artifact-agnostic
runner used for v4/v1.1/NH3). Per observable channel it reports a four-way median:
prior-predictive / SBI posterior-predictive / MCMC posterior-predictive / observed.
All pushforwards reuse the exact `MCMCRunner.generate_sbi_dataset(theta_override=...)`
theta→x loop (validated exact, max|fwd−stored|=0); NOISELESS (obs_noise=False).
MCMC-pp median is importance-WEIGHTED by the reference weights.

**Preregistered flag statistic (fixed before running, never moved):**
a channel is FLAGGED when `|median_pp(SBI) − median_pp(MCMC)| > 0.5 * sigma_obs`
(FLAG_K = 0.5, the module constant). Same statistic validated on NH3 (flags both
k2 channels, neither gravity) — license persisted at
`validation_reports/titan_freegrav_nh3_1m/ppc/ppc_pushforward_report.json`.

**Scope — BASELINES ONLY (deliberate).** Each campaign has three arms
(baseline 11D + noinduction + nok2), but a matching reference MCMC exists ONLY for
the three baselines (`Test52_seawater_v5`, `Test53_seawater_v6`, `Test54_seawater_v7`).
The ablation arms condition on a DIFFERENT observable subset (noinduction drops the
14 Bind channels; nok2 drops Re/Im k2+h2), so there is no reference posterior
conditioned on the same data — pushing the baseline reference against an ablation
SBI would be an apples-to-oranges comparison, and the flag statistic
(|SBI−MCMC|) is undefined without a same-conditioning reference. The arms would
need their own reference MCMC before a PPC; that is not in this batch and is
NOT required for the flow-fidelity question the manager posed (the deployed and
candidate-default artifacts are the 11D baselines).

---

## Result: all three baselines clean — 0 channels flagged, and non-trivially so

| Artifact | channels | flagged | max \|SBI−MCMC\| | strongest informative channel (MCMC update) |
|---|---|---|---|---|
| v5 geodesy 11D | 21 | **0** | 0.12σ (C20) | C20 **22.6σ** → tracked 0.12σ |
| v6 freegrav 11D | 20 | **0** | 0.07σ (Bind_synodic_x_real) | Bind_synodic_x_imag 3.14σ → tracked 0.04σ |
| v7 open-\|Ae\| 11D | 21 | **0** | 0.27σ (Bind_synodic_x_imag) | Bind_synodic_x_real **52.9σ** → tracked 0.18σ |

Reports:
- `validation_reports/europa_clipper_v5_baseline_1m/ppc/ppc_pushforward_report.json`
- `validation_reports/europa_clipper_v6_baseline_1m/ppc/ppc_pushforward_report.json`
- `validation_reports/europa_clipper_v7_openae_1m/ppc/ppc_pushforward_report.json`

**These are STRONG clean results, not structurally-insensitive zeros.** Recall the
reviewer's v4/v1.1 caveat: the flag can only fire on a channel the MCMC reference
*materially updates* off the prior-predictive (an *informative* channel). On v4 the
clean verdict rested on a single informative tidal datum (Re_k2). Here the 11D
baselines carry MANY strongly-informative channels, and the flow tracks the MCMC
update on every one:

- **v5:** C20 updates **22.6σ**, C22 16.6σ, two Bind_synodic_x channels (3.06σ,
  3.31σ), CMR2 2.8σ, Re_k2 1.6σ, Re_h2 1.0σ — ALL tracked to ≤0.12σ.
- **v6:** the free-gravity design deliberately loosens C20/C22 (update only 0.4–0.5σ,
  as intended — gravity is agnostic), so Re_k2 (1.66σ) and the induction sector
  (Bind_synodic_x 2.8–3.1σ) carry the information; both tracked to ≤0.07σ.
- **v7:** the most stringent test in the whole PPC program — Bind_synodic_x_real
  updates **52.9σ** off the prior, Bind_synodic_y_imag 16.9σ, C20 23.0σ; the flow
  tracks each to ≤0.18σ, and the batch-wide maximum gap (0.27σ, Bind_synodic_x_imag
  on a 4.3σ update) is still comfortably under the 0.5σ flag. v7's observable-space
  gaps are the batch-largest (0.10–0.27σ, ~2–3× v5/v6); the most likely cause is
  Monte-Carlo noise from the diffuse open-|Ae| posterior, but this is unconfirmed
  and is carried forward as a flow-side input to B5 (do not smooth it away).

Across all three, the informative tidal channel Re_k2 (datum 0.23, MCMC-updated
1.35–1.66σ off the prior) is reproduced to ≤0.05σ — the same behavior seen on
deployed v4. On these 20–21-channel Europa vectors the flow faithfully reproduces
the reference-MCMC pushforward even on the most informative channels. ("Reproduces
the reference," not "gets the right answer" — the statistic establishes only
SBI-pp ≈ MCMC-pp; for v7 the reference is itself the pre-B3 object under
investigation, see below.)

---

## Relevance to the under-update diagnosis (obs-vector-width hypothesis)

The manager's countersign (1d644be7) refined the diagnosis: the v4 Re_k2 result
weakened tail-sparsity (#2) as a standalone mechanism, and pointed at
**obs-vector width** — with 21 channels the k2 pair may retain relative weight that
a 4-channel NH3 vector loses. **This batch adds three same-width observations
consistent with, but not controlling for, that read:**

- Every clean Europa PPC here is a **20- or 21-channel** vector, and every one
  reproduces its informative channels (including Re_k2) faithfully.
- The ONLY campaign that exhibits the k2 under-update (NH3) is a **4-channel**
  vector {C20, C22, Re_k2, Im_k2}.
- So the under-update DIFFERS between 4ch and 20–21ch: 20–21ch → no under-update;
  4ch → under-update. Note this batch adds three points at the SAME width level
  (20–21ch) — it does NOT vary width, so it is not a controlled width experiment;
  the only width contrast remains Europa vs NH3, which is confounded by body,
  composition, sampler, AND the NH3 joint no-ocean+ocean mixture. (v1.1, the one
  intermediate 5-channel point, is excluded because its tidal channels are
  non-informative by construction and cannot exhibit the pathology.) So this does
  not *prove* width is the mechanism, but it is consistent with hypotheses #1
  (noise-augmentation swamping — fewer channels means each noisy channel is a
  larger fraction of the conditioning signal) and #3 (x-norm scale), and it points
  at the genuinely controlled probe: **the noinduction/nok2 ablation ARMS** — a
  6-channel v6-noinduction and a 4-channel-analog would test width at fixed body +
  physics + sampler. Running those arms' PPC (needs each arm's own matched reference
  MCMC first) is the cheapest next diagnostic, ahead of any retraining, and must NOT
  be substituted by the same-width evidence above.

## Model-adequacy caveat (carried over from v4/v1.1, not a warning trigger)

`x_obs interior on all: False` on all three — same out-of-envelope defaults noted
for v4/v1.1: Re_k2 = 0.23 sits at the 0.0–0.6 prior-predictive percentile (below
the model's achievable tail; the posterior is pulled to the achievable edge ~0.251
and cannot fully reach 0.23) and Im_h2 = 0 is physically unreachable
(dissipation > 0). On both channels the SBI–MCMC agreement is ≤0.05σ — both samplers
agree at the model edge, so this is a model-adequacy note about the deployed
defaults, NOT a flow-fidelity flag.

## Verdict

**All three v5/v6/v7 baselines pass the PPC with no flagged channels**, and — unlike
the deployed v4/v1.1 where the clean verdict leaned on one or two informative
channels — these 11D baselines faithfully reproduce the reference-MCMC pushforward
across a large set of strongly-informative gravity, tidal, and induction channels
(up to a 52.9σ MCMC update tracked to 0.18σ). No NH3-style datum-local under-update
appears in any Europa baseline.

**This is a flow-fidelity diagnostic only.** It does NOT overturn or affect the
v5/v6/v7 gate adjudication (v5 NOT ratifiable on D_iceIh shape-excess; v6
conditionally ratifiable pending powered SBC; v7 blocked on reference-wander
diagnostics B3–B5). The PPC clean result is one more piece of evidence that the
FLOWS reproduce their references where they can be tested; the ratification blocks
are about SBC power, reference-MCMC wander, and a shape-excess in a science
parameter — separate questions this check does not address.

**Caveat on the v7 (and v5) references — orthogonal to the wander block.** The
|SBI−MCMC| statistic here uses the Test52 (v5) and Test54 (v7) reference MCMCs —
precisely the runs B3 will re-examine for wander (v5's log_Z_err 0.281 flags it as
the more suspect). This does NOT corroborate or refute the wander block: the
adjudication found the wander is a *parameter-space* shift (D_iceIh, ~1.06 km) that
carries ≈0 observable-space posterior-mass difference, so an observable-space PPC is
orthogonal to it. When B3 re-runs those references on fresh seeds, re-confirm these
observable-space PPC medians are stable (expected, since the wander is
observable-orthogonal) as a cheap corroboration that this check is insensitive to
the wander.

## Provenance
- v5 reference MCMC: `Test52_seawater_v5/europa_clipper_v5_reference_result.pkl`
  (4343 samples). v6: `Test53_seawater_v6/europa_clipper_v6_reference_result.pkl`
  (4349). v7: `Test54_seawater_v7/europa_clipper_v7_reference_result.pkl` (4352).
- Artifacts: `europa_clipper_{v5_geodesy_11D,v6_freegrav_11D,v7_openae_11D}_posterior_1m.pt`.
- Training `.npz` were `/tmp`-only; prior-predictive column regenerated from each
  config prior through the identical forward loop (n=10000; v5/v7 show the expected
  high support-rejection rate for their wide extended-Tb / open-|Ae| support — a
  warning, not a failure, consistent with the adjudication finding).
- Seeds: v5 55, v6 56, v7 57; n_post 4000 each.
- NH3 license (statistic sanity): `titan_freegrav_nh3_1m/ppc/ppc_pushforward_report.json`
  (flags both k2, neither gravity).
- Scientific-reviewer: **PASS WITH CONCERNS** (2026-08-04). Headline updates/gaps
  reproduce exactly from the JSONs; no flag-insensitivity conflation. Corrections
  folded in: MODERATE-1 ("direct controlled evidence" → three same-width
  observations, not a controlled width experiment); MODERATE-2 (v7 gap gradient
  0.10–0.27σ noted + carried to B5; v5/v7 references are pre-B3 and this PPC is
  orthogonal to the parameter-space wander); MINOR (v5 two Bind_synodic_x channels
  not three; "differs between 4ch/20–21ch" not "monotone"; v1.1 5ch exclusion
  stated; "reproduce the reference" not "assimilate correctly").
