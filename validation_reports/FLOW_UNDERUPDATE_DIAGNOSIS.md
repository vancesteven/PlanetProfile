# Flow under-update diagnosis — NH3 tidal-sector k2 miss (2026-08-04, in progress)

**Question.** The Titan NH3 joint SBI artifact's posterior-predictive under-updates
the tidal sector: SBI-pp |Im_k2| median 0.042 vs MCMC-pp reference 0.100 (obs 0.135;
JSON `mcmc_pp_median`), and the PPC flags both k2 channels (Re 0.94σ, Im 1.64σ).
SBC-PASS + per-param crosscheck coexist with this datum-local miss. Why does the NH3
flow under-update k2 when every deployed/candidate Europa flow assimilates its
informative k2 cleanly?

**Manager guidance (countersign 1d644be7).** Lead with hypothesis #1
(noise-augmentation swamping) and #3 (x-norm scale); #2 (tail sparsity) is weakened
as a standalone mechanism by v4's clean assimilation of a MORE extreme-tail Re_k2
(0.3 pctile) than NH3's (86th). Look for what DIFFERS between the NH3 (4-channel)
and Europa (20–21-channel) training setups — obs-vector width and the NH3 joint
no-ocean+ocean mixture are the two candidates.

Ranked hypotheses (MACHINE-B-HANDOFF §0.6): #1 noise-augmentation swamping,
#2 tail sparsity, #3 x-normalization scale, #4 capacity/embedding.

---

## Hypothesis #3 — x-normalization scale: RULED OUT at the channel level

**Mechanism.** sbi's `posterior_nn` uses `z_score='independent'`: each x channel is
standardized by its TRAINING-SET std (the prior-predictive spread, ≈ √(signal²+σ_obs²)),
stored in `artifact['x_norm']['std']`. A σ_obs-scale shift in x becomes
(σ_obs / train_std) of a z-unit. If train_std ≫ σ_obs the conditional must resolve
very fine z-differences to update at the datum — a candidate under-update mechanism.

**Test.** Pure artifact inspection (no retraining): per-channel `train_std / σ_obs`
and the backed-out signal-to-noise `signal_std/σ_obs = √(train_std²−σ_obs²)/σ_obs`,
across NH3 (under-updates), the Titan no-ocean control, and the clean Europa flows.

| channel | train_std | σ_obs | train_std/σ_obs | signal/noise |
|---|---|---|---|---|
| **NH3 Re_k2** | 0.1626 | 0.048 | **3.39** | 3.24 |
| **NH3 Im_k2** | 0.1176 | 0.035 | **3.36** | 3.21 |
| Titan-noocean Re_k2 | 0.2270 | 0.048 | 4.73 | 4.62 |
| Titan-noocean Im_k2 | 0.0759 | 0.035 | 2.17 | 1.93 |
| v4 Re_k2 (clean) | 0.0761 | 0.015 | 5.07 | 4.97 |
| v4 Im_k2 (clean) | 0.0229 | 0.015 | 1.53 | 1.15 |
| v6 Re_k2 (clean) | 0.0767 | 0.015 | 5.11 | 5.01 |
| v1.1 Re_k2 (clean) | 0.0901 | 0.050 | 1.80 | 1.50 |
| v1.1 Im_k2 (clean) | 0.0529 | 0.050 | 1.06 | 0.34 |

All-channel `train_std/σ_obs`: NH3 min 3.4 / median 14.0 / max 46.7;
v4 min 1.0 / median 1.6 / max 80.1; v6 min 1.0 / median 1.4 / max 12.7.

**Reading — the scale story does NOT hold:**
- **NH3's k2 channels are the LEAST compressed in the NH3 vector** (ratio 3.4, vs the
  all-channel median 14). If z-scoring were starving k2 of dynamic range, k2 would
  have the HIGHEST ratio, not the lowest. The opposite of the #3 prediction.
- **NH3's k2 signal-to-noise (3.2) is as good as or better than the clean Europa k2
  channels** it should be compared against (v4/v6 Im_k2 1.15, v1.1 Im_k2 0.34,
  v1.1 Re_k2 1.50). A channel that Europa assimilates fine at S/N ~1 cannot be
  under-updated in NH3 at S/N ~3 by a scale/noise argument.
- **The Titan no-ocean control has the SAME σ_obs and comparable ratios** (Re 4.7,
  Im 2.2) — so if scale/noise drove the miss, no-ocean would miss too. (Its PPC is
  the decisive control — see below.)

**Verdict #3: RULED OUT as the mechanism.** The k2 channels are not scale-compressed
relative to channels the same architecture assimilates cleanly elsewhere.

**Corollary for #1 (noise-augmentation swamping).** #1 turns on the same
signal-to-noise quantity: the flow rationally down-weights a channel whose signal
spread is small relative to the noise it was trained against. NH3 k2 S/N (3.2) is
NOT small and is comparable to clean cases — so the per-channel form of #1 is also
weakened by this inspection. #1's pilot-retrain test (reduced/zero k2 noise) is still
worth running IF the mixture control below comes back clean-but-inconclusive, but the
scale inspection already argues the mechanism is not a simple per-channel S/N deficit.

---

## The discriminating control — Titan no-ocean vs NH3 (rules out width; isolates the ocean-admitting apparatus)

The manager's two candidates for "what differs" are **obs-vector width** (4ch vs
20–21ch) and the **joint no-ocean+ocean mixture**. These are confounded in the
Europa-vs-NH3 contrast (different body/composition/sampler). But the repo already
holds a near-ideal control for the WIDTH axis: the **Titan free-gravity NO-OCEAN**
artifact (`titan_freegrav_noocean_posterior_1m.pt`) with its own matched reference
MCMC (`titan_freegrav_noocean_1m/reference/`), conditioned on the SAME 4 channels
{C20,C22,Re_k2,Im_k2}, SAME body, SAME σ_obs, SAME real Petricca k2 datum, SAME
free-gravity design.

**The control is NOT "NH3 minus the mixture only."** The no-ocean and NH3 setups
differ by more than the bimodality (confirmed from the two configs + artifacts):
- **Parameter count:** no-ocean = 12 params; NH3 = 13 params. NH3 carries the extra
  `log10_wOcean_ppt` salinity axis; no-ocean does not.
- **Phase guard:** no-ocean sets `phase_stability = {'enforce':'no_ocean_Ih',
  'margin_K':0.1}`; NH3 sets `phase_stability = None`.

The salinity axis and the ocean/no-ocean bimodality co-vary (no liquid layer → no
salinity dimension and no ocean branch), so they are NOT physically separable in
this pair. The control therefore isolates the **entire ocean-admitting apparatus**
— the joint no-ocean+ocean mixture **plus** the added (near-degenerate, strong Tb↔w)
salinity dimension **plus** the removed phase guard — as a bundle, not the mixture
alone.

- If no-ocean **assimilates k2 cleanly (0 flagged)** → obs-vector width (4ch) is
  cleared as the *causal* story, and the mechanism lives inside the ocean-admitting
  apparatus (mixture and/or salinity axis). This kills the width reading.
- If no-ocean **also under-updates k2** → the 4-channel width (or something shared
  by both Titan 4ch flows) is implicated, independent of the apparatus.

### Result (seed 58): no-ocean assimilates k2 CLEANLY — 0/4 flagged

| channel | obs | prior-pp | MCMC-pp | SBI-pp | MCMC update | \|SBI−MCMC\| |
|---|---|---|---|---|---|---|
| C20 | −3.351e-5 | −3.236e-5 | −3.351e-5 | −3.350e-5 | 2.44σ | 0.04σ |
| C22 | +1.011e-5 | +9.70e-6 | +1.011e-5 | +1.011e-5 | 6.65σ | 0.03σ |
| **Re_k2** | 0.608 | 0.375 | 0.563 | 0.561 | **3.92σ** | **0.05σ** |
| Im_k2 | 0.135 | 0.093 | 0.096 | 0.095 | 0.10σ | 0.03σ |

Report: `/tmp/ppc_titan_noocean/ppc_pushforward_report.json` (copied to
`validation_reports/titan_freegrav_noocean_1m/ppc/`).

**The decisive channel is Re_k2 — and ONLY Re_k2.** In the no-ocean control, Re_k2
is strongly informative: MCMC-pp updates 3.92σ off the prior-predictive and the flow
tracks it to 0.05σ. That is a clean positive assimilation test on an informative
channel that a 4-channel vector must resolve. But **Im_k2 is non-informative in the
control**: obs 0.135 lies above the no-ocean reachable range, so the MCMC-pp barely
moves (prior-pp 0.0926 → mcmc-pp 0.0960 = 0.10σ update) and there is nothing to
assimilate. The control positively tests clean k2 assimilation for **Re_k2 only**.

**Verdict: obs-vector width RULED OUT; the ocean-admitting apparatus isolated.**
- A 4-channel vector is NOT the cause: the no-ocean 4-channel flow assimilates its
  informative Re_k2 cleanly. This kills the obs-vector-width reading (which the
  v5/v6/v7 batch already flagged as "supporting, not controlled") as a *causal*
  story. **Positively established.**
- The mechanism lives inside the **ocean-admitting apparatus** — the joint mixture
  and/or the added near-degenerate salinity dimension (they co-vary and cannot be
  separated in this control pair). This is narrower than "the mixture ONLY," but
  still confirms the manager's alternate candidate direction (countersign 1d644be7:
  "the joint no-ocean mixture may bimodalize the conditional").

**Candidate mechanism (bimodal-conditional mode-assignment) — NOT yet positively
established.** The plausible story: the joint posterior is bimodal in structure
space — a frozen no-ocean branch (low tidal response; the reachability check gave
frozen Re_k2≈0.113) and an ocean branch (higher k2). At the NH3 datum the reference
MCMC concentrates on the high-k2 ocean branch; the amortized flow's conditional
retains mass on BOTH modes, dragging its pushforward median toward the low-k2 frozen
mode → the under-update. **What the control does NOT show:** it rules out ONE
alternative (width); it does not positively evidence bimodality. Two caveats weaken
the two-mode reading as currently evidenced:
- NH3 SBI Re_k2 pushforward median is 0.541 — close to the ocean value, far from the
  frozen 0.113. A retained-frozen-mode drag would sit lower; 0.541 reads at least as
  easily as a generic failure-to-concentrate.
- NH3 SBI Im_k2 median 0.042 sits near the *prior*-predictive median (0.026), which
  again reads as failure-to-concentrate rather than a specific two-mode median drag.
- The added near-degenerate salinity dimension (documented strong Tb↔w degeneracy,
  corr −0.986) is a live alternative contributor to the conditional smearing and is
  NOT excluded by this control. To separate salinity from bimodality would require an
  ocean-only artifact that still carries `log10_wOcean_ppt` but no frozen branch (or
  the joint mixture with salinity fixed) — not run.

Positively establishing bimodal mode-assignment would require showing the NH3 SBI k2
pushforward is actually bimodal with retained mass on a frozen mode ≈0.11. Recorded
as a **candidate mechanism**, not established.

**Bearing on #1 (noise-augmentation swamping) — cleared for Re_k2 only.** The
no-ocean control shares NH3's σ_obs, noise convention, and datum and assimilates
**Re_k2** cleanly, so noise convention is not the discriminator *for Re_k2*. But the
control is non-informative for **Im_k2** (0.10σ), and NH3 flags Im_k2 HARDER (1.64σ)
than Re_k2 (0.94σ). Im_k2 is exactly the channel where the abs-fold + additive-noise
convention is most relevant (signal near zero, noised value can cross zero). So the
"#1 not on the critical path" conclusion holds for Re_k2 but is **untested for the
dominant Im_k2 miss** — the #1 pilot (reduced/zero k2-noise retrain) remains a live
test for Im_k2 specifically, not off the table.

---

## Consequences

1. **MgSO4/NaCl inheritance is CONDITIONAL on the mechanism.** IF the driver is the
   joint mixture (or the co-varying salinity axis), then MgSO4/NaCl — planned as the
   same joint no-ocean+ocean build (`build_tbw_grid_cache(..., retry_frozen_as_no_ocean
   =True)`, a `cache_builder.py` argument, NOT a sampler-config field — the config's
   `retry_frozen_as_no_ocean` is null in both NH3 and no-ocean) — inherit the same
   deficiency by construction. But if the true driver is instead an Im_k2-specific
   reachability/noise-fold gap between the frozen and ocean branches, a different salt
   with a different phase diagram and different reachable k2 range need NOT inherit
   identically. Either way the conservative posture is the same: the split-status
   discipline adopted for NH3 (gravity/structure VERIFIED, tidal sector MCMC-
   authoritative) and the first-class pushforward-observable gate (must flag the
   known NH3 k2 miss) apply. **Verify the joint-build flag from the actual cache
   metadata for each MgSO4/NaCl build** rather than assuming it — do not assert
   inheritance as a certainty.

2. **The likely remedy is apparatus-aware, not "more epochs" / "less noise."** If the
   candidate mechanism holds, a fix would target the bimodal conditional / salinity
   degeneracy: e.g. a mixture-aware embedding, a has_ocean-labelled conditional, or a
   sequential/truncated round near the datum. But because the mechanism is not yet
   positively established (and #1 is still open for Im_k2), a reduced-noise pilot is
   NOT excluded as a cheaper first probe of the Im_k2 miss. None of these is authorized
   compute — surface to the manager; do NOT fold into a production build without a
   frozen config.

3. **MgSO4/NaCl compute remains HELD** — the diagnosis explains the NH3 miss at the
   width axis (ruled out) and localizes the mechanism to the ocean-admitting apparatus,
   but has NOT positively established which sub-cause (bimodal mixture vs salinity
   degeneracy vs Im_k2 noise-fold) drives it. Training a third/fourth artifact that
   likely inherits the SAME deficiency, absent a chosen remedy, is exactly the waste
   the manager's sequencing guards against. Decision to proceed (split-status as-is)
   vs. run the two cheap separators first (ocean-only-with-salinity PPC; #1 reduced-
   noise pilot for Im_k2) vs. remediate is a manager/user call.

## Status
- #3 (x-norm scale): **RULED OUT** at the channel level (artifact inspection —
  reviewer-reproduced, all ratios exact).
- Obs-vector width (4ch): **RULED OUT as causal** — no-ocean 4ch flow assimilates
  informative Re_k2 (3.92σ MCMC update) to 0.05σ. Positively established.
- Mechanism: localized to the **ocean-admitting apparatus** (joint mixture and/or the
  co-varying salinity axis + removed phase guard) — NOT the mixture alone (the control
  bundles a 13th salinity param + phase-guard removal). Bimodal mode-assignment is a
  **candidate**, not positively established (would need the NH3 SBI k2 pushforward
  shown bimodal with mass on the frozen mode ≈0.11).
- #1 (noise-augmentation swamping): **cleared for Re_k2** by the control; **untested
  for the dominant Im_k2 miss** (control non-informative for Im_k2, 0.10σ). Pilot
  remains a live cheap probe for Im_k2.
- #2 (tail sparsity), #4 (capacity): not reached.

## Reviewer disposition (2026-08-04)
Scientific-reviewer PASS WITH CONCERNS. Upheld: #3 falsification (ratios reproduced
exactly; k2 lowest ratio, gravity highest — prediction inverted), the variance-
subtraction decomposition under the abs-fold (noise added post-fold, independent of
θ, not re-folded → Var(x)=Var(|Im|)+σ² exact), and "width ruled out." Required
corrections folded above: (1) control isolates the ocean-admitting apparatus, not the
mixture alone (12 vs 13 params; phase_stability enforce vs None); (2) bimodal
mode-assignment downgraded to candidate; (3) clean-assimilation restricted to Re_k2,
Im_k2 miss untested by the control; (4) MgSO4/NaCl inheritance made conditional +
verify the joint-build flag from cache metadata. Abstract Im_k2 reference reconciled
to the JSON `mcmc_pp_median` 0.100.
