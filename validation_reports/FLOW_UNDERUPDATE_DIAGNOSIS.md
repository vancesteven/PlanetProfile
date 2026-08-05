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

## Separators S1/S2 + capped anchor (2026-08-04, complete)
All four flows: nsf, seed 72, max_num_epochs=60 (all ran to cap, epochs_trained=61,
no early stop), n_post=4000, PPC noiseless. No new forward sims.

| flow | Re_k2 pp_med | Im_k2 pp_med | (obs Re 0.608 / Im 0.135; MCMC-pp Im 0.100) |
|---|---|---|---|
| **ANCHOR** capped full joint | 0.5446 | 0.0431 | reproduces deployed converged 0.042 |
| **S1** ocean-only (mixture off) | 0.5408 | 0.0392 | — |
| **S2** k2 σ/4 | 0.5381 | 0.0461 | — |
| **S2** k2 zero-noise | 0.5400 | 0.0469 | — |

**Findings (scientific-reviewer PASS WITH CONCERNS, 2026-08-04):**
- **Cap validity — CONFIRMED.** Capped anchor (0.0431) reproduces the externally
  converged deployed value (0.042), so 60 epochs suffice for the Im_k2 statistic and
  cap-vs-cap pilot readings are trustworthy. Depends on the deployed 0.042 being
  genuinely converged. `best_validation_log_prob` is null in all manifests (loss
  trajectory not independently inspectable — minor).
- **Mixture (S1) — CONFIRMED not the driver, direction strengthens it.** Removing the
  frozen no-ocean branch moved Im_k2 the WRONG way (0.0431→0.0392): a bimodal
  low-k2-drag mode would have RAISED the median toward obs. Robust to boundary-row
  misclassification (effect would have to reverse sign, not shrink).
- **Noise-swamping #1 (S2) — CONFIRMED closed for Im_k2.** Zero-noise (0.0469) and
  σ/4 (0.0461) both ≈ anchor despite fully removing k2 training noise. The zero-noise
  *prior*-pp median did drop (0.0369→0.0265), confirming the abs-fold+additive-noise
  convention inflates the prior Im — but the posterior-pp still cannot reach even
  0.10. Strongest possible falsification of #1 for Im_k2.
- **B3 reference-wander — CONFIRMED artifact.** Matched n_eff=2000 D_iceIh gap
  −0.19±0.22 km vs between-seed floor √(0.15²+0.09²)=0.18 km (|mean|/std=0.88, <1σ).
  Legacy 1.06 km shrinks ~5.6×; residual dominated by one seed (303 −0.44) inside the
  scatter. Neither v5 nor v7 is an outlier. ~0.18 km is a defensible 1σ floor for the
  §0.7 step-3 shape-excess re-eval **only at matched n_eff=2000**; reviewer advises a
  conservative ~2× (0.36 km) gate.

**Residual mechanism (by elimination, PLAUSIBLE not positively established):** the
amortized conditional fails to concentrate on the weakly-identified high-k2
ocean-branch tail (SBI-pp Im_k2 ≈ prior-pp ~0.04 vs MCMC-pp 0.100); the near-degenerate
salinity axis (Tb↔w −0.986) is the most probable specific contributor because it is
the one apparatus element S1 left in. S1/S2 do NOT rule out: (a) the salinity axis,
(b) the abs-fold *as a representation* (S2 changed only noise; x_clean is still the
folded |Im|), (c) the k2 support guard, (d) intrinsic sub-ceiling identifiability
(obs 0.135 > MCMC-pp 0.100, so part of the miss-vs-obs is model-datum tension).

**Reviewer note (important):** all four PPCs compare SBI-pp to obs (0.135), not to the
reachable MCMC-pp ceiling (0.100); the true SBI-vs-MCMC gap (~0.043 vs 0.100) was not
re-measured this pass.

**Reviewer-required follow-up separators before MgSO4/NaCl proceed:**
1. **Salinity-fixed (or sharply-narrowed-w) ocean-only retrain** at same cap/seed/arch
   — the one axis S1 could not touch; decisive cut for salinity-vs-fold/support.
2. **Re-measure matched MCMC-pp** for the anchor (ideally ocean-only) so the reported
   target is SBI-pp vs MCMC-pp, not vs obs.
3. **Plot the anchor SBI Im_k2 pushforward distribution** (not just median):
   unimodal-under-concentrated near 0.04 vs a mode near 0.10–0.13 — distinguishes
   "fails to concentrate" from "concentrates on the wrong mode."

Artifacts: `validation_reports/nh3_diagnosis/{capped_full_joint_anchor,s1_ocean_only,
s2_reduced_noise}/`; `validation_reports/b3_reference_wander/`. Drivers:
`plans/scripts/nh3_diag_*.py`, `plans/scripts/b3_reference_wander.py`.

**Manager (Machine A) call pending:** proceed to MgSO4/NaCl vs. run reviewer
follow-up #1 (salinity-fixed) first. MgSO4/NaCl remain HELD per §0.8 until the manager
adjudicates.

---

## Follow-up #3 EXECUTED — pushforward-shape diagnostic (2026-08-05)

Manager §0.9 (MACHINE-B-HANDOFF, 5ae7d17c) reordered the reviewer's follow-ups
to **#3 first** (FREE, no new MCMC), then #2 (matched-resolution MCMC re-measure),
then #1 (salinity-fixed retrain) only if a material gap survives #2. Follow-up #3
is done. Driver: `plans/scripts/nh3_diag_pushforward_plot.py`. Method: sample
N=4000 flow draws at x_obs, push theta back through the byte-identical forward
physics (`generate_sbi_dataset(theta_override, obs_noise=False)`), take |Im_k2|;
prior-predictive is the clean-recovered training column (one-shot additive noise
subtracted, min≥0 confirmed); MCMC reference is the weighted tidal posterior.

Run on BOTH the 60-epoch capped anchor AND the deployed 1M flow (reviewer-required
transfer check — the manager's conclusion is about the deployed flow, and the cap
alone cannot separate undertraining from non-identifiability).

| |Im k2\| | prior-pp | capped anchor | deployed 1M | MCMC ref (wtd) |
|---|---|---|---|---|
| p5 | 0.0021 | 0.0037 | 0.0034 | 0.0308 |
| median | 0.0265 | 0.0431 | 0.0423 | 0.0999 |
| p95 | 0.263 | 0.340 | 0.321 | 0.165 |
| iqr90 | 0.261 | 0.337 | 0.317 | 0.134 |
| concentration ratio vs prior | 1.0 | **1.289** | **1.215** | — |
| frac ≥ 0.100 (ceiling) | 0.176 | 0.250 | 0.233 | — |
| frac ≥ 0.135 (obs) | 0.128 | 0.173 | 0.161 | — |

**Verdict — concentration-failure CONFIRMED, wrong-mode EXCLUDED (scientific-reviewer
PASS, 2026-08-05).** Both concerns the reviewer raised on the first #3 run are closed:

- **Transfer (concern 1) — CLOSED.** The deployed 1M flow reproduces the capped
  anchor's shape: concentration ratio 1.215 vs 1.289 (both **>1** — the posterior-
  predictive is not narrower than the prior, it is marginally *broader*), high-k2 tail
  retained (p95 0.321 exceeds even the prior p95 0.263), near-zero pile retained
  (p5 0.0034). frac≥obs 0.161 (anchor 0.173), median 0.0423 (anchor 0.0431). The
  shape is a property of the trained flow, not of the cap. concentration ratio >1 is
  *stronger* than "barely narrows": the flow fails to condition this observable at all.
- **Nonfinite-drop bias (concern 2) — CLOSED on magnitude.** Only 7/3842 (0.18%)
  deployed draws drop for nonfinite forward output (anchor 9/3827, 0.22%). The dropped
  subset does lean *high-dissipation* (deployed: log10_zeta −1.64 vs kept −1.07;
  log10_eta_Ih 11.97 vs 12.89 — both directions raise |Im k2|), so the drop is
  *conservative*: including those draws would nudge SBI-pp mass toward obs, weakening
  the apparent under-update, never manufacturing it. Max conceivable shift to frac≥obs
  ~0.002 vs an under-update gap ~0.34. Closed on magnitude + conservative direction,
  NOT on kept/dropped similarity (reviewer correction — the drivers are not neutral).

**Ceiling-INDEPENDENT evidence (replaces the earlier "returns the prior" framing).**
obs=0.135 is reachable under the prior (12.8% of prior draws already exceed it). A
posterior correctly conditioned on 0.135±0.035 should heavily up-weight that mass —
qualitatively toward frac≥obs ~0.5 (heuristic target for a symmetric noiseless
pushforward centred on obs, not an exact identity). The flow moved frac≥obs only
0.128→0.161 (deployed) / →0.173 (anchor). This under-update is measured **entirely
against the prior** and does NOT depend on where the MCMC ceiling (0.100) sits — so it
survives whatever follow-up #2 finds.

**#3 does NOT obviate #2.** The dichotomy #3 tests is concentration-failure vs
wrong-mode only; it does not rule out reference/forward-model issues in the MCMC
ceiling itself. #2 (matched-resolution NH3 reference MCMC at n_eff=2000, re-measure the
0.100 ceiling — measured at the B3-discredited n_eff=500) still gates MgSO4/NaCl:
if the matched-resolution SBI-vs-MCMC-pp gap < 0.5 σ_obs → MgSO4/NaCl PROCEED (standard
gates + pushforward gate); if it survives → run #1 (salinity-fixed retrain) with
reviewer + user sign-off before any MgSO4/NaCl compute.

Artifacts: `validation_reports/nh3_diagnosis/pushforward_shape/` (anchor) and
`.../pushforward_shape_deployed/` (deployed 1M) — each has `pushforward_shape_report.json`,
`pushforward_arrays.npz` (raw prior/sbi/mcmc |Im k2| for re-plot without recompute),
`anchor_imk2_pushforward.pdf`/`.png`.

## Follow-up #2 EXECUTED — matched-resolution NH3 reference MCMC (2026-08-05)

**Question.** The 0.100 |Im_k2| MCMC-pp "ceiling" the flow was judged against was
measured at pocoMC n_effective=500 — the resolution class B3 discredited (the
v5/v7 1.06 km reference wander shrank ~5.6× at n_eff=2000). Re-measure the target
at matched resolution before spending a retrain chasing it.

**Method.** NH3 joint reference MCMC re-run at n_effective=2000 / n_active=1024
(regime ratio 1.953, byte-identical to the legacy 500/256 ratio, so flow-train
cadence + annealing preserved), seeds 72 + 172. The **tracked config was NOT
mutated** — config_hash `1611b65fff3f06c9` (referenced by the deployed flow and
the committed n_eff=500 reference) is intact; the knobs are set on the runner
instance only, exactly as B3. Both seeds annealed to β=1 (7187 samples each).
Driver: `plans/scripts/nh3_diag_matched_reference.py`. Reviewer-verified PASS
(numbers reproduced exactly, plumbing confirmed).

**Result — the ceiling is NOT a resolution artifact (decisive contrast with B3).**

| quantity | value |
|---|---|
| n_eff=500 ceiling (legacy) | 0.0999 |
| n_eff=2000 pooled median | **0.1037** (seed 72: 0.1043, seed 172: 0.1031) |
| move | +0.0038 = +0.11σ_obs, and *away* from the SBI value |
| between-seed std / range | 0.00086 / 0.0012 |

Unlike B3 (where the target statistic collapsed under the same resolution jump),
the ceiling held. So the n_eff=500 ceiling was sound.

**Four-way table (|Im_k2|).**

| quantity | value | vs obs (σ 0.035) |
|---|---|---|
| deployed SBI posterior-predictive median | 0.0423 | −3.2σ |
| matched-res MCMC posterior-predictive median | 0.1037 | −0.89σ |
| observed | 0.135 | — |

**Gate outcome — the gap SURVIVES.** SBI-pp 0.0423 vs matched MCMC-pp 0.1037 →
gap +0.0614 = **+1.76σ_obs ≫ 0.5σ_obs** threshold (robust to seed: worst per-seed
still 1.74σ). Combined with #3 (concentration-failure confirmed, wrong-mode
excluded), the under-update is a **genuine flow-training deficiency of the tidal
sector**, not a comparison/reference artifact and not mode collapse.

**Reviewer clarification (recorded).** Two distinct offsets — only one is the
remediation target:
- (a) matched MCMC-pp sits 0.89σ below obs = ordinary model/data tension (the
  physics grid cannot fully reach 0.135; expected, NOT a target — a retrain
  should not and will not force agreement with obs).
- (b) SBI-pp sits 1.76σ_obs below MCMC-pp = the flow deficiency, the actual
  target of #1.

**Decision (manager §0.9 preregistered).** Gap survived → **run #1
(salinity-fixed ocean-only retrain) with reviewer + user sign-off before any
MgSO4/NaCl compute**; MgSO4/NaCl PROCEED is falsified and those campaigns stay
HELD.

Artifacts: `validation_reports/nh3_diagnosis/matched_reference/matched_reference_report.json`,
per-seed pickles `nh3_reference_seed{72,172}_neff2000.pkl` + `prog_seed*.jsonl`.
Does NOT replace the committed n_eff=500 reference pkl (that stays the
ratification-time artifact).

## Follow-up #1 EXECUTED — salinity-fixed ocean-only pilot + matched-N control (2026-08-05)

**Manager AUTHORIZED (§0.9, 2026-08-06):** #2's gap-survives outcome satisfied the
preregistered branch, and #1 is a pilot retrain on the existing 1M dataset with no
artifact-design change — no further sign-off needed to RUN it (reviewer + user
sign-off attaches to REMEDY selection after #1).

**Question #1 separates.** S1 (ocean-only, salinity still varying) left the strong
Tb↔w salinity degeneracy (corr −0.986) intact. #1 removes it: fix salinity at the
reference posterior weighted median (log10_w = 1.1007 = 12.61 ppt) and retrain an
ocean-only pilot on that slice. Preregistered reading: Im_k2 gap COLLAPSES →
salinity axis is the driver; PERSISTS → capacity/embedding (#4).

**Design review (scientific-reviewer, PASS WITH CONCERNS 2026-08-05).** Required
before the reading could be acted on:
- MAJOR — the fixed-salinity band keeps only ~9% of rows (60,039) vs the capped
  anchor's ~690k, an ~11× cut. Under-resourced training generically biases a flow
  toward under-concentration (the same Im_k2≈0.04 signature under diagnosis), so
  a PERSIST read against the 690k anchor alone would confound "salinity removed"
  with "11× less data." Fix: also train a **sample-size-matched, salinity-VARYING
  control** (same 60,039 rows, seed 72, nsf, 60-epoch cap). Read banded-vs-control
  at matched N.
- MODERATE — the "nearest w-node IS the median node" rationale is false (the band
  is not the node's Voronoi cell, and the forward model blends bilinearly). Physics
  unbiased; re-documented as a symmetric ±0.084-dex (~±21%) salinity slice about
  12.6 ppt centered on the true median.
- MINOR — record epochs_trained + early-stop parity per pilot.
All three folded in; both pilots trained to the 60-epoch cap (61 ep, no early-stop),
config_hash e596574d1e81567c matches S1/anchor (no artifact-design drift). Banded
kept-log10_w 5/50/95 = [1.026, 1.102, 1.176] confirms centering on the true median.

**Results (PPC posterior-predictive, noiseless, n_post=4000):**

| pilot | N_train | salinity | Im_k2 pp-median | dev vs obs | frac ≤3σ | n_finite |
|---|---|---|---|---|---|---|
| banded (fixed ~12.6 ppt) | 60,039 | fixed | 0.03132 | 3.16σ | 0.429 | 3788 |
| control (varying, matched N) | 60,039 | varying | 0.03213 | 3.15σ | 0.436 | 3768 |
| S1 (varying) | 642,558 | varying | 0.03918 | 2.98σ | 0.504 | 3787 |
| capped full-joint anchor | ~690k | joint | 0.043 | — | — | — |
| matched MCMC-pp ceiling (#2) | — | — | 0.1037 | — | — | — |
| observed | — | — | 0.135 | — | — | — |

**Key deltas (σ_obs = 0.035):**
- **banded − control = +0.023σ_obs** — salinity axis IN vs OUT at matched N is
  essentially zero, and in the WRONG SIGN for the salinity hypothesis (fixing
  salinity should have raised Im_k2, not lowered it).
- control(60k) − S1(643k) = −0.201σ_obs — the pure size effect: more data moves
  TOWARD obs, exactly the under-concentration-from-data-starvation the MAJOR
  control was added to isolate. Discharges the confound.
- banded − ceiling(0.1037) = +2.07σ_obs — the still-open flow gap.

**Outcome — PERSIST. Reviewer PASS (2026-08-05).** Fixing salinity at matched N
produced no recovery of the Im_k2 update; the pure size effect is small, monotone,
and toward obs (cannot rescue salinity); the SBI-pp remains ~2σ_obs below the
matched-MCMC ceiling. **The salinity axis/degeneracy is ELIMINATED as the driver
of the ocean-branch under-update; the remaining candidate is flow capacity /
embedding (#4).**

**Two-gap scoping to surface (reviewer).** #1 addresses only gap (a):
- Gap (a): SBI-pp (0.039–0.043) vs MCMC-pp ceiling (0.1037), +1.76–2.07σ_obs — the
  flow deficiency; salinity now cleared → #4 (capacity/embedding) is the right
  remaining candidate. The banded pp 5–95 = [0.0023, 0.340] already spans past
  obs → an expressiveness/concentration signature, not a support failure.
- Gap (b): MCMC-pp ceiling (0.1037) vs obs (0.135), +0.95σ_obs — forward-model /
  prior-support / obs tension; #1 neither tests nor targets this.

The size gain is strongly sublinear — 11× data bought +0.007 in Im_k2, while
closing gap (a) to the ceiling needs +0.065 (~9× that gain) — so more ocean-only
data will not plausibly close gap (a). (Two size points, so treat the saturation
shape as heuristic; direction/magnitude robust.)

**STOPPED per protocol.** Remedy (#4) selection is a manager + reviewer + user
decision — NOT Machine B's call; MgSO4/NaCl stay HELD. Optional #4 hardening
(non-blocking, if #4 is later authorized): 2–3 additional training seeds for
banded + control to attach a formal training-noise band to the 0.023σ null; one
single-salinity-node retrain to confirm the ±21% band residual does not carry the
null.

Artifacts: `plans/scripts/nh3_diag_1_salinity_fixed.py` (driver);
`validation_reports/nh3_diagnosis/f1_salinity_fixed/f1_train_manifest.json`;
`validation_reports/nh3_diagnosis/f1_salinity_fixed/{banded,control_varying}/ppc_interior_report.json`.
Pilot .pt artifacts under `/tmp/nh3_diag/f1/` (not committed; diagnostic pilots,
not deployment artifacts).
