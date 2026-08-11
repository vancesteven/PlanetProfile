# Track 2a — information-gain-per-observable-channel: estimator design

**Status: RATIFIED (scientific-reviewer, 2026-08-10, verdict APPROVE-WITH-CHANGES,
agent a872477da30620eb0). The additive-Gaussian-noise reduction is accepted as
EXACT for this dataset; the five required changes below are folded in. The
2b/2c CANCEL decision is gated on the NH3 CMI (regen required); NaCl/MgSO4 are
reported as free context only.** (CLAUDE.md: scientific adjudication is never
self-made — this section records the reviewer's ruling verbatim in intent.)

## Reviewer ruling — required changes (folded into the estimator)

1. **Estimand is the EFFECTIVE-PRIOR CMI**, I_{p'}(θ; Im_k2 | C20,C22,Re_k2),
   where p'(θ) ∝ p(θ)·1[θ∈A] is the support-guarded, non-finite-dropped prior.
   Verified load-bearing fact: EVERY rejection (k2-support box on the NOISELESS
   forward output, no-ocean/density-inversion guards, drop_nonfinite) is a
   deterministic function of θ ALONE — independent of the noise realization ε.
   Therefore x_S | θ = N(g_S(θ), diag σ_S²) survives rejection intact, so
   H(x_S|θ) = analytic noise entropy still holds; rejection only redefines the
   marginal to p', which is exactly the distribution any flow/reference sees. NO
   prior-predictive reweighting is needed or wanted.
2. **Rank-Gaussian transform the THREE CONDITIONING channels (C20,C22,Re_k2)
   only** — H(Im|rest) is invariant under any invertible reparameterization of
   the conditioners, so this tames k-NN boundary bias (esp. Re_k2's heavy right
   tail, p999=1.40, max 1.63) WITHOUT changing the estimand. Do NOT nonlinearly
   transform Im: it must stay in linear units for the noise-entropy subtraction
   (linear z-score of Im is fine; subtract ½ln(2πe(σ_Im/s_Im)²) with σ_Im=0.035
   EXACT from stats.obs_noise.sigma). Use the SAME k for H(x₄) and H(x₃) (the
   difference does NOT cancel k-NN bias — d differs). Enforce structural
   CMI ≥ 0 (negative flags bias).
3. **CI by m-out-of-n subsampling WITHOUT replacement** (20 draws of 200k), NOT
   with-replacement bootstrap (duplicated rows create zero neighbor distances →
   downward entropy bias). Combined uncertainty band = variance CI ⊕ k-sweep
   systematic (bias) band. CANCEL 2b/2c only if the UPPER bound of the combined
   band < 0.1 nat (require positive evidence AGAINST information before
   cancelling costly downstream work).
4. **NH3 gates the decision** — composition changes g(θ), hence H(Im|rest) and
   the CMI (NaCl Im<0 frac 0.173 vs MgSO4 0.243; NH3 restructures the column
   again). Regenerate the NH3 dataset (~1.7 CPU-h) and gate on its CMI; run
   NaCl/MgSO4 now as free robustness context, NOT as a substitute.
5. **Report the Gaussian baseline as a BRACKET, not a bound** (no inequality
   forces the true CMI above the moment-matched Gaussian CMI). State that the
   estimand is an UPPER bound on the *tidal-parameter-specific* recoverable
   information (it counts Im's info about ANY of the 13 θ dims): CMI < 0.1 ⇒
   safe cancel; CMI > 0.1 is necessary-not-sufficient for 2b/2c to help.

Reviewer sanity magnitude (heuristic, NOT the gate): std(Im)=0.1075, σ_Im=0.035,
Im↔each conditioner corr < 0.06 ⇒ CMI_gauss ≈ ½ln(0.01156/0.001225) ≈ 1.1 nat
≫ 0.1. The decision likely lands far above threshold, so k-NN bias near 0.1 nat
is unlikely to be decision-critical — but the proper estimator must still run.

Required validations before accepting a number (reviewer): (1) confirm CMI ≥ 0;
(2) report H(x₄),H(x₃) across k∈{3,5,10,20}, raw vs transformed, and confirm the
z-score constant cancels in the difference to the analytic value; (3) verify the
subtracted noise entropy uses σ_Im=0.035 exactly in the same units as the
entropy; (4) CI by m-out-of-n subsampling, apply upper-bound<0.1 to the combined
band; (5) gate on NH3, report salts as context.

## What §0.18 Phase-1 item 4 asks

> Track 2a (information-gain in nats) on spare cores — FREE, and now
> decision-theoretic: report achievable-KL per observable channel; if Im k2
> carries **<~0.1 nat beyond {C20,C22,Re_k2}**, CANCEL 2b and 2c.

So the decision quantity is the **incremental information Im_k2 carries about
the interior θ beyond the other three channels** — i.e. the conditional mutual
information

    CMI = I(θ ; Im_k2 | C20, C22, Re_k2)   [nats]

with threshold ~0.1 nat. (Track 2's purpose, §0.12/plan: "quantify how much
signal a perfect posterior has — an upper bound on any flow's update." A CMI
below threshold means no amount of training/transform can recover a tidal
update, so 2b/2c are pointless.)

## Data available (this session)

- `/tmp/titanG_build/datasets/titanG_nacl_1m.npz` — θ (631214, 13), x (631214, 4)
- `/tmp/titanG_build/datasets/titanG_mgso4_1m.npz` — θ (631214? , 13), x (…, 4)
- **NH3 npz was CLEARED** — regen (~1.7 CPU-h) needed to run NH3's 2a; flag.

x columns are the config observable order {C20, C22, Re_k2, Im_k2}. Verified
from the NaCl gen manifest: the stored x is **NOISED** —
`obs_noise: gaussian_diagonal`, fixed σ = {C20 4.71e-7, C22 6.2e-8,
Re_k2 0.048, Im_k2 0.035}, `correlations: null`, `imag_convention: abs`,
`refold_im: false`. The abs fold is applied to the CLEAN forward Im_k2, THEN
independent Gaussian noise is added (this is why stored Im_k2 has negative
values, min −0.154). So x = g(θ) + ε, ε ~ N(0, diag σ²), ε ⟂ θ, where
g(θ) = (C20, C22, Re_k2, |Im_k2|)(θ) is a fixed deterministic map.

## Proposed estimand + estimator (RECOMMENDED — for the reviewer to accept/replace)

### The additive-Gaussian-noise reduction (exact under this noise model)

Because ε is additive, diagonal-Gaussian, and θ-independent, the conditional
law x|θ is N(g(θ), diag σ²) whose differential entropy is the CONSTANT noise
entropy regardless of θ. Hence for any channel subset S:

    I(θ ; x_S) = H(x_S) − H(x_S | θ) = H(x_S) − Σ_{c∈S} ½ ln(2πe σ_c²)

(the conditional entropy collapses to the analytic noise entropy). By the chain
rule the decision quantity is

    CMI = I(θ; all 4) − I(θ; {C20,C22,Re_k2})
        = [H(x₄) − H(x₃)] − ½ ln(2πe σ_Im²)
        = H(Im_k2 | C20,C22,Re_k2)  −  ½ ln(2πe σ_Im²)

So the whole thing collapses to a **≤4-dim differential-entropy estimation on
the noised-observable marginal** — NOT a 13+4-dim MI problem. The 13-dim θ
never enters the estimator; it is integrated out exactly by the noise-model
identity. This is the load-bearing claim the reviewer must check.

### Estimator: Kozachenko–Leonenko k-NN differential entropy

- Estimate H(x₄) and H(x₃) with the k-NN (KL) differential-entropy estimator on
  z-scored channels (z-scoring adds a known analytic constant that cancels in
  the H₄−H₃ difference; report both raw and z-scored to confirm cancellation).
- Sweep k ∈ {3,5,10,20} for the k-NN bias check; report the plateau.
- Subtract the analytic noise entropy ½ ln(2πe σ_Im²) (σ in the SAME z-scored
  units used for H, so the subtraction is consistent).
- Report per-channel: (a) marginal I(θ; x_c) for each of the 4 channels, and
  (b) the incremental CMI of each channel given the other three (Im_k2 is the
  decision one; the other three reported for context/ranking).
- Uncertainty: block-bootstrap over rows (say 20 resamples of 200k) → CI on the
  CMI; the 0.1-nat decision is made on the CI, not the point estimate.

### Cross-checks (independent estimator + analytic floor)

1. **Gaussian-copula / linear-Gaussian analytic baseline.** Fit a joint Gaussian
   to the 4 z-scored channels; the CMI has a closed form
   (½ ln[ (1) / (1 − ρ²_{Im·rest}) ] style partial-correlation expression).
   This is a LOWER-fidelity floor (ignores non-Gaussian structure the §0.16
   pushforwards showed is real) but a useful sanity bracket: k-NN CMI should be
   ≥ the Gaussian CMI if the extra Im structure is informative.
2. **KSG CMI** estimated directly on (θ_reduced, Im_k2 | rest) as a redundant
   check IF the reviewer wants the θ side kept explicit rather than trusting the
   noise-model reduction — but θ is 13-dim; would need a sufficient reduction
   (e.g. the tidal sub-block log10_eta_*, log10_zeta, alpha), which introduces
   its own approximation. The reduction-based KL estimator above AVOIDS this,
   which is its main advantage.

## Questions for the reviewer

1. Is the additive-Gaussian-noise reduction accepted as EXACT for this dataset
   (so the estimand is purely the low-dim noised-observable entropy difference)?
   The abs-fold-before-noise wrinkle: does it threaten the reduction? (Claim: no
   — folding changes g(θ), not the noise model, so I(θ;x_S)=H(x_S)−H_noise still
   holds. Confirm.)
2. Is Kozachenko–Leonenko k-NN the right primary estimator at n~6e5, d≤4, or is
   a KDE/other estimator preferred? Any bias concern near the 0.1-nat threshold?
3. Is the 0.1-nat threshold applied to the CMI **point estimate** or the
   **bootstrap CI** (my proposal: the decision is CANCEL 2b/2c only if the CI
   UPPER bound is < 0.1 nat, i.e. conservatively require evidence AGAINST
   information before cancelling costly downstream work)?
4. Should this run per-composition (NaCl + MgSO4 now; NH3 after ~1.7 CPU-h
   regen), or is the NH3 channel — where the quarantine actually lives — the
   only one whose 2a verdict gates 2b/2c? (NH3 is the quarantined comp; salts
   are split-status. Regen NH3 or accept salt CMI as a proxy?)
5. Any objection to reporting the OTHER three channels' incremental CMI too
   (cheap, and it contextualizes whether Im_k2 is uniquely uninformative vs the
   whole vector being weak)?

## RESULTS — salt context (NaCl, MgSO4), 2026-08-10

Estimator: `plans/scripts/track2a_infogain.py`; report
`validation_reports/track2a_infogain/track2a_cmi.json`. Run on the two
datasets present (NH3 npz cleared → NH3 gate pending regen).

| comp | n_kept | CMI k-sweep [k3..k20] (nat) | combined band (nat) | Gaussian bracket | CMI≥0 | CANCEL 2b/2c |
|---|---|---|---|---|---|---|
| NaCl  | 631214 | [0.6249, 0.6309] | [0.616, 0.640] | 1.120 | yes | **No** |
| MgSO4 | 691075 | [0.3877, 0.3944] | [0.380, 0.402] | 0.737 | yes | **No** |

Cross-checks (reviewer required validations, all satisfied):
- **CMI ≥ 0** structural check passes for both (min-over-k 0.625 / 0.388).
- **Raw (z-only) vs rank-transformed CMI agree to < 0.01 nat** (NaCl raw
  0.635–0.637 vs transformed 0.625–0.631; MgSO4 raw 0.397–0.398 vs
  0.388–0.394) — the conditioner reparameterization does not move the
  estimand, as required.
- **k-sweep stability**: H(x₄), H(x₃) each move < 0.01 nat across
  k∈{3,5,10,20}; m-out-of-n subsample std ≈ 0.002–0.003 nat.
- **Noise entropy** subtracted with σ_Im = 0.035 exact (stats.obs_noise),
  placed in the z-scored Im units (σ_Im/s_Im).
- **Gaussian value reported as a BRACKET** — the k-NN CMI sits below the
  moment-matched Gaussian, consistent (Gaussianization inflates each
  entropy block); NOT treated as a violated bound.

Interpretation (salts are CONTEXT, not the gate — reviewer ruling 4):
Im_k2 carries **substantial** incremental information about the interior
beyond {C20,C22,Re_k2} for both salt compositions — ~4–6× the 0.1-nat
cancel threshold. On this evidence the transform-retrain (2b) and
Re/Im-only ablation (2c) are **NOT cancellable** for the salts. The 2b/2c
CANCEL decision is gated on the **NH3** CMI (composition changes g(θ) hence
H(Im|rest)); NH3 requires the ~1.7 CPU-h npz regen before its CMI can be
computed. Given the salts land far above threshold, the expected NH3
outcome is also above threshold (i.e. 2b/2c stay on the table), but that is
NOT decided until the NH3 CMI is measured. Reviewer caveat: this CMI upper-
bounds the *any-θ-dimension* recoverable info, so > 0.1 nat is
necessary-not-sufficient for 2b/2c to help; it forbids CANCEL, it does not
by itself promise a recoverable tidal update.

## Non-gate framing

Per §0.18: 2a is FREE and runs on spare cores; it does NOT touch the manager
gate (C3/C5.3/C10/C16) and reads NO corrected/IS-reweighted quantity — only the
training (θ, x) pairs. It cannot advance or block Phase-1 Track-1 compute; its
only consequence is the 2b/2c CANCEL decision. Safe to run while the NH3 seed
272 reference MCMC anneals.
