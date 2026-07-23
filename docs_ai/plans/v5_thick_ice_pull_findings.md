# v5 reference-MCMC findings: thick-ice pull + degeneracy sign-flip

**Author:** genai session, 2026-07-20
**Source:** scientific-reviewer agent ae91a71114b49ff13 (definitive verdict), run
against the pinned v5 cache + `europa_clipper_v5_reference_result.pkl`
(4441 samples, ESS 4125, r̂ 1.0).
**Verdict:** **NON-BLOCKING — CHARACTERIZE-AND-REPORT.** Proceed with the three
1M SBI trainings. The four items below MUST appear in each v5 gate report.

---

## Summary

Two posterior features flagged during the reference-MCMC review turned out to be
scientifically sound, not defects:

1. **Thick-ice pull** — D_iceIh posterior median ≈ 50.8–51.5 km (≈ +2.25σ into the
   truncated-Gaussian prior tail; prior μ=29, σ=10, [5,60]). This is a **genuine,
   likelihood-driven data preference dominated by the tightened Mazarico k2
   constraint** (σ 0.06 → 0.015, 4×). It is **not** a cache-edge artifact: at the
   posterior's own salinity range (w 1–99 pct = 16.3–32.7 ppt) the cache reaches
   D_max ≈ 96–101 km, far above the 59.4 km support wall. **The support wall is the
   prior's 60 km truncation, not the cache's D_max(w) envelope.**

2. **Degeneracy sign flip** — corr(D_iceIh, log10_w) = **+0.863** (v4 had
   corr(Tb, w) = −0.986). This is the **algebraically-correct** consequence of the
   reparameterization, not new information and not mis-wiring. The v4 ridge is
   preserved: recovered corr(Tb_derived, w) = −0.935, corr(Tb_derived, D_iceIh) =
   −0.984. Since D ∝ −Tb and Tb ∝ −w ⟹ D ∝ +w. Positive is the expected sign.

---

## Root cause of the thick-ice pull (physics)

The forward model's **minimum reachable Re_k2 anywhere in the 11D posterior is
0.241**, which is **+0.72σ above the 0.23 target**. Thinner ice gives *higher* k2
(D < 40 km → Re_k2 ≥ 0.252). The tightened Mazarico k2 pulls toward the smallest
achievable k2, which the model delivers only with the thickest ice. The target
k2 = 0.23 is **unreachable**, so the posterior D is set by where the model's k2
floor meets the prior tail — a soft model–data tension, **prior-truncation-
regularized, not data-localized**.

## Decisive diagnostic — profile likelihood vs D (nuisance-maximized, prior-independent)

```
 D (km):   37    41    45    49    53    57    59
 maxLL: -12.4  -8.5  -5.3  -3.5  -2.3  -1.4  -1.2
```

Rises monotonically by ~11 log-units 37→59 km and is **still climbing at the
truncation wall** (slope has not flattened). Signature of a data preference the
model cannot satisfy within the prior box — the likelihood wants D > 60 km. A
prior-volume artifact would instead show a FLAT profile likelihood. → genuine
likelihood pull.

## Flat-prior reweight (how much the prior tail is doing)

Removing the truncated-Gaussian tilt (reweight by 1/prior_pdf) on the existing 4441
samples: likelihood-only median moves 51.5 → **54.4 km**, 97.5-pct → 58.8 km. The
informative prior pulls the median *down* ~3 km toward its 29 km mean; the
likelihood alone wants even thicker ice. The prior's upper tail is **fighting** the
likelihood and losing — it is not manufacturing the pull.

---

## Four MANDATORY gate-report entries (reviewer-required, all characterize-and-report)

1. **k2 model–data tension.** Report that min reachable Re_k2 = 0.241 exceeds the
   0.23 target by +0.72σ, that Re_k2 decreases monotonically with ice thickness,
   and that THIS — not a genuine thick-ice measurement — drives the D pull. State
   plainly that the posterior D upper edge is prior-truncation-regularized, not
   data-localized (profile likelihood still rising at 60 km).

2. **Prior-truncation sensitivity.** The 60 km upper truncation is load-bearing.
   Report the flat-prior reweight (likelihood-only median 54.4 km, 97.5% 58.8 km)
   and state the reported D median (50.8–51.5 km) is a prior-mean-vs-likelihood-tail
   compromise that would shift upward under a looser truncation. (Cheap: reweight
   existing samples, no new compute.)

3. **k2 normalization confirmation (one line).** That Re_k2 cannot reach the
   ocean-consistent 0.23 anywhere in an 11D posterior warrants confirming the k2
   forward map and the 0.23 target share the same normalization/definition
   (Re/Im, sign, unnormalized-Love). Cross-check against a known Europa k2 baseline
   in the PP test suite. Pre-existing forward-model property (inherited from the
   cache), NOT a v5 regression → does not block. BUT: if a future finding shows a
   k2 normalization mismatch, the entire thick-ice pull would be reinterpreted —
   flag explicitly.

4. **Degeneracy deliverable.** Report corr(D, w) = +0.863 **together with**
   recovered corr(Tb_derived, w) = −0.935, so readers see the ridge is preserved and
   the sign flip is a reparameterization artifact, not newly-broken degeneracy. Do
   NOT present +0.863 as evidence of a new degeneracy.

---

## Verified sound (no action needed)
- Cache monotonicity + D_max(w) envelope (D_max ≈ 96–101 km at posterior salinities;
  support wall = prior 60 km, not cache).
- Sampled-vs-derived D consistency (< 1e-6 km vs `D_iceIh_results`).
- Tb↔w ridge preservation across reparameterization (−0.935 vs v4 −0.986).
- Algebraic consistency of the +0.863 sign flip.
- Convergence (ESS 4125, weighted, r̂ 1.0).
- Likelihood/observable channels vary smoothly + monotonically in D (no
  discontinuities, no NaN pathology in retained samples).
