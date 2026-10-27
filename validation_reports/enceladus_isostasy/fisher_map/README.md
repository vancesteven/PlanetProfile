# Fisher at the MAP (Fisher-checkpoint condition 7)

Run on 2026-10-27 with `plans/scripts/enceladus_fisher_identifiability.py`, unchanged from the Fisher checkpoint, against the Mac production cache. It evaluates two points of the run (a) posterior:
- **maxL_a**: the maximum-likelihood posterior sample, with log L −0.52 (χ² ≈ 1.0 for 4 data). Its coordinates are zb 26.473, log10 w 0.276, ρ_ice 920.57, ρ_s 2402.2, C2 0.8122 and sys_frac −0.0074.
- **postmean_a**: the posterior mean.

Raw output is in `fisher_identifiability.{md,json}`. The headline comparison uses variant C: data Fisher, plus the sys_frac Gaussian prior, plus uniform-box width proxies. That is the like-for-like linearised counterpart of the posterior under the declared priors.

| param | posterior sd (run a) | Fisher σ at max-L (post/Fisher) | Fisher σ at post-mean (post/Fisher) |
|---|---|---|---|
| zb_km | 0.9408 | 0.9036 (1.04) | 0.9277 (1.01) |
| log10_wOcean_ppt | 0.8687 | 0.866 (1.00) | 0.8659 (1.00) |
| rho_ice_kgm3 | 5.78 | 5.767 (1.00) | 5.766 (1.00) |
| rho_rock_kgm3 | 40.21 | 37.15 (1.08) | 37.97 (1.06) |
| compensation_C2 | 0.09099 | 0.08957 (1.02) | 0.08922 (1.02) |
| libration_sys_frac | 0.003879 | 0.004 (0.97) | 0.004 (0.97) |
| C/MR² (derived) | 0.001676 | 0.001611 (1.04) | 0.001618 (1.04) |

| quantity | Fisher at max-L | Fisher at post-mean | posterior |
|---|---|---|---|
| corr(ρ_s, C2) | −0.886 | −0.884 | −0.877 |
| corr(zb, log10 w) | 0.04 | 0.07 | 0.20 |

**Finite-difference step stability** (relative change of the whitened Jacobian columns between ×1 and ×2 steps): ≤ 0.6% for zb, ρ_ice, ρ_s, C2 and sys_frac, and 2.4–3.2% for log10 w. The larger log10 w figure comes from the cache's log-w node spacing; it is not a defect.

## Reading

- **Widths.** At the posterior's own location, the linearised widths reproduce every posterior sd to within 0–8%. The pre-MCMC Fisher at fixed points gave 1.05–1.19×, so the residual non-Gaussianity at the MAP is small.
- **ρ_s is the widest-ratio parameter.** It is about 6–8% wider in the posterior, which is the expected curvature of the R_s ∝ (ρ_s − ρ_o)^(−1/3) closure.
- **ρ_s–C2 degeneracy.** It is linear: corr −0.88 in both Fisher and posterior.
- **zb–w coupling.** It is weaker in the local Fisher (0.04–0.07) than in the posterior (0.20). This is expected and consistent with run (g): the zb shift with w is concentrated at log10 w ≳ 1.5, away from the MAP, where the local Jacobian cannot see it. The global coupling is the w-dependent dzb/dw, which the posterior samples and a point Fisher does not.
- **Conclusion.** The posterior is close to Gaussian about its mode, and the sampler widths are confirmed by an independent linearised calculation. Condition 7 is met.

Status: `verified` (artifacts above).
