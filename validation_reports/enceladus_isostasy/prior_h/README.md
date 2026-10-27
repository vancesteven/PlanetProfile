# Run (h): prior-only pushforward, by direct prior Monte Carlo (Fisher-checkpoint condition 3)

- **Driver:** `plans/scripts/enceladus_prior_pushforward.py` (committed 2026-10-27).
- **Inputs:** the run (a) config (`ocean6D_rhorock`) and the Mac production cache.
- **Draws:** 20,000 prior draws, numpy `default_rng(75)`, about 6 ms per draw, 2 min wall.
- **Pipeline:** every draw goes through the same runner stages as the likelihood: cache lookup, then isostatic gravity with mass closure, then libration including the model budget.
- **Classifier check:** on 600 draws the classification matched `log_likelihood_fn`'s finite/reject outcome for every draw.
- **Outputs:** `prior_pushforward.json` (summary) and `prior_pushforward_draws.npz` (all draws, 2.0 MB, committed).

## Support: the hole and physical rejects, tallied separately

| category | draws | fraction |
|---|---|---|
| ok (all observables finite) | 17,857 | 89.3 % |
| **hole**: cache None-corner (a numerical artefact, not physics) | 2,143 | **10.7 %** |
| physical reject (declared support edges, incl. R_s < R_o) | **0** | 0 % |
| exception | 0 | 0 % |

**The hole.** It sits entirely at zb 5.0–13.0 km and gets steeply worse the thinner the shell. The hole fraction by zb band:

| zb band (km) | [5,7) | [7,9) | [9,11) | [11,13) | [13,15) |
|---|---|---|---|---|---|
| hole fraction | 0.82 | 0.68 | 0.42 | 0.23 | 0.00 |

These are the 121 brentq bracket-failure nodes of the r7 rebuild, where the Tb root-find has no sign change at thin shells. They are a known model-reachability restriction (`ocean_reachability_restriction`), not a physical edge.

**The R_s < R_o edge (an ocean must exist) never binds.**
- Zero rejects in 20k draws.
- The thinnest ocean drawn is 2.7 km, and only 0.34 % of draws have an ocean under 5 km.
- This matches the cache-wide scan (3399/3399 feasible at ρ_s = 2200 and 2600). The effective joint prior therefore equals the nominal prior with the hole removed.
- The hole is a region in (zb, w), bounded by the reachability limit zb_min(w) and quantised to grid cells. It is not a cut in zb alone.
- At zb 5–7 km the hole fraction is 0.56, 1.00 and 0.89 for log10 w in [−1, 0), [0, 1) and [1, 2). At zb 11–13 km it is 0.26, 0.36 and 0.07.

**Effective vs nominal prior.**
- Only zb changes: mean 25.0 → 27.0 km, sd 11.5 → 10.5. That change comes entirely from the hole removing thin shells.
- The marginals of ρ_ice, ρ_s, C2 and sys_frac are unchanged to ≤ 0.01 sd.
- The w marginal shifts by 0.024 sd, and the hole tilts the joint (zb, w) prior.
- **Acceptance figures.** This run's servable fraction is 0.893 ± 0.002 (binomial, n = 20,000). Run (a)'s 0.884 is its 1,000-draw preflight (`posterior_a/preflight.json`, 884/1000, ± 0.010). The two agree within 1σ of the preflight's own error, so there is no discrepancy. Use 0.893 as the measurement.

## Prior pushforward of the observables vs data uncertainty

| observable | prior pushforward 2.5 / 50 / 97.5 % | observed ± σ | prior sd / σ_obs |
|---|---|---|---|
| C20 | −6.18e-3 / −5.73e-3 / −5.17e-3 | −5.477e-3 ± 3.7e-5 | 7.2 |
| C22 | 1.468e-3 / 1.551e-3 / 1.634e-3 | 1.518e-3 ± 1.5e-5 | 3.2 |
| C30 | 1.12e-4 / 2.83e-4 / 3.82e-4 | 1.78e-4 ± 3.3e-5 | 2.2 |
| libration (°) | 0.061 / 0.092 / 0.374 | 0.092 ± 0.003 | 33 |

- All four data are well inside the prior pushforward, so there is no prior–data conflict.
  - The observed values sit at prior-pushforward percentiles of 0.80, 0.27, 0.12 and 0.50 (C20, C22, C30, libration).
  - The posterior's best χ² is 1.0 for 4 data, so there is no joint conflict either.
- All four are informative. Libration is the most informative, which is consistent with zb narrowing 11× from prior to posterior.
- The libration pushforward is strongly skewed. Its robust 16–84% half-width is 0.054°, or 18× σ_obs, against 33× from the sd.
- C30 is the least informative, at 2.2× σ_obs.

## Derived C/MR²: answering the reviewer's condition

The reviewer's instruction was "Do not present C/MR² as a measurement until (h) shows the prior pushforward is much wider."

- **Prior pushforward:** C/MR² = 0.333 ± 0.0052 (95 % [0.3251, 0.3424]).
- **Posterior (run (a)):** C/MR² = 0.3324 ± 0.0017 (95 % [0.3288, 0.3353]).
- **The data narrow C/MR² 3.1×.**

**That information is carried by ρ_s.**
- In the prior, corr(ρ_s, C/MR²) = −0.98. At fixed ρ_s (10 kg/m³ bins) the prior C/MR² sd is only ~0.0009.
- ρ_s itself narrows 2.9× (116 → 40 kg/m³).
- So about 96% of the C/MR² variance is the ρ_s constraint, expressed through the mass closure. The spread at fixed ρ_s also narrows, from 0.0009 to 0.0003, through zb.
- It is not an independent moment-of-inertia measurement.
- It also inherits the uniform-rock assumption and the shared Cassini J2/C22 origin of the literature MoI values (see `moi_literature_2026-10-25`).

**Wording (reviewer-adjudicated 2026-10-27):**

> Under the uniform-rock, Airy-isostatic H&M (2019) model with ρ_s ~ U[2200, 2600] kg m⁻³, the derived C/MR² is 0.3324 ± 0.0017. That is 3× narrower than its prior pushforward of 0.333 ± 0.005. The narrowing comes almost entirely from the rock-density constraint through mass closure: corr(ρ_s, C/MR²) = −0.98, and ρ_s accounts for about 96% of the C/MR² variance. It is therefore a model-dependent derived quantity, not an independent moment-of-inertia measurement, and it should not be combined with or compared to literature MoI values as independent evidence. Those values derive from the same Cassini degree-2 gravity field.

## Consequences

- **For the posterior (run (a)):** none. The likelihood is negligible above the hole (thin check: best loglike −318 at zb 13–15 km), and restricting the pushforward to zb ≥ 13 km leaves C/MR² unchanged (median 0.33286 vs 0.33305, same sd).
- **For SBI training (reviewer-adjudicated 2026-10-27): declare the effective prior; do not fill the corners.**
  - The effective prior is the nominal prior restricted to the servable (zb, w) region, meaning cells with all four corners built.
  - That region is recorded alongside the training set.
  - Prior draws in the hole are discarded, not given a reject label.
- **log Z (reviewer-adjudicated 2026-10-27).** Run (a)'s log Z is computed over the nominal prior with L = 0 in the hole.
  - L in the hole is negligible: the thin check gives best log L −318 at zb 13–15 km. So this already equals the evidence under the nominal prior, and it is **not** renormalised. It is not a bias.
  - The two-run evidence comparison uses nominal-prior log Z for both branches.
  - log Z under the effective prior would be larger by ln(1/f_servable): 0.113 nat for (a), 0.132 for g1 (hole 12.4%) and 0.115 for g2 (10.9%). Apply that factor only if an effective-prior evidence is ever wanted.

## Status

- Run (h) prior pushforward and support tally: `verified`. The artefacts are above, and the classifier was cross-checked against the likelihood.
- The C/MR² wording, SBI-prior handling and log Z handling were reviewer-adjudicated on 2026-10-27 (ACCEPT WITH EDITS; all edits applied).
- The reviewer's second independent check reproduced every number above. Its classifier cross-check on 400 fresh draws (seed 12345) gave 359 ok and 41 hole, with 0 mismatches.
- **Remaining before scientific acceptance of (a)/(g)/(h):**
  - injections, including one truth at log10 w ≈ 1.8, where the hole and the w-dependence interact;
  - ρ_ice prior variants;
  - Fisher at the MAP;
  - the Park erratum libration source into papers/;
  - the 2200 ρ_s floor stated as an open assumption. It does not bind: the posterior 95% range is [2345, 2501].
