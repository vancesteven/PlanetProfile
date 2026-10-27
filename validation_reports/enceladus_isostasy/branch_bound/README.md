# Enceladus branch comparison: frozen vs ocean, by an evidence UPPER BOUND (reviewer GO-with-bound, 2026-10-27)

- **Driver:** `plans/scripts/enceladus_frozen_branch_bound.py`, run on the Mac production cache.
- **Record:** `frozen_branch_bound.json`.
- **Rulings applied:** the reviewer pre-run ruling (2026-10-27) and the run (h) normalisation ruling.

## Method

The frozen prior integrates to 1, so **Z_frozen ≤ max_θ L(θ)** over the declared frozen prior. The maximum was found with the PRODUCTION likelihood: the ocean6D_rhorock config reduced to the frozen parameters, so the observables, σ, sigma_model_add, the C20–C22 correlation of 0.47 and the 1.008 libration factor are identical. Two searches were run:
- a dense grid of 161 ρ_rock × 21 ρ_ice × 13 sys_frac (43,953 points; 26 rejects);
- a 12-start bounded Nelder–Mead.

**Frozen construction** (reviewer-verified):
- θ with ρ_rock and no zb is dispatched to the frozen path.
- zb comes from mass via `frozen_zb_from_mass` (M = 1.08022e20 kg, R = 252.1 km).
- The production cache's 39 frozen nodes span 46.5–65.5 km. All are no-ocean, and all have exact mass.
- Gravity is rigid-supported; libration uses `ocean = False`.
- Rock model and MoI treatment are equalised with the ocean branch: uniform ρ_s ~ U[2200, 2600], exact mass closure, derived C/MR², and no MoI term in the likelihood.

**Declared priors (nominal; neither branch is renormalised):**
- Ocean: zb U[5, 45], log10 w U[−1, 2], ρ_ice U[915, 935], ρ_rock U[2200, 2600], C2 U[0, 1], and sys_frac truncated normal (σ 0.004, bounded ±0.012).
- Frozen: ρ_rock U[2200, 2600] × ρ_ice U[915, 935] × the same sys_frac prior.

## Result

| quantity | value |
|---|---|
| max log L, frozen (grid = multistart) | **−287.29** at the prior corner ρ_rock 2600, ρ_ice 915, sys_frac +0.012 |
| residuals at the maximum (σ_infl) | C20 −9.0, C22 −1.1, C30 +6.1, **libration −21.1** |
| max libration over the frozen box | 0.0312°, against observed 0.092 ± 0.0030 (inflated) |
| libration-only bound on log L | **−203.86** (support-model independent) |
| frozen prior mass outside the cache axis | 0.015 %, counted as L = 0, which keeps the bound conservative |
| log Z_ocean (run (a), nominal prior) | −7.096 ± 0.009 |
| **log B = log Z_frozen − log Z_ocean** | **≤ −280.2** (full likelihood); **≤ −196.8** (libration-only) |

Ocean-branch reporting lines, which are separate mechanisms and are not summed:
- servable fraction 0.893, so the hole fraction is 0.107 (a cache artefact; see run (h));
- MoI-window exclusions ≈ 0 (non-binding);
- effective-prior alternative for log Z_ocean: +0.113 nat. It is not applied.

## Verdict (one-sided)

**The frozen, no-ocean branch is excluded.** log B ≤ −280 under the full likelihood. The headline-robust statement is that the libration channel alone gives log B ≤ −197, which holds whether the frozen shell is rigidly supported or Airy-compensated, because libration does not use the gravity support model.

- **No ratio or probability is quoted.** e^−280 is representable in double precision, so `branch_bayes_factor`'s underflow-based `one_sided` flag is False here. The reason not to quote a ratio is **model dependence**: switching the frozen gravity support from rigid to Airy was already measured to move χ² from 1105 to 655, about 225 nat.
- **Caveats:**
  - The bound is conditional on the declared ρ_rock ≤ 2600, since the maximum is pinned at that corner.
  - The 1.008 libration factor is carried onto the frozen branch (immaterial).
  - Branch-construction asymmetry (3), the salinity dimension, is by design.
  - The two shared assumptions are unchanged: the uniform-rock interior and the borrowed 2200 floor.

Status: `verified`. Artefacts are above. The reviewer's independent probe (5103-point grid plus multistart) gave −287.29 and −287.93, matching this run.
