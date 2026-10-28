# Can a no-ocean Enceladus with low-viscosity ice reach the observed libration?

Asked by the user on 2026-10-28, prompted by Downey & Nimmo (2025, *Sci. Adv.* 11, eadl4741), `papers/downey2025titans.pdf`. Their shell–interior coupling parameter K/C_s is the quantity used here.

- **Driver:** `plans/scripts/enceladus_noocean_libration_bound.py`, using PlanetProfile's own `Gravity/Librations.py`.
- **Record:** `noocean_libration_bound.json`.
- **Target:** the observed libration of 0.092 ± 0.003°, i.e. ≥ 0.083° (observed − 3σ).

## (A) Whole body, no decoupled layer

Van Hoolst et al. (2013) oceanless branch, `librations(ocean=False, rigid=False, rek2)`, with Re(k₂) scanned over [0, 1.5]. The cases are ρ_rock ∈ {2200, 2400, 2600} and zb ∈ {20, 45, 65} km.

| | value |
|---|---|
| rigid amplitude | 0.0281–0.0306° |
| maximum over all Re(k₂) | the rigid value |

Tidal deformation only lowers the amplitude: it shrinks the effective B − A, falling about 0.0012° per 0.1 in k₂. The analytic estimate 6e(B − A)/C ≈ 0.030° agrees.

**A deformable frozen body peaks at 0.031°, about a third of the observed value.**

## (B) Decoupled ice shell over a soft (viscous or Maxwell) ice layer

**Stack.** Three layers: uniform rock core, soft ice layer of thickness h, ice shell of thickness zb, mass-closed to Enceladus.

**Physics.** The Van Hoolst et al. (2008) rigid-shell system, plus a frequency-domain coupling torque transmitted by the soft layer. That torque is the exact Stokes torque between concentric rotating spheres, K̃ = 8π iη_eff r₁³r₂³/(r₂³ − r₁³), with a Maxwell rheology η_eff = ημ/(μ + inη) and μ = 3.3 GPa.

**Check.** At zero coupling the solver reproduces `librations(ocean=True, rigid=True)` to a relative difference of 5 × 10⁻¹⁶.

**Coverage.** 38 feasible cases: zb 5–45 km, h 1–100 km, ρ_rock 2200–2600, η 10³–10¹⁸ Pa·s.

| result | value |
|---|---|
| cases reaching ≥ 0.083° at ANY η | 15 of 38, all with zb ≤ 20 km and h ≤ 10 km |
| largest viscosity that reaches it (η*) | **8 × 10⁶ Pa·s** |
| maximum amplitude at η = 10¹³–10¹⁴ Pa·s (ice near its melting point) | **0.058°** |
| near-melting ice relative to η* | **6–7 orders of magnitude too viscous** |
| at zb = 27 km (the headline thickness), any h, any η | ≤ 0.078°, never reaching 0.083° |

Viscosities as low as 10⁶–10⁷ Pa·s are reached only by liquid water or slush with a high melt fraction, not by solid ice. That regime is an ocean or mush model, not "no ocean".

## Biases (each one makes the soft-layer amplitude too LARGE, so the bound is conservative)

- **The shell is rigid.** Shell elasticity lowers the amplitude (Van Hoolst et al. 2013).
- **The soft layer's own inertia is not carried.** In the strong-coupling limit this gives up to 0.055°, against the true whole-body 0.029°. It inflates the η ≳ 10¹² amplitudes; it does not deflate them.
- **The soft layer has ice density (925 kg/m³).** An actual ocean (about 1005 kg/m³) has larger pressure coupling. Even in the inviscid limit, a 925 kg/m³ layer at zb = 27 km gives only 0.078°. The observed value needs the ocean–ice density contrast, consistent with the ocean-branch fit.

## Conclusion

**No.** A frozen Enceladus librates as one body at about 0.03° whatever its viscoelastic state. Decoupling the shell needs a sub-surface layer with an effective viscosity below about 10⁷ Pa·s, which is only possible for a liquid or a slush with high melt fraction. Ice at its melting point (10¹³–10¹⁴ Pa·s) couples the shell to the interior on the 1.37-day forcing period and caps the amplitude at ≤ 0.058° even with every bias in its favour.

This agrees with Downey & Nimmo's framing: liquid-ocean boundary layers give K/C_s ≈ 10⁻¹³–10⁻¹⁰ s⁻¹ ≪ n, while a solid layer gives K/C_s ≳ 1 s⁻¹ ≫ n. It strengthens the branch-comparison verdict (`../branch_bound/README.md`): that verdict's libration-only bound assumed a rigid frozen body, and this shows that allowing viscoelastic or low-viscosity ice does not rescue it.

**Status:** `verified` for the computation (check above). Not yet reviewed by the scientific-reviewer.
