# C3 root cause: Enceladus ocean-node mass deficit (2026-10-21)

- **Method:** an independent Opus investigation using instrumented single-node builds.
- **Environment:** Linux aarch64 container, Perple_X tables from GitHub. The builds reproduce the cached residuals exactly.
- **Files:** the script is `diag.py`; the logs are `b70_x*.log` (node (70,18): zb 39.5, w 2.424) and `b77_x*.log` (node (77,10)). The logs record the per-candidate numbers.
- **Status:** a diagnosis, not a ruling. The scientific-reviewer adjudication follows.

## (A) PP candidate selection, a modelling convention

- With Fe_CORE False, `FindInnerWithMoIAndEOS` (LayerPropagators.py ~1990-2060) sweeps porosity (nPoros=8, so 9 φ values).
- Per φ it keeps the FIRST candidate with Mtot ≤ M (Silicates.py:95-118).
- All 9 pass the ±0.08 CMR2 window, and the selection is argmin|CMR2 − 0.335| (~2097), which takes no account of mass.
- At node 70 x1 the nine candidate deficits are −0.085, −0.411, −0.358, −0.378, −0.464, −0.146, −0.347, −0.131 and −0.427%. The pick was −0.464% (φ 0.337).
- A candidate only +0.0011% over M was discarded by the ≤ M rule.
- The step between candidates is ~0.465-0.47%.

## (B) cache_builder mass convention, the floor

- `Mtot_achieved` (cache_builder.py ~576-590) integrates `Planet.Gravity.ALMAModel`, the REDUCED model.
- `GetReducedLayers` (reducedPlanetModel.py:178-188) linearly resamples the silicate layers (120 → 50) and the ocean (→ 50) using outer-edge density. Rock density rises inward (2236 → 2397 kg/m³), so mass is lost systematically.
- The loss is −0.044% at x1 and −0.042% at x4, independent of deltaP.

## Why C20/C22 move and libration does not

- `_derive_gravity_isostatic` (mcmc_runner.py:2240+) is dominated by a hydrostatic Tricarico 3-layer term normalized by the achieved stack mass.
- σ_obs as a fraction of the observed value: C20 0.66%, C22 0.94%, libration 4.6%.
- Closing node 70's mass by rescaling densities moves C20 +0.46σ (uniform) to +0.90σ (rock only), and C22 −0.34σ to −0.64σ. Libration moves ≤0.06σ.

## Shift between builds

- Node 70, x1 → x2: C20 moves +1.23σ. About 0.6σ is mass closure and about 0.6σ is a discrete candidate jump (φ 0.337 → 0.322, R_sil 197.93 → 197.19 km, ρ_sil 2301 → 2327).
- Node 77, x2 → x4, same φ: the deficit fell by 0.119 points but C20 moved only 0.044σ.
- C20 and C22 are therefore a jagged function of (zb, w) through the discrete φ/R_sil choice, and the C3 study's extrapolation to zero deficit is not valid.

## Remedy options (unadjudicated)

1. Make the reduction mass-conserving (volume-weighted), or use full-resolution arrays for Mtot_achieved and the gravity stack. Then do a full rebuild.
2. Make candidate selection mass-aware. This changes PP selection semantics.
3. deltaP/4.
4. More nPoros.
5. Add model σ for C20/C22.
6. Normalize gravity by Bulk.M instead of the stack mass.

## Scientific-reviewer ruling (2026-10-21): PASS WITH CONCERNS

The reviewer checked code paths only and did not re-run the numbers.

1. **(A) holds.** Silicates.py:110 and :118-121 keep the first candidate with Mtot ≤ M; LayerPropagators.py:2135-2137 then picks argmin|CMR2 − Cmeasured|, which ignores mass. This is a PP MODELLING CONVENTION. Remedy 2 (mass-aware selection) crosses the "fundamental assumptions" line and **needs explicit user approval**.
   **(B) holds.** cache_builder.py:274-283 and :577-579 integrate the reduced ALMAModel, and reducedPlanetModel.py:196-202 resamples r and ρ linearly using outer-edge ρ. This is a numerical defect. A fix inside the Inference builder is allowed; editing GetReducedLayers changes ALMA Love numbers for all PP users and needs approval.
   **Normalization confirmed:** mcmc_runner.py:2373 sets M_body to the reduced 3-layer stack mass, used at :2410 and :2420.
2. **Remedy (medium confidence).**
   - **Primary:** close the mass of the Enceladus 3-layer stack in the runner by rescaling rock density, as `derive_silicate_density` already does for the Clairaut campaigns. This enforces PP's own target M = Bulk.M and needs no rebuild.
   - **Then:** add σ_model for C20/C22 equal to the measured jitter.
   - **Escalate to nPoros/deltaP only if the gate fails.**
   - Remedy 1 becomes optional and diagnostic-only. Remedy 6 alone is REJECTED as inconsistent.
   **Gate that replaces C3:**
   - Take 8-10 fixed nodes, including 70 and 77. Rebuild each over deltaP {x1, x2, x4} × nPoros {8, 12, 16}, with mass closure applied.
   - PASS requires all of:
     - the C20/C22 spread is ≤ 0.5 σ_obs, or σ_model = std is ≤ 0.5 σ_obs in quadrature;
     - the libration spread is ≤ 0.1 σ;
     - no node departs from a local quadratic fit in (zb, w) by more than 3× the jitter std;
     - posterior medians shift by < 0.2 posterior σ between σ_model on and off.
3. **Blast radius.** Titan NH3, Titan NaCl/MgSO4 and Europa v5/v6 use `clairaut_hydrostatic`, normalized by Mtot_kg = Bulk.M_kg (cache_builder.py:608; mcmc_runner.py:1591, 2564/2578, 2598). Their deficit goes into ρ_sil, not into q, so they are **NOT exposed**, and **the NH3 flag flip is NOT blocked**. Residual exposures:
   - The discrete R_oceanbot enters through the hydrosphere layers. It is bounded by deltaP and was present at ratification.
   - k2 comes from the reduced ALMA model; its effect is estimated ≪ σ but unmeasured.
   - Suggested follow-up: a Europa spot-check of the C20/C22 spread across deltaP x1/x2, against σ/obs of 0.19% / 0.15%.
4. **May proceed:** libration-only Enceladus analyses; Titan and Europa work, including the NH3 flag; the two-run reference protocol and the mass-closure implementation, with Enceladus C20/C22 labelled provisional.
   **Frozen:** C3 acceptance; any Enceladus posterior conditioned on C20/C22; extrapolated corrections; edits to Silicates.py, LayerPropagators.py selection or reducedPlanetModel.py without user approval.
