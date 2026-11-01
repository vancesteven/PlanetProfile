# Enceladus forced libration vs ice-shell rheology (2026-10-31)

User question (2026-10-30): "libration depends on viscosity though?" — and
the k2 direction: future missions will measure Enceladus k2
(papers/bagheri2025exploring.pdf: k2 to < 0.002 with X-band DSN, < 0.0002
with Ka-band).

Status: computation `verified` (script + JSON below; reproduced by the reviewer). The BOUND is NOT accepted: Opus scientific-reviewer 2026-10-31 — claim 1 (ALMA units) CONFIRMED, claim 2 (bound) PARTLY: the Maxwell-only bound is not conservative (see Review).

## Method

`plans/scripts/enceladus_libration_rheology_bound.py`. Production (zb, w)
cache nodes at w = 3.46 ppt (posterior median salinity) and zb = 20, 26.7,
35 km. Deformable libration = Van Hoolst et al. 2013 branch
(`librations(rigid=False, y1=...)`) with y1 from a pyALMA3 viscoelastic
solve at the 1.37 d forcing, the recipe of
`tests/librations_test.py::_cache_node_diagnostics`. The ice shell's shear
modulus mu and Maxwell viscosity eta are scanned; results are RELATIVE to the
as-built elastic node (as-built ice mu 3.5–4.7 GPa). The ratio never uses the
rigid branch, so it is independent of the open K_int defect — but it ASSUMES
the deformable branch's mu/eta dependence is right, which the wrong-sign
deformable/rigid ratio does not establish. Only Re(y1) enters the libration,
so for eta ≲ 1e13 (large Im k2) the response phase is ignored. zb shift =
f · L / (dL/dzb), slope from a 7-node least-squares fit of the rigid
libration (dlnL/dzb −5.9, −3.3, −2.5 %/km); the rigid branch enters only
here. Sign: a positive dzb means a softer shell than modelled makes the
rigid-model zb too THICK by that amount.

## Result

Fractional libration change vs the as-built elastic shell (zb shift, km):

| shell | zb 20 km | zb 26.7 km | zb 35 km |
|---|---|---|---|
| elastic, mu 4.0 GPa | −0.15 % (+0.03) | −0.10 % (+0.03) | −0.07 % (+0.03) |
| elastic, mu 3.3 GPa | −0.75 % (+0.13) | −0.52 % (+0.16) | −0.37 % (+0.15) |
| elastic, mu 2.0 GPa | −2.9 % (+0.50) | −2.0 % (+0.63) | −1.5 % (+0.58) |
| elastic, mu 1.0 GPa | −8.0 % (+1.36) | −5.7 % (+1.75) | −4.1 % (+1.64) |
| Maxwell mu 3.3 GPa, eta ≥ 1e14 | = elastic 3.3 GPa | = | = |
| Maxwell mu 3.3 GPa, eta 3e13 | −1.1 % (+0.19) | −0.7 % (+0.22) | −0.5 % (+0.18) |
| Maxwell mu 3.3 GPa, eta 1e13 | −4.3 % (+0.73) | −2.3 % (+0.72) | −1.3 % (+0.50) |
| Maxwell, eta 1e12 (whole shell) | −73 % | −58 % | −37 % |

- Viscosity matters only when the WHOLE shell is near or below the Maxwell
  time (eta ≲ 3e13 Pa s at mu 3.3 GPa; tau_M = P at 3.9e14). For realistic
  shells (cold, elastic lid; at most a warm base) the effect is set by mu.
- MAXWELL ONLY: for mu 2–4 GPa and eta ≥ 3e13 Pa s the libration moves by
  ≤ ~3 % and zb by ≤ ~0.7 km. This is NOT a conservative bound (reviewer):
  at node 1220, Andrade alpha 0.3 gives −2.3 % (1.25e14) and −3.4 % (3e13);
  Andrade alpha 0.2, mu 2 GPa, 2.5e13 gives −6.4 % (~1.3–2 km in zb); an
  elastic lid over a 1e12 Pa s base gives −3.3 % (50 % lid) and −1.8 %
  (70 % lid). Realistic bound: up to ~6 % / ~2 km — comparable to the
  posterior 1σ (~1 km) and well above the frozen ±0.4 % libration_sys_frac
  residual.
- The eta = 1e12 whole-shell rows are not a physical shell (a 27 km shell
  that soft would not sustain the observed figure); listed as the limit.
- Qualitatively consistent with Bagheri et al. 2025 (p. 6, App. C: low
  sensitivity to shell viscosity; far more sensitive to thickness than
  rigidity). Quantitatively their ~10 % at eta ~1.25e14 is an ANDRADE
  number; Maxwell gives ~0 % beyond elastic there, Andrade here gives
  2–6 %. Bagheri also calls a uniform 1e14 shell unrealistic.
- k2 IS rheology-sensitive where libration is not: |k2| 0.012 (elastic
  3.3 GPa) → 0.039 (1 GPa) → 0.075 (Maxwell 1e13) at zb 26.7 km, with Im k2
  up to −0.07. Bagheri et al. predict |k2| 0.024–0.032, between these.

## Defects found (not fixed here; reported)

1. **Units in PlanetProfile's ALMA gravity path.**
   `Gravity.ComputeGravityObservations` calls
   `pyALMA3Updated.build_model(..., norm=False)` with the period in kyr
   (`ALMA_TIME_UNIT_S`) but the viscosity column in Pa·s
   (`Gravity.py:797-806`, `defineStructs.py:572-573` convert only VP, VS,
   GS). In periodic mode the Maxwell rigidity mu·s/(s + mu/eta) then mixes
   kyr⁻¹ with s⁻¹: every Maxwell/Andrade ice layer is ~3.2e10× too viscous,
   i.e. effectively elastic, so PP's ALMA k2 is the elastic k2 regardless of
   the configured ice viscosity. Demonstrated: at a production node, k2 is
   0.01344 for eta from 1e10 to 1e22 Pa·s as-is; converting eta to Pa·kyr
   gives the expected Maxwell curve (|k2| 0.0134 → 0.082 at 1e13). This
   affects any PP run that relies on ALMA k2 with a viscous ice layer.
   Reviewer CONFIRMED (norm=False with eta in Pa·kyr equals norm=True with
   eta in Pa·s to 1e-7; upstream alma uses the same t0 = 1 kyr). Andrade
   layers are affected too (anelastic term ∝ (s·eta/mu)^-α, suppressed
   ~7e-4 at α 0.3); defaults make Ih_conv/III_conv/V_conv/Clath_conv
   Andrade and II–VI Maxwell, all effectively elastic in this path.
   BLAST RADIUS (default backend, Params.Gravity.backend != 'tidalpy',
   Gravity.py:106-111): Planet.Gravity.h/l/k/delta/y/y1, libration_m (the
   deformable branch, via ALMA y1 and Re k2), kAmp/hAmp and phases, saved
   gravity files, Exploreogram/Monte Carlo kLoveAmp/hLoveAmp, Inversion
   kLoveAmp gridding (Inversion.py:309), and any heating from ALMA Im k2 —
   Im k2 is underestimated by orders of magnitude. NOT affected: the TidalPy
   backend and forward_models.forward_model_k2_flexible (eta in Pa·s, SI
   frequency).
   **Fix hazard:** do NOT divide every layer's eta. The Newtonian ocean
   (eta 0.0019 Pa·s) works only by accident as a small regularizing
   rigidity; converting it too gives k2 ≈ −0.85, and norm=True fails the
   same way (that was the "nonsense k2"). Convert only viscoelastic layers,
   or give fluid layers a reviewed floor, and show k2 converges in ocean eta.
   Before acceptance: k2 vs an independent analytic or TidalPy Maxwell
   result, an ocean-eta convergence scan, and a regression test that k2 is
   NOT invariant in eta.
2. **Rigid-branch K_int normalization** (known, pinned by
   `tests/librations_test.py`): deformable/rigid = 1.008–1.010 here, the
   wrong sign. Relative results above cancel it.

## Implication for an Enceladus k2 campaign

Adding k2 (user direction 2026-10-31) re-introduces the rheology
parameters B6 dropped, and needs the tidal solve. With k2 as an observable,
libration should move to the deformable branch so that libration and k2
share one rheology — which requires defect 2's reviewed repair first, and a
consistent viscosity unit convention (defect 1) if ALMA is used.
