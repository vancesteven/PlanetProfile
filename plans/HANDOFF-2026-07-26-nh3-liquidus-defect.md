# HANDOFF → Machine A: NH3 mu-liquidus solver defect (Phase-1 blocker)

**Date:** 2026-07-26
**From:** Machine B (Titan NH3 ocean campaign, Task #51 Phase B)
**Owner (per A/B split):** Machine A owns NH3-CoolProp thermodynamics
**Status of NH3 Phase B on B:** **HELD** pending this fix. MgSO4/NaCl continue on B in parallel.
**Severity:** CRITICAL — blocks the Titan NH3 Phase 1 2D cache AND corrupts any deployed
`cache_builder.build_single_structure` run that uses `comp='NH3'`.

---

## TL;DR

`PlanetProfile/Thermodynamics/NH3/NH3Props.py::muLiquidusCurve_K` (commit 9fb55fb7) returns a
**non-monotonic, partly Tb-independent liquidus**. The pressure-by-pressure melt-point depression
`dTliq(P) = Tm_pure(P) − Tliq(P)` — which should be smooth and mildly increasing with P — instead
has ~30% of pressures wrong: isolated **warm** points pinned at the L-K fallback value and isolated
**cold spikes of 15–22 K**. Downstream, `GetPfreeze` latches onto a spurious "liquid" cell at a
Tb-independent pressure, pinning the ice-Ih basal pressure and producing (a) flat `D_iceIh(Tb)` and
(b) physically impossible interleaved `liquid–iceV–liquid` structures.

This is **not** a scan-harness bug (see §5) and **not** the mu-equality *method* — it is the
*numerical solver* (bracketing + silent fallback) inside `muLiquidusCurve_K`.

---

## 1. Reproduction (cheap — no full PP build)

```bash
mamba run -n PPcl env PYTHONPATH=. KMP_DUPLICATE_LIB_OK=TRUE python -c "
import numpy as np, logging; logging.disable(logging.WARNING)
from PlanetProfile.Thermodynamics.NH3.NH3Props import muLiquidusCurve_K
P = np.arange(10.0, 220.0, 4.0)
Tliq, Tmpure = muLiquidusCurve_K(P, 30.0)
dep = Tmpure - Tliq
for p,d in zip(P, dep): print(f'{p:6.0f}  dep={d:6.2f} K')
"
```

**Observed at w=30 ppt** (should be a smooth ~3.2 K curve):

| symptom | pressures | value |
|---|---|---|
| baseline (correct) | most | ~3.2 K, smooth |
| spurious WARM (L-K fallback) | 70, 82, 106, 138, 178, 190, 194, 198, 202, 206, 214, 218 | pinned at **1.6315 K** = `NH3liquidusShift_K(30)` |
| spurious COLD spike | 94 | **22.0 K** |
| spurious COLD spike | 134 | **15.3 K** |

16 of 53 pressures deviate >0.5 K from their neighbours. The runtime log confirms the scale of the
failure: **"mu-based NH3 liquidus: L-K fallback at 135/298 pressures"** (~45% of pressures fail to
bracket and silently fall back to the Lee-Kesler polynomial).

The 1.6315 K value is exactly `NH3liquidusShift_K(30)` — i.e. these are the `except` fallback branch
(`NH3Props.py:201-203`), not real mu-equality solutions.

---

## 2. Root cause (in `muLiquidusCurve_K`, NH3Props.py:135-208)

1. **Silent bracketing failure → L-K fallback.** The sign-change scan runs on a **coarse 24-point
   T grid** (`Tscan = np.linspace(Tlo, min(T0+0.5,380), 24)`, line 165). At pressures where `_ln_aw`
   returns a noisy/NaN activity for some grid points, either no sign change is found (`raise
   RuntimeError('no bracket')`, line 185) or the *first* sign change is a spurious pair. Any exception
   drops to `Tliq[i] = T0 − NH3liquidusShift_K(w)` (line 202) — the warm 1.63 K points. There is **no
   log per-pressure**, only an aggregate count, so the non-monotonicity is invisible without probing.

2. **`_ln_aw` flash noise (NH3Props.py:119-132).** `_ln_aw` builds a fresh `AbstractState('HEOS',
   'Ammonia&Water')`, forces `iphase_liquid`, and reads `chemical_potential(1)` at each (P,T). At some
   (P,T) the forced-liquid flash of the estimated (linear simple-mixing-rule) binary returns a bad
   activity — that's the source of both the failed brackets and the 15–22 K cold spikes (a wrong
   sign-change bracket bisects to a wildly wrong root).

3. **No monotonicity/outlier guard.** Nothing enforces that `dTliq(P)` is smooth in P or monotone in
   w. A single bad `_ln_aw` at one pressure silently corrupts that pressure's root with no downstream
   check.

---

## 3. Downstream consequence (why this blocks Phase 1)

- `GetPfreeze` (`HydroEOS.py:1255-1266`) scans P upward for the first Ih→liquid sign change in the
  melt-EOS phase function. With phantom liquid cells sprinkled through the ice-Ih field, it latches
  onto the **first spurious cell (~56 MPa at w=30)**. That pressure is **Tb-independent**, so a full
  PP build gives `PbI_MPa = 56.41` identically at Tb=241/250/262 K → **flat D_iceIh** (42.2 km across
  the whole 241–263 K range at w=30; 48.5 km at w=50). The physically correct PbI at w=30/Tb=250
  should be ≈195 MPa (where Tliq(w=30) ≈ 250 K), a much thicker shell.
- Built structures acquire **multiple disconnected liquid layers**, e.g. v2 w=100/Tb=241:
  `Sil | 0 | V | 0 | V | 0 | V | 0 | Ih | Clath` (four "oceans" interleaved with ice V). Hydrostatically
  impossible at monotonically increasing P,T.
- **Deployed exposure:** `build_single_structure` uses the identical PP melt path (no PfreezeRes
  override → 0.05 default), so this is *not* probe-specific. Any Phase 1 `build_tbw_grid_cache` with
  `comp='NH3'` bakes in these wrong structures, wrong D_ocean, wrong basal P, wrong conductivity
  layering — and `region_phases`-keyed bilinear interpolation mis-fires because the bogus multi-liquid
  signatures differ node-to-node.

---

## 4. Required fix + validation (reviewer aad02da53, 2026-07-26)

**Fix `muLiquidusCurve_K` so the returned liquidus is physical:**
- Robust bracketing: refine the T grid where sign changes are ambiguous; do not accept the first
  sign-change if neighbouring pressures disagree by >~0.3 K.
- Guard `_ln_aw`: detect/repair CoolProp flash noise (metastable-liquid mis-flash); reject brackets
  built on non-finite or outlier activity points.
- Add a **monotonicity/outlier filter** on `dTliq(P)`: it must be smooth in P and increasing in w.
  Where a pressure fails, interpolate across P from good neighbours rather than silently dropping to
  L-K (or at minimum log per-pressure so the fallback is visible).
- Reduce reliance on the silent aggregate-only fallback; a 45% fallback rate should be a hard warning.

**Validation gates before Phase 1 can proceed (all on B or A, whoever runs first):**
1. `dTliq(P)` smooth in P, increasing in w; no isolated point >~0.3 K off neighbours. **Target table
   (reviewer-verified against spec):** depression ≈ **1.1 / 3 / 5 / 15 K at w = 10 / 30 / 50 / 150 ppt**.
2. Melt-EOS phase function along T=Tb is **monotone non-decreasing in phase index vs P** (no liquid
   cell above an ice-Ih cell at higher P until the true ocean).
3. 3-point Tb probe (241/250/262 K) at w=30: **PbI_MPa and D_iceIh vary monotonically with Tb**
   (expect a ~150–200 MPa PbI swing, not the current 0).
4. No built structure has >1 contiguous liquid layer across the proposed domain.
5. Tighten `PfreezeRes` to **≤0.5 MPa** for the production cache and re-verify <1% Pb quantization
   (D_iceIh maps steeply to Pb). NOTE: the current `PfreezeRes=2.0` "<1% error" claim is
   **unverifiable while the liquidus oscillates** — coarsening the grid changes *which* spurious cell
   GetPfreeze hits. Only meaningful post-fix.
6. Cross-check `probe_one` (sums on `Planet.phase`/`Planet.r_m` primary grid) vs
   `build_single_structure` (sums on padded Gravity/ALMA `model` grid) D_iceIh on one node — they can
   disagree for multi-liquid stacks (MINOR).

---

## 5. What is NOT the problem (already cleared on B)

- **Not the scan-harness caching.** B initially reused the melt EOS across the Tb axis (perf
  optimization); that was a *separate* real bug (poisoned by the Tb-centred melt-EOS T-window,
  `SetupInit.py:319/327`, so reuse extrapolates the phase interpolator out-of-range at other Tb). It
  was reverted to per-node `EOSlist.loaded.clear()` (matching deployed `cache_builder.py:122`). The
  corrected v2 scan (per-node clear) still shows the liquidus defect → it is upstream of the harness.
- **Not `PfreezeRes=2.0` alone.** Entangled with the defect (see gate 5) but not the root cause.
- **Not process-global state.** `_MIX_READY`/`_ensure_mixture` (NH3Props.py:51,211) is a
  comp/w-independent one-time binary-pair registration — safe to persist across `EOSlist.loaded.clear()`.
- **Not the mu-equality method.** The *mean* behaviour (dilute limit, monotone-in-w depression)
  matches the spec; the defect is isolated per-pressure numerical noise, not a systematic method bias.
  **Do not change the chemical-potential-equality formulation** — only the solver robustness.

---

## 6. Artifacts

- Contaminated (caching bug, v1): `/tmp/titan_nh3_phase0_results.json`,
  `/tmp/titan_nh3_phase0_v1_contaminated.log`
- Correct harness / defective liquidus (v2): `/tmp/titan_nh3_phase0_results_v2.json`,
  `/tmp/titan_nh3_phase0_v2.log`
- Caching-bug control: `/tmp/nh3_cache_control.py`
- Scan driver: `plans/scripts/titan_nh3_phase0_scan.py` (per-node clear, PfreezeRes=2.0)
- Reduced-stack probe: `plans/scripts/titan_tb_probe.py::probe_one`
- Campaign spec: `plans/active/titan-nh3-ocean-campaign-spec.md`
- Full reviewer report: Task #51 notes (agent aad02da53feef2ff7, 2026-07-26)

## 7. Post-fix Phase 1 domain (reviewer recommendation)

Once gates 1–4 pass: **Tb ∈ [241, 259] K × w ∈ [30, 100] ppt**. Exclude w=10 (fails Tb<253,
HP-heavy, sparse) and w=150 (fails Tb>256) from the first cache; span them only in a separate
campaign. Do NOT span w=10→150 in one bilinear cache — adjacent corners rarely share `region_phases`,
so bilinear would fall back to nearest-neighbour almost everywhere.
