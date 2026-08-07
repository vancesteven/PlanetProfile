"""Build a Titan ocean-bearing 2D (Tb x w) JOINT no-ocean+ocean structure cache
for the MgSO4 / NaCl free-gravity SBI campaign (Phase B of
plans/fluffy-snacking-fountain.md).

This is the standalone joint-build driver the Titan campaign was missing: the
NH3 joint cache was built by a since-deleted /tmp orchestrator, so MgSO4/NaCl
needed a committed, reviewer-traceable driver that clones
`cache_builder.build_tbw_grid_cache(...)` directly with the ratified settings
(acc725ea R1-R5 + the 2026-08-06 onset-probe/reviewer sign-off, §0.13 of
plans/MACHINE-B-HANDOFF.md).

Reviewer-ratified per-composition settings (2026-08-06, agent a3c9ed8ae24664527):

  NaCl (CLEARED):
    - Tb union-sweep [233, 272] K, fine ~1 K across the whole TILTED diagonal
      (the ocean stripe slides ~40 K down as salinity rises; every Tb is a
      boundary for some w, so uniform-fine is required -- NOT NH3's single band).
    - w log-spaced [1, 290] ppt, capped at log10(290)=2.4624 (R3), including the
      ~45-55 ppt nodes R4 requires (worst sub-step 0.377 sigma when present).
    - extrap_ocean=True: a no-op (== clamp) because NaCl Pmax is truncated at
      1000.1 MPa unconditionally + flat-extrapolating spline (reviewer Item-1).
    - Monotonicity sanity PASS (/tmp/nacl_monotonicity_check.json), retry-corner
      discriminator PASS (/tmp/nacl_corner_discriminator.json): the cold corner
      retries to a physical frozen no-ocean interior; the warm/high-w corner
      either builds a real ocean or the NO_OCEAN retry raises a distinct
      no-ice-lid ValueError -> stored None (excluded), never a phantom no-ocean.

  MgSO4 (extrap_ocean HELD pending reviewer verdict -- see --mgso4-extrap):
    - Genuine monotone onset ~[248, 255] K (nearly flat): compact fine band
      [248, 258] K. The w=194 low-Tb ocean "island" (ocean Tb=240, frozen
      242-249) is a PATHOLOGICAL Margules-lookup melting-monotonicity violation
      (reviewer 2026-08-06); EXCLUDE by construction -> NO Tb nodes below 248 K
      at w>=180 ppt. Post-build invariant asserts no cached node with w>=180 AND
      Tb<248 carries has_ocean=True (the retry canNOT catch these -- they build
      as phantom oceans, they do not raise NoIceLiquidTransitionError).
    - w log-spaced [1, 194] ppt (2-molal cap = log10(194)=2.2878), finer in the
      150-194 corner where boundary curvature drives ocean fraction.
    - MgSO4 phaseType='lookup' (Margules, uncapped Pmax=inf; R2). MgSO4 density
      table Pmax=800 MPa but deep w=194 ocean reaches ~1100 MPa, so extrap_ocean
      is NOT a no-op here (3.26% rho divergence at 1100 MPa):
      /tmp/mgso4_extrap_divergence.json. Setting deferred to the reviewer.

Env (mandatory): run under PPcl with thread pins; PP generation (numba) and any
torch training MUST be separate processes (libomp hazard).

    mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
      KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
      MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
      python plans/scripts/titanG_build_ocean_cache.py --comp NaCl \
        --out /tmp/titan_nacl_joint_structure_grid_2d.pkl

Use --validate for a coarse diagonal subset (fast geometry/timing check before
the full multi-hour build), then re-run without it for production.

Author: PlanetProfile team (Machine B), 2026-08-06.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

# Thread pins + numba cache before any PP/TidalPy import.
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/pp_numba_cache")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import logging
import warnings
for _n in ["PlanetProfile", "SeaFreeze", "tidalpy", "TidalPy", "lbftd"]:
    logging.getLogger(_n).setLevel(logging.ERROR)
    logging.getLogger(_n).propagate = False
warnings.filterwarnings("ignore")

import numpy as np

MODULE = "PlanetProfile.Default.Titan.PPTitan"
CUNCERTAINTY = 0.06


# ---------------------------------------------------------------------------
# Grid definitions (reviewer-ratified, 2026-08-06)
# ---------------------------------------------------------------------------
def nacl_grids(validate: bool):
    """NaCl: fine ~1 K Tb across the full [233,272] tilted diagonal; log w to 290."""
    if validate:
        # Coarse diagonal subset: geometry + per-node timing sanity only.
        Tb = [233.0, 238.0, 244.0, 250.0, 256.0, 262.0, 268.0, 272.0]
        w = [1.0, 10.0, 55.0, 100.0, 200.0, 290.0]
        return np.array(Tb), np.array(w)
    # Production: 1 K over [233, 272] -> 40 nodes (reviewer "~30-40").
    Tb = np.arange(233.0, 272.0 + 0.5, 1.0)
    # Log-spaced w [1, 290] with the R4-required 45/55 ppt nodes inserted.
    w_log = np.logspace(np.log10(1.0), np.log10(290.0), 13)
    w = np.unique(np.round(np.concatenate([w_log, [45.0, 55.0]]), 4))
    return Tb, w


def mgso4_grids(validate: bool):
    """MgSO4: compact fine band across the flat ~250 K monotone onset; NO nodes
    below 248 K at w>=180 (island exclusion). Log w to 194 (2-molal cap), finer
    in the 150-194 corner."""
    if validate:
        Tb = [248.0, 250.0, 252.0, 254.0, 258.0]
        w = [1.0, 10.0, 55.0, 100.0, 150.0, 194.0]
        return np.array(Tb), np.array(w)
    # Production: fine 0.5 K across [248, 254], coarser 1 K up to 258.
    Tb = np.unique(np.round(np.concatenate([
        np.arange(248.0, 254.0 + 0.25, 0.5),
        np.arange(255.0, 258.0 + 0.5, 1.0),
    ]), 3))
    # Log w [1, 194] with extra density in the 150-194 curvature corner.
    w_log = np.logspace(np.log10(1.0), np.log10(194.0), 11)
    w = np.unique(np.round(np.concatenate([w_log, [45.0, 55.0, 160.0, 175.0, 185.0]]), 4))
    return Tb, w


CONFIGS = {"NaCl": nacl_grids, "MgSO4": mgso4_grids}

# MgSO4 density-table pressure ceiling (EOS2_MgSO4_planetary_smaller) and the
# reviewer-set extrapolation ceiling (2026-08-06, agent a3c9ed8ae24664527).
MGSO4_TABLE_PMAX_MPA = 800.0
# Ceiling = 1400 MPa = production PPTitan deepest w=194 node (~1371 MPa) + ~29 MPa
# margin; linear extrap verified monotone/stiffening to 1500 MPa, projection
# error ~1.5% at the ceiling (pure-water-benchmarked). Reviewer raised this from
# the initial 1200 (which was sized against the coreless PPTest50 probe,
# P_basal<=1109, and would have gutted the entire w=194 salinity cap). Reject
# beyond 1400 to fence genuine runaway (a node at ~1800 MPa = interior problem).
MGSO4_EXTRAP_CEILING_MPA = 1400.0
MGSO4_STRATIFY_PMA = 1200.0  # >this = "deep" extrap band (~1-1.5% error); flag distinctly
# Physical dρ/dP band for the linear extrapolant above the table (reviewer
# condition #3): the in-table edge slope is ~0.157; extrap sits ~0.15-0.157 and
# gently stiffens. Reject any column whose above-800 dρ/dP leaves this band —
# catches a defective w-interpolation that plain positivity would miss.
MGSO4_DRHODP_BAND = (0.14, 0.17)  # kg/m^3/MPa
MGSO4_DEEP_T_MAX_K = 300.0  # reviewer condition #5: spot-flag ocean-base T above this


def _ocean_extrap_diagnostics(struct):
    """Return diagnostics for an MgSO4 ocean struct's deep (>table Pmax) liquid.

    Returns dict with:
      P_basal_MPa   : max P over liquid (phase 0) layers
      T_basal_K     : ocean-base temperature (max T over liquid layers)
      drhodp        : mean dρ/dP across liquid layers with P>800 MPa (or None)
      band_ok       : drhodp within MGSO4_DRHODP_BAND (True if <2 deep pts)
    """
    out = {"P_basal_MPa": None, "T_basal_K": None, "drhodp": None, "band_ok": True}
    try:
        P = np.asarray(struct["P_MPa"], float)
        rho = np.asarray(struct["rho"], float)
        T = np.asarray(struct["T_K"], float)
        ph = np.asarray(struct["phases"]).astype(int)
    except Exception:
        return out
    n = min(len(P), len(rho), len(T), len(ph))
    P, rho, T, ph = P[:n], rho[:n], T[:n], ph[:n]
    liq = ph == 0
    if not np.any(liq):
        return out
    out["P_basal_MPa"] = float(np.max(P[liq]))
    out["T_basal_K"] = float(np.max(T[liq]))
    mask = liq & (P > MGSO4_TABLE_PMAX_MPA)
    if np.count_nonzero(mask) >= 2:
        Pm = P[mask]; rm = rho[mask]
        order = np.argsort(Pm)
        Ps = Pm[order]; rs = rm[order]
        slope = float((rs[-1] - rs[0]) / (Ps[-1] - Ps[0])) if Ps[-1] > Ps[0] else None
        out["drhodp"] = slope
        if slope is not None:
            out["band_ok"] = bool(MGSO4_DRHODP_BAND[0] <= slope <= MGSO4_DRHODP_BAND[1])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--comp", required=True, choices=list(CONFIGS))
    ap.add_argument("--out", required=True, help="output pickle path (build in /tmp)")
    ap.add_argument("--validate", action="store_true",
                    help="coarse diagonal subset for geometry/timing (not production)")
    ap.add_argument("--mgso4-extrap", choices=["true", "false"], default=None,
                    help="MgSO4 extrap_ocean setting (HELD on reviewer; required for MgSO4)")
    ap.add_argument("--pfreeze-res", type=float, default=None,
                    help="override Planet.PfreezeRes_MPa (default: production fine default)")
    args = ap.parse_args()

    comp = args.comp
    Tb, w = CONFIGS[comp](args.validate)

    # extrap_ocean per composition (reviewer sign-off).
    if comp == "NaCl":
        extrap_ocean = True  # no-op == clamp (reviewer Item-1)
    else:  # MgSO4
        if args.mgso4_extrap is None:
            print("ERROR: MgSO4 extrap_ocean is HELD on scientific-reviewer adjudication "
                  "(table Pmax=800 MPa < ~1100 MPa basal P at w=194; NOT a no-op). "
                  "Re-run with --mgso4-extrap true|false once the verdict lands.",
                  file=sys.stderr)
            raise SystemExit(2)
        extrap_ocean = (args.mgso4_extrap == "true")

    ocean_overrides = {"comp": comp}
    if comp == "MgSO4":
        # R2: Margules (uncapped) path; PPTitan already uses 'lookup', set
        # explicitly for provenance clarity.
        ocean_overrides["phaseType"] = "lookup"

    print(f"=== Titan {comp} joint cache build "
          f"({'VALIDATION subset' if args.validate else 'PRODUCTION'}) ===", flush=True)
    print(f"  Tb ({len(Tb)}): {np.round(Tb, 2).tolist()}", flush=True)
    print(f"  w  ({len(w)}): {np.round(w, 3).tolist()}", flush=True)
    print(f"  extrap_ocean={extrap_ocean}  retry_frozen_as_no_ocean=True  "
          f"n_nodes={len(Tb) * len(w)}", flush=True)
    print(f"  out={args.out}", flush=True)

    from PlanetProfile.Inference import cache_builder as _cb
    # MANDATORY float coercion of PfreezeUpper_MPa (Machine B 2026-08-07 root-cause;
    # scientific-reviewer ab2d3748d0b00c78c PROCEED). PPTitan ships
    # Planet.PfreezeUpper_MPa = 230 as a Python INT. The melt-EOS grid
    # np.arange(PfreezeLower, PfreezeUpper, PfreezeRes) has Pmax=229.96 < 230, so
    # GetPfreeze's upper search endpoint P=230 always exceeds Pmax and is clamped
    # by ResetNearestExtrap (DataManip.py:15-35). Its size-1 branch (line 26,
    # outVar1=np.array(var1)) preserves the int64 dtype, so assigning the float
    # bound 229.96 into that int array TRUNCATES to 229 -> phase evaluates liquid
    # instead of ice III -> GetPfreeze routes to a failing GetZero bracket ->
    # NoIceLiquidTransitionError -> the onset node one step into the ocean region
    # is mislabeled frozen (MgSO4) / None (NaCl). Passing PfreezeUpper as a float
    # forces the clamp to 229.96 (ice III) -> the brute-force scan finds the real
    # Ih->liquid transition ~200 MPa -> a genuine ocean. VERIFIED to repair EXACTLY
    # the two Tb=252 K defect nodes (MgSO4 w=4.857 -> OCEAN D=15.60; NaCl w=4.127 ->
    # OCEAN D=31.06) with every other node byte-identical and the genuine frozen
    # onset preserved. This is a SYMPTOMATIC cache-build fix; the underlying
    # ResetNearestExtrap int-truncation bug (which affects ANY scalar-int P/T query
    # beyond the EOS domain) is referred to Machine A as a shared-thermodynamics
    # correction (dtype=float in the size-1 branch). See
    # validation_reports/titan_saltcaches/tb252_rootcause_2026_08_07.json.
    planet_overrides = {"PfreezeUpper_MPa": 230.0}
    # Optionally override the ice-Ih/liquid boundary search resolution.
    # Production uses the fine default (reviewer launch-condition 3); the probe's
    # 2 MPa is scan-only. build_tbw_grid_cache pulls PfreezeRes_MPa from the
    # template, so an override is only needed for a deliberate coarse test.
    if args.pfreeze_res is not None:
        print(f"  PfreezeRes_MPa override -> {args.pfreeze_res}", flush=True)
        planet_overrides["PfreezeRes_MPa"] = float(args.pfreeze_res)
    print(f"  planet_overrides -> {planet_overrides}  (PfreezeUpper float-coerced)",
          flush=True)

    t0 = time.time()
    cache = _cb.build_tbw_grid_cache(
        MODULE, list(Tb), list(w), args.out, progress=True,
        ocean_overrides=ocean_overrides,
        bulk_overrides={"Cuncertainty": CUNCERTAINTY},
        planet_overrides=planet_overrides,
        extrap_ocean=extrap_ocean,
        retry_frozen_as_no_ocean=True,
    )
    wall_min = (time.time() - t0) / 60.0

    # -------- has_ocean map + invariants --------
    structs = cache["structures"]
    n_w = len(w)
    n_none = sum(1 for s in structs if s is None)
    n_no_ocean = sum(1 for s in structs if s is not None and s.get("has_ocean") is False)
    n_ocean = sum(1 for s in structs if s is not None and s.get("has_ocean") is True)

    rows = []
    for it in range(len(Tb)):
        r = ""
        for iw in range(n_w):
            s = structs[it * n_w + iw]
            if s is None:
                r += "N"
            elif s.get("has_ocean") is True:
                r += "O"
            else:
                r += "F"
        rows.append(r)

    print(f"\n=== {comp} build done in {wall_min:.1f} wall-min ===", flush=True)
    print(f"  nodes={len(structs)} ocean={n_ocean} no_ocean={n_no_ocean} none={n_none}", flush=True)
    print("  has_ocean map (rows=Tb asc, cols=w asc; O=ocean F=frozen-no-ocean N=None):", flush=True)
    for it, r in enumerate(rows):
        print(f"    Tb={Tb[it]:6.1f}  {r}", flush=True)

    # --- Sandwich invariant (reviewer ab2d3748d0b00c78c, 2026-08-07) ---
    # No None or frozen (has_ocean=False) node may lie between two ocean nodes in
    # a single w-column (Tb ascending). This catches the Tb=252 K int/float
    # ResetNearestExtrap truncation artifact and any future endpoint-clamp
    # regression: a genuine freeze onset is a single monotone F/N -> O transition
    # per column, so an ocean ABOVE and BELOW a non-ocean node is unphysical.
    sandwich_violations = []
    for iw in range(n_w):
        col_cls = []  # (iT, 'O'/'F'/'N')
        for it in range(len(Tb)):
            s = structs[it * n_w + iw]
            c = "N" if s is None else ("O" if s.get("has_ocean") is True else "F")
            col_cls.append(c)
        ocean_iTs = [i for i, c in enumerate(col_cls) if c == "O"]
        if len(ocean_iTs) >= 2:
            lo, hi = min(ocean_iTs), max(ocean_iTs)
            for it in range(lo + 1, hi):
                if col_cls[it] != "O":
                    sandwich_violations.append(
                        (float(Tb[it]), float(w[iw]), col_cls[it]))
    if sandwich_violations:
        print(f"  *** SANDWICH INVARIANT VIOLATED at {sandwich_violations} "
              f"(non-ocean node between two ocean nodes in a w-column — "
              f"int/float clamp artifact or onset regression) ***", flush=True)
    else:
        print("  sandwich invariant OK (no non-ocean node between two ocean "
              "nodes in any w-column).", flush=True)

    # MgSO4 island invariant: no cached node with w>=180 AND Tb<248 may be ocean.
    invariant_violations = []
    eos_extrap_nodes = []       # (Tb,w,P_basal) for 800<P<=1200 -> flag (mild)
    eos_extrap_deep_nodes = []  # (Tb,w,P_basal) for 1200<P<=1400 -> flag (deep ~1-1.5%)
    extrap_ceiling_rejects = []  # (Tb,w,P_basal) for P>1400 -> reject (store None)
    band_violations = []        # (Tb,w,drhodp) outside physical dρ/dP band
    hot_base_nodes = []         # (Tb,w,T_basal) with ocean-base T > 300 K (spot-flag)
    if comp == "MgSO4":
        for it in range(len(Tb)):
            for iw in range(n_w):
                flat = it * n_w + iw
                s = structs[flat]
                # (1) island exclusion invariant
                if w[iw] >= 180.0 and Tb[it] < 248.0:
                    if s is not None and s.get("has_ocean") is True:
                        invariant_violations.append((float(Tb[it]), float(w[iw])))
                if s is None or s.get("has_ocean") is not True:
                    continue
                d = _ocean_extrap_diagnostics(s)
                P_basal = d["P_basal_MPa"]
                if P_basal is None:
                    continue
                s["P_basal_ocean_MPa"] = round(P_basal, 2)
                s["T_basal_ocean_K"] = round(d["T_basal_K"], 2) if d["T_basal_K"] else None
                if d["drhodp"] is not None:
                    s["ocean_drhodp_above800"] = round(d["drhodp"], 5)
                # (reviewer #5) ocean-base T spot-flag
                if d["T_basal_K"] is not None and d["T_basal_K"] > MGSO4_DEEP_T_MAX_K:
                    hot_base_nodes.append((float(Tb[it]), float(w[iw]), round(d["T_basal_K"], 1)))
                # (reviewer #2 tightened) ceiling reject beyond 1400
                if P_basal > MGSO4_EXTRAP_CEILING_MPA:
                    extrap_ceiling_rejects.append((float(Tb[it]), float(w[iw]), round(P_basal, 1)))
                    structs[flat] = None
                    continue
                # (reviewer #2) stratified eos_extrapolated flag
                if P_basal > MGSO4_STRATIFY_PMA:
                    s["eos_extrapolated"] = "deep"   # >1200 MPa, ~1-1.5% error
                    eos_extrap_deep_nodes.append((float(Tb[it]), float(w[iw]), round(P_basal, 1)))
                elif P_basal > MGSO4_TABLE_PMAX_MPA:
                    s["eos_extrapolated"] = "mild"   # 800-1200 MPa, <0.5% error
                    eos_extrap_nodes.append((float(Tb[it]), float(w[iw]), round(P_basal, 1)))
                else:
                    s["eos_extrapolated"] = False
                # (reviewer #3 strengthened) dρ/dP physical-band check
                if d["drhodp"] is not None and not d["band_ok"]:
                    band_violations.append((float(Tb[it]), float(w[iw]), round(d["drhodp"], 5)))

        if invariant_violations:
            print(f"  *** ISLAND INVARIANT VIOLATED at {invariant_violations} "
                  f"(phantom ocean below 248 K at w>=180) ***", flush=True)
        else:
            print("  island invariant OK (no has_ocean=True node with w>=180 & Tb<248).",
                  flush=True)
        print(f"  EOS-extrapolated mild (800<P<=1200): {len(eos_extrap_nodes)}  |  "
              f"deep (1200<P<=1400): {len(eos_extrap_deep_nodes)}", flush=True)
        if extrap_ceiling_rejects:
            print(f"  extrap-ceiling REJECTS (P_basal>1400 MPa, set None): "
                  f"{extrap_ceiling_rejects}", flush=True)
        else:
            print("  extrap-ceiling OK (no ocean node exceeds 1400 MPa basal P).", flush=True)
        if band_violations:
            print(f"  *** dρ/dP OUT OF BAND {MGSO4_DRHODP_BAND} at {band_violations} "
                  f"(defective extrap column — investigate) ***", flush=True)
        else:
            print(f"  dρ/dP within physical band {MGSO4_DRHODP_BAND} in all extrap columns.",
                  flush=True)
        if hot_base_nodes:
            print(f"  NOTE ocean-base T>300 K at {hot_base_nodes} (reviewer spot-check: "
                  f"extrap verified to 285 K; confirm these).", flush=True)
        else:
            print("  ocean-base T <= 300 K in all ocean nodes (extrap T-range OK).", flush=True)
        # Recompute counts after ceiling rejects
        n_none = sum(1 for s in structs if s is None)
        n_no_ocean = sum(1 for s in structs if s is not None and s.get("has_ocean") is False)
        n_ocean = sum(1 for s in structs if s is not None and s.get("has_ocean") is True)

    # MgSO4 guards mutate structs in place (flag + ceiling-reject); the builder
    # already wrote the pickle, so re-persist the guarded cache. Recompute the
    # has_ocean map to reflect any ceiling rejects.
    if comp == "MgSO4":
        import pickle
        rows = []
        for it in range(len(Tb)):
            r = ""
            for iw in range(n_w):
                s = structs[it * n_w + iw]
                r += "N" if s is None else ("O" if s.get("has_ocean") is True else "F")
            rows.append(r)
        tmp = args.out + ".tmp"
        with open(tmp, "wb") as f:
            pickle.dump(cache, f)
        os.replace(tmp, args.out)
        print("  re-persisted guarded MgSO4 cache (eos_extrapolated flags + "
              "ceiling rejects applied).", flush=True)

    manifest = {
        "comp": comp, "validate": args.validate, "out": args.out,
        "Tb_grid": [float(x) for x in Tb], "w_grid_ppt": [float(x) for x in w],
        "n_nodes": len(structs), "n_ocean": n_ocean, "n_no_ocean": n_no_ocean,
        "n_none": n_none, "extrap_ocean": extrap_ocean,
        "has_ocean_map_Tb_rows": rows, "build_wall_min": round(wall_min, 1),
        "pfreeze_upper_float_coerced": True,
        "pfreeze_upper_MPa": 230.0,
        "sandwich_invariant_violations": sandwich_violations,
        "tb252_int_float_fix_provenance": (
            "PfreezeUpper_MPa float-coerced (was int 230 in PPTitan template) to "
            "avoid ResetNearestExtrap size-1 int-dtype truncation (229.96->229) that "
            "mislabeled the Tb=252 K onset node frozen/None. Machine B root-cause + "
            "scientific-reviewer ab2d3748d0b00c78c PROCEED, 2026-08-07. See "
            "validation_reports/titan_saltcaches/tb252_rootcause_2026_08_07.json."
        ),
        "island_invariant_violations": invariant_violations if comp == "MgSO4" else None,
        "eos_extrapolated_mild_nodes": eos_extrap_nodes if comp == "MgSO4" else None,
        "eos_extrapolated_deep_nodes": eos_extrap_deep_nodes if comp == "MgSO4" else None,
        "extrap_ceiling_rejects": extrap_ceiling_rejects if comp == "MgSO4" else None,
        "drhodp_band_violations": band_violations if comp == "MgSO4" else None,
        "hot_ocean_base_nodes": hot_base_nodes if comp == "MgSO4" else None,
        "mgso4_table_Pmax_MPa": MGSO4_TABLE_PMAX_MPA if comp == "MgSO4" else None,
        "mgso4_extrap_ceiling_MPa": MGSO4_EXTRAP_CEILING_MPA if comp == "MgSO4" else None,
        "mgso4_drhodp_band": list(MGSO4_DRHODP_BAND) if comp == "MgSO4" else None,
        "mgso4_ceiling_provenance": (
            "1400 MPa = production PPTitan deepest w=194 node (~1371) + margin; "
            "linear extrap verified monotone/stiffening to 1500 MPa, projection "
            "error ~1.5% at ceiling (water-benchmarked); reviewer a3c9ed8ae24664527 "
            "2026-08-06 raised from 1200 (probe-sized) to preserve the w=194 cap."
        ) if comp == "MgSO4" else None,
    }
    mpath = args.out.replace(".pkl", "_manifest.json")
    with open(mpath, "w") as f:
        json.dump(manifest, f, indent=1)
    print(f"\nManifest -> {mpath}", flush=True)


if __name__ == "__main__":
    main()
