"""
Probe Tb_K window for Titan ocean-bearing SBI artifacts (NaCl and NH3).

Empirically determines the Tb_K range over which a Titan ocean EXISTS at both
low and high salinity ends of each composition axis. Runs a coarse grid of
single PP model builds — NOT a production cache build.

Hard EOS limits (from HydroEOS.py):
  NaCl: T[229, 501] K, Pmax 1000.1 MPa, wMax 290.3 ppt
  NH3:  T[241, 399.2] K, Pmax 2228.4 MPa, wMax 290.1 ppt
  NH3 Tb floor: >= 241 K (hard EOS lower bound)
  NaCl wMax floor: <= 290 ppt (log10 2.462)

Template: PlanetProfile.Test.PPTest50 (no-ocean Titan analog with
CONSTANT_INNER_DENSITY mechanism; we switch it to ocean mode by clearing
NO_OCEAN_EXCEPT_INNER_ICES, setting ConstantProps['Inner']=True, and
overriding comp/w/Tb).

The probe calls PlanetProfile directly (like build_single_structure does
internally) so we get full PP output including D_ocean_km, D_iceIh_km, etc.

Usage (run from repo root):
    mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
        KMP_DUPLICATE_LIB_OK=TRUE \
        python plans/scripts/titan_tb_probe.py

Author: PlanetProfile team, 2026-07-25
"""
from __future__ import annotations

import logging
import os
import sys
import time
import traceback
from copy import deepcopy
from typing import Any, Dict, List, Optional

# Ensure NUMBA cache dir is set before any TidalPy import
if not os.environ.get("NUMBA_CACHE_DIR"):
    _nd = "/tmp/pp_numba_cache"
    os.makedirs(_nd, exist_ok=True)
    os.environ["NUMBA_CACHE_DIR"] = _nd
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np

logging.basicConfig(
    level=logging.WARNING,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
# Quieten PP's very verbose EOS loading output
for _noisy in [
    "PlanetProfile",
    "SeaFreeze",
    "tidalpy",
    "TidalPy",
    "lbftd",
]:
    logging.getLogger(_noisy).setLevel(logging.ERROR)

log = logging.getLogger("titan_tb_probe")

# ---------------------------------------------------------------------------
# Grid definition
# ---------------------------------------------------------------------------
# Template: PPTest50 (no-ocean Titan analog; Sil.rhoSilWithCore_kgm3 set to match
# a realistic silicate density for an ocean-bearing Titan).
TEMPLATE = "PlanetProfile.Test.PPTest50"

# Silicate density for ocean-bearing Titan models:
# From the working NaCl-100 profile (rhoSilWithCore_kgm3=3539 in PPTitanNaCl100ppt)
RHO_SIL_KGM3 = 3539.0

# Probe grid — designed to BRACKET the ocean-existing band:
#
# NaCl melting point at Titan ice-ocean pressures (~50-500 MPa):
#   w=1 ppt (dilute):  Tf ~ 270-273 K (nearly pure water)
#   w=100 ppt:         Tf ~ 243-246 K (known working range)
#   w=290 ppt (brine): Tf ~ 220-235 K (very depressed)
#
# Strategy: probe two separate Tb grids — one targeting dilute and one
# targeting concentrated ends. The production 2D grid spans this wedge;
# a single Tb axis is shared, so the grid covers only the OVERLAP of
# ocean-existing regions across the full salinity range.
#
# NH3 melting point at these conditions:
#   w=1 ppt:  Tf ~ 271-273 K
#   w=70 ppt: Tf ~ ~250-260 K (estimated; depends on P)
#   EOS lower bound: 241 K (hard floor for NH3)
#
# For Phase A: we probe a REPRESENTATIVE set to find the shared Tb band.
# Per-salinity probes at key w values reveal the Tb window per composition.

# Shared Tb grid (wide sweep); NH3 floor >= 241 K
# For NaCl dilute (w=1): expect ocean at Tb close to 271 K
# For NaCl brine (w=290): expect ocean at Tb < ~241 K (where ice cap exists)
# For NaCl mid (w=100): expect ocean at Tb ~ 243-246 K

# NaCl: cover full range from brine-friendly to dilute-friendly.
# Reviewer (2026-08-06) computed the ocean-onset diagonal: onset ~251 K at
# w=1, dropping to ~220 K at w=290 (a ~31 K sweep). NaCl EOS T-floor is 229 K,
# so the low-Tb end of that diagonal (w>~150, onset<233) is UNREACHABLE — the
# achievable boundary crossing sits in ~[233, 252] K. Probe densely there and
# add a couple of near-floor nodes to expose the truncation empirically.
TB_GRID_NACL = [233.0, 235.0, 238.0, 241.0, 244.0, 247.0, 249.0, 250.0, 251.0,
                252.0, 255.0, 260.0, 265.0, 270.0, 272.0]

# MgSO4: reviewer onset ~251.5 K (w=1) -> ~249.5 K (w=100-150) -> ~240 K
# (w=194); a ~11 K sweep that is nearly flat to ~150 ppt then bends down at the
# 194 ppt cap. Probe densely over ~[238, 255] K.
TB_GRID_MGSO4 = [238.0, 240.0, 242.0, 244.0, 246.0, 248.0, 249.0, 250.0, 251.0,
                 252.0, 253.0, 255.0]

# NH3: cover range respecting 241 K floor (already deployed; kept for parity /
# re-runs, not run by default).
TB_GRID_NH3  = [241.0, 243.0, 246.0, 249.0, 252.0, 255.0, 258.0, 260.0, 262.0, 265.0, 268.0]

# Key salinity probe points — chosen to sample the extremes of each axis
# NaCl: test w=1 (dilute), w=100 (mid/prior-tested), w=290 (max brine)
W_NACL = [1.0, 100.0, 290.0]
# MgSO4: w=1 (dilute), w=100 (mid), w=194 (2 molal cap)
W_MGSO4 = [1.0, 100.0, 194.0]
# NH3: test w=1 (dilute), w=35 (mid), w=70 (max)
W_NH3  = [1.0, 35.0, 70.0]

# Scan-resolution ice-Ih/liquid boundary search: 2.0 MPa is reviewer-endorsed
# for LOCATING the boundary (<1% Pb error) and keeps NaCl/MgSO4 nodes fast. The
# production CACHE build must use the fine default (reviewer launch-condition 3).
PROBE_PFREEZE_RES_MPA = 2.0
PROBE_PFREEZE_UPPER_MPA = 600.0

# Widen Cuncertainty to permit off-default Tb (same as Phase-C1 ocean builds)
CUNCERTAINTY = 0.060


# ---------------------------------------------------------------------------
# Single probe call — direct PP invocation
# ---------------------------------------------------------------------------
def probe_one(
    Tb_K: float,
    comp: str,
    w_ppt: float,
    clear_eos: bool = True,
    PfreezeRes_MPa: Optional[float] = None,
    PfreezeUpper_MPa: Optional[float] = None,
) -> Dict[str, Any]:
    """Run PP once for (Tb_K, comp, w_ppt) and return a summary dict.

    clear_eos: when True (default), wipe EOSlist.loaded before building.
        Set False to REUSE cached EOS across Tb nodes that share (comp, w):
        Tb does not change the ocean/melt EOS, only the integration start
        condition, so clearing per-Tb needlessly rebuilds the expensive
        NH3 mu-liquidus meltEOS every node.
    PfreezeRes_MPa: when set, overrides the ice-ocean-boundary search
        pressure resolution (default 0.05 MPa → 12,000 meltEOS P-pts).
        For a feasibility scan, 2.0 MPa (300 P-pts) cuts the NH3 meltEOS
        build from ~68 min to ~100 s with <1% Pb error.
    """
    import importlib

    from PlanetProfile.GetConfig import Params as configParams
    from PlanetProfile.Gravity.Gravity import SetupGravity
    from PlanetProfile.Main import PlanetProfile as RunPP
    from PlanetProfile.Utilities.defineStructs import Constants, EOSlist

    if clear_eos:
        EOSlist.loaded.clear()

    # Reload and deepcopy template
    if TEMPLATE in sys.modules:
        importlib.reload(sys.modules[TEMPLATE])
    else:
        importlib.import_module(TEMPLATE)
    mod = sys.modules[TEMPLATE]
    Planet = deepcopy(mod.Planet)

    # --- Switch from no-ocean to ocean-bearing mode ---
    Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES = False
    Planet.Do.NO_OCEAN = False

    # ConstantProps['Inner'] = True + fixed silicate density: this is the
    # correct mechanism (equivalent to AssignPlanetVal('rhoSilInput_kgm3', ...)).
    Planet.Do.ConstantProps['Inner'] = True
    Planet.Sil.rhoSilWithCore_kgm3 = RHO_SIL_KGM3
    Planet.Do.POROUS_ROCK = False
    Planet.Do.Fe_CORE = False

    # Disable Yao2014 spherical convection (that model is only valid for the
    # thick no-ocean ice shell; with an ocean, the ice is thin and the
    # ConvectionYao2014 call raises TypeError on the GetTfreeze API).
    Planet.Do.SPHERICAL_CONVECTION = False
    Planet.Do.ARRHENIUS_VISCOSITY = False
    Planet.Do.KALOUSOVA_CONVECTION = False
    Planet.Do.NO_ICE_CONVECTION = True

    # Hydrosphere settings matching PPTitanNaCl100ppt.
    # PfreezeUpper_MPa controls the upper search bound for the Ih-liquid
    # phase transition. PPTitanNaCl100ppt sets it to 300 MPa; for the probe
    # we use 600 MPa to also cover high-salinity brines where the Ih-L
    # boundary occurs at higher pressures (NaCl 290ppt, NH3 high-w).
    Planet.PfreezeUpper_MPa = 600.0
    if PfreezeUpper_MPa is not None:
        Planet.PfreezeUpper_MPa = float(PfreezeUpper_MPa)
    # meltEOS grid resolution for the ice-ocean boundary search. The NH3
    # mu-liquidus costs ~340 ms/pressure, so the 0.05 MPa default (12,000
    # P-pts → ~68 min) is untenable for a scan. 2.0 MPa (300 P-pts, ~100 s)
    # keeps Pb error <1% of a 50-200 MPa Titan ocean base.
    if PfreezeRes_MPa is not None:
        Planet.PfreezeRes_MPa = float(PfreezeRes_MPa)
    Planet.Ocean.PHydroMax_MPa = 1800.0
    Planet.Ocean.THydroMax_K = 350.0
    Planet.Ocean.deltaP = 8.0
    Planet.Ocean.deltaT = 0.5
    Planet.Ocean.phaseType = 'lookup'

    # Apply overrides
    Planet.Bulk.Tb_K = Tb_K
    Planet.Bulk.Cuncertainty = CUNCERTAINTY
    Planet.Ocean.comp = comp
    Planet.Ocean.wOcean_ppt = w_ppt

    # PP config
    configParams.Gravity.backend = "tidalpy"
    configParams.CALC_NEW = True
    configParams.CALC_NEW_GRAVITY = False
    configParams.NO_SAVEFILE = True
    configParams.SKIP_PLOTS = True

    t0 = time.time()
    try:
        Planet, Params = RunPP(Planet, configParams)
        Params.CALC_NEW_GRAVITY = True
        Planet, Params = SetupGravity(Planet, Params)
        elapsed = time.time() - t0

        # Compute D_ocean_km directly from liquid layers (phase=0)
        # Use Planet data structures to extract layer thicknesses
        # Planet.r_m and Planet.phase_m / Planet.zTop_km etc.
        # Best approach: iterate over phases using PlanetProfile's layer info

        D_ocean_km = 0.0
        D_iceIh_km = 0.0
        D_iceIII_km = 0.0
        D_iceV_km = 0.0
        D_iceVI_km = 0.0
        D_clath_km = 0.0
        P_ocean_basal_MPa = 0.0

        # Use Planet.phase and Planet.r_m arrays (same as cache_builder)
        try:
            from PlanetProfile.Utilities.defineStructs import Constants as C
            phases_raw = np.asarray(Planet.phase)
            r_raw = np.asarray(Planet.r_m)
            P_raw = np.asarray(Planet.P_MPa)
            n = min(len(phases_raw), len(r_raw), len(P_raw))
            phases_raw = phases_raw[:n]
            r_raw = r_raw[:n]
            P_raw = P_raw[:n]

            # Group consecutive same-phase runs
            if n > 1:
                boundaries = np.where(np.diff(phases_raw.astype(int)) != 0)[0] + 1
                boundaries = np.concatenate(([0], boundaries, [n]))
                for i in range(len(boundaries) - 1):
                    s, e = boundaries[i], boundaries[i+1]
                    ph = int(phases_raw[s])
                    thick = abs(r_raw[e-1] - r_raw[s]) / 1e3
                    if ph == 0:
                        D_ocean_km += thick
                        if len(P_raw[s:e]) > 0:
                            P_ocean_basal_MPa = max(P_ocean_basal_MPa, float(np.max(P_raw[s:e])))
                    elif ph == 1:
                        D_iceIh_km += thick
                    elif ph == 3:
                        D_iceIII_km += thick
                    elif ph == 5:
                        D_iceV_km += thick
                    elif ph == 6:
                        D_iceVI_km += thick
                    elif C.phaseClath <= ph < C.phaseClath + 10:
                        D_clath_km += thick
        except Exception as exc:
            log.warning(f"Phase extraction failed: {exc}")

        # Fallback: use Planet.D_H2O_km - Planet.zIce_m type attributes if above fails
        # (belt-and-suspenders; the phase loop above is the canonical approach)

        # Get CMR2
        cmr2 = np.nan
        try:
            cmr2 = float(Planet.CMR2mean)
        except (AttributeError, TypeError):
            try:
                cmr2 = float(Planet.Bulk.Cmeasured)
            except Exception:
                pass

        # Region phases from Reduced
        region_phases = []
        try:
            from PlanetProfile.Utilities.Indexing import PhaseConv
            ci = np.max(Planet.Reduced.changeIndices) - np.flipud(Planet.Reduced.changeIndices)
            n_layers = len(ci) - 1
            ph_arr = None
            # Get phases from Gravity model if available
            try:
                cols = Planet.Gravity.columns
                model = Planet.Gravity.ALMAModel["model"]
                pIdx = cols.index("phase")
                ph_arr = model[:, pIdx]
            except Exception:
                pass
            if ph_arr is not None:
                iConv = np.flipud(Planet.Reduced.iConv)
                for i in range(n_layers):
                    ph = int(ph_arr[ci[i]])
                    if C.phaseClath <= ph < C.phaseClath + 10:
                        ph = C.phaseClath
                    pstr = PhaseConv(ph, liq="0")
                    if iConv[ci[i]]:
                        pstr += "_conv"
                    region_phases.append(pstr)
        except Exception:
            pass

        return {
            "comp": comp,
            "Tb_K": Tb_K,
            "w_ppt": w_ppt,
            "status": "built",
            "D_ocean_km": round(D_ocean_km, 2),
            "D_iceIh_km": round(D_iceIh_km, 2),
            "D_iceIII_km": round(D_iceIII_km, 2),
            "D_iceV_km": round(D_iceV_km, 2),
            "D_iceVI_km": round(D_iceVI_km, 2),
            "D_clath_km": round(D_clath_km, 2),
            "P_ocean_basal_MPa": round(P_ocean_basal_MPa, 1),
            # Ice-Ih basal pressure (PbI) is the gate-3 diagnostic: under the
            # liquidus defect it was pinned Tb-independently (~56 MPa at w=30);
            # a healthy liquidus makes PbI swing ~150-200 MPa across Tb.
            "PbI_MPa": (round(float(Planet.PbI_MPa), 2)
                        if getattr(Planet, "PbI_MPa", None) is not None else None),
            "Pb_MPa": (round(float(Planet.Pb_MPa), 2)
                       if getattr(Planet, "Pb_MPa", None) is not None else None),
            "zb_km": (round(float(Planet.zb_km), 2)
                      if getattr(Planet, "zb_km", None) is not None else None),
            "region_phases": region_phases,
            "CMR2": float(cmr2) if np.isfinite(cmr2) else None,
            "elapsed_s": round(elapsed, 1),
        }
    except Exception as exc:
        elapsed = time.time() - t0
        tb_str = traceback.format_exc()
        return {
            "comp": comp,
            "Tb_K": Tb_K,
            "w_ppt": w_ppt,
            "status": f"FAILED: {type(exc).__name__}: {exc}",
            "traceback": tb_str[-400:],
            "D_ocean_km": None,
            "D_iceIh_km": None,
            "D_iceIII_km": None,
            "D_iceV_km": None,
            "D_iceVI_km": None,
            "D_clath_km": None,
            "P_ocean_basal_MPa": None,
            "region_phases": None,
            "CMR2": None,
            "elapsed_s": round(elapsed, 1),
        }


# ---------------------------------------------------------------------------
# Coverage table printer
# ---------------------------------------------------------------------------
def print_coverage_table(rows: List[Dict[str, Any]], comp: str) -> None:
    w_vals = sorted(set(r["w_ppt"] for r in rows))
    tb_vals = sorted(set(r["Tb_K"] for r in rows))

    print(f"\n{'='*70}")
    print(f"Coverage table: {comp}")
    print(f"  D_ocean_km (0.0 = frozen/no-ocean; FAILED = exception)")
    print(f"  (*) = HP ice present (D_iceIII+D_iceV+D_iceVI > 0.1 km)")
    w_header = "  ".join(f"w={w:6.1f}" for w in w_vals)
    print(f"{'Tb_K':>8}  {w_header}  notes")
    print("-" * 70)

    for Tb in tb_vals:
        row_strs = []
        notes = []
        both_ocean = True
        for w in w_vals:
            match = [r for r in rows if r["Tb_K"] == Tb and r["w_ppt"] == w]
            if not match:
                row_strs.append(f"{'---':>13}")
                both_ocean = False
                continue
            r = match[0]
            if r["status"] != "built" or r["D_ocean_km"] is None:
                row_strs.append(f"{'FAILED':>13}")
                both_ocean = False
            else:
                d = r["D_ocean_km"]
                hp = (r["D_iceIII_km"] or 0) + (r["D_iceV_km"] or 0) + (r["D_iceVI_km"] or 0)
                tag = "(*)" if hp > 0.1 else "   "
                row_strs.append(f"{d:8.1f}{tag}")
                if d <= 0.0:
                    both_ocean = False
                if hp > 0.1:
                    notes.append(f"w={w}: HP_ice {hp:.1f}km")
        flag = " <-- OCEAN both ends" if both_ocean else ""
        row_line = "  ".join(row_strs)
        note_str = "; ".join(notes)
        print(f"{Tb:8.1f}  {row_line}{flag}  {note_str}")

    print()


def print_detailed(rows: List[Dict[str, Any]], comp: str) -> None:
    print(f"\n{'='*70}")
    print(f"Detailed results: {comp}")
    header = (f"{'Tb_K':>7} {'w_ppt':>6} {'D_ocn':>7} {'D_Ih':>7} "
              f"{'D_III':>6} {'D_V':>6} {'D_VI':>6} {'Pcl_bot':>9} "
              f"{'CMR2':>7} {'t_s':>5}  region_phases")
    print(header)
    print("-" * 100)
    for r in sorted(rows, key=lambda x: (x["Tb_K"], x["w_ppt"])):
        if r["status"] != "built":
            short = r["status"][:80]
            print(f"{r['Tb_K']:7.1f} {r['w_ppt']:6.1f}  FAILED: {short}")
            continue
        rph = "|".join(r["region_phases"]) if r["region_phases"] else "?"
        cmr2 = f"{r['CMR2']:.4f}" if r["CMR2"] is not None else "  NaN "
        print(
            f"{r['Tb_K']:7.1f} {r['w_ppt']:6.1f} "
            f"{r['D_ocean_km']:7.1f} {r['D_iceIh_km']:7.1f} "
            f"{r['D_iceIII_km']:6.1f} {r['D_iceV_km']:6.1f} {r['D_iceVI_km']:6.1f} "
            f"{r['P_ocean_basal_MPa']:9.1f} "
            f"{cmr2:>7} {r['elapsed_s']:5.1f}  {rph}"
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    import json

    all_rows: List[Dict[str, Any]] = []

    # Default run targets the two compositions that still need onset tables
    # (NaCl + MgSO4). NH3 is already deployed; include it only via --nh3.
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--comps", nargs="+", default=["NaCl", "MgSO4"],
                    help="subset of NaCl MgSO4 NH3 to probe")
    args = ap.parse_args()

    ALL_CONFIGS = {
        "NaCl":  (TB_GRID_NACL,  W_NACL),
        "MgSO4": (TB_GRID_MGSO4, W_MGSO4),
        "NH3":   (TB_GRID_NH3,   W_NH3),
    }
    configs = [(c, ALL_CONFIGS[c][0], ALL_CONFIGS[c][1])
               for c in args.comps if c in ALL_CONFIGS]

    for comp, tb_grid, w_list in configs:
        print(f"\n{'#'*70}")
        print(f"# Probing {comp}  ({len(tb_grid)} Tb values x {len(w_list)} w values = "
              f"{len(tb_grid)*len(w_list)} builds)")
        print(f"{'#'*70}")
        sys.stdout.flush()
        for Tb in tb_grid:
            for w in w_list:
                print(f"  -> {comp} Tb={Tb:.1f} K  w={w:.1f} ppt ... ", end="", flush=True)
                r = probe_one(Tb, comp, w,
                              PfreezeRes_MPa=PROBE_PFREEZE_RES_MPA,
                              PfreezeUpper_MPa=PROBE_PFREEZE_UPPER_MPA)
                if r["status"] == "built":
                    hp = (r["D_iceIII_km"] or 0) + (r["D_iceV_km"] or 0) + (r["D_iceVI_km"] or 0)
                    print(
                        f"D_ocean={r['D_ocean_km']:.1f} km  "
                        f"D_Ih={r['D_iceIh_km']:.1f} km  "
                        f"HP={hp:.1f} km  "
                        f"t={r['elapsed_s']:.0f}s"
                    )
                else:
                    print(f"FAILED ({r['status'][:60]})  t={r['elapsed_s']:.0f}s")
                sys.stdout.flush()
                all_rows.append(r)

    # Summary tables
    for comp, _, _ in configs:
        comp_rows = [r for r in all_rows if r["comp"] == comp]
        print_coverage_table(comp_rows, comp)
        print_detailed(comp_rows, comp)

    # Recommendations
    print(f"\n{'='*70}")
    print("RECOMMENDED Tb_K BOUNDS (both salinity ends have D_ocean_km > 0)")
    print("="*70)
    for comp, _, w_list in configs:
        comp_rows = [r for r in all_rows if r["comp"] == comp]
        tb_vals = sorted(set(r["Tb_K"] for r in comp_rows))
        ocean_band = []
        for Tb in tb_vals:
            tb_rows = [r for r in comp_rows if r["Tb_K"] == Tb]
            all_built = all(r["status"] == "built" for r in tb_rows)
            all_ocean = all(
                (r["D_ocean_km"] is not None and r["D_ocean_km"] > 0)
                for r in tb_rows
            )
            if all_built and all_ocean:
                ocean_band.append(Tb)
        if ocean_band:
            print(f"  {comp}: Tb_K in [{min(ocean_band):.1f}, {max(ocean_band):.1f}] K "
                  f"({len(ocean_band)} coarse grid points)")
            # Check for HP ice in any node of the ocean band
            hp_tbs = []
            for Tb in ocean_band:
                for r in comp_rows:
                    if r["Tb_K"] == Tb and r["status"] == "built":
                        hp = (r["D_iceIII_km"] or 0) + (r["D_iceV_km"] or 0) + (r["D_iceVI_km"] or 0)
                        if hp > 0.1:
                            hp_tbs.append((Tb, r["w_ppt"], hp))
            if hp_tbs:
                print(f"    WARNING: HP ice present in ocean band:")
                for Tb, w, hp in hp_tbs:
                    print(f"      Tb={Tb:.1f} K, w={w:.1f} ppt: {hp:.1f} km HP ice")
            else:
                print(f"    No HP ice detected in ocean band nodes.")
        else:
            print(f"  {comp}: NO Tb value found with ocean at both salinity ends.")
            for Tb in tb_vals:
                tb_rows = [r for r in comp_rows if r["Tb_K"] == Tb]
                ocean_ws = [r["w_ppt"] for r in tb_rows
                            if r["status"] == "built" and r["D_ocean_km"] is not None and r["D_ocean_km"] > 0]
                frozen_ws = [r["w_ppt"] for r in tb_rows
                             if r["status"] == "built" and r["D_ocean_km"] is not None and r["D_ocean_km"] <= 0]
                fail_ws = [r["w_ppt"] for r in tb_rows if r["status"] != "built"]
                print(f"    Tb={Tb:.1f}: ocean@w={ocean_ws}  frozen@w={frozen_ws}  fail@w={fail_ws}")

    # Save JSON
    out_json = "/tmp/titan_tb_probe_results.json"

    def _serial(obj: Any) -> Any:
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return str(obj)

    with open(out_json, "w") as f:
        json.dump(all_rows, f, indent=2, default=_serial)
    print(f"\nFull results saved to {out_json}")


if __name__ == "__main__":
    main()
