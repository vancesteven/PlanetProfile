"""
Titan Phase B SBI Cache — Binding Diagnostic Gates P1 and P2.

P1: NaCl region-boundary smoothness gate.
    Build a coarse NaCl Tb x w grid, compute {C20, C22, Re_k2, Im_k2} at
    each node, finite-difference adjacent pairs, flag any step at a
    region_phases boundary that exceeds the published sigma.

P2: MgSO4 density-envelope error gate.
    Build a coarse MgSO4 Tb x w grid, identify nodes where the ocean base
    exceeds 800 MPa (the MgSO4 EOS density-table cap), and estimate the
    clamping error on C20/C22.

Usage:
    mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \\
        KMP_DUPLICATE_LIB_OK=TRUE \\
        python plans/scripts/titan_smoothness_gate.py

Author: PlanetProfile team, 2026-07-25
"""
from __future__ import annotations

import json
import logging
import os
import sys
import time
import traceback
from copy import deepcopy
from typing import Any, Dict, List, Optional, Tuple

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
for _noisy in ["PlanetProfile", "SeaFreeze", "tidalpy", "TidalPy", "lbftd"]:
    logging.getLogger(_noisy).setLevel(logging.ERROR)

log = logging.getLogger("titan_smoothness_gate")

# ---------------------------------------------------------------------------
# Published observables and sigmas (from titan_freegrav_noocean.json)
# ---------------------------------------------------------------------------
OBS_SIGMA = {
    "C20":    4.71e-7,
    "C22":    6.2e-8,
    "Re_k2":  0.048,
    "Im_k2":  0.035,
}

# Gravity config for Titan (from config metadata)
GRAVITY_REF_RADIUS_M = 2575000.0
GRAVITY_J2_OVER_C22  = 3.3333333333

# Fiducial Andrade rheology parameters (midpoints of Phase-A config bounds)
FIDUCIAL_THETA = {
    "alpha":         0.30,
    "log10_zeta":   -0.50,
    "log10_eta_Ih":  13.0,
    "log10_eta_III": 13.0,
    "log10_eta_V":   13.0,
    "log10_eta_VI":  13.0,
    "log10_eta_sil": 20.0,
}

# Arrhenius params from Phase-A config
ARRHENIUS_PARAMS = {
    "activation_energy_kJ_mol": {"Ih": 60.0},
    "R_J_mol_K": 8.314462,
}

# Template module
TEMPLATE = "PlanetProfile.Test.PPTest50"
RHO_SIL_KGM3 = 3539.0

# MgSO4 density table max (reviewer-confirmed)
MGSО4_EOS_PMAX_MPA = 800.0

# ---------------------------------------------------------------------------
# Gate grids (per task brief)
# ---------------------------------------------------------------------------
# P1 NaCl: Tb spans the HP-ice on/off transition (probe shows transition ~265 K)
TB_GRID_P1 = [244, 247, 250, 255, 260, 263, 265, 267]
W_GRID_P1  = [3.0, 30.0, 100.0]

# P2 MgSO4: salinity range 1-194 ppt; same Tb set
TB_GRID_P2 = [244, 247, 250, 255, 260, 263, 265, 267]
W_GRID_P2  = [3.0, 50.0, 150.0]

# MgSO4 max salinity that EOS supports (safe to request)
MGSО4_W_MAX = 194.0


# ---------------------------------------------------------------------------
# Build one structure and extract everything needed
# ---------------------------------------------------------------------------
def build_one(
    Tb_K: float,
    comp: str,
    w_ppt: float,
    extrap_ocean: bool = False,
) -> Dict[str, Any]:
    """
    Build a single Titan ocean structure for (Tb_K, comp, w_ppt).

    Returns a dict with:
      status, D_ocean_km, D_ice{Ih/III/V/VI}_km, P_ocean_basal_MPa,
      region_phases, r_m (array), rho (array), R_body_m, Mtot_kg,
      omega, K_Pa (array), mu_Pa (array), eta_Pa_base (array),
      phases (array), changeIndices (array), n_layers (int),
      layer_upper_radii (tuple), layer_types (tuple),
      CMR2, Tb_K, w_ppt.

    If extrap_ocean=True, sets Params.EXTRAP_OCEAN=True before the run
    (used for the P2 EXTRAP comparison).
    """
    import importlib
    from PlanetProfile.GetConfig import Params as configParams
    from PlanetProfile.Gravity.Gravity import SetupGravity
    from PlanetProfile.Main import PlanetProfile as RunPP
    from PlanetProfile.Utilities.defineStructs import Constants, EOSlist
    from PlanetProfile.Utilities.Indexing import PhaseConv

    EOSlist.loaded.clear()

    # Reload template
    if TEMPLATE in sys.modules:
        importlib.reload(sys.modules[TEMPLATE])
    else:
        importlib.import_module(TEMPLATE)
    mod = sys.modules[TEMPLATE]
    Planet = deepcopy(mod.Planet)

    # Switch to ocean-bearing mode (matching probe_one / PPTest_NaClOcean)
    Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES = False
    Planet.Do.NO_OCEAN = False
    Planet.Do.ConstantProps = getattr(Planet.Do, 'ConstantProps', {})
    Planet.Do.ConstantProps['Inner'] = True
    Planet.Sil.rhoSilWithCore_kgm3 = RHO_SIL_KGM3
    Planet.Do.POROUS_ROCK = False
    Planet.Do.Fe_CORE = False
    Planet.Do.NONHYDROSTATIC = False

    # Disable thick-shell convection (raises on thin ocean shells)
    Planet.Do.SPHERICAL_CONVECTION   = False
    Planet.Do.ARRHENIUS_VISCOSITY    = False
    Planet.Do.KALOUSOVA_CONVECTION   = False
    Planet.Do.NO_ICE_CONVECTION      = True

    # Hydrosphere settings from working template
    Planet.PfreezeUpper_MPa      = 300.0
    Planet.Ocean.PHydroMax_MPa   = 1800.0
    Planet.Ocean.THydroMax_K     = 350.0
    Planet.Ocean.deltaP          = 8.0
    Planet.Ocean.deltaT          = 0.5
    Planet.Ocean.phaseType       = 'lookup'

    # Grid parameter overrides
    Planet.Bulk.Tb_K         = Tb_K
    Planet.Bulk.Cuncertainty = 0.060
    Planet.Ocean.comp        = comp
    Planet.Ocean.wOcean_ppt  = w_ppt

    # PP config
    configParams.Gravity.backend  = "tidalpy"
    configParams.CALC_NEW         = True
    configParams.CALC_NEW_GRAVITY = False
    configParams.NO_SAVEFILE      = True
    configParams.SKIP_PLOTS       = True

    if extrap_ocean:
        configParams.EXTRAP_OCEAN = True
    else:
        # ensure default (False)
        configParams.EXTRAP_OCEAN = False

    t0 = time.time()
    try:
        Planet, Params = RunPP(Planet, configParams)
        Params.CALC_NEW_GRAVITY = True
        Planet, Params = SetupGravity(Planet, Params)
        elapsed = time.time() - t0
    except Exception as exc:
        elapsed = time.time() - t0
        tb_str = traceback.format_exc()
        return {
            "comp": comp, "Tb_K": Tb_K, "w_ppt": w_ppt,
            "status": f"FAILED: {type(exc).__name__}: {exc}",
            "traceback": tb_str[-500:],
            "elapsed_s": round(elapsed, 1),
        }

    # ---- Extract from GRAVITY model (same as cache_builder.build_single_structure) ----
    try:
        model = Planet.Gravity.ALMAModel["model"]
        cols = Planet.Gravity.columns
        rIndex   = cols.index("r")
        rhoIndex = cols.index("rho")
        VPIndex  = cols.index("VP")
        GSIndex  = cols.index("GS")
        etaIndex = cols.index("eta")
        pIndex   = cols.index("phase")

        r_m          = model[:, rIndex].astype(np.float64)
        rho_arr      = model[:, rhoIndex].astype(np.float64)
        mu_Pa_arr    = model[:, GSIndex].astype(np.float64)
        VP_ms        = model[:, VPIndex].astype(np.float64)
        eta_Pa_base  = model[:, etaIndex].astype(np.float64)
        phases_arr   = model[:, pIndex]

        K_Pa_arr = rho_arr * VP_ms ** 2 - (4.0 / 3.0) * mu_Pa_arr
        nan_mask = ~np.isfinite(K_Pa_arr) | (K_Pa_arr <= 0)
        if np.any(nan_mask):
            for i in np.where(nan_mask)[0]:
                ph = int(phases_arr[i])
                if 50 <= ph < 100:
                    nu = 0.25
                elif ph >= 100:
                    nu = 0.29
                else:
                    nu = 0.33
                K_Pa_arr[i] = 2.0 * mu_Pa_arr[i] * (1.0 + nu) / (3.0 * (1.0 - 2.0 * nu))
        K_Pa_arr = np.maximum(K_Pa_arr, 1e6)

        changeIndices = np.max(Planet.Reduced.changeIndices) - np.flipud(Planet.Reduced.changeIndices)
        n_layers = len(changeIndices) - 1

        _orig_iConv = np.flipud(Planet.Reduced.iConv)
        region_phases = []
        for i_layer in range(n_layers):
            start = changeIndices[i_layer]
            phase = phases_arr[start]
            if phase >= Constants.phaseClath and phase < Constants.phaseClath + 10:
                phase = Constants.phaseClath
            convection = _orig_iConv[start]
            phase_str = PhaseConv(phase, liq="0")
            if convection:
                phase_str += "_conv"
            region_phases.append(phase_str)

        layer_upper_radii = []
        layer_types = []
        for i_layer in range(n_layers):
            end = changeIndices[i_layer + 1]
            layer_upper_radii.append(r_m[end - 1])
            layer_types.append("liquid" if phases_arr[changeIndices[i_layer]] == 0 else "solid")

        omega = float(Planet.Bulk.meanMotion_radps)
        ecc   = float(Planet.Bulk.eccentricity)
        R_body_m  = float(Planet.Bulk.R_m)
        Mtot_kg   = float(Planet.Bulk.M_kg)

    except Exception as exc:
        elapsed = time.time() - t0
        return {
            "comp": comp, "Tb_K": Tb_K, "w_ppt": w_ppt,
            "status": f"GRAVITY_EXTRACT_FAILED: {exc}",
            "elapsed_s": round(elapsed, 1),
        }

    # ---- Layer thicknesses and ocean base pressure from Planet arrays ----
    D_ocean_km = D_iceIh_km = D_iceIII_km = D_iceV_km = D_iceVI_km = D_clath_km = 0.0
    P_ocean_basal_MPa = 0.0
    try:
        phases_raw = np.asarray(Planet.phase)
        r_raw      = np.asarray(Planet.r_m)
        P_raw      = np.asarray(Planet.P_MPa)
        n = min(len(phases_raw), len(r_raw), len(P_raw))
        phases_raw, r_raw, P_raw = phases_raw[:n], r_raw[:n], P_raw[:n]
        if n > 1:
            boundaries = np.where(np.diff(phases_raw.astype(int)) != 0)[0] + 1
            boundaries = np.concatenate(([0], boundaries, [n]))
            for i in range(len(boundaries) - 1):
                s, e = boundaries[i], boundaries[i + 1]
                ph   = int(phases_raw[s])
                thick = abs(float(r_raw[e - 1]) - float(r_raw[s])) / 1e3
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
                elif Constants.phaseClath <= ph < Constants.phaseClath + 10:
                    D_clath_km += thick
    except Exception as exc:
        log.warning(f"Phase extraction failed: {exc}")

    CMR2 = np.nan
    try:
        CMR2 = float(Planet.CMR2mean)
    except (AttributeError, TypeError):
        try:
            CMR2 = float(Planet.Bulk.Cmeasured)
        except Exception:
            pass

    return {
        "comp": comp,
        "Tb_K": Tb_K,
        "w_ppt": w_ppt,
        "extrap_ocean": extrap_ocean,
        "status": "built",
        "D_ocean_km":          round(D_ocean_km,       2),
        "D_iceIh_km":          round(D_iceIh_km,       2),
        "D_iceIII_km":         round(D_iceIII_km,      2),
        "D_iceV_km":           round(D_iceV_km,        2),
        "D_iceVI_km":          round(D_iceVI_km,       2),
        "D_clath_km":          round(D_clath_km,       2),
        "P_ocean_basal_MPa":   round(P_ocean_basal_MPa, 1),
        "region_phases": region_phases,
        "CMR2": float(CMR2) if np.isfinite(CMR2) else None,
        "r_m":           r_m.tolist(),
        "rho":           rho_arr.tolist(),
        "K_Pa":          K_Pa_arr.tolist(),
        "mu_Pa":         mu_Pa_arr.tolist(),
        "eta_Pa_base":   eta_Pa_base.tolist(),
        "phases":        phases_arr.astype(int).tolist(),
        "changeIndices": changeIndices.tolist(),
        "n_layers":      n_layers,
        "layer_upper_radii": list(layer_upper_radii),
        "layer_types":       list(layer_types),
        "R_body_m":  R_body_m,
        "Mtot_kg":   Mtot_kg,
        "omega":     omega,
        "eccentricity": ecc,
        "elapsed_s": round(elapsed, 1),
    }


# ---------------------------------------------------------------------------
# Gravity observable: C20, C22 from the built structure
# ---------------------------------------------------------------------------
def derive_gravity(node: Dict[str, Any]) -> Tuple[Optional[float], Optional[float]]:
    """Return (C20, C22) from a built node using Clairaut k_f + published Titan config."""
    from PlanetProfile.Inference.gravity_obs import clairaut_kf, hydrostatic_c20_c22
    try:
        r_arr   = np.asarray(node["r_m"],  dtype=float)
        rho_arr = np.asarray(node["rho"],  dtype=float)
        omega   = float(node["omega"])
        R_m     = float(node["R_body_m"])
        M_kg    = float(node["Mtot_kg"])
        kf = clairaut_kf(r_arr, rho_arr)
        c20, c22 = hydrostatic_c20_c22(
            kf, omega, R_m, M_kg,
            R_ref_m=GRAVITY_REF_RADIUS_M,
            j2_over_c22=GRAVITY_J2_OVER_C22,
        )
        return float(c20), float(c22)
    except Exception as exc:
        log.warning(f"derive_gravity failed for Tb={node['Tb_K']} w={node['w_ppt']}: {exc}")
        return None, None


# ---------------------------------------------------------------------------
# k2 observable: Re_k2, Im_k2 at fiducial rheology
# ---------------------------------------------------------------------------
def derive_k2(node: Dict[str, Any]) -> Tuple[Optional[float], Optional[float]]:
    """
    Compute (Re_k2, Im_k2) at fiducial Andrade rheology from the built node's
    r_m / rho / K_Pa / mu_Pa / eta arrays via the production forward model path
    (forward_models.forward_model_k2_flexible).

    IMPORTANT: forward_model_k2_flexible fires the 'Tb_K' parameter hook
    (apply_bottom_temperature), which expects structure_data to be a grid-cache
    format, NOT a bare structure dict. We use Format 2 (Tb_K_grid + structures
    list with a single entry) to satisfy the hook without building a full grid.
    This is the correct production path — the hook selects the single entry and
    passes it to the TidalPy radial solver.

    The node dict must have the full structural arrays from build_one().
    """
    from PlanetProfile.Inference.forward_models import forward_model_k2_flexible

    try:
        # Single-node inner structure dict (same fields as cache_builder output)
        inner_structure = {
            "r_m":               np.asarray(node["r_m"],          dtype=float),
            "rho":               np.asarray(node["rho"],          dtype=float),
            "K_Pa":              np.asarray(node["K_Pa"],         dtype=float),
            "mu_Pa":             np.asarray(node["mu_Pa"],        dtype=float),
            "eta_Pa_base":       np.asarray(node["eta_Pa_base"],  dtype=float),
            "bulk_visc":         np.zeros(len(node["r_m"]),       dtype=float),
            "phases":            np.asarray(node["phases"]),
            "changeIndices":     np.asarray(node["changeIndices"],dtype=int),
            "n_layers":          int(node["n_layers"]),
            "layer_upper_radii": tuple(node["layer_upper_radii"]),
            "layer_types":       tuple(node["layer_types"]),
            "region_phases":     list(node["region_phases"]),
            "omega":             float(node["omega"]),
            "eccentricity":      float(node["eccentricity"]),
            "host_mass":         1.0e26,  # not used by TidalPy radial_solver
            "a_m":               1.0e9,
            "R_body_m":          float(node["R_body_m"]),
            "Mtot_kg":           float(node["Mtot_kg"]),
            "T_K":               np.full(len(node["r_m"]), node["Tb_K"], dtype=float),
            "Tb_K":              float(node["Tb_K"]),   # stored so Arrhenius uses correct Tb
        }

        # Wrap in Format 2 (Tb_K_grid + structures) to satisfy the apply_bottom_temperature hook
        structure_data = {
            "Tb_K_grid":  np.array([float(node["Tb_K"])]),
            "structures": [inner_structure],
        }

        # Build fiducial theta dict
        theta_dict = dict(FIDUCIAL_THETA)
        theta_dict["Tb_K"] = float(node["Tb_K"])

        Re_k2, Im_k2, _, _, _ = forward_model_k2_flexible(
            theta_dict,
            structure_data,
            return_heating=False,
            arrhenius_params=ARRHENIUS_PARAMS,
        )
        return (float(Re_k2)  if np.isfinite(Re_k2)  else None,
                float(Im_k2) if np.isfinite(Im_k2) else None)
    except Exception as exc:
        log.warning(f"derive_k2 failed for Tb={node['Tb_K']} w={node['w_ppt']}: {exc}")
        return None, None


# ---------------------------------------------------------------------------
# Build grid and attach observables
# ---------------------------------------------------------------------------
def build_grid(
    comp: str,
    tb_grid: List[float],
    w_grid: List[float],
    extrap_ocean: bool = False,
    label: str = "",
) -> List[Dict[str, Any]]:
    """Build a Tb x w grid, compute all observables at each built node."""
    nodes = []
    total = len(tb_grid) * len(w_grid)
    print(f"\n{'='*70}")
    print(f"Building {label or comp} grid: {len(tb_grid)} Tb x {len(w_grid)} w = {total} nodes")
    print(f"{'='*70}")
    sys.stdout.flush()

    for Tb in tb_grid:
        for w in w_grid:
            print(f"  {comp} Tb={Tb:.0f} K  w={w:.0f} ppt  extrap={extrap_ocean} ... ",
                  end="", flush=True)
            node = build_one(Tb, comp, w, extrap_ocean=extrap_ocean)

            if node["status"] != "built" or node.get("D_ocean_km", 0) == 0:
                print(f"  {node.get('status','?')[:80]}")
                # Keep node with minimal fields for reporting
                node["C20"] = None; node["C22"] = None
                node["Re_k2"] = None; node["Im_k2"] = None
                nodes.append(node)
                continue

            # Gravity
            C20, C22 = derive_gravity(node)
            node["C20"] = C20
            node["C22"] = C22

            # k2
            Re_k2, Im_k2 = derive_k2(node)
            node["Re_k2"] = Re_k2
            node["Im_k2"] = Im_k2

            hp = node["D_iceIII_km"] + node["D_iceV_km"] + node["D_iceVI_km"]
            print(
                f"  D_ocean={node['D_ocean_km']:.0f} km  HP={hp:.0f} km  "
                f"P_base={node['P_ocean_basal_MPa']:.0f} MPa  "
                f"C22={C22:.3e}  Re_k2={Re_k2}  "
                f"t={node['elapsed_s']:.0f}s"
            )
            nodes.append(node)

    return nodes


# ---------------------------------------------------------------------------
# P1: Region-boundary smoothness analysis
# ---------------------------------------------------------------------------
def analyze_p1(nodes: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    For every adjacent pair (same w, consecutive Tb), finite-difference all
    observables. Identify pairs where region_phases DIFFERS, report steps
    vs sigma.
    """
    # Index valid (built, ocean>0) nodes by (Tb, w)
    valid = {(n["Tb_K"], n["w_ppt"]): n
             for n in nodes
             if n.get("status") == "built" and (n.get("D_ocean_km") or 0) > 0}

    tb_vals = sorted(set(k[0] for k in valid))
    w_vals  = sorted(set(k[1] for k in valid))

    obs_fields = ["C20", "C22", "Re_k2", "Im_k2"]
    boundaries = []
    all_pairs  = []

    def _region_key(n):
        return tuple(n.get("region_phases") or [])

    # Adjacent Tb pairs (same w)
    for w in w_vals:
        for i in range(len(tb_vals) - 1):
            Tb_lo, Tb_hi = tb_vals[i], tb_vals[i + 1]
            n_lo = valid.get((Tb_lo, w))
            n_hi = valid.get((Tb_hi, w))
            if n_lo is None or n_hi is None:
                continue
            steps = {}
            for obs in obs_fields:
                v_lo = n_lo.get(obs)
                v_hi = n_hi.get(obs)
                if v_lo is not None and v_hi is not None:
                    steps[obs] = abs(v_hi - v_lo)
                else:
                    steps[obs] = None

            rk_lo = _region_key(n_lo)
            rk_hi = _region_key(n_hi)
            across = (rk_lo != rk_hi)

            pair = {
                "axis": "Tb",
                "w_ppt": w,
                "Tb_lo": Tb_lo, "Tb_hi": Tb_hi,
                "region_lo": list(rk_lo),
                "region_hi": list(rk_hi),
                "across_boundary": across,
                "steps": steps,
                "step_over_sigma": {
                    obs: (steps[obs] / OBS_SIGMA[obs]
                          if steps[obs] is not None else None)
                    for obs in obs_fields
                },
            }
            all_pairs.append(pair)
            if across:
                boundaries.append(pair)

    # Adjacent w pairs (same Tb)
    for Tb in tb_vals:
        for j in range(len(w_vals) - 1):
            w_lo, w_hi = w_vals[j], w_vals[j + 1]
            n_lo = valid.get((Tb, w_lo))
            n_hi = valid.get((Tb, w_hi))
            if n_lo is None or n_hi is None:
                continue
            steps = {}
            for obs in obs_fields:
                v_lo = n_lo.get(obs)
                v_hi = n_hi.get(obs)
                if v_lo is not None and v_hi is not None:
                    steps[obs] = abs(v_hi - v_lo)
                else:
                    steps[obs] = None

            rk_lo = _region_key(n_lo)
            rk_hi = _region_key(n_hi)
            across = (rk_lo != rk_hi)

            pair = {
                "axis": "w",
                "Tb_K": Tb,
                "w_lo": w_lo, "w_hi": w_hi,
                "region_lo": list(rk_lo),
                "region_hi": list(rk_hi),
                "across_boundary": across,
                "steps": steps,
                "step_over_sigma": {
                    obs: (steps[obs] / OBS_SIGMA[obs]
                          if steps[obs] is not None else None)
                    for obs in obs_fields
                },
            }
            all_pairs.append(pair)
            if across:
                boundaries.append(pair)

    # Verdict: check observable-wise
    # Separate computeable vs uncomputed observables
    COMPUTABLE = ["C20", "C22"]  # definitely computed via derive_gravity
    # k2 may be None if TidalPy solver fails

    fails = []
    k2_computed = any(
        b["steps"].get("Re_k2") is not None
        for b in boundaries
    )

    for b in boundaries:
        for obs in obs_fields:
            step = b["steps"].get(obs)
            sigma = OBS_SIGMA[obs]
            if step is None:
                continue
            if step >= sigma:
                fails.append({
                    "obs": obs,
                    "boundary_axis": b["axis"],
                    "Tb_lo": b.get("Tb_lo") or b.get("Tb_K"),
                    "Tb_hi": b.get("Tb_hi") or b.get("Tb_K"),
                    "w": b.get("w_ppt") or b.get("Tb_K"),
                    "w_lo": b.get("w_lo"),
                    "w_hi": b.get("w_hi"),
                    "step": step,
                    "sigma": sigma,
                    "step_over_sigma": step / sigma,
                    "region_lo": b["region_lo"],
                    "region_hi": b["region_hi"],
                })

    pass_c20_c22 = all(
        (b["steps"].get(obs) is None or b["steps"].get(obs) < OBS_SIGMA[obs])
        for b in boundaries
        for obs in ["C20", "C22"]
    )
    pass_k2 = (not k2_computed) or all(
        (b["steps"].get(obs) is None or b["steps"].get(obs) < OBS_SIGMA[obs])
        for b in boundaries
        for obs in ["Re_k2", "Im_k2"]
    )

    return {
        "n_valid_nodes": len(valid),
        "n_region_boundaries": len(boundaries),
        "n_all_pairs": len(all_pairs),
        "k2_computed": k2_computed,
        "all_pairs": all_pairs,
        "boundary_pairs": boundaries,
        "fails": fails,
        "pass_C20_C22": pass_c20_c22,
        "pass_k2": pass_k2,
        "verdict": ("PASS" if (pass_c20_c22 and pass_k2)
                    else ("PASS_gravity_UNVERIFIED_k2" if (pass_c20_c22 and not k2_computed)
                          else "FAIL")),
    }


# ---------------------------------------------------------------------------
# P2: MgSO4 density-envelope error analysis
# ---------------------------------------------------------------------------
def analyze_p2(
    nodes_clamp:  List[Dict[str, Any]],
    nodes_extrap: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Compare clamped vs extrapolated nodes for MgSO4.
    Report ocean-base pressure distribution, fraction >800 MPa,
    and the C20/C22 error from clamping.
    """
    valid_clamp = {(n["Tb_K"], n["w_ppt"]): n
                   for n in nodes_clamp
                   if n.get("status") == "built" and (n.get("D_ocean_km") or 0) > 0}
    valid_extrap = {(n["Tb_K"], n["w_ppt"]): n
                    for n in nodes_extrap
                    if n.get("status") == "built" and (n.get("D_ocean_km") or 0) > 0}

    n_total  = len(valid_clamp)
    keys_over = [k for k, n in valid_clamp.items()
                 if (n.get("P_ocean_basal_MPa") or 0) > MGSО4_EOS_PMAX_MPA]
    n_over = len(keys_over)
    frac_over = n_over / n_total if n_total > 0 else 0.0

    all_P = sorted(
        [(n.get("P_ocean_basal_MPa") or 0) for n in valid_clamp.values()],
        reverse=True
    )
    max_P = all_P[0] if all_P else 0.0

    # For nodes where extrap run succeeded, compute clamping error = |C20_clamp - C20_extrap|
    clamp_errors = []
    for key in keys_over:
        nc = valid_clamp.get(key)
        ne = valid_extrap.get(key)
        if nc is None or ne is None:
            continue
        c20_clamp = nc.get("C20")
        c22_clamp = nc.get("C22")
        c20_extrap = ne.get("C20")
        c22_extrap = ne.get("C22")
        rek2_clamp  = nc.get("Re_k2")
        rek2_extrap = ne.get("Re_k2")

        err_C20 = abs(c20_clamp - c20_extrap) if (c20_clamp is not None and c20_extrap is not None) else None
        err_C22 = abs(c22_clamp - c22_extrap) if (c22_clamp is not None and c22_extrap is not None) else None
        err_Rek2 = abs(rek2_clamp - rek2_extrap) if (rek2_clamp is not None and rek2_extrap is not None) else None

        clamp_errors.append({
            "Tb_K": key[0], "w_ppt": key[1],
            "P_ocean_basal_MPa": nc.get("P_ocean_basal_MPa"),
            "err_C20":  err_C20,  "err_C22":  err_C22,  "err_Re_k2": err_Rek2,
            "err_C20_over_sigma":  (err_C20  / OBS_SIGMA["C20"])  if err_C20  is not None else None,
            "err_C22_over_sigma":  (err_C22  / OBS_SIGMA["C22"])  if err_C22  is not None else None,
            "err_Rek2_over_sigma": (err_Rek2 / OBS_SIGMA["Re_k2"]) if err_Rek2 is not None else None,
        })

    # Verdict
    if clamp_errors:
        max_C20_err = max((e["err_C20"] or 0) for e in clamp_errors)
        max_C22_err = max((e["err_C22"] or 0) for e in clamp_errors)
        max_Rek2_err = max((e["err_Re_k2"] or 0) for e in clamp_errors)
        pass_p2 = (
            max_C20_err  < OBS_SIGMA["C20"]  and
            max_C22_err  < OBS_SIGMA["C22"]  and
            max_Rek2_err < OBS_SIGMA["Re_k2"]
        )
    else:
        max_C20_err = max_C22_err = max_Rek2_err = None
        pass_p2 = None  # can't determine without extrap comparison

    # Define safe sub-grid (P_base <= 800 MPa everywhere)
    safe_keys = [k for k, n in valid_clamp.items()
                 if (n.get("P_ocean_basal_MPa") or 0) <= MGSО4_EOS_PMAX_MPA]
    safe_Tb = sorted(set(k[0] for k in safe_keys))
    safe_w  = sorted(set(k[1] for k in safe_keys))
    # Only include Tb where ALL tested w values are safe
    safe_Tb_full = [Tb for Tb in safe_Tb
                    if all((Tb, w) in safe_keys for w in safe_w)]

    return {
        "n_valid_nodes": n_total,
        "n_above_800MPa": n_over,
        "frac_above_800MPa": round(frac_over, 3),
        "max_P_ocean_basal_MPa": max_P,
        "keys_above_800MPa": [list(k) for k in sorted(keys_over)],
        "clamp_errors": clamp_errors,
        "max_C20_clamp_err": max_C20_err,
        "max_C22_clamp_err": max_C22_err,
        "max_Re_k2_clamp_err": max_Rek2_err,
        "verdict": (
            "PASS" if (n_over == 0)
            else ("PASS (clamping err < sigma)" if pass_p2 is True
                  else ("FAIL (clamping err >= sigma)" if pass_p2 is False
                        else "INDETERMINATE (extrap comparison unavailable)"))
        ),
        "safe_subgrid": {
            "Tb_K": safe_Tb_full,
            "w_ppt": safe_w,
            "description": "Tb values where ALL tested w have P_ocean_base <= 800 MPa",
        },
    }


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def _serial(obj):
    if isinstance(obj, (np.integer,)):   return int(obj)
    if isinstance(obj, (np.floating,)):  return float(obj)
    if isinstance(obj, np.ndarray):      return obj.tolist()
    return str(obj)


def strip_arrays(node: Dict[str, Any]) -> Dict[str, Any]:
    """Return a copy of a node without the large array fields (for the summary JSON)."""
    SKIP = {"r_m", "rho", "K_Pa", "mu_Pa", "eta_Pa_base", "phases",
            "changeIndices", "traceback"}
    return {k: v for k, v in node.items() if k not in SKIP}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    out_dir = "/tmp"
    raw_json  = f"{out_dir}/titan_phaseB_P1P2_raw.json"
    summ_json = f"{out_dir}/titan_phaseB_P1P2_results.json"

    print(f"\n{'#'*70}")
    print(f"# TITAN PHASE B BINDING GATE — P1 (NaCl smoothness) + P2 (MgSO4 envelope)")
    print(f"{'#'*70}")

    # ----------------------------------------------------------------
    # P1: NaCl grid
    # ----------------------------------------------------------------
    print("\n>>> TASK P1: Building NaCl grid for smoothness gate")
    nacl_nodes = build_grid("NaCl", TB_GRID_P1, W_GRID_P1, label="NaCl P1")

    p1_result = analyze_p1(nacl_nodes)

    print(f"\n--- P1 SUMMARY ---")
    print(f"  Valid nodes (built, ocean>0): {p1_result['n_valid_nodes']}")
    print(f"  Region-boundary pairs found:  {p1_result['n_region_boundaries']}")
    print(f"  k2 successfully computed:     {p1_result['k2_computed']}")
    print(f"\n  Boundary pairs:")
    for b in p1_result["boundary_pairs"]:
        axis = b["axis"]
        if axis == "Tb":
            loc = f"Tb {b['Tb_lo']}->{b['Tb_hi']} K, w={b['w_ppt']} ppt"
        else:
            loc = f"w {b['w_lo']}->{b['w_hi']} ppt, Tb={b['Tb_K']} K"
        print(f"    {loc}:")
        print(f"      region_lo = {b['region_lo']}")
        print(f"      region_hi = {b['region_hi']}")
        for obs in ["C20", "C22", "Re_k2", "Im_k2"]:
            step = b["steps"].get(obs)
            sov  = b["step_over_sigma"].get(obs)
            if step is not None:
                flag = " <-- EXCEEDS SIGMA" if step >= OBS_SIGMA[obs] else ""
                print(f"      |Δ{obs}| = {step:.3e}  sigma={OBS_SIGMA[obs]:.3e}  "
                      f"ratio={sov:.2f}{flag}")
            else:
                print(f"      |Δ{obs}| = NOT COMPUTED")

    if p1_result["fails"]:
        print(f"\n  P1 FAIL — {len(p1_result['fails'])} observable steps exceed sigma:")
        for f in p1_result["fails"]:
            print(f"    obs={f['obs']}  step={f['step']:.3e}  sigma={f['sigma']:.3e}  "
                  f"ratio={f['step_over_sigma']:.2f}")
    print(f"\n  P1 VERDICT: {p1_result['verdict']}")
    if not p1_result["k2_computed"]:
        print("  NOTE: k2 was NOT computed at any boundary — "
              "verdict covers C20/C22 only; k2 smoothness is UNVERIFIED.")

    # ----------------------------------------------------------------
    # P2: MgSO4 grid (clamp + extrap)
    # ----------------------------------------------------------------
    print("\n>>> TASK P2: Building MgSO4 grid for density-envelope gate")
    # Clamp run (EXTRAP_OCEAN=False, the default)
    mgsо4_clamp_nodes = build_grid("MgSO4", TB_GRID_P2, W_GRID_P2,
                                    extrap_ocean=False, label="MgSO4 P2 clamp")

    # Identify which nodes have P_base > 800 MPa (need extrap comparison)
    over_800 = [
        (n["Tb_K"], n["w_ppt"])
        for n in mgsо4_clamp_nodes
        if n.get("status") == "built"
        and (n.get("D_ocean_km") or 0) > 0
        and (n.get("P_ocean_basal_MPa") or 0) > MGSО4_EOS_PMAX_MPA
    ]
    print(f"\n  Nodes with P_ocean_base > 800 MPa: {len(over_800)}")
    for Tb, w in over_800:
        n = next(x for x in mgsо4_clamp_nodes if x["Tb_K"] == Tb and x["w_ppt"] == w)
        print(f"    Tb={Tb} K, w={w} ppt, P_base={n['P_ocean_basal_MPa']} MPa")

    # Extrap run only for those nodes (keep total build count low)
    mgsо4_extrap_nodes = []
    if over_800:
        print(f"\n  Building MgSO4 EXTRAP nodes ({len(over_800)} nodes) for clamping error...")
        for Tb, w in over_800:
            print(f"    MgSO4 extrap Tb={Tb} K  w={w} ppt ...", end="", flush=True)
            n = build_one(Tb, "MgSO4", w, extrap_ocean=True)
            if n["status"] == "built" and (n.get("D_ocean_km") or 0) > 0:
                C20, C22 = derive_gravity(n)
                n["C20"] = C20; n["C22"] = C22
                Re_k2, Im_k2 = derive_k2(n)
                n["Re_k2"] = Re_k2; n["Im_k2"] = Im_k2
                print(f"  C22={C22:.3e}  Re_k2={Re_k2}  t={n['elapsed_s']:.0f}s")
            else:
                print(f"  {n.get('status','?')[:80]}")
            mgsо4_extrap_nodes.append(n)
    else:
        print("  No >800 MPa nodes — extrap comparison not needed.")

    p2_result = analyze_p2(mgsо4_clamp_nodes, mgsо4_extrap_nodes)

    print(f"\n--- P2 SUMMARY ---")
    print(f"  Valid MgSO4 nodes (built, ocean>0): {p2_result['n_valid_nodes']}")
    print(f"  Nodes with P_ocean_base > 800 MPa: {p2_result['n_above_800MPa']} "
          f"({p2_result['frac_above_800MPa']*100:.1f}%)")
    print(f"  Max ocean-base pressure: {p2_result['max_P_ocean_basal_MPa']:.1f} MPa")
    if p2_result["clamp_errors"]:
        print(f"  Clamping errors at >800 MPa nodes:")
        for e in p2_result["clamp_errors"]:
            print(f"    Tb={e['Tb_K']} K, w={e['w_ppt']} ppt, P_base={e['P_ocean_basal_MPa']:.1f} MPa:")
            if e["err_C20"] is not None:
                print(f"      |ΔC20|  = {e['err_C20']:.3e}  (sigma={OBS_SIGMA['C20']:.3e}  "
                      f"ratio={e['err_C20_over_sigma']:.2f})")
            else:
                print(f"      |ΔC20|  = NOT COMPUTED")
            if e["err_C22"] is not None:
                print(f"      |ΔC22|  = {e['err_C22']:.3e}  (sigma={OBS_SIGMA['C22']:.3e}  "
                      f"ratio={e['err_C22_over_sigma']:.2f})")
            else:
                print(f"      |ΔC22|  = NOT COMPUTED")
            if e["err_Re_k2"] is not None:
                print(f"      |ΔRe_k2| = {e['err_Re_k2']:.3e}  (sigma={OBS_SIGMA['Re_k2']:.3e}  "
                      f"ratio={e['err_Rek2_over_sigma']:.2f})")
            else:
                print(f"      |ΔRe_k2| = NOT COMPUTED")
    else:
        if p2_result["n_above_800MPa"] > 0:
            print("  Extrap comparison nodes failed to build — clamping error INDETERMINATE.")

    sg = p2_result["safe_subgrid"]
    if sg["Tb_K"]:
        print(f"  Safe MgSO4 sub-grid (P_base <= 800 MPa for all tested w):")
        print(f"    Tb_K in {sg['Tb_K']}")
        print(f"    w_ppt in {sg['w_ppt']}")
    print(f"\n  P2 VERDICT: {p2_result['verdict']}")

    # ----------------------------------------------------------------
    # Save JSON outputs
    # ----------------------------------------------------------------
    # Raw: all nodes with arrays stripped for size
    raw_out = {
        "p1_nacl_nodes": [strip_arrays(n) for n in nacl_nodes],
        "p2_mgsо4_clamp_nodes": [strip_arrays(n) for n in mgsо4_clamp_nodes],
        "p2_mgsо4_extrap_nodes": [strip_arrays(n) for n in mgsо4_extrap_nodes],
    }
    with open(raw_json, "w") as f:
        json.dump(raw_out, f, indent=2, default=_serial)
    print(f"\nRaw node data saved to {raw_json}")

    # Summary: gate results
    summary = {
        "metadata": {
            "script": "plans/scripts/titan_smoothness_gate.py",
            "date": "2026-07-25",
            "TB_GRID_P1": TB_GRID_P1,
            "W_GRID_P1":  W_GRID_P1,
            "TB_GRID_P2": TB_GRID_P2,
            "W_GRID_P2":  W_GRID_P2,
            "sigmas": OBS_SIGMA,
            "fiducial_theta": FIDUCIAL_THETA,
        },
        "P1": {
            "verdict": p1_result["verdict"],
            "k2_computed": p1_result["k2_computed"],
            "n_valid_nodes": p1_result["n_valid_nodes"],
            "n_region_boundaries": p1_result["n_region_boundaries"],
            "boundary_details": p1_result["boundary_pairs"],
            "fails": p1_result["fails"],
        },
        "P2": p2_result,
    }
    with open(summ_json, "w") as f:
        json.dump(summary, f, indent=2, default=_serial)
    print(f"Summary saved to {summ_json}")

    print(f"\n{'='*70}")
    print(f"FINAL GATE VERDICTS:")
    print(f"  P1 (NaCl smoothness):         {p1_result['verdict']}")
    print(f"  P2 (MgSO4 density envelope):  {p2_result['verdict']}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
