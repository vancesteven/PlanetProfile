"""
Forward models for Bayesian inference.

Provides fast k2 and heating calculations given parameter samples and cached
Planet structure. Supports both Andrade and Maxwell rheologies with optional
Arrhenius temperature-dependent viscosity.

Designed for HPC-scale MCMC/SBI: ~0.1s per evaluation.

Author: PlanetProfile Team
Date: 2026-04-29
"""
import numpy as np
from typing import Dict, Tuple, Optional, Any, List
import logging

from PlanetProfile.Inference.grid_interp_2d import (
    is_2d_cache as _is_2d_cache,
    bilinear_weights as _bilinear_weights,
    resolve_none_corners as _resolve_none_corners,
    wOcean_ppt_from_theta as _wOcean_ppt_from_theta,
)

log = logging.getLogger('PlanetProfile')


class UnservableSampleError(ValueError):
    """A sample falls in an unbuilt region of the 2D (Tb × w) cache.

    Raised by :func:`_apply_bottom_temperature_2d` when all four bilinear
    corners are ``None`` (the tilted valid band leaves two opposite corners
    unbuilt). This is NOT a configuration error — the prior samples the full
    rectangular box including the unbuilt corners, so the likelihood and the
    SBI dataset generator must catch this SPECIFIC subclass and treat it as a
    hard reject (``-1e30`` / support-reject), exactly like the induction
    support cut. A dedicated subclass (rather than a bare ``ValueError``)
    keeps the genuine ``ValueError``s raised by the rheology hooks
    (misconfigured Andrade params, etc.) surfacing as real errors.
    """

# Eager TidalPy import (was lazy inside forward_model_k2_flexible — that
# triggered a numba "no locator available" RuntimeError when pocoMC's eval
# loop hit the first sample, because partial_melt.melting_models's
# @njit(cacheable=True) decoration could not resolve __file__ from inside
# the sampler's invocation context. Importing at module load time keeps
# numba's caching state clean.)
try:
    from TidalPy.rheology import Andrade as _TPAndrade  # noqa: F401
    from TidalPy.rheology import Maxwell as _TPMaxwell  # noqa: F401
    from TidalPy.rheology import Elastic as _TPElastic  # noqa: F401
    from TidalPy.rheology import Newton as _TPNewton  # noqa: F401
    from TidalPy.RadialSolver import radial_solver as _TPradial_solver  # noqa: F401
    from TidalPy.RadialSolver import build_rs_input_from_data as _TPbuild_rs_input  # noqa: F401
    from TidalPy.tides.multilayer.heating import (
        calc_radial_volumetric_tidal_heating_from_rs_solution as _TPheating,  # noqa: F401
    )
    _TIDALPY_OK = True
    _TIDALPY_ERR = None
except ImportError as _exc:
    _TIDALPY_OK = False
    _TIDALPY_ERR = _exc
    log.error(f"TidalPy not available: {_exc}")


# ============================================================================
# Parameter Hook System for Extensibility
# ============================================================================

_PARAMETER_HOOKS = {}


def parameter_hook(*param_ids: str):
    """
    Decorator to register parameter application function.

    Multiple param_ids can be specified for grouped parameters that should be
    applied together (e.g., 'alpha' and 'log10_zeta' for Andrade).

    The decorated function receives theta_dict and structure_data, and must
    return modified structure_data.

    Example:
        @parameter_hook('alpha', 'log10_zeta')
        def apply_andrade_params(theta_dict, structure_data):
            # Reads both alpha and log10_zeta from theta_dict
            return modified_structure_data
    """
    def decorator(func):
        for param_id in param_ids:
            _PARAMETER_HOOKS[param_id] = func
        return func
    return decorator


def apply_parameters(theta_dict: Dict[str, float], structure_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Apply all parameters to structure data via registered hooks.

    Deduplicates hook calls: each function is only called once even if it's
    registered for multiple parameters (e.g., Andrade hook for alpha + log10_zeta).

    Args:
        theta_dict: Parameter values keyed by parameter ID
        structure_data: Cached structural arrays

    Returns:
        Modified structure_data with parameter overrides applied
    """
    modified_structure = structure_data.copy()
    called_functions = set()

    # Structure-swapping hooks (Tb_K) must run before viscosity/rheology hooks
    ordered_params = sorted(theta_dict.keys(), key=lambda p: 0 if p == 'Tb_K' else 1)
    for param_id in ordered_params:
        if param_id in _PARAMETER_HOOKS:
            hook_func = _PARAMETER_HOOKS[param_id]

            # Deduplicate: only call each function once
            if hook_func not in called_functions:
                modified_structure = hook_func(theta_dict, modified_structure)
                called_functions.add(hook_func)

    return modified_structure


# ============================================================================
# Parameter Hook Implementations
# ============================================================================

@parameter_hook('alpha', 'log10_zeta', 'log10_zeta_Ih', 'log10_zeta_HP', 'log10_zeta_sil')
def apply_andrade_params(theta_dict: Dict[str, float], structure_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Apply Andrade rheology parameters (alpha, per-phase zeta).

    PlanetProfile uses different zeta values per ice phase:
    - Ice Ih: zeta=0.005 (default)
    - HP ice (III, V, VI): zeta=0.05 (default, ~10x more compliant)
    - Silicate: zeta=0.005 (default)

    The per-phase zeta is critical for reproducing observed k2 ~ 0.6 for Titan.
    """
    if 'alpha' not in theta_dict:
        raise ValueError("Andrade rheology requires 'alpha' parameter")

    has_single_zeta = 'log10_zeta' in theta_dict
    has_phase_zeta = any(k.startswith('log10_zeta_') for k in theta_dict)
    if not has_single_zeta and not has_phase_zeta:
        raise ValueError("Andrade rheology requires 'log10_zeta' or at least one log10_zeta_* parameter")

    alpha = theta_dict['alpha']

    def _zeta_to_tp(log10_zeta, alpha):
        zeta_pa = 10 ** log10_zeta
        return zeta_pa ** (1.0 / alpha) if zeta_pa != 1.0 else 1.0

    structure_data['rheology_type'] = 'andrade'
    structure_data['andrade_alpha'] = alpha

    # Per-phase zeta_tp values
    zeta_tp = {}

    # Single log10_zeta applies to all solid phases (backward-compatible with
    # Test48/Test50 single-zeta parameterisation)
    if has_single_zeta:
        val = _zeta_to_tp(theta_dict['log10_zeta'], alpha)
        zeta_tp['Ih'] = val
        zeta_tp['III'] = val
        zeta_tp['V'] = val
        zeta_tp['VI'] = val
        zeta_tp['Sil'] = val

    # Per-phase overrides take precedence over single log10_zeta
    if 'log10_zeta_Ih' in theta_dict:
        zeta_tp['Ih'] = _zeta_to_tp(theta_dict['log10_zeta_Ih'], alpha)
    if 'log10_zeta_HP' in theta_dict:
        val = _zeta_to_tp(theta_dict['log10_zeta_HP'], alpha)
        zeta_tp['III'] = val
        zeta_tp['V'] = val
        zeta_tp['VI'] = val
    if 'log10_zeta_sil' in theta_dict:
        zeta_tp['Sil'] = _zeta_to_tp(theta_dict['log10_zeta_sil'], alpha)

    # Fallback for phases not explicitly set
    zeta_tp.setdefault('Fe', 1.0)
    zeta_tp.setdefault('Clath', 1.0)
    structure_data['andrade_zeta_tp'] = zeta_tp

    return structure_data


@parameter_hook('log10_eta_Ih', 'log10_eta_III', 'log10_eta_V', 'log10_eta_VI',
                'log10_eta_HP', 'log10_eta_sil')
def apply_viscosity_params(theta_dict: Dict[str, float], structure_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Apply viscosity overrides for ice and silicate layers.

    Modifies eta_Pa_base array based on phase IDs.

    Supports two modes:
    - Lumped HP: ``log10_eta_HP`` applies to all HP ices (III, V, VI).
    - Per-phase: ``log10_eta_III``, ``log10_eta_V``, ``log10_eta_VI`` applied
      independently. Per-phase values override ``log10_eta_HP`` if both are
      present for the same phase.
    """
    eta_mod = structure_data['eta_Pa_base'].copy()
    phases = structure_data['phases']
    changeIndices = structure_data['changeIndices']
    n_layers = structure_data['n_layers']

    # Convert log viscosities
    eta_Ih = 10 ** theta_dict['log10_eta_Ih'] if 'log10_eta_Ih' in theta_dict else None
    eta_HP = 10 ** theta_dict['log10_eta_HP'] if 'log10_eta_HP' in theta_dict else None
    eta_III = 10 ** theta_dict['log10_eta_III'] if 'log10_eta_III' in theta_dict else eta_HP
    eta_V   = 10 ** theta_dict['log10_eta_V']   if 'log10_eta_V'   in theta_dict else eta_HP
    eta_VI  = 10 ** theta_dict['log10_eta_VI']  if 'log10_eta_VI'  in theta_dict else eta_HP
    eta_sil = 10 ** theta_dict['log10_eta_sil'] if 'log10_eta_sil' in theta_dict else None

    _phase_eta = {1: eta_Ih, 3: eta_III, 5: eta_V, 6: eta_VI}

    # Apply per-phase
    for i_layer in range(n_layers):
        s, e = changeIndices[i_layer], changeIndices[i_layer + 1]
        ph = int(phases[s])

        if ph in _phase_eta and _phase_eta[ph] is not None:
            eta_mod[s:e] = _phase_eta[ph]
        elif ph >= 50 and ph < 100 and eta_sil is not None:  # Silicate
            eta_mod[s:e] = eta_sil
        # Fe core, clathrates, liquid: keep baseline

    structure_data['eta_Pa'] = eta_mod
    return structure_data


@parameter_hook('log10_mu_Ih')
def apply_shear_modulus(theta_dict: Dict[str, float], structure_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Apply shear modulus override for Maxwell rheology.

    Typically used for Ice Ih layers in Maxwell viscoelastic model.
    """
    if 'log10_mu_Ih' in theta_dict:
        mu_Ih = 10 ** theta_dict['log10_mu_Ih']

        # Store for Maxwell rheology setup
        structure_data['rheology_type'] = 'maxwell'
        structure_data['mu_Ih_override'] = mu_Ih

        # Apply to Ice Ih layers
        mu_mod = structure_data['mu_Pa'].copy()
        phases = structure_data['phases']
        changeIndices = structure_data['changeIndices']
        n_layers = structure_data['n_layers']

        for i_layer in range(n_layers):
            s, e = changeIndices[i_layer], changeIndices[i_layer + 1]
            ph = int(phases[s])
            if ph == 1:  # Ice Ih
                mu_mod[s:e] = mu_Ih

        structure_data['mu_Pa'] = mu_mod

    return structure_data


# ----------------------------------------------------------------------
# B-layered blend (per-layer boundary interpolation + intra-layer resample).
# Supersedes the legacy element-wise blend, which was unsafe whenever
# (a) cell counts differed between adjacent grid points, or (b) layer
# boundary radii moved between adjacent Tb (the universal case for any
# body with an ocean — see the cache_builder transition-detection logic).
# ----------------------------------------------------------------------

_BLEND_CONT_FIELDS = ('rho', 'K_Pa', 'mu_Pa', 'eta_Pa_base',
                      'T_K', 'P_MPa', 'bulk_visc', 'sigma_Sm')
_BLEND_SCALAR_FIELDS = ('Tb_K', 'CMR2', 'J2', 'C22',
                        'rhoSil_kgm3', 'D_hsphere_km',
                        'D_iceIh_km', 'D_iceIII_km', 'D_iceV_km', 'D_iceVI_km',
                        # Ocean thickness (cache_builder writes this as
                        # D_hsphere − sum(D_ice_phases); needed by the runner's
                        # per-sample D_ocean_results recomputation and by the
                        # corner / layers_vs_docean plotters)
                        'D_ocean_km')
_BLEND_META_FIELDS = ('omega', 'eccentricity', 'host_mass', 'a_m',
                      'R_body_m', 'Mtot_kg', 'phases', 'changeIndices',
                      'n_layers', 'layer_types', 'region_phases',
                      # Phase-code → label dict required by compute_heating;
                      # identical across all per-Tb structures, so copy from s0.
                      'phase_map',
                      # Induction body-fixed props copied from s0
                      'Texc_hr')
# Per-region arrays from PP's induction setup. Element-wise blend if
# s0 and s1 have matching lengths (same number of conducting regions —
# guaranteed by the matching-region_phases precondition); copy from s0
# otherwise.
_BLEND_PER_REGION_ARRAY_FIELDS = ('rSigChange_m', 'sigmaLayers_Sm')


def _regions_match(s0: Dict[str, Any], s1: Dict[str, Any]) -> bool:
    """True if two structures share the same layer set."""
    return (
        s0.get('region_phases') == s1.get('region_phases')
        and s0.get('n_layers') == s1.get('n_layers')
    )


def _blend_b_layered(s0: Dict[str, Any], s1: Dict[str, Any], w: float) -> Dict[str, Any]:
    """B-layered blend — per-layer boundary motion + intra-layer resample.

    Preconditions (caller's responsibility): ``_regions_match(s0, s1)`` is
    True, both structures carry ``r_m``, ``changeIndices``, and at least
    one of the continuous fields.

    Algorithm: for each layer index k, blend the layer's boundary radii
    (`r_lo`, `r_hi`) linearly in `w`. Build the target intra-layer cell
    grid by taking s0's normalised cell positions and rescaling to the
    blended boundaries. Resample s1's layer values onto that target grid
    via `np.interp`, then linear-blend with s0's values cell-wise.

    Phase identity, cell count, layer count, and layer types come from
    s0 (identical to s1's by the matching-regions precondition). Scalar
    body parameters (``omega``, ``R_body_m``, …) are Tb-invariant; copied
    from s0. Layer boundary radii (``layer_upper_radii``) are blended.
    Per-Tb scalars (``CMR2``, ``rhoSil_kgm3``, ``D_iceIh_km``, …) are
    linear-blended.
    """
    n_layers = int(s0['n_layers'])
    s0_ci = np.asarray(s0['changeIndices'], dtype=int)
    s1_ci = np.asarray(s1['changeIndices'], dtype=int)
    s0_r = np.asarray(s0['r_m'], dtype=float)
    s1_r = np.asarray(s1['r_m'], dtype=float)

    # Precondition: every per-cell field in either structure must have
    # the same length as that structure's r_m. Catches cache-builder
    # padding bugs surgically — without this, np.interp deep inside
    # the loop raises "fp and xp are not of the same length" with no
    # hint which field is at fault.
    for label, struct, r_arr in (('s0', s0, s0_r), ('s1', s1, s1_r)):
        n = len(r_arr)
        if int(struct.get('changeIndices', [0])[-1]) != n:
            raise ValueError(
                f"_blend_b_layered: {label}.changeIndices[-1] "
                f"({int(struct['changeIndices'][-1])}) != len(r_m) ({n}); "
                "cache is internally inconsistent."
            )
        for field in _BLEND_CONT_FIELDS:
            if field in struct:
                m = len(np.asarray(struct[field]))
                if m != n:
                    raise ValueError(
                        f"_blend_b_layered: {label}['{field}'] has length "
                        f"{m} but len(r_m)={n}. Cache-builder padding likely "
                        f"missed this field. Re-run cache build after fixing "
                        f"build_single_structure's per-cell padding pass."
                    )

    n_cells = len(s0_r)
    out_r = np.empty(n_cells, dtype=float)
    out_cont: Dict[str, np.ndarray] = {}
    for field in _BLEND_CONT_FIELDS:
        if field in s0 and field in s1:
            out_cont[field] = np.empty(n_cells, dtype=float)

    for k in range(n_layers):
        s0_lo, s0_hi = int(s0_ci[k]), int(s0_ci[k + 1])
        s1_lo, s1_hi = int(s1_ci[k]), int(s1_ci[k + 1])
        if s0_hi <= s0_lo or s1_hi <= s1_lo:
            continue  # degenerate layer; skip (shouldn't happen post-MIN_POINTS pad)

        # Layer boundary radii in s0 / s1
        r0_lo_b = float(s0_r[s0_lo])
        r0_hi_b = float(s0_r[s0_hi - 1])
        r1_lo_b = float(s1_r[s1_lo])
        r1_hi_b = float(s1_r[s1_hi - 1])
        # Blended boundaries
        r_lo = (1.0 - w) * r0_lo_b + w * r1_lo_b
        r_hi = (1.0 - w) * r0_hi_b + w * r1_hi_b

        # Target intra-layer cell grid: s0's normalised positions, rescaled.
        if r0_hi_b > r0_lo_b:
            t0 = (s0_r[s0_lo:s0_hi] - r0_lo_b) / (r0_hi_b - r0_lo_b)
        else:
            t0 = np.linspace(0.0, 1.0, s0_hi - s0_lo)
        out_r[s0_lo:s0_hi] = r_lo + t0 * (r_hi - r_lo)

        # s1's normalised positions for resampling.
        if r1_hi_b > r1_lo_b:
            t1_src = (s1_r[s1_lo:s1_hi] - r1_lo_b) / (r1_hi_b - r1_lo_b)
        else:
            t1_src = np.linspace(0.0, 1.0, s1_hi - s1_lo)
        # np.interp wants xp ascending — both t arrays are ascending by
        # construction (r_m is ascending within each layer).

        for field, out_arr in out_cont.items():
            s0_val = np.asarray(s0[field], dtype=float)[s0_lo:s0_hi]
            s1_val = np.asarray(s1[field], dtype=float)[s1_lo:s1_hi]
            s1_resampled = np.interp(t0, t1_src, s1_val)
            out_arr[s0_lo:s0_hi] = (1.0 - w) * s0_val + w * s1_resampled

    out: Dict[str, Any] = {'r_m': out_r}
    out.update(out_cont)

    # Layer-set / phase / metadata fields: take from s0 (identical to s1)
    for k in _BLEND_META_FIELDS:
        if k in s0:
            v = s0[k]
            out[k] = v.copy() if isinstance(v, np.ndarray) else v

    # Per-Tb scalar fields: linear blend
    for k in _BLEND_SCALAR_FIELDS:
        if k in s0 and k in s1:
            try:
                out[k] = (1.0 - w) * float(s0[k]) + w * float(s1[k])
            except (TypeError, ValueError):
                out[k] = s0[k]

    # layer_upper_radii: linear-blend the tuple element-wise
    if 'layer_upper_radii' in s0 and 'layer_upper_radii' in s1:
        a, b = s0['layer_upper_radii'], s1['layer_upper_radii']
        if len(a) == len(b):
            out['layer_upper_radii'] = tuple(
                (1.0 - w) * float(av) + w * float(bv) for av, bv in zip(a, b)
            )
        else:
            out['layer_upper_radii'] = a  # fallback, shouldn't happen post-precondition

    # Per-region induction arrays (rSigChange_m, sigmaLayers_Sm): same-region
    # bracket → element-wise linear blend; mismatched length → copy from s0.
    for field in _BLEND_PER_REGION_ARRAY_FIELDS:
        a = s0.get(field)
        b = s1.get(field)
        if a is None and b is None:
            continue
        if a is None or b is None:
            out[field] = a if a is not None else b
            continue
        a_arr = np.asarray(a, dtype=float)
        b_arr = np.asarray(b, dtype=float)
        if a_arr.shape == b_arr.shape:
            out[field] = (1.0 - w) * a_arr + w * b_arr
        else:
            out[field] = a_arr  # length mismatch — use s0 (defensive)

    return out


def _copy_struct(s: Dict[str, Any]) -> Dict[str, Any]:
    """Shallow copy with deep-copy of ndarray fields."""
    return {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in s.items()}


def _apply_bottom_temperature_2d(
    theta_dict: Dict[str, float], structure_data: Dict[str, Any]
) -> Dict[str, Any]:
    """Bilinear (Tb, log10 w) structural interpolation for the v3.0 2D cache.

    Uses the shared :mod:`grid_interp_2d` interpolant: the four bracketing
    corners and bilinear weights are computed once; ``None`` corners (the
    tilted valid band leaves two opposite corners unbuilt) are dropped and the
    weights renormalized. The surviving corners are combined as nested 1-D
    ``_blend_b_layered`` calls (blend the two w-corners at each Tb row, then
    blend the two rows), falling back to the highest-weight valid corner when
    region-sets differ (the same across-transition nearest-neighbour rule the
    1D path uses). If NO corner is valid the sample is unservable → raise so
    the likelihood rejects it (mirrors the support-guard REJECT).
    """
    Tb_grid = np.asarray(structure_data['Tb_K_grid'], dtype=float)
    w_grid = np.asarray(structure_data['wOcean_ppt_grid'], dtype=float)
    structs = structure_data['structures']
    w_ppt = _wOcean_ppt_from_theta(theta_dict)

    corners, weights = _bilinear_weights(Tb_grid, w_grid, Tb_K_of(theta_dict), w_ppt)
    valid = [structs[c] is not None for c in corners]
    resolved = _resolve_none_corners(corners, weights, valid)
    if resolved is None:
        raise UnservableSampleError(
            f"2D cache: all bilinear corners are None at "
            f"(Tb={theta_dict.get('Tb_K')}, w={w_ppt:.4g} ppt); "
            "sample sits in an unbuilt corner — reject."
        )
    kept_corners, kept_w = resolved

    def _blend_pair(cA, wA, cB, wB):
        """Blend two corners by their (already-normalized-within-pair) weights."""
        sA, sB = structs[cA], structs[cB]
        tot = wA + wB
        if tot <= 0:
            return _copy_struct(sA)
        frac_B = wB / tot
        if _regions_match(sA, sB):
            return _blend_b_layered(sA, sB, frac_B)
        # across-transition: nearest (higher-weight) corner
        return _copy_struct(sA if wA >= wB else sB)

    if len(kept_corners) == 1:
        out = _copy_struct(structs[kept_corners[0]])
    elif len(kept_corners) == 2:
        out = _blend_pair(kept_corners[0], kept_w[0],
                          kept_corners[1], kept_w[1])
    else:
        # 3 or 4 corners: blend within each Tb row (corners share i_Tb when
        # their flat index // n_w matches), then blend the row results by
        # their summed weights. Grouping by i_Tb keeps w-blends (same region
        # set, monotic freezing) separate from the Tb-blend.
        n_w = w_grid.size
        rows: Dict[int, list] = {}
        for c, wt in zip(kept_corners, kept_w):
            rows.setdefault(c // n_w, []).append((c, wt))
        row_structs = []
        row_weights = []
        for i_Tb, members in rows.items():
            if len(members) == 1:
                (c, wt), = members
                row_structs.append(_copy_struct(structs[c]))
                row_weights.append(wt)
            else:
                (cA, wA), (cB, wB) = members[0], members[1]
                row_structs.append(_blend_pair(cA, wA, cB, wB))
                row_weights.append(wA + wB)
        # blend the (≤2) row results across Tb
        if len(row_structs) == 1:
            out = row_structs[0]
        else:
            sA, sB = row_structs[0], row_structs[1]
            wA, wB = row_weights[0], row_weights[1]
            tot = wA + wB
            frac_B = 0.0 if tot <= 0 else wB / tot
            if _regions_match(sA, sB):
                out = _blend_b_layered(sA, sB, frac_B)
            else:
                out = sA if wA >= wB else sB

    # Carry grid references + sampled values for downstream hooks.
    out['Tb_K_grid'] = Tb_grid
    out['wOcean_ppt_grid'] = w_grid
    out['structures'] = structs
    out['Tb_K_sampled'] = theta_dict['Tb_K']
    out['wOcean_ppt_sampled'] = w_ppt
    return out


def Tb_K_of(theta_dict: Dict[str, float]) -> float:
    """Tb_K accessor (small shim so _apply_bottom_temperature_2d reads clean)."""
    return float(theta_dict['Tb_K'])


@parameter_hook('Tb_K')
def apply_bottom_temperature(theta_dict: Dict[str, float], structure_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Apply bottom temperature by selecting/interpolating from a grid cache.

    Supports two grid cache formats:

    1. **Dict-keyed (MCMCRunner default)**: ``structure_data`` contains
       ``'grid_cache'`` (dict[float -> structure_dict]) and
       ``'grid_Tb_values'`` (sorted 1-D array). Nearest-neighbour selection.

    2. **List format (Test50 + v2.1 transition-aware)**: ``structure_data``
       contains ``'Tb_K_grid'`` (1-D array) and ``'structures'`` (list of
       structure dicts). Bracketing grid points are blended via the
       B-layered scheme (per-layer boundary motion + intra-layer
       resample) when their layer sets match. When the bracket straddles
       a layer-set transition (different ``region_phases``), the runner
       falls back to nearest-neighbour selection — these intervals are
       narrow (cache_builder pins them to ε_T = 0.01 K) so the
       discontinuity sits well below typical posterior width.
    """
    Tb_K = theta_dict['Tb_K']

    # --- Format 3: v3.0 2D (Tb × wOcean_ppt) grid ---
    # Salinity sampled as log10_wOcean_ppt; interpolate bilinearly in
    # (Tb, log10 w) using the SHARED interpolant so this structural blend and
    # the induction Ae blend (mcmc_runner) agree on corners/weights/None-policy.
    if _is_2d_cache(structure_data):
        return _apply_bottom_temperature_2d(theta_dict, structure_data)

    # --- Format 2: Test50 list grid (now v2.1 transition-aware) ---
    if 'Tb_K_grid' in structure_data and 'structures' in structure_data:
        Tb_grid = structure_data['Tb_K_grid']
        structs = structure_data['structures']
        Tb_clamped = float(np.clip(Tb_K, Tb_grid[0], Tb_grid[-1]))
        j = int(np.searchsorted(Tb_grid, Tb_clamped))
        if j <= 0:
            out = _copy_struct(structs[0])
        elif j >= len(Tb_grid):
            out = _copy_struct(structs[-1])
        else:
            t0_v, t1_v = Tb_grid[j - 1], Tb_grid[j]
            w = 0.0 if t1_v == t0_v else (Tb_clamped - t0_v) / (t1_v - t0_v)
            s0, s1 = structs[j - 1], structs[j]
            if _regions_match(s0, s1):
                out = _blend_b_layered(s0, s1, w)
            else:
                # Across-transition bracket: pick the nearer grid point.
                # cache_builder narrows transitions to ε_T = 0.01 K so the
                # mid-bracket discontinuity is much smaller than typical
                # Tb posterior widths.
                nearer = s0 if (Tb_clamped - t0_v) <= (t1_v - Tb_clamped) else s1
                out = _copy_struct(nearer)
        # Carry grid reference so downstream hooks can still access it
        out['Tb_K_grid'] = Tb_grid
        out['structures'] = structs
        if 'transitions' in structure_data:
            out['transitions'] = structure_data['transitions']
        # Store sampled Tb as reference for Arrhenius hook
        out['Tb_K_sampled'] = Tb_K
        return out

    # --- Format 1: dict-keyed grid ---
    if 'grid_cache' not in structure_data:
        raise ValueError(
            "Parameter 'Tb_K' requires structure grid cache. "
            "Use MCMCRunner._load_grid_cache() or provide 'Tb_K_grid'+'structures' format."
        )
    grid_cache = structure_data['grid_cache']
    grid_Tb_values = structure_data['grid_Tb_values']
    idx = np.argmin(np.abs(grid_Tb_values - Tb_K))
    nearest_Tb = grid_Tb_values[idx]
    selected = grid_cache[nearest_Tb].copy()
    selected['grid_cache'] = grid_cache
    selected['grid_Tb_values'] = grid_Tb_values
    selected['Tb_K_sampled'] = Tb_K
    return selected


# ============================================================================
# Flexible Forward Model (Dict-based Parameters)
# ============================================================================

def forward_model_k2_flexible(
    theta_dict: Dict[str, float],
    structure_data: Dict[str, Any],
    return_heating: bool = False,
    arrhenius_params: Optional[Dict[str, Any]] = None
) -> Tuple[float, float, float, float, Optional[Dict[str, float]]]:
    """
    Flexible forward model using dict-based parameters.

    This is the extensible entry point that uses the parameter hook system.
    New parameters can be added by registering hooks - no changes needed here.

    Args:
        theta_dict: Parameters as dict (e.g., {'alpha': 0.3, 'log10_eta_Ih': 14.0})
        structure_data: Cached structural arrays from load_structure_cache()
        return_heating: If True, compute per-phase tidal heating (slower)
        arrhenius_params: Optional Arrhenius temperature-dependent viscosity params

    Returns:
        ``(Re_k2, Im_k2, Re_h2, Im_h2, perPhase_W)`` where ``perPhase_W`` is
        ``None`` if ``return_heating=False``. Returns ``(nan, nan, nan, nan,
        None)`` if the TidalPy solver fails. Both Love numbers are extracted
        from the same ``radial_solver`` result; the h₂ extraction is
        essentially free.
    """
    # TidalPy is imported at module top (eager) — see header comment.
    if not _TIDALPY_OK:
        log.error(f"TidalPy not available: {_TIDALPY_ERR}")
        return np.nan, np.nan, np.nan, np.nan, None
    Andrade, Maxwell, Elastic, Newton = _TPAndrade, _TPMaxwell, _TPElastic, _TPNewton
    radial_solver = _TPradial_solver
    build_rs_input_from_data = _TPbuild_rs_input
    calc_radial_volumetric_tidal_heating_from_rs_solution = _TPheating

    # Apply parameter hooks
    modified_structure = apply_parameters(theta_dict, structure_data)

    # Extract modified arrays
    eta_mod = modified_structure.get('eta_Pa', modified_structure['eta_Pa_base'].copy())
    phases = modified_structure['phases']
    changeIndices = modified_structure['changeIndices']
    n_layers = modified_structure['n_layers']

    # Apply Arrhenius temperature-dependence if requested
    if arrhenius_params is not None:
        # When Tb_K is sampled, use it as the Arrhenius reference temperature so
        # that the sampled eta_Ih is interpreted as the basal (bottom) viscosity —
        # identical to Test50's inline Arrhenius logic.  The sampled Tb_K is
        # stored in structure by the Tb_K parameter hook as 'Tb_K_sampled'.
        arrhenius_params_effective = dict(arrhenius_params)
        if 'Tb_K_sampled' in modified_structure and 'reference_temp_K' not in arrhenius_params:
            arrhenius_params_effective['reference_temp_K'] = float(
                modified_structure['Tb_K_sampled'])
        eta_mod = apply_arrhenius_viscosity(
            eta_mod,
            phases,
            changeIndices,
            n_layers,
            modified_structure.get('T_K'),
            arrhenius_params_effective
        )

    # Build rheology models per layer based on applied parameters
    rheology_type = modified_structure.get('rheology_type', 'maxwell')

    # Per-layer rheology assignment must mirror PP's Params.Gravity.rheology_models
    # (defined in PlanetProfile/Gravity/defaultConfigGravity.py:26): Fe core is
    # ELASTIC, ocean is NEWTON (viscous fluid), clathrate is ELASTIC, ices/silicate
    # are ANDRADE/MAXWELL. Treating Fe as Andrade (the previous behaviour) caused
    # TidalPy radial_solver to fail at "layer 0" with "step size less than spacing
    # between numbers" because Andrade applied to a metal core produces a near-
    # singular complex shear modulus the integrator can't traverse. This bug was
    # hidden by Ganymede/Callisto's higher tolerance and surfaced as 100% k2-NaN
    # for Europa during C2 smoke MCMC (2026-05-22).
    _ELASTIC_PHASES = frozenset({'Fe', 'Clath'})
    _NEWTON_PHASES = frozenset({'0'})  # liquid ocean

    if rheology_type == 'andrade':
        alpha = modified_structure['andrade_alpha']
        zeta_tp_map = modified_structure['andrade_zeta_tp']
        shear_models = []
        for rp in modified_structure['region_phases']:
            base_phase = rp.replace('_conv', '')
            if base_phase in _ELASTIC_PHASES:
                shear_models.append(Elastic())
            elif base_phase in _NEWTON_PHASES:
                shear_models.append(Newton())
            else:
                phase_zeta = zeta_tp_map.get(base_phase, 1.0)
                shear_models.append(Andrade(args=(alpha, phase_zeta)))

    elif rheology_type == 'maxwell':
        shear_models = []
        for rp in modified_structure['region_phases']:
            base_phase = rp.replace('_conv', '')
            if base_phase in _ELASTIC_PHASES:
                shear_models.append(Elastic())
            elif base_phase in _NEWTON_PHASES:
                shear_models.append(Newton())
            else:
                shear_models.append(Maxwell())
    else:
        raise ValueError(f"Unknown rheology type '{rheology_type}'. Must be 'andrade' or 'maxwell'")

    # Bulk rheology: always elastic
    bulk_models = [Elastic() for _ in shear_models]

    # Run TidalPy radial solver
    try:
        build_data = build_rs_input_from_data(
            modified_structure['omega'],
            modified_structure['r_m'],
            modified_structure['rho'],
            modified_structure['K_Pa'],
            modified_structure['mu_Pa'],
            modified_structure['bulk_visc'],
            np.ascontiguousarray(eta_mod),
            modified_structure['layer_upper_radii'],
            modified_structure['layer_types'],
            tuple([False] * n_layers),  # is_static_bylayer
            tuple([False] * n_layers),  # is_incompressible_bylayer
            tuple(shear_models),
            tuple(bulk_models),
            perform_checks=False,
            warnings=False
        )

        result = radial_solver(
            build_data.radius_array,
            build_data.density_array,
            build_data.complex_bulk_modulus_array,
            build_data.complex_shear_modulus_array,
            build_data.frequency,
            build_data.planet_bulk_density,
            build_data.layer_types,
            build_data.is_static_bylayer,
            build_data.is_incompressible_bylayer,
            build_data.upper_radius_bylayer_array,
            degree_l=2,
            solve_for=('tidal',),
            verbose=False,
            raise_on_fail=False
        )

        if not result.success:
            return np.nan, np.nan, np.nan, np.nan, None

        # Extract k2 and h2 from the same RadialSolver result.
        k2 = complex(result.k)
        Re_k2 = k2.real
        Im_k2 = k2.imag
        try:
            h2 = complex(result.h)
            Re_h2 = h2.real
            Im_h2 = h2.imag
        except AttributeError:
            log.debug("RadialSolver result has no .h attribute; h2 unavailable.")
            Re_h2 = np.nan
            Im_h2 = np.nan

        # Compute per-phase heating if requested. Guard separately from the
        # outer except so a heating-only failure (a) does not silently nuke
        # the k2/h2 result and (b) is LOUD — a broken heating import chain
        # previously vanished into the outer log.debug and returned empty
        # heating for every sample with no visible symptom.
        perPhase_W = None
        if return_heating and modified_structure['eccentricity'] > 0:
            try:
                perPhase_W = compute_heating(result, modified_structure)
            except Exception as heat_exc:
                log.warning(f"Per-phase heating computation failed "
                            f"(k2/h2 unaffected): {heat_exc}")

        return Re_k2, Im_k2, Re_h2, Im_h2, perPhase_W

    except Exception as e:
        log.debug(f"Forward model failed: {e}")
        return np.nan, np.nan, np.nan, np.nan, None


# ============================================================================
# Legacy Array-based Forward Model (Backward Compatible)
# ============================================================================

def forward_model_k2(
    theta: np.ndarray,
    structure_data: Dict[str, Any],
    rheology: str = 'andrade',
    return_heating: bool = False,
    arrhenius_params: Optional[Dict[str, Any]] = None
) -> Tuple[float, float, Optional[Dict[str, float]]]:
    """
    Compute tidal Love number k2 (and optionally heating) for given parameters.

    LEGACY INTERFACE: Backward compatible array-based parameter format.
    Internally wraps forward_model_k2_flexible() with dict conversion.

    Args:
        theta: Parameter vector. Format depends on rheology:
            - 'andrade' (7-param): [alpha, log10_zeta_Ih, log10_zeta_HP, log10_zeta_sil,
                                     log10_eta_Ih, log10_eta_HP, log10_eta_sil]
            - 'andrade' (5-param, no HP ice): [alpha, log10_zeta_Ih, log10_zeta_sil,
                                                log10_eta_Ih, log10_eta_sil]
            - 'maxwell' (3-param): [log10_eta_Ih, log10_eta_HP, log10_eta_sil]
            - 'maxwell' (4-param, with shear modulus): [log10_mu_Ih, log10_eta_Ih, log10_eta_HP, log10_eta_sil]
        structure_data: Cached structural arrays from load_structure_cache()
        rheology: 'andrade' or 'maxwell'
        return_heating: If True, compute per-phase tidal heating (slower)
        arrhenius_params: Optional dict with {'activation_energy_kJ_mol': Dict[phase: E],
                                             'reference_temp_K': float, ...}

    Returns:
        (Re_k2, Im_k2, perPhase_W) where perPhase_W is None if return_heating=False
        Returns (np.nan, np.nan, None) if TidalPy solver fails
    """
    # Convert theta array to dict based on rheology
    if rheology == 'andrade':
        if len(theta) == 7:
            theta_dict = {
                'alpha': theta[0],
                'log10_zeta_Ih': theta[1],
                'log10_zeta_HP': theta[2],
                'log10_zeta_sil': theta[3],
                'log10_eta_Ih': theta[4],
                'log10_eta_HP': theta[5],
                'log10_eta_sil': theta[6]
            }
        elif len(theta) == 5:
            # Europa-style: no HP ice layer (alpha, zeta_Ih, zeta_sil, eta_Ih, eta_sil)
            theta_dict = {
                'alpha': theta[0],
                'log10_zeta_Ih': theta[1],
                'log10_zeta_sil': theta[2],
                'log10_eta_Ih': theta[3],
                'log10_eta_sil': theta[4]
            }
        else:
            raise ValueError(f"Andrade rheology expects 5 or 7 parameters, got {len(theta)}")

    elif rheology == 'maxwell':
        if len(theta) == 4:
            # With shear modulus override: [log10_mu_Ih, log10_eta_Ih, log10_eta_HP, log10_eta_sil]
            theta_dict = {
                'log10_mu_Ih': theta[0],
                'log10_eta_Ih': theta[1],
                'log10_eta_HP': theta[2],
                'log10_eta_sil': theta[3]
            }
        elif len(theta) == 3:
            theta_dict = {
                'log10_eta_Ih': theta[0],
                'log10_eta_HP': theta[1],
                'log10_eta_sil': theta[2]
            }
        else:
            raise ValueError(f"Maxwell rheology expects 3 or 4 parameters, got {len(theta)}")
    else:
        raise ValueError(f"Unknown rheology '{rheology}'. Must be 'andrade' or 'maxwell'")

    # Call flexible version
    return forward_model_k2_flexible(
        theta_dict,
        structure_data,
        return_heating=return_heating,
        arrhenius_params=arrhenius_params
    )


def compute_heating(
    rs_solution,
    structure_data: Dict[str, Any]
) -> Dict[str, float]:
    """
    Compute per-phase tidal heating from TidalPy radial solver solution.

    Args:
        rs_solution: TidalPy RadialSolver result object
        structure_data: Cached structural arrays

    Returns:
        Dict mapping phase names (e.g., 'Ih', 'III', 'V', 'VI', 'Sil') to total
        heating in that phase (Watts)
    """
    from TidalPy.tides.multilayer.heating import calc_radial_volumetric_tidal_heating_from_rs_solution
    from PlanetProfile.Thermodynamics.Electrical import PhaseConv

    # Compute volumetric heating profile
    heating_profile = calc_radial_volumetric_tidal_heating_from_rs_solution(
        structure_data['eccentricity'],
        structure_data['omega'],
        structure_data['a_m'],
        structure_data['host_mass'],
        rs_solution,
        perform_checks=False
    )

    result_radii = np.asarray(rs_solution.radius_array)
    r_m_local = structure_data['r_m']
    phases = structure_data['phases']
    changeIndices = structure_data['changeIndices']
    n_layers = structure_data['n_layers']
    phase_map = structure_data.get('phase_map', {})

    # Integrate heating over each layer
    perPhase_W = {}
    for i_layer in range(n_layers):
        s_idx = changeIndices[i_layer]
        e_idx = changeIndices[i_layer + 1]
        layer_phase = int(phases[s_idx])

        # Map numeric phase to string
        phase_str = phase_map.get(layer_phase, PhaseConv(layer_phase, liq='0'))

        # Find radii in this layer
        r_lo = r_m_local[s_idx]
        r_hi = r_m_local[e_idx - 1]
        mask = (result_radii >= r_lo - 1.0) & (result_radii <= r_hi + 1.0)

        if np.any(mask):
            lr = result_radii[mask]
            lh = heating_profile[mask]

            # Integrate heating: ∫ heating * 4πr² dr
            if len(lr) > 1:
                total_power = np.trapezoid(lh * 4.0 * np.pi * lr**2, lr)
            else:
                # Single point: use layer volume
                V = (4.0/3.0) * np.pi * (r_hi**3 - r_lo**3)
                total_power = lh[0] * V

            # Accumulate (same phase may appear in multiple layers)
            perPhase_W[phase_str] = perPhase_W.get(phase_str, 0) + total_power

    # Aggregate into Titan reservoirs for heat partitioning tracking
    # Groupings follow standard ocean world definitions:
    # - Silicate: Silicate matrix + Fe core (typically negligible in tides, but included)
    # - HP Ices: Ice III, V, VI and their mixed clathrate variants
    # - Ice Ih: Ice Ih shell and its clathrate/mixed variants
    hp_phases = {'III', 'V', 'VI', 'MixedClathrateIII', 'MixedClathrateV', 'MixedClathrateVI'}
    ih_phases = {'Ih', 'Clath', 'MixedClathrateIh'}

    perPhase_W['Silicate_W'] = perPhase_W.get('Sil', 0.0) + perPhase_W.get('Fe', 0.0)
    perPhase_W['HP_Ice_W'] = sum(perPhase_W.get(p, 0.0) for p in hp_phases)
    perPhase_W['Ice_Ih_W'] = sum(perPhase_W.get(p, 0.0) for p in ih_phases)

    return perPhase_W


def apply_arrhenius_viscosity(
    eta_base: np.ndarray,
    phases: np.ndarray,
    changeIndices: np.ndarray,
    n_layers: int,
    T_K: Optional[np.ndarray],
    arrhenius_params: Dict[str, Any]
) -> np.ndarray:
    """
    Apply Arrhenius temperature-dependent viscosity scaling.

    η(T) = η_ref * exp(E/R * (1/T - 1/T_ref))

    Args:
        eta_base: Base viscosity array (Pa·s)
        phases: Phase ID array
        changeIndices: Layer boundary indices
        n_layers: Number of layers
        T_K: Temperature profile (K). If None, returns eta_base unchanged
        arrhenius_params: Dict with:
            - 'activation_energy_kJ_mol': Dict[phase_name: E_kJ_mol]
            - 'reference_temp_K': float (default 273.15 K)
            - 'R_J_mol_K': float (default 8.314 J/mol/K)

    Returns:
        Modified viscosity array with Arrhenius scaling applied
    """
    if T_K is None:
        log.warning("Temperature profile not available, Arrhenius viscosity not applied")
        return eta_base

    eta_mod = eta_base.copy()
    activation_energies = arrhenius_params.get('activation_energy_kJ_mol', {})
    T_ref = arrhenius_params.get('reference_temp_K', 273.15)
    R = arrhenius_params.get('R_J_mol_K', 8.314)

    # Phase name mapping (numeric -> string)
    phase_names = {1: 'Ih', 3: 'III', 5: 'V', 6: 'VI'}

    for i_layer in range(n_layers):
        s, e = changeIndices[i_layer], changeIndices[i_layer + 1]
        ph = int(phases[s])

        # Check if this phase has Arrhenius parameters
        phase_name = phase_names.get(ph)
        if phase_name is None:
            continue  # Not an ice phase with Arrhenius

        E_kJ_mol = activation_energies.get(phase_name)
        if E_kJ_mol is None:
            continue  # No activation energy specified for this phase

        # Apply Arrhenius scaling
        E_J_mol = E_kJ_mol * 1000.0
        T_layer = T_K[s:e]

        # Avoid division by zero or negative temperatures
        T_layer = np.maximum(T_layer, 1.0)

        # η(T) = η_ref * exp(E/R * (1/T - 1/T_ref))
        exponent = (E_J_mol / R) * (1.0 / T_layer - 1.0 / T_ref)
        eta_mod[s:e] *= np.exp(exponent)

    return eta_mod


def evaluate_heating_on_posterior(
    samples: np.ndarray,
    structure_data: Dict[str, Any],
    rheology: str = 'andrade',
    n_eval: Optional[int] = None,
    random_state: int = 42,
    arrhenius_params: Optional[Dict[str, Any]] = None
) -> Tuple[np.ndarray, np.ndarray, List[Dict[str, float]]]:
    """
    Re-evaluate forward model on posterior samples to get k2 and heating.

    This is called after MCMC/SBI completes to compute heating for a subset
    of posterior samples (since heating calculation is slow).

    Args:
        samples: Posterior samples (n_samples, n_params)
        structure_data: Cached structural arrays
        rheology: 'andrade' or 'maxwell'
        n_eval: Number of samples to evaluate (default: all)
        random_state: Random seed for subsampling
        arrhenius_params: Optional Arrhenius parameters

    Returns:
        (sample_indices, k2_results, heating_results)
        - sample_indices: Indices of evaluated samples (n_eval,)
        - k2_results: Array of (Re_k2, Im_k2) for each sample (n_eval, 2)
        - heating_results: List of per-phase heating dicts (n_eval,)
    """
    n_samples = len(samples)
    if n_eval is None or n_eval > n_samples:
        n_eval = n_samples

    # Subsample posterior
    rng = np.random.RandomState(random_state)
    if n_eval < n_samples:
        idx = rng.choice(n_samples, n_eval, replace=False)
        idx.sort()
    else:
        idx = np.arange(n_samples)

    log.info(f"Evaluating {n_eval} posterior samples for heating...")

    k2_results = []
    heating_results = []

    for i, si in enumerate(idx):
        theta = samples[si]
        Re_k2, Im_k2, perPhase_W = forward_model_k2(
            theta,
            structure_data,
            rheology=rheology,
            return_heating=True,
            arrhenius_params=arrhenius_params
        )
        k2_results.append((Re_k2, Im_k2))
        heating_results.append(perPhase_W if perPhase_W is not None else {})

        if (i + 1) % 100 == 0:
            log.info(f"  {i+1}/{n_eval} samples evaluated")

    log.info(f"Heating evaluation complete: {n_eval} samples")

    return idx, np.array(k2_results), heating_results


# ============================================================================
# Induction Forward Model
# ============================================================================

def forward_model_induction(
    structure_data: Dict[str, Any],
    freq_dict: Dict[str, float],
    nn: int = 1,
    do_parallel: bool = False
) -> Optional[Dict[str, complex]]:
    """
    Compute complex induction amplitude Ae for a set of excitation frequencies.

    For Path A (discrete reference structures) the induction response depends
    only on the interior conductivity profile, not on rheological MCMC
    parameters.  Call this once per reference structure and cache the result.

    Args:
        structure_data: Cached structure dict from extract_structure_from_planet().
            Must contain 'rSigChange_m', 'sigmaLayers_Sm', 'R_body_m'.
        freq_dict: Excitation periods in hours, keyed by canonical PP label.
            Matches PlanetProfile.MagneticInduction.Moments.ExcitationsList.
            Priority order for Europa: synodic (11.23 hr), orbital (85.21 hr),
            synodic 2nd (5.62 hr), true anomaly (84.64 hr).
            Example: {'synodic': 11.23, 'orbital': 85.21, 'synodic 2nd': 5.62}
        nn: Spherical harmonic degree (1 = dipole, the dominant term).
        do_parallel: Use multiprocessing for multiple frequencies.  Default
            False; set True for batch precomputation outside the MCMC sampler.

    Returns:
        Dict[str, complex] mapping label → complex Ae, or None if induction
        data is not available in structure_data.
    """
    rSigChange_m = structure_data.get('rSigChange_m')
    sigmaLayers_Sm = structure_data.get('sigmaLayers_Sm')
    R_body_m = structure_data.get('R_body_m')

    if rSigChange_m is None or sigmaLayers_Sm is None or R_body_m is None:
        log.debug(
            "Induction data not in structure_data. "
            "Rebuild cache with CALC_CONDUCT=True."
        )
        return None

    try:
        from PlanetProfile.MagneticInduction.MoonMag.symmetry_funcs import (
            InducedAeList,
        )
    except ImportError as e:
        log.warning(f"MoonMag unavailable for induction calculation: {e}")
        return None

    labels = list(freq_dict.keys())
    T_hr = np.array([freq_dict[k] for k in labels], dtype=np.float64)
    omegas = 2.0 * np.pi / (T_hr * 3600.0)

    # rscale_moments = 1/R_body (units: 1/m) — MoonMag convention for the
    # dimensionless scaling factor applied to the outer boundary radius.
    rscale_moments = 1.0 / R_body_m

    try:
        Aes, _, _ = InducedAeList(
            rSigChange_m, sigmaLayers_Sm, omegas, rscale_moments,
            nn=nn, writeout=False, do_parallel=do_parallel,
        )
    except Exception as e:
        log.warning(f"InducedAeList failed: {e}")
        return None

    return {label: complex(ae) for label, ae in zip(labels, Aes)}


def create_log_likelihood(
    observables: Dict[str, Tuple[float, float]],
    structure_data: Dict[str, Any],
    rheology: str = 'andrade',
    arrhenius_params: Optional[Dict[str, Any]] = None
) -> callable:
    """
    Create log-likelihood function for MCMC/SBI sampling.

    Returns a function that accepts theta and returns log-likelihood value.

    Args:
        observables: Dict of {observable_name: (value, uncertainty)}
            Supported observables: 'Re_k2', 'Im_k2', 'k2', 'abs_Im_k2'
        structure_data: Cached structural arrays
        rheology: 'andrade' or 'maxwell'
        arrhenius_params: Optional Arrhenius parameters

    Returns:
        Callable log_likelihood(theta) -> float
    """
    def log_likelihood(theta):
        """Gaussian log-likelihood on tidal observables."""
        # Only compute heating if specifically requested as an observable
        # (it is ~5-10x slower than k2-only)
        needs_heating = 'total_heating_W' in observables

        Re_k2, Im_k2, perPhase_W = forward_model_k2(
            theta,
            structure_data,
            rheology=rheology,
            return_heating=needs_heating,
            arrhenius_params=arrhenius_params
        )

        if np.isnan(Re_k2):
            return -1e30  # TidalPy solver failed

        # Compute chi-squared
        chi2 = 0.0

        if 'Re_k2' in observables:
            obs_val, obs_err = observables['Re_k2']
            chi2 += ((Re_k2 - obs_val) / obs_err)**2

        if 'Im_k2' in observables:
            # Note: TidalPy returns negative Im(k2), so use abs() to match positive observable
            obs_val, obs_err = observables['Im_k2']
            chi2 += ((abs(Im_k2) - obs_val) / obs_err)**2

        if 'abs_Im_k2' in observables:
            obs_val, obs_err = observables['abs_Im_k2']
            chi2 += ((abs(Im_k2) - obs_val) / obs_err)**2

        if 'k2' in observables:
            # Magnitude constraint: |k2| = sqrt(Re^2 + Im^2)
            obs_val, obs_err = observables['k2']
            k2_mag = np.sqrt(Re_k2**2 + Im_k2**2)
            chi2 += ((k2_mag - obs_val) / obs_err)**2

        if 'total_heating_W' in observables and perPhase_W is not None:
            obs_val, obs_err = observables['total_heating_W']
            # Total is the sum of all individual phase contributions
            # (ignoring the aggregate _W keys)
            agg_keys = {'Silicate_W', 'HP_Ice_W', 'Ice_Ih_W'}
            total_W = sum(v for k, v in perPhase_W.items() if k not in agg_keys)
            chi2 += ((total_W - obs_val) / obs_err)**2

        return -0.5 * chi2

    return log_likelihood
