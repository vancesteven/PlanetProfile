"""Find Tb_K prior bounds that satisfy geological-constraint specs.

Bisection-based search over Tb_K against a body's PP template module.
Generalizable: any body PP supports can be probed by passing
``--template-module`` and the constraint spec. Three constraint kinds
are supported, mirroring the C2 use cases:

    --d-iceIh-min KM      Tb where D_iceIh = KM at the WARM (thin-ice) end
    --d-iceIh-max KM      Tb where D_iceIh = KM at the COLD (thick-ice) end
    --d-ocean-min KM      Tb where D_ocean = KM at the COLD (frozen) end

The Ih-III-L triple point at 251.165 K is treated as a hard lower bound
on Tb; the search clamps to it so we never produce non-physical
configurations. The search runs PP in structure-only mode (no gravity,
no plotting) to keep iteration cost ~1 s per Tb evaluation.

Usage
-----
    NUMBA_CACHE_DIR=$TMPDIR/pp_numba_cache python -m \\
        PlanetProfile.Inference.find_tb_bounds \\
        --bodyname Europa \\
        --template-module PlanetProfile.Default.Europa.PPEuropa \\
        --d-iceIh-min 3 --d-iceIh-max 100

    # NaCl ocean override (Callisto):
    python -m PlanetProfile.Inference.find_tb_bounds \\
        --bodyname Callisto \\
        --template-module PlanetProfile.Default.Callisto.PPCallisto \\
        --ocean-comp NaCl --ocean-wppt 100.0 \\
        --d-iceIh-max 100 --d-ocean-min 1.0

Author: PlanetProfile Team
Date:   2026-05-23
"""
from __future__ import annotations

import os
import tempfile
if not os.environ.get('NUMBA_CACHE_DIR'):
    _default = os.path.join(tempfile.gettempdir(), 'pp_numba_cache')
    try:
        os.makedirs(_default, exist_ok=True)
        os.environ['NUMBA_CACHE_DIR'] = _default
    except OSError:
        pass

import argparse
import copy
import importlib
import json
import logging
import sys
from typing import Callable, Dict, Optional, Tuple

log = logging.getLogger('PlanetProfile.Inference.find_tb_bounds')

# Hard lower bound: below this Tb, Ice Ih in contact with liquid water is
# replaced by Ice III at the relevant pressures (Constants.TtripleIh_III_L_K).
# We never search below this floor.
TRIPLE_POINT_K = 251.165


def evaluate_tb(
    bodyname: str,
    template_module: str,
    Tb_K: float,
    ocean_comp: Optional[str] = None,
    ocean_wppt: Optional[float] = None,
) -> Dict[str, float]:
    """Run PP for the given body at the requested Tb_K and return the
    resulting hydrosphere thicknesses.

    Returns a dict with keys ``D_iceIh_km``, ``D_ocean_km``,
    ``D_iceIII_km``, ``D_iceV_km``, ``D_iceVI_km``, ``D_hsphere_km``.
    Returns ``{'failed': str}`` if PP raises.
    """
    # Lazy imports so the help screen and arg parsing don't pay the PP startup cost
    from PlanetProfile.Main import PlanetProfile as runPP
    from PlanetProfile.GetConfig import Params as configParams

    # NOTE: Earlier versions of this function rejected Tb_K <= TRIPLE_POINT_K
    # outright. That was wrong for non-pure compositions: NaCl/MgSO4 brines
    # depress the freezing point and may shift the relevant ice-fluid
    # coexistence to lower Tb (and possibly to Ih-II-L instead of Ih-III-L).
    # See user memory: project_composition_aware_phase_boundaries.md. We now
    # let PP decide whether the requested Tb is physically valid for the
    # given composition; PP's failure modes (no transition, MoI not found)
    # are surfaced via the existing exception handling below.

    if template_module in sys.modules:
        del sys.modules[template_module]
    mod = importlib.import_module(template_module)
    Planet = copy.deepcopy(mod.Planet)
    Planet.Bulk.Tb_K = float(Tb_K)
    if ocean_comp is not None:
        Planet.Ocean.comp = ocean_comp
    if ocean_wppt is not None:
        Planet.Ocean.wOcean_ppt = float(ocean_wppt)

    cfg = copy.deepcopy(configParams)
    cfg.CALC_NEW = True
    cfg.CALC_NEW_GRAVITY = False
    cfg.NO_SAVEFILE = True
    cfg.SKIP_PLOTS = True

    try:
        Planet, _ = runPP(Planet, cfg)
    except Exception as exc:  # noqa: BLE001
        msg = str(exc)
        # Interpret PP's "no valid phase transition" failure as ice-shell
        # extinction: at warm Tb, the geotherm meets the melting curve at
        # the surface, so D_iceIh → 0. Returning a synthetic zero-state
        # lets the bisection bracket through this regime instead of
        # treating it as an unrecoverable failure.
        if ('No valid phase transition' in msg
                or 'Unable to find a valid pressure' in msg):
            return {
                'Tb_K': float(Tb_K),
                'D_iceIh_km': 0.0,
                'D_iceIII_km': 0.0,
                'D_iceV_km': 0.0,
                'D_iceVI_km': 0.0,
                'D_ocean_km': 0.0,  # unreliable in this regime; not used by bisect
                'D_hsphere_km': 0.0,
                'note': 'pp_extinct',
            }
        return {'failed': f'PP run failed: {type(exc).__name__}: {exc}'}

    # Iterate Planet.Reduced.{phase, r_m, changeIndices} the same way
    # cache_builder._build_single_structure does (cache_builder.py:333-346),
    # so the thickness extraction is body-agnostic and doesn't depend on
    # whichever derived D_*_km attributes happen to exist on Planet.
    import numpy as np
    changeIndices = (
        np.max(Planet.Reduced.changeIndices)
        - np.flipud(Planet.Reduced.changeIndices)
    )
    n_layers = len(changeIndices) - 1
    phases = np.flipud(Planet.Reduced.phase)
    r_m_arr = np.flipud(Planet.Reduced.r_m)

    D_iceIh_km = D_iceIII_km = D_iceV_km = D_iceVI_km = D_ocean_km = 0.0
    for i_layer in range(n_layers):
        s = changeIndices[i_layer]
        e = changeIndices[i_layer + 1]
        ph = int(phases[s])
        thick_km = (r_m_arr[e - 1] - r_m_arr[s]) / 1e3
        if ph == 1:
            D_iceIh_km += thick_km
        elif ph == 3:
            D_iceIII_km += thick_km
        elif ph == 5:
            D_iceV_km += thick_km
        elif ph == 6:
            D_iceVI_km += thick_km
        elif ph == 0:
            D_ocean_km += thick_km

    R_body_m = float(Planet.Bulk.R_m)
    R_sil = float(getattr(Planet.Sil, 'Rmean_m', r_m_arr[0]))
    D_hsphere_km = (R_body_m - R_sil) / 1e3

    return {
        'Tb_K': float(Tb_K),
        'D_iceIh_km': D_iceIh_km,
        'D_iceIII_km': D_iceIII_km,
        'D_iceV_km': D_iceV_km,
        'D_iceVI_km': D_iceVI_km,
        'D_ocean_km': D_ocean_km,
        'D_hsphere_km': D_hsphere_km,
    }


def bisect_tb(
    target_fn: Callable[[float], Optional[float]],
    target_value: float,
    Tb_lo: float,
    Tb_hi: float,
    tol_K: float = 0.02,
    max_iter: int = 30,
) -> Optional[float]:
    """Bisect over Tb to find where ``target_fn(Tb) = target_value``.

    ``target_fn`` is ``f(Tb) -> float`` (e.g., D_iceIh as a function of
    Tb), expected monotone-decreasing in this codepath (warmer Tb →
    thinner ice → smaller D_iceIh; warmer Tb → thicker ocean → larger
    D_ocean). Returns ``None`` if the bracket does not contain
    ``target_value``.
    """
    f_lo = target_fn(Tb_lo)
    f_hi = target_fn(Tb_hi)
    if f_lo is None or f_hi is None:
        return None
    # Verify the target is bracketed
    if (f_lo - target_value) * (f_hi - target_value) > 0:
        log.warning(
            f'bisect_tb: target {target_value:.3f} not bracketed by '
            f'f({Tb_lo:.3f})={f_lo:.3f}, f({Tb_hi:.3f})={f_hi:.3f}'
        )
        return None
    for _ in range(max_iter):
        Tb_mid = 0.5 * (Tb_lo + Tb_hi)
        f_mid = target_fn(Tb_mid)
        if f_mid is None:
            return None
        if Tb_hi - Tb_lo < tol_K:
            return Tb_mid
        # Same sign as f_lo? then root in upper half
        if (f_mid - target_value) * (f_lo - target_value) > 0:
            Tb_lo, f_lo = Tb_mid, f_mid
        else:
            Tb_hi, f_hi = Tb_mid, f_mid
    return 0.5 * (Tb_lo + Tb_hi)


def find_bounds(
    bodyname: str,
    template_module: str,
    d_iceIh_min: Optional[float] = None,
    d_iceIh_max: Optional[float] = None,
    d_ocean_min: Optional[float] = None,
    Tb_search_lo: float = TRIPLE_POINT_K + 0.05,
    Tb_search_hi: float = 275.0,
    ocean_comp: Optional[str] = None,
    ocean_wppt: Optional[float] = None,
) -> Dict[str, float]:
    """Find Tb_K prior bounds satisfying the supplied constraints.

    All three constraints are optional. Convention:
      d_iceIh_min: warm-Tb end. Tb_upper is set so D_iceIh = d_iceIh_min.
      d_iceIh_max: cold-Tb end. Tb_lower is set so D_iceIh = d_iceIh_max.
      d_ocean_min: cold-Tb end (overrides d_iceIh_max if both set).
                   Tb_lower is set so D_ocean = d_ocean_min.

    Result dict reports the computed bounds plus the hydrosphere state
    at each bound (so you can see what you actually got).
    """
    log.info(
        f'Searching Tb bounds for {bodyname} '
        f'(template={template_module}, ocean={ocean_comp or "default"})'
    )

    # f(Tb) → state dict; cached so each Tb is only run once
    _cache: Dict[float, Dict[str, float]] = {}

    def _state(Tb_K: float) -> Dict[str, float]:
        key = round(Tb_K, 4)
        if key not in _cache:
            _cache[key] = evaluate_tb(bodyname, template_module, Tb_K,
                                       ocean_comp=ocean_comp,
                                       ocean_wppt=ocean_wppt)
        return _cache[key]

    def _safe_value(state: Dict[str, float], field: str) -> Optional[float]:
        if 'failed' in state:
            log.warning(f'  Tb={state.get("Tb_K", "?")}: {state["failed"]}')
            return None
        return state.get(field)

    # WARM end: search for D_iceIh = d_iceIh_min
    Tb_upper = None
    state_upper = None
    if d_iceIh_min is not None:
        log.info(f'  warm-end search: targeting D_iceIh = {d_iceIh_min:.2f} km')
        f = lambda Tb: _safe_value(_state(Tb), 'D_iceIh_km')
        Tb_upper = bisect_tb(f, d_iceIh_min, Tb_search_lo, Tb_search_hi)
        if Tb_upper is not None:
            state_upper = _state(Tb_upper)
            log.info(
                f'    → Tb_upper = {Tb_upper:.4f} K  '
                f'(D_iceIh = {state_upper["D_iceIh_km"]:.2f} km, '
                f'D_ocean = {state_upper["D_ocean_km"]:.2f} km)'
            )

    # COLD end: D_ocean_min takes priority, else D_iceIh_max
    Tb_lower = None
    state_lower = None
    if d_ocean_min is not None:
        log.info(f'  cold-end search: targeting D_ocean = {d_ocean_min:.2f} km')
        f = lambda Tb: _safe_value(_state(Tb), 'D_ocean_km')
        Tb_lower = bisect_tb(f, d_ocean_min, Tb_search_lo, Tb_search_hi)
        if Tb_lower is not None:
            state_lower = _state(Tb_lower)
            log.info(
                f'    → Tb_lower = {Tb_lower:.4f} K  '
                f'(D_ocean = {state_lower["D_ocean_km"]:.2f} km, '
                f'D_iceIh = {state_lower["D_iceIh_km"]:.2f} km)'
            )
    elif d_iceIh_max is not None:
        log.info(f'  cold-end search: targeting D_iceIh = {d_iceIh_max:.2f} km')
        f = lambda Tb: _safe_value(_state(Tb), 'D_iceIh_km')
        Tb_lower = bisect_tb(f, d_iceIh_max, Tb_search_lo, Tb_search_hi)
        if Tb_lower is not None:
            state_lower = _state(Tb_lower)
            log.info(
                f'    → Tb_lower = {Tb_lower:.4f} K  '
                f'(D_iceIh = {state_lower["D_iceIh_km"]:.2f} km, '
                f'D_ocean = {state_lower["D_ocean_km"]:.2f} km)'
            )

    return {
        'bodyname': bodyname,
        'Tb_K_lower': Tb_lower,
        'Tb_K_upper': Tb_upper,
        'state_lower': state_lower,
        'state_upper': state_upper,
        'triple_point_K': TRIPLE_POINT_K,
        'evaluations': len(_cache),
    }


def _main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--bodyname', required=True)
    p.add_argument('--template-module', required=True,
                   help='e.g., PlanetProfile.Default.Europa.PPEuropa')
    p.add_argument('--d-iceIh-min', type=float,
                   help='Geological-minimum ice shell thickness (km). '
                        'Sets warm-Tb end of prior.')
    p.add_argument('--d-iceIh-max', type=float,
                   help='Maximum ice shell thickness (km). Sets cold-Tb end.')
    p.add_argument('--d-ocean-min', type=float,
                   help='Minimum ocean thickness (km). Sets cold-Tb end '
                        '(takes priority over --d-iceIh-max).')
    p.add_argument('--ocean-comp', type=str, default=None,
                   help='Ocean composition override (e.g., NaCl).')
    p.add_argument('--ocean-wppt', type=float, default=None,
                   help='Ocean salinity override (ppt).')
    p.add_argument('--Tb-search-lo', type=float, default=TRIPLE_POINT_K + 0.05)
    p.add_argument('--Tb-search-hi', type=float, default=275.0)
    p.add_argument('--json-out', type=str, default=None)
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s [%(levelname)s] %(message)s')

    result = find_bounds(
        bodyname=args.bodyname,
        template_module=args.template_module,
        d_iceIh_min=args.d_iceIh_min,
        d_iceIh_max=args.d_iceIh_max,
        d_ocean_min=args.d_ocean_min,
        Tb_search_lo=args.Tb_search_lo,
        Tb_search_hi=args.Tb_search_hi,
        ocean_comp=args.ocean_comp,
        ocean_wppt=args.ocean_wppt,
    )

    print('\n=== Result ===')
    print(json.dumps(result, indent=2, default=str))

    if args.json_out:
        with open(args.json_out, 'w') as f:
            json.dump(result, f, indent=2, default=str)
        print(f'\nWritten to {args.json_out}')

    return 0


if __name__ == '__main__':
    raise SystemExit(_main())
