"""TidalPy gravity backend readiness (scientific-reviewer 2026-11-01, items
1-5 before it can become the default): libration computed from TidalPy's own
y1 + Re k2, mean-motion fallback from Bulk.Torb_s, one shared Andrade
exponent, and agreement with ALMA on an elastic Enceladus case."""
import copy
import pickle
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'tests'))
CACHE = REPO / ('PlanetProfile/Test/mcmc_results/Enceladus/Cassini_isostasy/'
                'enceladus_seawater_zbw_grid_2d.pkl')


def test_andrade_exponent_single_source():
    """Reviewer 2026-11-01: both attributes default to 0.2, so an override
    of EITHER must win (PPTest47/50 set andradExponent = 0.3; the exploreogram
    AndradeAlpha sweep sets andradAlpha); only two conflicting explicit
    overrides raise."""
    from PlanetProfile.Gravity.Gravity import AndradeExponent
    from PlanetProfile.Utilities.defineStructs import GravitySubstruct
    def g(**kw):
        G = GravitySubstruct()
        for k, v in kw.items():
            setattr(G, k, v)
        return SimpleNamespace(Gravity=G)
    assert AndradeExponent(g()) == 0.2                         # struct defaults
    assert AndradeExponent(g(andradExponent=0.3)) == 0.3       # PPTest47/50 path
    assert AndradeExponent(g(andradAlpha=0.35)) == 0.35        # exploreogram sweep
    assert AndradeExponent(g(andradAlpha=0.3, andradExponent=0.3)) == 0.3
    assert AndradeExponent(g(andradAlpha=None, andradExponent=None)) == 0.2
    with pytest.raises(ValueError):
        AndradeExponent(g(andradAlpha=0.25, andradExponent=0.3))


@pytest.mark.skipif(not CACHE.exists(), reason='cache absent')
def test_tidalpy_libration_matches_alma_on_elastic_enceladus_node():
    """Reviewer gate: <= 1 % on an elastic Enceladus case (measured ~0.2 %).
    Exercises _tidalpy_libration_m's interface-aware y1 sampling."""
    from alma_viscosity_units_test import _arrays
    from PlanetProfile.Gravity import pyALMA3Updated
    from PlanetProfile.Gravity.Gravity import (
        ALMA_TIME_UNIT_S, AlmaViscosity, LibrationModelInputs,
        _tidalpy_libration_m)
    from PlanetProfile.Gravity.Librations import librations
    from PlanetProfile.Utilities.defineStructs import Constants
    from TidalPy.RadialSolver import radial_solver
    from TidalPy.RadialSolver.helpers import build_rs_input_from_data
    from TidalPy.rheology import Elastic, Newton
    s = pickle.load(open(CACHE, 'rb'))['structures'][1220]
    r, rho, mu, eta, rheo, ph = _arrays(s, 1e14, ice='elastic')
    om = float(s['omega'])
    M, R = float(s['Mtot_kg']), float(s['R_body_m'])
    lu = np.asarray(s['layer_upper_radii'], float)
    lt = tuple(s['layer_types'])
    bd = build_rs_input_from_data(
        om, r, rho, np.asarray(s['K_Pa'], float), mu, np.zeros_like(eta), eta,
        tuple(lu), lt, tuple(False for _ in lt), tuple(False for _ in lt),
        tuple(Newton() if t == 'liquid' else Elastic() for t in lt),
        tuple(Elastic() for _ in lt), perform_checks=True, warnings=False)
    res = radial_solver(
        bd.radius_array, bd.density_array, bd.complex_bulk_modulus_array,
        bd.complex_shear_modulus_array, bd.frequency, bd.planet_bulk_density,
        bd.layer_types, bd.is_static_bylayer, bd.is_incompressible_bylayer,
        bd.upper_radius_bylayer_array, degree_l=2, solve_for=('tidal',),
        raise_on_fail=False)
    planet = SimpleNamespace(Bulk=SimpleNamespace(eccentricity=float(s['eccentricity'])),
                             Gravity=SimpleNamespace(k=complex(res.k)))
    change = np.r_[0, np.where(np.diff(ph) != 0)[0] + 1, len(ph)]
    lib_tp = _tidalpy_libration_m(planet, res, r, rho, ph, change, om)

    m = pyALMA3Updated.build_model(r, rho, mu, AlmaViscosity(eta, rheo, ph), rheo,
                                   np.zeros((r.size, 2)), norm=False)
    _, _, k, y = pyALMA3Updated.love_numbers(
        m, [2], [(2 * np.pi / om) / ALMA_TIME_UNIT_S], 'tidal', 'periodic', 8,
        verbose=False, internal=True)
    g = Constants.G * M / R ** 2
    rr, rrho, oi, ry = LibrationModelInputs(r, ph, rho, np.real(y[:, 0]) * M / R / g)
    lib_alma = librations(rr, rrho, om, float(s['eccentricity']), rigid=False,
                          ocean=True, ocean_idx=oi, y1=ry,
                          rek2=float(np.real(k[0, 0, -1])))
    assert np.isfinite(lib_tp)
    assert lib_tp == pytest.approx(lib_alma, rel=0.01)


def _run(body, backend):
    import logging
    from importlib import import_module
    from PlanetProfile.GetConfig import Params as P0
    from PlanetProfile.Main import PlanetProfile
    logging.disable(logging.WARNING)
    planet = copy.deepcopy(import_module(f'PlanetProfile.Default.{body}.PP{body}').Planet)
    params = copy.deepcopy(P0)
    params.CALC_NEW = True
    params.CALC_NEW_GRAVITY = True
    params.CALC_NEW_INDUCT = False
    params.NO_SAVEFILE = True
    params.SKIP_PLOTS = True
    params.QUIET_ALMA = True
    params.Gravity.backend = backend
    planet, params = PlanetProfile(planet, params)
    logging.disable(logging.NOTSET)
    return planet


@pytest.mark.slow
def test_titan_tidalpy_backend_end_to_end():
    """Default Titan sets no meanMotion_radps: the backend must fall back to
    Bulk.Torb_s, return a finite libration and per-phase heating, and agree
    with ALMA on Re k2 within the solvers' known 2-7 % offset."""
    tp = _run('Titan', 'tidalpy')
    alma = _run('Titan', 'pyalma')
    assert getattr(tp.Bulk, 'meanMotion_radps', None) is None
    assert np.isfinite(tp.Gravity.libration_m) and tp.Gravity.libration_m > 0
    assert tp.Gravity.tidalpy_Htidal_perPhase_W
    assert tp.Gravity.k.real == pytest.approx(alma.Gravity.k.real, rel=0.07)
    assert np.isfinite(tp.Gravity.kAmp) and np.isfinite(tp.Gravity.deltaAmp)
