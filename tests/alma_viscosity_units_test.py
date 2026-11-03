"""PlanetProfile's ALMA gravity path must give a viscosity-dependent k2
(scientific-reviewer 2026-10-31: with eta in Pa*s and the period in kyr every
Maxwell/Andrade layer was effectively elastic, k2 invariant from 1e10 to 1e22
Pa*s). Checks Gravity.AlmaViscosity on a production Enceladus node:
non-invariance, the elastic limit, agreement with TidalPy's independent radial
solver, and ocean-viscosity convergence (the Newtonian ocean must NOT be
converted)."""
import pickle
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
CACHE = REPO / ('PlanetProfile/Test/mcmc_results/Enceladus/Cassini_isostasy/'
                'enceladus_seawater_zbw_grid_2d.pkl')
pytestmark = pytest.mark.skipif(not CACHE.exists(), reason='cache absent')


@pytest.fixture(scope='module')
def node():
    return pickle.load(open(CACHE, 'rb'))['structures'][1220]


def _arrays(s, eta_ice, ice='maxwell', eta_ocean=None):
    r = np.asarray(s['r_m'], float)
    ph = np.asarray(s['phases'])
    mu = np.array(s['mu_Pa'], float)
    eta = np.array(s['eta_Pa_base'], float)
    mu[ph == 1] = 3.3e9
    eta[ph == 1] = eta_ice
    if eta_ocean is not None:
        eta[ph == 0] = eta_ocean
    rheo = np.where(ph == 0, 'newton', np.where(ph == 1, ice, 'elastic')).tolist()
    return r, np.asarray(s['rho'], float), mu, eta, rheo, ph


def _alma_k2(s, eta_ice, ice='maxwell', eta_ocean=None, convert=True):
    from PlanetProfile.Gravity import pyALMA3Updated
    from PlanetProfile.Gravity.Gravity import ALMA_TIME_UNIT_S, AlmaViscosity
    r, rho, mu, eta, rheo, ph = _arrays(s, eta_ice, ice, eta_ocean)
    params = np.zeros((r.size, 2))
    if ice == 'andrade':
        params[:, 0] = 0.3
    vis = AlmaViscosity(eta, rheo, ph) if convert else eta
    m = pyALMA3Updated.build_model(r, rho, mu, vis, rheo, params, norm=False)
    per = (2 * np.pi / float(s['omega'])) / ALMA_TIME_UNIT_S
    _, _, k, _ = pyALMA3Updated.love_numbers(m, [2], [per], 'tidal', 'periodic',
                                             8, verbose=False, internal=True)
    return complex(k[0, 0, -1])


def _tidalpy_k2(s, eta_ice):
    from TidalPy.RadialSolver import radial_solver
    from TidalPy.RadialSolver.helpers import build_rs_input_from_data
    from TidalPy.rheology import Elastic, Maxwell, Newton
    r, rho, mu, eta, _, ph = _arrays(s, eta_ice)
    lu = np.asarray(s['layer_upper_radii'], float)
    lt = tuple(s['layer_types'])
    shear = tuple(Newton() if t == 'liquid' else
                  (Maxwell() if ph[r <= ub][-1] == 1 else Elastic())
                  for t, ub in zip(lt, lu))
    bd = build_rs_input_from_data(
        float(s['omega']), r, rho, np.asarray(s['K_Pa'], float), mu,
        np.zeros_like(eta), eta, tuple(lu), lt, tuple(False for _ in lt),
        tuple(False for _ in lt), shear, tuple(Elastic() for _ in lt),
        perform_checks=True, warnings=False)
    res = radial_solver(
        bd.radius_array, bd.density_array, bd.complex_bulk_modulus_array,
        bd.complex_shear_modulus_array, bd.frequency, bd.planet_bulk_density,
        bd.layer_types, bd.is_static_bylayer, bd.is_incompressible_bylayer,
        bd.upper_radius_bylayer_array, degree_l=2, solve_for=('tidal',),
        raise_on_fail=False)
    return complex(res.k)


def test_k2_depends_on_ice_viscosity(node):
    k_soft = _alma_k2(node, 1e13)
    k_stiff = _alma_k2(node, 1e20)
    assert abs(k_soft) > 3 * abs(k_stiff)
    assert -k_soft.imag > 100 * -k_stiff.imag


def test_unconverted_path_reproduces_the_defect(node):
    """Guards the diagnosis: without conversion k2 is eta-invariant."""
    a = _alma_k2(node, 1e13, convert=False)
    b = _alma_k2(node, 1e20, convert=False)
    # Re k2 is eta-invariant to <1e-9 on any platform. Im k2 is ~3e-5 here
    # and dominated by ALMA's ocean regularization, which differs by ~6e-4
    # relative between macOS and Linux (CI 2026-11-03), so compare it
    # absolutely: the converted path moves Im k2 by ~8e-2 at 1e13 Pa s, four
    # orders of magnitude above this tolerance.
    assert a.real == pytest.approx(b.real, rel=1e-6)
    assert a.imag == pytest.approx(b.imag, abs=1e-7)


def test_stiff_maxwell_equals_elastic(node):
    assert _alma_k2(node, 1e22) == pytest.approx(
        _alma_k2(node, 1e22, ice='elastic'), rel=1e-3, abs=1e-6)


@pytest.mark.parametrize('eta', [1e13, 1e14, 3.9e14])
def test_imag_k2_matches_tidalpy(eta, node):
    """Independent solver. Im k2 (dissipation) agrees to 5 %; Re k2 carries
    a ~6-7 % ALMA-vs-TidalPy offset that is present even for an elastic
    shell (compressibility/bulk-modulus treatment), so it is not a
    viscosity-unit effect and is checked separately."""
    ka, kt = _alma_k2(node, eta), _tidalpy_k2(node, eta)
    assert ka.imag == pytest.approx(kt.imag, rel=0.05)


def test_re_k2_offset_is_viscosity_independent(node):
    off = [(_alma_k2(node, e).real / _tidalpy_k2(node, e).real)
           for e in (1e13, 1e14, 1e20)]
    assert max(off) - min(off) < 0.05
    assert all(0.85 < o < 1.05 for o in off)


def test_ocean_viscosity_convergence(node):
    """The Newtonian ocean is left unconverted. Its eta*s (s in kyr^-1) is a
    ~kPa regularizing rigidity, not physics: k2 sits on a plateau for ocean
    eta 1e-3..1e-5 (Re/Im within 0.5 %), degrades above ~1e-2 (spurious ocean
    dissipation) and collapses below ~1e-7 (propagator ill-conditioning). The
    cache default 1.9e-3 is on the plateau."""
    ref = _alma_k2(node, 1e14, eta_ocean=1e-4)
    for e in (1.9e-3, 1e-3, 1e-5):
        k = _alma_k2(node, 1e14, eta_ocean=e)
        assert k.real == pytest.approx(ref.real, rel=5e-3)
        assert k.imag == pytest.approx(ref.imag, rel=5e-3)


def test_andrade_is_converted_too(node):
    soft = _alma_k2(node, 1e14, ice='andrade')
    stiff = _alma_k2(node, 1e22, ice='andrade')
    assert -soft.imag > 10 * -stiff.imag


def test_solid_newtonian_clathrate_is_converted_ocean_is_not():
    """Reviewer 2026-11-01 MAJOR: the default 'Clath': 'newton' is a SOLID
    creeping layer and must be converted; the phase-0 ocean must not be."""
    from PlanetProfile.Gravity.Gravity import ALMA_TIME_UNIT_S, AlmaViscosity
    eta = np.array([1e20, 1e-3, 1e15, 1e14])
    rheo = ['elastic', 'newton', 'newton', 'maxwell']
    phase = np.array([50, 0, 30, 1])
    out = AlmaViscosity(eta, rheo, phase)
    assert out[0] == 1e20 and out[1] == 1e-3
    assert out[2] == pytest.approx(1e15 / ALMA_TIME_UNIT_S)
    assert out[3] == pytest.approx(1e14 / ALMA_TIME_UNIT_S)


def test_integer_rheology_codes_are_rejected():
    from PlanetProfile.Gravity.Gravity import AlmaViscosity
    with pytest.raises(TypeError):
        AlmaViscosity(np.ones(2), [1, 2], np.array([50, 1]))


@pytest.mark.xfail(strict=True, reason=(
    'ALMA ocean regularization: an inviscid ocean between elastic layers '
    'should give Im k2 = 0 (TidalPy does), ALMA gives ~-3e-5 here and '
    '-0.0029 for default Europa; reviewer 2026-11-01. Quote Im k2 from '
    'TidalPy.'))
def test_all_elastic_with_ocean_has_no_dissipation(node):
    assert abs(_alma_k2(node, 1e14, ice='elastic').imag) < 1e-6
