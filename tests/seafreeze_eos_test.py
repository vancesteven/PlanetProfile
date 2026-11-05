"""SeaFreeze EOS: direct phase finder and per-property table builds.

* ``SFphase`` (``Ocean.phaseType = 'calc'``) returns a scalar for scalar
  input, as the lookup-table phase finder does. ``GetPfreeze``'s brentq
  multiplies two results, and a 1-element array there crashed default
  Ganymede with ``phaseType = 'calc'``.
* The EOS tables request only the SeaFreeze properties they use
  (``SF_OCEAN_PROPS``, ``SF_ICE_PROPS``). Those values must be identical
  to the all-property call.

Run: python -m pytest tests/seafreeze_eos_test.py
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import seafreeze.seafreeze as sfz  # noqa: E402

from PlanetProfile.Thermodynamics.HydroEOS import (  # noqa: E402
    SF_OCEAN_PROPS, SFphase, sfPTgrid, sfPTmGrid)
from PlanetProfile.Thermodynamics.Seafreeze.SeafreezeProps import (  # noqa: E402
    SF_ICE_PROPS, IceSeaFreezeProps)
from PlanetProfile.Utilities.defineStructs import EOSlist  # noqa: E402
from PlanetProfile.Utilities.DataManip import Nearest2DInterpolator  # noqa: E402


# --- direct phase finder ---------------------------------------------------

def test_scalar_input_returns_scalar():
    finder = SFphase(0.0, 'PureH2O')
    for P in (50.0, np.float64(50.0), np.array(50.0)):
        out = finder(P, 260.0)
        assert np.ndim(out) == 0 and int(out) == 1      # ice Ih
        # brentq-style use: the product of two results is a Python scalar
        assert float(finder(P, 260.0)) * float(finder(300.0, 260.0)) >= 0


def test_vector_input_unchanged():
    finder = SFphase(0.0, 'PureH2O')
    out = finder(np.array([50.0, 300.0, 300.0, 400.0]),
                 np.array([260.0, 260.0, 250.0, 250.0]))
    # SeaFreeze water1: Ih; liquid (300 MPa is above Ih's limit and the
    # liquidus there is below 260 K); ice III; ice V.
    assert out.shape == (4,) and list(out) == [1, 0, 3, 5]


def test_matches_lookup_away_from_boundaries():
    Pg = np.arange(0.1, 600, 1.0)
    Tg = np.arange(240, 300, 0.5)
    lookup = Nearest2DInterpolator(
        Pg, Tg, sfz.whichphase(sfPTgrid(Pg, Tg)).astype(np.int8))
    rng = np.random.default_rng(0)
    P, T = rng.uniform(1, 590, 2000), rng.uniform(241, 299, 2000)
    finder = SFphase(0.0, 'PureH2O')
    direct = finder(P, T)
    # Disagreement only within half a grid step of a phase boundary.
    bad = np.where(lookup(P, T) != direct)[0]
    for i in bad:
        alongP = finder(np.array([P[i] - 0.5, P[i] + 0.5]), np.array([T[i]] * 2))
        alongT = finder(np.array([P[i]] * 2), np.array([T[i] - 0.25, T[i] + 0.25]))
        assert alongP[0] != alongP[1] or alongT[0] != alongT[1]
    assert len(bad) / len(P) < 0.01


# --- per-property table builds ----------------------------------------------

_OCEAN_P = np.linspace(0.1, 800, 60)
_OCEAN_T = np.linspace(240, 350, 50)
_CASES = [
    ('water1', sfPTgrid(_OCEAN_P, _OCEAN_T), SF_OCEAN_PROPS),
    ('NaClaq', sfPTmGrid(_OCEAN_P, _OCEAN_T, 0.5), SF_OCEAN_PROPS),
    ('Ih', sfPTgrid(np.linspace(0.1, 200, 40), np.linspace(100, 272, 40)),
     SF_ICE_PROPS),
] + [(phase, sfPTgrid(np.linspace(200, 900, 40), np.linspace(150, 350, 40)),
      SF_ICE_PROPS) for phase in ('II', 'III', 'V', 'VI')]


@pytest.mark.parametrize('phase, PT, props', _CASES,
                         ids=[case[0] for case in _CASES])
def test_requested_properties_equal_the_all_property_call(phase, PT, props):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        full = sfz.getProp(PT.copy(), phase)
        part = sfz.getProp(PT.copy(), phase, sfz.defpath, *props)
    for name in props:
        np.testing.assert_array_equal(
            np.asarray(getattr(part, name)), np.asarray(getattr(full, name)),
            err_msg=f'{phase}.{name}')


def test_ice_cache_reuse_returns_the_cached_grid():
    """A coarser request inside an already-built finer grid reuses the
    cached table; the P and T it returns must match that table's shape,
    now read from rho (G is no longer computed)."""
    EOSlist.loaded.pop('II', None)
    EOSlist.ranges.pop('II', None)
    fine = sfPTgrid(np.linspace(200, 600, 81), np.linspace(150, 250, 51))
    out_fine, P_fine, T_fine = IceSeaFreezeProps(fine, 'II')
    coarse = sfPTgrid(np.linspace(250, 550, 31), np.linspace(160, 240, 21))
    out, P, T = IceSeaFreezeProps(coarse, 'II')
    assert out is out_fine
    assert (P.size, T.size) == out.rho.shape == (81, 51)
    np.testing.assert_array_equal(P, P_fine)
    np.testing.assert_array_equal(T, T_fine)
    assert not hasattr(out, 'G')
