"""Direct SeaFreeze phase finder (Ocean.phaseType = 'calc', HydroEOS.SFphase):
scalar in -> scalar out (GetPfreeze's brentq multiplies two results), and
agreement with the lookup table away from boundaries."""
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))


def test_scalar_input_returns_scalar():
    from PlanetProfile.Thermodynamics.HydroEOS import SFphase
    f = SFphase(0.0, 'PureH2O')
    for P in (50.0, np.float64(50.0), np.array(50.0)):
        out = f(P, 260.0)
        assert np.ndim(out) == 0 and int(out) == 1      # ice Ih
        assert float(f(P, 260.0)) * float(f(300.0, 260.0)) >= 0  # brentq-style use


def test_vector_input_unchanged():
    from PlanetProfile.Thermodynamics.HydroEOS import SFphase
    f = SFphase(0.0, 'PureH2O')
    out = f(np.array([50.0, 300.0, 300.0, 400.0]),
            np.array([260.0, 260.0, 250.0, 250.0]))
    # SeaFreeze water1: Ih; liquid (300 MPa is above Ih's limit and the
    # liquidus there is below 260 K); ice III; ice V.
    assert out.shape == (4,) and list(out) == [1, 0, 3, 5]


def test_matches_lookup_away_from_boundaries():
    import seafreeze.seafreeze as sf
    from PlanetProfile.Thermodynamics.HydroEOS import SFphase, sfPTgrid
    from PlanetProfile.Utilities.DataManip import Nearest2DInterpolator
    Pg = np.arange(0.1, 600, 1.0)
    Tg = np.arange(240, 300, 0.5)
    lut = Nearest2DInterpolator(Pg, Tg, sf.whichphase(sfPTgrid(Pg, Tg)).astype(np.int8))
    rng = np.random.default_rng(0)
    P, T = rng.uniform(1, 590, 2000), rng.uniform(241, 299, 2000)
    direct = SFphase(0.0, 'PureH2O')(P, T)
    # disagreement only within half a grid step of a phase boundary
    bad = np.where(lut(P, T) != direct)[0]
    for i in bad:
        nb = SFphase(0.0, 'PureH2O')(np.array([P[i] - 0.5, P[i] + 0.5]), np.array([T[i]] * 2))
        nbT = SFphase(0.0, 'PureH2O')(np.array([P[i]] * 2), np.array([T[i] - 0.25, T[i] + 0.25]))
        assert nb[0] != nb[1] or nbT[0] != nbT[1]
    assert len(bad) / len(P) < 0.01
