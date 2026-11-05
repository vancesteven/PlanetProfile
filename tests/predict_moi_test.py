"""Tests for Do.PREDICT_MOI (LayerPropagators.PredictMoI).

With the flag, the C/MR^2 of the single mass-closed candidate of a waterless body with an iron
core is reported as a prediction and compared with a range; it no longer selects or invalidates
the structure. Uses a minimal stand-in for Planet -- no table or full run required.
"""
from types import SimpleNamespace

import numpy as np
import pytest

from PlanetProfile.Thermodynamics.LayerPropagators import PredictMoI


def _planet(NO_H2O=True, Fe_CORE=True, comparison=None):
    return SimpleNamespace(Do=SimpleNamespace(NO_H2O=NO_H2O, Fe_CORE=Fe_CORE),
                           Bulk=SimpleNamespace(Cmeasured=0.3931, CuncertaintyLower=0.0002,
                                                CuncertaintyUpper=0.0002, CMR2comparison=comparison),
                           CMR2predicted=None, CMR2withinComparison=None)


def test_single_candidate_is_kept_and_reported():
    P = _planet(comparison=(0.39311, 0.39338))
    inds = PredictMoI(P, np.array([0.393165, 0.0]), np.array([0]))
    assert inds == [0]
    assert P.CMR2predicted == pytest.approx(0.393165)
    assert P.CMR2withinComparison is True


def test_value_outside_range_is_kept_but_flagged():
    P = _planet(comparison=(0.39311, 0.39338))
    inds = PredictMoI(P, np.array([0.393080]), np.array([0]))
    assert inds == [0]
    assert P.CMR2withinComparison is False


def test_default_comparison_is_the_cmeasured_band():
    P = _planet()
    PredictMoI(P, np.array([0.39335]), np.array([0]))
    assert P.CMR2withinComparison is False  # outside 0.3931 +/- 0.0002
    P = _planet()
    PredictMoI(P, np.array([0.39320]), np.array([0]))
    assert P.CMR2withinComparison is True


def test_no_valid_candidate_returns_empty():
    P = _planet()
    assert PredictMoI(P, np.array([0.0]), np.array([0])) == []
    assert np.isnan(P.CMR2predicted) and P.CMR2withinComparison is False


@pytest.mark.parametrize('NO_H2O,Fe_CORE', [(False, True), (True, False)])
def test_requires_waterless_body_with_iron_core(NO_H2O, Fe_CORE):
    with pytest.raises(ValueError):
        PredictMoI(_planet(NO_H2O, Fe_CORE), np.array([0.3931]), np.array([0]))
