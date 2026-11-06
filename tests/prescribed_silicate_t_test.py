"""Tests for the externally prescribed silicate temperature mode.

The mode must replace ONLY the silicate temperature assignment: the
structural recursion and the core thermodynamics stay untouched, and no
query may silently extrapolate beyond the supplied thermal data.
"""
import os

import numpy as np
import pytest

from PlanetProfile import _ROOT

from PlanetProfile.Thermodynamics.PrescribedT import (
    LoadPrescribedTProfileCSV,
    PrescribedTemperatureProfile,
    PrescribedTemperatureRangeError,
)

CSV = os.path.join(_ROOT, 'Default', 'Luna', 'santangelo_preferred_radial_temperature.csv')


@pytest.fixture(scope='module')
def prof():
    return LoadPrescribedTProfileCSV(CSV, source='Santangelo preferred, global mean')


# --- construction / validation -----------------------------------------

def test_rejects_non_monotonic_depth():
    with pytest.raises(ValueError, match='strictly increasing'):
        PrescribedTemperatureProfile([0.0, 100.0, 50.0], [250.0, 300.0, 400.0])


def test_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        PrescribedTemperatureProfile([0.0, 100.0], [250.0])


def test_rejects_nonfinite():
    with pytest.raises(ValueError, match='non-finite'):
        PrescribedTemperatureProfile([0.0, 100.0], [250.0, np.nan])


def test_rejects_spline_request():
    """Linear is deliberate: splines ring across sharp thermal gradients."""
    with pytest.raises(ValueError, match='only "linear"'):
        PrescribedTemperatureProfile([0.0, 1.0], [1.0, 2.0], method='cubic')


# --- exact node recovery ------------------------------------------------

def test_exact_recovery_at_every_source_node(prof):
    got = prof(prof.depth_m)
    assert np.allclose(got, prof.T_K, rtol=0, atol=1e-12), \
        'prescribed T must reproduce source values exactly at tabulated depths'


def test_surface_value_matches_source(prof):
    assert float(prof(0.0)) == pytest.approx(250.0, abs=1e-9)


# --- interpolation ------------------------------------------------------

def test_midpoint_is_linear_interpolation(prof):
    i = 20
    z0, z1 = prof.depth_m[i], prof.depth_m[i + 1]
    T0, T1 = prof.T_K[i], prof.T_K[i + 1]
    zm = 0.5 * (z0 + z1)
    assert float(prof(zm)) == pytest.approx(0.5 * (T0 + T1), rel=1e-12)


def test_interpolation_is_monotone_between_increasing_nodes(prof):
    z = np.linspace(prof.zMin_m, prof.zMax_m, 5000)
    T = prof(z)
    # The lunar geotherm is monotonically increasing with depth; linear
    # interpolation must not introduce overshoot.
    assert np.all(np.diff(T) >= -1e-9)
    assert T.min() >= prof.TMin_K - 1e-9
    assert T.max() <= prof.TMax_K + 1e-9


# --- range enforcement: no silent extrapolation -------------------------

def test_below_range_raises(prof):
    with pytest.raises(PrescribedTemperatureRangeError, match='outside its data range'):
        prof(-1.0e3)


def test_above_range_raises(prof):
    with pytest.raises(PrescribedTemperatureRangeError, match='outside its data range'):
        prof(prof.zMax_m + 1.0e3)


def test_range_error_reports_how_far_out(prof):
    with pytest.raises(PrescribedTemperatureRangeError) as e:
        prof(prof.zMax_m + 57.0e3, label='mantle')
    msg = str(e.value)
    assert 'mantle' in msg
    assert 'deeper than the profile' in msg
    assert '57.0' in msg, 'exceedance magnitude should be reported'


def test_array_with_one_bad_entry_raises(prof):
    z = np.array([0.0, 100e3, prof.zMax_m + 1e3])
    with pytest.raises(PrescribedTemperatureRangeError):
        prof(z)


# --- multiple candidate profiles (the important case) -------------------

def test_evaluates_correctly_for_every_candidate_radius_grid(prof):
    """PlanetProfile solves many candidate mantle/core configurations at once,
    each with its OWN radius grid of shape (nProfiles, nSilMax+1). T must be
    correct for every candidate independently, not just the first."""
    R = 1737.15e3
    nProfiles, nShells = 7, 40
    # Each candidate spans the surface down to a different inner radius.
    rInner = np.linspace(500e3, 420e3, nProfiles)
    rSil = np.array([np.linspace(R, ri, nShells + 1) for ri in rInner])

    T = np.empty((nProfiles, nShells + 1))
    for j in range(nShells + 1):
        T[:, j] = prof(R - rSil[:, j])

    # Row-by-row independent check against a direct per-candidate evaluation.
    for i in range(nProfiles):
        expect = prof(R - rSil[i, :])
        assert np.allclose(T[i, :], expect, rtol=0, atol=1e-12), \
            f'candidate {i} temperature mismatch'

    # Candidates must not be identical (otherwise the test is vacuous), and
    # each must be monotonic in depth.
    assert not np.allclose(T[0, :], T[-1, :])
    for i in range(nProfiles):
        assert np.all(np.diff(T[i, :]) >= -1e-9)


def test_deepest_candidate_beyond_data_raises(prof):
    """A candidate whose mantle extends past the thermal data must fail, not
    silently extrapolate -- this is how an over-deep mantle is caught."""
    R = 1737.15e3
    rSil = np.linspace(R, R - prof.zMax_m - 50e3, 30)
    with pytest.raises(PrescribedTemperatureRangeError):
        prof(R - rSil)


# --- Moho continuity ----------------------------------------------------

def test_temperature_continuous_across_moho(prof):
    """Crust and mantle draw from ONE profile, so T is continuous at the Moho
    by construction. Property discontinuities there are handled separately by
    the material mask."""
    zMoho = 40e3
    eps = 1.0
    below = float(prof(zMoho - eps))
    at = float(prof(zMoho))
    above = float(prof(zMoho + eps))
    assert below < at < above
    assert at - below == pytest.approx(above - at, rel=0.05), \
        'no artificial kink should appear at the Moho'


def test_moho_temperature_matches_project_value(prof):
    # ~475 K at 40 km, the value used in this project's shallow-mantle work.
    assert float(prof(40e3)) == pytest.approx(475.0, abs=1.0)
