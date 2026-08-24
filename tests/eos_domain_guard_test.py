"""Tests for the EOS-domain policy guard.

Motivation: with EXTRAP=False the historical code silently replaced any
out-of-domain (P, T) with the nearest table edge, so an unvalidated state
returned plausible numbers instead of an error. These tests pin the three
supported policies.

Uses a small synthetic domain -- no Perple_X table or full run required.
"""
import logging

import numpy as np
import pytest

from PlanetProfile.Utilities.DataManip import (
    EnforceEOSDomain,
    EOSDomainError,
    EOSDomainExceedance,
    ResetNearestExtrap,
)

# Synthetic stand-in for the unified Warren stx24 domain
# (P 1000-50000 bar = 100-5000 MPa, T 300-2200 K).
PMIN, PMAX, TMIN, TMAX = 100.0, 5000.0, 300.0, 2200.0


def _enforce(P, T, policy, CLAMP=True):
    return EnforceEOSDomain(np.atleast_1d(np.array(P, dtype=float)),
                            np.atleast_1d(np.array(T, dtype=float)),
                            PMIN, PMAX, TMIN, TMAX,
                            policy=policy, CLAMP=CLAMP, EOSname='test EOS')


# --- In-domain ----------------------------------------------------------

@pytest.mark.parametrize('policy', ['legacy_clamp', 'warn', 'error'])
def test_in_domain_values_unchanged_and_silent(policy, caplog):
    P = np.array([100.0, 2500.0, 5000.0])
    T = np.array([300.0, 1500.0, 2200.0])
    with caplog.at_level(logging.WARNING, logger='PlanetProfile'):
        Pout, Tout = _enforce(P, T, policy)
    assert np.array_equal(Pout, P)
    assert np.array_equal(Tout, T)
    assert caplog.text == '', 'in-domain query must not warn'


def test_in_domain_reports_no_exceedance():
    r = EOSDomainExceedance(np.array([200.0]), np.array([800.0]), PMIN, PMAX, TMIN, TMAX)
    assert r['nOutside'] == 0
    assert r['max1Above'] == 0.0 and r['max2Above'] == 0.0


# --- Individual violations ---------------------------------------------

def test_lower_pressure_violation_detected():
    r = EOSDomainExceedance(np.array([50.0]), np.array([800.0]), PMIN, PMAX, TMIN, TMAX)
    assert r['n1Low'] == 1 and r['nOutside'] == 1
    assert r['max1Below'] == pytest.approx(50.0)


def test_upper_pressure_violation_detected():
    r = EOSDomainExceedance(np.array([9000.0]), np.array([800.0]), PMIN, PMAX, TMIN, TMAX)
    assert r['n1High'] == 1
    assert r['max1Above'] == pytest.approx(4000.0)


def test_lower_temperature_violation_detected():
    r = EOSDomainExceedance(np.array([500.0]), np.array([250.0]), PMIN, PMAX, TMIN, TMAX)
    assert r['n2Low'] == 1
    assert r['max2Below'] == pytest.approx(50.0)


def test_upper_temperature_violation_detected():
    r = EOSDomainExceedance(np.array([500.0]), np.array([2600.0]), PMIN, PMAX, TMIN, TMAX)
    assert r['n2High'] == 1
    assert r['max2Above'] == pytest.approx(400.0)


# --- Mixed array: only violating entries reported -----------------------

def test_mixed_array_counts_only_violating_nodes():
    P = np.array([50.0, 500.0, 2500.0, 9000.0])     # 2 violations (low, high)
    T = np.array([800.0, 800.0, 800.0, 800.0])      # all in domain
    r = EOSDomainExceedance(P, T, PMIN, PMAX, TMIN, TMAX)
    assert r['nTotal'] == 4
    assert r['nOutside'] == 2
    assert r['n1Low'] == 1 and r['n1High'] == 1
    assert r['n2Low'] == 0 and r['n2High'] == 0


def test_warn_message_names_eos_counts_and_exceedance(caplog):
    P = np.array([500.0, 9000.0])
    T = np.array([800.0, 2600.0])
    with caplog.at_level(logging.WARNING, logger='PlanetProfile'):
        _enforce(P, T, 'warn')
    msg = caplog.text
    assert 'test EOS' in msg
    assert '1 of 2 node(s)' in msg
    assert '4000' in msg, 'max pressure exceedance not reported'
    assert '400' in msg, 'max temperature exceedance not reported'


# --- Legacy compatibility ----------------------------------------------

def test_legacy_mode_reproduces_reset_nearest_extrap_exactly():
    P = np.array([50.0, 500.0, 9000.0])
    T = np.array([250.0, 800.0, 2600.0])
    expectP, expectT = ResetNearestExtrap(P.copy(), T.copy(), PMIN, PMAX, TMIN, TMAX)
    gotP, gotT = _enforce(P, T, 'legacy_clamp')
    assert np.array_equal(gotP, expectP)
    assert np.array_equal(gotT, expectT)


def test_legacy_mode_is_silent(caplog):
    with caplog.at_level(logging.WARNING, logger='PlanetProfile'):
        _enforce(np.array([9000.0]), np.array([2600.0]), 'legacy_clamp')
    assert caplog.text == '', 'legacy mode must stay silent for backwards compatibility'


def test_warn_mode_clamps_identically_to_legacy():
    P = np.array([9000.0]); T = np.array([2600.0])
    legacyP, legacyT = _enforce(P, T, 'legacy_clamp')
    warnP, warnT = _enforce(P, T, 'warn')
    assert np.array_equal(warnP, legacyP) and np.array_equal(warnT, legacyT)
    assert warnP[0] == PMAX and warnT[0] == TMAX


# --- Strict mode --------------------------------------------------------

def test_strict_mode_raises_before_evaluation():
    with pytest.raises(EOSDomainError, match='outside the tabulated EOS domain'):
        _enforce(np.array([9000.0]), np.array([800.0]), 'error')


def test_strict_mode_passes_in_domain():
    P = np.array([2500.0]); T = np.array([1500.0])
    Pout, Tout = _enforce(P, T, 'error')
    assert Pout[0] == 2500.0 and Tout[0] == 1500.0


def test_strict_mode_raises_even_when_extrapolation_allowed():
    """EXTRAP=True must not make a deliberate out-of-domain query silent."""
    with pytest.raises(EOSDomainError):
        _enforce(np.array([9000.0]), np.array([800.0]), 'error', CLAMP=False)


def test_unknown_policy_rejected():
    with pytest.raises(ValueError, match='Unrecognized EOS domain policy'):
        _enforce(np.array([500.0]), np.array([800.0]), 'clamp_quietly_please')
