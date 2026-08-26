"""EOS-domain strictness must apply to physically active nodes only.

PlanetProfile propagates every candidate silicate profile toward r = 0 as
scratch storage for the core-radius search. Those unused nodes routinely exceed
the EOS pressure ceiling. Aborting on them would reject sound models; ignoring
a violation on a node that IS part of a candidate structure would hide a real
one. These tests pin both halves.
"""
import numpy as np
import pytest

from PlanetProfile.Thermodynamics.InnerEOS import (
    EOSCandidatePropagation,
    EOSCandidatePropagationActive,
    ValidateCandidateEOSDomains,
)
from PlanetProfile.Utilities.DataManip import EOSDomainError
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist

PMIN, PMAX, TMIN, TMAX = 100.0, 5000.0, 300.0, 2200.0


class _FakeEOS:
    Pmin, Pmax, Tmin, Tmax = PMIN, PMAX, TMIN, TMAX
    comp = 'fake'
    Fe_EOS = False


class _Slot:
    def __init__(self, key): self.key = key


class _Crust:
    def __init__(self, enabled=False, mask=None):
        self.enabled = enabled
        self.material_mask = mask


class _Planet:
    """Minimal stand-in exposing only what the validator reads."""
    def __init__(self, crust=None):
        self.Sil = type('S', (), {'EOS': _Slot('fakeSil')})()
        self.Core = type('C', (), {'EOS': _Slot('fakeCore')})()
        self.Crust = crust if crust is not None else _Crust()


@pytest.fixture(autouse=True)
def _register(monkeypatch):
    monkeypatch.setitem(EOSlist.loaded, 'fakeSil', _FakeEOS())
    monkeypatch.setitem(EOSlist.loaded, 'fakeCore', _FakeEOS())
    yield


def _grid(nCand=1, nShell=10, Pact=1000.0, Tact=900.0, Pscratch=9000.0, nAct=6):
    """Active nodes in-domain by default; scratch nodes far out of domain."""
    P = np.full((nCand, nShell), Pscratch)
    T = np.full((nCand, nShell), 1500.0)
    P[:, :nAct] = Pact
    T[:, :nAct] = Tact
    return P, T, np.full(nCand, nAct), np.arange(nCand)


# --- scratch violations must NOT abort -----------------------------------

def test_scratch_node_out_of_domain_does_not_abort():
    P, T, nSilFinal, inds = _grid(Pscratch=9000.0)
    rep = ValidateCandidateEOSDomains(_Planet(), P, T, nSilFinal, inds, policy='error')
    assert rep['nBadNodes'] == 0, 'scratch storage must not count as a violation'


def test_scratch_thermal_violation_does_not_abort():
    P, T, nSilFinal, inds = _grid()
    T[:, 6:] = 4000.0          # scratch far above the T ceiling
    rep = ValidateCandidateEOSDomains(_Planet(), P, T, nSilFinal, inds, policy='error')
    assert rep['nBadNodes'] == 0


# --- active violations MUST abort ----------------------------------------

def test_active_pressure_violation_aborts():
    P, T, nSilFinal, inds = _grid()
    P[0, 3] = 9000.0           # inside the active extent
    with pytest.raises(EOSDomainError, match='PHYSICALLY ACTIVE'):
        ValidateCandidateEOSDomains(_Planet(), P, T, nSilFinal, inds, policy='error')


def test_active_thermal_violation_aborts():
    P, T, nSilFinal, inds = _grid()
    T[0, 2] = 2600.0
    with pytest.raises(EOSDomainError, match='PHYSICALLY ACTIVE'):
        ValidateCandidateEOSDomains(_Planet(), P, T, nSilFinal, inds, policy='error')


def test_active_violation_in_second_candidate_also_aborts():
    P, T, nSilFinal, inds = _grid(nCand=3)
    P[2, 1] = 20.0             # below Pmin, active extent of candidate 2
    with pytest.raises(EOSDomainError):
        ValidateCandidateEOSDomains(_Planet(), P, T, nSilFinal, inds, policy='error')


# --- fully valid candidate passes ----------------------------------------

def test_fully_in_domain_candidate_passes():
    P, T, nSilFinal, inds = _grid(nCand=4)
    rep = ValidateCandidateEOSDomains(_Planet(), P, T, nSilFinal, inds, policy='error')
    assert rep['nBadNodes'] == 0 and rep['nBadCandidates'] == 0
    assert rep['nCandidates'] == 4


def test_core_nodes_are_always_checked():
    P, T, nSilFinal, inds = _grid()
    Pc = np.array([[4000.0, 9000.0]]); Tc = np.array([[1700.0, 1700.0]])
    with pytest.raises(EOSDomainError, match='core'):
        ValidateCandidateEOSDomains(_Planet(), P, T, nSilFinal, inds,
                                    Pcore_MPa=Pc, Tcore_K=Tc, policy='error')


# --- crust nodes are excluded from the MANTLE EOS check ------------------

def test_crust_nodes_excluded_from_mantle_domain_check():
    """An explicit crust supplies prescribed constants, so its shallow nodes
    legitimately sit below the mantle table's P and T floors."""
    nShell, nAct = 10, 6
    P, T, nSilFinal, inds = _grid(nShell=nShell, nAct=nAct)
    P[0, :2] = 0.0             # crust: below Pmin
    T[0, :2] = 250.0           # crust: below Tmin
    mask = np.full((1, nShell), Constants.materialMantle, dtype=np.int8)
    mask[0, :2] = Constants.materialCrust
    planet = _Planet(crust=_Crust(enabled=True, mask=mask))
    rep = ValidateCandidateEOSDomains(planet, P, T, nSilFinal, inds, policy='error')
    assert rep['nBadNodes'] == 0
    assert rep['nCrustExcluded'] == 2


def test_mantle_violation_still_caught_when_crust_present():
    nShell, nAct = 10, 6
    P, T, nSilFinal, inds = _grid(nShell=nShell, nAct=nAct)
    P[0, :2] = 0.0; T[0, :2] = 250.0
    P[0, 4] = 9000.0           # a genuine MANTLE violation
    mask = np.full((1, nShell), Constants.materialMantle, dtype=np.int8)
    mask[0, :2] = Constants.materialCrust
    planet = _Planet(crust=_Crust(enabled=True, mask=mask))
    with pytest.raises(EOSDomainError):
        ValidateCandidateEOSDomains(planet, P, T, nSilFinal, inds, policy='error')


# --- policy passthrough ---------------------------------------------------

def test_warn_policy_reports_without_raising():
    P, T, nSilFinal, inds = _grid()
    P[0, 3] = 9000.0
    rep = ValidateCandidateEOSDomains(_Planet(), P, T, nSilFinal, inds, policy='warn')
    assert rep['nBadNodes'] == 1


def test_legacy_policy_reports_without_raising():
    P, T, nSilFinal, inds = _grid()
    P[0, 3] = 9000.0
    rep = ValidateCandidateEOSDomains(_Planet(), P, T, nSilFinal, inds, policy='legacy_clamp')
    assert rep['nBadNodes'] == 1


# --- propagation phase ----------------------------------------------------

def test_candidate_propagation_context_is_scoped():
    assert not EOSCandidatePropagationActive()
    with EOSCandidatePropagation():
        assert EOSCandidatePropagationActive()
    assert not EOSCandidatePropagationActive()


def test_propagation_phase_downgrades_only_error_policy():
    """'error' is relaxed during candidate propagation; 'warn' is untouched so
    diagnostics stay visible."""
    from PlanetProfile.Thermodynamics.InnerEOS import PerplexEOSStruct
    import logging

    class Stub:
        Pmin, Pmax, Tmin, Tmax = PMIN, PMAX, TMIN, TMAX
        EXTRAP = False; comp = 'stub'; Fe_EOS = False
        _EnforceDomain = PerplexEOSStruct._EnforceDomain

    s = Stub(); s.EOSdomainPolicy = 'error'
    with EOSCandidatePropagation():
        P, _ = s._EnforceDomain(np.array([9000.0]), np.array([1500.0]))
    assert P[0] == PMAX, 'strict policy should clamp, not raise, during propagation'

    with pytest.raises(EOSDomainError):
        s._EnforceDomain(np.array([9000.0]), np.array([1500.0]))
