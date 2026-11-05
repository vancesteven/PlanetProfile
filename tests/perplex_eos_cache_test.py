"""Tests for direct construction of a PerplexEOSStruct on a table already in the cache.

Motivation: PerplexEOSStruct.__init__ used to skip every attribute when its table
label was already in EOSlist.loaded, so a second direct constructor call returned
an object without domain limits or property functions. Callers outside
GetInnerEOS (scripts that evaluate a table directly) then had to memoise objects
themselves.

Uses a Perple_X table from the PlanetProfile cache; skipped if it is absent.
"""
import os

import numpy as np
import pytest

from PlanetProfile import _PERPLEXCACHE
from PlanetProfile.Thermodynamics.InnerEOS import PerplexEOSStruct
from PlanetProfile.Utilities.defineStructs import EOSlist

TABLE = 'moon_w05_preferred_stx24_unified.tab'
pytestmark = pytest.mark.skipif(not os.path.exists(os.path.join(_PERPLEXCACHE, TABLE)),
                                reason=f'{TABLE} not in the Perple_X cache')


def _fresh_pair(policy2='legacy_clamp'):
    EOSlist.loaded.clear(); EOSlist.ranges.clear()
    first = PerplexEOSStruct(TABLE, EOSinterpMethod='nearest', Fe_EOS=False)
    EOSlist.loaded[first.EOSlabel] = first
    second = PerplexEOSStruct(TABLE, EOSinterpMethod='nearest', Fe_EOS=False,
                              EOSdomainPolicy=policy2)
    return first, second


def test_second_construction_is_flagged_as_cached():
    first, second = _fresh_pair()
    assert not first.ALREADY_LOADED
    assert second.ALREADY_LOADED


def test_second_construction_keeps_domain_limits():
    first, second = _fresh_pair()
    for k in ('Pmin', 'Pmax', 'Tmin', 'Tmax'):
        assert getattr(second, k) == getattr(first, k)


def test_second_construction_gives_identical_properties():
    first, second = _fresh_pair()
    P = np.array([500.0, 2000.0, 4500.0]); T = np.array([800.0, 1300.0, 1700.0])
    assert np.array_equal(first.fn_rho_kgm3(P, T), second.fn_rho_kgm3(P, T))
    for a, b in zip(first.fn_Seismic(P, T), second.fn_Seismic(P, T)):
        assert np.array_equal(a, b)


def test_second_construction_keeps_its_own_domain_policy():
    first, second = _fresh_pair(policy2='error')
    assert first.EOSdomainPolicy == 'legacy_clamp'
    assert second.EOSdomainPolicy == 'error'
