"""Tests for InnerEOS.EnsurePerplexTableInCache: on a fresh install the lunar tables reach the cache.

The mantle table shipped in PlanetProfile/Thermodynamics/EOStables/Perple_X is copied to an empty
cache, and the traceable Fe-S core table is built there from the werami tables of the repository
and reproduces its reference semantic SHA-256.
"""
import hashlib
import importlib.util
import os

import pytest

import PlanetProfile.Thermodynamics.InnerEOS as InnerEOS
from PlanetProfile import _ROOT

MANTLE = 'moon_w05_preferred_stx24_unified.tab'
MANTLE_SHA256 = 'eb254ac5d454cd59f5ca1c29220214eb8d036cc21b44251b6db558c9df2d8348'


@pytest.fixture
def empty_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(InnerEOS, '_PERPLEXCACHE', str(tmp_path))
    return tmp_path


def _sha256(fpath):
    h = hashlib.sha256()
    with open(fpath, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def test_shipped_mantle_table_is_copied_to_empty_cache(empty_cache):
    fpath = InnerEOS.EnsurePerplexTableInCache(MANTLE)
    assert fpath == os.path.join(str(empty_cache), MANTLE)
    assert _sha256(fpath) == MANTLE_SHA256
    assert sorted(os.listdir(empty_cache)) == [MANTLE]


def test_existing_cache_entry_is_left_untouched(empty_cache):
    (empty_cache / MANTLE).write_text('user table')
    InnerEOS.EnsurePerplexTableInCache(MANTLE)
    assert (empty_cache / MANTLE).read_text() == 'user table'


def test_traceable_core_table_is_built_and_matches_reference(empty_cache):
    from scipy.io import loadmat
    fpath = InnerEOS.EnsurePerplexTableInCache(InnerEOS.TRACEABLE_FES_NAME)
    builder = os.path.join(os.path.dirname(_ROOT), 'Thermodynamics', 'Perple_X', 'stackFeSTraceable.py')
    spec = importlib.util.spec_from_file_location('stackFeSTraceable', builder)
    stack = importlib.util.module_from_spec(spec); spec.loader.exec_module(stack)
    m = loadmat(fpath)
    out = {k: m[k] for k in stack.PROPS}
    out.update({k: m[k][0] for k in stack.FIELDS[:4]})
    assert stack.SemanticSHA256(out) == stack.SEMANTIC_SHA256
    assert sorted(os.listdir(empty_cache)) == [InnerEOS.TRACEABLE_FES_NAME]


def test_unknown_missing_table_is_left_to_the_loader(empty_cache):
    fpath = InnerEOS.EnsurePerplexTableInCache('not_a_table.tab')
    assert not os.path.exists(fpath)
