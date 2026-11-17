"""Ae induction sidecar keying: a runner must never overwrite a sidecar whose
key differs from its own (two configs with different excitation-label sets
share the Test51 cache; alternating them used to recompute and rewrite a
tracked file from the MCMCRunner constructor).

All I/O happens on a COPY of the Test51 cache + committed sidecar in tmp_path.
``forward_model_induction`` is replaced by a deterministic fake so the test is
fast; the keying/caching plumbing under test is independent of the physics.
"""
import json
import os
import pickle
import shutil
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
T51 = REPO / 'PlanetProfile/Test/mcmc_results/Europa/Test51_seawater'
CACHE_NAME = 'europa_seawater_structure_grid.pkl'
CFG1 = REPO / 'PlanetProfile/Inference/configs/europa_galileo_v1p1_8D.json'
CFG3 = REPO / 'PlanetProfile/Inference/configs/europa_seawater_andrade_clipper_v2.json'

pytestmark = pytest.mark.skipif(
    not ((T51 / CACHE_NAME).exists()
         and (T51 / (CACHE_NAME + '.ae_sidecar.pkl')).exists()),
    reason='Test51 cache/sidecar not present')


@pytest.fixture
def env(tmp_path, monkeypatch):
    for name in (CACHE_NAME, CACHE_NAME + '.ae_sidecar.pkl'):
        shutil.copyfile(T51 / name, tmp_path / name)
    cache = tmp_path / CACHE_NAME
    legacy = tmp_path / (CACHE_NAME + '.ae_sidecar.pkl')

    from PlanetProfile.Inference import forward_models, mcmc_runner
    calls = {'compute': 0}

    def fake_induction(struct, freq_dict, nn=1, do_parallel=False):
        return {lbl: complex(sum(map(ord, lbl)), len(freq_dict))
                for lbl in freq_dict}
    monkeypatch.setattr(forward_models, 'forward_model_induction',
                        fake_induction)
    real_compute = mcmc_runner.MCMCRunner._compute_ae_grid

    def spy(self, structures, labels):
        calls['compute'] += 1
        return real_compute(self, structures, labels)
    monkeypatch.setattr(mcmc_runner.MCMCRunner, '_compute_ae_grid', spy)

    def make(cfg_path):
        from PlanetProfile.Inference.inference_core import InferenceConfig
        raw = json.loads(Path(cfg_path).read_text())
        raw['structure_cache_path'] = str(cache)
        return mcmc_runner.MCMCRunner(InferenceConfig(**raw))
    return {'dir': tmp_path, 'legacy': legacy, 'calls': calls, 'make': make}


def _siblings(d):
    return sorted(p for p in d.iterdir()
                  if '.ae_sidecar.' in p.name and p.name.endswith('.pkl')
                  and not p.name.endswith('.ae_sidecar.pkl'))


def test_a_one_label_config_uses_legacy_no_recompute(env):
    before = env['legacy'].read_bytes()
    r = env['make'](CFG1)
    assert env['calls']['compute'] == 0
    assert env['legacy'].read_bytes() == before
    assert _siblings(env['dir']) == []
    assert len(r._ae_grid_cache) > 0


def test_b_three_label_config_writes_keyed_sibling_once(env):
    before = env['legacy'].read_bytes()
    r1 = env['make'](CFG3)
    assert env['calls']['compute'] == 1
    assert env['legacy'].read_bytes() == before, 'legacy sidecar was modified'
    sib = _siblings(env['dir'])
    assert len(sib) == 1
    assert not list(env['dir'].glob('*.tmp'))
    r2 = env['make'](CFG3)
    assert env['calls']['compute'] == 1, 'second construction recomputed'
    assert env['legacy'].read_bytes() == before
    # alternating back to the 1-label config is still instant and clean
    env['make'](CFG1)
    assert env['calls']['compute'] == 1
    assert env['legacy'].read_bytes() == before
    assert r2._ae_grid_cache == r1._ae_grid_cache


def test_c_keyed_sibling_equals_fresh_in_memory_compute(env):
    r1 = env['make'](CFG3)
    r2 = env['make'](CFG3)             # loaded from the keyed sibling
    assert env['calls']['compute'] == 1
    sd = r2.structure_data
    labels = {'orbital', 'synodic', 'synodic 2nd'}
    r2._ae_grid_cache = {}
    r2._compute_ae_grid(sd['structures'], labels)   # fresh in-memory
    fresh = r2._ae_grid_cache
    with open(_siblings(env['dir'])[0], 'rb') as f:
        stored = pickle.load(f)['ae_grid']
    assert stored == fresh
    assert r1._ae_grid_cache == fresh


@pytest.mark.skipif(os.geteuid() == 0, reason='root ignores directory modes')
def test_d_read_only_directory_falls_back_to_memory(env):
    d = env['dir']
    before = env['legacy'].read_bytes()
    d.chmod(0o555)
    try:
        r = env['make'](CFG3)
        assert env['calls']['compute'] == 1
        assert len(r._ae_grid_cache) > 0
        assert _siblings(d) == []
        assert env['legacy'].read_bytes() == before
    finally:
        d.chmod(0o755)


def test_keyed_path_helpers():
    from PlanetProfile.Inference.mcmc_runner import (
        _ae_sidecar_keyed_path, _ae_sidecar_atomic_write,
        _ae_sidecar_try_load)
    k1 = {'labels': ('a',), 'n': 1, 'cache_size': 2}
    k2 = {'labels': ('a', 'b'), 'n': 1, 'cache_size': 2}
    p = '/x/c.pkl.ae_sidecar.pkl'
    a, b = _ae_sidecar_keyed_path(p, k1), _ae_sidecar_keyed_path(p, k2)
    assert a != b and a != p and a.startswith('/x/c.pkl.ae_sidecar.')
    assert a == _ae_sidecar_keyed_path(p, dict(k1))   # deterministic
