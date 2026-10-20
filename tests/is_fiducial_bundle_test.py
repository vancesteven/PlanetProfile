"""Tests for the shipped fiducial IS bundle loader (conditions B1-B4,
validation_reports/titan_freegrav_nh3_1m/RATIFICATION.md 2026-10-20)."""
import copy
import hashlib
import json
import sys
import types
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Inference import is_fiducial_bundle as fb  # noqa: E402
from PlanetProfile.Inference.is_correction import (  # noqa: E402
    compute_is_correction)

N = 4000
X_OBS = {'C20': -3.3511e-5, 'C22': 1.0107e-5, 'Re_k2': 0.608, 'Im_k2': 0.135}
SIG = {'C20': 4.71e-7, 'C22': 6.2e-8, 'Re_k2': 0.048, 'Im_k2': 0.035}


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _write_bundle(tmp, log_w_scale=0.5, c53_pass=True, mutate_manifest=None):
    rng = np.random.default_rng(0)
    flp = rng.normal(size=N)
    lw = rng.normal(scale=log_w_scale, size=N)
    ll = lw + flp
    rej = rng.random(N) < 0.01
    ll[rej] = -1e30
    log_w = np.where(rej, -np.inf, ll - flp)
    d_ocean = np.where(rng.random(N) < 0.9, 100.0, 0.0)
    npz = tmp / 'b.npz'
    np.savez_compressed(
        npz, samples=rng.normal(size=(N, 3)), log_w=log_w,
        log_likelihoods=ll, flow_log_prob=flp, D_ocean_results=d_ocean,
        k2_results=np.column_stack([np.full(N, .6), np.full(N, .1)]),
        param_names=np.array(['a', 'b', 'c']))
    art = tmp / 'art.pt'; art.write_bytes(b'artifact')
    cache = tmp / 'cache.pkl'; cache.write_bytes(b'cache')
    corr = fb.recompute_diagnostics(log_w)
    m = {
        'npz_file': 'b.npz', 'npz_sha256': _sha(npz),
        'artifact_sha256': _sha(art), 'structure_cache_sha256': _sha(cache),
        'config_hash': 'abc123', 'config_file': 'nonexistent.json',
        'x_obs': X_OBS, 'obs_sigma': SIG, 'observables': sorted(X_OBS),
        'n_samples': N,
        'diagnostics': {'ess': corr.ess, 'pareto_k': corr.pareto_k,
                        'w_max_frac': corr.w_max_frac,
                        'frac_rejected': corr.frac_rejected},
        'c5_3': {'pass': c53_pass},
    }
    if mutate_manifest:
        mutate_manifest(m)
    mp = tmp / 'b.json'
    mp.write_text(json.dumps(m))
    return mp, art, cache, npz


def _kw(art, cache, **over):
    kw = dict(artifact_path=art, config='abc123', cache_path=cache,
              x_obs=dict(X_OBS), obs_sigma=dict(SIG), n_post=N,
              n_derived=None)
    kw.update(over)
    return kw


def test_exact_match_activates(tmp_path):
    mp, art, cache, _ = _write_bundle(tmp_path)
    r = fb.load_fiducial_correction(mp, **_kw(art, cache))
    assert isinstance(r, fb.CorrectedBundle), getattr(r, 'reason', '')
    assert r.samples.shape == (N, 3)
    assert abs(r.weights.sum() - 1) < 1e-12
    assert 'D_ocean_results' in r.derived and 'k2_results' in r.derived
    # config may also be an object exposing generate_hash()
    cfg = types.SimpleNamespace(generate_hash=lambda: 'abc123')
    assert isinstance(fb.load_fiducial_correction(
        mp, **_kw(art, cache, config=cfg)), fb.CorrectedBundle)
    # n_post=None accepted
    assert isinstance(fb.load_fiducial_correction(
        mp, **_kw(art, cache, n_post=None)), fb.CorrectedBundle)


def test_feature_flag_default_off():
    assert fb.FEATURE_ENABLED is False


@pytest.mark.parametrize('key', sorted(X_OBS))
def test_x_obs_value_change(tmp_path, key):
    mp, art, cache, _ = _write_bundle(tmp_path)
    x = dict(X_OBS); x[key] *= (1 + 1e-6)
    r = fb.load_fiducial_correction(mp, **_kw(art, cache, x_obs=x))
    assert isinstance(r, fb.Fallback) and r.reason.startswith(fb.R_MODIFIED)
    # sub-tolerance perturbation (<1e-12 rel) still matches
    x = dict(X_OBS); x[key] *= (1 + 1e-14)
    assert isinstance(fb.load_fiducial_correction(
        mp, **_kw(art, cache, x_obs=x)), fb.CorrectedBundle)


def test_x_obs_keyset_change(tmp_path):
    mp, art, cache, _ = _write_bundle(tmp_path)
    x = dict(X_OBS); x['extra'] = 1.0
    r = fb.load_fiducial_correction(mp, **_kw(art, cache, x_obs=x))
    assert r.reason.startswith(fb.R_MODIFIED)
    x = dict(X_OBS); x.pop('C22')
    r = fb.load_fiducial_correction(mp, **_kw(art, cache, x_obs=x))
    assert r.reason.startswith(fb.R_MODIFIED)


@pytest.mark.parametrize('key', sorted(SIG))
def test_sigma_change(tmp_path, key):
    mp, art, cache, _ = _write_bundle(tmp_path)
    s = dict(SIG); s[key] *= 1.01
    r = fb.load_fiducial_correction(mp, **_kw(art, cache, obs_sigma=s))
    assert r.reason.startswith(fb.R_MODIFIED)


def test_hash_changes(tmp_path):
    mp, art, cache, _ = _write_bundle(tmp_path)
    r = fb.load_fiducial_correction(mp, **_kw(art, cache, config='other'))
    assert r.reason.startswith(fb.R_MODIFIED) and 'config' in r.reason
    art2 = tmp_path / 'art2.pt'; art2.write_bytes(b'artifact2')
    r = fb.load_fiducial_correction(mp, **_kw(art2, cache))
    assert r.reason.startswith(fb.R_MODIFIED) and 'artifact' in r.reason
    cache2 = tmp_path / 'cache2.pkl'; cache2.write_bytes(b'cache2')
    r = fb.load_fiducial_correction(mp, **_kw(art, cache2))
    assert r.reason.startswith(fb.R_MODIFIED) and 'cache' in r.reason


@pytest.mark.parametrize('nd', [1000, 1, N - 1])
def test_n_derived_truncated_falls_back(tmp_path, nd):
    mp, art, cache, _ = _write_bundle(tmp_path)
    r = fb.load_fiducial_correction(mp, **_kw(art, cache, n_derived=nd))
    assert isinstance(r, fb.Fallback) and r.reason.startswith(fb.R_FULLN)


def test_n_post_mismatch_falls_back(tmp_path):
    mp, art, cache, _ = _write_bundle(tmp_path)
    r = fb.load_fiducial_correction(mp, **_kw(art, cache, n_post=N + 1))
    assert isinstance(r, fb.Fallback)


def test_tampered_npz(tmp_path):
    mp, art, cache, npz = _write_bundle(tmp_path)
    npz.write_bytes(npz.read_bytes() + b'\0')
    r = fb.load_fiducial_correction(mp, **_kw(art, cache))
    assert r.reason.startswith(fb.R_INTEGRITY) and 'sha256' in r.reason
    # missing npz
    npz.unlink()
    r = fb.load_fiducial_correction(mp, **_kw(art, cache))
    assert r.reason.startswith(fb.R_INTEGRITY)


def test_inconsistent_logw_with_matching_sha(tmp_path):
    """Tamper log_w and re-sign the sha: the log_w == ll - flp cross-check
    still catches it."""
    mp, art, cache, npz = _write_bundle(tmp_path)
    with np.load(npz) as z:
        a = {k: z[k] for k in z.files}
    a['log_w'][5] += 1.0
    np.savez_compressed(npz, **a)
    m = json.loads(mp.read_text()); m['npz_sha256'] = _sha(npz)
    mp.write_text(json.dumps(m))
    r = fb.load_fiducial_correction(mp, **_kw(art, cache))
    assert r.reason.startswith(fb.R_INTEGRITY)


def test_c53_absent_or_failing(tmp_path):
    mp, art, cache, _ = _write_bundle(tmp_path, c53_pass=False)
    assert fb.load_fiducial_correction(
        mp, **_kw(art, cache)).reason.startswith(fb.R_INTEGRITY)
    d = tmp_path / 'd'; d.mkdir()
    mp, art, cache, _ = _write_bundle(
        d, mutate_manifest=lambda m: m.pop('c5_3'))
    assert fb.load_fiducial_correction(
        mp, **_kw(art, cache)).reason.startswith(fb.R_INTEGRITY)


def test_manifest_diag_mismatch(tmp_path):
    def mut(m):
        m['diagnostics']['ess'] *= 1.1
    mp, art, cache, _ = _write_bundle(tmp_path, mutate_manifest=mut)
    r = fb.load_fiducial_correction(mp, **_kw(art, cache))
    assert r.reason.startswith(fb.R_INTEGRITY) and 'ess' in r.reason


def test_low_ess_falls_back(tmp_path):
    mp, art, cache, _ = _write_bundle(tmp_path, log_w_scale=2.0)
    r = fb.load_fiducial_correction(mp, **_kw(art, cache))
    assert isinstance(r, fb.Fallback)
    assert r.reason.split(':')[0] in (fb.R_ESS, fb.R_PARETO, fb.R_WMAX)


def test_missing_manifest_never_raises(tmp_path):
    r = fb.load_fiducial_correction(tmp_path / 'nope.json',
                                    **_kw(tmp_path / 'a', tmp_path / 'c'))
    assert isinstance(r, fb.Fallback) and r.reason.startswith(fb.R_INTEGRITY)


def test_diagnostics_recompute_equals_manifest(tmp_path):
    mp, art, cache, _ = _write_bundle(tmp_path)
    r = fb.load_fiducial_correction(mp, **_kw(art, cache))
    md = r.manifest['diagnostics']
    assert np.isclose(r.diagnostics['ess'], md['ess'], rtol=1e-12)
    assert np.isclose(r.diagnostics['pareto_k'], md['pareto_k'], rtol=1e-12)
    assert np.isclose(r.diagnostics['w_max_frac'], md['w_max_frac'],
                      rtol=1e-12)
    assert np.isclose(r.diagnostics['frac_rejected'], md['frac_rejected'])


def test_compute_is_correction_raises_on_nan_padding():
    """B4 backstop: a truncated-derived result (NaN-padded likelihoods)
    must raise, not yield silent wrong weights."""
    ll = np.random.default_rng(0).normal(size=100)
    ll[50:] = np.nan
    res = types.SimpleNamespace(
        log_likelihoods=ll, metadata={'flow_log_prob': np.zeros(100)},
        D_ocean_results=np.array([]), k2_results=np.array([]))
    with pytest.raises(ValueError, match='n_derived=None'):
        compute_is_correction(res)


# ---- shipped bundle (skipped if not built / inputs absent) --------------
SHIPPED = REPO / 'PlanetProfile/Inference/sbi_artifacts/' \
    'titan_freegrav_nh3_is_fiducial.json'


@pytest.mark.skipif(not SHIPPED.exists(), reason='shipped bundle absent')
def test_shipped_bundle_manifest_and_recompute():
    m = json.loads(SHIPPED.read_text())
    npz = SHIPPED.parent / m['npz_file']
    assert _sha(npz) == m['npz_sha256']
    with np.load(npz) as z:
        log_w = z['log_w']
    assert len(log_w) == m['n_samples'] == 20000
    committed = np.load(REPO / 'validation_reports/titan_freegrav_nh3_1m/'
                        'is_correction/log_weights.npy')
    assert np.array_equal(np.isneginf(log_w), np.isneginf(committed))
    f = np.isfinite(committed)
    assert np.allclose(log_w[f], committed[f], rtol=1e-10, atol=0)
    corr = fb.recompute_diagnostics(log_w)
    d = m['diagnostics']
    assert np.isclose(corr.ess, d['ess'], rtol=1e-9)
    assert np.isclose(corr.pareto_k, d['pareto_k'], rtol=1e-9)
    assert np.isclose(corr.w_max_frac, d['w_max_frac'], rtol=1e-9)
    assert corr.ess >= 1000 and corr.pareto_k <= 0.7
    assert m['c5_3']['pass'] is True
    assert m['x_obs']['Re_k2'] == 0.608


@pytest.mark.skipif(not SHIPPED.exists(), reason='shipped bundle absent')
def test_shipped_bundle_activation_and_modified_datum():
    m = json.loads(SHIPPED.read_text())
    art = REPO / 'PlanetProfile/Inference/sbi_artifacts' / m['artifact_name']
    cache = REPO / 'PlanetProfile/Test/mcmc_results/Titan/Test54_nh3_ocean/' \
        'titan_nh3_joint_structure_grid_2d.pkl'
    if not (art.exists() and cache.exists()):
        pytest.skip('artifact/cache absent')
    kw = dict(artifact_path=art, config=m['config_hash'], cache_path=cache,
              x_obs=m['x_obs'], obs_sigma=m['obs_sigma'], n_post=20000,
              n_derived=None)
    r = fb.load_fiducial_correction(SHIPPED, **kw)
    assert isinstance(r, fb.CorrectedBundle), getattr(r, 'reason', '')
    assert r.samples.shape == (20000, 13)
    x = dict(m['x_obs']); x['Re_k2'] += 0.01
    kw2 = dict(kw, x_obs=x)
    assert fb.load_fiducial_correction(SHIPPED, **kw2).reason.startswith(
        fb.R_MODIFIED)
    assert fb.load_fiducial_correction(
        SHIPPED, **dict(kw, n_derived=1000)).reason.startswith(fb.R_FULLN)
