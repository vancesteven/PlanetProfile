"""End-to-end exercise of the Enceladus ocean6D_rhorock config on a
production-schema zbw cache (nodes WITHOUT ``eta_Pa_base``).

Regression for: ``run_inference_cli`` crashed on the first likelihood call with
``KeyError: 'eta_Pa_base'`` because ``log_likelihood`` unconditionally solved
tidal k2 (and rejected on NaN) although the config observes only
C20/C22/C30/libration. Gate: ``MCMCRunner._tidal_solve_needed``.

The cache here is SYNTHETIC (mass-consistent 3-layer nodes, same schema keys
as ``cache_builder.build_zbw_grid_cache`` nodes minus the tidal-solver
inputs production nodes also lack).
"""
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402
from PlanetProfile.Inference.mcmc_runner import MCMCRunner  # noqa: E402

CFG = (REPO / 'PlanetProfile' / 'Inference' / 'configs'
       / 'enceladus_cassini_isostasy_ocean6D_rhorock.json')

R_T = 252.22e3
M = 1.08022e20
OMEGA = 5.307e-5
ECC = 0.0047
V = 4.0 * np.pi / 3.0
ZB_GRID = np.array([5.0, 15.0, 25.0, 35.0, 45.0])
W_GRID = np.array([0.1, 1.0, 10.0, 100.0])


def _node(zb_km, w_ppt, d_ocean_km=40.0):
    """Mass-consistent rock | ocean | ice-Ih node, 4 cells per layer.
    Deliberately has NO ``eta_Pa_base`` (nor any other tidal-solver input)."""
    rho_i = 925.0
    rho_o = 1000.0 + 0.7 * w_ppt          # salinity-dependent, ppt -> kg/m3
    R_b = R_T - zb_km * 1e3
    R_s = R_b - d_ocean_km * 1e3
    rho_s = (M - rho_o * V * (R_b ** 3 - R_s ** 3)
             - rho_i * V * (R_T ** 3 - R_b ** 3)) / (V * R_s ** 3)
    r = np.concatenate([np.linspace(R_s / 4, R_s, 4),
                        np.linspace(R_s, R_b, 5)[1:],
                        np.linspace(R_b, R_T, 5)[1:]])
    rho = np.concatenate([np.full(4, rho_s), np.full(4, rho_o),
                          np.full(4, rho_i)])
    ph = np.concatenate([np.full(4, 50), np.zeros(4, int), np.ones(4, int)])
    return {
        'r_m': r, 'rho': rho, 'phases': ph,
        'K_Pa': np.full(12, 5e10), 'mu_Pa': np.full(12, 3e10),
        'T_K': np.linspace(260.0, 273.0, 12), 'P_MPa': np.linspace(0, 20, 12),
        'omega': OMEGA, 'eccentricity': ECC, 'Mtot_kg': M, 'R_body_m': R_T,
        'n_layers': 3, 'changeIndices': np.array([0, 4, 8, 12]),
        'layer_types': ['rock', 'ocean', 'ice'],
        'region_phases': ['Sil', 'Ocean', 'Ih'],
        'Tb_K': 270.0, 'D_hsphere_km': zb_km + d_ocean_km,
        'D_iceIh_km': zb_km, 'D_ocean_km': d_ocean_km,
        'CMR2': np.nan,
    }


def _write_cache(path):
    structs = [_node(zb, w) for zb in ZB_GRID for w in W_GRID]
    cache = {'schema_version': 'v3.1-zbw', 'zb_km_grid': ZB_GRID,
             'wOcean_ppt_grid': W_GRID, 'structures': structs}
    with open(path, 'wb') as f:
        pickle.dump(cache, f)
    return cache


def _make_cfg_dict(cache_path, **sampler):
    d = json.load(open(CFG))
    d['structure_cache_path'] = str(cache_path)
    d['sampler_settings'] = dict(d['sampler_settings'])
    d['sampler_settings'].update(sampler)
    return d


@pytest.fixture(scope='module')
def cache_path(tmp_path_factory):
    p = tmp_path_factory.mktemp('zbwcache') / 'syn_zbw.pkl'
    _write_cache(p)
    return p


def _runner(cache_path, cfg_dict=None, tmp_path=None):
    d = cfg_dict or _make_cfg_dict(cache_path)
    cp = Path(tmp_path or cache_path.parent) / 'cfg.json'
    cp.write_text(json.dumps(d))
    return MCMCRunner(InferenceConfig.from_json(str(cp)))


def _theta(runner):
    """Interior theta (in sampler order) built from the config midpoints."""
    vals = {'zb_km': 22.0, 'log10_wOcean_ppt': 0.3, 'rho_ice_kgm3': 925.0,
            'rho_rock_kgm3': 2400.0, 'compensation_C2': 0.6,
            'libration_sys_frac': 0.0}
    return np.array([vals[n] for n in runner.param_names])


def test_gate_definition(cache_path):
    r = _runner(cache_path)
    assert r._tidal_solve_needed() is False
    d = _make_cfg_dict(cache_path)
    d['observables']['Re_k2'] = [0.02, 0.01]
    r2 = _runner(cache_path, d)
    assert r2._tidal_solve_needed() is True
    d = _make_cfg_dict(cache_path)
    d['gravity_forward_model'] = 'clairaut_hydrostatic'
    r3 = _runner(cache_path, d)
    assert r3._tidal_solve_needed() is True


def test_a_loglike_finite_and_matches_independent_chi2(cache_path):
    r = _runner(cache_path)
    th = _theta(r)
    ll = r.log_likelihood_fn(th)
    assert np.isfinite(ll) and ll > -1e29

    td = r._expand_theta(th)
    g = r._derive_gravity_isostatic(td)
    lib = r._derive_libration_deg(td)
    assert g is not None and lib is not None
    obs = r.config.observables
    sm = r.config.metadata['sigma_model_add']
    s = {k: np.hypot(obs[k][1], sm[k]) for k in obs}
    rho = float(r.config.metadata['observable_correlations']['C20,C22'])
    pred = {'C20': g[(2, 0)], 'C22': g[(2, 2)], 'C30': g[(3, 0)],
            'libration_deg': lib}
    res = np.array([pred[k] - obs[k][0]
                    for k in ('C20', 'C22', 'C30', 'libration_deg')])
    sg = np.array([s[k] for k in ('C20', 'C22', 'C30', 'libration_deg')])
    cov = np.diag(sg ** 2)
    cov[0, 1] = cov[1, 0] = rho * sg[0] * sg[1]
    chi2 = float(res @ np.linalg.solve(cov, res))
    # Runner convention: log-likelihood = -0.5 * chi2 (no normalisation const)
    assert ll == pytest.approx(-0.5 * chi2, rel=1e-9)


def test_b_unpatched_path_would_raise_keyerror(cache_path):
    """The k2 solver on a production-schema zbw node raises KeyError
    ('eta_Pa_base') -- this is what the unpatched likelihood hit."""
    from PlanetProfile.Inference.forward_models import forward_model_k2_flexible
    r = _runner(cache_path)
    td = r._expand_theta(_theta(r))
    with pytest.raises(KeyError):
        forward_model_k2_flexible(td, r.structure_data, return_heating=False,
                                  arrhenius_params=r.arrhenius_params)


def test_d_gate_negative_k2_observable_still_calls_solver(cache_path, monkeypatch):
    import PlanetProfile.Inference.forward_models as fm
    calls = []

    def fake(theta_dict, structure_data, return_heating=False,
             arrhenius_params=None):
        calls.append(1)
        return 0.02, -0.001, 0.1, -0.01, None

    monkeypatch.setattr(fm, 'forward_model_k2_flexible', fake)
    d = _make_cfg_dict(cache_path)
    d['observables']['Re_k2'] = [0.02, 0.01]
    r = _runner(cache_path, d)
    ll = r.log_likelihood_fn(_theta(r))
    assert calls, 'k2 solver must still be called when a tidal observable exists'
    assert np.isfinite(ll)
    # NaN k2 must still hard-reject for such a config
    monkeypatch.setattr(fm, 'forward_model_k2_flexible',
                        lambda *a, **k: (np.nan, np.nan, np.nan, np.nan, None))
    r = _runner(cache_path, d)   # closure binds the solver at build time
    assert r.log_likelihood_fn(_theta(r)) == -1e30
    # and for the no-tidal config the solver is never called
    calls.clear()
    monkeypatch.setattr(fm, 'forward_model_k2_flexible', fake)
    r0 = _runner(cache_path)
    assert r0.log_likelihood_fn(_theta(r0)) > -1e29
    assert not calls


def test_sbi_dataset_generator_skips_k2(cache_path, monkeypatch):
    import PlanetProfile.Inference.forward_models as fm
    monkeypatch.setattr(fm, 'forward_model_k2_flexible',
                        lambda *a, **k: (_ for _ in ()).throw(
                            AssertionError('k2 solver called')))
    r = _runner(cache_path)
    theta, x = r.generate_sbi_dataset(
        n_samples=40, seed=3, apply_support_guard=True,
        imag_convention='abs', drop_nonfinite=True)
    assert x.shape[1] == 4 and len(x) > 0
    assert np.all(np.isfinite(x))
    assert r.last_sbi_dataset_stats['n_requested'] == 40


def test_c_end_to_end_run_inference_cli(cache_path, tmp_path, monkeypatch):
    import pocomc as pc
    from PlanetProfile.Inference import run_inference_cli
    # MCMCRunner.run hard-codes sampler.run(n_total=4096); at n_active=50
    # that is ~80 flow-training iterations. Cap it for the test (sampler
    # budget only; likelihood/derived-quantity paths are unchanged).
    _orig_run = pc.Sampler.run

    def _capped_run(self, *a, **k):
        k['n_total'] = min(int(k.get('n_total', 300)), 300)
        return _orig_run(self, *a, **k)

    monkeypatch.setattr(pc.Sampler, 'run', _capped_run)
    d = _make_cfg_dict(cache_path, n_effective=100, n_active=50,
                       n_posterior_samples=64, n_reeval=4)
    # Pipeline exercise, not a physics posterior: inflate the observation
    # sigmas so the toy posterior is broad and pocoMC finishes in seconds
    # (production sigmas need many tempering steps on this toy cache).
    d['observables'] = {k: [v[0], 100.0 * v[1]]
                        for k, v in d['observables'].items()}
    cfg = tmp_path / 'cfg.json'
    cfg.write_text(json.dumps(d))
    out = tmp_path / 'res.pkl'
    argv = ['run_inference_cli', '--config', str(cfg), '--output', str(out)]
    old = sys.argv
    sys.argv = argv
    t0 = time.time()
    try:
        rc = run_inference_cli.main()
    finally:
        sys.argv = old
    dt = time.time() - t0
    assert rc in (0, None)
    with open(out, 'rb') as f:
        res = pickle.load(f)
    n = len(res.samples)
    assert n >= 32
    assert res.k2_results.shape == (n, 2) and np.all(np.isnan(res.k2_results))
    assert np.all(np.isnan(res.h2_results))
    assert np.all(np.isfinite(res.log_likelihoods))
    iso = res.metadata['isostatic_derived']
    for k in ('C20', 'C22', 'C30', 'cmr2_derived', 'R_s_m', 'rho_s'):
        assert iso[k].shape == (n,)
        assert np.all(np.isfinite(iso[k])), k
    assert np.all((iso['cmr2_derived'] > 0.1) & (iso['cmr2_derived'] < 0.4))
    assert len(res.heating_results) == 0
    assert res.metadata['tidal_solve_skipped'] is True
    assert np.all(np.abs(iso['mass_residual_after']) < 1e-12)

    # Per-draw alignment: -0.5*chi2 rebuilt from the STORED diagnostics
    # (+ libration at the stored theta) must equal the stored log-like.
    r = _runner(cache_path, d, tmp_path)
    obs = r.config.observables
    sm = r.config.metadata['sigma_model_add']
    names = ('C20', 'C22', 'C30', 'libration_deg')
    sg = np.array([np.hypot(obs[k][1], sm[k]) for k in names])
    rho = float(r.config.metadata['observable_correlations']['C20,C22'])
    cov = np.diag(sg ** 2)
    cov[0, 1] = cov[1, 0] = rho * sg[0] * sg[1]
    for j in range(n):
        lib = r._derive_libration_deg(r._expand_theta(res.samples[j]))
        pred = np.array([iso['C20'][j], iso['C22'][j], iso['C30'][j], lib])
        resid = pred - np.array([obs[k][0] for k in names])
        chi2 = float(resid @ np.linalg.solve(cov, resid))
        assert -0.5 * chi2 == pytest.approx(res.log_likelihoods[j], rel=1e-9)
    print(f'e2e wall time {dt:.1f}s for {n} draws')
