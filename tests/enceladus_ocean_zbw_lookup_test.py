"""E2 (r6 training blocker): the OCEAN-branch (zb, w) cache lookup.

The zb x w cache (schema ``v3.1-zbw`` / ``v3.2-zbw-joint``) got a FROZEN-axis
lookup when A8 landed. Its OCEAN axis never got one, so
``MCMCRunner._struct_for_hydrosphere`` fell through every ocean sample to
``return sd`` and handed back the WHOLE CACHE DICT as if it were a structure.
Nothing downstream finds ``r_m`` in that dict, so every ocean sample would
have rejected to ``-inf``: no training run could have produced a posterior.
r6's D2 check missed it because it ran on a Tb-keyed smoke cache, which is
served by an entirely different branch of the dispatch.

The repair mirrors ``_frozen_struct_for_theta``: bilinear in
``(zb, log10 w)``, and ANY ``None`` bracketing corner rejects the sample with
no fallback.

The last test here is the load-bearing one: the D3-D6 ruled zb bracket
(25.99 / 27.34 km) must still reproduce THROUGH THE PRODUCTION DISPATCH when
the structures arrive through the new lookup rather than as a bare struct.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import brentq

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402
from PlanetProfile.Inference.mcmc_runner import MCMCRunner  # noqa: E402

CANDIDATE_CONFIG = (REPO / 'PlanetProfile' / 'Inference' / 'configs'
                    / 'enceladus_cassini_isostasy_7D.json')

# B2' fiducial, identical to tests/enceladus_libration_wiring_test.py.
R_T = 252.22e3
M_ENC = 1.08022e20
OMEGA = 5.307e-5
ECC = 0.0047
RHO_ICE = 925.0
RHO_OCEAN = 1005.0
D_OCEAN_KM = 36.0
SIGMA_LIB_DEG = 0.003
LIB_OBS_DEG = 0.092
H_OBS = {(2, 0): -3510.0, (2, 2): 857.0, (3, 0): 420.0}


def _b2prime_struct(zb_km):
    """Mass-conserved 3-layer B2' structure as a CACHE NODE."""
    R_b = R_T - zb_km * 1e3
    R_c = R_b - D_OCEAN_KM * 1e3
    V = 4.0 * np.pi / 3.0
    m_shell = RHO_ICE * V * (R_T ** 3 - R_b ** 3)
    m_ocean = RHO_OCEAN * V * (R_b ** 3 - R_c ** 3)
    rho_core = (M_ENC - m_shell - m_ocean) / (V * R_c ** 3)
    return {
        'r_m': np.array([R_c, R_b, R_T]),
        'rho': np.array([rho_core, RHO_OCEAN, RHO_ICE]),
        'phases': np.array([50, 0, 1]),
        'omega': OMEGA,
        'eccentricity': ECC,
        'Mtot_kg': M_ENC,
        'zb_km_node': float(zb_km),
        # The layer bookkeeping a real cache node carries, so the shared
        # _blend_b_layered interpolant runs exactly as it does in production.
        'changeIndices': np.array([0, 1, 2, 3]),
        'n_layers': 3,
        'region_phases': ['Sil', 'water', 'Ih'],
        'layer_types': ['solid', 'liquid', 'solid'],
        'layer_upper_radii': [R_c, R_b, R_T],
        'R_body_m': R_T,
    }


def _zbw_cache(zb_grid, w_grid, holes=()):
    """A v3.1-zbw cache of B2' nodes. ``holes`` are (i_zb, i_w) build rejects."""
    structs = []
    for i, zb in enumerate(zb_grid):
        for j, _w in enumerate(w_grid):
            structs.append(None if (i, j) in holes else _b2prime_struct(zb))
    return {
        'zb_km_grid': np.asarray(zb_grid, dtype=float),
        'wOcean_ppt_grid': np.asarray(w_grid, dtype=float),
        'structures': structs,
        'schema_version': 'v3.1-zbw',
        'ocean_comp': 'Seawater',
        'zb_tol_km': 0.125,
        'n_zb_placement_rejected': 0,
        'frozen_branch_supported': False,
    }


def _config_libration_block(**overrides):
    with open(CANDIDATE_CONFIG) as fh:
        blk = dict(json.load(fh)['metadata']['libration_model_correction'])
    blk.update(overrides)
    return blk


_UNCORRECTED = _config_libration_block(multiplicative_factor=1.0)


def _runner(structure_data, libration_block=None):
    metadata = {'isostasy': {
        'H_obs_lm_m': {'2,0': H_OBS[(2, 0)], '2,2': H_OBS[(2, 2)],
                       '3,0': H_OBS[(3, 0)]},
        'shape_ref_radius_m': R_T,
        'gravity_ref_radius_m': R_T,
        'finite_amplitude': True,
    }}
    if libration_block is not None:
        metadata['libration_model_correction'] = libration_block
    config = InferenceConfig(
        mode='mcmc', bodyname='Enceladus',
        param_space={'compensation_C2': {'prior_type': 'uniform',
                                         'bounds': [0.0, 1.0]}},
        observables={'libration_deg': (LIB_OBS_DEG, SIGMA_LIB_DEG)},
        sampler_settings={},
        gravity_forward_model='isostatic_hm2019',
        metadata=metadata,
    )
    runner = object.__new__(MCMCRunner)
    runner.config = config
    runner.structure_data = structure_data
    return runner


# Production-shaped axes: the config's finest zb segment step (0.25 km) over
# the posterior region, and the declared 40-point log10 w axis restricted to
# the same spacing (0.0769 dex).
ZB_GRID = list(np.arange(22.0, 30.0001, 0.25))
W_GRID = list(10.0 ** np.linspace(-1.0, 2.0, 40))


# ---------------------------------------------------------------------------
# THE DEFECT
# ---------------------------------------------------------------------------

def test_the_cache_dict_is_not_a_structure():
    """What the fall-through returned. Pinned so the defect cannot come back
    disguised as 'it returns something'."""
    cache = _zbw_cache(ZB_GRID, W_GRID)
    assert 'r_m' not in cache
    assert 'structures' in cache


def test_an_ocean_theta_gets_a_structure_not_the_cache():
    runner = _runner(_zbw_cache(ZB_GRID, W_GRID))
    out = runner._struct_for_hydrosphere(
        {'zb_km': 25.6, 'log10_wOcean_ppt': 0.7})
    assert out is not None
    assert isinstance(out, dict)
    assert 'r_m' in out and 'rho' in out and 'phases' in out
    assert np.asarray(out['r_m']).size == 3
    # and it is NOT the cache
    assert 'structures' not in out
    assert 'zb_km_grid' not in out
    assert out['zb_km_sampled'] == pytest.approx(25.6)
    assert out['wOcean_ppt_sampled'] == pytest.approx(10.0 ** 0.7)


def test_derive_libration_is_finite_through_the_new_lookup():
    runner = _runner(_zbw_cache(ZB_GRID, W_GRID), _UNCORRECTED)
    lib = runner._derive_libration_deg(
        {'zb_km': 25.6, 'log10_wOcean_ppt': 0.7, 'compensation_C2': 0.5})
    assert lib is not None
    assert np.isfinite(lib)
    assert 0.05 < lib < 0.15


# ---------------------------------------------------------------------------
# Interpolation behaviour
# ---------------------------------------------------------------------------

def test_an_exact_node_is_served_bit_for_bit():
    """At a node the bilinear weights are exactly 0/1, so no blend runs."""
    cache = _zbw_cache(ZB_GRID, W_GRID)
    runner = _runner(cache)
    i, j = 12, 20
    out = runner._struct_for_hydrosphere(
        {'zb_km': ZB_GRID[i], 'log10_wOcean_ppt': float(np.log10(W_GRID[j]))})
    node = cache['structures'][i * len(W_GRID) + j]
    assert np.array_equal(out['r_m'], node['r_m'])
    assert np.array_equal(out['rho'], node['rho'])


def test_interpolation_is_monotone_and_brackets_its_nodes():
    cache = _zbw_cache(ZB_GRID, W_GRID)
    runner = _runner(cache)

    def _R_b(zb):
        return runner._struct_for_hydrosphere(
            {'zb_km': zb, 'log10_wOcean_ppt': 0.7})['r_m'][1]

    lo, mid, hi = _R_b(25.0), _R_b(25.125), _R_b(25.25)
    assert lo > mid > hi                       # thicker shell -> smaller R_b
    assert mid == pytest.approx(0.5 * (lo + hi), rel=1e-12)


def test_w_axis_is_weighted_in_log10_not_ppt(monkeypatch):
    """The prior is uniform in log10 w and the grid is log-spaced, so the
    interpolant must weight in log10 w. Probed through the corner weights."""
    from PlanetProfile.Inference import grid_interp_2d
    w_grid = np.array([1.0, 100.0])
    corners, weights = grid_interp_2d.bilinear_weights(
        np.array([25.0, 26.0]), w_grid, 25.0, 10.0)
    # 10 ppt is the log-midpoint of [1, 100] and the ppt-value 50.5 is not.
    by_corner = {}
    for c, wt in zip(corners, weights):
        by_corner[c] = by_corner.get(c, 0.0) + wt
    assert by_corner[0] == pytest.approx(0.5)
    assert by_corner[1] == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# Reject rules -- no fallback
# ---------------------------------------------------------------------------

def test_a_single_none_corner_rejects_the_sample():
    """The ruled rule for this campaign: an unbuilt bracketing node means the
    sample's own coordinates were never built. Serving its neighbours would
    be a shell-thickness snap on the sole thickness channel."""
    cache = _zbw_cache(ZB_GRID, W_GRID, holes={(12, 20)})
    runner = _runner(cache)
    # a query bracketed by the hole
    assert runner._struct_for_hydrosphere(
        {'zb_km': ZB_GRID[12] + 0.1,
         'log10_wOcean_ppt': float(np.log10(W_GRID[20])) + 0.01}) is None
    # a query whose four corners are all built is unaffected
    assert runner._struct_for_hydrosphere(
        {'zb_km': ZB_GRID[2] + 0.1, 'log10_wOcean_ppt': 0.0}) is not None


def test_out_of_support_rejects_on_both_axes():
    runner = _runner(_zbw_cache(ZB_GRID, W_GRID))
    assert runner._struct_for_hydrosphere(
        {'zb_km': 21.0, 'log10_wOcean_ppt': 0.7}) is None      # below zb grid
    assert runner._struct_for_hydrosphere(
        {'zb_km': 31.0, 'log10_wOcean_ppt': 0.7}) is None      # above zb grid
    assert runner._struct_for_hydrosphere(
        {'zb_km': 25.0, 'log10_wOcean_ppt': 2.5}) is None      # above w grid
    assert runner._struct_for_hydrosphere(
        {'zb_km': 25.0, 'log10_wOcean_ppt': -1.5}) is None     # below w grid


def test_box_edges_are_servable():
    """Out-of-support must not eat the prior's own edges on float noise."""
    runner = _runner(_zbw_cache(ZB_GRID, W_GRID))
    for zb in (ZB_GRID[0], ZB_GRID[-1]):
        for lw in (-1.0, 2.0):
            assert runner._struct_for_hydrosphere(
                {'zb_km': zb, 'log10_wOcean_ppt': lw}) is not None


def test_a_frozen_theta_is_not_served_by_the_ocean_lookup():
    """Branch dispatch order: rho_rock_kgm3 routes to the frozen lookup even
    on a joint cache, and never falls into the ocean interpolant."""
    cache = _zbw_cache(ZB_GRID, W_GRID)
    cache['frozen_branch_supported'] = True
    cache['frozen_zb_km_grid'] = np.array([])
    cache['frozen_structures'] = []
    runner = _runner(cache)
    assert runner._struct_for_hydrosphere(
        {'rho_rock_kgm3': 2400.0, 'rho_ice_kgm3': 925.0}) is None


def test_a_tb_keyed_cache_is_untouched():
    """The v3.0 (Tb, w) path must not be captured by the new branch."""
    runner = _runner({'Tb_K_grid': np.array([270.0, 271.0]),
                      'wOcean_ppt_grid': np.array([1.0, 10.0]),
                      'structures': [_b2prime_struct(25.0)] * 4})
    out = runner._struct_for_hydrosphere(
        {'Tb_K': 270.5, 'log10_wOcean_ppt': 0.5})
    assert out is not None and 'r_m' in out


# ---------------------------------------------------------------------------
# THE LOAD-BEARING CHECK: D3-D6 still reproduces through the new lookup
# ---------------------------------------------------------------------------

def _zb_matching_observed_through_cache(C2, cache):
    """Shell thickness reproducing the Park libration, solved through
    _struct_for_hydrosphere -> _derive_libration_deg on the zb x w cache."""
    runner = _runner(cache, _UNCORRECTED)

    def _resid(zb_km):
        lib = runner._derive_libration_deg(
            {'zb_km': zb_km, 'log10_wOcean_ppt': 0.7,
             'compensation_C2': C2})
        assert lib is not None, f'sample rejected at zb={zb_km}'
        return lib - LIB_OBS_DEG

    return brentq(_resid, ZB_GRID[0] + 1e-9, ZB_GRID[-1] - 1e-9, xtol=1e-6)


@pytest.mark.parametrize('C2,target_zb_km', [(0.0, 27.34), (1.0, 25.99)])
def test_ruled_bracket_reproduces_through_the_zbw_lookup(C2, target_zb_km):
    """The D3-D6 ruled bracket, with the structures now arriving through the
    E2 lookup instead of as a bare struct. The residual against the
    bare-struct answer is INTERPOLATION error on a 0.25 km grid -- the very
    quantity B5 must report -- so it is asserted small, not zero."""
    cache = _zbw_cache(ZB_GRID, W_GRID)
    zb = _zb_matching_observed_through_cache(C2, cache)
    assert zb == pytest.approx(target_zb_km, abs=0.02)


def test_bracket_ordering_survives_the_lookup():
    cache = _zbw_cache(ZB_GRID, W_GRID)
    zb_c2_0 = _zb_matching_observed_through_cache(0.0, cache)
    zb_c2_1 = _zb_matching_observed_through_cache(1.0, cache)
    assert zb_c2_0 > zb_c2_1                      # C2=0 edge is the thicker
    assert 1.2 < zb_c2_0 - zb_c2_1 < 1.5          # ruled span 27.34 - 25.99


def test_interpolation_error_on_the_production_grid_is_reported():
    """B5's quantity, measured here on the 0.25 km segment: the zb the cache
    lookup returns vs the same solve on continuously-evaluated structures."""
    cache = _zbw_cache(ZB_GRID, W_GRID)
    zb_cache = _zb_matching_observed_through_cache(1.0, cache)

    runner_exact = _runner(_b2prime_struct(25.0), _UNCORRECTED)

    def _resid(zb_km):
        runner_exact.structure_data = _b2prime_struct(zb_km)
        return runner_exact._derive_libration_deg(
            {'compensation_C2': 1.0}) - LIB_OBS_DEG

    zb_exact = brentq(_resid, 22.0, 30.0, xtol=1e-6)
    err_km = abs(zb_cache - zb_exact)
    # 0.25 km nodes, ~1.1 sigma_obs per km: the interpolation error must stay
    # well inside a tenth of a node.
    assert err_km < 0.02, f'interpolation error {err_km:.4f} km'
