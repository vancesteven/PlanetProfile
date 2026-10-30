"""
Smoke tests for CMR2 (moment-of-inertia) observable in MCMC likelihood.

v1 path (cached CMR² lookup, Titan Test50):
  1. CMR2 chi² term fires correctly with a mock structure dict
  2. Graceful NaN fallback when CMR2 key is absent (old cache)
  3. Grid Tb_K nearest-neighbor lookup returns the correct CMR2
  4. cmr2_results wired into InferenceResult

v2 path (Phase C1 mass-conservation: derived rho_sil + sampled core,
recomputed CMR² per sample):
  5. Mass-conservation roundtrip recovers ρ_sil and matches direct CMR²
  6. ρ_sil out of bounds → log-likelihood floored at -1e30 (rejection)
  7. Missing R_core_km / rho_core_kgm3 in theta_dict → rejection
"""
import sys
import os
import numpy as np
import pytest

# Ensure project root on path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))


# ---------------------------------------------------------------------------
# Helpers — build minimal structure dicts without touching PlanetProfile
# ---------------------------------------------------------------------------

def _minimal_structure(cmr2=0.343):
    """Flat single-layer structure sufficient for v1 likelihood tests."""
    n = 20
    return {
        'r_m': np.linspace(0, 1e6, n),
        'rho': np.ones(n) * 1000.0,
        'K_Pa': np.ones(n) * 1e9,
        'mu_Pa': np.ones(n) * 1e9,
        'eta_Pa_base': np.ones(n) * 1e14,
        'phases': np.zeros(n),
        'bulk_visc': np.zeros(n),
        'changeIndices': np.array([0, n]),
        'n_layers': 1,
        'layer_upper_radii': (1e6,),
        'layer_types': ('solid',),
        'region_phases': ['Ih'],
        'omega': 4.56e-6,
        'eccentricity': 0.0292,
        'host_mass': 5.683e26,
        'a_m': 1.22e9,
        'phase_map': {0: '0', 1: 'Ih'},
        'T_K': None,
        'bodyname': 'TestBody',
        'CMR2': cmr2,
    }


def _grid_structure(tb_values, cmr2_values):
    """Grid cache dict: Tb_K -> structure."""
    grid = {}
    for tb, cmr2 in zip(tb_values, cmr2_values):
        s = _minimal_structure(cmr2=cmr2)
        grid[float(tb)] = s
    grid_Tb_values = np.array(sorted(grid.keys()))
    return {'grid_cache': grid, 'grid_Tb_values': grid_Tb_values}


def _v2_two_layer_structure(R_body_m, R_oceanbot_m, rho_hydro=1100.0, n_hydro=8):
    """Build a hydrosphere-only cache for v2 tests.

    Layout: ``n_hydro`` cells of constant rho_hydro from R_oceanbot_m to
    R_body_m, with one inner sentinel silicate cell (phase=50) below
    R_oceanbot_m so the hydrosphere mask has a clean lower edge to
    detect.
    """
    # Inner silicate sentinel runs from 0 to R_oceanbot_m (one cell).
    # Hydrosphere cells run from R_oceanbot_m to R_body_m (n_hydro cells).
    r_outer_hydro = np.linspace(R_oceanbot_m, R_body_m, n_hydro + 1)[1:]
    r_m = np.concatenate(([R_oceanbot_m], r_outer_hydro))  # ascending
    # Mark phase=50 for the silicate cell, phase=1 (Ih) for hydrosphere cells.
    phases = np.concatenate(([50], np.ones(n_hydro, dtype=int)))
    rho = np.concatenate(([3300.0], np.full(n_hydro, rho_hydro)))
    return {
        'r_m': r_m,
        'rho': rho,
        'phases': phases,
        'R_body_m': R_body_m,
        # Filled in by callers based on the v2 mass-conservation target
        'Mtot_kg': np.nan,
    }


class _StubConfig:
    """Minimal substitute for ``InferenceConfig`` — only the attributes the
    log-likelihood closure actually reads."""

    def __init__(self, sampler_settings=None, derived_params=None):
        self.sampler_settings = sampler_settings or {}
        self.derived_params = derived_params or {}


def _make_ll(
    param_names,
    observables,
    structure_data,
    *,
    param_groups=None,
    fixed_params=None,
    sampler_settings=None,
    derived_params=None,
):
    """Build the log-likelihood closure with a stub that mimics MCMCRunner's
    attribute surface as of Phase C1 Stage 2.

    The runner reads ``self.param_names``, ``self.param_groups``,
    ``self.fixed_params``, ``self.config.sampler_settings``, and
    ``self.config.derived_params`` during closure construction.
    """
    from PlanetProfile.Inference.mcmc_runner import MCMCRunner

    class _Stub:
        pass

    stub = _Stub()
    stub.param_names = param_names
    stub.param_groups = param_groups or {}
    stub.fixed_params = fixed_params or {}
    stub.config = _StubConfig(
        sampler_settings=sampler_settings,
        derived_params=derived_params,
    )
    return MCMCRunner._make_flexible_log_likelihood(stub, observables, structure_data)


@pytest.fixture
def mock_forward_model():
    """Patch forward_model_k2_flexible to return constant (Re_k2, Im_k2).

    Used so v1/v2 likelihood tests can isolate the CMR² chi² term.
    """
    import PlanetProfile.Inference.forward_models as fm
    _orig = fm.forward_model_k2_flexible

    def _mock_fwd(theta_dict, structure_data, **kw):
        return 0.608, -0.135, None

    fm.forward_model_k2_flexible = _mock_fwd
    try:
        yield _mock_fwd
    finally:
        fm.forward_model_k2_flexible = _orig


# ---------------------------------------------------------------------------
# v1 path: cached CMR² lookup
# ---------------------------------------------------------------------------

def test_cmr2_chi2_flat_structure(mock_forward_model):
    """Likelihood decreases when CMR2 deviates from observed value."""
    obs_val, obs_err = 0.343, 0.001
    observables = {'CMR2': (obs_val, obs_err)}

    structure_at_match = _minimal_structure(cmr2=obs_val)
    structure_off = _minimal_structure(cmr2=obs_val + 5 * obs_err)

    fn_match = _make_ll([], observables, structure_at_match)
    fn_off = _make_ll([], observables, structure_off)

    ll_match = fn_match(np.array([]))
    ll_off = fn_off(np.array([]))

    assert ll_match > ll_off, (
        f"Expected higher likelihood at CMR2={obs_val} than at {obs_val + 5*obs_err}, "
        f"got {ll_match} vs {ll_off}"
    )


def test_cmr2_nan_fallback(mock_forward_model):
    """No crash and no chi² contribution when structure cache lacks CMR2."""
    structure_no_cmr2 = _minimal_structure()
    del structure_no_cmr2['CMR2']   # simulate old cache

    observables = {'CMR2': (0.343, 0.001)}
    fn = _make_ll([], observables, structure_no_cmr2)
    ll = fn(np.array([]))

    # With only CMR2 observable and NaN CMR2, chi2=0 → ll = 0.0
    assert ll == 0.0, f"Expected 0.0 log-likelihood with NaN CMR2, got {ll}"


def test_cmr2_grid_lookup(mock_forward_model):
    """Grid lookup selects CMR2 from the nearest Tb_K grid point."""
    tb_values = [250.0, 260.0, 270.0]
    cmr2_values = [0.340, 0.343, 0.346]
    structure_data = _grid_structure(tb_values, cmr2_values)

    observables = {'CMR2': (0.343, 0.001)}
    fn = _make_ll(['Tb_K'], observables, structure_data)

    # Tb_K=261 → nearest grid point is 260 → CMR2=0.343 → chi2≈0
    ll_near_match = fn(np.array([261.0]))
    # Tb_K=251 → nearest grid point is 250 → CMR2=0.340 → chi2>0
    ll_near_off = fn(np.array([251.0]))

    assert ll_near_match > ll_near_off, (
        f"Expected higher LL near CMR2 match, got {ll_near_match} vs {ll_near_off}"
    )


# ---------------------------------------------------------------------------
# v2 path: derived ρ_sil + sampled core, recomputed CMR² per sample
# ---------------------------------------------------------------------------


def _v2_target_mass(R_body_m, R_oceanbot_m, rho_hydro, R_core_m, rho_core, rho_sil):
    """Total mass of an idealized two-shell-interior + uniform-hydrosphere body."""
    M_core = (4.0 / 3.0) * np.pi * R_core_m ** 3 * rho_core
    M_sil = (
        (4.0 / 3.0) * np.pi * (R_oceanbot_m ** 3 - R_core_m ** 3) * rho_sil
    )
    M_hydro = (
        (4.0 / 3.0) * np.pi * (R_body_m ** 3 - R_oceanbot_m ** 3) * rho_hydro
    )
    return M_core + M_sil + M_hydro


def test_v2_mass_conservation_recovers_cmr2(mock_forward_model):
    """End-to-end v2 path: build a body with known geometry, set Mtot_kg
    accordingly, sample (R_core, rho_core), expect CMR² to match the direct
    integral over the assembled body."""
    R_body_m = 2.5e6
    R_oceanbot_m = 0.7 * R_body_m
    rho_hydro = 1100.0
    R_core_km = 1000.0          # 1000 km core
    rho_core = 7500.0
    rho_sil_true = 3200.0

    M_total = _v2_target_mass(
        R_body_m, R_oceanbot_m, rho_hydro,
        R_core_km * 1000.0, rho_core, rho_sil_true,
    )

    structure = _v2_two_layer_structure(R_body_m, R_oceanbot_m, rho_hydro=rho_hydro)
    structure['Mtot_kg'] = M_total

    # Direct CMR² for the full body (the value the v2 path must reproduce)
    from PlanetProfile.Inference.structure_derivation import compute_cmr2
    direct_layers = [
        (0.0, R_core_km * 1000.0, rho_core),
        (R_core_km * 1000.0, R_oceanbot_m, rho_sil_true),
        (R_oceanbot_m, R_body_m, rho_hydro),
    ]
    cmr2_direct = compute_cmr2(direct_layers, R_body_m=R_body_m, M_body_kg=M_total)

    observables = {'CMR2': (cmr2_direct, 1e-6)}
    derived_params = {
        'rho_sil_kgm3': {
            'derivation': 'mass_conservation',
            'bounds': [2200.0, 3500.0],
            'reject_if_outside_bounds': True,
        },
    }

    fn = _make_ll(
        ['R_core_km', 'rho_core_kgm3'],
        observables, structure,
        derived_params=derived_params,
    )
    ll = fn(np.array([R_core_km, rho_core]))

    # CMR² match → chi² ≈ 0 → ll ≈ 0
    assert ll > -0.5, (
        f"v2 mass-conservation should recover CMR² to obs precision; "
        f"got ll={ll!r} (expected ≈ 0)."
    )


def test_v2_rho_sil_out_of_bounds_rejects(mock_forward_model):
    """Sample with implausible rho_sil → log-likelihood floored at -1e30."""
    R_body_m = 2.5e6
    R_oceanbot_m = 0.7 * R_body_m
    rho_hydro = 1100.0

    # Construct M_total assuming rho_sil = 3200 (in-bounds), but the SAMPLED
    # core will displace mass such that the derived rho_sil falls below 2200.
    rho_sil_true = 3200.0
    R_core_in_bounds_km = 500.0
    rho_core = 7500.0
    M_total = _v2_target_mass(
        R_body_m, R_oceanbot_m, rho_hydro,
        R_core_in_bounds_km * 1000.0, rho_core, rho_sil_true,
    )

    structure = _v2_two_layer_structure(R_body_m, R_oceanbot_m, rho_hydro=rho_hydro)
    structure['Mtot_kg'] = M_total

    observables = {'CMR2': (0.34, 0.005)}
    derived_params = {
        'rho_sil_kgm3': {
            'derivation': 'mass_conservation',
            'bounds': [2200.0, 3500.0],
            'reject_if_outside_bounds': True,
        },
    }

    fn = _make_ll(
        ['R_core_km', 'rho_core_kgm3'],
        observables, structure,
        derived_params=derived_params,
    )
    # Sample with an absurdly large dense core so mass-conservation drives
    # rho_sil far below 2200 (rejected).
    ll = fn(np.array([1500.0, 8000.0]))
    assert ll == -1e30, f"Expected -1e30 on bounds rejection, got {ll}"


def test_v2_rejects_when_core_params_missing_in_theta(mock_forward_model):
    """v2 path requires both R_core_km and rho_core_kgm3 in theta_dict."""
    R_body_m = 2.5e6
    R_oceanbot_m = 0.7 * R_body_m
    structure = _v2_two_layer_structure(R_body_m, R_oceanbot_m)
    structure['Mtot_kg'] = 1.0e23

    observables = {'CMR2': (0.34, 0.005)}
    derived_params = {
        'rho_sil_kgm3': {
            'derivation': 'mass_conservation',
            'bounds': [2200.0, 3500.0],
            'reject_if_outside_bounds': True,
        },
    }

    # Only R_core_km in param_names — rho_core_kgm3 is missing
    fn = _make_ll(
        ['R_core_km'],
        observables, structure,
        derived_params=derived_params,
    )
    ll = fn(np.array([500.0]))
    assert ll == -1e30, (
        f"Expected -1e30 when rho_core_kgm3 is absent from theta, got {ll}"
    )


# ---------------------------------------------------------------------------
# Test 4: cmr2_results field exists in InferenceResult and is correct shape
# ---------------------------------------------------------------------------

def test_cmr2_results_in_inference_result():
    """InferenceResult accepts and stores cmr2_results array."""
    from PlanetProfile.Inference.inference_core import InferenceResult, InferenceConfig

    config = InferenceConfig(
        mode='mcmc',
        bodyname='TestBody',
        param_space={'log10_eta_Ih': {'prior_type': 'uniform', 'bounds': [12, 16]}},
        observables={'Re_k2': (0.608, 0.048)},
        sampler_settings={'n_effective': 100},
    )

    n = 10
    samples = np.random.rand(n, 1)
    log_likes = np.random.rand(n)
    k2_results = np.random.rand(n, 2)
    cmr2_results = np.full(n, 0.343)

    result = InferenceResult(
        config=config,
        samples=samples,
        log_likelihoods=log_likes,
        param_names=['log10_eta_Ih'],
        param_labels=[r'$\log_{10}(\eta_{\rm Ih})$'],
        k2_results=k2_results,
        cmr2_results=cmr2_results,
    )

    assert result.cmr2_results is not None
    assert result.cmr2_results.shape == (n,)
    np.testing.assert_array_equal(result.cmr2_results, cmr2_results)


# ---------------------------------------------------------------------------
# Allow running directly
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    test_cmr2_chi2_flat_structure()
    print("Test 1 passed: CMR2 chi² fires correctly")

    test_cmr2_nan_fallback()
    print("Test 2 passed: NaN fallback graceful")

    test_cmr2_grid_lookup()
    print("Test 3 passed: Grid Tb_K nearest-neighbor correct")

    test_cmr2_results_in_inference_result()
    print("Test 4 passed: cmr2_results wired into InferenceResult")

    print("\nAll CMR2 smoke tests passed.")
