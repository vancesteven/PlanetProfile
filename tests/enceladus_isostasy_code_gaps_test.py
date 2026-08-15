"""Tests closing the Enceladus isostasy campaign's code_gaps list.

Mechanical wiring only -- no physics changes, no scientific judgment
calls (see validation_reports/enceladus_isostasy/RESUME_NOTE.md and
PlanetProfile/Inference/configs/enceladus_cassini_isostasy_7D.json,
metadata.blockers_open.code_gaps):

  1. parameter_registry round-trip for zb_km / compensation_C2 /
     rho_ice_kgm3 / libration_sys_frac.
  2. C30 observable dispatch (mcmc_runner.log_likelihood) returns finite
     values on a 3-layer synthetic structure, through the
     gravity_forward_model='isostatic_hm2019' path.
  3. isostatic_hm2019 dispatch (MCMCRunner._derive_gravity_isostatic)
     reproduces direct isostasy-module calls assembled independently
     from the same primitives, and is checked (loosely -- see the test
     docstring) against tests/isostasy_hm2019_gate_test.py::predict().
  4. B4 zb-axis cache-builder mode (cache_builder.build_zbw_grid_cache)
     smoke build on a tiny 3x2 grid (marked slow: real PlanetProfile
     runs, ~2-3 min).

Explicitly NOT covered here (out of scope; see RESUME_NOTE.md and the
docstrings in mcmc_runner.py / cache_builder.py): the frozen-branch
isostatic gravity / zb-derivation path, Sigma_model likelihood
inflation, and the MAJOR-1 per-node mass-conservation invariant (that
fix is reviewer-blocked production physics and must not land
unilaterally).
"""
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Inference.parameter_registry import (  # noqa: E402
    PARAMETER_REGISTRY, validate_parameter_combination)
from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402
from PlanetProfile.Inference.mcmc_runner import MCMCRunner  # noqa: E402
from PlanetProfile.Gravity.isostasy import (  # noqa: E402
    isostatic_gravity, finite_amplitude_coeffs, interface_gravity_coeff,
    eccentricities_to_axes, triaxial_to_H2m)
from PlanetProfile.Gravity.Librations import compute_eccentricities  # noqa: E402


# ===========================================================================
# Task 1 -- parameter registry round-trip
# ===========================================================================

def test_registry_zb_km_roundtrip():
    p = PARAMETER_REGISTRY['zb_km']
    assert p.category == 'structure'
    assert p.default_prior == 'uniform'
    assert p.default_bounds == [5.0, 45.0]
    assert p.units == 'km'
    assert p.requires_structure_rebuild is True


def test_registry_compensation_c2_roundtrip():
    p = PARAMETER_REGISTRY['compensation_C2']
    assert p.category == 'gravity'
    assert p.default_prior == 'uniform'
    assert p.default_bounds == [0.0, 1.0]
    assert p.requires_structure_rebuild is False


def test_registry_rho_ice_kgm3_roundtrip():
    p = PARAMETER_REGISTRY['rho_ice_kgm3']
    assert p.default_prior == 'uniform'
    assert p.default_bounds == [915.0, 935.0]
    assert p.units == 'kg/m^3'
    assert p.requires_structure_rebuild is False


def test_registry_libration_sys_frac_roundtrip():
    p = PARAMETER_REGISTRY['libration_sys_frac']
    assert p.default_prior == 'truncated_gaussian'
    assert p.default_bounds == [-0.012, 0.012]
    assert p.default_mean == 0.0
    assert p.default_std == 0.004
    assert p.requires_structure_rebuild is False


def test_registry_enceladus_5d_combination_valid():
    result = validate_parameter_combination(
        ['zb_km', 'compensation_C2', 'rho_ice_kgm3', 'libration_sys_frac'])
    assert result['valid'] is True
    assert result['requires_rebuild'] is True  # zb_km alone triggers it


# ===========================================================================
# Tasks 2/3 -- C30 dispatch + isostatic_hm2019 forward model
# ===========================================================================

# H&M 2019 Tajeddine-block inputs (unnormalized), matching both
# tests/isostasy_hm2019_gate_test.py and isostasy_test.py's fiducials.
R_REF = 252.22e3
M_ENC = 1.08022e20
OMEGA = 5.307e-5
RHO_ICE = 925.0
RHO_OCEAN = 1020.0
H_OBS = {(2, 0): -3510.0, (2, 2): 857.0, (3, 0): 420.0}
SIG_C = {(2, 0): 35.8e-6, (2, 2): 15.9e-6, (3, 0): 23.8e-6}


def _synthetic_ocean_struct(D_shell_km, D_ocean_km):
    """Idealized 3-layer (interior / ocean / shell) Enceladus-class
    structure, mass-conserved against M_ENC -- the same construction as
    tests/isostasy_hm2019_gate_test.py::predict()."""
    R_t = R_REF
    R_b = R_t - D_shell_km * 1e3
    R_c = R_b - D_ocean_km * 1e3
    V = 4.0 * np.pi / 3.0
    m_shell = RHO_ICE * V * (R_t ** 3 - R_b ** 3)
    m_ocean = RHO_OCEAN * V * (R_b ** 3 - R_c ** 3)
    rho_core = (M_ENC - m_shell - m_ocean) / (V * R_c ** 3)
    return {
        'r_m': np.array([R_c, R_b, R_t]),
        'rho': np.array([rho_core, RHO_OCEAN, RHO_ICE]),
        'phases': np.array([50, 0, 1]),  # interior(porous rock) / ocean / ice Ih
        'omega': OMEGA,
    }


def _isostatic_config(**metadata_extra):
    iso = {
        'H_obs_lm_m': {'2,0': H_OBS[(2, 0)], '2,2': H_OBS[(2, 2)],
                       '3,0': H_OBS[(3, 0)]},
        'shape_ref_radius_m': R_REF,
        'gravity_ref_radius_m': R_REF,  # no B14 rescale in this fixture
        'finite_amplitude': True,
    }
    metadata = {'isostasy': iso}
    metadata.update(metadata_extra)
    return InferenceConfig(
        mode='mcmc', bodyname='Enceladus',
        param_space={'compensation_C2': {'prior_type': 'uniform',
                                         'bounds': [0.0, 1.0]}},
        observables={'C20': (-5521e-6, 35.8e-6), 'C22': (1574e-6, 15.9e-6),
                     'C30': (118e-6, 23.8e-6)},
        sampler_settings={},
        gravity_forward_model='isostatic_hm2019',
        metadata=metadata,
    )


def _bare_runner(config, structure_data):
    """MCMCRunner instance bypassing __init__ (no disk cache, no prior
    build) -- _derive_gravity_isostatic only touches self.config and
    self.structure_data."""
    runner = object.__new__(MCMCRunner)
    runner.config = config
    runner.structure_data = structure_data
    return runner


def _direct_isostasy_reference(D_shell_km, D_ocean_km, C2=1.0):
    """Reference computed by calling isostasy.py's primitives directly
    (independent of mcmc_runner), assembling the SAME hydrostatic-sum +
    isostatic_gravity decomposition _derive_gravity_isostatic uses."""
    struct = _synthetic_ocean_struct(D_shell_km, D_ocean_km)
    r, rho = struct['r_m'], struct['rho']
    ep, eq = compute_eccentricities(r, rho, OMEGA)
    H_hyd = []
    for i, R_i in enumerate(r):
        a, b, c = eccentricities_to_axes(float(R_i), float(ep[i]), float(eq[i]))
        H20, H22 = triaxial_to_H2m(a, b, c)
        H_hyd.append({(2, 0): H20, (2, 2): H22})
    n = 3
    C_hyd = {lm: 0.0 for lm in H_OBS}
    for i in range(n):
        drho = float(rho[i] - (rho[i + 1] if i + 1 < n else 0.0))
        H_i = finite_amplitude_coeffs(H_hyd[i], float(r[i]))
        for lm in C_hyd:
            C_hyd[lm] += interface_gravity_coeff(
                H_i.get(lm, 0.0), drho, float(r[i]), lm[0], M_ENC, R_REF)
    C_nh = isostatic_gravity(
        H_obs_lm=H_OBS, H_hyd_t_lm=H_hyd[-1], R_t=float(r[-1]),
        R_b=float(r[1]), rho_ice=float(rho[2]), rho_ocean=float(rho[1]),
        M_body=M_ENC, R_ref=R_REF, C2=C2, finite_amplitude=True)
    return {lm: C_hyd[lm] + C_nh.get(lm, 0.0) for lm in H_OBS}


def test_derive_gravity_isostatic_finite_on_3layer_synthetic():
    """C20/C22/C30 dispatch: finite values on a 3-layer synthetic ocean
    structure (isostatic_hm2019 forward model)."""
    config = _isostatic_config()
    struct = _synthetic_ocean_struct(20.0, 39.0)
    runner = _bare_runner(config, struct)
    assert runner._gravity_isostatic_active()
    out = runner._derive_gravity_isostatic({'compensation_C2': 1.0})
    assert out is not None
    for lm in ((2, 0), (2, 2), (3, 0)):
        assert np.isfinite(out[lm]), lm


def test_derive_gravity_isostatic_matches_direct_isostasy_calls():
    """The dispatch must reproduce a hand-assembled call into
    isostasy.py's own primitives (isostatic_gravity + the hydrostatic
    Tricarico interface sum) to numerical precision -- this is what
    verifies the WIRING, independent of any finite-amplitude convention
    question (see the gate-comparison test below for that)."""
    config = _isostatic_config()
    struct = _synthetic_ocean_struct(20.0, 39.0)
    runner = _bare_runner(config, struct)
    out = runner._derive_gravity_isostatic({'compensation_C2': 1.0})
    ref = _direct_isostasy_reference(20.0, 39.0, C2=1.0)
    for lm in ref:
        assert out[lm] == pytest.approx(ref[lm], rel=1e-9), lm


@pytest.mark.parametrize('C2', [0.0, 0.5, 1.0])
def test_derive_gravity_isostatic_c2_sweep_matches_direct(C2):
    config = _isostatic_config()
    struct = _synthetic_ocean_struct(22.0, 36.0)
    runner = _bare_runner(config, struct)
    out = runner._derive_gravity_isostatic({'compensation_C2': C2})
    ref = _direct_isostasy_reference(22.0, 36.0, C2=C2)
    for lm in ref:
        assert out[lm] == pytest.approx(ref[lm], rel=1e-9), (C2, lm)


def _predict_reference(D_shell_km, D_ocean_km):
    """Trimmed reimplementation of
    tests/isostasy_hm2019_gate_test.py::predict() -- kept inline (rather
    than imported) so this file has no cross-test-module coupling."""
    from PlanetProfile.Gravity.isostasy import airy_root, gravity_g
    struct = _synthetic_ocean_struct(D_shell_km, D_ocean_km)
    r, rho = struct['r_m'], struct['rho']
    R_c, R_b, R_t = r
    rho_core, rho_ocean, rho_ice = rho
    ep, eq = compute_eccentricities(r, rho, OMEGA)
    H_hyd = []
    for i, R_i in enumerate(r):
        a, b, c = eccentricities_to_axes(float(R_i), float(ep[i]), float(eq[i]))
        H20, H22 = triaxial_to_H2m(a, b, c)
        H_hyd.append({(2, 0): H20, (2, 2): H22})
    H_nh_t = {lm: H_OBS.get(lm, 0.0)
             - (H_hyd[2].get(lm, 0.0) if lm[0] == 2 else 0.0) for lm in H_OBS}
    m_shell = rho_ice * (4 * np.pi / 3) * (R_t ** 3 - R_b ** 3)
    M_below_shell = M_ENC - m_shell
    g_t = gravity_g(M_ENC, R_t)
    g_b = gravity_g(M_below_shell, R_b)
    H_nh_b = airy_root(H_nh_t, rho_ice, rho_ocean, g_t, g_b, C2=1.0)
    interfaces = (
        (R_t, rho_ice, {lm: H_hyd[2].get(lm, 0.0) + H_nh_t.get(lm, 0.0)
                        for lm in H_OBS}),
        (R_b, rho_ocean - rho_ice,
         {lm: H_hyd[1].get(lm, 0.0) + H_nh_b.get(lm, 0.0) for lm in H_OBS}),
        (R_c, rho_core - rho_ocean, {lm: H_hyd[0].get(lm, 0.0) for lm in H_OBS}),
    )
    C_pred = {lm: 0.0 for lm in H_OBS}
    for R_i, drho, H_lm in interfaces:
        H_adj = finite_amplitude_coeffs(H_lm, float(R_i), n_theta=48, n_phi=96)
        for lm in C_pred:
            C_pred[lm] += interface_gravity_coeff(
                H_adj[lm], drho, float(R_i), lm[0], M_ENC, R_REF)
    return C_pred


def test_derive_gravity_isostatic_reproduces_hm2019_gate_predict():
    """Looser cross-check against the B13 gate's predict() -- the
    config's own required_unit_test
    (gravity_forward_model_contract.required_unit_test). The two
    computations differ by a genuine nonlinear finite-amplitude
    cross-term: predict() applies Wieczorek & Phillips (H&M Eqs. 5-7) to
    the COMBINED (hydrostatic + non-hydrostatic) topography at each
    interface, while this dispatch (module spec Physics map item 4, "C_lm
    = C_lm_hyd + C_lm_nh(isostatic)") applies it to each part SEPARATELY
    and sums -- matching isostatic_gravity()'s own already-tested,
    already-ratified (B11) convention of FA-ing the non-hydrostatic part
    alone. Measured discrepancy at shell=20/ocean=39 km (near the H&M
    solution): C20 +0.038%, C22 -0.121%, C30 +0.948% (rel), all < 0.06
    sigma_obs. Both a relative and a sigma_obs-normalized bound are
    asserted."""
    config = _isostatic_config()
    struct = _synthetic_ocean_struct(20.0, 39.0)
    runner = _bare_runner(config, struct)
    out = runner._derive_gravity_isostatic({'compensation_C2': 1.0})
    C_pred = _predict_reference(20.0, 39.0)
    for lm in C_pred:
        rel = abs(out[lm] - C_pred[lm]) / abs(C_pred[lm])
        sigma_frac = abs(out[lm] - C_pred[lm]) / SIG_C[lm]
        assert rel < 0.02, (lm, rel)
        assert sigma_frac < 0.2, (lm, sigma_frac)


def test_derive_gravity_isostatic_frozen_branch_returns_none():
    """Ocean-branch-only scope (module spec: 'implement the ocean branch
    first'): a structure with no phase-0 layer must be rejected, not
    silently treated as rigid support."""
    config = _isostatic_config()
    struct = {
        'r_m': np.array([100e3, 252.1e3]),
        'rho': np.array([2400.0, 925.0]),
        'phases': np.array([50, 1]),
        'omega': OMEGA,
    }
    runner = _bare_runner(config, struct)
    assert runner._derive_gravity_isostatic({'compensation_C2': 1.0}) is None


def test_derive_gravity_isostatic_rho_ice_override_is_mass_neutral():
    """rho_ice_kgm3 nuisance (registered above) must flow through the
    same mass_neutral_shell_density path _derive_libration_deg uses --
    the shell mass delta is absorbed by the interior, not left to drift."""
    config = _isostatic_config()
    struct = _synthetic_ocean_struct(25.0, 36.0)
    runner = _bare_runner(config, struct)
    baseline = runner._derive_gravity_isostatic(
        {'compensation_C2': 1.0, 'rho_ice_kgm3': float(struct['rho'][2])})
    identity_override = runner._derive_gravity_isostatic(
        {'compensation_C2': 1.0})
    for lm in baseline:
        assert baseline[lm] == pytest.approx(identity_override[lm], rel=1e-9)
    perturbed = runner._derive_gravity_isostatic(
        {'compensation_C2': 1.0, 'rho_ice_kgm3': 935.0})
    assert perturbed is not None
    assert perturbed[(2, 0)] != pytest.approx(baseline[(2, 0)])


# ===========================================================================
# Task 2 (cont'd) -- full log_likelihood integration on the real Enceladus
# smoke cache (existing, committed; not modified by this task).
# ===========================================================================

SMOKE_CACHE = (REPO / 'PlanetProfile' / 'Test' / 'mcmc_results' / 'Enceladus'
              / 'Cassini_smoke' / 'enceladus_seawater_smoke_grid.pkl')


def _smoke_isostatic_config():
    return InferenceConfig(
        mode='mcmc', bodyname='Enceladus',
        param_space={
            'Tb_K': {'prior_type': 'uniform', 'bounds': [271.8, 272.46]},
            'compensation_C2': {'prior_type': 'uniform', 'bounds': [0.0, 1.0]},
        },
        observables={'C20': (-5521e-6, 3.58e-4), 'C22': (1574e-6, 1.59e-4),
                     'C30': (118e-6, 2.38e-4)},
        sampler_settings={},
        structure_cache_path=str(SMOKE_CACHE),
        gravity_forward_model='isostatic_hm2019',
        planet_template_module='PlanetProfile.Default.Enceladus.PPEnceladus',
        metadata={'isostasy': {
            'H_obs_lm_m': {'2,0': -3510.0, '2,2': 857.0, '3,0': 420.0},
            'shape_ref_radius_m': 252220.0,
            'gravity_ref_radius_m': 256600.0,
            'finite_amplitude': True,
        }},
    )


@pytest.mark.skipif(not SMOKE_CACHE.exists(),
                    reason='Enceladus smoke structure cache not present')
def test_log_likelihood_isostatic_dispatch_on_real_smoke_cache():
    """Full log_likelihood dispatch (C20/C22/C30 via isostatic_hm2019)
    against the existing, committed 3-node Enceladus smoke cache: finite
    at the one ocean node (Tb=272.46, node 2 per RESUME_NOTE.md), hard
    -1e30 at a frozen node (Tb=271.8, node 0) -- ocean-branch-only
    scope."""
    runner = MCMCRunner(_smoke_isostatic_config())
    ll_ocean = runner.log_likelihood_fn(np.array([272.46, 1.0]))
    assert np.isfinite(ll_ocean)
    assert ll_ocean > -1e29

    ll_frozen = runner.log_likelihood_fn(np.array([271.8, 1.0]))
    assert ll_frozen == -1e30


# ===========================================================================
# Task 4 -- B4 zb-axis cache-builder mode (real PlanetProfile builds).
# ===========================================================================

@pytest.mark.slow
def test_build_zbw_grid_cache_smoke_3x2(tmp_path):
    """Tiny (3 zb x 2 w = 6 node) real-PlanetProfile smoke build of the
    B4 zb-axis cache mode. Minutes, not hours. Verifies: the two-of-three
    zb/Tb pairing actually drives the build (solved Tb_K differs node to
    node and is recorded, not echoed from a placeholder), the zb-
    placement invariant rejects out-of-tolerance nodes to None rather
    than silently mis-recording them, and every surviving node has an
    ocean layer (ocean-branch-only scope)."""
    from PlanetProfile.Inference.cache_builder import build_zbw_grid_cache

    out_path = tmp_path / 'enceladus_zbw_smoke.pkl'
    cache = build_zbw_grid_cache(
        'PlanetProfile.Default.Enceladus.PPEnceladus',
        zb_km_grid=[20.0, 25.0, 30.0],
        wOcean_ppt_grid=[5.0, 20.0],
        output_path=str(out_path),
        ocean_overrides={'comp': 'Seawater', 'deltaT': 0.002},
        # Widen the MoI acceptance window -- same accommodation the
        # frozen-config candidate's own bulk_overrides note documents is
        # needed for non-default Enceladus builds (config metadata.
        # template_radius_and_MoI_window).
        bulk_overrides={'Cuncertainty': 0.08},
        zb_tol_km=5.0,
        progress=False,
    )
    assert out_path.exists()
    assert cache['schema_version'] == 'v3.1-zbw'
    assert cache['frozen_branch_supported'] is False
    assert len(cache['structures']) == 6

    n_ok = sum(s is not None for s in cache['structures'])
    assert n_ok >= 1, 'every node was rejected -- zb-driven build not working'

    solved_Tb = set()
    for i, s in enumerate(cache['structures']):
        if s is None:
            continue
        # Every surviving node is on the ocean branch.
        assert bool(np.any(np.asarray(s['phases']) == 0)) is True
        assert s['has_ocean'] is True
        # zb-placement invariant actually holds for anything that survived.
        assert abs(s['zb_km_actual'] - s['zb_km_node']) < cache['zb_tol_km']
        # The solved Tb_K must NOT be the unused placeholder (272.0 K) --
        # confirms the two-of-three zb/Tb pairing drove the build, not an
        # echoed input.
        assert s['Tb_K'] != pytest.approx(272.0, abs=1e-6)
        solved_Tb.add(round(float(s['Tb_K']), 4))

    # Different zb targets should not all solve to the identical Tb_K.
    assert len(solved_Tb) >= 1
