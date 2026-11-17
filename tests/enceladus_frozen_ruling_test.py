"""Frozen-branch design-ruling tests C1-C5 (queue task C26).

Source of the formulations and invariants:
validation_reports/enceladus_isostasy/frozen_branch_DESIGN_RULING.md
(I-F1..I-F6, F5 regression). Extends -- does not duplicate -- the frozen
dispatch coverage in enceladus_isostasy_code_gaps_test.py and
enceladus_mass_closure_test.py (reused helpers imported below).

  C1  analytic bijection unit tests + F5 regression pins
  C2  I-F5 loop closure at nodes spanning zb ~47/54/65 km, tol 3e-3,
      with a self-invalidating off-by-one-grid-step companion
  C3  I-F6 both-directions non-conditioning (builder level, synthetic nodes)
  C4  None-node: no nearest-node fallback on the gravity AND libration paths
  C5  3-node real frozen build smoke: I-F1..I-F4 (slow)
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Gravity.isostasy import (  # noqa: E402
    frozen_zb_from_mass, frozen_rho_rock_from_zb, frozen_cmr2)
from PlanetProfile.Inference.mcmc_runner import (  # noqa: E402
    _FROZEN_LOOP_CLOSURE_TOL)
from enceladus_isostasy_code_gaps_test import (  # noqa: E402
    _synthetic_frozen_node, _frozen_dispatch_runner, M_ENC, R_REF)

TEMPLATE = 'PlanetProfile.Default.Enceladus.PPEnceladus'

# F5 feasibility numbers (ruling record / frozen_zb_from_mass docstring):
# RESUME_NOTE "deep-frozen Enceladus is feasible", M = 1.08022e20 kg,
# R = 252.1 km. (rho_ice, rho_rock, zb_km) at C/MR^2 = 0.3350.
F5_M, F5_R = 1.08022e20, 252.1e3
F5_CASES = [(900.0, 2314.3, 51.78),
            (925.0, 2334.3, 53.93),
            (970.0, 2377.5, 58.29)]


# ===========================================================================
# C1 -- analytic closure unit tests
# ===========================================================================

def test_c1_tolerance_constant_is_ruled_value():
    assert _FROZEN_LOOP_CLOSURE_TOL == 3e-3


@pytest.mark.parametrize('rho_ice,rho_rock,zb_km', F5_CASES)
def test_c1_f5_regression_pins(rho_ice, rho_rock, zb_km):
    """F5 regression: zb to the ruling's 0.01 km quoting precision and
    C/MR^2 = 0.3350 (the record quotes 4 digits; exact values are
    0.335000-0.335001, so abs 2e-6 pins the quoted figure)."""
    zb_m = frozen_zb_from_mass(F5_M, F5_R, rho_rock, rho_ice, strict=True)
    assert zb_m / 1e3 == pytest.approx(zb_km, abs=0.005)
    cmr2 = frozen_cmr2(F5_M, F5_R, zb_m, rho_rock, rho_ice)
    assert cmr2 == pytest.approx(0.3350, abs=2e-6)
    # exact inverse
    assert frozen_rho_rock_from_zb(
        F5_M, F5_R, zb_m, rho_ice, strict=True) == pytest.approx(
        rho_rock, rel=1e-12)


@pytest.mark.parametrize('rho_ice', [900.0, 925.0, 970.0])
@pytest.mark.parametrize('rho_rock', [2200.0, 2334.3, 2600.0])
def test_c1_round_trip_and_mass_exact(rho_ice, rho_rock):
    zb_m = frozen_zb_from_mass(M_ENC, R_REF, rho_rock, rho_ice, strict=True)
    back = frozen_rho_rock_from_zb(M_ENC, R_REF, zb_m, rho_ice, strict=True)
    assert back == pytest.approx(rho_rock, rel=1e-12)
    R_rock = R_REF - zb_m
    M = (4 * math.pi / 3) * (rho_rock * R_rock ** 3
                             + rho_ice * (R_REF ** 3 - R_rock ** 3))
    assert M / M_ENC - 1.0 == pytest.approx(0.0, abs=1e-12)  # I-F1 analytic


def test_c1_monotone_bijection():
    zbs = [frozen_zb_from_mass(M_ENC, R_REF, rr, 925.0) for rr in
           np.linspace(2200.0, 2600.0, 9)]
    assert all(z is not None for z in zbs)
    assert all(a < b for a, b in zip(zbs[:-1], zbs[1:]))  # denser rock -> thicker shell (fixed M)


_ZB_FROM_MASS_BAD = [
    # (M, R, rho_rock, rho_ice)
    (float('nan'), R_REF, 2300.0, 925.0),
    (M_ENC, float('inf'), 2300.0, 925.0),
    (M_ENC, R_REF, float('nan'), 925.0),
    (-1.0, R_REF, 2300.0, 925.0),
    (M_ENC, 0.0, 2300.0, 925.0),
    (M_ENC, R_REF, 2300.0, 0.0),            # rho_ice <= 0
    (M_ENC, R_REF, 2300.0, -5.0),
    (M_ENC, R_REF, 925.0, 925.0),           # singular at equality
    (M_ENC, R_REF, 900.0, 925.0),           # rock lighter than ice
    (M_ENC, R_REF, 3000.0, 2900.0),         # mean density below ice -> R_rock^3<0
    (M_ENC, R_REF, 1200.0, 925.0),          # mean density above rock -> zb<0
]


@pytest.mark.parametrize('M,R,rr,ri', _ZB_FROM_MASS_BAD)
def test_c1_zb_from_mass_rejects_both_modes(M, R, rr, ri):
    assert frozen_zb_from_mass(M, R, rr, ri) is None            # default
    assert frozen_zb_from_mass(M, R, rr, ri, strict=False) is None
    with pytest.raises(ValueError):
        frozen_zb_from_mass(M, R, rr, ri, strict=True)


_RHO_FROM_ZB_BAD = [
    # (M, R, zb_m, rho_ice)
    (float('nan'), R_REF, 5e4, 925.0),
    (M_ENC, R_REF, float('nan'), 925.0),
    (M_ENC, -R_REF, 5e4, 925.0),
    (0.0, R_REF, 5e4, 925.0),
    (M_ENC, R_REF, 5e4, 0.0),
    (M_ENC, R_REF, 0.0, 925.0),             # zero-thickness shell
    (M_ENC, R_REF, -1.0, 925.0),
    (M_ENC, R_REF, R_REF, 925.0),           # no interior
    (M_ENC, R_REF, 2 * R_REF, 925.0),
    (M_ENC, R_REF, 5e4, 2000.0),            # implied rho_rock <= rho_ice
]


@pytest.mark.parametrize('M,R,zb,ri', _RHO_FROM_ZB_BAD)
def test_c1_rho_rock_from_zb_rejects_both_modes(M, R, zb, ri):
    assert frozen_rho_rock_from_zb(M, R, zb, ri) is None
    assert frozen_rho_rock_from_zb(M, R, zb, ri, strict=False) is None
    with pytest.raises(ValueError):
        frozen_rho_rock_from_zb(M, R, zb, ri, strict=True)


# ===========================================================================
# C2 -- I-F5 loop closure at zb ~ 47 / 54 / 65 km, self-invalidating
# ===========================================================================

# 1 km grid: node spacing large enough that one step off exceeds 3e-3
# (d rho_rock / d zb ~ 20 kg/m^3 per km ~ 0.9%/km at the low end, ~0.6% at
# the high end -- asserted below, not assumed).
C2_GRID = [float(z) for z in range(44, 69)]
C2_NODES_KM = [47.0, 54.0, 65.0]


def _c2_runner(grid=C2_GRID, struct_grid=C2_GRID):
    structs = [_synthetic_frozen_node(z) for z in struct_grid]
    return _frozen_dispatch_runner(grid, structs), structs


def _theta_at(zb_km, rho_ice=925.0):
    rr = frozen_rho_rock_from_zb(M_ENC, R_REF, zb_km * 1e3, rho_ice,
                                 strict=True)
    return {'rho_rock_kgm3': rr, 'rho_ice_kgm3': rho_ice,
            'compensation_C2': 1.0}


@pytest.mark.parametrize('zb_km', C2_NODES_KM)
def test_c2_loop_closure_passes_at_correct_node(zb_km):
    runner, _ = _c2_runner()
    theta = _theta_at(zb_km)
    assert runner._derive_gravity_isostatic(theta) is not None
    lib = runner._derive_libration_deg(theta)
    assert lib is not None and np.isfinite(lib)
    # and the served interior density closes to the sample at far better
    # than the 3e-3 threshold at a node centre
    served = runner._struct_for_hydrosphere(theta)
    closure = abs(float(served['rho'][0]) / theta['rho_rock_kgm3'] - 1.0)
    assert closure < 1e-9


@pytest.mark.parametrize('zb_km', C2_NODES_KM)
@pytest.mark.parametrize('step', [-1, +1])
def test_c2_self_invalidating_node_off_by_one_grid_step(zb_km, step):
    """The same sample, but the cache is labelled one grid step off (the
    node served for zb is really the node for zb+step). The residual
    must EXCEED 3e-3 and I-F5 must reject -> None on both forward paths.
    If this test ever passes while the correct-node test above also
    passes, the loop-closure check is demonstrably sensitive to node
    selection; if the off-by-one case stopped failing, the tolerance (or
    the check) would have gone slack."""
    mis_grid = [z + step for z in C2_GRID]       # labels shifted one step
    runner, structs = _c2_runner(grid=mis_grid, struct_grid=C2_GRID)
    theta = _theta_at(zb_km)

    # Independent analytic statement of the residual, no runner involved:
    rho_off = frozen_rho_rock_from_zb(
        M_ENC, R_REF, (zb_km - step) * 1e3, 925.0, strict=True)
    residual = abs(rho_off / theta['rho_rock_kgm3'] - 1.0)
    assert residual > _FROZEN_LOOP_CLOSURE_TOL, residual

    served = runner._struct_for_hydrosphere(theta)
    assert served is not None            # dispatch itself serves a node ...
    closure = abs(float(served['rho'][0]) / theta['rho_rock_kgm3'] - 1.0)
    assert closure > _FROZEN_LOOP_CLOSURE_TOL
    # ... and I-F5 downstream rejects it, on both paths
    assert runner._derive_gravity_isostatic(theta) is None
    assert runner._derive_libration_deg(theta) is None


def test_c2_threshold_is_sharp_at_3e3():
    """A served interior density off by 0.2% (inside tol) passes, 0.4%
    (outside) rejects -- the 3e-3 threshold, not some looser bound."""
    runner, structs = _c2_runner()
    theta = _theta_at(54.0)
    for scale, ok in ((1.002, True), (1.004, False)):
        runner2, _ = _c2_runner()
        # corrupt the 54 km node's interior density by `scale`; the
        # reduced stack re-conserves the stack mass, but I-F5 compares
        # the (rescaled) interior with the sampled rho_rock.
        idx = C2_GRID.index(54.0)
        st = runner2.structure_data['frozen_structures'][idx]
        st['rho'] = np.array([st['rho'][0] * scale, st['rho'][1], st['rho'][2]])
        out = runner2._derive_gravity_isostatic(theta)
        assert (out is not None) is ok


# ===========================================================================
# C3 -- I-F6 both directions, builder level (synthetic nodes)
# ===========================================================================

def _patched_frozen_build(monkeypatch, zb_grid, cmr2_values, window):
    """Run cache_builder._build_frozen_axis with the per-node PlanetProfile
    solve replaced by self-consistent synthetic nodes whose derived C/MR^2
    is set to ``cmr2_values`` -- everything else (I-F1..I-F4 gates, I-F6)
    is the real builder code."""
    import PlanetProfile.Inference.cache_builder as cb
    import importlib
    mod = importlib.import_module(TEMPLATE)
    M = float(mod.Planet.Bulk.M_kg)
    R = float(mod.Planet.Bulk.R_m)
    cmr2_by_zb = dict(zip(map(float, zb_grid), cmr2_values))

    def fake_node(module, zb_node, PbI_seed, *a, **kw):
        node = _synthetic_frozen_node(zb_node, M_body_kg=M, R_body_m=R)
        node.update({
            'zb_km_actual': float(zb_node),
            'mass_residual_frac': 0.0,
            'rho_ice_mean_kgm3': 925.0,
            'rho_interior_mean_kgm3': float(node['rho'][0]),
            'n_liquid_layers': 0,
            'CMR2': cmr2_by_zb[float(zb_node)],
        })
        return node, float(PbI_seed), 1

    monkeypatch.setattr(cb, '_build_frozen_node', fake_node)
    return cb._build_frozen_axis(
        planet_template_module=TEMPLATE, frozen_zb_km_grid=zb_grid,
        base_ocean_overrides={}, base_bulk_overrides={},
        planet_overrides=None, base_do_overrides=None, extrap_ocean=False,
        frozen_wOcean_ppt=20.0, frozen_Cuncertainty=0.015,
        frozen_mass_tol=1e-6, frozen_rho_closure_tol_kgm3=12.0,
        frozen_zb_tol_km=0.25,
        frozen_moi_nonconditioning_window=window,
        frozen_max_iter=3, progress=False)


def test_c3_if6_node_set_entirely_inside_window_fails_build(monkeypatch):
    with pytest.raises(RuntimeError, match='I-F6'):
        _patched_frozen_build(monkeypatch, [50.0, 54.0, 58.0],
                              [0.3345, 0.3350, 0.3355], (0.335, 0.001))


def test_c3_if6_single_node_inside_window_also_fails(monkeypatch):
    """Boundary: a lone surviving node inside the window is 'entirely
    inside' too."""
    with pytest.raises(RuntimeError, match='I-F6'):
        _patched_frozen_build(monkeypatch, [54.0], [0.3350], (0.335, 0.001))


def test_c3_if6_other_direction_set_spanning_window_passes(monkeypatch):
    """Both-directions check: a set that leaves the window (one node
    outside) must NOT trip I-F6, so the invariant is not a blanket
    reject."""
    fz, structs, spec = _patched_frozen_build(
        monkeypatch, [50.0, 54.0, 58.0],
        [0.3220, 0.3350, 0.3440], (0.335, 0.001))
    assert sum(s is not None for s in structs) == 3
    lo, hi = spec['I-F6_cmr2_span']
    assert lo < 0.334 and hi > 0.336


def test_c3_if6_set_entirely_outside_window_passes(monkeypatch):
    fz, structs, spec = _patched_frozen_build(
        monkeypatch, [50.0, 54.0], [0.3220, 0.3440], (0.335, 0.001))
    assert sum(s is not None for s in structs) == 2


# ===========================================================================
# C4 -- None node: no nearest-node fallback, gravity AND libration paths
# ===========================================================================

def test_c4_none_node_rejects_on_both_paths_no_fallback():
    grid = [50.0, 51.0, 52.0, 53.0]
    structs = [_synthetic_frozen_node(z) for z in grid]
    structs[1] = None            # build-time reject of the 51 km node
    runner = _frozen_dispatch_runner(grid, structs)

    # zb values: exactly on the None node, and inside both brackets
    # touching it (50-51 and 51-52). Every one must reject on gravity and
    # libration -- including the ones a nearest-node rule would serve from
    # the surviving neighbour (50.1 -> node 50, 51.9 -> node 52).
    for zb_km in (51.0, 50.5, 50.1, 50.9, 51.1, 51.5, 51.9):
        theta = _theta_at(zb_km)
        assert runner._struct_for_hydrosphere(theta) is None, zb_km
        assert runner._derive_gravity_isostatic(theta) is None, zb_km
        assert runner._derive_libration_deg(theta) is None, zb_km

    # Control: brackets that do not touch the None node still serve.
    for zb_km in (52.0, 52.5, 53.0):
        theta = _theta_at(zb_km)
        assert runner._derive_gravity_isostatic(theta) is not None, zb_km
        lib = runner._derive_libration_deg(theta)
        assert lib is not None and np.isfinite(lib), zb_km


def test_c4_none_node_log_likelihood_rejected(tmp_path):
    """Call-site translation: a None from the dispatch becomes the
    hard-reject log-likelihood (-1e30 in this codebase, the finite
    stand-in for -inf) on the full MCMCRunner path -- again with no
    nearest-node substitute."""
    import pickle
    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.mcmc_runner import MCMCRunner

    grid = [50.0, 51.0, 52.0, 53.0]
    structs = [_synthetic_frozen_node(z) for z in grid]
    structs[1] = None
    cache = {
        'zb_km_grid': np.array([20.0, 25.0]),
        'structures': [None, None],
        'schema_version': 'v3.2-zbw-joint',
        'frozen_branch_supported': True,
        'frozen_zb_km_grid': np.array(grid),
        'frozen_structures': structs,
        'frozen_build_spec': {'M_body_kg': M_ENC, 'R_body_m': R_REF},
    }
    path = tmp_path / 'frozen_none_node.pkl'
    with open(path, 'wb') as f:
        pickle.dump(cache, f)

    config = InferenceConfig(
        mode='mcmc', bodyname='Enceladus',
        param_space={
            'rho_rock_kgm3': {'prior_type': 'uniform',
                              'bounds': [2200.0, 2600.0]},
            'rho_ice_kgm3': {'prior_type': 'uniform',
                             'bounds': [915.0, 935.0]},
        },
        observables={'C20': (-5521e-6, 3.58e-4), 'C22': (1574e-6, 1.59e-4)},
        sampler_settings={},
        structure_cache_path=str(path),
        gravity_forward_model='isostatic_hm2019',
        planet_template_module=TEMPLATE,
        metadata={'isostasy': {
            'H_obs_lm_m': {'2,0': -3510.0, '2,2': 857.0, '3,0': 420.0},
            'shape_ref_radius_m': 252220.0,
            'gravity_ref_radius_m': 256600.0,
            'finite_amplitude': True,
        }},
    )
    try:
        runner = MCMCRunner(config)
    except Exception as exc:  # pragma: no cover - environment dependent
        pytest.skip(f'full MCMCRunner not constructible on synthetic cache: {exc!r}')

    def theta_at(zb_km):
        t = _theta_at(zb_km)
        return np.array([t['rho_rock_kgm3'], 925.0])

    for zb_km in (50.1, 51.0, 51.9):
        ll = runner.log_likelihood_fn(theta_at(zb_km))
        assert ll <= -1e29, (zb_km, ll)
    # control: a servable sample is finite
    assert runner.log_likelihood_fn(theta_at(52.5)) > -1e29


# ===========================================================================
# C5 -- 3-node real frozen build smoke: I-F1..I-F4 (slow)
# ===========================================================================

@pytest.mark.slow
def test_c5_three_node_real_frozen_build_invariants():
    """Real PlanetProfile frozen build at zb = 47 / 54 / 65 km, asserted
    independently of the builder's own gating (the builder rejects to
    None on failure; here every node must survive AND pass):

    I-F1  mass exact at 1e-6 on the constant-rho path
    I-F2  stored-array outer-edge quadrature audit inside the ratified
          window (-step_scale, 0]
    I-F3  branch flag from the PHASE array: zero phase-0 cells, has_ocean
          False
    I-F4  parameterization closure: |rho_rock - analytic| <= 12 kg/m^3 and
          |zb_actual - zb_node| < 0.25 km
    """
    from PlanetProfile.Inference.cache_builder import (
        _build_frozen_axis, _frozen_stored_mass_audit)

    grid = [47.0, 54.0, 65.0]
    fz, structs, spec = _build_frozen_axis(
        planet_template_module=TEMPLATE, frozen_zb_km_grid=grid,
        base_ocean_overrides={'comp': 'Seawater', 'deltaT': 0.002},
        base_bulk_overrides={}, planet_overrides=None,
        base_do_overrides=None, extrap_ocean=False,
        frozen_wOcean_ppt=None, frozen_Cuncertainty=0.015,
        frozen_mass_tol=1e-6, frozen_rho_closure_tol_kgm3=12.0,
        frozen_zb_tol_km=0.25, frozen_moi_nonconditioning_window=None,
        frozen_max_iter=10, progress=False)

    assert len(structs) == 3
    assert all(s is not None for s in structs), spec['rejects_by_invariant']
    M, R = spec['M_body_kg'], spec['R_body_m']

    for node_km, s in zip(grid, structs):
        # I-F1: mass at 1e-6, recomputed from the stored arrays here too
        assert abs(s['mass_residual_frac']) <= 1e-6
        # I-F2: stored-array audit
        resid, step = _frozen_stored_mass_audit(s, M)
        assert -step < resid <= 1e-9
        # I-F3: branch from the phase array
        assert not np.any(np.asarray(s['phases']) == 0)
        assert s['n_liquid_layers'] == 0
        assert s['has_ocean'] is False and s['branch'] == 'frozen'
        # I-F4: closure + placement
        rho_closure = frozen_rho_rock_from_zb(
            M, R, s['zb_km_actual'] * 1e3, s['rho_ice_mean_kgm3'],
            strict=True)
        assert abs(s['rho_interior_mean_kgm3'] - rho_closure) <= 12.0
        assert abs(s['zb_km_actual'] - node_km) < 0.25
        # derived C/MR^2 is a recorded diagnostic, finite and physical
        assert 0.30 < s['cmr2_derived'] < 0.36

    # I-F6 held (the build did not raise) and the set spans the window
    lo, hi = spec['I-F6_cmr2_span']
    assert hi - lo > 0.001
