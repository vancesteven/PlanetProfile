"""Enceladus 3-layer stack mass closure (2026-10-21 reviewer ruling).

validation_reports/enceladus_isostasy/c3_root_cause_2026-10-21/README.md:
the reduced 3-layer stack loses ~0.04-0.46% of M_target (candidate pick +
lossy ALMA reduction). ``metadata.isostasy.stack_mass_closure ==
'rock_density'`` rescales the interior (rock) density so the stack mass is
exactly ``struct['Mtot_kg']`` in BOTH ``_derive_gravity_isostatic`` and
``_derive_libration_deg``. Absent flag == historical behaviour, bit-identical.
"""
import inspect
import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402
from PlanetProfile.Inference.mcmc_runner import MCMCRunner  # noqa: E402

CFG_DIR = REPO / 'PlanetProfile' / 'Inference' / 'configs'
CANDIDATE = CFG_DIR / 'enceladus_cassini_isostasy_7D.json'

R_T = 252.22e3
OMEGA = 5.307e-5
ECC = 0.0047
V = 4.0 * np.pi / 3.0
H_OBS = {'2,0': -3510.0, '2,2': 857.0, '3,0': 420.0}


def _struct(deficit_frac=0.0, rho_rock=None):
    """Cache-style structure: rock | ocean | shell, many cells per layer.
    Mtot_kg is the TARGET; the cell densities carry `deficit_frac` mass
    deficit in the rock (rock density low by construction)."""
    R_b, R_c = R_T - 20e3, R_T - 20e3 - 40e3
    R_sil = R_c - 0.0
    r_rock = np.linspace(R_sil / 20, R_sil, 20)
    r_oc = np.linspace(R_sil, R_b, 10)[1:]
    r_sh = np.linspace(R_b, R_T, 6)[1:]
    r = np.concatenate([r_rock, r_oc, r_sh])
    rho_o, rho_i = 1005.0, 925.0
    m_o = rho_o * V * (R_b ** 3 - R_sil ** 3)
    m_i = rho_i * V * (R_T ** 3 - R_b ** 3)
    M = 1.08e20
    rr = (M - m_o - m_i) / (V * R_sil ** 3)
    if rho_rock is None:
        rho_rock = rr * (1.0 - deficit_frac)
    rho = np.concatenate([np.full(r_rock.size, rho_rock),
                          np.full(r_oc.size, rho_o),
                          np.full(r_sh.size, rho_i)])
    ph = np.concatenate([np.full(r_rock.size, 3), np.zeros(r_oc.size, int),
                         np.ones(r_sh.size, int)])
    return {'r_m': r, 'rho': rho, 'phases': ph, 'omega': OMEGA,
            'eccentricity': ECC, 'Mtot_kg': M}, rr


def _runner(struct, closure=None):
    iso = {'H_obs_lm_m': dict(H_OBS), 'shape_ref_radius_m': R_T,
           'gravity_ref_radius_m': R_T, 'finite_amplitude': True}
    if closure is not None:
        iso['stack_mass_closure'] = closure
    cfg = InferenceConfig(
        mode='mcmc', bodyname='Enceladus',
        param_space={'compensation_C2': {'prior_type': 'uniform',
                                         'bounds': [0.0, 1.0]}},
        observables={'libration_deg': (0.092, 0.003)}, sampler_settings={},
        gravity_forward_model='isostatic_hm2019',
        metadata={'isostasy': iso})
    r = object.__new__(MCMCRunner)
    r.config = cfg
    r.structure_data = struct
    return r


THETA = {'compensation_C2': 0.5}


def _grav(struct, closure=None):
    return _runner(struct, closure)._derive_gravity_isostatic(THETA)


def test_closure_makes_stack_mass_equal_target():
    s, _ = _struct(deficit_frac=0.0046)
    r = _runner(s, 'rock_density')
    rr = np.array([100e3, 192.22e3, R_T])
    rin = np.concatenate(([0.0], rr[:-1]))
    rho = np.array([2300.0, 1005.0, 925.0])
    s = dict(s, Mtot_kg=float(np.sum(rho * V * (rr ** 3 - rin ** 3))) * 1.0046)
    out = r._close_reduced_stack_mass(rr, rho, s)
    M = np.sum(out * V * (rr ** 3 - rin ** 3))
    assert abs(M / s['Mtot_kg'] - 1.0) < 1e-12
    assert out[1] == rho[1] and out[2] == rho[2]      # only rock moved
    d = r.last_mass_closure
    assert d['rock_density_factor'] == pytest.approx(out[0] / rho[0])
    M0 = np.sum(rho * V * (rr ** 3 - rin ** 3))
    assert d['pre_closure_residual_frac'] == pytest.approx(
        (M0 - s['Mtot_kg']) / s['Mtot_kg'])


def test_closed_deficit_node_equals_consistent_node():
    """A node with a 0.46% rock-mass deficit, closure ON, must give the
    same C20/C22/C30 and libration as the mass-consistent node (closure
    recovers exactly the density the deficit removed from the rock)."""
    s_def, rr = _struct(deficit_frac=0.0046)
    s_ok, _ = _struct(deficit_frac=0.0)
    g_on = _grav(s_def, 'rock_density')
    g_ok = _grav(s_ok)
    for lm in g_on:
        assert g_on[lm] == pytest.approx(g_ok[lm], rel=1e-9)
    l_on = _runner(s_def, 'rock_density')._derive_libration_deg(THETA)
    l_ok = _runner(s_ok)._derive_libration_deg(THETA)
    assert l_on == pytest.approx(l_ok, rel=1e-9)
    r = _runner(s_def, 'rock_density')
    r._derive_gravity_isostatic(THETA)
    assert r.last_mass_closure['rock_density_factor'] > 1.0
    assert r.last_mass_closure['pre_closure_residual_frac'] < 0.0


def test_deficit_node_c20_moves_documented_direction():
    """README: closing a mass deficit moves C20 up (+) and C22 down (-)."""
    s_def, _ = _struct(deficit_frac=0.0046)
    off, on = _grav(s_def), _grav(s_def, 'rock_density')
    assert on[(2, 0)] > off[(2, 0)]
    assert on[(2, 2)] < off[(2, 2)]
    # magnitude sanity: a 0.46% deficit moves C20 by a sub-percent fraction
    assert abs(on[(2, 0)] / off[(2, 0)] - 1.0) < 0.05


@pytest.mark.parametrize('flag', [None, 'none', 'off', ''])
def test_flag_off_is_bit_identical_and_ignores_target_mass(flag):
    s, _ = _struct(deficit_frac=0.0046)
    base = _grav(s)
    assert _grav(s, flag) == base                      # exact equality
    s2 = dict(s)
    s2['Mtot_kg'] = s['Mtot_kg'] * 1.01                # target must not matter
    assert _grav(s2, flag) == base
    lib = _runner(s)._derive_libration_deg(THETA)
    assert _runner(s2, flag)._derive_libration_deg(THETA) == lib


def test_closure_on_rejects_missing_target():
    s, _ = _struct(0.001)
    s.pop('Mtot_kg')
    assert _grav(s, 'rock_density') is None
    assert _grav(s) is not None


def test_rho_ice_override_still_conserves_closed_mass():
    s, _ = _struct(0.0046)
    th = dict(THETA, rho_ice_kgm3=920.0)
    r = _runner(s, 'rock_density')
    g = r._derive_gravity_isostatic(th)
    assert g is not None and all(np.isfinite(v) for v in g.values())
    # identity override == no override (mass-neutral about the closed stack)
    g0 = _grav(s, 'rock_density')
    g1 = _runner(s, 'rock_density')._derive_gravity_isostatic(
        dict(THETA, rho_ice_kgm3=925.0))
    for lm in g0:
        assert g1[lm] == pytest.approx(g0[lm], rel=1e-9)


def test_enceladus_config_enables_flag_others_do_not():
    cfg = InferenceConfig.from_json(str(CANDIDATE))
    r = object.__new__(MCMCRunner)
    r.config = cfg
    assert r._stack_mass_closure_mode() == 'rock_density'
    # Every other campaign config (Titan / Europa clairaut, Callisto, ...)
    # must not carry the flag.
    for p in sorted(CFG_DIR.glob('*.json')):
        if p == CANDIDATE:
            continue
        txt = p.read_text()
        assert 'stack_mass_closure' not in txt, p.name


def test_sigma_model_add_carries_c20_c22():
    """The existing B15 wiring takes {C20, C22} model sigma; the jitter
    gate's measured sigma is added in quadrature (hypot) by _inflate_sigma."""
    meta = json.loads(CANDIDATE.read_text())['metadata']
    assert {'C20', 'C22'} <= set(meta['sigma_model_add'])
    src = inspect.getsource(MCMCRunner)
    assert 'sigma_model_add' in src


def test_clairaut_paths_do_not_use_closure():
    """Titan/Europa (clairaut_hydrostatic) close mass via
    derive_silicate_density; they must never touch the stack closure."""
    src = inspect.getsource(MCMCRunner)
    users = set()
    for m in re.finditer(r'\n    def (\w+)\(', src):
        users.add((m.start(), m.group(1)))
    order = sorted(users)
    calls = [m.start() for m in re.finditer(
        r'self\._close_reduced_stack_mass\(', src)]
    callers = set()
    for c in calls:
        callers.add([n for p, n in order if p < c][-1])
    assert callers == {'_derive_gravity_isostatic', '_derive_libration_deg'}
    for name in ('_derive_gravity_pair', '_derive_cmr2_via_mass_conservation'):
        body = inspect.getsource(getattr(MCMCRunner, name))
        assert '_close_reduced_stack_mass' not in body
        assert '_stack_mass_closure_mode' not in body


def test_closure_rejects_broken_stack_beyond_1pct():
    """Reviewer 2026-10-21: |pre-closure residual| > 1e-2 is a hard reject."""
    rr = np.array([100e3, 192.22e3, R_T])
    rho = np.array([2300.0, 1005.0, 925.0])
    rin = np.concatenate(([0.0], rr[:-1]))
    m = rho * V * (rr ** 3 - rin ** 3)
    for frac, ok in ((-0.02, False), (+0.02, False), (-0.005, True)):
        # target such that the stack residual (M_stack - M)/M = frac
        M = float(m.sum()) / (1.0 + frac)
        s = {'Mtot_kg': M}
        r = _runner(s, 'rock_density')
        out = r._close_reduced_stack_mass(rr, rho, s)
        assert (out is not None) == ok
        if ok:
            assert abs(np.sum(out * V * (rr ** 3 - rin ** 3)) / M - 1) < 1e-12
