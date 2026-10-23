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


def _runner(struct, closure=None, **iso_extra):
    iso = {'H_obs_lm_m': dict(H_OBS), 'shape_ref_radius_m': R_T,
           'gravity_ref_radius_m': R_T, 'finite_amplitude': True}
    if closure is not None:
        iso['stack_mass_closure'] = closure
    iso.update(iso_extra)
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
    assert r._stack_mass_closure_mode() == 'mass_moi'
    assert r._stack_moi_target() == pytest.approx(0.335)
    iso = json.loads(CANDIDATE.read_text())['metadata']['isostasy']
    assert 'MoI-equality assumption, user-approved 2026-10-22' in \
        iso['moi_target_provenance']
    # Every other campaign config (Titan / Europa clairaut, Callisto, ...)
    # must not carry the flag.
    for p in sorted(CFG_DIR.glob('*.json')):
        if p in (CANDIDATE, CFG_DIR / 'enceladus_cassini_isostasy_ocean6D_rhorock.json'):
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


# ---------------------------------------------------------------------------
# "mass_moi" closure (Option 4, user-approved 2026-10-22)
# ---------------------------------------------------------------------------
V5 = 8.0 * np.pi / 15.0
M_TRUE = 1.08e20


def _stack_M_C(r, rho):
    rin = np.concatenate(([0.0], r[:-1]))
    M = np.sum(rho * V * (r ** 3 - rin ** 3))
    C = np.sum(rho * V5 * (r ** 5 - rin ** 5))
    return M, C / (M * r[-1] ** 2)


def _true_stack(R_s=95e3, R_o=192.22e3, rho_o=1005.0, rho_i=925.0):
    """Rock density chosen so the stack mass is exactly M_TRUE."""
    r = np.array([R_s, R_o, R_T])
    a = 3 * M_TRUE / (4 * np.pi) - rho_o * R_o ** 3 - rho_i * (R_T ** 3 - R_o ** 3)
    rho_s = rho_o + a / R_s ** 3
    return r, np.array([rho_s, rho_o, rho_i])


def test_mass_moi_recovers_mass_and_cmr2_to_1e12():
    r_t, rho_t = _true_stack()
    M, C_t = _stack_M_C(r_t, rho_t)
    assert abs(M / M_TRUE - 1) < 1e-14
    # "PP" stack: wrong R_s (+1.5 km) and a rock density that leaves a 0.3%
    # mass deficit, i.e. inside the 1% guard.
    r_pp = np.array([r_t[0] + 1.5e3, r_t[1], r_t[2]])
    rho_pp = rho_t.copy()
    Mpp = _stack_M_C(r_pp, rho_pp)[0]
    rho_pp[0] *= 1.0 + ((M_TRUE * 0.997 - Mpp)
                        / (V * (r_pp[0] ** 3))) / rho_pp[0]
    s = {'Mtot_kg': M_TRUE}
    r = _runner(s, 'mass_moi', moi_target=C_t)
    out_r, out_rho = r._close_reduced_stack_mass_moi(r_pp, rho_pp, s)
    M2, C2 = _stack_M_C(out_r, out_rho)
    assert abs(M2 / M_TRUE - 1) < 1e-12
    assert abs(C2 / C_t - 1) < 1e-12
    # unique solution == the stack the target was generated from
    assert out_r[0] == pytest.approx(r_t[0], rel=1e-10)
    assert out_rho[0] == pytest.approx(rho_t[0], rel=1e-10)
    assert out_r[1] == r_pp[1] and out_r[2] == r_pp[2]
    assert out_rho[1] == rho_pp[1] and out_rho[2] == rho_pp[2]
    d = r.last_mass_closure
    assert d['feasible'] and d['reason'] is None
    assert d['dR_s_vs_pp_m'] == pytest.approx(-1.5e3, rel=1e-8)
    assert d['drho_s_vs_pp'] == pytest.approx(out_rho[0] - rho_pp[0])
    assert d['moi_stack_after'] == pytest.approx(C_t, rel=1e-12)
    assert d['moi_stack_before'] != pytest.approx(C_t, rel=1e-6)
    assert d['pre_closure_residual_frac'] == pytest.approx(-0.003, abs=1e-9)


def test_mass_moi_ice_override_is_held_exactly():
    r_t, rho_t = _true_stack()
    _, C_t = _stack_M_C(r_t, rho_t)
    s = {'Mtot_kg': M_TRUE}
    r = _runner(s, 'mass_moi', moi_target=C_t)
    out_r, out_rho = r._close_reduced_stack_mass_moi(r_t, rho_t, s, 900.0)
    assert out_rho[2] == 900.0
    M2, C2 = _stack_M_C(out_r, out_rho)
    assert abs(M2 / M_TRUE - 1) < 1e-12 and abs(C2 / C_t - 1) < 1e-12


@pytest.mark.parametrize('case', [
    'a_not_positive', 'b_not_positive', 'R_s_not_in_0_Ro',
    'pre_closure_mass_residual_gt_1pct', 'missing_mass_target',
    'missing_moi_target'])
def test_mass_moi_rejects_every_infeasible_case(case):
    # NB rho_s > rho_o is equivalent to a > 0 (a = (rho_s-rho_o) R_s^3), so
    # the rock-lighter-than-ocean case is the 'a_not_positive' case.
    r_t, rho_t = _true_stack()
    _, C_t = _stack_M_C(r_t, rho_t)
    s = {'Mtot_kg': M_TRUE}
    kw = {'moi_target': C_t}
    rr, rho = r_t.copy(), rho_t.copy()
    if case == 'a_not_positive':
        rho = np.array([500.0, 1005.0, 925.0])     # rock lighter than ocean
        s = {'Mtot_kg': _stack_M_C(rr, rho)[0]}    # 1% guard passes
    elif case == 'b_not_positive':
        kw['moi_target'] = 0.05                    # impossibly condensed
    elif case == 'R_s_not_in_0_Ro':
        kw['moi_target'] = 0.39                    # b/a > R_o^2
    elif case == 'pre_closure_mass_residual_gt_1pct':
        s = {'Mtot_kg': M_TRUE * 1.02}
    elif case == 'missing_mass_target':
        s = {}
    elif case == 'missing_moi_target':
        kw = {}
    r = _runner(s, 'mass_moi', **kw)
    assert r._close_reduced_stack_mass_moi(rr, rho, s) is None
    assert r.last_mass_closure['feasible'] is False
    assert r.last_mass_closure['reason'] == case


def test_mass_moi_infeasible_node_rejected_by_both_paths():
    s, _ = _struct(0.0)
    s = dict(s, Mtot_kg=s['Mtot_kg'] * 1.02)      # broken node (>1%)
    r = _runner(s, 'mass_moi', moi_target=0.335)
    assert r._derive_gravity_isostatic(THETA) is None
    assert r.last_mass_closure['reason'] == 'pre_closure_mass_residual_gt_1pct'
    r = _runner(s, 'mass_moi', moi_target=0.335)
    assert r._derive_libration_deg(THETA) is None


def test_mass_moi_gravity_and_libration_use_closed_stack():
    s, _ = _struct(0.0)
    # moi target the synthetic stack already satisfies -> closure is the
    # identity and mass_moi must reproduce the unclosed result exactly.
    first_r = np.array([R_T - 60e3, R_T - 20e3, R_T])
    r0 = _runner(s)
    # reduced stack of _struct: R_sil=R_T-60km, uniform layers
    rho = s['rho']
    ph = s['phases']
    rho_red = np.array([rho[ph == 3][0], rho[ph == 0][0], rho[ph == 1][0]])
    c_self = _stack_M_C(first_r, rho_red)[1]
    calls = []
    for fn in ('_derive_gravity_isostatic', '_derive_libration_deg'):
        r = _runner(s, 'mass_moi', moi_target=c_self)
        orig = r._close_reduced_stack_mass_moi

        def spy(*a, _orig=orig, **k):
            out = _orig(*a, **k)
            calls.append(out)
            return out
        r._close_reduced_stack_mass_moi = spy
        val = getattr(r, fn)(THETA)
        assert val is not None and calls[-1] is not None
        assert r.last_mass_closure['feasible']
        assert r.last_mass_closure['moi_stack_after'] == pytest.approx(
            c_self, rel=1e-12)
        base = getattr(r0, fn)(THETA)
        if isinstance(val, dict):
            for lm in val:
                assert val[lm] == pytest.approx(base[lm], rel=1e-9)
        else:
            assert val == pytest.approx(base, rel=1e-9)
    # and a target that moves R_s changes both channels
    r = _runner(s, 'mass_moi', moi_target=c_self * 1.01)
    g = r._derive_gravity_isostatic(THETA)
    assert abs(r.last_mass_closure['dR_s_vs_pp_m']) > 100.0
    base = r0._derive_gravity_isostatic(THETA)
    assert g[(2, 0)] != pytest.approx(base[(2, 0)], rel=1e-6)
    lib = _runner(s, 'mass_moi', moi_target=c_self * 1.01)._derive_libration_deg(THETA)
    assert lib != pytest.approx(r0._derive_libration_deg(THETA), rel=1e-6)


def test_rock_density_mode_unchanged_by_mass_moi_addition():
    s, _ = _struct(deficit_frac=0.0046)
    r = _runner(s, 'rock_density', moi_target=0.335)   # target ignored
    g = r._derive_gravity_isostatic(THETA)
    assert 'mode' not in r.last_mass_closure
    assert set(r.last_mass_closure) == {'pre_closure_residual_frac',
                                        'rock_density_factor'}
    assert g is not None


def _struct_rsil(R_sil, mass_frac=0.997):
    """Cache-style structure with a chosen R_sil; rock density leaves a
    (1 - mass_frac) deficit vs the common target M (inside the 1% guard)."""
    R_b = R_T - 20e3
    r_rock = np.linspace(R_sil / 20, R_sil, 20)
    r_oc = np.linspace(R_sil, R_b, 10)[1:]
    r_sh = np.linspace(R_b, R_T, 6)[1:]
    r = np.concatenate([r_rock, r_oc, r_sh])
    rho_o, rho_i, M = 1005.0, 925.0, 1.08e20
    m_o = rho_o * V * (R_b ** 3 - R_sil ** 3)
    m_i = rho_i * V * (R_T ** 3 - R_b ** 3)
    rho_r = (M * mass_frac - m_o - m_i) / (V * R_sil ** 3)
    rho = np.concatenate([np.full(r_rock.size, rho_r),
                          np.full(r_oc.size, rho_o),
                          np.full(r_sh.size, rho_i)])
    ph = np.concatenate([np.full(r_rock.size, 3), np.zeros(r_oc.size, int),
                         np.ones(r_sh.size, int)])
    return {'r_m': r, 'rho': rho, 'phases': ph, 'omega': OMEGA,
            'eccentricity': ECC, 'Mtot_kg': M}


def test_mass_moi_independent_of_candidate_rsil_rhosil():
    """Two PP candidates differing only in (R_sil, rho_sil) -- same M, R_o,
    R, rho_o, rho_i -- close to the SAME rock layer and give identical
    C20/C22/C30 and libration (the point of Option 4)."""
    sA = _struct_rsil(R_T - 60e3, 0.997)
    sB = _struct_rsil(R_T - 61.5e3, 0.9985)
    kw = dict(moi_target=0.335)
    rA, rB = _runner(sA, 'mass_moi', **kw), _runner(sB, 'mass_moi', **kw)
    gA, gB = rA._derive_gravity_isostatic(THETA), rB._derive_gravity_isostatic(THETA)
    assert gA is not None and gB is not None
    for lm in gA:
        assert gA[lm] == pytest.approx(gB[lm], rel=1e-10)
    assert rA.last_mass_closure['R_s_m'] == pytest.approx(
        rB.last_mass_closure['R_s_m'], rel=1e-12)
    assert rA.last_mass_closure['dR_s_vs_pp_m'] != pytest.approx(
        rB.last_mass_closure['dR_s_vs_pp_m'], abs=100.0)
    lA = rA._derive_libration_deg(THETA)
    lB = rB._derive_libration_deg(THETA)
    assert lA == pytest.approx(lB, rel=1e-10)
    # control: without the closure the two candidates DO differ
    oA = _runner(sA)._derive_gravity_isostatic(THETA)
    oB = _runner(sB)._derive_gravity_isostatic(THETA)
    assert oA[(2, 0)] != pytest.approx(oB[(2, 0)], rel=1e-6)


def test_moi_target_consistency_guard():
    s, _ = _struct(0.0)
    tmpl = 'PlanetProfile.Default.Enceladus.PPEnceladus'
    cfg = InferenceConfig.from_json(str(CANDIDATE))
    assert cfg.planet_template_module == tmpl
    r = object.__new__(MCMCRunner)
    r.config = cfg
    r.structure_data = s
    assert r._stack_moi_target() == pytest.approx(0.335)   # consistent
    # mismatched target vs template Cmeasured -> raises
    r = _runner(s, 'mass_moi', moi_target=0.336)
    r.config.planet_template_module = tmpl
    with pytest.raises(ValueError, match='differs from'):
        r._stack_moi_target()
    # bulk_overrides Cmeasured takes precedence over the template
    r = _runner(s, 'mass_moi', moi_target=0.336)
    r.config.planet_template_module = tmpl
    r.config.metadata['bulk_overrides'] = {'Cmeasured': 0.336}
    assert r._stack_moi_target() == pytest.approx(0.336)
    r = _runner(s, 'mass_moi', moi_target=0.335)
    r.config.metadata['bulk_overrides'] = {'Cmeasured': 0.336}
    with pytest.raises(ValueError):
        r._stack_moi_target()
    # raised through the forward model too (not swallowed into None)
    r = _runner(s, 'mass_moi', moi_target=0.336)
    r.config.planet_template_module = tmpl
    with pytest.raises(ValueError):
        r._derive_gravity_isostatic(THETA)


def test_mass_moi_loose_rejects_grain_max_and_dR_10km():
    r_t, rho_t = _true_stack()
    _, C_t = _stack_M_C(r_t, rho_t)
    s = {'Mtot_kg': M_TRUE}
    # baseline feasible (no grain key -> check skipped)
    r = _runner(s, 'mass_moi', moi_target=C_t)
    assert r._close_reduced_stack_mass_moi(r_t, rho_t, s) is not None
    # rho_s above the grain maximum
    rho_s = r.last_mass_closure['rho_s']
    r = _runner(s, 'mass_moi', moi_target=C_t,
                rho_grain_max_kgm3=rho_s - 1.0)
    assert r._close_reduced_stack_mass_moi(r_t, rho_t, s) is None
    assert r.last_mass_closure['reason'] == 'rho_s_above_grain_max'
    assert r.last_mass_closure['feasible'] is False
    r = _runner(s, 'mass_moi', moi_target=C_t,
                rho_grain_max_kgm3=rho_s + 1.0)
    assert r._close_reduced_stack_mass_moi(r_t, rho_t, s) is not None
    # |R_s - R_sil(PP)| > 10 km: PP stack with R_sil 11 km off the closed one
    r_pp = np.array([r_t[0] + 11e3, r_t[1], r_t[2]])
    rho_pp = rho_t.copy()
    Mpp = _stack_M_C(r_pp, rho_pp)[0]
    rho_pp[0] += (M_TRUE - Mpp) / (V * r_pp[0] ** 3)    # mass-consistent
    r = _runner(s, 'mass_moi', moi_target=C_t)
    assert r._close_reduced_stack_mass_moi(r_pp, rho_pp, s) is None
    assert r.last_mass_closure['reason'] == 'dR_s_gt_10km'
    assert abs(r.last_mass_closure['dR_s_vs_pp_m']) > 10e3
    # 9 km is accepted
    r_pp = np.array([r_t[0] + 9e3, r_t[1], r_t[2]])
    rho_pp = rho_t.copy()
    Mpp = _stack_M_C(r_pp, rho_pp)[0]
    rho_pp[0] += (M_TRUE - Mpp) / (V * r_pp[0] ** 3)
    assert r._close_reduced_stack_mass_moi(r_pp, rho_pp, s) is not None


def test_config_carries_grain_max_and_provenance():
    iso = json.loads(CANDIDATE.read_text())['metadata']['isostasy']
    assert 2000.0 < iso['rho_grain_max_kgm3'] < 4000.0
    assert 'Comet_67P' in iso['rho_grain_max_provenance']
    for k in ('mass_moi_gate_result_2026_10_23',
              'sigma_model_numerical_jitter_dropped',
              'mass_moi_ablation_statement',
              'moi_target_uncertainty_treatment'):
        assert k in iso


# ---------------------------------------------------------------------------
# "mass_rho_rock" closure (user decision 2026-10-23): rho_s = sampled
# rho_rock_kgm3 ~ U[2200, 2600]; R_s solved from mass; C/MR^2 DERIVED.
# ---------------------------------------------------------------------------
CAND_8D = CFG_DIR / 'enceladus_cassini_isostasy_ocean6D_rhorock.json'
EDGES = (2200.0, 2600.0)


def _enc_stack(R_o=R_T - 25e3, rho_o=1005.0, rho_i=925.0, R_sil=None):
    """Enceladus-like PP stack consistent with M_TRUE (rock R_sil chosen
    so rho_sil ~ 2450)."""
    if R_sil is None:
        R_sil = 192.0e3
    r = np.array([R_sil, R_o, R_T])
    a = 3 * M_TRUE / (4 * np.pi) - rho_o * R_o ** 3 - rho_i * (R_T ** 3 - R_o ** 3)
    return r, np.array([rho_o + a / R_sil ** 3, rho_o, rho_i])


def _rr_runner(s, **kw):
    return _runner(s, 'mass_rho_rock', **kw)


def test_mass_rho_rock_recovers_mass_and_derived_cmr2_matches_shell_sum():
    r_pp, rho_pp = _enc_stack()
    rho_pp = rho_pp.copy()
    rho_pp[0] *= 0.997                       # 0.3% deficit, inside the guard
    s = {'Mtot_kg': M_TRUE}
    for rho_s in (2200.0, 2400.0, 2600.0):
        r = _rr_runner(s)
        out_r, out_rho = r._close_reduced_stack_mass_rho_rock(
            r_pp, rho_pp, s, None, rho_s)
        M, C = _stack_M_C(out_r, out_rho)          # independent shell sum
        assert abs(M / M_TRUE - 1) < 1e-12
        assert out_rho[0] == rho_s
        assert out_r[1] == r_pp[1] and out_r[2] == r_pp[2]
        assert out_rho[1] == rho_pp[1] and out_rho[2] == rho_pp[2]
        d = r.last_mass_closure
        assert d['feasible'] and d['reason'] is None
        assert d['mode'] == 'mass_rho_rock'
        assert d['rho_s'] == rho_s
        assert d['R_s_m'] == pytest.approx(out_r[0], rel=1e-14)
        assert d['cmr2_derived'] == pytest.approx(C, rel=1e-12)
        assert d['dR_s_vs_pp_m'] == pytest.approx(out_r[0] - r_pp[0])
        assert d['ocean_thickness_m'] == pytest.approx(r_pp[1] - out_r[0])
        assert abs(d['mass_residual_after']) < 1e-12
        # closed form: (rho_s - rho_o) R_s^3 = a
        a = 3 * M_TRUE / (4 * np.pi) - rho_pp[1] * r_pp[1] ** 3 \
            - rho_pp[2] * (R_T ** 3 - r_pp[1] ** 3)
        assert (rho_s - rho_pp[1]) * out_r[0] ** 3 == pytest.approx(a, rel=1e-12)


def test_mass_rho_rock_ice_override_held_exactly():
    r_pp, rho_pp = _enc_stack()
    s = {'Mtot_kg': M_TRUE}
    r = _rr_runner(s)
    out_r, out_rho = r._close_reduced_stack_mass_rho_rock(
        r_pp, rho_pp, s, 935.0, 2400.0)
    assert out_rho[2] == 935.0
    assert abs(_stack_M_C(out_r, out_rho)[0] / M_TRUE - 1) < 1e-12


def test_mass_rho_rock_derived_cmr2_monotone_decreasing_in_rho_s():
    r_pp, rho_pp = _enc_stack()
    s = {'Mtot_kg': M_TRUE}
    cs, rs = [], []
    for rho_s in np.linspace(2200.0, 2600.0, 9):
        r = _rr_runner(s)
        assert r._close_reduced_stack_mass_rho_rock(
            r_pp, rho_pp, s, None, float(rho_s)) is not None
        cs.append(r.last_mass_closure['cmr2_derived'])
        rs.append(r.last_mass_closure['R_s_m'])
    assert all(np.diff(cs) < 0)              # denser rock -> more central
    assert all(np.diff(rs) < 0)              # and smaller rock core


def test_mass_rho_rock_feasible_at_both_box_edges_synthetic_enceladus():
    """Enceladus-like node (zb ~ 25 km, ocean ~ 36 km thick at the PP
    candidate): both prior-box edges feasible; ocean thickness > 0."""
    r_pp, rho_pp = _enc_stack(R_o=R_T - 25e3, R_sil=R_T - 25e3 - 36e3)
    s = {'Mtot_kg': M_TRUE}
    for rho_s in EDGES:
        r = _rr_runner(s, rho_grain_max_kgm3=2858.1098)
        assert r._close_reduced_stack_mass_rho_rock(
            r_pp, rho_pp, s, 925.0, rho_s) is not None, r.last_mass_closure
        d = r.last_mass_closure
        assert d['feasible'] and d['ocean_thickness_m'] > 0.0
        assert 0.30 < d['cmr2_derived'] < 0.37


@pytest.mark.parametrize('case', [
    'missing_mass_target', 'missing_rho_rock', 'nonfinite_rho_rock',
    'pre_closure_mass_residual_gt_1pct', 'a_not_positive',
    'rho_s_not_gt_rho_o', 'R_s_not_in_0_Ro', 'rho_s_above_grain_max'])
def test_mass_rho_rock_rejects_every_infeasible_case_with_reason(case):
    r_pp, rho_pp = _enc_stack()
    s = {'Mtot_kg': M_TRUE}
    kw, rho_s, rr, rho = {}, 2400.0, r_pp.copy(), rho_pp.copy()
    if case == 'missing_mass_target':
        s = {}
    elif case == 'missing_rho_rock':
        rho_s = None
    elif case == 'nonfinite_rho_rock':
        rho_s = float('nan')
    elif case == 'pre_closure_mass_residual_gt_1pct':
        s = {'Mtot_kg': M_TRUE * 1.02}
    elif case == 'a_not_positive':
        # dense ocean over a large volume: ice+ocean outweigh M (a <= 0);
        # M consistent with the PP stack so the 1% guard passes.
        rho = np.array([100.0, 4500.0, 925.0])
        rr = np.array([100e3, R_T - 25e3, R_T])
        s = {'Mtot_kg': _stack_M_C(rr, rho)[0]}
    elif case == 'rho_s_not_gt_rho_o':
        rho_s = 1000.0
    elif case == 'R_s_not_in_0_Ro':
        rho_s = 1020.0                       # barely denser than ocean
    elif case == 'rho_s_above_grain_max':
        kw = {'rho_grain_max_kgm3': 2300.0}
    r = _rr_runner(s, **kw)
    assert r._close_reduced_stack_mass_rho_rock(rr, rho, s, None, rho_s) is None
    assert r.last_mass_closure['feasible'] is False
    assert r.last_mass_closure['reason'] == (
        'missing_rho_rock' if case == 'nonfinite_rho_rock' else case)


def test_mass_rho_rock_does_not_apply_10km_reject_but_records_it():
    """Pending ruling: dR_s is recorded, never rejected on, in this mode."""
    r_pp, rho_pp = _enc_stack()
    s = {'Mtot_kg': M_TRUE}
    r_far = r_pp.copy()
    r_far[0] -= 14e3                          # PP R_sil 14 km below the closure
    rho_far = rho_pp.copy()
    rho_far[0] += (M_TRUE - _stack_M_C(r_far, rho_far)[0]) / (V * r_far[0] ** 3)
    r = _rr_runner(s)
    out = r._close_reduced_stack_mass_rho_rock(r_far, rho_far, s, None, 2400.0)
    assert out is not None
    assert abs(r.last_mass_closure['dR_s_vs_pp_m']) > 10e3
    # same stack under mass_moi-style 10 km rule would reject: control
    from PlanetProfile.Inference import mcmc_runner as mr
    assert mr._MASS_RHO_ROCK_APPLY_DR_S_REJECT is False


def test_mass_rho_rock_independent_of_candidate_rsil_rhosil():
    """Two PP candidates differing only in (R_sil, rho_sil): identical
    C20/C22/C30, libration and closed stack, for each rho_s."""
    sA = _struct_rsil(R_T - 60e3, 0.997)
    sB = _struct_rsil(R_T - 61.5e3, 0.9985)
    for rho_s in (2200.0, 2400.0, 2600.0):
        th = dict(THETA, rho_rock_kgm3=rho_s)
        rA, rB = _rr_runner(sA), _rr_runner(sB)
        gA, gB = rA._derive_gravity_isostatic(th), rB._derive_gravity_isostatic(th)
        assert gA is not None and gB is not None
        for lm in gA:
            assert gA[lm] == pytest.approx(gB[lm], rel=1e-10)
        assert rA.last_mass_closure['R_s_m'] == pytest.approx(
            rB.last_mass_closure['R_s_m'], rel=1e-12)
        assert rA.last_mass_closure['dR_s_vs_pp_m'] != pytest.approx(
            rB.last_mass_closure['dR_s_vs_pp_m'], abs=100.0)
        assert rA._derive_libration_deg(th) == pytest.approx(
            rB._derive_libration_deg(th), rel=1e-10)


def test_mass_rho_rock_honoured_identically_by_gravity_and_libration():
    s, _ = _struct(0.0)
    for rho_s in EDGES + (2400.0,):
        th = dict(THETA, rho_rock_kgm3=rho_s, rho_ice_kgm3=930.0)
        rg = _rr_runner(s)
        assert rg._derive_gravity_isostatic(th) is not None
        dg = dict(rg.last_mass_closure)
        rl = _rr_runner(s)
        assert rl._derive_libration_deg(th) is not None
        dl = dict(rl.last_mass_closure)
        assert dg['feasible'] and dl['feasible']
        assert dg['rho_s'] == dl['rho_s'] == rho_s
        assert dg['R_s_m'] == dl['R_s_m']
        assert dg['cmr2_derived'] == dl['cmr2_derived']
        # the closure really moves both channels with rho_s
    gs = [_rr_runner(s)._derive_gravity_isostatic(
        dict(THETA, rho_rock_kgm3=x)) for x in EDGES]
    ls = [_rr_runner(s)._derive_libration_deg(
        dict(THETA, rho_rock_kgm3=x)) for x in EDGES]
    assert gs[0][(2, 0)] != pytest.approx(gs[1][(2, 0)], rel=1e-4)
    assert ls[0] != pytest.approx(ls[1], rel=1e-4)


def test_mass_rho_rock_node_rejected_by_both_paths_no_fallback():
    s, _ = _struct(0.0)
    th = dict(THETA, rho_rock_kgm3=2400.0)
    s_bad = dict(s, Mtot_kg=s['Mtot_kg'] * 1.02)      # >1% broken node
    r = _rr_runner(s_bad)
    assert r._derive_gravity_isostatic(th) is None
    assert r.last_mass_closure['reason'] == 'pre_closure_mass_residual_gt_1pct'
    assert _rr_runner(s_bad)._derive_libration_deg(th) is None
    # missing rho_rock in theta -> reject, never a default/fallback value
    r = _rr_runner(s)
    assert r._derive_gravity_isostatic(THETA) is None
    assert r.last_mass_closure['reason'] == 'missing_rho_rock'
    assert _rr_runner(s)._derive_libration_deg(THETA) is None
    # infeasible rho_s -> None through the public paths
    r = _rr_runner(s)
    assert r._derive_gravity_isostatic(dict(THETA, rho_rock_kgm3=900.0)) is None
    assert r.last_mass_closure['reason'] == 'rho_s_not_gt_rho_o'


def test_mass_rho_rock_guard_moi_target_and_mode_both_live_raises():
    s, _ = _struct(0.0)
    r = _runner(s, 'mass_rho_rock', moi_target=0.335)
    with pytest.raises(ValueError, match='double-count'):
        r._stack_mass_closure_mode()
    with pytest.raises(ValueError):
        r._derive_gravity_isostatic(dict(THETA, rho_rock_kgm3=2400.0))
    # the Cmeasured consistency guard must not fire in the new mode
    r = _rr_runner(s)
    r.config.planet_template_module = 'PlanetProfile.Default.Enceladus.PPEnceladus'
    r.config.metadata['bulk_overrides'] = {'Cmeasured': 0.999}
    assert r._derive_gravity_isostatic(dict(THETA, rho_rock_kgm3=2400.0)) is not None
    # mass_moi still honours (and checks) moi_target
    r = _runner(s, 'mass_moi', moi_target=0.335)
    assert r._stack_mass_closure_mode() == 'mass_moi'


def test_ocean_theta_with_rho_rock_requires_mass_rho_rock_mode():
    s, _ = _struct(0.0)
    ps = {'zb_km': {'prior_type': 'uniform', 'bounds': [5.0, 45.0]},
          'rho_rock_kgm3': {'prior_type': 'uniform', 'bounds': [2200.0, 2600.0]}}
    for mode, ok in ((None, False), ('mass_moi', False),
                     ('rock_density', False), ('mass_rho_rock', True)):
        r = _runner(s, mode)
        r.config.param_space = ps
        r.param_names, r.fixed_params = list(ps), {}
        if ok:
            r._check_ocean_rho_rock_config()
        else:
            with pytest.raises(ValueError, match='silently ignored'):
                r._check_ocean_rho_rock_config()
    # frozen-style theta (rho_rock without zb_km) is unaffected
    r = _runner(s, None)
    r.param_names, r.fixed_params = ['rho_rock_kgm3', 'rho_ice_kgm3'], {}
    r._check_ocean_rho_rock_config()


def test_dispatch_ocean_theta_with_rho_rock_is_not_routed_to_frozen(monkeypatch):
    from PlanetProfile.Inference import mcmc_runner as mr
    s, _ = _struct(0.0)
    r = _rr_runner(s)
    r.structure_data = {'frozen_branch_supported': True}
    monkeypatch.setattr(mr, 'is_zbw_cache', lambda sd: True)
    r._frozen_struct_for_theta = lambda t: 'FROZEN'
    r._ocean_zbw_struct_for_theta = lambda t: 'OCEAN'
    assert r._struct_for_hydrosphere(
        {'zb_km': 20.0, 'log10_wOcean_ppt': 0.5, 'rho_rock_kgm3': 2400.0,
         'rho_ice_kgm3': 925.0}) == 'OCEAN'
    assert r._struct_for_hydrosphere(
        {'rho_rock_kgm3': 2400.0, 'rho_ice_kgm3': 925.0}) == 'FROZEN'


def _frozen_runner(mode, **kw):
    from PlanetProfile.Gravity.isostasy import frozen_rho_rock_from_zb
    zb_m, rho_ice = 55e3, 925.0
    rho_rock = frozen_rho_rock_from_zb(M_TRUE, R_T, zb_m, rho_ice)
    R_rock = R_T - zb_m
    st = {'r_m': np.array([R_rock, R_rock + 0.5 * zb_m, R_T]),
          'rho': np.array([rho_rock, rho_ice, rho_ice]),
          'phases': np.array([50, 1, 1]), 'omega': OMEGA, 'eccentricity': ECC,
          'Mtot_kg': M_TRUE}
    r = _runner(st, mode, **kw)
    return r, st


def test_frozen_branch_unchanged_by_mass_rho_rock_mode():
    """Frozen-branch forward predictions (the only mode-dependent inputs of
    its log L) are bit-identical for mode none / mass_rho_rock / mass_moi,
    and a sampled rho_rock_kgm3 keeps its existing frozen meaning (it
    selects the node via the zb bijection, checked by I-F5, not by the
    closure)."""
    r0, st = _frozen_runner(None)
    base_g = r0._derive_gravity_isostatic(THETA)
    base_l = r0._derive_libration_deg(THETA)
    assert base_g is not None and base_l is not None
    for mode, kw in (('mass_rho_rock', {}), ('mass_moi', {'moi_target': 0.335})):
        r, _ = _frozen_runner(mode, **kw)
        assert r._derive_gravity_isostatic(THETA) == base_g
        assert r._derive_libration_deg(THETA) == base_l
        # sampled rho_rock keeps its existing frozen meaning (I-F5 loop
        # closure against the served stack): the stack's own interior
        # density is accepted, a different one is rejected -- in every mode.
        own = float(st['rho'][0])
        th = dict(THETA, rho_rock_kgm3=own)
        assert r._derive_gravity_isostatic(th) == base_g
        assert r._derive_libration_deg(th) == base_l
        th = dict(THETA, rho_rock_kgm3=own * 1.05)
        assert r._derive_gravity_isostatic(th) is None
        assert r0._derive_gravity_isostatic(th) is None


def test_8d_config_wiring_and_7d_untouched():
    c7 = json.loads(CANDIDATE.read_text())
    c8 = json.loads(CAND_8D.read_text())
    assert c7['metadata']['isostasy']['stack_mass_closure'] == 'mass_moi'
    assert 'rho_rock_kgm3' not in c7['param_space']
    iso = c8['metadata']['isostasy']
    assert iso['stack_mass_closure'] == 'mass_rho_rock'
    assert 'moi_target' not in iso
    assert iso['moi_iess2014_reference_display_only'] == 0.335
    p = c8['param_space']['rho_rock_kgm3']
    assert p['prior_type'] == 'uniform' and p['bounds'] == [2200.0, 2600.0]
    from PlanetProfile.Inference.parameter_registry import PARAMETER_REGISTRY
    assert (list(PARAMETER_REGISTRY['rho_rock_kgm3'].default_bounds)
            == p['bounds'])
    assert 'rho_rock_kgm3' in c8['metadata']['branch_model']['ocean_branch'][
        'sampled']
    # same keys as 7D except the documented ones
    assert set(c7['param_space']) | {'rho_rock_kgm3'} == set(c8['param_space'])
    txt = json.dumps(c8['metadata'])
    assert 'ADDRESSED BY THIS CHANGE, PENDING VERIFICATION' in txt
    assert 'NOT independent confirmation' in iso['stack_mass_closure_provenance'] \
        or 'not independent confirmation' in \
        iso['stack_mass_closure_provenance'].lower()
    assert '2026-10-23' in iso['stack_mass_closure_provenance']
    cfg = InferenceConfig.from_json(str(CAND_8D))
    r = object.__new__(MCMCRunner)
    r.config = cfg
    assert r._stack_mass_closure_mode() == 'mass_rho_rock'
    r.param_names, r.fixed_params = list(cfg.param_space), {}
    r._check_ocean_rho_rock_config()
