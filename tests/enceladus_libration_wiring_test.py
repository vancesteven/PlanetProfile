"""D3/D4/D5/D6 -- the libration channel, wired at the PRODUCTION call sites.

r5 likelihood defects (validation_reports/enceladus_isostasy/
r5_ADJUDICATION.md), all on the campaign's SOLE shell-thickness channel:

  D3  the ruled B2' Delta_rho-consistent Eq.-12 figure treatment existed in
      ``Gravity/Librations.py`` but was wired at NEITHER production call
      site, so the likelihood ran the historical hydrostatic convention
      (1.3-2.6 km on the headline deliverable; dlibration/dC2 identically
      zero as shipped).
  D4  ``libration_model_correction.multiplicative_factor`` (1.008, the
      deterministic +0.8% B2 formulation-bias correction) was never applied.
  D5  ``libration_sys_frac`` was sampled and registered but never entered
      chi2 -- a dead parameter by the campaign's own CRITICAL-1 criterion.
  D6  ``rho_ice_kgm3`` never reached the OCEAN branch's libration
      (``rho_shell_override`` was unreachable from the likelihood) while it
      already reached the frozen branch's libration and the gravity channel.

EVERYTHING HERE GOES THROUGH ``MCMCRunner._derive_libration_deg`` (and, for
the likelihood deltas, through ``log_likelihood_fn``) -- never through
``librations()`` directly. That is the point: ``tests/librations_test.py``
already pins the same bracket at the module level, and it passed the whole
time the production dispatch was running a different convention.

The pinned targets are the adjudicated values in
``b2prime_ADJUDICATED_drho_weighting.json`` CRITICAL_2.reviewer_measurement,
at that file's own fiducial (zb = 25 km, D_ocean = 36 km, rho_ocean = 1005,
H22_obs = 857 m).
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

# B2' fiducial (b2prime_ADJUDICATED_drho_weighting.json).
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


def _b2prime_struct(zb_km, D_ocean_km=D_OCEAN_KM, rho_ocean=RHO_OCEAN):
    """Mass-conserved 3-layer structure, identical construction to
    tests/librations_test.py::_enceladus_b2prime_structure -- but returned
    as a CACHE STRUCTURE DICT so it is served through the production
    dispatch (_struct_for_hydrosphere -> _derive_libration_deg)."""
    R_b = R_T - zb_km * 1e3
    R_c = R_b - D_ocean_km * 1e3
    V = 4.0 * np.pi / 3.0
    m_shell = RHO_ICE * V * (R_T ** 3 - R_b ** 3)
    m_ocean = rho_ocean * V * (R_b ** 3 - R_c ** 3)
    rho_core = (M_ENC - m_shell - m_ocean) / (V * R_c ** 3)
    return {
        'r_m': np.array([R_c, R_b, R_T]),
        'rho': np.array([rho_core, rho_ocean, RHO_ICE]),
        'phases': np.array([50, 0, 1]),
        'omega': OMEGA,
        'eccentricity': ECC,
        'Mtot_kg': M_ENC,
    }


def _config_libration_block(**overrides):
    """The RULED block, read from the committed config so this test cannot
    drift from the record, with targeted overrides."""
    with open(CANDIDATE_CONFIG) as fh:
        blk = dict(json.load(fh)['metadata']['libration_model_correction'])
    blk.update(overrides)
    return blk


def _runner(struct, libration_block=None, observables=None):
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
        observables=observables or {'libration_deg': (LIB_OBS_DEG,
                                                      SIGMA_LIB_DEG)},
        sampler_settings={},
        gravity_forward_model='isostatic_hm2019',
        metadata=metadata,
    )
    runner = object.__new__(MCMCRunner)
    runner.config = config
    runner.structure_data = struct
    return runner


def _lib_deg(zb_km, theta, libration_block=None, **struct_kwargs):
    """Libration [deg] THROUGH THE PRODUCTION DISPATCH."""
    return _runner(_b2prime_struct(zb_km, **struct_kwargs),
                   libration_block)._derive_libration_deg(theta)


# Uncorrected model (multiplicative_factor forced to 1.0) -- the adjudicated
# bracket is a statement about the MODEL libration, before the deterministic
# +0.8% budget line D4 adds on top.
_UNCORRECTED = _config_libration_block(multiplicative_factor=1.0)


# ===========================================================================
# D3 -- the ruled figure treatment reaches the production call site.
# ===========================================================================

def test_config_declares_the_ruled_figure_convention():
    blk = _config_libration_block()
    assert blk['figure_convention'] == 'drho_consistent_eq12'
    assert blk['H22_obs_m'] == 857.0
    assert blk['multiplicative_factor'] == 1.008


@pytest.mark.parametrize('C2,target_deg,target_sigma', [
    (0.0, 0.10054, 3.25),
    (0.5, 0.09814, 2.46),
    (1.0, 0.09576, 1.66),
])
def test_production_dispatch_reproduces_the_adjudicated_c2_sweep(
        C2, target_deg, target_sigma):
    """b2prime_ADJUDICATED_drho_weighting.json
    CRITICAL_2.corrected_Eq12_C2_sweep_at_zb25_D36_rho1005, reproduced
    through _derive_libration_deg rather than through librations()."""
    theta = {'compensation_C2': C2}
    lib_hyd = _lib_deg(25.0, theta, None)          # historical convention
    lib_ruled = _lib_deg(25.0, theta, _UNCORRECTED)  # ruled treatment
    assert lib_ruled == pytest.approx(target_deg, abs=5e-6)
    assert ((lib_ruled - lib_hyd) / SIGMA_LIB_DEG
            == pytest.approx(target_sigma, abs=5e-3))


def test_shipped_hydrostatic_path_is_what_the_defect_left_running():
    """The defect, made explicit: with no block declared the dispatch
    returns the hydrostatic value, which is 1.66-3.25 sigma away from the
    ruled treatment and is NOT any of the adjudicated numbers."""
    lib_hyd = _lib_deg(25.0, {'compensation_C2': 1.0}, None)
    for adjudicated in (0.10054, 0.09814, 0.09576):
        assert abs(lib_hyd - adjudicated) > 4e-3 * SIGMA_LIB_DEG * 100


def _zb_matching_observed(C2, libration_block=_UNCORRECTED):
    """Shell thickness reproducing the Park libration, solved THROUGH the
    production dispatch."""
    def _resid(zb_km):
        return _lib_deg(zb_km, {'compensation_C2': C2},
                        libration_block) - LIB_OBS_DEG
    return brentq(_resid, 20.0, 35.0, xtol=1e-6)


@pytest.mark.parametrize('C2,target_zb_km', [(0.0, 27.34), (1.0, 25.99)])
def test_production_dispatch_reproduces_the_adjudicated_zb_bracket(
        C2, target_zb_km):
    """THE RULED BRACKET: shell thickness matching the conditioned Park
    libration 0.092 deg at D_ocean = 36 km is bracketed by C2 between
    27.34 km (C2 = 0, identical to surface-only) and 25.99 km (C2 = 1)
    (config metadata.libration_model_correction.effect)."""
    assert _zb_matching_observed(C2) == pytest.approx(target_zb_km, abs=0.005)


def test_bracket_is_a_real_shift_off_the_hydrostatic_answer():
    """The hydrostatic convention the defect left running answers 24.70 km
    -- 1.3 to 2.6 km thinner than the ruled bracket."""
    zb_hyd = _zb_matching_observed(1.0, libration_block=None)
    assert zb_hyd == pytest.approx(24.70, abs=0.01)
    assert 25.99 - zb_hyd > 1.28
    assert 27.34 - zb_hyd > 2.63


def test_dlibration_dC2_is_nonzero_through_the_dispatch():
    """As shipped, dlibration/dC2 was IDENTICALLY zero (the Airy root never
    reached the torque). Under the ruled treatment C2 is shared between the
    gravity and libration channels, at a REAL BUT SMALL 1.6 sigma across
    C2 in [0, 1] (the ~30 sigma span once recorded is retracted)."""
    theta = lambda c2: {'compensation_C2': c2}  # noqa: E731
    # The defect: exactly zero, not merely small.
    assert (_lib_deg(26.0, theta(0.0), None)
            == _lib_deg(26.0, theta(1.0), None))
    lo = _lib_deg(26.0, theta(0.0), _UNCORRECTED)
    hi = _lib_deg(26.0, theta(1.0), _UNCORRECTED)
    assert lo != hi
    span_sigma = (lo - hi) / SIGMA_LIB_DEG
    assert 1.4 < span_sigma < 1.8, span_sigma
    # Local derivative, nonzero and negative (more compensation -> smaller
    # libration).
    d = (_lib_deg(26.0, theta(0.55), _UNCORRECTED)
         - _lib_deg(26.0, theta(0.45), _UNCORRECTED)) / 0.1
    assert d < -1e-4


# ===========================================================================
# D4 -- the deterministic +0.8% correction.
# ===========================================================================

def test_libration_model_correction_is_the_declared_multiplicative_factor():
    zb = _zb_matching_observed(1.0)
    theta = {'compensation_C2': 1.0}
    uncorrected = _lib_deg(zb, theta, _UNCORRECTED)
    corrected = _lib_deg(zb, theta, _config_libration_block())
    assert corrected / uncorrected == pytest.approx(1.008, rel=1e-12)


def test_correction_is_plus_0245_sigma_at_the_solved_bracket():
    """+0.8% of a prediction sitting at the conditioned 0.092 deg is
    +0.245 sigma_obs -- the number the config's B2 provenance budgets."""
    zb = _zb_matching_observed(1.0)
    theta = {'compensation_C2': 1.0}
    uncorrected = _lib_deg(zb, theta, _UNCORRECTED)
    corrected = _lib_deg(zb, theta, _config_libration_block())
    assert uncorrected == pytest.approx(LIB_OBS_DEG, abs=1e-6)
    shift_sigma = (corrected - uncorrected) / SIGMA_LIB_DEG
    assert shift_sigma == pytest.approx(0.24533, abs=1e-4)


def test_correction_shows_in_a_likelihood_delta():
    """Same shift, seen where it matters: the log-likelihood."""
    zb = _zb_matching_observed(1.0)
    struct = _b2prime_struct(zb)
    theta = {'compensation_C2': 1.0}

    def _ll(block):
        runner = _runner(struct, block)
        pred = runner._derive_libration_deg(theta)
        return -0.5 * ((pred - LIB_OBS_DEG) / SIGMA_LIB_DEG) ** 2

    ll_uncorrected = _ll(_UNCORRECTED)
    ll_corrected = _ll(_config_libration_block())
    # Uncorrected sits exactly on the observed value -> chi2 = 0; the
    # correction moves it +0.245 sigma off, so the likelihood must DROP by
    # 0.5 * 0.245**2.
    assert ll_uncorrected == pytest.approx(0.0, abs=1e-9)
    assert ll_corrected == pytest.approx(-0.5 * 0.24533 ** 2, abs=1e-4)
    assert ll_corrected < ll_uncorrected


# ===========================================================================
# D5 -- libration_sys_frac reaches chi2.
# ===========================================================================

@pytest.mark.parametrize('s', [-0.012, -0.004, 0.004, 0.012])
def test_libration_sys_frac_scales_the_prediction(s):
    zb = 26.0
    base = _lib_deg(zb, {'compensation_C2': 0.5}, _config_libration_block())
    with_s = _lib_deg(zb, {'compensation_C2': 0.5, 'libration_sys_frac': s},
                      _config_libration_block())
    assert with_s / base == pytest.approx(1.0 + s, rel=1e-12)


def test_libration_sys_frac_is_not_a_dead_parameter_in_chi2():
    """The CRITICAL-1 criterion: a sampled parameter that cannot move the
    likelihood is dead. At its declared 1-sigma (0.004) the nuisance must
    move the libration channel by ~0.12 sigma_obs."""
    zb = _zb_matching_observed(1.0)
    struct = _b2prime_struct(zb)
    blk = _config_libration_block(multiplicative_factor=1.0)

    def _pred(s):
        return _runner(struct, blk)._derive_libration_deg(
            {'compensation_C2': 1.0, 'libration_sys_frac': s})

    p0, pp, pm = _pred(0.0), _pred(0.004), _pred(-0.004)
    assert pp != p0 and pm != p0
    assert (pp - p0) / SIGMA_LIB_DEG == pytest.approx(0.1227, abs=1e-3)
    assert (pm - p0) / SIGMA_LIB_DEG == pytest.approx(-0.1227, abs=1e-3)
    # ... and that difference survives into chi2.
    chi2 = lambda p: ((p - LIB_OBS_DEG) / SIGMA_LIB_DEG) ** 2  # noqa: E731
    assert chi2(pp) > chi2(p0)
    assert chi2(pm) > chi2(p0)


# ===========================================================================
# D6 -- rho_ice_kgm3 reaches the OCEAN branch's libration.
# ===========================================================================

def test_rho_ice_moves_the_ocean_branch_libration():
    """The third branch asymmetry: rho_ice already moved the frozen branch's
    libration and the ocean branch's gravity, but not the ocean branch's
    libration."""
    zb = 26.0
    theta = {'compensation_C2': 0.5}
    blk = _config_libration_block()
    base = _lib_deg(zb, theta, blk)
    lo = _lib_deg(zb, {**theta, 'rho_ice_kgm3': 915.0}, blk)
    hi = _lib_deg(zb, {**theta, 'rho_ice_kgm3': 935.0}, blk)
    assert lo != base and hi != base
    assert lo != hi
    # Monotone across the declared U[915, 935] prior, and a real but
    # sub-sigma systematic on the channel.
    span_sigma = abs(hi - lo) / SIGMA_LIB_DEG
    assert 0.05 < span_sigma < 1.0, span_sigma


def test_rho_ice_identity_override_is_still_a_noop():
    """Mass-neutrality invariant: handing the stack back its OWN shell
    density must not move the answer (the bug that once injected -0.19
    sigma merely by enabling the nuisance)."""
    zb = 26.0
    theta = {'compensation_C2': 0.5}
    blk = _config_libration_block()
    struct = _b2prime_struct(zb)
    base = _runner(struct, blk)._derive_libration_deg(theta)
    ident = _runner(struct, blk)._derive_libration_deg(
        {**theta, 'rho_ice_kgm3': RHO_ICE})
    assert ident == pytest.approx(base, rel=1e-12)


def test_explicit_rho_shell_override_still_wins():
    """The kwarg remains the authority for the mass-neutrality unit tests
    that drive it directly."""
    zb = 26.0
    struct = _b2prime_struct(zb)
    blk = _config_libration_block()
    from_theta = _runner(struct, blk)._derive_libration_deg(
        {'compensation_C2': 0.5, 'rho_ice_kgm3': 935.0})
    from_kwarg = _runner(struct, blk)._derive_libration_deg(
        {'compensation_C2': 0.5, 'rho_ice_kgm3': 915.0},
        rho_shell_override=935.0)
    assert from_kwarg == pytest.approx(from_theta, rel=1e-12)


# ===========================================================================
# The no-op invariant: configs without the block are unchanged.
# ===========================================================================

@pytest.mark.parametrize('C2', [0.0, 0.5, 1.0])
def test_absent_block_is_byte_identical_to_the_historical_path(C2):
    """A config that declares no libration_model_correction block and does
    not sample the nuisance keeps PlanetProfile's historical behaviour
    exactly -- which is what protects every other campaign from this
    repair."""
    from PlanetProfile.Gravity.Librations import librations
    zb = 26.0
    struct = _b2prime_struct(zb)
    pred = _runner(struct, None)._derive_libration_deg(
        {'compensation_C2': C2})
    lib_m = librations(struct['r_m'], struct['rho'], OMEGA, ECC,
                       rigid=True, ocean=True, ocean_idx=1)
    assert pred == float(np.degrees(lib_m / struct['r_m'][-1]))


def test_sys_frac_is_independent_of_the_correction_block():
    """libration_sys_frac is a registered SAMPLED parameter in its own
    right, so it applies whenever a sample carries it -- the correction
    block only supplies the deterministic factor. A config that samples the
    nuisance and gets no chi2 response from it is the D5 defect."""
    zb = 26.0
    struct = _b2prime_struct(zb)
    base = _runner(struct, None)._derive_libration_deg(
        {'compensation_C2': 0.5})
    with_s = _runner(struct, None)._derive_libration_deg(
        {'compensation_C2': 0.5, 'libration_sys_frac': 0.01})
    assert with_s / base == pytest.approx(1.01, rel=1e-12)


def test_unknown_figure_convention_does_not_silently_select_a_treatment():
    """The RESCINDED whole-difference treatment and the ruled one are
    different physics; an unrecognized declaration must fall back to the
    hydrostatic path, never guess."""
    zb = 26.0
    theta = {'compensation_C2': 1.0}
    hyd = _lib_deg(zb, theta, None)
    bogus = _lib_deg(zb, theta,
                     _config_libration_block(figure_convention='whole_difference',
                                             multiplicative_factor=1.0))
    assert bogus == hyd
