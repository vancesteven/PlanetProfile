"""Sigma_model additive-variance likelihood inflation (H&M 2019 Eq. 22/24 --
freeze blocker B15 wiring).

Mechanical wiring only -- no physics changes, no scientific judgment calls.
A config MAY declare ``metadata['sigma_model_add'] = {observable_name:
sigma_model}``; ``MCMCRunner._make_flexible_log_likelihood`` then replaces
each named observable's sigma_obs with ``sqrt(sigma_obs**2 +
sigma_model**2)`` before forming chi2 (H&M's Eq. 24 combined-sigma
construction). Configs that do not declare the block -- i.e. every
pre-existing campaign -- get an empty dict back from
``metadata.get('sigma_model_add', {})``, so ``_inflate_sigma`` returns
sigma_obs completely unchanged and the likelihood is bit-identical to the
pre-B15 formula. That is the critical invariant this file protects.

Uses the committed real Enceladus smoke structure cache (3 Tb-grid nodes,
2 frozen + 1 ocean) via the same fixture pattern as
tests/enceladus_isostasy_code_gaps_test.py's
``test_log_likelihood_isostatic_dispatch_on_real_smoke_cache``.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402
from PlanetProfile.Inference.mcmc_runner import MCMCRunner  # noqa: E402

SMOKE_CACHE = (REPO / 'PlanetProfile' / 'Test' / 'mcmc_results' / 'Enceladus'
              / 'Cassini_smoke' / 'enceladus_seawater_smoke_grid.pkl')

C20_OBS = (-5521e-6, 3.58e-4)
C22_OBS = (1574e-6, 1.59e-4)
C30_OBS = (118e-6, 2.38e-4)
LIB_OBS = (0.092, 0.003)

OCEAN_TB = 272.46  # node 2 (the only phase-0 ocean node) per RESUME_NOTE.md


def _config(sigma_model_add=None, include_libration=False):
    metadata = {'isostasy': {
        'H_obs_lm_m': {'2,0': -3510.0, '2,2': 857.0, '3,0': 420.0},
        'shape_ref_radius_m': 252220.0,
        'gravity_ref_radius_m': 256600.0,
        'finite_amplitude': True,
    }}
    if sigma_model_add is not None:
        metadata['sigma_model_add'] = sigma_model_add
    observables = {'C20': C20_OBS, 'C22': C22_OBS, 'C30': C30_OBS}
    if include_libration:
        observables['libration_deg'] = LIB_OBS
    return InferenceConfig(
        mode='mcmc', bodyname='Enceladus',
        param_space={
            'Tb_K': {'prior_type': 'uniform', 'bounds': [271.8, 272.46]},
            'compensation_C2': {'prior_type': 'uniform', 'bounds': [0.0, 1.0]},
        },
        observables=observables,
        sampler_settings={},
        structure_cache_path=str(SMOKE_CACHE),
        gravity_forward_model='isostatic_hm2019',
        planet_template_module='PlanetProfile.Default.Enceladus.PPEnceladus',
        metadata=metadata,
    )


pytestmark = pytest.mark.skipif(
    not SMOKE_CACHE.exists(), reason='Enceladus smoke structure cache not present')


# ===========================================================================
# No-op branch: absent / empty block must be bit-identical to the pre-B15
# formula.
# ===========================================================================

@pytest.mark.parametrize('c2', [0.0, 0.35, 1.0])
def test_absent_and_empty_block_are_bit_identical(c2):
    runner_absent = MCMCRunner(_config(sigma_model_add=None))
    runner_empty = MCMCRunner(_config(sigma_model_add={}))
    theta = np.array([OCEAN_TB, c2])
    ll_absent = runner_absent.log_likelihood_fn(theta)
    ll_empty = runner_empty.log_likelihood_fn(theta)
    assert ll_absent == ll_empty


@pytest.mark.parametrize('c2', [0.0, 0.35, 0.6, 1.0])
def test_absent_block_matches_hand_computed_pre_b15_chi2(c2):
    """The core regression assertion: with no sigma_model_add block, the
    likelihood must equal exactly the old (pre-B15) formula -- independent
    Gaussian terms built directly from the raw observable sigmas, with no
    inflation whatsoever."""
    runner = MCMCRunner(_config(sigma_model_add=None))
    theta_dict = {'Tb_K': OCEAN_TB, 'compensation_C2': c2}
    theta = np.array([OCEAN_TB, c2])

    ll = runner.log_likelihood_fn(theta)
    assert np.isfinite(ll)

    pred = runner._derive_gravity_isostatic(theta_dict)
    assert pred is not None

    chi2 = 0.0
    for name, obs, lm in (('C20', C20_OBS, (2, 0)),
                          ('C22', C22_OBS, (2, 2)),
                          ('C30', C30_OBS, (3, 0))):
        v, s = obs
        chi2 += ((pred[lm] - v) / s) ** 2

    assert ll == pytest.approx(-0.5 * chi2, rel=1e-12, abs=1e-12)


# ===========================================================================
# Inflation branch: declared sigma_model_add must apply
# sqrt(sigma_obs**2 + sigma_model**2) per named observable.
# ===========================================================================

SIGMA_MODEL_GRAVITY = {'C20': 5.3e-06, 'C22': 1.7e-06, 'C30': 4.4e-06}


@pytest.mark.parametrize('c2', [0.0, 0.35, 0.6, 1.0])
def test_sigma_model_add_inflates_exactly(c2):
    runner = MCMCRunner(_config(sigma_model_add=SIGMA_MODEL_GRAVITY))
    theta_dict = {'Tb_K': OCEAN_TB, 'compensation_C2': c2}
    theta = np.array([OCEAN_TB, c2])

    ll = runner.log_likelihood_fn(theta)
    assert np.isfinite(ll)

    pred = runner._derive_gravity_isostatic(theta_dict)
    assert pred is not None

    chi2 = 0.0
    for name, obs, lm in (('C20', C20_OBS, (2, 0)),
                          ('C22', C22_OBS, (2, 2)),
                          ('C30', C30_OBS, (3, 0))):
        v, s = obs
        s_eff = float(np.sqrt(s ** 2 + SIGMA_MODEL_GRAVITY[name] ** 2))
        chi2 += ((pred[lm] - v) / s_eff) ** 2

    assert ll == pytest.approx(-0.5 * chi2, rel=1e-12, abs=1e-12)


def test_sigma_model_add_reduces_chi2_magnitude():
    """Inflating sigma can only shrink |residual / sigma|, so the inflated
    log-likelihood must be >= the uninflated one at the same theta (given a
    nonzero residual)."""
    theta = np.array([OCEAN_TB, 0.6])
    ll_plain = MCMCRunner(_config(sigma_model_add=None)).log_likelihood_fn(theta)
    ll_inflated = MCMCRunner(
        _config(sigma_model_add=SIGMA_MODEL_GRAVITY)).log_likelihood_fn(theta)
    assert ll_inflated > ll_plain


def test_sigma_model_add_only_affects_named_observables():
    """Declaring inflation for C20 only must leave C22/C30 sigma
    untouched -- the mechanism is keyed by observable name, not
    all-or-nothing."""
    runner = MCMCRunner(_config(sigma_model_add={'C20': 5.3e-06}))
    theta_dict = {'Tb_K': OCEAN_TB, 'compensation_C2': 0.4}
    theta = np.array([OCEAN_TB, 0.4])
    ll = runner.log_likelihood_fn(theta)
    pred = runner._derive_gravity_isostatic(theta_dict)

    v20, s20 = C20_OBS
    v22, s22 = C22_OBS
    v30, s30 = C30_OBS
    s20_eff = float(np.sqrt(s20 ** 2 + 5.3e-06 ** 2))
    chi2 = (((pred[(2, 0)] - v20) / s20_eff) ** 2
            + ((pred[(2, 2)] - v22) / s22) ** 2
            + ((pred[(3, 0)] - v30) / s30) ** 2)
    assert ll == pytest.approx(-0.5 * chi2, rel=1e-12, abs=1e-12)


# ===========================================================================
# libration_deg channel (the fourth H&M Eq. 22 observable).
# ===========================================================================

def test_libration_deg_absent_block_matches_hand_computed_chi2():
    runner = MCMCRunner(_config(sigma_model_add=None, include_libration=True))
    theta_dict = {'Tb_K': OCEAN_TB, 'compensation_C2': 0.5}
    theta = np.array([OCEAN_TB, 0.5])
    pred_lib = runner._derive_libration_deg(theta_dict)
    if pred_lib is None:
        pytest.skip('libration forward model unservable on this smoke node')
    ll = runner.log_likelihood_fn(theta)
    assert np.isfinite(ll)

    pred = runner._derive_gravity_isostatic(theta_dict)
    v, s = LIB_OBS
    chi2_lib = ((pred_lib - v) / s) ** 2
    chi2_grav = sum(((pred[lm] - obs[0]) / obs[1]) ** 2
                    for obs, lm in ((C20_OBS, (2, 0)), (C22_OBS, (2, 2)),
                                    (C30_OBS, (3, 0))))
    assert ll == pytest.approx(-0.5 * (chi2_grav + chi2_lib), rel=1e-12, abs=1e-12)


def test_libration_deg_sigma_model_add_inflates():
    sigma_model = dict(SIGMA_MODEL_GRAVITY)
    sigma_model['libration_deg'] = 0.00025
    runner = MCMCRunner(_config(sigma_model_add=sigma_model, include_libration=True))
    theta_dict = {'Tb_K': OCEAN_TB, 'compensation_C2': 0.5}
    theta = np.array([OCEAN_TB, 0.5])
    pred_lib = runner._derive_libration_deg(theta_dict)
    if pred_lib is None:
        pytest.skip('libration forward model unservable on this smoke node')
    ll = runner.log_likelihood_fn(theta)

    pred = runner._derive_gravity_isostatic(theta_dict)
    v, s = LIB_OBS
    s_eff = float(np.sqrt(s ** 2 + 0.00025 ** 2))
    chi2_lib = ((pred_lib - v) / s_eff) ** 2
    chi2_grav = 0.0
    for name, obs, lm in (('C20', C20_OBS, (2, 0)),
                          ('C22', C22_OBS, (2, 2)),
                          ('C30', C30_OBS, (3, 0))):
        v_, s_ = obs
        s_eff_ = float(np.sqrt(s_ ** 2 + SIGMA_MODEL_GRAVITY[name] ** 2))
        chi2_grav += ((pred[lm] - v_) / s_eff_) ** 2
    assert ll == pytest.approx(-0.5 * (chi2_grav + chi2_lib), rel=1e-12, abs=1e-12)


# ===========================================================================
# The actual candidate config's declared sigma_model_add block (post config-
# schema reconciliation) parses and produces finite log-likelihoods.
# ===========================================================================

CANDIDATE_CONFIG = (REPO / 'PlanetProfile' / 'Inference' / 'configs'
                    / 'enceladus_cassini_isostasy_7D.json')


@pytest.mark.skipif(not CANDIDATE_CONFIG.exists(), reason='candidate config not present')
def test_candidate_config_sigma_model_add_reachable_and_wired():
    cfg = InferenceConfig.from_json(str(CANDIDATE_CONFIG))
    sigma_model_add = cfg.metadata.get('sigma_model_add')
    assert sigma_model_add == {
        'C20': 5.3e-06, 'C22': 1.7e-06, 'C30': 4.4e-06, 'libration_deg': 0.00025,
    }
    # And these are exactly the already-ratified (B15 CLOSED) isostasy
    # values, just re-keyed by observable name.
    iso = cfg.metadata['isostasy']
    assert sigma_model_add['C20'] == iso['sigma_model_lm']['2,0']
    assert sigma_model_add['C22'] == iso['sigma_model_lm']['2,2']
    assert sigma_model_add['C30'] == iso['sigma_model_lm']['3,0']
    assert sigma_model_add['libration_deg'] == iso['sigma_model_libration_deg']
