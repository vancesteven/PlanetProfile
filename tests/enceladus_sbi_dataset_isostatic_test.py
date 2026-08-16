"""D2 -- generate_sbi_dataset can produce the observable vector under
``gravity_forward_model='isostatic_hm2019'``.

r5 BUILD BLOCKER D2 (validation_reports/enceladus_isostasy/
r5_ADJUDICATION.md): "generate_sbi_dataset cannot produce the observable
vector under isostatic_hm2019 (gravity_pair None -> C20/C22 fall through;
NO C30 arm; drop_nonfinite rejects everything -> empty dataset). C30 was
wired for the likelihood only."

Three separate holes, all closed by the one isostatic arm:

  * ``gravity_pair`` was only computed under ``_gravity_clairaut_active()``,
    so under the isostatic model it stayed ``None`` and the ``C20``/``C22``
    arms' ``and gravity_pair is not None`` guards failed;
  * there was no ``C30`` arm at all, in either forward model;
  * so every row hit the ``else: xi.append(np.nan)`` default and
    ``drop_nonfinite=True`` -- the setting an SBI training run uses --
    rejected 100% of them.

Uses the committed real Enceladus smoke structure cache (3 Tb nodes:
2 frozen, 1 ocean), the same fixture tests/sigma_model_inflation_test.py
uses.
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

C20_OBS = (-5477.45e-6, 36.99e-6)
C22_OBS = (1517.90e-6, 14.70e-6)
C30_OBS = (177.82e-6, 33.42e-6)
LIB_OBS = (0.092, 0.003)

OBS_NAMES = ['C20', 'C22', 'C30', 'libration_deg']

pytestmark = pytest.mark.skipif(
    not SMOKE_CACHE.exists(),
    reason='Enceladus smoke structure cache not present')


def _config():
    """The campaign's observable vector and forward model, on the smoke
    cache's own Tb axis (the production zb axis needs the 87x40 build)."""
    return InferenceConfig(
        mode='mcmc', bodyname='Enceladus',
        param_space={
            'Tb_K': {'prior_type': 'uniform', 'bounds': [271.8, 272.46]},
            'compensation_C2': {'prior_type': 'uniform', 'bounds': [0.0, 1.0]},
            'rho_ice_kgm3': {'prior_type': 'uniform',
                             'bounds': [915.0, 935.0]},
            'libration_sys_frac': {'prior_type': 'truncated_gaussian',
                                   'mean': 0.0, 'std': 0.004,
                                   'bounds': [-0.012, 0.012]},
        },
        observables={'C20': C20_OBS, 'C22': C22_OBS, 'C30': C30_OBS,
                     'libration_deg': LIB_OBS},
        sampler_settings={},
        structure_cache_path=str(SMOKE_CACHE),
        gravity_forward_model='isostatic_hm2019',
        planet_template_module='PlanetProfile.Default.Enceladus.PPEnceladus',
        metadata={
            'isostasy': {
                'H_obs_lm_m': {'2,0': -3510.0, '2,2': 857.0, '3,0': 420.0},
                'shape_ref_radius_m': 252220.0,
                'gravity_ref_radius_m': 256600.0,
                'finite_amplitude': True,
            },
            'libration_model_correction': {
                'multiplicative_factor': 1.008,
                'figure_convention': 'drho_consistent_eq12',
                'H22_obs_m': 857.0,
                'rigid': True,
            },
        },
    )


@pytest.fixture(scope='module')
def mini_dataset():
    runner = MCMCRunner(_config())
    theta, x = runner.generate_sbi_dataset(
        n_samples=50, drop_nonfinite=True, seed=75)
    return theta, x


def test_mini_generation_is_not_empty(mini_dataset):
    """THE D2 PIN: before the isostatic arm this returned zero rows."""
    theta, x = mini_dataset
    assert len(x) > 0, (
        'empty dataset -- drop_nonfinite rejected every row, which is '
        'exactly r5 blocker D2')
    assert 0.0 < len(x) / 50.0 <= 1.0


def test_every_kept_row_is_a_finite_four_channel_vector(mini_dataset):
    theta, x = mini_dataset
    assert x.shape[1] == len(OBS_NAMES)
    assert np.all(np.isfinite(x)), 'a kept row carries a non-finite channel'
    assert theta.shape[0] == x.shape[0]


def test_all_four_channels_actually_vary(mini_dataset):
    """A column that is constant across the kept rows would mean the arm is
    emitting a placeholder rather than the forward model."""
    theta, x = mini_dataset
    for j, name in enumerate(OBS_NAMES):
        assert np.ptp(x[:, j]) > 0.0, f'{name} column is constant'


def test_channels_land_in_physically_sane_ranges(mini_dataset):
    """Order-of-magnitude sanity against the Park 2024 observables -- this
    catches an arm wired to the wrong (l, m) key, which a finiteness check
    alone would pass."""
    theta, x = mini_dataset
    c20, c22, c30, lib = (x[:, j] for j in range(4))
    assert np.all(c20 < 0) and np.all(np.abs(c20) < 5e-2)
    assert np.all(c22 > 0) and np.all(c22 < 5e-2)
    assert np.all(np.abs(c30) < 5e-2)
    assert np.all(lib > 0) and np.all(lib < 1.0)


def test_c30_has_its_own_arm_and_is_not_a_copy_of_c20_or_c22(mini_dataset):
    """C30 was 'wired for the likelihood only'. It must now come from the
    degree-3 key, distinct from both degree-2 channels."""
    theta, x = mini_dataset
    c20, c22, c30 = x[:, 0], x[:, 1], x[:, 2]
    assert not np.allclose(c30, c20)
    assert not np.allclose(c30, c22)


def test_dataset_rows_match_the_likelihood_forward_model(mini_dataset):
    """The training generative model and the MCMC target must be the same
    forward model -- re-derive one kept row's channels through the
    likelihood's own dispatch methods."""
    theta, x = mini_dataset
    runner = MCMCRunner(_config())
    names = list(runner.param_names)
    row = {n: float(v) for n, v in zip(names, theta[0])}
    iso = runner._derive_gravity_isostatic(row)
    assert iso is not None
    assert x[0, 0] == pytest.approx(iso[(2, 0)], rel=1e-12)
    assert x[0, 1] == pytest.approx(iso[(2, 2)], rel=1e-12)
    assert x[0, 2] == pytest.approx(iso[(3, 0)], rel=1e-12)
    assert x[0, 3] == pytest.approx(runner._derive_libration_deg(row),
                                    rel=1e-12)


def test_clairaut_configs_are_untouched():
    """The isostatic arm must not intercept a Clairaut config's C20/C22."""
    cfg = _config()
    cfg.gravity_forward_model = 'clairaut_hydrostatic'
    runner = MCMCRunner(cfg)
    assert not runner._gravity_isostatic_active()
    assert runner._gravity_clairaut_active()
