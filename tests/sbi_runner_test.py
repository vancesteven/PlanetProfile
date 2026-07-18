"""
Tests for the SBI (amortized NPE) pipeline: PlanetProfile/Inference/sbi_runner.py
and the SBI-specific opt-in kwargs of MCMCRunner.generate_sbi_dataset.

Phase C of plans/monte-carlo-sbi-implementation-plan.md.

Fixture strategy
-----------------
The Test50 structure-cache pickle referenced by
``configs/test50_titan_noocean_andrade_8D.json``
(``PlanetProfile/Test/mcmc_results/Titan/Test50_andrade_noocean_yao2014/
titan_allice_yao2014_structure_grid.pkl``) does not exist in this checkout, so
tests that need a real ``MCMCRunner`` (dataset-generation tests) build a tiny
synthetic 2-point Tb_K grid cache in a tmp dir and point
``structure_cache_path`` at it. The synthetic grid deliberately gives its two
grid points different ``region_phases`` tuples so ``apply_bottom_temperature``
always takes its nearest-neighbour fallback (never the b-layered blend), which
keeps the fixture's required fields to the small set ``_check_no_ocean`` and
``apply_viscosity_params``/``apply_andrade_params`` actually touch.

The real forward model (``forward_model_k2_flexible``, which calls into
TidalPy's radial_solver) is monkeypatched with a small deterministic stub for
every test that exercises ``generate_sbi_dataset``/``generate_training_set``,
since the fixture grid is not a physically valid PlanetProfile structure.
Everything *around* the forward-model call (prior sampling, parameter
expansion/hooks, the no-ocean support guard, non-finite row dropping, seed
handling, and .npz schema) is exercised for real.

Artifact/training tests (7, 8) test sbi_runner.py's train/save/load/sample
plumbing directly against a cheap synthetic BoxUniform + linear-Gaussian toy
simulator -- no PlanetProfile forward model or structure cache involved.
"""
import json
import pickle
import tempfile
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.mcmc_runner import MCMCRunner

try:
    import torch  # noqa: F401
    import sbi  # noqa: F401
    _HAVE_SBI = True
except ImportError:
    _HAVE_SBI = False

if _HAVE_SBI:
    from PlanetProfile.Inference.sbi_runner import SBIRunner


REPO_ROOT = Path(__file__).resolve().parents[1]
BASE_CONFIG_PATH = (
    REPO_ROOT / 'PlanetProfile' / 'Inference' / 'configs'
    / 'test50_titan_noocean_andrade_8D.json'
)

TB_LOW = 249.0
TB_HIGH = 250.965
TB_MID = 0.5 * (TB_LOW + TB_HIGH)


# ---------------------------------------------------------------------------
# Config / structure-cache fixtures
# ---------------------------------------------------------------------------

def _load_base_config_dict():
    with open(BASE_CONFIG_PATH) as f:
        return json.load(f)


def _make_synthetic_structure(tag, T_Ih):
    """Minimal structure dict covering only the fields exercised by the
    Tb_K/andrade/viscosity parameter hooks and the no-ocean guard."""
    phases = np.array([1, 1, 50, 50, 3, 3], dtype=float)
    P_MPa = np.array([0.1, 5.0, 50.0, 100.0, 150.0, 200.0])
    T_K = np.array([T_Ih[0], T_Ih[1], 245.0, 255.0, 265.0, 270.0])
    return {
        'phases': phases,
        'P_MPa': P_MPa,
        'T_K': T_K,
        'eta_Pa_base': np.full(6, 1e14),
        'changeIndices': np.array([0, 2, 4, 6]),
        'n_layers': 3,
        # Distinct region_phases per grid point -> apply_bottom_temperature's
        # _regions_match() is always False -> nearest-neighbour fallback,
        # never the b-layered blend (which needs r_m and more fields).
        'region_phases': (tag,),
        'CMR2': 0.31,
        'Mtot_kg': 1.3452e23,
        'D_ocean_km': 0.0 if tag == 'lowT' else 60.0,
        'D_iceIh_km': 190.0,
    }


def _write_fake_grid_cache(path):
    """Write a 2-point Tb_K grid cache (Test50 list format) to ``path``.

    Ice-Ih temperatures are chosen so the low-Tb grid point is well below the
    melting curve (no-ocean guard: accept) and the high-Tb grid point is
    above it (no-ocean guard: reject). Tm_Ih_lin(P) = 273.16 - 0.068*P; at
    P=0.1 MPa that's ~273.09 K, at P=5.0 MPa ~272.82 K.
    """
    s_low = _make_synthetic_structure('lowT', T_Ih=(230.0, 235.0))
    s_high = _make_synthetic_structure('highT', T_Ih=(273.05, 273.4))
    grid = {'Tb_K_grid': np.array([TB_LOW, TB_HIGH]), 'structures': [s_low, s_high]}
    with open(path, 'wb') as f:
        pickle.dump(grid, f)


def _build_mcmc_runner(tmp_path, observables=None):
    """Construct a real MCMCRunner against the Test50 8D config, with
    structure_cache_path overridden to a synthetic tmp fixture (the real
    Test50 .pkl is not present in this checkout)."""
    data = _load_base_config_dict()
    cache_path = tmp_path / 'fake_grid.pkl'
    _write_fake_grid_cache(cache_path)
    data['structure_cache_path'] = str(cache_path)
    data['sampler_settings'] = dict(data['sampler_settings'])
    data['sampler_settings']['checkpoint_dir'] = str(tmp_path / 'checkpoints')
    if observables is not None:
        data['observables'] = observables
    config = InferenceConfig.from_dict(data)
    return MCMCRunner(config)


_FORWARD_MODEL_TARGET = 'PlanetProfile.Inference.forward_models.forward_model_k2_flexible'


def _stub_forward_model_basic(theta_dict, structure_data, return_heating=False,
                               arrhenius_params=None):
    """Deterministic, always-finite stand-in for the real (TidalPy-backed)
    forward model. Re_k2/Im_k2 are simple functions of theta only, so runs
    with identical theta are bit-reproducible regardless of structure_data."""
    alpha = float(theta_dict['alpha'])
    log10_eta_Ih = float(theta_dict['log10_eta_Ih'])
    Re_k2 = 0.5 + 0.2 * alpha
    Im_k2 = -(0.05 + 0.01 * (log10_eta_Ih - 10.0))  # always negative
    perPhase = {'Ih': 1.0e9} if return_heating else None
    return Re_k2, Im_k2, np.nan, np.nan, perPhase


def _make_nan_above_alpha_stub(threshold):
    def _stub(theta_dict, structure_data, return_heating=False, arrhenius_params=None):
        if float(theta_dict['alpha']) > threshold:
            return np.nan, np.nan, np.nan, np.nan, None
        return _stub_forward_model_basic(theta_dict, structure_data, return_heating,
                                          arrhenius_params)
    return _stub


def _make_toy_sbi_config(bounds=((-2.0, 2.0), (-1.0, 3.0)), seed=42):
    """Cheap 2D uniform-prior SBI config for the artifact/training tests --
    no PlanetProfile forward model or structure cache needed."""
    param_names = ['theta1', 'theta2']
    param_space = {
        name: {'prior_type': 'uniform', 'bounds': list(b)}
        for name, b in zip(param_names, bounds)
    }
    data = {
        'mode': 'sbi',
        'bodyname': 'ToyBody',
        'param_space': param_space,
        'observables': {'obs_a': [0.0, 1.0], 'obs_b': [0.0, 1.0]},
        'sampler_settings': {
            'num_simulations': 200, 'density_estimator': 'maf',
            'n_posterior_samples': 50, 'n_reeval': 10,
        },
        'structure_cache_path': 'unused_nonexistent.pkl',
        'random_state': seed,
    }
    return InferenceConfig.from_dict(data)


class SBIPriorConsistencyTests(unittest.TestCase):
    """Task 1: SBIRunner.build_prior() bounds vs MCMCRunner._build_prior()."""

    @unittest.skipUnless(_HAVE_SBI, 'sbi/torch not installed')
    def test_sbi_prior_matches_mcmc_bounds(self):
        data = _load_base_config_dict()
        mcmc_config = InferenceConfig.from_dict(data)
        sbi_data = dict(data)
        sbi_data['mode'] = 'sbi'
        sbi_config = InferenceConfig.from_dict(sbi_data)

        # MCMCRunner._build_prior only touches self.config/self.param_names;
        # call it unbound to avoid constructing a full MCMCRunner (and thus
        # needing a structure cache) just for a prior-bounds comparison.
        fake_self = SimpleNamespace(
            config=mcmc_config, param_names=list(mcmc_config.param_space.keys())
        )
        mcmc_prior = MCMCRunner._build_prior(fake_self)
        mcmc_bounds = np.asarray(mcmc_prior.bounds, dtype=float)

        runner = SBIRunner(sbi_config)
        sbi_prior = runner.build_prior()
        sbi_bounds = np.stack(
            [sbi_prior.base_dist.low.numpy(), sbi_prior.base_dist.high.numpy()],
            axis=1,
        )

        # SBIRunner.build_prior() casts bounds to torch.float32 (BoxUniform
        # requirement); MCMCRunner's bounds stay float64. Tolerance accounts
        # for float32 round-trip only, not a substantive bounds mismatch.
        np.testing.assert_allclose(sbi_bounds, mcmc_bounds, rtol=1e-6, atol=1e-4)


class SBIDatasetGenerationTests(unittest.TestCase):
    """Tasks 2-6: MCMCRunner.generate_sbi_dataset opt-in kwargs."""

    def test_dataset_shapes_and_schema(self):
        with tempfile.TemporaryDirectory() as td:
            tmp_path = Path(td)
            runner = _build_mcmc_runner(tmp_path)
            n_params = len(runner.param_names)
            n_obs = len(runner.config.observables)

            with mock.patch(_FORWARD_MODEL_TARGET, side_effect=_stub_forward_model_basic):
                # Default call: exactly the original 4 .npz keys.
                out_default = tmp_path / 'default.npz'
                theta, x = runner.generate_sbi_dataset(
                    n_samples=20, output_path=str(out_default)
                )
                self.assertEqual(theta.shape, (20, n_params))
                self.assertEqual(x.shape, (20, n_obs))
                with np.load(out_default, allow_pickle=True) as npz:
                    self.assertEqual(
                        set(npz.files), {'theta', 'x', 'param_names', 'obs_names'}
                    )

                # Extended call: additive keys present alongside the original 4.
                out_extended = tmp_path / 'extended.npz'
                theta2, x2 = runner.generate_sbi_dataset(
                    n_samples=20, output_path=str(out_extended), seed=7,
                    imag_convention='abs', drop_nonfinite=True,
                    provenance={'note': 'unit-test'},
                )
                self.assertEqual(theta2.shape[1], n_params)
                self.assertEqual(x2.shape[1], n_obs)
                with np.load(out_extended, allow_pickle=True) as npz2:
                    expected_additive = {
                        'theta', 'x', 'param_names', 'obs_names',
                        'param_bounds', 'param_units', 'config_hash', 'git_sha',
                        'seed', 'n_requested', 'n_rejected_nonfinite',
                        'n_rejected_support', 'imag_convention', 'provenance',
                    }
                    self.assertEqual(set(npz2.files), expected_additive)

    def test_imag_convention_abs(self):
        with tempfile.TemporaryDirectory() as td:
            tmp_path = Path(td)
            runner = _build_mcmc_runner(
                tmp_path, observables={'Re_k2': [0.608, 0.048], 'Im_k2': [0.135, 0.035]}
            )
            im_idx = list(runner.config.observables.keys()).index('Im_k2')

            with mock.patch(_FORWARD_MODEL_TARGET, side_effect=_stub_forward_model_basic):
                _, x_signed = runner.generate_sbi_dataset(
                    n_samples=25, seed=17, imag_convention='signed'
                )
                _, x_abs = runner.generate_sbi_dataset(
                    n_samples=25, seed=17, imag_convention='abs'
                )

            self.assertTrue(np.all(x_abs[:, im_idx] >= 0))
            np.testing.assert_allclose(x_abs[:, im_idx], np.abs(x_signed[:, im_idx]))

            with self.assertRaises(ValueError):
                runner.generate_sbi_dataset(n_samples=5, imag_convention='bogus')

            # abs_Im_k2 observable name: 'abs' convention fills it (not NaN);
            # 'signed' convention NaN-fills it (documented legacy behavior).
            runner_alias = _build_mcmc_runner(
                tmp_path, observables={'Re_k2': [0.608, 0.048], 'abs_Im_k2': [0.135, 0.035]}
            )
            alias_idx = list(runner_alias.config.observables.keys()).index('abs_Im_k2')
            with mock.patch(_FORWARD_MODEL_TARGET, side_effect=_stub_forward_model_basic):
                _, x_alias_abs = runner_alias.generate_sbi_dataset(
                    n_samples=10, seed=3, imag_convention='abs'
                )
                _, x_alias_signed = runner_alias.generate_sbi_dataset(
                    n_samples=10, seed=3, imag_convention='signed'
                )
            self.assertTrue(np.all(np.isfinite(x_alias_abs[:, alias_idx])))
            self.assertTrue(np.all(np.isnan(x_alias_signed[:, alias_idx])))

    def test_nonfinite_rows_dropped_and_counted(self):
        with tempfile.TemporaryDirectory() as td:
            tmp_path = Path(td)
            runner = _build_mcmc_runner(tmp_path)
            n_samples = 200
            seed = 123
            threshold = 0.40  # alpha in [0.15, 0.45]

            np.random.seed(seed)
            theta_check = runner.prior.rvs(n_samples)
            alpha_idx = runner.param_names.index('alpha')
            expected_nan = int(np.sum(theta_check[:, alpha_idx] > threshold))
            self.assertGreater(expected_nan, 0)
            self.assertLess(expected_nan, n_samples)

            with mock.patch(
                _FORWARD_MODEL_TARGET,
                side_effect=_make_nan_above_alpha_stub(threshold),
            ):
                theta, x = runner.generate_sbi_dataset(
                    n_samples=n_samples, seed=seed, drop_nonfinite=True,
                    apply_support_guard=False,
                )

            self.assertEqual(theta.shape[0], n_samples - expected_nan)
            self.assertEqual(x.shape[0], n_samples - expected_nan)
            self.assertTrue(np.all(np.isfinite(x)))
            self.assertEqual(
                runner.last_sbi_dataset_stats['n_rejected_nonfinite'], expected_nan
            )
            self.assertEqual(runner.last_sbi_dataset_stats['n_rejected_support'], 0)

    def test_no_ocean_guard_applied(self):
        with tempfile.TemporaryDirectory() as td:
            tmp_path = Path(td)
            runner = _build_mcmc_runner(tmp_path)
            n_samples = 200
            seed = 55

            np.random.seed(seed)
            theta_check = runner.prior.rvs(n_samples)
            tb_idx = runner.param_names.index('Tb_K')
            expected_rejected = int(np.sum(theta_check[:, tb_idx] > TB_MID))
            self.assertGreater(expected_rejected, 0)
            self.assertLess(expected_rejected, n_samples)

            with mock.patch(_FORWARD_MODEL_TARGET, side_effect=_stub_forward_model_basic):
                theta_guard, x_guard = runner.generate_sbi_dataset(
                    n_samples=n_samples, seed=seed, apply_support_guard=True
                )
                stats_guard = dict(runner.last_sbi_dataset_stats)

                theta_noguard, x_noguard = runner.generate_sbi_dataset(
                    n_samples=n_samples, seed=seed, apply_support_guard=False
                )
                stats_noguard = dict(runner.last_sbi_dataset_stats)

            self.assertEqual(stats_guard['n_rejected_support'], expected_rejected)
            self.assertEqual(theta_guard.shape[0], n_samples - expected_rejected)
            self.assertEqual(x_guard.shape[0], n_samples - expected_rejected)

            self.assertEqual(stats_noguard['n_rejected_support'], 0)
            self.assertEqual(theta_noguard.shape[0], n_samples)
            self.assertEqual(x_noguard.shape[0], n_samples)

    def test_seed_reproducibility(self):
        with tempfile.TemporaryDirectory() as td:
            tmp_path = Path(td)
            runner = _build_mcmc_runner(tmp_path)

            with mock.patch(_FORWARD_MODEL_TARGET, side_effect=_stub_forward_model_basic):
                theta_a, x_a = runner.generate_sbi_dataset(n_samples=30, seed=99)
                theta_b, x_b = runner.generate_sbi_dataset(n_samples=30, seed=99)
                theta_c, x_c = runner.generate_sbi_dataset(n_samples=30, seed=100)

            np.testing.assert_array_equal(theta_a, theta_b)
            np.testing.assert_array_equal(x_a, x_b)
            self.assertFalse(np.array_equal(theta_a, theta_c))

    def test_obs_noise_injection(self):
        with tempfile.TemporaryDirectory() as td:
            tmp_path = Path(td)
            runner = _build_mcmc_runner(tmp_path)
            obs_names = list(runner.config.observables.keys())
            sigmas = np.array([float(runner.config.observables[n][1])
                               for n in obs_names])

            with mock.patch(_FORWARD_MODEL_TARGET, side_effect=_stub_forward_model_basic):
                theta_c, x_clean = runner.generate_sbi_dataset(
                    n_samples=200, seed=7, imag_convention='abs')
                theta_n, x_noisy = runner.generate_sbi_dataset(
                    n_samples=200, seed=7, imag_convention='abs',
                    obs_noise=True, noise_seed=11)
                theta_n2, x_noisy2 = runner.generate_sbi_dataset(
                    n_samples=200, seed=7, imag_convention='abs',
                    obs_noise=True, noise_seed=11)

            # Prior draws unaffected by the noise Generator; noise reproducible.
            np.testing.assert_array_equal(theta_c, theta_n)
            np.testing.assert_array_equal(x_noisy, x_noisy2)

            # Residuals are Gaussian with the config sigmas (loose n=200 band).
            resid = x_noisy - x_clean
            self.assertFalse(np.allclose(resid, 0.0))
            for s_emp, s_true in zip(resid.std(axis=0), sigmas):
                self.assertGreater(s_emp, 0.6 * s_true)
                self.assertLess(s_emp, 1.4 * s_true)

            # Provenance recorded; likelihood convention (no re-fold).
            noise_meta = runner.last_sbi_dataset_stats['obs_noise']
            self.assertEqual(noise_meta['noise_seed'], 11)
            self.assertFalse(noise_meta['refold_im'])
            self.assertEqual(noise_meta['sigma'],
                             {n: float(s) for n, s in zip(obs_names, sigmas)})

            # Signed-Im convention is ambiguous under the noise model: rejected.
            with self.assertRaises(ValueError):
                runner.generate_sbi_dataset(n_samples=5, obs_noise=True)

    def test_induction_channels_computed_from_ae_grid_cache(self):
        """Induction observables mirror the likelihood's grid-cache lookup
        (2026-07-10 fix; previously silently NaN-filled in SBI datasets)."""
        obs = {
            'Re_k2': [0.608, 0.048],
            'Im_k2': [0.135, 0.035],
            'Ae_synodic_real': [0.9, 0.05],
            'Ae_synodic_imag': [-0.2, 0.05],
            'BiAmp_synodic': [0.92, 0.05],
            'BiPhase_synodic_deg': [-12.5, 5.0],
        }
        with tempfile.TemporaryDirectory() as td:
            runner = _build_mcmc_runner(Path(td), observables=obs)
            # Fixture structures carry no Texc_hr, so inject the Ae cache the
            # likelihood would have precomputed (one complex Ae per Tb point).
            ae_low, ae_high = 0.90 - 0.20j, 0.80 - 0.30j
            runner._ae_grid_cache = {0: {'synodic': ae_low},
                                     1: {'synodic': ae_high}}

            with mock.patch(_FORWARD_MODEL_TARGET, side_effect=_stub_forward_model_basic):
                theta, x = runner.generate_sbi_dataset(
                    n_samples=40, seed=13, imag_convention='abs',
                    drop_nonfinite=True)

            self.assertEqual(x.shape[1], len(obs))
            names = list(obs.keys())
            tb = theta[:, runner.param_names.index('Tb_K')]
            grid = np.array([TB_LOW, TB_HIGH])
            expected_ae = np.where(
                np.abs(tb - TB_LOW) <= np.abs(tb - TB_HIGH), ae_low, ae_high)
            np.testing.assert_allclose(
                x[:, names.index('Ae_synodic_real')], expected_ae.real, atol=1e-12)
            np.testing.assert_allclose(
                x[:, names.index('Ae_synodic_imag')], expected_ae.imag, atol=1e-12)
            np.testing.assert_allclose(
                x[:, names.index('BiAmp_synodic')], np.abs(expected_ae), atol=1e-12)
            np.testing.assert_allclose(
                x[:, names.index('BiPhase_synodic_deg')],
                np.degrees(np.angle(expected_ae)), atol=1e-12)

            # Unavailable cache -> NaN -> rows rejected under drop_nonfinite.
            runner._ae_grid_cache = {}
            with mock.patch(_FORWARD_MODEL_TARGET, side_effect=_stub_forward_model_basic):
                theta2, x2 = runner.generate_sbi_dataset(
                    n_samples=10, seed=13, imag_convention='abs',
                    drop_nonfinite=True)
            self.assertEqual(len(theta2), 0)
            self.assertEqual(
                runner.last_sbi_dataset_stats['n_rejected_nonfinite'], 10)

    def test_collect_posterior_ae_metadata(self):
        """_collect_posterior_Ae: JSON-safe per-sample complex Ae per label,
        nearest-grid lookup identical to the likelihood's; None without an
        Ae cache. Feeds metadata['induction_Ae'] in BOTH runners."""
        import json as _json
        with tempfile.TemporaryDirectory() as td:
            runner = _build_mcmc_runner(Path(td))
            ae_low, ae_high = 0.90 - 0.20j, 0.80 - 0.30j
            runner._ae_grid_cache = {0: {'synodic': ae_low},
                                     1: {'synodic': ae_high}}
            n = 25
            rng = np.random.RandomState(3)
            samples = np.column_stack([
                rng.uniform(*runner.config.param_space[p]['bounds'], size=n)
                for p in runner.param_names])
            out = runner._collect_posterior_Ae(samples)
            self.assertEqual(set(out), {'synodic'})
            self.assertEqual(len(out['synodic']['re']), n)
            tb = samples[:, runner.param_names.index('Tb_K')]
            expected = np.where(
                np.abs(tb - TB_LOW) <= np.abs(tb - TB_HIGH), ae_low, ae_high)
            np.testing.assert_allclose(out['synodic']['re'], expected.real,
                                       atol=1e-12)
            np.testing.assert_allclose(out['synodic']['im'], expected.imag,
                                       atol=1e-12)
            # JSON-safe (saved-run round trip)
            _json.dumps(out)

            # Be metadata: None without Bind_ observables; JSON-safe when set.
            self.assertIsNone(runner._be_excitation_metadata())
            runner._be_excitation = {'synodic': {'x': 128.5 - 170.1j}}
            be = runner._be_excitation_metadata()
            self.assertEqual(be['synodic']['x'], [128.5, -170.1])
            _json.dumps(be)

            # No cache -> None
            runner._ae_grid_cache = {}
            self.assertIsNone(runner._collect_posterior_Ae(samples))

    def test_induction_bounds_support_rejection(self):
        """induction_bounds (ratified 2026-07-12): one-sided Ae support cuts
        reject dataset rows under apply_support_guard and hard-reject in the
        likelihood; both use the same Tb-grid Ae lookup."""
        obs = {'Re_k2': [0.608, 0.048], 'Im_k2': [0.135, 0.035]}
        with tempfile.TemporaryDirectory() as td:
            runner = _build_mcmc_runner(Path(td), observables=obs)
            # Isolate the induction bound: the fixture's high-Tb grid point
            # is built to violate the no-ocean guard, which would otherwise
            # reject the very rows this test needs to survive (Europa-like
            # ocean-bearing configs carry no phase_stability block anyway).
            runner.config.sampler_settings = dict(runner.config.sampler_settings)
            runner.config.sampler_settings.pop('phase_stability', None)
            runner.config.induction_bounds = {
                'synodic': {'amp_min': 0.7, 'im_abs_max': 0.4}}
            # Low-Tb grid point violates amp_min; high-Tb point passes both.
            ae_low, ae_high = (0.10 - 0.05j), (0.90 - 0.05j)
            runner._ae_grid_cache = {0: {'synodic': ae_low},
                                     1: {'synodic': ae_high}}

            with mock.patch(_FORWARD_MODEL_TARGET, side_effect=_stub_forward_model_basic):
                # Rebuild the likelihood INSIDE the patch: the builder
                # imports forward_model_k2_flexible at build time, so the
                # closure must capture the stub.
                runner.log_likelihood_fn = runner._make_flexible_log_likelihood(
                    runner.config.observables, runner.structure_data,
                    arrhenius_params=runner.arrhenius_params)
                theta, x = runner.generate_sbi_dataset(
                    n_samples=60, seed=21, imag_convention='abs',
                    apply_support_guard=True, drop_nonfinite=True)
            # Every surviving row must sit at the passing (high-Tb) grid point.
            tb = theta[:, runner.param_names.index('Tb_K')]
            self.assertGreater(len(theta), 0)
            self.assertTrue(np.all(np.abs(tb - TB_HIGH) < np.abs(tb - TB_LOW)))
            self.assertGreater(
                runner.last_sbi_dataset_stats['n_rejected_support'], 0)

            # Likelihood: hard reject at the violating grid point, finite at
            # the passing one. theta order: [alpha, log10_zeta, eta_Ih,
            # eta_III, eta_V, eta_VI, eta_sil, Tb_K].
            base = [0.3, -1.0, 12.0, 12.0, 12.0, 12.0, 20.0, TB_LOW]
            with mock.patch(_FORWARD_MODEL_TARGET, side_effect=_stub_forward_model_basic):
                ll_low = runner.log_likelihood_fn(np.array(base))
                base[-1] = TB_HIGH
                ll_high = runner.log_likelihood_fn(np.array(base))
            self.assertEqual(ll_low, -1e30)
            self.assertTrue(np.isfinite(ll_high) and ll_high > -1e29)

            # im_abs_max bites independently of amp_min.
            runner._ae_grid_cache = {0: {'synodic': ae_low},
                                     1: {'synodic': (0.9 - 0.6j)}}
            with mock.patch(_FORWARD_MODEL_TARGET, side_effect=_stub_forward_model_basic):
                runner.log_likelihood_fn = runner._make_flexible_log_likelihood(
                    runner.config.observables, runner.structure_data,
                    arrhenius_params=runner.arrhenius_params)
                ll_high_im = runner.log_likelihood_fn(np.array(base))
            self.assertEqual(ll_high_im, -1e30)


@unittest.skipUnless(_HAVE_SBI, 'sbi/torch not installed')
class InferFromArtifactTests(unittest.TestCase):
    """infer_from_artifact: compatibility validation + truncation bounds
    (packaging path is covered by the real-artifact smoke below and the GUI
    verification; these tests need no forward model)."""

    @classmethod
    def setUpClass(cls):
        cls.config = _make_toy_sbi_config()
        runner = SBIRunner(cls.config)
        theta, x = _toy_train_pair(runner, n=300, seed=0)
        runner.train(theta, x, seed=0, density_estimator='maf', max_num_epochs=40)
        cls.td = tempfile.TemporaryDirectory()
        cls.artifact_path = Path(cls.td.name) / 'toy_artifact.pt'
        runner.save_artifact(cls.artifact_path)

    @classmethod
    def tearDownClass(cls):
        cls.td.cleanup()

    def test_rejects_mismatched_prior_box(self):
        cfg = _make_toy_sbi_config(bounds=((-1.0, 1.0), (-1.0, 3.0)))
        runner = SBIRunner(cfg)
        with self.assertRaisesRegex(ValueError, 'prior box'):
            runner.infer_from_artifact(self.artifact_path)

    def test_rejects_mismatched_params(self):
        cfg = _make_toy_sbi_config()
        d = cfg.to_dict()
        d['param_space'] = {'other1': d['param_space']['theta1'],
                            'other2': d['param_space']['theta2']}
        from PlanetProfile.Inference.inference_core import InferenceConfig
        runner = SBIRunner(InferenceConfig.from_dict(d))
        with self.assertRaisesRegex(ValueError, 'parameters'):
            runner.infer_from_artifact(self.artifact_path)

    def test_truncation_bounds_validated(self):
        runner = SBIRunner(self.config)
        artifact = runner._load_artifact_dict(self.artifact_path)
        runner._posterior = artifact['posterior']
        with self.assertRaisesRegex(ValueError, 'inside the trained prior box'):
            runner._validate_truncation({'theta1': (-5.0, 0.0)})
        with self.assertRaisesRegex(ValueError, 'unknown parameter'):
            runner._validate_truncation({'nope': (0.0, 1.0)})
        # In-box truncation validates cleanly.
        runner._validate_truncation({'theta1': (-1.0, 1.0)})


def _toy_train_pair(runner, n, seed):
    prior = runner.build_prior()
    torch.manual_seed(seed)
    theta_t = prior.sample((n,))
    x_t = theta_t + 0.05 * torch.randn_like(theta_t)
    return theta_t.numpy().astype(np.float64), x_t.numpy().astype(np.float64)


@unittest.skipUnless(_HAVE_SBI, 'sbi/torch not installed')
class SBIArtifactAndTrainingTests(unittest.TestCase):
    """Tasks 7-8: sbi_runner.py train/save/load/sample plumbing against a
    cheap synthetic BoxUniform + linear-Gaussian toy simulator."""

    def _toy_training_data(self, runner, n, seed):
        prior = runner.build_prior()
        torch.manual_seed(seed)
        theta_t = prior.sample((n,))
        x_t = theta_t + 0.05 * torch.randn_like(theta_t)
        return (
            theta_t.numpy().astype(np.float64),
            x_t.numpy().astype(np.float64),
        )

    def test_artifact_roundtrip_metadata(self):
        config = _make_toy_sbi_config()
        runner = SBIRunner(config)
        theta, x = self._toy_training_data(runner, n=300, seed=0)

        runner.train(theta, x, seed=0, density_estimator='maf', max_num_epochs=40)

        with tempfile.TemporaryDirectory() as td:
            artifact_path = Path(td) / 'toy_artifact.pt'
            runner.save_artifact(artifact_path)

            loaded = SBIRunner.load_artifact(artifact_path)
            self.assertEqual(loaded.param_names, runner.param_names)
            self.assertEqual(loaded.imag_convention, 'abs')
            self.assertEqual(loaded.artifact_meta['config_hash'], config.generate_hash())
            self.assertEqual(loaded.artifact_meta['prior_spec'], runner._prior_spec())

            lows, highs = runner._param_bounds()
            np.testing.assert_allclose(
                loaded.artifact_meta['param_bounds'],
                [[lo, hi] for lo, hi in zip(lows, highs)],
            )

            samples = loaded.sample_posterior(
                {'obs_a': 0.2, 'obs_b': -0.3}, n_samples=64, seed=1
            )
            self.assertEqual(samples.shape, (64, 2))
            self.assertTrue(np.all(np.isfinite(samples)))
            for d, (lo, hi) in enumerate(zip(lows, highs)):
                self.assertTrue(np.all(samples[:, d] >= lo - 1e-4))
                self.assertTrue(np.all(samples[:, d] <= hi + 1e-4))

            # Tamper sbi_version -> loud RuntimeWarning on load (schema v1
            # pickles a full DirectPosterior; version mismatches must not be
            # silent per sbi_runner.py's _load_artifact_dict).
            raw = torch.load(artifact_path, map_location='cpu', weights_only=False)
            raw['sbi_version'] = '0.0.0-tampered'
            tampered_path = Path(td) / 'tampered.pt'
            torch.save(raw, tampered_path)
            with self.assertWarns(RuntimeWarning):
                SBIRunner.load_artifact(tampered_path)

            # Local build suffixes are NOT a version mismatch: an artifact
            # saved with '2.8.0' loading under '2.8.0+cpu' (cloud CPU wheel)
            # must neither warn nor defeat validated_version_pairs.
            raw2 = torch.load(artifact_path, map_location='cpu', weights_only=False)
            raw2['torch_version'] = torch.__version__.split('+')[0] + '+fakebuild'
            suffixed_path = Path(td) / 'suffixed.pt'
            torch.save(raw2, suffixed_path)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                SBIRunner.load_artifact(suffixed_path)
            self.assertFalse(
                any('misbehave' in str(x.message) for x in w),
                'build-suffix-only difference must not raise the version warning')

            # A genuinely different release still warns even with a suffix.
            raw3 = torch.load(artifact_path, map_location='cpu', weights_only=False)
            raw3['torch_version'] = '0.0.1+cpu'
            diff_path = Path(td) / 'diffver.pt'
            torch.save(raw3, diff_path)
            with self.assertWarns(RuntimeWarning):
                SBIRunner.load_artifact(diff_path)

    def test_tiny_train_and_sample_smoke(self):
        config = _make_toy_sbi_config()
        runner = SBIRunner(config)
        theta, x = self._toy_training_data(runner, n=200, seed=7)

        posterior = runner.train(theta, x, seed=7, density_estimator='maf',
                                  max_num_epochs=40)
        self.assertIsNotNone(posterior)

        samples = runner.sample_posterior(
            {'obs_a': 0.5, 'obs_b': 0.5}, n_samples=128, seed=2
        )
        self.assertEqual(samples.shape, (128, 2))
        self.assertTrue(np.all(np.isfinite(samples)))

        lows, highs = runner._param_bounds()
        for d, (lo, hi) in enumerate(zip(lows, highs)):
            self.assertTrue(np.all(samples[:, d] >= lo - 1e-4))
            self.assertTrue(np.all(samples[:, d] <= hi + 1e-4))


def _make_spline_assertion(marker_ok=True):
    """Return an AssertionError whose traceback frame filename does (or does
    not) contain the nflows rational_quadratic marker, so the guard's
    frame-filename detector can be exercised deterministically without a real
    nsf flow (whose float-roundoff failure is platform dependent)."""
    if marker_ok:
        # Fabricate a traceback frame whose filename carries the marker the
        # guard matches on (the real assert lives in .../splines/
        # rational_quadratic.py and has no message text).
        code = compile(
            "def f():\n assert False\n",
            '/fake/site-packages/nflows/transforms/splines/'
            'rational_quadratic.py', 'exec')
        ns = {}
        exec(code, ns)
        try:
            ns['f']()
        except AssertionError as e:
            return e
    else:
        try:
            assert False, 'some unrelated assertion'
        except AssertionError as e:
            return e


@unittest.skipUnless(_HAVE_SBI, 'sbi/torch not installed')
class SBISplineResampleGuardTests(unittest.TestCase):
    """The nflows spline-inversion resample guard in sample_posterior.

    The underlying `assert (discriminant >= 0).all()` fires only on a
    floating-point knife-edge that is platform/BLAS dependent (it does not
    reproduce naturally on every machine), so these tests drive the guard's
    control flow deterministically by patching the posterior's .sample to
    raise the marker assertion on the first N calls and succeed afterwards.
    """

    def _trained_toy_runner(self):
        config = _make_toy_sbi_config()
        runner = SBIRunner(config)
        theta, x = _toy_train_pair(runner, n=200, seed=7)
        runner.train(theta, x, seed=7, density_estimator='maf',
                     max_num_epochs=40)
        return runner

    def test_detector_matches_only_spline_assert(self):
        from PlanetProfile.Inference.sbi_runner import (
            _is_nflows_spline_assertion)
        self.assertTrue(
            _is_nflows_spline_assertion(_make_spline_assertion(marker_ok=True)))
        self.assertFalse(
            _is_nflows_spline_assertion(_make_spline_assertion(marker_ok=False)))
        self.assertFalse(_is_nflows_spline_assertion(ValueError('nope')))

    def test_resamples_then_succeeds(self):
        runner = self._trained_toy_runner()
        real_sample = runner._posterior.sample
        calls = {'n': 0}

        def flaky(*args, **kwargs):
            calls['n'] += 1
            if calls['n'] <= 2:  # fail twice, then succeed
                raise _make_spline_assertion(marker_ok=True)
            return real_sample(*args, **kwargs)

        with mock.patch.object(runner._posterior, 'sample', side_effect=flaky):
            samples = runner.sample_posterior(
                {'obs_a': 0.5, 'obs_b': 0.5}, n_samples=64, seed=3)
        self.assertEqual(calls['n'], 3)  # 2 failed + 1 success
        self.assertEqual(samples.shape, (64, 2))
        self.assertTrue(np.all(np.isfinite(samples)))

    def test_persistent_failure_raises_runtimeerror(self):
        runner = self._trained_toy_runner()

        def always_fail(*args, **kwargs):
            raise _make_spline_assertion(marker_ok=True)

        with mock.patch.object(runner._posterior, 'sample',
                               side_effect=always_fail):
            with self.assertRaisesRegex(RuntimeError, 'resample attempts'):
                runner.sample_posterior(
                    {'obs_a': 0.5, 'obs_b': 0.5}, n_samples=64, seed=3)

    def test_non_spline_assertion_propagates(self):
        runner = self._trained_toy_runner()

        def wrong_assert(*args, **kwargs):
            raise _make_spline_assertion(marker_ok=False)

        with mock.patch.object(runner._posterior, 'sample',
                               side_effect=wrong_assert):
            with self.assertRaises(AssertionError):
                runner.sample_posterior(
                    {'obs_a': 0.5, 'obs_b': 0.5}, n_samples=64, seed=3)

    def test_diagnostic_path_not_guarded(self):
        """reject_outside_prior=False is the diagnostic sweep path; the guard
        must NOT resample there — the assertion must propagate on the first
        hit so validate_sbi's containment check measures leakage honestly."""
        runner = self._trained_toy_runner()
        calls = {'n': 0}

        def fail_once(*args, **kwargs):
            calls['n'] += 1
            raise _make_spline_assertion(marker_ok=True)

        with mock.patch.object(runner._posterior, 'sample',
                               side_effect=fail_once):
            with self.assertRaises(AssertionError):
                runner.sample_posterior(
                    {'obs_a': 0.5, 'obs_b': 0.5}, n_samples=64, seed=3,
                    reject_outside_prior=False)
        self.assertEqual(calls['n'], 1)  # no resample in the diagnostic path

    def test_guard_preserves_valid_draw_statistics(self):
        """A run that never trips the assert must be byte-identical to the
        pre-guard behaviour: same seed -> same draws (the guard's reseed only
        fires on attempt>0, which a clean run never reaches)."""
        runner = self._trained_toy_runner()
        a = runner.sample_posterior(
            {'obs_a': 0.5, 'obs_b': 0.5}, n_samples=128, seed=5)
        b = runner.sample_posterior(
            {'obs_a': 0.5, 'obs_b': 0.5}, n_samples=128, seed=5)
        np.testing.assert_array_equal(a, b)


if __name__ == '__main__':
    unittest.main()
