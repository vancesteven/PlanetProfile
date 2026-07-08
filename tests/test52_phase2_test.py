"""
Regression tests for Test52 Phase 2 (plans/test52-differentiated-titan-plan.md).

Phase 2 adds two things to
``MCMCRunner._derive_cmr2_via_mass_conservation`` (the single shared CMR2
mass-conservation derivation used by the likelihood, ``run()``'s
``cmr2_results``, and ``generate_sbi_dataset``'s ``CMR2`` column):

1. **CMR2 discretization-offset anchor** (plan decision 1): if a sidecar
   JSON exists next to ``config.structure_cache_path`` (named
   ``<cache_stem>_offsets.json``, produced by the Phase 1 cache build), its
   per-Tb ``CMR2_offset`` is linearly interpolated in ``Tb_K`` and ADDED to
   the reconstructed CMR2 before it is returned. No sidecar (every
   pre-Test52 config, including all current Callisto configs) -> no
   offset, numerics unchanged.
2. **Density-inversion guard** (plan decision 3): if
   ``derived_params.rho_sil_kgm3.density_inversion_guard`` is truthy, a
   sample whose derived ``rho_sil`` exceeds the sampled ``rho_core_kgm3``
   is hard-rejected (``None``) exactly like the existing
   ``reject_if_outside_bounds`` check. Absent/false for every config that
   predates this key (including Callisto) -> no-op there.

Physical note (verified numerically while writing these tests, worth
disclosing): under configs/test52_titan_noocean_andrade_10D.json's actual
prior (``rho_core_kgm3`` in [2591, 3600]), the derived silicate density
never exceeds ~2554 kg/m^3 (its no-core-limit maximum, reached as
``R_core_km -> 0``) across the whole Tb grid and R_core domain -- i.e.
strictly below the 2591 kg/m^3 floor. The guard is therefore DORMANT under
Test52's own prior as currently configured (by construction: 2591 kg/m^3
was chosen as the Petricca et al. 2025 bulk-rock density specifically as a
floor above the achievable silicate density). The inversion-guard tests
below exercise the guard mechanism directly with an out-of-prior
``rho_core_kgm3`` value (2500 kg/m^3, i.e. calling
``_derive_cmr2_via_mass_conservation``/``generate_sbi_dataset`` with theta
values a real MCMC run would never draw) to prove the code path is wired
correctly -- this is a defense-in-depth / future-prior-widening guard, not
one that fires during ordinary Test52 sampling today.
"""
import copy
import json
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.mcmc_runner import MCMCRunner
from PlanetProfile.Inference.structure_derivation import (
    compute_cmr2,
    derive_silicate_density,
    extract_hydrosphere_layers,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
TEST52_CONFIG_PATH = (
    REPO_ROOT / 'PlanetProfile' / 'Inference' / 'configs'
    / 'test52_titan_noocean_andrade_10D.json'
)
CALLISTO_CONFIG_PATH = (
    REPO_ROOT / 'PlanetProfile' / 'Inference' / 'configs'
    / 'callisto_nacl_andrade_8D.json'
)


def _load_config_dict(path):
    with open(path) as f:
        return json.load(f)


def _test52_cache_path():
    data = _load_config_dict(TEST52_CONFIG_PATH)
    return REPO_ROOT / data['structure_cache_path']


def _test52_sidecar_path():
    cache_path = _test52_cache_path()
    return cache_path.with_name(cache_path.stem + '_offsets.json')


def _test52_cache_exists():
    return _test52_cache_path().exists() and _test52_sidecar_path().exists()


def _callisto_cache_path():
    data = _load_config_dict(CALLISTO_CONFIG_PATH)
    return REPO_ROOT / data['structure_cache_path']


def _callisto_cache_exists():
    return _callisto_cache_path().exists()


def _reference_cmr2_pre_phase2(runner, theta_dict):
    """CMR2 mass-conservation derivation exactly as
    ``MCMCRunner._derive_cmr2_via_mass_conservation`` read (verbatim) BEFORE
    the Test52 Phase 2 edit -- no offset-anchor addition, no
    density-inversion guard. Used only to prove, for configs where both
    Phase 2 additions are inert by construction (no sidecar file;
    ``density_inversion_guard`` absent/false), that the current method's
    numeric output has not changed.
    """
    derived_params_cfg = getattr(runner.config, 'derived_params', {}) or {}
    rho_sil_cfg = derived_params_cfg.get('rho_sil_kgm3', {}) or {}
    rho_sil_bounds = tuple(rho_sil_cfg.get('bounds', (2200.0, 3500.0)))
    rho_sil_reject = bool(rho_sil_cfg.get('reject_if_outside_bounds', True))

    structure_data = runner.structure_data
    Tb_sample = theta_dict.get('Tb_K')
    if (Tb_sample is not None
            and 'Tb_K_grid' in structure_data
            and 'structures' in structure_data):
        grid_Tb_values = np.asarray(structure_data['Tb_K_grid'])
        idx = int(np.argmin(np.abs(grid_Tb_values - Tb_sample)))
        struct_for_cmr2 = structure_data['structures'][idx]
    elif 'grid_cache' in structure_data and Tb_sample is not None:
        grid_Tb_values = structure_data['grid_Tb_values']
        idx = int(np.argmin(np.abs(grid_Tb_values - Tb_sample)))
        struct_for_cmr2 = structure_data['grid_cache'][grid_Tb_values[idx]]
    else:
        struct_for_cmr2 = structure_data

    R_body_m = struct_for_cmr2.get('R_body_m', np.nan)
    M_total_kg = struct_for_cmr2.get('Mtot_kg', np.nan)
    R_core_km = theta_dict.get('R_core_km')
    rho_core_kgm3 = theta_dict.get('rho_core_kgm3')
    if not (np.isfinite(R_body_m) and np.isfinite(M_total_kg)
            and R_core_km is not None and rho_core_kgm3 is not None):
        return None
    R_core_m = float(R_core_km) * 1000.0
    try:
        hydro_layers, R_oceanbot_m, M_hydro_kg = (
            extract_hydrosphere_layers(struct_for_cmr2)
        )
    except (KeyError, ValueError):
        return None
    if not hydro_layers:
        return None
    rho_sil, accepted = derive_silicate_density(
        M_total_kg=float(M_total_kg),
        M_hydrosphere_kg=M_hydro_kg,
        R_oceanbot_m=R_oceanbot_m,
        R_core_m=R_core_m,
        rho_core_kgm3=float(rho_core_kgm3),
        bounds=rho_sil_bounds,
    )
    if rho_sil_reject and not accepted:
        return None
    assembled = []
    if R_core_m > 0.0:
        assembled.append((0.0, R_core_m, float(rho_core_kgm3)))
    assembled.append((R_core_m, R_oceanbot_m, rho_sil))
    assembled.extend(hydro_layers)
    try:
        cmr2_val = compute_cmr2(
            assembled, R_body_m=float(R_body_m), M_body_kg=float(M_total_kg),
        )
    except ValueError:
        return None
    return cmr2_val


class _RNGIsolatedTestCase(unittest.TestCase):
    """Snapshot/restore the global numpy RNG state around every test.

    Several tests below call ``np.random.seed(...)`` directly, or invoke
    ``generate_sbi_dataset(seed=...)`` (which reseeds the global numpy RNG
    internally per its docstring -- pocoMC's installed ``Prior.rvs()`` has
    no seed argument of its own). Without this isolation, those calls would
    permanently mutate the process-global RNG state and leak into whichever
    test module happens to run next (e.g. ``tests/cmr2_reporting_test.py``'s
    ``test_reported_cmr2_varies_with_core_params``, which draws from the
    prior without seeding its own draws and is therefore sensitive to
    leftover global state) when run together via
    ``python -m unittest test52_phase2_test cmr2_reporting_test ...``.
    """

    def setUp(self):
        self._np_random_state = np.random.get_state()

    def tearDown(self):
        np.random.set_state(self._np_random_state)


def _base_test52_theta(Tb_K=249.0, R_core_km=0.0, rho_core_kgm3=2591.0):
    return {
        'alpha': 0.25, 'log10_zeta': 0.0,
        'log10_eta_Ih': 13.0, 'log10_eta_III': 13.0, 'log10_eta_V': 13.0,
        'log10_eta_VI': 13.0, 'log10_eta_sil': 20.0,
        'Tb_K': Tb_K, 'R_core_km': R_core_km, 'rho_core_kgm3': rho_core_kgm3,
    }


@unittest.skipUnless(
    _test52_cache_exists(),
    'Test52 structure cache/sidecar not present in this checkout',
)
class OffsetAnchorTests(unittest.TestCase):
    """CMR2 discretization-offset anchor (plan decision 1)."""

    @classmethod
    def setUpClass(cls):
        config = InferenceConfig.from_dict(_load_config_dict(TEST52_CONFIG_PATH))
        cls.runner = MCMCRunner(config)
        with open(_test52_sidecar_path()) as f:
            cls.sidecar_json = json.load(f)

    def test_sidecar_found_and_applied(self):
        sidecar = self.runner._load_cmr2_offset_sidecar()
        self.assertIsNotNone(sidecar)
        np.testing.assert_allclose(
            sidecar['Tb_K_grid'], self.sidecar_json['Tb_K_grid']
        )
        np.testing.assert_allclose(
            sidecar['offsets'], self.sidecar_json['CMR2_offset_per_point']
        )

    def test_no_core_cmr2_matches_pp_native_anchor(self):
        """No-core draw's reported CMR2 should land in the PP-native
        0.34297-0.34300 band (reconstructed ~0.342 + offset ~0.00095),
        per the plan's physical sanity anchor (no-core CMR2 = 0.34299)."""
        theta_dict = _base_test52_theta(Tb_K=249.0)
        val = self.runner._compute_model_cmr2(theta_dict)
        self.assertTrue(np.isfinite(val))
        self.assertGreaterEqual(val, 0.34290)
        self.assertLessEqual(val, 0.34305)

    def test_offset_added_matches_sidecar_at_grid_points(self):
        """At exact grid Tb values, (post-edit CMR2) - (reference,
        no-offset CMR2) must equal the sidecar's recorded offset to within
        1e-6."""
        runner = self.runner
        grid_Tb = np.asarray(self.sidecar_json['Tb_K_grid'])
        offsets = np.asarray(self.sidecar_json['CMR2_offset_per_point'])
        for Tb, expected_offset in zip(grid_Tb, offsets):
            theta_dict = _base_test52_theta(Tb_K=float(Tb))
            post = runner._derive_cmr2_via_mass_conservation(theta_dict)
            ref = _reference_cmr2_pre_phase2(runner, theta_dict)
            self.assertIsNotNone(post)
            self.assertIsNotNone(ref)
            self.assertAlmostEqual(post - ref, expected_offset, delta=1e-6)

    def test_offset_interpolation_is_linear_in_tb(self):
        """Off-grid Tb: the applied offset must equal linear interpolation
        of the sidecar's per-point offsets in Tb_K."""
        runner = self.runner
        grid_Tb = np.asarray(self.sidecar_json['Tb_K_grid'])
        offsets = np.asarray(self.sidecar_json['CMR2_offset_per_point'])
        Tb_mid = float(0.5 * (grid_Tb[0] + grid_Tb[1]))
        expected_offset = float(np.interp(Tb_mid, grid_Tb, offsets))

        theta_dict = _base_test52_theta(Tb_K=Tb_mid)
        post = runner._derive_cmr2_via_mass_conservation(theta_dict)
        ref = _reference_cmr2_pre_phase2(runner, theta_dict)
        self.assertAlmostEqual(post - ref, expected_offset, delta=1e-6)


class OffsetPinnedTests(unittest.TestCase):
    """Regression pin on the Phase 1 sidecar's recorded offset values
    themselves (independent of any live model run)."""

    @unittest.skipUnless(
        _test52_sidecar_path().exists(),
        'Test52 offset sidecar not present in this checkout',
    )
    def test_offset_mean_pinned(self):
        with open(_test52_sidecar_path()) as f:
            data = json.load(f)
        self.assertAlmostEqual(data['offset_mean'], 0.00095, delta=1e-4)
        offsets = np.asarray(data['CMR2_offset_per_point'])
        self.assertAlmostEqual(float(np.mean(offsets)), 0.00095, delta=1e-4)


@unittest.skipUnless(
    _callisto_cache_exists(),
    'Callisto structure cache not present in this checkout',
)
class CallistoUnchangedTests(_RNGIsolatedTestCase):
    """No-sidecar / no-guard-key configs (Callisto) must be byte-identical
    to the pre-Phase-2 derivation, proven both "by construction" (no
    sidecar file; density_inversion_guard absent) and numerically on
    seeded prior draws."""

    @classmethod
    def setUpClass(cls):
        config = InferenceConfig.from_dict(_load_config_dict(CALLISTO_CONFIG_PATH))
        cls.runner = MCMCRunner(config)

    def test_no_sidecar_found_for_callisto(self):
        self.assertIsNone(self.runner._load_cmr2_offset_sidecar())
        self.assertFalse(_callisto_cache_path().with_name(
            _callisto_cache_path().stem + '_offsets.json').exists())

    def test_density_inversion_guard_key_absent_for_callisto(self):
        rho_sil_cfg = (self.runner.config.derived_params or {}).get('rho_sil_kgm3', {})
        self.assertFalse(bool(rho_sil_cfg.get('density_inversion_guard', False)))

    def test_derived_cmr2_matches_pre_phase2_reference_on_seeded_draws(self):
        runner = self.runner
        np.random.seed(20260708)
        theta_draws = runner.prior.rvs(10)
        for t in theta_draws:
            td = runner._expand_theta(t)
            post = runner._derive_cmr2_via_mass_conservation(td)
            ref = _reference_cmr2_pre_phase2(runner, td)
            if post is None or ref is None:
                self.assertEqual(post, ref)
            else:
                self.assertAlmostEqual(post, ref, places=12)

    def test_log_likelihood_unchanged_on_seeded_draws(self):
        """Directly compares log_likelihood_fn(theta) computed with the
        real (post-edit) method against the same call with
        _derive_cmr2_via_mass_conservation monkeypatched to the verbatim
        pre-Phase-2 reference implementation. Equality here proves Phase 2
        made no numeric difference to the Callisto likelihood (its only
        possible effect funnels through that one method, and Callisto has
        neither a sidecar nor the guard key)."""
        runner = self.runner
        np.random.seed(20260708)
        theta_draws = runner.prior.rvs(10)

        actual_lls = [runner.log_likelihood_fn(t) for t in theta_draws]

        def _bound_reference(self_runner, theta_dict):
            return _reference_cmr2_pre_phase2(self_runner, theta_dict)

        with mock.patch.object(
            MCMCRunner, '_derive_cmr2_via_mass_conservation',
            autospec=True, side_effect=_bound_reference,
        ):
            reference_lls = [runner.log_likelihood_fn(t) for t in theta_draws]

        for actual, reference in zip(actual_lls, reference_lls):
            if np.isnan(actual) or np.isnan(reference):
                self.assertTrue(np.isnan(actual) and np.isnan(reference))
            else:
                self.assertEqual(actual, reference)


@unittest.skipUnless(
    _test52_cache_exists(),
    'Test52 structure cache/sidecar not present in this checkout',
)
class DensityInversionGuardTests(_RNGIsolatedTestCase):
    """Guard mechanism itself (plan decision 3), exercised with an
    out-of-prior rho_core_kgm3=2500 (below the achievable no-core
    rho_sil~2554 -- see module docstring) to force a concrete inversion."""

    @classmethod
    def setUpClass(cls):
        cls.config_dict = _load_config_dict(TEST52_CONFIG_PATH)
        # Verified numerically: R_core_km=5, rho_core_kgm3=2500 at Tb=249.0
        # gives derived rho_sil ~= 2554.1 kg/m^3 > rho_core -> inversion.
        cls.inversion_theta = _base_test52_theta(
            Tb_K=249.0, R_core_km=5.0, rho_core_kgm3=2500.0
        )

    def test_inversion_theta_actually_inverts(self):
        """Sanity precondition: confirm the constructed theta really is a
        density inversion before testing the guard's response to it."""
        config = InferenceConfig.from_dict(copy.deepcopy(self.config_dict))
        runner = MCMCRunner(config)
        structure_data = runner.structure_data
        grid_Tb_values = np.asarray(structure_data['Tb_K_grid'])
        idx = int(np.argmin(np.abs(grid_Tb_values - self.inversion_theta['Tb_K'])))
        struct = structure_data['structures'][idx]
        hydro_layers, R_oceanbot_m, M_hydro_kg = extract_hydrosphere_layers(struct)
        rho_sil, accepted = derive_silicate_density(
            M_total_kg=float(struct['Mtot_kg']),
            M_hydrosphere_kg=M_hydro_kg,
            R_oceanbot_m=R_oceanbot_m,
            R_core_m=self.inversion_theta['R_core_km'] * 1000.0,
            rho_core_kgm3=self.inversion_theta['rho_core_kgm3'],
            bounds=(1800.0, 3500.0),
        )
        self.assertTrue(accepted)
        self.assertGreater(rho_sil, self.inversion_theta['rho_core_kgm3'])

    def test_guard_on_rejects_via_none_and_likelihood_minus_1e30(self):
        config = InferenceConfig.from_dict(copy.deepcopy(self.config_dict))
        self.assertTrue(
            config.derived_params['rho_sil_kgm3']['density_inversion_guard']
        )
        runner = MCMCRunner(config)

        derived = runner._derive_cmr2_via_mass_conservation(self.inversion_theta)
        self.assertIsNone(derived)
        self.assertEqual(runner._last_cmr2_reject_reason, 'density_inversion')

        theta_arr = np.array(
            [self.inversion_theta[name] for name in runner.param_names]
        )
        ll = runner.log_likelihood_fn(theta_arr)
        self.assertEqual(ll, -1e30)

    def test_guard_off_does_not_reject_same_theta(self):
        data = copy.deepcopy(self.config_dict)
        data['derived_params']['rho_sil_kgm3']['density_inversion_guard'] = False
        config = InferenceConfig.from_dict(data)
        runner = MCMCRunner(config)

        derived = runner._derive_cmr2_via_mass_conservation(self.inversion_theta)
        self.assertIsNotNone(derived)
        self.assertTrue(np.isfinite(derived))
        self.assertIsNone(runner._last_cmr2_reject_reason)

    def test_sbi_dataset_guard_counts_as_support_rejection_when_enabled(self):
        """generate_sbi_dataset: with apply_support_guard=True and a
        rho_core prior narrowed below the achievable rho_sil (forcing
        near-certain inversion), rejected rows must be dropped and counted
        under n_rejected_support -- not silently NaN'd."""
        data = copy.deepcopy(self.config_dict)
        data['param_space']['rho_core_kgm3']['bounds'] = [2000.0, 2100.0]
        config = InferenceConfig.from_dict(data)
        runner = MCMCRunner(config)

        theta, x = runner.generate_sbi_dataset(
            n_samples=15, apply_support_guard=True, seed=20260708
        )
        stats = runner.last_sbi_dataset_stats
        self.assertGreater(stats['n_rejected_support'], 0)
        self.assertEqual(
            theta.shape[0], 15 - stats['n_rejected_support'] - stats['n_rejected_nonfinite']
        )

    def test_sbi_dataset_guard_inert_when_support_guard_disabled(self):
        """Documented behavior: the density-inversion guard only affects
        generate_sbi_dataset's row-dropping/counting when
        apply_support_guard=True. With the default False, rejected samples
        still flow through as NaN in the CMR2 column (same as any other
        mass-conservation reject cause), and no rows are dropped/counted."""
        data = copy.deepcopy(self.config_dict)
        data['param_space']['rho_core_kgm3']['bounds'] = [2000.0, 2100.0]
        config = InferenceConfig.from_dict(data)
        runner = MCMCRunner(config)

        theta, x = runner.generate_sbi_dataset(
            n_samples=15, apply_support_guard=False, seed=20260708
        )
        stats = runner.last_sbi_dataset_stats
        self.assertEqual(stats['n_rejected_support'], 0)
        self.assertEqual(theta.shape[0], 15)
        obs_names = list(runner.config.observables.keys())
        cmr2_idx = obs_names.index('CMR2')
        self.assertTrue(np.any(np.isnan(x[:, cmr2_idx])))


if __name__ == '__main__':
    unittest.main()
