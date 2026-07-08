"""
Regression test for the CMR2 reporting-path bug (2026-07 fix).

Bug (diagnosed, fixed 2026-07): ``InferenceResult.cmr2_results`` (filled in
``MCMCRunner.run()``) and the SBI dataset's ``CMR2`` column (filled in
``MCMCRunner.generate_sbi_dataset``) both called
``MCMCRunner._get_cache_scalar(theta_dict, 'CMR2')`` unconditionally. For v2
(mass-conservation) configs -- e.g.
``configs/callisto_nacl_andrade_8D.json`` -- this ignored the sampled core
parameters (``R_core_km``, ``rho_core_kgm3``) entirely and reported the
cache-builder's no-core placeholder, which only tracks ``Tb_K``. The
likelihood itself was never affected: it already re-derived CMR2 via mass
conservation for its own chi-squared term
(``_make_flexible_log_likelihood``'s ``'CMR2'`` observable branch). Only the
two reporting paths were wrong, causing e.g. Callisto runs to report CMR2 ~=
0.3412 (-3.3 sigma vs. the observation) while the likelihood the sampler
actually evaluated saw CMR2 ~= 0.3558 (+0.22 sigma).

Fix: ``MCMCRunner._compute_model_cmr2(theta_dict)`` is now the single shared
call site for both reporting paths (``run()``'s ``cmr2_results`` loop and
``generate_sbi_dataset``'s ``CMR2`` branch). It dispatches on
``config.derived_params.rho_sil_kgm3.derivation``:

- ``'mass_conservation'``: calls ``_derive_cmr2_via_mass_conservation``, the
  same derivation the likelihood uses (``extract_hydrosphere_layers`` +
  ``derive_silicate_density`` + ``compute_cmr2``), factored out of the
  likelihood as a pure restructure (verified bit-identical old-vs-new
  likelihood values on 20 random Callisto draws during development; not
  re-run here since it required a throwaway copy of the pre-fix module).
- absent / anything else: falls back to
  ``_get_cache_scalar(theta_dict, 'CMR2')`` exactly as before -- v1
  (non-derived) configs, e.g. Titan Test50, are unaffected by this fix.

This test:

1. Uses the real Callisto v2 config and its real structure cache (present
   in this checkout at
   ``PlanetProfile/Test/mcmc_results/Callisto/C2_andrade/
   callisto_nacl_structure_grid.pkl``) to assert the reported CMR2 varies
   with ``R_core_km`` across ~100 prior draws (Pearson |r| > 0.3), and that
   the old buggy path (``_get_cache_scalar``) does *not* vary with core
   params on the same draws -- confirming both the bug and the fix in one
   test. A second, more direct test perturbs a single accepted draw's
   ``R_core_km`` alone and shows the reported CMR2 changes while the old
   path's value does not.
2. Uses the Test50 Titan (non-derived) config with a synthetic grid-cache
   fixture (same fixture pattern as ``tests/sbi_runner_test.py`` -- the
   real Test50 ``.pkl`` is not present in this checkout) to assert
   ``_compute_model_cmr2`` is identical to ``_get_cache_scalar`` for a
   non-derived config, i.e. the fix is a no-op there.
"""
import json
import pickle
import tempfile
import unittest
from pathlib import Path

import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.mcmc_runner import MCMCRunner

REPO_ROOT = Path(__file__).resolve().parents[1]
CALLISTO_CONFIG_PATH = (
    REPO_ROOT / 'PlanetProfile' / 'Inference' / 'configs'
    / 'callisto_nacl_andrade_8D.json'
)
TITAN_CONFIG_PATH = (
    REPO_ROOT / 'PlanetProfile' / 'Inference' / 'configs'
    / 'test50_titan_noocean_andrade_8D.json'
)

TB_LOW = 249.0
TB_HIGH = 250.965


def _load_config_dict(path):
    with open(path) as f:
        return json.load(f)


def _callisto_cache_path():
    data = _load_config_dict(CALLISTO_CONFIG_PATH)
    return REPO_ROOT / data['structure_cache_path']


def _callisto_cache_exists():
    return _callisto_cache_path().exists()


# --- Titan synthetic fixture (mirrors tests/sbi_runner_test.py) -----------

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
        'region_phases': (tag,),
        'CMR2': 0.31,
        'Mtot_kg': 1.3452e23,
        'D_ocean_km': 0.0 if tag == 'lowT' else 60.0,
        'D_iceIh_km': 190.0,
    }


def _write_fake_grid_cache(path):
    s_low = _make_synthetic_structure('lowT', T_Ih=(230.0, 235.0))
    s_high = _make_synthetic_structure('highT', T_Ih=(273.05, 273.4))
    grid = {'Tb_K_grid': np.array([TB_LOW, TB_HIGH]), 'structures': [s_low, s_high]}
    with open(path, 'wb') as f:
        pickle.dump(grid, f)


def _build_titan_runner(tmp_path):
    data = _load_config_dict(TITAN_CONFIG_PATH)
    cache_path = tmp_path / 'fake_grid.pkl'
    _write_fake_grid_cache(cache_path)
    data['structure_cache_path'] = str(cache_path)
    data['sampler_settings'] = dict(data['sampler_settings'])
    data['sampler_settings']['checkpoint_dir'] = str(tmp_path / 'checkpoints')
    config = InferenceConfig.from_dict(data)
    return MCMCRunner(config)


@unittest.skipUnless(
    _callisto_cache_exists(),
    'Callisto structure cache not present in this checkout',
)
class CMR2ReportingMassConservationTests(unittest.TestCase):
    """v2 (mass-conservation) config: reporting CMR2 must be core-sensitive."""

    @classmethod
    def setUpClass(cls):
        data = _load_config_dict(CALLISTO_CONFIG_PATH)
        config = InferenceConfig.from_dict(data)
        cls.runner = MCMCRunner(config)

    def test_reported_cmr2_varies_with_core_params(self):
        runner = self.runner
        rng_theta = runner.prior.rvs(100)
        r_core_idx = runner.param_names.index('R_core_km')

        theta_dicts = [runner._expand_theta(t) for t in rng_theta]
        reported = np.array(
            [runner._compute_model_cmr2(td) for td in theta_dicts]
        )
        old_buggy = np.array(
            [runner._get_cache_scalar(td, 'CMR2') for td in theta_dicts]
        )
        r_core_vals = rng_theta[:, r_core_idx]

        finite_mask = np.isfinite(reported)
        n_finite = int(finite_mask.sum())
        self.assertGreater(
            n_finite, 20,
            f"too few mass-conservation-accepted draws ({n_finite}/100) "
            f"to test correlation -- check rho_sil bounds/config"
        )

        pearson_r = np.corrcoef(r_core_vals[finite_mask], reported[finite_mask])[0, 1]
        self.assertGreater(
            abs(pearson_r), 0.3,
            f"fixed reporting path (_compute_model_cmr2) should correlate "
            f"with R_core_km, got r={pearson_r}"
        )

        # Sanity check on the bug itself: the OLD reporting path
        # (_get_cache_scalar) is core-blind and must show ~zero correlation
        # with R_core_km on the exact same draws.
        old_finite_mask = np.isfinite(old_buggy)
        self.assertGreater(int(old_finite_mask.sum()), 20)
        old_pearson_r = np.corrcoef(
            r_core_vals[old_finite_mask], old_buggy[old_finite_mask]
        )[0, 1]
        self.assertLess(
            abs(old_pearson_r), 0.05,
            "sanity check failed: _get_cache_scalar should be core-blind "
            "(no correlation with R_core_km) -- if this fails, the bug "
            "characterization itself needs re-examination"
        )

    def test_perturbing_core_radius_changes_reported_cmr2_not_old_path(self):
        runner = self.runner
        draws = runner.prior.rvs(30)

        base_dict = None
        base_cmr2 = np.nan
        for t in draws:
            td = runner._expand_theta(t)
            val = runner._compute_model_cmr2(td)
            if np.isfinite(val):
                base_dict, base_cmr2 = td, val
                break
        self.assertIsNotNone(
            base_dict,
            "could not find a single mass-conservation-accepted draw in "
            "30 tries -- check rho_sil bounds/config"
        )

        perturbed_dict = dict(base_dict)
        perturbed_dict['R_core_km'] = base_dict['R_core_km'] * 0.5
        perturbed_cmr2 = runner._compute_model_cmr2(perturbed_dict)

        # The perturbed sample may itself be mass-conservation-rejected
        # (NaN) -- that is still a different outcome from the finite base
        # value, so treat NaN-vs-finite as a pass too.
        if np.isfinite(perturbed_cmr2):
            self.assertNotEqual(
                base_cmr2, perturbed_cmr2,
                "reported CMR2 did not change when R_core_km was halved"
            )

        # Old buggy path is blind to R_core_km -> identical for both.
        old_base = runner._get_cache_scalar(base_dict, 'CMR2')
        old_perturbed = runner._get_cache_scalar(perturbed_dict, 'CMR2')
        self.assertEqual(
            old_base, old_perturbed,
            "sanity check failed: _get_cache_scalar should not depend on "
            "R_core_km"
        )


class CMR2ReportingNonDerivedFallbackTests(unittest.TestCase):
    """v1 (non-derived) config: reporting CMR2 must be unchanged (cache path)."""

    def test_cache_path_used_unchanged_for_non_derived_config(self):
        with tempfile.TemporaryDirectory() as td:
            tmp_path = Path(td)
            runner = _build_titan_runner(tmp_path)
            self.assertEqual(runner.config.derived_params, {})

            theta = runner.prior.rvs(20)
            for t in theta:
                theta_dict = runner._expand_theta(t)
                reported = runner._compute_model_cmr2(theta_dict)
                cache_val = runner._get_cache_scalar(theta_dict, 'CMR2')
                if np.isfinite(cache_val):
                    self.assertEqual(reported, cache_val)
                else:
                    self.assertTrue(np.isnan(reported))


if __name__ == '__main__':
    unittest.main()
