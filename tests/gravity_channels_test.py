"""Tests for the v4 geodesy gravity channels (C20/C22 via Clairaut k_f).

Covers (plans/europa-clipper-v4-geodesy-plan.md):
- _derive_gravity_pair runs the Clairaut integration on the SAME
  composite (core + mass-conserved silicate + hydrosphere) profile the
  CMR2 derivation assembles, at the Tricarico ratio -C20/C22 = 3.324.
- Sampled nuisance offsets dC20_nh/dC22_nh shift the pair additively.
- Likelihood dispatch: gravity_forward_model='clairaut_hydrostatic'
  computes the pair; legacy configs keep the cached-scalar read.
- Config-hash stability: the new field changes NOTHING when unset.

Uses the real Callisto v2 (mass-conservation) config + committed cache,
same fixture pattern as tests/cmr2_reporting_test.py.
Run: python -m unittest tests.gravity_channels_test (env PP/PPcl).
"""
import json
import sys
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402
from PlanetProfile.Inference.mcmc_runner import MCMCRunner  # noqa: E402
from PlanetProfile.Inference.gravity_obs import (  # noqa: E402
    clairaut_kf, hydrostatic_c20_c22, radau_darwin_kf, J2_OVER_C22,
)

CALLISTO_CONFIG_PATH = (
    REPO_ROOT / 'PlanetProfile' / 'Inference' / 'configs'
    / 'callisto_nacl_andrade_8D.json'
)


def _load_config_dict():
    with open(CALLISTO_CONFIG_PATH) as f:
        return json.load(f)


def _cache_exists():
    data = _load_config_dict()
    return (REPO_ROOT / data['structure_cache_path']).exists()


def _gravity_config_dict(with_nuisance=True):
    """Callisto v2 config dict adapted to the v4 gravity forward model."""
    data = _load_config_dict()
    data['gravity_forward_model'] = 'clairaut_hydrostatic'
    obs = dict(data['observables'])
    # Placeholder centrals; per-test code conditions on model predictions.
    obs['C20'] = [-1.0e-4, 8.5e-7]
    obs['C22'] = [3.0e-5, 2.0e-7]
    data['observables'] = obs
    if with_nuisance:
        ps = dict(data['param_space'])
        for nm in ('dC20_nh', 'dC22_nh'):
            ps[nm] = {'prior_type': 'uniform', 'bounds': [-2.0e-5, 2.0e-5]}
        data['param_space'] = ps
    return data


def _first_accepted_theta(runner, n=200, seed=20260719):
    np.random.seed(seed)
    for theta in runner.prior.rvs(n):
        td = runner._expand_theta(theta)
        if runner._derive_gravity_pair(td) is not None:
            return td
    raise AssertionError('no mass-conservation-accepted draw found')


@unittest.skipUnless(_cache_exists(), 'Callisto structure cache not present')
class DeriveGravityPairTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.runner = MCMCRunner(
            InferenceConfig.from_dict(_gravity_config_dict()))
        cls.td = _first_accepted_theta(cls.runner)

    def test_pair_matches_composite_profile_clairaut(self):
        """The pair must equal gravity_obs applied to the EXACT composite
        layer stack the CMR2 derivation cached (reviewer-binding)."""
        runner, td = self.runner, dict(self.td)
        td['dC20_nh'] = 0.0
        td['dC22_nh'] = 0.0
        pair = runner._derive_gravity_pair(td)
        self.assertIsNotNone(pair)
        c20, c22 = pair
        assembled, R_body_m, M_total_kg = runner._last_composite_layers
        r_arr = np.array([ro for (_, ro, _) in assembled])
        rho_arr = np.array([rh for (_, _, rh) in assembled])
        struct = runner._struct_for_hydrosphere(td)
        kf = clairaut_kf(r_arr, rho_arr)
        c20_ref, c22_ref = hydrostatic_c20_c22(
            kf, float(struct['omega']), R_body_m, M_total_kg)
        self.assertAlmostEqual(c20, c20_ref, places=12)
        self.assertAlmostEqual(c22, c22_ref, places=12)
        # Hydrostatic ratio at the Tricarico value, exact by construction
        self.assertAlmostEqual(-c20 / c22, J2_OVER_C22, places=9)
        # Clairaut k_f agrees with Radau-Darwin from the same composite's
        # CMR2 to the documented ~1-2% systematic level.
        cmr2 = runner._derive_cmr2_via_mass_conservation(dict(self.td))
        kf_rd = radau_darwin_kf(cmr2)
        self.assertLess(abs(kf - kf_rd) / kf_rd, 0.02,
                        f'Clairaut {kf:.4f} vs RD {kf_rd:.4f}')

    def test_nuisance_offsets_are_additive(self):
        runner, td0 = self.runner, dict(self.td)
        td0['dC20_nh'] = 0.0
        td0['dC22_nh'] = 0.0
        c20_0, c22_0 = runner._derive_gravity_pair(td0)
        td1 = dict(td0)
        td1['dC20_nh'] = -7.0e-6
        td1['dC22_nh'] = 1.0e-5
        c20_1, c22_1 = runner._derive_gravity_pair(td1)
        self.assertAlmostEqual(c20_1 - c20_0, -7.0e-6, places=12)
        self.assertAlmostEqual(c22_1 - c22_0, 1.0e-5, places=12)

    def test_core_sensitivity(self):
        """A different sampled core must move the hydrostatic pair (the
        Clairaut integration sees the composite, not the raw cache)."""
        runner, td = self.runner, dict(self.td)
        td['dC20_nh'] = td['dC22_nh'] = 0.0
        _, c22_a = runner._derive_gravity_pair(td)
        td2 = dict(td)
        td2['R_core_km'] = td['R_core_km'] * 0.7
        pair2 = runner._derive_gravity_pair(td2)
        if pair2 is None:
            self.skipTest('perturbed core rejected by mass conservation')
        self.assertNotAlmostEqual(c22_a, pair2[1], places=8)


@unittest.skipUnless(_cache_exists(), 'Callisto structure cache not present')
class LikelihoodDispatchTests(unittest.TestCase):
    def test_gravity_likelihood_conditions_on_computed_pair(self):
        """Log-likelihood must be maximal when the C20/C22 observation
        equals the model prediction and drop when the observation moves."""
        runner = MCMCRunner(InferenceConfig.from_dict(_gravity_config_dict()))
        td = _first_accepted_theta(runner)
        td['dC20_nh'] = td['dC22_nh'] = 0.0
        c20_m, c22_m = runner._derive_gravity_pair(td)

        def _ll_with_obs(c20_obs, c22_obs):
            data = _gravity_config_dict()
            data['observables']['C20'] = [c20_obs, 8.5e-7]
            data['observables']['C22'] = [c22_obs, 2.0e-7]
            r = MCMCRunner(InferenceConfig.from_dict(data))
            theta = [td[p] for p in r.param_names]
            return float(r.log_likelihood_fn(theta))

        ll_match = _ll_with_obs(c20_m, c22_m)
        ll_off = _ll_with_obs(c20_m, c22_m + 5 * 2.0e-7)
        self.assertGreater(ll_match, -1e29, 'valid draw hard-rejected')
        self.assertAlmostEqual(ll_match - ll_off, 0.5 * 25.0, places=5,
                               msg='5-sigma C22 offset must cost 12.5 in ll')

    def test_legacy_config_unchanged(self):
        """Without gravity_forward_model, a C22 observable still reads the
        cached scalar (NaN in this cache -> hard reject), and the config
        hash is unchanged by the new field's existence."""
        data = _load_config_dict()
        h_before = InferenceConfig.from_dict(data).generate_hash()
        data2 = dict(data)
        data2['gravity_forward_model'] = None
        h_none = InferenceConfig.from_dict(data2).generate_hash()
        self.assertEqual(h_before, h_none,
                         'unset gravity_forward_model must not change hash')
        data3 = dict(data)
        data3['gravity_forward_model'] = 'clairaut_hydrostatic'
        h_set = InferenceConfig.from_dict(data3).generate_hash()
        self.assertNotEqual(h_before, h_set)

        # Legacy cached-scalar path: C22 in observables without the flag
        # (cache C22 is NaN for this body) -> -1e30 rejection, as before.
        data4 = _load_config_dict()
        obs = dict(data4['observables'])
        obs['C22'] = [1.0e-5, 1.0e-6]
        data4['observables'] = obs
        r = MCMCRunner(InferenceConfig.from_dict(data4))
        np.random.seed(20260719)
        theta = r.prior.rvs(1)[0]
        self.assertLessEqual(float(r.log_likelihood_fn(theta)), -1e29)


if __name__ == '__main__':
    unittest.main()
