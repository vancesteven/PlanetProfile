import unittest

import numpy as np

from PlanetProfile.Inference.forward_models import (
    compute_tidal_log_likelihood,
    select_structure_from_grid,
)


class InferenceDispatchTests(unittest.TestCase):
    def test_select_structure_from_grid_uses_nearest_tb(self):
        grid_cache = {
            255.0: {'CMR2': 0.33, 'marker': 'low'},
            260.0: {'CMR2': 0.34, 'marker': 'mid'},
            270.0: {'CMR2': 0.35, 'marker': 'high'},
        }
        structure_data = {
            'grid_cache': grid_cache,
            'grid_Tb_values': np.array(sorted(grid_cache)),
        }

        selected, nearest = select_structure_from_grid(structure_data, 258.2)

        self.assertEqual(nearest, 260.0)
        self.assertEqual(selected['marker'], 'mid')
        self.assertIn('grid_cache', selected)
        self.assertIn('grid_Tb_values', selected)
        self.assertEqual(selected['selected_Tb_K'], 260.0)

    def test_tidal_likelihood_supports_abs_imaginary_constraint(self):
        observables = {
            'Re_k2': (0.608, 0.048),
            'abs_Im_k2': (0.135, 0.035),
        }

        log_like = compute_tidal_log_likelihood(0.608, -0.135, observables)

        self.assertEqual(log_like, 0.0)

    def test_tidal_likelihood_accumulates_legacy_imaginary_key(self):
        observables = {
            'Im_k2': (0.100, 0.010),
            'k2': (0.5, 0.1),
        }

        log_like = compute_tidal_log_likelihood(0.3, -0.2, observables)
        expected_chi2 = ((0.2 - 0.1) / 0.01) ** 2
        expected_chi2 += ((np.sqrt(0.3**2 + 0.2**2) - 0.5) / 0.1) ** 2

        self.assertAlmostEqual(log_like, -0.5 * expected_chi2)


if __name__ == '__main__':
    unittest.main()
