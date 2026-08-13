"""Standard-run parity tests for hydrostatic degree-2 gravity output."""
import sys
import unittest
from copy import deepcopy
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from PlanetProfile.Default.Europa.PPEuropa import Planet as EUROPA  # noqa: E402
from PlanetProfile.Default.Titan.PPTitan import Planet as TITAN  # noqa: E402
from PlanetProfile.Inference.gravity_obs import (  # noqa: E402
    clairaut_kf,
    hydrostatic_c20_c22,
)
from PlanetProfile.Main import PostProcessingProfile  # noqa: E402


class StandardGravityOutputTests(unittest.TestCase):
    @staticmethod
    def _completed_profile(template):
        planet = deepcopy(template)
        radius = planet.Bulk.R_m
        planet.r_m = np.array([
            radius,
            0.94 * radius,
            0.76 * radius,
            0.34 * radius,
            0.0,
        ])
        planet.rho_kgm3 = np.array([1000.0, 1500.0, 3400.0, 6500.0])
        return planet

    def test_default_flag_leaves_gravity_outputs_unset(self):
        planet = self._completed_profile(EUROPA)
        self.assertFalse(planet.Do.CALC_C20_C22)

        result = PostProcessingProfile(planet, None)

        self.assertTrue(np.isnan(result.Gravity.kf))
        self.assertTrue(np.isnan(result.Gravity.C20))
        self.assertTrue(np.isnan(result.Gravity.C22))

    def test_europa_and_titan_match_direct_gravity_observable_map(self):
        for template in (EUROPA, TITAN):
            with self.subTest(body=template.name):
                planet = self._completed_profile(template)
                planet.Do.CALC_C20_C22 = True
                shell_radii = planet.r_m[:-1]
                expected_kf = clairaut_kf(shell_radii, planet.rho_kgm3)
                expected_c20, expected_c22 = hydrostatic_c20_c22(
                    expected_kf,
                    2.0 * np.pi / planet.Bulk.Torb_s,
                    planet.Bulk.R_m,
                    planet.Bulk.M_kg,
                    R_ref_m=planet.Gravity.Rref_m,
                    j2_over_c22=planet.Gravity.J2overC22,
                )

                with self.assertLogs('PlanetProfile', level='INFO') as logs:
                    result = PostProcessingProfile(planet, None)

                self.assertEqual(result.Gravity.kf, expected_kf)
                self.assertEqual(result.Gravity.C20, expected_c20)
                self.assertEqual(result.Gravity.C22, expected_c22)
                printout = '\n'.join(logs.output)
                self.assertIn('Hydrostatic degree-2 gravity', printout)
                self.assertIn('C20 =', printout)
                self.assertIn('C22 =', printout)


if __name__ == '__main__':
    unittest.main()
