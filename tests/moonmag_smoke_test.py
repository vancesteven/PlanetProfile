"""Smoke tests for the vendored, in-tree MoonMag package.

MoonMag is maintained as first-party PlanetProfile code (see
PlanetProfile/MagneticInduction/MoonMag/README.md). These tests exist to catch
dependency API drift early: the SciPy 1.17 removal of scipy.special.sph_harm
broke every full forward-model run via the GetConfig -> asymmetry_funcs import
chain and was only discovered during a structure-cache build. Import and
evaluate the affected paths directly so the next break fails here instead.
"""

import unittest

import numpy as np


class MoonMagImportTests(unittest.TestCase):
    def test_asymmetry_funcs_imports(self):
        # The module that broke under scipy 1.17 (sph_harm removal).
        from PlanetProfile.MagneticInduction.MoonMag import asymmetry_funcs  # noqa: F401

    def test_getconfig_import_chain(self):
        # GetConfig pulls in MoonMag at PlanetProfile import time; this is the
        # chain every forward-model/cache-build run exercises.
        import PlanetProfile.GetConfig  # noqa: F401

    def test_symmetry_funcs_imports(self):
        from PlanetProfile.MagneticInduction.MoonMag import symmetry_funcs  # noqa: F401


class SphericalHarmonicEvalTests(unittest.TestCase):
    def test_eval_dev_finite_and_stable(self):
        """eval_dev wraps the sph_harm_y call migrated from sph_harm.

        Reference value computed 2026-07-08 with scipy 1.17.1 immediately after
        the migration, which was itself verified equivalent to the pre-1.17
        sph_harm to ~1e-15 over 500 random (degree<=9) evaluations.
        """
        from PlanetProfile.MagneticInduction.MoonMag.asymmetry_funcs import eval_dev

        val = eval_dev(3, 2, 1.0 + 0.5j, np.array([0.7]), np.array([1.2]), (1,))
        self.assertTrue(np.all(np.isfinite(val)))
        self.assertAlmostEqual(float(val[0]), -0.34877145, places=6)

    def test_eval_dev_zero_coefficient_short_circuit(self):
        from PlanetProfile.MagneticInduction.MoonMag.asymmetry_funcs import eval_dev

        val = eval_dev(3, 2, 0, np.array([0.7]), np.array([1.2]), (1,))
        self.assertTrue(np.array_equal(val, np.zeros((1,))))

    def test_orthonormality_spot_check(self):
        """Y_2^1 integrates to unit norm over the sphere (quadrature check).

        Guards against a future harmonic-convention regression (normalization
        or theta/phi swap), which a single point value might miss.
        """
        from scipy.special import sph_harm_y

        n_t, n_p = 200, 400
        tht = np.linspace(0, np.pi, n_t)
        phi = np.linspace(0, 2 * np.pi, n_p)
        TT, PP = np.meshgrid(tht, phi, indexing='ij')
        Y = sph_harm_y(2, 1, TT, PP)
        integrand = np.abs(Y) ** 2 * np.sin(TT)
        norm = np.trapezoid(np.trapezoid(integrand, phi, axis=1), tht)
        self.assertAlmostEqual(norm, 1.0, places=3)


if __name__ == '__main__':
    unittest.main()
