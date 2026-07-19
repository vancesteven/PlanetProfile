"""Unit tests for PlanetProfile.Inference.gravity_obs (Clipper v4 plan).

Run: python -m unittest tests/gravity_obs_test.py (env PP/PPcl).
"""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from PlanetProfile.Inference.gravity_obs import (  # noqa: E402
    clairaut_kf, radau_darwin_kf, radau_darwin_cmr2, hydrostatic_c20_c22,
    cmr2_from_c22_rd, rotation_parameter, J2_OVER_C22,
)

# Europa constants (cache values)
OMEGA = 2.0478e-5   # rad/s, synchronous spin
R_M = 1560.8e3
M_KG = 4.7998e22


class TestClairautKf(unittest.TestCase):
    def test_homogeneous_body(self):
        """Uniform density: eta_s = 0 exactly -> k_f = 3/2."""
        r = np.linspace(1.0, R_M, 200)
        rho = np.full_like(r, 3013.0)
        kf = clairaut_kf(r, rho)
        self.assertAlmostEqual(kf, 1.5, places=4)

    def test_roche_point_mass_limit(self):
        """Tiny ultra-dense core + near-massless envelope: k_f -> 0."""
        r = np.concatenate(([0.005 * R_M], np.linspace(0.01 * R_M, R_M, 200)))
        rho = np.concatenate(([1.0e9], np.full(200, 1.0)))
        kf = clairaut_kf(r, rho)
        self.assertLess(kf, 0.02)

    def test_radau_darwin_crosscheck_realistic(self):
        """3-layer Europa-like profile: Clairaut k_f agrees with the
        Radau-Darwin estimate from the profile's own C/MR^2 to ~1%
        (the documented RD systematic scale)."""
        # metal core 600 km rho 5500, rock to 1420 km rho 3300,
        # H2O shell to surface rho 1050 — bulk density ~ Europa's.
        r = np.concatenate([
            np.linspace(10e3, 600e3, 60),
            np.linspace(605e3, 1420e3, 80),
            np.linspace(1425e3, R_M, 40),
        ])
        rho = np.concatenate([
            np.full(60, 5500.0), np.full(80, 3300.0), np.full(40, 1050.0),
        ])
        kf = clairaut_kf(r, rho)
        # C/MR^2 of the same piecewise profile
        r_edges = np.concatenate(([0.0], r))
        shell_m = (4/3) * np.pi * rho * (r_edges[1:]**3 - r_edges[:-1]**3)
        M = shell_m.sum()
        shell_I = (8/15) * np.pi * rho * (r_edges[1:]**5 - r_edges[:-1]**5)
        cmr2 = shell_I.sum() / (M * r[-1]**2)
        kf_rd = radau_darwin_kf(cmr2)
        self.assertLess(abs(kf - kf_rd) / kf_rd, 0.015,
                        f"Clairaut {kf:.4f} vs RD {kf_rd:.4f} at C={cmr2:.4f}")

    def test_condensation_monotonicity(self):
        """More centrally condensed -> smaller k_f."""
        kfs = []
        for rho_core in (3013.0, 4500.0, 7000.0):
            r = np.concatenate([np.linspace(10e3, 700e3, 60),
                                np.linspace(705e3, R_M, 80)])
            # hold bulk mass roughly fixed is not required for monotonicity
            rho = np.concatenate([np.full(60, rho_core),
                                  np.full(80, 1500.0)])
            kfs.append(clairaut_kf(r, rho))
        self.assertGreater(kfs[0], kfs[1])
        self.assertGreater(kfs[1], kfs[2])


class TestRadauDarwin(unittest.TestCase):
    def test_anderson_roundtrip(self):
        """C = 0.346 -> k_f = 1.044 -> C22_h = 130e-6, J2_h = 433e-6
        (plan-verified reproduction of Anderson's hydrostatic pair;
        J2 via the Tricarico-corrected 3.324 ratio, and evaluated at
        the physical radius to match the classical numbers)."""
        kf = radau_darwin_kf(0.346)
        self.assertAlmostEqual(kf, 1.044, places=3)
        c20, c22 = hydrostatic_c20_c22(kf, OMEGA, R_M, M_KG, R_ref_m=R_M)
        self.assertAlmostEqual(c22 * 1e6, 130.0, delta=1.5)
        self.assertAlmostEqual(-c20 / c22, J2_OVER_C22, places=9)
        self.assertAlmostEqual(-c20 * 1e6, 3.324 * c22 * 1e6, delta=0.01)

    def test_inverse_consistency(self):
        for cmr2 in (0.30, 0.3405, 0.3475, 0.3547, 0.40):
            kf = radau_darwin_kf(cmr2)
            self.assertAlmostEqual(radau_darwin_cmr2(kf), cmr2, places=12)

    def test_homogeneous_limit(self):
        """C/MR^2 = 0.4 (homogeneous) -> k_f = 3/2 exactly."""
        self.assertAlmostEqual(radau_darwin_kf(0.4), 1.5, places=12)

    def test_kf_monotone_in_cmr2(self):
        grid = np.linspace(0.30, 0.40, 50)
        kfs = [radau_darwin_kf(c) for c in grid]
        self.assertTrue(np.all(np.diff(kfs) > 0))


class TestObservableMap(unittest.TestCase):
    def test_reference_radius_scaling(self):
        """(R/R_ref)^2 rescale: ~0.5% smaller C22 at 1565 km than at
        the physical radius — 3x sigma(C22), load-bearing."""
        kf = 1.044
        _, c22_ref = hydrostatic_c20_c22(kf, OMEGA, R_M, M_KG)
        _, c22_phys = hydrostatic_c20_c22(kf, OMEGA, R_M, M_KG, R_ref_m=R_M)
        self.assertAlmostEqual(c22_ref / c22_phys, (R_M / 1565.0e3) ** 2,
                               places=12)
        self.assertGreater(c22_phys, c22_ref)

    def test_gui_inverse_readout(self):
        """cmr2_from_c22_rd inverts the forward map (same convention)."""
        for cmr2 in (0.3405, 0.3475, 0.3547):
            kf = radau_darwin_kf(cmr2)
            _, c22 = hydrostatic_c20_c22(kf, OMEGA, R_M, M_KG)
            self.assertAlmostEqual(
                cmr2_from_c22_rd(c22, OMEGA, R_M, M_KG), cmr2, places=10)

    def test_rotation_parameter_europa(self):
        q = rotation_parameter(OMEGA, R_M, M_KG)
        self.assertAlmostEqual(q * 1e4, 4.98, delta=0.05)


if __name__ == '__main__':
    unittest.main()
