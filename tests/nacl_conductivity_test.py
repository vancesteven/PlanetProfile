"""Validation of the Pan et al. (2021) NaCl(aq) conductivity model against
the authors' own regression dataset (Mendeley Data
https://doi.org/10.17632/g43xkvm3gx.6, 69 experimental rows spanning
1/5/10 wt%, 255-295 K, 212-817 MPa), committed at
tests/data/pan2021_validation.csv with columns
molality_molkg, T_K, P_MPa, rho_gcm3, sigma_Sm (rho = SeaFreeze water
density as used by the paper; sigma = measured).

The equation test uses the dataset's own rho (isolates the regression
formula); the class test exercises the SeaFreeze-density path.
"""
import csv
import sys
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from PlanetProfile.Thermodynamics.NaCl.NaClProps import (
    NaClConductPan2021, Pan2021_sigma_Sm,
)

DATA = REPO_ROOT / 'tests' / 'data' / 'pan2021_validation.csv'


def _load_rows():
    with open(DATA) as f:
        return [tuple(float(v) for v in row.values())
                for row in csv.DictReader(f)]


class Pan2021EquationTests(unittest.TestCase):

    def test_reproduces_authors_model_exactly(self):
        """PRIMARY formula test: our implementation must reproduce the
        authors' own regression-sheet sigma_calc column. Most rows match to
        machine precision; a handful differ at the third decimal from
        rounding in the cached spreadsheet values."""
        rows = _load_rows()
        self.assertGreaterEqual(len(rows), 60)
        resid = np.abs([np.log10(Pan2021_sigma_Sm(c, T, rho))
                        - np.log10(sig_calc)
                        for c, T, P, rho, _, sig_calc in rows])
        self.assertLess(np.median(resid), 1e-10,
                        f'median |dlog10| vs authors sigma_calc {np.median(resid):.2e}')
        self.assertLess(np.max(resid), 0.03,
                        f'max |dlog10| vs authors sigma_calc {np.max(resid):.3f}')

    def test_scatter_vs_measured_matches_paper(self):
        """Model vs MEASURED sigma: the residual scatter must match the
        paper's own fit quality (R^2 = 0.993, ~10% geometry error = 0.041
        in log10; tails to ~0.2 near freezing per their Fig. 3)."""
        rows = _load_rows()
        resid = np.abs([np.log10(Pan2021_sigma_Sm(c, T, rho))
                        - np.log10(sigma_meas)
                        for c, T, P, rho, sigma_meas, _ in rows])
        self.assertLess(np.median(resid), 0.05,
                        f'median |log10 residual| {np.median(resid):.4f}')
        self.assertGreater(np.mean(resid <= 0.1), 0.80,
                           f'only {np.mean(resid <= 0.1):.0%} within 0.1')
        self.assertLess(np.max(resid), 0.30,
                        f'max |log10 residual| {np.max(resid):.3f}')

    def test_physical_trends(self):
        """Conductivity falls on cooling and rises with concentration."""
        rho = 1.08
        sig_cold = Pan2021_sigma_Sm(0.9, 260.0, rho)
        sig_warm = Pan2021_sigma_Sm(0.9, 290.0, rho)
        self.assertLess(sig_cold, sig_warm)
        sig_dilute = Pan2021_sigma_Sm(0.17, 280.0, rho)
        sig_conc = Pan2021_sigma_Sm(1.9, 280.0, rho)
        self.assertLess(sig_dilute, sig_conc)


class NaClConductClassTests(unittest.TestCase):

    def test_seafreeze_density_path(self):
        """Full class (SeaFreeze water density): a few dataset rows within a
        slightly looser tolerance (their rho and ours may differ at the
        third decimal)."""
        model = NaClConductPan2021(100.0)  # 100 ppt = 1.901 mol/kg
        self.assertAlmostEqual(model.molality_molkg, 1.901, places=2)
        rows = [r for r in _load_rows() if abs(r[0] - 1.901) < 0.01][:5]
        self.assertGreaterEqual(len(rows), 3)
        for c, T, P, rho, sigma_meas, _ in rows:
            sigma = model(P, T)
            self.assertTrue(np.isfinite(sigma) and sigma > 0)
            self.assertLess(abs(np.log10(sigma) - np.log10(sigma_meas)), 0.15,
                            f'P={P} T={T}: model {sigma:.3f} vs meas {sigma_meas:.3f}')

    def test_ocean_magnitude(self):
        """A 100 ppt NaCl Callisto-like ocean is conductive (order 1-10 S/m),
        not the 1e-5 S/m placeholder it used to get."""
        model = NaClConductPan2021(100.0)
        sigma = model(300.0, 260.0)
        self.assertGreater(sigma, 0.5)
        self.assertLess(sigma, 30.0)


if __name__ == '__main__':
    unittest.main()
