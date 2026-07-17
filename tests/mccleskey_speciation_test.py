"""McCleskey conductivity must use Reaktoro equilibrium FREE-ION speciation.

User directive (2026-07-16): neutral complexes / ion pairs (e.g. MgSO4(aq))
remove ions from solution, so McCleskey conductivity computed from
equilibrium speciation must fall below the nominal fully-dissociated value
— the same way McCleskey (2012) is applied with WATEQ4F speciation in
terrestrial work. A dilute NaCl control (weak pairing) guards against a
broken species mapping masquerading as a speciation effect.

Requires reaktoro (conda-forge); skipped otherwise.
"""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    import reaktoro  # noqa: F401
    HAVE_RKT = True
except ImportError:
    HAVE_RKT = False


@unittest.skipUnless(HAVE_RKT, "reaktoro not installed")
class McCleskeySpeciationTests(unittest.TestCase):

    @staticmethod
    def _make_conduct(comp_string, w_ppt, use_equilibrium):
        from PlanetProfile.Thermodynamics.Reaktoro.reaktoroProps import (
            SpeciesParser, RktHydroSpecies, RktConduct)
        aq, ratio, solids, _label = SpeciesParser(comp_string, w_ppt)
        fn_species = RktHydroSpecies(aq, ratio, solids)
        return RktConduct(aq, ratio, solids, fn_species,
                          use_equilibrium_speciation=use_equilibrium)

    def test_mgso4_pairing_reduces_conductivity(self):
        comp = "CustomSolutionMgSO4Test = Mg+2: 0.1, SO4-2: 0.1"
        P = np.array([10.0, 50.0])
        T = np.array([275.0, 285.0])

        sig_nom = self._make_conduct(comp, None, False)(P, T)
        sig_eq = self._make_conduct(comp, None, True)(P, T)

        self.assertTrue(np.all(np.isfinite(sig_nom)))
        self.assertTrue(np.all(np.isfinite(sig_eq)))
        self.assertTrue(np.all(sig_nom > 0))
        # MgSO4(aq) ion pairing removes free Mg+2/SO4-2: equilibrium-path
        # conductivity must be strictly below the fully-dissociated value.
        self.assertTrue(
            np.all(sig_eq < sig_nom),
            f"expected sigma_eq < sigma_nominal everywhere; got "
            f"eq={sig_eq}, nom={sig_nom}")
        # And materially so for 0.1 mol/kg MgSO4 (pairing is strong).
        self.assertTrue(np.all(sig_eq < 0.95 * sig_nom))

    def test_dilute_nacl_control_paths_agree(self):
        comp = "CustomSolutionNaClDilute = Na+: 0.01, Cl-: 0.01"
        P = np.array([10.0, 50.0])
        T = np.array([275.0, 285.0])

        sig_nom = self._make_conduct(comp, None, False)(P, T)
        sig_eq = self._make_conduct(comp, None, True)(P, T)

        # Weak pairing: the two paths agree closely — a large gap here
        # would indicate a broken speciated-name mapping, not physics.
        rel = np.abs(sig_eq - sig_nom) / sig_nom
        self.assertTrue(np.all(rel < 0.10),
                        f"dilute NaCl paths diverged: rel diff {rel}")


if __name__ == '__main__':
    unittest.main()
