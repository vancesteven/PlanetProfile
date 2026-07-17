"""Tests for the CustomSolution composition-string builder feeding
Planet.Ocean.comp (format consumed verbatim by SpeciesParser)."""
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'PlanetProfileApp'))

from Utilities.custom_solution import build_custom_solution_comp


class BuildCompTests(unittest.TestCase):

    def test_format_matches_legacy(self):
        comp, kept, warns = build_custom_solution_comp(
            'NaClMgTest', [('Na+', 0.5), ('Cl-', 0.5), ('Mg+2', 0.05)])
        self.assertEqual(
            comp, 'CustomSolutionNaClMgTest = Na+: 0.5, Cl-: 0.5, Mg+2: 0.05')
        self.assertEqual(len(kept), 3)
        self.assertEqual(warns, [])

    def test_null_zero_nan_rows_mean_absent(self):
        comp, kept, _ = build_custom_solution_comp(
            'X', [('Na+', 0.5), ('Cl-', 0.0), ('K+', None),
                  ('SO4-2', math.nan), ('', 1.0), (None, 1.0)])
        self.assertEqual(comp, 'CustomSolutionX = Na+: 0.5')
        self.assertEqual(kept, [('Na+', 0.5)])

    def test_duplicates_keep_last_with_warning(self):
        comp, kept, warns = build_custom_solution_comp(
            'X', [('Na+', 0.1), ('Na+', 0.9)])
        self.assertEqual(comp, 'CustomSolutionX = Na+: 0.9')
        self.assertEqual(len(warns), 1)

    def test_no_name_or_no_rows_gives_none(self):
        comp, kept, _ = build_custom_solution_comp('', [('Na+', 0.5)])
        self.assertIsNone(comp)
        self.assertEqual(kept, [('Na+', 0.5)])
        comp2, kept2, _ = build_custom_solution_comp('X', [('Na+', 0.0)])
        self.assertIsNone(comp2)
        self.assertEqual(kept2, [])


if __name__ == '__main__':
    unittest.main()
