"""Acceptance + guard tests for build_cmr2_offset_sidecar.

1. Regenerating the Test52 sidecar reproduces the committed JSON to <1e-6
   per grid point (the committed values were produced by the original
   session-scratch Phase 1 code; this tool must match them exactly).
2. The validity guard refuses the pre-existing Callisto C2_andrade cache,
   whose structured silicate interior makes native-recon a physical
   artifact (-0.035), not a discretization offset (opus design review,
   2026-07-10).

Both tests skip when the (gitignored) cache pkl is absent on this machine.
"""
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from PlanetProfile.Inference import build_cmr2_offset_sidecar as builder

TEST52_CACHE = (REPO_ROOT / 'PlanetProfile' / 'Test' / 'mcmc_results' / 'Titan'
                / 'Test52_andrade_noocean_diff'
                / 'titan_diff_noocean_structure_grid.pkl')
TEST52_SIDECAR = TEST52_CACHE.with_name(TEST52_CACHE.stem + '_offsets.json')
CALLISTO_CACHE = (REPO_ROOT / 'PlanetProfile' / 'Test' / 'mcmc_results'
                  / 'Callisto' / 'C2_andrade' / 'callisto_nacl_structure_grid.pkl')


class SidecarAcceptanceTests(unittest.TestCase):

    @unittest.skipUnless(TEST52_CACHE.exists() and TEST52_SIDECAR.exists(),
                         'Test52 cache/sidecar not present on this machine')
    def test_reproduces_committed_test52_sidecar(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / 'sidecar.json'
            rc = builder.main(['--cache', str(TEST52_CACHE), '--out', str(out)])
            self.assertEqual(rc, 0)
            with open(out) as f:
                built = json.load(f)
        with open(TEST52_SIDECAR) as f:
            committed = json.load(f)

        np.testing.assert_allclose(
            built['Tb_K_grid'], committed['Tb_K_grid'], rtol=0, atol=1e-12)
        np.testing.assert_allclose(
            built['CMR2_offset_per_point'],
            committed['CMR2_offset_per_point'], rtol=0, atol=1e-6)
        self.assertEqual(built['definition'], committed['definition'])

    @unittest.skipUnless(TEST52_CACHE.exists(), 'Test52 cache not present')
    def test_refuses_overwrite_without_force(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / 'sidecar.json'
            out.write_text('{}')
            rc = builder.main(['--cache', str(TEST52_CACHE), '--out', str(out)])
            self.assertEqual(rc, 1)
            rc = builder.main(['--cache', str(TEST52_CACHE), '--out', str(out),
                               '--force'])
            self.assertEqual(rc, 0)


class SidecarValidityGuardTests(unittest.TestCase):

    @unittest.skipUnless(CALLISTO_CACHE.exists(),
                         'Callisto C2_andrade cache not present on this machine')
    def test_refuses_structured_silicate_callisto_cache(self):
        """The physical structured-vs-uniform silicate difference (-0.035)
        must never be written as a 'discretization offset'."""
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / 'sidecar.json'
            rc = builder.main(['--cache', str(CALLISTO_CACHE), '--out', str(out)])
            self.assertEqual(rc, 3)
            self.assertFalse(out.exists())


if __name__ == '__main__':
    unittest.main()
