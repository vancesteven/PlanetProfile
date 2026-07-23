"""Config-load guards for the Europa Clipper v6 free-gravity trio.

Guards the reviewer's block-level correction (2026-07-22): GC21 Table 2 is
already UNNORMALIZED at R_ref=1565 km, so C22 must be entered directly as
~138.6e-6 with NO sqrt(5)/sqrt(10/24) fully-normalized conversion. A C22 near
89e-6 would signal the double-conversion bug. Also pins the v6 invariants:
CMR2 dropped, C20/C22 = SOL-A, and the widened agnostic offset boxes.

Run: python -m unittest tests/v6_config_test.py (env PPcl).
"""
import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

_CFGDIR = Path(__file__).resolve().parents[1] / "PlanetProfile/Inference/configs"
_ARMS = {
    "europa_clipper_v6_freegrav_11D.json": 20,
    "europa_clipper_v6_freegrav_noinduction_6obs.json": 6,
    "europa_clipper_v6_freegrav_nok2_16obs.json": 16,
}

# GC21 Table 2 SOL-A, unnormalized at 1565 km.
_C20 = -4.3759e-4
_C22 = 1.3862e-4
_C22_DOUBLE_CONV_BUG = 8.95e-5  # what sqrt(10/24)*C22 would wrongly give (~35% low)


class TestV6Configs(unittest.TestCase):
    def _load(self, fname):
        with open(_CFGDIR / fname) as f:
            return json.load(f)

    def test_all_arms_present(self):
        for fname in _ARMS:
            self.assertTrue((_CFGDIR / fname).exists(), f"missing {fname}")

    def test_cmr2_dropped(self):
        for fname in _ARMS:
            obs = self._load(fname)["observables"]
            self.assertNotIn("CMR2", obs, f"{fname}: CMR2 must not be an observable in v6")

    def test_obs_count(self):
        for fname, n in _ARMS.items():
            obs = self._load(fname)["observables"]
            self.assertEqual(len(obs), n, f"{fname}: expected {n} obs, got {len(obs)}")

    def test_c20_c22_sola_unnormalized(self):
        """Reviewer block-level guard: C22 ~138.6e-6, NOT ~89e-6 (conversion bug)."""
        for fname in _ARMS:
            obs = self._load(fname)["observables"]
            self.assertAlmostEqual(obs["C22"][0], _C22, places=8,
                                   msg=f"{fname}: C22 must be SOL-A unnormalized")
            self.assertAlmostEqual(obs["C20"][0], _C20, places=8,
                                   msg=f"{fname}: C20 must be SOL-A unnormalized")
            # explicit conversion-bug tripwire
            self.assertGreater(obs["C22"][0], 1.2e-4,
                               f"{fname}: C22={obs['C22'][0]:.3e} looks fully-normalized "
                               f"(~{_C22_DOUBLE_CONV_BUG:.2e}) -- double-conversion bug")

    def test_sola_ratio_nonhydrostatic(self):
        """-C20/C22 = 3.157 (SOL-A), NOT the hydrostatic 3.324 -- the free signal."""
        for fname in _ARMS:
            obs = self._load(fname)["observables"]
            ratio = -obs["C20"][0] / obs["C22"][0]
            self.assertAlmostEqual(ratio, 3.157, places=2, msg=f"{fname}: ratio {ratio}")

    def test_agnostic_offset_boxes(self):
        """Both offset boxes widened so gravity carries ~zero interior C/MR2 info."""
        for fname in _ARMS:
            ps = self._load(fname)["param_space"]
            self.assertEqual(ps["dC20_nh"]["bounds"], [-3.9e-4, 3.9e-4], f"{fname} dC20_nh")
            self.assertEqual(ps["dC22_nh"]["bounds"], [-5.0e-5, 5.0e-5], f"{fname} dC22_nh")

    def test_loads_through_inference_config(self):
        from PlanetProfile.Inference.inference_core import InferenceConfig
        for fname in _ARMS:
            c = InferenceConfig.from_json(str(_CFGDIR / fname))
            self.assertNotIn("CMR2", c.observables)
            self.assertIn("C20", c.observables)
            self.assertIn("C22", c.observables)


if __name__ == "__main__":
    unittest.main()
