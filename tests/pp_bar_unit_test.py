"""Regression tests for the bar-unit correction in PlanetProfile.

Verifies the smallest safe refactor documented in
lunar-bsm-interior-models/docs/PLANETPROFILE_BAR_UNIT_AUDIT.md, § 5:

- Constants.bar2MPa = 1.0e-1 (exact SI, 1 bar = 0.1 MPa)
- Constants.bar2GPa = 1.0e-4 (exact SI, 1 bar = 1e-4 GPa)
- Constants.P0 remains 0.101325 MPa (one standard atmosphere)
- Thermodynamics/Seawater/SwProps.MPa2seaPressure subtracts P0
  (not bar2MPa) for the atmospheric offset.

Runnable with the stdlib unittest runner (no pytest dependency):

    python -m unittest tests.pp_bar_unit_test

The Perple_X-loader synthetic-table tests require numpy and scipy in
the runtime environment; they are skipped automatically if those
imports fail.
"""

from __future__ import annotations

import os
import sys
import tempfile
import unittest

# Prefer running from the PP repo root so the package resolves without
# environment surgery, but allow running from anywhere.
_HERE = os.path.dirname(os.path.abspath(__file__))
_PP_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
if _PP_ROOT not in sys.path:
    sys.path.insert(0, _PP_ROOT)

try:
    import numpy as np
    _HAS_NUMPY = True
except Exception:
    _HAS_NUMPY = False

try:
    import scipy  # noqa: F401
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False


class TestExactConstants(unittest.TestCase):
    """Phase-4 test 1: exact scalar conversions."""

    def test_bar2MPa_is_exact_physical(self):
        from PlanetProfile.Utilities.defineStructs import Constants
        self.assertEqual(Constants.bar2MPa, 1.0e-1)

    def test_bar2GPa_is_exact_physical(self):
        from PlanetProfile.Utilities.defineStructs import Constants
        self.assertEqual(Constants.bar2GPa, 1.0e-4)

    def test_1000_bar_equals_100_MPa(self):
        from PlanetProfile.Utilities.defineStructs import Constants
        self.assertAlmostEqual(1000.0 * Constants.bar2MPa, 100.0, places=12)

    def test_one_million_bar_equals_100_GPa(self):
        from PlanetProfile.Utilities.defineStructs import Constants
        self.assertAlmostEqual(1.0e6 * Constants.bar2GPa, 100.0, places=9)

    def test_P0_remains_one_standard_atmosphere_in_MPa(self):
        from PlanetProfile.Utilities.defineStructs import Constants
        # 1 atm = 101 325 Pa = 0.101325 MPa.
        self.assertEqual(Constants.P0, 0.101325)

    def test_bar2MPa_and_P0_are_distinct(self):
        """The atmospheric offset and the bar-to-MPa conversion are
        different physical quantities and must not be numerically equal
        after the correction."""
        from PlanetProfile.Utilities.defineStructs import Constants
        self.assertNotAlmostEqual(Constants.bar2MPa, Constants.P0, places=6)


class TestSeawaterConversion(unittest.TestCase):
    """Phase-4 test 2: seawater atmospheric offset and dbar conversion."""

    def test_sea_pressure_at_one_atm_absolute_is_zero_dbar(self):
        from PlanetProfile.Utilities.defineStructs import Constants
        from PlanetProfile.Thermodynamics.Seawater.SwProps import MPa2seaPressure
        # Absolute pressure = 1 atm exactly => sea pressure = 0 dbar.
        self.assertAlmostEqual(
            float(MPa2seaPressure(Constants.P0)), 0.0, places=10
        )

    def test_one_MPa_above_one_atm_is_100_dbar(self):
        from PlanetProfile.Utilities.defineStructs import Constants
        from PlanetProfile.Thermodynamics.Seawater.SwProps import MPa2seaPressure
        # Absolute pressure = 1 atm + 1 MPa = 0.101325 + 1.0 = 1.101325 MPa
        # Sea pressure = 1 MPa = 10 bar = 100 dbar exactly (physical).
        result = float(MPa2seaPressure(Constants.P0 + 1.0))
        self.assertAlmostEqual(result, 100.0, places=10)

    def test_atmospheric_offset_is_101325_Pa_not_100000_Pa(self):
        """After the refactor MPa2seaPressure(P0) must be exactly zero,
        which is only true if the code subtracts P0 (0.101325 MPa)
        rather than bar2MPa (0.1 MPa)."""
        from PlanetProfile.Utilities.defineStructs import Constants
        from PlanetProfile.Thermodynamics.Seawater.SwProps import MPa2seaPressure
        # Subtracting one bar (0.1 MPa) from P0 = 0.101325 MPa would
        # leave 0.001325 MPa = 0.1325 dbar of residual sea pressure.
        # After the correction the residual must be zero.
        residual_dbar = float(MPa2seaPressure(Constants.P0))
        self.assertAlmostEqual(residual_dbar, 0.0, places=10)
        # And to confirm the refactor did *not* accidentally leave the
        # atm offset at 1 bar, the result at absolute = 1 bar exactly
        # (0.1 MPa) must be negative (sub-atm) and equal to -0.1325 dbar.
        sub_atm = float(MPa2seaPressure(Constants.bar2MPa))
        self.assertAlmostEqual(sub_atm, -0.1325, places=6)


class TestPerpleXNativeLoader(unittest.TestCase):
    """Phase-4 test 3: PerplexEOSStruct loads a synthetic native table.

    Verifies that PP's Perple_X-native loader applies the corrected
    physical conversions exactly:

    - Pmin = 1000 bar loads as 100.0 MPa (not 101.325 MPa).
    - Ks_bar = 1e6 loads as 100.0 GPa (not 101.325 GPa).
    - Gs_bar = 5e5 loads as  50.0 GPa (not  50.6625 GPa).
    - rho, Vp, Vs, alpha, Cp behave correctly (Cp is converted from
      J/K/m3 to J/kg/K by dividing by rho on load).
    - EXTRAP=False clamps out-of-domain queries to the nearest boundary
      and does not mutate the caller's P, T arrays.
    """

    def _write_synthetic_table(self, path, *,
                               P_bar_min, P_bar_step, P_n,
                               T_K_min, T_K_step, T_n,
                               rho_kgm3, VP_kms, VS_kms,
                               Cp_JKm3, alpha_pK, KS_bar, GS_bar,
                               first_axis="T"):
        """Emit a Perple_X-native 13-line-header tab file with constant
        property values so the interpolator returns them exactly at any
        interior query."""
        lines = [
            "|6.6.6",
            "synthetic_bar_unit_test.tab" + " " * 60,
            f"{2:12d}",
        ]
        # Axis blocks in the requested order.
        axes = [("T(K)", T_K_min, T_K_step, T_n),
                ("P(bar)", P_bar_min, P_bar_step, P_n)]
        if first_axis == "P":
            axes = list(reversed(axes))
        for label, v0, step, n in axes:
            lines.append(f"{label:<8s}")
            lines.append(f"{v0:>24.15g}     ")
            lines.append(f"{step:>24.15g}     ")
            lines.append(f"{n:12d}")
        # Column-count and column-header line.
        lines.append(f"{9:12d}")
        cols = ["T(K)", "P(bar)", "rho,kg/m3", "vp,km/s", "vs,km/s",
                "cp,J/K/m3", "alpha,1/K", "Ks,bar", "Gs,bar"]
        if first_axis == "P":
            cols = ["P(bar)", "T(K)"] + cols[2:]
        lines.append("".join(f"{c:<15s}" for c in cols).rstrip())
        # Emit data rows.
        Ts = [T_K_min + i * T_K_step for i in range(T_n)]
        Ps = [P_bar_min + j * P_bar_step for j in range(P_n)]
        # Native writer order for first_axis="T" is: for each P, iterate T
        # (i.e., T is the fastest-changing column). This matches the study
        # exporter's T-first convention.
        if first_axis == "T":
            for P in Ps:
                for T in Ts:
                    row = [T, P, rho_kgm3, VP_kms, VS_kms,
                           Cp_JKm3, alpha_pK, KS_bar, GS_bar]
                    lines.append(" ".join(f"{v:>14.6g}" for v in row))
        else:
            for T in Ts:
                for P in Ps:
                    row = [P, T, rho_kgm3, VP_kms, VS_kms,
                           Cp_JKm3, alpha_pK, KS_bar, GS_bar]
                    lines.append(" ".join(f"{v:>14.6g}" for v in row))
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")

    _synth_counter = 0

    @classmethod
    def _synth_table_name(cls):
        cls._synth_counter += 1
        return f"PP_BAR_UNIT_TEST_SYNTH_{os.getpid()}_{cls._synth_counter}.tab"

    def _make_and_load(self, tmpdir):
        # Small regular grid that is still large enough for
        # RectBivariateSpline's default order (kx=ky=3).
        from PlanetProfile import _PERPLEXCACHE
        table_name = self._synth_table_name()
        cached_path = os.path.join(_PERPLEXCACHE, table_name)
        # Constant per-cell values so the interpolator returns them
        # exactly at any interior sample.
        params = dict(
            P_bar_min=1000.0, P_bar_step=100.0, P_n=8,     # 1000..1700 bar (100..170 MPa)
            T_K_min=800.0,    T_K_step=50.0,   T_n=8,     # 800..1150 K
            rho_kgm3=3300.0,
            VP_kms=8.0,
            VS_kms=4.5,
            Cp_JKm3=3.3e6,     # will load as Cp_JkgK = 3.3e6 / 3300 = 1000
            alpha_pK=3.0e-5,
            KS_bar=1.0e6,      # exact bar values; expected 100 GPa physical
            GS_bar=5.0e5,      # exact bar values; expected  50 GPa physical
        )
        self._write_synthetic_table(cached_path, first_axis="T", **params)
        from PlanetProfile.Thermodynamics.InnerEOS import PerplexEOSStruct
        eos = PerplexEOSStruct(
            table_name, EOSinterpMethod="nearest", nHeaders=13,
            EXTRAP=False, phiTop_frac=0, Pclosure_MPa=350, phiMin_frac=None,
        )
        return eos, cached_path

    def _cleanup(self, cached_path):
        try:
            os.remove(cached_path)
        except FileNotFoundError:
            pass

    @unittest.skipUnless(_HAS_NUMPY and _HAS_SCIPY,
                         "numpy and scipy are required for the loader test.")
    def test_perplex_native_loader_uses_physical_conversions(self):
        with tempfile.TemporaryDirectory(prefix="pp_bar_unit_test_") as tmp:
            # Run from a temp cwd so any UserConfigs/ that PP creates lands in /tmp.
            prev_cwd = os.getcwd()
            os.chdir(tmp)
            try:
                eos, cached_path = self._make_and_load(tmp)
                try:
                    # Domain: 1000..1700 bar => 100..170 MPa physical.
                    self.assertAlmostEqual(eos.Pmin, 100.0, places=9)
                    self.assertAlmostEqual(eos.Pmax, 170.0, places=9)
                    self.assertAlmostEqual(eos.Tmin, 800.0, places=9)
                    self.assertAlmostEqual(eos.Tmax, 1150.0, places=9)
                    # Interior spot values at a mid-grid point.
                    Pq, Tq = 130.0, 975.0
                    self.assertAlmostEqual(
                        float(eos.fn_rho_kgm3(Pq, Tq, grid=False)), 3300.0, places=6)
                    # Ks: 1e6 bar * 1e-4 GPa/bar = 100 GPa
                    Vp, Vs, Ks, Gs = eos.fn_Seismic(Pq, Tq, grid=False)
                    self.assertAlmostEqual(float(Ks), 100.0, places=6)
                    self.assertAlmostEqual(float(Gs),  50.0, places=6)
                    self.assertAlmostEqual(float(Vp),   8.0, places=6)
                    self.assertAlmostEqual(float(Vs),   4.5, places=6)
                    # Cp: 3.3e6 J/K/m3 / 3300 kg/m3 = 1000 J/kg/K
                    Cp = float(eos.fn_Cp_JkgK(Pq, Tq, grid=False))
                    self.assertAlmostEqual(Cp, 1000.0, places=6)
                    # alpha unchanged
                    alp = float(eos.fn_alpha_pK(Pq, Tq, grid=False))
                    self.assertAlmostEqual(alp, 3.0e-5, places=9)
                finally:
                    self._cleanup(cached_path)
            finally:
                os.chdir(prev_cwd)

    @unittest.skipUnless(_HAS_NUMPY and _HAS_SCIPY,
                         "numpy and scipy are required for the loader test.")
    def test_extrap_false_clamps_and_does_not_mutate_caller(self):
        with tempfile.TemporaryDirectory(prefix="pp_bar_unit_test_") as tmp:
            prev_cwd = os.getcwd()
            os.chdir(tmp)
            try:
                eos, cached_path = self._make_and_load(tmp)
                try:
                    P_query = np.array([50.0, 100.0, 200.0])  # 50 below Pmin
                    T_query = np.array([800.0, 800.0, 800.0])
                    P_before = P_query.copy()
                    T_before = T_query.copy()
                    rho_50   = float(eos.fn_rho_kgm3(np.array([50.0]),
                                                    np.array([800.0]),
                                                    grid=False))
                    rho_pmin = float(eos.fn_rho_kgm3(np.array([eos.Pmin]),
                                                    np.array([800.0]),
                                                    grid=False))
                    self.assertAlmostEqual(rho_50, rho_pmin, places=9,
                        msg="EXTRAP=False should silently clamp P<Pmin to Pmin.")
                    _ = eos.fn_rho_kgm3(P_query, T_query, grid=False)
                    self.assertTrue(np.array_equal(P_query, P_before),
                        msg="Caller's P array must not be mutated.")
                    self.assertTrue(np.array_equal(T_query, T_before),
                        msg="Caller's T array must not be mutated.")
                finally:
                    self._cleanup(cached_path)
            finally:
                os.chdir(prev_cwd)


class TestElasticIdentityConsistency(unittest.TestCase):
    """Phase-4 test 4: source vs post-loader elastic identities.

    Uses a synthetic table whose (rho, K_S, G_S, V_P, V_S) are
    physically self-consistent — i.e., the values were chosen so that
    V_S = sqrt(G_S / rho) and V_P = sqrt((K_S + 4 G_S / 3) / rho)
    hold within source precision. The test asserts that after PP loads
    them, the identities still hold (i.e., the loader did not scale K_S
    or G_S away from the source velocities).
    """

    @unittest.skipUnless(_HAS_NUMPY and _HAS_SCIPY,
                         "numpy and scipy are required.")
    def test_post_loader_elastic_identities(self):
        from PlanetProfile import _PERPLEXCACHE
        from PlanetProfile.Thermodynamics.InnerEOS import PerplexEOSStruct

        # Build a consistent source cell:
        rho = 3300.0                           # kg/m^3
        Gs_GPa = 50.0                          # target physical
        Ks_GPa = 120.0                         # target physical
        Vs_kms = np.sqrt(Gs_GPa * 1e9 / rho) / 1000.0
        Vp_kms = np.sqrt((Ks_GPa + 4.0 / 3.0 * Gs_GPa) * 1e9 / rho) / 1000.0
        Ks_bar = Ks_GPa / 1.0e-4               # physical inverse: GPa -> bar
        Gs_bar = Gs_GPa / 1.0e-4               # physical inverse: GPa -> bar

        # Sanity — the source values themselves must satisfy the identity to
        # the precision of double-precision arithmetic. Test A of § 6a.
        Vs_inferred_source = np.sqrt(Gs_GPa * 1e9 / rho)
        Vp_inferred_source = np.sqrt((Ks_GPa + 4.0/3.0*Gs_GPa) * 1e9 / rho)
        self.assertAlmostEqual(Vs_kms * 1000, Vs_inferred_source, places=6)
        self.assertAlmostEqual(Vp_kms * 1000, Vp_inferred_source, places=6)

        table_name = f"PP_BAR_UNIT_TEST_ELASTIC_{os.getpid()}.tab"
        cached_path = os.path.join(_PERPLEXCACHE, table_name)
        # Reuse the writer from TestPerpleXNativeLoader.
        writer = TestPerpleXNativeLoader()
        writer._write_synthetic_table(
            cached_path,
            P_bar_min=1000.0, P_bar_step=100.0, P_n=8,
            T_K_min=800.0,    T_K_step=50.0,   T_n=8,
            rho_kgm3=rho, VP_kms=Vp_kms, VS_kms=Vs_kms,
            Cp_JKm3=3.3e6, alpha_pK=3.0e-5,
            KS_bar=Ks_bar, GS_bar=Gs_bar,
            first_axis="T",
        )

        prev_cwd = os.getcwd()
        with tempfile.TemporaryDirectory(prefix="pp_bar_unit_test_elastic_") as tmp:
            os.chdir(tmp)
            try:
                eos = PerplexEOSStruct(
                    table_name, EOSinterpMethod="nearest", nHeaders=13,
                    EXTRAP=False, phiTop_frac=0, Pclosure_MPa=350, phiMin_frac=None,
                )
                # Test B of § 6a: identities on PP-loaded values.
                Vp_pp, Vs_pp, Ks_pp, Gs_pp = eos.fn_Seismic(130.0, 975.0, grid=False)
                rho_pp = float(eos.fn_rho_kgm3(130.0, 975.0, grid=False))
                Vp_pp = float(Vp_pp); Vs_pp = float(Vs_pp)
                Ks_pp = float(Ks_pp); Gs_pp = float(Gs_pp)

                # Loaded moduli must match physical inputs.
                self.assertAlmostEqual(Ks_pp, Ks_GPa, places=6)
                self.assertAlmostEqual(Gs_pp, Gs_GPa, places=6)

                Vs_from_loaded = np.sqrt(Gs_pp * 1e9 / rho_pp) / 1000.0
                Vp_from_loaded = np.sqrt((Ks_pp + 4.0/3.0*Gs_pp) * 1e9 / rho_pp) / 1000.0
                # After the correction, |inferred − loaded| must be at
                # source precision (well below 1 m/s).
                self.assertAlmostEqual(Vs_from_loaded, Vs_pp, places=5)
                self.assertAlmostEqual(Vp_from_loaded, Vp_pp, places=5)
            finally:
                try:
                    os.remove(cached_path)
                except FileNotFoundError:
                    pass
                os.chdir(prev_cwd)


class TestPREMPressureConversion(unittest.TestCase):
    """Phase-4 test 6a: PREM lookup uses corrected physical conversion.

    The PREM table in PlanetProfile/Thermodynamics/EOStables/PREMtable.txt
    stores pressures in kbar (Dziewonski & Anderson 1981). The loader
    at InnerEOS.py:331 computes
        PPREM_MPa = PPREM_kbar * 1e3 * Constants.bar2MPa
    which must equal PPREM_kbar * 100 exactly under the correction.
    """

    def test_kbar_to_MPa_via_bar2MPa(self):
        from PlanetProfile.Utilities.defineStructs import Constants
        # 1 kbar = 1000 bar = 100 MPa exactly.
        self.assertAlmostEqual(1.0 * 1e3 * Constants.bar2MPa, 100.0, places=12)
        self.assertAlmostEqual(3.370 * 1e3 * Constants.bar2MPa, 337.0, places=9)


class TestPlottingMultiplier(unittest.TestCase):
    """Phase-4 test 6b: PmultHydro (a plot-axis multiplier used to convert
    Planet.P_MPa to true bar for readable plot labels) must equal 10.0
    after the correction so that P_MPa * PmultHydro yields the same number
    of physical bar.
    """

    def test_PmultHydro_equals_10(self):
        from PlanetProfile.Utilities.defineStructs import Constants
        # PmultHydro is defined by 1/bar2MPa inside the plotting labels code path.
        self.assertAlmostEqual(1.0 / Constants.bar2MPa, 10.0, places=12)


if __name__ == "__main__":
    unittest.main()
