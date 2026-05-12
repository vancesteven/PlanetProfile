import unittest

import numpy as np

from PlanetProfile.Thermodynamics.HydroEOS import (
    GetIceArrheniusViscosityKwargs,
    ViscIceArrhenius_Pas,
)
from PlanetProfile.Thermodynamics.LayerPropagators import (
    _ConvectionDeschampsSotinHPIceDiagnostic,
    _FixedPhaseEOS,
    _SetHPIceDiagnosticFields,
)
from PlanetProfile.Thermodynamics.ThermalProfiles.ThermalProfiles import (
    ConvectionKalousova2018,
)
from PlanetProfile.Utilities.defineStructs import Constants, OceanSubstruct


class _Namespace:
    pass


def _planet_with_arrhenius_flags():
    planet = _Namespace()
    planet.Do = _Namespace()
    for flag in (
        "ARRHENIUS_VISCOSITY",
        "ARRHENIUS_VISCOSITY_Ih",
        "ARRHENIUS_VISCOSITY_III",
        "ARRHENIUS_VISCOSITY_V",
        "ARRHENIUS_VISCOSITY_VI",
    ):
        setattr(planet.Do, flag, False)

    planet.Ocean = _Namespace()
    planet.Ocean.Eact_kJmol = {"Ih": np.nan, "III": np.nan, "V": np.nan, "VI": np.nan}

    planet.Bulk = _Namespace()
    planet.Bulk.TbIII_K = 253.0
    planet.Bulk.TbV_K = 263.0
    return planet


class _SimpleIceEOS:
    phaseID = 5
    POROUS = False

    def fn_rho_kgm3(self, P_MPa, T_K):
        return 1300.0

    def fn_Cp_JkgK(self, P_MPa, T_K):
        return 2000.0

    def fn_alpha_pK(self, P_MPa, T_K):
        return 1.0e-4

    def fn_kTherm_WmK(self, P_MPa, T_K):
        return 3.0


class Test01ArrheniusViscosity(unittest.TestCase):

    def test_arrhenius_kwargs_empty_when_flags_disabled(self):
        planet = _planet_with_arrhenius_flags()
        self.assertEqual(GetIceArrheniusViscosityKwargs(planet, "Ih"), {})
        self.assertEqual(GetIceArrheniusViscosityKwargs(planet, "V"), {})

    def test_arrhenius_kwargs_return_phase_parameters_when_enabled(self):
        planet = _planet_with_arrhenius_flags()
        planet.Do.ARRHENIUS_VISCOSITY_Ih = True

        kwargs = GetIceArrheniusViscosityKwargs(planet, "Ih")

        self.assertTrue(kwargs["ARRHENIUS_VISCOSITY"])
        self.assertGreater(kwargs["etaMelt_Pas"], 0.0)
        self.assertGreater(kwargs["Eact_Jmol"], 0.0)
        self.assertGreater(kwargs["Tmelt_K"], 0.0)

    def test_hp_arrhenius_phase_flags_do_not_crash(self):
        for phase in ("III", "V", "VI"):
            with self.subTest(phase=phase):
                planet = _planet_with_arrhenius_flags()
                setattr(planet.Do, f"ARRHENIUS_VISCOSITY_{phase}", True)
                kwargs = GetIceArrheniusViscosityKwargs(planet, phase)
                self.assertTrue(kwargs["ARRHENIUS_VISCOSITY"])
                self.assertGreater(kwargs["etaMelt_Pas"], 0.0)

    def test_arrhenius_viscosity_decreases_toward_melt(self):
        viscosity = ViscIceArrhenius_Pas(
            etaMelt_Pas=1.0e13,
            Eact_Jmol=60.0e3,
            Tmelt_K=273.16,
        )

        eta_cold, eta_warm = viscosity(0.0, np.array([250.0, 273.16]))

        self.assertGreater(eta_cold, eta_warm)
        np.testing.assert_allclose(eta_warm, 1.0e13)

    def test_arrhenius_update_convection_viscosity_is_noop(self):
        viscosity = ViscIceArrhenius_Pas(
            etaMelt_Pas=1.0e13,
            Eact_Jmol=60.0e3,
            Tmelt_K=273.16,
        )
        before = viscosity(0.0, np.array([250.0, 260.0, 270.0]))

        viscosity.updateConvectionViscosity(1.0e99, 100.0)
        after = viscosity(0.0, np.array([250.0, 260.0, 270.0]))

        np.testing.assert_allclose(after, before)


class Test03HPIceConvectionDiagnostics(unittest.TestCase):

    def test_fixed_phase_eos_reports_existing_hp_phase(self):
        ocean_eos = _Namespace()
        ocean_eos.Tmin = 200.0
        phase_eos = _FixedPhaseEOS(ocean_eos, 5)

        phases = phase_eos.fn_phase(600.0, np.array([250.0, 260.0]))

        self.assertEqual(phase_eos.Tmin, 200.0)
        np.testing.assert_array_equal(phases, np.array([5, 5]))

    def test_hp_diagnostic_writer_keeps_top_level_and_dict_in_sync(self):
        planet = _Namespace()
        planet.HPIceDiagnostics = {}

        _SetHPIceDiagnosticFields(
            planet,
            "V",
            status="computed",
            thickness_m=2.0e5,
            Tconv_K=264.0,
            etaConv_Pas=1.0e15,
            etaMelt_Pas=5.0e14,
            eLid_m=1.0e4,
            Dconv_m=1.8e5,
            deltaTBL_m=1.0e4,
            Ra=1.0e6,
            RaCrit=1.0e4,
            Qbot_W=1.0e12,
            Pmid_MPa=600.0,
            method="Deschamps and Sotin (2001)",
        )

        self.assertEqual(planet.TconvV_K, planet.HPIceDiagnostics["V"]["Tconv_K"])
        self.assertEqual(planet.etaConvV_Pas, planet.HPIceDiagnostics["V"]["etaConv_Pas"])
        self.assertEqual(planet.DconvV_m, planet.HPIceDiagnostics["V"]["Dconv_m"])
        self.assertEqual(planet.RaConvectV, planet.HPIceDiagnostics["V"]["RaConvect"])
        self.assertEqual(planet.HPIceDiagnostics["V"]["status"], "computed")

    def test_ds2001_hp_diagnostic_returns_finite_values(self):
        diagnostics = _ConvectionDeschampsSotinHPIceDiagnostic(
            Ttop_K=250.0,
            rTop_m=2.4e6,
            kTop_WmK=3.0,
            Tb_K=270.0,
            zb_m=2.0e5,
            gtop_ms2=1.0,
            Pmid_MPa=600.0,
            iceEOS=_SimpleIceEOS(),
            phaseID=5,
            EQUIL_Q=True,
            Eact_kJmol={"V": np.nan},
        )

        Tconv_K, etaConv_Pas, eLid_m, Dconv_m, deltaTBL_m, Qbot_W, Ra, RaCrit = diagnostics
        self.assertTrue(np.all(np.isfinite(diagnostics)))
        self.assertGreater(Tconv_K, 250.0)
        self.assertGreater(etaConv_Pas, 0.0)
        self.assertGreaterEqual(eLid_m, 0.0)
        self.assertGreaterEqual(Dconv_m, 0.0)
        self.assertGreaterEqual(deltaTBL_m, 0.0)
        self.assertGreater(Qbot_W, 0.0)
        self.assertGreater(Ra, RaCrit)


class Test04KalousovaHPIceDiagnostics(unittest.TestCase):

    def test_kalousova_subcritical_branch_has_no_melt_fraction(self):
        (
            Tconv_K,
            etaConv_Pas,
            eLid_m,
            Dconv_m,
            deltaTBL_m,
            Qbot_W,
            Ra,
            RaCrit,
        ) = ConvectionKalousova2018(
            Ttop_K=250.0,
            rTop_m=2.4e6,
            kTop_WmK=3.0,
            Tb_K=270.0,
            zb_m=5.0e5,
            gtop_ms2=1.0,
            Pmid_MPa=1000.0,
            oceanEOS=None,
            iceEOS=_SimpleIceEOS(),
            phaseBot=5,
            EQUIL_Q=True,
            Eact_kJmol={},
            qBot_Wm2=0.1,
            etaMelt_Pas=1.0e14,
        )

        self.assertEqual(Tconv_K, 250.0)
        self.assertEqual(etaConv_Pas, 1.0e14)
        self.assertEqual(eLid_m, 5.0e5)
        self.assertEqual(Dconv_m, 0.0)
        self.assertEqual(deltaTBL_m, 0.0)
        self.assertGreater(Qbot_W, 0.0)
        self.assertLess(Ra, RaCrit)

    def test_kalousova_ice_v_viscosity_is_local_to_kalousova(self):
        ocean = OceanSubstruct()
        ice_v_phase_id = 5

        self.assertEqual(ocean.etaMeltKalousova_Pas["V"], 2.8e14)
        self.assertEqual(Constants.etaMelt_Pas[ice_v_phase_id], 5e14)

    def test_kalousova_supercritical_branch_returns_finite_diagnostics(self):
        diagnostics = ConvectionKalousova2018(
            Ttop_K=250.0,
            rTop_m=2.4e6,
            kTop_WmK=3.0,
            Tb_K=270.0,
            zb_m=5.0e5,
            gtop_ms2=1.0,
            Pmid_MPa=1000.0,
            oceanEOS=None,
            iceEOS=_SimpleIceEOS(),
            phaseBot=5,
            EQUIL_Q=True,
            Eact_kJmol={},
            qBot_Wm2=1.0e-3,
            etaMelt_Pas=1.0e14,
        )

        Tconv_K, etaConv_Pas, eLid_m, Dconv_m, deltaTBL_m, Qbot_W, Ra, RaCrit = diagnostics
        self.assertEqual(Tconv_K, 250.0)
        self.assertEqual(etaConv_Pas, 1.0e14)
        self.assertGreaterEqual(eLid_m, 0.0)
        self.assertGreaterEqual(Dconv_m, 0.0)
        self.assertGreaterEqual(deltaTBL_m, 0.0)
        self.assertGreater(Qbot_W, 0.0)
        self.assertGreater(Ra, RaCrit)
        self.assertTrue(np.all(np.isfinite(diagnostics)))

    def test_kalousova_scaling_formulas_match_genai_validation(self):
        diagnostics = ConvectionKalousova2018(
            Ttop_K=250.0,
            rTop_m=2.4e6,
            kTop_WmK=3.0,
            Tb_K=270.0,
            zb_m=5.0e5,
            gtop_ms2=1.0,
            Pmid_MPa=1000.0,
            oceanEOS=None,
            iceEOS=_SimpleIceEOS(),
            phaseBot=5,
            EQUIL_Q=True,
            Eact_kJmol={},
            qBot_Wm2=1.0e-3,
            etaMelt_Pas=1.0e14,
        )

        _, _, eLid_m, _, deltaTBL_m, _, Ra, RaCrit = diagnostics
        qs_mWm2 = 1.0
        expected_RaCrit = 19.965e3 * qs_mWm2**3.690
        expected_eLid_m = ((0.145e-3 * qs_mWm2 + 0.015) * (1.0e14**0.21)) * 1.0e3
        expected_deltaTBL_m = (2.746 * Ra**(-0.271)) * 5.0e5

        np.testing.assert_allclose(RaCrit, expected_RaCrit)
        np.testing.assert_allclose(eLid_m, expected_eLid_m)
        np.testing.assert_allclose(deltaTBL_m, expected_deltaTBL_m)


if __name__ == "__main__":
    unittest.main()
