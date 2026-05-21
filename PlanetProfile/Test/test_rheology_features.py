import unittest
from unittest import mock

import numpy as np

from PlanetProfile.Thermodynamics.HydroEOS import (
    GetIceArrheniusViscosityKwargs,
    GetIceEOS,
    ViscIceArrhenius_Pas,
)
from PlanetProfile.Thermodynamics.LayerPropagators import (
    ApplyActiveIceVIRollbackPolicy,
    BuildActiveIceVIProductionCandidateCopy,
    BuildActiveIceVIThermalUpdateCandidate,
    EvaluateActiveIceVICandidateResiduals,
    EvaluatePosthocIceVIProductionCandidate,
    _ConvectionDeschampsSotinHPIceDiagnostic,
    _FixedPhaseEOS,
    _GetIceVIMeltCurveCandidateChecks,
    _SetHPIceDiagnosticFields,
    HPIceConvectionDiagnostics,
)
from PlanetProfile.Thermodynamics.ThermalProfiles.ThermalProfiles import (
    ConvectionKalousova2018,
)
from PlanetProfile.Utilities.defineStructs import (
    Constants,
    DoSubstruct,
    EOSlist,
    HPIcePhaseState,
    OceanSubstruct,
    ResetMutableModelState,
    ResolveHPIceConvectionModel,
    _EvaluateIceVIProductionAcceptance,
)


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


def _reset_ice_eos_state(eta_ice=None, tvisc_ice=None):
    ResetMutableModelState(reset_eos=True)
    if eta_ice is not None:
        Constants.etaIce_Pas[:] = list(eta_ice)
    if tvisc_ice is not None:
        Constants.TviscIce_K[:] = list(tvisc_ice)


class _SimpleIceEOS:
    POROUS = False

    def __init__(self, phaseID=5):
        self.phaseID = phaseID

    def fn_rho_kgm3(self, P_MPa, T_K):
        return 1300.0

    def fn_Cp_JkgK(self, P_MPa, T_K):
        return 2000.0

    def fn_alpha_pK(self, P_MPa, T_K):
        return 1.0e-4

    def fn_kTherm_WmK(self, P_MPa, T_K):
        return 3.0


class _SyntheticIceVIMeltEOS:
    Pmin = 700.0
    Pmax = 1300.0
    Tmin = 250.0
    Tmax = 340.0

    @staticmethod
    def tmelt_K(P_MPa):
        return 292.0 + 0.01 * (P_MPa - 800.0)

    def fn_phase(self, P_MPa, T_K):
        phase = np.where(np.asarray(T_K) < self.tmelt_K(P_MPa), 6, 0)
        return phase


class _SyntheticWrongPhaseEOS(_SyntheticIceVIMeltEOS):

    def fn_phase(self, P_MPa, T_K):
        return np.zeros_like(np.asarray(T_K), dtype=int)


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

    def test_arrhenius_update_convection_viscosity_preserves_arrhenius_curve(self):
        eta_ice_original = list(Constants.etaIce_Pas)
        tvisc_ice_original = list(Constants.TviscIce_K)
        viscosity = ViscIceArrhenius_Pas(
            etaMelt_Pas=1.0e13,
            Eact_Jmol=60.0e3,
            Tmelt_K=273.16,
        )
        before = viscosity(0.0, np.array([250.0, 260.0, 270.0]))

        try:
            viscosity.updateConvectionViscosity(1.0e16, 260.0)
            after = viscosity(0.0, np.array([250.0, 260.0, 270.0]))

            np.testing.assert_allclose(after, before)
            self.assertEqual(Constants.etaIce_Pas[-1], 1.0e16)
            self.assertEqual(Constants.TviscIce_K[-1], 260.0)
        finally:
            _reset_ice_eos_state(eta_ice_original, tvisc_ice_original)

    def test_ih_arrhenius_changes_only_ih_in_synthetic_hp_profile(self):
        eta_ice_original = list(Constants.etaIce_Pas)
        tvisc_ice_original = list(Constants.TviscIce_K)
        phase = np.array([1, 1, 5, 6])
        P_MPa = np.array([5.0, 50.0, 500.0, 1000.0])
        T_K = np.array([220.0, 260.0, 270.0, 290.0])
        eta_conv_Pas = 6.0e14
        Tconv_K = 255.0

        def synthetic_eta(ih_arrhenius):
            _reset_ice_eos_state(eta_ice_original, tvisc_ice_original)
            ih_kwargs = {}
            if ih_arrhenius:
                ih_kwargs = {
                    "ARRHENIUS_VISCOSITY": True,
                    "etaMelt_Pas": Constants.etaMelt_Pas[1],
                    "Eact_Jmol": Constants.Eact_kJmol[1] * 1.0e3,
                    "Tmelt_K": Constants.T0,
                }

            ih_eos = GetIceEOS(
                np.linspace(0.1, 100.0, 5),
                np.linspace(200.0, 270.0, 5),
                "Ih",
                **ih_kwargs,
            )
            ih_eos.updateConvectionViscosity(eta_conv_Pas, Tconv_K)
            v_eos = GetIceEOS(
                np.linspace(350.0, 700.0, 5),
                np.linspace(260.0, 290.0, 5),
                "V",
            )
            vi_eos = GetIceEOS(
                np.linspace(700.0, 1500.0, 5),
                np.linspace(270.0, 320.0, 5),
                "VI",
            )

            eta_Pas = np.zeros_like(T_K, dtype=float)
            eta_Pas[phase == 1] = ih_eos.fn_eta_Pas(P_MPa[phase == 1], T_K[phase == 1])
            eta_Pas[phase == 5] = v_eos.fn_eta_Pas(P_MPa[phase == 5], T_K[phase == 5])
            eta_Pas[phase == 6] = vi_eos.fn_eta_Pas(P_MPa[phase == 6], T_K[phase == 6])
            return eta_Pas

        try:
            default_eta_Pas = synthetic_eta(ih_arrhenius=False)
            ih_arrhenius_eta_Pas = synthetic_eta(ih_arrhenius=True)

            self.assertTrue(np.any(default_eta_Pas[phase == 1] != ih_arrhenius_eta_Pas[phase == 1]))
            np.testing.assert_allclose(
                ih_arrhenius_eta_Pas[phase == 5],
                default_eta_Pas[phase == 5],
            )
            np.testing.assert_allclose(
                ih_arrhenius_eta_Pas[phase == 6],
                default_eta_Pas[phase == 6],
            )
        finally:
            _reset_ice_eos_state(eta_ice_original, tvisc_ice_original)

    def test_reset_mutable_model_state_restores_viscosity_defaults_and_clears_eos_cache(self):
        _reset_ice_eos_state()
        ih_eos = GetIceEOS(
            np.linspace(0.1, 100.0, 5),
            np.linspace(200.0, 270.0, 5),
            "Ih",
        )
        ih_eos.updateConvectionViscosity(2.0e16, 260.0)
        EOSlist.loaded["synthetic_run_local_cache"] = object()
        EOSlist.ranges["synthetic_run_local_cache"] = "changed"

        ResetMutableModelState(reset_eos=True)

        self.assertEqual(Constants.etaIce_Pas, [1.0e19, 1.0e15])
        self.assertEqual(Constants.TviscIce_K, [241])
        self.assertNotIn("synthetic_run_local_cache", EOSlist.loaded)
        self.assertEqual(set(EOSlist.loaded), {"CustomSolutionEOS", "ReaktoroDatabases"})
        self.assertEqual(EOSlist.ranges, {})

    def test_reset_preserves_within_run_uniform_viscosity_handoff(self):
        _reset_ice_eos_state()
        ih_eos = GetIceEOS(
            np.linspace(0.1, 100.0, 5),
            np.linspace(200.0, 270.0, 5),
            "Ih",
        )
        ih_eos.updateConvectionViscosity(4.0e16, 258.0)

        v_eos = GetIceEOS(
            np.linspace(350.0, 700.0, 5),
            np.linspace(250.0, 290.0, 5),
            "V",
        )

        self.assertEqual(Constants.etaIce_Pas[-1], 4.0e16)
        self.assertEqual(Constants.TviscIce_K[-1], 258.0)
        np.testing.assert_allclose(
            v_eos.fn_eta_Pas(np.array([500.0]), np.array([270.0])),
            np.array([4.0e16]),
        )


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


class Test05HPIceConvectionModelSelector(unittest.TestCase):

    def test_default_selector_is_none_and_diagnostics_disabled(self):
        do = DoSubstruct()

        self.assertEqual(do.HP_ICE_CONVECTION_MODEL, "none")
        self.assertFalse(do.HP_ICE_CONVECTION_DIAGNOSTICS)
        self.assertFalse(do.KALOUSOVA_CONVECTION)
        self.assertEqual(ResolveHPIceConvectionModel(do), "none")

    def test_legacy_hp_diagnostics_flag_maps_to_ds2001_diagnostic(self):
        do = DoSubstruct()
        do.HP_ICE_CONVECTION_DIAGNOSTICS = True

        self.assertEqual(ResolveHPIceConvectionModel(do), "DS2001_diagnostic")

    def test_legacy_kalousova_flag_maps_to_kalousova_diagnostic(self):
        do = DoSubstruct()
        do.KALOUSOVA_CONVECTION = True

        self.assertEqual(ResolveHPIceConvectionModel(do), "Kalousova2018_diagnostic")

    def test_explicit_selector_wins_over_legacy_flags(self):
        do = DoSubstruct()
        do.HP_ICE_CONVECTION_MODEL = "DS2001_diagnostic"
        do.KALOUSOVA_CONVECTION = True

        self.assertEqual(ResolveHPIceConvectionModel(do), "DS2001_diagnostic")

    def test_invalid_selector_raises_clear_error(self):
        do = DoSubstruct()
        do.HP_ICE_CONVECTION_MODEL = "not_a_model"

        with self.assertRaisesRegex(ValueError, "HP_ICE_CONVECTION_MODEL"):
            ResolveHPIceConvectionModel(do)

    def test_production_selector_requires_experimental_gate(self):
        do = DoSubstruct()
        do.HP_ICE_CONVECTION_MODEL = "Kalousova2018_production_experimental"

        with self.assertRaisesRegex(ValueError, "ALLOW_EXPERIMENTAL_HP_KALOUSOVA_PRODUCTION=True"):
            ResolveHPIceConvectionModel(do)

    def test_production_selector_with_gate_is_reserved_not_silent(self):
        do = DoSubstruct()
        do.HP_ICE_CONVECTION_MODEL = "Kalousova2018_production_experimental"
        do.ALLOW_EXPERIMENTAL_HP_KALOUSOVA_PRODUCTION = True

        self.assertEqual(ResolveHPIceConvectionModel(do), "Kalousova2018_production_experimental")


class Test06HPIcePhaseState(unittest.TestCase):

    def test_absent_phase_state_is_consistent(self):
        state = HPIcePhaseState("VI", phaseID=6, status="absent")

        self.assertFalse(state.present)
        self.assertEqual(state.validityStatus, "absent")
        self.assertEqual(state.skipReason, "absent")

    def test_valid_phase_state_stores_diagnostic_fields(self):
        state = HPIcePhaseState(
            "V",
            phaseID=5,
            present=True,
            status="computed",
            thickness_m=2.0e5,
            Ttop_K=260.0,
            Tbot_K=270.0,
            Tconv_K=264.0,
            etaConv_Pas=1.0e15,
            etaMelt_Pas=5.0e14,
            eLid_m=1.0e4,
            Dconv_m=1.8e5,
            deltaTBL_m=1.0e4,
            RaConvect=1.0e6,
            RaCrit=1.0e4,
            meltFraction=0.05,
            DO_HP_MELT=True,
        )

        self.assertEqual(state.validityStatus, "ok")
        self.assertTrue(state.DO_HP_MELT)
        self.assertEqual(state.as_dict()["etaConv_Pas"], 1.0e15)

    def test_diagnostic_writer_keeps_phase_state_synchronized(self):
        planet = _Namespace()
        planet.HPIceDiagnostics = {}

        _SetHPIceDiagnosticFields(
            planet,
            "V",
            status="computed",
            phaseID=5,
            iTop=4,
            iBot=8,
            thickness_m=2.0e5,
            Ttop_K=260.0,
            Tbot_K=270.0,
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
            meltFraction=0.05,
        )

        state = planet.HPIceDiagnostics["V"]
        self.assertEqual(planet.TconvV_K, state["Tconv_K"])
        self.assertEqual(planet.etaConvV_Pas, state["etaConv_Pas"])
        self.assertEqual(planet.meltFractionV, state["meltFraction"])
        self.assertEqual(state["phaseName"], "V")
        self.assertEqual(state["phaseID"], 5)
        self.assertTrue(state["present"])
        self.assertEqual(state["validityStatus"], "ok")
        self.assertTrue(state["DO_HP_MELT"])

    def test_negative_thickness_is_flagged(self):
        state = HPIcePhaseState(
            "III",
            phaseID=3,
            present=True,
            status="computed",
            thickness_m=-1.0,
            Tconv_K=250.0,
            etaConv_Pas=1.0e14,
            eLid_m=0.0,
            Dconv_m=0.0,
            deltaTBL_m=0.0,
            RaConvect=1.0e6,
            RaCrit=1.0e4,
        )

        self.assertEqual(state.validityStatus, "invalid_geometry")

    def test_nonfinite_computed_state_is_flagged(self):
        state = HPIcePhaseState(
            "V",
            phaseID=5,
            present=True,
            status="computed",
            thickness_m=2.0e5,
            Tconv_K=np.nan,
            etaConv_Pas=1.0e15,
            eLid_m=1.0e4,
            Dconv_m=1.8e5,
            deltaTBL_m=1.0e4,
            RaConvect=1.0e6,
            RaCrit=1.0e4,
        )

        self.assertEqual(state.validityStatus, "nonfinite")

    def test_melt_fraction_outside_bounds_is_flagged(self):
        state = HPIcePhaseState(
            "VI",
            phaseID=6,
            present=True,
            status="computed",
            thickness_m=2.0e5,
            Tconv_K=286.0,
            etaConv_Pas=5.0e14,
            eLid_m=1.0e4,
            Dconv_m=1.8e5,
            deltaTBL_m=1.0e4,
            RaConvect=1.0e6,
            RaCrit=1.0e4,
            meltFraction=1.5,
        )

        self.assertEqual(state.validityStatus, "invalid_melt_fraction")

    def test_selector_behavior_is_unchanged_by_phase_state_scaffold(self):
        do = DoSubstruct()
        do.HP_ICE_CONVECTION_DIAGNOSTICS = True

        self.assertEqual(ResolveHPIceConvectionModel(do), "DS2001_diagnostic")


class Test07HPIceHeatBookkeeping(unittest.TestCase):

    def test_valid_state_defaults_to_conserved_heat_when_no_sources_active(self):
        state = HPIcePhaseState(
            "VI",
            phaseID=6,
            present=True,
            status="computed",
            rTop_m=2.3e6,
            rBot_m=2.1e6,
            thickness_m=2.0e5,
            Ttop_K=280.0,
            Tbot_K=300.0,
            Tconv_K=286.0,
            etaConv_Pas=5.0e14,
            eLid_m=1.0e4,
            Dconv_m=1.8e5,
            deltaTBL_m=1.0e4,
            RaConvect=1.0e6,
            RaCrit=1.0e4,
            Q_in_W=1.0e12,
        )

        self.assertEqual(state.Q_out_W, 1.0e12)
        self.assertEqual(state.Q_internal_W, 0.0)
        self.assertEqual(state.Q_latent_W, 0.0)
        self.assertEqual(state.energyResidual_W, 0.0)
        self.assertEqual(state.energyResidual_frac, 0.0)
        self.assertEqual(state.heatBookkeepingStatus, "ok")

    def test_zero_contrast_preserves_through_heat_diagnostic(self):
        state = HPIcePhaseState(
            "V",
            phaseID=5,
            present=True,
            status="computed",
            rTop_m=2.3e6,
            rBot_m=2.1e6,
            thickness_m=2.0e5,
            Ttop_K=270.0,
            Tbot_K=270.0,
            Tconv_K=270.0,
            etaConv_Pas=5.0e14,
            eLid_m=0.0,
            Dconv_m=0.0,
            deltaTBL_m=0.0,
            RaConvect=0.0,
            RaCrit=np.inf,
            Q_in_W=2.0e12,
        )

        self.assertEqual(state.Q_out_W, 2.0e12)
        self.assertEqual(state.energyResidual_W, 0.0)
        self.assertEqual(state.validityStatus, "zero_contrast")

    def test_absent_phase_heat_fields_are_unavailable(self):
        state = HPIcePhaseState("III", phaseID=3, status="absent")

        self.assertTrue(np.isnan(state.Q_in_W))
        self.assertTrue(np.isnan(state.Q_out_W))
        self.assertEqual(state.heatBookkeepingStatus, "absent")

    def test_too_thin_phase_does_not_invent_heat_values(self):
        state = HPIcePhaseState(
            "III",
            phaseID=3,
            present=True,
            status="too thin",
            thickness_m=500.0,
        )

        self.assertTrue(np.isnan(state.Q_in_W))
        self.assertTrue(np.isnan(state.Q_out_W))
        self.assertEqual(state.validityStatus, "too_thin")
        self.assertEqual(state.heatBookkeepingStatus, "unavailable")

    def test_energy_residual_sign_convention(self):
        state = HPIcePhaseState(
            "VI",
            phaseID=6,
            present=True,
            status="computed",
            thickness_m=2.0e5,
            Ttop_K=280.0,
            Tbot_K=300.0,
            Tconv_K=286.0,
            etaConv_Pas=5.0e14,
            eLid_m=1.0e4,
            Dconv_m=1.8e5,
            deltaTBL_m=1.0e4,
            RaConvect=1.0e6,
            RaCrit=1.0e4,
            Q_in_W=100.0,
            Q_out_W=110.0,
            Q_internal_W=5.0,
            Q_latent_W=2.0,
        )

        self.assertEqual(state.energyResidual_W, 7.0)
        self.assertEqual(state.energyResidual_frac, 7.0 / 110.0)

    def test_heat_bookkeeping_does_not_mutate_profile_arrays(self):
        planet = _Namespace()
        planet.HPIceDiagnostics = {}
        planet.T_K = np.array([250.0, 260.0, 270.0])
        planet.P_MPa = np.array([100.0, 200.0, 300.0])
        planet.phase = np.array([0, 5, 5])
        planet.rho_kgm3 = np.array([1000.0, 1200.0, 1300.0])
        before = {
            "T_K": planet.T_K.copy(),
            "P_MPa": planet.P_MPa.copy(),
            "phase": planet.phase.copy(),
            "rho_kgm3": planet.rho_kgm3.copy(),
        }

        _SetHPIceDiagnosticFields(
            planet,
            "V",
            status="computed",
            phaseID=5,
            thickness_m=2.0e5,
            Ttop_K=260.0,
            Tbot_K=270.0,
            Tconv_K=264.0,
            etaConv_Pas=1.0e15,
            etaMelt_Pas=5.0e14,
            eLid_m=1.0e4,
            Dconv_m=1.8e5,
            deltaTBL_m=1.0e4,
            Ra=1.0e6,
            RaCrit=1.0e4,
            Q_in_W=1.0e12,
            Q_out_W=1.0e12,
        )

        np.testing.assert_array_equal(planet.T_K, before["T_K"])
        np.testing.assert_array_equal(planet.P_MPa, before["P_MPa"])
        np.testing.assert_array_equal(planet.phase, before["phase"])
        np.testing.assert_array_equal(planet.rho_kgm3, before["rho_kgm3"])
        self.assertEqual(planet.HPIceDiagnostics["V"]["energyResidual_W"], 0.0)


class Test08HPIceProductionDryRun(unittest.TestCase):

    def _synthetic_hp_planet(self, qBot_Wm2=1.0e-3, vi_deltaT_K=20.0, vi_thickness_m=2.0e5):
        planet = _Namespace()
        planet.Do = DoSubstruct()
        planet.Do.VALID = True
        planet.Do.NO_H2O = False
        planet.Do.NO_OCEAN = False
        planet.Do.NO_OCEAN_EXCEPT_INNER_ICES = False
        planet.Do.BOTTOM_ICEV = False
        planet.Do.BOTTOM_ICEIII = False
        planet.Do.NO_ICE_CONVECTION = False
        planet.Do.EQUIL_Q = True
        planet.Do.HP_ICE_CONVECTION_MODEL = "Kalousova2018_production_experimental"
        planet.Do.ALLOW_EXPERIMENTAL_HP_KALOUSOVA_PRODUCTION = True

        planet.Steps = _Namespace()
        planet.Steps.nSurfIce = 0
        planet.Steps.nOceanMax = 6

        z_m = np.array([0.0, 2.0e5, 2.5e5, 4.5e5, 5.0e5, 5.0e5 + vi_thickness_m])
        r_m = 3.0e6 - z_m
        vi_top_T_K = 280.0
        planet.T_K = np.array([250.0, 270.0, 260.0, 280.0, vi_top_T_K, vi_top_T_K + vi_deltaT_K])
        planet.P_MPa = np.array([100.0, 200.0, 300.0, 500.0, 800.0, 1100.0])
        planet.phase = np.array([3, 3, 5, 5, 6, 6])
        planet.r_m = r_m
        planet.z_m = z_m
        planet.kTherm_WmK = np.zeros(6) + 3.0
        planet.g_ms2 = np.zeros(6) + 1.0
        planet.rho_kgm3 = np.array([1200.0, 1210.0, 1220.0, 1230.0, 1300.0, 1310.0])
        planet.Cp_JkgK = np.zeros(6) + 2200.0
        planet.alpha_pK = np.zeros(6) + 1.0e-4
        planet.qSurf_Wm2 = 0.02
        planet.qCon_Wm2 = 0.02
        planet.Mtot_kg = 1.0e23
        planet.CMR2mean = 0.33
        planet.eta_Pas = np.zeros(6) + 1.0e15

        planet.Ocean = _Namespace()
        planet.Ocean.QfromMantle_W = qBot_Wm2 * 4.0 * np.pi * r_m[-1]**2
        planet.Ocean.iceEOS = {
            "III": _SimpleIceEOS(3),
            "V": _SimpleIceEOS(5),
            "VI": _SimpleIceEOS(6),
        }
        planet.Ocean.EOS = _Namespace()
        planet.Ocean.Eact_kJmol = {"III": np.nan, "V": np.nan, "VI": np.nan}
        planet.Ocean.HtidalIce_Wm3 = 0.0
        planet.Ocean.etaMeltKalousova_Pas = {"III": 5.0e12, "V": 2.8e14, "VI": 5.0e14}
        planet.Ocean.phiPercolationKalousova_frac = 0.05
        planet.HPIceDiagnostics = {}
        return planet

    def _profile_snapshot(self, planet):
        return {
            "T_K": planet.T_K.copy(),
            "P_MPa": planet.P_MPa.copy(),
            "phase": planet.phase.copy(),
            "rho_kgm3": planet.rho_kgm3.copy(),
            "qSurf_Wm2": planet.qSurf_Wm2,
            "qCon_Wm2": planet.qCon_Wm2,
            "Mtot_kg": planet.Mtot_kg,
            "CMR2mean": planet.CMR2mean,
            "eta_Pas": planet.eta_Pas.copy(),
        }

    def test_production_selector_without_gate_raises(self):
        do = DoSubstruct()
        do.HP_ICE_CONVECTION_MODEL = "Kalousova2018_production_experimental"

        with self.assertRaisesRegex(ValueError, "ALLOW_EXPERIMENTAL_HP_KALOUSOVA_PRODUCTION=True"):
            ResolveHPIceConvectionModel(do)

    def test_gated_dry_run_records_ice_vi_rejection_without_profile_mutation(self):
        planet = self._synthetic_hp_planet(qBot_Wm2=1.0e-3)
        before = self._profile_snapshot(planet)

        HPIceConvectionDiagnostics(planet, _Namespace(), ResolveHPIceConvectionModel(planet))

        state = planet.HPIceDiagnostics["VI"]
        self.assertFalse(state["productionCandidate"])
        self.assertEqual(state["productionMode"], "Kalousova2018_production_experimental")
        self.assertEqual(state["candidateStatus"], "missing_melt_curve_rejected")
        self.assertFalse(state["updateAccepted"])
        self.assertEqual(state["candidateReason"], "missing_TmeltVI_P")
        self.assertEqual(state["massResidual_kg"], 0.0)
        self.assertEqual(state["massResidual_frac"], 0.0)
        for key, value in before.items():
            if isinstance(value, np.ndarray):
                np.testing.assert_array_equal(getattr(planet, key), value)
            else:
                self.assertEqual(getattr(planet, key), value)

    def test_ice_vi_subcritical_case_is_rejected(self):
        state = HPIcePhaseState(
            "VI",
            phaseID=6,
            present=True,
            status="computed",
            productionMode="Kalousova2018_production_experimental",
            rTop_m=2.3e6,
            rBot_m=2.1e6,
            zTop_m=1.0e5,
            zBot_m=3.0e5,
            thickness_m=2.0e5,
            Ttop_K=280.0,
            Tbot_K=300.0,
            rho_kgm3=1300.0,
            Cp_JkgK=2200.0,
            alpha_pK=1.0e-4,
            kTherm_WmK=3.0,
            Tconv_K=280.0,
            etaConv_Pas=5.0e14,
            etaMelt_Pas=5.0e14,
            eLid_m=2.0e5,
            Dconv_m=0.0,
            deltaTBL_m=0.0,
            RaConvect=1.0e4,
            RaCrit=1.0e6,
            Q_in_W=1.0e12,
        )

        self.assertFalse(state.productionCandidate)
        self.assertEqual(state.candidateStatus, "subcritical_rejected")
        self.assertFalse(state.updateAccepted)

    def test_ice_iii_and_v_remain_diagnostic_only(self):
        planet = self._synthetic_hp_planet(qBot_Wm2=1.0e-3)

        HPIceConvectionDiagnostics(planet, _Namespace(), ResolveHPIceConvectionModel(planet))

        self.assertEqual(planet.HPIceDiagnostics["III"]["candidateStatus"], "diagnostic_only_extrapolative")
        self.assertEqual(planet.HPIceDiagnostics["V"]["candidateStatus"], "diagnostic_only_extrapolative")
        self.assertFalse(planet.HPIceDiagnostics["III"]["updateAccepted"])
        self.assertFalse(planet.HPIceDiagnostics["V"]["updateAccepted"])

    def test_zero_contrast_rejected_without_collapsing_heat_to_zero(self):
        planet = self._synthetic_hp_planet(qBot_Wm2=1.0e-3, vi_deltaT_K=0.0)

        HPIceConvectionDiagnostics(planet, _Namespace(), ResolveHPIceConvectionModel(planet))

        state = planet.HPIceDiagnostics["VI"]
        self.assertEqual(state["candidateStatus"], "zero_contrast_rejected")
        self.assertGreater(state["Q_in_W"], 0.0)
        self.assertGreater(state["Q_out_W"], 0.0)
        self.assertEqual(state["energyResidual_W"], 0.0)
        self.assertFalse(state["updateAccepted"])

    def test_too_thin_ice_vi_case_is_rejected(self):
        planet = self._synthetic_hp_planet(qBot_Wm2=1.0e-3, vi_thickness_m=500.0)

        HPIceConvectionDiagnostics(planet, _Namespace(), ResolveHPIceConvectionModel(planet))

        state = planet.HPIceDiagnostics["VI"]
        self.assertEqual(state["candidateStatus"], "too_thin_rejected")
        self.assertFalse(state["productionCandidate"])
        self.assertFalse(state["updateAccepted"])

    def test_no_mass_transport_fields_are_active_outputs(self):
        state = HPIcePhaseState(
            "VI",
            phaseID=6,
            present=True,
            status="computed",
            productionMode="Kalousova2018_production_experimental",
            thickness_m=2.0e5,
            Ttop_K=280.0,
            Tbot_K=300.0,
            Tconv_K=280.0,
            etaConv_Pas=5.0e14,
            eLid_m=1.0e4,
            Dconv_m=1.8e5,
            deltaTBL_m=1.0e4,
            RaConvect=1.0e8,
            RaCrit=1.0e4,
            Q_in_W=1.0e12,
        )

        self.assertEqual(state.massResidual_kg, 0.0)
        self.assertEqual(state.massResidual_frac, 0.0)
        self.assertFalse(state.updateAccepted)


class Test09IceVIProductionAcceptanceCriteria(unittest.TestCase):

    def _accepted_ice_vi_state(self, **overrides):
        kwargs = dict(
            phaseName="VI",
            phaseID=6,
            present=True,
            status="computed",
            productionMode="Kalousova2018_production_experimental",
            iTop=4,
            iBot=8,
            rTop_m=2.3e6,
            rBot_m=2.1e6,
            zTop_m=1.0e5,
            zBot_m=3.0e5,
            thickness_m=2.0e5,
            Ttop_K=280.0,
            Tbot_K=300.0,
            Ptop_MPa=800.0,
            Pmid_MPa=950.0,
            Pbot_MPa=1100.0,
            rho_kgm3=1300.0,
            Cp_JkgK=2200.0,
            alpha_pK=1.0e-4,
            kTherm_WmK=3.0,
            Q_in_W=1.0e12,
            Q_out_W=1.0e12,
            Tconv_K=285.0,
            etaConv_Pas=5.0e14,
            etaMelt_Pas=5.0e14,
            eLid_m=1.0e4,
            Dconv_m=1.8e5,
            deltaTBL_m=1.0e4,
            RaConvect=1.0e8,
            RaCrit=1.0e4,
            meltFraction=0.05,
            phaseBoundaryResidual_K=0.01,
            Tmelt_top_K=292.0,
            Tmelt_mid_K=295.0,
            Tmelt_bot_K=298.0,
            viscositySource="fixed_reference",
        )
        kwargs.update(overrides)
        return HPIcePhaseState(**kwargs)

    def test_synthetic_ice_vi_candidate_can_be_accepted_without_profile_mutation(self):
        planet = _Namespace()
        planet.HPIceDiagnostics = {}
        planet.T_K = np.array([280.0, 300.0])
        planet.P_MPa = np.array([800.0, 1100.0])
        planet.phase = np.array([6, 6])
        planet.rho_kgm3 = np.array([1300.0, 1310.0])
        planet.eta_Pas = np.array([5.0e14, 5.0e14])
        before = {
            "T_K": planet.T_K.copy(),
            "P_MPa": planet.P_MPa.copy(),
            "phase": planet.phase.copy(),
            "rho_kgm3": planet.rho_kgm3.copy(),
            "eta_Pas": planet.eta_Pas.copy(),
        }

        _SetHPIceDiagnosticFields(
            planet,
            "VI",
            status="computed",
            phaseID=6,
            rTop_m=2.3e6,
            rBot_m=2.1e6,
            zTop_m=1.0e5,
            zBot_m=3.0e5,
            thickness_m=2.0e5,
            Ttop_K=280.0,
            Tbot_K=300.0,
            Ptop_MPa=800.0,
            Pmid_MPa=950.0,
            Pbot_MPa=1100.0,
            rho_kgm3=1300.0,
            Cp_JkgK=2200.0,
            alpha_pK=1.0e-4,
            kTherm_WmK=3.0,
            Q_in_W=1.0e12,
            Q_out_W=1.0e12,
            Tconv_K=285.0,
            etaConv_Pas=5.0e14,
            etaMelt_Pas=5.0e14,
            eLid_m=1.0e4,
            Dconv_m=1.8e5,
            deltaTBL_m=1.0e4,
            Ra=1.0e8,
            RaCrit=1.0e4,
            meltFraction=0.05,
            productionMode="Kalousova2018_production_experimental",
            phaseBoundaryResidual_K=0.01,
            Tmelt_top_K=292.0,
            Tmelt_mid_K=295.0,
            Tmelt_bot_K=298.0,
            viscositySource="fixed_reference",
        )

        state = planet.HPIceDiagnostics["VI"]
        self.assertEqual(state["candidateStatus"], "accepted_candidate_state_only")
        self.assertTrue(state["acceptanceCriteriaPassed"])
        self.assertTrue(state["updateAccepted"])
        self.assertEqual(state["acceptanceBlockers"], ())
        for key, value in before.items():
            np.testing.assert_array_equal(getattr(planet, key), value)

    def test_helper_returns_state_and_preserves_accepted_candidate(self):
        state = self._accepted_ice_vi_state()

        self.assertIs(_EvaluateIceVIProductionAcceptance(state), state)
        self.assertTrue(state.updateAccepted)
        self.assertEqual(state.candidateStatus, "accepted_candidate_state_only")

    def test_real_body_like_state_without_melt_curve_remains_rejected(self):
        state = self._accepted_ice_vi_state(Tmelt_top_K=np.nan, Tmelt_mid_K=np.nan, Tmelt_bot_K=np.nan)

        self.assertFalse(state.updateAccepted)
        self.assertFalse(state.acceptanceCriteriaPassed)
        self.assertEqual(state.candidateStatus, "missing_melt_curve_rejected")
        self.assertIn("missing_melt_curve_rejected", state.acceptanceBlockers)

    def test_ice_iii_and_v_remain_diagnostic_only_extrapolative(self):
        for phase_name, phase_id in (("III", 3), ("V", 5)):
            with self.subTest(phase=phase_name):
                state = self._accepted_ice_vi_state(phaseName=phase_name, phaseID=phase_id)

                self.assertFalse(state.updateAccepted)
                self.assertEqual(state.candidateStatus, "diagnostic_only_extrapolative")

    def test_phase_boundary_residual_too_large_rejects(self):
        state = self._accepted_ice_vi_state(phaseBoundaryResidual_K=0.2)

        self.assertFalse(state.updateAccepted)
        self.assertEqual(state.candidateStatus, "phase_boundary_rejected")

    def test_energy_residual_too_large_rejects(self):
        state = self._accepted_ice_vi_state(Q_in_W=100.0, Q_out_W=101.0)

        self.assertFalse(state.updateAccepted)
        self.assertEqual(state.candidateStatus, "energy_residual_rejected")

    def test_nonzero_mass_residual_rejects(self):
        state = self._accepted_ice_vi_state(massResidual_kg=1.0, massResidual_frac=1.0e-12)

        self.assertFalse(state.updateAccepted)
        self.assertEqual(state.candidateStatus, "mass_residual_rejected")

    def test_subcritical_ice_vi_rejects(self):
        state = self._accepted_ice_vi_state(RaConvect=1.0e4, RaCrit=1.0e6)

        self.assertFalse(state.updateAccepted)
        self.assertEqual(state.candidateStatus, "subcritical_rejected")

    def test_zero_contrast_ice_vi_rejects(self):
        state = self._accepted_ice_vi_state(Tbot_K=280.0)

        self.assertFalse(state.updateAccepted)
        self.assertEqual(state.candidateStatus, "zero_contrast_rejected")

    def test_too_thin_ice_vi_rejects(self):
        state = self._accepted_ice_vi_state(
            zBot_m=1.005e5,
            rBot_m=2.2995e6,
            thickness_m=500.0,
            eLid_m=100.0,
            Dconv_m=300.0,
            deltaTBL_m=100.0,
        )

        self.assertFalse(state.updateAccepted)
        self.assertEqual(state.candidateStatus, "too_thin_rejected")

    def test_boundary_layer_exceeds_layer_rejects(self):
        state = self._accepted_ice_vi_state(eLid_m=5.0e4, Dconv_m=1.8e5, deltaTBL_m=1.0e4)

        self.assertFalse(state.updateAccepted)
        self.assertEqual(state.candidateStatus, "boundary_layer_exceeds_layer_rejected")

    def test_invalid_geometry_and_nonfinite_values_reject(self):
        invalid_geometry = self._accepted_ice_vi_state(rBot_m=2.4e6)
        nonfinite_value = self._accepted_ice_vi_state(etaConv_Pas=np.nan)

        self.assertEqual(invalid_geometry.candidateStatus, "invalid_geometry_rejected")
        self.assertEqual(nonfinite_value.candidateStatus, "nonfinite_rejected")
        self.assertFalse(invalid_geometry.updateAccepted)
        self.assertFalse(nonfinite_value.updateAccepted)


class Test10IceVIMeltCurveCandidateChecks(unittest.TestCase):

    def _production_planet(self):
        planet = _Namespace()
        planet.Ocean = _Namespace()
        planet.Ocean.EOS = _SyntheticIceVIMeltEOS()
        planet.Ocean.THydroMax_K = 340.0
        planet.TfreezeUpper_K = 340.0
        planet.TfreezeRes_K = 0.01
        return planet

    def _ice_vi_state(self, **overrides):
        kwargs = dict(
            phaseName="VI",
            phaseID=6,
            present=True,
            status="computed",
            productionMode="Kalousova2018_production_experimental",
            rTop_m=2.3e6,
            rBot_m=2.1e6,
            zTop_m=1.0e5,
            zBot_m=3.0e5,
            thickness_m=2.0e5,
            Ttop_K=285.0,
            Tbot_K=287.0,
            Ptop_MPa=800.0,
            Pmid_MPa=900.0,
            Pbot_MPa=1000.0,
            rho_kgm3=1300.0,
            Cp_JkgK=2200.0,
            alpha_pK=1.0e-4,
            kTherm_WmK=3.0,
            Q_in_W=1.0e12,
            Q_out_W=1.0e12,
            Tconv_K=286.0,
            etaConv_Pas=5.0e14,
            etaMelt_Pas=5.0e14,
            eLid_m=1.0e4,
            Dconv_m=1.8e5,
            deltaTBL_m=1.0e4,
            RaConvect=1.0e8,
            RaCrit=1.0e4,
            meltFraction=0.05,
            phaseBoundaryResidual_K=0.0,
            Tmelt_top_K=292.01,
            Tmelt_mid_K=293.01,
            Tmelt_bot_K=294.01,
            TmeltSource="GetTfreeze",
            meltCurveStatus="melt_curve_ok",
            phaseBoundaryStatus="phase_boundary_ok",
            viscositySource="fixed_reference",
        )
        kwargs.update(overrides)
        return HPIcePhaseState(**kwargs)

    def test_missing_melt_curve_keeps_ice_vi_rejected(self):
        state = self._ice_vi_state(Tmelt_top_K=np.nan, Tmelt_mid_K=np.nan, Tmelt_bot_K=np.nan)

        self.assertFalse(state.updateAccepted)
        self.assertEqual(state.candidateStatus, "missing_melt_curve_rejected")
        self.assertEqual(state.meltCurveStatus, "missing_melt_curve_rejected")

    def test_nonfinite_melt_curve_values_reject(self):
        state = self._ice_vi_state(Tmelt_mid_K=np.inf)

        self.assertFalse(state.updateAccepted)
        self.assertEqual(state.candidateStatus, "melt_curve_nonfinite_rejected")
        self.assertEqual(state.meltCurveStatus, "melt_curve_nonfinite_rejected")

    def test_valid_synthetic_melt_curve_records_three_temperatures(self):
        planet = self._production_planet()

        fields = _GetIceVIMeltCurveCandidateChecks(
            planet,
            phaseID=6,
            Ptop_MPa=800.0,
            Pmid_MPa=900.0,
            Pbot_MPa=1000.0,
            Ttop_K=285.0,
            Tconv_K=286.0,
            Tbot_K=287.0,
            productionMode="Kalousova2018_production_experimental",
        )

        self.assertEqual(fields["TmeltSource"], "GetTfreeze")
        self.assertEqual(fields["meltCurveStatus"], "melt_curve_ok")
        self.assertEqual(fields["phaseBoundaryStatus"], "phase_boundary_ok")
        self.assertEqual(fields["eosBoundaryContext"], "in_run_candidate")
        self.assertEqual(fields["eosBoundaryStatus"], "in_run_candidate_boundary_in_domain")
        self.assertEqual(fields["candidateBoundarySource"], "in_run_candidate_bounds")
        self.assertEqual(fields["finalProfileCoverageStatus"], "final_profile_not_evaluated")
        np.testing.assert_allclose(
            [fields["Tmelt_top_K"], fields["Tmelt_mid_K"], fields["Tmelt_bot_K"]],
            [292.01, 293.01, 294.01],
            atol=2.0e-2,
        )

    def test_outside_eos_domain_records_in_run_boundary_context(self):
        planet = self._production_planet()

        fields = _GetIceVIMeltCurveCandidateChecks(
            planet,
            phaseID=6,
            Ptop_MPa=800.0,
            Pmid_MPa=900.0,
            Pbot_MPa=1400.0,
            Ttop_K=285.0,
            Tconv_K=286.0,
            Tbot_K=287.0,
            productionMode="Kalousova2018_production_experimental",
        )
        state = self._ice_vi_state(**fields)

        self.assertEqual(fields["meltCurveStatus"], "outside_eos_domain_rejected")
        self.assertEqual(fields["eosBoundaryContext"], "in_run_candidate")
        self.assertEqual(fields["eosBoundaryStatus"], "in_run_candidate_boundary_outside_eos_domain")
        self.assertEqual(fields["eosBoundaryReason"], "candidate_boundary_outside_declared_eos_domain")
        self.assertEqual(fields["candidateBoundarySource"], "in_run_candidate_bounds")
        self.assertFalse(state.updateAccepted)
        self.assertEqual(state.candidateStatus, "outside_eos_domain_rejected")
        self.assertEqual(state.candidateReason, "in_run_candidate_boundary_outside_eos_domain")

    def test_final_profile_coverage_status_does_not_change_acceptance(self):
        state = self._ice_vi_state(
            Tmelt_top_K=np.nan,
            Tmelt_mid_K=np.nan,
            Tmelt_bot_K=np.nan,
            TmeltSource=None,
            meltCurveStatus="outside_eos_domain_rejected",
            phaseBoundaryStatus="phase_boundary_unavailable_rejected",
            eosBoundaryContext="in_run_candidate",
            eosBoundaryStatus="in_run_candidate_boundary_outside_eos_domain",
            eosBoundaryReason="candidate_boundary_outside_declared_eos_domain",
            candidateBoundarySource="in_run_candidate_bounds",
            finalProfileCoverageStatus="final_profile_nodes_in_domain",
        )

        self.assertFalse(state.updateAccepted)
        self.assertEqual(state.candidateStatus, "outside_eos_domain_rejected")
        self.assertEqual(state.candidateReason, "in_run_candidate_boundary_outside_eos_domain")
        self.assertEqual(state.finalProfileCoverageStatus, "final_profile_nodes_in_domain")

    def test_phase_boundary_residual_can_pass_candidate_state_only(self):
        planet = self._production_planet()
        fields = _GetIceVIMeltCurveCandidateChecks(
            planet, 6, 800.0, 900.0, 1000.0, 285.0, 286.0, 287.0,
            "Kalousova2018_production_experimental",
        )
        state = self._ice_vi_state(**fields)

        self.assertTrue(state.updateAccepted)
        self.assertTrue(state.acceptanceCriteriaPassed)
        self.assertEqual(state.candidateStatus, "accepted_candidate_state_only")

    def test_phase_boundary_residual_above_threshold_rejects(self):
        state = self._ice_vi_state(phaseBoundaryResidual_K=0.2, phaseBoundaryStatus="phase_boundary_rejected")

        self.assertFalse(state.updateAccepted)
        self.assertEqual(state.candidateStatus, "phase_boundary_rejected")

    def test_real_body_like_state_without_verified_melt_curve_remains_rejected(self):
        state = self._ice_vi_state(
            Tmelt_top_K=np.nan,
            Tmelt_mid_K=np.nan,
            Tmelt_bot_K=np.nan,
            TmeltSource=None,
            meltCurveStatus="not_evaluated",
            phaseBoundaryResidual_K=np.nan,
            phaseBoundaryStatus="not_evaluated",
        )

        self.assertFalse(state.updateAccepted)
        self.assertEqual(state.candidateStatus, "missing_melt_curve_rejected")

    def test_fixed_diagnostic_fallback_is_not_allowed_for_acceptance(self):
        state = self._ice_vi_state(
            Tmelt_top_K=290.0,
            Tmelt_mid_K=290.0,
            Tmelt_bot_K=290.0,
            TmeltSource="fixed_diagnostic_fallback",
        )

        self.assertFalse(state.updateAccepted)
        self.assertEqual(state.candidateStatus, "missing_melt_curve_rejected")
        self.assertEqual(state.candidateReason, "fixed_diagnostic_fallback_not_allowed")

    def test_ice_iii_and_v_remain_diagnostic_only_extrapolative(self):
        for phase_name, phase_id in (("III", 3), ("V", 5)):
            with self.subTest(phase=phase_name):
                state = self._ice_vi_state(phaseName=phase_name, phaseID=phase_id)

                self.assertFalse(state.updateAccepted)
                self.assertEqual(state.candidateStatus, "diagnostic_only_extrapolative")

    def test_candidate_melt_curve_check_does_not_mutate_profile_arrays(self):
        planet = self._production_planet()
        planet.HPIceDiagnostics = {}
        planet.T_K = np.array([285.0, 287.0])
        planet.P_MPa = np.array([800.0, 1000.0])
        planet.phase = np.array([6, 6])
        planet.rho_kgm3 = np.array([1300.0, 1310.0])
        planet.eta_Pas = np.array([5.0e14, 5.0e14])
        before = {
            "T_K": planet.T_K.copy(),
            "P_MPa": planet.P_MPa.copy(),
            "phase": planet.phase.copy(),
            "rho_kgm3": planet.rho_kgm3.copy(),
            "eta_Pas": planet.eta_Pas.copy(),
        }
        fields = _GetIceVIMeltCurveCandidateChecks(
            planet, 6, 800.0, 900.0, 1000.0, 285.0, 286.0, 287.0,
            "Kalousova2018_production_experimental",
        )

        _SetHPIceDiagnosticFields(
            planet,
            "VI",
            status="computed",
            phaseID=6,
            rTop_m=2.3e6,
            rBot_m=2.1e6,
            zTop_m=1.0e5,
            zBot_m=3.0e5,
            thickness_m=2.0e5,
            Ttop_K=285.0,
            Tbot_K=287.0,
            Ptop_MPa=800.0,
            Pmid_MPa=900.0,
            Pbot_MPa=1000.0,
            rho_kgm3=1300.0,
            Cp_JkgK=2200.0,
            alpha_pK=1.0e-4,
            kTherm_WmK=3.0,
            Q_in_W=1.0e12,
            Q_out_W=1.0e12,
            Tconv_K=286.0,
            etaConv_Pas=5.0e14,
            etaMelt_Pas=5.0e14,
            eLid_m=1.0e4,
            Dconv_m=1.8e5,
            deltaTBL_m=1.0e4,
            Ra=1.0e8,
            RaCrit=1.0e4,
            meltFraction=0.05,
            productionMode="Kalousova2018_production_experimental",
            **fields,
        )

        self.assertTrue(planet.HPIceDiagnostics["VI"]["updateAccepted"])
        self.assertEqual(
            planet.HPIceDiagnostics["VI"]["eosBoundaryStatus"],
            "in_run_candidate_boundary_in_domain",
        )
        for key, value in before.items():
            np.testing.assert_array_equal(getattr(planet, key), value)


class Test11PosthocIceVIProductionCandidate(unittest.TestCase):

    def _planet(self, phase=None, P_MPa=None, T_K=None, eos=None, diagnostics=None):
        planet = _Namespace()
        planet.Ocean = _Namespace()
        planet.Ocean.EOS = eos or _SyntheticIceVIMeltEOS()
        planet.Ocean.THydroMax_K = 340.0
        planet.TfreezeUpper_K = 340.0
        planet.TfreezeRes_K = 0.01
        planet.phase = np.array([0, 6, 6, 6, 0] if phase is None else phase)
        planet.P_MPa = np.array([100.0, 800.0, 900.0, 1000.0, 1100.0] if P_MPa is None else P_MPa, dtype=float)
        planet.T_K = np.array([275.0, 285.0, 286.0, 287.0, 288.0] if T_K is None else T_K, dtype=float)
        n_nodes = planet.phase.size
        planet.r_m = np.linspace(2.5e6, 2.1e6, n_nodes)
        planet.z_m = 2.5e6 - planet.r_m
        planet.rho_kgm3 = np.zeros(n_nodes) + 1300.0
        planet.eta_Pas = np.zeros(n_nodes) + 5.0e14
        planet.kTherm_WmK = np.zeros(n_nodes) + 3.0
        planet.qSurf_Wm2 = 0.02
        planet.qCon_Wm2 = 0.02
        planet.Mtot_kg = 1.0e23
        planet.CMR2mean = 0.33
        if diagnostics is None:
            q_in_Wm2 = 1.0e12 / (4.0 * np.pi * planet.r_m[3]**2)
            q_out_Wm2 = 1.0e12 / (4.0 * np.pi * planet.r_m[1]**2)
            diagnostics = {
                "VI": {
                    "energyResidual_W": 0.0,
                    "energyResidual_frac": 0.0,
                    "Q_in_W": 1.0e12,
                    "Q_out_W": 1.0e12,
                    "q_in_Wm2": q_in_Wm2,
                    "q_out_Wm2": q_out_Wm2,
                    "massResidual_kg": 0.0,
                    "massResidual_frac": 0.0,
                    "etaConv_Pas": 5.0e14,
                    "etaMelt_Pas": 5.0e14,
                    "RaConvect": 1.0e8,
                    "RaCrit": 1.0e4,
                    "layerClosureResidual_m": 0.0,
                    "kTherm_WmK": 3.0,
                },
                "III": {"candidateStatus": "diagnostic_only_extrapolative"},
                "V": {"candidateStatus": "diagnostic_only_extrapolative"},
            }
        planet.HPIceDiagnostics = diagnostics
        return planet

    @staticmethod
    def _snapshot(planet):
        return {
            "phase": planet.phase.copy(),
            "P_MPa": planet.P_MPa.copy(),
            "T_K": planet.T_K.copy(),
            "rho_kgm3": planet.rho_kgm3.copy(),
            "eta_Pas": planet.eta_Pas.copy(),
            "qSurf_Wm2": planet.qSurf_Wm2,
            "qCon_Wm2": planet.qCon_Wm2,
            "Mtot_kg": planet.Mtot_kg,
            "CMR2mean": planet.CMR2mean,
        }

    def test_finalized_in_domain_ice_vi_can_pass_posthoc_candidate_state(self):
        planet = self._planet()
        before = self._snapshot(planet)

        result = EvaluatePosthocIceVIProductionCandidate(
            planet, productionMode="Kalousova2018_production_experimental",
        )

        self.assertTrue(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocCandidateStatus"], "posthoc_candidate_passed")
        self.assertEqual(result["posthocCandidateReason"], "posthoc_candidate_state_only")
        self.assertEqual(result["posthocBoundarySource"], "finalized_profile_nodes")
        self.assertEqual(result["posthocMeltCurveStatus"], "melt_curve_ok")
        self.assertEqual(result["posthocPhaseBoundaryStatus"], "phase_boundary_ok")
        self.assertEqual(result["finalProfileCoverageStatus"], "final_profile_nodes_in_domain")
        self.assertEqual(result["posthocSensitivityRiskStatus"], "nominal")
        self.assertEqual(
            planet.HPIceDiagnostics["VI"]["posthocProductionCandidate"],
            result,
        )
        self.assertEqual(
            planet.HPIceDiagnostics["III"]["posthocProductionCandidate"]["posthocCandidateStatus"],
            "diagnostic_only_extrapolative",
        )
        self.assertEqual(
            planet.HPIceDiagnostics["V"]["posthocProductionCandidate"]["posthocCandidateStatus"],
            "diagnostic_only_extrapolative",
        )
        for key, value in before.items():
            if isinstance(value, np.ndarray):
                np.testing.assert_array_equal(getattr(planet, key), value)
            else:
                self.assertEqual(getattr(planet, key), value)

    def test_outside_eos_finalized_bounds_reject(self):
        planet = self._planet(P_MPa=[100.0, 800.0, 900.0, 1400.0, 1500.0])

        result = EvaluatePosthocIceVIProductionCandidate(
            planet, productionMode="Kalousova2018_production_experimental",
        )

        self.assertFalse(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocCandidateStatus"], "posthoc_outside_eos_domain_rejected")
        self.assertEqual(result["posthocCandidateReason"], "posthoc_all_nodes_outside_eos_domain")

    def test_wrong_eos_phase_rejects(self):
        planet = self._planet(eos=_SyntheticWrongPhaseEOS())

        result = EvaluatePosthocIceVIProductionCandidate(
            planet, productionMode="Kalousova2018_production_experimental",
        )

        self.assertFalse(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocCandidateStatus"], "posthoc_outside_eos_domain_rejected")
        self.assertEqual(result["posthocCandidateReason"], "posthoc_all_nodes_not_ice_vi")

    def test_missing_ice_vi_rejects(self):
        planet = self._planet(phase=[0, 5, 5, 5, 0])

        result = EvaluatePosthocIceVIProductionCandidate(
            planet, productionMode="Kalousova2018_production_experimental",
        )

        self.assertFalse(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocCandidateStatus"], "posthoc_missing_ice_vi")

    def test_nonfinite_melt_curve_rejects(self):
        planet = self._planet()

        with mock.patch("PlanetProfile.Thermodynamics.LayerPropagators.GetTfreeze", return_value=np.nan):
            result = EvaluatePosthocIceVIProductionCandidate(
                planet, productionMode="Kalousova2018_production_experimental",
            )

        self.assertFalse(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocCandidateStatus"], "posthoc_nonfinite_rejected")
        self.assertEqual(result["posthocMeltCurveStatus"], "melt_curve_nonfinite_rejected")

    def test_phase_boundary_residual_failure_rejects(self):
        planet = self._planet()

        with mock.patch("PlanetProfile.Thermodynamics.LayerPropagators.GetTfreeze", return_value=280.0):
            result = EvaluatePosthocIceVIProductionCandidate(
                planet, productionMode="Kalousova2018_production_experimental",
            )

        self.assertFalse(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocCandidateStatus"], "posthoc_phase_boundary_rejected")
        self.assertEqual(result["posthocPhaseBoundaryStatus"], "phase_boundary_rejected")

    def test_energy_residual_failure_rejects(self):
        planet = self._planet()
        planet.HPIceDiagnostics["VI"]["energyResidual_W"] = 1.0e9
        planet.HPIceDiagnostics["VI"]["energyResidual_frac"] = 1.0e-3

        result = EvaluatePosthocIceVIProductionCandidate(
            planet, productionMode="Kalousova2018_production_experimental",
        )

        self.assertFalse(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocCandidateStatus"], "posthoc_energy_residual_rejected")

    def test_mass_residual_failure_rejects(self):
        planet = self._planet()
        planet.HPIceDiagnostics["VI"]["massResidual_kg"] = 1.0

        result = EvaluatePosthocIceVIProductionCandidate(
            planet, productionMode="Kalousova2018_production_experimental",
        )

        self.assertFalse(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocCandidateStatus"], "posthoc_mass_residual_rejected")


class Test12PosthocIceVIMarginRiskDiagnostics(unittest.TestCase):

    def _planet(self, **kwargs):
        return Test11PosthocIceVIProductionCandidate()._planet(**kwargs)

    @staticmethod
    def _evaluate(planet):
        return EvaluatePosthocIceVIProductionCandidate(
            planet, productionMode="Kalousova2018_production_experimental",
        )

    def test_nominal_synthetic_ice_vi_records_margins_and_nominal_risk(self):
        planet = self._planet()

        result = self._evaluate(planet)

        self.assertTrue(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocSensitivityRiskStatus"], "nominal")
        self.assertEqual(result["posthocRiskReasons"], ())
        self.assertTrue(result["posthocAllNodesInEOSDomain"])
        self.assertTrue(result["posthocAllNodesIceVI"])
        self.assertGreater(result["posthocEOSPressureMargin_MPa"], 0.0)
        self.assertGreater(result["posthocEOSTemperatureMargin_K"], 0.0)
        self.assertGreater(result["posthocMinPhaseBoundaryMargin_K"], 0.0)

    def test_near_eos_boundary_records_risk(self):
        planet = self._planet(P_MPa=[100.0, 800.0, 900.0, 1299.9, 1100.0])

        result = self._evaluate(planet)

        self.assertTrue(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocSensitivityRiskStatus"], "near_eos_boundary")
        self.assertIn("near_eos_boundary", result["posthocRiskReasons"])
        self.assertAlmostEqual(result["posthocEOSPressureMargin_MPa"], 0.1)

    def test_near_phase_boundary_records_risk(self):
        planet = self._planet(T_K=[275.0, 285.0, 286.0, 293.98, 288.0])

        result = self._evaluate(planet)

        self.assertTrue(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocSensitivityRiskStatus"], "near_phase_boundary")
        self.assertIn("near_phase_boundary", result["posthocRiskReasons"])
        self.assertGreaterEqual(result["posthocMinPhaseBoundaryMargin_K"], 0.0)
        self.assertLessEqual(result["posthocMinPhaseBoundaryMargin_K"], 0.1)

    def test_zero_contrast_is_high_risk_rejected(self):
        planet = self._planet(T_K=[275.0, 285.0, 285.0, 285.0, 288.0])

        result = self._evaluate(planet)

        self.assertFalse(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocCandidateStatus"], "posthoc_high_risk_rejected")
        self.assertEqual(result["posthocTemperatureContrastStatus"], "invalid_contrast")
        self.assertIn("invalid_contrast", result["posthocRiskReasons"])

    def test_subcritical_rayleigh_is_high_risk_rejected(self):
        planet = self._planet()
        planet.HPIceDiagnostics["VI"]["RaConvect"] = 1.0e3
        planet.HPIceDiagnostics["VI"]["RaCrit"] = 1.0e8

        result = self._evaluate(planet)

        self.assertFalse(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocCandidateStatus"], "posthoc_high_risk_rejected")
        self.assertEqual(result["posthocRayleighRegimeStatus"], "subcritical")
        self.assertIn("subcritical", result["posthocRiskReasons"])

    def test_too_thin_geometry_is_high_risk_rejected(self):
        planet = self._planet()
        planet.z_m[1:4] = np.array([100000.0, 100400.0, 100800.0])

        result = self._evaluate(planet)

        self.assertFalse(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocCandidateStatus"], "posthoc_high_risk_rejected")
        self.assertEqual(result["posthocThicknessStatus"], "too_thin")
        self.assertIn("too_thin", result["posthocRiskReasons"])

    def test_boundary_layer_exceeds_layer_is_high_risk_rejected(self):
        planet = self._planet()
        planet.HPIceDiagnostics["VI"]["layerClosureResidual_m"] = 1.0

        result = self._evaluate(planet)

        self.assertFalse(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocCandidateStatus"], "posthoc_high_risk_rejected")
        self.assertEqual(result["posthocLayerClosureStatus"], "boundary_layer_exceeds_layer")
        self.assertIn("boundary_layer_exceeds_layer", result["posthocRiskReasons"])

    def test_nonfinite_viscosity_is_high_risk_rejected(self):
        planet = self._planet()
        planet.HPIceDiagnostics["VI"]["etaConv_Pas"] = np.nan

        result = self._evaluate(planet)

        self.assertFalse(result["posthocUpdateAccepted"])
        self.assertEqual(result["posthocCandidateStatus"], "posthoc_high_risk_rejected")
        self.assertEqual(result["posthocViscosityStatus"], "invalid_viscosity")
        self.assertIn("invalid_viscosity", result["posthocRiskReasons"])

    def test_all_node_eos_and_phase_checks_are_recorded(self):
        outside_domain = self._planet(
            phase=[0, 6, 6, 6, 6, 6, 0],
            P_MPa=[100.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1100.0],
            T_K=[275.0, 285.0, 400.0, 286.0, 287.0, 287.0, 288.0],
        )

        result = self._evaluate(outside_domain)

        self.assertFalse(result["posthocUpdateAccepted"])
        self.assertFalse(result["posthocAllNodesInEOSDomain"])
        self.assertEqual(result["posthocCandidateReason"], "posthoc_all_nodes_outside_eos_domain")

        wrong_phase = self._planet(
            phase=[0, 6, 6, 6, 6, 6, 0],
            P_MPa=[100.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1100.0],
            T_K=[275.0, 285.0, 295.0, 286.0, 287.0, 287.0, 288.0],
        )

        result = self._evaluate(wrong_phase)

        self.assertFalse(result["posthocUpdateAccepted"])
        self.assertTrue(result["posthocAllNodesInEOSDomain"])
        self.assertFalse(result["posthocAllNodesIceVI"])
        self.assertEqual(result["posthocCandidateReason"], "posthoc_all_nodes_not_ice_vi")

    def test_no_profile_mutation_from_margin_checks(self):
        planet = self._planet()
        before = Test11PosthocIceVIProductionCandidate._snapshot(planet)

        self._evaluate(planet)

        for key, value in before.items():
            if isinstance(value, np.ndarray):
                np.testing.assert_array_equal(getattr(planet, key), value)
            else:
                self.assertEqual(getattr(planet, key), value)


class Test13ActiveIceVICandidateProfileCopy(unittest.TestCase):

    def _planet(self, **kwargs):
        planet = Test11PosthocIceVIProductionCandidate()._planet(**kwargs)
        EvaluatePosthocIceVIProductionCandidate(
            planet, productionMode="Kalousova2018_production_experimental",
        )
        return planet

    @staticmethod
    def _copy(planet):
        return BuildActiveIceVIProductionCandidateCopy(
            planet, productionMode="Kalousova2018_production_experimental",
        )

    def test_candidate_copy_created_for_valid_finalized_posthoc_state(self):
        planet = self._planet()

        result = self._copy(planet)

        self.assertTrue(result["candidateCopyCreated"])
        self.assertEqual(result["candidateCopySource"], "finalized_posthoc_ice_vi_nodes")
        self.assertEqual(result["candidatePhase"], "VI")
        self.assertEqual(result["candidateNodeCount"], 3)
        self.assertEqual(result["candidateIndexStart"], 1)
        self.assertEqual(result["candidateIndexEnd"], 3)
        self.assertEqual(result["candidateStatus"], "candidate_copy_created")
        self.assertEqual(result["candidateReason"], "candidate_copy_state_only")
        self.assertFalse(result["candidateAppliedToProfile"])
        self.assertFalse(result["candidateAccepted"])
        self.assertTrue(result["protectedFieldsUnchanged"])
        self.assertIs(
            planet.HPIceDiagnostics["VI"]["activeProductionCandidate"],
            result,
        )

    def test_candidate_copy_contains_only_phase_6_nodes(self):
        planet = self._planet()

        result = self._copy(planet)

        np.testing.assert_array_equal(result["candidatePhaseArray"], np.array([6, 6, 6]))
        np.testing.assert_allclose(result["candidateP_MPa"], planet.P_MPa[1:4])
        np.testing.assert_allclose(result["candidateT_K"], planet.T_K[1:4])

    def test_candidate_arrays_are_independent_copies(self):
        planet = self._planet()

        result = self._copy(planet)

        self.assertFalse(np.shares_memory(result["candidateP_MPa"], planet.P_MPa))
        self.assertFalse(np.shares_memory(result["candidateT_K"], planet.T_K))
        self.assertFalse(np.shares_memory(result["candidatePhaseArray"], planet.phase))
        self.assertFalse(np.shares_memory(result["candidateR_m"], planet.r_m))
        self.assertFalse(np.shares_memory(result["candidateZ_m"], planet.z_m))
        self.assertFalse(np.shares_memory(result["candidateEta_Pas"], planet.eta_Pas))

    def test_modifying_candidate_copy_does_not_modify_planet_fields(self):
        planet = self._planet()
        before = Test11PosthocIceVIProductionCandidate._snapshot(planet)

        result = self._copy(planet)
        result["candidateP_MPa"][0] += 10.0
        result["candidateT_K"][0] += 10.0
        result["candidatePhaseArray"][0] = 0
        result["candidateEta_Pas"][0] *= 2.0

        for key, value in before.items():
            if isinstance(value, np.ndarray):
                np.testing.assert_array_equal(getattr(planet, key), value)
            else:
                self.assertEqual(getattr(planet, key), value)

    def test_missing_ice_vi_rejects(self):
        planet = self._planet(phase=[0, 5, 5, 5, 0])

        result = self._copy(planet)

        self.assertFalse(result["candidateCopyCreated"])
        self.assertFalse(result["candidateAccepted"])
        self.assertEqual(result["candidateStatus"], "posthoc_candidate_not_passed_rejected")
        self.assertEqual(result["candidateReason"], "requires_at_least_two_finalized_ice_vi_nodes")

    def test_high_risk_posthoc_candidate_rejects(self):
        planet = self._planet(T_K=[275.0, 285.0, 285.0, 285.0, 288.0])

        result = self._copy(planet)

        self.assertFalse(result["candidateCopyCreated"])
        self.assertFalse(result["candidateAccepted"])
        self.assertEqual(result["candidateStatus"], "high_risk_posthoc_rejected")
        self.assertEqual(result["candidateReason"], "posthoc_candidate_is_high_risk")

    def test_ice_iii_and_v_remain_diagnostic_only_extrapolative(self):
        planet = self._planet()

        self._copy(planet)

        self.assertEqual(
            planet.HPIceDiagnostics["III"]["activeProductionCandidate"]["candidateStatus"],
            "diagnostic_only_extrapolative",
        )
        self.assertEqual(
            planet.HPIceDiagnostics["V"]["activeProductionCandidate"]["candidateStatus"],
            "diagnostic_only_extrapolative",
        )

    def test_candidate_applied_and_accepted_remain_false(self):
        planet = self._planet()

        result = self._copy(planet)

        self.assertFalse(result["candidateAppliedToProfile"])
        self.assertFalse(result["candidateAccepted"])
        self.assertFalse(result["rollbackRequired"])
        self.assertFalse(result["rollbackApplied"])


class Test14ActiveIceVICandidateResidualEvaluator(unittest.TestCase):

    def _planet_with_copy(self):
        planet = Test13ActiveIceVICandidateProfileCopy()._planet()
        BuildActiveIceVIProductionCandidateCopy(
            planet, productionMode="Kalousova2018_production_experimental",
        )
        return planet

    @staticmethod
    def _evaluate(planet):
        return EvaluateActiveIceVICandidateResiduals(
            planet, productionMode="Kalousova2018_production_experimental",
        )

    @staticmethod
    def _candidate(planet):
        return planet.HPIceDiagnostics["VI"]["activeProductionCandidate"]

    def test_residual_evaluator_requires_candidate_copy(self):
        planet = Test13ActiveIceVICandidateProfileCopy()._planet()

        result = self._evaluate(planet)

        self.assertFalse(result["candidateResidualsPassed"])
        self.assertEqual(result["candidateResidualStatus"], "candidate_missing_required_field_rejected")
        self.assertIn("requires_active_candidate_copy", result["candidateResidualReasons"])
        self.assertFalse(result["candidateAccepted"])
        self.assertFalse(result["candidateAppliedToProfile"])

    def test_nominal_candidate_copy_records_residuals(self):
        planet = self._planet_with_copy()

        result = self._evaluate(planet)

        self.assertTrue(result["candidateResidualsPassed"])
        self.assertEqual(result["candidateResidualStatus"], "candidate_residuals_passed")
        self.assertEqual(result["candidateResidualReasons"], ())
        self.assertEqual(result["candidateEnergyResidual_W"], 0.0)
        self.assertEqual(result["candidateEnergyResidual_frac"], 0.0)
        self.assertEqual(result["candidateHeatPowerResidual_W"], 0.0)
        self.assertEqual(result["candidateSphericalFluxScalingStatus"], "spherical_area_scaled")
        self.assertAlmostEqual(
            result["candidateHeatFluxResidual_Wm2"],
            result["candidateExpected_q_out_Wm2"] - result["candidateExpected_q_in_Wm2"],
        )
        self.assertEqual(result["candidateMassResidual_kg"], 0.0)
        self.assertEqual(result["candidateMassResidual_frac"], 0.0)
        self.assertGreater(result["candidateRaOverRaCrit"], 1.0)
        self.assertFalse(result["candidateAccepted"])
        self.assertFalse(result["candidateAppliedToProfile"])

    def test_energy_residual_failure_rejects(self):
        planet = self._planet_with_copy()
        planet.HPIceDiagnostics["VI"]["energyResidual_W"] = 1.0e9
        planet.HPIceDiagnostics["VI"]["energyResidual_frac"] = 1.0e-3

        result = self._evaluate(planet)

        self.assertFalse(result["candidateResidualsPassed"])
        self.assertEqual(result["candidateResidualStatus"], "candidate_energy_residual_rejected")
        self.assertIn("energy_residual_exceeds_tolerance", result["candidateResidualReasons"])

    def test_heat_flux_residual_failure_rejects(self):
        planet = self._planet_with_copy()
        candidate = self._candidate(planet)
        candidate["candidateq_out_Wm2"] = 0.021

        result = self._evaluate(planet)

        self.assertFalse(result["candidateResidualsPassed"])
        self.assertEqual(result["candidateResidualStatus"], "candidate_heat_flux_residual_rejected")
        self.assertIn("spherical_flux_scaling_mismatch", result["candidateResidualReasons"])

    def test_nonzero_mass_residual_rejects(self):
        planet = self._planet_with_copy()
        planet.HPIceDiagnostics["VI"]["massResidual_kg"] = 1.0

        result = self._evaluate(planet)

        self.assertFalse(result["candidateResidualsPassed"])
        self.assertEqual(result["candidateResidualStatus"], "candidate_mass_residual_rejected")
        self.assertIn("mass_residual_nonzero", result["candidateResidualReasons"])

    def test_phase_boundary_residual_failure_rejects(self):
        planet = self._planet_with_copy()
        planet.HPIceDiagnostics["VI"]["posthocProductionCandidate"]["posthocPhaseBoundaryResidual_K"] = 0.2

        result = self._evaluate(planet)

        self.assertFalse(result["candidateResidualsPassed"])
        self.assertEqual(result["candidateResidualStatus"], "candidate_phase_boundary_rejected")
        self.assertIn("phase_boundary_residual_exceeds_tolerance", result["candidateResidualReasons"])

    def test_layer_closure_residual_failure_rejects(self):
        planet = self._planet_with_copy()
        planet.HPIceDiagnostics["VI"]["layerClosureResidual_m"] = 1.0

        result = self._evaluate(planet)

        self.assertFalse(result["candidateResidualsPassed"])
        self.assertEqual(result["candidateResidualStatus"], "candidate_layer_closure_rejected")
        self.assertIn("layer_closure_residual_exceeds_tolerance", result["candidateResidualReasons"])

    def test_eos_margin_failure_rejects(self):
        planet = self._planet_with_copy()
        planet.HPIceDiagnostics["VI"]["posthocProductionCandidate"]["posthocEOSPressureMargin_MPa"] = 0.0

        result = self._evaluate(planet)

        self.assertFalse(result["candidateResidualsPassed"])
        self.assertEqual(result["candidateResidualStatus"], "candidate_outside_eos_domain_rejected")
        self.assertIn("eos_margin_nonpositive", result["candidateResidualReasons"])

    def test_zero_temperature_contrast_rejects(self):
        planet = self._planet_with_copy()
        candidate = self._candidate(planet)
        candidate["candidateT_K"] = np.array([285.0, 285.0, 285.0])

        result = self._evaluate(planet)

        self.assertFalse(result["candidateResidualsPassed"])
        self.assertEqual(result["candidateResidualStatus"], "candidate_invalid_contrast_rejected")
        self.assertIn("candidate_temperature_contrast_not_positive", result["candidateResidualReasons"])

    def test_subcritical_rayleigh_rejects(self):
        planet = self._planet_with_copy()
        planet.HPIceDiagnostics["VI"]["RaConvect"] = 1.0e3
        planet.HPIceDiagnostics["VI"]["RaCrit"] = 1.0e8

        result = self._evaluate(planet)

        self.assertFalse(result["candidateResidualsPassed"])
        self.assertEqual(result["candidateResidualStatus"], "candidate_subcritical_rejected")
        self.assertIn("rayleigh_ratio_not_supercritical", result["candidateResidualReasons"])

    def test_nonfinite_viscosity_rejects(self):
        planet = self._planet_with_copy()
        candidate = self._candidate(planet)
        candidate["candidateEta_Pas"][0] = np.nan

        result = self._evaluate(planet)

        self.assertFalse(result["candidateResidualsPassed"])
        self.assertEqual(result["candidateResidualStatus"], "candidate_invalid_viscosity_rejected")
        self.assertIn("candidate_profile_viscosity_invalid", result["candidateResidualReasons"])

    def test_missing_required_field_rejects_with_reason(self):
        planet = self._planet_with_copy()
        candidate = self._candidate(planet)
        del candidate["candidateP_MPa"]

        result = self._evaluate(planet)

        self.assertFalse(result["candidateResidualsPassed"])
        self.assertEqual(result["candidateResidualStatus"], "candidate_missing_required_field_rejected")
        self.assertIn("missing_candidateP_MPa", result["candidateResidualReasons"])

    def test_candidate_booleans_remain_false(self):
        planet = self._planet_with_copy()
        candidate = self._candidate(planet)
        candidate["candidateAccepted"] = True
        candidate["candidateAppliedToProfile"] = True

        result = self._evaluate(planet)

        self.assertFalse(result["candidateAccepted"])
        self.assertFalse(result["candidateAppliedToProfile"])

    def test_residual_evaluation_does_not_modify_planet_fields(self):
        planet = self._planet_with_copy()
        before = Test11PosthocIceVIProductionCandidate._snapshot(planet)

        self._evaluate(planet)

        for key, value in before.items():
            if isinstance(value, np.ndarray):
                np.testing.assert_array_equal(getattr(planet, key), value)
            else:
                self.assertEqual(getattr(planet, key), value)


class Test15ActiveIceVICandidateRollbackPolicy(unittest.TestCase):

    def _planet_with_copy(self):
        return Test14ActiveIceVICandidateResidualEvaluator()._planet_with_copy()

    @staticmethod
    def _evaluate_residuals(planet):
        return EvaluateActiveIceVICandidateResiduals(
            planet, productionMode="Kalousova2018_production_experimental",
        )

    @staticmethod
    def _apply(planet):
        return ApplyActiveIceVIRollbackPolicy(
            planet, productionMode="Kalousova2018_production_experimental",
        )

    @staticmethod
    def _candidate(planet):
        return planet.HPIceDiagnostics["VI"]["activeProductionCandidate"]

    def _planet_with_failed_residuals(self, modifier):
        planet = self._planet_with_copy()
        modifier(planet)
        self._evaluate_residuals(planet)
        return planet

    def _full_chain(self, planet, modifier=None, skip_copy=False, skip_residuals=False):
        EvaluatePosthocIceVIProductionCandidate(
            planet, productionMode="Kalousova2018_production_experimental",
        )
        if not skip_copy:
            BuildActiveIceVIProductionCandidateCopy(
                planet, productionMode="Kalousova2018_production_experimental",
            )
        if modifier is not None:
            modifier(planet)
        if not skip_residuals:
            self._evaluate_residuals(planet)
        return self._apply(planet)

    def _assert_discarded(self, result, residual_status, reason):
        self.assertTrue(result["rollbackRequired"])
        self.assertTrue(result["candidateDiscarded"])
        self.assertTrue(result["rollbackApplied"])
        self.assertEqual(result["rollbackStatus"], "rollback_candidate_discarded")
        self.assertEqual(result["rollbackReason"], residual_status)
        self.assertIn(reason, result["rollbackReasons"])
        self.assertFalse(result["candidateAccepted"])
        self.assertFalse(result["candidateAppliedToProfile"])

    def test_rollback_helper_handles_missing_candidate_copy(self):
        planet = Test13ActiveIceVICandidateProfileCopy()._planet()

        result = self._apply(planet)

        self.assertEqual(result["rollbackStatus"], "rollback_missing_candidate")
        self.assertTrue(result["rollbackRequired"])
        self.assertTrue(result["candidateDiscarded"])
        self.assertFalse(result["rollbackApplied"])
        self.assertEqual(result["rollbackReason"], "missing_candidate_copy")
        self.assertIn("missing_candidate_copy", result["rollbackReasons"])
        self.assertFalse(result["candidateAccepted"])
        self.assertFalse(result["candidateAppliedToProfile"])

    def test_clean_residuals_do_not_require_rollback(self):
        planet = self._planet_with_copy()
        self._evaluate_residuals(planet)

        result = self._apply(planet)

        self.assertEqual(result["rollbackStatus"], "rollback_not_required")
        self.assertFalse(result["rollbackRequired"])
        self.assertFalse(result["candidateDiscarded"])
        self.assertFalse(result["rollbackApplied"])
        self.assertIsNone(result["rollbackReason"])
        self.assertEqual(result["rollbackReasons"], ())
        self.assertFalse(result["candidateAccepted"])
        self.assertFalse(result["candidateAppliedToProfile"])

    def test_energy_residual_failure_discards_candidate(self):
        planet = self._planet_with_failed_residuals(
            lambda p: p.HPIceDiagnostics["VI"].update({
                "energyResidual_W": 1.0e9,
                "energyResidual_frac": 1.0e-3,
            })
        )

        result = self._apply(planet)

        self._assert_discarded(
            result, "candidate_energy_residual_rejected",
            "energy_residual_exceeds_tolerance",
        )

    def test_heat_flux_residual_failure_discards_candidate(self):
        def modifier(planet):
            self._candidate(planet)["candidateq_out_Wm2"] = 0.021

        planet = self._planet_with_failed_residuals(modifier)

        result = self._apply(planet)

        self._assert_discarded(
            result, "candidate_heat_flux_residual_rejected",
            "spherical_flux_scaling_mismatch",
        )

    def test_mass_residual_failure_discards_candidate(self):
        planet = self._planet_with_failed_residuals(
            lambda p: p.HPIceDiagnostics["VI"].update({"massResidual_kg": 1.0})
        )

        result = self._apply(planet)

        self._assert_discarded(
            result, "candidate_mass_residual_rejected",
            "mass_residual_nonzero",
        )

    def test_phase_boundary_failure_discards_candidate(self):
        def modifier(planet):
            planet.HPIceDiagnostics["VI"]["posthocProductionCandidate"]["posthocPhaseBoundaryResidual_K"] = 0.2

        planet = self._planet_with_failed_residuals(modifier)

        result = self._apply(planet)

        self._assert_discarded(
            result, "candidate_phase_boundary_rejected",
            "phase_boundary_residual_exceeds_tolerance",
        )

    def test_layer_closure_failure_discards_candidate(self):
        planet = self._planet_with_failed_residuals(
            lambda p: p.HPIceDiagnostics["VI"].update({"layerClosureResidual_m": 1.0})
        )

        result = self._apply(planet)

        self._assert_discarded(
            result, "candidate_layer_closure_rejected",
            "layer_closure_residual_exceeds_tolerance",
        )

    def test_eos_domain_failure_discards_candidate(self):
        def modifier(planet):
            planet.HPIceDiagnostics["VI"]["posthocProductionCandidate"]["posthocEOSPressureMargin_MPa"] = 0.0

        planet = self._planet_with_failed_residuals(modifier)

        result = self._apply(planet)

        self._assert_discarded(
            result, "candidate_outside_eos_domain_rejected",
            "eos_margin_nonpositive",
        )

    def test_invalid_contrast_discards_candidate(self):
        def modifier(planet):
            self._candidate(planet)["candidateT_K"] = np.array([285.0, 285.0, 285.0])

        planet = self._planet_with_failed_residuals(modifier)

        result = self._apply(planet)

        self._assert_discarded(
            result, "candidate_invalid_contrast_rejected",
            "candidate_temperature_contrast_not_positive",
        )

    def test_subcritical_rayleigh_discards_candidate(self):
        planet = self._planet_with_failed_residuals(
            lambda p: p.HPIceDiagnostics["VI"].update({"RaConvect": 1.0e3, "RaCrit": 1.0e8})
        )

        result = self._apply(planet)

        self._assert_discarded(
            result, "candidate_subcritical_rejected",
            "rayleigh_ratio_not_supercritical",
        )

    def test_invalid_viscosity_discards_candidate(self):
        def modifier(planet):
            self._candidate(planet)["candidateEta_Pas"][0] = np.nan

        planet = self._planet_with_failed_residuals(modifier)

        result = self._apply(planet)

        self._assert_discarded(
            result, "candidate_invalid_viscosity_rejected",
            "candidate_profile_viscosity_invalid",
        )

    def test_missing_required_field_discards_candidate(self):
        def modifier(planet):
            del self._candidate(planet)["candidateP_MPa"]

        planet = self._planet_with_failed_residuals(modifier)

        result = self._apply(planet)

        self._assert_discarded(
            result, "candidate_missing_required_field_rejected",
            "missing_candidateP_MPa",
        )

    def test_protected_fields_changed_blocks_rollback_application(self):
        planet = self._planet_with_copy()
        self._evaluate_residuals(planet)
        self._candidate(planet)["protectedFieldsUnchanged"] = False

        result = self._apply(planet)

        self.assertEqual(result["rollbackStatus"], "rollback_protected_fields_changed")
        self.assertTrue(result["rollbackRequired"])
        self.assertTrue(result["candidateDiscarded"])
        self.assertFalse(result["rollbackApplied"])
        self.assertEqual(result["rollbackReason"], "protected_fields_changed")
        self.assertIn("protected_fields_changed", result["rollbackReasons"])
        self.assertFalse(result["candidateAccepted"])
        self.assertFalse(result["candidateAppliedToProfile"])

    def test_repeated_rollback_evaluation_is_idempotent(self):
        planet = self._planet_with_failed_residuals(
            lambda p: p.HPIceDiagnostics["VI"].update({"massResidual_kg": 1.0})
        )

        first = self._apply(planet).copy()
        second = self._apply(planet)

        for key in (
            "candidateDiscarded",
            "candidateDiscardReason",
            "rollbackStatus",
            "rollbackRequired",
            "rollbackApplied",
            "rollbackReason",
            "rollbackReasons",
            "candidateAccepted",
            "candidateAppliedToProfile",
        ):
            self.assertEqual(second[key], first[key])

    def test_discarded_residual_failures_are_terminal_on_full_chain_rerun(self):
        cases = (
            (
                "energy",
                lambda p: p.HPIceDiagnostics["VI"].update({
                    "energyResidual_W": 1.0e9,
                    "energyResidual_frac": 1.0e-3,
                }),
                "candidate_energy_residual_rejected",
                "energy_residual_exceeds_tolerance",
            ),
            (
                "mass",
                lambda p: p.HPIceDiagnostics["VI"].update({"massResidual_kg": 1.0}),
                "candidate_mass_residual_rejected",
                "mass_residual_nonzero",
            ),
            (
                "layer_closure",
                lambda p: p.HPIceDiagnostics["VI"].update({"layerClosureResidual_m": 1.0}),
                "candidate_layer_closure_rejected",
                "layer_closure_residual_exceeds_tolerance",
            ),
            (
                "subcritical",
                lambda p: p.HPIceDiagnostics["VI"].update({"RaConvect": 1.0e3, "RaCrit": 1.0e8}),
                "candidate_subcritical_rejected",
                "rayleigh_ratio_not_supercritical",
            ),
        )
        for label, modifier, residual_status, reason in cases:
            with self.subTest(label=label):
                planet = Test13ActiveIceVICandidateProfileCopy()._planet()
                before = Test11PosthocIceVIProductionCandidate._snapshot(planet)

                first = self._full_chain(planet, modifier)
                second = self._full_chain(planet, modifier)

                self.assertEqual(first["candidateResidualStatus"], residual_status)
                self.assertEqual(second["candidateResidualStatus"], residual_status)
                self.assertEqual(second["rollbackStatus"], "rollback_candidate_discarded")
                self.assertEqual(second["rollbackReasons"], first["rollbackReasons"])
                self.assertIn(reason, second["rollbackReasons"])
                self.assertFalse(second["candidateAccepted"])
                self.assertFalse(second["candidateAppliedToProfile"])
                for key, value in before.items():
                    if isinstance(value, np.ndarray):
                        np.testing.assert_array_equal(getattr(planet, key), value)
                    else:
                        self.assertEqual(getattr(planet, key), value)

    def test_missing_candidate_copy_is_terminal_on_full_chain_rerun(self):
        planet = Test13ActiveIceVICandidateProfileCopy()._planet()
        before = Test11PosthocIceVIProductionCandidate._snapshot(planet)

        first = self._full_chain(planet, skip_copy=True, skip_residuals=True)
        second = self._full_chain(planet, skip_copy=True, skip_residuals=True)

        self.assertEqual(first["rollbackStatus"], "rollback_missing_candidate")
        self.assertEqual(second["rollbackStatus"], "rollback_missing_candidate")
        self.assertEqual(second["rollbackReasons"], first["rollbackReasons"])
        self.assertFalse(second["candidateAccepted"])
        self.assertFalse(second["candidateAppliedToProfile"])
        for key, value in before.items():
            if isinstance(value, np.ndarray):
                np.testing.assert_array_equal(getattr(planet, key), value)
            else:
                self.assertEqual(getattr(planet, key), value)

    def test_candidate_booleans_remain_false(self):
        planet = self._planet_with_copy()
        self._evaluate_residuals(planet)
        candidate = self._candidate(planet)
        candidate["candidateAccepted"] = True
        candidate["candidateAppliedToProfile"] = True

        result = self._apply(planet)

        self.assertFalse(result["candidateAccepted"])
        self.assertFalse(result["candidateAppliedToProfile"])

    def test_rollback_policy_does_not_modify_planet_fields(self):
        planet = self._planet_with_failed_residuals(
            lambda p: p.HPIceDiagnostics["VI"].update({"massResidual_kg": 1.0})
        )
        before = Test11PosthocIceVIProductionCandidate._snapshot(planet)

        self._apply(planet)

        for key, value in before.items():
            if isinstance(value, np.ndarray):
                np.testing.assert_array_equal(getattr(planet, key), value)
            else:
                self.assertEqual(getattr(planet, key), value)


class Test16ActiveIceVIThermalUpdateCandidate(unittest.TestCase):

    def _planet_with_residuals(self):
        planet = Test14ActiveIceVICandidateResidualEvaluator()._planet_with_copy()
        EvaluateActiveIceVICandidateResiduals(
            planet, productionMode="Kalousova2018_production_experimental",
        )
        return planet

    @staticmethod
    def _build(planet):
        return BuildActiveIceVIThermalUpdateCandidate(
            planet, productionMode="Kalousova2018_production_experimental",
        )

    @staticmethod
    def _candidate(planet):
        return planet.HPIceDiagnostics["VI"]["activeProductionCandidate"]

    def test_nominal_candidate_creates_independent_temperature_copy(self):
        planet = self._planet_with_residuals()
        candidate = self._candidate(planet)

        result = self._build(planet)

        self.assertEqual(result["candidateThermalUpdateStatus"], "candidate_thermal_update_no_op_reconstruction")
        self.assertTrue(result["candidateThermalUpdateAppliedToCopy"])
        self.assertFalse(result["candidateThermalUpdateAppliedToPlanet"])
        self.assertFalse(result["candidateThermalUpdateAccepted"])
        np.testing.assert_allclose(result["candidateUpdatedT_K"], candidate["candidateT_K"])
        self.assertFalse(np.shares_memory(result["candidateUpdatedT_K"], candidate["candidateT_K"]))
        self.assertFalse(np.shares_memory(result["candidateUpdatedT_K"], planet.T_K))

    def test_updated_heat_power_and_flux_are_spherical_copy_metadata(self):
        planet = self._planet_with_residuals()
        candidate = self._candidate(planet)

        result = self._build(planet)

        expected_Q_W = np.full(candidate["candidateT_K"].shape, candidate["candidateQ_in_W"])
        expected_q_Wm2 = expected_Q_W / (4.0 * np.pi * candidate["candidateR_m"]**2)
        np.testing.assert_allclose(result["candidateUpdatedQ_W"], expected_Q_W)
        np.testing.assert_allclose(result["candidateUpdatedq_Wm2"], expected_q_Wm2)
        self.assertEqual(result["candidateThermalHeatPowerResidual_W"], 0.0)

    def test_missing_candidate_rejects(self):
        planet = Test13ActiveIceVICandidateProfileCopy()._planet()

        result = self._build(planet)

        self.assertEqual(result["candidateThermalUpdateStatus"], "candidate_thermal_update_missing_candidate_rejected")
        self.assertFalse(result["candidateThermalUpdateAppliedToCopy"])
        self.assertFalse(result["candidateThermalUpdateAppliedToPlanet"])
        self.assertFalse(result["candidateThermalUpdateAccepted"])

    def test_discarded_candidate_rejects(self):
        planet = self._planet_with_residuals()
        candidate = self._candidate(planet)
        candidate["candidateResidualsPassed"] = False
        candidate["candidateResidualStatus"] = "candidate_mass_residual_rejected"
        candidate["candidateResidualReasons"] = ("mass_residual_nonzero",)
        ApplyActiveIceVIRollbackPolicy(
            planet, productionMode="Kalousova2018_production_experimental",
        )

        result = self._build(planet)

        self.assertEqual(result["candidateThermalUpdateStatus"], "candidate_thermal_update_discarded_candidate_rejected")
        self.assertFalse(result["candidateThermalUpdateAppliedToCopy"])
        self.assertFalse(result["candidateThermalUpdateAppliedToPlanet"])
        self.assertFalse(result["candidateThermalUpdateAccepted"])
        self.assertEqual(candidate["rollbackStatus"], "rollback_candidate_discarded")

    def test_residuals_not_passed_candidate_rejects(self):
        planet = Test14ActiveIceVICandidateResidualEvaluator()._planet_with_copy()

        result = self._build(planet)

        self.assertEqual(result["candidateThermalUpdateStatus"], "candidate_thermal_update_residuals_not_passed_rejected")
        self.assertFalse(result["candidateThermalUpdateAppliedToCopy"])

    def test_missing_radius_rejects(self):
        planet = self._planet_with_residuals()
        candidate = self._candidate(planet)
        candidate["candidateR_m"] = np.array([], dtype=float)

        result = self._build(planet)

        self.assertEqual(result["candidateThermalUpdateStatus"], "candidate_thermal_update_missing_radius_rejected")
        self.assertFalse(result["candidateThermalUpdateAppliedToCopy"])
        self.assertIn(
            "candidate_radius_required_for_spherical_flux_reconstruction",
            result["candidateThermalUpdateReasons"],
        )

    def test_nonfinite_candidate_temperature_rejects(self):
        planet = self._planet_with_residuals()
        candidate = self._candidate(planet)
        candidate["candidateT_K"][0] = np.nan

        result = self._build(planet)

        self.assertEqual(result["candidateThermalUpdateStatus"], "candidate_thermal_update_nonfinite_rejected")
        self.assertFalse(result["candidateThermalUpdateAppliedToCopy"])

    def test_update_and_candidate_booleans_remain_false(self):
        planet = self._planet_with_residuals()
        candidate = self._candidate(planet)
        candidate["candidateAccepted"] = True
        candidate["candidateAppliedToProfile"] = True

        result = self._build(planet)

        self.assertFalse(result["candidateThermalUpdateAppliedToPlanet"])
        self.assertFalse(result["candidateThermalUpdateAccepted"])
        self.assertFalse(candidate["candidateAccepted"])
        self.assertFalse(candidate["candidateAppliedToProfile"])

    def test_thermal_update_does_not_modify_planet_fields(self):
        planet = self._planet_with_residuals()
        before = Test11PosthocIceVIProductionCandidate._snapshot(planet)

        self._build(planet)

        for key, value in before.items():
            if isinstance(value, np.ndarray):
                np.testing.assert_array_equal(getattr(planet, key), value)
            else:
                self.assertEqual(getattr(planet, key), value)


class Test17ActiveIceVILinearThermalUpdateCandidate(unittest.TestCase):

    LINEAR_STRATEGY = "linear_conservative_reconstruction"

    def _planet_with_residuals(self):
        return Test16ActiveIceVIThermalUpdateCandidate()._planet_with_residuals()

    @staticmethod
    def _build(planet):
        return BuildActiveIceVIThermalUpdateCandidate(
            planet,
            productionMode="Kalousova2018_production_experimental",
            strategy=Test17ActiveIceVILinearThermalUpdateCandidate.LINEAR_STRATEGY,
        )

    @staticmethod
    def _candidate(planet):
        return planet.HPIceDiagnostics["VI"]["activeProductionCandidate"]

    def _safe_linear_candidate(self):
        planet = self._planet_with_residuals()
        candidate = self._candidate(planet)
        candidate["candidateQ_in_W"] = 1.0e7
        candidate["candidateQ_out_W"] = 1.0e7
        return planet

    def test_no_op_strategy_still_works(self):
        planet = self._planet_with_residuals()

        result = BuildActiveIceVIThermalUpdateCandidate(
            planet, productionMode="Kalousova2018_production_experimental",
        )

        self.assertEqual(result["candidateThermalUpdateStrategy"], "no_op_reconstruction")
        self.assertEqual(result["candidateThermalUpdateStatus"], "candidate_thermal_update_no_op_reconstruction")
        self.assertTrue(result["candidateThermalUpdateAppliedToCopy"])

    def test_linear_conservative_strategy_reconstructs_copy_only_profile(self):
        planet = self._safe_linear_candidate()
        candidate = self._candidate(planet)

        result = self._build(planet)

        self.assertEqual(result["candidateThermalUpdateStrategy"], self.LINEAR_STRATEGY)
        self.assertEqual(
            result["candidateThermalUpdateStatus"],
            "candidate_thermal_update_linear_conservative_reconstruction",
        )
        self.assertTrue(result["candidateThermalUpdateAppliedToCopy"])
        self.assertFalse(result["candidateThermalUpdateAppliedToPlanet"])
        self.assertFalse(result["candidateThermalUpdateAccepted"])
        self.assertFalse(candidate["candidateAccepted"])
        self.assertFalse(candidate["candidateAppliedToProfile"])
        self.assertFalse(np.shares_memory(result["candidateUpdatedT_K"], candidate["candidateT_K"]))
        self.assertFalse(np.shares_memory(result["candidateUpdatedT_K"], planet.T_K))
        self.assertTrue(np.all(np.isfinite(result["candidateUpdatedT_K"])))
        np.testing.assert_array_equal(result["candidateUpdatedPhaseArray"], np.array([6, 6, 6]))

    def test_linear_strategy_conserves_heat_power_and_uses_spherical_flux(self):
        planet = self._safe_linear_candidate()
        candidate = self._candidate(planet)

        result = self._build(planet)

        expected_Q_W = np.full(candidate["candidateT_K"].shape, 1.0e7)
        expected_q_Wm2 = expected_Q_W / (4.0 * np.pi * candidate["candidateR_m"]**2)
        np.testing.assert_allclose(result["candidateUpdatedQ_W"], expected_Q_W)
        np.testing.assert_allclose(result["candidateUpdatedq_Wm2"], expected_q_Wm2)
        self.assertEqual(result["candidateThermalHeatPowerResidual_W"], 0.0)

    def test_linear_strategy_rejects_missing_conductivity(self):
        planet = self._safe_linear_candidate()
        candidate = self._candidate(planet)
        candidate["candidatekTherm_WmK"] = np.array([], dtype=float)

        result = self._build(planet)

        self.assertEqual(
            result["candidateThermalUpdateStatus"],
            "candidate_thermal_update_missing_conductivity_rejected",
        )
        self.assertFalse(result["candidateThermalUpdateAppliedToCopy"])
        self.assertIn(
            "candidate_thermal_conductivity_and_geometry_required_for_linear_reconstruction",
            result["candidateThermalUpdateReasons"],
        )

    def test_linear_strategy_rejects_phase_boundary_crossing(self):
        planet = self._planet_with_residuals()

        result = self._build(planet)

        self.assertEqual(
            result["candidateThermalUpdateStatus"],
            "candidate_thermal_update_phase_boundary_rejected",
        )
        self.assertFalse(result["candidateThermalUpdateAppliedToCopy"])
        self.assertFalse(result["candidateThermalUpdateAppliedToPlanet"])
        self.assertFalse(result["candidateThermalUpdateAccepted"])

    def test_discarded_candidate_rejects_linear_strategy(self):
        planet = self._safe_linear_candidate()
        candidate = self._candidate(planet)
        candidate["candidateResidualsPassed"] = False
        candidate["candidateResidualStatus"] = "candidate_mass_residual_rejected"
        candidate["candidateResidualReasons"] = ("mass_residual_nonzero",)
        ApplyActiveIceVIRollbackPolicy(
            planet, productionMode="Kalousova2018_production_experimental",
        )

        result = self._build(planet)

        self.assertEqual(
            result["candidateThermalUpdateStatus"],
            "candidate_thermal_update_discarded_candidate_rejected",
        )
        self.assertFalse(result["candidateThermalUpdateAppliedToCopy"])
        self.assertFalse(result["candidateThermalUpdateAppliedToPlanet"])

    def test_linear_strategy_does_not_modify_planet_fields(self):
        planet = self._safe_linear_candidate()
        before = Test11PosthocIceVIProductionCandidate._snapshot(planet)

        self._build(planet)

        for key, value in before.items():
            if isinstance(value, np.ndarray):
                np.testing.assert_array_equal(getattr(planet, key), value)
            else:
                self.assertEqual(getattr(planet, key), value)


class Test18ActiveIceVIOriginalBoundaryThermalReconstruction(unittest.TestCase):

    ORIGINAL_BOUNDARY_STRATEGY = "original_boundary_reconstruction"

    def _planet_with_residuals(self):
        return Test16ActiveIceVIThermalUpdateCandidate()._planet_with_residuals()

    @staticmethod
    def _build(planet):
        return BuildActiveIceVIThermalUpdateCandidate(
            planet,
            productionMode="Kalousova2018_production_experimental",
            strategy=Test18ActiveIceVIOriginalBoundaryThermalReconstruction.ORIGINAL_BOUNDARY_STRATEGY,
        )

    @staticmethod
    def _candidate(planet):
        return planet.HPIceDiagnostics["VI"]["activeProductionCandidate"]

    def test_no_op_strategy_still_works(self):
        planet = self._planet_with_residuals()

        result = BuildActiveIceVIThermalUpdateCandidate(
            planet, productionMode="Kalousova2018_production_experimental",
        )

        self.assertEqual(result["candidateThermalUpdateStatus"], "candidate_thermal_update_no_op_reconstruction")
        self.assertTrue(result["candidateThermalUpdateAppliedToCopy"])

    def test_linear_strategy_remains_available_and_fail_closed(self):
        planet = self._planet_with_residuals()

        result = BuildActiveIceVIThermalUpdateCandidate(
            planet,
            productionMode="Kalousova2018_production_experimental",
            strategy="linear_conservative_reconstruction",
        )

        self.assertEqual(result["candidateThermalUpdateStatus"], "candidate_thermal_update_phase_boundary_rejected")
        self.assertFalse(result["candidateThermalUpdateAppliedToCopy"])
        self.assertFalse(result["candidateThermalUpdateAppliedToPlanet"])

    def test_original_boundary_preserves_top_and_bottom_temperatures(self):
        planet = self._planet_with_residuals()
        candidate = self._candidate(planet)

        result = self._build(planet)

        self.assertEqual(
            result["candidateThermalUpdateStatus"],
            "candidate_thermal_update_original_boundary_reconstruction",
        )
        self.assertTrue(result["candidateThermalUpdateAppliedToCopy"])
        self.assertEqual(result["candidateThermalBoundaryCondition"], "preserve_finalized_top_bottom_temperature")
        self.assertEqual(result["candidateThermalInterpolationCoordinate"], "candidateZ_m")
        self.assertAlmostEqual(result["candidateUpdatedT_K"][0], candidate["candidateT_K"][0])
        self.assertAlmostEqual(result["candidateUpdatedT_K"][-1], candidate["candidateT_K"][-1])

    def test_original_boundary_updated_temperature_is_independent_copy(self):
        planet = self._planet_with_residuals()
        candidate = self._candidate(planet)

        result = self._build(planet)

        self.assertFalse(np.shares_memory(result["candidateUpdatedT_K"], candidate["candidateT_K"]))
        self.assertFalse(np.shares_memory(result["candidateUpdatedT_K"], planet.T_K))
        result["candidateUpdatedT_K"][0] += 1.0
        self.assertNotEqual(result["candidateUpdatedT_K"][0], candidate["candidateT_K"][0])
        self.assertNotEqual(result["candidateUpdatedT_K"][0], planet.T_K[1])

    def test_original_boundary_interpolates_with_depth_grid(self):
        planet = self._planet_with_residuals()
        candidate = self._candidate(planet)
        candidate["candidateT_K"] = np.array([285.0, 285.5, 287.0])
        expected = np.interp(
            candidate["candidateZ_m"],
            (candidate["candidateZ_m"][0], candidate["candidateZ_m"][-1]),
            (candidate["candidateT_K"][0], candidate["candidateT_K"][-1]),
        )

        result = self._build(planet)

        np.testing.assert_allclose(result["candidateUpdatedT_K"], expected)
        self.assertEqual(result["candidateThermalInterpolationCoordinate"], "candidateZ_m")

    def test_original_boundary_missing_or_invalid_depth_rejects(self):
        planet = self._planet_with_residuals()
        candidate = self._candidate(planet)
        candidate["candidateZ_m"] = np.array([1.0, 2.0, 1.5])

        result = self._build(planet)

        self.assertEqual(
            result["candidateThermalUpdateStatus"],
            "candidate_thermal_update_original_boundary_missing_depth_rejected",
        )
        self.assertFalse(result["candidateThermalUpdateAppliedToCopy"])
        self.assertIn("candidate_depth_grid_must_increase_downward", result["candidateThermalUpdateReasons"])

    def test_original_boundary_q_metadata_uses_spherical_scaling(self):
        planet = self._planet_with_residuals()
        candidate = self._candidate(planet)

        result = self._build(planet)

        expected_Q_W = np.full(candidate["candidateT_K"].shape, candidate["candidateQ_in_W"])
        expected_q_Wm2 = expected_Q_W / (4.0 * np.pi * candidate["candidateR_m"]**2)
        np.testing.assert_allclose(result["candidateUpdatedQ_W"], expected_Q_W)
        np.testing.assert_allclose(result["candidateUpdatedq_Wm2"], expected_q_Wm2)
        self.assertEqual(result["candidateThermalHeatPowerResidual_W"], 0.0)

    def test_original_boundary_eos_and_phase_boundary_checks_pass_for_valid_candidate(self):
        planet = self._planet_with_residuals()

        result = self._build(planet)

        self.assertEqual(result["candidateThermalMeltCurveStatus"], "melt_curve_ok")
        self.assertEqual(result["candidateThermalPhaseBoundaryStatus"], "phase_boundary_ok")
        self.assertGreater(result["candidateThermalMinPhaseBoundaryMargin_K"], 0.0)

    def test_original_boundary_phase_crossing_rejects(self):
        planet = self._planet_with_residuals()
        candidate = self._candidate(planet)
        candidate["candidateT_K"] = np.array([285.0, 286.0, 300.0])

        result = self._build(planet)

        self.assertEqual(
            result["candidateThermalUpdateStatus"],
            "candidate_thermal_update_original_boundary_wrong_phase_rejected",
        )
        self.assertFalse(result["candidateThermalUpdateAppliedToCopy"])
        self.assertFalse(result["candidateThermalUpdateAppliedToPlanet"])

    def test_original_boundary_missing_radius_rejects(self):
        planet = self._planet_with_residuals()
        candidate = self._candidate(planet)
        candidate["candidateR_m"] = np.array([], dtype=float)

        result = self._build(planet)

        self.assertEqual(result["candidateThermalUpdateStatus"], "candidate_thermal_update_missing_radius_rejected")
        self.assertFalse(result["candidateThermalUpdateAppliedToCopy"])

    def test_original_boundary_booleans_remain_false(self):
        planet = self._planet_with_residuals()
        candidate = self._candidate(planet)
        candidate["candidateAccepted"] = True
        candidate["candidateAppliedToProfile"] = True

        result = self._build(planet)

        self.assertFalse(result["candidateThermalUpdateAppliedToPlanet"])
        self.assertFalse(result["candidateThermalUpdateAccepted"])
        self.assertFalse(candidate["candidateAccepted"])
        self.assertFalse(candidate["candidateAppliedToProfile"])

    def test_original_boundary_does_not_modify_planet_fields(self):
        planet = self._planet_with_residuals()
        before = Test11PosthocIceVIProductionCandidate._snapshot(planet)

        self._build(planet)

        for key, value in before.items():
            if isinstance(value, np.ndarray):
                np.testing.assert_array_equal(getattr(planet, key), value)
            else:
                self.assertEqual(getattr(planet, key), value)


if __name__ == "__main__":
    unittest.main()
