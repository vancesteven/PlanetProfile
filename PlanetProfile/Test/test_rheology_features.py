import unittest

import numpy as np

from PlanetProfile.Thermodynamics.HydroEOS import (
    GetIceArrheniusViscosityKwargs,
    GetIceEOS,
    ViscIceArrhenius_Pas,
)
from PlanetProfile.Thermodynamics.LayerPropagators import (
    _ConvectionDeschampsSotinHPIceDiagnostic,
    _FixedPhaseEOS,
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

    def test_gated_dry_run_records_ice_vi_candidate_without_profile_mutation(self):
        planet = self._synthetic_hp_planet(qBot_Wm2=1.0e-3)
        before = self._profile_snapshot(planet)

        HPIceConvectionDiagnostics(planet, _Namespace(), ResolveHPIceConvectionModel(planet))

        state = planet.HPIceDiagnostics["VI"]
        self.assertTrue(state["productionCandidate"])
        self.assertEqual(state["productionMode"], "Kalousova2018_production_experimental")
        self.assertEqual(state["candidateStatus"], "ice_vi_candidate")
        self.assertFalse(state["updateAccepted"])
        self.assertEqual(state["candidateReason"], "dry_run_only")
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
            thickness_m=2.0e5,
            Ttop_K=280.0,
            Tbot_K=300.0,
            Tconv_K=280.0,
            etaConv_Pas=5.0e14,
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


if __name__ == "__main__":
    unittest.main()
