import numpy as np
import logging as log
from Utilities.defineStructs import Constants
from Thermodynamics.FromLiterature.HydroEOS import PhaseConv
from Thermodynamics.FromLiterature.ThermalProfiles import ConductiveTemperature

def EvalLayerProperties(Planet, Params, iStart, iEnd, EOS, P_MPa, T_K):

    Planet.rhoMatrix_kgm3[iStart:iEnd] = EOS.fn_rho_kgm3(  P_MPa, T_K, grid=False)
    Planet.Cp_JkgK[iStart:iEnd] =        EOS.fn_Cp_JkgK(   P_MPa, T_K, grid=False)
    Planet.alpha_pK[iStart:iEnd] =       EOS.fn_alpha_pK(  P_MPa, T_K, grid=False)
    Planet.kTherm_WmK[iStart:iEnd] =     EOS.fn_kTherm_WmK(P_MPa, T_K, grid=False)

    return Planet


def PorosityCorrectionVacIce(Planet, Params, iStart, iEnd, EOS, P_MPa, T_K):

    Planet.phi_frac[iStart:iEnd] = EOS.fn_phi_frac(P_MPa, T_K)
    # Only evaluate at places where porosity is non-negligible, so that we
    # can limit calling EOS functions with invalid parameters, which can give NaNs
    iEval = np.where(Planet.phi_frac[iStart:iEnd] >= Planet.Ocean.phiMin_frac)[0]
    iSolid = np.where(Planet.phi_frac[iStart:iEnd] < Planet.Ocean.phiMin_frac)[0]

    Planet.rho_kgm3[iStart:iEnd][iSolid] = Planet.rhoMatrix_kgm3[iStart:iEnd][iSolid]
    Planet.rho_kgm3[iStart:iEnd][iEval] = EOS.fn_porosCorrect(
        Planet.rhoMatrix_kgm3[iStart:iEnd][iEval], 0,
        Planet.phi_frac[iStart:iEnd][iEval], Planet.Ocean.Jrho)
    Planet.Cp_JkgK[iStart:iEnd][iEval] = EOS.fn_porosCorrect(
        Planet.Cp_JkgK[iStart:iEnd][iEval] * Planet.rhoMatrix_kgm3[iStart:iEnd][iEval], 0,
        Planet.phi_frac[iStart:iEnd][iEval], Planet.Ocean.JCp) / Planet.rho_kgm3[iStart:iEnd][iEval]
    Planet.alpha_pK[iStart:iEnd][iEval] = EOS.fn_porosCorrect(
        Planet.alpha_pK[iStart:iEnd][iEval], 0,
        Planet.phi_frac[iStart:iEnd][iEval], Planet.Ocean.Jalpha)
    Planet.kTherm_WmK[iStart:iEnd][iEval] = EOS.fn_porosCorrect(
        Planet.kTherm_WmK[iStart:iEnd][iEval], 0,
        Planet.phi_frac[iStart:iEnd][iEval], Planet.Ocean.JkTherm)

    return Planet


def PorosityCorrectionFilledIce(Planet, Params, iStart, iEnd, EOS, EOSpore):

    alphaPeff = Planet.Ocean.alphaPeff[EOS.phaseStr]
    DeltaPpore_MPa = 0.0
    if iStart == 0:
        # Record last entry in Ppore, and set it as needed to ensure
        # Ppore[i-1] reference will work correctly in loop
        PporeSave = Planet.Ppore_MPa[-1] + 0.0
        Planet.Ppore_MPa[-1] = Planet.P_MPa[0] + 0.0
    elif Planet.Ppore_MPa[iStart-1] == 0:
        # Make sure Ppore[i-1] has been assigned, and set it so Ppore[iStart] = P_MPa[iStart],
        # because we are at the top of a porous layer, above which is liquid or vacuum
        Planet.Ppore_MPa[iStart-1] = Planet.P_MPa[iStart] + 0.0

    for i in range(iStart, iEnd):
        Planet.Ppore_MPa[i] = Planet.Ppore_MPa[i-1] + DeltaPpore_MPa
        Peff_MPa = Planet.P_MPa[i] - alphaPeff * Planet.Ppore_MPa[i]
        Planet.phi_frac[i] = EOS.fn_phi_frac(Peff_MPa, Planet.T_K[i])
        if Planet.phi_frac[i] >= Planet.Ocean.phiMin_frac:
            # In ices, assume pores are filled only with liquid.
            rhoPore_kgm3 =   EOSpore.fn_rho_kgm3(  Planet.Ppore_MPa[i], Planet.T_K[i], grid=False)
            CpPore_JkgK =    EOSpore.fn_Cp_JkgK(   Planet.Ppore_MPa[i], Planet.T_K[i], grid=False)
            alphaPore_pK =   EOSpore.fn_alpha_pK(  Planet.Ppore_MPa[i], Planet.T_K[i], grid=False)
            kThermPore_WmK = EOSpore.fn_kTherm_WmK(Planet.Ppore_MPa[i], Planet.T_K[i], grid=False)

            Planet.rho_kgm3[i] = EOS.fn_porosCorrect(Planet.rhoMatrix_kgm3[i], rhoPore_kgm3,
                                                       Planet.phi_frac[i], Planet.Ocean.Jrho)
            Planet.Cp_JkgK[i] = EOS.fn_porosCorrect(Planet.Cp_JkgK[i] * Planet.rhoMatrix_kgm3[i], CpPore_JkgK * rhoPore_kgm3,
                                                       Planet.phi_frac[i], Planet.Ocean.JCp) / Planet.rho_kgm3[i]
            Planet.alpha_pK[i] = EOS.fn_porosCorrect(Planet.alpha_pK[i], alphaPore_pK,
                                                       Planet.phi_frac[i], Planet.Ocean.Jalpha)
            Planet.kTherm_WmK[i] = EOS.fn_porosCorrect(Planet.kTherm_WmK[i], kThermPore_WmK,
                                                       Planet.phi_frac[i], Planet.Ocean.JkTherm)
            # Get approximate pore pressure change, using the previous depth change. This value
            # should be close to the correct value -- we do it this way to avoid array indexing
            # issues at the start/end points and the need to calculate the depth change twice
            DeltaPpore_MPa = 1e-6 * rhoPore_kgm3 * Planet.g_ms2[i] * (Planet.z_m[i] - Planet.z_m[i-1])
        else:
            Planet.phi_frac[i] = 0.0
            Planet.rho_kgm3[i] = Planet.rhoMatrix_kgm3[i] + 0.0
            Planet.Ppore_MPa[i] = Planet.P_MPa[i] + 0.0
            DeltaPpore_MPa = Planet.P_MPa[i] - Planet.P_MPa[i-1]

    # Restore final Ppore entry if it was beyond those affected by the loop,
    # because we used it as the -1 index to start the loop
    if iStart == 0 and iEnd < np.size(Planet.Ppore_MPa):
        Planet.Ppore_MPa[-1] = PporeSave

    return Planet, phasePore


def PropagateConduction(Planet, Params, iStart, iEnd):

    thisMAbove_kg = np.sum(Planet.MLayer_kg[:iStart])
    for i in range(iStart+1, iEnd+1):
        # Increment depth based on change in pressure, combined with gravity and density
        Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1]) * 1e6 / Planet.g_ms2[i-1] / \
                        Planet.rho_kgm3[i-1]
        Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
        Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1] ** 3 - Planet.r_m[i] ** 3)
        thisMAbove_kg += Planet.MLayer_kg[i-1]
        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i] ** 2
        log.debug(f'il: {i:d}; P_MPa: {Planet.P_MPa[i]:.3f}; T_K: {Planet.T_K[i]:.3f}; phase: {Planet.phase[i]:d}')

    return Planet


def PropagateAdiabaticSolid(Planet, Params, iStart, iEnd, EOS):

    thisMAbove_kg = np.sum(Planet.MLayer_kg[:iStart-1])
    for i in range(iStart, iEnd):
        # Increment depth based on change in pressure, combined with gravity and density
        Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1]) * 1e6 / Planet.g_ms2[i-1] / \
                        Planet.rho_kgm3[i-1]
        Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
        Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1] ** 3 - Planet.r_m[i] ** 3)
        thisMAbove_kg += Planet.MLayer_kg[i-1]
        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i] ** 2

        # Propagate adiabatic thermal profile
        Planet.T_K[i] = Planet.T_K[i-1] + Planet.T_K[i-1] * Planet.alpha_pK[i-1] / \
                        Planet.Cp_JkgK[i-1] / Planet.rho_kgm3[i-1] * (Planet.P_MPa[i] - Planet.P_MPa[i-1]) * 1e6
        Planet.rhoMatrix_kgm3[i] = EOS.fn_rho_kgm3(  Planet.P_MPa[i], Planet.T_K[i], grid=False)
        Planet.Cp_JkgK[i] =        EOS.fn_Cp_JkgK(   Planet.P_MPa[i], Planet.T_K[i], grid=False)
        Planet.alpha_pK[i] =       EOS.fn_alpha_pK(  Planet.P_MPa[i], Planet.T_K[i], grid=False)
        Planet.kTherm_WmK[i] =     EOS.fn_kTherm_WmK(Planet.P_MPa[i], Planet.T_K[i], grid=False)
        Planet.rho_kgm3[i] = Planet.rhoMatrix_kgm3[i] + 0.0

        log.debug(f'il: {i:d}; P_MPa: {Planet.P_MPa[i]:.3f}; T_K: {Planet.T_K[i]:.3f}; phase: {Planet.phase[i]:d}')

    return Planet


def PropagateAdiabaticPorousVacIce(Planet, Params, iStart, iEnd, EOS):

    if iStart == 0:
        raise RuntimeError('Adiabatic calculations rely on overlying layer properties to begin recursion.')
    thisMAbove_kg = np.sum(Planet.MLayer_kg[:iStart-1])

    for i in range(iStart, iEnd):
            # Increment depth based on change in pressure, combined with gravity and density
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2

            Planet.T_K[i] = Planet.T_K[i-1] + Planet.T_K[i-1] * Planet.alpha_pK[i-1] / \
                            Planet.Cp_JkgK[i-1] / Planet.rho_kgm3[i-1] * (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6
            Planet.rhoMatrix_kgm3[i] = EOS.fn_rho_kgm3(  Planet.P_MPa[i], Planet.T_K[i], grid=False)
            Planet.Cp_JkgK[i] =        EOS.fn_Cp_JkgK(   Planet.P_MPa[i], Planet.T_K[i], grid=False)
            Planet.alpha_pK[i] =       EOS.fn_alpha_pK(  Planet.P_MPa[i], Planet.T_K[i], grid=False)
            Planet.kTherm_WmK[i] =     EOS.fn_kTherm_WmK(Planet.P_MPa[i], Planet.T_K[i], grid=False)

            # Correct for porosity in ice I layers
            Planet.phi_frac[i] = EOS.fn_phi_frac(Planet.P_MPa[i], Planet.T_K[i])
            if Planet.phi_frac[i] >= Planet.Ocean.phiMin_frac:
                Planet.rho_kgm3[i] = EOS.fn_porosCorrect(Planet.rhoMatrix_kgm3[i], 0,
                                                           Planet.phi_frac[i], Planet.Ocean.Jrho)
                Planet.Cp_JkgK[i] = EOS.fn_porosCorrect(Planet.Cp_JkgK[i] * Planet.rhoMatrix_kgm3[i], 0,
                                                           Planet.phi_frac[i], Planet.Ocean.JCp) / Planet.rho_kgm3[i]
                Planet.alpha_pK[i] = EOS.fn_porosCorrect(Planet.alpha_pK[i], 0,
                                                           Planet.phi_frac[i], Planet.Ocean.Jalpha)
                Planet.kTherm_WmK[i] = EOS.fn_porosCorrect(Planet.kTherm_WmK[i], 0,
                                                           Planet.phi_frac[i], Planet.Ocean.JkTherm)
            else:
                Planet.phi_frac[i] = 0.0
                Planet.rho_kgm3[i] = Planet.rhoMatrix_kgm3[i] + 0.0

            log.debug(f'il: {i:d}; P_MPa: {Planet.P_MPa[i]:.3f}; T_K: {Planet.T_K[i]:.3f}; phase: {Planet.phase[i]:d}')

    return Planet


def PropagateAdiabaticPorousFilledIce(Planet, Params, iStart, iEnd, EOS, EOSpore):
    # For use when the phase of all matrix layers is known, e.g. in convection
    # calculations when we're modeling the middle convecting layers

    alphaPeff = Planet.Ocean.alphaPeff[EOS.phaseStr]

    phasePore = np.zeros(iEnd - iStart, dtype=np.int_)
    DeltaPpore_MPa = 0.0

    if iStart == 0:
        raise RuntimeError('Adiabatic calculations rely on overlying layer properties to begin recursion.')
    thisMAbove_kg = np.sum(Planet.MLayer_kg[:iStart-1])
    if Planet.Ppore_MPa[iStart-1] == 0:
        Planet.Ppore_MPa[iStart-1] = Planet.P_MPa[i-1] + 0.0

    for i in range(iStart, iEnd):
        # Increment depth based on change in pressure, combined with gravity and density
        Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
        Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
        Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
        thisMAbove_kg += Planet.MLayer_kg[i-1]
        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2

        Planet.T_K[i] = Planet.T_K[i-1] + Planet.T_K[i-1] * Planet.alpha_pK[i-1] / \
                        Planet.Cp_JkgK[i-1] / Planet.rho_kgm3[i-1] * (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6
        Planet.rhoMatrix_kgm3[i] = EOS.fn_rho_kgm3(  Planet.P_MPa[i], Planet.T_K[i], grid=False)
        Planet.Cp_JkgK[i] =        EOS.fn_Cp_JkgK(   Planet.P_MPa[i], Planet.T_K[i], grid=False)
        Planet.alpha_pK[i] =       EOS.fn_alpha_pK(  Planet.P_MPa[i], Planet.T_K[i], grid=False)
        Planet.kTherm_WmK[i] =     EOS.fn_kTherm_WmK(Planet.P_MPa[i], Planet.T_K[i], grid=False)

        # Correct for porosity in ice I layers
        Planet.Ppore_MPa[i] = Planet.Ppore_MPa[i-1] + DeltaPpore_MPa
        Peff_MPa = Planet.P_MPa[i] - alphaPeff * Planet.Ppore_MPa[i]
        Planet.phi_frac[i] = EOS.fn_phi_frac(Peff_MPa, Planet.T_K[i])
        if Planet.phi_frac[i] >= Planet.Ocean.phiMin_frac:
            # In ices, assume pores are filled only with liquid.

            rhoPore_kgm3 =   EOSpore.fn_rho_kgm3(  Planet.Ppore_MPa[i], Planet.T_K[i], grid=False)
            CpPore_JkgK =    EOSpore.fn_Cp_JkgK(   Planet.Ppore_MPa[i], Planet.T_K[i], grid=False)
            alphaPore_pK =   EOSpore.fn_alpha_pK(  Planet.Ppore_MPa[i], Planet.T_K[i], grid=False)
            kThermPore_WmK = EOSpore.fn_kTherm_WmK(Planet.Ppore_MPa[i], Planet.T_K[i], grid=False)

            Planet.rho_kgm3[i] = EOS.fn_porosCorrect(Planet.rhoMatrix_kgm3[i], rhoPore_kgm3,
                                                       Planet.phi_frac[i], Planet.Ocean.Jrho)
            Planet.Cp_JkgK[i] = EOS.fn_porosCorrect(Planet.Cp_JkgK[i] * Planet.rhoMatrix_kgm3[i], CpPore_JkgK * rhoPore_kgm3,
                                                       Planet.phi_frac[i], Planet.Ocean.JCp) / Planet.rho_kgm3[i]
            Planet.alpha_pK[i] = EOS.fn_porosCorrect(Planet.alpha_pK[i], alphaPore_pK,
                                                       Planet.phi_frac[i], Planet.Ocean.Jalpha)
            Planet.kTherm_WmK[i] = EOS.fn_porosCorrect(Planet.kTherm_WmK[i], kThermPore_WmK,
                                                       Planet.phi_frac[i], Planet.Ocean.JkTherm)
            # Get approximate pore pressure change, using the previous depth change. This value
            # should be close to the correct value -- we do it this way to avoid array indexing
            # issues at the start/end points and the need to calculate the depth change twice
            DeltaPpore_MPa = 1e-6 * rhoPore_kgm3 * Planet.g_ms2[i] * (Planet.z_m[i] - Planet.z_m[i-1])
        else:
            Planet.phi_frac[i] = 0.0
            Planet.rho_kgm3[i] = Planet.rhoMatrix_kgm3[i] + 0.0
            Planet.Ppore_MPa[i] = Planet.P_MPa[i] + 0.0
            DeltaPpore_MPa = Planet.P_MPa[i] - Planet.P_MPa[i-1]

        log.debug(f'il: {i:d}; P_MPa: {Planet.P_MPa[i]:.3f}; T_K: {Planet.T_K[i]:.3f}; phase: {Planet.phase[i]:d}')

    return Planet, phasePore


def PropagateConductionProfilesSolid(Planet, Params, nProfiles, profRange, rSilEnd_m):

    # Initialize output arrays and working arrays
    Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg, \
    MHydro_kg, gSil_ms2, phiSil_frac, HtidalSil_Wm3, KSsil_GPa, GSsil_GPa, \
    Ppore_MPa, rhoPore_kgm3, phasePore, qTop_Wm2, fn_g_ms2 \
        = InitSil(Planet, Params, nProfiles, profRange, rSilEnd_m)

    # Calculate initial values based on matrix properties since we're not modeling porosity
    MLayerSil_kg[:,0] = rhoSil_kgm3[:,0] * 4/3*np.pi*(rSil_m[:,0]**3 - rSil_m[:,1]**3)
    HtidalSil_Wm3[:,0] = Planet.Sil.fn_Htidal_Wm3(rhoSil_kgm3[:,0], gSil_ms2[:,0], KSsil_GPa, GSsil_GPa)

    # Layer recursion across profiles
    Psil_MPa, Tsil_K, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, \
    HtidalSil_Wm3, kThermSil_WmK = SilRecursionSolid(Planet, Params,
        Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg,
        gSil_ms2, HtidalSil_Wm3, KSsil_GPa, GSsil_GPa, qTop_Wm2, fn_g_ms2)

    return Planet, Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, \
           gSil_ms2, phiSil_frac, HtidalSil_Wm3, kThermSil_WmK, Ppore_MPa, rhoSil_kgm3, \
           rhoPore_kgm3, phasePore


def PropagateConductionProfilesPorous(Planet, Params, nProfiles, profRange, rSilEnd_m):

    # Initialize output arrays and working arrays in common with non-porous modeling
    Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg, \
    MHydro_kg, gSil_ms2, phiSil_frac, HtidalSil_Wm3, KSsil_GPa, GSsil_GPa, \
    Ppore_MPa, rhoPore_kgm3, phasePore, qTop_Wm2, fn_g_ms2 \
        = InitSil(Planet, Params, nProfiles, profRange, rSilEnd_m)

    # Initialize porosity-specific arrays
    kThermPore_WmK, KSpore_GPa, GSpore_GPa, DeltaPpore_MPa, \
    kThermTot_WmK, KStot_GPa, GStot_GPa, rhoTot_kgm3, \
    MLayerSil_kg[:,0], HtidalSil_Wm3[:,0], rhoPore_kgm3[:,0], phasePore[:,0], \
    Ppore_MPa[:,0], phiSil_frac[:,0] = InitPorous(Planet, Params, nProfiles,
        rSil_m[:,0], rSil_m[:,1], Psil_MPa[:,0], Tsil_K[:,0], gSil_ms2[:,0],
        kThermSil_WmK[:,0], rhoSil_kgm3[:,0], KSsil_GPa[:,0], GSsil_GPa[:,0])

    # Layer recursion across profiles
    Psil_MPa, Tsil_K, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, \
    HtidalSil_Wm3, kThermSil_WmK, rhoTot_kgm3, phiSil_frac, kThermTot_WmK, \
    Ppore_MPa, rhoPore_kgm3, phasePore = SilRecursionPorous(Planet, Params,
        Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg,
        gSil_ms2, HtidalSil_Wm3, qTop_Wm2, fn_g_ms2,
        KSsil_GPa, GSsil_GPa, phiSil_frac, Ppore_MPa, phasePore,
        rhoPore_kgm3, kThermPore_WmK, KSpore_GPa, GSpore_GPa, DeltaPpore_MPa,
        rhoTot_kgm3, kThermTot_WmK, KStot_GPa, GStot_GPa)

    return Planet, Psil_MPa, Tsil_K, rSil_m, rhoTot_kgm3, MLayerSil_kg, MAboveSil_kg, \
           gSil_ms2, phiSil_frac, HtidalSil_Wm3, kThermTot_WmK, Ppore_MPa, rhoSil_kgm3, \
           rhoPore_kgm3, phasePore


def InitSil(Planet, Params, nProfiles, profRange, rSilEnd_m):

    Psil_MPa, Tsil_K, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg, gSil_ms2, \
    phiSil_frac, HtidalSil_Wm3, Ppore_MPa, rhoPore_kgm3, KSsil_GPa, GSsil_GPa \
        = (np.zeros((nProfiles, Planet.Steps.nSilMax)) for _ in range(13))
    phasePore = np.zeros((nProfiles, Planet.Steps.nSilMax), dtype=np.int_)

    # Assign ocean-bottom options as silicate-top values
    rSil_m = np.array([np.linspace(Planet.r_m[i+Planet.Steps.iSilStart], rSilEnd_m, Planet.Steps.nSilMax+1) for i in profRange])
    Psil_MPa[:,0] = [Planet.P_MPa[i+Planet.Steps.iSilStart] for i in profRange]
    Tsil_K[:,0] = [Planet.T_K[i+Planet.Steps.iSilStart] for i in profRange]
    rhoSil_kgm3[:,0] = Planet.Sil.EOS.fn_rho_kgm3(Psil_MPa[:,0], Tsil_K[:,0], grid=False)
    kThermSil_WmK[:,0] = Planet.Sil.EOS.fn_kTherm_WmK(Psil_MPa[:,0], Tsil_K[:,0], grid=False)
    KSsil_GPa[:,0] = Planet.Sil.EOS.fn_KS_GPa(Psil_MPa[:,0], Tsil_K[:,0], grid=False)
    GSsil_GPa[:,0] = Planet.Sil.EOS.fn_GS_GPa(Psil_MPa[:,0], Tsil_K[:,0], grid=False)
    gSil_ms2[:,0] = [Planet.g_ms2[i+Planet.Steps.iSilStart] for i in profRange]

    MHydro_kg = np.array([np.sum(Planet.MLayer_kg[:i]) for i in range(Planet.Steps.iSilStart, Planet.Steps.iSilStart + nProfiles)])
    # Initialize MAbove_kg to 0th silicate layer, so that the hydrosphere mass is equal to the mass above the silicates.
    MAboveSil_kg[:,0] = MHydro_kg + 0.0
    # Initialize qTop_WmK, the heat flux leaving the top of each layer.
    # For the top layer, this is equal to the heat flux warming the hydrosphere.
    qTop_Wm2 = Planet.Ocean.QfromMantle_W / 4/np.pi/rSil_m[:,0]**2

    # Adjust gravity model behavior in case these layers extend to the origin.
    # This implementation may need to be changed to account for varying silicate
    # physical properties with depth, esp. porosity, that will affect the validity
    # of the constant-density approximation giving linear gravity.
    if rSilEnd_m < 1:
        # Constant density approximation yields usable gravity values that
        # introduce less inaccuracy than gravity values exploding because of
        # bulk mass misfit from MoI searching
        fn_g_ms2 = lambda MAboveLayer_kg, rLayerBot_m: \
            gSil_ms2[:,0] * rLayerBot_m / rSil_m[:,0]
    else:
        fn_g_ms2 = lambda MAboveLayer_kg, rLayerBot_m: \
            Constants.G * np.abs(Planet.Bulk.M_kg - MAboveLayer_kg) / rLayerBot_m**2

    return Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg, \
           MHydro_kg, gSil_ms2, phiSil_frac, HtidalSil_Wm3, KSsil_GPa, GSsil_GPa, \
           Ppore_MPa, rhoPore_kgm3, phasePore, qTop_Wm2, fn_g_ms2


def InitPorous(Planet, Params, nProfiles, rSil_m0, rSil_m1, Psil_MPa0, Tsil_K0, gSil_ms20, kThermSil_WmK0,
    rhoSil_kgm30, KSsil_GPa0, GSsil_GPa0):

    kThermPore_WmK, KSpore_GPa, GSpore_GPa, DeltaPpore_MPa, \
        kThermTot_WmK, KStot_GPa, GStot_GPa, rhoTot_kgm3 \
        = (np.zeros((nProfiles, Planet.Steps.nSilMax)) for _ in range(8))
    Ppore_MPa0 = Psil_MPa0 + 0.0
    Peff_MPa = Psil_MPa0 - Planet.Sil.alphaPeff * Ppore_MPa0
    phiSil_frac0 = Planet.Sil.fn_phi_frac(Peff_MPa, Tsil_K0)
    # Get pore fluid properties
    # First check for HP ice phases
    phasePore0 = Planet.Ocean.EOS.fn_phase(Ppore_MPa0, Tsil_K0)
    liqP = phasePore0 == 0
    # If all the pores are filled with liquid, do everything together
    if np.all(liqP):
        rhoPore_kgm30 = Planet.Ocean.EOS.fn_rho_kgm3(Ppore_MPa0, Tsil_K0, grid=False)
        kThermPore_WmK[:,0] = Planet.Ocean.EOS.fn_kTherm_WmK(Ppore_MPa0, Tsil_K0, grid=False)
        _, KSpore_GPa[:,0] = Planet.Ocean.EOS.fn_Seismic(Ppore_MPa0, Tsil_K0)
        GSpore_GPa[:,0] = np.zeros_like(KSpore_GPa[:,0])
        DeltaPpore_MPa[:,0] = 1e-6 * rhoPore_kgm30 * gSil_ms20 * (rSil_m0 - rSil_m1)
    else:
        rhoPore_kgm30 = np.zeros_like(Psil_MPa0)
        # For the pores filled with liquid, use ocean liquid EOS
        rhoPore_kgm30[liqP] = Planet.Ocean.EOS.fn_rho_kgm3(Ppore_MPa0[liqP], Tsil_K0[liqP], grid=False)
        kThermPore_WmK[liqP,0] = Planet.Ocean.EOS.fn_kTherm_WmK(Ppore_MPa0[liqP], Tsil_K0[liqP], grid=False)
        _, KSpore_GPa[liqP,0] = Planet.Ocean.EOS.fn_Seismic(Ppore_MPa0[liqP], Tsil_K0[liqP])
        GSpore_GPa[liqP,0] = np.zeros_like(KSpore_GPa[liqP,0])
        DeltaPpore_MPa[liqP,0] = 1e-6 * rhoPore_kgm30[liqP] * gSil_ms20[liqP] * (rSil_m0[liqP] - rSil_m1[liqP])
        phases = np.unique(phasePore0)
        # For each HP ice phase represented, perform the same calculations
        for phase in phases[phases != 0]:
            if phase == 1:
                # We only get here if the MoI is consistent with a fully frozen ocean
                thisIceEOS = Planet.Ocean.surfIceEOS['Ih']
            else:
                thisIceEOS = Planet.Ocean.iceEOS[PhaseConv(phase)]
            iP = np.where(phasePore0 == phase)[0]
            rhoPore_kgm30[iP] = thisIceEOS.fn_rho_kgm3(Ppore_MPa0[iP], Tsil_K0[iP], grid=False)
            kThermPore_WmK[iP,0] = thisIceEOS.fn_kTherm_WmK(Ppore_MPa0[iP], Tsil_K0[iP], grid=False)
            GSpore_GPa[iP,0], KSpore_GPa[iP,0], _, _ = thisIceEOS.fn_Seismic(Ppore_MPa0[iP], Tsil_K0[iP])
            # Keep DeltaPpore calculation separate for liquid and ice in case of future updates to adjust matrix coupling/pressure scaling
            DeltaPpore_MPa[iP,0] = 1e-6 * rhoPore_kgm30[iP] * gSil_ms20[iP] * (rSil_m0[iP] - rSil_m1[iP])
    # Combine properties using rules defined in EOS functions for porosity
    rhoTot_kgm3[:,0] = Planet.Sil.EOS.fn_porosCorrect(rhoSil_kgm30, rhoPore_kgm30, phiSil_frac0, Planet.Sil.Jrho)
    kThermTot_WmK[:,0] = Planet.Sil.EOS.fn_porosCorrect(kThermSil_WmK0, kThermPore_WmK[:,0], phiSil_frac0, Planet.Sil.JkTherm)
    KStot_GPa[:,0] = Planet.Sil.EOS.fn_porosCorrect(KSsil_GPa0, KSpore_GPa[:,0], phiSil_frac0, Planet.Sil.JKS)
    GStot_GPa[:,0] = Planet.Sil.EOS.fn_porosCorrect(GSsil_GPa0, GSpore_GPa[:,0], phiSil_frac0, Planet.Sil.JGS)

    MLayerSil_kg0 = rhoTot_kgm3[:,0] * 4/3*np.pi*(rSil_m0**3 - rSil_m1**3)
    HtidalSil_Wm30 = Planet.Sil.fn_Htidal_Wm3(rhoTot_kgm3[:,0], gSil_ms20, KStot_GPa[:,0], GStot_GPa[:,0])

    return kThermPore_WmK, KSpore_GPa, GSpore_GPa, DeltaPpore_MPa, \
           kThermTot_WmK, KStot_GPa, GStot_GPa, rhoTot_kgm3, \
           MLayerSil_kg0, HtidalSil_Wm30, rhoPore_kgm30, phasePore0, \
           Ppore_MPa0, phiSil_frac0


def SilRecursionSolid(Planet, Params,
    Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg,
    gSil_ms2, HtidalSil_Wm3, KSsil_GPa, GSsil_GPa, qTop_Wm2, fn_g_ms2):

    for j in range(1, Planet.Steps.nSilMax):
        MAboveSil_kg[:,j] = MAboveSil_kg[:,j-1] + MLayerSil_kg[:,j-1]
        Psil_MPa[:,j] = Psil_MPa[:,j-1] + 1e-6 * MLayerSil_kg[:,j-1] * gSil_ms2[:,j-1] / (4*np.pi*rSil_m[:,j]**2)
        Tsil_K[:,j], qTop_Wm2 = ConductiveTemperature(Tsil_K[:,j-1], rSil_m[:,j-1], rSil_m[:,j],
                    kThermSil_WmK[:,j-1], rhoSil_kgm3[:,j-1], Planet.Sil.Qrad_Wkg, HtidalSil_Wm3[:,j-1],
                    qTop_Wm2)
        rhoSil_kgm3[:,j] = Planet.Sil.EOS.fn_rho_kgm3(Psil_MPa[:,j], Tsil_K[:,j], grid=False)
        kThermSil_WmK[:,j] = Planet.Sil.EOS.fn_kTherm_WmK(Psil_MPa[:,j], Tsil_K[:,j], grid=False)
        # Get KS and GS now as they are needed for Htidal calculation;
        # we will calculate them again later along with other seismic calcs
        KSsil_GPa[:,j] = Planet.Sil.EOS.fn_KS_GPa(Psil_MPa[:,j], Tsil_K[:,j], grid=False)
        GSsil_GPa[:,j] = Planet.Sil.EOS.fn_GS_GPa(Psil_MPa[:,j], Tsil_K[:,j], grid=False)
        # Calculate gravity using absolute values, as we will use MAboveSil to check for exceeding body mass later.
        gSil_ms2[:,j] = fn_g_ms2(MAboveSil_kg[:,j], rSil_m[:,j])

        MLayerSil_kg[:,j] = rhoSil_kgm3[:,j] * 4/3*np.pi*(rSil_m[:,j]**3 - rSil_m[:,j+1]**3)
        HtidalSil_Wm3[:,j] = Planet.Sil.fn_Htidal_Wm3(rhoSil_kgm3[:,j], gSil_ms2[:,j], KSsil_GPa[:,j], GSsil_GPa[:,j])


    return Psil_MPa, Tsil_K, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, \
           HtidalSil_Wm3, kThermSil_WmK


def SilRecursionPorous(Planet, Params,
    Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg,
    gSil_ms2, HtidalSil_Wm3, qTop_Wm2, fn_g_ms2,
    KSsil_GPa, GSsil_GPa, phiSil_frac, Ppore_MPa, phasePore,
    rhoPore_kgm3, kThermPore_WmK, KSpore_GPa, GSpore_GPa, DeltaPpore_MPa,
    rhoTot_kgm3, kThermTot_WmK, KStot_GPa, GStot_GPa):

    for j in range(1, Planet.Steps.nSilMax):
        MAboveSil_kg[:,j] = MAboveSil_kg[:,j-1] + MLayerSil_kg[:,j-1]
        Psil_MPa[:,j] = Psil_MPa[:,j-1] + 1e-6 * MLayerSil_kg[:,j-1] * gSil_ms2[:,j-1] / (4*np.pi*rSil_m[:,j]**2)
        Tsil_K[:,j], qTop_Wm2 = ConductiveTemperature(Tsil_K[:,j-1], rSil_m[:,j-1], rSil_m[:,j],
                    kThermSil_WmK[:,j-1], rhoSil_kgm3[:,j-1], Planet.Sil.Qrad_Wkg, HtidalSil_Wm3[:,j-1],
                    qTop_Wm2)
        rhoSil_kgm3[:,j] = Planet.Sil.EOS.fn_rho_kgm3(Psil_MPa[:,j], Tsil_K[:,j], grid=False)
        kThermSil_WmK[:,j] = Planet.Sil.EOS.fn_kTherm_WmK(Psil_MPa[:,j], Tsil_K[:,j], grid=False)
        # Get KS and GS now as they are needed for Htidal calculation;
        # we will calculate them again later along with other seismic calcs
        KSsil_GPa[:,j] = Planet.Sil.EOS.fn_KS_GPa(Psil_MPa[:,j], Tsil_K[:,j], grid=False)
        GSsil_GPa[:,j] = Planet.Sil.EOS.fn_GS_GPa(Psil_MPa[:,j], Tsil_K[:,j], grid=False)
        # Calculate gravity using absolute values, as we will use MAboveSil to check for exceeding body mass later.
        gSil_ms2[:,j] = fn_g_ms2(MAboveSil_kg[:,j], rSil_m[:,j])

        # Adjust for properties of pore material
        Ppore_MPa[:,j] = Ppore_MPa[:,j-1] + DeltaPpore_MPa[:,j-1]
        Peff_MPa = Psil_MPa[:,j] - Planet.Sil.alphaPeff * Ppore_MPa[:,j]
        phiSil_frac[:,j] = Planet.Sil.fn_phi_frac(Peff_MPa, Tsil_K[:,j])
        # Get pore fluid properties
        # First check for HP ice phases
        phasePore[:,j] = Planet.Ocean.EOS.fn_phase(Ppore_MPa[:,j], Tsil_K[:,j])
        liqP = phasePore[:,j] == 0
        # If all the pores are filled with liquid, do everything together
        if np.all(liqP):
            rhoPore_kgm3[:,j] = Planet.Ocean.EOS.fn_rho_kgm3(Ppore_MPa[:,j], Tsil_K[:,j], grid=False)
            kThermPore_WmK[:,j] = Planet.Ocean.EOS.fn_kTherm_WmK(Ppore_MPa[:,j], Tsil_K[:,j], grid=False)
            _, KSpore_GPa[:,j] = Planet.Ocean.EOS.fn_Seismic(Ppore_MPa[:,j], Tsil_K[:,j])
            GSpore_GPa[:,j] = np.zeros_like(KSpore_GPa[:,j])
            DeltaPpore_MPa[:,j] = 1e-6 * rhoPore_kgm3[:,j] * gSil_ms2[:,j] * (rSil_m[:,j] - rSil_m[:,j+1])
        else:
            # For the pores filled with liquid, use ocean liquid EOS
            rhoPore_kgm3[liqP,j] = Planet.Ocean.EOS.fn_rho_kgm3(Ppore_MPa[liqP,j], Tsil_K[liqP,j], grid=False)
            kThermPore_WmK[liqP,j] = Planet.Ocean.EOS.fn_kTherm_WmK(Ppore_MPa[liqP,j], Tsil_K[liqP,j], grid=False)
            _, KSpore_GPa[liqP,j] = Planet.Ocean.EOS.fn_Seismic(Ppore_MPa[liqP,j], Tsil_K[liqP,j])
            GSpore_GPa[liqP,j] = np.zeros_like(KSpore_GPa[liqP,j])
            DeltaPpore_MPa[liqP,j] = 1e-6 * rhoPore_kgm3[liqP,j] * gSil_ms2[liqP,j] * (rSil_m[liqP,j] - rSil_m[liqP,j+1])
            phases = np.unique(phasePore[:,j])
            # For each HP ice phase represented, perform the same calculations
            for phase in phases[phases != 0]:
                if phase == 1:
                    # We only get here if the MoI is consistent with a fully frozen ocean
                    thisIceEOS = Planet.Ocean.surfIceEOS['Ih']
                else:
                    thisIceEOS = Planet.Ocean.iceEOS[PhaseConv(phase)]
                iP = np.where(phasePore[:,j] == phase)[0]
                rhoPore_kgm3[iP,j] = thisIceEOS.fn_rho_kgm3(Ppore_MPa[iP,j], Tsil_K[iP,j], grid=False)
                kThermPore_WmK[iP,j] = thisIceEOS.fn_kTherm_WmK(Ppore_MPa[iP,j], Tsil_K[iP,j], grid=False)
                GSpore_GPa[iP,j], KSpore_GPa[iP,j], _, _ = thisIceEOS.fn_Seismic(Ppore_MPa[iP,j], Tsil_K[iP,j])
                # Keep DeltaPpore calculation separate for liquid and ice in case of future updates to adjust matrix coupling/pressure scaling
                DeltaPpore_MPa[iP,j] = 1e-6 * rhoPore_kgm3[iP,j] * gSil_ms2[iP,j] * (rSil_m[iP,j] - rSil_m[iP,j+1])
        # Combine properties using rule defined in EOS function for porosity as in
        # Yu et al. (2016): http://dx.doi.org/10.1016/j.jrmge.2015.07.004
        rhoTot_kgm3[:,j] = Planet.Sil.EOS.fn_porosCorrect(rhoSil_kgm3[:,j], rhoPore_kgm3[:,j], phiSil_frac[:,j], Planet.Sil.Jrho)
        kThermTot_WmK[:,j] = Planet.Sil.EOS.fn_porosCorrect(kThermSil_WmK[:,j], kThermPore_WmK[:,j], phiSil_frac[:,j], Planet.Sil.JkTherm)
        KStot_GPa[:,j] = Planet.Sil.EOS.fn_porosCorrect(KSsil_GPa[:,j], KSpore_GPa[:,j], phiSil_frac[:,j], Planet.Sil.JKS)
        GStot_GPa[:,j] = Planet.Sil.EOS.fn_porosCorrect(GSsil_GPa[:,j], GSpore_GPa[:,j], phiSil_frac[:,j], Planet.Sil.JGS)

        MLayerSil_kg[:,j] = rhoTot_kgm3[:,j] * 4/3*np.pi*(rSil_m[:,j]**3 - rSil_m[:,j+1]**3)
        HtidalSil_Wm3[:,j] = Planet.Sil.fn_Htidal_Wm3(rhoTot_kgm3[:,j], gSil_ms2[:,j], KStot_GPa[:,j], GStot_GPa[:,j])

        if j % 10 == 0: log.debug(f'il_Sil: {Planet.Steps.iSilStart+j}')

    return Psil_MPa, Tsil_K, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, \
           HtidalSil_Wm3, kThermSil_WmK, rhoTot_kgm3, phiSil_frac, kThermTot_WmK, \
           Ppore_MPa, rhoPore_kgm3, phasePore


