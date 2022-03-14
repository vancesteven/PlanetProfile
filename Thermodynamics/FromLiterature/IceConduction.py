import numpy as np
import logging as log
import scipy.interpolate as spi
from Utilities.dataStructs import Constants
from Thermodynamics.FromLiterature.HydroEOS import PhaseConv, IceEOSStruct
from Thermodynamics.FromLiterature.ThermalProfiles import GetPbConduct

def IceIWholeConduct(Planet, Params):
    """ Calculate conductive profile in ice Ih layers based on a linear pressure profile
        from the surface to the melting pressure that coincides with the assumed bottom
        melting temperature, Tb_K. Use when Do.CLATHRATE = False and when
        Bulk.clathType = 'whole'.

        Assigns Planet attributes:
            All physical layer arrays
    """
    icePhase = PhaseConv(Planet.phase[0])

    # Set linear P and adiabatic T in ice I layers. Include 1 extra for P and T to assign next phase to the values
    # at the phase transition
    PIceI_MPa = np.linspace(Planet.P_MPa[0], Planet.PbI_MPa, Planet.Steps.nIbottom+1)
    Pratios = (PIceI_MPa - Planet.P_MPa[0]) / (Planet.PbI_MPa - Planet.P_MPa[0])
    TIceI_K = Planet.Bulk.Tb_K**(Pratios) * Planet.T_K[0]**(1 - PIceI_MPa/Planet.PbI_MPa)
    Planet.P_MPa[:Planet.Steps.nIbottom+1] = PIceI_MPa
    Planet.T_K[:Planet.Steps.nIbottom+1] = TIceI_K

    # Get ice EOS
    Planet.Ocean.surfIceEOS[icePhase] = IceEOSStruct(PIceI_MPa, TIceI_K, icePhase,
                                                     porosType=Planet.Ocean.porosType[icePhase],
                                                     phiTop_frac=Planet.Ocean.phiMax_frac[icePhase],
                                                     Pclosure_MPa=Planet.Ocean.Pclosure_MPa[icePhase])

    # Evaluate thermodynamic properties of uppermost ice
    Planet.rhoMatrix_kgm3[:Planet.Steps.nIbottom] \
        = Planet.Ocean.surfIceEOS[icePhase].fn_rho_kgm3(PIceI_MPa[:-1], TIceI_K[:-1], grid=False)
    Planet.Cp_JkgK[:Planet.Steps.nIbottom] \
        = Planet.Ocean.surfIceEOS[icePhase].fn_Cp_JkgK(PIceI_MPa[:-1], TIceI_K[:-1], grid=False)
    Planet.alpha_pK[:Planet.Steps.nIbottom] \
        = Planet.Ocean.surfIceEOS[icePhase].fn_alpha_pK(PIceI_MPa[:-1], TIceI_K[:-1], grid=False)
    Planet.kTherm_WmK[:Planet.Steps.nIbottom] \
        = Planet.Ocean.surfIceEOS[icePhase].fn_kTherm_WmK(PIceI_MPa[:-1], TIceI_K[:-1], grid=False)

    # Correct for porosity in ice, assuming pores contain only vacuum
    if Planet.Do.POROUS_ICE:
        Planet.phi_frac[:Planet.Steps.nIbottom] = Planet.Ocean.surfIceEOS[icePhase].fn_phi_frac(PIceI_MPa[:-1], TIceI_K[:-1], grid=False)
        Planet.rho_kgm3[:Planet.Steps.nIbottom] = Planet.Ocean.surfIceEOS[icePhase].fn_porosCorrect(
            Planet.rhoMatrix_kgm3[:Planet.Steps.nIbottom], 0, Planet.phi_frac[:Planet.Steps.nIbottom], Planet.Ocean.Jrho)
        Planet.Cp_JkgK[:Planet.Steps.nIbottom] = Planet.Ocean.surfIceEOS[icePhase].fn_porosCorrect(
            Planet.Cp_JkgK[:Planet.Steps.nIbottom] * Planet.rhoMatrix_kgm3[:Planet.Steps.nIbottom],
            0, Planet.phi_frac[:Planet.Steps.nIbottom], Planet.Ocean.JCp) / Planet.rho_kgm3[:Planet.Steps.nIbottom]
        Planet.alpha_pK[:Planet.Steps.nIbottom] = Planet.Ocean.surfIceEOS[icePhase].fn_porosCorrect(
            Planet.alpha_pK[:Planet.Steps.nIbottom], 0, Planet.phi_frac[:Planet.Steps.nIbottom], Planet.Ocean.Jalpha)
        Planet.kTherm_WmK[:Planet.Steps.nIbottom] = Planet.Ocean.surfIceEOS[icePhase].fn_porosCorrect(
            Planet.kTherm_WmK[:Planet.Steps.nIbottom], 0, Planet.phi_frac[:Planet.Steps.nIbottom], Planet.Ocean.JkTherm)
    else:
        Planet.rho_kgm3[:Planet.Steps.nIbottom] = Planet.rhoMatrix_kgm3[:Planet.Steps.nIbottom] + 0.0


    # Calculate remaining physical properties of upper ice I
    thisMAbove_kg = 0
    for i in range(1, Planet.Steps.nIbottom+1):
        # Increment depth based on change in pressure, combined with gravity and density
        Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
        Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
        Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
        thisMAbove_kg += Planet.MLayer_kg[i-1]
        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
        log.debug(f'il: {i:d}; P_MPa: {Planet.P_MPa[i]:.3f}; T_K: {Planet.T_K[i]:.3f}; phase: {Planet.phase[i]:d}')

    return Planet


def IceIConductClathLid(Planet, Params):
    """ Calculate conductive profile in ice Ih layers beneath a clathrate lid based
        on a linear pressure profile from the surface to the melting pressure that
        coincides with the assumed bottom melting temperature, Tb_K.

        Assigns Planet attributes:
            All physical layer arrays
    """
    # Assign phases for clathrates, as number of layers is fixed in this case
    Planet.phase[:Planet.Steps.nClath] = Constants.phaseClath
    # Set linear P and T in ice I layers to use in surfIce.EOS functions
    Plin_MPa = np.linspace(Planet.Bulk.Psurf_MPa, Planet.PbI_MPa, Planet.Steps.nIbottom+1)
    Tlin_K = np.linspace(Planet.Bulk.Tsurf_K, Planet.Bulk.Tb_K, Planet.Steps.nIbottom+1)

    # Get ice Ih EOS
    Planet.Ocean.surfIceEOS['Ih'] = IceEOSStruct(Plin_MPa, Tlin_K, 'Ih',
                                                 porosType=Planet.Ocean.porosType['Ih'],
                                                 phiTop_frac=Planet.Ocean.phiMax_frac['Ih'],
                                                 Pclosure_MPa=Planet.Ocean.Pclosure_MPa['Ih'])
    # Get clathrate EOS
    Planet.Ocean.surfIceEOS['Clath'] = IceEOSStruct(Plin_MPa, Tlin_K, 'Clath',
                                                    porosType=Planet.Ocean.porosType['Clath'],
                                                    phiTop_frac=Planet.Ocean.phiMax_frac['Clath'],
                                                    Pclosure_MPa=Planet.Ocean.Pclosure_MPa['Clath'])

    # Get approximate temperature of transition between clathrates and ice I
    zbApprox_m = (Planet.PbI_MPa - Planet.Bulk.Psurf_MPa) * 1e6 / Planet.g_ms2[0] / \
                 Planet.Ocean.surfIceEOS['Ih'].fn_rho_kgm3(Planet.PbI_MPa, Planet.Bulk.Tb_K, grid=False)
    if Planet.Bulk.clathMaxThick_m > 0.9*zbApprox_m:
        raise ValueError('Bulk.clathMaxThick_m is greater than 90% of the approximate ice shell thickness, ' +
                         'but for this model type it should be substantially less. Reduce Bulk.clathMaxThick_m ' +
                         'or decrease Tb_K to address this issue.')
    Pmid_MPa = (Planet.PbI_MPa - Planet.Bulk.Psurf_MPa) / 2
    Tmid_K = (Planet.Bulk.Tb_K - Planet.Bulk.Tsurf_K) / 2
    kImid_WmK = Planet.Ocean.surfIceEOS['Ih'].fn_kTherm_WmK(Pmid_MPa, Tmid_K, grid=False)
    kClathTop_WmK = Planet.Ocean.surfIceEOS['Clath'].fn_kTherm_WmK(Pmid_MPa, Tmid_K, grid=False)
    # Adjust thermal conductivities in the case of modeling porosity
    if Planet.Do.POROUS_ICE:
        kImid_WmK = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(kImid_WmK, 0,
                    Planet.Ocean.surfIceEOS['Ih'].fn_phi_frac(Pmid_MPa, Tmid_K, grid=False), Planet.Ocean.JkTherm)
        kClathTop_WmK = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(kClathTop_WmK, 0,
                        Planet.Ocean.surfIceEOS['Clath'].fn_phi_frac(Pmid_MPa, Tmid_K, grid=False), Planet.Ocean.JkTherm)
    thickRatio = (zbApprox_m - Planet.Bulk.clathMaxThick_m)/Planet.Bulk.clathMaxThick_m
    kRatio = kImid_WmK / kClathTop_WmK
    Planet.TclathTrans_K = (Planet.Bulk.Tsurf_K * thickRatio + Planet.Bulk.Tb_K * kRatio) / (thickRatio + kRatio)
    qTop_Wm2 = (Planet.TclathTrans_K - Planet.Bulk.Tsurf_K) / Planet.Bulk.clathMaxThick_m * kClathTop_WmK
    Planet.Ocean.QfromMantle_W = qTop_Wm2 * 4*np.pi*Planet.Bulk.R_m**2
    # Get P at the transition depth
    Planet.PbClathMax_MPa = GetPbConduct(Planet.Bulk.Tsurf_K, Planet.TclathTrans_K, Planet.Bulk.R_m, Planet.Bulk.Psurf_MPa,
                               Planet.g_ms2[0], qTop_Wm2, Planet.Ocean.surfIceEOS['Clath'],
                               rRes_m=Planet.Bulk.clathMaxThick_m/Planet.Steps.nClath/2)
    log.debug('Found approx. bottom pressure of clathrate layer at ' +
             f'{Planet.PbClathMax_MPa:.3f} MPa for Bulk.clathMaxThick = ' +
             f'{Planet.Bulk.clathMaxThick_m/1e3:.3f} km and DeltaTclath = ' +
             f'{Planet.TclathTrans_K - Planet.Bulk.Tsurf_K:.3f} K.')
    if Planet.PbClathMax_MPa > Planet.PbI_MPa:
        raise ValueError(f'Clathrate transition pressure of {Planet.PbClathMax_MPa:.3f} MPa is less ' +
                         f'than ice I melting pressure of {Planet.PbI_MPa:.3f} MPa. Bulk.clathMaxThick ' +
                         f'= {Planet.Bulk.clathMaxThick_m/1e3:.3f} km is likely too large for Tb = ' +
                         f'{Planet.Bulk.Tb_K:.3f} K.')

    # Reset thermal and pressure profiles to be linear in P and adiabatic in T
    # Include 1 extra for P and T to assign next phase to the values at the
    # ice I bottom phase transition
    Pclath_MPa = np.linspace(Planet.Bulk.Psurf_MPa, Planet.PbClathMax_MPa, Planet.Steps.nClath+1)[:-1]
    PIceI_MPa = np.linspace(Planet.PbClathMax_MPa, Planet.PbI_MPa, Planet.Steps.nIceI+1)
    PratiosClath = (Pclath_MPa - Planet.P_MPa[0]) / (Planet.PbClathMax_MPa - Planet.P_MPa[0])
    Tclath_K = Planet.TclathTrans_K**(PratiosClath) * Planet.Bulk.Tsurf_K**(1 - PratiosClath)
    PratiosIceI = (PIceI_MPa - Planet.PbClathMax_MPa) / (Planet.PbI_MPa - Planet.PbClathMax_MPa)
    TIceI_K = Planet.Bulk.Tb_K**(PratiosIceI) * Planet.TclathTrans_K**(1 - PratiosIceI)
    Planet.P_MPa[:Planet.Steps.nClath] = Pclath_MPa
    Planet.P_MPa[Planet.Steps.nClath:Planet.Steps.nIbottom+1] = PIceI_MPa
    Planet.T_K[:Planet.Steps.nClath] = Tclath_K
    Planet.T_K[Planet.Steps.nClath:Planet.Steps.nIbottom+1] = TIceI_K

    # Evaluate thermodynamic properties of clathrate lid
    Planet.rhoMatrix_kgm3[:Planet.Steps.nClath] \
        = Planet.Ocean.surfIceEOS['Clath'].fn_rho_kgm3(Pclath_MPa, Tclath_K, grid=False)
    Planet.Cp_JkgK[:Planet.Steps.nClath] \
        = Planet.Ocean.surfIceEOS['Clath'].fn_Cp_JkgK(Pclath_MPa, Tclath_K, grid=False)
    Planet.alpha_pK[:Planet.Steps.nClath] \
        = Planet.Ocean.surfIceEOS['Clath'].fn_alpha_pK(Pclath_MPa, Tclath_K, grid=False)
    Planet.kTherm_WmK[:Planet.Steps.nClath] \
        = Planet.Ocean.surfIceEOS['Clath'].fn_kTherm_WmK(Pclath_MPa, Tclath_K, grid=False)
    # Evaluate thermodynamic properties of uppermost ice I
    Planet.rhoMatrix_kgm3[Planet.Steps.nClath:Planet.Steps.nIbottom] \
        = Planet.Ocean.surfIceEOS['Ih'].fn_rho_kgm3(PIceI_MPa[:-1], TIceI_K[:-1], grid=False)
    Planet.Cp_JkgK[Planet.Steps.nClath:Planet.Steps.nIbottom] \
        = Planet.Ocean.surfIceEOS['Ih'].fn_Cp_JkgK(PIceI_MPa[:-1], TIceI_K[:-1], grid=False)
    Planet.alpha_pK[Planet.Steps.nClath:Planet.Steps.nIbottom] \
        = Planet.Ocean.surfIceEOS['Ih'].fn_alpha_pK(PIceI_MPa[:-1], TIceI_K[:-1], grid=False)
    Planet.kTherm_WmK[Planet.Steps.nClath:Planet.Steps.nIbottom] \
        = Planet.Ocean.surfIceEOS['Ih'].fn_kTherm_WmK(PIceI_MPa[:-1], TIceI_K[:-1], grid=False)

    # Correct for porosity in ice I and clathrates, assuming pores contain only vacuum
    if Planet.Do.POROUS_ICE:
        Planet.phi_frac[:Planet.Steps.nClath] = Planet.Ocean.surfIceEOS['Clath'].fn_phi_frac(Pclath_MPa, Tclath_K, grid=False)
        Planet.rho_kgm3[:Planet.Steps.nClath] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
            Planet.rhoMatrix_kgm3[:Planet.Steps.nClath], 0, Planet.phi_frac[:Planet.Steps.nClath], Planet.Ocean.Jrho)
        Planet.Cp_JkgK[:Planet.Steps.nClath] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
            Planet.Cp_JkgK[:Planet.Steps.nClath] * Planet.rhoMatrix_kgm3[:Planet.Steps.nClath],
            0, Planet.phi_frac[:Planet.Steps.nClath], Planet.Ocean.JCp) / Planet.rho_kgm3[:Planet.Steps.nClath]
        Planet.alpha_pK[:Planet.Steps.nClath] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
            Planet.alpha_pK[:Planet.Steps.nClath], 0, Planet.phi_frac[:Planet.Steps.nClath], Planet.Ocean.Jalpha)
        Planet.kTherm_WmK[:Planet.Steps.nClath] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
            Planet.kTherm_WmK[:Planet.Steps.nClath], 0, Planet.phi_frac[:Planet.Steps.nClath], Planet.Ocean.JkTherm)

        Planet.phi_frac[Planet.Steps.nClath:Planet.Steps.nIbottom] = Planet.Ocean.surfIceEOS['Ih'].fn_phi_frac(PIceI_MPa[:-1], TIceI_K[:-1], grid=False)
        Planet.rho_kgm3[Planet.Steps.nClath:Planet.Steps.nIbottom] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.rhoMatrix_kgm3[Planet.Steps.nClath:Planet.Steps.nIbottom], 0, Planet.phi_frac[Planet.Steps.nClath:Planet.Steps.nIbottom], Planet.Ocean.Jrho)
        Planet.Cp_JkgK[Planet.Steps.nClath:Planet.Steps.nIbottom] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Cp_JkgK[Planet.Steps.nClath:Planet.Steps.nIbottom] * Planet.rhoMatrix_kgm3[Planet.Steps.nClath:Planet.Steps.nIbottom],
            0, Planet.phi_frac[Planet.Steps.nClath:Planet.Steps.nIbottom], Planet.Ocean.JCp) / Planet.rho_kgm3[Planet.Steps.nClath:Planet.Steps.nIbottom]
        Planet.alpha_pK[Planet.Steps.nClath:Planet.Steps.nIbottom] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.alpha_pK[Planet.Steps.nClath:Planet.Steps.nIbottom], 0, Planet.phi_frac[Planet.Steps.nClath:Planet.Steps.nIbottom], Planet.Ocean.Jalpha)
        Planet.kTherm_WmK[Planet.Steps.nClath:Planet.Steps.nIbottom] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.kTherm_WmK[Planet.Steps.nClath:Planet.Steps.nIbottom], 0, Planet.phi_frac[Planet.Steps.nClath:Planet.Steps.nIbottom], Planet.Ocean.JkTherm)
    else:
        Planet.rho_kgm3[:Planet.Steps.nIbottom] = Planet.rhoMatrix_kgm3[:Planet.Steps.nIbottom] + 0.0

    # Calculate remaining physical properties of upper ice I and clathrates
    thisMAbove_kg = 0
    for i in range(1, Planet.Steps.nIbottom+1):
        # Increment depth based on change in pressure, combined with gravity and density
        Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
        Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
        Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
        thisMAbove_kg += Planet.MLayer_kg[i-1]
        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
        log.debug(f'il: {i:d}; P_MPa: {Planet.P_MPa[i]:.3f}; T_K: {Planet.T_K[i]:.3f}; phase: {Planet.phase[i]:d}')

    Planet.zClath_m = Planet.z_m[Planet.Steps.nClath]

    return Planet


def IceIConductClathUnderplate(Planet, Params):
    """ Calculate conductive profiles in ice Ih and clathrate layers based on an assumed
        thickness for a clathrate underplate layer. Use when Do.CLATHRATE = True and
        Bulk.clathType = 'bottom'.

        Assigns Planet attributes:
            All physical layer arrays
    """
    # Get clathrate EOS
    PIceFull_MPa = np.linspace(Planet.P_MPa[0], Planet.PbI_MPa, Planet.Steps.nIbottom+1)
    TIceFull_K = np.linspace(Planet.T_K[0], Planet.Bulk.Tb_K, Planet.Steps.nIbottom)
    if Planet.Do.POROUS_ICE and Planet.Ocean.phiMax_frac['Clath'] is not None:
        log.warning('Ocean.phiMaxClath_frac is set, but for clathrate underplating (Bulk.clathType = "bottom"), porosity is not modeled. Resetting Ocean.phiMaxClath_frac to None.')
        Planet.Ocean.phiMax_frac['Clath'] = None
    Planet.Ocean.surfIceEOS['Clath'] = IceEOSStruct(PIceFull_MPa, TIceFull_K, 'Clath')

    rhoBot_kgm3 = Planet.Ocean.surfIceEOS['Clath'].fn_rho_kgm3(Planet.PbI_MPa, Planet.Bulk.Tb_K, grid=False)
    # Get approximate pressure change across clathrate layer (assuming surface gravity)
    DeltaPClath_MPa = Planet.Bulk.clathMaxThick_m * Planet.g_ms2[0] * rhoBot_kgm3 / 1e6
    # Get pressure at ice Ih-clathrate transition
    PbTrans_MPa = Planet.PbI_MPa - DeltaPClath_MPa
    PIceI_MPa = np.linspace(Planet.P_MPa[0], PbTrans_MPa, Planet.Steps.nIceI+1)

    # Get ice I EOS
    Planet.Ocean.surfIceEOS['Ih'] = IceEOSStruct(PIceI_MPa, TIceFull_K, 'Ih',
                                                 porosType=Planet.Ocean.porosType['Ih'],
                                                 phiTop_frac=Planet.Ocean.phiMax_frac['Ih'],
                                                 Pclosure_MPa=Planet.Ocean.Pclosure_MPa['Ih'])

    # Get approximate temperature at top of clathrate layer based on assumed surface heat flux
    # Need approx. depth first for curvature change to heat flux
    rhoTransIceI_kgm3 = Planet.Ocean.surfIceEOS['Ih'].fn_rho_kgm3(PbTrans_MPa, Planet.Bulk.Tb_K, grid=False)
    if Planet.Do.POROUS_ICE:
        rhoTransIceI_kgm3 = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(rhoTransIceI_kgm3, 0,
                            Planet.Ocean.surfIceEOS['Ih'].fn_phi_frac(PbTrans_MPa, Planet.Bulk.Tb_K, grid=False),
                            Planet.Ocean.Jrho)
    zbApprox_m = (PbTrans_MPa - Planet.Bulk.Psurf_MPa) * 1e6 / Planet.g_ms2[0] / \
                 rhoTransIceI_kgm3
    rTopApprox_m = Planet.r_m[0] - zbApprox_m
    # Get minimum heat flux through clathrate layer (assuming it's relatively thin)
    qClathMin_Wm2 = Planet.Bulk.qSurf_Wm2 * Planet.r_m[0]**2 / rTopApprox_m**2
    kBot_WmK = Planet.Ocean.surfIceEOS['Clath'].fn_kTherm_WmK(Planet.PbI_MPa, Planet.Bulk.Tb_K, grid=False)
    Ttop_K = Planet.Bulk.Tb_K - qClathMin_Wm2 * Planet.Bulk.clathMaxThick_m / kBot_WmK
    if Ttop_K < Planet.Bulk.Tsurf_K:
        raise RuntimeError('Bulk.qSurf_Wm2 is set too high for the values set for Bulk.Tb_K, Bulk.Tsurf_K, and ' +
                           'Bulk.clathMaxThick_m. Try decreasing Bulk.qSurf_Wm2.')

    # Propagate from the surface down to the temperature at the top of the clathrate layer
    # using a spatial resolution comparable to what we expect for the end result in the shell
    rhoTop_kgm3 = Planet.Ocean.surfIceEOS['Ih'].fn_rho_kgm3(Planet.P_MPa[0], Planet.T_K[0], grid=False)
    if Planet.Do.POROUS_ICE:
        rhoTop_kgm3 = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(rhoTop_kgm3, 0,
                      Planet.Ocean.surfIceEOS['Ih'].fn_phi_frac(Planet.P_MPa[0], Planet.T_K[0], grid=False),
                      Planet.Ocean.Jrho)
    rStep_m = (Planet.PbI_MPa - Planet.P_MPa[0])*1e6 / Planet.g_ms2[0] / rhoTop_kgm3 / Planet.Steps.nIbottom

    PbTrans_MPa = GetPbConduct(Planet.T_K[0], Ttop_K, Planet.r_m[0], Planet.P_MPa[0], Planet.g_ms2[0],
                               Planet.Bulk.qSurf_Wm2, Planet.Ocean.surfIceEOS['Ih'], rRes_m=rStep_m,
                               Qrad_Wkg=0, Htidal_Wm3=0)

    # Now we have our (P,T) profile for the ice shell.
    # Now set linear P and adiabatic T for the ice I shell
    # This is not self-consistent, because we calculated PbTrans_MPa using the thermal conductivity
    # of the ice as we progressed through the ice I shell, but this is consistent with the current
    # implementation of the thermal profile of conductive ice used in other PlanetProfile models.
    PIceI_MPa = np.linspace(Planet.P_MPa[0], PbTrans_MPa, Planet.Steps.nIceI+1)[:-1]
    PIceIratios = (PIceI_MPa - Planet.P_MPa[0]) / (PbTrans_MPa - Planet.P_MPa[0])
    TIceI_K = Ttop_K**(PIceIratios) * Planet.T_K[0]**(1 - PIceIratios)
    # Next, set linear P and adiabatic T for the clathrate layer
    Pclath_MPa = np.linspace(PbTrans_MPa, Planet.PbI_MPa, Planet.Steps.nClath+1)
    PclathRatios = (Pclath_MPa - Pclath_MPa[0]) / (Planet.PbI_MPa - Pclath_MPa[0])
    Tclath_K = Planet.Bulk.Tb_K**(PclathRatios) * Ttop_K**(1 - PclathRatios)

    Planet.P_MPa[:Planet.Steps.nIceI] = PIceI_MPa
    Planet.T_K[:Planet.Steps.nIceI] = TIceI_K
    Planet.P_MPa[Planet.Steps.nIceI:Planet.Steps.nIbottom+1] = Pclath_MPa
    Planet.T_K[Planet.Steps.nIceI:Planet.Steps.nIbottom+1] = Tclath_K

    # Evaluate thermodynamic properties of uppermost ice
    Planet.rhoMatrix_kgm3[:Planet.Steps.nIceI] \
        = Planet.Ocean.surfIceEOS['Ih'].fn_rho_kgm3(PIceI_MPa, TIceI_K, grid=False)
    Planet.Cp_JkgK[:Planet.Steps.nIceI] \
        = Planet.Ocean.surfIceEOS['Ih'].fn_Cp_JkgK(PIceI_MPa, TIceI_K, grid=False)
    Planet.alpha_pK[:Planet.Steps.nIceI] \
        = Planet.Ocean.surfIceEOS['Ih'].fn_alpha_pK(PIceI_MPa, TIceI_K, grid=False)
    Planet.kTherm_WmK[:Planet.Steps.nIceI] \
        = Planet.Ocean.surfIceEOS['Ih'].fn_kTherm_WmK(PIceI_MPa, TIceI_K, grid=False)

    if Planet.Do.POROUS_ICE:
        Planet.phi_frac[:Planet.Steps.nIceI] = Planet.Ocean.surfIceEOS['Ih'].fn_phi_frac(PIceI_MPa, TIceI_K, grid=False)
        Planet.rho_kgm3[:Planet.Steps.nIceI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.rhoMatrix_kgm3[:Planet.Steps.nIceI], 0, Planet.phi_frac[:Planet.Steps.nIceI], Planet.Ocean.Jrho)
        Planet.Cp_JkgK[:Planet.Steps.nIceI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Cp_JkgK[:Planet.Steps.nIceI] * Planet.rhoMatrix_kgm3[:Planet.Steps.nIceI],
            0, Planet.phi_frac[:Planet.Steps.nIceI], Planet.Ocean.JCp) / Planet.rho_kgm3[:Planet.Steps.nIceI]
        Planet.alpha_pK[:Planet.Steps.nIceI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.alpha_pK[:Planet.Steps.nIceI], 0, Planet.phi_frac[:Planet.Steps.nIceI], Planet.Ocean.Jalpha)
        Planet.kTherm_WmK[:Planet.Steps.nIceI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.kTherm_WmK[:Planet.Steps.nIceI], 0, Planet.phi_frac[:Planet.Steps.nIceI], Planet.Ocean.JkTherm)
    else:
        Planet.rho_kgm3[:Planet.Steps.nIceI] = Planet.rhoMatrix_kgm3[:Planet.Steps.nIceI] + 0.0

    # Evaluate thermodynamic properties of clathrate underplate layer
    Planet.rho_kgm3[Planet.Steps.nIceI:Planet.Steps.nIbottom] \
        = Planet.Ocean.surfIceEOS['Clath'].fn_rho_kgm3(Pclath_MPa[:-1], Tclath_K[:-1], grid=False)
    Planet.Cp_JkgK[Planet.Steps.nIceI:Planet.Steps.nIbottom] \
        = Planet.Ocean.surfIceEOS['Clath'].fn_Cp_JkgK(Pclath_MPa[:-1], Tclath_K[:-1], grid=False)
    Planet.alpha_pK[Planet.Steps.nIceI:Planet.Steps.nIbottom] \
        = Planet.Ocean.surfIceEOS['Clath'].fn_alpha_pK(Pclath_MPa[:-1], Tclath_K[:-1], grid=False)
    Planet.kTherm_WmK[Planet.Steps.nIceI:Planet.Steps.nIbottom] \
        = Planet.Ocean.surfIceEOS['Clath'].fn_kTherm_WmK(Pclath_MPa[:-1], Tclath_K[:-1], grid=False)
    Planet.rhoMatrix_kgm3[Planet.Steps.nIceI:Planet.Steps.nIbottom] = Planet.rho_kgm3[Planet.Steps.nIceI:Planet.Steps.nIbottom] + 0.0

    # Calculate remaining physical properties of upper ice and clathrates
    thisMAbove_kg = 0
    for i in range(1, Planet.Steps.nIbottom+1):
        # Increment depth based on change in pressure, combined with gravity and density
        Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
        Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
        Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
        thisMAbove_kg += Planet.MLayer_kg[i-1]
        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
        log.debug(f'il: {i:d}; P_MPa: {Planet.P_MPa[i]:.3f}; T_K: {Planet.T_K[i]:.3f}; phase: {Planet.phase[i]:d}')

    # Set actual thickness of clathrate underplate layer
    Planet.zClath_m = Planet.z_m[Planet.Steps.nIbottom] - Planet.z_m[Planet.Steps.nIceI]

    return Planet


def IceIIIConduct(Planet, Params):
    """ Calculate conductive profile for ice III layers
    """

    # Set linear P and adiabatic T in ice III layers
    PIceIII_MPa = np.linspace(Planet.P_MPa[Planet.Steps.nIbottom], Planet.PbIII_MPa, Planet.Steps.nIceIIILitho+1)
    PIceIIIratios = (PIceIII_MPa - PIceIII_MPa[0]) / (Planet.PbIII_MPa - PIceIII_MPa[0])
    TIceIII_K = Planet.Bulk.TbIII_K**(PIceIIIratios) * \
                Planet.T_K[Planet.Steps.nIbottom]**(1 - PIceIIIratios)
    Planet.P_MPa[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom+1] = PIceIII_MPa
    Planet.T_K[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom+1] = TIceIII_K

    # Get ice III EOS
    Planet.Ocean.surfIceEOS['III'] = IceEOSStruct(PIceIII_MPa, TIceIII_K, 'III')

    # Evaluate thermodynamic properties of upper ice III
    Planet.rho_kgm3[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] \
        = Planet.Ocean.surfIceEOS['III'].fn_rho_kgm3(PIceIII_MPa[:-1], TIceIII_K[:-1], grid=False)
    Planet.Cp_JkgK[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] \
        = Planet.Ocean.surfIceEOS['III'].fn_Cp_JkgK(PIceIII_MPa[:-1], TIceIII_K[:-1], grid=False)
    Planet.alpha_pK[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] \
        = Planet.Ocean.surfIceEOS['III'].fn_alpha_pK(PIceIII_MPa[:-1], TIceIII_K[:-1], grid=False)
    Planet.kTherm_WmK[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] \
        = Planet.Ocean.surfIceEOS['III'].fn_kTherm_WmK(PIceIII_MPa[:-1], TIceIII_K[:-1], grid=False)

    # Calculate remaining physical properties of upper ice III
    thisMAbove_kg = np.sum(Planet.MLayer_kg[:Planet.Steps.nIbottom])
    for i in range(Planet.Steps.nIbottom+1, Planet.Steps.nIIIbottom+1):
        # Increment depth based on change in pressure, combined with gravity and density
        Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
        Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
        Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
        thisMAbove_kg += Planet.MLayer_kg[i-1]
        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
        log.debug(f'il: {i:d}; P_MPa: {Planet.P_MPa[i]:.3f}; T_K: {Planet.T_K[i]:.3f}; phase: {Planet.phase[i]:d}')

    return Planet


def IceVConduct(Planet, Params):
    """ Calculate conductive profile for ice V layers
    """

    # Set linear P and adiabatic T in ice V layers
    PIceV_MPa = np.linspace(Planet.P_MPa[Planet.Steps.nIIIbottom], Planet.PbV_MPa, Planet.Steps.nIceVLitho+1)
    PIceVratios = (PIceV_MPa - PIceV_MPa[0]) / (Planet.PbV_MPa - PIceV_MPa[0])
    TIceV_K = Planet.Bulk.TbV_K**(PIceVratios) * \
                Planet.T_K[Planet.Steps.nIIIbottom]**(1 - PIceVratios)
    Planet.P_MPa[Planet.Steps.nIIIbottom:Planet.Steps.nSurfIce+1] = PIceV_MPa
    Planet.T_K[Planet.Steps.nIIIbottom:Planet.Steps.nSurfIce+1] = TIceV_K

    # Get ice V EOS
    Planet.Ocean.surfIceEOS['V'] = IceEOSStruct(PIceV_MPa, TIceV_K, 'V')

    # Evaluate thermodynamic properties of upper ice V
    Planet.rho_kgm3[Planet.Steps.nIIIbottom:Planet.Steps.nSurfIce] \
        = Planet.Ocean.surfIceEOS['V'].fn_rho_kgm3(PIceV_MPa[:-1], TIceV_K[:-1], grid=False)
    Planet.Cp_JkgK[Planet.Steps.nIIIbottom:Planet.Steps.nSurfIce] \
        = Planet.Ocean.surfIceEOS['V'].fn_Cp_JkgK(PIceV_MPa[:-1], TIceV_K[:-1], grid=False)
    Planet.alpha_pK[Planet.Steps.nIIIbottom:Planet.Steps.nSurfIce] \
        = Planet.Ocean.surfIceEOS['V'].fn_alpha_pK(PIceV_MPa[:-1], TIceV_K[:-1], grid=False)
    Planet.kTherm_WmK[Planet.Steps.nIIIbottom:Planet.Steps.nSurfIce] \
        = Planet.Ocean.surfIceEOS['V'].fn_kTherm_WmK(PIceV_MPa[:-1], TIceV_K[:-1], grid=False)

    # Calculate remaining physical properties of upper ice V
    thisMAbove_kg = np.sum(Planet.MLayer_kg[:Planet.Steps.nIIIbottom])
    for i in range(Planet.Steps.nIIIbottom+1, Planet.Steps.nSurfIce+1):
        # Increment depth based on change in pressure, combined with gravity and density
        Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
        Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
        Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
        thisMAbove_kg += Planet.MLayer_kg[i-1]
        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
        log.debug(f'il: {i:d}; P_MPa: {Planet.P_MPa[i]:.3f}; T_K: {Planet.T_K[i]:.3f}; phase: {Planet.phase[i]:d}')

    return Planet
