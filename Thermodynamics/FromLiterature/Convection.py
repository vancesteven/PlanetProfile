import numpy as np
from Utilities.dataStructs import Constants
from Thermodynamics.FromLiterature.HydroEOS import GetIceThermo
from Thermodynamics.FromLiterature.ThermalProfiles import ConductionClathLid, ConvectionDeschampsSotin2001, \
    kThermIsobaricAnderssonIbari2005, kThermHobbs1974, kThermMelinder2007

def IceIConvect(Planet, Params):
    """ Apply convection models from literature to determine thermal profile and
        state variables for possibly convecting ice layers

        Assigns Planet attributes:
            Tconv_K, etaConv_Pas, eLid_m, deltaTBL_m, QfromMantle_W, all physical layer arrays
    """

    if Params.VERBOSE: print('Applying solid-state convection to surface ice I based on Deschamps and Sotin (2001).')
    zbI_m = Planet.z_m[Planet.Steps.nIbottom-1]
    # Get "middle" pressure
    Pmid_MPa = (Planet.PbI_MPa - Planet.PbClath_MPa) / 2

    # Run calculations to get convection layer parameters
    Planet.Tconv_K, Planet.etaConv_Pas, Planet.eLid_m, Planet.deltaTBL_m, Planet.Ocean.QfromMantle_W, Planet.RaConvect = \
        ConvectionDeschampsSotin2001(Planet.Bulk.Tsurf_K, Planet.Bulk.R_m, Planet.kTherm_WmK[0], Planet.Bulk.Tb_K,
                                     zbI_m, Planet.g_ms2[0], Pmid_MPa, Planet.Ocean.comp, Planet.Ocean.wOcean_ppt,
                                     1, Planet.Do.EQUIL_Q)

    if Params.VERBOSE: print('Ice I convection parameters:\nT_convect = ' + str(round(Planet.Tconv_K,3)) + ' K,\n' +
                             'Viscosity etaConvect = {:.3e} Pa*s,\n'.format(Planet.etaConv_Pas) +
                             'Conductive lid thickness eLid_m = ' + str(round(Planet.eLid_m/1e3,1)) + ' km,\n' +
                             'Lower TBL thickness deltaTBL_m = ' + str(round(Planet.deltaTBL_m/1e3,1)) + ' km,\n' +
                             'Rayleigh number Ra = {:.3e}.'.format(Planet.RaConvect))

    if(Planet.zClath_m > Planet.eLid_m):
        Planet.Bulk.clathMaxDepth_m = Planet.eLid_m
        print('WARNING: Clathrate lid thickness was greater than the conductive lid thickness.' +
              'Planet.Bulk.clathMaxDepth_m has been reduced to be equal to the conductive lid thickness.')

    # Check for whole-lid conduction
    if(zbI_m <= Planet.eLid_m + Planet.deltaTBL_m):
        if Params.VERBOSE: print('Ice shell thickness ('+str(round(zbI_m/1e3,1))+' km) is less than that of the thermal'+
                                 ' boundary layers, convection is absent.\nApplying whole-shell conductive profile.')

        # Recalculate heat flux, as it will be too high for conduction-only:
        qSurf_Wm2 = (Planet.T_K[1] - Planet.Bulk.Tsurf_K) / (Planet.Bulk.R_m - Planet.r_m[1]) * Planet.kTherm_WmK[0]
        Planet.Ocean.QfromMantle_W = qSurf_Wm2 * 4*np.pi * Planet.Bulk.R_m**2

        # We leave the remaining quantities as initially assigned,
        # as we find the initial profile assuming conduction only.
    else:
        # Now model conductive + convective layers
        # Get layer transition indices
        nConduct = next(i[0] for i,val in np.ndenumerate(Planet.z_m) if val > Planet.eLid_m)
        nConvect = next(i[0] for i,val in np.ndenumerate(Planet.z_m) if val > zbI_m - Planet.deltaTBL_m) - nConduct
        indsTBL = range(nConduct + nConvect, Planet.Steps.nIbottom+1)
        if (Planet.Steps.nClath > nConduct):
            raise ValueError('Planet.Steps.nClath is greater than the fraction of surface ice layers used for the ' +
                             'upper conductive portion of the shell. Either increase Planet.Steps.eTBL_frac or ' +
                             'decrease Planet.Steps.nClath and re-run to solve this problem.')

        # Get pressure at the convecting transition
        PconvTop_MPa = Planet.P_MPa[nConduct]
        # Reset profile of upper layers, keeping pressure values fixed
        if Params.VERBOSE: print('Modeling ice I conduction in stagnant lid...')
        if Planet.Do.CLATHRATE:
            if Params.VERBOSE: print('Evaluating clathrate layers in stagnant lid.')
            ClathrateLayers(Planet, Params)
        else:
            thisMAbove_kg = 0
            # Reassign conductive profile with new bottom temperature for conductive layer
            PlidRatios = Planet.P_MPa[:nConduct+1]/PconvTop_MPa
            Planet.T_K[:nConduct+1] = Planet.Tconv_K**(PlidRatios) * Planet.T_K[0]**(1 - PlidRatios)

        # Get physical properties of upper conducting layer, and include 1 layer of convective layer for next step
        rhoCond_kgm3, CpCond_JkgK, alphaCond_pK, kThermCond_WmK = GetIceThermo(Planet.P_MPa[:nConduct+1],
                                                        Planet.T_K[:nConduct+1], Planet.phase[:nConduct+1])
        Planet.rho_kgm3[:nConduct+1] = rhoCond_kgm3
        Planet.Cp_JkgK[:nConduct+1] = CpCond_JkgK
        Planet.alpha_pK[:nConduct+1] = alphaCond_pK
        Planet.kTherm_WmK[:nConduct+1] = kThermCond_WmK

        for i in range(1, nConduct):
            # Increment depth based on change in pressure, combined with gravity and density
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

        if Params.VERBOSE: print('Stagnant lid conductive profile complete. Modeling ice I convecting layer...')

        for i in range(nConduct, nConduct + nConvect):
            # Increment depth based on change in pressure, combined with gravity and density
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2

            # Propagate adiabatic thermal profile
            Planet.T_K[i] = Planet.T_K[i-1] + Planet.T_K[i-1] * Planet.alpha_pK[i-1] / \
                            Planet.Cp_JkgK[i-1] / Planet.rho_kgm3[i-1] * (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6
            Planet.rho_kgm3[i], Planet.Cp_JkgK[i], Planet.alpha_pK[i], Planet.kTherm_WmK[i] \
                = GetIceThermo([Planet.P_MPa[i]], [Planet.T_K[i]], [Planet.phase[i]])
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

        if Params.VERBOSE: print('Convective profile complete. Modeling conduction in lower thermal boundary layer...')

        # Reassign conductive profile with new top temperature for conductive layer
        PTBLratios = Planet.P_MPa[indsTBL] / Planet.PbI_MPa
        Planet.T_K[indsTBL] = Planet.Bulk.Tb_K**(PTBLratios) * Planet.T_K[nConduct+nConvect-1]**(1 - PTBLratios)

        # Get physical properties of thermal boundary layer
        rhoTBL_kgm3, CpTBL_JkgK, alphaTBL_pK, kThermTBL_WmK = GetIceThermo(Planet.P_MPa[indsTBL[:-1]],
                                                Planet.T_K[indsTBL[:-1]], Planet.phase[indsTBL[:-1]])
        Planet.rho_kgm3[indsTBL[:-1]] = rhoTBL_kgm3
        Planet.Cp_JkgK[indsTBL[:-1]] = CpTBL_JkgK
        Planet.alpha_pK[indsTBL[:-1]] = alphaTBL_pK
        Planet.kTherm_WmK[indsTBL[:-1]] = kThermTBL_WmK

        for i in indsTBL:
            # Increment depth based on change in pressure, combined with gravity and density
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

    if Params.VERBOSE: print('Ice I convection calculations complete.')

    return Planet


def IceIIIConvect(Planet, Params):
    """ Apply convection models from literature to determine thermal profile and
        state variables for possibly convecting ice layers

        Assigns Planet attributes:
            TconvIII_K, etaConvIII_Pas, eLidIII_m, deltaTBLIII_m, QfromMantle_W, all physical layer arrays
    """

    if Params.VERBOSE: print('Applying solid-state convection to surface ice III based on Deschamps and Sotin (2001).')
    zbIII_m = Planet.z_m[Planet.Steps.nIIIbottom-1] - Planet.z_m[Planet.Steps.nIbottom-1]
    # Get "middle" pressure
    PmidIII_MPa = (Planet.PbIII_MPa - Planet.PbI_MPa) / 2

    # Run calculations to get convection layer parameters
    Planet.TconvIII_K, Planet.etaConvIII_Pas, Planet.eLidIII_m, Planet.deltaTBLIII_m, Planet.Ocean.QfromMantle_W,\
        Planet.RaConvectIII = ConvectionDeschampsSotin2001(Planet.Bulk.Tb_K, Planet.r_m[Planet.Steps.nIbottom],
                                     Planet.kTherm_WmK[Planet.Steps.nIbottom], Planet.Bulk.TbIII_K, zbIII_m,
                                     Planet.g_ms2[Planet.Steps.nIbottom], PmidIII_MPa, Planet.Ocean.comp,
                                     Planet.Ocean.wOcean_ppt, 3, Planet.Do.EQUIL_Q)

    if Params.VERBOSE: print('Ice III convection parameters:\nT_convectIII = ' + str(round(Planet.TconvIII_K,3)) + ' K,\n' +
                             'Viscosity etaConvectIII = {:.3e} Pa*s,\n'.format(Planet.etaConvIII_Pas) +
                             'Conductive lid thickness eLidIII_m = ' + str(round(Planet.eLidIII_m/1e3,1)) + ' km,\n' +
                             'Lower TBL thickness deltaTBLIII_m = ' + str(round(Planet.deltaTBLIII_m/1e3,1)) + ' km,\n' +
                             'Rayleigh number RaIII = {:.3e}.'.format(Planet.RaConvectIII))

    # Check for whole-lid conduction
    if(zbIII_m <= Planet.eLidIII_m + Planet.deltaTBLIII_m):
        if Params.VERBOSE: print('Underplate ice III thickness ('+str(round(zbIII_m/1e3,1))+' km) is less than that of the thermal'+
                                 ' boundary layers, convection is absent.\nApplying whole-layer conductive profile.')

        # We leave the remaining quantities as initially assigned,
        # as we find the initial profile assuming conduction only.
    else:
        # Now model conductive + convective layers
        # Get layer transition indices
        iConductEnd = next(i[0] for i,val in np.ndenumerate(Planet.z_m) if val > Planet.z_m[Planet.Steps.nIbottom] + Planet.eLidIII_m)
        iConvectEnd = next(i[0] for i,val in np.ndenumerate(Planet.z_m) if val > Planet.z_m[Planet.Steps.nIbottom] + zbIII_m - Planet.deltaTBLIII_m)
        indsTBL = range(iConvectEnd, Planet.Steps.nIIIbottom+1)

        # Get pressure at the convecting transition
        PconvTopIII_MPa = Planet.P_MPa[iConductEnd]
        # Reset profile of upper layers, keeping pressure values fixed
        if Params.VERBOSE: print('Modeling ice III conduction in stagnant lid...')

        thisMAbove_kg = np.sum(Planet.MLayer_kg[:Planet.Steps.nIbottom])
        # Reassign conductive profile with new bottom temperature for conductive layer
        PlidRatios = Planet.P_MPa[Planet.Steps.nIbottom:iConductEnd+1]/PconvTopIII_MPa
        Planet.T_K[Planet.Steps.nIbottom:iConductEnd+1] = Planet.TconvIII_K**(PlidRatios) * Planet.T_K[Planet.Steps.nIbottom]**(1 - PlidRatios)

        # Get physical properties of upper conducting layer, and include 1 layer of convective layer for next step
        rhoCondIII_kgm3, CpCondIII_JkgK, alphaCondIII_pK, kThermCondIII_WmK \
            = GetIceThermo(Planet.P_MPa[Planet.Steps.nIbottom:iConductEnd+1],
                           Planet.T_K[Planet.Steps.nIbottom:iConductEnd+1],
                           Planet.phase[Planet.Steps.nIbottom:iConductEnd+1])
        Planet.rho_kgm3[Planet.Steps.nIbottom:iConductEnd+1] = rhoCondIII_kgm3
        Planet.Cp_JkgK[Planet.Steps.nIbottom:iConductEnd+1] = CpCondIII_JkgK
        Planet.alpha_pK[Planet.Steps.nIbottom:iConductEnd+1] = alphaCondIII_pK
        Planet.kTherm_WmK[Planet.Steps.nIbottom:iConductEnd+1] = kThermCondIII_WmK

        for i in range(Planet.Steps.nIbottom+1, iConductEnd):
            # Increment depth based on change in pressure, combined with gravity and density
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

        if Params.VERBOSE: print('Stagnant lid conductive profile complete. Modeling ice III convecting layer...')

        for i in range(iConductEnd, iConvectEnd):
            # Increment depth based on change in pressure, combined with gravity and density
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2

            # Propagate adiabatic thermal profile
            Planet.T_K[i] = Planet.T_K[i-1] + Planet.T_K[i-1] * Planet.alpha_pK[i-1] / \
                            Planet.Cp_JkgK[i-1] / Planet.rho_kgm3[i-1] * (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6
            Planet.rho_kgm3[i], Planet.Cp_JkgK[i], Planet.alpha_pK[i], Planet.kTherm_WmK[i] \
                = GetIceThermo([Planet.P_MPa[i]], [Planet.T_K[i]], [Planet.phase[i]])
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

        if Params.VERBOSE: print('Convective profile complete. Modeling conduction in lower thermal boundary layer...')

        if(Planet.T_K[iConvectEnd-1] > Planet.Bulk.TbIII_K):
            raise ValueError('Ice III bottom temperature of ' + str(round(Planet.Bulk.TbIII_K,3)) +
                             ' K is less than the temperature at the lower TBL transition of ' +
                             str(round(Planet.T_K[iConvectEnd-1],3)) + ' K. Try increasing TbIII_K ' +
                             'to create a more realistic thermal profile.')

        # Reassign conductive profile with new top temperature for conductive layer
        PTBLratios = Planet.P_MPa[indsTBL] / Planet.PbIII_MPa
        Planet.T_K[indsTBL] = Planet.Bulk.TbIII_K**(PTBLratios) * Planet.T_K[iConvectEnd-1]**(1 - PTBLratios)

        # Get physical properties of thermal boundary layer
        rhoTBLIII_kgm3, CpTBLIII_JkgK, alphaTBLIII_pK, kThermTBLIII_WmK \
            = GetIceThermo(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], Planet.phase[indsTBL[:-1]])
        Planet.rho_kgm3[indsTBL[:-1]] = rhoTBLIII_kgm3
        Planet.Cp_JkgK[indsTBL[:-1]] = CpTBLIII_JkgK
        Planet.alpha_pK[indsTBL[:-1]] = alphaTBLIII_pK
        Planet.kTherm_WmK[indsTBL[:-1]] = kThermTBLIII_WmK

        for i in indsTBL:
            # Increment depth based on change in pressure, combined with gravity and density
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

    if Params.VERBOSE: print('Ice III convection calculations complete.')

    return Planet


def IceVConvect(Planet, Params):
    """ Apply convection models from literature to determine thermal profile and
        state variables for possibly convecting ice layers

        Assigns Planet attributes:
            Tconv_K, etaConv_Pas, eLid_m, deltaTBL_m, QfromMantle_W, all physical layer arrays
    """

    if Params.VERBOSE: print('Applying solid-state convection to surface ice I based on Deschamps and Sotin (2001).')
    zbI_m = Planet.z_m[Planet.Steps.nIbottom-1]
    # Get "middle" pressure
    Pmid_MPa = (Planet.PbI_MPa - Planet.PbClath_MPa) / 2
    # Get thermal conductivity at surface
    if Planet.Do.CLATHRATE:
        Planet.kTherm_WmK[0] = Constants.kClath_WmK
    else:
        #Planet.kTherm_WmK[0] = kThermIsobaricAnderssonIbari2005(Planet.Bulk.Tsurf_K, 1)
        #Planet.kTherm_WmK[0] = kThermHobbs1974(Planet.Bulk.Tsurf_K)
        Planet.kTherm_WmK[0] = kThermMelinder2007(Planet.Bulk.Tsurf_K, Constants.T0)

    # Run calculations to get convection layer parameters
    Planet.Tconv_K, Planet.etaConv_Pas, Planet.eLid_m, Planet.deltaTBL_m, Planet.Ocean.QfromMantle_W, Planet.RaConvect = \
        ConvectionDeschampsSotin2001(Planet.Bulk.Tsurf_K, Planet.Bulk.R_m, Planet.kTherm_WmK[0],
                                     Planet.Bulk.Tb_K, zbI_m, Planet.g_ms2[0], Pmid_MPa, 1, Planet.Do.EQUIL_Q)

    if Params.VERBOSE: print('Convection parameters:\nT_convect = ' + str(round(Planet.Tconv_K,3)) + ' K,\n' +
                             'Viscosity etaConvect = {:.3e} Pa*s,\n'.format(Planet.etaConv_Pas) +
                             'Conductive lid thickness eLid_m = ' + str(round(Planet.eLid_m/1e3,1)) + ' km,\n' +
                             'Lower TBL thickness deltaTBL_m = ' + str(round(Planet.deltaTBL_m/1e3,1)) + ' km,\n' +
                             'Rayleigh number Ra = {:.3e}.'.format(Planet.RaConvect))

    if(Planet.zClath_m > Planet.eLid_m):
        Planet.Bulk.clathMaxDepth_m = Planet.eLid_m
        print('WARNING: Clathrate lid thickness was greater than the conductive lid thickness.' +
              'Planet.Bulk.clathMaxDepth_m has been reduced to be equal to the conductive lid thickness.')

    # Check for whole-lid conduction
    if(zbI_m <= Planet.eLid_m + Planet.deltaTBL_m):
        if Params.VERBOSE: print('Ice shell thickness ('+str(round(zbI_m/1e3,1))+' km) is less than that of the thermal'+
                                 ' boundary layers, convection is absent.\nApplying whole-shell conductive profile.')

        # Recalculate heat flux, as it will be too high for conduction-only:
        qSurf_Wm2 = (Planet.T_K[1] - Planet.Bulk.Tsurf_K) / (Planet.Bulk.R_m - Planet.r_m[1]) * Planet.kTherm_WmK[0]
        Planet.Ocean.QfromMantle_W = qSurf_Wm2 * 4*np.pi * Planet.Bulk.R_m**2

        # We again need to do clathrate layers first if there are any
        if Planet.Do.CLATHRATE:
            if Params.VERBOSE: print('Evaluating clathrate layers.')
            ClathrateLayers(Planet, Params)
        else:
            thisMAbove_kg = 0

        Planet.z_m[Planet.Steps.nClath] = Planet.zClath_m + 0.0
        Planet.r_m[Planet.Steps.nClath] = Planet.Bulk.R_m - Planet.z_m[Planet.Steps.nClath]
        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        Planet.g_ms2[Planet.Steps.nClath] = Constants.G * thisMBelow_kg / Planet.r_m[Planet.Steps.nClath]**2
        qTop_Wm2 = Planet.Ocean.QfromMantle_W / 4/np.pi/Planet.r_m[Planet.Steps.nClath]**2

        # Evaluate thermodynamic properties at top of ice I with SeaFreeze
        Planet.rho_kgm3[Planet.Steps.nClath], Planet.Cp_JkgK[Planet.Steps.nClath], Planet.alpha_pK[Planet.Steps.nClath] \
            = GetIceThermo([Planet.P_MPa[Planet.Steps.nClath]], [Planet.T_K[Planet.Steps.nClath]], [Planet.phase[Planet.Steps.nClath]])

        for i in range(Planet.Steps.nClath+1, Planet.Steps.nIbottom):
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
            Planet.kTherm_WmK[i] = Planet.kTherm_WmK[i-1]  # Placeholder until a self-consistent determination is implemented

            Planet.T_K[i], qTop_Wm2 = ConductiveTemperature(Planet.T_K[i-1], Planet.r_m[i-1], Planet.r_m[i],
                                                            Planet.kTherm_WmK[i], Planet.rho_kgm3[i], 0, 0, qTop_Wm2)
            Planet.rho_kgm3[i], Planet.Cp_JkgK[i], Planet.alpha_pK[i] = GetIceThermo([Planet.P_MPa[i]], [Planet.T_K[i]],
                                                                                      [Planet.phase[i]])
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

        # Force final temperature value to be the assigned value
        Planet.T_K[Planet.Steps.nIbottom] = Planet.Bulk.Tb_K

    else:
        # Now model conductive + convective layers
        # Get layer transition indices
        nConduct = next(i[0] for i,val in np.ndenumerate(Planet.z_m) if val > Planet.eLid_m)
        nConvect = next(i[0] for i,val in np.ndenumerate(Planet.z_m) if val > zbI_m - Planet.deltaTBL_m) - nConduct
        indsTBL = range(nConduct + nConvect, Planet.Steps.nIbottom)
        if (Planet.Steps.nClath > nConduct):
            raise ValueError('Planet.Steps.nClath is greater than the fraction of surface ice layers used for the ' +
                             'upper conductive portion of the shell. Either increase Planet.Steps.eTBL_frac or ' +
                             'decrease Planet.Steps.nClath and re-run to solve this problem.')

        # Get pressure at the convecting transition
        PconvTop_MPa = Planet.P_MPa[nConduct]
        # Reset profile of upper layers, keeping pressure values fixed
        if Params.VERBOSE: print('Modeling ice I conduction in stagnant lid...')
        if Planet.Do.CLATHRATE:
            if Params.VERBOSE: print('Evaluating clathrate layers in stagnant lid.')
            ClathrateLayers(Planet, Params)
        else:
            thisMAbove_kg = 0
            # Reassign conductive profile with new bottom temperature for conductive layer
            PlidRatios = Planet.P_MPa[:nConduct+1]/PconvTop_MPa
            Planet.T_K[:nConduct+1] = Planet.Tconv_K**(PlidRatios) * Planet.T_K[0]**(1 - PlidRatios)

        # Get physical properties of upper conducting layer, and include 1 layer of convective layer for next step
        rhoCond_kgm3, CpCond_JkgK, alphaCond_pK = GetIceThermo(Planet.P_MPa[:nConduct+1], Planet.T_K[:nConduct+1],
                                                               Planet.phase[:nConduct+1])
        Planet.rho_kgm3[:nConduct+1] = rhoCond_kgm3
        Planet.Cp_JkgK[:nConduct+1] = CpCond_JkgK
        Planet.alpha_pK[:nConduct+1] = alphaCond_pK

        for i in range(1, nConduct):
            # Increment depth based on change in pressure, combined with gravity and density
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

        if Params.VERBOSE: print('Stagnant lid conductive profile complete. Modeling ice I convecting layer...')

        for i in range(nConduct + 1, nConduct + nConvect):
            # Increment depth based on change in pressure, combined with gravity and density
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2

            # Propagate adiabatic thermal profile
            Planet.T_K[i] = Planet.T_K[i-1] + Planet.T_K[i-1] * Planet.alpha_pK[i-1] / \
                            Planet.Cp_JkgK[i-1] / Planet.rho_kgm3[i-1] * (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6
            Planet.rho_kgm3[i], Planet.Cp_JkgK[i], Planet.alpha_pK[i] \
                = GetIceThermo([Planet.P_MPa[i]], [Planet.T_K[i]], [Planet.phase[i]])
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

        if Params.VERBOSE: print('Convective profile complete. Modeling conduction in lower thermal boundary layer...')

        # Reassign conductive profile with new top temperature for conductive layer
        PTBLratios = Planet.P_MPa[indsTBL] / Planet.PbI_MPa
        Planet.T_K[indsTBL] = Planet.Bulk.Tb_K**(PTBLratios) * Planet.T_K[nConduct+nConvect-1]**(1 - PTBLratios)

        # Get physical properties of thermal boundary layer
        rhoTBL_kgm3, CpTBL_JkgK, alphaTBL_pK = GetIceThermo(Planet.P_MPa[indsTBL], Planet.T_K[indsTBL],
                                                            Planet.phase[indsTBL])
        Planet.rho_kgm3[indsTBL] = rhoTBL_kgm3
        Planet.Cp_JkgK[indsTBL] = CpTBL_JkgK
        Planet.alpha_pK[indsTBL] = alphaTBL_pK

        for i in indsTBL:
            # Increment depth based on change in pressure, combined with gravity and density
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

    if Params.VERBOSE: print('Convection calculations complete.')

    return Planet


