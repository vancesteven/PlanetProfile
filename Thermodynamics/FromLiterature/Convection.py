import numpy as np
from Utilities.dataStructs import Constants
from Thermodynamics.FromLiterature.HydroEOS import PhaseConv
from Thermodynamics.FromLiterature.ThermalProfiles import ConvectionDeschampsSotin2001, \
    kThermIsobaricAnderssonIbari2005, kThermHobbs1974, kThermMelinder2007

def IceIConvect(Planet, Params):
    """ Apply convection models from literature to determine thermal profile and
        state variables for possibly convecting ice I layers

        Assigns Planet attributes:
            Tconv_K, etaConv_Pas, eLid_m, deltaTBL_m, QfromMantle_W, all physical layer arrays
    """

    if Params.VERBOSE: print('Applying solid-state convection to surface ice I based on Deschamps and Sotin (2001).')
    zbI_m = Planet.z_m[Planet.Steps.nIbottom-1]
    # Get "middle" pressure
    Pmid_MPa = (Planet.PbI_MPa + Planet.P_MPa[0]) / 2
    # Get a lower estimate for thermal conductivity if there's no clathrate lid, to be more consistent
    # with the Matlab implementation of Deschamps and Sotin (2001). The thermal conductivity of clathrates
    # is essentially fixed at 0.5 W/(m K), so doing this operation when the top layer is clathrate
    # gives us what we want in that case, too.
    Planet.kTherm_WmK[0] = Planet.Ocean.surfIceEOS[PhaseConv(Planet.phase[0])].fn_kTherm_WmK(Pmid_MPa, Planet.Bulk.Tb_K)

    # Run calculations to get convection layer parameters
    Planet.Tconv_K, Planet.etaConv_Pas, Planet.eLid_m, Planet.deltaTBL_m, Planet.Ocean.QfromMantle_W, Planet.RaConvect = \
        ConvectionDeschampsSotin2001(Planet.T_K[0], Planet.r_m[0], Planet.kTherm_WmK[0], Planet.Bulk.Tb_K,
                                     zbI_m, Planet.g_ms2[0], Pmid_MPa, Planet.Ocean.EOS,
                                     Planet.Ocean.surfIceEOS['Ih'], 1, Planet.Do.EQUIL_Q)

    if Params.VERBOSE: print('Ice I convection parameters:\nT_convect = ' + str(round(Planet.Tconv_K,3)) + ' K,\n' +
                             'Viscosity etaConvect = {:.3e} Pa*s,\n'.format(Planet.etaConv_Pas) +
                             'Conductive lid thickness eLid_m = ' + str(round(Planet.eLid_m/1e3,1)) + ' km,\n' +
                             'Lower TBL thickness deltaTBL_m = ' + str(round(Planet.deltaTBL_m/1e3,1)) + ' km,\n' +
                             'Rayleigh number Ra = {:.3e}.'.format(Planet.RaConvect))

    # Check for whole-lid conduction
    if(zbI_m <= Planet.eLid_m + Planet.deltaTBL_m):
        if Params.VERBOSE: print('Ice shell thickness ('+str(round(zbI_m/1e3,1))+' km) is less than that of the thermal'+
                                 ' boundary layers, convection is absent.\nApplying whole-shell conductive profile.')

        # Recalculate heat flux, as it will be too high for conduction-only:
        qSurf_Wm2 = (Planet.T_K[1] - Planet.T_K[0]) / (Planet.r_m[0] - Planet.r_m[1]) * Planet.kTherm_WmK[0]
        Planet.Ocean.QfromMantle_W = qSurf_Wm2 * 4*np.pi * Planet.Bulk.R_m**2

        # We leave the remaining quantities as initially assigned,
        # as we find the initial profile assuming conduction only.
    else:
        # Now model conductive + convective layers
        # Get layer transition indices from previous profile
        try:
            nConduct = next(i[0] for i,val in np.ndenumerate(Planet.z_m) if val > Planet.eLid_m)
        except StopIteration:
            raise RuntimeError('Failed to find any depth indices for upper TBL of ice III. Try increasing Steps.nIceVLitho.')
        try:
            nConvect = next(i[0] for i,val in np.ndenumerate(Planet.z_m) if val > zbI_m - Planet.deltaTBL_m) - nConduct
        except StopIteration:
            raise RuntimeError('Failed to find any depth indices for lower TBL of ice III. Try increasing Steps.nIceVLitho.')
        indsTBL = range(nConduct + nConvect, Planet.Steps.nIbottom+1)
        # Get pressure at the convecting transition
        PconvTop_MPa = Planet.P_MPa[nConduct]

        # Reset profile of upper layers, keeping pressure values fixed
        if Planet.Do.CLATHRATE:
            if Params.VERBOSE: print('Evaluating clathrate layers in stagnant lid.')

            if Planet.Bulk.clathType == 'top':
                if (Planet.eLid_m < Planet.Bulk.zClath_m):
                    Planet.Bulk.clathMaxDepth_m = Planet.eLid_m
                    if Params.VERBOSE: print('Clathrate lid thickness was greater than the conductive lid thickness. ' +
                                             'Planet.Bulk.clathMaxDepth_m has been reduced to be equal to the conductive lid thickness.')
                if Planet.PbClathMax_MPa > PconvTop_MPa:
                    Planet.PbClathMax_MPa = PconvTop_MPa
                    Planet.TclathTrans_K = Planet.Tconv_K
                    Planet.P_MPa[:Planet.Steps.nClath] = np.linspace(Planet.P_MPa[0], PconvTop_MPa, Planet.Steps.nClath+1)[:-1]
                    Planet.P_MPa[Planet.Steps.nClath:Planet.Steps.nIbottom+1] = \
                        np.linspace(PconvTop_MPa, Planet.PbI_MPa, Planet.Steps.nIceI+1)
                    # Reassign T profile to be consistent with conduction
                    PlidRatios = Planet.P_MPa[:Planet.Steps.nClath+1] / PconvTop_MPa
                    Planet.T_K[:Planet.Steps.nClath+1] = Planet.Tconv_K**(PlidRatios) * Planet.T_K[0]**(1 - PlidRatios)

                    # Reset nConduct/nConvect/indsTBL to account for the index shift moving clathrates to be above the
                    # transition to ice I
                    nConduct = Planet.Steps.nClath
                    nConvect = next(i[0] for i,val in np.ndenumerate(Planet.z_m) if val > zbI_m - Planet.deltaTBL_m) - nConduct
                    indsTBL = range(nConduct + nConvect, Planet.Steps.nIbottom+1)
                else:
                    print('WARNING: Max. clathrate layer thickness is less than that of the stagnant lid. Lid thickness ' +
                          'was calculated using a constant thermal conductivity equal to that of clathrates at the ' +
                          'surface, so the properties of the ice I conductive layer between the clathrate lid and ' +
                          'convective region are likely to be physically inconsistent (i.e., an unrealistically low ' +
                          'thermal conductivity). Increase Bulk.clathMaxDepth_m to greater than ' +
                          str(round(Planet.eLid_m/1e3,3)) + ' km for these run settings to avoid this problem.')
                    # Keep existing transition (P,T) from clathrates to ice I
                    Planet.P_MPa[:Planet.Steps.nClath] = np.linspace(Planet.P_MPa[0], Planet.PbClathMax_MPa, Planet.Steps.nClath+1)[:-1]
                    Planet.P_MPa[Planet.Steps.nClath:Planet.Steps.nIbottom+1] = \
                        np.linspace(Planet.PbClathMax_MPa, Planet.PbI_MPa, Planet.Steps.nIceI+1)
                    # Model conduction in ice I between clathrate lid and convective region
                    PlidRatiosClath = (Planet.P_MPa[:Planet.Steps.nClath+1] - Planet.P_MPa[0]) / (Planet.PbClathMax_MPa - Planet.P_MPa[0])
                    Planet.T_K[:Planet.Steps.nClath+1] = Planet.TclathTrans_K**(PlidRatiosClath) * Planet.T_K[0]**(1 - PlidRatiosClath)
                    PlidRatiosIceI = (Planet.P_MPa[Planet.Steps.nClath:nConduct+1] - Planet.PbClathMax_MPa) / (PconvTop_MPa - Planet.PbClathMax_MPa)
                    Planet.T_K[Planet.Steps.nClath:nConduct+1] = Planet.Tconv_K**(PlidRatiosIceI) * Planet.TclathTrans_K**(1 - PlidRatiosIceI)

                # Get physical properties of clathrate lid
                Planet.rho_kgm3[:Planet.Steps.nClath] \
                    = Planet.Ocean.surfIceEOS['Clath'].fn_rho_kgm3(Planet.P_MPa[:Planet.Steps.nClath],
                                                                   Planet.T_K[:Planet.Steps.nClath], grid=False)
                Planet.Cp_JkgK[:Planet.Steps.nClath] \
                    = Planet.Ocean.surfIceEOS['Clath'].fn_Cp_JkgK(Planet.P_MPa[:Planet.Steps.nClath],
                                                                  Planet.T_K[:Planet.Steps.nClath], grid=False)
                Planet.alpha_pK[:Planet.Steps.nClath] \
                    = Planet.Ocean.surfIceEOS['Clath'].fn_alpha_pK(Planet.P_MPa[:Planet.Steps.nClath],
                                                                   Planet.T_K[:Planet.Steps.nClath], grid=False)
                Planet.kTherm_WmK[:Planet.Steps.nClath] \
                    = Planet.Ocean.surfIceEOS['Clath'].fn_kTherm_WmK(Planet.P_MPa[:Planet.Steps.nClath],
                                                                     Planet.T_K[:Planet.Steps.nClath], grid=False)
            else:
                raise ValueError('IceIConvect behavior is not defined for Bulk.clathType "' +
                                 Planet.Bulk.clathType + '".')
        else:
            if Params.VERBOSE: print('Modeling ice I conduction in stagnant lid...')
            # Reassign conductive profile with new bottom temperature for conductive layer
            PlidRatios = (Planet.P_MPa[:nConduct+1] - Planet.P_MPa[0]) / (PconvTop_MPa - Planet.P_MPa[0])
            Planet.T_K[:nConduct+1] = Planet.Tconv_K**(PlidRatios) * Planet.T_K[0]**(1 - PlidRatios)

        # Get physical properties of upper conducting layer, and include 1 layer of convective layer for next step
        Planet.rho_kgm3[Planet.Steps.nClath:nConduct+1] \
            = Planet.Ocean.surfIceEOS['Ih'].fn_rho_kgm3(Planet.P_MPa[Planet.Steps.nClath:nConduct+1],
                                                        Planet.T_K[Planet.Steps.nClath:nConduct+1], grid=False)
        Planet.Cp_JkgK[Planet.Steps.nClath:nConduct+1] \
            = Planet.Ocean.surfIceEOS['Ih'].fn_Cp_JkgK(Planet.P_MPa[Planet.Steps.nClath:nConduct+1],
                                                       Planet.T_K[Planet.Steps.nClath:nConduct+1], grid=False)
        Planet.alpha_pK[Planet.Steps.nClath:nConduct+1] \
            = Planet.Ocean.surfIceEOS['Ih'].fn_alpha_pK(Planet.P_MPa[Planet.Steps.nClath:nConduct+1],
                                                        Planet.T_K[Planet.Steps.nClath:nConduct+1], grid=False)
        Planet.kTherm_WmK[Planet.Steps.nClath:nConduct+1] \
            = Planet.Ocean.surfIceEOS['Ih'].fn_kTherm_WmK(Planet.P_MPa[Planet.Steps.nClath:nConduct+1],
                                                          Planet.T_K[Planet.Steps.nClath:nConduct+1], grid=False)

        thisMAbove_kg = 0
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
            Planet.rho_kgm3[i] = Planet.Ocean.surfIceEOS['Ih'].fn_rho_kgm3(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            Planet.Cp_JkgK[i] = Planet.Ocean.surfIceEOS['Ih'].fn_Cp_JkgK(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            Planet.alpha_pK[i] = Planet.Ocean.surfIceEOS['Ih'].fn_alpha_pK(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            Planet.kTherm_WmK[i] = Planet.Ocean.surfIceEOS['Ih'].fn_kTherm_WmK(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

        if Params.VERBOSE: print('Convective profile complete. Modeling conduction in lower thermal boundary layer...')

        # Reassign conductive profile with new top temperature for conductive layer
        PTBLratios = (Planet.P_MPa[indsTBL] - Planet.P_MPa[nConduct+nConvect-1]) / (Planet.PbI_MPa - Planet.P_MPa[nConduct+nConvect-1])
        Planet.T_K[indsTBL] = Planet.Bulk.Tb_K**(PTBLratios) * Planet.T_K[nConduct+nConvect-1]**(1 - PTBLratios)

        # Get physical properties of thermal boundary layer
        Planet.rho_kgm3[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['Ih'].fn_rho_kgm3(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)
        Planet.Cp_JkgK[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['Ih'].fn_Cp_JkgK(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)
        Planet.alpha_pK[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['Ih'].fn_alpha_pK(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)
        Planet.kTherm_WmK[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['Ih'].fn_kTherm_WmK(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)

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
    zbIII_m = Planet.z_m[Planet.Steps.nIIIbottom-1] - Planet.z_m[Planet.Steps.nIbottom]
    # Get "middle" pressure
    PmidIII_MPa = (Planet.PbIII_MPa + Planet.PbI_MPa) / 2

    # Run calculations to get convection layer parameters
    Planet.TconvIII_K, Planet.etaConvIII_Pas, Planet.eLidIII_m, Planet.deltaTBLIII_m, Planet.Ocean.QfromMantle_W,\
        Planet.RaConvectIII = ConvectionDeschampsSotin2001(Planet.Bulk.Tb_K, Planet.r_m[Planet.Steps.nIbottom],
                                     Planet.kTherm_WmK[Planet.Steps.nIbottom], Planet.Bulk.TbIII_K, zbIII_m,
                                     Planet.g_ms2[Planet.Steps.nIbottom], PmidIII_MPa, Planet.Ocean.EOS,
                                     Planet.Ocean.surfIceEOS['III'], 3, Planet.Do.EQUIL_Q)

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
        PlidRatios = (Planet.P_MPa[Planet.Steps.nIbottom:iConductEnd+1] - Planet.P_MPa[Planet.Steps.nIbottom]) / \
                     (PconvTopIII_MPa - Planet.P_MPa[Planet.Steps.nIbottom])
        Planet.T_K[Planet.Steps.nIbottom:iConductEnd+1] = Planet.TconvIII_K**(PlidRatios) * Planet.T_K[Planet.Steps.nIbottom]**(1 - PlidRatios)

        # Get physical properties of upper conducting layer, and include 1 layer of convective layer for next step
        Planet.rho_kgm3[Planet.Steps.nIbottom:iConductEnd+1] \
            = Planet.Ocean.surfIceEOS['III'].fn_rho_kgm3(Planet.P_MPa[Planet.Steps.nIbottom:iConductEnd+1],
                                                          Planet.T_K[Planet.Steps.nIbottom:iConductEnd+1], grid=False)
        Planet.Cp_JkgK[Planet.Steps.nIbottom:iConductEnd+1] \
            = Planet.Ocean.surfIceEOS['III'].fn_Cp_JkgK(Planet.P_MPa[Planet.Steps.nIbottom:iConductEnd+1],
                                                          Planet.T_K[Planet.Steps.nIbottom:iConductEnd+1], grid=False)
        Planet.alpha_pK[Planet.Steps.nIbottom:iConductEnd+1] \
            = Planet.Ocean.surfIceEOS['III'].fn_alpha_pK(Planet.P_MPa[Planet.Steps.nIbottom:iConductEnd+1],
                                                          Planet.T_K[Planet.Steps.nIbottom:iConductEnd+1], grid=False)
        Planet.kTherm_WmK[Planet.Steps.nIbottom:iConductEnd+1] \
            = Planet.Ocean.surfIceEOS['III'].fn_kTherm_WmK(Planet.P_MPa[Planet.Steps.nIbottom:iConductEnd+1],
                                                          Planet.T_K[Planet.Steps.nIbottom:iConductEnd+1], grid=False)

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
            Planet.rho_kgm3[i] = Planet.Ocean.surfIceEOS['III'].fn_rho_kgm3(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            Planet.Cp_JkgK[i] = Planet.Ocean.surfIceEOS['III'].fn_Cp_JkgK(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            Planet.alpha_pK[i] = Planet.Ocean.surfIceEOS['III'].fn_alpha_pK(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            Planet.kTherm_WmK[i] = Planet.Ocean.surfIceEOS['III'].fn_kTherm_WmK(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

        if Params.VERBOSE: print('Convective profile complete. Modeling conduction in lower thermal boundary layer...')

        if(Planet.T_K[iConvectEnd-1] > Planet.Bulk.TbIII_K):
            raise ValueError('Ice III bottom temperature of ' + str(round(Planet.Bulk.TbIII_K,3)) +
                             ' K is less than the temperature at the lower TBL transition of ' +
                             str(round(Planet.T_K[iConvectEnd-1],3)) + ' K. Try increasing Bulk.TbIII_K ' +
                             ' or decreasing Bulk.Tb_K to create a more realistic thermal profile.')

        # Reassign conductive profile with new top temperature for conductive layer

        PTBLratios = (Planet.P_MPa[indsTBL] - Planet.P_MPa[iConvectEnd-1]) / (Planet.PbIII_MPa - Planet.P_MPa[iConvectEnd-1])
        Planet.T_K[indsTBL] = Planet.Bulk.TbIII_K**(PTBLratios) \
                              * Planet.T_K[iConvectEnd-1]**(1 - PTBLratios)

        # Get physical properties of thermal boundary layer
        Planet.rho_kgm3[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['III'].fn_rho_kgm3(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)
        Planet.Cp_JkgK[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['III'].fn_Cp_JkgK(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)
        Planet.alpha_pK[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['III'].fn_alpha_pK(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)
        Planet.kTherm_WmK[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['III'].fn_kTherm_WmK(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)

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
            TconvV_K, etaConvV_Pas, eLidV_m, deltaTBLV_m, QfromMantle_W, all physical layer arrays
    """

    if Params.VERBOSE: print('Applying solid-state convection to surface ice V based on Deschamps and Sotin (2001).')
    zbV_m = Planet.z_m[Planet.Steps.nSurfIce-1] - Planet.z_m[Planet.Steps.nIIIbottom]
    # Get "middle" pressure
    PmidV_MPa = (Planet.PbV_MPa + Planet.PbIII_MPa) / 2

    # Run calculations to get convection layer parameters
    Planet.TconvV_K, Planet.etaConvV_Pas, Planet.eLidV_m, Planet.deltaTBLV_m, Planet.Ocean.QfromMantle_W,\
        Planet.RaConvectV = ConvectionDeschampsSotin2001(Planet.Bulk.TbIII_K, Planet.r_m[Planet.Steps.nIIIbottom],
                                     Planet.kTherm_WmK[Planet.Steps.nIIIbottom], Planet.Bulk.TbV_K, zbV_m,
                                     Planet.g_ms2[Planet.Steps.nIIIbottom], PmidV_MPa, Planet.Ocean.EOS,
                                     Planet.Ocean.surfIceEOS['V'], 5, Planet.Do.EQUIL_Q)

    if Params.VERBOSE: print('Ice V convection parameters:\nT_convectV = ' + str(round(Planet.TconvV_K,3)) + ' K,\n' +
                             'Viscosity etaConvectV = {:.3e} Pa*s,\n'.format(Planet.etaConvV_Pas) +
                             'Conductive lid thickness eLidV_m = ' + str(round(Planet.eLidV_m/1e3,1)) + ' km,\n' +
                             'Lower TBL thickness deltaTBLV_m = ' + str(round(Planet.deltaTBLV_m/1e3,1)) + ' km,\n' +
                             'Rayleigh number RaV = {:.3e}.'.format(Planet.RaConvectV))

    # Check for whole-lid conduction
    if(zbV_m <= Planet.eLidV_m + Planet.deltaTBLV_m):
        if Params.VERBOSE: print('Underplate ice V thickness ('+str(round(zbV_m/1e3,1))+' km) is less than that of the thermal'+
                                 ' boundary layers, convection is absent.\nApplying whole-layer conductive profile.')

        # We leave the remaining quantities as initially assigned,
        # as we find the initial profile assuming conduction only.
    else:
        # Now model conductive + convective layers
        # Get layer transition indices
        try:
            iConductEnd = next(i[0] for i,val in np.ndenumerate(Planet.z_m) if val > Planet.z_m[Planet.Steps.nIIIbottom] + Planet.eLidV_m)
        except StopIteration:
            raise RuntimeError('Failed to find any depth indices for upper TBL of ice V. Try increasing Steps.nIceVLitho.')
        try:
            iConvectEnd = next(i[0] for i,val in np.ndenumerate(Planet.z_m) if val > Planet.z_m[Planet.Steps.nIIIbottom] + zbV_m - Planet.deltaTBLV_m)
        except StopIteration:
            raise RuntimeError('Failed to find any depth indices for lower TBL of ice V. Try increasing Steps.nIceVLitho.')
        indsTBL = range(iConvectEnd, Planet.Steps.nSurfIce+1)

        # Get pressure at the convecting transition
        PconvTopV_MPa = Planet.P_MPa[iConductEnd]
        # Reset profile of upper layers, keeping pressure values fixed
        if Params.VERBOSE: print('Modeling ice V conduction in stagnant lid...')

        thisMAbove_kg = np.sum(Planet.MLayer_kg[:Planet.Steps.nIIIbottom])
        # Reassign conductive profile with new bottom temperature for conductive layer
        PlidRatios = (Planet.P_MPa[Planet.Steps.nIIIbottom:iConductEnd+1] - Planet.P_MPa[Planet.Steps.nIIIbottom]) / \
                     (PconvTopV_MPa - Planet.P_MPa[Planet.Steps.nIIIbottom])
        Planet.T_K[Planet.Steps.nIIIbottom:iConductEnd+1] = Planet.TconvV_K**(PlidRatios) \
                                                            * Planet.T_K[Planet.Steps.nIIIbottom]**(1 - PlidRatios)

        # Get physical properties of upper conducting layer, and include 1 layer of convective layer for next step
        Planet.rho_kgm3[Planet.Steps.nIIIbottom:iConductEnd+1] \
            = Planet.Ocean.surfIceEOS['V'].fn_rho_kgm3(Planet.P_MPa[Planet.Steps.nIIIbottom:iConductEnd+1],
                                                        Planet.T_K[Planet.Steps.nIIIbottom:iConductEnd+1], grid=False)
        Planet.Cp_JkgK[Planet.Steps.nIIIbottom:iConductEnd+1] \
            = Planet.Ocean.surfIceEOS['V'].fn_Cp_JkgK(Planet.P_MPa[Planet.Steps.nIIIbottom:iConductEnd+1],
                                                        Planet.T_K[Planet.Steps.nIIIbottom:iConductEnd+1], grid=False)
        Planet.alpha_pK[Planet.Steps.nIIIbottom:iConductEnd+1] \
            = Planet.Ocean.surfIceEOS['V'].fn_alpha_pK(Planet.P_MPa[Planet.Steps.nIIIbottom:iConductEnd+1],
                                                        Planet.T_K[Planet.Steps.nIIIbottom:iConductEnd+1], grid=False)
        Planet.kTherm_WmK[Planet.Steps.nIIIbottom:iConductEnd+1] \
            = Planet.Ocean.surfIceEOS['V'].fn_kTherm_WmK(Planet.P_MPa[Planet.Steps.nIIIbottom:iConductEnd+1],
                                                        Planet.T_K[Planet.Steps.nIIIbottom:iConductEnd+1], grid=False)

        for i in range(Planet.Steps.nIIIbottom+1, iConductEnd):
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

        if Params.VERBOSE: print('Stagnant lid conductive profile complete. Modeling ice V convecting layer...')

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
            Planet.rho_kgm3[i] = Planet.Ocean.surfIceEOS['V'].fn_rho_kgm3(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            Planet.Cp_JkgK[i] = Planet.Ocean.surfIceEOS['V'].fn_Cp_JkgK(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            Planet.alpha_pK[i] = Planet.Ocean.surfIceEOS['V'].fn_alpha_pK(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            Planet.kTherm_WmK[i] = Planet.Ocean.surfIceEOS['V'].fn_kTherm_WmK(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

        if Params.VERBOSE: print('Convective profile complete. Modeling conduction in lower thermal boundary layer...')

        if(Planet.T_K[iConvectEnd-1] > Planet.Bulk.TbV_K):
            raise ValueError('Ice V bottom temperature of ' + str(round(Planet.Bulk.TbV_K,3)) +
                             ' K is less than the temperature at the lower TBL transition of ' +
                             str(round(Planet.T_K[iConvectEnd-1],3)) + ' K. Try increasing TbV_K ' +
                             'to create a more realistic thermal profile.')

        # Reassign conductive profile with new top temperature for conductive layer
        PTBLratios = (Planet.P_MPa[indsTBL] - Planet.P_MPa[iConvectEnd-1]) / (Planet.PbV_MPa - Planet.P_MPa[iConvectEnd-1])
        Planet.T_K[indsTBL] = Planet.Bulk.TbV_K**(PTBLratios) * Planet.T_K[iConvectEnd-1]**(1 - PTBLratios)

        # Get physical properties of thermal boundary layer
        Planet.rho_kgm3[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['V'].fn_rho_kgm3(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)
        Planet.Cp_JkgK[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['V'].fn_Cp_JkgK(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)
        Planet.alpha_pK[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['V'].fn_alpha_pK(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)
        Planet.kTherm_WmK[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['V'].fn_kTherm_WmK(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)

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

    if Params.VERBOSE: print('Ice V convection calculations complete.')

    return Planet


def ClathShellConvect(Planet, Params):
    """ Apply convection models from literature to determine thermal profile and
        state variables for possibly convecting clathrate layers when the whole
        shell is made of clathrates.

        Assigns Planet attributes:
            Tconv_K, etaConv_Pas, eLid_m, deltaTBL_m, QfromMantle_W, all physical layer arrays
    """

    if Params.VERBOSE: print('Applying solid-state convection to surface clathrates based on Deschamps and Sotin (2001).')
    zbI_m = Planet.z_m[Planet.Steps.nIbottom-1]
    # Get "middle" pressure
    Pmid_MPa = (Planet.PbI_MPa - Planet.P_MPa[0]) / 2
    Planet.kTherm_WmK[0] = Planet.Ocean.surfIceEOS['Clath'].fn_kTherm_WmK(Planet.P_MPa[0], Planet.Bulk.Tb_K)

    # Run calculations to get convection layer parameters
    Planet.Tconv_K, Planet.etaConv_Pas, Planet.eLid_m, Planet.deltaTBL_m, Planet.Ocean.QfromMantle_W, Planet.RaConvect = \
        ConvectionDeschampsSotin2001(Planet.T_K[0], Planet.r_m[0], Planet.kTherm_WmK[0], Planet.Bulk.Tb_K,
                                     zbI_m, Planet.g_ms2[0], Pmid_MPa, Planet.Ocean.surfIceEOS['Clath'],
                                     Planet.Ocean.surfIceEOS['Clath'], Constants.phaseClath, Planet.Do.EQUIL_Q)

    if Params.VERBOSE: print('Clathrate shell convection parameters:\nT_convect = ' + str(round(Planet.Tconv_K,3)) + ' K,\n' +
                             'Viscosity etaConvect = {:.3e} Pa*s,\n'.format(Planet.etaConv_Pas) +
                             'Conductive lid thickness eLid_m = ' + str(round(Planet.eLid_m/1e3,1)) + ' km,\n' +
                             'Lower TBL thickness deltaTBL_m = ' + str(round(Planet.deltaTBL_m/1e3,1)) + ' km,\n' +
                             'Rayleigh number Ra = {:.3e}.'.format(Planet.RaConvect))

    # Check for whole-lid conduction
    if(zbI_m <= Planet.eLid_m + Planet.deltaTBL_m):
        if Params.VERBOSE: print('Ice shell thickness ('+str(round(zbI_m/1e3,1))+' km) is less than that of the thermal'+
                                 ' boundary layers, convection is absent.\nApplying whole-shell conductive profile.')

        # Recalculate heat flux, as it will be too high for conduction-only:
        qSurf_Wm2 = (Planet.T_K[1] - Planet.T_K[0]) / (Planet.r_m[0] - Planet.r_m[1]) * Planet.kTherm_WmK[0]
        Planet.Ocean.QfromMantle_W = qSurf_Wm2 * 4*np.pi * Planet.Bulk.R_m**2

        # We leave the remaining quantities as initially assigned,
        # as we find the initial profile assuming conduction only.
    else:
        # Now model conductive + convective layers
        # Get layer transition indices from previous profile
        nConduct = next(i[0] for i,val in np.ndenumerate(Planet.z_m) if val > Planet.eLid_m)
        nConvect = next(i[0] for i,val in np.ndenumerate(Planet.z_m) if val > zbI_m - Planet.deltaTBL_m) - nConduct
        indsTBL = range(nConduct + nConvect, Planet.Steps.nIbottom+1)
        # Get pressure at the convecting transition
        PconvTop_MPa = Planet.P_MPa[nConduct]

        # Reset profile of upper layers, keeping pressure values fixed
        if Params.VERBOSE: print('Modeling clathrate conduction in stagnant lid...')
        # Reassign conductive profile with new bottom temperature for conductive layer
        PlidRatios = (Planet.P_MPa[:nConduct+1] - Planet.P_MPa[0]) / (PconvTop_MPa - Planet.P_MPa[0])
        Planet.T_K[:nConduct+1] = Planet.Tconv_K**(PlidRatios) * Planet.T_K[0]**(1 - PlidRatios)

        # Get physical properties of upper conducting layer, and include 1 layer of convective layer for next step
        Planet.rho_kgm3[:nConduct+1] \
            = Planet.Ocean.surfIceEOS['Clath'].fn_rho_kgm3(Planet.P_MPa[:nConduct+1],
                                                           Planet.T_K[:nConduct+1], grid=False)
        Planet.Cp_JkgK[:nConduct+1] \
            = Planet.Ocean.surfIceEOS['Clath'].fn_Cp_JkgK(Planet.P_MPa[:nConduct+1],
                                                          Planet.T_K[:nConduct+1], grid=False)
        Planet.alpha_pK[:nConduct+1] \
            = Planet.Ocean.surfIceEOS['Clath'].fn_alpha_pK(Planet.P_MPa[:nConduct+1],
                                                           Planet.T_K[:nConduct+1], grid=False)
        Planet.kTherm_WmK[:nConduct+1] \
            = Planet.Ocean.surfIceEOS['Clath'].fn_kTherm_WmK(Planet.P_MPa[:nConduct+1],
                                                             Planet.T_K[:nConduct+1], grid=False)

        thisMAbove_kg = 0
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
            Planet.rho_kgm3[i] = Planet.Ocean.surfIceEOS['Clath'].fn_rho_kgm3(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            Planet.Cp_JkgK[i] = Planet.Ocean.surfIceEOS['Clath'].fn_Cp_JkgK(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            Planet.alpha_pK[i] = Planet.Ocean.surfIceEOS['Clath'].fn_alpha_pK(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            Planet.kTherm_WmK[i] = Planet.Ocean.surfIceEOS['Clath'].fn_kTherm_WmK(Planet.P_MPa[i], Planet.T_K[i], grid=False)
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

        if Params.VERBOSE: print('Convective profile complete. Modeling conduction in lower thermal boundary layer...')

        # Reassign conductive profile with new top temperature for conductive layer
        PTBLratios = (Planet.P_MPa[indsTBL] - Planet.P_MPa[nConduct+nConvect-1]) / (Planet.PbI_MPa - Planet.P_MPa[nConduct+nConvect-1])
        Planet.T_K[indsTBL] = Planet.Bulk.Tb_K**(PTBLratios) * Planet.T_K[nConduct+nConvect-1]**(1 - PTBLratios)

        # Get physical properties of thermal boundary layer
        Planet.rho_kgm3[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['Clath'].fn_rho_kgm3(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)
        Planet.Cp_JkgK[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['Clath'].fn_Cp_JkgK(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)
        Planet.alpha_pK[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['Clath'].fn_alpha_pK(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)
        Planet.kTherm_WmK[indsTBL[:-1]] \
            = Planet.Ocean.surfIceEOS['Clath'].fn_kTherm_WmK(Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]], grid=False)

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

    if Params.VERBOSE: print('Clathrate convection calculations complete.')

    return Planet
