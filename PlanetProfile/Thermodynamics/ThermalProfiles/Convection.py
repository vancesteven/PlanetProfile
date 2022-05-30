import numpy as np
import logging
from PlanetProfile.Thermodynamics.Geophysical import PropagateConduction, EvalLayerProperties, \
    PorosityCorrectionVacIce, PorosityCorrectionFilledIce, PropagateAdiabaticSolid, \
    PropagateAdiabaticPorousVacIce, PropagateAdiabaticPorousFilledIce
from PlanetProfile.Thermodynamics.HydroEOS import PhaseConv
from PlanetProfile.Thermodynamics.ThermalProfiles.ThermalProfiles import ConvectionDeschampsSotin2001, \
    kThermIsobaricAnderssonIbari2005, kThermHobbs1974, kThermMelinder2007
from PlanetProfile.Utilities.defineStructs import Constants

# Assign logger
log = logging.getLogger('PlanetProfile')

def IceIConvectSolid(Planet, Params):
    """ Apply convection models from literature to determine thermal profile and
        state variables for possibly convecting ice I layers

        Assigns Planet attributes:
            Tconv_K, etaConv_Pas, eLid_m, deltaTBL_m, QfromMantle_W, all physical layer arrays
    """

    log.debug('Applying solid-state convection to surface ice I based on Deschamps and Sotin (2001).')
    zbI_m = Planet.z_m[Planet.Steps.nIbottom]
    # Get "middle" pressure
    Pmid_MPa = (Planet.PbI_MPa + Planet.P_MPa[0]) / 2
    # Get a lower estimate for thermal conductivity if there's no clathrate lid, to be more consistent
    # with the Matlab implementation of Deschamps and Sotin (2001). The thermal conductivity of clathrates
    # is essentially fixed at 0.5 W/(m K), so doing this operation when the top layer is clathrate
    # gives us what we want in that case, too.
    phaseTop = PhaseConv(Planet.phase[0])
    Planet.kTherm_WmK[0] = Planet.Ocean.surfIceEOS[phaseTop].fn_kTherm_WmK(Pmid_MPa, Planet.Bulk.Tb_K)

    # Run calculations to get convection layer parameters
    Planet.Tconv_K, Planet.etaConv_Pas, Planet.eLid_m, Planet.Dconv_m, Planet.deltaTBL_m, Planet.Ocean.QfromMantle_W, \
        Planet.RaConvect, Planet.RaCrit = \
        ConvectionDeschampsSotin2001(Planet.T_K[0], Planet.r_m[0], Planet.kTherm_WmK[0], Planet.Bulk.Tb_K,
                                     zbI_m, Planet.g_ms2[0], Pmid_MPa, Planet.Ocean.EOS,
                                     Planet.Ocean.surfIceEOS['Ih'], 1, Planet.Do.EQUIL_Q)

    log.debug(f'Ice I convection parameters:\n    T_convect = {Planet.Tconv_K:.3f} K,\n' +
              f'    Viscosity etaConvect = {Planet.etaConv_Pas:.3e} Pa*s,\n' +
              f'    Conductive lid thickness eLid = {Planet.eLid_m/1e3:.1f} km,\n' +
              f'    Convecting layer thickness Dconv = {Planet.Dconv_m/1e3:.1f} km,\n' +
              f'    Lower TBL thickness deltaTBL = {Planet.deltaTBL_m/1e3:.1f} km,\n' +
              f'    Rayleigh number Ra = {Planet.RaConvect:.3e}.')

    # Check for whole-lid conduction
    if(zbI_m <= Planet.eLid_m + Planet.deltaTBL_m):
        log.info(f'Ice shell thickness ({zbI_m/1e3:.1f} km) is less than that of the thermal ' +
                  'boundary layers--convection is absent. Applying whole-shell conductive profile.')
        Planet.eLid_m = zbI_m
        Planet.Dconv_m = 0.0
        Planet.deltaTBL_m = 0.0

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
            log.debug('Evaluating clathrate layers in conductive lid.')

            if Planet.Bulk.clathType == 'top':
                if (Planet.eLid_m < Planet.zClath_m):
                    Planet.Bulk.clathMaxDepth_m = Planet.eLid_m
                    log.debug('Clathrate lid thickness was greater than the conductive lid thickness. ' +
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
                    log.warning('Max. clathrate layer thickness is less than that of the stagnant lid. Lid thickness ' +
                                'was calculated using a constant thermal conductivity equal to that of clathrates at the ' +
                                'surface, so the properties of the ice I conductive layer between the clathrate lid and ' +
                                'convective region are likely to be physically inconsistent (i.e., an unrealistically low ' +
                                'thermal conductivity). Increase Bulk.clathMaxDepth_m to greater than ' +
                               f'{Planet.eLid_m/1e3:.3f} km for these run settings to avoid this problem.')
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
                Planet = EvalLayerProperties(Planet, Params, 0, Planet.Steps.nClath, Planet.Ocean.surfIceEOS['Clath'],
                                                  Planet.P_MPa[:Planet.Steps.nClath], Planet.T_K[:Planet.Steps.nClath])

                Planet.rho_kgm3[:Planet.Steps.nClath] = Planet.rhoMatrix_kgm3[:Planet.Steps.nClath] + 0.0

            else:
                raise ValueError(f'IceIConvect behavior is not defined for Bulk.clathType "{Planet.Bulk.clathType}".')
        else:
            log.debug('Modeling ice I conduction in stagnant lid...')
            # Reassign conductive profile with new bottom temperature for conductive layer
            PlidRatios = (Planet.P_MPa[:nConduct+1] - Planet.P_MPa[0]) / (PconvTop_MPa - Planet.P_MPa[0])
            Planet.T_K[:nConduct+1] = Planet.Tconv_K**(PlidRatios) * Planet.T_K[0]**(1 - PlidRatios)

        # Get physical properties of upper conducting layer, and include 1 layer of convective layer for next step
        Planet = EvalLayerProperties(Planet, Params, Planet.Steps.nClath, nConduct+1, Planet.Ocean.surfIceEOS['Ih'],
                                          Planet.P_MPa[Planet.Steps.nClath:nConduct+1], Planet.T_K[Planet.Steps.nClath:nConduct+1])
        Planet.rho_kgm3[Planet.Steps.nClath:nConduct+1] = Planet.rhoMatrix_kgm3[Planet.Steps.nClath:nConduct+1] + 0.0

        Planet = PropagateConduction(Planet, Params, 0, nConduct-1)
        log.debug('Stagnant lid conductive profile complete. Modeling ice I convecting layer...')

        Planet = PropagateAdiabaticSolid(Planet, Params, nConduct, nConduct+nConvect, Planet.Ocean.surfIceEOS['Ih'])
        log.debug('Convective profile complete. Modeling conduction in lower thermal boundary layer...')

        # Reassign conductive profile with new top temperature for conductive layer
        PTBLratios = (Planet.P_MPa[indsTBL] - Planet.P_MPa[nConduct+nConvect-1]) / (Planet.PbI_MPa - Planet.P_MPa[nConduct+nConvect-1])
        Planet.T_K[indsTBL] = Planet.Bulk.Tb_K**(PTBLratios) * Planet.T_K[nConduct+nConvect-1]**(1 - PTBLratios)

        # Get physical properties of thermal boundary layer
        Planet = EvalLayerProperties(Planet, Params, indsTBL[0], indsTBL[-1]+1,
                                     Planet.Ocean.surfIceEOS['Ih'],
                                     Planet.P_MPa[indsTBL], Planet.T_K[indsTBL])
        Planet.rho_kgm3[indsTBL] = Planet.rhoMatrix_kgm3[indsTBL] + 0.0

        # Apply conductive profile to lower TBL
        Planet = PropagateConduction(Planet, Params, indsTBL[0]-1, indsTBL[-1])

    log.debug('Ice I convection calculations complete.')

    return Planet


def IceIConvectPorous(Planet, Params):
    """ Apply convection models from literature to determine thermal profile and
        state variables for possibly convecting ice I layers

        Assigns Planet attributes:
            Tconv_K, etaConv_Pas, eLid_m, deltaTBL_m, QfromMantle_W, all physical layer arrays
    """

    log.debug('Applying solid-state convection to surface ice I based on Deschamps and Sotin (2001).')
    zbI_m = Planet.z_m[Planet.Steps.nIbottom]
    # Get "middle" pressure
    Pmid_MPa = (Planet.PbI_MPa + Planet.P_MPa[0]) / 2
    # Get a lower estimate for thermal conductivity if there's no clathrate lid, to be more consistent
    # with the Matlab implementation of Deschamps and Sotin (2001). The thermal conductivity of clathrates
    # is essentially fixed at 0.5 W/(m K), so doing this operation when the top layer is clathrate
    # gives us what we want in that case, too.
    phaseTop = PhaseConv(Planet.phase[0])
    Planet.kTherm_WmK[0] = Planet.Ocean.surfIceEOS[phaseTop].fn_kTherm_WmK(Pmid_MPa, Planet.Bulk.Tb_K)
    Planet.kTherm_WmK[0] = Planet.Ocean.surfIceEOS[phaseTop].fn_porosCorrect(Planet.kTherm_WmK[0], 0,
                           Planet.Ocean.surfIceEOS[phaseTop].fn_phi_frac(Pmid_MPa, Planet.Bulk.Tb_K),
                           Planet.Ocean.JkTherm)

    # Run calculations to get convection layer parameters
    Planet.Tconv_K, Planet.etaConv_Pas, Planet.eLid_m, Planet.Dconv_m, Planet.deltaTBL_m, Planet.Ocean.QfromMantle_W, \
        Planet.RaConvect, Planet.RaCrit = \
        ConvectionDeschampsSotin2001(Planet.T_K[0], Planet.r_m[0], Planet.kTherm_WmK[0], Planet.Bulk.Tb_K,
                                     zbI_m, Planet.g_ms2[0], Pmid_MPa, Planet.Ocean.EOS,
                                     Planet.Ocean.surfIceEOS['Ih'], 1, Planet.Do.EQUIL_Q)

    log.debug(f'Ice I convection parameters:\n    T_convect = {Planet.Tconv_K:.3f} K,\n' +
              f'    Viscosity etaConvect = {Planet.etaConv_Pas:.3e} Pa*s,\n' +
              f'    Conductive lid thickness eLid = {Planet.eLid_m/1e3:.1f} km,\n' +
              f'    Convecting layer thickness Dconv = {Planet.Dconv_m/1e3:.1f} km,\n' +
              f'    Lower TBL thickness deltaTBL = {Planet.deltaTBL_m/1e3:.1f} km,\n' +
              f'    Rayleigh number Ra = {Planet.RaConvect:.3e}.')

    # Check for whole-lid conduction
    if(zbI_m <= Planet.eLid_m + Planet.deltaTBL_m):
        log.info(f'Ice shell thickness ({zbI_m/1e3:.1f} km) is less than that of the thermal ' +
                  'boundary layers--convection is absent. Applying whole-shell conductive profile.')
        Planet.eLid_m = zbI_m
        Planet.Dconv_m = 0.0
        Planet.deltaTBL_m = 0.0
        indsTBL = [1, 0]

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
            log.debug('Evaluating clathrate layers in stagnant lid.')

            if Planet.Bulk.clathType == 'top':
                if (Planet.eLid_m < Planet.zClath_m):
                    Planet.Bulk.clathMaxDepth_m = Planet.eLid_m
                    log.debug('Clathrate lid thickness was greater than the conductive lid thickness. ' +
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
                    log.warning('Max. clathrate layer thickness is less than that of the stagnant lid. Lid thickness ' +
                                'was calculated using a constant thermal conductivity equal to that of clathrates at the ' +
                                'surface, so the properties of the ice I conductive layer between the clathrate lid and ' +
                                'convective region are likely to be physically inconsistent (i.e., an unrealistically low ' +
                                'thermal conductivity). Increase Bulk.clathMaxDepth_m to greater than ' +
                               f'{Planet.eLid_m/1e3:.3f} km for these run settings to avoid this problem.')
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
                Planet = EvalLayerProperties(Planet, Params, 0, Planet.Steps.nClath,
                                             Planet.Ocean.surfIceEOS['Clath'],
                                             Planet.P_MPa[:Planet.Steps.nClath],
                                             Planet.T_K[:Planet.Steps.nClath])
                # Correct for porosity in clathrate layers
                Planet = PorosityCorrectionVacIce(Planet, Params, 0, Planet.Steps.nClath,
                                                  Planet.Ocean.surfIceEOS['Clath'],
                                                  Planet.P_MPa[:Planet.Steps.nClath],
                                                  Planet.T_K[:Planet.Steps.nClath])

            else:
                raise ValueError(f'IceIConvect behavior is not defined for Bulk.clathType "{Planet.Bulk.clathType}".')
        else:
            log.debug('Modeling ice I conduction in stagnant lid...')
            # Reassign conductive profile with new bottom temperature for conductive layer
            PlidRatios = (Planet.P_MPa[:nConduct+1] - Planet.P_MPa[0]) / (PconvTop_MPa - Planet.P_MPa[0])
            Planet.T_K[:nConduct+1] = Planet.Tconv_K**(PlidRatios) * Planet.T_K[0]**(1 - PlidRatios)

        # Get physical properties of upper conducting layer, and include 1 layer of convective layer for next step
        Planet = EvalLayerProperties(Planet, Params, Planet.Steps.nClath, nConduct+1,
                                     Planet.Ocean.surfIceEOS['Ih'],
                                     Planet.P_MPa[Planet.Steps.nClath:nConduct+1],
                                     Planet.T_K[Planet.Steps.nClath:nConduct+1])
        # Correct for porosity in ice I layers
        Planet = PorosityCorrectionVacIce(Planet, Params, Planet.Steps.nClath, nConduct+1,
                                          Planet.Ocean.surfIceEOS['Ih'],
                                          Planet.P_MPa[Planet.Steps.nClath:nConduct+1],
                                          Planet.T_K[Planet.Steps.nClath:nConduct+1])

        Planet = PropagateConduction(Planet, Params, 0, nConduct-1)
        log.debug('Stagnant lid conductive profile complete. Modeling ice I convecting layer...')

        # Propagate adiabatic thermal profile
        Planet = PropagateAdiabaticPorousVacIce(Planet, Params, nConduct, nConduct+nConvect,
                                                Planet.Ocean.surfIceEOS['Ih'])
        log.debug('Convective profile complete. Modeling conduction in lower thermal boundary layer...')

        # Reassign conductive profile with new top temperature for conductive layer
        PTBLratios = (Planet.P_MPa[indsTBL] - Planet.P_MPa[nConduct+nConvect-1]) / (Planet.PbI_MPa - Planet.P_MPa[nConduct+nConvect-1])
        Planet.T_K[indsTBL] = Planet.Bulk.Tb_K**(PTBLratios) * Planet.T_K[nConduct+nConvect-1]**(1 - PTBLratios)

        # Get physical properties of thermal boundary layer
        Planet = EvalLayerProperties(Planet, Params, indsTBL[0], indsTBL[-1],
                                     Planet.Ocean.surfIceEOS['Ih'],
                                     Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]])

        # Correct for porosity in ice I layers, assuming pores may be treated as vacuum

        Planet = PorosityCorrectionVacIce(Planet, Params, indsTBL[0], indsTBL[-1]+1,
                                          Planet.Ocean.surfIceEOS['Ih'],
                                          Planet.P_MPa[indsTBL], Planet.T_K[indsTBL])

    # Apply conductive profile to lower TBL
    Planet = PropagateConduction(Planet, Params, indsTBL[0]-1, indsTBL[-1])
    log.debug('Ice I convection calculations complete.')

    return Planet


def IceIIIConvectSolid(Planet, Params):
    """ Apply convection models from literature to determine thermal profile and
        state variables for possibly convecting ice layers

        Assigns Planet attributes:
            TconvIII_K, etaConvIII_Pas, eLidIII_m, deltaTBLIII_m, QfromMantle_W, all physical layer arrays
    """

    log.debug('Applying solid-state convection to surface ice III based on Deschamps and Sotin (2001).')
    zbIII_m = Planet.z_m[Planet.Steps.nIIIbottom] - Planet.z_m[Planet.Steps.nIbottom]
    # Get "middle" pressure
    PmidIII_MPa = (Planet.PbIII_MPa + Planet.PbI_MPa) / 2

    # Run calculations to get convection layer parameters
    Planet.TconvIII_K, Planet.etaConvIII_Pas, Planet.eLidIII_m, Planet.DconvIII_m, Planet.deltaTBLIII_m, Planet.Ocean.QfromMantle_W, \
        Planet.RaConvectIII, Planet.RaCritIII = ConvectionDeschampsSotin2001(Planet.Bulk.Tb_K, Planet.r_m[Planet.Steps.nIbottom],
                                     Planet.kTherm_WmK[Planet.Steps.nIbottom], Planet.Bulk.TbIII_K, zbIII_m,
                                     Planet.g_ms2[Planet.Steps.nIbottom], PmidIII_MPa, Planet.Ocean.EOS,
                                     Planet.Ocean.surfIceEOS['III'], 3, Planet.Do.EQUIL_Q)

    log.debug(f'Ice III convection parameters:\n    T_convectIII = {Planet.TconvIII_K:.3f} K,\n' +
              f'    Viscosity etaConvectIII = {Planet.etaConvIII_Pas:.3e} Pa*s,\n' +
              f'    Conductive lid thickness eLidIII = {Planet.eLidIII_m/1e3:.1f} km,\n' +
              f'    Convecting layer thickness DconvIII = {Planet.DconvIII_m/1e3:.1f} km,\n' +
              f'    Lower TBL thickness deltaTBLIII = {Planet.deltaTBLIII_m/1e3:.1f} km,\n' +
              f'    Rayleigh number RaIII = {Planet.RaConvectIII:.3e}.')

    # Check for whole-lid conduction
    if(zbIII_m <= Planet.eLidIII_m + Planet.deltaTBLIII_m):
        log.info(f'Underplate ice III thickness ({zbIII_m/1e3:.1f} km) is less than that of the thermal ' +
                  'boundary layers--convection is absent. Applying whole-layer conductive profile.')
        Planet.eLidIII_m = zbIII_m
        Planet.DconvIII_m = 0.0
        Planet.deltaTBLIII_m = 0.0

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
        log.debug('Modeling ice III conduction in stagnant lid...')

        thisMAbove_kg = np.sum(Planet.MLayer_kg[:Planet.Steps.nIbottom])
        # Reassign conductive profile with new bottom temperature for conductive layer
        PlidRatios = (Planet.P_MPa[Planet.Steps.nIbottom:iConductEnd+1] - Planet.P_MPa[Planet.Steps.nIbottom]) / \
                     (PconvTopIII_MPa - Planet.P_MPa[Planet.Steps.nIbottom])
        Planet.T_K[Planet.Steps.nIbottom:iConductEnd+1] = Planet.TconvIII_K**(PlidRatios) * Planet.T_K[Planet.Steps.nIbottom]**(1 - PlidRatios)

        # Get physical properties of upper conducting layer, and include 1 layer of convective layer for next step
        Planet = EvalLayerProperties(Planet, Params, Planet.Steps.nIbottom, iConductEnd+1,
                                    Planet.Ocean.surfIceEOS['III'],
                                    Planet.P_MPa[Planet.Steps.nIbottom:iConductEnd+1],
                                    Planet.T_K[Planet.Steps.nIbottom:iConductEnd+1])

        Planet = PropagateConduction(Planet, Params, Planet.Steps.nIbottom, iConductEnd-1)

        log.debug('Stagnant lid conductive profile complete. Modeling ice III convecting layer...')

        Planet = PropagateAdiabaticSolid(Planet, Params, iConductEnd, iConvectEnd,
                                         Planet.Ocean.surfIceEOS['III'])

        log.debug('Convective profile complete. Modeling conduction in lower thermal boundary layer...')

        if(Planet.T_K[iConvectEnd-1] > Planet.Bulk.TbIII_K):
            raise ValueError(f'Ice III bottom temperature of {Planet.Bulk.TbIII_K:.3f} K ' +
                              'is less than the temperature at the lower TBL transition of ' +
                             f'{Planet.T_K[iConvectEnd-1]:.3f} K. Try increasing Bulk.TbIII_K ' +
                              'or decreasing Bulk.Tb_K to create a more realistic thermal profile.')

        # Reassign conductive profile with new top temperature for conductive layer

        PTBLratios = (Planet.P_MPa[indsTBL] - Planet.P_MPa[iConvectEnd-1]) / (Planet.PbIII_MPa - Planet.P_MPa[iConvectEnd-1])
        Planet.T_K[indsTBL] = Planet.Bulk.TbIII_K**(PTBLratios) \
                              * Planet.T_K[iConvectEnd-1]**(1 - PTBLratios)

        # Get physical properties of thermal boundary layer
        Planet = EvalLayerProperties(Planet, Params, indsTBL[0], indsTBL[-1]+1,
                                    Planet.Ocean.surfIceEOS['III'],
                                    Planet.P_MPa[indsTBL], Planet.T_K[indsTBL])

        Planet = PropagateConduction(Planet, Params, indsTBL[0]-1, indsTBL[-1])

    log.debug('Ice III convection calculations complete.')

    return Planet


def IceIIIConvectPorous(Planet, Params):
    """ Apply convection models from literature to determine thermal profile and
        state variables for possibly convecting ice layers

        Assigns Planet attributes:
            TconvIII_K, etaConvIII_Pas, eLidIII_m, deltaTBLIII_m, QfromMantle_W, all physical layer arrays
    """

    log.debug('Applying solid-state convection to surface ice III based on Deschamps and Sotin (2001).')
    zbIII_m = Planet.z_m[Planet.Steps.nIIIbottom] - Planet.z_m[Planet.Steps.nIbottom]
    # Get "middle" pressure
    PmidIII_MPa = (Planet.PbIII_MPa + Planet.PbI_MPa) / 2
    # Porosity is unlikely to change the rough kTherm estimate we use here, so we do not
    # correct for porosity, unlike in the ice Ih/clathrate functions.

    # Run calculations to get convection layer parameters
    Planet.TconvIII_K, Planet.etaConvIII_Pas, Planet.eLidIII_m, Planet.DconvIII_m, Planet.deltaTBLIII_m, Planet.Ocean.QfromMantle_W, \
        Planet.RaConvectIII, Planet.RaCritIII = ConvectionDeschampsSotin2001(Planet.Bulk.Tb_K, Planet.r_m[Planet.Steps.nIbottom],
                                     Planet.kTherm_WmK[Planet.Steps.nIbottom], Planet.Bulk.TbIII_K, zbIII_m,
                                     Planet.g_ms2[Planet.Steps.nIbottom], PmidIII_MPa, Planet.Ocean.EOS,
                                     Planet.Ocean.surfIceEOS['III'], 3, Planet.Do.EQUIL_Q)

    log.debug(f'Ice III convection parameters:\n    T_convectIII = {Planet.TconvIII_K:.3f} K,\n' +
              f'    Viscosity etaConvectIII = {Planet.etaConvIII_Pas:.3e} Pa*s,\n' +
              f'    Conductive lid thickness eLidIII = {Planet.eLidIII_m/1e3:.1f} km,\n' +
              f'    Convecting layer thickness DconvIII = {Planet.DconvIII_m/1e3:.1f} km,\n' +
              f'    Lower TBL thickness deltaTBLIII = {Planet.deltaTBLIII_m/1e3:.1f} km,\n' +
              f'    Rayleigh number RaIII = {Planet.RaConvectIII:.3e}.')

    # Check for whole-lid conduction
    if(zbIII_m <= Planet.eLidIII_m + Planet.deltaTBLIII_m):
        log.info(f'Underplate ice III thickness ({zbIII_m/1e3:.1f} km) is less than that of the thermal ' +
                  'boundary layers--convection is absent. Applying whole-layer conductive profile.')
        Planet.eLidIII_m = zbIII_m
        Planet.DconvIII_m = 0.0
        Planet.deltaTBLIII_m = 0.0

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
        log.debug('Modeling ice III conduction in stagnant lid...')

        thisMAbove_kg = np.sum(Planet.MLayer_kg[:Planet.Steps.nIbottom])
        # Reassign conductive profile with new bottom temperature for conductive layer
        PlidRatios = (Planet.P_MPa[Planet.Steps.nIbottom:iConductEnd+1] - Planet.P_MPa[Planet.Steps.nIbottom]) / \
                     (PconvTopIII_MPa - Planet.P_MPa[Planet.Steps.nIbottom])
        Planet.T_K[Planet.Steps.nIbottom:iConductEnd+1] = Planet.TconvIII_K**(PlidRatios) * Planet.T_K[Planet.Steps.nIbottom]**(1 - PlidRatios)

        # Get physical properties of upper conducting layer, and include 1 layer of convective layer for next step
        Planet = EvalLayerProperties(Planet, Params, Planet.Steps.nIbottom, iConductEnd+1,
                                     Planet.Ocean.surfIceEOS['III'],
                                     Planet.P_MPa[Planet.Steps.nIbottom:iConductEnd+1],
                                     Planet.T_K[Planet.Steps.nIbottom:iConductEnd+1])

        Planet = PorosityCorrectionVacIce(Planet, Params, Planet.Steps.nIbottom, iConductEnd+1,
                                          Planet.Ocean.surfIceEOS['III'],
                                          Planet.P_MPa[Planet.Steps.nIbottom:iConductEnd+1],
                                          Planet.T_K[Planet.Steps.nIbottom:iConductEnd+1])

        Planet = PropagateConduction(Planet, Params, Planet.Steps.nIbottom, iConductEnd-1)

        log.debug('Stagnant lid conductive profile complete. Modeling ice III convecting layer...')

        Planet = PropagateAdiabaticPorousVacIce(Planet, Params, iConductEnd, iConvectEnd,
                                                Planet.Ocean.surfIceEOS['III'])

        log.debug('Convective profile complete. Modeling conduction in lower thermal boundary layer...')

        if(Planet.T_K[iConvectEnd-1] > Planet.Bulk.TbIII_K):
            raise ValueError(f'Ice III bottom temperature of {Planet.Bulk.TbIII_K:.3f} K ' +
                              'is less than the temperature at the lower TBL transition of ' +
                             f'{Planet.T_K[iConvectEnd-1]:.3f} K. Try increasing Bulk.TbIII_K ' +
                              'or decreasing Bulk.Tb_K to create a more realistic thermal profile.')

        # Reassign conductive profile with new top temperature for conductive layer

        PTBLratios = (Planet.P_MPa[indsTBL] - Planet.P_MPa[iConvectEnd-1]) / (Planet.PbIII_MPa - Planet.P_MPa[iConvectEnd-1])
        Planet.T_K[indsTBL] = Planet.Bulk.TbIII_K**(PTBLratios) \
                              * Planet.T_K[iConvectEnd-1]**(1 - PTBLratios)

        # Get physical properties of thermal boundary layer
        Planet = EvalLayerProperties(Planet, Params, indsTBL[0], indsTBL[-1],
                                     Planet.Ocean.surfIceEOS['III'],
                                     Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]])

        Planet = PorosityCorrectionVacIce(Planet, Params, indsTBL[0], indsTBL[-1]+1,
                                          Planet.Ocean.surfIceEOS['III'],
                                          Planet.P_MPa[indsTBL], Planet.T_K[indsTBL])

        Planet = PropagateConduction(Planet, Params, indsTBL[0]-1, indsTBL[-1])

    log.debug('Ice III convection calculations complete.')

    return Planet


def IceVConvectSolid(Planet, Params):
    """ Apply convection models from literature to determine thermal profile and
        state variables for possibly convecting ice layers

        Assigns Planet attributes:
            TconvV_K, etaConvV_Pas, eLidV_m, deltaTBLV_m, QfromMantle_W, all physical layer arrays
    """

    log.debug('Applying solid-state convection to surface ice V based on Deschamps and Sotin (2001).')
    zbV_m = Planet.z_m[Planet.Steps.nSurfIce-1] - Planet.z_m[Planet.Steps.nIIIbottom]
    # Get "middle" pressure
    PmidV_MPa = (Planet.PbV_MPa + Planet.PbIII_MPa) / 2

    # Run calculations to get convection layer parameters
    Planet.TconvV_K, Planet.etaConvV_Pas, Planet.eLidV_m, Planet.DconvV_m, Planet.deltaTBLV_m, Planet.Ocean.QfromMantle_W, \
        Planet.RaConvectV, Planet.RaCritV = ConvectionDeschampsSotin2001(Planet.Bulk.TbIII_K, Planet.r_m[Planet.Steps.nIIIbottom],
                                     Planet.kTherm_WmK[Planet.Steps.nIIIbottom], Planet.Bulk.TbV_K, zbV_m,
                                     Planet.g_ms2[Planet.Steps.nIIIbottom], PmidV_MPa, Planet.Ocean.EOS,
                                     Planet.Ocean.surfIceEOS['V'], 5, Planet.Do.EQUIL_Q)

    log.debug(f'Ice V convection parameters:\n    T_convectV = {Planet.TconvV_K:.3f} K,\n' +
              f'    Viscosity etaConvectV = {Planet.etaConvV_Pas:.3e} Pa*s,\n' +
              f'    Conductive lid thickness eLidV = {Planet.eLidV_m/1e3:.1f} km,\n' +
              f'    Convecting layer thickness DconvV = {Planet.DconvV_m/1e3:.1f} km,\n' +
              f'    Lower TBL thickness deltaTBLV = {Planet.deltaTBLV_m/1e3:.1f} km,\n' +
              f'    Rayleigh number RaV = {Planet.RaConvectV:.3e}.')

    # Check for whole-lid conduction
    if(zbV_m <= Planet.eLidV_m + Planet.deltaTBLV_m):
        log.info(f'Underplate ice V thickness ({zbV_m/1e3:.1f} km) is less than that of the thermal ' +
                  'boundary layers--convection is absent. Applying whole-layer conductive profile.')
        Planet.eLidV_m = zbV_m
        Planet.DconvV_m = 0.0
        Planet.deltaTBLV_m = 0.0

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
        log.debug('Modeling ice V conduction in stagnant lid...')

        thisMAbove_kg = np.sum(Planet.MLayer_kg[:Planet.Steps.nIIIbottom])
        # Reassign conductive profile with new bottom temperature for conductive layer
        PlidRatios = (Planet.P_MPa[Planet.Steps.nIIIbottom:iConductEnd+1] - Planet.P_MPa[Planet.Steps.nIIIbottom]) / \
                     (PconvTopV_MPa - Planet.P_MPa[Planet.Steps.nIIIbottom])
        Planet.T_K[Planet.Steps.nIIIbottom:iConductEnd+1] = Planet.TconvV_K**(PlidRatios) \
                                                            * Planet.T_K[Planet.Steps.nIIIbottom]**(1 - PlidRatios)

        # Get physical properties of upper conducting layer, and include 1 layer of convective layer for next step
        Planet = EvalLayerProperties(Planet, Params, Planet.Steps.nIIIbottom, iConductEnd+1,
                                     Planet.Ocean.surfIceEOS['V'],
                                     Planet.P_MPa[Planet.Steps.nIIIbottom:iConductEnd+1],
                                     Planet.T_K[Planet.Steps.nIIIbottom:iConductEnd+1])

        Planet = PropagateConduction(Planet, Params, Planet.Steps.nIIIbottom, iConductEnd-1)
        log.debug('Stagnant lid conductive profile complete. Modeling ice V convecting layer...')

        Planet = PropagateAdiabaticSolid(Planet, Params, iConductEnd, iConvectEnd,
                                         Planet.Ocean.surfIceEOS['V'])

        log.debug('Convective profile complete. Modeling conduction in lower thermal boundary layer...')

        if(Planet.T_K[iConvectEnd-1] > Planet.Bulk.TbV_K):
            raise ValueError(f'Ice V bottom temperature of {Planet.Bulk.TbV_K:.3f} K ' +
                              'is less than the temperature at the lower TBL transition of ' +
                             f'{Planet.T_K[iConvectEnd-1]:.3f} K. Try increasing TbV_K ' +
                              'to create a more realistic thermal profile.')

        # Reassign conductive profile with new top temperature for conductive layer
        PTBLratios = (Planet.P_MPa[indsTBL] - Planet.P_MPa[iConvectEnd-1]) / (Planet.PbV_MPa - Planet.P_MPa[iConvectEnd-1])
        Planet.T_K[indsTBL] = Planet.Bulk.TbV_K**(PTBLratios) * Planet.T_K[iConvectEnd-1]**(1 - PTBLratios)

        # Get physical properties of thermal boundary layer
        Planet = EvalLayerProperties(Planet, Params, indsTBL[0], indsTBL[-1]+1,
                                    Planet.Ocean.surfIceEOS['V'],
                                    Planet.P_MPa[indsTBL], Planet.T_K[indsTBL])

        Planet = PropagateConduction(Planet, Params, indsTBL[0]-1, indsTBL[-1])

    log.debug('Ice V convection calculations complete.')

    return Planet


def IceVConvectPorous(Planet, Params):
    """ Apply convection models from literature to determine thermal profile and
        state variables for possibly convecting ice layers

        Assigns Planet attributes:
            TconvV_K, etaConvV_Pas, eLidV_m, deltaTBLV_m, QfromMantle_W, all physical layer arrays
    """

    log.debug('Applying solid-state convection to surface ice V based on Deschamps and Sotin (2001).')
    zbV_m = Planet.z_m[Planet.Steps.nSurfIce-1] - Planet.z_m[Planet.Steps.nIIIbottom]
    # Get "middle" pressure
    PmidV_MPa = (Planet.PbV_MPa + Planet.PbIII_MPa) / 2

    # Run calculations to get convection layer parameters
    Planet.TconvV_K, Planet.etaConvV_Pas, Planet.eLidV_m, Planet.DconvV_m, Planet.deltaTBLV_m, Planet.Ocean.QfromMantle_W, \
        Planet.RaConvectV, Planet.RaCritV = ConvectionDeschampsSotin2001(Planet.Bulk.TbIII_K, Planet.r_m[Planet.Steps.nIIIbottom],
                                     Planet.kTherm_WmK[Planet.Steps.nIIIbottom], Planet.Bulk.TbV_K, zbV_m,
                                     Planet.g_ms2[Planet.Steps.nIIIbottom], PmidV_MPa, Planet.Ocean.EOS,
                                     Planet.Ocean.surfIceEOS['V'], 5, Planet.Do.EQUIL_Q)

    log.debug(f'Ice V convection parameters:\n    T_convectV = {Planet.TconvV_K:.3f} K,\n' +
              f'    Viscosity etaConvectV = {Planet.etaConvV_Pas:.3e} Pa*s,\n' +
              f'    Conductive lid thickness eLidV = {Planet.eLidV_m/1e3:.1f} km,\n' +
              f'    Convecting layer thickness DconvV = {Planet.DconvV_m/1e3:.1f} km,\n' +
              f'    Lower TBL thickness deltaTBLV = {Planet.deltaTBLV_m/1e3:.1f} km,\n' +
              f'    Rayleigh number RaV = {Planet.RaConvectV:.3e}.')

    # Check for whole-lid conduction
    if(zbV_m <= Planet.eLidV_m + Planet.deltaTBLV_m):
        log.info(f'Underplate ice V thickness ({zbV_m/1e3:.1f} km) is less than that of the thermal ' +
                  'boundary layers--convection is absent. Applying whole-layer conductive profile.')
        Planet.eLidV_m = zbV_m
        Planet.DconvV_m = 0.0
        Planet.deltaTBLV_m = 0.0

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
        log.debug('Modeling ice V conduction in stagnant lid...')

        thisMAbove_kg = np.sum(Planet.MLayer_kg[:Planet.Steps.nIIIbottom])
        # Reassign conductive profile with new bottom temperature for conductive layer
        PlidRatios = (Planet.P_MPa[Planet.Steps.nIIIbottom:iConductEnd+1] - Planet.P_MPa[Planet.Steps.nIIIbottom]) / \
                     (PconvTopV_MPa - Planet.P_MPa[Planet.Steps.nIIIbottom])
        Planet.T_K[Planet.Steps.nIIIbottom:iConductEnd+1] = Planet.TconvV_K**(PlidRatios) \
                                                            * Planet.T_K[Planet.Steps.nIIIbottom]**(1 - PlidRatios)

        # Get physical properties of upper conducting layer, and include 1 layer of convective layer for next step
        Planet = EvalLayerProperties(Planet, Params, Planet.Steps.nIIIbottom, iConductEnd+1,
                                     Planet.Ocean.surfIceEOS['V'],
                                     Planet.P_MPa[Planet.Steps.nIIIbottom:iConductEnd+1],
                                     Planet.T_K[Planet.Steps.nIIIbottom:iConductEnd+1])

        Planet = PorosityCorrectionVacIce(Planet, Params, Planet.Steps.nIIIbottom, iConductEnd+1,
                                          Planet.Ocean.surfIceEOS['V'],
                                          Planet.P_MPa[Planet.Steps.nIIIbottom:iConductEnd+1],
                                          Planet.T_K[Planet.Steps.nIIIbottom:iConductEnd+1])

        Planet = PropagateConduction(Planet, Params, Planet.Steps.nIIIbottom, iConductEnd-1)
        log.debug('Stagnant lid conductive profile complete. Modeling ice V convecting layer...')

        Planet = PropagateAdiabaticPorousVacIce(Planet, Params, iConductEnd, iConvectEnd,
                                                Planet.Ocean.surfIceEOS['V'])

        log.debug('Convective profile complete. Modeling conduction in lower thermal boundary layer...')

        if(Planet.T_K[iConvectEnd-1] > Planet.Bulk.TbV_K):
            raise ValueError(f'Ice V bottom temperature of {Planet.Bulk.TbV_K:.3f} K ' +
                              'is less than the temperature at the lower TBL transition of ' +
                             f'{Planet.T_K[iConvectEnd-1]:.3f} K. Try increasing TbV_K ' +
                              'to create a more realistic thermal profile.')

        # Reassign conductive profile with new top temperature for conductive layer
        PTBLratios = (Planet.P_MPa[indsTBL] - Planet.P_MPa[iConvectEnd-1]) / (Planet.PbV_MPa - Planet.P_MPa[iConvectEnd-1])
        Planet.T_K[indsTBL] = Planet.Bulk.TbV_K**(PTBLratios) * Planet.T_K[iConvectEnd-1]**(1 - PTBLratios)

        # Get physical properties of thermal boundary layer
        Planet = EvalLayerProperties(Planet, Params, indsTBL[0], indsTBL[-1],
                                    Planet.Ocean.surfIceEOS['V'],
                                    Planet.P_MPa[indsTBL[:-1]], Planet.T_K[indsTBL[:-1]])

        Planet = PorosityCorrectionVacIce(Planet, Params, indsTBL[0], indsTBL[-1]+1,
                                          Planet.Ocean.surfIceEOS['V'],
                                          Planet.P_MPa[indsTBL], Planet.T_K[indsTBL])

        Planet = PropagateConduction(Planet, Params, indsTBL[0]-1, indsTBL[-1])

    log.debug('Ice V convection calculations complete.')

    return Planet


def ClathShellConvectSolid(Planet, Params):
    """ Apply convection models from literature to determine thermal profile and
        state variables for possibly convecting clathrate layers when the whole
        shell is made of clathrates.

        Assigns Planet attributes:
            Tconv_K, etaConv_Pas, eLid_m, deltaTBL_m, QfromMantle_W, all physical layer arrays
    """

    log.debug('Applying solid-state convection to surface clathrates based on Deschamps and Sotin (2001).')
    zbI_m = Planet.z_m[Planet.Steps.nIbottom]
    # Get "middle" pressure
    Pmid_MPa = (Planet.PbI_MPa - Planet.P_MPa[0]) / 2
    Planet.kTherm_WmK[0] = Planet.Ocean.surfIceEOS['Clath'].fn_kTherm_WmK(Planet.P_MPa[0], Planet.Bulk.Tb_K)

    # Run calculations to get convection layer parameters
    Planet.Tconv_K, Planet.etaConv_Pas, Planet.eLid_m, Planet.Dconv_m, Planet.deltaTBL_m, Planet.Ocean.QfromMantle_W, \
        Planet.RaConvect, Planet.RaCrit = \
        ConvectionDeschampsSotin2001(Planet.T_K[0], Planet.r_m[0], Planet.kTherm_WmK[0], Planet.Bulk.Tb_K,
                                     zbI_m, Planet.g_ms2[0], Pmid_MPa, Planet.Ocean.surfIceEOS['Clath'],
                                     Planet.Ocean.surfIceEOS['Clath'], Constants.phaseClath, Planet.Do.EQUIL_Q)

    log.debug(f'Clathrate shell convection parameters:\n    T_convect = {Planet.Tconv_K:.3f} K,\n' +
              f'    Viscosity etaConvect = {Planet.etaConv_Pas:.3e} Pa*s,\n' +
              f'    Conductive lid thickness eLid = {Planet.eLid_m/1e3:.1f} km,\n' +
              f'    Convecting layer thickness Dconv = {Planet.Dconv_m/1e3:.1f} km,\n' +
              f'    Lower TBL thickness deltaTBL = {Planet.deltaTBL_m/1e3:.1f} km,\n' +
              f'    Rayleigh number Ra = {Planet.RaConvect:.3e}.')

    # Check for whole-lid conduction
    if(zbI_m <= Planet.eLid_m + Planet.deltaTBL_m):
        log.info(f'Ice shell thickness ({zbI_m/1e3:.1f} km) is less than that of the thermal ' +
                  'boundary layers--convection is absent. Applying whole-shell conductive profile.')
        Planet.eLid_m = zbI_m
        Planet.Dconv_m = 0.0
        Planet.deltaTBL_m = 0.0

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
        log.debug('Modeling clathrate conduction in stagnant lid...')
        # Reassign conductive profile with new bottom temperature for conductive layer
        PlidRatios = (Planet.P_MPa[:nConduct+1] - Planet.P_MPa[0]) / (PconvTop_MPa - Planet.P_MPa[0])
        Planet.T_K[:nConduct+1] = Planet.Tconv_K**(PlidRatios) * Planet.T_K[0]**(1 - PlidRatios)

        # Get physical properties of upper conducting layer, and include 1 layer of convective layer for next step
        Planet = EvalLayerProperties(Planet, Params, 0, nConduct+1,
                                     Planet.Ocean.surfIceEOS['Clath'],
                                     Planet.P_MPa[:nConduct+1], Planet.T_K[:nConduct+1])

        Planet = PropagateConduction(Planet, Params, 0, nConduct-1)

        log.debug('Stagnant lid conductive profile complete. Modeling ice I convecting layer...')

        Planet = PropagateAdiabaticSolid(Planet, Params, nConduct, nConduct + nConvect,
                                         Planet.Ocean.surfIceEOS['Clath'])

        log.debug('Convective profile complete. Modeling conduction in lower thermal boundary layer...')

        # Reassign conductive profile with new top temperature for conductive layer
        PTBLratios = (Planet.P_MPa[indsTBL] - Planet.P_MPa[nConduct+nConvect-1]) / (Planet.PbI_MPa - Planet.P_MPa[nConduct+nConvect-1])
        Planet.T_K[indsTBL] = Planet.Bulk.Tb_K**(PTBLratios) * Planet.T_K[nConduct+nConvect-1]**(1 - PTBLratios)

        # Get physical properties of thermal boundary layer
        Planet = EvalLayerProperties(Planet, Params, indsTBL[0], indsTBL[-1]+1,
                                    Planet.Ocean.surfIceEOS['Clath'],
                                    Planet.P_MPa[indsTBL], Planet.T_K[indsTBL])

        Planet = PropagateConduction(Planet, Params, indsTBL[0]-1, indsTBL[-1])

    Planet.zClath_m = Planet.z_m[Planet.Steps.nIbottom]

    log.debug('Clathrate convection calculations complete.')

    return Planet


def ClathShellConvectPorous(Planet, Params):
    """ Apply convection models from literature to determine thermal profile and
        state variables for possibly convecting clathrate layers when the whole
        shell is made of clathrates.

        Assigns Planet attributes:
            Tconv_K, etaConv_Pas, eLid_m, deltaTBL_m, QfromMantle_W, all physical layer arrays
    """

    log.debug('Applying solid-state convection to surface clathrates based on Deschamps and Sotin (2001).')
    zbI_m = Planet.z_m[Planet.Steps.nIbottom]
    # Get "middle" pressure
    Pmid_MPa = (Planet.PbI_MPa - Planet.P_MPa[0]) / 2
    Planet.kTherm_WmK[0] = Planet.Ocean.surfIceEOS['Clath'].fn_kTherm_WmK(Planet.P_MPa[0], Planet.Bulk.Tb_K)

    # Run calculations to get convection layer parameters
    Planet.Tconv_K, Planet.etaConv_Pas, Planet.eLid_m, Planet.Dconv_m, Planet.deltaTBL_m, Planet.Ocean.QfromMantle_W, \
        Planet.RaConvect, Planet.RaCrit = \
        ConvectionDeschampsSotin2001(Planet.T_K[0], Planet.r_m[0], Planet.kTherm_WmK[0], Planet.Bulk.Tb_K,
                                     zbI_m, Planet.g_ms2[0], Pmid_MPa, Planet.Ocean.surfIceEOS['Clath'],
                                     Planet.Ocean.surfIceEOS['Clath'], Constants.phaseClath, Planet.Do.EQUIL_Q)

    log.debug(f'Clathrate shell convection parameters:\n    T_convect = {Planet.Tconv_K:.3f} K,\n' +
              f'    Viscosity etaConvect = {Planet.etaConv_Pas:.3e} Pa*s,\n' +
              f'    Conductive lid thickness eLid = {Planet.eLid_m/1e3:.1f} km,\n' +
              f'    Convecting layer thickness Dconv = {Planet.Dconv_m/1e3:.1f} km,\n' +
              f'    Lower TBL thickness deltaTBL = {Planet.deltaTBL_m/1e3:.1f} km,\n' +
              f'    Rayleigh number Ra = {Planet.RaConvect:.3e}.')

    # Check for whole-lid conduction
    if(zbI_m <= Planet.eLid_m + Planet.deltaTBL_m):
        log.info(f'Ice shell thickness ({zbI_m/1e3:.1f} km) is less than that of the thermal ' +
                  'boundary layers--convection is absent. Applying whole-shell conductive profile.')
        Planet.eLid_m = zbI_m
        Planet.Dconv_m = 0.0
        Planet.deltaTBL_m = 0.0

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
        log.debug('Modeling clathrate conduction in stagnant lid...')
        # Reassign conductive profile with new bottom temperature for conductive layer
        PlidRatios = (Planet.P_MPa[:nConduct+1] - Planet.P_MPa[0]) / (PconvTop_MPa - Planet.P_MPa[0])
        Planet.T_K[:nConduct+1] = Planet.Tconv_K**(PlidRatios) * Planet.T_K[0]**(1 - PlidRatios)

        # Get physical properties of upper conducting layer, and include 1 layer of convective layer for next step
        Planet = EvalLayerProperties(Planet, Params, 0, nConduct+1,
                                     Planet.Ocean.surfIceEOS['Clath'],
                                     Planet.P_MPa[:nConduct+1], Planet.T_K[:nConduct+1])

        Planet = PorosityCorrectionVacIce(Planet, Params, 0, nConduct+1,
                                          Planet.Ocean.surfIceEOS['Clath'],
                                          Planet.P_MPa[:nConduct+1], Planet.T_K[:nConduct+1])

        Planet = PropagateConduction(Planet, Params, 0, nConduct-1)

        log.debug('Stagnant lid conductive profile complete. Modeling ice I convecting layer...')

        Planet = PropagateAdiabaticPorousVacIce(Planet, Params, nConduct, nConduct + nConvect,
                                                Planet.Ocean.surfIceEOS['Clath'])

        log.debug('Convective profile complete. Modeling conduction in lower thermal boundary layer...')

        # Reassign conductive profile with new top temperature for conductive layer
        PTBLratios = (Planet.P_MPa[indsTBL] - Planet.P_MPa[nConduct+nConvect-1]) / (Planet.PbI_MPa - Planet.P_MPa[nConduct+nConvect-1])
        Planet.T_K[indsTBL] = Planet.Bulk.Tb_K**(PTBLratios) * Planet.T_K[nConduct+nConvect-1]**(1 - PTBLratios)

        # Get physical properties of thermal boundary layer
        Planet = EvalLayerProperties(Planet, Params, indsTBL[0], indsTBL[-1]+1,
                                    Planet.Ocean.surfIceEOS['Clath'],
                                    Planet.P_MPa[indsTBL], Planet.T_K[indsTBL])

        Planet = PorosityCorrectionVacIce(Planet, Params, indsTBL[0], indsTBL[-1]+1,
                                          Planet.Ocean.surfIceEOS['Clath'],
                                          Planet.P_MPa[indsTBL], Planet.T_K[indsTBL])

        Planet = PropagateConduction(Planet, Params, indsTBL[0]-1, indsTBL[-1])

    Planet.zClath_m = Planet.z_m[Planet.Steps.nIbottom]

    log.debug('Clathrate convection calculations complete.')

    return Planet
