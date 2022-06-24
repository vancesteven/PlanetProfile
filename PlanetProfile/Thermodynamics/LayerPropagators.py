import numpy as np
import logging

from PlanetProfile.Thermodynamics.IronCore import IronCoreLayers
from PlanetProfile.Thermodynamics.HydroEOS import GetPfreeze, GetTfreeze, \
    PhaseConv, GetPhaseIndices, GetIceEOS, GetOceanEOS
from PlanetProfile.Thermodynamics.InnerEOS import GetHtidalFunc, GetphiCalc
from PlanetProfile.Thermodynamics.Silicates import SilicateLayers
from PlanetProfile.Thermodynamics.ThermalProfiles.Convection import IceIConvectSolid, IceIConvectPorous, \
    IceIIIConvectSolid, IceIIIConvectPorous, IceVConvectSolid, IceVConvectPorous, \
    ClathShellConvectSolid, ClathShellConvectPorous
from PlanetProfile.Thermodynamics.ThermalProfiles.IceConduction import IceIWholeConductSolid, IceIWholeConductPorous, \
    IceIConductClathLidSolid, IceIConductClathLidPorous, IceIConductClathUnderplateSolid, IceIConductClathUnderplatePorous, \
    IceIIIConductSolid, IceIIIConductPorous, IceVConductSolid, IceVConductPorous
from PlanetProfile.Thermodynamics.ThermalProfiles.ThermalProfiles import ConductiveTemperature, ConvectionDeschampsSotin2001
from PlanetProfile.Utilities.defineStructs import Constants

# Assign logger
log = logging.getLogger('PlanetProfile')

def IceLayers(Planet, Params):
    """ Layer propagation from surface downward through the ice using geophysics.
        Iteratively sets up the thermal profile (the density and temperature)
        of the layer with each pressure step for all ice layers including
        ice Ih, ice III, ice V by calling the necessary subfunctions

        Assigns Planet attributes:
            Steps.nSurfIce, phase, r_m, z_m, g_ms2, T_K, P_MPa, rho_kgm3, Cp_JkgK, alpha_pK, MLayer_kg, PbI_MPa, Pb_MPa
    """
    Planet.Steps.nIbottom = Planet.Steps.nClath + Planet.Steps.nIceI
    Planet.Steps.nIIIbottom = Planet.Steps.nIbottom + Planet.Steps.nIceIIILitho
    Planet.Steps.nSurfIce = Planet.Steps.nIIIbottom + Planet.Steps.nIceVLitho
    # Assign phase values for near-surface ices
    Planet.phase[:Planet.Steps.nIbottom] = 1  # Ice Ih layers (some will be reassigned if Do.CLATHRATE = True)
    Planet.phase[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] = 3  # Ice III layers
    Planet.phase[Planet.Steps.nIIIbottom:Planet.Steps.nSurfIce] = 5  # Ice V layers

    # Get the pressure consistent with the bottom of the surface ice layer that is
    # consistent with the choice of Tb_K we suppose for this model
    Planet.PbI_MPa = GetPfreeze(Planet.Ocean.meltEOS, 1, Planet.Bulk.Tb_K,
                                PLower_MPa=Planet.PfreezeLower_MPa, PUpper_MPa=Planet.PfreezeUpper_MPa,
                                PRes_MPa=Planet.PfreezeRes_MPa, UNDERPLATE=(Planet.Do.BOTTOM_ICEIII or Planet.Do.BOTTOM_ICEV),
                                ALLOW_BROKEN_MODELS=Params.ALLOW_BROKEN_MODELS, DO_EXPLOREOGRAM=Params.DO_EXPLOREOGRAM)
    if(Planet.Do.CLATHRATE and
            (Planet.Bulk.clathType == 'bottom' or
             Planet.Bulk.clathType == 'whole')):
        PbClath_MPa = Planet.Ocean.ClathDissoc.PbClath_MPa()
        if not np.isnan(PbClath_MPa):
            log.debug(f'Clathrate dissociation pressure: {PbClath_MPa:.3f} MPa.')
            if PbClath_MPa < Planet.PbI_MPa:
                raise ValueError('Dissociation pressure for clathrates is lower than the ice Ih ' +
                                 'melting pressure consistent with Bulk.Tb_K. This means ice Ih layers ' +
                                 'will be found underneath the clathrate layers, inconsistent with the ' +
                                 'assumption that clathrates are in contact with the ocean. Increase ' +
                                 f'Bulk.Tb_K until PbClath_MPa ({PbClath_MPa:.3f} MPa) ' +
                                 f'exceeds PbI_MPa ({Planet.PbI_MPa:.3f} MPa).')
        Planet.PbI_MPa = PbClath_MPa
    else:
        if np.isnan(Planet.PbI_MPa):
            msg = f'No valid phase transition was found for Tb_K = {Planet.Bulk.Tb_K:.3f} K for P in the range ' + \
                  f'[{Planet.PfreezeLower_MPa:.1f} MPa, {Planet.PfreezeUpper_MPa:.1f} MPa]. ' + \
                  'This likely means Tb_K is too high and the phase at the lower end of this range matches ' + \
                  'the phase at the upper end. Try decreasing Tb_K. The ice shell will be set to zero thickness.'
            if not Params.DO_EXPLOREOGRAM:
                if Planet.Bulk.Tb_K > 271:
                    log.warning(msg)
                else:
                    raise ValueError(msg)
            Planet.PbI_MPa = 0.0
        log.debug(f'Ice Ih transition pressure: {Planet.PbI_MPa:.3f} MPa.')

    if Planet.PbI_MPa > 0:
        # Now do the same for HP ices, if present, to make sure we have a possible configuration before continuing
        if Planet.Do.BOTTOM_ICEV:
            Planet.PbIII_MPa = GetPfreeze(Planet.Ocean.meltEOS, 3, Planet.Bulk.TbIII_K,
                       PLower_MPa=Planet.PbI_MPa, PUpper_MPa=Planet.Ocean.PHydroMax_MPa,
                       PRes_MPa=Planet.PfreezeRes_MPa, UNDERPLATE=True,
                       ALLOW_BROKEN_MODELS=Params.ALLOW_BROKEN_MODELS, DO_EXPLOREOGRAM=Params.DO_EXPLOREOGRAM)
            if(Planet.PbIII_MPa <= Planet.PbI_MPa) or np.isnan(Planet.PbIII_MPa):
                msg = 'Ice III bottom pressure is not greater than ice I bottom pressure. ' + \
                      'This likely indicates TbIII_K is too high for the corresponding Tb_K.' + \
                      f'\nPbI_MPa = {Planet.PbI_MPa:.3f}' + \
                      f', Tb_K = {Planet.Bulk.Tb_K:.3f}' + \
                      f'\nPbIII_MPa = {Planet.PbIII_MPa:.3f}' + \
                      f', TbIII_K = {Planet.Bulk.TbIII_K:.3f}'
                if Params.ALLOW_BROKEN_MODELS:
                    Planet.PbIII_MPa = np.nan
                    if Params.DO_EXPLOREOGRAM:
                        log.info(msg)
                    else:
                        log.error(msg)
                    Planet.Do.VALID = False
                else:
                    raise ValueError(msg)
            if not np.isnan(Planet.PbIII_MPa):
                Planet.PbV_MPa = GetPfreeze(Planet.Ocean.meltEOS, 5, Planet.Bulk.TbV_K,
                                              PLower_MPa=Planet.PbIII_MPa, PUpper_MPa=Planet.Ocean.PHydroMax_MPa,
                                              PRes_MPa=Planet.PfreezeRes_MPa, UNDERPLATE=False,
                                              ALLOW_BROKEN_MODELS=Params.ALLOW_BROKEN_MODELS,
                                              DO_EXPLOREOGRAM=Params.DO_EXPLOREOGRAM)
            else:
                Planet.PbV_MPa = np.nan
            Planet.Pb_MPa = Planet.PbV_MPa
            if(Planet.PbV_MPa <= Planet.PbIII_MPa) or np.isnan(Planet.PbV_MPa):
                msg = 'Ice V bottom pressure is not greater than ice III bottom pressure. ' + \
                      'This likely indicates TbV_K is too high for the corresponding TbIII_K.' + \
                      f'\nPbIII_MPa = {Planet.PbIII_MPa:.3f}' + \
                      f', TbIII_K = {Planet.Bulk.TbIII_K:.3f}' + \
                      f'\nPbV_MPa = {Planet.PbV_MPa:.3f}' + \
                      f', TbV_K = {Planet.Bulk.TbV_K:.3f}'
                if Params.ALLOW_BROKEN_MODELS:
                    Planet.PbIII_MPa = np.nan
                    if Params.DO_EXPLOREOGRAM:
                        log.info(msg)
                    else:
                        log.error(msg)
                    Planet.Do.VALID = False
                else:
                    raise ValueError(msg)
        elif Planet.Do.BOTTOM_ICEIII:
            Planet.PbIII_MPa = GetPfreeze(Planet.Ocean.meltEOS, 3, Planet.Bulk.TbIII_K,
                                          PLower_MPa=Planet.PbI_MPa, PUpper_MPa=Planet.Ocean.PHydroMax_MPa,
                                          PRes_MPa=Planet.PfreezeRes_MPa, UNDERPLATE=False,
                                          ALLOW_BROKEN_MODELS=Params.ALLOW_BROKEN_MODELS,
                                          DO_EXPLOREOGRAM=Params.DO_EXPLOREOGRAM)
            if(Planet.PbIII_MPa <= Planet.PbI_MPa) or np.isnan(Planet.PbIII_MPa):
                msg = 'Ice III bottom pressure is not greater than ice I bottom pressure. ' + \
                      'This likely indicates TbIII_K is too high for the corresponding Tb_K.' + \
                      f'\nPbI_MPa = {Planet.PbI_MPa:.3f}' + \
                      f', Tb_K = {Planet.Bulk.Tb_K:.3f}' + \
                      f'\nPbIII_MPa = {Planet.PbIII_MPa:.3f}' + \
                      f', TbIII_K = {Planet.Bulk.TbIII_K:.3f}'
                if Params.ALLOW_BROKEN_MODELS:
                    Planet.PbIII_MPa = np.nan
                    if Params.DO_EXPLOREOGRAM:
                        log.info(msg)
                    else:
                        log.error(msg)
                    Planet.Do.VALID = False
                else:
                    raise ValueError(msg)
            Planet.Pb_MPa = Planet.PbIII_MPa
        else:
            Planet.Pb_MPa = Planet.PbI_MPa

    elif Planet.Pb_MPa == 0 and Planet.Bulk.Tsurf_K == Planet.Bulk.Tb_K:
        # This config needs to be caught in SetupInit.
        pass

    else:
        Planet.Pb_MPa = np.nan
        Planet.Do.VALID = False
        if not Params.ALLOW_BROKEN_MODELS:
            raise RuntimeError('Unable to find a valid pressure corresponding to Bulk.TbX_K values. ' +
                               f'This is usually because Bulk.Tb_K (currently {Planet.Bulk.Tb_K:.3f}) ' +
                               'is set too high. Try decreasing Bulk.Tb_K before running again.')

    # Now, we want to check for a convective profile. First, we need to get zb_km, so we need to suppose
    # a whole-layer conductive profile. The densities will change slightly, so we depart from self-consistency
    # here. Repeated applications of IceConvect will get closer to self-consistency.

    if Planet.Pb_MPa > 0:
        if Planet.Do.CLATHRATE:
            """ For ice shells insulated by a layer of clathrate at the surface or against the bottom
                Calculates state variables of the layer with each pressure step
                Applies different behavior based on Bulk.clathType:
                    'top': Models a conductive lid of clathrates limited to Bulk.clathMaxThick_m or eLid_m
                        (calculated for convection), whichever is less
                    'bottom': Models a clathrate layer at the ice-ocean interface with a fixed thickness
                        equal to Bulk.clathMaxThick_m. Assumes a purely conductive lid, as justified in
                        Kamata et al. (2019) for Pluto: https://doi.org/10.1038/s41561-019-0369-8
                    'whole': Models clathrates as present throughout the outer ice shell, checking for
                        convection, and assumes no ice I is present in the shell. This option is handled in IceLayers.
            """
            if Planet.Bulk.clathType == 'top':
                log.debug('Applying clathrate lid conduction.')
                Planet.phase[:Planet.Steps.nClath] = Constants.phaseClath
                if Planet.Do.POROUS_ICE:
                    Planet = IceIConductClathLidPorous(Planet, Params)
                else:
                    Planet = IceIConductClathLidSolid(Planet, Params)
            elif Planet.Bulk.clathType == 'bottom':
                log.debug('Applying clathrate underplating to ice I shell.')
                Planet.phase[Planet.Steps.nIceI:Planet.Steps.nIbottom] = Constants.phaseClath
                if Planet.Do.POROUS_ICE:
                    Planet = IceIConductClathUnderplatePorous(Planet, Params)
                else:
                    Planet = IceIConductClathUnderplateSolid(Planet, Params)

            elif Planet.Bulk.clathType == 'whole':
                log.debug('Applying whole-shell clathrate modeling with possible convection.')
                Planet.phase[:Planet.Steps.nIbottom] = Constants.phaseClath
                if Planet.Do.POROUS_ICE:
                    Planet = IceIWholeConductPorous(Planet, Params)
                else:
                    Planet = IceIWholeConductSolid(Planet, Params)
            else:
                raise ValueError(f'Bulk.clathType option "{Planet.Bulk.clathType}" is not supported. ' +
                                 'Options are "top", "bottom", and "whole".')
        else:
            if Planet.Do.POROUS_ICE:
                Planet = IceIWholeConductPorous(Planet, Params)
            else:
                Planet = IceIWholeConductSolid(Planet, Params)

        log.debug('Upper ice initial conductive profile complete.')

        if not Planet.Do.NO_ICE_CONVECTION and not Planet.Bulk.clathType == 'bottom':
            # Record zb_m to see if it gets adjusted significantly
            zbOld_m = Planet.z_m[Planet.Steps.nIbottom-1] + 0.0
            # Now check for convective region and get dimensions if present
            if Planet.Do.CLATHRATE and Planet.Bulk.clathType == 'whole':
                if Planet.Do.POROUS_ICE:
                    Planet = ClathShellConvectPorous(Planet, Params)
                else:
                    Planet = ClathShellConvectSolid(Planet, Params)
            else:
                if Planet.Do.POROUS_ICE:
                    Planet = IceIConvectPorous(Planet, Params)
                else:
                    Planet = IceIConvectSolid(Planet, Params)
            # Run IceIConvect a second time if zbI_m changed by more than a set tolerance
            if(np.abs(Planet.z_m[Planet.Steps.nIbottom-1] - zbOld_m)/Planet.z_m[Planet.Steps.nIbottom-1] > Planet.Bulk.zbChangeTol_frac):
                log.debug('The bottom depth of surface ice I changed by ' +
                        f'{(Planet.z_m[Planet.Steps.nIbottom-1] - zbOld_m)/1e3:.2f} km from IceIConvect, which is greater than ' +
                        f'{Planet.Bulk.zbChangeTol_frac * 100:.0f}%. running IceIConvect a second time...')
                if Planet.Do.CLATHRATE and Planet.Bulk.clathType == 'whole':
                    if Planet.Do.POROUS_ICE:
                        Planet = ClathShellConvectPorous(Planet, Params)
                    else:
                        Planet = ClathShellConvectSolid(Planet, Params)
                else:
                    if Planet.Do.POROUS_ICE:
                        Planet = IceIConvectPorous(Planet, Params)
                    else:
                        Planet = IceIConvectSolid(Planet, Params)
        else:
            if Planet.Do.NO_ICE_CONVECTION:
                log.debug('NO_ICE_CONVECTION is True -- skipping ice I convection calculations.')
                Planet.RaConvect = np.nan
                Planet.RaCrit = np.nan
                Planet.Tconv_K = np.nan
                Planet.etaConv_Pas = np.nan
            else:
                Pbot_MPa = np.arange(Planet.PfreezeLower_MPa, Planet.PfreezeUpper_MPa, Planet.PfreezeRes_MPa)
                Tbot_K = np.arange(Planet.T_K[Planet.Steps.nIceI], np.maximum(273, Planet.T_K[Planet.Steps.nIceI]+20), 0.5)
                iceImeltEOS = GetOceanEOS('PureH2O', 0.0, Pbot_MPa, Tbot_K, None,
                                          phaseType=Planet.Ocean.phaseType, FORCE_NEW=True)
                Planet.Tconv_K, Planet.etaConv_Pas, _, _, _, _, Planet.RaConvect, Planet.RaCrit \
                    = ConvectionDeschampsSotin2001(Planet.T_K[0], Planet.r_m[0], Planet.kTherm_WmK[0],
                                                   Planet.T_K[Planet.Steps.nIceI], Planet.z_m[Planet.Steps.nIceI],
                                                   Planet.g_ms2[0], np.mean([Planet.P_MPa[Planet.Steps.nIceI], Planet.P_MPa[0]]),
                                                   iceImeltEOS, Planet.Ocean.surfIceEOS['Ih'], 1, Planet.Do.EQUIL_Q)
            Planet.eLid_m = Planet.z_m[Planet.Steps.nSurfIce]
            Planet.Dconv_m = 0.0
            Planet.deltaTBL_m = 0.0
            # Find the surface heat flux from the conductive profile. This assumes there is no tidal heating!
            Planet.Ocean.QfromMantle_W = Planet.kTherm_WmK[Planet.Steps.nIbottom-2] * Planet.T_K[Planet.Steps.nIbottom-2] / \
                                         (Planet.z_m[Planet.Steps.nIbottom-1] - Planet.z_m[Planet.Steps.nIbottom-2]) \
                                         * np.log(Planet.T_K[Planet.Steps.nIbottom-1]/Planet.T_K[Planet.Steps.nIbottom-2]) \
                                         * 4*np.pi*(Planet.Bulk.R_m - Planet.z_m[Planet.Steps.nIbottom-1])**2

        # Additional adiabats + possible convection in ice III and/or V underplate layers --
        # for thick, cold ice shells and saline oceans
        if Planet.Do.BOTTOM_ICEV:
            log.debug('Modeling ice III and V underplating...')
            Planet = IceIIIUnderplate(Planet, Params)
            Planet = IceVUnderplate(Planet, Params)
        elif Planet.Do.BOTTOM_ICEIII:
            log.debug('Modeling ice III underplating...')
            Planet = IceIIIUnderplate(Planet, Params)

        # Print and save transition pressure and upper ice thickness
        Planet.zb_km = Planet.z_m[Planet.Steps.nSurfIce] / 1e3
        log.info(f'Upper ice transition pressure: {Planet.Pb_MPa:.3f} MPa, ' +
                 f'thickness: {Planet.zb_km:.3f} km.')

        # Set surface HP ice layers to have negative phase ID to differentiate from in-ocean HP ices
        indsHP = np.where(np.logical_and(abs(Planet.phase[:Planet.Steps.nSurfIce]) > 1,
                                         abs(Planet.phase[:Planet.Steps.nSurfIce]) <= 6))[0]
        Planet.phase[:Planet.Steps.nSurfIce][indsHP] = -Planet.phase[:Planet.Steps.nSurfIce][indsHP]

        # Get heat flux out of the possibly convecting region
        Planet.qCon_Wm2 = Planet.Ocean.QfromMantle_W / (4*np.pi * (Planet.Bulk.R_m - Planet.z_m[Planet.Steps.nSurfIce])**2)
        # Get heat flux at the surface, assuming Htidal = Qrad = 0 throughout the entire hydrosphere.
        Planet.qSurf_Wm2 = Planet.Ocean.QfromMantle_W / (4*np.pi * Planet.Bulk.R_m**2)

    elif Planet.Pb_MPa == 0 and Planet.Bulk.Tsurf_K == Planet.Bulk.Tb_K:
        Planet.zb_km = 0
        # This configuration should be accounted for in SetupInit.
    else:
        # Set necessary empty variables for when we have an invalid profile
        Planet.Do.VALID = False
        Planet.zb_km, Planet.PbClathMax_MPa, Planet.PbIII_MPa, Planet.PbV_MPa, Planet.RaConvect, \
        Planet.RaCrit, Planet.Tconv_K, Planet.TconvIII_K, Planet.TconvV_K, Planet.etaConv_Pas, \
        Planet.etaConvIII_Pas, Planet.etaConvV_Pas, Planet.eLid_m, Planet.Dconv_m, Planet.deltaTBL_m, \
        Planet.qCon_Wm2, Planet.qSurf_Wm2, Planet.TclathTrans_K, Planet.Ocean.QfromMantle_W \
            = (np.nan for _ in range(19))

    return Planet


def IceIIIUnderplate(Planet, Params):
    """ Conductive and convective profile calculations for ice III layers between
        the ocean/underplate ice V and surface ice I layer.

        Assigns Planet attributes:
            PbIII_MPa, all physical layer arrays
    """

    log.debug(f'Ice III bottom phase transition pressure: {Planet.PbIII_MPa:.3f} MPa ' +
             f'at TbIII_K = {Planet.Bulk.TbIII_K:.3f} K.')

    if Planet.Do.POROUS_ICE:
        Planet = IceIIIConductPorous(Planet, Params)
    else:
        Planet = IceIIIConductSolid(Planet, Params)

    if not Planet.Do.NO_ICE_CONVECTION:
        # Record zbIII_m to see if it gets adjusted significantly
        zbIIIold_m = Planet.z_m[Planet.Steps.nIIIbottom-1] + 0.0
        # Now check for convective region and get dimensions if present
        if Planet.Do.POROUS_ICE:
            Planet = IceIIIConvectPorous(Planet, Params)
        else:
            Planet = IceIIIConvectSolid(Planet, Params)
        # Run IceIIIConvect a second time if zbIII_m changed by more than a set tolerance
        if(np.abs(Planet.z_m[Planet.Steps.nIIIbottom-1] - zbIIIold_m)/Planet.z_m[Planet.Steps.nIIIbottom-1] > Planet.Bulk.zbChangeTol_frac):
            log.debug('The bottom depth of underplate ice III changed by ' +
                    f'{(Planet.z_m[Planet.Steps.nIIIbottom-1] - zbIIIold_m)/1e3:.2f} km from IceIIIConvect, which is greater than ' +
                    f'{Planet.Bulk.zbChangeTol_frac * 100:.0f}%. running IceIIIConvect a second time...')
            if Planet.Do.POROUS_ICE:
                Planet = IceIIIConvectPorous(Planet, Params)
            else:
                Planet = IceIIIConvectSolid(Planet, Params)
    else:
        log.debug('NO_ICE_CONVECTION is True -- skipping ice III convection calculations.')
        Planet.eLidIII_m = Planet.Planet.z_m[Planet.Steps.nIIIbottom-1]
        Planet.DconvIII_m = 0.0
        Planet.deltaTBLIII_m = 0.0

    return Planet


def IceVUnderplate(Planet, Params):
    """ Conductive and convective profile calculations for ice V layers between
        the ocean and surface ice III layer.

        Assigns Planet attributes:
            PbV_MPa, all physical layer arrays
    """
    log.debug(f'Ice V bottom phase transition pressure: {Planet.PbV_MPa:.3f} MPa ' +
                             f'at TbV_K = {Planet.Bulk.TbV_K:.3f} K.')

    if Planet.Do.POROUS_ICE:
        Planet = IceVConductPorous(Planet, Params)
    else:
        Planet = IceVConductSolid(Planet, Params)

    if not Planet.Do.NO_ICE_CONVECTION:
        # Record zbV_m to see if it gets adjusted significantly
        zbVold_m = Planet.z_m[Planet.Steps.nSurfIce-1] + 0.0
        # Now check for convective region and get dimensions if present
        if Planet.Do.POROUS_ICE:
            Planet = IceVConvectPorous(Planet, Params)
        else:
            Planet = IceVConvectSolid(Planet, Params)
        # Run IceVConvect a second time if zbV_m changed by more than a set tolerance
        if(np.abs(Planet.z_m[Planet.Steps.nSurfIce-1] - zbVold_m)/Planet.z_m[Planet.Steps.nSurfIce-1] > Planet.Bulk.zbChangeTol_frac):
            log.debug('The bottom depth of underplate ice V changed by ' +
                     f'{(Planet.z_m[Planet.Steps.nSurfIce-1] - zbVold_m)/1e3:.2f} km from IceVConvect, which is greater than ' +
                     f'{Planet.Bulk.zbChangeTol_frac * 100:.0f}%. running IceVConvect a second time...')
            if Planet.Do.POROUS_ICE:
                Planet = IceVConvectPorous(Planet, Params)
            else:
                Planet = IceVConvectSolid(Planet, Params)
    else:
        log.debug('NO_ICE_CONVECTION is True -- skipping ice V convection calculations.')
        Planet.eLidV_m = Planet.Planet.z_m[Planet.Steps.nSurfIce-1]
        Planet.DconvV_m = 0.0
        Planet.deltaTBLV_m = 0.0

    return Planet


def OceanLayers(Planet, Params):
    """ Geophysical and thermodynamic calculations for ocean layer
        Calculates state variables of the layer with each pressure step

        Assigns Planet attributes:
            phase, r_m, z_m, g_ms2, T_K, P_MPa, rho_kgm3, Cp_JkgK, alpha_pK, MLayer_kg
    """
    if Planet.Do.VALID:
        log.debug('Evaluating ocean layers.')

        # Confirm that we haven't made mistakes in phase assignment in IceLayers()
        # Note that this assignment is only temporary, as we re-check the phase of
        # this layer with the ocean EOS momentarily--it could also be ice II or III.
        if Planet.phase[Planet.Steps.nSurfIce] != 0:
            raise ValueError('Phase of first "ocean" layer is not zero.')

        # Start ocean pressure at the melting pressure and assign linear pressure profile for layers
        POcean_MPa = np.arange(Planet.Pb_MPa, Planet.Ocean.PHydroMax_MPa, Planet.Ocean.deltaP)
        Planet.Steps.nOceanMax = np.size(POcean_MPa)

        # Initialize remaining local arrays
        TOcean_K, rhoOcean_kgm3, CpOcean_JkgK, alphaOcean_pK, kThermOcean_WmK \
            = (np.zeros(Planet.Steps.nOceanMax) for _ in range(5))
        TOcean_K = np.insert(TOcean_K, 0, Planet.T_K[Planet.Steps.nSurfIce])

        # Add HP ices to iceEOS if needed
        PHydroMax_MPa = np.maximum(Planet.Ocean.PHydroMax_MPa, Planet.Sil.PHydroMax_MPa)
        if PHydroMax_MPa > Planet.Ocean.PHydroMax_MPa:
            PHPices_MPa = np.arange(POcean_MPa[0], PHydroMax_MPa, Planet.Ocean.deltaP)
        else:
            PHPices_MPa = POcean_MPa
        if PHydroMax_MPa > Constants.PminHPices_MPa:
            GetOceanHPIceEOS(Planet, Params, PHPices_MPa)

        # Do initial ocean step separately in order to catch potential Melosh layer--
        # see Melosh et al. (2004): https://doi.org/10.1016/j.icarus.2003.11.026
        log.debug(f'il: {Planet.Steps.nSurfIce:d}; P_MPa: {POcean_MPa[0]:.3f}; ' +
                  f'T_K: {TOcean_K[0]:.3f}; phase: {Planet.phase[Planet.Steps.nSurfIce]:d}')
        rhoOcean_kgm3[0] = Planet.Ocean.EOS.fn_rho_kgm3(POcean_MPa[0], TOcean_K[0])
        CpOcean_JkgK[0] = Planet.Ocean.EOS.fn_Cp_JkgK(POcean_MPa[0], TOcean_K[0])
        alphaOcean_pK[0] = Planet.Ocean.EOS.fn_alpha_pK(POcean_MPa[0], TOcean_K[0])
        kThermOcean_WmK[0] = Planet.Ocean.EOS.fn_kTherm_WmK(POcean_MPa[0], TOcean_K[0])
        if alphaOcean_pK[0] < 0:
            log.info(f'Thermal expansivity alpha at the ice-ocean interface is negative. Modeling Melosh et al. conductive layer.')
            # Layer should be thin, so we just use a fixed dT/dz value
            dTdz = Planet.Ocean.QfromMantle_W / (4*np.pi * Planet.r_m[Planet.Steps.nSurfIce]**2) / kThermOcean_WmK[0]
            i = 0
            # Use a smaller pressure step to make sure we don't overshoot by a lot
            deltaPMelosh = Planet.Ocean.deltaP / 100
            # Initialize for while loop
            alphaMelosh_pK = alphaOcean_pK[0] + 0
            rhoMelosh_kgm3 = rhoOcean_kgm3[0] + 0
            gMelosh_ms2 = Planet.g_ms2[Planet.Steps.nSurfIce] + 0
            TMelosh_K = TOcean_K[0] + 0
            deltaPtop = 0
            zMelosh = 0
            while alphaMelosh_pK < 0:
                dz = deltaPMelosh*1e6 / gMelosh_ms2 / rhoMelosh_kgm3
                zMelosh += dz
                deltaPtop += deltaPMelosh
                thisP_MPa = deltaPtop + POcean_MPa[0]
                # Model temperature as linear and conductive in this layer
                TMelosh_K += dTdz * dz

                rhoMelosh_kgm3 = Planet.Ocean.EOS.fn_rho_kgm3(thisP_MPa, TMelosh_K)
                alphaMelosh_pK = Planet.Ocean.EOS.fn_alpha_pK(thisP_MPa, TMelosh_K)

                if alphaMelosh_pK > 0 or deltaPtop >= (i+1)*Planet.Ocean.deltaP:
                    i += 1
                    POcean_MPa[i] = thisP_MPa + 0
                    TOcean_K[i] = TMelosh_K + 0
                    rhoOcean_kgm3[i] = Planet.Ocean.EOS.fn_rho_kgm3(POcean_MPa[i], TOcean_K[i])
                    CpOcean_JkgK[i] = Planet.Ocean.EOS.fn_Cp_JkgK(POcean_MPa[i], TOcean_K[i])
                    alphaOcean_pK[i] = Planet.Ocean.EOS.fn_alpha_pK(POcean_MPa[i], TOcean_K[i])
                    kThermOcean_WmK[i] = Planet.Ocean.EOS.fn_kTherm_WmK(POcean_MPa[i], TOcean_K[i])
                    Planet.phase[Planet.Steps.nSurfIce+i] = Planet.Ocean.EOS.fn_phase(POcean_MPa[i], TOcean_K[i]).astype(np.int_)
                    log.debug(f'il: {Planet.Steps.nSurfIce+i:d}; P_MPa: {POcean_MPa[i]:.3f}; ' +
                              f'T_K: {TOcean_K[i]:.3f}; phase: {Planet.phase[Planet.Steps.nSurfIce+i]:d}')
            iStart = i
            # Reset pressure profile to use standard pressure step below Melosh layer bottom
            POcean_MPa[i+1:] = np.linspace(POcean_MPa[i], POcean_MPa[-1], Planet.Steps.nOceanMax - i)[1:]
            log.info(f'Melosh et al. layer complete, thickness {zMelosh:.1f} m.')

        else:
            # Now use the present layer's properties to calculate an adiabatic thermal profile for layers below
            TOcean_K[1] = TOcean_K[0] + alphaOcean_pK[0] * TOcean_K[0] / \
                              CpOcean_JkgK[0] / rhoOcean_kgm3[0] * Planet.Ocean.deltaP*1e6
            iStart = 1

        for i in range(iStart, Planet.Steps.nOceanMax):
            Planet.phase[Planet.Steps.nSurfIce+i] = Planet.Ocean.EOS.fn_phase(POcean_MPa[i], TOcean_K[i]).astype(np.int_)
            log.debug(f'il: {Planet.Steps.nSurfIce+i:d}; P_MPa: {POcean_MPa[i]:.3f}; ' +
                      f'T_K: {TOcean_K[i]:.3f}; phase: {Planet.phase[Planet.Steps.nSurfIce+i]:d}')
            if Planet.phase[Planet.Steps.nSurfIce+i] < 2:
                # Liquid water layers -- get fluid properties for the present layer but with the
                # overlaying layer's temperature. Note that we include ice Ih in these layers because
                # ice Ih layers result only from instabilities in phase diagram calculations. There should
                # not be any ice Ih below the ice--ocean interface at Tb.
                rhoOcean_kgm3[i] = Planet.Ocean.EOS.fn_rho_kgm3(POcean_MPa[i], TOcean_K[i])
                CpOcean_JkgK[i] = Planet.Ocean.EOS.fn_Cp_JkgK(POcean_MPa[i], TOcean_K[i])
                alphaOcean_pK[i] = Planet.Ocean.EOS.fn_alpha_pK(POcean_MPa[i], TOcean_K[i])
                kThermOcean_WmK[i] = Planet.Ocean.EOS.fn_kTherm_WmK(POcean_MPa[i], TOcean_K[i])
                # Now use the present layer's properties to calculate an adiabatic thermal profile for layers below
                TOcean_K[i+1] = TOcean_K[i] + alphaOcean_pK[i] * TOcean_K[i] / \
                                CpOcean_JkgK[i] / rhoOcean_kgm3[i] * Planet.Ocean.deltaP*1e6
            else:
                # Undersea high-pressure ices -- we use GetTfreeze here to propagate the layer temperatures.
                # This is based on an assumption that the undersea HP ices are vigorously mixed by
                # two-phase convection, such that each layer is in local equilibrium with the liquid,
                # meaning each layer's temperature is equal to the melting temperature.
                # We implement this by averaging the upper layer temp with the melting temp minus a small offset,
                # to step more gently and avoid overshooting that causes phase oscillations.
                thisPhase = PhaseConv(Planet.phase[Planet.Steps.nSurfIce+i])
                rhoOcean_kgm3[i] = Planet.Ocean.iceEOS[thisPhase].fn_rho_kgm3(POcean_MPa[i], TOcean_K[i])
                CpOcean_JkgK[i] = Planet.Ocean.iceEOS[thisPhase].fn_Cp_JkgK(POcean_MPa[i], TOcean_K[i])
                alphaOcean_pK[i] = Planet.Ocean.iceEOS[thisPhase].fn_alpha_pK(POcean_MPa[i], TOcean_K[i])
                kThermOcean_WmK[i] = Planet.Ocean.iceEOS[thisPhase].fn_kTherm_WmK(POcean_MPa[i], TOcean_K[i])
                TOcean_K[i+1] = np.mean([GetTfreeze(Planet.Ocean.EOS, POcean_MPa[i], TOcean_K[i])
                                         - Planet.Ocean.TfreezeOffset_K, TOcean_K[i]])

        Planet.P_MPa[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = POcean_MPa
        Planet.T_K[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = TOcean_K[:-1]
        Planet.rho_kgm3[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = rhoOcean_kgm3
        Planet.Cp_JkgK[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = CpOcean_JkgK
        Planet.alpha_pK[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = alphaOcean_pK
        Planet.kTherm_WmK[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = kThermOcean_WmK

        # Evaluate remaining physical quantities for ocean layers
        MAbove_kg = np.sum(Planet.MLayer_kg[:Planet.Steps.nSurfIce])
        for i in range(Planet.Steps.nSurfIce, Planet.Steps.nSurfIce + Planet.Steps.nOceanMax):
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1]) * 1e6 / Planet.g_ms2[i-1] / \
                            Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            MAbove_kg += Planet.MLayer_kg[i-1]
            MBelow_kg = Planet.Bulk.M_kg - MAbove_kg
            Planet.g_ms2[i] = Constants.G * MBelow_kg / Planet.r_m[i]**2

        log.info(f'Ocean layers complete. zMax: {Planet.z_m[Planet.Steps.nSurfIce + Planet.Steps.nOceanMax - 1]/1e3:.1f} km, ' +
                 f'upper ice thickness zb: {Planet.zb_km:.3f} km.')

    return Planet


def GetOceanHPIceEOS(Planet, Params, POcean_MPa):
    """ Assign EOS functions for possible high-pressure ices expected in the ocean
        based on the min/max temperatures and pressures we plan to model.

        Args:
            POcean_MPa (float, shape N): Range of pressures expected in the ocean, from
                Pb_MPa to PHydroMax_MPa
        Assigns Planet attributes:
            Ocean.iceEOS
    """

    # Generate linear sets of outside possibilities of P,T combinations for the ocean
    POceanHPices_MPa = POcean_MPa[POcean_MPa>Constants.PminHPices_MPa]
    TOceanHPices_K = np.arange(Planet.Bulk.Tb_K, Planet.Ocean.THydroMax_K, Planet.Ocean.deltaT)
    PHPicesLin_MPa = np.array([P for P in POceanHPices_MPa for _ in TOceanHPices_K])
    THPicesLin_K = np.array([T for _ in POceanHPices_MPa for T in TOceanHPices_K])

    # Stopgap measure to avoid MgSO4 calcs taking ages with the slow Margules formulation phase calcs
    # Remove this if/else block (just do the "else") when a faster phase calculation is implemented!
    if Planet.Ocean.comp == 'MgSO4' or Planet.Sil.poreComp == 'MgSO4':
        # Just load all HP ice phases in case we need them. This part is way faster than Margules phase calcs
        Planet.Ocean.iceEOS['II'] = GetIceEOS(POceanHPices_MPa, TOceanHPices_K, 'II',
                                              porosType=Planet.Ocean.porosType['II'],
                                              phiTop_frac=Planet.Ocean.phiMax_frac['II'],
                                              Pclosure_MPa=Planet.Ocean.Pclosure_MPa['II'],
                                              phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['II'])
        Planet.Ocean.iceEOS['III'] = GetIceEOS(POceanHPices_MPa, TOceanHPices_K, 'III',
                                               porosType=Planet.Ocean.porosType['III'],
                                               phiTop_frac=Planet.Ocean.phiMax_frac['III'],
                                               Pclosure_MPa=Planet.Ocean.Pclosure_MPa['III'],
                                               phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['III'])
        Planet.Ocean.iceEOS['V'] = GetIceEOS(POceanHPices_MPa, TOceanHPices_K, 'V',
                                             porosType=Planet.Ocean.porosType['V'],
                                             phiTop_frac=Planet.Ocean.phiMax_frac['V'],
                                             Pclosure_MPa=Planet.Ocean.Pclosure_MPa['V'],
                                             phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['V'])
        Planet.Ocean.iceEOS['VI'] = GetIceEOS(POceanHPices_MPa, TOceanHPices_K, 'VI',
                                              porosType=Planet.Ocean.porosType['VI'],
                                              phiTop_frac=Planet.Ocean.phiMax_frac['VI'],
                                              Pclosure_MPa=Planet.Ocean.Pclosure_MPa['VI'],
                                              phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['VI'])
    else:
        # Get phase of each P,T combination
        expandPhases = Planet.Ocean.EOS.fn_phase(PHPicesLin_MPa, THPicesLin_K)
        # Check if any of them are not liquid
        if np.any(expandPhases != 0):
            _, _, _, indsIceII, _, indsIceIII, _, indsIceV, _, indsIceVI, _, _, _, _, _, _, _, _, _, _, _ \
                = GetPhaseIndices(expandPhases)

            if(np.size(indsIceII) != 0):
                log.debug('Loading ice II EOS functions for ocean layers...')
                PiceIImin_MPa = PHPicesLin_MPa[indsIceII[0]]
                PiceIImax_MPa = PHPicesLin_MPa[indsIceII[-1]]
                TiceIImin_K = np.min(THPicesLin_K[indsIceII])
                TiceIImax_K = np.max(THPicesLin_K[indsIceII])
                Planet.Ocean.iceEOS['II'] = GetIceEOS(np.linspace(PiceIImin_MPa, PiceIImax_MPa, Planet.Steps.nPsHP),
                                                      np.linspace(TiceIImin_K, TiceIImax_K, Planet.Steps.nTsHP), 'II',
                                                      porosType=Planet.Ocean.porosType['II'],
                                                      phiTop_frac=Planet.Ocean.phiMax_frac['II'],
                                                      Pclosure_MPa=Planet.Ocean.Pclosure_MPa['II'],
                                                      phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['II'])
            if(np.size(indsIceIII) != 0):
                log.debug('Loading ice III EOS functions for ocean layers...')
                PiceIIImin_MPa = PHPicesLin_MPa[indsIceIII[0]]
                PiceIIImax_MPa = PHPicesLin_MPa[indsIceIII[-1]]
                TiceIIImin_K = np.min(THPicesLin_K[indsIceIII])
                TiceIIImax_K = np.max(THPicesLin_K[indsIceIII])
                Planet.Ocean.iceEOS['III'] = GetIceEOS(np.linspace(PiceIIImin_MPa, PiceIIImax_MPa, Planet.Steps.nPsHP),
                                                       np.linspace(TiceIIImin_K, TiceIIImax_K, Planet.Steps.nTsHP), 'III',
                                                       porosType=Planet.Ocean.porosType['III'],
                                                       phiTop_frac=Planet.Ocean.phiMax_frac['III'],
                                                       Pclosure_MPa=Planet.Ocean.Pclosure_MPa['III'],
                                                       phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['III'])
            if(np.size(indsIceV) != 0):
                log.debug('Loading ice V EOS functions for ocean layers...')
                PiceVmin_MPa = PHPicesLin_MPa[indsIceV[0]]
                PiceVmax_MPa = PHPicesLin_MPa[indsIceV[-1]]
                TiceVmin_K = np.min(THPicesLin_K[indsIceV])
                TiceVmax_K = np.max(THPicesLin_K[indsIceV])
                Planet.Ocean.iceEOS['V'] = GetIceEOS(np.linspace(PiceVmin_MPa, PiceVmax_MPa, Planet.Steps.nPsHP),
                                                     np.linspace(TiceVmin_K, TiceVmax_K, Planet.Steps.nTsHP), 'V',
                                                     porosType=Planet.Ocean.porosType['V'],
                                                     phiTop_frac=Planet.Ocean.phiMax_frac['V'],
                                                     Pclosure_MPa=Planet.Ocean.Pclosure_MPa['V'],
                                                     phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['V'])
            if(np.size(indsIceVI) != 0):
                log.debug('Loading ice VI EOS functions for ocean layers...')
                PiceVImin_MPa = PHPicesLin_MPa[indsIceVI[0]]
                PiceVImax_MPa = PHPicesLin_MPa[indsIceVI[-1]]
                TiceVImin_K = np.min(THPicesLin_K[indsIceVI])
                TiceVImax_K = np.max(THPicesLin_K[indsIceVI])
                Planet.Ocean.iceEOS['VI'] = GetIceEOS(np.linspace(PiceVImin_MPa, PiceVImax_MPa, Planet.Steps.nPsHP),
                                                      np.linspace(TiceVImin_K, TiceVImax_K, Planet.Steps.nTsHP), 'VI',
                                                      porosType=Planet.Ocean.porosType['VI'],
                                                      phiTop_frac=Planet.Ocean.phiMax_frac['VI'],
                                                      Pclosure_MPa=Planet.Ocean.Pclosure_MPa['VI'],
                                                      phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['VI'])
        else:
            log.debug('No high-pressure ices found in ocean layers.')

    return Planet, Params


def InnerLayers(Planet, Params):
    """ Geophysical and thermodynamic calculations for silicate and core layers
        Calculates state variables of the layer with each pressure step

        Assigns Planet attributes:
            Steps.nTotal, all layer arrays
    """
    if Planet.Do.VALID:
        if Planet.Do.CONSTANT_INNER_DENSITY or Params.SKIP_INNER:
            Planet, mantleProps, coreProps = CalcMoIConstantRho(Planet, Params)
        else:
            Planet, mantleProps, coreProps = CalcMoIWithEOS(Planet, Params)

        if(Planet.Steps.nHydro <= Planet.Steps.nSurfIce) and (not Planet.Do.NO_H2O):
            log.warning('For these run settings, the hydrosphere is entirely frozen and contains only surface ice.')
        Planet.Steps.nTotal = Planet.Steps.nHydro + Planet.Steps.nSil + Planet.Steps.nCore

        log.debug('Evaluating remaining quantities for layer arrays...')
        # Extend Planet layer arrays to make space for silicate and possible core layers
        extend = np.zeros(Planet.Steps.nSil + Planet.Steps.nCore)
        Planet.P_MPa = np.concatenate((Planet.P_MPa[:Planet.Steps.nHydro], extend))
        Planet.T_K = np.concatenate((Planet.T_K[:Planet.Steps.nHydro], extend))
        Planet.r_m = np.concatenate((Planet.r_m[:Planet.Steps.nHydro], extend, [0.0]))
        Planet.phase = np.concatenate((Planet.phase[:Planet.Steps.nHydro], extend.astype(np.int_)))
        Planet.rho_kgm3 = np.concatenate((Planet.rho_kgm3[:Planet.Steps.nHydro], extend))
        Planet.Cp_JkgK = np.concatenate((Planet.Cp_JkgK[:Planet.Steps.nHydro], extend))
        Planet.alpha_pK = np.concatenate((Planet.alpha_pK[:Planet.Steps.nHydro], extend))
        Planet.kTherm_WmK = np.concatenate((Planet.kTherm_WmK[:Planet.Steps.nHydro], extend))
        Planet.g_ms2 = np.concatenate((Planet.g_ms2[:Planet.Steps.nHydro], extend))
        Planet.phi_frac = np.concatenate((Planet.phi_frac[:Planet.Steps.nHydro], extend))
        Planet.Htidal_Wm3 = np.concatenate((Planet.Htidal_Wm3[:Planet.Steps.nHydro], extend))
        Planet.z_m = np.concatenate((Planet.z_m[:Planet.Steps.nHydro], extend, [0.0]))
        Planet.Ppore_MPa = np.concatenate((Planet.Ppore_MPa[:Planet.Steps.nHydro], extend))
        Planet.rhoMatrix_kgm3 = np.concatenate((Planet.rhoMatrix_kgm3[:Planet.Steps.nHydro], extend))
        Planet.rhoPore_kgm3 = np.concatenate((Planet.rhoPore_kgm3[:Planet.Steps.nHydro], extend))
        Planet.MLayer_kg = np.concatenate((Planet.MLayer_kg[:Planet.Steps.nHydro], extend))

        # Unpack results from MoI calculations
        iOS = Planet.Steps.nHydro
        iSC = Planet.Steps.nHydro + Planet.Steps.nSil
        Planet.P_MPa[iOS:iSC], Planet.T_K[iOS:iSC], Planet.r_m[iOS:iSC], Planet.rho_kgm3[iOS:iSC], \
        Planet.g_ms2[iOS:iSC], Planet.phi_frac[iOS:iSC], Planet.Htidal_Wm3[iOS:iSC], Planet.kTherm_WmK[iOS:iSC], \
        Planet.Ppore_MPa[iOS:iSC], Planet.rhoMatrix_kgm3[iOS:iSC], Planet.rhoPore_kgm3[iOS:iSC], \
        Planet.MLayer_kg[iOS:iSC], phasePore \
            = mantleProps

        iCC = Planet.Steps.nTotal
        if Planet.Do.Fe_CORE:
            # Unpack results from MoI calculations
            Planet.P_MPa[iSC:iCC], Planet.T_K[iSC:iCC], Planet.r_m[iSC:iCC], Planet.g_ms2[iSC:iCC], Planet.rho_kgm3[iSC:iCC], \
            Planet.Cp_JkgK[iSC:iCC], Planet.alpha_pK[iSC:iCC], Planet.kTherm_WmK[iSC:iCC], Planet.MLayer_kg[iSC:iCC] \
                = coreProps

        Planet.z_m[iOS:iCC+1] = Planet.Bulk.R_m - Planet.r_m[iOS:iCC+1]

        # Assign phase values for silicates and core
        Planet.phase[iOS:iSC] = Constants.phaseSil + phasePore
        Planet.phase[iSC:iCC] = Constants.phaseFe

        # Record ocean layer thickness
        if Planet.Do.NO_H2O or not np.any(Planet.phase == 0):
            Planet.D_km = 0
        else:
            # Get first index of phase changing from 0 to something different ---
            # this is the bottom of the contiguous ocean layer.
            iOceanBot = next(i for i,phase in enumerate(Planet.phase[:Planet.Steps.nHydro]) if phase == 0 and phase != Planet.phase[i+1])
            Planet.D_km = Planet.z_m[iOceanBot + 1]/1e3 - Planet.zb_km

        # Calculate total salt and water masses
        Planet.Mcore_kg = np.sum(Planet.MLayer_kg[iSC:iCC])
        Planet.Mrock_kg = 4/3*np.pi * np.sum((Planet.r_m[iOS:iSC]**3 - Planet.r_m[iOS+1:iSC+1]**3)
                                             * Planet.rhoMatrix_kgm3[iOS:iSC] * (1 - Planet.phi_frac[iOS:iSC]))
        # Next, fetch the phase IDs of the silicate layers, which are incremented when
        # they contain non-liquid phases.
        silPhases = Planet.phase[iOS:iSC]
        # Add matrix density for ice phases; rhoMatrix is not set for phase == 0, so we
        # can safely include those layers in the sum. We also must include ices in the
        # pore space of silicates.
        Planet.Mice_kg = 4/3*np.pi * (np.sum((Planet.r_m[0:iOS]**3 - Planet.r_m[1:iOS+1]**3)
                                             * Planet.rhoMatrix_kgm3[0:iOS] * (1 - Planet.phi_frac[0:iOS]))
                                    + np.sum((Planet.r_m[iOS:iSC][silPhases != Constants.phaseSil]**3
                                            - Planet.r_m[iOS+1:iSC+1][silPhases != Constants.phaseSil]**3)
                                             * Planet.rhoPore_kgm3[iOS:iSC][silPhases != Constants.phaseSil]
                                             * Planet.phi_frac[iOS:iSC][silPhases != Constants.phaseSil]))
        # The remainder is the ocean fluids, including H2O and salts and pore spaces.
        Planet.Mfluid_kg = Planet.Mtot_kg - Planet.Mcore_kg - Planet.Mrock_kg - Planet.Mice_kg
        # Get the mass contained in clathrate layers and just the trapped gas
        Planet.Mclath_kg = np.sum(Planet.MLayer_kg[abs(Planet.phase) == Constants.phaseClath])
        Planet.MclathGas_kg = Planet.Mclath_kg * Constants.clathGasFrac_ppt / 1e3
        # The portion just in the ocean is simple to evaluate:
        Planet.Mocean_kg = np.sum(Planet.MLayer_kg[Planet.phase == 0])
        # The difference is the amount contained in the pore space:
        Planet.MporeFluid_kg = Planet.Mfluid_kg - Planet.Mocean_kg
        # Multiply mass concentration of solute to get total mass of salt in ocean
        Planet.MoceanSalt_kg = Planet.Mfluid_kg * Planet.Ocean.wOcean_ppt / 1e3
        # Record the mass of salt in the pore space in case we want to track it separately
        Planet.MporeSalt_kg = Planet.MporeFluid_kg * Planet.Sil.wPore_ppt / 1e3
        # Combine these to get total salt content
        Planet.Msalt_kg = Planet.MoceanSalt_kg + Planet.MporeSalt_kg
        # The remainder, plus ice mass excluding gasses in clathrates,
        # is the total mass contained in water molecules for the body
        Planet.MH2O_kg = Planet.Mfluid_kg - Planet.Msalt_kg + Planet.Mice_kg - Planet.MclathGas_kg

        # Get the mean density of ocean layers and conducting/convecting upper ice layers
        Planet.VLayer_m3 = 4/3*np.pi * (Planet.r_m[:-1]**3 - Planet.r_m[1:]**3)
        Planet.Ocean.Vtot_m3 = np.sum(Planet.VLayer_m3[Planet.phase == 0])
        if Planet.Do.NO_H2O:
            Planet.Ocean.rhoMean_kgm3 = np.nan
            Planet.Ocean.Tmean_K = np.nan
        else:
            Planet.Ocean.rhoMean_kgm3 = Planet.Mocean_kg / Planet.Ocean.Vtot_m3
            # Get average ocean temperature by summing the total heat energy in the
            # ocean and dividing by the total heat storage capacity
            oceanHeat_pK = Planet.Cp_JkgK[Planet.phase == 0] * Planet.MLayer_kg[Planet.phase == 0]
            Planet.Ocean.Tmean_K = np.sum(Planet.T_K[Planet.phase == 0] * oceanHeat_pK) / np.sum(oceanHeat_pK)

        # Get mean tidal heating in silicate layers
        Planet.Sil.HtidalMean_Wm3 = np.mean(Planet.Htidal_Wm3[iOS:iSC])

        # Check for any negative temperature gradient (indicates non-equilibrium conditions)
        gradTneg = np.where(np.diff(Planet.T_K) < 0)
        if np.any(gradTneg) and not Params.SKIP_INNER:
            log.warning(f'Negative temperature gradient starting at index {gradTneg[0][0]:d}. This indicates that ' +
                         'internal heating parameters Qrad_Wkg and/or Htidal_Wm3 are likely set too high to be consistent ' +
                         'with the heat flux into the ocean. The current configuration represents a ' +
                         'non-equilibrium state.')

    else:
        # Set remaining quantities that are still None if the profile is invalid:
        Planet.Ocean.Tmean_K, Planet.Sil.HtidalMean_Wm3, Planet.Ocean.rhoMean_kgm3, Planet.Ocean.Vtot_m3, \
        Planet.D_km, Planet.CMR2mean, Planet.CMR2less, Planet.CMR2more, Planet.MH2O_kg, Planet.MclathGas_kg, \
        Planet.Mclath_kg, Planet.Mcore_kg, Planet.Mfluid_kg, Planet.Mice_kg, Planet.MoceanSalt_kg, \
        Planet.Mocean_kg, Planet.MporeFluid_kg, Planet.MporeSalt_kg, Planet.Mrock_kg, Planet.Msalt_kg, \
        Planet.Mtot_kg, Planet.Sil.Rmean_m, Planet.Sil.Rrange_m, Planet.Core.Rmean_m, Planet.Core.Rrange_m, \
        Planet.Sil.rhoMean_kgm3, Planet.Core.rhoMean_kgm3, Planet.Sil.GSmean_GPa, Planet.Core.GSmean_GPa \
            = (np.nan for _ in range(29))
        Planet.Sil.Rtrade_m, Planet.Core.Rtrade_m, Planet.Sil.rhoTrade_kgm3 \
            = (np.array([np.nan]) for _ in range(3))
        Planet.Steps.nHydro = 1
        Planet.Steps.nSil = 0
        Planet.Steps.nTotal = 1

    return Planet


def CalcMoIConstantRho(Planet, Params):
    """ Find the relative sizes of silicate, core, and hydrosphere layers that are
        consistent with the measured moment of inertia, based on calculated hydrosphere
        properties and assumptions about the silicate and possible core layers.

        Assigns Planet attributes:
            CMR2mean, Sil.RsilMean_m, Sil.RsilRange_m, Core.RFeMean_m, Core.RFeRange_m, Steps.nHydro
    """
    log.debug('Finding MoI consistent with measured value for constant-density inner layers...')
    # Get MR^2 -- we will need to divide each C by this later.
    MR2_kgm2 = Planet.Bulk.M_kg * Planet.Bulk.R_m**2

    # Get final number of layers modeled in "overshoot" hydrosphere
    if Planet.Do.NO_H2O:
        Planet.Steps.iSilStart = 0
        nHydroActual = 2
    else:
        nHydroActual = Planet.Steps.nSurfIce + Planet.Steps.nOceanMax
    # Find contribution to axial moment of inertia C from each ocean layer
    dChydro_kgm2 = 8*np.pi/15 * Planet.rho_kgm3[:-1] * (Planet.r_m[:-1]**5 - Planet.r_m[1:]**5)
    # Find total mass contained above each hydrosphere layer
    MAbove_kg = np.array([np.sum(Planet.MLayer_kg[:i]) for i in range(nHydroActual)])
    # Find volume of a full sphere of silicate corresponding to each valid layer
    VsilSphere_m3 = 4/3*np.pi * Planet.r_m[Planet.Steps.iSilStart:]**3

    if Planet.Do.Fe_CORE:
        # Find core bulk density based on assumed sulfide content (Eq X of Styczinski et al., 2022)
        rhoCore_kgm3 = Planet.Core.rhoFeS_kgm3 * Planet.Core.rhoFe_kgm3 \
            * (Planet.Core.xFeS * Constants.mFeS_gmol + (1 - Planet.Core.xFeS) * Constants.mFe_gmol ) \
            / (Planet.Core.xFeS * Constants.mFeS_gmol * Planet.Core.rhoFe_kgm3 + (1 - Planet.Core.xFeS) * Constants.mFe_gmol * Planet.Core.rhoFeS_kgm3)
            # / (Planet.Core.xFeS * (Planet.Core.rhoFe_kgm3 - Planet.Core.rhoFeS_kgm3) + Planet.Core.rhoFeS_kgm3)  # Vance et al. (2014) Eq. 10
        # Calculate core volume for a silicate layer with outer radius equal to bottom of each hydrosphere layer
        # and inner radius equal to the core radius
        VCore_m3 = np.array([(Planet.Bulk.M_kg - MAbove_kg[i] - VsilSphere_m3[i-Planet.Steps.iSilStart] *
                             Planet.Sil.rhoSilWithCore_kgm3) / (rhoCore_kgm3 - Planet.Sil.rhoSilWithCore_kgm3)
                             for i in range(Planet.Steps.iSilStart, nHydroActual-1)])
        # Find values for which the silicate radius is too large
        try:
            nTooBig = next((i[0] for i, val in np.ndenumerate(VCore_m3) if val>0))
        except StopIteration:
            msg = f'Failed to find a core size consistent with rhoSil = {Planet.Sil.rhoSilWithCore_kgm3:.1f} kg/m3 ' + \
                  f'and xFeS = {Planet.Core.xFeS:.3f} for PHydroMax_MPa = {Planet.Ocean.PHydroMax_MPa:.1f}. ' + \
                   'Core size will be set to zero.'
            if Params.DO_EXPLOREOGRAM:
                log.debug(msg)
            else:
                log.warning(msg)
            nTooBig = 0
            rCore_m = np.zeros_like(VCore_m3)
        else:
            # Calculate corresponding core radii based on above density
            rCore_m = (VCore_m3[nTooBig:]*3/4/np.pi)**(1/3)

        # Assign fixed density to an array for dual-use code looking for compatible C/MR^2
        rhoSil_kgm3 = np.ones_like(rCore_m) * Planet.Sil.rhoSilWithCore_kgm3
    else:
        # Find silicate density consistent with observed bulk mass for each radius
        rhoSil_kgm3 = np.array([(Planet.Bulk.M_kg - MAbove_kg[i]) / VsilSphere_m3[i-Planet.Steps.iSilStart]
                              for i in range(Planet.Steps.iSilStart, nHydroActual-1)])
        # Density of silicates is scaled to fit the total mass, so there is no nTooBig in this case.
        nTooBig = 0
        # Set core radius and density to zero so calculations can proceed
        rCore_m = np.zeros(nHydroActual-1 - Planet.Steps.iSilStart)
        rhoCore_kgm3 = 0

    # Calculate C for a mantle extending up to each hydrosphere layer in turn
    C_kgm2 = np.zeros(nHydroActual - 1)
    C_kgm2[Planet.Steps.iSilStart + nTooBig:] = [np.sum(dChydro_kgm2[:i + Planet.Steps.iSilStart + nTooBig + 1]) +
            8*np.pi/15 * rhoSil_kgm3[i] * (Planet.r_m[i + Planet.Steps.iSilStart + nTooBig]**5 - rCore_m[i]**5) +
            8*np.pi/15 * rhoCore_kgm3 * rCore_m[i]**5
            for i in range(nHydroActual - Planet.Steps.iSilStart - nTooBig - 1)]
    CMR2 = C_kgm2 / MR2_kgm2

    CMR2inds = [i[0] for i, valCMR2 in np.ndenumerate(CMR2)
                 if valCMR2 > Planet.Bulk.Cmeasured - Planet.Bulk.Cuncertainty
                and valCMR2 < Planet.Bulk.Cmeasured + Planet.Bulk.Cuncertainty]

    if len(CMR2inds) == 0:
        if Planet.Do.NO_H2O:
            suggestion = '\nTry adjusting properties of silicates and core to get C/MR^2 values in range.'
        else:
            suggestion = '\nTry adjusting properties of silicates and core to get C/MR^2 values in range. ' + \
                         'Increasing PHydroMax_MPa can also lower C/MR^2 values.'
        msg = f'No MoI found matching C/MR^2 = {Planet.Bulk.Cmeasured:.3f}±{Planet.Bulk.Cuncertainty:.3f}. ' + \
                  f'Min: {np.min(CMR2[CMR2>0]):.3f}, Max: {np.max(CMR2):.3f}.'
        if Params.ALLOW_BROKEN_MODELS:
            fullMsg = msg + suggestion + ' Params.ALLOW_BROKEN_MODELS is True, so calculations will proceed with many values set to nan.'
            if Params.DO_EXPLOREOGRAM:
                log.info(msg)
            else:
                log.error(fullMsg)
            Planet.Do.VALID = False
        else:
            raise RuntimeError(msg + suggestion)

        nans = np.array([np.nan])
        Planet.CMR2mean = np.nan
        Planet.CMR2less = Planet.CMR2mean
        Planet.CMR2more = Planet.CMR2mean
        Planet.Sil.rhoTrade_kgm3 = nans
        Planet.Sil.Rmean_m = np.nan
        Planet.Sil.Rtrade_m = nans
        Planet.Sil.Rrange_m = np.nan
        Planet.Core.Rmean_m = np.nan
        Planet.Core.Rtrade_m = nans
        Planet.Core.Rrange_m = np.nan
        Planet.Steps.nSil = Planet.Steps.nSilMax
        # Use Rset_m to indicate that we have already determined the core size in using SilicateLayers
        Planet.Core.Rset_m = np.nan
        iCMR2inner = 0
        if Planet.Do.NO_H2O:
            Planet.Steps.nHydro = 0
        else:
            Planet.Steps.nHydro = Planet.Steps.nOceanMax

    else:
        # Find the C/MR^2 value most closely matching the measured value
        CMR2diff = np.abs(CMR2[CMR2inds] - Planet.Bulk.Cmeasured)
        # Get index of closest match in CMR2inds
        iCMR2ind = np.argmin(CMR2diff)
        # Find Planet array index corresponding to closest matching value
        iCMR2 = CMR2inds[iCMR2ind]
        iCMR2inner = iCMR2 - Planet.Steps.iSilStart - nTooBig
        CMR2indsInner = [ind - Planet.Steps.iSilStart - nTooBig for ind in CMR2inds]
        # Record the best-match C/MR^2 value
        Planet.CMR2mean = CMR2[iCMR2]
        # We don't have neighboring values because we used the MoI to calculate properties
        Planet.CMR2less = Planet.CMR2mean
        Planet.CMR2more = Planet.CMR2mean
        # Record interior sizes
        Planet.Sil.rhoTrade_kgm3 = rhoSil_kgm3[CMR2indsInner]
        Planet.Sil.Rmean_m = Planet.r_m[iCMR2]
        Planet.Sil.Rtrade_m = Planet.r_m[CMR2inds]
        Planet.Sil.Rrange_m = np.max(Planet.Sil.Rtrade_m) - np.min(Planet.Sil.Rtrade_m)
        Planet.Core.Rmean_m = rCore_m[iCMR2inner]
        Planet.Core.Rtrade_m = rCore_m[CMR2indsInner]
        Planet.Core.Rrange_m = np.max(Planet.Core.Rtrade_m) - np.min(Planet.Core.Rtrade_m)
        # Now we finally know how many layers there are in the hydrosphere
        Planet.Steps.nHydro = iCMR2
        # Number of steps in the silicate layer is fixed for the constant-density approach
        Planet.Steps.nSil = Planet.Steps.nSilMax
        # Use Rset_m to indicate that we have already determined the core size in using SilicateLayers
        Planet.Core.Rset_m = Planet.Core.Rmean_m + 0.0

    if not Params.SKIP_INNER:
        Planet.Sil.fn_phi_frac = GetphiCalc(Planet.Sil.phiRockMax_frac, Planet.Sil.EOS.fn_phi_frac, Planet.Sil.phiMin_frac)
        Planet.Sil.fn_Htidal_Wm3 = GetHtidalFunc(Planet.Sil.Htidal_Wm3)  # Placeholder until we implement a self-consistent calc
        # Evaluate the silicate EOS for each layer
        indsSilValid, nProfiles, Psil_MPa, Tsil_K, rSil_m, rhoSilEOS_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, \
        phiSil_frac, HtidalSil_Wm3, kThermSil_WmK, Ppore_MPa, rhoSilMatrix_kgm3, rhoPore_kgm3, phasePore \
            = SilicateLayers(Planet, Params)
        nSilTooBig = nProfiles - np.size(indsSilValid)

        # Fill core/mantle trade arrays and set mean values consistent with MoI
        MtotSil_kg = np.sum(MLayerSil_kg)
        Planet.Sil.rhoMean_kgm3 = MtotSil_kg / (4/3*np.pi * (rSil_m[0,0]**3 - rSil_m[0,-1]**3))

        if Planet.Do.Fe_CORE:
            # Evaluate the core EOS for each layer
            _, Pcore_MPa, Tcore_K, rCoreEOS_m, rhoCore_kgm3, MLayerCore_kg, gCore_ms2, CpCore_JkgK, alphaCore_pK, \
                kThermCore_WmK = IronCoreLayers(Planet, Params,
                               indsSilValid, nSilTooBig, nProfiles, Psil_MPa, Tsil_K, rSil_m, MAboveSil_kg, gSil_ms2)

            MtotCore_kg = np.sum(MLayerCore_kg)
            Planet.Core.rhoMean_kgm3 = MtotCore_kg / VCore_m3[nTooBig:][iCMR2inner]

            coreProps = (Pcore_MPa, Tcore_K, rCoreEOS_m[0,:-1], gCore_ms2, rhoCore_kgm3, CpCore_JkgK, alphaCore_pK,
                         kThermCore_WmK, MLayerCore_kg)
        else:
            MtotCore_kg = 0
            Planet.Core.rhoMean_kgm3 = 0
            Planet.Core.Rtrade_m = np.zeros_like(Planet.Sil.Rtrade_m)
            Planet.Core.Rrange_m = 0
            coreProps = None

        Planet.Mtot_kg = np.sum(Planet.MLayer_kg[:iCMR2]) + MtotSil_kg + MtotCore_kg
        if not np.isnan(Planet.CMR2mean):
            log.info(f'Found matching MoI of {Planet.CMR2mean:.4f} ' +
                     f'(C/MR^2 = {Planet.Bulk.Cmeasured:.4f}±{Planet.Bulk.Cuncertainty:.4f}) for ' +
                     f'R_sil = {Planet.Sil.Rmean_m / Planet.Bulk.R_m:.2f} R_{Planet.name[0]}, ' +
                     f'R_core = {Planet.Core.Rmean_m / Planet.Bulk.R_m:.2f} R_{Planet.name[0]}, ' +
                     f'rho_sil (found) = {rhoSil_kgm3[iCMR2inner]:.0f} kg/m^3, ' +
                     f'rho_sil (actual) = {Planet.Sil.rhoMean_kgm3:.0f} kg/m^3, ' +
                     f'M_tot = {Planet.Mtot_kg/Planet.Bulk.M_kg:.4f} M_{Planet.name[0]}.')
            log.warning('Because silicate and core properties were determined from the EOS after finding their ' +
                        'sizes by assuming constant densities, the body mass may not match the measured value.')

    else:
        if not np.isnan(Planet.CMR2mean):
            log.info(f'Found matching MoI of {Planet.CMR2mean:.4f} ' +
                     f'(C/MR^2 = {Planet.Bulk.Cmeasured:.4f}±{Planet.Bulk.Cuncertainty:.4f}) for ' +
                     f'R_sil = {Planet.Sil.Rmean_m / Planet.Bulk.R_m:.2f} R_{Planet.name[0]}, ' +
                     f'R_core = {Planet.Core.Rmean_m / Planet.Bulk.R_m:.2f} R_{Planet.name[0]}, ' +
                     f'rho_sil = {rhoSil_kgm3[iCMR2inner]:.0f} kg/m^3, ' +
                     f'M_tot = {1.0:.5f} M_{Planet.name[0]} (fixed).')
        if not Params.DO_EXPLOREOGRAM:
            log.debug('Params.SKIP_INNER is True, assigning interior properties to 0.')
        Psil_MPa, Tsil_K, rhoSilEOS_kgm3, gSil_ms2, phiSil_frac, kThermSil_WmK, Ppore_MPa, rhoSilMatrix_kgm3, \
            rhoPore_kgm3, HtidalSil_Wm3, MLayerSil_kg \
            = (np.zeros(Planet.Steps.nSil) for _ in range(11))
        phasePore = np.zeros(Planet.Steps.nSil, dtype=np.int_)
        rSil_m = np.zeros((1, Planet.Steps.nSil+1))
        rSil_m[0,0] = Planet.Sil.Rmean_m
        coreProps = (np.zeros(Planet.Steps.nCore) for _ in range(9))
        Planet.Sil.rhoMean_kgm3 = rhoSil_kgm3[iCMR2inner]
        Planet.Core.rhoMean_kgm3 = rhoCore_kgm3
        Planet.Mtot_kg = Planet.Bulk.M_kg

    mantleProps = (Psil_MPa, Tsil_K, rSil_m[0,:-1], rhoSilEOS_kgm3, gSil_ms2, phiSil_frac, HtidalSil_Wm3, kThermSil_WmK,
                   Ppore_MPa, rhoSilMatrix_kgm3, rhoPore_kgm3, MLayerSil_kg, phasePore)

    return Planet, mantleProps, coreProps


def CalcMoIWithEOS(Planet, Params):
    """ Find the relative sizes of silicate, core, and hydrosphere layers that are
        consistent with the measured moment of inertia, based on calculated hydrosphere
        properties and EOS data for assumed mantle and core compositions output by Perple_X.

        For models with a core, consistency with the MoI and body mass is determined with
        Sil.EOS and Core.EOS, Ocean.QfromMantle_W, Sil.Qrad_Wkg, and Sil.Htidal_Wm3. The radii
        of the silicate and core layers are treated as the free variables in 2 equations to match
        the 2 "unknowns", MoI and body mass.
        For models without a core, consistency with the MoI and body mass is determined with
        Sil.EOS, Ocean.QfromMantle_W, and Sil.Qrad_Wkg. Sil.Htidal_Wm3 and the radius of the
        silicate layer are treated as the free variables in 2 equations to match MoI and body mass.
        Sil.Htidal_Wm3 is started at 0 and increased until the thermal profile extends beyond the
        T domain of the Perplex_X EOS.

        Assigns Planet attributes:
            CMR2mean, Sil.RsilMean_m, Sil.RsilRange_m, Core.RFeMean_m, Core.RFeRange_m, Steps.nHydro, Steps.nSil,
            all layer arrays
    """
    log.debug('Finding MoI consistent with measured value...')
    # Get MR^2 -- we will need to divide each C by this later.
    MR2_kgm2 = Planet.Bulk.M_kg * Planet.Bulk.R_m**2

    # Find contribution to axial moment of inertia C from each ocean layer
    dCfromH2O_kgm2 = 8*np.pi/15 * Planet.rho_kgm3[:-1] * (Planet.r_m[:-1]**5 - Planet.r_m[1:]**5)

    if Planet.Do.Fe_CORE:
        Planet.Sil.fn_Htidal_Wm3 = GetHtidalFunc(Planet.Sil.Htidal_Wm3)  # Placeholder until we implement a self-consistent calc
        Planet.Sil.fn_phi_frac = GetphiCalc(Planet.Sil.phiRockMax_frac, Planet.Sil.EOS.fn_phi_frac, Planet.Sil.phiMin_frac)
        # Propagate the silicate EOS from each hydrosphere layer to the center of the body
        log.debug(f'Propagating silicate EOS for each possible mantle size...')
        indsSilValid, nProfiles, Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, \
        phiSil_frac, HtidalSil_Wm3, kThermSil_WmK, PsilPore_MPa, rhoSilMatrix_kgm3, rhoSilPore_kgm3, phaseSilPore \
            = SilicateLayers(Planet, Params)
        nSilTooBig = nProfiles - np.size(indsSilValid)
        # Propagate the core EOS from each silicate layer at the max core radius to the center of the body
        nSilFinal, Pcore_MPa, Tcore_K, rCore_m, rhoCore_kgm3, MLayerCore_kg, gCore_ms2, CpCore_JkgK, alphaCore_pK, \
            kThermCore_WmK = IronCoreLayers(Planet, Params,
                           indsSilValid, nSilTooBig, nProfiles, Psil_MPa, Tsil_K, rSil_m, MAboveSil_kg, gSil_ms2)

        dCfromCore_kgm2 = 8*np.pi/15 * rhoCore_kgm3 * (rCore_m[:,:-1]**5 - rCore_m[:,1:]**5)
        Ccore_kgm2 = np.array([np.sum(dCfromCore_kgm2[i,:]) for i in range(nProfiles - nSilTooBig)])

        # Get indices of valid silicate portions of the layer profile
        iValid = np.array([Planet.Steps.iSilStart + i for i in indsSilValid]).astype(np.int_)
        C_kgm2 = np.zeros(nProfiles + Planet.Steps.iSilStart)
    else:
        # Propagate the silicate EOS from each hydrosphere layer to the center of the body

        if Planet.Do.POROUS_ROCK and not Planet.Do.FIXED_POROSITY:
            # In this case, we use the user-specified phiTop_frac value as a "middle" option,
            # and vary between 1/Planet.Sil.phiRangeMult times the difference from this middle
            # value to the endpoints, 0 and 1, to get the range of porosities to model.
            # For example, the user specifies 0.22 and phiRangeMult = 5. phiTop will then
            # take values from a min value of 0.22 - 0.22 / 5 = 0.176 to a max value of
            # 0.22 + (1 - 0.22) / 5 = 0.376.
            thisHtidal_Wm3 = Planet.Sil.Htidal_Wm3
            phiMin_frac = Planet.Sil.phiRockMax_frac - Planet.Sil.phiRockMax_frac / Planet.Sil.phiRangeMult
            phiMax_frac = Planet.Sil.phiRockMax_frac + (1 - Planet.Sil.phiRockMax_frac) / Planet.Sil.phiRangeMult
            multphi_frac = (phiMax_frac/phiMin_frac)**(1/Planet.Steps.nPoros)
            log.debug(f'Propagating silicate EOS for each possible mantle size and porosity from phiVac = {phiMin_frac:.3f} to {phiMax_frac:.3f}...')
        else:
            # In this case, we will use Sil.HtidalMin_Wm3 and Sil.deltaHtidal_logUnits to get
            # a valid set of profiles.
            HtidalStart_Wm3 = Planet.Sil.HtidalMin_Wm3
            multHtidal_Wm3 = 10**Planet.Sil.deltaHtidal_logUnits
            thisHtidal_Wm3 = 0
            log.debug(f'Propagating silicate EOS for each possible mantle size and heating from Htidal = {thisHtidal_Wm3:.2e} to {Planet.Sil.HtidalMax_Wm3:.2e} W/m^3...')
            phiMin_frac = Planet.Sil.phiRockMax_frac
            phiMax_frac = phiMin_frac

        nProfiles = 0

        # Start at minimum tidal heating and initialize arrays
        thisphiTop_frac = phiMin_frac + 0.0
        Planet.Sil.fn_Htidal_Wm3 = GetHtidalFunc(thisHtidal_Wm3)  # Placeholder until we implement a self-consistent calc
        Planet.Sil.fn_phi_frac = GetphiCalc(Planet.Sil.phiRockMax_frac, Planet.Sil.EOS.fn_phi_frac, Planet.Sil.phiMin_frac)
        indsSilValidTemp, _, PsilTemp_MPa, TsilTemp_K, rSilTemp_m, rhoSilTemp_kgm3, MLayerSilTemp_kg, MAboveSilTemp_kg, gSilTemp_ms2, \
        phiSilTemp_frac, HtidalSilTemp_Wm3, kThermSilTemp_WmK, PporeTemp_MPa, rhoSilMatrixTemp_kgm3, rhoSilPoreTemp_kgm3, phaseSilPoreTemp \
            = SilicateLayers(Planet, Params)
        if (np.size(indsSilValidTemp) == 0):
            raise RuntimeError('No silicate mantle size was less than the total body mass for the initialization ' +
                               f'setting of {thisHtidal_Wm3} W/m^3 tidal heating and the expected maximum porosity ' +
                               f' of {thisphiTop_frac}. Try adjusting run settings that affect mantle density, ' +
                               'like porosity, silicate composition, and radiogenic heat flux.')

        # Initialize empty arrays for silicate layer properties
        Psil_MPa, Tsil_K, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, phiSil_frac, HtidalSil_Wm3, \
            kThermSil_WmK, PsilPore_MPa, rhoSilMatrix_kgm3, rhoSilPore_kgm3, phaseSilPore, phiTop_frac \
            = (np.empty((0, Planet.Steps.nSilMax)) for _ in range(14))
        rSil_m = np.empty((0, Planet.Steps.nSilMax+1))
        iValid = [indsSilValidTemp + Planet.Steps.iSilStart]

        while np.size(indsSilValidTemp) != 0 and (thisHtidal_Wm3 <= Planet.Sil.HtidalMax_Wm3 and thisphiTop_frac <= phiMax_frac):
            nProfiles += 1
            # Append the newest mass matches alongside existing ones --
            # do this in the loop so that we don't append the final,
            # invalid profile
            Psil_MPa = np.vstack([Psil_MPa, PsilTemp_MPa[indsSilValidTemp,:]])
            Tsil_K = np.vstack([Tsil_K, TsilTemp_K[indsSilValidTemp,:]])
            rSil_m = np.vstack([rSil_m, rSilTemp_m[indsSilValidTemp,:]])
            rhoSil_kgm3 = np.vstack([rhoSil_kgm3, rhoSilTemp_kgm3[indsSilValidTemp,:]])
            MLayerSil_kg = np.vstack([MLayerSil_kg, MLayerSilTemp_kg[indsSilValidTemp,:]])
            MAboveSil_kg = np.vstack([MAboveSil_kg, MAboveSilTemp_kg[indsSilValidTemp,:]])
            gSil_ms2 = np.vstack([gSil_ms2, gSilTemp_ms2[indsSilValidTemp,:]])
            phiSil_frac = np.vstack([phiSil_frac, phiSilTemp_frac[indsSilValidTemp,:]])
            HtidalSil_Wm3 = np.vstack([HtidalSil_Wm3, HtidalSilTemp_Wm3[indsSilValidTemp,:]])
            kThermSil_WmK = np.vstack([kThermSil_WmK, kThermSilTemp_WmK[indsSilValidTemp,:]])
            PsilPore_MPa = np.vstack([PsilPore_MPa, PporeTemp_MPa[indsSilValidTemp,:]])
            rhoSilMatrix_kgm3 = np.vstack([rhoSilMatrix_kgm3, rhoSilMatrixTemp_kgm3[indsSilValidTemp,:]])
            rhoSilPore_kgm3 = np.vstack([rhoSilPore_kgm3, rhoSilPoreTemp_kgm3[indsSilValidTemp,:]])
            phaseSilPore = np.vstack([phaseSilPore, phaseSilPoreTemp[indsSilValidTemp,:]])
            phiTop_frac = np.append(phiTop_frac, thisphiTop_frac)

            if Planet.Do.POROUS_ROCK and not Planet.Do.FIXED_POROSITY:
                rSilOuter_m = rSilTemp_m[indsSilValidTemp, 0]
                log.debug(f'Silicate match for phiTop = {thisphiTop_frac:.3f} with ' +
                          f'rSil = {rSilOuter_m / 1e3:.1f} km ({rSilOuter_m / Planet.Bulk.R_m:.3f} R_{Planet.name[0]}).')
                thisphiTop_frac = phiMin_frac * multphi_frac**nProfiles
                Planet.Sil.fn_phi_frac.update(thisphiTop_frac)
            else:
                rSilOuter_m = rSilTemp_m[indsSilValidTemp, 0]
                log.debug(f'Silicate match for Htidal = {thisHtidal_Wm3:.3e} W/m^3 with ' +
                          f'rSil = {rSilOuter_m / 1e3:.1f} km ({rSilOuter_m / Planet.Bulk.R_m:.3f} R_{Planet.name[0]}).')
                thisHtidal_Wm3 = HtidalStart_Wm3 * multHtidal_Wm3**nProfiles
                Planet.Sil.fn_Htidal_Wm3 = GetHtidalFunc(thisHtidal_Wm3)  # Placeholder until we implement a self-consistent calc

            # Perform a check to avoid the final evaluation of SilicateLayers if we already know we won't use the outputs
            if (thisHtidal_Wm3 <= Planet.Sil.HtidalMax_Wm3 and thisphiTop_frac <= phiMax_frac):
                indsSilValidTemp, _, PsilTemp_MPa, TsilTemp_K, rSilTemp_m, rhoSilTemp_kgm3, MLayerSilTemp_kg, MAboveSilTemp_kg, gSilTemp_ms2, \
                phiSilTemp_frac, HtidalSilTemp_Wm3, kThermSilTemp_WmK, PporeTemp_MPa, rhoSilMatrixTemp_kgm3, rhoPoreTemp_kgm3, phaseSilPoreTemp \
                    = SilicateLayers(Planet, Params)
                # Record hydrosphere indices to include along with each profile
            iValid = np.append(iValid, indsSilValidTemp + Planet.Steps.iSilStart)

        # Mark all mass-matching profiles as valid
        indsSilValid = range(nProfiles)
        # Final profile is the one that violates the loop condition (or final value repeated if SilicateLayers is skipped),
        # so it is not valid -- pop it off
        iValid = iValid[:-1]
        # Now fill in values we're missing from not doing core calculations
        # We need copies of Steps.nSilMax for each possible result from the C/MR^2 search
        # to be compatible with the same infrastructure for when we have a core.
        nSilFinal = Planet.Steps.nSilMax * np.ones_like(indsSilValid).astype(np.int_)
        Ccore_kgm2 = 0

        C_kgm2 = np.zeros(Planet.Steps.nHydroMax)

    # Find contribution to axial moment of inertia C from each silicate layer
    dCfromSil_kgm2 = 8*np.pi/15 * rhoSil_kgm3 * (rSil_m[:,:-1]**5 - rSil_m[:,1:]**5)

    # Calculate C for a mantle extending up to each hydrosphere layer in turn
    Chydro_kgm2 = np.array([np.sum(dCfromH2O_kgm2[:i+1]) for i in iValid])
    Csil_kgm2 = np.array([np.sum(dCfromSil_kgm2[i,:nSilFinal[i]]) for i in indsSilValid])
    C_kgm2[iValid] = Chydro_kgm2 + Csil_kgm2 + Ccore_kgm2
    CMR2 = C_kgm2 / MR2_kgm2

    CMR2inds = [i[0] for i, valCMR2 in np.ndenumerate(CMR2)
                 if valCMR2 > Planet.Bulk.Cmeasured - Planet.Bulk.Cuncertainty
                and valCMR2 < Planet.Bulk.Cmeasured + Planet.Bulk.Cuncertainty]

    if len(CMR2inds) == 0:
        if (np.min(CMR2[iValid]) < Planet.Bulk.Cmeasured) and (np.max(CMR2[iValid]) > Planet.Bulk.Cmeasured):
            if Planet.Do.Fe_CORE:
                tweakable = 'increasing Steps.nSil or decreasing Ocean.deltaP'
            elif Planet.Do.POROUS_ROCK and not Planet.Do.FIXED_POROSITY:
                tweakable = 'increasing Steps.nPoros'
            else:
                tweakable = 'decreasing Sil.deltaHtidal'
            suggestion = f'\nTry increasing the resolution in C/MR^2 by {tweakable}.'
        elif(np.max(CMR2) > Planet.Bulk.Cmeasured):
            suggestion = '\nTry decreasing Tb_K or adjusting silicate/core composition to get lower C/MR^2 values.'
        else:
            suggestion = '\nTry adjusting properties of silicates and core to get higher C/MR^2 values.'

        msg = f'No MoI found matching C/MR^2 = {Planet.Bulk.Cmeasured:.4f}±{Planet.Bulk.Cuncertainty:.4f}.\n' + \
              f'Min: {np.min(CMR2[CMR2>0]):.4f}, Max: {np.max(CMR2):.4f}. '
        if Params.ALLOW_BROKEN_MODELS:
            if Params.DO_EXPLOREOGRAM:
                log.info(msg)
            else:
                log.error(msg + suggestion + ' Params.ALLOW_BROKEN_MODELS is True, so calculations will proceed ' +
                                             'with many values set to nan.')
            Planet.Do.VALID = False
        else:
            raise ValueError(msg + suggestion)

        nans = np.array([np.nan])
        Planet.CMR2mean = np.nan
        Planet.CMR2less = Planet.CMR2mean
        Planet.CMR2more = Planet.CMR2mean
        Planet.Sil.rhoTrade_kgm3 = nans
        Planet.Sil.Rmean_m = np.nan
        Planet.Sil.Rtrade_m = nans
        Planet.Sil.Rrange_m = np.nan
        Planet.Core.Rmean_m = np.nan
        Planet.Core.Rtrade_m = nans
        Planet.Core.Rrange_m = np.nan
        Planet.Steps.nSil = Planet.Steps.nSilMax
        iCMR2sil = 0
        if Planet.Do.NO_H2O:
            Planet.Steps.nHydro = 0
        else:
            Planet.Steps.nHydro = Planet.Steps.nOceanMax

        Psil_MPa, Tsil_K, rhoSilEOS_kgm3, gSil_ms2, phiSil_frac, kThermSil_WmK, Ppore_MPa, rhoSilMatrix_kgm3, \
        rhoPore_kgm3, HtidalSil_Wm3, MLayerSil_kg \
            = (np.zeros(Planet.Steps.nSil) for _ in range(11))
        phaseSilPore = np.zeros(Planet.Steps.nSil, dtype=np.int_)
        rSil_m = np.zeros((1, Planet.Steps.nSil+1))
        coreProps = (np.zeros(Planet.Steps.nCore) for _ in range(9))
        Planet.Sil.rhoMean_kgm3 = np.nan
        Planet.Core.rhoMean_kgm3 = np.nan
        Planet.Mtot_kg = Planet.Bulk.M_kg
        RcoreOrHtidalLine = ''

    else:
        # Find the C/MR^2 value most closely matching the measured value
        CMR2diff = np.abs(CMR2[CMR2inds] - Planet.Bulk.Cmeasured)
        # Get index of closest match in CMR2inds
        iCMR2ind = np.argmin(CMR2diff)
        # Check that we have only one best-fit option
        if np.size(iCMR2ind) > 1:
            log.warning('Multiple mass-matching profiles had the same MoI value. Using the first one. Try increasing the ' +
                     'resolution in the hydrosphere (by decreasing Ocean.deltaP) to avoid this problem.')
            iCMR2ind = iCMR2ind[0]
        # Find Planet array index corresponding to closest matching value
        iCMR2 = CMR2inds[iCMR2ind]
        # Get indices for silicate layer arrays
        if Planet.Do.Fe_CORE:
            iCMR2sil = iCMR2 - Planet.Steps.iSilStart
            CMR2indsSil = [ind - Planet.Steps.iSilStart for ind in CMR2inds]
        else:
            iCMR2sil = np.where(iCMR2 == iValid)[0][0]
            CMR2indsSil = [np.where(ind == iValid)[0][0] for ind in CMR2inds]
        # Record the best-match C/MR^2 value and neighboring values
        Planet.CMR2mean = CMR2[iCMR2]
        if Planet.Do.NO_H2O:
            log.warning('Only one profile is calculated for waterless bodies. Another method ' +
                        'of giving a trade space between total mass and MoI should replace this ' +
                        'implementation, varying surface heat flux, tidal heating, and/or porosity.')
            Planet.CMR2less = np.nan
            Planet.CMR2more = np.nan
        else:
            CMR2validSorted = np.sort(CMR2[CMR2inds])
            iSortedCMR2mean = np.where(CMR2validSorted == Planet.CMR2mean)[0][0]
            if iSortedCMR2mean == 0:
                Planet.CMR2less = Planet.CMR2mean
            else:
                Planet.CMR2less = CMR2validSorted[iSortedCMR2mean-1]
            if iSortedCMR2mean == np.size(CMR2validSorted) - 1:
                Planet.CMR2more = Planet.CMR2mean
            else:
                Planet.CMR2more = CMR2validSorted[iSortedCMR2mean+1]
        # Now we finally know how many layers there are in the hydrosphere and silicates
        Planet.Steps.nHydro = iCMR2
        Planet.Steps.nSil = nSilFinal[iCMR2sil]

        # Fill core/mantle trade arrays and set mean values consistent with MoI
        MtotSil_kg = np.sum(MLayerSil_kg[iCMR2sil,:nSilFinal[iCMR2sil]])
        VtotSil_m3 = 4/3*np.pi * (rSil_m[iCMR2sil,0]**3 - rSil_m[iCMR2sil,nSilFinal[iCMR2sil]]**3)
        Planet.Sil.rhoMean_kgm3 = MtotSil_kg / VtotSil_m3
        Planet.Sil.rhoTrade_kgm3 = np.array([np.sum(MLayerSil_kg[i,:nSilFinal[i]]) / (4/3*np.pi * (rSil_m[i,0]**3 - rSil_m[i,nSilFinal[i]-1]**3)) for i in CMR2indsSil])
        Planet.Sil.Rmean_m = Planet.r_m[iCMR2]
        Planet.Sil.Rtrade_m = Planet.r_m[CMR2inds]
        Planet.Sil.Rrange_m = np.max(Planet.Sil.Rtrade_m) - np.min(Planet.Sil.Rtrade_m)
        HtotSil_W = np.sum(HtidalSil_Wm3[iCMR2sil,:nSilFinal[iCMR2sil]] * 4/3*np.pi *
                           (rSil_m[iCMR2sil,:nSilFinal[iCMR2sil]]**3 - rSil_m[iCMR2sil,1:nSilFinal[iCMR2sil]+1]**3))
        Planet.Sil.Htidal_Wm3 = HtotSil_W / VtotSil_m3

        if Planet.Do.Fe_CORE:
            # Get indices for iron layer arrays
            iCMR2core = np.where(indsSilValid == iCMR2sil)[0][0]
            CMR2indsCore = [np.where(indsSilValid == i)[0][0] for i in CMR2indsSil]
            MtotCore_kg = np.sum(MLayerCore_kg[iCMR2core,:])
            Planet.Core.rhoMean_kgm3 = MtotCore_kg / (4/3*np.pi * rCore_m[iCMR2core,0]**3)
            Planet.Core.Rmean_m = rCore_m[iCMR2core,0]
            Planet.Core.Rtrade_m = rCore_m[CMR2indsCore,0]
            Planet.Core.Rrange_m = np.max(Planet.Core.Rtrade_m) - np.min(Planet.Core.Rtrade_m)
            RcoreOrHtidalLine = f'R_core = {Planet.Core.Rmean_m / Planet.Bulk.R_m:.2f} R, '

            # Package up core properties for returning
            coreProps = (Pcore_MPa[iCMR2core,:], Tcore_K[iCMR2core,:], rCore_m[iCMR2core,:-1],
                         gCore_ms2[iCMR2core, :], rhoCore_kgm3[iCMR2core,:], CpCore_JkgK[iCMR2core,:],
                         alphaCore_pK[iCMR2core,:], kThermCore_WmK[iCMR2core,:],
                         MLayerCore_kg[iCMR2core,:])
        else:
            if Planet.Do.POROUS_ROCK and not Planet.Do.FIXED_POROSITY:
                RcoreOrHtidalLine = f'phi_vac = {phiTop_frac[iCMR2sil]:.3f}, '
            else:
                RcoreOrHtidalLine = f'H_tidal = {Planet.Sil.Htidal_Wm3:.2e} W/m^3, '
            MtotCore_kg = 0
            Planet.Core.rhoMean_kgm3 = 0
            Planet.Core.Rmean_m = 0
            Planet.Core.Rtrade_m = np.zeros_like(Planet.Sil.Rtrade_m)
            Planet.Core.Rrange_m = 0
            coreProps = None

        Planet.Mtot_kg = np.sum(Planet.MLayer_kg[:iCMR2]) + MtotSil_kg + MtotCore_kg

    if not np.isnan(Planet.CMR2mean):
        log.info(f'Found matching MoI of {Planet.CMR2mean:.4f} ' +
                 f'(C/MR^2 = {Planet.Bulk.Cmeasured:.4f}±{Planet.Bulk.Cuncertainty:.4f}) for ' +
                 f'rho_sil = {Planet.Sil.rhoMean_kgm3:.0f} kg/m^3, ' +
                 f'R_sil = {Planet.Sil.Rmean_m / Planet.Bulk.R_m:.3f} R_{Planet.name[0]}, ' +
                 RcoreOrHtidalLine +
                 f'M_tot = {Planet.Mtot_kg/Planet.Bulk.M_kg:.4f} M_{Planet.name[0]}.')

    mantleProps = (Psil_MPa[iCMR2sil,:nSilFinal[iCMR2sil]], Tsil_K[iCMR2sil,:nSilFinal[iCMR2sil]],
                   rSil_m[iCMR2sil,:nSilFinal[iCMR2sil]], rhoSil_kgm3[iCMR2sil,:nSilFinal[iCMR2sil]],
                   gSil_ms2[iCMR2sil,:nSilFinal[iCMR2sil]], phiSil_frac[iCMR2sil,:nSilFinal[iCMR2sil]],
                   HtidalSil_Wm3[iCMR2sil,:nSilFinal[iCMR2sil]], kThermSil_WmK[iCMR2sil,:nSilFinal[iCMR2sil]],
                   PsilPore_MPa[iCMR2sil,:nSilFinal[iCMR2sil]], rhoSilMatrix_kgm3[iCMR2sil,:nSilFinal[iCMR2sil]],
                   rhoSilPore_kgm3[iCMR2sil,:nSilFinal[iCMR2sil]], MLayerSil_kg[iCMR2sil,:nSilFinal[iCMR2sil]],
                   phaseSilPore[iCMR2sil,:nSilFinal[iCMR2sil]])

    return Planet, mantleProps, coreProps
