import numpy as np
import logging as log
from collections.abc import Iterable
from Utilities.dataStructs import Constants
from Thermodynamics.FromLiterature.HydroEOS import GetPfreeze, GetTfreeze, \
    PhaseConv, GetPhaseIndices, IceEOSStruct, GetPbClath
from Thermodynamics.FromLiterature.ThermalProfiles import ConductiveTemperature
from Thermodynamics.FromLiterature.IceConduction import IceIWholeConduct, IceIConductClathLid, \
    IceIConductClathUnderplate, IceIIIConduct, IceVConduct
from Thermodynamics.FromLiterature.Convection import IceIConvect, IceIIIConvect, IceVConvect, ClathShellConvect
from Thermodynamics.FromLiterature.InnerEOS import PerplexEOSStruct

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
    Planet.PbI_MPa = GetPfreeze(Planet.Ocean.EOS, 1, Planet.Bulk.Tb_K,
                                PLower_MPa=Planet.PfreezeLower_MPa, PUpper_MPa=Planet.PfreezeUpper_MPa,
                                PRes_MPa=Planet.PfreezeRes_MPa, UNDERPLATE=(Planet.Do.BOTTOM_ICEIII or Planet.Do.BOTTOM_ICEV))
    if(Planet.Do.CLATHRATE and
            (Planet.Bulk.clathType == 'bottom' or
             Planet.Bulk.clathType == 'whole')):
        PbClath_MPa = GetPbClath(Planet.Bulk.Tb_K)
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
            raise ValueError(
                f'No valid phase transition was found for Tb_K = {Planet.Bulk.Tb_K:.3f} K for P in the range ' +
                f'[{Planet.PfreezeLower_MPa:.1f} MPa, {Planet.PfreezeUpper_MPa:.1f} MPa]. ' +
                'This likely means Tb_K is too high and the phase at the lower end of this range matches ' +
                'the phase at the upper end. Try decreasing Tb_K.')
        log.debug(f'Ice Ih transition pressure: {Planet.PbI_MPa:.3f} MPa.')

    # Now do the same for HP ices, if present, to make sure we have a possible configuration before continuing
    if Planet.Do.BOTTOM_ICEV:
        Planet.PbIII_MPa = GetPfreeze(Planet.Ocean.EOS, 3, Planet.Bulk.TbIII_K,
                   PLower_MPa=Planet.PbI_MPa, PUpper_MPa=Planet.Ocean.PHydroMax_MPa,
                   PRes_MPa=Planet.PfreezeRes_MPa, UNDERPLATE=True)
        if(Planet.PbIII_MPa <= Planet.PbI_MPa):
            raise ValueError('Ice III bottom pressure is not greater than ice I bottom pressure. ' +
                             'This likely indicates TbIII_K is too high for the corresponding Tb_K.' +
                             f'\nPbI_MPa = {Planet.PbI_MPa:.3f}' +
                             f', Tb_K = {Planet.Bulk.Tb_K:.3f}' +
                             f'\nPbIII_MPa = {Planet.PbIII_MPa:.3f}' +
                             f', TbIII_K = {Planet.Bulk.TbIII_K:.3f}')
        Planet.PbV_MPa = GetPfreeze(Planet.Ocean.EOS, 5, Planet.Bulk.TbV_K,
                                      PLower_MPa=Planet.PbIII_MPa, PUpper_MPa=Planet.Ocean.PHydroMax_MPa,
                                      PRes_MPa=Planet.PfreezeRes_MPa, UNDERPLATE=False)
        Planet.Pb_MPa = Planet.PbV_MPa
        if(Planet.PbV_MPa <= Planet.PbIII_MPa):
            raise ValueError('Ice V bottom pressure is not greater than ice III bottom pressure. ' +
                             'This likely indicates TbV_K is too high for the corresponding TbIII_K.' +
                             f'\nPbIII_MPa = {Planet.PbIII_MPa:.3f}' +
                             f', TbIII_K = {Planet.Bulk.TbIII_K:.3f}' +
                             f'\nPbV_MPa = {Planet.PbV_MPa:.3f}' +
                             f', TbV_K = {Planet.Bulk.TbV_K:.3f}')
    elif Planet.Do.BOTTOM_ICEIII:
        Planet.PbIII_MPa = GetPfreeze(Planet.Ocean.EOS, 3, Planet.Bulk.TbIII_K,
                                      PLower_MPa=Planet.PbI_MPa, PUpper_MPa=Planet.Ocean.PHydroMax_MPa,
                                      PRes_MPa=Planet.PfreezeRes_MPa, UNDERPLATE=False)
        if(Planet.PbIII_MPa <= Planet.PbI_MPa):
            raise ValueError('Ice III bottom pressure is not greater than ice I bottom pressure. ' +
                             'This likely indicates TbIII_K is too high for the corresponding Tb_K.' +
                             f'\nPbI_MPa = {Planet.PbI_MPa:.3f}' +
                             f', Tb_K = {Planet.Bulk.Tb_K:.3f}' +
                             f'\nPbIII_MPa = {Planet.PbIII_MPa:.3f}' +
                             f', TbIII_K = {Planet.Bulk.TbIII_K:.3f}')
        Planet.Pb_MPa = Planet.PbIII_MPa
    else:
        Planet.Pb_MPa = Planet.PbI_MPa

    # Now, we want to check for a convective profile. First, we need to get zb_km, so we need to suppose
    # a whole-layer conductive profile. The densities will change slightly, so we depart from self-consistency
    # here. Repeated applications of IceConvect will get closer to self-consistency.

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
            Planet = IceIConductClathLid(Planet, Params)
        elif Planet.Bulk.clathType == 'bottom':
            log.debug('Applying clathrate underplating to ice I shell.')
            Planet.phase[Planet.Steps.nIceI:Planet.Steps.nIbottom] = Constants.phaseClath
            Planet = IceIConductClathUnderplate(Planet, Params)
        elif Planet.Bulk.clathType == 'whole':
            log.debug('Applying whole-shell clathrate modeling with possible convection.')
            Planet.phase[:Planet.Steps.nIbottom] = Constants.phaseClath
            Planet = IceIWholeConduct(Planet, Params)
        else:
            raise ValueError(f'Bulk.clathType option "{Planet.Bulk.clathType}" is not supported. ' +
                             'Options are "top", "bottom", and "whole".')
    else:
        Planet = IceIWholeConduct(Planet, Params)

    log.debug('Upper ice initial conductive profile complete.')

    if not Planet.Do.NO_ICE_CONVECTION and not Planet.Bulk.clathType == 'bottom':
        # Record zb_m to see if it gets adjusted significantly
        zbOld_m = Planet.z_m[Planet.Steps.nIbottom-1] + 0.0
        # Now check for convective region and get dimensions if present
        if Planet.Do.CLATHRATE and Planet.Bulk.clathType == 'whole':
            Planet = ClathShellConvect(Planet, Params)
        else:
            Planet = IceIConvect(Planet, Params)
        # Run IceIConvect a second time if zbI_m changed by more than a set tolerance
        if(np.abs(Planet.z_m[Planet.Steps.nIbottom-1] - zbOld_m)/Planet.z_m[Planet.Steps.nIbottom-1] > Planet.Bulk.zbChangeTol_frac):
            log.debug('The bottom depth of surface ice I changed by ' +
                    f'{(Planet.z_m[Planet.Steps.nIbottom-1] - zbOld_m)/1e3:.2f} km from IceIConvect, which is greater than ' +
                    f'{Planet.Bulk.zbChangeTol_frac * 100:.0f}%. running IceIConvect a second time...')
            if Planet.Do.CLATHRATE and Planet.Bulk.clathType == 'whole':
                Planet = ClathShellConvect(Planet, Params)
            else:
                Planet = IceIConvect(Planet, Params)
    else:
        if Planet.Do.NO_ICE_CONVECTION: log.debug('NO_ICE_CONVECTION is True -- skipping ice I convection calculations.')
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

    Planet.zb_km = Planet.z_m[Planet.Steps.nSurfIce-1] / 1e3
    log.info(f'Upper ice transition pressure: {Planet.Pb_MPa:.3f} MPa, ' +
             f'thickness: {Planet.zb_km:.3f} km.')

    return Planet


def IceIIIUnderplate(Planet, Params):
    """ Conductive and convective profile calculations for ice III layers between
        the ocean/underplate ice V and surface ice I layer.

        Assigns Planet attributes:
            PbIII_MPa, all physical layer arrays
    """

    log.debug(f'Ice III bottom phase transition pressure: {Planet.PbIII_MPa:.3f} MPa ' +
             f'at TbIII_K = {Planet.Bulk.TbIII_K:.3f} K.')

    Planet = IceIIIConduct(Planet, Params)

    if not Planet.Do.NO_ICE_CONVECTION:
        # Record zbIII_m to see if it gets adjusted significantly
        zbIIIold_m = Planet.z_m[Planet.Steps.nIIIbottom-1] + 0.0
        # Now check for convective region and get dimensions if present
        Planet = IceIIIConvect(Planet, Params)
        # Run IceIIIConvect a second time if zbIII_m changed by more than a set tolerance
        if(np.abs(Planet.z_m[Planet.Steps.nIIIbottom-1] - zbIIIold_m)/Planet.z_m[Planet.Steps.nIIIbottom-1] > Planet.Bulk.zbChangeTol_frac):
            log.debug('The bottom depth of underplate ice III changed by ' +
                    f'{(Planet.z_m[Planet.Steps.nIIIbottom-1] - zbIIIold_m)/1e3:.2f} km from IceIIIConvect, which is greater than ' +
                    f'{Planet.Bulk.zbChangeTol_frac * 100:.0f}%. running IceIIIConvect a second time...')
            Planet = IceIIIConvect(Planet, Params)
    else:
        log.debug('NO_ICE_CONVECTION is True -- skipping ice III convection calculations.')

    return Planet


def IceVUnderplate(Planet, Params):
    """ Conductive and convective profile calculations for ice V layers between
        the ocean and surface ice III layer.

        Assigns Planet attributes:
            PbV_MPa, all physical layer arrays
    """
    log.debug(f'Ice V bottom phase transition pressure: {Planet.PbV_MPa:.3f} MPa ' +
                             f'at TbV_K = {Planet.Bulk.TbV_K:.3f} K.')

    Planet = IceVConduct(Planet, Params)

    if not Planet.Do.NO_ICE_CONVECTION:
        # Record zbV_m to see if it gets adjusted significantly
        zbVold_m = Planet.z_m[Planet.Steps.nSurfIce-1] + 0.0
        # Now check for convective region and get dimensions if present
        Planet = IceVConvect(Planet, Params)
        # Run IceVConvect a second time if zbV_m changed by more than a set tolerance
        if(np.abs(Planet.z_m[Planet.Steps.nSurfIce-1] - zbVold_m)/Planet.z_m[Planet.Steps.nSurfIce-1] > Planet.Bulk.zbChangeTol_frac):
            log.debug('The bottom depth of underplate ice V changed by ' +
                     f'{(Planet.z_m[Planet.Steps.nSurfIce-1] - zbVold_m)/1e3:.2f} km from IceVConvect, which is greater than ' +
                     f'{Planet.Bulk.zbChangeTol_frac * 100:.0f}%. running IceVConvect a second time...')
            Planet = IceVConvect(Planet, Params)
    else:
        log.debug('NO_ICE_CONVECTION is True -- skipping ice V convection calculations.')

    return Planet


def OceanLayers(Planet, Params):
    """ Geophysical and thermodynamic calculations for ocean layer
        Calculates state variables of the layer with each pressure step

        Assigns Planet attributes:
            phase, r_m, z_m, g_ms2, T_K, P_MPa, rho_kgm3, Cp_JkgK, alpha_pK, MLayer_kg
    """
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
    if(Planet.Ocean.PHydroMax_MPa > Constants.PminHPices_MPa):
        GetOceanHPIceEOS(Planet, Params, POcean_MPa)

    # Do initial ocean step separately in order to catch potential Melosh layer--
    # see Melosh et al. (2004): https://doi.org/10.1016/j.icarus.2003.11.026
    log.debug(f'il: {Planet.Steps.nSurfIce:d}; P_MPa: {POcean_MPa[0]:.3f}; ' +
              f'T_K: {TOcean_K[0]:.3f}; phase: {Planet.phase[Planet.Steps.nSurfIce]:d}')
    rhoOcean_kgm3[0] = Planet.Ocean.EOS.fn_rho_kgm3(POcean_MPa[0], TOcean_K[0], grid=False)
    CpOcean_JkgK[0] = Planet.Ocean.EOS.fn_Cp_JkgK(POcean_MPa[0], TOcean_K[0], grid=False)
    alphaOcean_pK[0] = Planet.Ocean.EOS.fn_alpha_pK(POcean_MPa[0], TOcean_K[0], grid=False)
    kThermOcean_WmK[0] = Planet.Ocean.EOS.fn_kTherm_WmK(POcean_MPa[0], TOcean_K[0], grid=False)
    if alphaOcean_pK[0] < 0:
        log.info(f'Thermal expansivity alpha at the ice-ocean interface is negative. Modeling Melosh et al. conductive layer.')
        # Layer should be thin, so we just use a fixed dT/dz value
        dTdz = Planet.Ocean.QfromMantle_W / (4*np.pi * Planet.r_m[Planet.Steps.nSurfIce]**2) / kThermOcean_WmK[0]
        i = 0
        while alphaOcean_pK[i] < 0:
            i += 1
            dz = Planet.Ocean.deltaP*1e6 / Planet.g_ms2[Planet.Steps.nSurfIce] / rhoOcean_kgm3[i-1]
            # Model temperature as linear and conductive in this layer
            TOcean_K[i] = TOcean_K[i-1] + dTdz * dz
            rhoOcean_kgm3[i] = Planet.Ocean.EOS.fn_rho_kgm3(POcean_MPa[i], TOcean_K[i], grid=False)
            CpOcean_JkgK[i] = Planet.Ocean.EOS.fn_Cp_JkgK(POcean_MPa[i], TOcean_K[i], grid=False)
            alphaOcean_pK[i] = Planet.Ocean.EOS.fn_alpha_pK(POcean_MPa[i], TOcean_K[i], grid=False)
            kThermOcean_WmK[i] = Planet.Ocean.EOS.fn_kTherm_WmK(POcean_MPa[i], TOcean_K[i], grid=False)
        iStart = i
    else:
        # Now use the present layer's properties to calculate an adiabatic thermal profile for layers below
        TOcean_K[1] = TOcean_K[0] + alphaOcean_pK[0] * TOcean_K[0] / \
                          CpOcean_JkgK[0] / rhoOcean_kgm3[0] * Planet.Ocean.deltaP*1e6
        iStart = 1

    for i in range(iStart, Planet.Steps.nOceanMax):
        Planet.phase[Planet.Steps.nSurfIce+i] = Planet.Ocean.EOS.fn_phase(POcean_MPa[i], TOcean_K[i]).astype(np.int_)
        log.debug(f'il: {Planet.Steps.nSurfIce+i:d}; P_MPa: {POcean_MPa[i]:.3f}; ' +
                  f'T_K: {TOcean_K[i]:.3f}; phase: {Planet.phase[Planet.Steps.nSurfIce+i]:d}')
        if Planet.phase[Planet.Steps.nSurfIce+i] == 0:
            # Liquid water layers -- get fluid properties for the present layer but with the
            # overlaying layer's temperature
            rhoOcean_kgm3[i] = Planet.Ocean.EOS.fn_rho_kgm3(POcean_MPa[i], TOcean_K[i], grid=False)
            CpOcean_JkgK[i] = Planet.Ocean.EOS.fn_Cp_JkgK(POcean_MPa[i], TOcean_K[i], grid=False)
            alphaOcean_pK[i] = Planet.Ocean.EOS.fn_alpha_pK(POcean_MPa[i], TOcean_K[i], grid=False)
            kThermOcean_WmK[i] = Planet.Ocean.EOS.fn_kTherm_WmK(POcean_MPa[i], TOcean_K[i], grid=False)
            # Now use the present layer's properties to calculate an adiabatic thermal profile for layers below
            TOcean_K[i+1] = TOcean_K[i] + alphaOcean_pK[i] * TOcean_K[i] / \
                            CpOcean_JkgK[i] / rhoOcean_kgm3[i] * Planet.Ocean.deltaP*1e6
        else:
            # Undersea high-pressure ices -- we use GetTfreeze here to propagate the layer temperatures.
            # This is based on an assumption that the undersea HP ices are vigorously mixed by
            # two-phase convection, such that each layer is in local equilibrium with the liquid,
            # meaning each layer's temperature is equal to the melting temperature.
            thisPhase = PhaseConv(Planet.phase[Planet.Steps.nSurfIce+i])
            rhoOcean_kgm3[i] = Planet.Ocean.iceEOS[thisPhase].fn_rho_kgm3(POcean_MPa[i], TOcean_K[i], grid=False)
            CpOcean_JkgK[i] = Planet.Ocean.iceEOS[thisPhase].fn_Cp_JkgK(POcean_MPa[i], TOcean_K[i], grid=False)
            alphaOcean_pK[i] = Planet.Ocean.iceEOS[thisPhase].fn_alpha_pK(POcean_MPa[i], TOcean_K[i], grid=False)
            kThermOcean_WmK[i] = Planet.Ocean.iceEOS[thisPhase].fn_kTherm_WmK(POcean_MPa[i], TOcean_K[i], grid=False)
            TOcean_K[i+1] = GetTfreeze(Planet.Ocean.EOS, POcean_MPa[i], TOcean_K[i])

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

    log.info('Ocean layers complete.')

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
    # Remove this if block when a faster phase calculation is implemented!
    if Planet.Ocean.comp == 'MgSO4':
        # Just load all HP ice phases in case we need them. This part is way faster than Margules phase calcs
        Planet.Ocean.iceEOS['II'] = IceEOSStruct(POceanHPices_MPa, TOceanHPices_K, 'II')
        Planet.Ocean.iceEOS['III'] = IceEOSStruct(POceanHPices_MPa, TOceanHPices_K, 'III')
        Planet.Ocean.iceEOS['V'] = IceEOSStruct(POceanHPices_MPa, TOceanHPices_K, 'V')
        Planet.Ocean.iceEOS['VI'] = IceEOSStruct(POceanHPices_MPa, TOceanHPices_K, 'VI')
    else:
        # Get phase of each P,T combination
        expandPhases = Planet.Ocean.EOS.fn_phase(PHPicesLin_MPa, THPicesLin_K)
        # Check if any of them are not liquid
        if np.any(expandPhases != 0):
            _, _, indsIceII, indsIceIII, indsIceV, indsIceVI, _, _, _ = GetPhaseIndices(expandPhases)

            if(np.size(indsIceII) != 0):
                log.debug('Loading ice II EOS functions for ocean layers...')
                PiceIImin_MPa = PHPicesLin_MPa[indsIceII[0]]
                PiceIImax_MPa = PHPicesLin_MPa[indsIceII[-1]]
                TiceIImin_K = np.min(THPicesLin_K[indsIceII])
                TiceIImax_K = np.max(THPicesLin_K[indsIceII])
                Planet.Ocean.iceEOS['II'] = IceEOSStruct(np.linspace(PiceIImin_MPa, PiceIImax_MPa, Planet.Steps.nPsHP),
                                                         np.linspace(TiceIImin_K, TiceIImax_K, Planet.Steps.nTsHP), 'II')
            if(np.size(indsIceIII) != 0):
                log.debug('Loading ice III EOS functions for ocean layers...')
                PiceIIImin_MPa = PHPicesLin_MPa[indsIceIII[0]]
                PiceIIImax_MPa = PHPicesLin_MPa[indsIceIII[-1]]
                TiceIIImin_K = np.min(THPicesLin_K[indsIceIII])
                TiceIIImax_K = np.max(THPicesLin_K[indsIceIII])
                Planet.Ocean.iceEOS['III'] = IceEOSStruct(np.linspace(PiceIIImin_MPa, PiceIIImax_MPa, Planet.Steps.nPsHP),
                                                          np.linspace(TiceIIImin_K, TiceIIImax_K, Planet.Steps.nTsHP), 'III')
            if(np.size(indsIceV) != 0):
                log.debug('Loading ice V EOS functions for ocean layers...')
                PiceVmin_MPa = PHPicesLin_MPa[indsIceV[0]]
                PiceVmax_MPa = PHPicesLin_MPa[indsIceV[-1]]
                TiceVmin_K = np.min(THPicesLin_K[indsIceV])
                TiceVmax_K = np.max(THPicesLin_K[indsIceV])
                Planet.Ocean.iceEOS['V'] = IceEOSStruct(np.linspace(PiceVmin_MPa, PiceVmax_MPa, Planet.Steps.nPsHP),
                                                        np.linspace(TiceVmin_K, TiceVmax_K, Planet.Steps.nTsHP), 'V')
            if(np.size(indsIceVI) != 0):
                log.debug('Loading ice VI EOS functions for ocean layers...')
                PiceVImin_MPa = PHPicesLin_MPa[indsIceVI[0]]
                PiceVImax_MPa = PHPicesLin_MPa[indsIceVI[-1]]
                TiceVImin_K = np.min(THPicesLin_K[indsIceVI])
                TiceVImax_K = np.max(THPicesLin_K[indsIceVI])
                Planet.Ocean.iceEOS['VI'] = IceEOSStruct(np.linspace(PiceVImin_MPa, PiceVImax_MPa, Planet.Steps.nPsHP),
                                                         np.linspace(TiceVImin_K, TiceVImax_K, Planet.Steps.nTsHP), 'VI')
        else:
            log.debug('No high-pressure ices found in ocean layers.')

    return Planet, Params


def InnerLayers(Planet, Params):
    """ Geophysical and thermodynamic calculations for silicate and core layers
        Calculates state variables of the layer with each pressure step

        Assigns Planet attributes:
            Steps.nTotal, all layer arrays
    """
    if Planet.Do.CONSTANT_INNER_DENSITY or Params.SKIP_INNER:
        # It will be better to calculate the layers always using the EOS, but no-core models
        # require more inputs to be physically reasonable.
        if not Planet.Do.CONSTANT_INNER_DENSITY:
            # Force this flag on in case Params.SKIP_INNER is True and this flag is not,
            # to avoid problems in SilicateLayers:
            Planet.Do.CONSTANT_INNER_DENSITY = True
            log.warning('Do.CONSTANT_INNER_DENSITY forced ON based on implementation in SilicateLayers.')
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
    Planet.r_m = np.concatenate((Planet.r_m[:Planet.Steps.nHydro], extend))
    Planet.phase = np.concatenate((Planet.phase[:Planet.Steps.nHydro], extend.astype(np.int_)))
    Planet.rho_kgm3 = np.concatenate((Planet.rho_kgm3[:Planet.Steps.nHydro], extend))
    Planet.Cp_JkgK = np.concatenate((Planet.Cp_JkgK[:Planet.Steps.nHydro], extend))
    Planet.alpha_pK = np.concatenate((Planet.alpha_pK[:Planet.Steps.nHydro], extend))
    Planet.kTherm_WmK = np.concatenate((Planet.kTherm_WmK[:Planet.Steps.nHydro], extend))
    Planet.g_ms2 = np.concatenate((Planet.g_ms2[:Planet.Steps.nHydro], extend))
    Planet.phi_frac = np.concatenate((Planet.phi_frac[:Planet.Steps.nHydro], extend))
    Planet.z_m = np.concatenate((Planet.z_m[:Planet.Steps.nHydro], extend))

    # Assign phase values for silicates and core
    Planet.phase[Planet.Steps.nHydro:Planet.Steps.nHydro + Planet.Steps.nSil] = Constants.phaseSil
    Planet.phase[Planet.Steps.nHydro + Planet.Steps.nSil:Planet.Steps.nTotal] = Constants.phaseFe

    # Unpack results from MoI calculations
    iOS = Planet.Steps.nHydro
    iSC = Planet.Steps.nHydro + Planet.Steps.nSil
    Planet.P_MPa[iOS:iSC], Planet.T_K[iOS:iSC], Planet.r_m[iOS:iSC], Planet.rho_kgm3[iOS:iSC], \
    Planet.g_ms2[iOS:iSC], Planet.phi_frac[iOS:iSC] = mantleProps

    iCC = Planet.Steps.nTotal
    if Planet.Do.Fe_CORE:
        # Unpack results from MoI calculations
        Planet.P_MPa[iSC:iCC], Planet.T_K[iSC:iCC], Planet.r_m[iSC:iCC], Planet.g_ms2[iSC:iCC], Planet.rho_kgm3[iSC:iCC], \
        Planet.Cp_JkgK[iSC:iCC], Planet.alpha_pK[iSC:iCC] = coreProps

    Planet.z_m[iOS:iCC] = Planet.Bulk.R_m - Planet.r_m[iOS:iCC]

    # Check for any negative temperature gradient (indicates non-equilibrium conditions)
    gradTneg = np.where(np.diff(Planet.T_K) < 0)
    if np.any(gradTneg):
        log.warning(f'Negative temperature gradient starting at index {gradTneg[0][0]:d}. This indicates that ' +
                     'internal heating parameters Qrad_Wkg and/or Htidal_Wm3 are set too high to be consistent ' +
                     'with the heat flux into the ocean. The current configuration represents a ' +
                     'non-equilibrium state.')

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
        # Find core bulk density based on assumed sulfide content (Eq 10 of Vance et al., 2014)
        rhoCore_kgm3 = Planet.Core.rhoFeS_kgm3 * Planet.Core.rhoFe_kgm3 / \
            (Planet.Core.xFeS * (Planet.Core.rhoFe_kgm3 - Planet.Core.rhoFeS_kgm3) + Planet.Core.rhoFeS_kgm3)
        # Calculate core volume for a silicate layer with outer radius equal to bottom of each hydrosphere layer
        # and inner radius equal to the core radius
        VCore_m3 = np.array([(Planet.Bulk.M_kg - MAbove_kg[i] - VsilSphere_m3[i-Planet.Steps.iSilStart]*
                             Planet.Sil.rhoSilWithCore_kgm3) / (rhoCore_kgm3 - Planet.Sil.rhoSilWithCore_kgm3)
                             for i in range(Planet.Steps.iSilStart, nHydroActual-1)])
        # Find values for which the silicate radius is too large
        nTooBig = next((i[0] for i, val in np.ndenumerate(VCore_m3) if val>0))
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
        if(Planet.Do.NO_H2O or np.max(CMR2) > Planet.Bulk.Cmeasured):
            suggestion = 'Try adjusting properties of silicates and core to get C/MR^2 values in range.'
        else:
            suggestion = 'Try increasing PHydroMax_MPa to get lower C/MR^2 values.'
        raise ValueError(f'No MoI found matching C/MR^2 = {Planet.Bulk.Cmeasured:.3f}±{Planet.Bulk.Cuncertainty:.3f}.\n' +
                         f'Min: {np.min(CMR2[CMR2>0]):.3f}, Max: {np.max(CMR2):.3f}.\n' + suggestion)

    # Find the C/MR^2 value most closely matching the measured value
    CMR2diff = np.abs(CMR2[CMR2inds] - Planet.Bulk.Cmeasured)
    # Get index of closest match in CMR2inds
    iCMR2inds = np.argmin(CMR2diff)
    # Find Planet array index corresponding to closest matching value
    iCMR2 = CMR2inds[iCMR2inds]
    iCMR2inner = iCMR2 - Planet.Steps.iSilStart - nTooBig
    CMR2indsInner = [ind - Planet.Steps.iSilStart - nTooBig for ind in CMR2inds]
    # Record the best-match C/MR^2 value
    Planet.CMR2mean = CMR2[iCMR2]
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

    if not Params.SKIP_INNER:
        # Load in Perple_X table for silicate properties
        log.debug('Loading silicate Perple_X table...')
        Planet.Sil.EOS = PerplexEOSStruct(Planet.Sil.mantleEOS, EOSinterpMethod=Params.interpMethod,
                                      kThermConst_WmK=Planet.Sil.kTherm_WmK)
        # Evaluate the silicate EOS for each layer
        indsSilValid, nProfiles, Psil_MPa, Tsil_K, rSil_m, rhoSilEOS_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, phiSil_frac \
            = SilicateLayers(Planet, Params)
        nSilTooBig = nProfiles - np.size(indsSilValid)

        # Fill core/mantle trade arrays and set mean values consistent with MoI
        MtotSil_kg = np.sum(MLayerSil_kg)
        Planet.Sil.rhoMean_kgm3 = MtotSil_kg / (4/3*np.pi * (rSil_m[0,0]**3 - rSil_m[0,-1]**3))

        if Planet.Do.Fe_CORE:
            # Load in Perple_X table for core properties
            log.debug('Loading core Perple_X table...')
            Planet.Core.EOS = PerplexEOSStruct(Planet.Core.coreEOS, EOSinterpMethod=Params.interpMethod, Fe_EOS=True,
                                      kThermConst_WmK=Planet.Core.kTherm_WmK)
            # Evaluate the core EOS for each layer
            _, Pcore_MPa, Tcore_K, rCoreEOS_m, rhoCore_kgm3, MLayerCore_kg, gCore_ms2, CpCore_JkgK, alphaCore_pK = \
                IronCoreLayers(Planet, Params,
                               indsSilValid, nSilTooBig, nProfiles, Psil_MPa, Tsil_K, rSil_m, MAboveSil_kg, gSil_ms2)

            MtotCore_kg = np.sum(MLayerCore_kg)
            Planet.Core.rhoMean_kgm3 = MtotCore_kg / VCore_m3[iCMR2inner]

            coreProps = (Pcore_MPa, Tcore_K, rCoreEOS_m[0,:-1], gCore_ms2, rhoCore_kgm3, CpCore_JkgK, alphaCore_pK)
        else:
            MtotCore_kg = 0
            Planet.Core.rhoMean_kgm3 = 0
            Planet.Core.Rtrade_m = np.zeros_like(Planet.Sil.Rtrade_m)
            Planet.Core.Rrange_m = 0
            coreProps = ()

        Mtot_kg = np.sum(Planet.MLayer_kg[:iCMR2]) + MtotSil_kg + MtotCore_kg
        log.info(f'Found matching MoI of {Planet.CMR2mean:.3f} ' +
                 f'(C/MR^2 = {Planet.Bulk.Cmeasured:.3f}±{Planet.Bulk.Cuncertainty:.3f}) for ' +
                 f'R_sil = {Planet.Sil.Rmean_m / Planet.Bulk.R_m:.2f} R, ' +
                 f'R_core = {Planet.Core.Rmean_m / Planet.Bulk.R_m:.2f} R, ' +
                 f'rho_sil (found) = {rhoSil_kgm3[iCMR2inner]:.0f} kg/m^3, ' +
                 f'rho_sil (actual) = {Planet.Sil.rhoMean_kgm3:.0f} kg/m^3, ' +
                 f'M_tot = {Mtot_kg/Planet.Bulk.M_kg:.4f} M_{Planet.name[0]}.')
        log.warning('Because silicate and core properties were determined from the EOS after finding their ' +
                    'sizes by assuming constant densities, the body mass may not match the measured value.')

    else:
        Psil_MPa, Tsil_K, rhoSilEOS_kgm3, gSil_ms2, phiSil_frac = (np.zeros(Planet.Steps.nSil) for _ in range(5))
        rSil_m = np.zeros((1, Planet.Steps.nSil+1))
        coreProps = (np.zeros(Planet.Steps.nCore) for _ in range(7))
        Planet.Sil.rhoMean_kgm3 = 0
        Planet.Core.rhoMean_kgm3 = 0

    mantleProps = (Psil_MPa, Tsil_K, rSil_m[0,:-1], rhoSilEOS_kgm3, gSil_ms2, phiSil_frac)

    return Planet, mantleProps, coreProps


def CalcMoIWithEOS(Planet, Params):
    """ Find the relative sizes of silicate, core, and hydrosphere layers that are
        consistent with the measured moment of inertia, based on calculated hydrosphere
        properties and EOS data for assumed mantle and core compositions output by Perple_X.

        Assigns Planet attributes:
            CMR2mean, Sil.RsilMean_m, Sil.RsilRange_m, Core.RFeMean_m, Core.RFeRange_m, Steps.nHydro, Steps.nSil,
            all layer arrays
    """
    log.debug('Finding MoI consistent with measured value...')
    # Get MR^2 -- we will need to divide each C by this later.
    MR2_kgm2 = Planet.Bulk.M_kg * Planet.Bulk.R_m**2

    # Load in Perple_X table for silicate properties
    log.debug('Loading silicate Perple_X table...')
    Planet.Sil.EOS = PerplexEOSStruct(Planet.Sil.mantleEOS, EOSinterpMethod=Params.interpMethod,
                                      kThermConst_WmK=Planet.Sil.kTherm_WmK)
    # Propagate the silicate EOS from each hydrosphere layer to the center of the body
    indsSilValid, nProfiles, Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, phiSil_frac \
        = SilicateLayers(Planet, Params)
    nSilTooBig = nProfiles - np.size(indsSilValid)
    iValid = np.array([Planet.Steps.iSilStart + i for i in indsSilValid]).astype(np.int_)

    # Find contribution to axial moment of inertia C from each ocean layer
    dCfromH2O_kgm2 = 8*np.pi/15 * Planet.rho_kgm3[:-1] * (Planet.r_m[:-1]**5 - Planet.r_m[1:]**5)
    # Same for silicate layers
    dCfromSil_kgm2 = 8*np.pi/15 * rhoSil_kgm3 * (rSil_m[:,:-1]**5 - rSil_m[:,1:]**5)

    if Planet.Do.Fe_CORE:
        # Load in Perple_X table for core properties
        log.debug('Loading core Perple_X table...')
        Planet.Core.EOS = PerplexEOSStruct(Planet.Core.coreEOS, EOSinterpMethod=Params.interpMethod, Fe_EOS=True,
                                      kThermConst_WmK=Planet.Core.kTherm_WmK)
        # Propagate the core EOS from each silicate layer at the max core radius to the center of the body
        nSilFinal, Pcore_MPa, Tcore_K, rCore_m, rhoCore_kgm3, MLayerCore_kg, gCore_ms2, CpCore_JkgK, alphaCore_pK \
            = IronCoreLayers(Planet, Params,
                           indsSilValid, nSilTooBig, nProfiles, Psil_MPa, Tsil_K, rSil_m, MAboveSil_kg, gSil_ms2)

        dCfromCore_kgm2 = 8*np.pi/15 * rhoCore_kgm3 * (rCore_m[:,:-1]**5 - rCore_m[:,1:]**5)
    else:
        # To implement silicates-only here, we need some other variable besides silicate size to tweak
        # that can permit us to match both the body mass and MoI. Composition or porosity would each work,
        # but there's no example to draw from for tweaking these to match the MoI yet.
        coreProps = None
        nSilFinal = Planet.Steps.nSilMax
        raise ValueError('Handling for MoI matching with silicate EOS and no core is not implemented yet. ' +
                         'Set Planet.Do.CONSTANT_INNER_DENSITY to True and run again.')

    # Calculate C for a mantle extending up to each hydrosphere layer in turn
    C_kgm2 = np.zeros(nProfiles + Planet.Steps.iSilStart)
    Chydro_kgm2 = np.array([np.sum(dCfromH2O_kgm2[:i + 1]) for i in iValid])
    Csil_kgm2 = np.array([np.sum(dCfromSil_kgm2[i,:nSilFinal[i]]) for i in indsSilValid])
    Ccore_kgm2 = np.array([np.sum(dCfromCore_kgm2[i,:]) for i in range(nProfiles - nSilTooBig)])
    C_kgm2[iValid] = Chydro_kgm2 + Csil_kgm2 + Ccore_kgm2
    CMR2 = C_kgm2 / MR2_kgm2

    CMR2inds = [i[0] for i, valCMR2 in np.ndenumerate(CMR2)
                 if valCMR2 > Planet.Bulk.Cmeasured - Planet.Bulk.Cuncertainty
                and valCMR2 < Planet.Bulk.Cmeasured + Planet.Bulk.Cuncertainty]

    if len(CMR2inds) == 0:
        if(np.max(CMR2) > Planet.Bulk.Cmeasured):
            suggestion = 'Try increasing PHydroMax_MPa to get lower C/MR^2 values.'
        else:
            suggestion = 'Try adjusting properties of silicates and core to get higher C/MR^2 values.'
        raise ValueError(f'No MoI found matching C/MR^2 = {Planet.Bulk.Cmeasured:.3f}±{Planet.Bulk.Cuncertainty:.3f}.\n' +
                         f'Min: {np.min(CMR2[CMR2>0]):.3f}, Max: {np.max(CMR2):.3f}.\n ' + suggestion)

    # Find the C/MR^2 value most closely matching the measured value
    CMR2diff = np.abs(CMR2[CMR2inds] - Planet.Bulk.Cmeasured)
    # Get index of closest match in CMR2inds
    iCMR2inds = np.argmin(CMR2diff)
    # Find Planet array index corresponding to closest matching value
    iCMR2 = CMR2inds[iCMR2inds]
    # Get indices for inner layer arrays
    iCMR2sil = iCMR2 - Planet.Steps.iSilStart
    iCMR2core = np.argwhere(indsSilValid == iCMR2sil)[0][0]
    CMR2indsSil = [ind - Planet.Steps.iSilStart for ind in CMR2inds]
    CMR2indsCore = [np.argwhere(indsSilValid == i)[0][0] for i in CMR2indsSil]
    # Record the best-match C/MR^2 value
    Planet.CMR2mean = CMR2[iCMR2]
    # Now we finally know how many layers there are in the hydrosphere and silicates
    Planet.Steps.nHydro = iCMR2
    Planet.Steps.nSil = nSilFinal[iCMR2sil]

    # Fill core/mantle trade arrays and set mean values consistent with MoI
    MtotSil_kg = np.sum(MLayerSil_kg[iCMR2sil,:nSilFinal[iCMR2sil]])
    Planet.Sil.rhoMean_kgm3 = MtotSil_kg / (4/3*np.pi * (rSil_m[iCMR2sil,0]**3 - rSil_m[iCMR2sil,nSilFinal[iCMR2sil]-1]**3))
    Planet.Sil.rhoTrade_kgm3 = np.array([np.sum(MLayerSil_kg[i,:nSilFinal[i]]) / (4/3*np.pi * (rSil_m[i,0]**3 - rSil_m[i,nSilFinal[i]-1]**3)) for i in CMR2indsSil])
    Planet.Sil.Rmean_m = Planet.r_m[iCMR2]
    Planet.Sil.Rtrade_m = Planet.r_m[CMR2inds]
    Planet.Sil.Rrange_m = np.max(Planet.Sil.Rtrade_m) - np.min(Planet.Sil.Rtrade_m)
    if Planet.Do.Fe_CORE:
        MtotCore_kg = np.sum(MLayerCore_kg[iCMR2core,:])
        Planet.Core.rhoMean_kgm3 = MtotCore_kg / (4/3*np.pi * rCore_m[iCMR2core,0]**3)
        Planet.Core.Rmean_m = rCore_m[iCMR2core,0]
        Planet.Core.Rtrade_m = rCore_m[CMR2indsCore,0]
        Planet.Core.Rrange_m = np.max(Planet.Core.Rtrade_m) - np.min(Planet.Core.Rtrade_m)

        # Package up core properties for returning
        coreProps = (Pcore_MPa[iCMR2core,:], Tcore_K[iCMR2core,:], rCore_m[iCMR2core,:-1],
                     gCore_ms2[iCMR2core, :], rhoCore_kgm3[iCMR2core,:], CpCore_JkgK[iCMR2core,:],
                     alphaCore_pK[iCMR2core,:])
    else:
        MtotCore_kg = 0
        Planet.Core.rhoMean_kgm3 = 0
        Planet.Core.Rmean_m = 0
        Planet.Core.Rtrade_m = 0
        Planet.Core.Rrange_m = 0
        coreProps = None

    Mtot_kg = np.sum(Planet.MLayer_kg[:iCMR2]) + MtotSil_kg + MtotCore_kg
    log.info(f'Found matching MoI of {Planet.CMR2mean:.3f} ' +
             f'(C/MR^2 = {Planet.Bulk.Cmeasured:.3f}±{Planet.Bulk.Cuncertainty:.3f}) for ' +
             f'rho_sil = {Planet.Sil.rhoMean_kgm3:.0f} kg/m^3, ' +
             f'R_sil = {Planet.Sil.Rmean_m / Planet.Bulk.R_m:.2f} R, ' +
             f'R_core = {Planet.Core.Rmean_m / Planet.Bulk.R_m:.2f} R, ' +
             f'M_tot = {Mtot_kg/Planet.Bulk.M_kg:.4f} M_{Planet.name[0]}.')

    mantleProps = (Psil_MPa[iCMR2sil,:nSilFinal[iCMR2sil]], Tsil_K[iCMR2sil,:nSilFinal[iCMR2sil]],
                   rSil_m[iCMR2sil,:nSilFinal[iCMR2sil]], rhoSil_kgm3[iCMR2sil,:nSilFinal[iCMR2sil]],
                   gSil_ms2[iCMR2sil,:nSilFinal[iCMR2sil]], phiSil_frac[iCMR2sil,:nSilFinal[iCMR2sil]])

    return Planet, mantleProps, coreProps


def SilicateLayers(Planet, Params):
    """ Determines properties of silicate layers based on input Perple_X table
        and seafloor properties.

        Returns:
            nSilTooBig (int): Number of silicate profiles that have a mass that exceeds the body mass
            nProfiles (int): Number of silicate profiles considered
            Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, phiSil_frac
                (float, shape nHydroMax-2): State variables needed to determine underlying core properties
                 to proceed with MoI calculations.
    """
    if Planet.Do.CONSTANT_INNER_DENSITY or Planet.Do.NO_H2O:
        # If CONSTANT_INNER_DENSITY is True, we have already done the C/MR^2 calculations
        # and now we are just evaluating the EOS for the winning silicate layer set
        # Similarly, if Do.NO_H2O is True, we have 0 ice layers so there is just 1 profile to run
        nProfiles = 1
        if Planet.Do.NO_H2O:
            Planet.Steps.iSilStart = 0
            Planet.Steps.nHydro = 0
        profRange = [Planet.Steps.nHydro - Planet.Steps.iSilStart]
    else:
        nProfiles = Planet.Steps.nSurfIce + Planet.Steps.nOceanMax - Planet.Steps.iSilStart - 1
        profRange = range(nProfiles)
    # Initialize output arrays and working arrays
    Psil_MPa, Tsil_K, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg, gSil_ms2, phiSil_frac, \
        = (np.zeros((nProfiles, Planet.Steps.nSilMax)) for _ in range(8))
    # Check if we set the core radius to 0 or a found C/MR^2 value (for constant-density approach)
    if Planet.Core.Rmean_m is not None:
        rSilEnd_m = Planet.Core.Rmean_m
    else:
        rSilEnd_m = 0

    rSil_m = np.array([np.linspace(Planet.r_m[i+Planet.Steps.iSilStart], rSilEnd_m, Planet.Steps.nSilMax+1) for i in profRange])
    Psil_MPa[:,0] = [Planet.P_MPa[i+Planet.Steps.iSilStart] for i in profRange]
    Tsil_K[:,0] = [Planet.T_K[i+Planet.Steps.iSilStart] for i in profRange]
    rhoSil_kgm3[:,0] = Planet.Sil.EOS.fn_rho_kgm3(Psil_MPa[:nProfiles,0], Tsil_K[:nProfiles,0], grid=False)
    kThermSil_WmK[:,0] = Planet.Sil.EOS.fn_kTherm_WmK(Psil_MPa[:nProfiles,0], Tsil_K[:nProfiles,0], grid=False)
    MLayerSil_kg[:,0] = [rhoSil_kgm3[i,0] * 4/3*np.pi*(rSil_m[i,0]**3 - rSil_m[i,1]**3) for i in range(nProfiles)]
    gSil_ms2[:,0] = [Planet.g_ms2[i+Planet.Steps.iSilStart] for i in profRange]
    if Planet.Do.POROUS_ROCK:
        log.warning('POROUS_ROCK not implemented yet. Only the top silicate layer will have porosity set.')
        phiSil_frac[:,0] = [Planet.Sil.phiRockMax_frac for _ in profRange]

    MHydro_kg = np.array([np.sum(Planet.MLayer_kg[:i]) for i in range(Planet.Steps.iSilStart, Planet.Steps.iSilStart + nProfiles)])
    # Initialize MAbove_kg to 0th silicate layer, so that the hydrosphere mass is equal to the mass above the silicates.
    MAboveSil_kg[:,0] = MHydro_kg + 0.0
    # Initialize qTop_WmK, the heat flux leaving the top of each layer.
    # For the top layer, this is equal to the heat flux warming the hydrosphere.
    qTop_Wm2 = Planet.Ocean.QfromMantle_W / 4/np.pi/rSil_m[:,0]**2

    log.debug('Propagating silicate EOS for each possible mantle size...')
    for j in range(1, Planet.Steps.nSilMax):
        MAboveSil_kg[:,j] = MAboveSil_kg[:,j-1] + MLayerSil_kg[:,j-1]
        Psil_MPa[:,j] = Psil_MPa[:,j-1] + 1e-6 * MLayerSil_kg[:,j-1] * gSil_ms2[:,j-1] / (4*np.pi*rSil_m[:,j]**2)
        Tsil_K[:,j], qTop_Wm2 = ConductiveTemperature(Tsil_K[:,j-1], rSil_m[:,j-1], rSil_m[:,j],
                    kThermSil_WmK[:,j-1], rhoSil_kgm3[:,j-1], Planet.Sil.Qrad_Wkg, Planet.Sil.Htidal_Wm3,
                    qTop_Wm2)
        rhoSil_kgm3[:,j] = Planet.Sil.EOS.fn_rho_kgm3(Psil_MPa[:nProfiles,j], Tsil_K[:nProfiles,j], grid=False)
        kThermSil_WmK[:,j] = Planet.Sil.EOS.fn_kTherm_WmK(Psil_MPa[:nProfiles,j], Tsil_K[:nProfiles,j], grid=False)
        MLayerSil_kg[:,j] = rhoSil_kgm3[:,j] * 4/3*np.pi*(rSil_m[:,j]**3 - rSil_m[:,j+1]**3)
        # Calculate gravity using absolute values, as we will use MAboveSil to check for exceeding body mass later.
        gSil_ms2[:,j] = Constants.G * np.abs(Planet.Bulk.M_kg - MAboveSil_kg[:,j]) / rSil_m[:,j]**2

    if np.any(Tsil_K < 0):
        raise RuntimeError('Negative temperatures encountered in silicates. This likely indicates Qrad_Wkg + Htidal_Wm3 ' +
                           'is too high to be consistent with the heat flow through the ice shell.')

    if Planet.Do.CONSTANT_INNER_DENSITY:
        # Include all indices for later calculations if we already found the desired C/MR^2 match
        indsSilValid = profRange
    else:
        # Get total mass for each possible silicate layer size
        Mtot_kg = MLayerSil_kg[:,-1] + MAboveSil_kg[:,-1]
        # Find silicate radii for which the total mass is too high so we can exclude them
        indsSilValid = np.where(Mtot_kg <= Planet.Bulk.M_kg)[0]
        if(np.size(indsSilValid) == 0):
            raise RuntimeError('No silicate mantle size was less than the total body mass.\n' +
                               f'Min mass: {np.min(Mtot_kg/Planet.Bulk.M_kg):.3f} M_{Planet.name[0]}, ' +
                               f'max mass: {np.max(Mtot_kg/Planet.Bulk.M_kg):.3f} M_{Planet.name[0]}. ' +
                               'Try adjusting run settings that affect mantle density, like silicate composition ' +
                               'and heat flux settings.')

    return indsSilValid, nProfiles, Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, \
           MLayerSil_kg, MAboveSil_kg, gSil_ms2, phiSil_frac


def IronCoreLayers(Planet, Params,
                   indsSilValid, nSilTooBig, nProfiles, Psil_MPa, Tsil_K, rSil_m, MAboveSil_kg, gSil_ms2):
    """ Determines properties of core layers based on input Perple_X table
        and seafloor properties.

        Args:
            nSilTooBig (int): Number of silicate profiles to skip past due to masses exceeding body mass.
            Psil_MPa, Tsil_K, rSil_m, MLayerSil_kg, MAboveSil_kg, gSil_ms2 (float, shape NxM): Outputs from
                SilicateLayers with layer properties for each silicate region size possibility.
        Returns:
            nSilFinal (int): Index in silicate profiles of core with a total mass just under body mass.
            Pcore_MPa, Tcore_K, rCore_m, rhoCore_kgm3, MLayerCore_kg, gCore_ms2 (float, shape Planet.Steps.nCore):
                Core properties needed to determine MoI.
    """
    # Initialize output arrays and working arrays
    Pcore_MPa, Tcore_K, rhoCore_kgm3, MLayerCore_kg, gCore_ms2, CpCore_JkgK, alphaCore_pK = \
        (np.zeros((nProfiles-nSilTooBig, Planet.Steps.nCore)) for _ in range(7))
    rCore_m = np.zeros((nProfiles-nSilTooBig, Planet.Steps.nCore+1))
    # Initialize matching indices as -1 as a flag for unfilled values
    iCoreMatch, nSilFinal = (-1 * np.ones(nProfiles).astype(np.int_) for _ in range(2))

    if Planet.Do.CONSTANT_INNER_DENSITY:
        iCoreStart = [-1]
        silEnd = 0
        indsSilValid = [0]
    else:
        # Calculate maximum core size based on minimum plausible density setting
        MCore_kg = Planet.Bulk.M_kg - MAboveSil_kg[indsSilValid,:]
        rCoreMax_m = (MCore_kg/Planet.Core.rhoMin_kgm3 * 3/4/np.pi)**(1/3)
        # Find first silicate layer smaller than the max core radius
        iCoreStart = [next(j[0] for j,val in np.ndenumerate(rSil_m[i,:-1]) if val < rCoreMax_m[i-nSilTooBig,j])
                      for i in indsSilValid]
        silEnd = Planet.Steps.nSilMax

    log.debug('Evaluating core EOS for possible configurations...')
    for iValid in range(nProfiles - nSilTooBig):
        # Get profile index among silicate layers from number of profiles
        iProf = indsSilValid[iValid]
        # Get index for which silicate layer gets replaced by core layers there and below
        thisCoreStart = iCoreStart[iValid]
        # Get number of remaining silicate layers to iterate over for possible core configs
        nSilRemain = silEnd - thisCoreStart
        #(Re-)initialize placeholder arrays for each core possibility (they change length for each)
        thisPcore_MPa, thisTcore_K, thisrhoCore_kgm3, thisMLayerCore_kg, thisgCore_ms2, thisCpCore_JkgK, \
        thisalphaCore_pK = (np.zeros((nSilRemain, Planet.Steps.nCore)) for _ in range(7))

        # Set starting core values for all possibilities to be equal to silicates at this transition radius
        thisrCore_m = np.array([np.linspace(rSil_m[iProf,thisCoreStart+j], 0, Planet.Steps.nCore+1) for j in range(nSilRemain)])
        thisPcore_MPa[:,0] = [Psil_MPa[iProf,thisCoreStart+j] for j in range(nSilRemain)]
        thisTcore_K[:,0] = [Tsil_K[iProf,thisCoreStart+j] for j in range(nSilRemain)]
        thisrhoCore_kgm3[:,0] = Planet.Core.EOS.fn_rho_kgm3(thisPcore_MPa[:nSilRemain,0], thisTcore_K[:nSilRemain,0], grid=False)
        thisCpCore_JkgK[:,0] = Planet.Core.EOS.fn_Cp_JkgK(thisPcore_MPa[:nSilRemain,0], thisTcore_K[:nSilRemain,0], grid=False)
        thisalphaCore_pK[:,0] = Planet.Core.EOS.fn_alpha_pK(thisPcore_MPa[:nSilRemain,0], thisTcore_K[:nSilRemain,0], grid=False)
        thisMLayerCore_kg[:,0] = [thisrhoCore_kgm3[j,0] * 4/3*np.pi*(thisrCore_m[j,0]**3 - thisrCore_m[j,1]**3) for j in range(nSilRemain)]
        thisgCore_ms2[:,0] = [gSil_ms2[iProf,thisCoreStart+j] for j in range(nSilRemain)]
        MAbove_kg = np.array([MAboveSil_kg[iProf,thisCoreStart+j] for j in range(nSilRemain)])

        for k in range(1, Planet.Steps.nCore):
            MAbove_kg += thisMLayerCore_kg[:,k-1]
            thisDeltaP = 1e-6 * thisMLayerCore_kg[:,k-1] * thisgCore_ms2[:,k-1] / (4*np.pi*thisrCore_m[:,k]**2)
            thisPcore_MPa[:,k] = thisPcore_MPa[:,k-1] + thisDeltaP
            thisTcore_K[:,k] = thisTcore_K[:,k-1] + thisalphaCore_pK[:,k-1]*thisTcore_K[:,k] / \
                           thisCpCore_JkgK[:,k-1] / thisrhoCore_kgm3[:,k-1] * thisDeltaP*1e6
            thisrhoCore_kgm3[:,k] = Planet.Core.EOS.fn_rho_kgm3(thisPcore_MPa[:nSilRemain,k], thisTcore_K[:nSilRemain,k], grid=False)
            thisCpCore_JkgK[:,k] = Planet.Core.EOS.fn_Cp_JkgK(thisPcore_MPa[:nSilRemain,k], thisTcore_K[:nSilRemain,k], grid=False)
            thisalphaCore_pK[:,k] = Planet.Core.EOS.fn_alpha_pK(thisPcore_MPa[:nSilRemain,k], thisTcore_K[:nSilRemain,k], grid=False)
            thisMLayerCore_kg[:,k] = thisrhoCore_kgm3[:,k] * 4/3*np.pi*(thisrCore_m[:,k]**3 - thisrCore_m[:,k+1]**3)
            # Approximate gravity as linear to avoid blowing up for total mass less than body mass (accurate for constant density only)
            thisgCore_ms2[:,k] = thisgCore_ms2[:,0] * thisrCore_m[:,k] / thisrCore_m[:,0]

        if not Planet.Do.CONSTANT_INNER_DENSITY:
            # Find the first core profile that has a mass just below the body mass
            Mtot_kg = MAbove_kg + thisMLayerCore_kg[:,-1]
            iCoreMatch[iProf] = next(ii[0] for ii,val in np.ndenumerate(Mtot_kg) if val < Planet.Bulk.M_kg)
            nSilFinal[iProf] = iCoreStart[iValid] + iCoreMatch[iProf]
            log.debug(f'Core match for iProf = {iProf:d} with Steps.nSil = {nSilFinal[iProf]:d} ' +
                      f'and M = {Mtot_kg[iCoreMatch[iProf]]/Planet.Bulk.M_kg:.4f} M_{Planet.name[0]}.')

        # Assign the values for the core profile with matching total mass to output arrays
        Pcore_MPa[iValid,:] = thisPcore_MPa[iCoreMatch[iProf],:]
        Tcore_K[iValid,:] = thisTcore_K[iCoreMatch[iProf],:]
        rCore_m[iValid,:] = thisrCore_m[iCoreMatch[iProf],:]
        rhoCore_kgm3[iValid,:] = thisrhoCore_kgm3[iCoreMatch[iProf],:]
        MLayerCore_kg[iValid,:] = thisMLayerCore_kg[iCoreMatch[iProf],:]
        gCore_ms2[iValid,:] = thisgCore_ms2[iCoreMatch[iProf],:]
        CpCore_JkgK[iValid,:] = thisCpCore_JkgK[iCoreMatch[iProf],:]
        alphaCore_pK[iValid,:] = thisalphaCore_pK[iCoreMatch[iProf],:]

    return nSilFinal, Pcore_MPa, Tcore_K, rCore_m, rhoCore_kgm3, MLayerCore_kg, gCore_ms2, CpCore_JkgK, alphaCore_pK
