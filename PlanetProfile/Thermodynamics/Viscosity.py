import numpy as np
import logging
import scipy.interpolate as spi
from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS
from PlanetProfile.Utilities.defineStructs import EOSlist, Timing
import time
from PlanetProfile.Utilities.Indexing import GetPhaseIndices

# Assign logger
log = logging.getLogger('PlanetProfile')

def ViscosityCalcs(Planet, Params):
    Timing.setFunctionTime(time.time())
    # Initialize outputs as NaN so that we get errors if we missed any layers
    Planet.eta_Pas = np.zeros(Planet.Steps.nTotal) * np.nan

    # Only perform calculations if this is a valid profile
    if Planet.Do.VALID or (Params.ALLOW_BROKEN_MODELS and Planet.Do.STILL_CALCULATE_BROKEN_PROPERTIES):
        # Identify which indices correspond to which phases
        indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund, \
            indsClath, indsClathWet, indsMixedClathrateIh, indsMixedClathrateII, indsMixedClathrateIII, indsMixedClathrateV, indsMixedClathrateVI, \
            indsMixedClathrateIhwet, indsMixedClathrateIIund, indsMixedClathrateIIIund, indsMixedClathrateVund, indsMixedClathrateVIund, \
            indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI, \
            indsFe = GetPhaseIndices(Planet.phase)

        if Params.CALC_VISCOSITY:
            # Make sure the necessary EOSs have been loaded (mainly only important in parallel ExploreOgram runs)
            if not (Planet.Do.NO_H2O or Planet.Do.NO_OCEAN) and Planet.Ocean.EOS.key not in EOSlist.loaded.keys():
                POcean_MPa = np.arange(Planet.PfreezeLower_MPa, Planet.Ocean.PHydroMax_MPa,
                                       Planet.Ocean.deltaP)
                TOcean_K = np.arange(Planet.Bulk.Tb_K, Planet.Ocean.THydroMax_K,
                                     Planet.Ocean.deltaT)
                Planet.Ocean.EOS = GetOceanEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt,
                                               POcean_MPa, TOcean_K,
                                               Planet.Ocean.MgSO4elecType,
                                               rhoType=Planet.Ocean.MgSO4rhoType,
                                               scalingType=Planet.Ocean.MgSO4scalingType,
                                               FORCE_NEW=Params.FORCE_EOS_RECALC,
                                               phaseType=Planet.Ocean.phaseType,
                                               EXTRAP=Params.EXTRAP_OCEAN,
                                               sigmaFixed_Sm=Planet.Ocean.sigmaFixed_Sm,
                                               etaFixed_Pas=None, kThermConst_WmK=Planet.Ocean.kThermWater_WmK,
                                               propsStepReductionFactor=Planet.Ocean.propsStepReductionFactor)

            if Planet.Do.POROUS_ICE:
                Planet = CalcViscPorIce(Planet, Params, indsLiq, indsI, indsIwet, indsII, indsIIund,
                                        indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund,
                                        indsClath, indsClathWet, indsMixedClathrateIh, indsMixedClathrateII, indsMixedClathrateIII, indsMixedClathrateV, indsMixedClathrateVI, \
                                        indsMixedClathrateIhwet, indsMixedClathrateIIund, indsMixedClathrateIIIund, indsMixedClathrateVund, indsMixedClathrateVIund)
            else:
                Planet = CalcViscSolidIce(Planet, Params, indsLiq, indsI, indsII, indsIIund,
                                          indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund,
                                          indsClath, indsMixedClathrateIh, indsMixedClathrateII, indsMixedClathrateIII, indsMixedClathrateV, indsMixedClathrateVI, \
                                          indsMixedClathrateIIund, indsMixedClathrateIIIund, indsMixedClathrateVund, indsMixedClathrateVIund)

            if not Params.SKIP_INNER:
                if Planet.Do.POROUS_ROCK:
                    Planet = CalcViscPorRock(Planet, Params, indsSil, indsSilLiq, indsSilI,
                                             indsSilII, indsSilIII, indsSilV, indsSilVI)
                else:
                    Planet = CalcViscSolidRock(Planet, Params, indsSil)

                if Planet.Do.Fe_CORE:
                    Planet = CalcViscCore(Planet, Params, indsFe)

    Timing.printFunctionTimeDifference('ViscosityCalcs()', time.time())
    return Planet


def CalcViscPorIce(Planet, Params, indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund,
                                   indsV, indsVund, indsVI, indsVIund, indsClath, indsClathWet,
                                   indsMixedClathrateIh, indsMixedClathrateII, indsMixedClathrateIII, indsMixedClathrateV, indsMixedClathrateVI, \
                                        indsMixedClathrateIhwet, indsMixedClathrateIIund, indsMixedClathrateIIIund, indsMixedClathrateVund, indsMixedClathrateVIund):

    # First do ocean (if present) and dry surface ice and clathrates
    if np.size(indsLiq) != 0:
        Planet.eta_Pas[indsLiq] = Planet.Ocean.EOS.fn_eta_Pas(Planet.P_MPa[indsLiq], Planet.T_K[indsLiq])
    if np.size(indsI) != 0:
        Planet.eta_Pas[indsI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Ocean.surfIceEOS['Ih'].fn_eta_Pas(Planet.P_MPa[indsI], Planet.T_K[indsI]), 0,
            Planet.phi_frac[indsI], Planet.Ocean.Jvisc)
    if np.size(indsClath) != 0:
        Planet.eta_Pas[indsClath] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
            Planet.Ocean.surfIceEOS['Clath'].fn_eta_Pas(Planet.P_MPa[indsClath], Planet.T_K[indsClath]), 0,
            Planet.phi_frac[indsClath], Planet.Ocean.Jvisc)
    # We use the negative underplate phase IDs for dry HP ices
    if np.size(indsIIund) != 0:
        Planet.eta_Pas[indsIIund] = Planet.Ocean.surfIceEOS['II'].fn_porosCorrect(
            Planet.Ocean.surfIceEOS['II'].fn_eta_Pas(Planet.P_MPa[indsIIund], Planet.T_K[indsIIund]), 0,
            Planet.phi_frac[indsIIund], Planet.Ocean.Jvisc)
    if np.size(indsIIIund) != 0:
        Planet.eta_Pas[indsIIIund] = Planet.Ocean.surfIceEOS['III'].fn_porosCorrect(
            Planet.Ocean.surfIceEOS['III'].fn_eta_Pas(Planet.P_MPa[indsIIIund], Planet.T_K[indsIIIund]), 0,
            Planet.phi_frac[indsIIIund], Planet.Ocean.Jvisc)
    if np.size(indsVund) != 0:
        Planet.eta_Pas[indsVund] = Planet.Ocean.surfIceEOS['V'].fn_porosCorrect(
            Planet.Ocean.surfIceEOS['V'].fn_eta_Pas(Planet.P_MPa[indsVund], Planet.T_K[indsVund]), 0,
            Planet.phi_frac[indsVund], Planet.Ocean.Jvisc)
    if np.size(indsVIund) != 0:
        Planet.eta_Pas[indsVIund] = Planet.Ocean.surfIceEOS['VI'].fn_porosCorrect(
            Planet.Ocean.surfIceEOS['VI'].fn_eta_Pas(Planet.P_MPa[indsVIund], Planet.T_K[indsVIund]), 0,
            Planet.phi_frac[indsVIund], Planet.Ocean.Jvisc)
    # Get all mixed clathrate phases that are not of size zero
    indsMixedClathrateAll = [indsMixedClathrateIh, indsMixedClathrateIIund, indsMixedClathrateIIIund, indsMixedClathrateVund, indsMixedClathrateVIund]
    mixedPhases = ['MixedClathrateIh', 'MixedClathrateII', 'MixedClathrateIII', 'MixedClathrateV', 'MixedClathrateVI']
    for indsMixedClathrate, mixedPhase in zip(indsMixedClathrateAll, mixedPhases):
        if np.size(indsMixedClathrate) != 0:
            Planet.eta_Pas[indsMixedClathrate] = Planet.Ocean.surfIceEOS[mixedPhase].fn_porosCorrect(
                Planet.Ocean.surfIceEOS[mixedPhase].fn_eta_Pas(Planet.P_MPa[indsMixedClathrate], Planet.T_K[indsMixedClathrate]), 0,
                Planet.phi_frac[indsMixedClathrate], Planet.Ocean.Jvisc)
            
    # Next, do liquid-filled ice and clathrate pores
    if np.size(indsIwet) != 0:
        etaFluid_Pas = Planet.Ocean.EOS.fn_eta_Pas(Planet.Ppore_MPa[indsIwet], Planet.T_K[indsIwet])
        Planet.eta_Pas[indsIwet] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Ocean.surfIceEOS['Ih'].fn_eta_Pas(Planet.P_MPa[indsIwet], Planet.T_K[indsIwet]),
            etaFluid_Pas, Planet.phi_frac[indsIwet], Planet.Ocean.Jvisc)
    if np.size(indsClathWet) != 0:
        etaFluid_Pas = Planet.Ocean.EOS.fn_eta_Pas(Planet.Ppore_MPa[indsClathWet], Planet.T_K[indsClathWet])
        Planet.eta_Pas[indsClathWet] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
            Planet.Ocean.surfIceEOS['Clath'].fn_eta_Pas(Planet.P_MPa[indsClathWet], Planet.T_K[indsClathWet]),
            etaFluid_Pas, Planet.phi_frac[indsClathWet], Planet.Ocean.Jvisc)
    if np.size(indsII) != 0:
        etaFluid_Pas = Planet.Ocean.EOS.fn_eta_Pas(Planet.Ppore_MPa[indsII], Planet.T_K[indsII])
        Planet.eta_Pas[indsII] = Planet.Ocean.iceEOS['II'].fn_porosCorrect(
            Planet.Ocean.iceEOS['II'].fn_eta_Pas(Planet.P_MPa[indsII], Planet.T_K[indsII]),
            etaFluid_Pas, Planet.phi_frac[indsII], Planet.Ocean.Jvisc)
    if np.size(indsIII) != 0:
        etaFluid_Pas = Planet.Ocean.EOS.fn_eta_Pas(Planet.Ppore_MPa[indsIII], Planet.T_K[indsIII])
        Planet.eta_Pas[indsIII] = Planet.Ocean.iceEOS['III'].fn_porosCorrect(
            Planet.Ocean.iceEOS['III'].fn_eta_Pas(Planet.P_MPa[indsIII], Planet.T_K[indsIII]),
            etaFluid_Pas, Planet.phi_frac[indsIII], Planet.Ocean.Jvisc)
    if np.size(indsV) != 0:
        etaFluid_Pas = Planet.Ocean.EOS.fn_eta_Pas(Planet.Ppore_MPa[indsV], Planet.T_K[indsV])
        Planet.eta_Pas[indsV] = Planet.Ocean.iceEOS['V'].fn_porosCorrect(
            Planet.Ocean.iceEOS['V'].fn_eta_Pas(Planet.P_MPa[indsV], Planet.T_K[indsV]),
            etaFluid_Pas, Planet.phi_frac[indsV], Planet.Ocean.Jvisc)
    if np.size(indsVI) != 0:
        etaFluid_Pas = Planet.Ocean.EOS.fn_eta_Pas(Planet.Ppore_MPa[indsVI], Planet.T_K[indsVI])
        Planet.eta_Pas[indsVI] = Planet.Ocean.iceEOS['VI'].fn_porosCorrect(
            Planet.Ocean.iceEOS['VI'].fn_eta_Pas(Planet.P_MPa[indsVI], Planet.T_K[indsVI]),
            etaFluid_Pas, Planet.phi_frac[indsVI], Planet.Ocean.Jvisc)
    # Get all mixed clathrate phases that are not of size zero
    indsMixedWetClathrateAll = [indsMixedClathrateIhwet, indsMixedClathrateII, indsMixedClathrateIII, indsMixedClathrateV, indsMixedClathrateVI]
    mixedPhases = ['MixedClathrateIh', 'MixedClathrateII', 'MixedClathrateIII', 'MixedClathrateV', 'MixedClathrateVI']
    for indsMixedClathrate, mixedPhase in zip(indsMixedWetClathrateAll, mixedPhases):
        if np.size(indsMixedClathrate) != 0:
            Planet.eta_Pas[indsMixedClathrate] = Planet.Ocean.iceEOS[mixedPhase].fn_porosCorrect(
                Planet.Ocean.iceEOS[mixedPhase].fn_eta_Pas(Planet.P_MPa[indsMixedClathrate], Planet.T_K[indsMixedClathrate]), 0,
                Planet.phi_frac[indsMixedClathrate], Planet.Ocean.Jvisc)

    return Planet


def CalcViscSolidIce(Planet, Params, indsLiq, indsI, indsII, indsIIund, indsIII, indsIIIund,
                                     indsV, indsVund, indsVI, indsVIund, indsClath, indsMixedClathrateIh, indsMixedClathrateII, indsMixedClathrateIII, indsMixedClathrateV, indsMixedClathrateVI, \
                                          indsMixedClathrateIIund, indsMixedClathrateIIIund, indsMixedClathrateVund, indsMixedClathrateVIund):

    # Calculate and/or assign conductivities for each phase type
    if np.size(indsLiq) != 0:
        Planet.eta_Pas[indsLiq] = Planet.Ocean.EOS.fn_eta_Pas(Planet.P_MPa[indsLiq], Planet.T_K[indsLiq])

    if np.size(indsI) != 0:
        Planet.eta_Pas[indsI] = Planet.Ocean.surfIceEOS['Ih'].fn_eta_Pas(Planet.P_MPa[indsI], Planet.T_K[indsI])

    if np.size(indsIIund) != 0:
        Planet.eta_Pas[indsIIund] = Planet.Ocean.surfIceEOS['II'].fn_eta_Pas(Planet.P_MPa[indsIIund], Planet.T_K[indsIIund])
    if np.size(indsII) != 0:
        Planet.eta_Pas[indsII] = Planet.Ocean.iceEOS['II'].fn_eta_Pas(Planet.P_MPa[indsII], Planet.T_K[indsII])

    if np.size(indsIIIund) != 0:
        Planet.eta_Pas[indsIIIund] = Planet.Ocean.surfIceEOS['III'].fn_eta_Pas(Planet.P_MPa[indsIIIund], Planet.T_K[indsIIIund])
    if np.size(indsIII) != 0:
        Planet.eta_Pas[indsIII] = Planet.Ocean.iceEOS['III'].fn_eta_Pas(Planet.P_MPa[indsIII], Planet.T_K[indsIII])

    if np.size(indsVund) != 0:
        Planet.eta_Pas[indsVund] = Planet.Ocean.surfIceEOS['V'].fn_eta_Pas(Planet.P_MPa[indsVund], Planet.T_K[indsVund])
    if np.size(indsV) != 0:
        Planet.eta_Pas[indsV] = Planet.Ocean.iceEOS['V'].fn_eta_Pas(Planet.P_MPa[indsV], Planet.T_K[indsV])

    if np.size(indsVIund) != 0:
        Planet.eta_Pas[indsVIund] = Planet.Ocean.surfIceEOS['VI'].fn_eta_Pas(Planet.P_MPa[indsVIund], Planet.T_K[indsVIund])
    if np.size(indsVI) != 0:
        Planet.eta_Pas[indsVI] = Planet.Ocean.iceEOS['VI'].fn_eta_Pas(Planet.P_MPa[indsVI], Planet.T_K[indsVI])

    if np.size(indsClath) != 0:
        Planet.eta_Pas[indsClath] = Planet.Ocean.surfIceEOS['Clath'].fn_eta_Pas(Planet.P_MPa[indsClath], Planet.T_K[indsClath])

    # Get all mixed clathrate phases that are not of size zero
    indsAllMixedClathDry = (indsMixedClathrateIh, indsMixedClathrateIIund, indsMixedClathrateIIIund, indsMixedClathrateVund, indsMixedClathrateVIund)
    clathPhase = ('MixedClathrateIh',  'MixedClathrateII',  'MixedClathrateIII', 'MixedClathrateV', 'MixedClathrateVI')
    for indsMixedClath, clathPhase in zip(indsAllMixedClathDry, clathPhase):
        if np.size(indsMixedClath) != 0:
            Planet.eta_Pas[indsMixedClath] = Planet.Ocean.surfIceEOS[clathPhase].fn_eta_Pas(Planet.P_MPa[indsMixedClath], Planet.T_K[indsMixedClath])
    
    indsAllMixedClathWet = (indsMixedClathrateII, indsMixedClathrateIII, indsMixedClathrateV, indsMixedClathrateVI)
    clathPhase = ('MixedClathrateII',  'MixedClathrateIII', 'MixedClathrateV', 'MixedClathrateVI')
    for indsMixedClath, clathPhase in zip(indsAllMixedClathWet, clathPhase):
        if np.size(indsMixedClath) != 0:
            Planet.eta_Pas[indsMixedClath] = Planet.Ocean.iceEOS[clathPhase].fn_eta_Pas(Planet.P_MPa[indsMixedClath], Planet.T_K[indsMixedClath])

    return Planet


def CalcViscPorRock(Planet, Params, indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI):

    # Initialize viscosity array for all pore materials
    etaPore_Pas = np.zeros_like(Planet.eta_Pas)

    # Fill with appropriate values based on ice phase within pore material
    if np.size(indsSilLiq) > 0:
        # Get pore fluid viscosity
        etaFluid_Pas = Planet.Sil.poreEOS.fn_eta_Pas(Planet.Ppore_MPa[indsSilLiq], Planet.T_K[indsSilLiq])
        # Account for possible errors in pore fluid conductivity calcs
        validEtas = np.logical_not(np.isnan(etaFluid_Pas))
        # Interpolate over NaNs to remove them
        etaPore_Pas[indsSilLiq] = spi.griddata(Planet.r_m[indsSilLiq][validEtas], etaFluid_Pas[validEtas], Planet.r_m[indsSilLiq])

    if np.size(indsSilI) > 0:
        etaPore_Pas[indsSilI] =   Planet.Ocean.surfIceEOS['Ih'].fn_eta_Pas(Planet.Ppore_MPa[indsSilI], Planet.T_K[indsSilI])
    if np.size(indsSilII) > 0:
        etaPore_Pas[indsSilII] =  Planet.Ocean.iceEOS['II'].fn_eta_Pas(Planet.Ppore_MPa[indsSilII], Planet.T_K[indsSilII])
    if np.size(indsSilIII) > 0:
        etaPore_Pas[indsSilIII] =  Planet.Ocean.iceEOS['III'].fn_eta_Pas(Planet.Ppore_MPa[indsSilIII], Planet.T_K[indsSilIII])
    if np.size(indsSilV) > 0:
        etaPore_Pas[indsSilV] =  Planet.Ocean.iceEOS['V'].fn_eta_Pas(Planet.Ppore_MPa[indsSilV], Planet.T_K[indsSilV])
    if np.size(indsSilVI) > 0:
        etaPore_Pas[indsSilVI] =  Planet.Ocean.iceEOS['VI'].fn_eta_Pas(Planet.Ppore_MPa[indsSilVI], Planet.T_K[indsSilVI])
    # Chop sigmaPore down to just the silicate indices
    etaPore_Pas = etaPore_Pas[indsSil]

    # Now combine with the viscosity of the rock matrix
    Planet.eta_Pas[indsSil] = Planet.Sil.EOS.fn_porosCorrect(
        Planet.Sil.EOS.fn_eta_Pas(Planet.P_MPa[indsSil], Planet.T_K[indsSil]),
        etaPore_Pas, Planet.phi_frac[indsSil], Planet.Sil.Jvisc)

    return Planet


def CalcViscSolidRock(Planet, Params, indsSil):

    Planet.eta_Pas[indsSil] = Planet.Sil.EOS.fn_eta_Pas(Planet.P_MPa[indsSil], Planet.T_K[indsSil])

    return Planet


def CalcViscCore(Planet, Params, indsFe):

    Planet.eta_Pas[indsFe] = Planet.Core.EOS.fn_eta_Pas(Planet.P_MPa[indsFe], Planet.T_K[indsFe])

    return Planet
