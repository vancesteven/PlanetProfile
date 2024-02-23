import numpy as np
import logging
import scipy.interpolate as spi
from PlanetProfile.Thermodynamics.HydroEOS import GetPhaseIndices, GetOceanEOS
from PlanetProfile.Thermodynamics.InnerEOS import GetInnerEOS
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist

# Assign logger
log = logging.getLogger('PlanetProfile')

def ViscosityCalcs(Planet, Params):
    # Initialize outputs as NaN so that we get errors if we missed any layers
    Planet.eta_Pas = np.zeros(Planet.Steps.nTotal) * np.nan

    # Only perform calculations if this is a valid profile
    if Planet.Do.VALID:
        # Identify which indices correspond to which phases
        indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund, \
            indsClath, indsClathWet, indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI, \
            indsFe = GetPhaseIndices(Planet.phase)

        if Params.CALC_VISCOSITY:
            # Make sure the necessary EOSs have been loaded (mainly only important in parallel ExploreOgram runs)
            if not Planet.Do.NO_H2O and Planet.Ocean.EOS.key not in EOSlist.loaded.keys():
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
                                               etaFixed_Pas=None)  # Causes ocean EOS to use default behavior for this comp

            if Planet.Do.POROUS_ICE:
                Planet = CalcViscPorIce(Planet, Params, indsLiq, indsI, indsIwet, indsII, indsIIund,
                                        indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund,
                                        indsClath, indsClathWet)
            else:
                Planet = CalcViscSolidIce(Planet, Params, indsLiq, indsI, indsII, indsIIund,
                                          indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund,
                                          indsClath)

            if not Params.SKIP_INNER:
                if Planet.Do.POROUS_ROCK:
                    Planet = CalcViscPorRock(Planet, Params, indsSil, indsSilLiq, indsSilI,
                                             indsSilII, indsSilIII, indsSilV, indsSilVI)
                else:
                    Planet = CalcViscSolidRock(Planet, Params, indsSil)

                Planet = CalcViscCore(Planet, Params, indsFe)

    return Planet


def CalcViscPorIce(Planet, Params, indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund,
                                   indsV, indsVund, indsVI, indsVIund, indsClath, indsClathWet):

    # First do ocean (if present) and dry surface ice and clathrates
    if np.size(indsLiq) != 0:
        Planet.eta_Pas[indsLiq] = Planet.Ocean.EOS.fn_eta_Pas(Planet.P_MPa[indsLiq], Planet.T_K[indsLiq])
    if np.size(indsI) != 0:
        Planet.eta_Pas[indsI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Ocean.surfIceEOS['Ih'].fn_eta_Pas(Planet.P_MPa[indsI], Planet.T_K[indsI]), 0,
            Planet.phi_frac[indsI], Planet.Ocean.Jvisc)
    if np.size(indsClath) != 0:
        Planet.eta_Pas[indsClath] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Ocean.surfIceEOS['Clath'].fn_eta_Pas(Planet.P_MPa[indsClath], Planet.T_K[indsClath]), 0,
            Planet.phi_frac[indsClath], Planet.Ocean.Jvisc)
    # We use the negative underplate phase IDs for dry HP ices
    if np.size(indsIIund) != 0:
        Planet.eta_Pas[indsIIund] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Ocean.surfIceEOS['II'].fn_eta_Pas(Planet.P_MPa[indsIIund], Planet.T_K[indsIIund]), 0,
            Planet.phi_frac[indsIIund], Planet.Ocean.Jvisc)
    if np.size(indsIIIund) != 0:
        Planet.eta_Pas[indsIIIund] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Ocean.surfIceEOS['III'].fn_eta_Pas(Planet.P_MPa[indsIIIund], Planet.T_K[indsIIIund]), 0,
            Planet.phi_frac[indsIIIund], Planet.Ocean.Jvisc)
    if np.size(indsVund) != 0:
        Planet.eta_Pas[indsVund] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Ocean.surfIceEOS['V'].fn_eta_Pas(Planet.P_MPa[indsVund], Planet.T_K[indsVund]), 0,
            Planet.phi_frac[indsVund], Planet.Ocean.Jvisc)
    if np.size(indsVIund) != 0:
        Planet.eta_Pas[indsVIund] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Ocean.surfIceEOS['VI'].fn_eta_Pas(Planet.P_MPa[indsVIund], Planet.T_K[indsVIund]), 0,
            Planet.phi_frac[indsVIund], Planet.Ocean.Jvisc)
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

    return Planet


def CalcElecSolidIce(Planet, Params, indsLiq, indsI, indsII, indsIIund, indsIII, indsIIIund,
                                     indsV, indsVund, indsVI, indsVIund, indsClath):

    # Calculate and/or assign conductivities for each phase type
    if np.size(indsLiq) != 0:
        Planet.eta_Pas[indsLiq] = Planet.Ocean.EOS.fn_eta_Pas(Planet.P_MPa[indsLiq], Planet.T_K[indsLiq])

    Planet.eta_Pas[indsI] = Planet.Ocean.surfIceEOS['Ih'].fn_eta_Pas(Planet.P_MPa[indsI], Planet.T_K[indsI])

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    # Below here has not been finished yet.
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    # Need to add catch for is inds are size 0, and if not, apply the correct surfIceEOS or iceEOS
    # evaluation of the EOS viscosity function.

    indsIIall = np.concatenate((indsII, indsIIund))
    indsIIIall = np.concatenate((indsIII, indsIIIund))
    indsVall = np.concatenate((indsV, indsVund))
    indsVIall = np.concatenate((indsVI, indsVIund))
    Planet.eta_Pas[indsIIall] =  Planet.Ocean.sigmaIce_Sm['II']
    Planet.eta_Pas[indsIIIall] = Planet.Ocean.sigmaIce_Sm['III']
    Planet.eta_Pas[indsVall] =   Planet.Ocean.sigmaIce_Sm['V']
    Planet.eta_Pas[indsVIall] =  Planet.Ocean.sigmaIce_Sm['VI']
    Planet.eta_Pas[indsClath] =  Planet.Ocean.sigmaIce_Sm['Clath']

    return Planet


def CalcViscPorRock(Planet, Params, indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI):

    # Need to convert this functionality to viscosity calcs and create a parallel one for solid
    # rocks and another for the core.

    # First, get pore fluid conductivity
    sigmaFluid_Sm = Planet.Sil.poreEOS.fn_sigma_Sm(Planet.Ppore_MPa[indsSilLiq], Planet.T_K[indsSilLiq])
    # Account for possible errors in pore fluid conductivity calcs
    validSigs = np.logical_not(np.isnan(sigmaFluid_Sm))
    if np.size(sigmaFluid_Sm) > 0:
        # Interpolate over NaNs to remove them
        sigmaFluid_Sm = spi.griddata(Planet.r_m[indsSilLiq][validSigs], sigmaFluid_Sm[validSigs], Planet.r_m[indsSilLiq])

    # Initialize conductivity array for all pore materials
    sigmaPore_Sm = np.zeros_like(Planet.sigma_Sm)
    # Fill with appropriate values based on ice phase within pore material
    sigmaPore_Sm[indsSilLiq] = sigmaFluid_Sm
    sigmaPore_Sm[indsSilI] =   Planet.Ocean.sigmaIce_Sm['Ih']
    sigmaPore_Sm[indsSilII] =  Planet.Ocean.sigmaIce_Sm['II']
    sigmaPore_Sm[indsSilIII] = Planet.Ocean.sigmaIce_Sm['III']
    sigmaPore_Sm[indsSilV] =   Planet.Ocean.sigmaIce_Sm['V']
    sigmaPore_Sm[indsSilVI] =  Planet.Ocean.sigmaIce_Sm['VI']
    # Chop sigmaPore down to just the silicate indices
    sigmaPore_Sm = sigmaPore_Sm[indsSil]

    # Next, assign varying dependence based on porosity threshold
    belowThresh = Planet.phi_frac[indsSil] <  Planet.Sil.poreConductThresh_frac
    aboveThresh = Planet.phi_frac[indsSil] >= Planet.Sil.poreConductThresh_frac
    sigmaPore_Sm[belowThresh] = Planet.Sil.poreConductPrefac * \
                                sigmaPore_Sm[belowThresh] * Planet.phi_frac[indsSil][belowThresh]**Planet.Sil.poreConductBelowExp
    sigmaPore_Sm[aboveThresh] = sigmaPore_Sm[aboveThresh] * Planet.phi_frac[indsSil][aboveThresh]**Planet.Sil.poreConductAboveExp
    # Now combine with the conductivity of the rock matrix
    Planet.sigma_Sm[indsSil] = Planet.Sil.EOS.fn_porosCorrect(Planet.Sil.sigmaSil_Sm,
        sigmaPore_Sm, Planet.phi_frac[indsSil], Planet.Sil.Jvisc)

    # Record the mean conductivities of the pore fluids and porous layers,
    # but include only the contributing layers in the calculation
    if np.size(indsSilLiq) != 0:
        aboveThreshLiq = Planet.phi_frac[indsSilLiq] >= Planet.Sil.poreConductThresh_frac
        if np.any(aboveThreshLiq):
            Planet.Sil.sigmaPoreMean_Sm = np.mean(sigmaFluid_Sm[aboveThreshLiq])
            Planet.Sil.sigmaPorousLayerMean_Sm = np.mean(Planet.sigma_Sm[indsSilLiq][aboveThreshLiq])
        else:
            Planet.Sil.sigmaPoreMean_Sm = Constants.sigmaDef_Sm
            Planet.Sil.sigmaPorousLayerMean_Sm = Constants.sigmaDef_Sm
    else:
        # No pores contain liquid (only HP ices), so just record the mean combined conductivity
        Planet.Sil.sigmaPoreMean_Sm = np.nan
        Planet.Sil.sigmaPorousLayerMean_Sm = np.mean(Planet.sigma_Sm[indsSil][aboveThresh])

    return Planet


class ViscOceanUniform_Pas:
    def __init__(self, etaSet_Pas=None, comp=None):
        if etaSet_Pas is None:
            if comp == 'Seawater':
                self.eta_Pas = Constants.etaSeawater_Pas
            else:
                self.eta_Pas = Constants.etaH2O_Pas
        else:
            self.eta_Pas = etaSet_Pas

    def __call__(self, P_MPa, T_K, grid=False):
        if grid:
            return np.zeros((np.size(P_MPa), np.size(T_K))) + self.eta_Pas
        else:
            return np.zeros_like(P_MPa) + self.eta_Pas


class ViscIceUniform_Pas:
    def __init__(self, etaSet_Pas=None, TviscTrans_K=None):
        if etaSet_Pas is None:
            self.eta_Pas = Constants.etaIce_Pas
        else:
            self.eta_Pas = etaSet_Pas

        if TviscTrans_K is None:
            self.TviscTrans_K = Constants.TviscIce_K
        else:
            self.TviscTrans_K = TviscTrans_K

    def __call__(self, P_MPa, T_K, grid=False):
        Ttrans_K = np.insert([0.0, np.inf], 1, self.TviscTrans_K)
        if grid:
            eta_Pas = np.zeros((np.size(P_MPa), np.size(T_K)))
            for Tlow_K, Tupp_K, etaConst_Pas in zip(Ttrans_K[:-1], Ttrans_K[1:], self.eta_Pas):
                eta_Pas[:, np.logical_and(T_K >= Tlow_K, T_K < Tupp_K)] = etaConst_Pas
        else:
            eta_Pas = np.zeros_like(P_MPa)
            for Tlow_K, Tupp_K, etaConst_Pas in zip(Ttrans_K[:-1], Ttrans_K[1:], self.eta_Pas):
                eta_Pas[:, np.logical_and(T_K >= Tlow_K, T_K < Tupp_K)] = etaConst_Pas

        return eta_Pas


class ViscRockUniform_Pas:
    def __init__(self, etaSet_Pas=None, TviscTrans_K=None):
        if etaSet_Pas is None:
            self.eta_Pas = Constants.etaRock_Pas
        else:
            self.eta_Pas = etaSet_Pas

        if TviscTrans_K is None:
            self.TviscTrans_K = Constants.TviscRock_K
        else:
            self.TviscTrans_K = TviscTrans_K

    def __call__(self, P_MPa, T_K, grid=False):
        Ttrans_K = np.insert([0.0, np.inf], 1, self.TviscTrans_K)
        if grid:
            eta_Pas = np.zeros((np.size(P_MPa), np.size(T_K)))
            for Tlow_K, Tupp_K, etaConst_Pas in zip(Ttrans_K[:-1], Ttrans_K[1:], self.eta_Pas):
                eta_Pas[:, np.logical_and(T_K >= Tlow_K, T_K < Tupp_K)] = etaConst_Pas
        else:
            eta_Pas = np.zeros_like(P_MPa)
            for Tlow_K, Tupp_K, etaConst_Pas in zip(Ttrans_K[:-1], Ttrans_K[1:], self.eta_Pas):
                eta_Pas[:, np.logical_and(T_K >= Tlow_K, T_K < Tupp_K)] = etaConst_Pas

        return eta_Pas


class ViscCoreUniform_Pas:
    def __init__(self, etaSet_Pas=None, TviscTrans_K=None):
        if etaSet_Pas is None:
            self.eta_Pas = [Constants.etaFeSolid_Pas, Constants.etaFeLiquid_Pas]
        else:
            self.eta_Pas = etaSet_Pas

        if TviscTrans_K is None:
            self.TviscTrans_K = Constants.TviscFe_K
        else:
            self.TviscTrans_K = TviscTrans_K

    def __call__(self, P_MPa, T_K, grid=False):
        Ttrans_K = np.insert([0.0, np.inf], 1, self.TviscTrans_K)
        if grid:
            eta_Pas = np.zeros((np.size(P_MPa), np.size(T_K)))
            for Tlow_K, Tupp_K, etaConst_Pas in zip(Ttrans_K[:-1], Ttrans_K[1:], self.eta_Pas):
                eta_Pas[:, np.logical_and(T_K >= Tlow_K, T_K < Tupp_K)] = etaConst_Pas
        else:
            eta_Pas = np.zeros_like(P_MPa)
            for Tlow_K, Tupp_K, etaConst_Pas in zip(Ttrans_K[:-1], Ttrans_K[1:], self.eta_Pas):
                eta_Pas[:, np.logical_and(T_K >= Tlow_K, T_K < Tupp_K)] = etaConst_Pas

        return eta_Pas
