import numpy as np
import logging
import scipy.interpolate as spi
from PlanetProfile.Thermodynamics.HydroEOS import GetPhaseIndices
from PlanetProfile.Utilities.defineStructs import Constants

# Assign logger
log = logging.getLogger('PlanetProfile')

def ElecConduct(Planet, Params):
    """ Calculate/assign electrical conductivities for each layer

        Assigns Planet attributes:
            sigma_Sm
    """
    # Initialize outputs as NaN so that we get errors if we missed any layers
    Planet.sigma_Sm = np.zeros(Planet.Steps.nTotal) * np.nan

    # Only perform calculations if this is a valid profile
    if Planet.Do.VALID:
        # Identify which indices correspond to which phases
        indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund, \
            indsClath, indsClathWet, indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI, \
            indsFe = GetPhaseIndices(Planet.phase)

        if Params.CALC_CONDUCT:
            if Planet.Do.POROUS_ICE:
                Planet = CalcElecPorIce(Planet, Params, indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund,
                                                        indsV, indsVund, indsVI, indsVIund, indsClath, indsClathWet)
            else:
                Planet = CalcElecSolidIce(Planet, Params, indsLiq, indsI, indsII, indsIIund, indsIII, indsIIIund,
                                                          indsV, indsVund, indsVI, indsVIund, indsClath)

            if not Params.SKIP_INNER:
                if Planet.Do.POROUS_ROCK:
                    Planet = CalcElecPorRock(Planet, Params, indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI)
                else:
                    Planet.sigma_Sm[indsSil] = Planet.Sil.sigmaSil_Sm
                    Planet.Sil.sigmaPoreMean_Sm = np.nan
                    Planet.Sil.sigmaPorousLayerMean_Sm = np.nan

                Planet.sigma_Sm[indsFe] = Planet.Core.sigmaCore_Sm
            else:
                Planet.Sil.sigmaPoreMean_Sm = np.nan
                Planet.Sil.sigmaPorousLayerMean_Sm = np.nan
        else:
            Planet.Sil.sigmaPoreMean_Sm = np.nan
            Planet.Sil.sigmaPorousLayerMean_Sm = np.nan

        if np.size(indsLiq) != 0:
            Planet.Ocean.sigmaMean_Sm = np.mean(Planet.sigma_Sm[indsLiq])
            Planet.Ocean.sigmaTop_Sm = Planet.sigma_Sm[indsLiq[0]]
        else:
            Planet.Ocean.sigmaMean_Sm = np.nan
            Planet.Ocean.sigmaTop_Sm = np.nan

    else:
        Planet.Sil.sigmaPoreMean_Sm = np.nan
        Planet.Sil.sigmaPorousLayerMean_Sm = np.nan
        Planet.Ocean.sigmaMean_Sm = np.nan
        Planet.Ocean.sigmaTop_Sm = np.nan

    return Planet


def CalcElecPorIce(Planet, Params, indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund,
                                   indsV, indsVund, indsVI, indsVIund, indsClath, indsClathWet):
    """ Calculation of electrical conductivities for porous ice, which is assumed to contain
        either evacuated pores or ocean fluids in pores.

        Args:
            indsLiq (int, shape Ni): Indices of Planet layer arrays (e.g. Planet.sigma_Sm)
                where the layers are ocean fluid.
            indsI, indsIIund, indsIIIund, indsVund, indsVIund, indsClath (int, shape Nj):
                Indices relating to (dry) surface ices, where pores contain vacuum.
            indsIwet, indsII, indsIII, indsV, indsVI, indsClathWet (int, shape Nk):
                Indices relating to ices found within the ocean or containing ocean fluids
                (near the bottom of the ice shell).
        Assigns Planet attributes:
            sigma_Sm
    """

    # First do ocean (if present) and dry surface ice and clathrates
    if np.size(indsLiq) != 0:
        Planet.sigma_Sm[indsLiq] = Planet.Ocean.EOS.fn_sigma_Sm(Planet.P_MPa[indsLiq], Planet.T_K[indsLiq])
    if np.size(indsI) != 0:
        Planet.sigma_Sm[indsI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(Planet.Ocean.sigmaIce_Sm['Ih'], 0,
                                                                               Planet.phi_frac[indsI],
                                                                               Planet.Ocean.Jsigma)
    if np.size(indsClath) != 0:
        Planet.sigma_Sm[indsClath] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(Planet.Ocean.sigmaIce_Sm['Clath'], 0,
                                                                                      Planet.phi_frac[indsClath],
                                                                                      Planet.Ocean.Jsigma)
    # We use the negative underplate phase IDs for dry HP ices
    if np.size(indsIIund) != 0:
        Planet.sigma_Sm[indsIIund] = Planet.Ocean.surfIceEOS['II'].fn_porosCorrect(Planet.Ocean.sigmaIce_Sm['II'], 0,
                                                                                   Planet.phi_frac[indsIIund],
                                                                                   Planet.Ocean.Jsigma)
    if np.size(indsIIIund) != 0:
        Planet.sigma_Sm[indsIIIund] = Planet.Ocean.surfIceEOS['III'].fn_porosCorrect(Planet.Ocean.sigmaIce_Sm['III'], 0,
                                                                                     Planet.phi_frac[indsIIIund],
                                                                                     Planet.Ocean.Jsigma)
    if np.size(indsVund) != 0:
        Planet.sigma_Sm[indsVund] = Planet.Ocean.surfIceEOS['V'].fn_porosCorrect(Planet.Ocean.sigmaIce_Sm['V'], 0,
                                                                                 Planet.phi_frac[indsVund],
                                                                                 Planet.Ocean.Jsigma)
    if np.size(indsVIund) != 0:
        Planet.sigma_Sm[indsVIund] = Planet.Ocean.surfIceEOS['VI'].fn_porosCorrect(Planet.Ocean.sigmaIce_Sm['VI'], 0,
                                                                                   Planet.phi_frac[indsVIund],
                                                                                   Planet.Ocean.Jsigma)
    # Next, do liquid-filled ice and clathrate pores
    if np.size(indsIwet) != 0:
        # First, get pore fluid conductivity
        sigmaFluid_Sm = Planet.Ocean.EOS.fn_sigma_Sm(Planet.Ppore_MPa[indsIwet], Planet.T_K[indsIwet])
        # Account for possible errors in pore fluid conductivity calcs
        validSigs = np.logical_not(np.isnan(sigmaFluid_Sm))
        # Interpolate over NaNs to remove them
        sigmaFluid_Sm = spi.griddata(Planet.r_m[indsIwet][validSigs], sigmaFluid_Sm[validSigs], Planet.r_m[indsIwet])
        # Finally, combine the conductivity of the porous ice with the pore-filling fluid
        Planet.sigma_Sm[indsIwet] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(Planet.Ocean.sigmaIce_Sm['Ih'], sigmaFluid_Sm,
                                                                                  Planet.phi_frac[indsIwet],
                                                                                  Planet.Ocean.Jsigma)
    if np.size(indsClathWet) != 0:
        # First, get pore fluid conductivity
        sigmaFluid_Sm = Planet.Ocean.EOS.fn_sigma_Sm(Planet.Ppore_MPa[indsClathWet], Planet.T_K[indsClathWet])
        # Account for possible errors in pore fluid conductivity calcs
        validSigs = np.logical_not(np.isnan(sigmaFluid_Sm))
        # Interpolate over NaNs to remove them
        sigmaFluid_Sm = spi.griddata(Planet.r_m[indsClathWet][validSigs], sigmaFluid_Sm[validSigs], Planet.r_m[indsClathWet])
        # Finally, combine the conductivity of the porous ice with the pore-filling fluid
        Planet.sigma_Sm[indsClathWet] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(Planet.Ocean.sigmaIce_Sm['Clath'], sigmaFluid_Sm,
                                                                                         Planet.phi_frac[indsClathWet],
                                                                                         Planet.Ocean.Jsigma)
    if np.size(indsII) != 0:
        # First, get pore fluid conductivity
        sigmaFluid_Sm = Planet.Ocean.EOS.fn_sigma_Sm(Planet.Ppore_MPa[indsII], Planet.T_K[indsII])
        # Account for possible errors in pore fluid conductivity calcs
        validSigs = np.logical_not(np.isnan(sigmaFluid_Sm))
        # Interpolate over NaNs to remove them
        sigmaFluid_Sm = spi.griddata(Planet.r_m[indsII][validSigs], sigmaFluid_Sm[validSigs], Planet.r_m[indsII])
        # Finally, combine the conductivity of the porous ice with the pore-filling fluid
        Planet.sigma_Sm[indsII] = Planet.Ocean.iceEOS['II'].fn_porosCorrect(Planet.Ocean.sigmaIce_Sm['II'], sigmaFluid_Sm,
                                                                                  Planet.phi_frac[indsII],
                                                                                  Planet.Ocean.Jsigma)
    if np.size(indsIII) != 0:
        # First, get pore fluid conductivity
        sigmaFluid_Sm = Planet.Ocean.EOS.fn_sigma_Sm(Planet.Ppore_MPa[indsIII], Planet.T_K[indsIII])
        # Account for possible errors in pore fluid conductivity calcs
        validSigs = np.logical_not(np.isnan(sigmaFluid_Sm))
        # Interpolate over NaNs to remove them
        sigmaFluid_Sm = spi.griddata(Planet.r_m[indsIII][validSigs], sigmaFluid_Sm[validSigs], Planet.r_m[indsIII])
        # Finally, combine the conductivity of the porous ice with the pore-filling fluid
        Planet.sigma_Sm[indsIII] = Planet.Ocean.iceEOS['III'].fn_porosCorrect(Planet.Ocean.sigmaIce_Sm['III'], sigmaFluid_Sm,
                                                                                  Planet.phi_frac[indsIII],
                                                                                  Planet.Ocean.Jsigma)
    if np.size(indsV) != 0:
        # First, get pore fluid conductivity
        sigmaFluid_Sm = Planet.Ocean.EOS.fn_sigma_Sm(Planet.Ppore_MPa[indsV], Planet.T_K[indsV])
        # Account for possible errors in pore fluid conductivity calcs
        validSigs = np.logical_not(np.isnan(sigmaFluid_Sm))
        # Interpolate over NaNs to remove them
        sigmaFluid_Sm = spi.griddata(Planet.r_m[indsV][validSigs], sigmaFluid_Sm[validSigs], Planet.r_m[indsV])
        # Finally, combine the conductivity of the porous ice with the pore-filling fluid
        Planet.sigma_Sm[indsV] = Planet.Ocean.iceEOS['V'].fn_porosCorrect(Planet.Ocean.sigmaIce_Sm['V'], sigmaFluid_Sm,
                                                                                  Planet.phi_frac[indsV],
                                                                                  Planet.Ocean.Jsigma)
    if np.size(indsVI) != 0:
        # First, get pore fluid conductivity
        sigmaFluid_Sm = Planet.Ocean.EOS.fn_sigma_Sm(Planet.Ppore_MPa[indsVI], Planet.T_K[indsVI])
        # Account for possible errors in pore fluid conductivity calcs
        validSigs = np.logical_not(np.isnan(sigmaFluid_Sm))
        # Interpolate over NaNs to remove them
        sigmaFluid_Sm = spi.griddata(Planet.r_m[indsVI][validSigs], sigmaFluid_Sm[validSigs], Planet.r_m[indsVI])
        # Finally, combine the conductivity of the porous ice with the pore-filling fluid
        Planet.sigma_Sm[indsVI] = Planet.Ocean.iceEOS['VI'].fn_porosCorrect(Planet.Ocean.sigmaIce_Sm['VI'], sigmaFluid_Sm,
                                                                                  Planet.phi_frac[indsVI],
                                                                                  Planet.Ocean.Jsigma)


    return Planet


def CalcElecSolidIce(Planet, Params, indsLiq, indsI, indsII, indsIIund, indsIII, indsIIIund,
                                     indsV, indsVund, indsVI, indsVIund, indsClath):
    """ Calculation of electrical conductivities for solid (non-porous) ice. All porosity
        considerations are ignored.

        Args:
            indsLiq (int, shape Ni): Indices of Planet layer arrays (e.g. Planet.sigma_Sm)
                where the layers are ocean fluid.
            indsI, indsIwet (int, shape Nj): Ice Ih indices.
            indsII, indsIIund (int, shape Nk): Ice II indices.
            indsIII, indsIIIund (int, shape Nl): Ice III indices.
            indsV, indsVund (int, shape Nm): Ice V indices.
            indsVI, indsVIund (int, shape Nn): Ice VI indices.
            indsClath, indsClathWet (int, shape No): Clathrate ice indices.
        Assigns Planet attributes:
            sigma_Sm
    """

    # Calculate and/or assign conductivities for each phase type
    if np.size(indsLiq) != 0:
        Planet.sigma_Sm[indsLiq] = Planet.Ocean.EOS.fn_sigma_Sm(Planet.P_MPa[indsLiq], Planet.T_K[indsLiq])

    Planet.sigma_Sm[indsI] = Planet.Ocean.sigmaIce_Sm['Ih']
    Planet.sigma_Sm[np.concatenate((indsII,  indsIIund))] =  Planet.Ocean.sigmaIce_Sm['II']
    Planet.sigma_Sm[np.concatenate((indsIII, indsIIIund))] = Planet.Ocean.sigmaIce_Sm['III']
    Planet.sigma_Sm[np.concatenate((indsV,   indsVund))] =   Planet.Ocean.sigmaIce_Sm['V']
    Planet.sigma_Sm[np.concatenate((indsVI,  indsVIund))] =  Planet.Ocean.sigmaIce_Sm['VI']
    Planet.sigma_Sm[indsClath] = Planet.Ocean.sigmaIce_Sm['Clath']

    return Planet


def CalcElecPorRock(Planet, Params, indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI):
    """ Calculate electrical conductivities of porous rock layers. Pores are assumed to be
        filled with ocean fluids or HP ices, based on feeding the EOS the pressure in the
        pore space and the temperature of the surrounding layer. Pores are assumed to be
        supported by the rock matrix, so the internal pressures are determined by overlying
        ocean/ice pore materials. The method for calculating conductivities is as follows:
        Apply Wong et al. (1984) model based on comparing an Archie's law formulation against
        measurements with glass beads: https://doi.org/10.1103/PhysRevB.30.6606
        This assumes the porosity dependence of conductivity has a sharp knee,
        and pore connectivity ~shuts off below some threshold (~5%).

        Args:
            indsSil (int, shape Planet.Steps.nSil): Indices of layer arrays (e.g. Planet.P_MPa) that
                correspond to silicates. Essentially a sorted concatenation of all other inputs.
            indsSilLiq (int, shape Nl): Indices of silicate layers that contain ocean liquid within
                the pore space, according to applying the ocean EOS to (Ppore_MPa, T_K).
            indsSilI, ... indsSilVI (int, shape Ni): Indices of ice phases found in silicate pores.
        Assigns Planet attributes:
            sigma_Sm
    """

    # First, get pore fluid conductivity
    sigmaFluid_Sm = Planet.Ocean.EOS.fn_sigma_Sm(Planet.Ppore_MPa[indsSilLiq], Planet.T_K[indsSilLiq])
    # Account for possible errors in pore fluid conductivity calcs
    validSigs = np.logical_not(np.isnan(sigmaFluid_Sm))
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
        sigmaPore_Sm, Planet.phi_frac[indsSil], Planet.Sil.Jsigma)

    # Record the mean conductivities of the pore fluids and porous layers,
    # but include only the contributing layers in the calculation
    if np.size(indsSilLiq) != 0:
        aboveThreshLiq = Planet.phi_frac[indsSilLiq] >= Planet.Sil.poreConductThresh_frac
        Planet.Sil.sigmaPoreMean_Sm = np.mean(sigmaFluid_Sm[aboveThreshLiq])
        Planet.Sil.sigmaPorousLayerMean_Sm = np.mean(Planet.sigma_Sm[indsSilLiq][aboveThreshLiq])
    else:
        # No pores contain liquid (only HP ices), so just record the mean combined conductivity
        Planet.Sil.sigmaPoreMean_Sm = np.nan
        Planet.Sil.sigmaPorousLayerMean_Sm = np.mean(Planet.sigma_Sm[indsSil][aboveThresh])

    return Planet
