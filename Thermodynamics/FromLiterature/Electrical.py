import numpy as np
import logging as log
from Thermodynamics.FromLiterature.HydroEOS import GetPhaseIndices
from Utilities.defineStructs import Constants

def ElecConduct(Planet, Params):
    """ Calculate/assign electrical conductivities for each layer

        Assigns Planet attributes:
            sigma_Sm
    """
    # Initialize outputs as NaN so that we get errors if we missed any layers
    Planet.sigma_Sm = np.zeros(Planet.Steps.nTotal) * np.nan

    if Params.CALC_CONDUCT:
        # Identify which indices correspond to which phases
        indsLiquid, indsI, indsII, indsIII, indsV, indsVI, indsClath, indsSil, indsFe = GetPhaseIndices(Planet.phase)

        # Calculate and/or assign conductivities for each phase type
        if np.size(indsI) != 0:
            if Planet.Do.POROUS_ICE:
                Planet.sigma_Sm[indsI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(Planet.Ocean.sigmaIce_Sm, 0,
                                                                        Planet.phi_frac[indsI], Planet.Ocean.Jsigma)
            else:
                Planet.sigma_Sm[indsI] = Planet.Ocean.sigmaIce_Sm
        if np.size(indsLiquid) != 0:
            Planet.sigma_Sm[indsLiquid] = Planet.Ocean.EOS.fn_sigma_Sm(Planet.P_MPa[indsLiquid], Planet.T_K[indsLiquid])
        if np.size(indsII) != 0:
            Planet.sigma_Sm[indsII] = Planet.Ocean.sigmaIce_Sm
        if np.size(indsIII) != 0:
            Planet.sigma_Sm[indsIII] = Planet.Ocean.sigmaIce_Sm
        if np.size(indsV) != 0:
            Planet.sigma_Sm[indsV] = Planet.Ocean.sigmaIce_Sm
        if np.size(indsVI) != 0:
            Planet.sigma_Sm[indsVI] = Planet.Ocean.sigmaIce_Sm
        if np.size(indsClath) != 0:
            if Planet.Do.POROUS_ICE:
                Planet.sigma_Sm[indsClath] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(Constants.sigmaClath_Sm, 0,
                                                                        Planet.phi_frac[indsClath], Planet.Ocean.Jsigma)
            else:
                Planet.sigma_Sm[indsClath] = Constants.sigmaClath_Sm
        if np.size(indsSil) != 0:
            if Planet.Do.POROUS_ROCK and not Params.SKIP_INNER:
                Planet.sigma_Sm[indsSil] = Planet.Sil.EOS.fn_porosCorrect(Planet.Sil.sigmaSil_Sm,
                    Planet.Ocean.EOS.fn_sigma_Sm(Planet.Ppore_MPa[indsSil], Planet.T_K[indsSil]),
                    Planet.phi_frac[indsSil], Planet.Sil.Jsigma)
            else:
                Planet.sigma_Sm[indsSil] = Planet.Sil.sigmaSil_Sm
        if np.size(indsFe) != 0:
            Planet.sigma_Sm[indsFe] = Planet.Core.sigmaCore_Sm

        if np.any(np.isnan(Planet.sigma_Sm)):
            raise ValueError('Some layer conductivities are NaN. There is likely an error with layer phase assignments.')

    return Planet
