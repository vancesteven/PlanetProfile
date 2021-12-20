import numpy as np
from Thermodynamics.HydroEOS import GetPhaseIndices

def ElecConduct(Planet, Params):
    """ Calculate/assign electrical conductivities for each layer

        Assigns Planet attributes:
            sigma_Sm
    """
    # Initialize outputs as NaN so that we get errors if we missed any odd phases
    sigma_Sm = np.zeros(Planet.Steps.nTotal) * np.nan
    # Identify which indices correspond to which phases
    indsLiquid, indsI, indsII, indsIII, indsV, indsVI, indsClath, indsSil, indsFe = GetPhaseIndices(Planet.phase)

    # Calculate and/or assign conductivities for each phase type
    if len(indsI) != 0:
        sigma_Sm[indsI] = Planet.Ocean.sigmaIce_Sm
    if len(indsLiquid) != 0:
        print('WARNING: Electrical conductivity calculations in the ocean are not yet implemented.')
        sigma_Sm[indsLiquid] = 0
    if len(indsII) != 0:
        sigma_Sm[indsII] = Planet.Ocean.sigmaIce_Sm
    if len(indsIII) != 0:
        sigma_Sm[indsIII] = Planet.Ocean.sigmaIce_Sm
    if len(indsV) != 0:
        sigma_Sm[indsV] = Planet.Ocean.sigmaIce_Sm
    if len(indsVI) != 0:
        sigma_Sm[indsVI] = Planet.Ocean.sigmaIce_Sm
    if len(indsClath) != 0:
        sigma_Sm[indsClath] = Planet.Ocean.sigmaIce_Sm
    if len(indsSil) != 0:
        sigma_Sm[indsSil] = Planet.Sil.sigmaSil_Sm
    if len(indsFe) != 0:
        sigma_Sm[indsFe] = Planet.Core.sigmaCore_Sm

    if np.any(np.isnan(sigma_Sm)):
        raise ValueError('Some layer conductivities are NaN. There is likely an error with layer phase assignments.')
    Planet.sigma_Sm = sigma_Sm

    return Planet