""" File to reduce the plnaet profile based on parameters. Currently reduced based on magnetic induction calculations"""
import os
import numpy as np
import logging
import scipy.interpolate as spi
from PlanetProfile.Utilities.defineStructs import Constants

# Assign logger
log = logging.getLogger('PlanetProfile')

def GetReducedPlanetProfile(Planet, Params):
    """ Generate a reduced planet profile to be used in magnetic induction and/or gravity calculations
    #TODO Only implemented for magnetic induction currently"""
    if Params.REDUCE_ACCORDING_TO == 'MagneticInduction':
        Planet, Params = GetMagneticReducedLayers(Planet, Params)
        # We need to ensure we calculated seismic to be able to interpolate these values
        if not Params.SKIP_GRAVITY and Params.CALC_SEISMIC:
            Planet.Reduced.rPhase = spi.interp1d(Planet.r_m[:-1], Planet.phase, kind=Params.Induct.oceanInterpMethod,
                                                 bounds_error=False)(Planet.Reduced.rLayers_m)
            Planet.Reduced.rRho_kgm3 = spi.interp1d(Planet.r_m[:-1], Planet.rho_kgm3, kind=Params.Induct.oceanInterpMethod,
                                                 bounds_error=False)(Planet.Reduced.rLayers_m)
            Planet.Reduced.rVS_kms = spi.interp1d(Planet.r_m[:-1], Planet.Seismic.VS_kms,
                                                    kind=Params.Induct.oceanInterpMethod, bounds_error=False)(Planet.Reduced.rLayers_m)
            Planet.Reduced.rVP_kms = spi.interp1d(Planet.r_m[:-1], Planet.Seismic.VP_kms,
                                                  kind=Params.Induct.oceanInterpMethod, bounds_error=False)(Planet.Reduced.rLayers_m)
            Planet.Reduced.rGS_GPa= spi.interp1d(Planet.r_m[:-1], Planet.Seismic.GS_GPa,
                                                   kind=Params.Induct.oceanInterpMethod, bounds_error=False)(Planet.Reduced.rLayers_m)
            Planet.Reduced.reta_Pas = spi.interp1d(Planet.r_m[:-1], Planet.eta_Pas,
                                                  kind=Params.Induct.oceanInterpMethod, bounds_error=False)(
                Planet.Reduced.rLayers_m)
    return Planet, Params


def GetMagneticReducedLayers(Planet, Params):
    """
    Get the reduced layers based on magnetic calculations to reduce computational overhead
    """
    # Lots of errors can happen if we haven't calculated the electrical conductivity,
    # so we make this contingent on having it.
    indsLiq = np.where((Planet.phase) == 0)[0]
    if Params.CALC_CONDUCT and Planet.Do.VALID:
        # Reconfigure layer conducting boundaries as needed.
        # For inductOtype == 'sigma', we have already set these arrays.
        if not (Params.Induct.inductOtype == 'sigma' and Params.DO_INDUCTOGRAM):
            if Params.CALC_NEW or not Params.DO_INDUCTOGRAM:
                # Obtain radius and induction of layers
                rLayers_m = Planet.r_m[:-1]
                sigmaInduct_Sm = Planet.sigma_Sm
                # Obtain other properties for reduced planet that are necessary
                # Eliminate NaN values and 0 values, assigning them to a default minimum
                sigmaInduct_Sm[np.logical_or(np.isnan(sigmaInduct_Sm), sigmaInduct_Sm == 0)] = Constants.sigmaDef_Sm
                # Set low conductivities to all be the same default value so we can shrink them down to single layers
                sigmaInduct_Sm[sigmaInduct_Sm < Constants.sigmaMin_Sm] = Constants.sigmaDef_Sm

                # Optionally, further reduce computational overhead by shrinking the number of ocean layers modeled
                if np.size(indsLiq) != 0 and Params.Sig.REDUCED_INDUCT and not Planet.Do.NO_H2O:
                    if not np.all(np.diff(indsLiq) == 1):
                        log.warning(
                            'HP ices found in ocean while REDUCED_INDUCT is True. They will be ignored ' + 'in the interpolation.')

                    # Get radius values from D/nIntL above the seafloor to the ice shell
                    rBot_m = Planet.Bulk.R_m - (Planet.zb_km + Planet.D_km) * 1e3
                    rTop_m = rLayers_m[indsLiq[0]]
                    rOcean_m = np.linspace(rTop_m, rBot_m, Params.Induct.nIntL + 1)[:-1]

                    # Interpolate the conductivities corresponding to those radii
                    if np.size(indsLiq) == 1:
                        log.warning(
                            f'Only 1 layer found in ocean, but number of layers to ' + f'interpolate over is {Params.Induct.nIntL}. Arbitrary layers will be introduced.')
                        rModel_m = np.concatenate((np.array([rBot_m]), rLayers_m[indsLiq]))
                        sigmaModel_Sm = np.concatenate((sigmaInduct_Sm[indsLiq] * 1.001, sigmaInduct_Sm[indsLiq]))
                    else:
                        rModel_m = rLayers_m[indsLiq]
                        sigmaModel_Sm = sigmaInduct_Sm[indsLiq]
                    sigmaOcean_Sm = spi.interp1d(rModel_m, sigmaModel_Sm, kind=Params.Induct.oceanInterpMethod,
                                                 bounds_error=False, fill_value=Constants.sigmaDef_Sm)(rOcean_m)
                    # Stitch together the r and sigma arrays with the new ocean values
                    rLayers_m = np.concatenate((rLayers_m[:indsLiq[0]], rOcean_m, rLayers_m[indsLiq[-1] + 1:]))
                    sigmaInduct_Sm = np.concatenate(
                        (sigmaInduct_Sm[:indsLiq[0]], sigmaOcean_Sm, sigmaInduct_Sm[indsLiq[-1] + 1:]))

                # Get the indices of layers just below where changes happen
                iChange = [i for i, sig in enumerate(sigmaInduct_Sm) if sig != np.append(sigmaInduct_Sm, np.nan)[i + 1]]
                Planet.Reduced.rLayers_m = rLayers_m[iChange]
                Planet.Reduced.rSigma_Sm = sigmaInduct_Sm[iChange]
    return Planet, Params