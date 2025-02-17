""" File to reduce the plnaet profile based on parameters. Currently reduced based on magnetic induction calculations"""
import os
import numpy as np
import logging
import scipy.interpolate as spi
from PlanetProfile.Utilities.defineStructs import Constants
from PlanetProfile.Utilities.Indexing import PhaseConv

# Assign logger
log = logging.getLogger('PlanetProfile')

def GetReducedPlanetProfile(Planet, Params):
    """ Generate a reduced planet profile to be used in magnetic induction and/or gravity calculations
    #TODO Only implemented for magnetic induction currently"""
    if Params.REDUCE_ACCORDING_TO == 'MagneticInduction':
        Planet, Params = GetMagneticReducedLayers(Planet, Params)
        # We need to ensure we calculated seismic to be able to interpolate these values
        if not Params.SKIP_GRAVITY and Params.CALC_SEISMIC:
            Planet.Reduced.phase = spi.interp1d(Planet.r_m[:-1], Planet.phase, kind=Params.Induct.oceanInterpMethod,
                                                 bounds_error=False)(Planet.Reduced.r_m)
            Planet.Reduced.rho_kgm3 = spi.interp1d(Planet.r_m[:-1], Planet.rho_kgm3, kind=Params.Induct.oceanInterpMethod,
                                                 bounds_error=False)(Planet.Reduced.r_m)
            Planet.Reduced.Seismic.VS_kms = spi.interp1d(Planet.r_m[:-1], Planet.Seismic.VS_kms,
                                                    kind=Params.Induct.oceanInterpMethod, bounds_error=False)(Planet.Reduced.r_m)
            Planet.Reduced.Seismic.VP_kms = spi.interp1d(Planet.r_m[:-1], Planet.Seismic.VP_kms,
                                                  kind=Params.Induct.oceanInterpMethod, bounds_error=False)(Planet.Reduced.r_m)
            Planet.Reduced.Seismic.GS_GPa = spi.interp1d(Planet.r_m[:-1], Planet.Seismic.GS_GPa,
                                                   kind=Params.Induct.oceanInterpMethod, bounds_error=False)(Planet.Reduced.r_m)
            Planet.Reduced.eta_Pas = spi.interp1d(Planet.r_m[:-1], Planet.eta_Pas,
                                                  kind=Params.Induct.oceanInterpMethod, bounds_error=False)(
                Planet.Reduced.r_m)
    elif Params.REDUCE_ACCORDING_TO == 'AveragedLayers':
        Planet, Params = GetAverageLayers(Planet, Params)
    elif Params.REDUCE_ACCORDING_TO == 'ReducedLayers':
        Planet, Params = GetReducedLayers(Planet, Params)
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
                Planet.Reduced.r_m = rLayers_m[iChange]
                Planet.Reduced.sigma_Sm = sigmaInduct_Sm[iChange]
    return Planet, Params


def GetReducedLayers(Planet, Params):
    """
    Get the reduced planet according to Params.reduceLayersAccordingTo
    """
    phase = Planet.phase.copy()  # Create a copy to avoid modifying the original attribute

    # Find change points and include start/end
    change_indices = np.concatenate(([0], np.where(np.diff(phase) != 0)[0] + 1, [len(phase)]))

    rPhase_m = Planet.r_m[:-1]

    # Define attributes to interpolate (separating seismic attributes)
    attributes_to_reduce = {
        "phase": Planet,
        "r_m": Planet,
        "sigma_Sm": Planet,
        "rho_kgm3": Planet,
        "VS_kms": Planet.Seismic,
        "VP_kms": Planet.Seismic,
        "GS_GPa": Planet.Seismic,
        "eta_Pas": Planet
    }

    # Initialize attributes in the correct substructure
    for attr, source in attributes_to_reduce.items():
        if source is Planet.Seismic:
            setattr(Planet.Reduced.Seismic, attr, [])  # Store in Reduced.Seismic
        else:
            setattr(Planet.Reduced, attr, [])

    # Iterate through each segment
    for start, end in zip(change_indices[:-1], change_indices[1:]):
        layer_phase = phase[start]
        layer_str = PhaseConv(layer_phase, liq='0')  # Convert phase to string
        r_layer_m = rPhase_m[start:end]  # Extract radii for this layer

        target_layers = min(Params.REDUCED_LAYERS_SIZE[layer_str], len(r_layer_m))  # Get target layer count

        # Generate original and reduced depth points
        original_points = r_layer_m
        reduced_points = np.linspace(r_layer_m[0], r_layer_m[-1], target_layers)

        for attr, source in attributes_to_reduce.items():
            original_values = getattr(source, attr)[start:end]  # Extract original data
            interpolator = spi.interp1d(original_points, original_values, kind='linear', fill_value="extrapolate")
            reduced_values = interpolator(reduced_points)  # Interpolated values

            # Store reduced values in the correct location
            if source is Planet.Seismic:
                getattr(Planet.Reduced.Seismic, attr).extend(reduced_values.tolist())  # Store in Seismic
            else:
                getattr(Planet.Reduced, attr).extend(reduced_values.tolist())  # Convert all attributes in Planet.Reduced and Planet.Reduced.Seismic to NumPy arrays
    for attr, source in attributes_to_reduce.items():
        if source is Planet.Seismic:
            setattr(Planet.Reduced.Seismic, attr, np.array(getattr(Planet.Reduced.Seismic, attr)))
        else:
            setattr(Planet.Reduced, attr, np.array(getattr(Planet.Reduced, attr)))

    return Planet, Params


def GetAverageLayers(Planet, Params):
    """Reduce planet into each of its differentiated layers, taking the average of all of its properties"""
    iChange = np.where(Planet.phase[:-1] != Planet.phase[1:])[0] + 1
    layerChangeIndices = np.insert(iChange, 0, 0)
    if Planet.DO.VALID:

        r_phases = []
        indsLiq = np.where((Planet.phase) == 0)[0]
        for i in range(iChange.size):
            layerIndices = layerChangeIndices[i, i+1]
            reduced_layer = Planet.phase[layerIndices]
            # Obtain radius and induction of layers
            rLayers_m = Planet.r_m[layerIndices]
            sigmaInduct_Sm = Planet.sigma_Sm[layerIndices]
            rho_kgm3 = Planet.rho_kgm3[layerIndices]
            vs_kms = Planet.Seismic.VS_kms[layerIndices]
            vp_kms = Planet.Seismic.VP_kms[layerIndices]
            gs_GPa = Planet.Seismic.GS_GPa[layerIndices]
            reta_Pas = Planet.eta_Pas[layerIndices]

            if reduced_layer[0] == 0:
                if Params.Sig.REDUCED_INDUCT:
                    # Get radius values from D/nIntL above the seafloor to the ice shell
                    rBot_m = Planet.Bulk.R_m - (Planet.zb_km + Planet.D_km) * 1e3
                    rTop_m = rLayers_m[indsLiq[0]]
                    rOcean_m = np.linspace(rTop_m, rBot_m, Params.Induct.nIntL + 1)[:-1]
                    # Interpolate the conductivities corresponding to those radii
                    if np.size(indsLiq) == 1:
                        log.warning(
                            f'Only 1 layer found in ocean, but number of layers to ' + f'interpolate over is {Params.Induct.nIntL}. Arbitrary layers will be introduced.')
                        rModel_m = np.concatenate((np.array([rBot_m]), rLayers_m[indsLiq]))
                    else:
                        rModel_m = rLayers_m[indsLiq]
                    # Stitch together the r and sigma arrays with the new ocean values
                    rLayers_m = np.concatenate((rLayers_m[:indsLiq[0]], rOcean_m, rLayers_m[indsLiq[-1] + 1:]))





