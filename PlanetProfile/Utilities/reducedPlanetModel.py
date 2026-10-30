""" File to reduce the plnaet profile based on parameters. Currently reduced based on magnetic induction calculations"""
import os
import numpy as np
import logging
import scipy.interpolate as spi
from PlanetProfile.Utilities.defineStructs import Constants
from PlanetProfile.Utilities.Indexing import PhaseConv

# Assign logger
log = logging.getLogger('PlanetProfile')

def GetReducedPlanet(Planet, Params):
    """ Generate a reduced planet profile to be used in magnetic induction and/or gravity calculations
    #TODO Only implemented for reduced layer currently, and not used for magnetic induction calculations until it can be further stress tested"""
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
    # Find where phase changes occur
    phase_transitions = np.diff(Planet.phase) != 0
    layer_boundaries = phase_transitions.copy()
    
    # Also include convection boundaries if ice convection is enabled
    if not Planet.Do.NO_ICE_CONVECTION:
        # Find phase transitions within the ice shell (surface ice layers only)
        ice_phase_transitions = np.diff(Planet.phase[:Planet.Steps.nSurfIce]) != 0
        
        # Find where convection state changes (convective vs conductive layers)
        convection_transitions = np.diff(Planet.Steps.iConv) != 0
        
        # Combine both types of boundaries in ice - either phase change OR convection change
        ice_layer_boundaries = np.logical_or(ice_phase_transitions, convection_transitions)
        
        # Update the overall layer boundaries array for ice shell portion only
        # (nSurfIce-1 because diff reduces array size by 1)
        layer_boundaries[:Planet.Steps.nSurfIce-1] = ice_layer_boundaries
    
    # Get the indices where layer boundaries occur (add 1 to account for diff offset)
    boundary_indices = np.where(layer_boundaries)[0] + 1

    # Include start and end
    change_indices = np.concatenate(([0], boundary_indices, [len(Planet.phase)]))

    rPhase_m = Planet.r_m[:-1]

    # Pre-allocate lists for reduced attributes - more efficient than repeated extend
    reduced_phase = []
    reduced_r_m = []
    reduced_rho_kgm3 = []
    reduced_changeIndices = [0]
    reduced_VP_kms = [] if Params.CALC_SEISMIC else None
    reduced_VS_kms = [] if Params.CALC_SEISMIC else None
    reduced_GS_GPa = [] if Params.CALC_SEISMIC else None
    reduced_eta_Pas = [] if Params.CALC_VISCOSITY else None
    reduced_sigma_Sm = [] if Params.CALC_CONDUCT else None

    # Initialize reduced attributes
    Planet.Reduced.phase = []
    Planet.Reduced.r_m = []
    Planet.Reduced.rho_kgm3 = []
    Planet.Reduced.iConv = []
    
    if Params.CALC_SEISMIC:
        Planet.Reduced.Seismic.VP_kms = []
        Planet.Reduced.Seismic.VS_kms = []
        Planet.Reduced.Seismic.GS_GPa = []
    
    if Params.CALC_VISCOSITY:
        Planet.Reduced.eta_Pas = []
    
    if Params.CALC_CONDUCT:
        Planet.Reduced.sigma_Sm = []

    # Iterate through each segment
    for start, end in zip(change_indices[:-1], change_indices[1:]):
        layer_phase = Planet.phase[start]
        
        if not Planet.Do.NO_ICE_CONVECTION and start < Planet.Steps.nSurfIce:
            convection = Planet.Steps.iConv[start]
        else:
            convection = False
        
            
        # Skip inner layers if requested
        if layer_phase >= 10 and Params.SKIP_INNER:
            continue
            
        layer_str = PhaseConv(layer_phase, PORE=Planet.Do.POROUS_ROCK, liq='0')  # Convert phase to string
        r_layer_m = rPhase_m[start:end]  # Extract radii for this layer
        n_layers = len(r_layer_m)

        
        # Get target layer count
        if layer_str in Params.REDUCED_LAYERS_SIZE:
            target_layers = min(Params.REDUCED_LAYERS_SIZE[layer_str], n_layers)
        else:
            target_layers = min(Constants.defaultReducedLayerSize, n_layers)
        # Determine convection state for this layer segment
        # Check if we're in a surface ice layer and convection modeling is enabled
        if not Planet.Do.NO_ICE_CONVECTION and start < Planet.Steps.nSurfIce:
            # Use the convection state from the original model at this layer's starting index
            convection = Planet.Steps.iConv[start]
        else:
            # No convection for non-ice layers or when convection modeling is disabled
            convection = False

        # Apply the same convection state to all reduced layers in this segment
        # This maintains consistency within each homogeneous layer
        Planet.Reduced.iConv.extend([convection] * target_layers)
        
        # Only do interpolation if we actually need to reduce layers
        if n_layers > target_layers:
            # Interpolate to reduce number of layers
            original_points = r_layer_m
            reduced_points = np.linspace(r_layer_m[0], r_layer_m[-1], target_layers)
            
            # Direct attribute access for interpolation - much faster than getattr/setattr
            # Core attributes (always present)
            phase_interp = spi.interp1d(original_points, Planet.phase[start:end], kind='nearest')
            r_m_interp = spi.interp1d(original_points, Planet.r_m[start:end], kind='linear')
            rho_interp = spi.interp1d(original_points, Planet.rho_kgm3[start:end], kind='linear')
            
            reduced_phase.extend(phase_interp(reduced_points).tolist())
            reduced_r_m.extend(r_m_interp(reduced_points).tolist())
            reduced_rho_kgm3.extend(rho_interp(reduced_points).tolist())
            
            # Optional attributes - only interpolate if needed
            if Params.CALC_SEISMIC:
                vp_interp = spi.interp1d(original_points, Planet.Seismic.VP_kms[start:end], kind='linear')
                vs_interp = spi.interp1d(original_points, Planet.Seismic.VS_kms[start:end], kind='linear')
                gs_interp = spi.interp1d(original_points, Planet.Seismic.GS_GPa[start:end], kind='linear')
                
                reduced_VP_kms.extend(vp_interp(reduced_points).tolist())
                reduced_VS_kms.extend(vs_interp(reduced_points).tolist())
                reduced_GS_GPa.extend(gs_interp(reduced_points).tolist())
            
            if Params.CALC_VISCOSITY:
                eta_interp = spi.interp1d(original_points, Planet.eta_Pas[start:end], kind='linear')
                reduced_eta_Pas.extend(eta_interp(reduced_points).tolist())
            
            if Params.CALC_CONDUCT:
                sigma_interp = spi.interp1d(original_points, Planet.sigma_Sm[start:end], kind='linear')
                reduced_sigma_Sm.extend(sigma_interp(reduced_points).tolist())
                
        else:
            # Already reduced enough, just copy the original values - no interpolation needed
            reduced_phase.extend(Planet.phase[start:end].tolist())
            reduced_r_m.extend(Planet.r_m[start:end].tolist())
            reduced_rho_kgm3.extend(Planet.rho_kgm3[start:end].tolist())
            
            if Params.CALC_SEISMIC:
                reduced_VP_kms.extend(Planet.Seismic.VP_kms[start:end].tolist())
                reduced_VS_kms.extend(Planet.Seismic.VS_kms[start:end].tolist())
                reduced_GS_GPa.extend(Planet.Seismic.GS_GPa[start:end].tolist())
            
            if Params.CALC_VISCOSITY:
                reduced_eta_Pas.extend(Planet.eta_Pas[start:end].tolist())
            
            if Params.CALC_CONDUCT:
                reduced_sigma_Sm.extend(Planet.sigma_Sm[start:end].tolist())

        # Set the reduced change index
        reduced_changeIndices.extend([len(reduced_phase)])
    # Assign all reduced attributes at once - more efficient than individual assignments
    Planet.Reduced.phase = reduced_phase
    Planet.Reduced.r_m = reduced_r_m
    Planet.Reduced.rho_kgm3 = reduced_rho_kgm3
    Planet.Reduced.changeIndices = reduced_changeIndices
    
    if Params.CALC_SEISMIC:
        Planet.Reduced.Seismic.VP_kms = reduced_VP_kms
        Planet.Reduced.Seismic.VS_kms = reduced_VS_kms
        Planet.Reduced.Seismic.GS_GPa = reduced_GS_GPa
    
    if Params.CALC_VISCOSITY:
        Planet.Reduced.eta_Pas = reduced_eta_Pas
    
    if Params.CALC_CONDUCT:
        Planet.Reduced.sigma_Sm = reduced_sigma_Sm

    return Planet, Params
