import numpy as np
from alma import infer_rheology_pp, build_model, love_numbers
import logging
import ast
import os
from PlanetProfile.Utilities.Indexing import PhaseConv

# Assign logger
log = logging.getLogger('PlanetProfile')

def GravityParameters(Planet, Params):
    """ Calculate induced gravity responses for the body and prints them to disk."""
    if Planet.Do.VALID and Params.CALC_NEW_GRAVITY and Params.CALC_VISCOSITY and Params.CALC_SEISMIC and not Params.SKIP_INNER:
        # Set Magnetic struct layer arrays as we need for induction calculations
        Planet, Params = SetupGravity(Planet, Params)
        rheology, params = infer_rheology_pp(Planet.Gravity.ALMAModel, structure=Params.Gravity.rheology_structure,
                                             layer_radius=Params.Gravity.layer_radius,
                                             layer_radius_index=Params.Gravity.layer_radius_index)
        model_params = build_model(Planet.Gravity.ALMAModel['model'][:, Planet.Gravity.ALMAModel['columns'].index('r')],
                                   Planet.Gravity.ALMAModel['model'][:, Planet.Gravity.ALMAModel['columns'].index('rho')],
                                   Planet.Gravity.ALMAModel['mu'],
                                   Planet.Gravity.ALMAModel['vis'],
                                   rheology, params, ndigits=Params.Gravity.num_digits, verbose=Params.Gravity.verbose)
        # Compute love numbers
        Planet.Gravity.h, Planet.Gravity.l, Planet.Gravity.k = love_numbers(Params.Gravity.harmonic_degrees, Params.Gravity.time_log_kyrs,
                                                                Params.Gravity.loading_type, Params.Gravity.time_history_function, Params.Gravity.tau, model_params,
                                                                Params.Gravity.output_type, Params.Gravity.gorder, verbose=Params.Gravity.verbose,
                                                                parallel=Params.Gravity.parallel and not (Params.INDUCTOGRAM_IN_PROGRESS or Params.DO_EXPLOREOGRAM))
        # Compute delta relation
        Planet.Gravity.delta = 1 + Planet.Gravity.k - Planet.Gravity.h
        # If our love numbers are 1x1 numpy array, let's convert to float - important for plotting and output
        if len(Planet.Gravity.time_log_kyrs) == 1 and len(Planet.Gravity.harmonic_degrees) == 1:
            Planet.Gravity.h = float(Planet.Gravity.h[0, 0])
            Planet.Gravity.l = float(Planet.Gravity.l[0, 0])
            Planet.Gravity.k = float(Planet.Gravity.k[0, 0])
            Planet.Gravity.delta = float(Planet.Gravity.delta[0, 0])
        if (not Params.NO_SAVEFILE) and (not Params.INVERSION_IN_PROGRESS) and (not Params.DO_EXPLOREOGRAM):
            Planet, Params = WriteGravityParameters(Planet, Params)
    elif Planet.Do.VALID:
        if os.path.isfile(Params.DataFiles.gravityParametersFile):
            # Reload gravity parameters from disk
            Planet, Params = ReloadGravityParameters(Planet, Params)
        else:
            log.warning(
                f'CALC_NEW_GRAVITY is False, but {Params.DataFiles.gravityParametersFile} was not found. ' + f'Skipping gravity parameter calculations.')
    return Planet, Params




def SetupGravity(Planet, Params):
    """ Reconfigure layer boundaries and gravity model into a format
        usable by gravity response calculation functions of PyALMA3.

        NOTE: Function models the read_model_pp of PyALMA3 so that we can ensure compatibility with the package's functions.

        Requires Planet attributes:
    """
    if Params.CALC_NEW_GRAVITY and Planet.Do.VALID and not Params.SKIP_INNER:
        """Combine data into model format that is required by PyALMA3"""
        # Note we have to use r_m[:-1] since r_m has one extra value than other arrays
        Planet.Gravity.model = np.vstack(
            [Planet.Reduced.r_m, Planet.Reduced.phase, Planet.Reduced.rho_kgm3, Planet.Reduced.Seismic.VP_kms, Planet.Reduced.Seismic.VS_kms, Planet.Reduced.Seismic.GS_GPa, Planet.Reduced.eta_Pas]).T
    
        # Convert parameter units to Pa and meters
        for index, (header, unit) in enumerate(zip(Planet.Gravity.columns, Planet.Gravity.units_PyALMA3)):
            if header in Planet.Gravity.parameters_to_convert:
                Planet.Gravity.model[:, index] = Planet.Gravity.model[:, index] * Planet.Gravity.parameters_to_convert[header]
        # Get indices of used properties
        rIndex  = Planet.Gravity.columns.index('r')
        pIndex = Planet.Gravity.columns.index('phase')
        rhoIndex = Planet.Gravity.columns.index('rho')
        VPIndex = Planet.Gravity.columns.index('VP')
        VSIndex = Planet.Gravity.columns.index('VS')
        GSIndex = Planet.Gravity.columns.index('GS')
        etaIndex = Planet.Gravity.columns.index('eta')
        # Round radius - In similar PyALMA3 function, it rounds radius so will do so here as well
        Planet.Gravity.model[:, rIndex] = np.round(Planet.Gravity.model[:, rIndex], 0)

        # Flip model: core at top and surface at bottom
        Planet.Gravity.model = np.flipud(Planet.Gravity.model)

        ## Calculate elastic parameters
        # LAMBDA = rho (Vp^2 - 2 * Vs^2)
        # 1st Lame parameter
        Planet.Gravity.LAMBDA_Pa = Planet.Gravity.model[:, rhoIndex] * (
                np.power(Planet.Gravity.model[:, VPIndex], 2) -
                2. * np.power(Planet.Gravity.model[:, VSIndex], 2))

        # shear modulus G or MU = rho Vs^2
        Planet.Gravity.MU_Pa = Planet.Gravity.model[:, GSIndex]

        # Poissons ratio sigma = lambda / 2*(lambda + mu)
        Planet.Gravity.SIGMA = Planet.Gravity.LAMBDA_Pa / (2 * Planet.Gravity.LAMBDA_Pa + 2 * Planet.Gravity.MU_Pa)

        # Youngs modulus Y = 2 * MU * (1 + sigma)
        Planet.Gravity.Y_Pa = 2. * Planet.Gravity.MU_Pa * (1 + Planet.Gravity.SIGMA)

        Planet.Gravity.VISCOSITY_kgms = Planet.Gravity.model[:, etaIndex]

        Planet.Gravity.ALMAModel = {'columns': Planet.Gravity.columns,
                'units': Planet.Gravity.units_PyALMA3,
                'model': Planet.Gravity.model,
                'lambda': Planet.Gravity.LAMBDA_Pa,
                'mu': Planet.Gravity.MU_Pa,
                'sigma': Planet.Gravity.SIGMA,
                'y': Planet.Gravity.Y_Pa,
                'vis': Planet.Gravity.VISCOSITY_kgms}

        # Set Planet time scale and harmonic degrees from Params
        Planet.Gravity.time_log_kyrs = Params.Gravity.time_log_kyrs
        Planet.Gravity.harmonic_degrees = Params.Gravity.harmonic_degrees

        # Finally, we must setup the rheology structure, from core to surface
        # Track phase changes while preserving order
        phases = Planet.Gravity.model[:, pIndex]
        # Initialize with the first phase
        phase_transitions = [phases[0]]
        # Iterate through the model to find phase transitions
        for i in range(1, len(phases)):
            if phases[i] != phases[i-1]:
                # We've found a phase transition
                phase_transitions.append(phases[i])
        # Map phases to rheology models
        rheology_structure = []
        for phase in phase_transitions:
            # Convert numerical phase to string representation for dictionary lookup
            phase_str = PhaseConv(phase, liq='0')
        
            # Try to use PhaseConv if direct lookup fails
            if phase_str not in Params.Gravity.rheology_models:
                raise ValueError(f"Phase {phase_str} not found in rheology models.")
            else:
                rheology_model = Params.Gravity.rheology_models[phase_str]
            rheology_structure.append(rheology_model)
        
        # Store the compiled structures in Planet.Gravity
        Params.Gravity.rheology_structure = rheology_structure
        
    # Return Planet and Params
    return Planet, Params



def WriteGravityParameters(Planet, Params):
    """Write the gravity parameters to disk"""
    headerLines = [
        f'Time range (log kyrs) = {Planet.Gravity.time_log_kyrs}',
        f'harmonic degrees = {Planet.Gravity.harmonic_degrees}',
        f'h love number = {Planet.Gravity.h}',
        f'l love number = {Planet.Gravity.l}',
        f'k love number = {Planet.Gravity.k}'
        ]
    # Write out data from core/mantle trade
    with open(Params.DataFiles.gravityParametersFile, 'w') as f:
        f.write('\n  '.join(headerLines) + '\n')
    log.debug(f'Saved induced moments to file: {Params.DataFiles.gravityParametersFile}')

    return Planet, Params


def ReloadGravityParameters(Planet, Params):
    """Reload gravity paramters from disk"""
    with open(Params.DataFiles.gravityParametersFile) as f:
        Planet.Gravity.time_log_kyrs = ast.literal_eval(f.readline().split('=')[-1])
        Planet.Gravity.harmonic_degrees = ast.literal_eval(f.readline().split('=')[-1])
        if len(Planet.Gravity.time_log_kyrs) == 1 and len(Planet.Gravity.harmonic_degrees) == 1:
            # Make these variables floats
            Planet.Gravity.h = float((f.readline().split('=')[-1]))
            Planet.Gravity.l = float((f.readline().split('=')[-1]))
            Planet.Gravity.k = float((f.readline().split('=')[-1]))
        else:
            Planet.Gravity.h = np.array(ast.literal_eval(f.readline().split('=')[-1]))
            Planet.Gravity.l = np.array(ast.literal_eval(f.readline().split('=')[-1]))
            Planet.Gravity.k = np.array(ast.literal_eval(f.readline().split('=')[-1]))
        Planet.Gravity.delta = 1 + Planet.Gravity.k - Planet.Gravity.h
    return Planet, Params
