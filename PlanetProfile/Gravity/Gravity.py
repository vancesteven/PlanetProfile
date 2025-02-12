import numpy as np
from alma import infer_rheology_pp, build_model, love_numbers

def GravityParameters(Planet, Params):
    """ Calculate induced gravity responses for the body and prints them to disk."""
    SKIP = False
    if Planet.Do.VALID and Params.CALC_NEW_GRAVITY and Params.CALC_VISCOSITY:
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
                                                                parallel=Params.Gravity.parallel)
        # TEMPROARY PRINTING LOVE NUMBERS FOR NOW
        print(f"h number: {Planet.Gravity.h}, l number: {Planet.Gravity.l}, k number: {Planet.Gravity.k}")
    return Planet, Params




def SetupGravity(Planet, Params):
    """ Reconfigure layer boundaries and gravity model into a format
        usable by gravity response calculation functions of PyALMA3.

        NOTE: Function models the read_model_pp of PyALMA3 so that we can ensure compatibility with the package's functions.

        Requires Planet attributes:
    """
    # Get liquid indices
    indsLiq = np.where(np.flip(Planet.phase) == 0)[0]
    if Params.CALC_NEW_GRAVITY and Planet.Do.VALID:
        """Combine data into model format that is required by PyALMA3"""
        # Note we have to use r_m[:-1] since r_m has one extra value than other arrays
        Planet.Gravity.model = np.vstack(
            [Planet.Reduced.rLayers_m, Planet.Reduced.rPhase, Planet.Reduced.rRho_kgm3, Planet.Reduced.rVP_kms, Planet.Reduced.rVS_kms, Planet.Reduced.rGS_GPa, Planet.Reduced.reta_Pas]).T

        # Convert parameter units to Pa and meters
        for index, (header, unit) in enumerate(zip(Planet.Gravity.columns, Planet.Gravity.units_PyALMA3)):
            if header in Planet.Gravity.parameters_to_convert:
                Planet.Gravity.model[:, index] = Planet.Gravity.model[:, index] * Planet.Gravity.parameters_to_convert[header]
        # Get indices of used properties
        rIndex  = Planet.Gravity.columns.index('r')
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
    # Return Planet and Params
    return Planet, Params
