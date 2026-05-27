""" Configuration settings specific to gravity response calculations and plots """
from PlanetProfile.Utilities.defineStructs import GravityParamsStruct

configGravityVersion = 4 # Integer number for config file version. Increment when new settings are added to the default config file.

def gravityAssign():
    GravityParams = GravityParamsStruct()
    # Backend: 'pyalma' (default, uses PyALMA3) or 'tidalpy' (per-layer tidal heating dissipation)
    GravityParams.backend = 'pyalma'

    # Verbose settings of PyALMA - #TODO: Need to update PYALMA to use logger
    GravityParams.verbose = True

    # Parallel computing
    GravityParams.parallel = False  # Use Parallel computing for PyALMA calculations. #TODO: Need to implement way to do this if Parallel already being used in Exploreogram

    # Rheology model for each layer
    # Options: 'andrade', 'maxwell', 'elastic', 'newton'
    # TidalPy backend supports all four; PyALMA3 backend supports andrade, maxwell, elastic, newton
    GravityParams.rheology_models = {'0': 'newton', 'Ih': 'andrade', 'Ih_conv': 'andrade',
                                     'II': 'andrade', 'III': 'andrade', 'III_conv': 'andrade',
                                     'IV': 'andrade', 'V': 'andrade', 'V_conv': 'andrade', 'VI': 'andrade',
                                     'Sil': 'andrade', 'Fe': 'elastic', 'Clath': 'elastic', 'Clath_conv': 'andrade'}  # Rheology structure model, where each model corresponds to a layer

    # General parameters
    GravityParams.num_digits = 128  # Set precision
    GravityParams.gorder = 4  # Order of Gaver method
    GravityParams.tau = 2  # TODO: FIGURE OUT WHAT TAU DOES
    GravityParams.loading_type = 'tidal'  # Loading type to calculate love numbers - 'tidal' or 'loading'

    # Harmonic degrees parameters
    GravityParams.harmonic_degrees = [2]  # List of harmonic degrees to calculate - not that for compatibility with PlanetProfile plotting, user should only specify one harmonic to calculate if they desire PlanetProfile's plotting functionality.
    GravityParams.time_log_kyrs = [1e-9]  # List of time range in log_kyrs
    GravityParams.time_history_function = 'step'  # Function to use for time - 'step' or 'ramp'
    GravityParams.ramp_function_length_kyrs = None  # Ramp length in kyrs

    # Output parameters
    GravityParams.output_type = 'complex'  # Output type - 'complex' or 'real'

    return GravityParams
