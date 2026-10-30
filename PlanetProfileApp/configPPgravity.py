""" Configuration settings specific to gravity response calculations and plots """
from PlanetProfile.Utilities.defineStructs import GravityParamsStruct

configGravityVersion = 3 # Integer number for config file version. Increment when new settings are added to the default config file.

def gravityAssign():
    GravityParams = GravityParamsStruct()
    # Verbose settings of PyALMA - #TODO: Need to update PYALMA to use logger
    GravityParams.verbose = True

    # Parallel computing
    GravityParams.parallel = False  # Use Parallel computing for PyALMA calculations. #TODO: Need to implement way to do this if Parallel already being used in Exploreogram

    # Parsing parameters
    GravityParams.rheology_models = {'0': 'newton', 'Ih': 'maxwell', 'Ih_conv': 'andrade','II': 'maxwell', 'III': 'maxwell', 'III_conv': 'andrade',
                                     'IV': 'maxwell', 'V': 'maxwell','V_conv': 'andrade', 'VI': 'maxwell',
                                     'Sil': 'elastic', 'Fe': 'elastic', 'Clath': 'newton', 'Clath_conv': 'andrade'}  # Rheology structure model, where each model corresponds to a layer

    # General parameters
    GravityParams.num_digits = 128  # Set precision
    GravityParams.gorder = 8  # Order of Gaver method
    GravityParams.tau = 0  # TODO: FIGURE OUT WHAT TAU DOES
    GravityParams.loading_type = 'tidal'  # Loading type to calculate love numbers - 'tidal' or 'loading'

    # Harmonic degrees parameters
    GravityParams.harmonic_degrees = [2]  # List of harmonic degrees to calculate - not that for compatibility with PlanetProfile plotting, user should only specify one harmonic to calculate if they desire PlanetProfile's plotting functionality.
    GravityParams.time_log_kyrs = [1e-9]  # List of time range in log_kyrs
    GravityParams.time_history_function = 'step'  # Function to use for time - 'step' or 'ramp'
    GravityParams.ramp_function_length_kyrs = None  # Ramp length in kyrs

    # Output parameters
    GravityParams.output_type = 'real'  # Output type - 'complex' or 'real'

    return GravityParams
