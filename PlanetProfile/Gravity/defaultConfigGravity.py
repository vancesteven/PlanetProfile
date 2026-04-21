""" Configuration settings specific to gravity response calculations and plots """
from PlanetProfile.Utilities.defineStructs import GravityParamsStruct

configGravityVersion = 4 # Integer number for config file version. Increment when new settings are added to the default config file.

def gravityAssign():
    GravityParams = GravityParamsStruct()
    # Verbose settings of PyALMA - #TODO: Need to update PYALMA to use logger
    GravityParams.verbose = True

    # Parallel computing
    GravityParams.parallel = False  # Use Parallel computing for PyALMA calculations. #TODO: Need to implement way to do this if Parallel already being used in Exploreogram

    # Julia ALMA backend (experimental, not yet implemented)
    GravityParams.USE_JULIA_ALMA = False  # Use Julia ALMA.jl backend for Love number calculations instead of PyALMA3.
                                           # Requires: Julia >= 1.6 installed, juliacall Python package (pip install juliacall), ALMA.jl source
                                           # Benefits: Potentially faster for multi-degree calculations and large exploreograms
                                           # Drawbacks: ~30s JIT compilation startup cost, additional dependency management
                                           # Status: NOT YET IMPLEMENTED - deferred pending benchmarking justification
                                           # See JULIA_ALMA3_ASSESSMENT.md for implementation approach

    # Parsing parameters
    GravityParams.rheology_models = {'0': 'newton', 'Ih': 'elastic', 'Ih_conv': 'andrade','II': 'maxwell', 'III': 'maxwell', 'III_conv': 'andrade',
                                     'IV': 'maxwell', 'V': 'maxwell','V_conv': 'andrade', 'VI': 'maxwell',
                                     'Sil': 'elastic', 'Fe': 'elastic', 'Clath': 'newton', 'Clath_conv': 'andrade'}  # Rheology structure model, where each model corresponds to a layer

    # General parameters
    GravityParams.num_digits = 128  # Set precision
    GravityParams.gorder = 4  # Order of Gaver method
    GravityParams.tau = 2  # TODO: FIGURE OUT WHAT TAU DOES
    GravityParams.loading_type = 'tidal'  # Loading type to calculate love numbers - 'tidal' or 'loading'

    # Harmonic degrees parameters
    GravityParams.harmonic_degrees = [2]  # List of harmonic degrees to calculate. Supports multiple degrees, e.g. [2, 3, 4] for higher-order Love numbers. Love number arrays will be 2D: (len(harmonic_degrees), len(time_log_kyrs)).
    GravityParams.time_log_kyrs = [1e-9]  # List of time range in log_kyrs
    GravityParams.time_history_function = 'step'  # Function to use for time - 'step' or 'ramp'
    GravityParams.ramp_function_length_kyrs = None  # Ramp length in kyrs

    # Output parameters
    GravityParams.output_type = 'complex'  # Output type - 'complex' or 'real'

    return GravityParams
