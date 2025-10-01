""" Configuration settings specific to inversion calculations """
import numpy as np
from PlanetProfile.Utilities.defineStructs import InversionParamsStruct

configInversionVersion = 1  # Integer number for config file version. Increment when new settings are added to the default config file.

def inversionAssign():
    """
    Assign inversion parameters for PlanetProfile runs.
    """
    InversionParams = InversionParamsStruct()
    InversionParams.spacecraftUncertainties = \
    {'Clipper': {'InductionResponseUncertainty_nT': 1.5, 'kLoveAmpUncertainity': 0.018, 'hLoveAmpUncertainity': 0.1}}
    InversionParams.setSpaceCraft('Clipper')
    return InversionParams
    
    