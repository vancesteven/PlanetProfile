import numpy as np
import logging as log

def MagneticInduction(Planet, Params):
    """ Calculate induced magnetic moments for the body and prints them to disk.

        Requires Planet attributes:
            Magnetic.inductType, Magnetic.peaks_Hz
        Sets Planet attributes:

    """

    # Set Magnetic struct layer arrays as we need for induction calculations
    Planet = InductionProfile(Planet, Params)

    # Calculate induced magnetic moments
    Planet = CalcInducedMoments(Planet, Params)

    return Planet


def CalcInducedMoments(Planet, Params):
    """ Calculate induced magnetic moments based on conductivity profile,
        possible asymmetric shape, and excitation moments for this body.

        Sets Planet attributes:
            Magnetic.Binm_nT, Magnetic.Benm_nT
    """

    return Planet


def ReloadInduction(Planet, Params):
    """ Reload induced magnetic moments that have been printed to disk.

        Sets Planet attributes:
            Magnetic.Binm_nT
    """
    Planet.Magnetic.Binm_nT = np.loadtxt(Params.inducedMomentsFile, skiprows=1, unpack=False)

    return Planet


def InductionProfile(Planet, Params):
    """ Reconfigure layer boundaries and conductivities into a format
        usable by magnetic induction calculation functions.
        Optionally also identify asymmetric shape information from gravity.

        Requires Planet attributes:
        Sets Planet attributes:
            Magnetic.inductionBounds_m, Magnetic.inductionSigmas_m, Magnetic.asymShape
    """

    if Params.INCLUDE_ASYM:
        pass

    return Planet
