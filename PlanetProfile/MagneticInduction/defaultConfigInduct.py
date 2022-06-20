""" Configuration settings specific to induction calculations and plots """
import numpy as np
from PlanetProfile.Utilities.defineStructs import InductOgramParamsStruct, \
    ExcitationSpectrumParamsStruct, ConductLayerParamsStruct, Constants

configInductVersion = 2  # Integer number for config file version. Increment when new settings are added to the default config file.
inductOtype = 'rho'  # Type of inductogram plot to make. Options are "Tb", "phi", "rho", "sigma", where the first 3 are vs. salinity, and sigma is vs. thickness. Sigma/D plot is not self-consistent.
testBody = 'Europa'  # Assign test profiles to use excitation moments for this body
dftC = 5  # Default number of contours to include in induct-o-grams

def GetInductParams(inductOtype, cLevels, dftC, cFmt):
    """ Lots of these settings depend on each other, so for convenience we
        include settings in this function so they can live at the top of
        this file.
    """
    SigParams = ConductLayerParamsStruct()
    ExcSpecParams = ExcitationSpectrumParamsStruct()
    InductParams = InductOgramParamsStruct(inductOtype, cLevels, dftC, cFmt)

    # Conducting layer parameter settings
    SigParams.REDUCED_INDUCT = True  # Whether to limit number of ocean layers for faster computation of layered induction
    SigParams.INCLUDE_ASYM = False  # Whether to include asymmetry in the induction conductivity profile based on J2 and C22 values
    SigParams.CONCENTRIC_ASYM = False  # Whether to map a single asymmetric shape to all layers, concentrically, scaling by their radii.
    SigParams.ALLOW_LOW_PMAX = False  # Whether to allow Magnetic.pMax to be set to an integer less than 2.
    SigParams.asymFstring = 'Shape_4piNormDepth'

    # Excitation spectrum settings
    ExcSpecParams.nOmegaPts = 100  # Resolution in log frequency space for magnetic excitation spectra
    ExcSpecParams.interpMethod = 'cubic'  # Interpolation method for complex response amplitudes in Fourier spectrum
    ExcSpecParams.Tmin_hr = 1  # Cutoff period to limit range of Fourier space shown

    # Inductogram calculation and plot settings
    InductParams.colorType = 'zb'  # What parameter to use for color of points in phase space plots. Options are "Tmean", "zb".
    InductParams.SPECIFIC_CLEVELS = True  # Whether to use the specific cLevels listed below (in GetClevels) or use default numbers
    InductParams.excSelectionCalc = {'synodic': True, 'orbital': True, 'true anomaly': True,  'synodic harmonic': True}  # Which magnetic excitations to include in calculations
    InductParams.excSelectionPlot = {'synodic': True, 'orbital': True, 'true anomaly': False, 'synodic harmonic': True}  # Which magnetic excitations to include in plotting
    InductParams.nwPts = 80  # Resolution for salinity values in ocean salinity vs. other plots
    InductParams.wMin = {'Europa': np.log10(1)}
    InductParams.wMax = {'Europa': np.log10(100)}
    InductParams.nTbPts = 60  # Resolution for Tb values in ocean salinity/Tb plots
    InductParams.TbMin = {'Europa': 262.0}
    InductParams.TbMax = {'Europa': 267.0}
    InductParams.nphiPts = 60  # Resolution for phiRockMax values in ocean salinity/phiMax plots
    InductParams.phiMin = {'Europa': np.log10(0.01)}
    InductParams.phiMax = {'Europa': np.log10(0.75)}
    InductParams.nrhoPts = 100  # Resolution for silicate density values in ocean salinity/rho plots
    InductParams.rhoMin = {'Europa': 3300}
    InductParams.rhoMax = {'Europa': 3700}
    InductParams.nSigmaPts = 80  # Resolution for conductivity values in ocean conductivity/thickness plots
    InductParams.sigmaMin = {'Europa': np.log10(1e-1)}
    InductParams.sigmaMax = {'Europa': np.log10(1e2)}
    InductParams.nDpts = 90  # Resolution for ocean thickness as for conductivity
    InductParams.Dmin = {'Europa': np.log10(1e0)}
    InductParams.Dmax = {'Europa': np.log10(2e2)}
    InductParams.zbFixed_km = {'Europa': 20}
    InductParams.EckhardtSolveMethod = 'RK45'  # Numerical solution method for scipy.integrate.solve_ivp. See https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
    InductParams.rMinODE = 1e3  # Minimum radius to use for numerical solution. Cannot be zero because of singularity at the origin.
    InductParams.oceanInterpMethod = 'linear'  # Interpolation method for determining ocean conductivities when REDUCED_INDUCT is True.
    InductParams.nIntL = 5  # Number of ocean layers to use when REDUCED_INDUCT = 1

    return SigParams, ExcSpecParams, InductParams


def GetClevels(inductOtype, dftC, testBody):
    # Specific contours for plots
    if inductOtype == 'sigma':
        cLevels = {
            'Europa': {
                'synodic':         {'Amp': dftC, 'Bx': [15, 50, 80, 170, 190, 200], 'By': dftC, 'Bz': dftC, 'phase': dftC+2},
                'orbital':         {'Amp': dftC, 'Bx': [1, 4, 11, 12.5, 13.3, 14.3], 'By': dftC, 'Bz': dftC, 'phase': dftC-2},
                'true anomaly':    {'Amp': dftC, 'Bx': dftC, 'By': dftC, 'Bz': dftC, 'phase': dftC-3},
                'synodic harmonic':{'Amp': dftC, 'Bx': [6, 9, 9.9], 'By': dftC, 'Bz': dftC, 'phase': dftC-2}
            }
        }
    else:
        cLevels = {
            'Europa': {
                'synodic': {'Amp': dftC, 'Bx': dftC, 'By': dftC, 'Bz': dftC, 'phase': dftC+2},
                'orbital': {'Amp': dftC, 'Bx': dftC, 'By': dftC, 'Bz': dftC, 'phase': dftC-2},
                'true anomaly': {'Amp': dftC, 'Bx': dftC, 'By': dftC, 'Bz': dftC, 'phase': dftC-3},
                'synodic harmonic': {'Amp': dftC, 'Bx': dftC, 'By': dftC, 'Bz': dftC, 'phase': dftC-2}
            }
        }

    cLevels['Test'] = cLevels[testBody]
    return cLevels


def GetContourFmt(testBody):
    deftFmt = '%1.1f'  # Default contour label format string
    deftPhi = '%1.0f'  # Default for phases in degrees (whole numbers)
    deftAmp = '%1.1f'  # Default for amplitudes (which are typically less than 1)
    cFmt = {
        'Europa': {
            'synodic':         {'Amp': deftAmp, 'Bx': '%1.0f', 'By': '%1.0f', 'Bz': deftFmt, 'phase': deftPhi},
            'orbital':         {'Amp': deftAmp, 'Bx': deftFmt, 'By': deftFmt, 'Bz': deftFmt, 'phase': deftPhi},
            'true anomaly':    {'Amp': deftAmp, 'Bx': '%1.2f', 'By': '%1.2f', 'Bz': deftFmt, 'phase': deftPhi},
            'synodic harmonic':{'Amp': deftAmp, 'Bx': deftFmt, 'By': deftFmt, 'Bz': '%1.2f', 'phase': deftPhi}
        }
    }

    cFmt['Test'] = cFmt[testBody]
    return cFmt


# Construct information for accessing by functions
cLevels = GetClevels(inductOtype, dftC, testBody)
cFmt = GetContourFmt(testBody)
SigParams, ExcSpecParams, InductParams = GetInductParams(inductOtype, cLevels, dftC, cFmt)
