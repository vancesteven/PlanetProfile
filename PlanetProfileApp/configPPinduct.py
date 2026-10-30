""" Configuration settings specific to induction calculations and plots """
import numpy as np
from PlanetProfile.Utilities.defineStructs import InductOgramParamsStruct, \
    ExcitationSpectrumParamsStruct, ConductLayerParamsStruct, Constants

configInductVersion = 6  # Integer number for config file version. Increment when new settings are added to the default config file.

def inductAssign():
    inductOtype = 'rho'  # Type of inductogram plot to make. Options are "Tb", "phi", "rho", "sigma", where the first 3 are vs. salinity, and sigma is vs. thickness. Sigma/D plot is not self-consistent.
    testBody = 'Europa'  # Assign test profiles to use excitation moments for this body
    dftC = 5  # Default number of contours to include in induct-o-grams

    # Construct information for accessing by functions
    cLevels = GetClevels(inductOtype, dftC, testBody)
    cfmt = GetContourFmt(testBody)
    SigParams, ExcSpecParams, InductParams = GetInductParams(inductOtype, cLevels, dftC, cfmt)

    return SigParams, ExcSpecParams, InductParams, testBody


def GetInductParams(inductOtype, cLevels, dftC, cfmt):
    """ Lots of these settings depend on each other, so for convenience we
        include settings in this function so they can live at the top of
        this file.
    """
    SigParams = ConductLayerParamsStruct()
    ExcSpecParams = ExcitationSpectrumParamsStruct()
    InductParams = InductOgramParamsStruct(inductOtype, cLevels, dftC, cfmt)

    # Conducting layer parameter settings
    SigParams.REDUCED_INDUCT = True  # Whether to limit number of ocean layers for faster computation of layered induction
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
    InductParams.excSelectionCalc = {  # Which magnetic excitations to include in calculations
        'orbital': True,  # Key excitation
        'orbital 2nd': True,
        'orbital 2nd-TA beat': True,
        'orbital 2nd-TA-year beat': True,
        'orbital 2nd-TA-half year beat': True,
        'orbital 4th-TA 2nd-half year beat': True,
        'orbital-half year beat': True,
        'orbital-synodic beat': True,
        'orbital-year beat': True,
        'orbital+year beat': True,
        'synodic': True,  # Key excitation
        'synodic 2nd': True,  # Key excitation
        'synodic 2nd+orbital beat': True,
        'synodic 2nd-TA beat': True,
        'synodic 3rd': True,
        'synodic 4th': True,
        'synodic 5th': True,
        'synodic-orbital beat': True,
        'synodic+orbital beat': True,
        'synodic-TA beat': True,
        'synodic+TA beat': True,
        'true anomaly': True,  # Key excitation
        'TA-year beat': True,
        'TA+year beat': True,
        'TA 2nd-orbital beat': True
    }
    InductParams.excSelectionPlot = {  # Which magnetic excitations to include in plotting
        'orbital': True,  # Key excitation
        'orbital 2nd': False,
        'orbital 2nd-TA beat': False,
        'orbital 2nd-TA-year beat': False,
        'orbital 2nd-TA-half year beat': False,
        'orbital 4th-TA 2nd-half year beat': False,
        'orbital-half year beat': False,
        'orbital-synodic beat': False,
        'orbital-year beat': False,
        'orbital+year beat': False,
        'synodic': True,  # Key excitation
        'synodic 2nd': True,  # Key excitation
        'synodic 2nd+orbital beat': False,
        'synodic 2nd-TA beat': False,
        'synodic 3rd': False,
        'synodic 4th': False,
        'synodic 5th': False,
        'synodic-orbital beat': False,
        'synodic+orbital beat': False,
        'synodic-TA beat': False,
        'synodic+TA beat': False,
        'true anomaly': True,  # Key excitation
        'TA-year beat': False,
        'TA+year beat': False,
        'TA 2nd-orbital beat': False
    }
    InductParams.nwPts = 40  # Resolution for salinity values in ocean salinity vs. other plots
    InductParams.wMin = {'Europa': np.log10(1), 'Enceladus': np.log10(0.1)}
    InductParams.wMax = {'Europa': np.log10(Constants.stdSeawater_ppt), 'Enceladus': np.log10(1.5*Constants.stdSeawater_ppt)}
    InductParams.nTbPts = 30  # Resolution for Tb values in ocean salinity/Tb plots
    InductParams.TbMin = {'Europa': 262.0, 'Enceladus': 269.8}
    InductParams.TbMax = {'Europa': 267.0, 'Enceladus': 273.1}
    InductParams.nphiPts = 30  # Resolution for phiRockMax values in ocean salinity/phiMax plots
    InductParams.phiMin = {'Europa': np.log10(0.01), 'Enceladus': np.log10(0.01)}
    InductParams.phiMax = {'Europa': np.log10(0.75), 'Enceladus': np.log10(0.75)}
    InductParams.nrhoPts = 50  # Resolution for silicate density values in ocean salinity/rho plots
    InductParams.rhoMin = {'Europa': 3300, 'Enceladus': 1900}
    InductParams.rhoMax = {'Europa': 3700, 'Enceladus': 3100}
    InductParams.nSigmaPts = 40  # Resolution for conductivity values in ocean conductivity/thickness plots
    InductParams.sigmaMin = {'Europa': np.log10(1e-1), 'Enceladus': np.log10(1e-2)}
    InductParams.sigmaMax = {'Europa': np.log10(1e2), 'Enceladus': np.log10(1e1)}
    InductParams.nDpts = 50  # Resolution for ocean thickness as for conductivity
    InductParams.Dmin = {'Europa': np.log10(1e0), 'Enceladus': np.log10(1e0)}
    InductParams.Dmax = {'Europa': np.log10(2e2), 'Enceladus': np.log10(1e2)}
    InductParams.zbFixed_km = {'Europa': 20, 'Enceladus': 23}
    InductParams.EckhardtSolveMethod = 'RK45'  # Numerical solution method for scipy.integrate.solve_ivp. See https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
    InductParams.rMinODE = 1e3  # Minimum radius to use for numerical solution. Cannot be zero because of singularity at the origin.
    InductParams.oceanInterpMethod = 'linear'  # Interpolation method for determining ocean conductivities when REDUCED_INDUCT is True.
    InductParams.nIntL = 5  # Number of ocean layers to use when REDUCED_INDUCT = 1
    InductParams.SUM_NEAR = False  # Whether to sum together closely-spaced periods. Accuracy of this approach decreases with time away from J2000.
    InductParams.USE_NAMED_EXC = True  # Whether to make use of named periods defined in PlanetProfile.MagneticInduction.Moments for excitation calcs
    InductParams.minBe_nT = 1.0  # Minimum value in nT to use for excitation moments when not using specific periods
    
    return SigParams, ExcSpecParams, InductParams


def GetClevels(inductOtype, dftC, testBody):
    # Specific contours for plots
    if inductOtype == 'sigma':
        cLevels = {
            'Europa': {
                'synodic':         {'Amp': dftC, 'Bx': [15, 50, 80, 170, 190, 200], 'By': dftC, 'Bz': dftC, 'phase': dftC+2},
                'orbital':         {'Amp': dftC, 'Bx': [1, 4, 11, 12.5, 13.3, 14.3], 'By': dftC, 'Bz': dftC, 'phase': dftC-2},
                'true anomaly':    {'Amp': dftC, 'Bx': dftC, 'By': dftC, 'Bz': dftC, 'phase': dftC-3},
                'synodic 2nd':{'Amp': dftC, 'Bx': [6, 9, 9.9], 'By': dftC, 'Bz': dftC, 'phase': dftC-2}
            }
        }
    else:
        cLevels = {
            'Europa': {
                'synodic': {'Amp': dftC, 'Bx': dftC, 'By': dftC, 'Bz': dftC, 'phase': dftC+2},
                'orbital': {'Amp': dftC, 'Bx': dftC, 'By': dftC, 'Bz': dftC, 'phase': dftC-2},
                'true anomaly': {'Amp': dftC, 'Bx': dftC, 'By': dftC, 'Bz': dftC, 'phase': dftC-3},
                'synodic 2nd': {'Amp': dftC, 'Bx': dftC, 'By': dftC, 'Bz': dftC, 'phase': dftC-2}
            }
        }

    cLevels['Test'] = cLevels[testBody]
    return cLevels


def GetContourFmt(testBody):
    deftFmt = '%1.1f'  # Default contour label format string
    deftPhi = '%1.0f'  # Default for phases in degrees (whole numbers)
    deftAmp = '%1.1f'  # Default for amplitudes (which are typically less than 1)
    cfmt = {
        'Europa': {
            'synodic':         {'Amp': deftAmp, 'Bx': '%1.0f', 'By': '%1.0f', 'Bz': deftFmt, 'phase': deftPhi},
            'orbital':         {'Amp': deftAmp, 'Bx': deftFmt, 'By': deftFmt, 'Bz': deftFmt, 'phase': deftPhi},
            'true anomaly':    {'Amp': deftAmp, 'Bx': '%1.2f', 'By': '%1.2f', 'Bz': deftFmt, 'phase': deftPhi},
            'synodic 2nd':{'Amp': deftAmp, 'Bx': deftFmt, 'By': deftFmt, 'Bz': '%1.2f', 'phase': deftPhi}
        }
    }

    cfmt['Test'] = cfmt[testBody]
    return cfmt
