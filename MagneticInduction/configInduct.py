""" Configuration settings specific to induction calculations and plots """
import numpy as np
from Utilities.defineStructs import Constants

class InductOgramParams:
    def __init__(self):
        # Type of InductOgram plot to make. Options are "Tb", "phi", "rho", "sigma",
        # where the first 3 are vs. salinity, and sigma is vs. thickness.
        # Sigma/D plot is not self-consistent.
        self.inductOtype = 'rho'

        self.excSelectionCalc = {'synodic': True, 'orbital': True, 'true anomaly': True,  'synodic harmonic': True}  # Which magnetic excitations to include in calculations
        self.excSelectionPlot = {'synodic': True, 'orbital': True, 'true anomaly': False, 'synodic harmonic': True}  # Which magnetic excitations to include in plotting
        # Force calculations to be done for each oscillation to be plotted
        for osc in self.excSelectionPlot:
            if self.excSelectionPlot[osc] and not self.excSelectionCalc[osc]:
                self.excSelectionCalc[osc] = True
        self.nwPts = 5  # Resolution for salinity values in ocean salinity vs. other plots
        self.wMin = {'Europa': np.log10(0.05 * Constants.stdSeawater_ppt)}
        self.wMax = {'Europa': np.log10(Constants.stdSeawater_ppt)}
        self.nTbPts = 60  # Resolution for Tb values in ocean salinity/Tb plots
        self.nphiPts = 60  # Resolution for phiRockMax values in ocean salinity/phiMax plots
        self.phiMin = {'Europa': np.log10(0.01)}
        self.phiMax = {'Europa': np.log10(0.75)}
        self.nrhoPts = 6  # Resolution for silicate density values in ocean salinity/rho plots
        self.rhoMin = {'Europa': 3500}
        self.rhoMax = {'Europa': 3700}
        self.nSigmaPts = 50  # Resolution for conductivity values in ocean conductivity/thickness plots
        self.sigmaMin = {'Europa': np.log10(1e-1)}
        self.sigmaMax = {'Europa': np.log10(1e2)}
        self.nDpts = 60  # Resolution for ocean thickness as for conductivity
        self.Dmin = {'Europa': np.log10(1e0)}
        self.Dmax = {'Europa': np.log10(2e2)}
        self.zbFixed_km = {'Europa': 20}
        self.EckhardtSolveMethod = 'RK45'  # Numerical solution method for scipy.integrate.solve_ivp. See https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
        self.rMinODE = 1e3  # Minimum radius to use for numerical solution. Cannot be zero because of singularity at the origin.
        self.oceanInterpMethod = 'linear'  # Interpolation method for determining ocean conductivities when REDUCED_INDUCT is True.
        self.nIntL = 5  # Number of ocean layers to use when REDUCED_INDUCT = 1
        #self.opts_odeParams = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep', 2e3,'InitialStep',1e-2)
        #To be implemented- eventually need some ODE numerical solution parameters
        #self.opts_odeLayers = odeset('RelTol',1e-8, 'AbsTol',1e-10,'MaxStep',10e3,'InitialStep',1e-2)


class AsymmetryStruct:
    def __init__(self):
        self.pMax = {
            'Europa': 2
        }
        # The below is just a placeholder and eventually will need a robust way to read these in and manipulate.
        self.shape = {
            'Europa': np.zeros((2, self.pMax['Europa']+1, self.pMax['Europa']+1), dtype=np.complex_)
        }

    def chipq(self, bodyname, nBds, pMax):
        return np.tile(self.shape[bodyname][:, :pMax+1, :pMax+1], (nBds, 1))


class ExcitationSpectrumParams:
    def __init__(self):
        self.nOmegaPts = 100  # Resolution in log frequency space for magnetic excitation spectra
        self.nOmegaFine = 1000  # Fine-spacing resolution for log frequency spectrum

InductParams = InductOgramParams()
Asymmetry = AsymmetryStruct()
ExcSpecParams = ExcitationSpectrumParams()

dftC = 5  # Default number of contours to include in induct-o-grams
cLevels = {
    'Europa':{
        'synodic':         {'Amp': dftC, 'Bx': dftC, 'By': dftC, 'Bz': dftC, 'phase': dftC+2},
        'orbital':         {'Amp': dftC, 'Bx': dftC, 'By': dftC, 'Bz': dftC, 'phase': dftC-2},
        'true anomaly':    {'Amp': dftC, 'Bx': dftC, 'By': dftC, 'Bz': dftC, 'phase': dftC-3},
        'synodic harmonic':{'Amp': dftC, 'Bx': dftC, 'By': dftC, 'Bz': dftC, 'phase': dftC-2}
    }
}

deftFmt = '%1.1f'  # Default contour label format string
deftPhi = '%1.0f'  # Default for phases in degrees (whole numbers)
deftAmp = '%1.1f'  # Default for amplitudes (which are typically less than 1)
cFmt = {
    'Europa':{
        'synodic':         {'Amp': deftAmp, 'Bx': '%1.0f', 'By': '%1.0f', 'Bz': deftFmt, 'phase': deftPhi},
        'orbital':         {'Amp': deftAmp, 'Bx': deftFmt, 'By': deftFmt, 'Bz': deftFmt, 'phase': deftPhi},
        'true anomaly':    {'Amp': deftAmp, 'Bx': deftFmt, 'By': deftFmt, 'Bz': deftFmt, 'phase': deftPhi},
        'synodic harmonic':{'Amp': deftAmp, 'Bx': deftFmt, 'By': deftFmt, 'Bz': deftFmt, 'phase': deftPhi}
    }
}
