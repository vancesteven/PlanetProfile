"""
ResultsStructs: New hierarchical data structures for PlanetProfile results

This module contains the new architecture for storing and managing results from
Monte Carlo, Exploration, and Inductogram analyses with a common base structure.

DESIGN PRINCIPLE: All data arrays are stored in 2D format for consistency:
- Monte Carlo: shape (1, N_samples) - single "row" of N samples  
- Exploration: shape (N_x, N_y) - actual 2D parameter grid
- Inductogram: shape (N_params, N_frequencies) - parameters vs frequencies

This consistent 2D approach enables:
1. Unified plotting functions that work across all analysis types
2. Simpler data extraction without special case handling  
3. Easy extensibility to new analysis types
4. Consistent flatten_for_matlab methods without complex reshaping
"""

import numpy as np
import pickle
from scipy.io import savemat
import logging

# Assign logger
log = logging.getLogger('PlanetProfile')



class BaseResultsStruct:
    """
    Common data structure for all analysis types containing shared results.
    
    This structure holds data that is common across Monte Carlo, Exploration,
    and Inductogram analyses, eliminating duplication.
    """
    
    def __init__(self):
        # Identity
        self.bodyname = None  # Name of body modeled
        
        # Validity
        self.VALID = None  # Whether this profile is physically possible
        self.invalidReason = None  # Explanation for why any invalid solution failed
        
        # Core results
        self.CMR2calc = None  # Calculated C/MR^2 value
        self.Mtot_kg = None  # Total body mass in kg
        self.D_km = None  # Ocean layer thickness in km
        self.zb_km = None  # Ice shell thickness in km
        self.Rcore_km = None  # Core radius in km
        
        # Densities
        self.rhoOceanMean_kgm3 = None  # Mean ocean density in kg/m^3
        self.rhoSilMean_kgm3 = None  # Mean silicate density in kg/m^3
        self.rhoCoreMean_kgm3 = None  # Mean core density in kg/m^3
        
        # Electrical properties
        self.sigmaMean_Sm = None  # Mean ocean conductivity in S/m
        self.Tmean_K = None  # Mean ocean temperature in K
        self.oceanComp = None  # Ocean composition

        
        # Love numbers
        self.kLoveComplex = None  # k2 Love number
        self.hLoveComplex = None  # h2 Love number
        self.lLoveComplex = None  # l2 Love number
        self.deltaLoveComplex = None  # delta Love number
        self.kLoveAmp = None  # k2 Love number
        self.hLoveAmp = None  # h2 Love number
        self.lLoveAmp = None  # l2 Love number
        self.deltaLoveAmp = None  # delta Love number
        self.kLovePhase = None  # k2 Love number phase
        self.hLovePhase = None  # h2 Love number phase
        self.lLovePhase = None  # l2 Love number phase
        self.deltaLovePhase = None  # delta Love number phase
        
        # Input parameters (commonly varied)
        self.wOcean_ppt = None  # Ocean salinity in g/kg
        self.Tb_K = None  # Bottom temperature of ice shell in K
        self.Pb_MPa = None  # Bottom pressure of ice shell in MPa
        self.rhoSil_kgm3 = None  # Input silicate density in kg/m^3
        self.rhoCore_kgm3 = None  # Input core density in kg/m^3
        
        # Additional common fields from WriteExploreOgram and WriteMonteCarloResults
        # Identity and basic info
        self.NO_H2O = None  # Whether this is a waterless body
        
        # Moment of inertia constraints
        self.CMR2str = None  # LaTeX-formatted string for moment of inertia
        self.Cmeasured = None  # Input moment of inertia to match
        self.Cupper = None  # Upper bound for valid MoI matches
        self.Clower = None  # Lower bound for valid MoI matches
        
        # Body structure
        self.R_m = None  # Body radius in m
        
        # Layer thicknesses (ice layers)
        self.zSeafloor_km = None  # Depth to bottom of ocean in km
        self.dzIceI_km = None  # Thickness of surface ice layer in km
        self.dzClath_km = None  # Thickness of clathrate layer in km
        self.dzIceIII_km = None  # Thickness of undersea ice III layer in km
        self.dzIceIIIund_km = None  # Thickness of underplate ice III layer in km
        self.dzIceV_km = None  # Thickness of undersea ice V layer in km
        self.dzIceVund_km = None  # Thickness of underplate ice V layer in km
        self.dzIceVI_km = None  # Thickness of undersea ice VI layer in km
        self.dzWetHPs_km = None  # Total thickness of all undersea high-pressure ices in km
        self.eLid_km = None  # Thickness of stagnant lid conductive ice layer in km
        self.Dconv_m = None  # Thickness of convective layer in m
        
        # Additional electrical properties
        self.sigmaTop_Sm = None  # Ocean top conductivity in S/m
        
        
        # Seafloor and geochemistry
        self.Pseafloor_MPa = None  # Pressure at seafloor in MPa
        self.phiSeafloor_frac = None  # Rock porosity at seafloor
        self.affinitySeafloor_kJ = None  # Available chemical energy at seafloor
        self.affinityMean_kJ = None  # Mean available chemical energy
        self.pHSeafloor = None  # pH at the seafloor
        self.pHTop = None  # pH at the top of the ocean
        self.affinityTop_kJ = None  # Available chemical energy at top of ocean
        
        # Porosity and rock properties
        self.silPhiCalc_frac = None  # Calculated rock porosity (best-match P=0 value)
        
        # Additional input parameters from exploration
        self.zb_approximate_km = None  # Approximate ice shell thickness input
        self.xFeS = None  # Core FeS mole fraction
        self.rhoSilInput_kgm3 = None  # Input silicate density (alternative name for rhoSil_kgm3)
        self.silPhi_frac = None  # Silicate porosity fraction
        self.icePhi_frac = None  # Ice porosity fraction
        self.silPclosure_MPa = None  # Silicate pore closure pressure
        self.icePclosure_MPa = None  # Ice pore closure pressure
        self.ionosTop_km = None  # Ionosphere top altitude
        self.sigmaIonos_Sm = None  # Ionosphere conductivity
        self.Htidal_Wm3 = None  # Tidal heating rate
        self.Qrad_Wkg = None  # Radiogenic heating rate
        
class InversionData:
    def __init__(self):
        self.gridWithinAllUncertainty = None
        self.gridWithinInductionResponseUncertainty = None
        self.gridWithinkLoveAmpUncertainty = None
        self.gridWithinkLovePhaseUncertainty = None
        self.gridWithinhLoveAmpUncertainty = None
        self.gridWithinhLovePhaseUncertainty = None

class InductionData:
    """
    Standardized induction results structure.
    
    This holds magnetic induction calculation results in a consistent format
    across all analysis types.
    """
    
    def __init__(self):
        self.Amp = None  # Amplitude of dipole response (modulus of complex response)
        self.phase = None  # Phase delay in degrees (positive)
        self.Bix_nT = None  # Induced magnetic field x-component in nT
        self.Biy_nT = None  # Induced magnetic field y-component in nT  
        self.Biz_nT = None  # Induced magnetic field z-component in nT
        self.Bi1x_nT = None  # Complex induced magnetic field x-component in nT
        self.Bi1y_nT = None  # Complex induced magnetic field y-component in nT
        self.Bi1z_nT = None  # Complex induced magnetic field z-component in nT
        self.rBi1x_nT = None  # Real part of complex induced magnetic field x-component in nT
        self.rBi1y_nT = None  # Real part of complex induced magnetic field y-component in nT
        self.rBi1z_nT = None  # Real part of complex induced magnetic field z-component in nT
        self.iBi1x_nT = None  # Imaginary part of complex induced magnetic field x-component in nT
        self.iBi1y_nT = None  # Imaginary part of complex induced magnetic field y-component in nT
        self.iBi1z_nT = None  # Imaginary part of complex induced magnetic field z-component in nT
        self.Bi1Tot_nT = None  # Total complex induced magnetic field in nT
        self.rBi1Tot_nT = None  # Real part of total complex induced magnetic field in nT
        self.iBi1Tot_nT = None  # Imaginary part of total complex induced magnetic field in nT
        self.Texc_hr = None  # Excitation period in hours
        self.freq_Hz = None  # Excitation frequency in Hz
        self.period_hr = None  # Excitation period in hours (same as Texc_hr for consistency)
        
        # Additional fields from Monte Carlo results
        self.calcedExc = None # List of excitations calculated
        self.nPeaks = None  # Number of peaks in magnetic response
        self.excSelectionCalc = None  # Excitation selection calculation results

class MonteCarloStatistics:
    """
    Monte Carlo specific statistical data.
    
    This holds information specific to Monte Carlo analyses including
    run statistics and parameter sampling information.
    """
    
    def __init__(self):
        # Run statistics
        self.nRuns = None  # Number of Monte Carlo runs performed
        self.nSuccess = None  # Number of successful runs
        self.successRate = None  # Success rate as a fraction
        self.seed = None  # Random seed used for reproducibility
        
        # Parameter information
        self.paramsToSearch = None  # List of parameter names that were varied
        self.paramsRanges = None  # Dictionary of parameter ranges used
        self.paramsDistributions = None  # Dictionary of distribution types used
        self.paramValues = None  # Dictionary of parameter name: array of sampled values
        
        # Timing information  
        self.totalTime_s = None  # Total execution time in seconds
        self.avgTime_s = None  # Average time per model in seconds


class MonteCarloResults:
    """
    Container for Monte Carlo analysis results.
    
    This combines base results with Monte Carlo-specific statistics and provides
    methods for data management and export.
    """
    
    def __init__(self):
        self.base = BaseResultsStruct()  # Common results data
        self.induction = InductionData()  # Induction results (if applicable)
        self.statistics = MonteCarloStatistics()  # Monte Carlo specific data
        
        # Additional MC-specific results arrays
        self.runtimePerModel_s = None  # Array of runtime per model in seconds
    


class ExplorationResults:
    """
    Container for Exploration analysis results.
    
    This combines base results with grid data and provides methods for 
    data management and export.
    """
    
    def __init__(self):
        self.base = BaseResultsStruct()  # Common results data
        self.induction = InductionData()  # Induction results (if applicable)
        self.inversion = InversionData()  # Inversion results (if applicable)
        
        # Grid data stored directly in ExplorationResults (no separate GridData class)
        self.xData = None  # 2D grid data for x-axis variable
        self.yData = None  # 2D grid data for y-axis variable
        self.xName = None  # Name of x-axis variable  
        self.yName = None  # Name of y-axis variable
        self.nx = None  # Size of x-axis variable
        self.ny = None  # Size of y-axis variable 
        self.xUnits = None  # Units for x-axis variable
        self.yUnits = None  # Units for y-axis variable
        self.xScale = 'linear'  # Scale type for x-axis ('linear' or 'log')
        self.yScale = 'linear'  # Scale type for y-axis ('linear' or 'log')
        
        # Additional exploration-specific results
        self.zName = None  # Name of z-axis variable
        self.CMR2str = None  # LaTeX-formatted string for moment of inertia
        self.Cmeasured = None  # Input moment of inertia to match
        self.Cupper = None  # Upper bound for valid MoI matches
        self.Clower = None  # Lower bound for valid MoI matches
        
        self.excName = None  # Name of excitation variable
    


class InductionResults:
    """
    Container for Inductogram analysis results.
    
    This combines base results with induction data and provides methods for
    data management and export.
    """
    
    def __init__(self):
        self.base = BaseResultsStruct()  # Common results data
        self.induction = InductionData()  # Induction results
        
        self.bodyname = None  # Name of body modeled.
        self.yName = None  # Name of variable along y axis. Options are "Tb", "phi", "rho", "sigma", where the first 3 are vs. salinity, and sigma is vs. thickness.
        self.Texc_hr = None  # Dict of excitation periods modeled.
        self.Amp = None  # Amplitude of dipole response (modulus of complex dipole response).
        self.phase = None  # (Positive) phase delay in degrees.
        self.Bix_nT = None  # Induced Bx dipole moments relative to body surface in nT for each excitation.
        self.Biy_nT = None  # Induced By dipole moments relative to body surface in nT for each excitation.
        self.Biz_nT = None  # Induced Bz dipole moments relative to body surface in nT for each excitation.
        self.wOcean_ppt = None  # Values of salinity used.
        self.oceanComp = None  # Ocean composition used.
        self.Tb_K = None  # Values of Bulk.Tb_K used.
        self.rhoSilMean_kgm3 = None  # Values of Sil.rhoMean_kgm3 resulted (also equal to those set for all but phi inductOtype).
        self.phiRockMax_frac = None  # Values of Sil.phiRockMax_frac set.
        self.Tmean_K = None  # Ocean mean temperature result in K.
        self.sigmaMean_Sm = None  # Mean ocean conductivity. Used to map plots vs. salinity onto D/sigma plots.
        self.sigmaTop_Sm = None  # Ocean top conductivity. Used to map plots vs. salinity onto D/sigma plots.
        self.D_km = None  # Ocean layer thickness in km. Used to map plots vs. salinity onto D/sigma plots.
        self.zb_km = None  # Upper ice shell thickness in km.
        self.R_m = None  # Body radius in m, used to scale amplitudes.
        self.rBds_m = None  # Conducting layer upper boundaries in m.
        self.sigmaLayers_Sm = None  # Conductivities below each boundary in S/m.
        self.zb_approximate_km = None  # Upper approximate ice shell thickness in km.
        self.oceanComp = None  # Ocean compositions
        

        self.x = None  # Variable to plot on x axis of inductogram plots
        self.y = None  # Variable to plot on y axis of inductogram plots
        self.nx = None  # Size of x-axis variable
        self.ny = None  # Size of y-axis variable
        self.compsList = None  # Linear list of compositions for each model point
        self.comps = None  # Minimal list of compositions, with 1 entry per comp
        self.SINGLE_COMP = None  # Boolean flag for tracking if all of the models have the same composition
    
    def SetAxes(self, inductOtype):
        # Set the x and y variables to plot in inductograms based on inductOtype
        if inductOtype == 'sigma':
            self.x = self.sigmaMean_Sm
            self.y = self.D_km
        elif inductOtype == 'oceanComp':
            self.x = self.oceanComp
            self.y = self.zb_approximate_km
        else:
            self.x = self.wOcean_ppt
            if inductOtype == 'Tb':
                self.y = self.Tb_K
            elif inductOtype == 'rho':
                self.y = self.rhoSilMean_kgm3
            elif inductOtype == 'phi':
                self.y = self.phiRockMax_frac
            else:
                raise ValueError(f'inductOtype {inductOtype} not recognized.')

    def SetComps(self, inductOtype):
        # Set some attributes pertaining to handling multiple ocean compositions in plots
        self.compsList = self.oceanComp.flatten()
        # Change any array with CustomSolution to just CustomSolution
        for i in range(self.compsList.size):
            if 'CustomSolution' in self.compsList[i]:
                self.compsList[i] = 'CustomSolution'
        if np.all(self.compsList == self.compsList[0]) and inductOtype != 'sigma':
            self.SINGLE_COMP = True
            self.comps = [self.compsList[0]]
        else:
            self.SINGLE_COMP = False
            self.comps = np.unique(self.compsList)