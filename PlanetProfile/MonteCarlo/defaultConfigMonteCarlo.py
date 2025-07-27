""" Configuration settings specific to Monte Carlo calculations """
import numpy as np
from PlanetProfile.Utilities.defineStructs import MonteCarloParamsStruct

configMonteCarloVersion = 1  # Integer number for config file version. Increment when new settings are added to the default config file.

def montecarloAssign():
    """
    Assign Monte Carlo parameters for PlanetProfile runs.
    """
    MonteCarloParams = MonteCarloParamsStruct()
    
    # General settings
    MonteCarloParams.nRuns = 1000  # Number of Monte Carlo runs to perform
    MonteCarloParams.seed = None  # Random seed for reproducibility (None for random)
    MonteCarloParams.useParallel = True  # Whether to use parallel processing
    
    # Parameters to search over for self-consistent models
    MonteCarloParams.paramsToSearchSelfConsistent = [
        'R_m', 'Tb_K', 'wOcean_ppt', 'xFeS', 'rhoSilInput_kgm3'
    ]
    
    # Parameters to search over for non-self-consistent models
    MonteCarloParams.paramsToSearchNonSelfConsistent = [
        'dzIceI_km', 'D_km', 'Core_R_km', 'rho_iceIhCond_kgm3', 'rho_iceIhConv_kgm3',
        'rho_ocean_kgm3', 'rho_sil_kgm3', 'rho_core_kgm3', 'GS_ice_GPa', 'GS_sil_GPa', 'GS_core_GPa'
    ]
    
    # Distribution types for each parameter (currently only 'uniform' is supported)
    MonteCarloParams.paramsDistributions = {
        # Self-consistent parameters
        'R_m': 'uniform',
        'Tb_K': 'uniform', 
        'wOcean_ppt': 'uniform',
        'xFeS': 'uniform',
        'rhoSilInput_kgm3': 'uniform',
        
        # Non-self-consistent parameters
        'dzIceI_km': 'uniform',
        'D_km': 'uniform',
        'Core_R_km': 'uniform',
        'rho_iceIhCond_kgm3': 'uniform',
        'rho_iceIhConv_kgm3': 'uniform',
        'rho_ocean_kgm3': 'uniform',
        'rho_sil_kgm3': 'uniform',
        'rho_core_kgm3': 'uniform',
        'GSConvIh_GPa': 'uniform',
        'GS_sil_GPa': 'uniform',
        'GS_core_GPa': 'uniform',
        'kThermWater_WmK': 'uniform',
        'kThermIceIh_WmK': 'uniform',
        'kThermCore_WmK': 'uniform',
        'etaSil_Pas': 'uniform',
        'etaMelt_Pas': 'uniform',
        'TSurf_K': 'uniform',
        'EactIceIh_kJmol': 'uniform',
        'AndradeExponent': 'uniform'
        
    }
    
    # Parameter ranges [min, max] for each parameter
    MonteCarloParams.paramsRanges = {
        # Self-consistent parameters
        'R_m': [1.5e6, 1.6e6],  # Body radius in m (Europa range)
        'Tb_K': [260, 272],  # Bottom temperature in K
        'wOcean_ppt': [0, 100],  # Ocean salinity in ppt
        'xFeS': [0, 1],  # Core FeS mole fraction
        'rhoSilInput_kgm3': [2500, 3500],  # Silicate density in kg/m^3
        
        # Non-self-consistent parameters
        'dzIceI_km': [10, 50],  # Ice shell thickness in km
        'D_km': [50, 150],  # Ocean thickness in km
        'Core_R_km': [0, 1000],  # Core radius in km
        'rho_iceIhCond_kgm3': [850, 950],  # Conductive ice density in kg/m^3
        'rho_iceIhConv_kgm3': [950, 1050],  # Convective ice density in kg/m^3
        'rho_ocean_kgm3': [950, 1050],  # Ocean density in kg/m^3
        'rho_sil_kgm3': [3000, 3500],  # Silicate density in kg/m^3
        'rho_core_kgm3': [5000, 6000],  # Core density in kg/m^3
        'GSConvIh_GPa': [50, 150],  # Ice shear modulus in GPa
        'GS_sil_GPa': [30, 70],  # Silicate shear modulus in GPa
        'GS_core_GPa': [80, 120],  # Core shear modulus in GPa
        'kThermWater_WmK': [0.5, 1.5],  # Water thermal conductivity in W/mK
        'kThermIceIh_WmK': [0.5, 1.5],  # Ice Ih thermal conductivity in W/mK
        'kThermCore_WmK': [0.5, 1.5],  # Core thermal conductivity in W/mK
        'etaSil_Pas': [1e18, 1e20],  # Silicate viscosity in Pa s
        'etaMelt_Pas': [1e18, 1e20],  # Melt viscosity in Pa s
        'TSurf_K': [260, 272],  # Surface temperature in K
        'EactIceIh_kJmol': [0, 100],  # Activation energy for ocean in kJ/mol
        'AndradeExponent': [0, 10]  # Andrade exponent
    }
    
    # Output settings
    MonteCarloParams.saveResults = True  # Whether to save results to file
    MonteCarloParams.showPlots = True  # Whether to display plots
    MonteCarloParams.plotDistributions = True  # Whether to plot parameter distributions
    MonteCarloParams.plotResults = True  # Whether to plot Monte Carlo results distributions
    MonteCarloParams.plotOceanComps = False  # Whether to plot results by ocean composition
    MonteCarloParams.plotCorrelations = True  # Whether to plot parameter correlations
    MonteCarloParams.plotScatter = False  # Whether to plot scatter plots of parameter pairs
    MonteCarloParams.scatterParams = [['k2_love', 'Amp'], ['h2_love', 'phase']]  # Parameter pairs to plot
    MonteCarloParams.excSelectionScatter = {  # Which magnetic excitations to include in scatter plots
        'orbital': True,  # Key excitation
        'synodic': True,  # Key excitation  
        'synodic 2nd': False,
        'true anomaly': False
    }
    
    return MonteCarloParams
