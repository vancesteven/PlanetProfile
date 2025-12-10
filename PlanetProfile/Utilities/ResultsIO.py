"""
ResultsIO: Input/Output functions for PlanetProfile hierarchical results structures
This module contains functions for writing and reloading results using the new
hierarchical data structures, supporting both pickle and MATLAB formats.
"""

import numpy as np
import os
import logging
import pickle
from scipy.io import savemat
from PlanetProfile.MagneticInduction.Moments import Excitations
from PlanetProfile.MagneticInduction.MagneticInduction import Benm2absBexyz
from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc

# Assign logger
log = logging.getLogger('PlanetProfile')

def WriteResults(Results, pklFilePath, saveMatlab = False, matlabFilePath=None):
    """
    Save results to disk in pickle format (always) and optionally MATLAB format.
    
    Args:
        results: MonteCarloResults, ExplorationResults, or InductionResults object
        filepath (str): Base filepath (without extension)
        save_matlab (bool): Whether to also save as .mat file
    """
    with open(pklFilePath, 'wb') as f:
        pickle.dump(Results, f)
    log.info(f"Results saved to pickle file: {pklFilePath}")
    
    # Optionally save as MATLAB format for compatibility
    if saveMatlab:
        flat_dict = flatten_dict_for_matlab(Results)
        savemat(matlabFilePath, flat_dict)
        log.info(f"Results saved to MATLAB file: {matlabFilePath}")


def ReloadResultsFromPickle(fName):
    """
    Reload results from pickle file.
    
    Args:
        fName (str): Path to pickle file
        
    Returns:
        MonteCarloResults, ExplorationResults, or InductionResults object
    """
    if not os.path.isfile(fName):
        raise FileNotFoundError(f"File {fName} not found.")
    with open(fName, 'rb') as f:
        results = pickle.load(f)
    
    log.info(f"Results loaded from pickle file: {fName}")
    return results


def ExtractResults(Results, PlanetGrid, Params):
    """
    Extract results from PlanetGrid.
    """
    Results.nx = PlanetGrid.shape[0]
    Results.ny = PlanetGrid.shape[1]
    Results.base = ExtractBasePlanetData(Results.base, PlanetGrid)
    Results.induction = ExtractInductionData(Results.induction, Results.base.bodyname, PlanetGrid, Params)
    # Set exploration-specific fields that don't fit the base pattern
    Results.CMR2str = f'$C/MR^2 = {PlanetGrid[0,0].CMR2str}$'
    Results.Cmeasured = PlanetGrid[0,0].Bulk.Cmeasured
    Results.Cupper = PlanetGrid[0,0].Bulk.CuncertaintyUpper  
    Results.Clower = PlanetGrid[0,0].Bulk.CuncertaintyLower
    
    # Set grid information
    Results.xName = Params.Explore.xName
    Results.yName = Params.Explore.yName
    Results.zName  = Params.Explore.zName
    Results.xData = getattr(Results.base, Params.Explore.xName)
    Results.yData = getattr(Results.base, Params.Explore.yName)
    Results.nx = Params.Explore.nx
    Results.ny = Params.Explore.ny
    Results.titleAddendum = FigLbl.titleAddendum
    return Results


def ExtractBasePlanetData(baseStruct, PlanetGrid):
    """
    Extract common base data from Planet objects.
    Always returns 2D arrays for consistency across analysis types.
    
    Args:
        PlanetArray: Array of Planet objects (1D for Monte Carlo, 2D for Exploration)
        
    Returns:
        dict: Dictionary of extracted base data arrays (always 2D)
    """
    # Force to 2D array if it's 1D (Monte Carlo case: reshape to 1 x N)
    if PlanetGrid.ndim == 1:
        PlanetGrid = PlanetGrid.reshape(1, -1)
    
    
    base_data = {
        # Identity and validity
        'bodyname': PlanetGrid[0,0].name,
        'VALID': np.array([[Planeti.Do.VALID if hasattr(Planeti.Do, 'VALID') else False for Planeti in line] for line in PlanetGrid]),
        'invalidReason': np.array([[getattr(Planeti, 'invalidReason', '') for Planeti in line] for line in PlanetGrid]),
        'NO_H2O': PlanetGrid[0,0].Do.NO_H2O,
        
        # Core results
        'CMR2mean': np.array([[getattr(Planeti, 'CMR2mean', np.nan) for Planeti in line] for line in PlanetGrid]),
        'Mtot_kg': np.array([[getattr(Planeti, 'Mtot_kg', np.nan) for Planeti in line] for line in PlanetGrid]),
        'D_km': np.array([[getattr(Planeti, 'D_km', np.nan) for Planeti in line] for line in PlanetGrid]),
        'zb_km': np.array([[getattr(Planeti, 'zb_km', np.nan) for Planeti in line] for line in PlanetGrid]),
        'Rcore_km': np.array([[getattr(Planeti.Core, 'Rmean_m', np.nan) / 1e3 if getattr(Planeti.Core, 'Rmean_m', np.nan) is not None else np.nan for Planeti in line] for line in PlanetGrid]),
        
        # Body structure  
        'R_m': np.array([[Planeti.Bulk.R_m for Planeti in line] for line in PlanetGrid]),
        'zSeafloor_km': np.array([[getattr(Planeti, 'D_km', 0) + getattr(Planeti, 'zb_km', 0) for Planeti in line] for line in PlanetGrid]),
        'dzIceI_km': np.array([[getattr(Planeti, 'dzIceI_km', np.nan) for Planeti in line] for line in PlanetGrid]),
        'dzClath_km': np.array([[getattr(Planeti, 'dzClath_km', np.nan) for Planeti in line] for line in PlanetGrid]),
        'dzIceIII_km': np.array([[getattr(Planeti, 'dzIceIII_km', np.nan) for Planeti in line] for line in PlanetGrid]),
        'dzIceIIIund_km': np.array([[getattr(Planeti, 'dzIceIIIund_km', np.nan) for Planeti in line] for line in PlanetGrid]),
        'dzIceV_km': np.array([[getattr(Planeti, 'dzIceV_km', np.nan) for Planeti in line] for line in PlanetGrid]),
        'dzIceVund_km': np.array([[getattr(Planeti, 'dzIceVund_km', np.nan) for Planeti in line] for line in PlanetGrid]),
        'dzIceVI_km': np.array([[getattr(Planeti, 'dzIceVI_km', np.nan) for Planeti in line] for line in PlanetGrid]),
        'Dconv_m': np.array([[getattr(Planeti, 'Dconv_m', np.nan) for Planeti in line] for line in PlanetGrid]),
        'dzWetHPs_km': np.array([[getattr(Planeti, 'dzWetHPs_km', np.nan) for Planeti in line] for line in PlanetGrid]),
        'eLid_km': np.array([[getattr(Planeti, 'eLid_m', np.nan) / 1e3 if getattr(Planeti, 'eLid_m', np.nan) is not None else np.nan for Planeti in line] for line in PlanetGrid]),
        
        # Densities
        'rhoOceanMean_kgm3': np.array([[getattr(Planeti.Ocean, 'rhoMean_kgm3', np.nan) for Planeti in line] for line in PlanetGrid]),
        'rhoSilMean_kgm3': np.array([[getattr(Planeti.Sil, 'rhoMean_kgm3', np.nan) for Planeti in line] for line in PlanetGrid]),
        'rhoCoreMean_kgm3': np.array([[getattr(Planeti.Core, 'rhoMean_kgm3', np.nan) for Planeti in line] for line in PlanetGrid]),
        
        # Electrical properties
        'sigmaMean_Sm': np.array([[getattr(Planeti.Ocean, 'sigmaMean_Sm', np.nan) for Planeti in line] for line in PlanetGrid]),
        'sigmaTop_Sm': np.array([[getattr(Planeti.Ocean, 'sigmaTop_Sm', np.nan) for Planeti in line] for line in PlanetGrid]),
        'Tmean_K': np.array([[getattr(Planeti.Ocean, 'Tmean_K', np.nan) for Planeti in line] for line in PlanetGrid]),
        'oceanComp': np.array([[getattr(Planeti.Ocean, 'comp', '') for Planeti in line] for line in PlanetGrid]),
        
        # Love numbers
        'kLoveComplex': np.array([[getattr(Planeti.Gravity, 'k', np.nan) if hasattr(Planeti, 'Gravity') else np.nan for Planeti in line] for line in PlanetGrid]),
        'hLoveComplex': np.array([[getattr(Planeti.Gravity, 'h', np.nan) if hasattr(Planeti, 'Gravity') else np.nan for Planeti in line] for line in PlanetGrid]),
        'lLoveComplex': np.array([[getattr(Planeti.Gravity, 'l', np.nan) if hasattr(Planeti, 'Gravity') else np.nan for Planeti in line] for line in PlanetGrid]),
        'deltaLoveComplex': np.array([[getattr(Planeti.Gravity, 'delta', np.nan) if hasattr(Planeti, 'Gravity') else np.nan for Planeti in line] for line in PlanetGrid]),
        'kLoveAmp': np.array([[getattr(Planeti.Gravity, 'kAmp', np.nan) if hasattr(Planeti, 'Gravity') else np.nan for Planeti in line] for line in PlanetGrid]),
        'kLovePhase': np.array([[getattr(Planeti.Gravity, 'kPhase', np.nan) if hasattr(Planeti, 'Gravity') else np.nan for Planeti in line] for line in PlanetGrid]),
        'hLoveAmp': np.array([[getattr(Planeti.Gravity, 'hAmp', np.nan) if hasattr(Planeti, 'Gravity') else np.nan for Planeti in line] for line in PlanetGrid]),
        'hLovePhase': np.array([[getattr(Planeti.Gravity, 'hPhase', np.nan) if hasattr(Planeti, 'Gravity') else np.nan for Planeti in line] for line in PlanetGrid]),
        'lLoveAmp': np.array([[getattr(Planeti.Gravity, 'lAmp', np.nan) if hasattr(Planeti, 'Gravity') else np.nan for Planeti in line] for line in PlanetGrid]),
        'lLovePhase': np.array([[getattr(Planeti.Gravity, 'lPhase', np.nan) if hasattr(Planeti, 'Gravity') else np.nan for Planeti in line] for line in PlanetGrid]),
        'deltaLoveAmp': np.array([[getattr(Planeti.Gravity, 'deltaAmp', np.nan) if hasattr(Planeti, 'Gravity') else np.nan for Planeti in line] for line in PlanetGrid]),
        'deltaLovePhase': np.array([[getattr(Planeti.Gravity, 'deltaPhase', np.nan) if hasattr(Planeti, 'Gravity') else np.nan for Planeti in line] for line in PlanetGrid]),
        

        # Seafloor and geochemistry
        'Pseafloor_MPa': np.array([[getattr(Planeti, 'Pseafloor_MPa', np.nan) for Planeti in line] for line in PlanetGrid]),
        'phiSeafloor_frac': np.array([[getattr(Planeti, 'phiSeafloor_frac', np.nan) for Planeti in line] for line in PlanetGrid]),
        'affinitySeafloor_kJ': np.array([[getattr(Planeti.Ocean, 'affinitySeafloor_kJ', np.nan) for Planeti in line] for line in PlanetGrid]),
        'affinityMean_kJ': np.array([[getattr(Planeti.Ocean, 'affinityMean_kJ', np.nan) for Planeti in line] for line in PlanetGrid]),
        'pHSeafloor': np.array([[getattr(Planeti.Ocean, 'pHSeafloor', np.nan) for Planeti in line] for line in PlanetGrid]),
        'pHTop': np.array([[getattr(Planeti.Ocean, 'pHTop', np.nan) for Planeti in line] for line in PlanetGrid]),
        'affinityTop_kJ': np.array([[getattr(Planeti.Ocean, 'affinityTop_kJ', np.nan) for Planeti in line] for line in PlanetGrid]),
        'speciesRatioToChange': np.array([[getattr(Planeti.Ocean.Reaction, 'speciesRatioToChange', np.nan) for Planeti in line] for line in PlanetGrid]),
        'mixingRatioToH2O': np.array([[getattr(Planeti.Ocean.Reaction, 'speciesToChangeMixingRatio', np.nan) for Planeti in line] for line in PlanetGrid]),
        
        # Porosity and rock properties  
        'silPhiCalc_frac': np.array([[getattr(Planeti.Sil, 'phiCalc_frac', np.nan) for Planeti in line] for line in PlanetGrid]),
        
        # Input parameters (commonly varied)
        'wOcean_ppt': np.array([[Planeti.Ocean.wOcean_ppt for Planeti in line] for line in PlanetGrid]),
        'Tb_K': np.array([[Planeti.Bulk.Tb_K for Planeti in line] for line in PlanetGrid]),
        'zb_approximate_km': np.array([[getattr(Planeti.Bulk, 'zb_approximate_km', np.nan) for Planeti in line] for line in PlanetGrid]),
        'Pb_MPa': np.array([[getattr(Planeti.Bulk, 'Pb_MPa', np.nan) for Planeti in line] for line in PlanetGrid]),
        'rhoSil_kgm3': np.array([[getattr(Planeti.Sil, 'rhoSilWithCore_kgm3', np.nan) for Planeti in line] for line in PlanetGrid]),
        'rhoCore_kgm3': np.array([[getattr(Planeti.Core, 'rhoFe_kgm3', np.nan) for Planeti in line] for line in PlanetGrid]),
        'xFeS': np.array([[getattr(Planeti.Core, 'xFeS', np.nan) for Planeti in line] for line in PlanetGrid]),
        'silPhi_frac': np.array([[getattr(Planeti.Sil, 'phiRockMax_frac', np.nan) for Planeti in line] for line in PlanetGrid]),
        'icePhi_frac': np.array([[Planeti.Ocean.phiMax_frac.get('Ih', np.nan) if hasattr(Planeti.Ocean, 'phiMax_frac') and Planeti.Ocean.phiMax_frac else np.nan for Planeti in line] for line in PlanetGrid]),
        'silPclosure_MPa': np.array([[getattr(Planeti.Sil, 'Pclosure_MPa', np.nan) for Planeti in line] for line in PlanetGrid]),
        'icePclosure_MPa': np.array([[Planeti.Ocean.Pclosure_MPa.get('Ih', np.nan) if hasattr(Planeti.Ocean, 'Pclosure_MPa') and Planeti.Ocean.Pclosure_MPa else np.nan for Planeti in line] for line in PlanetGrid]),
        'ionosTop_km': np.array([[Planeti.Magnetic.ionosBounds_m[-1] / 1e3 if hasattr(Planeti, 'Magnetic') and hasattr(Planeti.Magnetic, 'ionosBounds_m') and Planeti.Magnetic.ionosBounds_m is not None else np.nan for Planeti in line] for line in PlanetGrid]),
        'sigmaIonos_Sm': np.array([[Planeti.Magnetic.sigmaIonosPedersen_Sm[-1] if hasattr(Planeti, 'Magnetic') and hasattr(Planeti.Magnetic, 'sigmaIonosPedersen_Sm') and Planeti.Magnetic.sigmaIonosPedersen_Sm is not None else np.nan for Planeti in line] for line in PlanetGrid]),
        'Htidal_Wm3': np.array([[getattr(Planeti.Sil, 'Htidal_Wm3', np.nan) for Planeti in line] for line in PlanetGrid]),
        'Qrad_Wkg': np.array([[getattr(Planeti.Sil, 'Qrad_Wkg', np.nan) for Planeti in line] for line in PlanetGrid]),
        
        # Additional common fields
        'qSurf_Wm2': np.array([[getattr(Planeti, 'qSurf_Wm2', np.nan) for Planeti in line] for line in PlanetGrid]),
    }

    # Ensure everything is set so things will play nicely with .mat saving and plotting functions
    nans = np.nan * base_data['R_m']
    for name, attr in base_data.items():
        if isinstance(attr, np.ndarray) and attr.ndim == 2 and np.all(attr == None):
            base_data[name] = nans
        setattr(baseStruct, name, base_data[name])
    for name, attr in baseStruct.__dict__.items():
        if attr is None:
            baseStruct.__dict__[name] = nans
    # Data is now always 2D - no need to flatten for Monte Carlo
    return baseStruct


def ExtractInductionData(InductionResults, bodyname, PlanetGrid, Params):
    """
    Extract induction-specific data from Planet objects.
    Always returns 2D arrays for consistency across analysis types.
    
    Args:
        PlanetArray: Array of Planet objects (1D for Monte Carlo, 2D for Exploration)
        Bex_nT, Bey_nT, Bez_nT: External field components
        
    Returns:
        dict: Dictionary of extracted induction data (always 2D)
    #TODO Make this data compatible with multiple n's
    """
    # Force to 2D array if it's 1D (Monte Carlo case: reshape to 1 x N)
    if PlanetGrid.ndim == 1:
        PlanetGrid = PlanetGrid.reshape(1, -1)
    
    BeList = Excitations(bodyname)
    eachT = np.logical_and([Params.Induct.excSelectionCalc[key] for key in BeList.keys()], [BeList[key] is not None for key in BeList.keys()])
    nPeaks = sum(eachT)
    # Extract magnetic induction results from the PlanetGrid
    Benm_nT = PlanetGrid[0, 0].Magnetic.Benm_nT
    # Organize data into a format that can be plotted/saved for plotting
    Bex_nT, Bey_nT, Bez_nT = Benm2absBexyz(Benm_nT)
    induction_data = {
        'nPeaks': nPeaks,
        'Amp': None,
        'Aen': None,
        'phase': None,
        'Bix_nT': None,
        'Biy_nT': None,
        'Biz_nT': None,
        'Bi1xyz_nT': None,
        'rBi1x_nT': None,
        'rBi1y_nT': None,
        'rBi1z_nT': None,
        'iBi1x_nT': None,
        'iBi1y_nT': None,
        'iBi1z_nT': None,
        'Bi1Tot_nT': None,
        'rBi1Tot_nT': None,
        'iBi1Tot_nT': None,
        'calcedExc': None,
        'Texc_hr': None,
    }
    induction_data['calcedExc'] = {}
    induction_data['Texc_hr'] = []
    if nPeaks > 0:
        # Extract amplitude and phase data as 3D arrays (nPeaks x rows x cols)
        # For 1D case, this will be nPeaks x 1 x N, then we'll flatten later
        Amp_3d = np.full((nPeaks, PlanetGrid.shape[0], PlanetGrid.shape[1]), np.nan)
        phase_3d = np.full((nPeaks, PlanetGrid.shape[0], PlanetGrid.shape[1]), np.nan)
        Aen_3d = np.full((nPeaks, PlanetGrid.shape[0], PlanetGrid.shape[1]), np.nan, dtype=np.complex_)
        # Create 3D arrays for Bi1xyz_nT of nPeaks x 2 (since complex number) x rows x cols
        Bi1x_nT_3D = np.full((nPeaks, PlanetGrid.shape[0], PlanetGrid.shape[1]), np.nan, dtype=np.complex_)
        Bi1y_nT_3D = np.full((nPeaks, PlanetGrid.shape[0], PlanetGrid.shape[1]), np.nan, dtype=np.complex_)
        Bi1z_nT_3D = np.full((nPeaks, PlanetGrid.shape[0], PlanetGrid.shape[1]), np.nan, dtype=np.complex_)
        Bi1Tot_3D = np.full((nPeaks, PlanetGrid.shape[0], PlanetGrid.shape[1]), np.nan, dtype=np.complex_)
        for i, line in enumerate(PlanetGrid):
            for j, Planet in enumerate(line):
                if hasattr(Planet, 'Magnetic') and hasattr(Planet.Magnetic, 'Amp') and Planet.Magnetic.Amp is not None:
                    planet_amp = np.array(Planet.Magnetic.Amp)
                    planet_phase = np.array(Planet.Magnetic.phase)
                    planet_Aen = np.array(Planet.Magnetic.Aen[:, 1])
                    Aen_3d[:nPeaks, i, j] = planet_Aen[:nPeaks]
                    Amp_3d[:nPeaks, i, j] = planet_amp[:nPeaks]
                    phase_3d[:nPeaks, i, j] = planet_phase[:nPeaks]
                    Bi1x_nT_3D[:nPeaks, i, j] = Planet.Magnetic.Bi1xyz_nT['x'][:]
                    Bi1y_nT_3D[:nPeaks, i, j] = Planet.Magnetic.Bi1xyz_nT['y'][:]
                    Bi1z_nT_3D[:nPeaks, i, j] = Planet.Magnetic.Bi1xyz_nT['z'][:]
                    Bi1Tot_3D[:nPeaks, i, j] = Planet.Magnetic.Bi1Tot_nT[:]
                    
                    # Set attributes that are consistent across all planets - in some cases where planets are invalid, they will not have this data, so this is why we only reset values for which their magnetic data has been calculated
                    induction_data['calcedExc'] = Planet.Magnetic.calcedExc
                    induction_data['Texc_hr'] = Planet.Magnetic.Texc_hr
                    
        
        induction_data['Amp'] = Amp_3d
        induction_data['phase'] = phase_3d
        induction_data['Aen'] = Aen_3d
        induction_data['Bi1x_nT'] = Bi1x_nT_3D
        induction_data['Bi1y_nT'] = Bi1y_nT_3D
        induction_data['Bi1z_nT'] = Bi1z_nT_3D
        # Calculate induced field components
        induction_data['Bix_nT'] = np.array([Amp_3d[i, ...] * Bex_nT[i] for i in range(nPeaks)])
        induction_data['Biy_nT'] = np.array([Amp_3d[i, ...] * Bey_nT[i] for i in range(nPeaks)])
        induction_data['Biz_nT'] = np.array([Amp_3d[i, ...] * Bez_nT[i] for i in range(nPeaks)])
        induction_data['rBi1x_nT'] = np.real(induction_data['Bi1x_nT'])
        induction_data['rBi1y_nT'] = np.real(induction_data['Bi1y_nT'])
        induction_data['rBi1z_nT'] = np.real(induction_data['Bi1z_nT'])
        induction_data['iBi1x_nT'] = np.imag(induction_data['Bi1x_nT'])
        induction_data['iBi1y_nT'] = np.imag(induction_data['Bi1y_nT'])
        induction_data['iBi1z_nT'] = np.imag(induction_data['Bi1z_nT'])
        induction_data['Bi1Tot_nT'] = Bi1Tot_3D
        induction_data['rBi1Tot_nT'] = np.real(induction_data['Bi1Tot_nT'])
        induction_data['iBi1Tot_nT'] = np.imag(induction_data['Bi1Tot_nT'])

    for key, value in induction_data.items():
        if hasattr(InductionResults, key):
            setattr(InductionResults, key, value)
    
    # Data is now always 2D - consistent across all analysis types
    return InductionResults


def InductionCalced(ExplorationList):
    """
    Check if induction data exists for all results in the list.
    """
    for Exploration in ExplorationList:
        if Exploration.induction.Amp is None:
            return False
    return True

def matlab_safe_key(key_name, max_length=31):
    """
    Create MATLAB-safe variable name by truncating if necessary.
    
    MATLAB variable names must be 31 characters or less. This function
    truncates long names while trying to preserve meaningful information.
    
    Args:
        key_name (str): Original variable name
        max_length (int): Maximum allowed length (default 31 for MATLAB)
        
    Returns:
        str: Truncated variable name safe for MATLAB
    """
    if len(key_name) <= max_length:
        return key_name
    
    # Try to preserve meaningful parts by removing common suffixes/prefixes
    truncated = key_name[:max_length]

    log.debug(f'Truncated MATLAB key "{key_name}" to "{truncated}"')
    return truncated


def flatten_dict_for_matlab(obj, prefix='', flat_dict=None, max_key_length=31):
    """
    Flatten an object's attributes into a MATLAB-compatible dictionary.
    
    This centralizes the flattening logic for all ResultsStruct objects.
    
    Args:
        obj: ResultsStruct object to flatten
        prefix (str): Prefix to add to attribute names (unused for top-level)
        flat_dict (dict): Dictionary to add flattened items to (created if None)
        max_key_length (int): Maximum key length for MATLAB compatibility
        
    Returns:
        dict: Flattened dictionary suitable for savemat()
    """
    if flat_dict is None:
        flat_dict = {}
    
    def _flatten_object_attributes(source_obj, key_prefix):
        """Helper to flatten attributes from a source object with given prefix."""
        for attr_name in dir(source_obj):
            if not attr_name.startswith('_') and not callable(getattr(source_obj, attr_name)):
                attr_value = getattr(source_obj, attr_name)
                if attr_value is not None:
                    full_key = f'{key_prefix}{attr_name}'
                    safe_key = matlab_safe_key(full_key, max_key_length)
                    
                    # Convert data types for MATLAB compatibility
                    if isinstance(attr_value, (bool, np.bool_)):
                        flat_dict[safe_key] = int(attr_value)
                    elif isinstance(attr_value, np.ndarray):
                        flat_dict[safe_key] = attr_value
                    elif isinstance(attr_value, (list, tuple)):
                        flat_dict[safe_key] = np.array(attr_value)
                    elif isinstance(attr_value, dict):
                        # Handle dictionaries by flattening their contents
                        for dict_key, dict_value in attr_value.items():
                            dict_full_key = f'{full_key}_{dict_key}'
                            dict_safe_key = matlab_safe_key(dict_full_key, max_key_length)
                            if isinstance(dict_value, (bool, np.bool_)):
                                flat_dict[dict_safe_key] = int(dict_value)
                            else:
                                flat_dict[dict_safe_key] = dict_value
                    else:
                        flat_dict[safe_key] = attr_value
    possibleNestedResultStructs = {'base', 'induction', 'inversion'}
    nestedResultsStructs = {}
    for nestedResultStruct in possibleNestedResultStructs:
        if hasattr(obj, nestedResultStruct):
            nestedResultsStructs[nestedResultStruct] = getattr(obj, nestedResultStruct)
        
    # Flatten nested objects (base and induction) with their respective prefixes
    for nestedResultStruct in nestedResultsStructs:
        _flatten_object_attributes(nestedResultsStructs[nestedResultStruct], f'{nestedResultStruct}_')
    
    # Handle top-level attributes that are not nested objects
    for attr_name in dir(obj):
        # Skip private attributes (starting with '_'), nested objects, and callable methods
        if (not attr_name.startswith('_') and 
            attr_name not in nestedResultsStructs.keys() and
            not callable(getattr(obj, attr_name))):
            attr_value = getattr(obj, attr_name)
            # Only process attributes that have actual values (not None)
            if attr_value is not None:
                safe_key = matlab_safe_key(attr_name, max_key_length)
                # Special handling for numpy arrays - flatten to 1D for MATLAB compatibility
                if isinstance(attr_value, np.ndarray):
                    flat_dict[safe_key] = attr_value.flatten()
                # Convert boolean values to integers (MATLAB doesn't handle booleans well)
                elif isinstance(attr_value, (bool, np.bool_)):
                    flat_dict[safe_key] = int(attr_value)
                # All other data types can be stored directly
                else:
                    flat_dict[safe_key] = attr_value
            
    return flat_dict
