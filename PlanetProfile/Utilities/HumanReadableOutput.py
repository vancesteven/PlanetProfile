"""
Functions for creating CSV output files from PlanetProfile results.

This module provides CSV export of profile data for easier analysis in
spreadsheet programs, complementing the full-precision .txt files.
"""

import pandas as pd
import numpy as np
import os


def WriteProfileCSV(Planet, Params, outputPath=None, precision=8):
    """
    Write profile data to CSV format with reasonable precision.

    This creates a spreadsheet-friendly version of the main profile data,
    with the same columns as the .txt file but in CSV format for easy
    import into Excel, Google Sheets, etc.

    Args:
        Planet: PlanetStruct object with results
        Params: ParamsStruct with configuration
        outputPath: Optional path override (defaults to DataFiles.saveFile with .csv extension)
        precision: Number of significant figures (default 8)

    Returns:
        Path to created CSV file
    """
    if outputPath is None:
        outputPath = Params.DataFiles.saveFile.replace('.txt', '.csv')

    # Create DataFrame with all profile data
    data = {
        'P_MPa': Planet.P_MPa[:Planet.Steps.nTotal],
        'T_K': Planet.T_K[:Planet.Steps.nTotal],
        'r_m': Planet.r_m[:Planet.Steps.nTotal],
        'phase': Planet.phase[:Planet.Steps.nTotal],
        'rho_kgm3': Planet.rho_kgm3[:Planet.Steps.nTotal],
        'Cp_JkgK': Planet.Cp_JkgK[:Planet.Steps.nTotal],
        'alpha_pK': Planet.alpha_pK[:Planet.Steps.nTotal],
        'g_ms2': Planet.g_ms2[:Planet.Steps.nTotal],
        'phi_frac': Planet.phi_frac[:Planet.Steps.nTotal],
        'sigma_Sm': Planet.sigma_Sm[:Planet.Steps.nTotal],
        'kTherm_WmK': Planet.kTherm_WmK[:Planet.Steps.nTotal],
        'VP_kms': Planet.Seismic.VP_kms[:Planet.Steps.nTotal],
        'VS_kms': Planet.Seismic.VS_kms[:Planet.Steps.nTotal],
        'QS': Planet.Seismic.QS[:Planet.Steps.nTotal],
        'KS_GPa': Planet.Seismic.KS_GPa[:Planet.Steps.nTotal],
        'GS_GPa': Planet.Seismic.GS_GPa[:Planet.Steps.nTotal],
        'Ppore_MPa': Planet.Ppore_MPa[:Planet.Steps.nTotal],
        'rhoMatrix_kgm3': Planet.rhoMatrix_kgm3[:Planet.Steps.nTotal],
        'rhoPore_kgm3': Planet.rhoPore_kgm3[:Planet.Steps.nTotal],
        'MLayer_kg': Planet.MLayer_kg[:Planet.Steps.nTotal],
        'VLayer_m3': Planet.VLayer_m3[:Planet.Steps.nTotal],
        'Htidal_Wm3': Planet.Htidal_Wm3[:Planet.Steps.nTotal],
        'eta_Pas': Planet.eta_Pas[:Planet.Steps.nTotal]
    }

    df = pd.DataFrame(data)

    # Round to reasonable precision (default 8 significant figures)
    # Use float_format in to_csv for consistent formatting
    df.to_csv(outputPath, index=False, float_format=f'%.{precision}g')

    return outputPath


def WriteProfileSpreadsheet(Planet, Params, outputPath=None):
    """
    Write profile data to Excel format (.xlsx) with metadata sheet.

    This creates an Excel workbook with two sheets:
    1. "Profile" - The profile data (same as CSV)
    2. "Metadata" - Header information from the model

    Args:
        Planet: PlanetStruct object with results
        Params: ParamsStruct with configuration
        outputPath: Optional path override

    Returns:
        Path to created Excel file
    """
    if outputPath is None:
        outputPath = Params.DataFiles.saveFile.replace('.txt', '.xlsx')

    try:
        import openpyxl
    except ImportError:
        raise ImportError('openpyxl required for Excel output. Install with: pip install openpyxl')

    # Create profile data
    data = {
        'P_MPa': Planet.P_MPa[:Planet.Steps.nTotal],
        'T_K': Planet.T_K[:Planet.Steps.nTotal],
        'r_m': Planet.r_m[:Planet.Steps.nTotal],
        'phase': Planet.phase[:Planet.Steps.nTotal],
        'rho_kgm3': Planet.rho_kgm3[:Planet.Steps.nTotal],
        'Cp_JkgK': Planet.Cp_JkgK[:Planet.Steps.nTotal],
        'alpha_pK': Planet.alpha_pK[:Planet.Steps.nTotal],
        'g_ms2': Planet.g_ms2[:Planet.Steps.nTotal],
        'phi_frac': Planet.phi_frac[:Planet.Steps.nTotal],
        'sigma_Sm': Planet.sigma_Sm[:Planet.Steps.nTotal],
        'kTherm_WmK': Planet.kTherm_WmK[:Planet.Steps.nTotal],
        'VP_kms': Planet.Seismic.VP_kms[:Planet.Steps.nTotal],
        'VS_kms': Planet.Seismic.VS_kms[:Planet.Steps.nTotal],
        'QS': Planet.Seismic.QS[:Planet.Steps.nTotal],
        'KS_GPa': Planet.Seismic.KS_GPa[:Planet.Steps.nTotal],
        'GS_GPa': Planet.Seismic.GS_GPa[:Planet.Steps.nTotal],
        'Ppore_MPa': Planet.Ppore_MPa[:Planet.Steps.nTotal],
        'rhoMatrix_kgm3': Planet.rhoMatrix_kgm3[:Planet.Steps.nTotal],
        'rhoPore_kgm3': Planet.rhoPore_kgm3[:Planet.Steps.nTotal],
        'MLayer_kg': Planet.MLayer_kg[:Planet.Steps.nTotal],
        'VLayer_m3': Planet.VLayer_m3[:Planet.Steps.nTotal],
        'Htidal_Wm3': Planet.Htidal_Wm3[:Planet.Steps.nTotal],
        'eta_Pas': Planet.eta_Pas[:Planet.Steps.nTotal]
    }
    df_profile = pd.DataFrame(data)

    # Create metadata sheet
    metadata = {
        'Parameter': [
            'Model Label',
            'PlanetProfile Version',
            'Iron Core',
            'Ocean Composition',
            'Salinity (ppt)',
            'Radius (m)',
            'Mass (kg)',
            'C/MR² measured',
            'C/MR² calculated',
            'Surface Temperature (K)',
            'Bottom Ice Temperature (K)',
            'Ice Shell Thickness (km)',
            'Ocean Thickness (km)',
            'Bottom Ice Pressure (MPa)',
            'Total Steps'
        ],
        'Value': [
            Planet.label,
            getattr(Params, 'version', 'N/A'),
            Planet.Do.Fe_CORE,
            Planet.Ocean.comp if not Planet.Do.NO_H2O else 'None',
            Planet.Ocean.wOcean_ppt if not Planet.Do.NO_H2O else 0,
            Planet.Bulk.R_m,
            Planet.Bulk.M_kg,
            Planet.Bulk.Cmeasured,
            Planet.CMR2mean,
            Planet.Bulk.Tsurf_K,
            Planet.Bulk.Tb_K,
            Planet.zb_km,
            Planet.D_km,
            Planet.Pb_MPa,
            Planet.Steps.nTotal
        ]
    }
    df_metadata = pd.DataFrame(metadata)

    # Write to Excel with multiple sheets
    with pd.ExcelWriter(outputPath, engine='openpyxl') as writer:
        df_metadata.to_excel(writer, sheet_name='Metadata', index=False)
        df_profile.to_excel(writer, sheet_name='Profile', index=False)

    return outputPath
