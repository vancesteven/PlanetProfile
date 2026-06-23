"""
Save and reload 3D lateral structure fields.

Stores lateral results as NPZ files alongside standard PlanetProfile
output, and provides utilities for reloading and inspecting 3D fields.
"""
import os
import pickle
import numpy as np
import logging

log = logging.getLogger('PlanetProfile.Lateral.LateralIO')


def SaveLateralResults(Planet, Params):
    """ Save 3D lateral structure results to an NPZ file.

        Output file is placed alongside the standard PlanetProfile output,
        with '_lateral3D' appended to the filename.

        Args:
            Planet: PlanetStruct with completed Lateral substruct.
            Params: ParamsStruct.
    """
    Lateral = Planet.Lateral
    if Lateral.dIce_m is None:
        log.debug('No lateral results to save.')
        return

    # Determine output directory
    if hasattr(Params, 'FigureFiles') and Params.FigureFiles is not None:
        outDir = Params.FigureFiles.path
    else:
        outDir = os.path.join(os.getcwd(), Planet.bodyname, 'lateral3D')

    os.makedirs(outDir, exist_ok=True)

    # Build output dict with numpy arrays
    results = {}

    # Pack scalar/string metadata into a 0-d object array
    meta_dict = {
        'bodyName': Planet.name,
        'gridType': Lateral.gridType,
        'nPix': Lateral.nPix,
        'Mtarget_kg': Lateral.Mtarget_kg,
        'Mactual_kg': Lateral.Mactual_kg,
        'massResidual_frac': Lateral.massResidual_frac,
        'repairLog': getattr(Lateral, 'repairLog', None),
        'maxRepairFrac': getattr(Lateral, 'maxRepairFrac', None),
    }
    if Lateral.dIce_Cpq_km is not None and hasattr(Lateral, 'dIce_pMax'):
        meta_dict['dIce_pMax'] = Lateral.dIce_pMax
    results['_meta'] = np.array([meta_dict], dtype=object)

    # Store numpy arrays (skip None values)
    array_fields = ['theta_rad', 'phi_rad', 'pixArea_sr', 'dIce_m', 'Tb_K',
                    'qSurf_Wm2', 'sigma_mean_Sm', 'kThermEff_WmK',
                    'HtidalIce_Wm3', 'fClath', 'dIce_Cpq_km', 'dIce_Spq_km',
                    'failedColumnMask', 'repairedColumnMask']
    for field in array_fields:
        if hasattr(Lateral, field):
            val = getattr(Lateral, field)
            if val is not None:
                results[field] = val

    outPath = os.path.join(outDir, f'{Planet.name}_lateral3D.npz')
    np.savez_compressed(outPath, **results)

    log.info(f'Lateral results saved to {outPath}')


def ReloadLateralResults(Planet, fPath):
    """ Reload 3D lateral structure results from an NPZ or legacy pickle file.

        Args:
            Planet: PlanetStruct to populate Lateral fields on.
            fPath: Path to the lateral3D .npz or .pkl file.

        Returns:
            Planet: Updated with loaded Lateral fields.
    """
    Lateral = Planet.Lateral

    # Detect file format
    if fPath.endswith('.pkl'):
        # Legacy pickle format
        log.warning('Loading legacy .pkl lateral file...')
        with open(fPath, 'rb') as f:
            results = pickle.load(f)

        for key in ['gridType', 'nPix', 'theta_rad', 'phi_rad', 'pixArea_sr',
                    'dIce_m', 'fClath', 'Tb_K', 'qSurf_Wm2', 'sigma_mean_Sm',
                    'kThermEff_WmK', 'HtidalIce_Wm3', 'Mtarget_kg', 'Mactual_kg',
                    'massResidual_frac', 'failedColumnMask', 'repairedColumnMask',
                    'repairLog', 'maxRepairFrac']:
            if key in results:
                setattr(Lateral, key, results[key])

        if 'dIce_Cpq_km' in results:
            Lateral.dIce_Cpq_km = results['dIce_Cpq_km']
            Lateral.dIce_Spq_km = results['dIce_Spq_km']
            Lateral.dIce_pMax = results['dIce_pMax']

    else:
        # NPZ format
        results = np.load(fPath, allow_pickle=True)

        # Restore metadata
        if '_meta' in results:
            meta_dict = results['_meta'].item()
            for key in ['gridType', 'nPix', 'Mtarget_kg', 'Mactual_kg',
                        'massResidual_frac', 'dIce_pMax', 'repairLog',
                        'maxRepairFrac']:
                if key in meta_dict:
                    setattr(Lateral, key, meta_dict[key])

        # Restore array fields
        for key in results.keys():
            if key != '_meta':
                setattr(Lateral, key, results[key])

    log.info(f'Lateral results loaded from {fPath}')
    return Planet
