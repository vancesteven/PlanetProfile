"""
Save and reload 3D lateral structure fields.

Stores lateral results as pickle files alongside standard PlanetProfile
output, and provides utilities for reloading and inspecting 3D fields.
"""
import os
import pickle
import numpy as np
import logging

log = logging.getLogger('PlanetProfile.Lateral.LateralIO')


def SaveLateralResults(Planet, Params):
    """ Save 3D lateral structure results to a pickle file.

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

    # Build output dict (only numpy/basic types for portability)
    results = {
        'bodyName': Planet.name,
        'gridType': Lateral.gridType,
        'nPix': Lateral.nPix,
        'theta_rad': Lateral.theta_rad,
        'phi_rad': Lateral.phi_rad,
        'pixArea_sr': Lateral.pixArea_sr,
        'dIce_m': Lateral.dIce_m,
        'fClath': Lateral.fClath,
        'Tb_K': Lateral.Tb_K,
        'qSurf_Wm2': Lateral.qSurf_Wm2,
        'sigma_mean_Sm': Lateral.sigma_mean_Sm,
        'kThermEff_WmK': Lateral.kThermEff_WmK,
        'HtidalIce_Wm3': Lateral.HtidalIce_Wm3,
        'Mtarget_kg': Lateral.Mtarget_kg,
        'Mactual_kg': Lateral.Mactual_kg,
        'massResidual_frac': Lateral.massResidual_frac,
    }

    # Include SH coefficients if available
    if Lateral.dIce_Cpq_km is not None:
        results['dIce_Cpq_km'] = Lateral.dIce_Cpq_km
        results['dIce_Spq_km'] = Lateral.dIce_Spq_km
        results['dIce_pMax'] = Lateral.dIce_pMax

    outPath = os.path.join(outDir, f'{Planet.name}_lateral3D.pkl')
    with open(outPath, 'wb') as f:
        pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)

    log.info(f'Lateral results saved to {outPath}')


def ReloadLateralResults(Planet, fPath):
    """ Reload 3D lateral structure results from a pickle file.

        Args:
            Planet: PlanetStruct to populate Lateral fields on.
            fPath: Path to the lateral3D pickle file.

        Returns:
            Planet: Updated with loaded Lateral fields.
    """
    with open(fPath, 'rb') as f:
        results = pickle.load(f)

    Lateral = Planet.Lateral

    for key in ['gridType', 'nPix', 'theta_rad', 'phi_rad', 'pixArea_sr',
                'dIce_m', 'fClath', 'Tb_K', 'qSurf_Wm2', 'sigma_mean_Sm',
                'kThermEff_WmK', 'HtidalIce_Wm3', 'Mtarget_kg', 'Mactual_kg',
                'massResidual_frac']:
        if key in results:
            setattr(Lateral, key, results[key])

    if 'dIce_Cpq_km' in results:
        Lateral.dIce_Cpq_km = results['dIce_Cpq_km']
        Lateral.dIce_Spq_km = results['dIce_Spq_km']
        Lateral.dIce_pMax = results['dIce_pMax']

    log.info(f'Lateral results loaded from {fPath}')
    return Planet
