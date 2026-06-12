"""
Regenerate structure-wedge plots for all MCMC tests that have existing pkl results.
Does NOT import TidalPy — safe to run in any environment.

Usage:
    python PlanetProfile/Test/replot_wedges.py
"""
import logging
import os
import pickle
import sys

import matplotlib
matplotlib.use('Agg')

_pp_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _pp_root)

logging.basicConfig(level=logging.INFO, format='%(name)s - %(message)s')
log = logging.getLogger('replot_wedges')

import importlib.util as _ilu
_spec = _ilu.spec_from_file_location(
    'mcmc_plots',
    os.path.join(_pp_root, 'PlanetProfile', 'Inference', 'mcmc_plots.py')
)
mp = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(mp)

_BASE = os.path.join(_pp_root, 'PlanetProfile', 'Test', 'mcmc_results')

TITAN_R_KM = 2575.5
EUROPA_R_KM = 1560.8
GANYMEDE_R_KM = 2634.1
CALLISTO_R_KM = 2410.3


def _load(path):
    with open(path, 'rb') as f:
        return pickle.load(f)


def _load_grid(path):
    """Load grid cache pkl, unwrapping the 'grid_cache' envelope if present."""
    d = _load(path)
    if isinstance(d, dict) and 'grid_cache' in d and isinstance(d['grid_cache'], dict):
        return d['grid_cache']
    return d


def replot_titan_test46():
    pkl    = os.path.join(_BASE, 'Titan', 'Test46_andrade_hybrid_hydro', 'hybrid_hydro_andrade_mcmc.pkl')
    grid   = os.path.join(_BASE, 'Titan', 'Test46_andrade_hybrid_hydro', 'titan_maxwell_hybrid_hydro_grid_cache.pkl')
    outdir = os.path.join(_BASE, 'Titan', 'Test46_andrade_hybrid_hydro')
    if not (os.path.exists(pkl) and os.path.exists(grid)):
        log.warning('Test46 andrade hybrid: pkls not found, skipping.')
        return
    r = _load(pkl)
    grid_cache = _load_grid(grid)
    mp.plot_structure_wedge(
        r['samples'], r['eval_idx'], grid_cache,
        output_path=os.path.join(outdir, 'hybrid_hydro_andrade_structure_profile.png'),
        R_body_km=TITAN_R_KM,
        body_name='Titan',
        param_indices={'Tb': 5, 'D_hydro': 6, 'f_core': 9},
    )


def replot_titan_test48():
    pkl    = os.path.join(_BASE, 'Titan', 'Test48_andrade_yao2014', 'hybrid_hydro_andrade_yao2014_mcmc.pkl')
    grid   = os.path.join(_BASE, 'Titan', 'Test48_andrade_yao2014', 'titan_yao2014_hybrid_hydro_grid_cache.pkl')
    outdir = os.path.join(_BASE, 'Titan', 'Test48_andrade_yao2014')
    if not (os.path.exists(pkl) and os.path.exists(grid)):
        log.warning('Test48 yao2014: pkls not found, skipping.')
        return
    r = _load(pkl)
    grid_cache = _load_grid(grid)
    mp.plot_structure_wedge(
        r['samples'], r['eval_idx'], grid_cache,
        output_path=os.path.join(outdir, 'hybrid_hydro_andrade_yao2014_structure_profile.png'),
        R_body_km=TITAN_R_KM,
        body_name='Titan',
        param_indices={'Tb': 5, 'D_hydro': 6, 'f_core': 9},
    )


def replot_titan_test49(tag, pkl_name, grid_name, outname):
    subdir = 'Test49_clathrate4km' if '4km' in tag else 'Test49_clathrate2km'
    dirpath = os.path.join(_BASE, 'Titan', subdir)
    pkl  = os.path.join(dirpath, pkl_name)
    grid = os.path.join(dirpath, grid_name)
    if not (os.path.exists(pkl) and os.path.exists(grid)):
        log.warning(f'{tag}: pkls not found, skipping.')
        return
    import numpy as np
    r = _load(pkl)
    grid_cache = _load_grid(grid)
    samples = r['samples']
    eval_idx = r.get('eval_idx', np.arange(len(samples)))
    mp.plot_structure_wedge(
        samples, eval_idx, grid_cache,
        output_path=os.path.join(dirpath, outname),
        R_body_km=TITAN_R_KM,
        body_name='Titan',
        param_indices={'Tb': 5, 'D_hydro': 6, 'f_core': 9},
    )


def replot_europa_test51():
    dirpath = os.path.join(_BASE, 'Europa', 'Test51_seawater')
    pkl  = os.path.join(dirpath, 'test51_europa_seawater_results.pkl')
    grid = os.path.join(dirpath, 'europa_seawater_structure_grid.pkl')
    if not os.path.exists(grid):
        log.warning('Test51 Europa: structure grid not found, skipping.')
        return
    if not os.path.exists(pkl):
        log.warning('Test51 Europa: results pkl not found, skipping.')
        return
    try:
        r = _load(pkl)
        grid_cache = _load_grid(grid)
    except Exception as e:
        log.warning(f'Test51 Europa: load failed ({e}), skipping.')
        return
    # Probe pkl structure — handle both dict and InferenceResult
    if isinstance(r, dict):
        samples  = r['samples']
        eval_idx = r.get('eval_idx', range(len(samples)))
        p_idx    = r.get('param_indices', {'Tb': 5, 'D_hydro': 6, 'f_core': 9})
    else:
        samples  = r.samples
        eval_idx = getattr(r, 'eval_idx', range(len(r.samples)))
        p_idx    = {'Tb': 5, 'D_hydro': 6, 'f_core': 9}
    mp.plot_structure_wedge(
        samples, eval_idx, grid_cache,
        output_path=os.path.join(dirpath, 'europa_test51_structure_wedge.png'),
        R_body_km=EUROPA_R_KM,
        body_name='Europa',
        param_indices=p_idx,
    )


if __name__ == '__main__':
    log.info('Regenerating structure-wedge plots with updated PP colors...')

    replot_titan_test46()
    replot_titan_test48()

    replot_titan_test49(
        tag='Test49_4km',
        pkl_name='test49_clathrate4km_mcmc_results.pkl',
        grid_name='titan_yao2014_clathrate4km_hybrid_hydro_grid_cache.pkl',
        outname='test49_clathrate4km_structure_wedge.png',
    )
    replot_titan_test49(
        tag='Test49_2km',
        pkl_name='test49_clathrate2km_mcmc_results.pkl',
        grid_name='titan_yao2014_clathrate2km_hybrid_hydro_grid_cache.pkl',
        outname='test49_clathrate2km_structure_wedge.png',
    )

    replot_europa_test51()

    log.info('Done.')
