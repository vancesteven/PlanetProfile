"""Pared-down 5D variant of Test50 MCMC.

Drops Tb_K from sampling (fixed at PPTest50 default = upper grid edge) and
collapses HP-ice viscosities (Ice III/V/VI) into a single log10_eta_HP.

Param order: [alpha, log10_zeta, log10_eta_Ih, log10_eta_HP, log10_eta_sil]

Reuses Test50's structure cache, forward_model, and OBS so the comparison
against the 8D run is apples-to-apples.  Drives the toolkit-generalization
exercise — if the 5D wrapper is clean, the GUI surface that exposes
"which params to sample / hold fixed" can mirror this pattern.
"""
import os
import pickle
import sys
import time

import numpy as np
import pocomc as pc
from scipy.stats import uniform

from PlanetProfile.Test.Test50_mcmc_andrade_noocean_yao2014 import (
    OBS,
    OUTPUT_DIR,
    RANDOM_STATE,
    STRUCTURE_TB_GRID,
    build_or_load_structure_grid,
    forward_model,
)
from PlanetProfile.Inference import mcmc_plots


# Fix Tb at PPTest50 default (upper edge of grid = TtripleIh_III_L_K - 0.2)
TB_FIXED_K = float(STRUCTURE_TB_GRID[-1])

PARAM_NAMES_5D = ['alpha', 'log10_zeta', 'log10_eta_Ih', 'log10_eta_HP',
                  'log10_eta_sil']
PARAM_LABELS_5D = [
    r'$\alpha$',
    r'$\log_{10}\zeta$',
    r'$\log_{10}\eta_\mathrm{Ih}$',
    r'$\log_{10}\eta_\mathrm{HP}$',
    r'$\log_{10}\eta_\mathrm{sil}$',
]
N_DIM_5D = 5
N_EFF_5D = 500


def theta5_to_theta8(theta5):
    """Expand 5D parameter vector to the 8D vector that forward_model expects."""
    alpha, log10_zeta, log10_eta_Ih, log10_eta_HP, log10_eta_sil = theta5
    return np.array([
        alpha, log10_zeta, log10_eta_Ih,
        log10_eta_HP, log10_eta_HP, log10_eta_HP,   # III, V, VI all share η_HP
        log10_eta_sil, TB_FIXED_K,
    ])


def log_likelihood_5D(theta5, structure_grid):
    Re_k2, Im_k2, _ = forward_model(theta5_to_theta8(theta5), structure_grid)
    if not np.isfinite(Re_k2):
        return -1e30
    chi2 = (
        ((Re_k2 - OBS['Re_k2']) / OBS['Re_k2_err']) ** 2
        + ((abs(Im_k2) - OBS['Im_k2']) / OBS['Im_k2_err']) ** 2
    )
    return -0.5 * chi2


def run_5D():
    grid = build_or_load_structure_grid()
    print(f'Tb fixed at {TB_FIXED_K:.3f} K (PPTest50 default)')
    print(f'Sampling {N_DIM_5D}D: {PARAM_NAMES_5D}')

    prior = pc.Prior([
        uniform(loc=0.15, scale=0.30),  # alpha:        [0.15, 0.45]
        uniform(loc=-3.0, scale=5.0),   # log10_zeta:   [-3,    2]
        uniform(loc=10.0, scale=6.0),   # log10_eta_Ih: [10,   16]
        uniform(loc=10.0, scale=6.0),   # log10_eta_HP: [10,   16]
        uniform(loc=18.0, scale=4.0),   # log10_eta_sil:[18,   22]
    ])

    sampler = pc.Sampler(
        prior=prior,
        likelihood=lambda theta: log_likelihood_5D(theta, grid),
        n_effective=N_EFF_5D,
        random_state=RANDOM_STATE,
    )
    t0 = time.time()
    sampler.run()
    elapsed = time.time() - t0
    print(f'\nMCMC completed in {elapsed/60:.1f} min')

    raw_x = sampler.results['x']
    raw_logl = sampler.results['logl']
    if raw_x.ndim == 3:
        samples = raw_x.reshape(-1, raw_x.shape[-1])
        log_prob = raw_logl.reshape(-1)
    else:
        samples = raw_x
        log_prob = raw_logl

    # Posterior summary
    print('\n=== 5D Posterior Summary (median [16th, 84th]) ===')
    qs = np.percentile(samples, [16, 50, 84], axis=0)
    for i, name in enumerate(PARAM_NAMES_5D):
        print(f'  {name:>15s}: {qs[1,i]:8.3f}  [{qs[0,i]:7.3f}, {qs[2,i]:7.3f}]')

    # Forward-evaluate a posterior subsample for k2 distribution
    n_eval = min(500, samples.shape[0])
    rng = np.random.default_rng(RANDOM_STATE)
    idx = rng.choice(samples.shape[0], size=n_eval, replace=False)
    Re_k2_arr, Im_k2_arr = [], []
    for j in idx:
        re, im, _ = forward_model(theta5_to_theta8(samples[j]), grid)
        if np.isfinite(re):
            Re_k2_arr.append(re)
            Im_k2_arr.append(abs(im))
    Re_k2_arr = np.array(Re_k2_arr)
    Im_k2_arr = np.array(Im_k2_arr)
    print(f'\n  Re(k2): {np.median(Re_k2_arr):.4f} '
          f'[{np.percentile(Re_k2_arr,16):.4f}, {np.percentile(Re_k2_arr,84):.4f}]')
    print(f'  |Im(k2)|: {np.median(Im_k2_arr):.4f} '
          f'[{np.percentile(Im_k2_arr,16):.4f}, {np.percentile(Im_k2_arr,84):.4f}]')

    # Save results
    pkl_path = os.path.join(OUTPUT_DIR, 'allice_yao2014_andrade_5D_results.pkl')
    results = {
        'samples': samples,
        'log_prob': log_prob,
        'param_names': PARAM_NAMES_5D,
        'obs': OBS,
        'Tb_fixed_K': TB_FIXED_K,
        'Re_k2_post': Re_k2_arr,
        'Im_k2_post': Im_k2_arr,
    }
    with open(pkl_path, 'wb') as f:
        pickle.dump(results, f)
    print(f'\nSaved {pkl_path}')

    # Corner plot
    corner_path = os.path.join(OUTPUT_DIR,
                               'allice_yao2014_andrade_5D_corner.png')
    mcmc_plots.plot_corner(
        samples=samples,
        labels=PARAM_LABELS_5D,
        title=f'Test50 5D (Andrade no-ocean Titan, Tb={TB_FIXED_K:.2f} K)',
        output_path=corner_path,
    )
    print(f'Saved {corner_path}')

    return results


if __name__ == '__main__':
    run_5D()
