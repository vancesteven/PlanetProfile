"""
MCMC exploration: Andrade rheology, no-ocean Titan (PPTest41).

Explores the posterior distribution of Andrade parameters and viscosities
that reproduce the observed tidal Love number k2 for Titan, using the
no-ocean model from PPTest41.

Parameter space (5D, from Petricca et al. 2025 Extended Data Table 2):
  alpha:         Andrade exponent           [0.2, 0.4]
  log10(zeta):   Andrade compliance         [-2, 2]
  log10(eta_Ih): Ice Ih viscosity (Pa s)    [12, 16]
  log10(eta_HP): HP ice viscosity (Pa s)    [10, 18]
  log10(eta_sil): Silicate viscosity (Pa s) [18, 22]

Observational constraints (Petricca et al. 2025):
  Re(k2) = 0.608 +/- 0.048
  |Im(k2)| = 0.135 +/- 0.035

Outputs:
  - MCMC chain and per-phase heating arrays saved to pickle for re-analysis
  - Corner plot of posterior parameters
  - Re(k2) vs Im(k2) scatter with error ellipse
  - Per-phase heating distribution across the posterior

Usage:
  python PlanetProfile/Test/Test41_mcmc_andrade_no_ocean.py
"""
import numpy as np
import os
import sys
import time
import pickle
import logging
import importlib

# --- Environment setup ---
import platformdirs
_pp_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
platformdirs.user_documents_dir = lambda: _pp_root
sys.path.insert(0, _pp_root)

logging.basicConfig(level=logging.WARNING, format='%(name)s - %(message)s')
log = logging.getLogger('PPTest41_MCMC')
log.setLevel(logging.INFO)

from PlanetProfile.Utilities.defineStructs import EOSlist, Constants
from PlanetProfile.Utilities.Indexing import PhaseConv
from PlanetProfile.Main import PlanetProfile as RunPP
from PlanetProfile.Gravity.Gravity import SetupGravity

# TidalPy imports
from TidalPy.RadialSolver import radial_solver, build_rs_input_from_data
from TidalPy.tides.multilayer.heating import calc_radial_volumetric_tidal_heating_from_rs_solution
from TidalPy.rheology import Andrade, Maxwell, Elastic, Newton

# ============================================================
# Configuration
# ============================================================
OUTPUT_DIR = os.path.join(_pp_root, 'PlanetProfile', 'Test', 'mcmc_results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Observational constraints (Petricca et al. 2025)
RE_K2_OBS = 0.608
RE_K2_ERR = 0.048
IM_K2_OBS = 0.135  # absolute value
IM_K2_ERR = 0.035

# MCMC settings
N_EFF = 500       # Target effective sample size for pocoMC
RANDOM_STATE = 42
N_REEVAL = 500    # Number of posterior samples to re-evaluate for heating

# Parameter labels
PARAM_NAMES = ['alpha', 'log10_zeta', 'log10_eta_Ih', 'log10_eta_HP', 'log10_eta_sil']
PARAM_LABELS = [r'$\alpha$', r'$\log_{10}\zeta$', r'$\log_{10}\eta_\mathrm{Ih}$',
                r'$\log_{10}\eta_\mathrm{HP}$', r'$\log_{10}\eta_\mathrm{sil}$']
N_DIM = 5

# ============================================================
# Step 1: Build structural model (one-time)
# ============================================================
def build_structural_model():
    """Run PlanetProfile once to get the fixed interior structure."""
    log.info('Building structural model from PPTest41...')
    t0 = time.time()

    EOSlist.loaded.clear()
    from PlanetProfile.GetConfig import Params as configParams

    test_module = 'PlanetProfile.Test.PPTest41'
    if test_module in sys.modules:
        importlib.reload(sys.modules[test_module])
    else:
        importlib.import_module(test_module)
    mod = sys.modules[test_module]
    Planet = mod.Planet

    configParams.Gravity.backend = 'tidalpy'
    configParams.CALC_NEW = True
    configParams.CALC_NEW_GRAVITY = True
    configParams.NO_SAVEFILE = True
    configParams.SKIP_PLOTS = True

    # Run PlanetProfile to get self-consistent structure
    Planet, Params = RunPP(Planet, configParams)
    Params.CALC_NEW_GRAVITY = True
    Planet, Params = SetupGravity(Planet, Params)

    # Extract and cache structural arrays
    model = Planet.Gravity.ALMAModel['model']
    cols = Planet.Gravity.columns
    rIndex = cols.index('r')
    rhoIndex = cols.index('rho')
    VPIndex = cols.index('VP')
    GSIndex = cols.index('GS')
    etaIndex = cols.index('eta')
    pIndex = cols.index('phase')

    r_m = model[:, rIndex].astype(np.float64)
    rho = model[:, rhoIndex].astype(np.float64)
    mu_Pa = model[:, GSIndex].astype(np.float64)
    VP_ms = model[:, VPIndex].astype(np.float64)
    eta_Pa_base = model[:, etaIndex].astype(np.float64)
    phases = model[:, pIndex]

    # Bulk modulus
    K_Pa = rho * VP_ms**2 - (4.0 / 3.0) * mu_Pa
    nan_mask = ~np.isfinite(K_Pa) | (K_Pa <= 0)
    if np.any(nan_mask):
        for i in np.where(nan_mask)[0]:
            ph = int(phases[i])
            if ph >= 50 and ph < 100:
                nu = 0.25
            elif ph >= 100:
                nu = 0.29
            else:
                nu = 0.33
            K_Pa[i] = 2.0 * mu_Pa[i] * (1.0 + nu) / (3.0 * (1.0 - 2.0 * nu))
    K_Pa = np.maximum(K_Pa, 1e6)

    # Layer boundaries
    changeIndices = np.max(Planet.Reduced.changeIndices) - np.flipud(Planet.Reduced.changeIndices)
    n_layers = len(changeIndices) - 1

    # Thin-layer padding (same logic as _run_tidalpy_backend)
    MIN_POINTS = 5
    needs_padding = any(
        changeIndices[i+1] - changeIndices[i] < MIN_POINTS
        for i in range(n_layers)
    )

    # Compute region_phases from original indices
    _orig_iConv = np.flipud(Planet.Reduced.iConv)
    _orig_changeIndices = changeIndices.copy()
    region_phases = []
    for i_layer in range(n_layers):
        start = _orig_changeIndices[i_layer]
        phase = phases[start]
        if phase >= Constants.phaseClath and phase < Constants.phaseClath + 10:
            phase = Constants.phaseClath
        convection = _orig_iConv[start]
        phase_str = PhaseConv(phase, liq='0')
        if convection:
            phase_str += '_conv'
        region_phases.append(phase_str)

    bulk_visc = np.zeros_like(eta_Pa_base)

    if needs_padding:
        new_r, new_rho, new_K, new_mu, new_eta, new_phases, new_bv = \
            [], [], [], [], [], [], []
        new_ci = [0]
        for i_layer in range(n_layers):
            s, e = changeIndices[i_layer], changeIndices[i_layer + 1]
            n_pts = e - s
            if n_pts < MIN_POINTS and n_pts >= 2:
                r_layer = r_m[s:e]
                r_interp = np.linspace(r_layer[0], r_layer[-1], MIN_POINTS)
                new_r.append(r_interp)
                new_rho.append(np.interp(r_interp, r_layer, rho[s:e]))
                new_K.append(np.interp(r_interp, r_layer, K_Pa[s:e]))
                new_mu.append(np.interp(r_interp, r_layer, mu_Pa[s:e]))
                new_eta.append(np.interp(r_interp, r_layer, eta_Pa_base[s:e]))
                new_phases.append(np.full(MIN_POINTS, phases[s]))
                new_bv.append(np.zeros(MIN_POINTS))
                new_ci.append(new_ci[-1] + MIN_POINTS)
            else:
                new_r.append(r_m[s:e])
                new_rho.append(rho[s:e])
                new_K.append(K_Pa[s:e])
                new_mu.append(mu_Pa[s:e])
                new_eta.append(eta_Pa_base[s:e])
                new_phases.append(phases[s:e])
                new_bv.append(bulk_visc[s:e])
                new_ci.append(new_ci[-1] + (e - s))

        r_m = np.concatenate(new_r)
        rho = np.concatenate(new_rho)
        K_Pa = np.concatenate(new_K)
        mu_Pa = np.concatenate(new_mu)
        eta_Pa_base = np.concatenate(new_eta)
        phases = np.concatenate(new_phases)
        bulk_visc = np.concatenate(new_bv)
        changeIndices = np.array(new_ci)

    # Layer metadata
    layer_upper_radii = []
    layer_types = []
    for i_layer in range(n_layers):
        end = changeIndices[i_layer + 1]
        layer_upper_radii.append(r_m[end - 1])
        layer_types.append('liquid' if phases[changeIndices[i_layer]] == 0 else 'solid')

    # Orbital parameters
    omega = Planet.Bulk.meanMotion_radps
    ecc = Planet.Bulk.eccentricity
    host_mass = Constants.parentMass_kg[Planet.parent]
    a_m = (Constants.G * host_mass / omega**2) ** (1.0 / 3.0)

    # Phase map for heating aggregation
    phase_map = {0: '0', 1: 'Ih', 2: 'II', 3: 'III', 5: 'V', 6: 'VI'}

    data = {
        'r_m': np.ascontiguousarray(r_m),
        'rho': np.ascontiguousarray(rho),
        'K_Pa': np.ascontiguousarray(K_Pa),
        'mu_Pa': np.ascontiguousarray(mu_Pa),
        'eta_Pa_base': eta_Pa_base,
        'phases': phases,
        'bulk_visc': np.ascontiguousarray(bulk_visc),
        'changeIndices': changeIndices,
        'n_layers': n_layers,
        'layer_upper_radii': tuple(layer_upper_radii),
        'layer_types': tuple(layer_types),
        'region_phases': region_phases,
        'omega': omega,
        'eccentricity': ecc,
        'host_mass': host_mass,
        'a_m': a_m,
        'phase_map': phase_map,
    }

    log.info(f'Structure cached: {n_layers} layers, {len(r_m)} radial points, '
             f'built in {time.time()-t0:.0f}s')
    return data


# ============================================================
# Step 2: Fast forward model
# ============================================================
def forward_model(theta, data):
    """Compute k2 and per-phase heating for given parameters.

    Args:
        theta: [alpha, log10_zeta, log10_eta_Ih, log10_eta_HP, log10_eta_sil]
        data: cached structural arrays from build_structural_model()

    Returns:
        (Re_k2, Im_k2, perPhase_W): Real and imaginary k2, per-phase heating dict (W)
    """
    alpha, log10_zeta, log10_eta_Ih, log10_eta_HP, log10_eta_sil = theta
    zeta_pa = 10**log10_zeta
    eta_Ih = 10**log10_eta_Ih
    eta_HP = 10**log10_eta_HP
    eta_sil = 10**log10_eta_sil

    # Modify viscosities per-phase
    eta_mod = data['eta_Pa_base'].copy()
    phases = data['phases']
    changeIndices = data['changeIndices']
    n_layers = data['n_layers']

    for i_layer in range(n_layers):
        s, e = changeIndices[i_layer], changeIndices[i_layer + 1]
        ph = int(phases[s])
        if ph == 1:  # Ice Ih
            eta_mod[s:e] = eta_Ih
        elif ph in (3, 5, 6):  # HP ices III, V, VI
            eta_mod[s:e] = eta_HP
        elif ph >= 50 and ph < 100:  # Silicate
            eta_mod[s:e] = eta_sil
        # Fe, Clath, liquid: keep baseline

    # Build rheology objects
    zeta_tp = zeta_pa ** (1.0 / alpha) if zeta_pa != 1.0 else 1.0
    shear_models = []
    for rp in data['region_phases']:
        base_phase = rp.replace('_conv', '')
        if base_phase in ('0',):
            shear_models.append(Elastic())
        elif base_phase in ('Clath',):
            shear_models.append(Elastic())
        else:
            shear_models.append(Andrade(args=(alpha, zeta_tp)))
    bulk_models = [Elastic() for _ in shear_models]

    try:
        build_data = build_rs_input_from_data(
            data['omega'],
            data['r_m'],
            data['rho'],
            data['K_Pa'],
            data['mu_Pa'],
            data['bulk_visc'],
            np.ascontiguousarray(eta_mod),
            data['layer_upper_radii'],
            data['layer_types'],
            tuple([False] * n_layers),
            tuple([False] * n_layers),
            tuple(shear_models),
            tuple(bulk_models),
            perform_checks=False,
            warnings=False
        )

        result = radial_solver(
            build_data.radius_array,
            build_data.density_array,
            build_data.complex_bulk_modulus_array,
            build_data.complex_shear_modulus_array,
            build_data.frequency,
            build_data.planet_bulk_density,
            build_data.layer_types,
            build_data.is_static_bylayer,
            build_data.is_incompressible_bylayer,
            build_data.upper_radius_bylayer_array,
            degree_l=2,
            solve_for=('tidal',),
            verbose=False,
            raise_on_fail=False
        )

        if not result.success:
            return np.nan, np.nan, {}

        k2 = complex(result.k)
        Re_k2 = k2.real
        Im_k2 = k2.imag

        # Per-phase heating
        perPhase_W = {}
        if data['eccentricity'] > 0:
            heating_profile = calc_radial_volumetric_tidal_heating_from_rs_solution(
                data['eccentricity'], data['omega'], data['a_m'],
                data['host_mass'], result, perform_checks=False
            )
            result_radii = np.asarray(result.radius_array)
            r_m_local = data['r_m']

            for i_layer in range(n_layers):
                s_idx = changeIndices[i_layer]
                e_idx = changeIndices[i_layer + 1]
                layer_phase = int(phases[s_idx])
                phase_str = data['phase_map'].get(layer_phase,
                    PhaseConv(layer_phase, liq='0'))

                r_lo = r_m_local[s_idx]
                r_hi = r_m_local[e_idx - 1]
                mask = (result_radii >= r_lo - 1.0) & (result_radii <= r_hi + 1.0)

                if np.any(mask):
                    lr = result_radii[mask]
                    lh = heating_profile[mask]
                    if len(lr) > 1:
                        total_power = np.trapz(lh * 4.0 * np.pi * lr**2, lr)
                    else:
                        V = (4.0/3.0) * np.pi * (r_hi**3 - r_lo**3)
                        total_power = lh[0] * V
                    perPhase_W[phase_str] = perPhase_W.get(phase_str, 0) + total_power

        return Re_k2, Im_k2, perPhase_W

    except Exception:
        return np.nan, np.nan, {}


# ============================================================
# Step 3: Likelihood and prior
# ============================================================
def log_likelihood(theta, data):
    """Gaussian log-likelihood on Re(k2) and |Im(k2)|."""
    Re_k2, Im_k2, _ = forward_model(theta, data)
    if np.isnan(Re_k2):
        return -1e30
    chi2 = ((Re_k2 - RE_K2_OBS) / RE_K2_ERR)**2 + \
           ((abs(Im_k2) - IM_K2_OBS) / IM_K2_ERR)**2
    return -0.5 * chi2


# ============================================================
# Step 4: Run MCMC
# ============================================================
def run_mcmc(data):
    """Run pocoMC sampler and return results."""
    import pocomc as pc
    from scipy.stats import uniform

    log.info(f'Starting pocoMC MCMC with {N_DIM} parameters, n_eff={N_EFF}')

    prior = pc.Prior([
        uniform(loc=0.2, scale=0.2),       # alpha: [0.2, 0.4]
        uniform(loc=-2.0, scale=4.0),      # log10_zeta: [-2, 2]
        uniform(loc=12.0, scale=4.0),      # log10_eta_Ih: [12, 16]
        uniform(loc=10.0, scale=8.0),      # log10_eta_HP: [10, 18]
        uniform(loc=18.0, scale=4.0),      # log10_eta_sil: [18, 22]
    ])

    def _log_like(theta):
        return log_likelihood(theta, data)

    t0 = time.time()
    sampler = pc.Sampler(
        prior=prior,
        likelihood=_log_like,
        n_effective=N_EFF,
        random_state=RANDOM_STATE,
    )
    sampler.run()
    elapsed = time.time() - t0

    # posterior() returns (samples, logl, logp, logw) with importance resampling
    samples, log_likes, logp, logw = sampler.posterior()

    log.info(f'MCMC completed in {elapsed/60:.1f} min, {len(samples)} samples')
    log.info(f'  Best log_like: {np.max(log_likes):.2f}')

    return samples, log_likes, sampler


# ============================================================
# Step 5: Re-evaluate heating on posterior
# ============================================================
def evaluate_heating(samples, data, n_eval=None):
    """Re-run forward model on posterior samples to get per-phase heating."""
    if n_eval is None:
        n_eval = min(len(samples), N_REEVAL)
    idx = np.random.choice(len(samples), n_eval, replace=False)
    idx.sort()

    heating_results = []
    k2_results = []
    log.info(f'Re-evaluating {n_eval} posterior samples for heating...')
    t0 = time.time()

    for i, si in enumerate(idx):
        theta = samples[si]
        Re_k2, Im_k2, perPhase_W = forward_model(theta, data)
        k2_results.append((Re_k2, Im_k2))
        heating_results.append(perPhase_W)
        if (i + 1) % 100 == 0:
            log.info(f'  {i+1}/{n_eval} done ({time.time()-t0:.0f}s)')

    log.info(f'Heating evaluation done in {time.time()-t0:.0f}s')
    return idx, np.array(k2_results), heating_results


# ============================================================
# Step 6: Plots
# ============================================================
def make_plots(samples, log_likes, k2_results, heating_results, eval_idx):
    """Generate corner plot and heating distribution plots."""
    import corner
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    # --- Corner plot ---
    fig = corner.corner(
        samples, labels=PARAM_LABELS,
        quantiles=[0.16, 0.5, 0.84], show_titles=True,
        title_fmt='.3f', title_kwargs={'fontsize': 10}
    )
    fig.suptitle('Andrade No-Ocean Titan: Posterior', fontsize=14, y=1.02)
    corner_path = os.path.join(OUTPUT_DIR, 'andrade_no_ocean_corner.png')
    fig.savefig(corner_path, dpi=150, bbox_inches='tight')
    log.info(f'Corner plot saved to {corner_path}')
    plt.close(fig)

    # --- k2 scatter with error ellipse ---
    fig, ax = plt.subplots(figsize=(8, 6))
    Re_arr = k2_results[:, 0]
    Im_arr = np.abs(k2_results[:, 1])

    # Collect heating fractions for coloring
    all_phases = set()
    for h in heating_results:
        all_phases.update(h.keys())
    ice_phases = [p for p in all_phases if p in ('Ih', 'II', 'III', 'V', 'VI')]

    f_sil = []
    for h in heating_results:
        total = sum(h.values()) if h else 1e-30
        sil = h.get('Sil', 0)
        f_sil.append(sil / total if total > 0 else 0)
    f_sil = np.array(f_sil)

    sc = ax.scatter(Re_arr, Im_arr, c=f_sil, cmap='RdYlBu_r', s=8, alpha=0.6,
                    vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label='Silicate heating fraction')

    # Error ellipse
    ell = Ellipse((RE_K2_OBS, IM_K2_OBS), 2*RE_K2_ERR, 2*IM_K2_ERR,
                  fill=False, color='red', linewidth=2, linestyle='--',
                  label=f'Petricca 1$\\sigma$')
    ax.add_patch(ell)
    ell2 = Ellipse((RE_K2_OBS, IM_K2_OBS), 4*RE_K2_ERR, 4*IM_K2_ERR,
                   fill=False, color='red', linewidth=1, linestyle=':',
                   label=f'Petricca 2$\\sigma$')
    ax.add_patch(ell2)
    ax.set_xlabel('Re(k$_2$)')
    ax.set_ylabel('|Im(k$_2$)|')
    ax.set_title('Andrade No-Ocean: k$_2$ Posterior')
    ax.legend()
    k2_path = os.path.join(OUTPUT_DIR, 'andrade_no_ocean_k2_scatter.png')
    fig.savefig(k2_path, dpi=150, bbox_inches='tight')
    log.info(f'k2 scatter saved to {k2_path}')
    plt.close(fig)

    # --- Heating distribution vs parameters ---
    eval_samples = samples[eval_idx]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()

    phase_colors = {'Ih': 'C0', 'III': 'C1', 'V': 'C2', 'VI': 'C3', 'Sil': 'C4', 'Clath': 'C5'}

    # Compute heating fractions per phase
    heating_fracs = {}
    for phase in ['Ih', 'III', 'V', 'VI', 'Sil']:
        fracs = []
        for h in heating_results:
            total = sum(h.values()) if h else 1e-30
            fracs.append(h.get(phase, 0) / total if total > 0 else 0)
        heating_fracs[phase] = np.array(fracs)

    for ip, (pname, plabel) in enumerate(zip(PARAM_NAMES, PARAM_LABELS)):
        ax = axes[ip]
        x = eval_samples[:, ip]
        for phase in ['Ih', 'III', 'V', 'VI', 'Sil']:
            if phase in heating_fracs:
                ax.scatter(x, heating_fracs[phase], s=4, alpha=0.3,
                          color=phase_colors.get(phase, 'gray'), label=phase)
        ax.set_xlabel(plabel)
        ax.set_ylabel('Heating fraction')
        ax.set_ylim(-0.05, 1.05)
        if ip == 0:
            ax.legend(fontsize=8, loc='best')

    # Stacked histogram in last panel
    ax = axes[5]
    phases_to_plot = ['Sil', 'V', 'VI', 'III', 'Ih']
    bottom = np.zeros(len(heating_results))
    for phase in phases_to_plot:
        if phase in heating_fracs:
            ax.bar(range(len(heating_results)), heating_fracs[phase], bottom=bottom,
                   color=phase_colors.get(phase, 'gray'), label=phase, width=1.0)
            bottom += heating_fracs[phase]
    ax.set_xlabel('Sample index (sorted by silicate fraction)')
    ax.set_ylabel('Heating fraction')
    ax.legend(fontsize=8)
    ax.set_title('Per-phase heating across posterior')

    fig.suptitle('Andrade No-Ocean: Heating Distribution', fontsize=14)
    fig.tight_layout()
    heat_path = os.path.join(OUTPUT_DIR, 'andrade_no_ocean_heating.png')
    fig.savefig(heat_path, dpi=150, bbox_inches='tight')
    log.info(f'Heating plot saved to {heat_path}')
    plt.close(fig)


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    # Build structure
    data = build_structural_model()

    # Quick sanity check
    theta_test = [0.2, np.log10(0.05), 12.0, 13.0, 19.0]
    Re, Im, hW = forward_model(theta_test, data)
    log.info(f'Sanity check: Re(k2)={Re:.4f}, Im(k2)={Im:.4f}')
    if np.isnan(Re):
        log.error('Forward model returned NaN — check structural model.')
        sys.exit(1)

    # Run MCMC
    samples, log_likes, sampler = run_mcmc(data)

    # Re-evaluate heating on posterior
    eval_idx, k2_results, heating_results = evaluate_heating(samples, data)

    # Save all results for re-analysis
    results = {
        'samples': samples,
        'log_likes': log_likes,
        'param_names': PARAM_NAMES,
        'param_labels': PARAM_LABELS,
        'eval_idx': eval_idx,
        'k2_results': k2_results,
        'heating_results': heating_results,
        'observational_constraints': {
            'Re_k2': (RE_K2_OBS, RE_K2_ERR),
            'Im_k2': (IM_K2_OBS, IM_K2_ERR),
        },
    }
    pkl_path = os.path.join(OUTPUT_DIR, 'andrade_no_ocean_mcmc.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump(results, f)
    log.info(f'Results saved to {pkl_path}')

    # Generate plots
    make_plots(samples, log_likes, k2_results, heating_results, eval_idx)

    # Print summary
    print('\n' + '=' * 60)
    print('ANDRADE NO-OCEAN MCMC SUMMARY')
    print('=' * 60)
    for ip, (pname, plabel) in enumerate(zip(PARAM_NAMES, PARAM_LABELS)):
        med = np.median(samples[:, ip])
        lo = np.percentile(samples[:, ip], 16)
        hi = np.percentile(samples[:, ip], 84)
        print(f'  {pname:15s}: {med:.3f}  [{lo:.3f}, {hi:.3f}]')

    Re_med = np.median(k2_results[:, 0])
    Im_med = np.median(np.abs(k2_results[:, 1]))
    print(f'\n  Re(k2) median: {Re_med:.4f}  (obs: {RE_K2_OBS} +/- {RE_K2_ERR})')
    print(f'  |Im(k2)| median: {Im_med:.4f}  (obs: {IM_K2_OBS} +/- {IM_K2_ERR})')
    print(f'\n  Output directory: {OUTPUT_DIR}')
