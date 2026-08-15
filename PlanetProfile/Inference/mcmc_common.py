"""
Body-agnostic MCMC helpers for PlanetProfile Bayesian inference.

Extracted from Test48_mcmc_andrade_yao2014.py to allow reuse across bodies
(Titan, Europa, Ganymede, ...) with consistent forward-model physics.

Public API:
    apply_arrhenius_ih        — T-dependent Ih viscosity (Yao convention)
    split_silicate_core       — recompute silicate/core density + moduli + layering
    build_andrade_shear_bulk  — per-layer Andrade shear + Elastic bulk models
    build_maxwell_shear_bulk  — per-layer Maxwell shear + Elastic bulk models
    compute_per_phase_heating — integrate radial dissipation by phase
    run_pocomc_sampler        — pocoMC wrapper with timing + logging
    evaluate_posterior        — re-evaluate a forward model on N posterior samples
    gaussian_chi2_terms       — helper to build Σ ((x - μ)/σ)² from OBS dict
"""
from __future__ import annotations

import logging
import time
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

log = logging.getLogger(__name__)


# ============================================================
# Arrhenius viscosity scaling (Ice Ih, Yao-style basal reference)
# ============================================================

def apply_arrhenius_ih(
    eta_mod: np.ndarray,
    phases: np.ndarray,
    ci: Sequence[int],
    n_layers: int,
    T_K_profile: Optional[np.ndarray],
    Tb_K: float,
    E_act_J_mol: float = 60e3,
    R_gas: float = 8.314462,
    ih_phase_id: int = 1,
) -> np.ndarray:
    """Apply Arrhenius T-dependence to Ice Ih layers in place.

    Sampled eta_Ih (already assigned uniformly to Ih layers) is interpreted as
    the BASAL viscosity at T=Tb (Yao 2014 convention).  Each node's viscosity
    is then scaled by exp(E/R · (1/T - 1/Tb)), so the cold stagnant lid
    becomes much stiffer and the warm convecting interior stays near the
    sampled value.

    Args:
        eta_mod: Viscosity array (modified in place).
        phases: Phase ID array, same length as eta_mod.
        ci: Layer boundary indices.
        n_layers: Number of layers (len(ci) - 1).
        T_K_profile: Per-node temperature profile.  If None or shape mismatch,
            the function is a no-op (returns eta_mod unchanged).
        Tb_K: Ice Ih basal temperature.
        E_act_J_mol: Activation energy (default 60 kJ/mol — Yao Table 5).
        R_gas: Gas constant.
        ih_phase_id: Phase ID for Ice Ih (default 1).

    Returns:
        eta_mod (same array, modified in place).
    """
    if T_K_profile is None:
        return eta_mod
    T_arr = np.asarray(T_K_profile)
    if T_arr.shape != eta_mod.shape:
        return eta_mod
    for i in range(n_layers):
        s, e = ci[i], ci[i + 1]
        if int(phases[min(s, len(phases) - 1)]) != ih_phase_id:
            continue
        T_layer = T_arr[s:e]
        if not (np.all(np.isfinite(T_layer)) and np.all(T_layer > 0)):
            continue
        exponent = (E_act_J_mol / R_gas) * (1.0 / T_layer - 1.0 / Tb_K)
        eta_mod[s:e] *= np.exp(exponent)
    return eta_mod


# ============================================================
# Silicate/core density + moduli + layering recomputation
# ============================================================

def split_silicate_core(
    r_m: np.ndarray,
    rho_mod: np.ndarray,
    K_Pa_mod: np.ndarray,
    mu_Pa_mod: np.ndarray,
    phases: np.ndarray,
    ci: List[int],
    n_layers: int,
    layer_upper_radii: List[float],
    layer_types: List[str],
    region_phases: List[str],
    rho_sil: float,
    rho_core: float,
    f_core: float,
    sil_phase_lo: int = 50,
    sil_phase_hi: int = 100,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[int], int,
           List[float], List[str], List[str]]:
    """Split the innermost silicate layer at R_core and apply two-layer density.

    Density/moduli update via K ∝ ρ, μ ∝ ρ at constant Vp, Vs.  If f_core < ~0,
    falls back to uniform ρ_sil through all silicate layers.

    Modifies rho_mod/K_Pa_mod/mu_Pa_mod in place and returns updated arrays
    plus the possibly-extended layering (ci, n_layers, layer_upper_radii,
    layer_types, region_phases).

    Args:
        r_m: Radial grid (m).
        rho_mod, K_Pa_mod, mu_Pa_mod: Arrays to modify (in place).
        phases: Phase IDs.
        ci, n_layers, layer_upper_radii, layer_types, region_phases:
            Layer metadata.  ci and the lists are NOT modified; copies with
            the new inner-core layer inserted are returned.
        rho_sil: Outer silicate density (kg/m³).
        rho_core: Inner core density (kg/m³).  Must be ≥ rho_sil for stability.
        f_core: R_core / r_sil_top in [0, 0.8].
        sil_phase_lo, sil_phase_hi: Phase ID range treated as silicate.

    Returns:
        (rho_mod, K_Pa_mod, mu_Pa_mod, ci_out, n_layers_out,
         layer_upper_radii_out, layer_types_out, region_phases_out)
    """
    # Local copies of lists that may be extended
    ci_out = list(ci)
    layer_upper_radii_out = list(layer_upper_radii)
    layer_types_out = list(layer_types)
    region_phases_out = list(region_phases)

    sil_layers = []
    for i in range(n_layers):
        s, e = ci_out[i], ci_out[i + 1]
        if sil_phase_lo <= int(phases[s]) < sil_phase_hi:
            sil_layers.append(i)

    if sil_layers and f_core > 0.001:
        # Split innermost silicate layer at R_core
        i_sil = sil_layers[0]
        s_sil, e_sil = ci_out[i_sil], ci_out[i_sil + 1]
        r_sil_layer = r_m[s_sil:e_sil]

        r_sil_top = layer_upper_radii_out[i_sil] if i_sil < len(layer_upper_radii_out) else r_m[e_sil - 1]
        R_core = f_core * r_sil_top

        idx_split = int(np.searchsorted(r_sil_layer, R_core))
        idx_split = max(2, min(idx_split, len(r_sil_layer) - 2))
        abs_split = s_sil + idx_split

        rho_cached = rho_mod[s_sil:e_sil].copy()
        rho_mod[s_sil:abs_split] = rho_core
        rho_mod[abs_split:e_sil] = rho_sil
        scale_core = np.where(rho_cached[:idx_split] > 0,
                              rho_core / rho_cached[:idx_split], 1.0)
        scale_mantle = np.where(rho_cached[idx_split:] > 0,
                                rho_sil / rho_cached[idx_split:], 1.0)
        K_Pa_mod[s_sil:abs_split] *= scale_core
        mu_Pa_mod[s_sil:abs_split] *= scale_core
        K_Pa_mod[abs_split:e_sil] *= scale_mantle
        mu_Pa_mod[abs_split:e_sil] *= scale_mantle

        # Additional silicate layers above the innermost use rho_sil
        for i_extra in sil_layers[1:]:
            s_e, e_e = ci_out[i_extra], ci_out[i_extra + 1]
            rho_cached_e = rho_mod[s_e:e_e].copy()
            rho_mod[s_e:e_e] = rho_sil
            scale_e = np.where(rho_cached_e > 0, rho_sil / rho_cached_e, 1.0)
            K_Pa_mod[s_e:e_e] *= scale_e
            mu_Pa_mod[s_e:e_e] *= scale_e

        # Insert boundary at the split point
        R_core_actual = float(r_m[abs_split])
        new_ci = ci_out[:i_sil + 1] + [abs_split] + ci_out[i_sil + 1:]
        layer_upper_radii_out.insert(i_sil, R_core_actual)
        layer_types_out.insert(i_sil, 'solid')
        region_phases_out.insert(i_sil, 'Sil')
        n_layers += 1
        ci_out = new_ci
    else:
        # Uniform ρ_sil across all silicate layers
        for i_s in sil_layers:
            s_s, e_s = ci_out[i_s], ci_out[i_s + 1]
            rho_cached_s = rho_mod[s_s:e_s].copy()
            rho_mod[s_s:e_s] = rho_sil
            scale_s = np.where(rho_cached_s > 0, rho_sil / rho_cached_s, 1.0)
            K_Pa_mod[s_s:e_s] *= scale_s
            mu_Pa_mod[s_s:e_s] *= scale_s

    return (rho_mod, K_Pa_mod, mu_Pa_mod,
            ci_out, n_layers, layer_upper_radii_out, layer_types_out, region_phases_out)


# ============================================================
# Rheology model builders (per-layer shear + bulk)
# ============================================================

def build_andrade_shear_bulk(
    region_phases: Sequence[str],
    alpha: float,
    log10_zeta: float,
    elastic_region_names: Tuple[str, ...] = ('0', 'Clath'),
) -> Tuple[List[Any], List[Any]]:
    """Per-layer shear (Andrade / Elastic) + bulk (Elastic) module lists.

    zeta_tp = (10**log10_zeta)**(1/alpha) is used as Andrade's τ-like parameter.

    Returns (shear_list, bulk_list) with the same length as region_phases.
    """
    from TidalPy.rheology import Andrade, Elastic  # local import

    zeta_tp = (10 ** log10_zeta) ** (1.0 / alpha)
    shear: List[Any] = []
    for rp in region_phases:
        base = rp.replace('_conv', '')
        if base in elastic_region_names:
            shear.append(Elastic())
        else:
            shear.append(Andrade(args=(alpha, zeta_tp)))
    bulk: List[Any] = [Elastic() for _ in shear]
    return shear, bulk


def build_maxwell_shear_bulk(
    region_phases: Sequence[str],
    elastic_region_names: Tuple[str, ...] = ('0', 'Clath'),
) -> Tuple[List[Any], List[Any]]:
    """Per-layer shear (Maxwell / Elastic) + bulk (Elastic) module lists."""
    from TidalPy.rheology import Maxwell, Elastic

    shear: List[Any] = []
    for rp in region_phases:
        base = rp.replace('_conv', '')
        if base in elastic_region_names:
            shear.append(Elastic())
        else:
            shear.append(Maxwell())
    bulk: List[Any] = [Elastic() for _ in shear]
    return shear, bulk


# ============================================================
# Per-phase heating integration (TidalPy radial solution → phase powers)
# ============================================================

def compute_per_phase_heating(
    result: Any,
    data: Dict[str, Any],
    phase_map: Optional[Dict[int, str]] = None,
) -> Dict[str, float]:
    """Integrate TidalPy volumetric dissipation over each layer by phase.

    Args:
        result: TidalPy radial_solver result (needs .radius_array).
        data: Structure dict with 'r_m', 'changeIndices', 'n_layers',
              'phases', 'eccentricity', 'omega', 'a_m', 'host_mass'.
        phase_map: Optional override for numeric→name mapping; defaults to
            data.get('phase_map', {0:'0', 1:'Ih', 2:'II', 3:'III', 5:'V', 6:'VI'}).

    Returns:
        Dict mapping phase name → total heating power (W).  Silicate layers
        (phase 50–99) use 'Sil'; core (phase ≥ 100) uses 'Core'.
    """
    from TidalPy.tides.multilayer.heating import (
        calc_radial_volumetric_tidal_heating_from_rs_solution,
    )
    from PlanetProfile.Utilities.Indexing import PhaseConv

    if phase_map is None:
        phase_map = data.get('phase_map',
                             {0: '0', 1: 'Ih', 2: 'II', 3: 'III', 5: 'V', 6: 'VI'})

    hp = calc_radial_volumetric_tidal_heating_from_rs_solution(
        data['eccentricity'], data['omega'], data['a_m'],
        data['host_mass'], result, perform_checks=False,
    )
    rr = np.asarray(result.radius_array)
    r_m = data['r_m']
    ci = data['changeIndices']
    n = data['n_layers']
    phases = data['phases']

    out: Dict[str, float] = {}
    for i in range(n):
        s, e = ci[i], ci[i + 1]
        ph = int(phases[s])
        if 50 <= ph < 100:
            pname = 'Sil'
        elif ph >= 100:
            pname = 'Core'
        else:
            pname = phase_map.get(ph, PhaseConv(ph, liq='0'))
        lo, hi = r_m[s], r_m[e - 1]
        mask = (rr >= lo - 1.0) & (rr <= hi + 1.0)
        if np.any(mask):
            lr, lh = rr[mask], hp[mask]
            if len(lr) > 1:
                power = float(np.trapezoid(lh * 4 * np.pi * lr**2, lr))
            else:
                power = float(lh[0] * (4.0 / 3.0) * np.pi * (hi**3 - lo**3))
            out[pname] = out.get(pname, 0.0) + power
    return out


# ============================================================
# pocoMC sampler wrapper
# ============================================================

def run_pocomc_sampler(
    prior: Any,
    log_likelihood_fn: Callable[[Sequence[float]], float],
    n_effective: int = 500,
    random_state: int = 42,
    return_sampler: bool = True,
) -> Tuple[np.ndarray, np.ndarray, Optional[Any]]:
    """Run pocoMC with standard settings and light logging.

    Returns (samples, log_likes, sampler).  sampler may be useful for
    evidence computation (sampler.logz).
    """
    import pocomc as pc  # local import keeps top-level light

    log.info(f'Starting pocoMC ({prior.dim}D, n_eff={n_effective})')
    t0 = time.time()
    sampler = pc.Sampler(
        prior=prior,
        likelihood=log_likelihood_fn,
        n_effective=n_effective,
        random_state=random_state,
    )
    sampler.run()
    elapsed = time.time() - t0
    samples, log_likes, _, _ = sampler.posterior()
    log.info(f'pocoMC done in {elapsed / 60:.1f} min, {len(samples)} samples '
             f'(best log-like: {np.max(log_likes):.2f})')
    return samples, log_likes, (sampler if return_sampler else None)


# ============================================================
# Posterior re-evaluation (slow: includes heating)
# ============================================================

def evaluate_posterior(
    samples: np.ndarray,
    forward_fn: Callable[[Sequence[float]], Tuple[Any, ...]],
    n_eval: Optional[int] = None,
    random_state: int = 42,
    log_every: int = 100,
) -> Tuple[np.ndarray, List[Tuple[Any, ...]]]:
    """Re-evaluate forward_fn on a subset of posterior samples.

    Returns (eval_idx, results_list) where results_list[i] is whatever
    forward_fn returned for samples[eval_idx[i]].  The caller is responsible
    for unpacking (e.g. (Re_k2, Im_k2, Mtot, CMR2, heating)).
    """
    n_total = len(samples)
    if n_eval is None:
        n_eval = min(500, n_total)
    n_eval = min(n_eval, n_total)

    rng = np.random.default_rng(random_state)
    idx = rng.choice(n_total, n_eval, replace=False)
    idx.sort()

    log.info(f'Re-evaluating {n_eval} posterior samples...')
    t0 = time.time()
    results: List[Tuple[Any, ...]] = []
    for i, si in enumerate(idx):
        results.append(forward_fn(samples[si]))
        if (i + 1) % log_every == 0:
            log.info(f'  {i + 1}/{n_eval} ({time.time() - t0:.0f}s)')
    log.info(f'Done in {time.time() - t0:.0f}s')
    return idx, results


# ============================================================
# Gaussian chi² builder
# ============================================================

def gaussian_chi2_terms(
    observed: Dict[str, Tuple[float, float]],
    predicted: Dict[str, float],
) -> float:
    """Accumulate Σ ((pred - obs)/sigma)² over matched observables.

    Args:
        observed: {name: (value, sigma)}.
        predicted: {name: value}.  Keys that are NaN or absent are skipped.

    Returns:
        Sum of squared normalized residuals.
    """
    chi2 = 0.0
    for name, (mu, sigma) in observed.items():
        val = predicted.get(name)
        if val is None or not np.isfinite(val):
            continue
        chi2 += ((val - mu) / sigma) ** 2
    return chi2
