"""Verify silicate shear/bulk moduli (and viscosity base) in the post-fix cache
are not contaminated by the diverging silicate T profile.

For Test50 with Do.CONSTANT_INNER_DENSITY=True + Do.Fe_CORE=False, silicate rho
is fixed at Sil.rhoNoCore_kgm3.  But K, G come from the Perple_X silicate EOS
queried at (P, T), so if T -> 50,000+ K they may extrapolate wildly.

Forward-model decisions depend on this:
  - If mu_Pa and K_Pa look reasonable O(10^10 Pa), Test50 MCMC is safe.
  - If either is NaN or huge, we need mitigation before running MCMC.

Also prints eta_Pa_base for silicate (overwritten by sampled eta_sil at MCMC
runtime, so only informational).
"""
import os, sys, pickle
import numpy as np

_pp_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
cache_path = os.path.join(_pp_root, 'PlanetProfile', 'Test', 'mcmc_results',
                          'titan_allice_yao2014_structure_grid.pkl')
with open(cache_path, 'rb') as f:
    grid = pickle.load(f)

s = grid['structures'][-1]  # upper-edge Tb
phases = np.asarray(s['phases']).astype(int)
r = np.asarray(s['r_m'])
rho = np.asarray(s['rho'])
mu = np.asarray(s['mu_Pa'])
K = np.asarray(s['K_Pa'])
eta_base = np.asarray(s['eta_Pa_base'])
T = np.asarray(s['T_K'])
P = np.asarray(s['P_MPa'])

# Sort top-down for clarity
idx = np.argsort(-r)
phases, r, rho, mu, K, eta_base, T, P = [a[idx] for a in (phases, r, rho, mu, K, eta_base, T, P)]

sil_m = (phases >= 50) & (phases < 100)
Ih_m = (phases == 1)

print(f"Cache: {cache_path}")
print(f"Inspecting Tb={s['Tb_K']:.4f} K structure")
print(f"Silicate cells: {sil_m.sum()}")
print()

# Silicate statistics
if sil_m.sum() > 0:
    si = np.where(sil_m)[0]
    print(f"Silicate range:")
    print(f"  r   = [{r[si].min()/1e3:.1f}, {r[si].max()/1e3:.1f}] km")
    print(f"  P   = [{P[si].min():.2f}, {P[si].max():.2f}] MPa")
    print(f"  T   = [{T[si].min():.2f}, {T[si].max():.2f}] K")
    print(f"  rho = [{rho[si].min():.2f}, {rho[si].max():.2f}] kg/m^3")
    print(f"  mu  = [{mu[si].min():.3e}, {mu[si].max():.3e}] Pa")
    print(f"  K   = [{K[si].min():.3e}, {K[si].max():.3e}] Pa")
    print(f"  eta = [{eta_base[si].min():.3e}, {eta_base[si].max():.3e}] Pa s (base; overwritten at MCMC)")

    n_nan_mu = int(np.sum(~np.isfinite(mu[si])))
    n_nan_K  = int(np.sum(~np.isfinite(K[si])))
    n_neg_mu = int(np.sum(mu[si] <= 0))
    n_neg_K  = int(np.sum(K[si] <= 0))
    print()
    print(f"Silicate NaN/non-positive counts:")
    print(f"  mu:  NaN={n_nan_mu}, non-pos={n_neg_mu}")
    print(f"  K:   NaN={n_nan_K}, non-pos={n_neg_K}")
    print()

    # Sample a few silicate cells spread across the stack
    print("Sample silicate cells (every 10th, top-down order):")
    print(f"  {'cell':>4}  {'r[km]':>8}  {'P[MPa]':>8}  {'T[K]':>10}  {'rho':>8}  {'mu[Pa]':>11}  {'K[Pa]':>11}")
    for i in si[::10]:
        print(f"  {i:>4d}  {r[i]/1e3:>8.1f}  {P[i]:>8.2f}  {T[i]:>10.1f}  {rho[i]:>8.1f}  "
              f"{mu[i]:>11.3e}  {K[i]:>11.3e}")

    # Sanity bracket: typical mantle silicate mu = 60-100 GPa = 6e10 - 1e11 Pa; K = 100-300 GPa
    mu_reasonable = ((1e10 <= mu[si]) & (mu[si] <= 5e11))
    K_reasonable = ((5e10 <= K[si]) & (K[si] <= 1e12))
    print(f"\nFraction of silicate cells with 'reasonable' moduli:")
    print(f"  mu in [1e10, 5e11] Pa: {mu_reasonable.sum() / sil_m.sum() * 100:.1f}%")
    print(f"  K  in [5e10, 1e12] Pa: {K_reasonable.sum() / sil_m.sum() * 100:.1f}%")

# Ice Ih reference for comparison
if Ih_m.sum() > 0:
    ii = np.where(Ih_m)[0]
    print(f"\nIce Ih reference (for comparison):")
    print(f"  mu = [{mu[ii].min():.3e}, {mu[ii].max():.3e}] Pa")
    print(f"  K  = [{K[ii].min():.3e}, {K[ii].max():.3e}] Pa")
