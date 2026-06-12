"""Ad-hoc diagnostic: inspect post-fix Test50 grid cache to verify T(r), P(r),
and implied heat flux q(r) through clathrate + ice Ih conductive lid.

Used in May-2026 session to validate the ConductiveTemperature c1 fix
(T&S 4.40 restored — factor-of-2 removed).  Expected outcomes:
  - Clathrate layer thickness now matches Bulk.clathMaxThick_m input (2 km)
  - Local q through clathrate ≈ Carnahan-matched qTop (~3 mW/m² at current Tb)
  - q approximately constant through the conductive lid (Fourier steady state)

Run from repo root:
    mamba activate PPcl
    python PlanetProfile/Test/scripts/probe_heat_flux.py
"""
import pickle
import sys
import os
import numpy as np

_pp_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, _pp_root)

cache_path = os.path.join(_pp_root, 'PlanetProfile', 'Test', 'mcmc_results',
                          'Titan', 'Test50_andrade_noocean_yao2014',
                          'titan_allice_yao2014_structure_grid.pkl')
with open(cache_path, 'rb') as f:
    grid = pickle.load(f)

# Upper-edge structure (Tb = 250.965 K)
s = grid['structures'][-1]

phase = np.asarray(s['phases']).astype(int)
P = np.asarray(s['P_MPa'])
T = np.asarray(s['T_K'])
r = np.asarray(s['r_m'])

# ALMA arrays are bottom-up; reverse to top-down for readable print order
idx = np.argsort(-r)
phase, P, T, r = phase[idx], P[idx], T[idx], r[idx]

print(f"Tb grid: {grid['Tb_K_grid']}")
print(f"Inspecting Tb = {s['Tb_K']:.4f} K  (grid upper edge)")
print(f"Cells: {len(phase)}")

clath_m = (phase >= 30) & (phase < 40)
Ih_m = (phase == 1)
III_m = (phase == 3)
V_m = (phase == 5)
VI_m = (phase == 6)
sil_m = (phase >= 50) & (phase < 100)

def summarize(name, m):
    if m.sum() == 0:
        return
    print(f'  {name:6s}: N={m.sum():4d}  T[K]=[{T[m].min():6.2f}, {T[m].max():6.2f}]  '
          f'P[MPa]=[{P[m].min():7.2f}, {P[m].max():7.2f}]  '
          f'r[km]=[{r[m].min()/1e3:7.1f}, {r[m].max()/1e3:7.1f}]')

print("\nPer-phase summary:")
for nm, m in [('Clath', clath_m), ('Ih', Ih_m), ('III', III_m), ('V', V_m), ('VI', VI_m), ('Sil', sil_m)]:
    summarize(nm, m)

# Local q = -k·dT/dr with approximate k(T, phase)
k = np.full_like(T, np.nan)
k[Ih_m] = 651.0 / T[Ih_m]           # Hobbs 1974 for Ice Ih
k[clath_m] = 0.5                     # Clathrate near-constant
for m in (III_m, V_m, VI_m):
    k[m] = 0.47 + 0.15 * (T[m] - 250.0) / 50.0  # Andersson & Inaba 2005 rough linearization

dT_dr = np.zeros_like(T)
dT_dr[1:-1] = (T[2:] - T[:-2]) / (r[2:] - r[:-2])
dT_dr[0] = (T[1] - T[0]) / (r[1] - r[0])
dT_dr[-1] = (T[-1] - T[-2]) / (r[-1] - r[-2])
q_Wm2 = -k * dT_dr

def print_cells(label, indices):
    print(f"\n{label}:")
    for i in indices:
        print(f'  cell {i:3d}: r={r[i]/1e3:8.2f} km, P={P[i]:7.2f} MPa, T={T[i]:7.2f} K, '
              f'k={k[i]:.3f} W/m/K, dT/dr={dT_dr[i]*1e3:.3f} K/km, q={q_Wm2[i]*1e3:.2f} mW/m²')

if clath_m.sum() > 0:
    clath_idx = np.where(clath_m)[0]
    print_cells('Clathrate — top 5 cells', clath_idx[:5])
    print_cells('Clathrate — base (last 3 cells)', clath_idx[-3:])

if Ih_m.sum() > 0:
    Ih_idx = np.where(Ih_m)[0]
    print_cells('Ice Ih — top 3 (just below clathrate)', Ih_idx[:3])
    print_cells('Ice Ih — deepest 3 (approaching TBL)', Ih_idx[-3:])

# Top-cell Fourier estimate for surface heat flux
if clath_m.sum() > 1:
    i0, i1 = clath_idx[0], clath_idx[1]
    q_surf_est = k[i0] * (T[i1] - T[i0]) / (r[i0] - r[i1])
    print(f'\nSurface q estimate (top 2 clathrate cells, Fourier): {q_surf_est*1e3:.3f} mW/m²')

# Clathrate base T relative to analytic expectation
if clath_m.sum() > 0 and Ih_m.sum() > 0:
    T_clath_base = T[clath_idx[-1]]
    H_clath_m = r[clath_idx[0]] - r[clath_idx[-1]]
    print(f'\nClathrate layer diagnostics:')
    print(f'  measured thickness = {H_clath_m/1e3:.3f} km')
    print(f'  T at clathrate base = {T_clath_base:.2f} K')
    print(f'  ΔT across clathrate  = {T_clath_base - T[clath_idx[0]]:.2f} K')
