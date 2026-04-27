"""
PPTest40 validation: Andrade vs Maxwell on ocean-bearing Titan.

Runs four configurations:
  A. Andrade + TidalPy (reference)
  B. Maxwell + TidalPy (default viscosities)
  C. Maxwell + TidalPy (tuned: eta_sil = 5.6e14 for k2 match)
  D. Maxwell + TidalPy (iterated: re-run C with TidalPy-derived Sil.Htidal_Wm3)

Compares Re(k2), Im(k2), Q, per-phase heating distribution.
Tests whether self-consistent per-phase heating redistribution works.
"""
import numpy as np
import importlib
import copy
import os
import sys
import platformdirs

# Set up paths relative to this script (PlanetProfile/Test/)
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_THIS_DIR, '..', '..'))
platformdirs.user_documents_dir = lambda: _REPO_ROOT
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import logging
logging.basicConfig(level=logging.INFO, format='%(name)s - %(message)s')

from PlanetProfile.Utilities.defineStructs import EOSlist, Constants
from PlanetProfile.Main import PlanetProfile as RunPP

MAXWELL_RHEOLOGY = {
    '0': 'newton', 'Ih': 'maxwell', 'Ih_conv': 'maxwell',
    'II': 'maxwell', 'III': 'maxwell', 'III_conv': 'maxwell',
    'IV': 'maxwell', 'V': 'maxwell', 'V_conv': 'maxwell', 'VI': 'maxwell',
    'Sil': 'maxwell', 'Fe': 'elastic', 'Clath': 'elastic', 'Clath_conv': 'maxwell'
}

results = {}

configs = [
    {
        'label': 'A. Andrade (reference)',
        'rheology_override': None,
        'eta_sil_override': None,
        'htidal_sil_override': None,
    },
    {
        'label': 'B. Maxwell (default eta)',
        'rheology_override': MAXWELL_RHEOLOGY,
        'eta_sil_override': None,
        'htidal_sil_override': None,
    },
    {
        'label': 'C. Maxwell (eta_sil=5.6e14)',
        'rheology_override': MAXWELL_RHEOLOGY,
        'eta_sil_override': 5.6e14,
        'htidal_sil_override': None,
    },
]

# Run configs A, B, C first; D depends on C's output
for cfg in configs:
    label = cfg['label']
    print("\n" + "=" * 70)
    print(f"  {label}")
    print("=" * 70)

    EOSlist.loaded.clear()
    from PlanetProfile.GetConfig import Params as configParams
    test_module = 'PlanetProfile.Test.PPTest40'
    importlib.reload(sys.modules[test_module] if test_module in sys.modules
                     else importlib.import_module(test_module))
    mod = sys.modules[test_module]
    Planet = mod.Planet

    configParams.Gravity.backend = 'tidalpy'
    configParams.CALC_NEW = True
    configParams.CALC_NEW_GRAVITY = True
    configParams.NO_SAVEFILE = True
    configParams.SKIP_PLOTS = True

    if cfg['rheology_override'] is not None:
        configParams.Gravity.rheology_models = cfg['rheology_override']

    if cfg['eta_sil_override'] is not None:
        Planet.Sil.etaRock_Pas = [cfg['eta_sil_override'], cfg['eta_sil_override']]

    if cfg['htidal_sil_override'] is not None:
        Planet.Sil.Htidal_Wm3 = cfg['htidal_sil_override']

    Planet, Params = RunPP(Planet, configParams)

    k2 = Planet.Gravity.k
    if k2 is not None and np.isfinite(k2):
        Im_k2 = abs(np.imag(k2))
        Re_k2 = np.real(k2)
        E_dot = (21./2.) * Planet.Bulk.meanMotion_radps**5 * Planet.Bulk.R_m**5 * Planet.Bulk.eccentricity**2 * Im_k2 / Constants.G
        Q = abs(Re_k2 / np.imag(k2)) if np.imag(k2) != 0 else float('inf')
    else:
        Im_k2 = Re_k2 = E_dot = 0
        Q = float('inf')

    perPhase_W = getattr(Planet.Gravity, 'tidalpy_Htidal_perPhase_W', None)
    perPhase_Wm3 = getattr(Planet.Gravity, 'tidalpy_Htidal_perPhase_Wm3', None)

    results[label] = {
        'k2': k2, 'Re_k2': Re_k2, 'Im_k2': Im_k2,
        'E_dot': E_dot, 'Q': Q,
        'perPhase_W': perPhase_W, 'perPhase_Wm3': perPhase_Wm3,
        'HtidalIce_Wm3': Planet.Ocean.HtidalIce_Wm3,
        'Htidal_sil_Wm3': Planet.Sil.Htidal_Wm3,
    }

    print(f"\n  k2 = {k2}")
    print(f"  Re(k2) = {Re_k2:.6f}, Im(k2) = {-Im_k2:.6f}, Q = {Q:.1f}")
    print(f"  E_dot = {E_dot/1e9:.1f} GW")
    print(f"  HtidalIce_Wm3 = {Planet.Ocean.HtidalIce_Wm3:.2e}")
    print(f"  Sil.Htidal_Wm3 = {Planet.Sil.Htidal_Wm3:.2e}")
    if perPhase_W:
        total = sum(perPhase_W.values())
        print(f"  Per-layer heating (total = {total/1e9:.1f} GW):")
        for ph in sorted(perPhase_W.keys()):
            pct = 100 * perPhase_W[ph] / total if total > 0 else 0
            print(f"    {ph:6s}: {perPhase_W[ph]/1e9:8.1f} GW  "
                  f"({perPhase_Wm3[ph]:.3e} W/m^3)  {pct:5.1f}%")


# Config D: re-run C with TidalPy-derived silicate heating
label_C = 'C. Maxwell (eta_sil=5.6e14)'
if label_C in results and results[label_C]['Htidal_sil_Wm3'] > 0:
    label = 'D. Maxwell (iterated Sil heating)'
    print("\n" + "=" * 70)
    print(f"  {label}")
    print(f"  Using Sil.Htidal_Wm3 = {results[label_C]['Htidal_sil_Wm3']:.2e} from config C")
    print("=" * 70)

    EOSlist.loaded.clear()
    from PlanetProfile.GetConfig import Params as configParams
    importlib.reload(sys.modules[test_module])
    mod = sys.modules[test_module]
    Planet = mod.Planet

    configParams.Gravity.backend = 'tidalpy'
    configParams.CALC_NEW = True
    configParams.CALC_NEW_GRAVITY = True
    configParams.NO_SAVEFILE = True
    configParams.SKIP_PLOTS = True
    configParams.Gravity.rheology_models = MAXWELL_RHEOLOGY

    Planet.Sil.etaRock_Pas = [5.6e14, 5.6e14]
    # Feed back TidalPy's silicate heating prediction from config C
    Planet.Sil.Htidal_Wm3 = results[label_C]['Htidal_sil_Wm3']
    # Also feed back ice heating
    Planet.Ocean.HtidalIce_Wm3 = results[label_C]['HtidalIce_Wm3']

    Planet, Params = RunPP(Planet, configParams)

    k2 = Planet.Gravity.k
    if k2 is not None and np.isfinite(k2):
        Im_k2 = abs(np.imag(k2))
        Re_k2 = np.real(k2)
        E_dot = (21./2.) * Planet.Bulk.meanMotion_radps**5 * Planet.Bulk.R_m**5 * Planet.Bulk.eccentricity**2 * Im_k2 / Constants.G
        Q = abs(Re_k2 / np.imag(k2)) if np.imag(k2) != 0 else float('inf')
    else:
        Im_k2 = Re_k2 = E_dot = 0
        Q = float('inf')

    perPhase_W = getattr(Planet.Gravity, 'tidalpy_Htidal_perPhase_W', None)
    perPhase_Wm3 = getattr(Planet.Gravity, 'tidalpy_Htidal_perPhase_Wm3', None)

    results[label] = {
        'k2': k2, 'Re_k2': Re_k2, 'Im_k2': Im_k2,
        'E_dot': E_dot, 'Q': Q,
        'perPhase_W': perPhase_W, 'perPhase_Wm3': perPhase_Wm3,
        'HtidalIce_Wm3': Planet.Ocean.HtidalIce_Wm3,
        'Htidal_sil_Wm3': Planet.Sil.Htidal_Wm3,
    }

    print(f"\n  k2 = {k2}")
    print(f"  Re(k2) = {Re_k2:.6f}, Im(k2) = {-Im_k2:.6f}, Q = {Q:.1f}")
    print(f"  E_dot = {E_dot/1e9:.1f} GW")
    print(f"  HtidalIce_Wm3 = {Planet.Ocean.HtidalIce_Wm3:.2e}")
    print(f"  Sil.Htidal_Wm3 = {Planet.Sil.Htidal_Wm3:.2e}")
    if perPhase_W:
        total = sum(perPhase_W.values())
        print(f"  Per-layer heating (total = {total/1e9:.1f} GW):")
        for ph in sorted(perPhase_W.keys()):
            pct = 100 * perPhase_W[ph] / total if total > 0 else 0
            print(f"    {ph:6s}: {perPhase_W[ph]/1e9:8.1f} GW  "
                  f"({perPhase_Wm3[ph]:.3e} W/m^3)  {pct:5.1f}%")


# ============================================================
# Summary comparison
# ============================================================
print("\n\n" + "=" * 70)
print("COMPARISON SUMMARY")
print("=" * 70)

print(f"\nPetricca et al. (2025):  Re(k2) ~ 0.6, Im(k2) ~ -0.12, Q ~ 5")
print(f"No-ocean Andrade ref:   Re(k2) = 0.578, Im(k2) = -0.123, Q = 4.7\n")

print(f"{'Configuration':<35} {'Re(k2)':>8} {'Im(k2)':>8} {'Q':>6} {'E_dot':>8} {'Ice GW':>8} {'Sil GW':>8}")
print("-" * 90)
for label in results:
    r = results[label]
    ice_W = 0
    sil_W = 0
    if r['perPhase_W']:
        for p in ['Ih', 'II', 'III', 'V', 'VI']:
            ice_W += r['perPhase_W'].get(p, 0)
        sil_W = r['perPhase_W'].get('Sil', 0)
    print(f"{label:<35} {r['Re_k2']:8.4f} {-r['Im_k2']:8.4f} {r['Q']:6.1f} "
          f"{r['E_dot']/1e9:8.1f} {ice_W/1e9:8.1f} {sil_W/1e9:8.1f}")

# Heating distribution comparison
print("\n--- Per-phase heating distribution ---")
all_phases = set()
for r in results.values():
    if r['perPhase_W']:
        all_phases.update(r['perPhase_W'].keys())
all_phases = sorted(all_phases)

header = f"{'Phase':<8}"
for label in results:
    short = label.split('.')[0] + '.'
    header += f" {short:>10}"
print(header)
print("-" * (8 + 11 * len(results)))

for ph in all_phases:
    row = f"{ph:<8}"
    for label in results:
        r = results[label]
        if r['perPhase_W'] and ph in r['perPhase_W']:
            row += f" {r['perPhase_W'][ph]/1e9:10.1f}"
        else:
            row += f" {'--':>10}"
    print(row)

# Total
row = f"{'TOTAL':<8}"
for label in results:
    r = results[label]
    if r['perPhase_W']:
        total = sum(r['perPhase_W'].values())
        row += f" {total/1e9:10.1f}"
    else:
        row += f" {'--':>10}"
print(row)

# Self-consistent heating values
print("\n--- Self-consistent heating values after run ---")
print(f"{'Configuration':<35} {'HtidalIce':>12} {'Htidal_sil':>12}")
print("-" * 60)
for label in results:
    r = results[label]
    print(f"{label:<35} {r['HtidalIce_Wm3']:12.2e} {r['Htidal_sil_Wm3']:12.2e}")

print("\n--- Key findings ---")
print("  1. With ocean, tidal response differs from no-ocean model.")
print("  2. Under Andrade, HP ices dominate dissipation (transient creep).")
print("  3. Under Maxwell (default eta), HP ices are sub-peak → low dissipation.")
print("  4. Under Maxwell (tuned eta_sil~5e14), silicate mantle dominates.")
print("  5. Self-consistent iteration (D) feeds TidalPy heating back to layers.")
