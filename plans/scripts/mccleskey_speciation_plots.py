"""Documentation plots for the McCleskey equilibrium-speciation change
(2026-07-16): for NaCl and MgSO4 oceans, contours over (T, molality) of

  1. conductivity WITHOUT speciation (nominal, fully dissociated),
  2. conductivity WITH equilibrium free-ion speciation,
  3. the amount of speciation (percent of the cation bound in complexes),
  4. the conductivity difference (percent reduction).

Output: PDF + PNG per salt under plans/figures/mccleskey_speciation/.
Run from repo root in the PP env. Smoke-scale (~minutes).
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath('.'))
from PlanetProfile.Thermodynamics.Reaktoro.reaktoroProps import (  # noqa: E402
    SpeciesParser, RktHydroSpecies, RktConduct)

# PlanetProfile's config import side-effects enable usetex, where a bare %
# is a comment char that truncates titles — plain mathtext here.
plt.rcParams['text.usetex'] = False

OUT_DIR = 'plans/figures/mccleskey_speciation'
P_MPA = 10.0            # McCleskey has no pressure dependence; fixed for the map
T_K = np.linspace(265.0, 320.0, 23)
MOLAL = np.logspace(-3, np.log10(1.5), 16)   # 0.001 - 1.5 mol/kg

SALTS = {
    'NaCl': {'cation': 'Na+', 'ions': lambda m: f"Na+: {m}, Cl-: {m}"},
    'MgSO4': {'cation': 'Mg+2', 'ions': lambda m: f"Mg+2: {m}, SO4-2: {m}"},
}


def conduct_pair(salt, m):
    """(sigma_nominal(T), sigma_equilibrium(T), free_cation_fraction(T))."""
    comp = f"CustomSolution{salt}Doc = {SALTS[salt]['ions'](m)}"
    aq, ratio, solids, _ = SpeciesParser(comp, None)

    rc_nom = RktConduct(aq, ratio, solids, RktHydroSpecies(aq, ratio, solids),
                        use_equilibrium_speciation=False)
    rc_eq = RktConduct(aq, ratio, solids, RktHydroSpecies(aq, ratio, solids),
                       use_equilibrium_speciation=True)

    P = np.full_like(T_K, P_MPA)
    sig_nom = np.asarray(rc_nom(P, T_K), float)
    sig_eq = np.asarray(rc_eq(P, T_K), float)

    # Speciation amount: fraction of the cation NOT free at equilibrium.
    pH, spec, names, aff = rc_eq.fn_species_conduct(P, T_K)
    cat = SALTS[salt]['cation']
    i_cat = list(names).index(cat)
    free = np.ravel(np.asarray(spec[i_cat], float))
    bound_frac = 1.0 - free / m
    return sig_nom, sig_eq, np.clip(bound_frac, 0, 1)


def build_maps(salt):
    nom = np.zeros((MOLAL.size, T_K.size))
    eq = np.zeros_like(nom)
    bound = np.zeros_like(nom)
    for i, m in enumerate(MOLAL):
        nom[i], eq[i], bound[i] = conduct_pair(salt, m)
        print(f"{salt}: molality {m:.4g} done", flush=True)
    return nom, eq, bound


def plot_salt(salt, nom, eq, bound):
    TT, MM = np.meshgrid(T_K, MOLAL)
    dsig_pct = 100.0 * (nom - eq) / np.where(nom > 0, nom, np.nan)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True, sharey=True)
    panels = [
        (nom, r'$\sigma$ without speciation (S/m)', 'viridis'),
        (eq, r'$\sigma$ with equilibrium speciation (S/m)', 'viridis'),
        (100 * bound, 'cation bound in complexes (%)', 'magma'),
        (dsig_pct, r'conductivity reduction $(\sigma_{nom}-\sigma_{eq})/\sigma_{nom}$ (%)', 'magma'),
    ]
    for ax, (Z, title, cmap) in zip(axes.ravel(), panels):
        cf = ax.contourf(TT, MM, Z, levels=14, cmap=cmap)
        cs = ax.contour(TT, MM, Z, levels=cf.levels[::2], colors='k',
                        linewidths=0.4)
        ax.clabel(cs, fontsize=6, fmt='%.2g')
        ax.set_yscale('log')
        ax.set_title(title, fontsize=10)
        fig.colorbar(cf, ax=ax)
    for ax in axes[1]:
        ax.set_xlabel('T (K)')
    for ax in axes[:, 0]:
        ax.set_ylabel('molality (mol/kg)')
    fig.suptitle(
        f'{salt}: McCleskey (2012) conductivity, nominal vs Reaktoro '
        f'equilibrium speciation (P = {P_MPA:g} MPa; supcrt16, 1:1 ion-pair '
        f'complexes)', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    os.makedirs(OUT_DIR, exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(f'{OUT_DIR}/mccleskey_speciation_{salt}.{ext}', dpi=170)
    plt.close(fig)
    print(f'{salt}: figures written to {OUT_DIR}/')


if __name__ == '__main__':
    for salt in SALTS:
        maps = build_maps(salt)
        plot_salt(salt, *maps)
