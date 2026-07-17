"""Predicted vs measured conductivity: Mahboub et al. (2026) comparisons.

Measured data: Table S2 of the supplement (papers/mahboub2026electricalSOM.pdf),
transcribed below — NaCl, MgSO4, NH4Cl, Na2CO3, and the 1:1:1
NaCl:MgSO4:Na2SO4 mixture, 263.15-298.15 K. (NOTE: Table S2 labels the
mixture Na2SO4 while the main-text Fig. 5 caption says Na2CO3 — we follow
the data table.)

Predictions: McCleskey (2012) via PlanetProfile's RktConduct in BOTH modes:
- nominal (fully dissociated — reproduces the paper's MC12 dashed curves),
- equilibrium speciation (2026-07 change: Reaktoro free-ion molalities,
  1:1 ion-pair complexes in the system).

Outputs (PDF+PNG) in plans/figures/mccleskey_speciation/:
- mahboub2026_fig3_analog: sigma vs molality per salt (3 temperatures),
  with percent-deviation-from-data panels for both models.
- mahboub2026_fig4_analog: sigma vs temperature per salt (2-3 concs).
- mahboub2026_mixture: the 1:1:1 mixture vs concentration.

Run from repo root in the PP env (reaktoro required). Smoke-scale.
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

plt.rcParams['text.usetex'] = False

OUT_DIR = 'plans/figures/mccleskey_speciation'
P_MPA = 0.1
T7 = np.array([263.15, 267.15, 270.15, 272.15, 278.15, 293.15, 298.15])
# Mixture rows carry only four values; temperature assignment inferred from
# the value progression vs the single-salt rows (263 -> ~272 -> 293 -> 298).
# TO VERIFY against the paper's supplemental script (MahboubEtAl2026.py).
T_MIX = np.array([263.15, 272.15, 293.15, 298.15])

# ---- Table S2 (sigma S/m, +/- absolute) --------------------------------
S2 = {
    'NaCl': {
        0.1711: [(0.7852, .0234), (0.8649, .0278), (0.9823, .0528), (0.9847, .0293), (1.1820, .0418), (1.6113, .0505), (1.7727, .0557)],
        0.5133: [(2.1743, .1021), (2.3533, .0810), (2.5050, .0964), (2.6990, .0951), (3.1343, .0959), (4.4307, .1323), (4.8737, .1455)],
        0.8556: [(3.1820, .1056), (3.5487, .1139), (3.9440, .1647), (4.1137, .1482), (4.8507, .1646), (6.8430, .2098), (7.5273, .2307)],
        1.2834: [(4.4507, .1768), (5.1140, .1648), (5.4793, .1638), (5.6953, .1748), (6.8403, .2168), (9.5800, .2871), (10.5380, .3158)],
        1.7112: [(5.7470, .2133), (6.3817, .2213), (7.0060, .2245), (7.1933, .2636), (8.6653, .2786), (12.0500, .3591), (13.2550, .3950)],
        2.5667: [(7.6640, .3239), (8.5290, .2586), (9.2313, .3059), (9.5450, .2906), (11.3400, .3618), (16.0033, .4805), (17.6037, .5285)],
    },
    'MgSO4': {
        0.0249: [(np.nan, np.nan), (0.1889, .0188), (0.1683, .0465), (0.1930, .0102), (0.2209, .0097), (0.3112, .0128), (0.3423, .0140)],
        0.3323: [(1.1090, .0330), (1.1723, .0458), (1.2797, .0469), (1.3787, .0580), (1.5993, .0498), (2.2943, .0710), (2.5237, .0781)],
        0.6646: [(1.5500, .0461), (1.8210, .0588), (1.9340, .0866), (2.1143, .0673), (2.5187, .0761), (3.6417, .1084), (4.0060, .1193)],
        0.9970: [(2.0370, .0606), (2.2263, .0805), (2.4713, .0758), (2.6320, .0805), (3.1110, .0960), (4.4403, .2288), (4.8847, .2515)],
        1.4124: [(2.1787, .0666), (2.4390, .0884), (2.7707, .1322), (2.9163, .0886), (3.5000, .1106), (5.1337, .1531), (5.6470, .1685)],
        1.6616: [(2.2353, .0774), (2.4713, .0955), (2.7767, .0830), (2.9483, .1109), (3.5430, .1082), (5.3050, .1740), (5.8353, .1915)],
    },
    'NH4Cl': {
        0.1870: [(0.9741, .1389), (1.1877, .0558), (1.2603, .0385), (1.3343, .0452), (1.5493, .0516), (2.1160, .0634), (2.3276, .0698)],
        0.9348: [(4.7517, .1697), (5.2710, .2081), (5.5857, .1879), (5.9827, .1975), (6.9087, .2165), (9.4327, .2826), (10.3759, .3109)],
        1.4021: [(6.9420, .2291), (7.6387, .2418), (8.2280, .2513), (8.6027, .2602), (9.9713, .2999), (13.5000, .4154), (14.8500, .4569)],
        1.8695: [(9.0340, .3579), (9.9220, .3127), (10.6400, .3579), (11.2633, .4601), (13.0533, .4011), (17.3867, .5288), (19.1253, .5817)],
    },
    'Na2CO3': {
        0.0472: [(0.2990, .0094), (0.3786, .0216), (0.6193, .0194), (0.6550, .0206), (0.5143, .0164), (0.7475, .0235), (0.8222, .0259)],
        0.0943: [(0.5329, .0161), (0.6548, .0197), (1.1038, .0335), (1.1667, .0354), (0.9145, .0307), (1.3323, .0403), (1.4656, .0444)],
        0.2359: [(1.1184, .0343), (1.4057, .0476), (2.3187, .0714), (2.4506, .0755), (1.8973, .0565), (2.7960, .0859), (3.0756, .0945)],
        0.3774: [(1.5959, .0517), (1.9567, .0636), (3.3086, .1075), (3.4932, .1134), (2.6770, .0806), (3.9897, .1292), (4.3886, .1421)],
        0.5189: [(1.9823, .0603), (2.4180, .0733), (4.1115, .1259), (4.3404, .1329), (3.3867, .1012), (4.9557, .1507), (5.4512, .1657)],
    },
    'Mixture': {   # NaCl:MgSO4:Na2SO4 (1:1:1); columns follow T_MIX
        0.05: [(0.6665, .0198), (0.8441, .0264), (1.3747, .0085), (1.51213, .0094)],
        0.10: [(1.1355, .0567), (1.5227, .0705), (2.4587, .0781), (2.70453, .0859)],
        0.40: [(3.2887, .1140), (4.0220, .1257), (6.6827, .0674), (7.35093, .0741)],
        0.70: [(3.0870, .2373), (3.5753, .1230), (5.8820, .1847), (6.47020, .2032)],
    },
}

IONS = {
    'NaCl': lambda m: f"Na+: {m}, Cl-: {m}",
    'MgSO4': lambda m: f"Mg+2: {m}, SO4-2: {m}",
    'NH4Cl': lambda m: f"NH4+: {m}, Cl-: {m}",
    'Na2CO3': lambda m: f"Na+: {2 * m}, CO3-2: {m}",
    # 1:1:1 NaCl + MgSO4 + Na2SO4, each at m:
    'Mixture': lambda m: (f"Na+: {3 * m}, Cl-: {m}, Mg+2: {m}, SO4-2: {2 * m}"),
}


def predict(salt, m, T_K):
    """(sigma_nominal(T), sigma_speciated(T)); NaN where a path fails."""
    comp = f"CustomSolution{salt}Cmp = {IONS[salt](m)}"
    aq, ratio, solids, _ = SpeciesParser(comp, None)
    P = np.full_like(np.asarray(T_K, float), P_MPA)
    out = []
    for eq in (False, True):
        try:
            rc = RktConduct(aq, dict(ratio), solids,
                            RktHydroSpecies(aq, dict(ratio), solids),
                            use_equilibrium_speciation=eq)
            out.append(np.asarray(rc(P, np.asarray(T_K, float)), float))
        except Exception as e:
            print(f"{salt} m={m} eq={eq}: FAILED ({e})")
            out.append(np.full(np.asarray(T_K).shape, np.nan))
    return out


def fig3_analog():
    """sigma vs molality per salt at three temperatures + deviation rows."""
    salts = ['NaCl', 'MgSO4', 'NH4Cl', 'Na2CO3']
    temps = [263.15, 278.15, 298.15]
    it = [int(np.where(T7 == t)[0][0]) for t in temps]
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(temps)))

    fig, axes = plt.subplots(2, len(salts), figsize=(4.2 * len(salts), 7.4),
                             sharex='col')
    for j, salt in enumerate(salts):
        mol = np.array(sorted(S2[salt]))
        meas = np.array([[S2[salt][m][k][0] for k in range(len(T7))]
                         for m in mol])
        errs = np.array([[S2[salt][m][k][1] for k in range(len(T7))]
                         for m in mol])
        nom = np.full((mol.size, len(temps)), np.nan)
        spc = np.full_like(nom, np.nan)
        for i, m in enumerate(mol):
            n, s = predict(salt, m, temps)
            nom[i], spc[i] = n, s
        ax, axd = axes[0, j], axes[1, j]
        for c, (t, k) in zip(colors, zip(temps, it)):
            ax.errorbar(mol, meas[:, k], yerr=errs[:, k], fmt='o', ms=4,
                        color=c, label=f'{t:.0f} K data')
            ki = temps.index(t)
            ax.plot(mol, nom[:, ki], '--', color=c, lw=1.2)
            ax.plot(mol, spc[:, ki], '-', color=c, lw=1.6)
            with np.errstate(invalid='ignore', divide='ignore'):
                axd.plot(mol, 100 * (nom[:, ki] - meas[:, k]) / meas[:, k],
                         '--', color=c, lw=1.2)
                axd.plot(mol, 100 * (spc[:, ki] - meas[:, k]) / meas[:, k],
                         '-', color=c, lw=1.6)
        axd.axhline(0, color='0.5', lw=0.6)
        ax.set_title(salt)
        ax.set_xscale('log')
        axd.set_xscale('log')
        axd.set_xlabel('molality (mol/kg)')
        if j == 0:
            ax.set_ylabel(r'$\sigma$ (S/m)')
            axd.set_ylabel('model − data (%)')
            ax.legend(fontsize=7)
    fig.suptitle('Mahboub et al. (2026) Table S2 vs McCleskey (2012): '
                 'dashed = nominal (fully dissociated, as in the paper), '
                 'solid = with Reaktoro equilibrium speciation', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig, 'mahboub2026_fig3_analog'


def fig4_analog():
    """sigma vs temperature per salt at 2 concentrations."""
    picks = {'NaCl': [0.5133, 2.5667], 'MgSO4': [0.3323, 1.6616],
             'NH4Cl': [0.9348, 1.8695], 'Na2CO3': [0.0943, 0.5189]}
    fig, axes = plt.subplots(1, 4, figsize=(16.5, 4.4))
    for ax, (salt, mols) in zip(axes, picks.items()):
        colors = plt.cm.plasma(np.linspace(0.2, 0.75, len(mols)))
        for c, m in zip(colors, mols):
            meas = np.array([v[0] for v in S2[salt][m]])
            errs = np.array([v[1] for v in S2[salt][m]])
            nom, spc = predict(salt, m, T7)
            ax.errorbar(T7, meas, yerr=errs, fmt='o', ms=4, color=c,
                        label=f'{m:g} mol/kg')
            ax.plot(T7, nom, '--', color=c, lw=1.2)
            ax.plot(T7, spc, '-', color=c, lw=1.6)
        ax.set_title(salt)
        ax.set_xlabel('T (K)')
        ax.legend(fontsize=7)
    axes[0].set_ylabel(r'$\sigma$ (S/m)')
    fig.suptitle('Conductivity vs temperature: data (points), nominal MC12 '
                 '(dashed), speciated MC12 (solid)', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig, 'mahboub2026_fig4_analog'


def fig_mixture():
    """1:1:1 NaCl:MgSO4:Na2SO4 mixture vs concentration."""
    mol = np.array(sorted(S2['Mixture']))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, T_MIX.size))
    fig, (ax, axd) = plt.subplots(2, 1, figsize=(6.2, 7.2), sharex=True)
    nom = np.full((mol.size, T_MIX.size), np.nan)
    spc = np.full_like(nom, np.nan)
    for i, m in enumerate(mol):
        nom[i], spc[i] = predict('Mixture', m, T_MIX)
    for k, (c, t) in enumerate(zip(colors, T_MIX)):
        meas = np.array([S2['Mixture'][m][k][0] for m in mol])
        errs = np.array([S2['Mixture'][m][k][1] for m in mol])
        ax.errorbar(mol, meas, yerr=errs, fmt='o', ms=4, color=c,
                    label=f'{t:.0f} K data')
        ax.plot(mol, nom[:, k], '--', color=c, lw=1.2)
        ax.plot(mol, spc[:, k], '-', color=c, lw=1.6)
        with np.errstate(invalid='ignore', divide='ignore'):
            axd.plot(mol, 100 * (nom[:, k] - meas) / meas, '--', color=c)
            axd.plot(mol, 100 * (spc[:, k] - meas) / meas, '-', color=c)
    axd.axhline(0, color='0.5', lw=0.6)
    ax.set_ylabel(r'$\sigma$ (S/m)')
    ax.legend(fontsize=8)
    axd.set_xlabel('concentration of each salt (mol/kg)')
    axd.set_ylabel('model − data (%)')
    fig.suptitle('1:1:1 NaCl:MgSO4:Na2SO4 (Table S2; Fig. 5 caption says '
                 'Na2CO3 — flagged): dashed = nominal, solid = speciated',
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig, 'mahboub2026_mixture'


if __name__ == '__main__':
    os.makedirs(OUT_DIR, exist_ok=True)
    for builder in (fig3_analog, fig4_analog, fig_mixture):
        fig, name = builder()
        for ext in ('pdf', 'png'):
            fig.savefig(f'{OUT_DIR}/{name}.{ext}', dpi=170)
        plt.close(fig)
        print(f'wrote {OUT_DIR}/{name}.(pdf|png)')
