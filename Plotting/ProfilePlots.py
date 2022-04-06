import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from matplotlib.patches import Wedge
from matplotlib.collections import PatchCollection
import os
from pathlib import Path
from Utilities.SetupInit import SetupFilenames
from Utilities.defineStructs import Constants

def GeneratePlots(Planet, Params):

    if Params.PLOT_GRAVITY: PlotGravPres(Planet, Params)
    if Params.PLOT_HYDROSPHERE and not Planet.Do.NO_H2O: PlotHydrosphereProps(Planet, Params)
    if Params.PLOT_TRADEOFF:
        if Planet.Do.Fe_CORE: PlotCoreTradeoff(Planet, Params)
        else: PlotSilTradeoff(Planet, Params)
    if Params.PLOT_WEDGE: PlotWedge(Planet, Params)
    return


def PlotGravPres(Planet, Params):
    data = {'radius': Planet.r_m/1000,
            'grav': Planet.g_ms2,
            'pressure': Planet.P_MPa/1000}
    fig, axes = plt.subplots(1, 2, figsize=Params.FigSize.vgrav)
    axes[0].plot('grav', 'radius', data=data)
    axes[0].set_xlabel('Gravity (m/s$^2$)')
    axes[0].set_ylabel('Radius (km)')

    axes[1].plot('pressure', 'radius', data = data)
    axes[1].set_xlabel('Pressure (GPa)')
    axes[1].set_ylabel('$r_\mathrm{' + Planet.name + '}$')

    fig.subplots_adjust(wspace=0.5)
    fig.suptitle('Gravity and pressure')
    fig.savefig(Params.FigureFiles.vgrav, format=Params.figFormat, dpi=300)
    plt.close()
    return


def PlotHydrosphereProps(Planet, Params):
    data = {'pressure': Planet.P_MPa[:Planet.Steps.nHydro],
            'density': Planet.rho_kgm3[:Planet.Steps.nHydro],
            'temp': Planet.T_K[:Planet.Steps.nHydro],
            'depth': Planet.z_m[:Planet.Steps.nHydro]/1000}
    fig, axes = plt.subplots(1, 2, figsize=Params.FigSize.vhydro)

    # Plot reference profiles first, so they plot on bottom
    if Params.PLOT_REF:
        # Get strings for referencing and labeling
        wList = [f'{w:.0f}' for w in Params.wRef_ppt[Planet.Ocean.comp]]
        # Take care to only plot the values consistent with layer solutions
        iPlot = Params.Pref_MPa < np.max(Planet.P_MPa[:Planet.Steps.nHydro])
        data['Pref'] = Params.Pref_MPa[iPlot]
        # Plot all reference melting curve densities
        for i in range(Params.nRef):
            data[wList[i]] = Params.rhoRef_kgm3[i,iPlot]
            axes[0].plot('Pref', f'{wList[i]}', data=data, color=Params.refColor,
                         lw=Params.refLW, ls=Params.refLS[Planet.Ocean.comp])

    # Plot density vs. pressure curve for hydrosphere
    axes[0].plot('pressure', 'density', data=data)
    axes[0].set_xlabel('Pressure (MPa)')
    axes[0].set_ylabel('Density (kg/m$^3$)')

    # Plot thermal profile vs. depth in hydrosphere
    axes[1].plot('temp', 'depth', data=data)
    axes[1].invert_yaxis()
    axes[1].set_xlabel('Temperature (K)')
    axes[1].set_ylabel('Depth (km)')

    fig.suptitle('Hydrosphere properties')
    fig.subplots_adjust(wspace=0.5)
    fig.savefig(Params.FigureFiles.vhydro, format=Params.figFormat, dpi=300)
    plt.close()
    return


def PlotCoreTradeoff(Planet, Params):
    data = {'Rsil': Planet.Sil.Rtrade_m/1000,
            'RFe': Planet.Core.Rtrade_m/1000}
    fig, axes = plt.subplots(1, 1, figsize=Params.FigSize.vcore)
    axes.plot('Rsil', 'RFe', data = data)
    axes.set_xlabel('Iron core outer radius (km)')
    axes.set_ylabel('Silicate layer outer radius (km)')
    fig.suptitle(r'With Fe core. $C/MR^2$: $' + f'{Planet.Bulk.Cmeasured:.3f}\pm{Planet.Bulk.Cuncertainty:.3f}' +
                 r'$; $w$: $0\,\mathrm{wt}\%$; $\rho_\mathrm{sil}$: $' + \
                 f'{Planet.Sil.rhoMean_kgm3:.0f}' + r'\,\mathrm{kg/m^3}$; $\rho_\mathrm{Fe}$: $' + \
                 f'{Planet.Core.rhoMean_kgm3:.0f}' + r'\,\mathrm{kg/m^3}$')
    fig.savefig(Params.FigureFiles.vcore, format=Params.figFormat, dpi=300)
    plt.close()
    return


def PlotSilTradeoff(Planet, Params):
    data = {'Rsil': Planet.Sil.Rtrade_m/1000,
            'rhoSil': Planet.Sil.rhoTrade_kgm3}
    fig, axes = plt.subplots(1, 1, figsize=Params.FigSize.vmant)
    axes.plot('rhoSil', 'Rsil', data = data)
    axes.set_xlabel('$\\rho_\mathrm{sil}$ (kg/m$^3$)')
    axes.set_ylabel('Silicate layer outer radius (km)')
    fig.suptitle(r'No Fe core. $C/MR^2$: $0.346\pm0.005$; $W$')
    fig.savefig(Params.FigureFiles.vmant, format=Params.figFormat, dpi=300)
    plt.close()
    return


def PlotWedge(Planet, Params):
    fig, ax = plt.subplots()
    width = (math.pi / 7)*180/math.pi  # angular width of wedge to be plotted
    patches = []  # for storing wedge objects
    colors = []  # colors for layers
    iPhaseTrans = 1+np.where(Planet.phase[1:] != Planet.phase[:-1])[0]  # finds indexes of transitions between layers
    iPhaseTrans = np.insert(iPhaseTrans, 0, 0) # this makes sure the ice phase is included
    for layerPhase in iPhaseTrans:

        if Planet.phase[layerPhase] == 0:
            colors.append(Params.Colors.OceanTop)
        elif Planet.phase[layerPhase] == 1:
            colors.append(Params.Colors.IceI)
        elif Planet.phase[layerPhase] == 2:
            colors.append(Params.Colors.IceII)
        elif Planet.phase[layerPhase] == 3:
            colors.append(Params.Colors.IceIII)
        elif Planet.phase[layerPhase] == 5:
            colors.append(Params.Colors.IceV)
        elif Planet.phase[layerPhase] == 6:
            colors.append(Params.Colors.IceVI)
        elif Planet.phase[layerPhase] == Constants.phaseClath:
            colors.append(Params.Colors.Clath)
        elif Planet.phase[layerPhase] == Constants.phaseSil:
            colors.append(Params.Colors.Rock)
        elif Planet.phase[layerPhase] == Constants.phaseFe:
            colors.append(Params.Colors.Core)
    phases = [Planet.phase[iShell] for iShell in iPhaseTrans]  # stores phase of particular layer
    radii = [Planet.r_m[iShell]/Planet.Bulk.R_m for iShell in iPhaseTrans]  # normalizes radii of layer transitions

    funNum = 1
    im = None

    for i, radius in enumerate(radii):
        iCol = i % np.size(colors)
        print(i, colors[iCol])

        patches.append(Wedge((0.5,0), radius, 90 - width, 90 + width, lw = 0.25, fc = "none" if i == funNum else colors[iCol], ec="k", zorder=i))  # creating wedges
        ax.add_patch(patches[-1])

        if i == funNum:
            print("Draw time!")
            delta = 0.025
            x = y = np.arange(0, 1.0, delta)
            X, Y = np.meshgrid(x, y)
            Z1 = np.exp(-(X-0.5) ** 2 - Y ** 2)
            Z2 = np.exp(-(X - 1.5) ** 2 - (Y - 1) ** 2)
            Z = (Z1 - Z2) * 2

            Z = ((X+0.5) ** 0.5 - Y ** 0.5)**2

            im = plt.imshow(Z, interpolation='bilinear', cmap=mpl.cm.bone,
                           origin='lower', extent=[0, 1, 0, 1],
                           clip_path=patches[-1], clip_on=True)
            im.set_clip_path(patches[-1])

    ax.set_aspect('equal')

    #fig.colorbar(p, ax = ax)
    if Planet.Ocean.comp == 'MgSO4':
        compstr = 'MgSO$_4$'
    elif Planet.Ocean.comp == 'PureH2O':
        compstr = 'Pure H$_2$O'
    else:
        compstr = Planet.Ocean.comp
    fig.suptitle(f'Interior wedge diagram\n$T_b = {Planet.Bulk.Tb_K}\,\mathrm{{K}}$, Composition = {compstr}, Salinity = ${Planet.Ocean.wOcean_ppt}\,\mathrm{{g/kg}}$')
    plt.margins(0.02)
    fig.savefig(Params.FigureFiles.vwedg, format=Params.figFormat, dpi=300)

    plt.close()
    return

