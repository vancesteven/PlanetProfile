import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from matplotlib.patches import Wedge
from matplotlib.collections import PatchCollection
import os
from pathlib import Path
from Utilities.SetupInit import SetupFilenames

def GeneratePlots(Planet, Params):

    if Params.PLOT_GRAVITY: PlotGravPres(Planet, Params)
    if Params.PLOT_HYDROSPHERE: PlotHydrosphereProps(Planet, Params)
    if Planet.Do.Fe_CORE:
        if Params.PLOT_TRADEOFF_WCORE: PlotCoreTradeoff(Planet, Params)
    else:
        if Params.PLOT_TRADEOFF_NOCORE: PlotSilTradeoff(Planet, Params)
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

    fig.subplots_adjust(wspace = 0.5)
    fig.suptitle('Gravity and pressure')
    fig.savefig(Params.FigureFiles.vgrav, format=Params.figFormat, dpi=300)
    plt.close()

def PlotHydrosphereProps(Planet, Params):
    data = {'pressure': Planet.P_MPa[:Planet.Steps.nHydro],
            'density': Planet.rho_kgm3[:Planet.Steps.nHydro],
            'temp': Planet.T_K[:Planet.Steps.nHydro],
            'depth': Planet.z_m[:Planet.Steps.nHydro]/1000}
    fig, axes = plt.subplots(1, 2, figsize=Params.FigSize.vhydro)
    axes[0].plot('pressure', 'density', data=data)
    axes[0].set_xlabel('Pressure (MPa)')
    axes[0].set_ylabel('Density (kg/m$^3$)')

    axes[1].plot('temp', 'depth', data = data)
    axes[1].invert_yaxis()
    axes[1].set_xlabel('Temperature (K)')
    axes[1].set_ylabel('Depth (km)')

    fig.suptitle('Hydrosphere properties')
    fig.subplots_adjust(wspace=0.5)
    fig.savefig(Params.FigureFiles.vhydro, format=Params.figFormat, dpi=300)
    plt.close()


def PlotCoreTradeoff(Planet, Params):
    data = {'Rsil': Planet.Sil.Rtrade_m/1000,
            'RFe': Planet.Core.Rtrade_m/1000}
    fig, axes = plt.subplots(1, 1, figsize=Params.FigSize.vcore)
    axes.plot('Rsil', 'RFe', data = data)
    axes.set_xlabel('Iron core outer radius (km)')
    axes.set_ylabel('Silicate layer outer radius (km)')
    fig.suptitle(r'With Fe core. $C/MR^2$: $0.346\pm0.005$; $w$: $0\,\mathrm{wt}\%$; $\rho_\mathrm{sil}$: $' \
                 + str(round(Planet.Sil.rhoMean_kgm3)) + r'\,\mathrm{kg/m^3}$; $\rho_\mathrm{Fe}$: $' \
                 + str(round(Planet.Core.rhoMean_kgm3)) + r'\,\mathrm{kg/m^3}$')
    fig.savefig(Params.FigureFiles.vcore, format=Params.figFormat, dpi=300)
    plt.close()


def PlotSilTradeoff(Planet, Params):
    data = {'Rsil': Planet.Sil.Rtrade_m/1000,
            'rhoSil': Planet.Sil.rhoTrade_kgm3}
    fig, axes = plt.subplots(1, 1, figsize=Params.FigSize.vmant)
    axes.plot('rhoSil', 'Rsil', data = data)
    axes.set_xlabel('$\\rho_\mathrm{sil}$ (kg/m$^3$)')
    axes.set_ylabel('$R_\mathrm{sil}$ (km)')
    fig.suptitle(r'No Fe core. $C/MR^2$: $0.346\pm0.005$; $W$')
    fig.savefig(Params.FigureFiles.vmant, format=Params.figFormat, dpi=300)
    plt.close()

def PlotWedge(Planet, Params):
    fig, ax = plt.subplots()
    width = (math.pi / 7)*180/math.pi
    patches = []
    colors = []
    iPhaseTrans = 1+np.argwhere(Planet.phase[1:] != Planet.phase[:-1])
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
        elif Planet.phase[layerPhase] == 30:
            colors.append(Params.Colors.Clath)
        elif Planet.phase[layerPhase] == 50:
            colors.append(Params.Colors.Rock)
        elif Planet.phase[layerPhase] == 100:
            colors.append(Params.Colors.Core)
    phases = [Planet.phase[iShell] for iShell in iPhaseTrans]
    radii = [Planet.r_m[iShell]/Planet.Bulk.R_m for iShell in iPhaseTrans]
    for i in range(len(radii)):
        patches += [Wedge((0.5,0), radii[i], 90 - width, 90 + width, lw = 0.25)]
    p = PatchCollection(patches)

    p.set_color(colors)
    p.set_edgecolor('k')
    ax.add_collection(p)
    ax.set_aspect('equal')

    #fig.colorbar(p, ax = ax)
    fig.suptitle('Interior wedge diagram\n$T_b = ' + str(Planet.Bulk.Tb_K) + '\,\mathrm{K}' + ', Composition = ' + str(Planet.Ocean.comp) + ', Salinity = ' + str(Planet.Ocean.wOcean_ppt) + '\,\mathrm{ppt}$')
    plt.margins(0.02)
    fig.savefig(Params.FigureFiles.vwedg, format=Params.figFormat, dpi=300)

    plt.close()


