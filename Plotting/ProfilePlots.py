import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from matplotlib.patches import Wedge
from matplotlib.collections import PatchCollection


def GeneratePlots(Planet, Params):

    if Params.PLOT_GRAVITY: PlotGravPres(Planet, Params)
    if Params.PLOT_HYDROSPHERE: PlotHydrosphereProps(Planet, Params)
    if Params.PLOT_TRADEOFF: PlotTradeoff(Planet, Params)
    if Params.PLOT_WEDGE: PlotWedge(Planet, Params)
    return


def PlotGravPres(Planet, Params):
    data = {'Radius': Planet.r_m/1000,
            'grav': Planet.g_ms2,
            'pressure': Planet.P_MPa}
    plt.subplot(1, 2, 1)
    plt.plot('grav', 'Radius', data=data)
    plt.xlabel('g (m/s$^2$)')
    plt.ylabel('Radius (km)')

    plt.subplot(1,2,2)
    plt.plot('pressure', 'Radius', data = data)
    plt.xlabel('Pressure (MPa)')
    plt.ylabel('$r_\mathrm{' + Planet.name + '}$')

    plt.subplots_adjust(wspace = 0.5)
    plt.suptitle('Gravity and pressure')
    plt.show()

def PlotHydrosphereProps(Planet, Params):
    data = {'pressure': Planet.P_MPa[:Planet.Steps.nHydro],
            'density': Planet.rho_kgm3[:Planet.Steps.nHydro],
            'temp': Planet.T_K[:Planet.Steps.nHydro],
            'depth': Planet.z_m[:Planet.Steps.nHydro]/1000}

    plt.subplot(1,2,1)
    plt.plot('pressure', 'density', data=data)
    plt.xlabel('Pressure (MPa)')
    plt.ylabel('Density (kg/m$^3$)')

    plt.subplot(2,2,2)
    plt.plot('temp', 'depth', data = data)
    plt.gca().invert_yaxis()
    plt.xlabel('Temperature (K)')
    plt.ylabel('Depth (km)')

    plt.suptitle('Hydrosphere properties')
    plt.subplots_adjust(wspace=0.5)
    plt.show()


def PlotTradeoff(Planet, Params):
    data = {'Rsil': Planet.Sil.Rtrade_m/1000,
            'Rfe': Planet.Core.Rtrade_m/1000}
    plt.plot('Rsil', 'Rfe', data = data)
    plt.xlabel('RFe (km)')
    plt.ylabel('Rsil (km)')
    plt.title(r'With Fe core. $C/MR^2$: $0.346\pm0.005$; $w$: $0\,\mathrm{wt}\%$; $\rho_\mathrm{sil}$: $' \
              + str(round(Planet.Sil.rhoMean_kgm3)) + r'\,\mathrm{kg/m^3}$; $\rho_\mathrm{Fe}$: $' \
              + str(round(Planet.Core.rhoMean_kgm3)) + r'\,\mathrm{kg/m^3}$')
    plt.show()


def PlotWedge(Planet, Params):
    fig, ax = plt.subplots()
    width = (math.pi / 7)*180/math.pi
    patches = []
    colors = []
    iPhaseTrans = 1+np.argwhere(Planet.phase[1:] != Planet.phase[:-1])
    iPhaseTrans = np.insert(iPhaseTrans, 0, 0) # this makes sure the ice phase is included
    #iPhaseTrans = np.insert(iPhaseTrans, int([Planet.Bulk.R_m]))  # this makes sure the ice phase is included
    print(iPhaseTrans)
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
    print(colors)
    phases = [Planet.phase[iShell] for iShell in iPhaseTrans]
    radii = [Planet.r_m[iShell]/Planet.Bulk.R_m for iShell in iPhaseTrans]
    #print(phases)
    print(radii)
    for i in range(len(radii)):
        patches += [Wedge((0.5,0), radii[i], 90 - width, 90 + width, lw = 0.25)]#, #1-radii[i-1] )]#color = colors[i])]
    print(patches)
    p = PatchCollection(patches)

    p.set_color(colors)
    p.set_edgecolor('k')
    ax.add_collection(p)
    ax.set_aspect('equal')

    #fig.colorbar(p, ax = ax)
    plt.suptitle('Interior wedge diagram')
    plt.title('$T_b = ' + str(Planet.Bulk.Tb_K) + '\,\mathrm{K}' + ', Composition = ' + str(Planet.Ocean.comp) + ', Salinity = ' + str(Planet.Ocean.wOcean_ppt) + '\,\mathrm{ppt}$')
    plt.margins(0.02)
    plt.show()




