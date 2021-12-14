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
    plt.suptitle('Gravity and Pressure')
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

    plt.suptitle('Hydrosphere Properties')
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
    data = {'r': Planet.r_m/1000,
            'phase': Planet.phase}
    fig, ax = plt.subplots()
    width = (math.pi / 7)*180/math.pi
    patches = []
    colors = []
    #length = len(Planet.r_m)
    #widths = []
    iPhaseTrans = 1+np.argwhere(Planet.phase[1:] != Planet.phase[:-1])
    iPhaseTrans = np.insert(iPhaseTrans, 0, 0) # this makes sure the ice phase is included
    #widths = iPhaseTrans/length
    phases = [Planet.phase[iShell] for iShell in iPhaseTrans]
    radii = [Planet.r_m[iShell]/Planet.Bulk.R_m for iShell in iPhaseTrans]
    #print(phases)
    #print(radii)
    patches = [Wedge((0.5,0), radii[1], 90 - width, 90 + width, 1-radii[0], color='b') ]
    #for i in range(len(radii)):
       # patches += [Wedge((0.5,0), radii[i], 90 - width, 90 + width, 1-radii[i-1])]
    print(patches)
    p = PatchCollection(patches)
    #colors = colors.append(Params.Colors.IceI)

    #p.set_array(colors)
    ax.add_collection(p)
    plt.title('Interior Wedge Diagram')
    plt.suptitle('$Tb_K = \mathrm{' + str(Planet.Bulk.Tb_K) + '}$')
    plt.show()
'''
    for layerPhase in iPhaseTrans:
        if layerPhase == 0:
            colors.append(Params.Colors.OceanTop)
        elif layerPhase == 1:
            colors.append(Params.Colors.IceI)
        elif layerPhase == 2:
            colors.append(Params.Colors.IceII)
        elif layerPhase == 3:
            colors.append(Params.Colors.IceIII)
        elif layerPhase == 5:
            colors.append(Params.Colors.IceV)
        elif layerPhase == 6:
            colors.append(Params.Colors.IceVI)
        elif layerPhase == 30:
            colors.append(Params.Colors.Clath)
        elif layerPhase == 50:
            colors.append(Params.Colors.Rock)
        elif layerPhase == 100:
            colors.append(Params.Colors.Core)
    print(colors)
'''




