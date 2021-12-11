import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

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

    plt.title('Interior Wedge Diagram')
    plt.suptitle('$Tb_K = \mathrm{' + str(Planet.Bulk.Tb_K) + '}$')
    plt.show()
    pass


def  GetWedgeColors():
    colors.IonosphereTop = [1, 0, 1] # matlab's magenta
    colors.Ionosphere = [1, 0, 1] # matlab's magenta
    colors.IonosphereBot = [1, 0, 1] # matlab's magenta
    colors.IceI = [150, 226, 241]/255
    colors.IceII = [3, 169, 252]/255
    colors.IceIII = [150, 226, 241]/255
    colors.IceV = [150, 226, 241]/255
    colors.IceVI =  [150, 226, 241]/255
    colors.IceVII =  [44, 115, 150]/255
    colors.OceanTop = [134, 149, 201]/255  # Darker and richer
    colors.OceanBot = [45, 55, 100]/255
    colors.Rock = [101, 46, 11]/255
    colors.Core = [141, 122, 121]/255

