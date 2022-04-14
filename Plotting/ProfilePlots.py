import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from matplotlib.patches import Wedge
from matplotlib.collections import PatchCollection
import os
import logging as log
from pathlib import Path
from Utilities.SetupInit import SetupFilenames
from Utilities.defineStructs import Constants
from Thermodynamics.RefProfiles.RefProfiles import CalcRefProfiles, ReloadRefProfiles

def GeneratePlots(PlanetList, Params):

    # Handle refprofiles first, so we can print log messages before silencing them
    if Params.PLOT_HYDROSPHERE and not PlanetList[0].Do.NO_H2O:
        if Params.CALC_NEW_REF:
            # Calculate reference profiles showing melting curves for
            # several salinities specified in config.py
            Params = CalcRefProfiles(PlanetList, Params)
        else:
            # Reload refprofiles for this composition
            Params = ReloadRefProfiles(PlanetList, Params)

    log.warning('Temporarily quieting INFO and DEBUG messages due to a high number of current latex errors.')
    saveLevel = log.getLogger().getEffectiveLevel()
    log.getLogger().setLevel(log.WARN)

    if Params.PLOT_GRAVITY: PlotGravPres(PlanetList, Params)
    if Params.PLOT_HYDROSPHERE and not PlanetList[0].Do.NO_H2O: PlotHydrosphereProps(PlanetList, Params)
    if Params.PLOT_TRADEOFF:
        if Planet.Do.Fe_CORE: PlotCoreTradeoff(PlanetList, Params)
        else: PlotSilTradeoff(PlanetList, Params)
    if Params.PLOT_WEDGE: PlotWedge(PlanetList, Params)

    log.getLogger().setLevel(saveLevel)

    return


def PlotGravPres(PlanetList, Params):
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
    fig.suptitle(f'{PlanetList[0].name} gravity and pressure')
    fig.savefig(Params.FigureFiles.vgrav, format=Params.figFormat, dpi=300)
    plt.close()
    return


def PlotHydrosphereProps(PlanetList, Params):
    # Generate canvas and add labels
    fig, axes = plt.subplots(1, 2, figsize=Params.FigSize.vhydro)
    axes[0].set_xlabel('Pressure (MPa)')
    axes[0].set_ylabel('Density (kg/m$^3$)')
    axes[1].invert_yaxis()
    axes[1].set_xlabel('Temperature (K)')
    axes[1].set_ylabel('Depth (km)')
    fig.subplots_adjust(wspace=0.5)
    fig.suptitle(f'{PlanetList[0].name} hydrosphere properties')


    # Plot reference profiles first, so they plot on bottom of everything
    if Params.PLOT_REF:
        # Keep track of which reference profiles have been plotted so that we do each only once
        comps = np.unique([Planet.Ocean.comp for Planet in PlanetList])
        newRef = {comp:True for comp in comps}

        # Get max pressure among all profiles so we know how far out to plot refs
        Plist = np.concatenate([Planet.P_MPa[:Planet.Steps.nHydro] for Planet in PlanetList])
        Pmax_MPa = np.max(Plist)

        for Planet in PlanetList:
            if newRef[Planet.Ocean.comp]:
                # Get strings for referencing and labeling
                wList = f'$\\rho_\mathrm{{melt}}$ \ce{{{Planet.Ocean.comp}}} \\{{'
                wList += ', '.join([f'{w:.0f}' for w in Params.wRef_ppt[Planet.Ocean.comp]])
                wList += '\}\,ppt'
                # Take care to only plot the values consistent with layer solutions
                iPlot = Params.Pref_MPa[Planet.Ocean.comp] < Pmax_MPa
                # Plot all reference melting curve densities
                for i in range(Params.nRef[Planet.Ocean.comp]):
                    thisRef, = axes[0].plot(Params.Pref_MPa[Planet.Ocean.comp][iPlot],
                                            Params.rhoRef_kgm3[Planet.Ocean.comp][i,iPlot],
                                            color=Params.refColor,
                                            lw=Params.refLW,
                                            ls=Params.refLS[Planet.Ocean.comp])
                    if Params.refsInLegend and i == 0: thisRef.set_label(wList)
                newRef[Planet.Ocean.comp] = False

    # Now plot all profiles together
    for Planet in PlanetList:

        # Plot density vs. pressure curve for hydrosphere
        axes[0].plot(Planet.P_MPa[:Planet.Steps.nHydro], Planet.rho_kgm3[:Planet.Steps.nHydro], label=Planet.label)
        # Plot thermal profile vs. depth in hydrosphere
        axes[1].plot(Planet.T_K[:Planet.Steps.nHydro], Planet.z_m[:Planet.Steps.nHydro]/1e3)

    if Params.LEGEND:
        box1 = axes[0].get_position()
        fig.legend(loc=Params.LegendPosition)
    fig.savefig(Params.FigureFiles.vhydro, format=Params.figFormat, dpi=300)
    plt.close()
    return


def PlotCoreTradeoff(PlanetList, Params):
    data = {'Rsil': Planet.Sil.Rtrade_m/1000,
            'RFe': Planet.Core.Rtrade_m/1000}
    fig, axes = plt.subplots(1, 1, figsize=Params.FigSize.vcore)
    axes.plot('Rsil', 'RFe', data = data)
    axes.set_xlabel('Iron core outer radius (km)')
    axes.set_ylabel('Silicate layer outer radius (km)')
    fig.suptitle(f'{PlanetList[0].name} with Fe core. $C/MR^2$: ${Planet.Bulk.Cmeasured:.3f}\pm{Planet.Bulk.Cuncertainty:.3f}' +
                 r'$; $w$: $0\,\mathrm{wt}\%$; $\rho_\mathrm{sil}$: $' + \
                 f'{Planet.Sil.rhoMean_kgm3:.0f}' + r'\,\mathrm{kg/m^3}$; $\rho_\mathrm{Fe}$: $' + \
                 f'{Planet.Core.rhoMean_kgm3:.0f}' + r'\,\mathrm{kg/m^3}$')
    fig.savefig(Params.FigureFiles.vcore, format=Params.figFormat, dpi=300)
    plt.close()
    return


def PlotSilTradeoff(PlanetList, Params):
    data = {'Rsil': Planet.Sil.Rtrade_m/1000,
            'rhoSil': Planet.Sil.rhoTrade_kgm3}
    fig, axes = plt.subplots(1, 1, figsize=Params.FigSize.vmant)
    axes.plot('rhoSil', 'Rsil', data = data)
    axes.set_xlabel('$\\rho_\mathrm{sil}$ (kg/m$^3$)')
    axes.set_ylabel('Silicate layer outer radius (km)')
    fig.suptitle(f'{PlanetList[0].name} no Fe core. $C/MR^2$: $0.346\pm0.005$; $W$')
    fig.savefig(Params.FigureFiles.vmant, format=Params.figFormat, dpi=300)
    plt.close()
    return


def PlotWedge(PlanetList, Params):
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
        compstr = '\ce{MgSO4}'
    elif Planet.Ocean.comp == 'PureH2O':
        compstr = 'Pure \ce{H_2O}'
    else:
        compstr = Planet.Ocean.comp
    fig.suptitle(f'{PlanetList[0].name} wedge diagram\n$T_b = {Planet.Bulk.Tb_K}\,\mathrm{{K}}$, Composition = {compstr}, Salinity = ${Planet.Ocean.wOcean_ppt}\,\mathrm{{g/kg}}$')
    plt.margins(0.02)
    fig.savefig(Params.FigureFiles.vwedg, format=Params.figFormat, dpi=300)

    plt.close()
    return


def PlotInductOgram(Induction, Params):
    log.warning('Temporarily quieting INFO and DEBUG messages due to a high number of current latex errors.')
    saveLevel = log.getLogger().getEffectiveLevel()
    log.getLogger().setLevel(log.WARN)
    inductionTitle = f'\\textbf{{{Induction.bodyname} induction response}}'

    
    for z, name, fLabel in zip([Induction.Amp, Induction.Bix_nT, Induction.Biy_nT, Induction.Biz_nT],
                          ['Amplitude $A$', '$B_x$ component', '$B_y$ component', '$B_z$ component'],
                               ['Amp',         'Bx',             'By',            'Bz']):

        # Generate canvas and add labels
        fig, axes = plt.subplots(1, 2, figsize=Params.FigSize.induct)
        fig.subplots_adjust(wspace=0.5)
        fig.suptitle(inductionTitle)
        [ax.set_xlabel('Mean conductivity $\overline{\sigma}$ ($\si{S/m}$)') for ax in axes]
        [ax.set_ylabel('Ocean thickness $D$ ($\si{km}$)') for ax in axes]
        [ax.set_xlim([10**Params.sigmaMin[Induction.bodyname], 10**Params.sigmaMax[Induction.bodyname]]) for ax in axes]
        [ax.set_ylim([10**Params.Dmin[Induction.bodyname], 10**Params.Dmax[Induction.bodyname]]) for ax in axes]
        [ax.set_xscale('log') for ax in axes]
        [ax.set_yscale('log') for ax in axes]

        zContours = [axes[0].contour(Induction.sigmaMean_Sm, Induction.D_km, z[i, ...],
                         colors=Params.Colors.Induction[T], linestyles=Params.LS_Induction[T],
                         linewidths=Params.LW_Induction[T], levels=Params.cLevels[Induction.bodyname][T][fLabel])
                         for i, T in enumerate(Induction.Texc_hr.keys())]
        phaseContours = [axes[1].contour(Induction.sigmaMean_Sm, Induction.D_km, Induction.phase[i, ...],
                         colors=Params.Colors.Induction[T], linestyles=Params.LS_Induction[T],
                         linewidths=Params.LW_Induction[T], levels=Params.cLevels[Induction.bodyname][T]['phase'])
                         for i, T in enumerate(Induction.Texc_hr.keys())]
        if Params.inductOtype == 'sigma':
            [axes[0].clabel(zContours[i], fmt=Params.cFmt[Induction.bodyname][T][fLabel],
                            fontsize=Params.cLabelSize, inline_spacing=Params.cLabelPad)
                            for i, T in enumerate(Induction.Texc_hr.keys())]
            [axes[1].clabel(phaseContours[i], fmt=Params.cFmt[Induction.bodyname][T]['phase'],
                            fontsize=Params.cLabelSize, inline_spacing=Params.cLabelPad)
                            for i, T in enumerate(Induction.Texc_hr.keys())]

        if Params.LEGEND:
            lines = np.array([contour.legend_elements()[0][0] for contour in phaseContours])
            labels = np.array([f'{T_h:.2f} h' for T_h in Induction.Texc_hr.values()])
            iSort = np.argsort(list(Induction.Texc_hr.values()))
            axes[1].legend(lines[iSort], labels[iSort], framealpha=Params.cLegendOpacity)

        fig.savefig(Params.FigureFiles.sigma[fLabel], format=Params.figFormat)
        plt.close()

        if Params.inductOtype != 'sigma':
            fig, axes = plt.subplots(1, 2, figsize=Params.FigSize.induct)
            fig.subplots_adjust(wspace=0.5)
            fig.suptitle(inductionTitle)
            axes[0].title.set_text(name)
            axes[1].title.set_text('Phase delay $\\upphi$ ($^\circ$)')
            [ax.set_xscale('log') for ax in axes]
            [ax.set_xlabel('Salinity $w$ ($\si{g/kg}$)') for ax in axes]
            if Params.inductOtype == 'Tb':
                [ax.set_ylabel('Ice bottom temp.\ $T_b$ ($\si{K}$)') for ax in axes]
            elif Params.inductOtype == 'rho':
                [ax.set_ylabel('Silicate density $\\rho_\mathrm{sil}$ ($\si{kg/m^3}$)') for ax in axes]
            elif Params.inductOtype == 'phi':
                [ax.set_ylabel('Seafloor porosity $\phi_\mathrm{sil}$ (void frac)') for ax in axes]
                [ax.set_yscale('log') for ax in axes]
            else:
                raise ValueError(f'inductOtype {Params.inductOtype} not recognized.')

            zContours = [axes[0].contour(Induction.x, Induction.y, z[i, ...],
                                 colors=Params.Colors.Induction[T], linestyles=Params.LS_Induction[T],
                                 linewidths=Params.LW_Induction[T], levels=Params.cLevels[Induction.bodyname][T][fLabel])
                                 for i, T in enumerate(Induction.Texc_hr.keys())]
            phaseContours = [axes[1].contour(Induction.x, Induction.y, Induction.phase[i, ...],
                                 colors=Params.Colors.Induction[T], linestyles=Params.LS_Induction[T],
                                 linewidths=Params.LW_Induction[T], levels=Params.cLevels[Induction.bodyname][T]['phase'])
                                 for i, T in enumerate(Induction.Texc_hr.keys())]
            [axes[0].clabel(zContours[i], fmt=Params.cFmt[Induction.bodyname][T][fLabel],
                            fontsize=Params.cLabelSize, inline_spacing=Params.cLabelPad)
                            for i, T in enumerate(Induction.Texc_hr.keys())]
            [axes[1].clabel(phaseContours[i], fmt=Params.cFmt[Induction.bodyname][T]['phase'],
                            fontsize=Params.cLabelSize, inline_spacing=Params.cLabelPad)
                            for i, T in enumerate(Induction.Texc_hr.keys())]
            
            if Params.LEGEND:
                lines = np.array([contour.legend_elements()[0][0] for contour in phaseContours])
                labels = np.array([f'{T_h:.2f} h' for T_h in Induction.Texc_hr.values()])
                iSort = np.argsort(list(Induction.Texc_hr.values()))
                axes[1].legend(lines[iSort], labels[iSort], framealpha=Params.cLegendOpacity)

            fig.savefig(Params.FigureFiles.induct[fLabel], format=Params.figFormat)
            plt.close()

    log.getLogger().setLevel(saveLevel)

    return
