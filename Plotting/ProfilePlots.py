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
from MagneticInduction.configInduct import cLevels as magClevels, cFmt as magCfmt
from Plotting.configPlots import FigSize, Color, Style, FigMisc

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
    fig, axes = plt.subplots(1, 2, figsize=FigSize.vgrav)
    axes[0].plot('grav', 'radius', data=data)
    axes[0].set_xlabel('Gravity (m/s$^2$)')
    axes[0].set_ylabel('Radius (km)')

    axes[1].plot('pressure', 'radius', data = data)
    axes[1].set_xlabel('Pressure (GPa)')
    axes[1].set_ylabel('$r_\mathrm{' + Planet.name + '}$')

    fig.subplots_adjust(wspace=0.5)
    fig.suptitle(f'{PlanetList[0].name} gravity and pressure')
    fig.savefig(Params.FigureFiles.vgrav, format=FigMisc.figFormat, dpi=FigMisc.dpi)
    plt.close()
    return


def PlotHydrosphereProps(PlanetList, Params):
    # Generate canvas and add labels
    fig, axes = plt.subplots(1, 2, figsize=FigSize.vhydro)
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
                                            color=Color.ref,
                                            lw=Style.LW_ref,
                                            ls=Style.LS_ref[Planet.Ocean.comp])
                    if FigMisc.refsInLegend and i == 0: thisRef.set_label(wList)
                newRef[Planet.Ocean.comp] = False

    # Now plot all profiles together
    for Planet in PlanetList:

        # Plot density vs. pressure curve for hydrosphere
        axes[0].plot(Planet.P_MPa[:Planet.Steps.nHydro], Planet.rho_kgm3[:Planet.Steps.nHydro], label=Planet.label)
        # Plot thermal profile vs. depth in hydrosphere
        axes[1].plot(Planet.T_K[:Planet.Steps.nHydro], Planet.z_m[:Planet.Steps.nHydro]/1e3)

    if FigMisc.LEGEND:
        box1 = axes[0].get_position()
        fig.legend(loc=FigMisc.LegendPosition)
    fig.savefig(Params.FigureFiles.vhydro, format=FigMisc.figFormat, dpi=FigMisc.dpi)
    plt.close()
    return


def PlotCoreTradeoff(PlanetList, Params):
    data = {'Rsil': Planet.Sil.Rtrade_m/1000,
            'RFe': Planet.Core.Rtrade_m/1000}
    fig, axes = plt.subplots(1, 1, figsize=FigSize.vcore)
    axes.plot('Rsil', 'RFe', data = data)
    axes.set_xlabel('Iron core outer radius (km)')
    axes.set_ylabel('Silicate layer outer radius (km)')
    fig.suptitle(f'{PlanetList[0].name} with Fe core. $C/MR^2$: ${Planet.Bulk.Cmeasured:.3f}\pm{Planet.Bulk.Cuncertainty:.3f}' +
                 r'$; $w$: $0\,\mathrm{wt}\%$; $\rho_\mathrm{sil}$: $' + \
                 f'{Planet.Sil.rhoMean_kgm3:.0f}' + r'\,\mathrm{kg/m^3}$; $\rho_\mathrm{Fe}$: $' + \
                 f'{Planet.Core.rhoMean_kgm3:.0f}' + r'\,\mathrm{kg/m^3}$')
    fig.savefig(Params.FigureFiles.vcore, format=FigMisc.figFormat, dpi=FigMisc.dpi)
    plt.close()
    return


def PlotSilTradeoff(PlanetList, Params):
    data = {'Rsil': Planet.Sil.Rtrade_m/1000,
            'rhoSil': Planet.Sil.rhoTrade_kgm3}
    fig, axes = plt.subplots(1, 1, figsize=FigSize.vmant)
    axes.plot('rhoSil', 'Rsil', data = data)
    axes.set_xlabel('$\\rho_\mathrm{sil}$ (kg/m$^3$)')
    axes.set_ylabel('Silicate layer outer radius (km)')
    fig.suptitle(f'{PlanetList[0].name} no Fe core. $C/MR^2$: $0.346\pm0.005$; $W$')
    fig.savefig(Params.FigureFiles.vmant, format=FigMisc.figFormat, dpi=FigMisc.dpi)
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
            colors.append(Color.OceanTop)
        elif Planet.phase[layerPhase] == 1:
            colors.append(Color.IceI)
        elif Planet.phase[layerPhase] == 2:
            colors.append(Color.IceII)
        elif Planet.phase[layerPhase] == 3:
            colors.append(Color.IceIII)
        elif Planet.phase[layerPhase] == 5:
            colors.append(Color.IceV)
        elif Planet.phase[layerPhase] == 6:
            colors.append(Color.IceVI)
        elif Planet.phase[layerPhase] == Constants.phaseClath:
            colors.append(Color.Clath)
        elif Planet.phase[layerPhase] == Constants.phaseSil:
            colors.append(Color.Rock)
        elif Planet.phase[layerPhase] == Constants.phaseFe:
            colors.append(Color.Core)
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
    fig.savefig(Params.FigureFiles.vwedg, format=FigMisc.figFormat, dpi=FigMisc.dpi)

    plt.close()
    return


def PlotInductOgram(Induction, Params):
    log.warning('Temporarily quieting INFO and DEBUG messages due to a high number of current latex errors.')
    saveLevel = log.getLogger().getEffectiveLevel()
    log.getLogger().setLevel(log.WARN)

    # Get all common labels and data for zipping
    zData = [Induction.Amp, Induction.Bix_nT, Induction.Biy_nT, Induction.Biz_nT]
    plotTitles = ['Amplitude $A$', '$B_x$ component', '$B_y$ component', '$B_z$ component']
    fLabels = ['Amp', 'Bx', 'By', 'Bz']
    inductionTitle = f'\\textbf{{{Induction.bodyname} induction response}}'
    phaseTitle = 'Phase delay $\\upphi$ ($^\circ$)'
    sigLabel = 'Mean conductivity $\overline{\sigma}$ ($\si{S/m}$)'
    Dlabel = 'Ocean thickness $D$ ($\si{km}$)'
    sigLims = [10**Params.Induct.sigmaMin[Induction.bodyname], 10**Params.Induct.sigmaMax[Induction.bodyname]]
    Dlims = [10**Params.Induct.Dmin[Induction.bodyname], 10**Params.Induct.Dmax[Induction.bodyname]]
    wLabel = 'Salinity $w$ ($\si{g/kg}$)'
    TbLabel = 'Ice bottom temp.\ $T_b$ ($\si{K}$)'
    rhoLabel = 'Silicate density $\\rho_\mathrm{sil}$ ($\si{kg/m^3}$)'
    phiLabel = 'Seafloor porosity $\phi_\mathrm{sil}$ (void frac)'
    legendLabels = np.array([f'{T_h:.2f} h' for T_h in Induction.Texc_hr.values()])
    iSort = np.argsort(list(Induction.Texc_hr.values()))

    if Params.Induct.inductOtype != 'sigma':
        x = Induction.w_ppt
        if Params.Induct.inductOtype == 'Tb':
            yLabel = TbLabel
            yScale = 'linear'
            y = Induction.Tb_K
        elif Params.Induct.inductOtype == 'rho':
            yLabel = rhoLabel
            yScale = 'linear'
            y = Induction.rhoSilMean_kgm3
        elif Params.Induct.inductOtype == 'phi':
            yLabel = phiLabel
            yScale = 'log'
            y = Induction.phiSilMax_frac
        else:
            raise ValueError(f'inductOtype {Params.Induct.inductOtype} not recognized.')

    if Params.COMBINE_BCOMPS:
        # Plot B components all together with phase. Amplitude is still separate
        # Generate canvas and add labels
        fig, axes = plt.subplots(2, 2, figsize=FigSize.inductCombo)
        allAxes = axes.flatten()
        fig.suptitle(inductionTitle)
        # Only label the bottom-left sides of axes
        [ax.set_xlabel(sigLabel) for ax in (axes[1,0], axes[1,1])]
        [ax.set_ylabel(Dlabel) for ax in (axes[0,0], axes[1,0])]
        [ax.set_xlim(sigLims) for ax in allAxes]
        [ax.set_ylim(Dlims) for ax in allAxes]
        [ax.set_xscale('log') for ax in allAxes]
        [ax.set_yscale('log') for ax in allAxes]
        coords = {'Bx': (0,0), 'By': (0,1), 'Bz': (1,0), 'phase': (1,1)}
        comboData = [Induction.Bix_nT, Induction.Biy_nT, Induction.Biz_nT, Induction.phase]
        comboTitles = np.append(plotTitles[1:], phaseTitle)
        comboLabels = list(coords.keys())

        for z, name, fLabel in zip(comboData, comboTitles, comboLabels):
            ax = axes[coords[fLabel]]
            ax.title.set_text(name)
            zContours = [ax.contour(Induction.sigmaMean_Sm, Induction.D_km, z[i, ...],
                           colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                           linewidths=Style.LW_Induction[T], levels=magClevels[Induction.bodyname][T][fLabel])
                           for i, T in enumerate(Induction.Texc_hr.keys())]
            if Params.Induct.inductOtype == 'sigma':
                [ax.clabel(zContours[i], fmt=magCfmt[Induction.bodyname][T][fLabel],
                           fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                           for i, T in enumerate(Induction.Texc_hr.keys())]

        if FigMisc.LEGEND:
            lines = np.array([contour.legend_elements()[0][0] for contour in zContours])
            axes[1,1].legend(lines[iSort], legendLabels[iSort], framealpha=FigMisc.cLegendOpacity)

        fig.savefig(Params.FigureFiles.sigma['Bcomps'], format=FigMisc.figFormat, dpi=FigMisc.dpi)
        plt.close()

        if Params.Induct.inductOtype != 'sigma':
            fig, axes = plt.subplots(2, 2, figsize=FigSize.inductCombo)
            allAxes = axes.flatten()
            fig.suptitle(inductionTitle)
            # Only label the bottom-left sides of axes
            [ax.set_xlabel(wLabel) for ax in (axes[1,0], axes[1,1])]
            [ax.set_ylabel(yLabel) for ax in (axes[0,0], axes[1,0])]
            [ax.set_xscale('log') for ax in allAxes]
            [ax.set_yscale(yScale) for ax in allAxes]

            for z, name, fLabel in zip(comboData, comboTitles, comboLabels):
                ax = axes[coords[fLabel]]
                ax.title.set_text(name)
                zContours = [ax.contour(x, y, z[i, ...],
                                        colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                                        linewidths=Style.LW_Induction[T],
                                        levels=magClevels[Induction.bodyname][T][fLabel])
                             for i, T in enumerate(Induction.Texc_hr.keys())]
                [ax.clabel(zContours[i], fmt=magCfmt[Induction.bodyname][T][fLabel],
                           fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                           for i, T in enumerate(Induction.Texc_hr.keys())]

            if FigMisc.LEGEND:
                lines = np.array([contour.legend_elements()[0][0] for contour in zContours])
                axes[1,1].legend(lines[iSort], legendLabels[iSort], framealpha=FigMisc.cLegendOpacity)

            fig.savefig(Params.FigureFiles.sigma['Bcomps'], format=FigMisc.figFormat, dpi=FigMisc.dpi)
            plt.close()

        # Set lists to just contain Amplitude now to reuse the remaining routines for that plot
        zData = [zData[0]]
        plotTitles = [plotTitles[0]]
        fLabels = [fLabels[0]]

    # Plot each component separately alongside phase
    for z, name, fLabel in zip(zData, plotTitles, fLabels):

        # Generate canvas and add labels
        fig, axes = plt.subplots(1, 2, figsize=FigSize.induct)
        fig.subplots_adjust(wspace=0.5)
        fig.suptitle(inductionTitle)
        axes[0].title.set_text(name)
        axes[1].title.set_text(phaseTitle)
        [ax.set_xlabel(sigLabel) for ax in axes]
        [ax.set_ylabel(Dlabel) for ax in axes]
        [ax.set_xlim(sigLims) for ax in axes]
        [ax.set_ylim(Dlims) for ax in axes]
        [ax.set_xscale('log') for ax in axes]
        [ax.set_yscale('log') for ax in axes]

        zContours = [axes[0].contour(Induction.sigmaMean_Sm, Induction.D_km, z[i, ...],
                         colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                         linewidths=Style.LW_Induction[T], levels=magClevels[Induction.bodyname][T][fLabel])
                         for i, T in enumerate(Induction.Texc_hr.keys())]
        phaseContours = [axes[1].contour(Induction.sigmaMean_Sm, Induction.D_km, Induction.phase[i, ...],
                         colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                         linewidths=Style.LW_Induction[T], levels=magClevels[Induction.bodyname][T]['phase'])
                         for i, T in enumerate(Induction.Texc_hr.keys())]
        if Params.Induct.inductOtype == 'sigma':
            [axes[0].clabel(zContours[i], fmt=magCfmt[Induction.bodyname][T][fLabel],
                            fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                            for i, T in enumerate(Induction.Texc_hr.keys())]
            [axes[1].clabel(phaseContours[i], fmt=magCfmt[Induction.bodyname][T]['phase'],
                            fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                            for i, T in enumerate(Induction.Texc_hr.keys())]
            fNameSigma = Params.FigureFiles.sigmaOnly[fLabel]
        else:
            fNameSigma = Params.FigureFiles.sigma[fLabel]

        if FigMisc.LEGEND:
            lines = np.array([contour.legend_elements()[0][0] for contour in phaseContours])
            axes[1].legend(lines[iSort], legendLabels[iSort], framealpha=FigMisc.cLegendOpacity)

        fig.savefig(fNameSigma, format=FigMisc.figFormat, dpi=FigMisc.dpi)
        plt.close()

        if Params.Induct.inductOtype != 'sigma':
            fig, axes = plt.subplots(1, 2, figsize=FigSize.induct)
            fig.subplots_adjust(wspace=0.5)
            fig.suptitle(inductionTitle)
            axes[0].title.set_text(name)
            axes[1].title.set_text(phaseTitle)
            [ax.set_xscale('log') for ax in axes]
            [ax.set_xlabel(wLabel) for ax in axes]
            if Params.Induct.inductOtype == 'Tb':
                [ax.set_ylabel(TbLabel) for ax in axes]
            elif Params.Induct.inductOtype == 'rho':
                [ax.set_ylabel(rhoLabel) for ax in axes]
            elif Params.Induct.inductOtype == 'phi':
                [ax.set_ylabel(phiLabel) for ax in axes]
                [ax.set_yscale('log') for ax in axes]
            else:
                raise ValueError(f'inductOtype {Params.Induct.inductOtype} not recognized.')

            zContours = [axes[0].contour(x, y, z[i, ...],
                             colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                             linewidths=Style.LW_Induction[T], levels=magClevels[Induction.bodyname][T][fLabel])
                             for i, T in enumerate(Induction.Texc_hr.keys())]
            phaseContours = [axes[1].contour(x, y, Induction.phase[i, ...],
                             colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                             linewidths=Style.LW_Induction[T], levels=magClevels[Induction.bodyname][T]['phase'])
                             for i, T in enumerate(Induction.Texc_hr.keys())]
            [axes[0].clabel(zContours[i], fmt=magCfmt[Induction.bodyname][T][fLabel],
                            fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                            for i, T in enumerate(Induction.Texc_hr.keys())]
            [axes[1].clabel(phaseContours[i], fmt=magCfmt[Induction.bodyname][T]['phase'],
                            fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                            for i, T in enumerate(Induction.Texc_hr.keys())]

            if FigMisc.LEGEND:
                lines = np.array([contour.legend_elements()[0][0] for contour in phaseContours])
                axes[1].legend(lines[iSort], legendLabels[iSort], framealpha=FigMisc.cLegendOpacity)

            fig.savefig(Params.FigureFiles.induct[fLabel], format=FigMisc.figFormat, dpi=FigMisc.dpi)
            plt.close()

    log.getLogger().setLevel(saveLevel)

    return
