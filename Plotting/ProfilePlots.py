import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from matplotlib.patches import Wedge
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import logging as log
from pathlib import Path
from scipy.interpolate import interp1d
from Utilities.SetupInit import SetupFilenames
from Utilities.defineStructs import Constants
from Thermodynamics.RefProfiles.RefProfiles import CalcRefProfiles, ReloadRefProfiles
from Plotting.configPlots import FigSize, Color, Style, FigMisc, FigLbl
from MagneticInduction.configInduct import InductParams as IndParams

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

    if Params.PLOT_GRAVITY: PlotGravPres(PlanetList, Params)
    if Params.PLOT_HYDROSPHERE and not PlanetList[0].Do.NO_H2O: PlotHydrosphereProps(PlanetList, Params)
    if Params.PLOT_TRADEOFF:
        if Planet.Do.Fe_CORE: PlotCoreTradeoff(PlanetList, Params)
        else: PlotSilTradeoff(PlanetList, Params)
    if Params.PLOT_WEDGE: PlotWedge(PlanetList, Params)

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
            if newRef[Planet.Ocean.comp] and Planet.Ocean.comp != 'none':
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
        # This is a hydrosphere-only plot, so skip waterless bodies
        if Planet.Ocean.comp != 'none':
            # Plot density vs. pressure curve for hydrosphere
            axes[0].plot(Planet.P_MPa[:Planet.Steps.nHydro], Planet.rho_kgm3[:Planet.Steps.nHydro], label=Planet.label)
            # Plot thermal profile vs. depth in hydrosphere
            axes[1].plot(Planet.T_K[:Planet.Steps.nHydro], Planet.z_m[:Planet.Steps.nHydro]/1e3)

    if FigMisc.LEGEND:
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


def PlotInductOgramPhaseSpace(InductionList, Params):
    """ For plotting points showing the various models used in making
        inductogram plots.
    """

    if InductionList[0].SINGLE_COMP:
        FigLbl.singleComp(InductionList[0].comps[0])
    FigLbl.setInduction(InductionList[0].bodyname, Params.Induct, InductionList[0].Texc_hr.values())

    sigma_Sm, D_km, ptColors = (np.empty_like(InductionList) for _ in range(3))
    for i,Induction in enumerate(InductionList):
        sigma_Sm[i] = Induction.sigmaMean_Sm.flatten()
        D_km[i] = Induction.D_km.flatten()
        if Params.Induct.inductOtype == 'sigma':
            # In this case, we likely don't have salinity and ocean temp information
            # so we need to set the colormap to use the info we do have
            sigmaNorm = sigma_Sm[i] / 10**Params.Induct.sigmaMax[Induction.bodyname]
            Dnorm = D_km[i] / np.max(D_km)
            ptColors[i] = Color.OceanCmap(Induction.compsList, sigmaNorm, Dnorm)
        else:
            w_ppt = Induction.x.flatten()
            Tmean_K = Induction.Tmean_K.flatten()
            if FigMisc.NORMALIZED_SALINITIES:
                wMax_ppt = np.array([Color.saturation[comp] for comp in Induction.compsList])
                w_normFrac = w_ppt / wMax_ppt
            else:
                w_normFrac = interp1d([np.min(w_ppt), np.max(w_ppt)], [0.0, 1.0])(w_ppt)
            if Params.Induct.colorType == 'Tmean':
                if FigMisc.NORMALIZED_TEMPERATURES:
                    Tmean_normFrac = Color.GetNormT(Tmean_K)
                else:
                    Tmean_normFrac = interp1d([np.min(Tmean_K), np.max(Tmean_K)], [0.0, 1.0])(Tmean_K)
                ptColors[i] = Color.OceanCmap(Induction.compsList, w_normFrac, Tmean_normFrac)
            elif Params.Induct.colorType == 'zb':
                zb_km = Induction.zb_km.flatten()
                zb_normFrac = interp1d([np.min(zb_km), np.max(zb_km)], [0.0, 1.0])(zb_km)
                ptColors[i] = Color.OceanCmap(Induction.compsList, w_normFrac, zb_normFrac)
            else:
                raise ValueError(f'Inductogram colortype {Params.Induct.colorType} not recognized.')

    if Params.Induct.inductOtype == 'sigma':
        fig, ax = plt.subplots(1, 1, figsize=FigSize.phaseSpaceSolo)
        axes = [ax]
        cbarUnits = InductionList[0].zb_km.flatten()
        cbarLabel = FigLbl.iceThickLbl
    else:
        w_ppt = InductionList[0].x.flatten()
        yFlat = InductionList[0].y.flatten()
        if Params.Induct.colorType == 'Tmean':
            cbarUnits = InductionList[0].Tmean_K.flatten()
            cbarLabel = FigLbl.oceanTempLbl
        elif Params.Induct.colorType == 'zb':
            cbarUnits = InductionList[0].zb_km.flatten()
            cbarLabel = FigLbl.iceThickLbl

        fig, axes = plt.subplots(1, 2, figsize=FigSize.phaseSpaceCombo)
        axes[1].set_xlabel(FigLbl.wLabel)
        axes[1].set_ylabel(FigLbl.yLabelInduct)
        axes[1].set_xscale(FigLbl.wScale)
        axes[1].set_yscale(FigLbl.yScaleInduct)
        axes[1].scatter(w_ppt, yFlat, s=Style.MW_Induction,
                        marker=Style.MS_Induction, c=ptColors[0])

    fig.suptitle(FigLbl.phaseSpaceTitle)
    axes[0].set_xlabel(FigLbl.sigLabel)
    axes[0].set_ylabel(FigLbl.Dlabel)
    axes[0].set_xlim(FigLbl.sigLims)
    axes[0].set_ylim(FigLbl.Dlims)
    axes[0].set_xscale(FigLbl.sigScale)
    axes[0].set_yscale(FigLbl.Dscale)

    pts = {}
    cbar = {}
    if Params.Induct.inductOtype == 'sigma':
        comps = ['Ice']
        divider = make_axes_locatable(axes[0])
        pts[comps[0]] = axes[0].scatter(sigma_Sm[0], D_km[0], s=Style.MW_Induction,
                              marker=Style.MS_Induction, c=ptColors[0])
        cbarAx = divider.new_horizontal(size=FigMisc.cbarSize, pad=FigMisc.cbarPad)
        cbar[comps[0]] = mpl.colorbar.ColorbarBase(cbarAx, cmap=Color.cmap[comps[0]],
                                               values=np.linspace(np.min(cbarUnits), np.max(cbarUnits), FigMisc.nCbarPts),
                                               format=FigMisc.cbarFmt, orientation='vertical')
        fig.add_axes(cbarAx)
    else:
        divider = make_axes_locatable(axes[1])
        extraPad = 0
        comps = np.unique(InductionList[0].comps)
        for comp in comps:
            thisComp = InductionList[0].compsList == comp
            pts[comp] = axes[0].scatter(sigma_Sm[0][thisComp], D_km[0][thisComp], s=Style.MW_Induction,
                            marker=Style.MS_Induction, c=ptColors[0][thisComp])
            cbarAx = divider.new_horizontal(size=FigMisc.cbarSize, pad=FigMisc.cbarPad + extraPad)
            extraPad = FigMisc.extraPad
            cbar[comp] = mpl.colorbar.ColorbarBase(cbarAx, cmap=Color.cmap[comp],
                                             values=np.linspace(np.min(cbarUnits[thisComp]), np.max(cbarUnits[thisComp]), FigMisc.nCbarPts),
                                             format=FigMisc.cbarFmt, orientation='vertical')
            fig.add_axes(cbarAx)
            cbarAx.set_title(f'\ce{{{comp}}}')

    cbar[comps[-1]].set_label(cbarLabel, size=12)
    fig.savefig(Params.FigureFiles.phaseSpace, format=FigMisc.figFormat, dpi=FigMisc.dpi)
    log.debug(f'Plot saved to file: {Params.FigureFiles.phaseSpace}')
    plt.close()

    # Plot combination
    if Params.COMPARE and np.size(InductionList) > 1 and Params.Induct.inductOtype != 'sigma':
        comps = np.unique(np.append([],[Induction.comps for Induction in InductionList]))
        nComps = np.size(comps)
        figWidth = FigSize.phaseSpaceSolo[0] + nComps * FigMisc.cbarSpace
        fig, ax = plt.subplots(1, 1, figsize=(figWidth, FigSize.phaseSpaceSolo[1]))

        fig.suptitle(FigLbl.phaseSpaceTitle)
        ax.set_xlabel(FigLbl.sigLabel)
        ax.set_ylabel(FigLbl.Dlabel)
        ax.set_xlim(FigLbl.sigLims)
        ax.set_ylim(FigLbl.Dlims)
        ax.set_xscale(FigLbl.sigScale)
        ax.set_yscale(FigLbl.Dscale)

        divider = make_axes_locatable(ax)
        extraPad = 0
        comboCompsList = np.concatenate(tuple(Induction.compsList for Induction in InductionList))
        comboSigma_Sm = np.concatenate(tuple(sigmai for sigmai in sigma_Sm))
        comboD_km = np.concatenate(tuple(Di for Di in D_km))
        comboColors = np.concatenate(tuple(ptColi for ptColi in ptColors))
        if Params.Induct.colorType == 'Tmean':
            comboCbarUnits = np.concatenate(tuple(Induction.Tmean_K.flatten() for Induction in InductionList))
        elif Params.Induct.colorType == 'zb':
            comboCbarUnits = np.concatenate(tuple(Induction.zb_km.flatten() for Induction in InductionList))
        pts = {}
        for comp in comps:
            thisComp = comboCompsList == comp
            pts[comp] = ax.scatter(comboSigma_Sm[thisComp], comboD_km[thisComp], s=Style.MW_Induction,
                                        marker=Style.MS_Induction, c=comboColors[thisComp])
            cbarAx = divider.new_horizontal(size=FigMisc.cbarSize, pad=FigMisc.cbarPad + extraPad)
            extraPad = FigMisc.extraPad
            cbar = mpl.colorbar.ColorbarBase(cbarAx, cmap=Color.cmap[comp],
                                             values=np.linspace(np.min(comboCbarUnits[thisComp]), np.max(comboCbarUnits[thisComp]), FigMisc.nCbarPts),
                                             format=FigMisc.cbarFmt, orientation='vertical')
            fig.add_axes(cbarAx)
            cbarAx.set_title(f'\ce{{{comp}}}')

        cbar.set_label(cbarLabel, size=12)
        fig.savefig(Params.FigureFiles.phaseSpaceCombo, format=FigMisc.figFormat, dpi=FigMisc.dpi)
        log.debug(f'Plot saved to file: {Params.FigureFiles.phaseSpaceCombo}')
        plt.close()

    return


def PlotInductOgram(Induction, Params):
    """ Plot contours showing magnetic induction responses for an array of models
    """

    # Get all common labels and data for zipping
    zData = [Induction.Amp, np.abs(Induction.Bix_nT),
             np.abs(Induction.Biy_nT), np.abs(Induction.Biz_nT)]
    if Induction.SINGLE_COMP:
        FigLbl.singleComp(Induction.comps[0])
    FigLbl.setInduction(Induction.bodyname, Params.Induct, Induction.Texc_hr.values())
    iSort = np.argsort(list(Induction.Texc_hr.values()))

    if Params.COMBINE_BCOMPS:
        # Plot B components all together with phase. Amplitude is still separate
        # Generate canvas and add labels
        fig, axes = plt.subplots(2, 2, figsize=FigSize.inductCombo)
        allAxes = axes.flatten()
        fig.suptitle(FigLbl.inductionTitle)
        # Only label the bottom-left sides of axes
        [ax.set_xlabel(FigLbl.sigLabel) for ax in (axes[1,0], axes[1,1])]
        [ax.set_ylabel(FigLbl.Dlabel) for ax in (axes[0,0], axes[1,0])]
        [ax.set_xlim(FigLbl.sigLims) for ax in allAxes]
        [ax.set_ylim(FigLbl.Dlims) for ax in allAxes]
        [ax.set_xscale(FigLbl.sigScale) for ax in allAxes]
        [ax.set_yscale(FigLbl.Dscale) for ax in allAxes]
        coords = {'Bx': (0,0), 'By': (0,1), 'Bz': (1,0), 'phase': (1,1)}
        comboData = [np.abs(Induction.Bix_nT), np.abs(Induction.Biy_nT),
                     np.abs(Induction.Biz_nT), Induction.phase]
        comboTitles = np.append(FigLbl.plotTitles[1:], FigLbl.phaseTitle)
        comboLabels = list(coords.keys())

        for z, name, fLabel in zip(comboData, comboTitles, comboLabels):
            ax = axes[coords[fLabel]]
            ax.title.set_text(name)
            zContours = [ax.contour(Induction.sigmaMean_Sm, Induction.D_km, z[i, ...],
                           colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                           linewidths=Style.LW_Induction[T], levels=IndParams.GetClevels(fLabel, T))
                           for i, T in enumerate(Induction.Texc_hr.keys())]
            if Params.Induct.inductOtype == 'sigma':
                [ax.clabel(zContours[i], fmt=IndParams.GetCfmt(fLabel, T),
                           fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                           for i, T in enumerate(Induction.Texc_hr.keys())]

        if FigMisc.LEGEND:
            lines = np.array([contour.legend_elements()[0][0] for contour in zContours])
            axes[1,1].legend(lines[iSort], FigLbl.legendTexc[iSort], framealpha=FigMisc.cLegendOpacity)

        if Params.Induct.inductOtype == 'sigma':
            fNameSigma = Params.FigureFiles.sigmaOnly['Bcomps']
        else:
            fNameSigma = Params.FigureFiles.sigma['Bcomps']
        fig.savefig(fNameSigma, format=FigMisc.figFormat, dpi=FigMisc.dpi)
        log.debug(f'Plot saved to file: {fNameSigma}')
        plt.close()

        if Params.Induct.inductOtype != 'sigma':
            fig, axes = plt.subplots(2, 2, figsize=FigSize.inductCombo)
            allAxes = axes.flatten()
            fig.suptitle(FigLbl.inductionTitle)
            # Only label the bottom-left sides of axes
            [ax.set_xlabel(FigLbl.wLabel) for ax in (axes[1,0], axes[1,1])]
            [ax.set_ylabel(FigLbl.yLabelInduct) for ax in (axes[0,0], axes[1,0])]
            [ax.set_xscale(FigLbl.wScale) for ax in allAxes]
            [ax.set_yscale(FigLbl.yScaleInduct) for ax in allAxes]

            for z, name, fLabel in zip(comboData, comboTitles, comboLabels):
                ax = axes[coords[fLabel]]
                ax.title.set_text(name)
                zContours = [ax.contour(Induction.x, Induction.y, z[i, ...],
                                        colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                                        linewidths=Style.LW_Induction[T],
                                        levels=IndParams.GetClevels(fLabel, T))
                             for i, T in enumerate(Induction.Texc_hr.keys())]
                [ax.clabel(zContours[i], fmt=IndParams.GetCfmt(fLabel, T),
                           fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                           for i, T in enumerate(Induction.Texc_hr.keys())]

            if FigMisc.LEGEND:
                lines = np.array([contour.legend_elements()[0][0] for contour in zContours])
                axes[1,1].legend(lines[iSort], FigLbl.legendTexc[iSort], framealpha=FigMisc.cLegendOpacity)

            fig.savefig(Params.FigureFiles.induct['Bcomps'], format=FigMisc.figFormat, dpi=FigMisc.dpi)
            log.debug(f'Plot saved to file: {Params.FigureFiles.induct["Bcomps"]}')
            plt.close()

            # Also plot a comparison of Bx, which is usually the strongest oscillation
            compChoice = 'Bx'
            fig, axes = plt.subplots(2, 2, figsize=FigSize.inductCombo)
            fig.subplots_adjust(wspace=0.25, hspace=0.35)
            allAxes = axes.flatten()
            fig.suptitle(FigLbl.inductCompareTitle)
            # Label all axes for clarity
            [ax.set_xlabel(FigLbl.wLabel) for ax in axes[0,:]]
            [ax.set_ylabel(FigLbl.yLabelInduct) for ax in axes[0,:]]
            [ax.set_xlabel(FigLbl.sigLabel) for ax in axes[1,:]]
            [ax.set_ylabel(FigLbl.Dlabel) for ax in axes[1,:]]
            [ax.set_xscale(FigLbl.sigScale) for ax in allAxes]
            [ax.set_yscale(FigLbl.Dscale) for ax in axes[1,:]]
            [ax.set_yscale(FigLbl.yScaleInduct) for ax in axes[0,:]]
            [ax.set_xlim(FigLbl.sigLims) for ax in axes[1,:]]
            [ax.set_ylim(FigLbl.Dlims) for ax in axes[1,:]]

            axes[0,0].title.set_text(comboTitles[0])
            axes[1,0].title.set_text(comboTitles[0])
            axes[0,1].title.set_text(comboTitles[-1])
            axes[1,1].title.set_text(comboTitles[-1])
            zContours = [axes[0,0].contour(Induction.x, Induction.y, comboData[0][i, ...],
                                    colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                                    linewidths=Style.LW_Induction[T],
                                    levels=IndParams.GetClevels(comboLabels[0], T))
                         for i, T in enumerate(Induction.Texc_hr.keys())]
            phaseContours = [axes[0,1].contour(Induction.x, Induction.y, comboData[-1][i, ...],
                                    colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                                    linewidths=Style.LW_Induction[T],
                                    levels=IndParams.GetClevels(comboLabels[-1], T))
                         for i, T in enumerate(Induction.Texc_hr.keys())]
            [axes[1,0].contour(Induction.sigmaMean_Sm, Induction.D_km, comboData[0][i, ...],
                                    colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                                    linewidths=Style.LW_Induction[T],
                                    levels=IndParams.GetClevels(comboLabels[0], T))
                         for i, T in enumerate(Induction.Texc_hr.keys())]
            [axes[1,1].contour(Induction.sigmaMean_Sm, Induction.D_km, comboData[-1][i, ...],
                                    colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                                    linewidths=Style.LW_Induction[T],
                                    levels=IndParams.GetClevels(comboLabels[-1], T))
                         for i, T in enumerate(Induction.Texc_hr.keys())]
            [axes[0,0].clabel(zContours[i], fmt=IndParams.GetCfmt(comboLabels[0], T),
                       fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                       for i, T in enumerate(Induction.Texc_hr.keys())]
            [axes[0,1].clabel(phaseContours[i], fmt=IndParams.GetCfmt(comboLabels[-1], T),
                       fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                       for i, T in enumerate(Induction.Texc_hr.keys())]

            if FigMisc.LEGEND:
                lines = np.array([contour.legend_elements()[0][0] for contour in zContours])
                axes[1,1].legend(lines[iSort], FigLbl.legendTexc[iSort], framealpha=FigMisc.cLegendOpacity)
            fig.savefig(Params.FigureFiles.inductCompare[compChoice], format=FigMisc.figFormat, dpi=FigMisc.dpi)
            log.debug(f'Plot saved to file: {Params.FigureFiles.inductCompare[compChoice]}')
            plt.close()

        # Set lists to just contain Amplitude now to reuse the remaining routines for that plot
        zData = [zData[0]]
        FigLbl.plotTitles = [FigLbl.plotTitles[0]]
        FigLbl.fLabels = [FigLbl.fLabels[0]]

    # Plot each component separately alongside phase
    for z, name, fLabel in zip(zData, FigLbl.plotTitles, FigLbl.fLabels):

        # Generate canvas and add labels
        fig, axes = plt.subplots(1, 2, figsize=FigSize.induct)
        fig.subplots_adjust(wspace=0.5)
        fig.suptitle(FigLbl.inductionTitle)
        axes[0].title.set_text(name)
        axes[1].title.set_text(FigLbl.phaseTitle)
        [ax.set_xlabel(FigLbl.sigLabel) for ax in axes]
        [ax.set_ylabel(FigLbl.Dlabel) for ax in axes]
        [ax.set_xlim(FigLbl.sigLims) for ax in axes]
        [ax.set_ylim(FigLbl.Dlims) for ax in axes]
        [ax.set_xscale(FigLbl.sigScale) for ax in axes]
        [ax.set_yscale(FigLbl.Dscale) for ax in axes]

        zContours = [axes[0].contour(Induction.sigmaMean_Sm, Induction.D_km, z[i, ...],
                         colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                         linewidths=Style.LW_Induction[T], levels=IndParams.GetClevels(fLabel, T))
                         for i, T in enumerate(Induction.Texc_hr.keys())]
        phaseContours = [axes[1].contour(Induction.sigmaMean_Sm, Induction.D_km, Induction.phase[i, ...],
                         colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                         linewidths=Style.LW_Induction[T], levels=IndParams.GetClevels('phase', T))
                         for i, T in enumerate(Induction.Texc_hr.keys())]
        if Params.Induct.inductOtype == 'sigma':
            [axes[0].clabel(zContours[i], fmt=IndParams.GetCfmt(fLabel, T),
                            fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                            for i, T in enumerate(Induction.Texc_hr.keys())]
            [axes[1].clabel(phaseContours[i], fmt=IndParams.GetCfmt('phase', T),
                            fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                            for i, T in enumerate(Induction.Texc_hr.keys())]
            fNameSigma = Params.FigureFiles.sigmaOnly[fLabel]
        else:
            fNameSigma = Params.FigureFiles.sigma[fLabel]

        if FigMisc.LEGEND:
            lines = np.array([contour.legend_elements()[0][0] for contour in phaseContours])
            axes[1].legend(lines[iSort], FigLbl.legendTexc[iSort], framealpha=FigMisc.cLegendOpacity)

        fig.savefig(fNameSigma, format=FigMisc.figFormat, dpi=FigMisc.dpi)
        log.debug(f'Plot saved to file: {fNameSigma}')
        plt.close()

        if Params.Induct.inductOtype != 'sigma':
            fig, axes = plt.subplots(1, 2, figsize=FigSize.induct)
            fig.subplots_adjust(wspace=0.5)
            fig.suptitle(FigLbl.inductionTitle)
            axes[0].title.set_text(name)
            axes[1].title.set_text(FigLbl.phaseTitle)
            [ax.set_xscale(FigLbl.wScale) for ax in axes]
            [ax.set_xlabel(FigLbl.wLabel) for ax in axes]
            [ax.set_ylabel(FigLbl.yLabelInduct) for ax in axes]
            [ax.set_yscale(FigLbl.yScaleInduct) for ax in axes]

            zContours = [axes[0].contour(Induction.x, Induction.y, z[i, ...],
                             colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                             linewidths=Style.LW_Induction[T], levels=IndParams.GetClevels(fLabel, T))
                             for i, T in enumerate(Induction.Texc_hr.keys())]
            phaseContours = [axes[1].contour(Induction.x, Induction.y, Induction.phase[i, ...],
                             colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                             linewidths=Style.LW_Induction[T], levels=IndParams.GetClevels('phase', T))
                             for i, T in enumerate(Induction.Texc_hr.keys())]
            [axes[0].clabel(zContours[i], fmt=IndParams.GetCfmt(fLabel, T),
                            fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                            for i, T in enumerate(Induction.Texc_hr.keys())]
            [axes[1].clabel(phaseContours[i], fmt=IndParams.GetCfmt('phase', T),
                            fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                            for i, T in enumerate(Induction.Texc_hr.keys())]

            if FigMisc.LEGEND:
                lines = np.array([contour.legend_elements()[0][0] for contour in phaseContours])
                axes[1].legend(lines[iSort], FigLbl.legendTexc[iSort], framealpha=FigMisc.cLegendOpacity)

            fig.savefig(Params.FigureFiles.induct[fLabel], format=FigMisc.figFormat, dpi=FigMisc.dpi)
            log.debug(f'Plot saved to file: {Params.FigureFiles.induct[fLabel]}')
            plt.close()

    return
