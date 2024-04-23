import numpy as np
import logging
from copy import deepcopy
from collections.abc import Iterable
import matplotlib.pyplot as plt
import matplotlib.colorbar as mcbar
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp1d
from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc
from PlanetProfile.Utilities.defineStructs import xyzComps, vecComps
from PlanetProfile.MagneticInduction.Moments import Excitations
from MoonMag.asymmetry_funcs import getMagSurf as GetMagSurf

# Assign logger
log = logging.getLogger('PlanetProfile')

def GenerateMagPlots(PlanetList, Params):
    
    # Catch if we get here and induction calcs have not been done
    if np.all([Planet.Magnetic.Binm_nT is not None for Planet in PlanetList]):

        if np.any([isinstance(Planet.Magnetic.Benm_nT, dict) for Planet in PlanetList]):
            log.warning('Multiple spacecraft eras are being modeled, which will break many ' +
                        'magnetic field plotting functions. General plots will be skipped.')

        else:
            # Remove latex styling from legend labels if Latex is not installed
            if not FigMisc.TEX_INSTALLED:
                for Planet in PlanetList:
                    Planet.label = FigLbl.StripLatexFromString(Planet.label)

            if Params.PLOT_BDIP:
                PlotComplexBdip(PlanetList, Params)
            if (Params.PLOT_MAG_SPECTRUM_COMBO or Params.PLOT_MAG_SPECTRUM) and np.any(
                    [Planet.Magnetic.FT_LOADED for Planet in PlanetList]):
                if Params.PLOT_MAG_SPECTRUM_COMBO:
                    PlotMagSpectrumInduced(PlanetList, Params)
                if Params.PLOT_MAG_SPECTRUM:
                    PlotMagSpectrum(PlanetList, Params)
            if Params.PLOT_BSURF:
                PlotMagSurface(PlanetList, Params)
            if Params.PLOT_ASYM and Params.CALC_ASYM:
                PlotAsym(PlanetList, Params)

            elif np.any([Planet.Magnetic.Binm_nT is not None for Planet in PlanetList]):
                log.warning('Induction calculations have only been completed for some bodies. Asymmetry contour plotting will be skipped.')
            else:
                log.warning('Induction calculations have not been completed for any bodies. Asymmetry contour plotting will be skipped. ' +
                            'Set Params.SKIP_INDUCTION = True in configPP.py to skip induction calculations.')

    return


def PlotInductOgramPhaseSpace(InductionList, Params):
    """ For plotting points showing the various models used in making
        inductogram plots.
    """

    if InductionList[0].SINGLE_COMP:
        FigLbl.singleComp(InductionList[0].comps[0])
    FigLbl.SetInduction(InductionList[0].bodyname, Params.Induct, InductionList[0].Texc_hr.values())
    if not FigMisc.TEX_INSTALLED:
        FigLbl.StripLatex()

    sigma_Sm, D_km, ptColors, w_ppt, yFlat, iValid, iValidFlat \
        = (np.empty_like(InductionList) for _ in range(7))
    for i, Induction in enumerate(InductionList):
        iValid[i] = np.where(np.isfinite(Induction.Amp[0,...]))
        iValidFlat[i] = np.where(np.isfinite(Induction.Amp[0,...].flatten()))[0]
        sigma_Sm[i] = Induction.sigmaMean_Sm[iValid[i]].flatten()
        D_km[i] = Induction.D_km[iValid[i]].flatten()
        if Params.Induct.inductOtype == 'sigma':
            # In this case, we likely don't have salinity and ocean temp information
            # so we need to set the colormap to use the info we do have
            sigmaNorm = sigma_Sm[i] / 10**Params.Induct.sigmaMax[Induction.bodyname]
            Dnorm = D_km[i] / np.max(D_km)
            ptColors[i] = Color.OceanCmap(Induction.compsList[iValidFlat[i]], sigmaNorm, Dnorm,
                                          DARKEN_SALINITIES=FigMisc.DARKEN_SALINITIES)
        else:
            w_ppt[i] = Induction.x[iValid[i]].flatten()
            yFlat[i] = Induction.y[iValid[i]].flatten()
            Tmean_K = Induction.Tmean_K[iValid[i]].flatten()
            if FigMisc.NORMALIZED_SALINITIES:
                wMax_ppt = np.array([Color.saturation[comp] for comp in Induction.compsList[iValidFlat[i]]])
                w_normFrac = w_ppt[i] / wMax_ppt
            else:
                w_normFrac = interp1d([np.min(w_ppt[i]), np.max(w_ppt[i])], [0.0, 1.0])(w_ppt[i])
            if Params.Induct.colorType == 'Tmean':
                if FigMisc.NORMALIZED_TEMPERATURES:
                    Tmean_normFrac = Color.GetNormT(Tmean_K)
                else:
                    Tmean_normFrac = interp1d([np.min(Tmean_K), np.max(Tmean_K)], [0.0, 1.0])(Tmean_K)
                ptColors[i] = Color.OceanCmap(Induction.compsList[iValidFlat[i]], w_normFrac, Tmean_normFrac,
                                              DARKEN_SALINITIES=FigMisc.DARKEN_SALINITIES)
            elif Params.Induct.colorType == 'zb':
                zb_km = Induction.zb_km[iValid[i]].flatten()
                zb_normFrac = interp1d([np.min(zb_km), np.max(zb_km)], [0.0, 1.0])(zb_km)
                ptColors[i] = Color.OceanCmap(Induction.compsList[iValidFlat[i]], w_normFrac, zb_normFrac,
                                              DARKEN_SALINITIES=FigMisc.DARKEN_SALINITIES)
            else:
                raise ValueError(f'Inductogram colortype {Params.Induct.colorType} not recognized.')

    widthPlot = 25
    widthCbar = 1
    if Params.Induct.inductOtype == 'sigma':
        comps = ['Ice']
        fig = plt.figure(figsize=FigSize.phaseSpaceSolo, constrained_layout=True)
        grid = GridSpec(1, 2, width_ratios=[widthPlot, widthCbar], figure=fig)
        axes = [fig.add_subplot(grid[0, 0])]
        if Style.GRIDS:
            axes[0].grid()
            axes[0].set_axisbelow(True)
        cbarAx = fig.add_subplot(grid[0, 1])
        cbarUnits = InductionList[0].zb_km[iValid[0]].flatten()
        cbarLabel = FigLbl.iceThickLbl
    else:
        comps = np.unique(InductionList[0].comps)
        if Params.Induct.colorType == 'Tmean':
            cbarUnits = InductionList[0].Tmean_K[iValid[0]].flatten()
            cbarLabel = FigLbl.oceanTempLbl
        elif Params.Induct.colorType == 'zb':
            cbarUnits = InductionList[0].zb_km[iValid[0]].flatten()
            cbarLabel = FigLbl.iceThickLbl

        fig = plt.figure(figsize=FigSize.phaseSpaceCombo, constrained_layout=True)
        nComps = np.size(comps)
        grid = GridSpec(1, 2 + nComps, width_ratios=np.append([widthPlot, widthPlot], [widthCbar for _ in range(nComps)]), figure=fig)
        axes = [fig.add_subplot(grid[0, i]) for i in range(2)]
        if Style.GRIDS:
            [ax.grid() for ax in axes]
            [ax.set_axisbelow(True) for ax in axes]
        cbarAxes = [fig.add_subplot(grid[0, i+2]) for i in range(nComps)]
        axes[1].set_xlabel(FigLbl.wLabel)
        axes[1].set_ylabel(FigLbl.yLabelInduct)
        axes[1].set_xscale(FigLbl.wScale)
        axes[1].set_yscale(FigLbl.yScaleInduct)
        axes[1].scatter(w_ppt[0], yFlat[0], s=Style.MW_Induction,
                        marker=Style.MS_Induction, c=ptColors[0])

    # Labels and titles
    if Params.TITLES:
        fig.suptitle(FigLbl.phaseSpaceTitle)
    axes[0].set_xlabel(FigLbl.sigMeanLabel)
    axes[0].set_ylabel(FigLbl.Dlabel)
    axes[0].set_xlim(FigLbl.sigLims)
    axes[0].set_ylim(FigLbl.Dlims)
    axes[0].set_xscale(FigLbl.sigScale)
    axes[0].set_yscale(FigLbl.Dscale)

    pts = {}
    cbar = {}
    if Params.Induct.inductOtype == 'sigma':
        pts[comps[0]] = axes[0].scatter(sigma_Sm[0], D_km[0], s=Style.MW_Induction,
                              marker=Style.MS_Induction, c=cbarUnits, cmap=Color.cmap[comps[0]])
        cbar[comps[0]] = fig.colorbar(pts[comps[0]], cax=cbarAx)
    else:
        for comp, cbarAx in zip(comps, cbarAxes):
            thisComp = InductionList[0].compsList[iValidFlat[0]] == comp
            pts[comp] = axes[0].scatter(sigma_Sm[0][thisComp], D_km[0][thisComp], s=Style.MW_Induction,
                            marker=Style.MS_Induction, c=cbarUnits[thisComp], cmap=Color.cmap[comp])
            cbar[comp] = fig.colorbar(pts[comp], cax=cbarAx)
            lbl = f'\ce{{{comp}}}'
            if not FigMisc.TEX_INSTALLED:
                lbl = FigLbl.StripLatexFromString(lbl)
            cbarAx.set_title(lbl, fontsize=FigMisc.cbarTitleSize)

            if FigMisc.MARK_INDUCT_BOUNDS:
                boundStyle = {'ls': Style.LS_BdipInset, 'lw': Style.LW_BdipInset, 'c': Color.BdipInset}
                _ = axes[0].plot(sigma_Sm[0][thisComp][w_ppt[0] == np.min(w_ppt[0])], D_km[0][thisComp][w_ppt[0] == np.min(w_ppt[0])], **boundStyle)
                _ = axes[0].plot(sigma_Sm[0][thisComp][w_ppt[0] == np.max(w_ppt[0])], D_km[0][thisComp][w_ppt[0] == np.max(w_ppt[0])], **boundStyle)
                _ = axes[0].plot(sigma_Sm[0][thisComp][yFlat[0] == np.min(yFlat[0])], D_km[0][thisComp][yFlat[0] == np.min(yFlat[0])], **boundStyle)
                _ = axes[0].plot(sigma_Sm[0][thisComp][yFlat[0] == np.max(yFlat[0])], D_km[0][thisComp][yFlat[0] == np.max(yFlat[0])], **boundStyle)

    cbar[comps[-1]].set_label(cbarLabel)
    fig.savefig(Params.FigureFiles.phaseSpace, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
    log.debug(f'InductOgram phase space plot saved to file: {Params.FigureFiles.phaseSpace}')
    plt.close()

    # Plot combination
    if Params.COMPARE and np.size(InductionList) > 1 and Params.Induct.inductOtype != 'sigma':
        comps = np.unique(np.append([], [Induction.comps for Induction in InductionList]))
        nComps = np.size(comps)
        figWidth = FigSize.phaseSpaceSolo[0] + nComps * FigMisc.cbarSpace
        fig, ax = plt.subplots(1, 1, figsize=(figWidth, FigSize.phaseSpaceSolo[1]))
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)

        # Labels and titles
        if Params.TITLES:
            fig.suptitle(FigLbl.phaseSpaceTitle)
        ax.set_xlabel(FigLbl.sigMeanLabel)
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
        combow_ppt = np.concatenate(tuple(wi for wi in w_ppt))
        comboyFlat = np.concatenate(tuple(yi for yi in yFlat))
        if Params.Induct.colorType == 'Tmean':
            comboCbarUnits = np.concatenate(tuple(Induction.Tmean_K[iValid[i]].flatten() for i,Induction in enumerate(InductionList)))
        elif Params.Induct.colorType == 'zb':
            comboCbarUnits = np.concatenate(tuple(Induction.zb_km[iValid[i]].flatten() for i,Induction in enumerate(InductionList)))
        pts = {}
        for comp in comps:
            thisComp = comboCompsList == comp
            pts[comp] = ax.scatter(comboSigma_Sm[thisComp], comboD_km[thisComp], s=Style.MW_Induction,
                                        marker=Style.MS_Induction, c=comboColors[thisComp])
            cbarAx = divider.new_horizontal(size=FigMisc.cbarSize, pad=FigMisc.cbarSpace + extraPad)
            extraPad = FigMisc.extraPad
            cbar = mcbar.ColorbarBase(cbarAx, cmap=Color.cmap[comp], format=FigMisc.cbarFmt,
                                      values=np.linspace(np.min(comboCbarUnits[thisComp]), np.max(comboCbarUnits[thisComp]), FigMisc.nCbarPts))
            fig.add_axes(cbarAx)
            lbl = f'\ce{{{comp}}}'
            if not FigMisc.TEX_INSTALLED:
                lbl = FigLbl.StripLatexFromString(lbl)
            cbarAx.set_title(lbl, fontsize=FigMisc.cbarTitleSize)

            if FigMisc.MARK_INDUCT_BOUNDS:
                thisw_ppt = combow_ppt[thisComp]
                thisyFlat = comboyFlat[thisComp]
                boundStyle = {'ls': Style.LS_BdipInset, 'lw': Style.LW_BdipInset, 'c': Color.BdipInset}
                _ = ax.plot(comboSigma_Sm[thisComp][thisw_ppt == np.min(thisw_ppt)], comboD_km[thisComp][thisw_ppt == np.min(thisw_ppt)], **boundStyle)
                _ = ax.plot(comboSigma_Sm[thisComp][thisw_ppt == np.max(thisw_ppt)], comboD_km[thisComp][thisw_ppt == np.max(thisw_ppt)], **boundStyle)
                _ = ax.plot(comboSigma_Sm[thisComp][thisyFlat == np.min(thisyFlat)], comboD_km[thisComp][thisyFlat == np.min(thisyFlat)], **boundStyle)
                _ = ax.plot(comboSigma_Sm[thisComp][thisyFlat == np.max(thisyFlat)], comboD_km[thisComp][thisyFlat == np.max(thisyFlat)], **boundStyle)

        cbar.set_label(cbarLabel, size=12)
        plt.tight_layout()
        fig.savefig(Params.FigureFiles.phaseSpaceCombo, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Plot saved to file: {Params.FigureFiles.phaseSpaceCombo}')
        plt.close()

    return


def PlotInductOgram(Induction, Params):
    """ Plot contours showing magnetic induction responses for an array of models
    """

    # Get indices for the oscillations that we can and want to plot
    excSelectionCalc = {key:(Texc is not None) for key, Texc in Induction.Texc_hr.items()}
    TexcCalcd = excSelectionCalc and Induction.Texc_hr
    whichTexc = excSelectionCalc and Params.Induct.excSelectionPlot
    allTexc_hr = np.fromiter(TexcCalcd.values(), dtype=np.float_)
    allAvailableTexc_hr = allTexc_hr[np.isfinite(allTexc_hr)]
    iTexc = [np.where(allTexc_hr == Texc)[0][0] for key, Texc in TexcCalcd.items() if whichTexc[key] and np.size(np.where(allTexc_hr == Texc)[0]) > 0]
    iTexcAvail = [np.where(allAvailableTexc_hr == Texc)[0][0] for key, Texc in TexcCalcd.items() if whichTexc[key] and np.size(np.where(allAvailableTexc_hr == Texc)[0]) > 0]
    TexcPlotNames = np.fromiter(Induction.Texc_hr.keys(), dtype='<U20')[iTexc]
    Texc_hr = allTexc_hr[iTexc]
    iSort = np.argsort(Texc_hr)
    if Induction.bodyname not in Params.Induct.sigmaMax.keys() or Induction.bodyname not in Params.Induct.sigmaMin.keys():
        sigmaMin = np.min(Induction.sigmaMean_Sm[np.isfinite(Induction.sigmaMean_Sm)])
        sigmaMax = np.max(Induction.sigmaMean_Sm[np.isfinite(Induction.sigmaMean_Sm)])
        Dmin = np.min(Induction.D_km[np.isfinite(Induction.D_km)])
        Dmax = np.max(Induction.D_km[np.isfinite(Induction.D_km)])
        log.warning(f'At least one of [sigmaMin, sigmaMax, Dmin, Dmax] are not set for {Induction.bodyname} in configPPinduct. ' +
                    f'Using min/max for this inductOgram run:\n' +
                    f'sigmaMin: {sigmaMin} S/m\n' +
                    f'sigmaMin: {sigmaMax} S/m\n' +
                    f'sigmaMin: {Dmin} S/m\n' +
                    f'sigmaMin: {Dmax} S/m\n')
        Params.Induct.sigmaMin[Induction.bodyname] = np.log10(sigmaMin)
        Params.Induct.sigmaMax[Induction.bodyname] = np.log10(sigmaMax)
        Params.Induct.Dmin[Induction.bodyname] = np.log10(Dmin)
        Params.Induct.Dmax[Induction.bodyname] = np.log10(Dmax)
        
    # Let matplotlib decide contour levels and format if we haven't explicitly set them
    if Induction.bodyname not in Params.Induct.cLevels.keys() or Induction.bodyname not in Params.Induct.cfmt.keys():
        Params.Induct.SPECIFIC_CLEVELS = False
        
    FigLbl.SetInduction(Induction.bodyname, Params.Induct, Texc_hr)
    if not FigMisc.TEX_INSTALLED:
        FigLbl.StripLatex()

    # Get all common labels and data for zipping
    zData = [Induction.Amp[iTexcAvail,...], np.abs(Induction.Bix_nT[iTexcAvail,...]),
             np.abs(Induction.Biy_nT[iTexcAvail,...]), np.abs(Induction.Biz_nT[iTexcAvail,...])]
    if Induction.SINGLE_COMP:
        FigLbl.singleComp(Induction.comps[0])

    # Adjust phi values in case we're plotting void volume % instead of void fraction
    if Params.Induct.inductOtype == 'phi':
        Induction.y = Induction.y * FigLbl.phiMult

    if Params.COMBINE_BCOMPS:
        # Plot B components all together with phase. Amplitude is still separate
        # Generate canvas and add labels
        fig = plt.figure(figsize=FigSize.inductCombo)
        grid = GridSpec(2, 2)
        axes = np.array([[fig.add_subplot(grid[i, j]) for j in range(2)] for i in range(2)])
        allAxes = axes.flatten()
        if Style.GRIDS:
            [ax.grid() for ax in allAxes]
            [ax.set_axisbelow(True) for ax in allAxes]

        if Params.TITLES:
            fig.suptitle(FigLbl.inductionTitle)
        # Only label the bottom-left sides of axes
        [ax.set_xlabel(FigLbl.sigMeanLabel) for ax in (axes[1,0], axes[1,1])]
        [ax.set_ylabel(FigLbl.Dlabel) for ax in (axes[0,0], axes[1,0])]
        [ax.set_xlim(FigLbl.sigLims) for ax in allAxes]
        [ax.set_ylim(FigLbl.Dlims) for ax in allAxes]
        [ax.set_xscale(FigLbl.sigScale) for ax in allAxes]
        [ax.set_yscale(FigLbl.Dscale) for ax in allAxes]
        coords = {'Bx': (0,0), 'By': (0,1), 'Bz': (1,0), 'phase': (1,1)}
        comboData = [np.abs(Induction.Bix_nT[iTexcAvail,...]), np.abs(Induction.Biy_nT[iTexcAvail,...]),
                     np.abs(Induction.Biz_nT[iTexcAvail,...]), Induction.phase[iTexcAvail,...]]
        comboTitles = np.append(FigLbl.plotTitles[1:], FigLbl.phaseTitle)
        comboLabels = list(coords.keys())

        for z, name, fLabel in zip(comboData, comboTitles, comboLabels):
            ax = axes[coords[fLabel]]
            ax.set_title(name)
            zContours = [ax.contour(Induction.sigmaMean_Sm, Induction.D_km, z[i, ...],
                           colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                           linewidths=Style.LW_Induction[T], levels=Params.Induct.GetClevels(fLabel, T))
                           for i, T in enumerate(TexcPlotNames)]
            if Params.Induct.inductOtype == 'sigma':
                [ax.clabel(zContours[i], fmt=Params.Induct.GetCfmt(fLabel, T),
                           fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                           for i, T in enumerate(TexcPlotNames)]
            elif FigMisc.MARK_INDUCT_BOUNDS:
                boundStyle = {'ls': Style.LS_BdipInset, 'lw': Style.LW_BdipInset, 'c': Color.BdipInset}
                _ = ax.plot(Induction.sigmaMean_Sm[Induction.x == np.min(Induction.x)], Induction.D_km[Induction.x == np.min(Induction.x)], **boundStyle)
                _ = ax.plot(Induction.sigmaMean_Sm[Induction.x == np.max(Induction.x)], Induction.D_km[Induction.x == np.max(Induction.x)], **boundStyle)
                _ = ax.plot(Induction.sigmaMean_Sm[Induction.y == np.min(Induction.y)], Induction.D_km[Induction.y == np.min(Induction.y)], **boundStyle)
                _ = ax.plot(Induction.sigmaMean_Sm[Induction.y == np.max(Induction.y)], Induction.D_km[Induction.y == np.max(Induction.y)], **boundStyle)

        if FigMisc.PLOT_V2021 and Induction.bodyname in ['Europa', 'Ganymede', 'Callisto']:
            AddV2021points(Params.Induct, Induction.bodyname, 'sigma', allAxes)

        if Params.LEGEND:
            lines = np.array([contour.legend_elements()[0][0] for contour in zContours])
            axes[1,1].legend(lines[iSort], FigLbl.legendTexc[iSort], framealpha=FigMisc.cLegendOpacity)

        if Params.Induct.inductOtype == 'sigma':
            fNameSigma = Params.FigureFiles.sigmaOnly['Bcomps']
        else:
            fNameSigma = Params.FigureFiles.sigma['Bcomps']

        plt.tight_layout()
        fig.savefig(fNameSigma, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Plot saved to file: {fNameSigma}')
        plt.close()

        if Params.Induct.inductOtype != 'sigma':
            fig = plt.figure(figsize=FigSize.inductCombo)
            grid = GridSpec(2, 2)
            axes = np.array([[fig.add_subplot(grid[i, j]) for j in range(2)] for i in range(2)])
            allAxes = axes.flatten()
            if Style.GRIDS:
                [ax.grid() for ax in allAxes]
                [ax.set_axisbelow(True) for ax in allAxes]

            if Params.TITLES:
                fig.suptitle(FigLbl.inductionTitle)
            # Only label the bottom-left sides of axes
            [ax.set_xlabel(FigLbl.wLabel) for ax in (axes[1,0], axes[1,1])]
            [ax.set_ylabel(FigLbl.yLabelInduct) for ax in (axes[0,0], axes[1,0])]
            [ax.set_xscale(FigLbl.wScale) for ax in allAxes]
            [ax.set_yscale(FigLbl.yScaleInduct) for ax in allAxes]

            for z, name, fLabel in zip(comboData, comboTitles, comboLabels):
                ax = axes[coords[fLabel]]
                ax.set_title(name)
                zContours = [ax.contour(Induction.x, Induction.y, z[i, ...],
                                        colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                                        linewidths=Style.LW_Induction[T],
                                        levels=Params.Induct.GetClevels(fLabel, T))
                             for i, T in enumerate(TexcPlotNames)]
                [ax.clabel(zContours[i], fmt=Params.Induct.GetCfmt(fLabel, T),
                           fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                           for i, T in enumerate(TexcPlotNames)]

            if FigMisc.PLOT_V2021 and Induction.bodyname in ['Europa', 'Ganymede', 'Callisto']:
                AddV2021points(Params.Induct, Induction.bodyname, Params.Induct.inductOtype, allAxes)

            if Params.LEGEND:
                lines = np.array([contour.legend_elements()[0][0] for contour in zContours])
                axes[1,1].legend(lines[iSort], FigLbl.legendTexc[iSort], framealpha=FigMisc.cLegendOpacity)

            plt.tight_layout()
            fig.savefig(Params.FigureFiles.induct['Bcomps'], format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
            log.debug(f'Plot saved to file: {Params.FigureFiles.induct["Bcomps"]}')
            plt.close()

        if Params.Induct.inductOtype != 'sigma':
            # Also plot a comparison of Bx, which is usually the strongest oscillation
            compChoice = 'Bx'
            fig = plt.figure(figsize=FigSize.inductCombo)
            grid = GridSpec(2, 2)
            axes = np.array([[fig.add_subplot(grid[i, j]) for j in range(2)] for i in range(2)])
            fig.subplots_adjust(wspace=0.25, hspace=0.35)
            allAxes = axes.flatten()
            if Style.GRIDS:
                [ax.grid() for ax in allAxes]
                [ax.set_axisbelow(True) for ax in allAxes]

            # Label all axes for clarity
            if Params.TITLES:
                fig.suptitle(FigLbl.inductCompareTitle)
            [ax.set_xlabel(FigLbl.wLabel) for ax in axes[0,:]]
            [ax.set_ylabel(FigLbl.yLabelInduct) for ax in axes[0,:]]
            [ax.set_xlabel(FigLbl.sigMeanLabel) for ax in axes[1,:]]
            [ax.set_ylabel(FigLbl.Dlabel) for ax in axes[1,:]]
            [ax.set_xscale(FigLbl.sigScale) for ax in allAxes]
            [ax.set_yscale(FigLbl.Dscale) for ax in axes[1,:]]
            [ax.set_yscale(FigLbl.yScaleInduct) for ax in axes[0,:]]
            [ax.set_xlim(FigLbl.sigLims) for ax in axes[1,:]]
            [ax.set_ylim(FigLbl.Dlims) for ax in axes[1,:]]

            axes[0,0].set_title(comboTitles[0])
            axes[1,0].set_title(comboTitles[0])
            axes[0,1].set_title(comboTitles[-1])
            axes[1,1].set_title(comboTitles[-1])
            zContours = [axes[0,0].contour(Induction.x, Induction.y, comboData[0][i, ...],
                                    colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                                    linewidths=Style.LW_Induction[T],
                                    levels=Params.Induct.GetClevels(comboLabels[0], T))
                         for i, T in enumerate(TexcPlotNames)]
            phaseContours = [axes[0,1].contour(Induction.x, Induction.y, comboData[-1][i, ...],
                                    colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                                    linewidths=Style.LW_Induction[T],
                                    levels=Params.Induct.GetClevels(comboLabels[-1], T))
                         for i, T in enumerate(TexcPlotNames)]
            [axes[1,0].contour(Induction.sigmaMean_Sm, Induction.D_km, comboData[0][i, ...],
                                    colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                                    linewidths=Style.LW_Induction[T],
                                    levels=Params.Induct.GetClevels(comboLabels[0], T))
                         for i, T in enumerate(TexcPlotNames)]
            [axes[1,1].contour(Induction.sigmaMean_Sm, Induction.D_km, comboData[-1][i, ...],
                                    colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                                    linewidths=Style.LW_Induction[T],
                                    levels=Params.Induct.GetClevels(comboLabels[-1], T))
                         for i, T in enumerate(TexcPlotNames)]
            [axes[0,0].clabel(zContours[i], fmt=Params.Induct.GetCfmt(comboLabels[0], T),
                       fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                       for i, T in enumerate(TexcPlotNames)]
            [axes[0,1].clabel(phaseContours[i], fmt=Params.Induct.GetCfmt(comboLabels[-1], T),
                       fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                       for i, T in enumerate(TexcPlotNames)]

            if FigMisc.PLOT_V2021 and Induction.bodyname in ['Europa', 'Ganymede', 'Callisto']:
                AddV2021points(Params.Induct, Induction.bodyname, Params.Induct.inductOtype, axes[0,:])
                AddV2021points(Params.Induct, Induction.bodyname, 'sigma', axes[1,:])

            if Params.LEGEND:
                lines = np.array([contour.legend_elements()[0][0] for contour in zContours])
                axes[1,1].legend(lines[iSort], FigLbl.legendTexc[iSort], framealpha=FigMisc.cLegendOpacity)

            if FigMisc.MARK_INDUCT_BOUNDS:
                boundStyle = {'ls': Style.LS_BdipInset, 'lw': Style.LW_BdipInset, 'c': Color.BdipInset}
                _ = [ax.plot(Induction.sigmaMean_Sm[Induction.x == np.min(Induction.x)], Induction.D_km[Induction.x == np.min(Induction.x)], **boundStyle) for ax in axes[1,:]]
                _ = [ax.plot(Induction.sigmaMean_Sm[Induction.x == np.max(Induction.x)], Induction.D_km[Induction.x == np.max(Induction.x)], **boundStyle) for ax in axes[1,:]]
                _ = [ax.plot(Induction.sigmaMean_Sm[Induction.y == np.min(Induction.y)], Induction.D_km[Induction.y == np.min(Induction.y)], **boundStyle) for ax in axes[1,:]]
                _ = [ax.plot(Induction.sigmaMean_Sm[Induction.y == np.max(Induction.y)], Induction.D_km[Induction.y == np.max(Induction.y)], **boundStyle) for ax in axes[1,:]]

            plt.tight_layout()
            fig.savefig(Params.FigureFiles.inductCompare[compChoice], format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
            log.debug(f'Plot saved to file: {Params.FigureFiles.inductCompare[compChoice]}')
            plt.close()

        # Set lists to just contain Amplitude now to reuse the remaining routines for that plot
        zData = [zData[0]]
        FigLbl.plotTitles = [FigLbl.plotTitles[0]]
        FigLbl.fLabels = [FigLbl.fLabels[0]]

    # Plot each component separately alongside phase
    for z, name, fLabel in zip(zData, FigLbl.plotTitles, FigLbl.fLabels):

        # Generate canvas and add labels
        fig = plt.figure(figsize=FigSize.induct)
        grid = GridSpec(1, 2)
        axes = [fig.add_subplot(grid[0, j]) for j in range(2)]
        if Style.GRIDS:
            [ax.grid() for ax in axes]
            [ax.set_axisbelow(True) for ax in axes]

        # Labels and titles
        if Params.TITLES:
            fig.suptitle(FigLbl.inductionTitle)
        axes[0].set_title(name)
        axes[1].set_title(FigLbl.phaseTitle)
        [ax.set_xlabel(FigLbl.sigMeanLabel) for ax in axes]
        [ax.set_ylabel(FigLbl.Dlabel) for ax in axes]
        [ax.set_xlim(FigLbl.sigLims) for ax in axes]
        [ax.set_ylim(FigLbl.Dlims) for ax in axes]
        [ax.set_xscale(FigLbl.sigScale) for ax in axes]
        [ax.set_yscale(FigLbl.Dscale) for ax in axes]

        zContours = [axes[0].contour(Induction.sigmaMean_Sm, Induction.D_km, z[i, ...],
                         colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                         linewidths=Style.LW_Induction[T], levels=Params.Induct.GetClevels(fLabel, T))
                         for i, T in enumerate(TexcPlotNames)]
        phaseContours = [axes[1].contour(Induction.sigmaMean_Sm, Induction.D_km, Induction.phase[i, ...],
                         colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                         linewidths=Style.LW_Induction[T], levels=Params.Induct.GetClevels('phase', T))
                         for i, T in enumerate(TexcPlotNames)]
        if Params.Induct.inductOtype == 'sigma':
            [axes[0].clabel(zContours[i], fmt=Params.Induct.GetCfmt(fLabel, T),
                            fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                            for i, T in enumerate(TexcPlotNames)]
            [axes[1].clabel(phaseContours[i], fmt=Params.Induct.GetCfmt('phase', T),
                            fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                            for i, T in enumerate(TexcPlotNames)]
            fNameSigma = Params.FigureFiles.sigmaOnly[fLabel]
        else:
            fNameSigma = Params.FigureFiles.sigma[fLabel]

            if FigMisc.MARK_INDUCT_BOUNDS:
                boundStyle = {'ls': Style.LS_BdipInset, 'lw': Style.LW_BdipInset, 'c': Color.BdipInset}
                _ = [ax.plot(Induction.sigmaMean_Sm[Induction.x == np.min(Induction.x)], Induction.D_km[Induction.x == np.min(Induction.x)], **boundStyle) for ax in axes]
                _ = [ax.plot(Induction.sigmaMean_Sm[Induction.x == np.max(Induction.x)], Induction.D_km[Induction.x == np.max(Induction.x)], **boundStyle) for ax in axes]
                _ = [ax.plot(Induction.sigmaMean_Sm[Induction.y == np.min(Induction.y)], Induction.D_km[Induction.y == np.min(Induction.y)], **boundStyle) for ax in axes]
                _ = [ax.plot(Induction.sigmaMean_Sm[Induction.y == np.max(Induction.y)], Induction.D_km[Induction.y == np.max(Induction.y)], **boundStyle) for ax in axes]

        if FigMisc.PLOT_V2021 and Induction.bodyname in ['Europa', 'Ganymede', 'Callisto']:
            AddV2021points(Params.Induct, Induction.bodyname, 'sigma', axes)

        if Params.LEGEND:
            lines = np.array([contour.legend_elements()[0][0] for contour in phaseContours])
            axes[1].legend(lines[iSort], FigLbl.legendTexc[iSort], framealpha=FigMisc.cLegendOpacity)

        plt.tight_layout()
        fig.savefig(fNameSigma, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Plot saved to file: {fNameSigma}')
        plt.close()

        if Params.Induct.inductOtype != 'sigma':
            fig = plt.figure(figsize=FigSize.induct)
            grid = GridSpec(1, 2)
            axes = [fig.add_subplot(grid[0, j]) for j in range(2)]
            if Style.GRIDS:
                [ax.grid() for ax in axes]
                [ax.set_axisbelow(True) for ax in axes]

            # Labels and titles
            if Params.TITLES:
                fig.suptitle(FigLbl.inductionTitle)
            axes[0].set_title(name)
            axes[1].set_title(FigLbl.phaseTitle)
            [ax.set_xscale(FigLbl.wScale) for ax in axes]
            [ax.set_xlabel(FigLbl.wLabel) for ax in axes]
            [ax.set_ylabel(FigLbl.yLabelInduct) for ax in axes]
            [ax.set_yscale(FigLbl.yScaleInduct) for ax in axes]

            zContours = [axes[0].contour(Induction.x, Induction.y, z[i, ...],
                             colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                             linewidths=Style.LW_Induction[T], levels=Params.Induct.GetClevels(fLabel, T))
                             for i, T in enumerate(TexcPlotNames)]
            phaseContours = [axes[1].contour(Induction.x, Induction.y, Induction.phase[i, ...],
                             colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                             linewidths=Style.LW_Induction[T], levels=Params.Induct.GetClevels('phase', T))
                             for i, T in enumerate(TexcPlotNames)]
            [axes[0].clabel(zContours[i], fmt=Params.Induct.GetCfmt(fLabel, T),
                            fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                            for i, T in enumerate(TexcPlotNames)]
            [axes[1].clabel(phaseContours[i], fmt=Params.Induct.GetCfmt('phase', T),
                            fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                            for i, T in enumerate(TexcPlotNames)]

            if FigMisc.PLOT_V2021 and Induction.bodyname in ['Europa', 'Ganymede', 'Callisto']:
                AddV2021points(Params.Induct, Induction.bodyname, Params.Induct.inductOtype, axes)

            if Params.LEGEND:
                lines = np.array([contour.legend_elements()[0][0] for contour in phaseContours])
                axes[1].legend(lines[iSort], FigLbl.legendTexc[iSort], framealpha=FigMisc.cLegendOpacity)

            plt.tight_layout()
            fig.savefig(Params.FigureFiles.induct[fLabel], format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
            log.debug(f'Plot saved to file: {Params.FigureFiles.induct[fLabel]}')
            plt.close()

    return


def AddV2021points(IndParams, bodyname, inductOtype, axes):
    if inductOtype == 'sigma' or inductOtype == 'Tb':
        log.debug('Adding Vance et al. (2021) markers to inductogram.')
        if inductOtype == 'sigma':
            xVals = IndParams.V2021_sigma_Sm[bodyname]
            yVals = IndParams.V2021_D_km[bodyname]
        else:
            xVals = IndParams.V2021_w_ppt[bodyname]
            yVals = IndParams.V2021_Tb_K[bodyname]
        UP = IndParams.V2021_zb_km[bodyname] == np.min(IndParams.V2021_zb_km[bodyname])
        DOWN = np.logical_not(UP)
        [ax.scatter(xVals[UP], yVals[UP], marker=IndParams.V2021_MS[bodyname][UP][0],
                    facecolor=IndParams.V2021_FC[bodyname][UP],
                    edgecolor=IndParams.V2021_EC[bodyname][UP],
                    label=r'Vance et al.\ (2021) models') for ax in axes]
        [ax.scatter(xVals[DOWN], yVals[DOWN], marker=IndParams.V2021_MS[bodyname][DOWN][0],
                    facecolor=IndParams.V2021_FC[bodyname][DOWN],
                    edgecolor=IndParams.V2021_EC[bodyname][DOWN]
                    ) for ax in axes]
    else:
        log.warning(f'FigMisc.PLOT_V2021 is True but inductOtype "{inductOtype}" is '
                    'not supported. Skipping.')

    return


def PlotComplexBdip(PlanetList, Params):
    if PlanetList[0].bodyname in Excitations.Texc_hr.keys():
        refs = {}
        nPeaksToPlot = np.sum([Tdo and Thave for Tdo, Thave, Texc in zip(Params.Induct.excSelectionPlot.values(),
                                                                         Params.Induct.excSelectionCalc.values(),
                                                                         Excitations.Texc_hr[
                                                                             PlanetList[0].bodyname].values()) if
                               Texc is not None])
        nPeaksCalced = np.sum([np.size(Planeti.Magnetic.calcedExc) for Planeti in PlanetList])
        if nPeaksToPlot >= 1 and nPeaksCalced >= 1:
            if nPeaksToPlot == 1:
                DO_ZOOM = False
                nCol = 1
                titleToUse = [FigLbl.BdipTitleNoZoom, FigLbl.BdipCompareTitleNoZoom]
            else:
                DO_ZOOM = True
                nCol = 2
                titleToUse = [FigLbl.BdipTitle, FigLbl.BdipCompareTitle]

            if FigMisc.MANUAL_HYDRO_COLORS:
                comps = np.unique([Planet.Ocean.comp for Planet in PlanetList if Planet.Ocean.comp != 'none'])
                wMinMax_ppt = {}
                TminMax_K = {}

                for comp in comps:
                    wAll_ppt = [Planet.Ocean.wOcean_ppt for Planet in PlanetList if Planet.Ocean.comp == comp]
                    wMinMax_ppt[comp] = [np.min(wAll_ppt), np.max(wAll_ppt)]
                    Tall_K = [Planet.Bulk.Tb_K for Planet in PlanetList if Planet.Ocean.comp == comp]
                    TminMax_K[comp] = [np.min(Tall_K), np.max(Tall_K)]
                    # Reset to default if all models are the same or if desired
                    if not FigMisc.RELATIVE_Tb_K or TminMax_K[comp][0] == TminMax_K[comp][1]:
                        TminMax_K[comp] = Color.Tbounds_K

            if Params.COMBINE_BCOMPS:
                if DO_ZOOM:
                    figSize = FigSize.BdipCombo
                else:
                    figSize = FigSize.BdipSoloCombo
                fig = plt.figure(figsize=figSize)
                grid = GridSpec(3, nCol)
                axes = np.array([[fig.add_subplot(grid[i, j]) for j in range(nCol)] for i in range(3)])
                axc = {vComp: row for vComp, row in zip(xyzComps, axes)}
                if Style.GRIDS:
                    axf = axes.flatten()
                    [ax.grid() for ax in axf]
                    [ax.set_axisbelow(True) for ax in axf]

                # Labels and titles
                [[ax.set_xlabel(FigLbl.BdipReLabel[axComp]) for ax in axs] for axComp, axs in axc.items()]
                [[ax.set_ylabel(FigLbl.BdipImLabel[axComp]) for ax in axs] for axComp, axs in axc.items()]
                if DO_ZOOM:
                    [axs[0].set_title(FigLbl.BdipZoomLabel[axComp]) for axComp, axs in axc.items()]
                    [axs[-1].set_title(FigLbl.BdipLabel[axComp]) for axComp, axs in axc.items()]
                else:
                    [axs[-1].set_title(FigLbl.BdipLabelNoZoom[axComp]) for axComp, axs in axc.items()]
                if Params.TITLES:
                    if Params.ALL_ONE_BODY:
                        fig.suptitle(f'{PlanetList[0].name}{titleToUse[0]}')
                    else:
                        fig.suptitle(titleToUse[1])

                insetx, insety = ({vComp: 0 for vComp in xyzComps} for _ in range(2))
                for iPlanet, Planet in enumerate(PlanetList):
                    if Planet.Magnetic.calcedExc is not None and np.size(Planet.Magnetic.calcedExc) >= 1:
                        # Set color options
                        if FigMisc.MANUAL_HYDRO_COLORS and not Planet.Do.NO_H2O:
                            if wMinMax_ppt[Planet.Ocean.comp][0] != wMinMax_ppt[Planet.Ocean.comp][1]:
                                thisAlpha = Style.GetMA(Planet.Ocean.wOcean_ppt, wMinMax_ppt[Planet.Ocean.comp])
                            else:
                                thisAlpha = Style.MAlims[-1]
                            Color.Tbounds_K = TminMax_K[Planet.Ocean.comp]
                            thisColor = Color.GetNormT(Planet.Bulk.Tb_K)
                            thisEdgeColor = Color.cmap[Planet.Ocean.comp](thisColor)
                            thisFaceColor = Color.cmap[Planet.Ocean.comp](thisColor, alpha=thisAlpha)
                        else:
                            thisEdgeColor = None
                            thisFaceColor = None

                        for iRow, axComp in enumerate(xyzComps):
                            xPlotted, yPlotted, absBiPlotted = (np.empty(0) for _ in range(3))
                            for iPeak, Tkey in enumerate(Planet.Magnetic.calcedExc):
                                if Params.Induct.excSelectionPlot[Tkey]:
                                    pts = [ax.scatter(np.real(Planet.Magnetic.Bi1xyz_nT[axComp][iPeak]),
                                                      np.imag(Planet.Magnetic.Bi1xyz_nT[axComp][iPeak]),
                                                      marker=Style.MS_dip[Tkey], s=Style.MW_dip[Tkey] ** 2,
                                                      facecolor=thisFaceColor, edgecolor=thisEdgeColor) for ax in
                                           axc[axComp]][0].get_offsets()[0, :]

                                    if iPlanet == 0 and iRow == 0:
                                        refs[Tkey] = axes[0, -1].scatter(0, 0, label=Tkey.capitalize(),
                                                                         marker=Style.MS_dip[Tkey],
                                                                         s=Style.MW_dip[Tkey] ** 2, facecolor=Color.ref,
                                                                         edgecolor=Color.ref)

                                    xPlotted = np.append(xPlotted, pts[0])
                                    yPlotted = np.append(yPlotted, pts[1])
                                    absBiPlotted = np.append(absBiPlotted,
                                                             np.abs(Planet.Magnetic.Bi1xyz_nT[axComp][iPeak]))

                            if DO_ZOOM:
                                iAllButLargest = np.argsort(absBiPlotted)[:-1]
                                secondLargestRe = np.max(xPlotted[iAllButLargest])
                                secondLargestIm = np.max(yPlotted[iAllButLargest])
                                insetx[axComp] = np.maximum(insetx[axComp], secondLargestRe * FigMisc.BdipZoomMult)
                                insety[axComp] = np.maximum(insety[axComp], secondLargestIm * FigMisc.BdipZoomMult)
                                axes[iRow, 0].set_xlim([0, insetx[axComp]])
                                axes[iRow, 0].set_ylim([0, insety[axComp]])

                if DO_ZOOM and FigMisc.SHOW_INSET:
                    refs[f'inset'] = \
                    axes[0, -1].plot([0, 1], [1, 0], color=Color.BdipInset, linewidth=Style.LW_BdipInset,
                                     linestyle=Style.LS_BdipInset, label='Inset region')[0]
                    for iRow, axComp in enumerate(xyzComps):
                        axes[iRow, -1].add_patch(
                            Rectangle((0, 0), insetx[axComp], insety[axComp], edgecolor=Color.BdipInset, zorder=-1,
                                      linewidth=Style.LW_BdipInset, linestyle=Style.LS_BdipInset, facecolor='None'))

                if Params.LEGEND:
                    axes[0, -1].legend()

                for refPt in refs.values():
                    refPt.remove()

                [axes[iRow, -1].set_xlim(left=0) for iRow in range(3)]
                [axes[iRow, -1].set_ylim(bottom=0) for iRow in range(3)]
                plt.tight_layout()
                fig.savefig(Params.FigureFiles.Bdip['all'], format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
                log.debug(f'Induced dipole surface strength plot saved to file: {Params.FigureFiles.Bdip["all"]}')
                plt.close()

            else:
                for axComp in xyzComps:
                    if DO_ZOOM:
                        figSize = FigSize.Bdip
                    else:
                        figSize = FigSize.BdipSolo
                    fig = plt.figure(figsize=figSize)
                    grid = GridSpec(1, nCol)
                    axes = np.array([fig.add_subplot(grid[0, j]) for j in range(nCol)])
                    if Style.GRIDS:
                        [ax.grid() for ax in axes]
                        [ax.set_axisbelow(True) for ax in axes]

                    # Labels and titles
                    [ax.set_xlabel(FigLbl.BdipReLabel[axComp]) for ax in axes]
                    [ax.set_ylabel(FigLbl.BdipImLabel[axComp]) for ax in axes]
                    if DO_ZOOM:
                        axes[0].set_title(FigLbl.BdipZoomLabel[axComp])
                        axes[-1].set_title(FigLbl.BdipLabel[axComp])
                    else:
                        axes[-1].set_title(FigLbl.BdipLabelNoZoom[axComp])
                    if Params.TITLES:
                        if Params.ALL_ONE_BODY:
                            fig.suptitle(f'{PlanetList[0].name}{titleToUse[0]}')
                        else:
                            fig.suptitle(titleToUse[1])

                    insetx, insety = (0, 0)
                    for iPlanet, Planet in enumerate(PlanetList):
                        if Planet.Magnetic.calcedExc is not None and np.size(Planet.Magnetic.calcedExc) >= 1:
                            # Set color options
                            if FigMisc.MANUAL_HYDRO_COLORS and not Planet.Do.NO_H2O:
                                if wMinMax_ppt[Planet.Ocean.comp][0] != wMinMax_ppt[Planet.Ocean.comp][1]:
                                    thisAlpha = Style.GetMA(Planet.Ocean.wOcean_ppt, wMinMax_ppt[Planet.Ocean.comp])
                                else:
                                    thisAlpha = Style.MAlims[-1]
                                Color.Tbounds_K = TminMax_K[Planet.Ocean.comp]
                                thisColor = Color.GetNormT(Planet.Bulk.Tb_K)
                                thisEdgeColor = Color.cmap[Planet.Ocean.comp](thisColor)
                                thisFaceColor = Color.cmap[Planet.Ocean.comp](thisColor, alpha=thisAlpha)
                            else:
                                thisEdgeColor = None
                                thisFaceColor = None

                            xPlotted, yPlotted, absBiPlotted = (np.empty(0) for _ in range(3))
                            for iPeak, Tkey in enumerate(Planet.Magnetic.calcedExc):
                                if Params.Induct.excSelectionPlot[Tkey]:
                                    pts = [ax.scatter(np.real(Planet.Magnetic.Bi1xyz_nT[axComp][iPeak]),
                                                      np.imag(Planet.Magnetic.Bi1xyz_nT[axComp][iPeak]),
                                                      marker=Style.MS_dip[Tkey], s=Style.MW_dip[Tkey] ** 2,
                                                      facecolor=thisFaceColor, edgecolor=thisEdgeColor) for ax in axes][
                                              0].get_offsets()[0, :]

                                    if iPlanet == 0:
                                        refs[Tkey] = axes[-1].scatter(0, 0, label=Tkey.capitalize(),
                                                                      marker=Style.MS_dip[Tkey],
                                                                      s=Style.MW_dip[Tkey] ** 2, facecolor=Color.ref,
                                                                      edgecolor=Color.ref)

                                    xPlotted = np.append(xPlotted, pts[0])
                                    yPlotted = np.append(yPlotted, pts[1])
                                    absBiPlotted = np.append(absBiPlotted,
                                                             np.abs(Planet.Magnetic.Bi1xyz_nT[axComp][iPeak]))

                            if np.size(absBiPlotted) == 0:
                                log.warning(f'No moments were plotted for {Planet.name}. Selected were: {Params.Induct.excSelectionPlot}.' +
                                            f'Calculated were: {Planet.Magnetic.calcedExc}. Zoom plotting will be skipped.')
                            elif DO_ZOOM:
                                iAllButLargest = np.argsort(absBiPlotted)[:-1]
                                secondLargestRe = np.max(xPlotted[iAllButLargest])
                                secondLargestIm = np.max(yPlotted[iAllButLargest])
                                insetx = np.maximum(insetx, secondLargestRe * FigMisc.BdipZoomMult)
                                insety = np.maximum(insety, secondLargestIm * FigMisc.BdipZoomMult)
                                axes[0].set_xlim([0, insetx])
                                axes[0].set_ylim([0, insety])

                    if DO_ZOOM and FigMisc.SHOW_INSET:
                        axes[-1].add_patch(Rectangle((0, 0), insetx, insety, edgecolor=Color.BdipInset, zorder=-1,
                                                     linewidth=Style.LW_BdipInset, linestyle=Style.LS_BdipInset,
                                                     facecolor='None'))
                        refs['inset'] = \
                        axes[-1].plot([0, 1], [1, 0], color=Color.BdipInset, linewidth=Style.LW_BdipInset,
                                      linestyle=Style.LS_BdipInset, label='Inset region')[0]

                    if Params.LEGEND:
                        axes[-1].legend()

                    for refPt in refs.values():
                        refPt.remove()

                    axes[-1].set_xlim(left=0)
                    axes[-1].set_ylim(bottom=0)
                    plt.tight_layout()
                    fig.savefig(Params.FigureFiles.Bdip[axComp], format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
                    log.debug(f'Induced dipole surface strength plot saved to file: {Params.FigureFiles.Bdip[axComp]}')
                    plt.close()

    return


def PlotMagSurface(PlanetList, Params):
    """ Plot magnetic field on a spherical surface of fixed radius using
        calculations from MoonMag.
    """
    PlanetListSubset = [Planet for Planet in PlanetList if Planet.bodyname in Excitations.Texc_hr.keys()]
    nts = np.size(FigMisc.tMagEval_s)
    BvecRe_nT, BvecReSym_nT, tMagEvalLbl, tMagEvalPrint, tFnameEnd, rMagEvalLbl, rMagEvalPrint \
        = (np.empty(nts, dtype=object) for _ in range(7))
    Bdiff_nT, BvecReDiff_nT = (None, None)

    # Confirm tMagEval and rMagEval are in a format we can use together
    if not isinstance(FigMisc.tMagEval_s, Iterable):
        FigMisc.tMagEval_s = np.array([FigMisc.tMagEval_s])
    if not isinstance(FigMisc.rMagEval_Rp, Iterable):
        FigMisc.rMagEval_Rp = np.array([FigMisc.rMagEval_Rp])
    if np.size(FigMisc.tMagEval_s) != np.size(FigMisc.rMagEval_Rp):
        if np.size(FigMisc.tMagEval_s) == 1:
            FigMisc.tMagEval_s = np.repeat(FigMisc.tMagEval_s[0], np.size(FigMisc.rMagEval_Rp))
        elif np.size(FigMisc.rMagEval_Rp) == 1:
            FigMisc.rMagEval_Rp = np.repeat(FigMisc.rMagEval_Rp[0], np.size(FigMisc.tMagEval_s))
        else:
            raise ValueError('Lists of tMagEval and rMagEval must have the same length or one must be length 1.')

    for iPlanet, Planet in enumerate(PlanetListSubset):
        phaseNow = np.zeros((nts, np.size(Planet.Magnetic.omegaExc_radps)), dtype=np.complex_)

        # Get plot text and axis labels
        if Planet.lonMap_deg is None:
            Planet.lonMap_deg = FigMisc.lonMap_deg
            Planet.nLonMap = np.size(FigMisc.lonMap_deg)
        if Planet.phiMap_rad is None:
            Planet.phiMap_rad = np.radians(FigMisc.lonMap_deg)
        if Planet.latMap_deg is None:
            Planet.latMap_deg = FigMisc.latMap_deg
            Planet.nLatMap = np.size(FigMisc.latMap_deg)
        if Planet.thetaMap_rad is None:
            Planet.thetaMap_rad = np.radians(90 - FigMisc.latMap_deg)
        if FigMisc.tMagLbl is None:
            tMagEvalLbl, tMagEvalPrint, tFnameEnd = FigLbl.tStr(FigMisc.tMagEval_s)
        else:
            tMagEvalLbl, tMagEvalPrint, tFnameEnd = FigLbl.tStrManual(FigMisc.tMagLbl, nts=nts)
        rMagEvalLbl, rMagEvalPrint = FigLbl.rStr(FigMisc.rMagEval_Rp, Planet.bodyname)

        if Params.CALC_ASYM:
            asymStr = ',3D'
            # Get separate spherically symmetric induced moments
            BinmSph_nT = np.zeros_like(Planet.Magnetic.Benm_nT)
            for n in range(1, Planet.Magnetic.nprmMax + 1):
                for iExc in range(Planet.Magnetic.nExc):
                    BinmSph_nT[iExc,:,n,:] = n/(n+1) * Planet.Magnetic.Aen[iExc,n] * Planet.Magnetic.Benm_nT[iExc,:,n,:]

        else:
            asymStr = ''
            BinmSph_nT = None

        # Evaluate magnetic field on the selected surface at the selected time
        for iEval, tEval_s in enumerate(FigMisc.tMagEval_s):
            BvecRe_nT[iEval] = {vComp: np.zeros((np.size(FigMisc.latMap_deg), np.size(FigMisc.lonMap_deg))) for vComp in vecComps}
            # Get complex phase factor for time dependence
            phaseNow[iEval,:] = np.exp(-1j * Planet.Magnetic.omegaExc_radps * tEval_s)

            # Evaluate field at selected time and radius
            log.debug(f'Evaluating induced magnetic field at {rMagEvalPrint[iEval]}, {tMagEvalPrint[iEval]}.')
            for iExc in range(Planet.Magnetic.nExc):
                BinmNow_nT = Planet.Magnetic.Binm_nT[iExc, ...] * phaseNow[iEval][iExc]
                Bx_nT, By_nT, Bz_nT = GetMagSurf(Planet.Magnetic.nLin, Planet.Magnetic.mLin, BinmNow_nT,
                                                 FigMisc.rMagEval_Rp[iEval], Planet.thetaMap_rad, Planet.phiMap_rad,
                                                 do_parallel=Params.DO_PARALLEL)
                BvecRe_nT[iEval]['x'] = BvecRe_nT[iEval]['x'] + np.real(Bx_nT)
                BvecRe_nT[iEval]['y'] = BvecRe_nT[iEval]['y'] + np.real(By_nT)
                BvecRe_nT[iEval]['z'] = BvecRe_nT[iEval]['z'] + np.real(Bz_nT)
                BvecRe_nT[iEval]['mag'] = np.sqrt(BvecRe_nT[iEval]['x']**2 + BvecRe_nT[iEval]['y']**2 + BvecRe_nT[iEval]['z']**2)

            if Params.CALC_ASYM:
                # Also evaluate field at selected time and radius for symmetric case for comparison
                log.debug(f'Evaluating symmetric induced magnetic field at {rMagEvalPrint[iEval]}, {tMagEvalPrint[iEval]}.')
                BvecReSym_nT[iEval] = {vComp: np.zeros((np.size(FigMisc.latMap_deg), np.size(FigMisc.lonMap_deg))) for vComp in vecComps}
                for iExc in range(Planet.Magnetic.nExc):
                    BinmSphNow_nT = BinmSph_nT[iExc, ...] * phaseNow[iEval][iExc]
                    Bx_nT, By_nT, Bz_nT = GetMagSurf(Planet.Magnetic.nprmLin, Planet.Magnetic.mprmLin, BinmSphNow_nT,
                                                     FigMisc.rMagEval_Rp[iEval], Planet.thetaMap_rad, Planet.phiMap_rad,
                                                     do_parallel=Params.DO_PARALLEL)
                    BvecReSym_nT[iEval]['x'] = BvecReSym_nT[iEval]['x'] + np.real(Bx_nT)
                    BvecReSym_nT[iEval]['y'] = BvecReSym_nT[iEval]['y'] + np.real(By_nT)
                    BvecReSym_nT[iEval]['z'] = BvecReSym_nT[iEval]['z'] + np.real(Bz_nT)
                    BvecReSym_nT[iEval]['mag'] = np.sqrt(BvecReSym_nT[iEval]['x']**2 + BvecReSym_nT[iEval]['y']**2 + BvecReSym_nT[iEval]['z']**2)

        if Params.CALC_ASYM:
            Bdiff_nT = np.array([{vComp: BvecRe_nT[iEval][vComp] - BvecReSym_nT[iEval][vComp]
                                  for vComp in vecComps} for iEval in range(nts)], dtype=object)

        if Params.COMPARE and np.size(PlanetListSubset) > 1:
            # Save first profile in list and take the difference for each subsequent profile
            if iPlanet == 0:
                BvecReBase_nT = deepcopy(BvecRe_nT)
            else:
                BvecReDiff_nT = np.array([{vComp: BvecRe_nT[iEval][vComp] - BvecReBase_nT[iEval][vComp]
                                 for vComp in vecComps} for iEval in range(nts)], dtype=object)
                # Set 'mag' comp of model comparison to be the magnitude of the difference vector, which is a better
                # estimate of signal size than the difference of the magnitudes
                for iEval in range(nts):
                    BvecReDiff_nT[iEval]['mag'] = np.sqrt(BvecReDiff_nT[iEval]['x']**2 + BvecReDiff_nT[iEval]['y']**2 + BvecReDiff_nT[iEval]['z']**2)

        if FigMisc.vCompMagSurf == 'all':
            vComps2plot = vecComps
        else:
            vComps2plot = [FigMisc.vCompMagSurf]

        for vCompMagSurf in vComps2plot:
            # Get colormap and contour settings for surface plots
            cmap, vmin, vmax, cLevels = GetZparams(BvecRe_nT, vCompMagSurf, nts, Bdiff_nT=Bdiff_nT, Bcomp_nT=BvecReDiff_nT)

            # Plot individual components on surface
            for iEval, tEval_s in enumerate(FigMisc.tMagEval_s):

                if np.size(PlanetListSubset) == 1:
                    if FigMisc.BASYM_WITH_SYM:
                        fig = plt.figure(figsize=FigSize.MagSurfCombo)
                        grid = GridSpec(1, 2)
                    else:
                        fig = plt.figure(figsize=FigSize.MagSurf)
                        grid = GridSpec(1, 1)

                    ax = fig.add_subplot(grid[0, -1])
                    SetMap(ax)
                    if FigMisc.LARGE_ADJUST:
                        if vCompMagSurf == 'mag':
                            compLabel = f'$B_\mathrm{{{vCompMagSurf}{asymStr}}}$'
                        else:
                            compLabel = f'$B_{{{vCompMagSurf}\mathrm{{{asymStr}}}}}$'
                        title = f'{FigLbl.MagSurfShortTitle} {compLabel} ' + \
                                f'at {rMagEvalLbl[iEval]} ($\si{{nT}}$)'
                    else:
                        if vCompMagSurf == 'mag':
                            compLabel = 'magnitude'
                        else:
                            compLabel = f'${vCompMagSurf}$ component'
                        title = f'{Planet.name} {FigLbl.MagSurfTitle} {compLabel}, ' + \
                                f'{rMagEvalLbl[iEval]} ($\si{{nT}}$), {tMagEvalLbl[iEval]}'

                    if Params.TITLES:
                        # Remove latex styling from labels if Latex is not installed
                        if not FigMisc.TEX_INSTALLED:
                            title = FigLbl.StripLatexFromString(title)
                        ax.set_title(title, size=FigMisc.mapTitleSize)

                    Bmap = ax.pcolormesh(FigMisc.lonMap_deg, FigMisc.latMap_deg, BvecRe_nT[iEval][vCompMagSurf],
                                         shading='auto', cmap=cmap['vec'], vmin=vmin['vec'][iEval],
                                         vmax=vmax['vec'][iEval], rasterized=FigMisc.PT_RASTER)
                    BmapContours = ax.contour(FigMisc.lonMap_deg, FigMisc.latMap_deg, BvecRe_nT[iEval][vCompMagSurf],
                                              levels=cLevels['vec'][iEval], colors='black')
                    ax.clabel(BmapContours, fmt=ticker.FuncFormatter(FigMisc.Cformat),
                              fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)

                    if not FigMisc.LARGE_ADJUST:
                        cbar = fig.colorbar(Bmap, ax=ax)
                        cbar.ax.set_title(FigLbl.MagSurfCbarTitle)

                    ax.set_aspect(1)
                    if not (Params.CALC_ASYM and FigMisc.BASYM_WITH_SYM):
                        plt.tight_layout()
                        fName = f'{Params.FigureFiles.MagSurf[vCompMagSurf]}{tFnameEnd[iEval]}{FigMisc.xtn}'
                        fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
                        log.debug(f'Induced field surface map saved to file: {fName}')
                        plt.close()

                    # Plot symmetric component and difference to asymmetric if applicable
                    if Params.CALC_ASYM:
                        if not FigMisc.BASYM_WITH_SYM:
                            fig = plt.figure(figsize=FigSize.MagSurf)
                            grid = GridSpec(1, 1)

                        ax = fig.add_subplot(grid[0, 0])
                        SetMap(ax)
                        if FigMisc.LARGE_ADJUST:
                            if vCompMagSurf == 'mag':
                                compLabel = f'$B_\mathrm{{{vCompMagSurf},sym}}$'
                            else:
                                compLabel = f'$B_{{{vCompMagSurf}\mathrm{{,sym}}}}$'
                            title = f'{FigLbl.MagSurfShortTitle} {compLabel} ' + \
                                    f'at {rMagEvalLbl[iEval]} ($\si{{nT}}$)'
                        else:
                            # compLabel able to be reused from asymmetric in this case
                            title = f'{Planet.name} {FigLbl.MagSurfSymTitle} {compLabel}, ' + \
                                    f'{rMagEvalLbl[iEval]} ($\si{{nT}}$), {tMagEvalLbl[iEval]}'
                        if Params.TITLES:
                            # Remove latex styling from labels if Latex is not installed
                            if not FigMisc.TEX_INSTALLED:
                                title = FigLbl.StripLatexFromString(title)
                            ax.set_title(title, size=FigMisc.mapTitleSize)

                        Bmap = ax.pcolormesh(FigMisc.lonMap_deg, FigMisc.latMap_deg, BvecReSym_nT[iEval][vCompMagSurf],
                                             shading='auto', cmap=cmap['vec'], vmin=vmin['vec'][iEval],
                                             vmax=vmax['vec'][iEval], rasterized=FigMisc.PT_RASTER)
                        BmapContours = ax.contour(FigMisc.lonMap_deg, FigMisc.latMap_deg, BvecReSym_nT[iEval][vCompMagSurf],
                                                  levels=cLevels['vec'][iEval], colors='black')
                        ax.clabel(BmapContours, fmt=ticker.FuncFormatter(FigMisc.Cformat),
                                  fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)

                        if not FigMisc.LARGE_ADJUST:
                            cbar = fig.colorbar(Bmap, ax=ax)
                            cbar.ax.set_title(FigLbl.MagSurfCbarTitle)

                        ax.set_aspect(1)
                        plt.tight_layout()
                        if FigMisc.BASYM_WITH_SYM:
                            fName = f'{Params.FigureFiles.MagSurfCombo[vCompMagSurf]}{tFnameEnd[iEval]}{FigMisc.xtn}'
                            fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
                            log.debug(f'Symmetric/asymmetric induced field surface maps saved to file: {fName}')
                        else:
                            fName = f'{Params.FigureFiles.MagSurfSym[vCompMagSurf]}{tFnameEnd[iEval]}{FigMisc.xtn}'
                            fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
                            log.debug(f'Symmetric induced field surface map saved to file: {fName}')

                        plt.close()

                        # Now plot the difference between asymmetric and symmetric for this component
                        fig = plt.figure(figsize=FigSize.MagSurf)
                        grid = GridSpec(1, 1)
                        ax = fig.add_subplot(grid[0, 0])
                        SetMap(ax)
                        if FigMisc.LARGE_ADJUST:
                            if vCompMagSurf == 'mag':
                                compLabel = f'$\Delta B_\mathrm{{{vCompMagSurf},sym}}$'
                            else:
                                compLabel = f'$\Delta B_{{{vCompMagSurf}\mathrm{{,sym}}}}$'
                            title = f'{FigLbl.MagSurfShortTitle} {compLabel} ' + \
                                    f'at {rMagEvalLbl[iEval]} ($\si{{nT}}$)'
                        else:
                            # compLabel able to be reused from above in this case
                            title = f'{Planet.name} {FigLbl.MagSurfTitle} {compLabel} {FigLbl.MagSurfDiffTitle}, ' + \
                                    f'{rMagEvalLbl[iEval]} ($\si{{nT}}$), {tMagEvalLbl[iEval]}'
                        if Params.TITLES:
                            # Remove latex styling from labels if Latex is not installed
                            if not FigMisc.TEX_INSTALLED:
                                title = FigLbl.StripLatexFromString(title)
                            ax.set_title(title, size=FigMisc.mapTitleSize)

                        Bmap = ax.pcolormesh(FigMisc.lonMap_deg, FigMisc.latMap_deg, Bdiff_nT[iEval][vCompMagSurf],
                                             shading='auto', cmap=cmap['diff'], vmin=vmin['diff'][iEval],
                                             vmax=vmax['diff'][iEval], rasterized=FigMisc.PT_RASTER)
                        BmapContours = ax.contour(FigMisc.lonMap_deg, FigMisc.latMap_deg, Bdiff_nT[iEval][vCompMagSurf],
                                                  levels=cLevels['diff'][iEval], colors='black')
                        ax.clabel(BmapContours, fmt=ticker.FuncFormatter(FigMisc.Cformat),
                                  fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)

                        if not FigMisc.LARGE_ADJUST:
                            fig.colorbar(Bmap, ax=ax, label=FigLbl.MagSurfCbarDiffLabel)

                        ax.set_aspect(1)
                        plt.tight_layout()
                        fName = f'{Params.FigureFiles.MagSurfDiff[vCompMagSurf]}{tFnameEnd[iEval]}{FigMisc.xtn}'
                        fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
                        log.debug(f'Induced field difference surface map saved to file: {fName}')
                        plt.close()

                # Plot difference map between each model we want to compare to the first in the list
                if Params.COMPARE and iPlanet != 0:

                    fig = plt.figure(figsize=FigSize.MagSurf)
                    grid = GridSpec(1, 1)
                    ax = fig.add_subplot(grid[0, 0])
                    SetMap(ax)
                    if FigMisc.LARGE_ADJUST:
                        if vCompMagSurf == 'mag':
                            compLabel = f'$|\Delta \mathbf{{B}}|$'
                        else:
                            compLabel = f'$\Delta B_{{{vCompMagSurf}}}$'
                        title = f'{FigLbl.MagSurfCompareShortTitle} {compLabel} ' + \
                                f'at {rMagEvalLbl[iEval]} ($\si{{nT}}$)'
                    else:
                        if vCompMagSurf == 'mag':
                            compLabel = 'magnitude'
                        else:
                            compLabel = f'${vCompMagSurf}$ component'
                        title = f'{Planet.name} {FigLbl.MagSurfCompareTitle} {compLabel}, ' + \
                                f'{rMagEvalLbl[iEval]} ($\si{{nT}}$), {tMagEvalLbl[iEval]}'

                    if Params.TITLES:
                        # Remove latex styling from labels if Latex is not installed
                        if not FigMisc.TEX_INSTALLED:
                            title = FigLbl.StripLatexFromString(title)
                        ax.set_title(title, size=FigMisc.mapTitleSize)

                    Bmap = ax.pcolormesh(FigMisc.lonMap_deg, FigMisc.latMap_deg, BvecReDiff_nT[iEval][vCompMagSurf],
                                         shading='auto', cmap=cmap['comp'], vmin=vmin['comp'][iEval],
                                         vmax=vmax['comp'][iEval], rasterized=FigMisc.PT_RASTER)
                    BmapContours = ax.contour(FigMisc.lonMap_deg, FigMisc.latMap_deg, BvecReDiff_nT[iEval][vCompMagSurf],
                                              levels=cLevels['comp'][iEval], colors='black')
                    ax.clabel(BmapContours, fmt=ticker.FuncFormatter(FigMisc.Cformat),
                              fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)

                    if not FigMisc.LARGE_ADJUST:
                        cbar = fig.colorbar(Bmap, ax=ax)
                        cbar.ax.set_title(FigLbl.MagSurfCbarTitle)

                    ax.set_aspect(1)
                    plt.tight_layout()
                    fName = f'{Params.FigureFiles.MagSurfComp[vCompMagSurf]}{tFnameEnd[iEval]}{FigMisc.xtn}'
                    fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
                    log.debug(f'Induced field surface map saved to file: {fName}')
                    plt.close()

    return


def SetMap(ax):
    ax.set_xticks(FigMisc.lonMapTicks_deg)
    ax.set_yticks(FigMisc.latMapTicks_deg)
    ax.tick_params(axis='both', which='major', labelsize=FigMisc.latlonSize)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(FigMisc.LonMapFormatter))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(FigMisc.LatMapFormatter))

    return


def GetZparams(Bvec_nT, vCompMagSurf, nts, Bdiff_nT=None, Bcomp_nT=None):
    """ Configure cmap, vmin, vmax, and levels kwargs for pcolormesh and contour plotting
    """

    vTypes = ['vec', 'diff', 'comp']
    vmin, vmax = ({vType: np.empty(nts) for vType in vTypes} for _ in range(2))
    cLevels = {vType: np.empty(nts, dtype=object) for vType in vTypes}
    cmap = {vType: None for vType in vTypes}

    # @@@@@@@@@@@@@@@@
    # @@ Main field @@
    # @@@@@@@@@@@@@@@@
    # Decide which colormap to use for main field plotting
    if vCompMagSurf == 'mag':
        cmap['vec'] = Color.cmap['BmapPos']
    else:
        cmap['vec'] = Color.cmap['BmapDiv']

    # Get min, max, and contour list for main field plotting
    if FigMisc.FIXED_COLORBAR:
        vmaxAll = {vComp: [np.max(np.abs(Bvec_nT[iEval][vComp])) for iEval in range(nts)] for vComp in vecComps}
        vminAll = {vComp: [np.min(np.abs(Bvec_nT[iEval][vComp])) for iEval in range(nts)] for vComp in vecComps}
        if FigMisc.FIXED_ALL_COMPS and not FigMisc.MAG_CBAR_SEPARATE:
            vAbsMax = np.max([vmaxAll[vComp] for vComp in vecComps])
            vAbsMin = np.min([vminAll[vComp] for vComp in vecComps])
        else:
            vAbsMax = np.max(vmaxAll[vCompMagSurf])
            vAbsMin = np.min(vminAll[vCompMagSurf])

        if vCompMagSurf != 'mag':
            vAbsMin = -vAbsMax

        vmax['vec'] = np.repeat(vAbsMax, nts)
        vmin['vec'] = np.repeat(vAbsMin, nts)

    else:
        if FigMisc.vmaxMagSurf_nT is None:
            vmax['vec'] = np.array([np.max(Bvec_nT[iEval][vCompMagSurf]) for iEval in range(nts)])
        else:
            vmax['vec'] = np.repeat(FigMisc.vmaxMagSurf_nT, nts)
        if FigMisc.vminMagSurf_nT is None:
            vmin['vec'] = np.array([np.min(Bvec_nT[iEval][vCompMagSurf]) for iEval in range(nts)])
        else:
            vmin['vec'] = np.repeat(FigMisc.vminMagSurf_nT, nts)

    if FigMisc.nMagContours is not None:
        for iEval in range(nts):
            cLevels['vec'][iEval] = np.linspace(vmin['vec'][iEval], vmax['vec'][iEval], FigMisc.nMagContours)

    # @@@@@@@@@@@@@@@@@@@
    # @@ Asym/sym diff @@
    # @@@@@@@@@@@@@@@@@@@
    if Bdiff_nT is not None:
        # Decide which colormap to use for difference plotting
        if np.all([Bdiff_nT[iEval][vCompMagSurf] >= 0 for iEval in range(nts)]):
            cmap['diff'] = Color.cmap['BmapPos']
            neg = 'none'
        elif np.all([Bdiff_nT[iEval][vCompMagSurf] < 0 for iEval in range(nts)]):
            cmap['diff'] = Color.cmap['BmapNeg']
            neg = 'all'
        else:
            cmap['diff'] = Color.cmap['BmapDiv']
            neg = 'some'

        # Get min, max, and contour list for difference plotting
        if FigMisc.FIXED_COLORBAR:
            vmaxAll = {vComp: [np.max(np.abs(Bdiff_nT[iEval][vComp])) for iEval in range(nts)] for vComp in vecComps}
            vminAll = {vComp: [np.min(np.abs(Bdiff_nT[iEval][vComp])) for iEval in range(nts)] for vComp in vecComps}
            if FigMisc.FIXED_ALL_COMPS and not FigMisc.MAG_CBAR_SEPARATE:
                vAbsMax = np.max([vmaxAll[vComp] for vComp in vecComps])
                vAbsMin = np.min([vminAll[vComp] for vComp in vecComps])
            else:
                vAbsMax = np.max(vmaxAll[vCompMagSurf])
                vAbsMin = np.min(vminAll[vCompMagSurf])

            if neg == 'all':
                vAbsMax = -vAbsMax
                vAbsMin = -vAbsMin
            elif neg == 'some':
                vAbsMin = -vAbsMax

            vmax['diff'] = np.repeat(vAbsMax, nts)
            vmin['diff'] = np.repeat(vAbsMin, nts)

        else:
            if FigMisc.vmaxMagSurf_nT is None:
                vmax['diff'] = np.array([np.max(Bdiff_nT[iEval][vCompMagSurf]) for iEval in range(nts)])
            else:
                vmax['diff'] = np.repeat(FigMisc.vmaxMagSurfDiff_nT, nts)
            if FigMisc.vminMagSurf_nT is None:
                vmin['diff'] = np.array([np.min(Bdiff_nT[iEval][vCompMagSurf]) for iEval in range(nts)])
            else:
                vmin['diff'] = np.repeat(FigMisc.vminMagSurfDiff_nT, nts)

        if FigMisc.nMagContours is not None:
            for iEval in range(nts):
                cLevels['diff'][iEval] = np.linspace(vmin['diff'][iEval], vmax['diff'][iEval], FigMisc.nMagContours)

    # @@@@@@@@@@@@@@@@@@@@@@
    # @@ Model comparison @@
    # @@@@@@@@@@@@@@@@@@@@@@
    if Bcomp_nT is not None:
        # Decide which colormap to use for model comparison plotting
        if np.all([Bcomp_nT[iEval][vCompMagSurf] >= 0 for iEval in range(nts)]):
            cmap['comp'] = Color.cmap['BmapPos']
            neg = 'none'
        elif np.all([Bcomp_nT[iEval][vCompMagSurf] < 0 for iEval in range(nts)]):
            cmap['comp'] = Color.cmap['BmapNeg']
            neg = 'all'
        else:
            cmap['comp'] = Color.cmap['BmapDiv']
            neg = 'some'

        # Get min, max, and contour list for model comparison plotting
        if FigMisc.FIXED_COLORBAR:
            vmaxAll = {vComp: [np.max(np.abs(Bcomp_nT[iEval][vComp])) for iEval in range(nts)] for vComp in vecComps}
            vminAll = {vComp: [np.min(np.abs(Bcomp_nT[iEval][vComp])) for iEval in range(nts)] for vComp in vecComps}
            if FigMisc.FIXED_ALL_COMPS and not FigMisc.MAG_CBAR_SEPARATE:
                vAbsMax = np.max([vmaxAll[vComp] for vComp in vecComps])
                vAbsMin = np.min([vminAll[vComp] for vComp in vecComps])
            else:
                vAbsMax = np.max(vmaxAll[vCompMagSurf])
                vAbsMin = np.min(vminAll[vCompMagSurf])

            if neg == 'all':
                vAbsMax = -vAbsMax
                vAbsMin = -vAbsMin
            elif neg == 'some':
                vAbsMin = -vAbsMax

            vmax['comp'] = np.repeat(vAbsMax, nts)
            vmin['comp'] = np.repeat(vAbsMin, nts)

        else:
            if FigMisc.vmaxMagSurf_nT is None:
                vmax['comp'] = np.array([np.max(Bcomp_nT[iEval][vCompMagSurf]) for iEval in range(nts)])
            else:
                vmax['comp'] = np.repeat(FigMisc.vmaxMagSurfDiff_nT, nts)
            if FigMisc.vminMagSurf_nT is None:
                vmin['comp'] = np.array([np.min(Bcomp_nT[iEval][vCompMagSurf]) for iEval in range(nts)])
            else:
                vmin['comp'] = np.repeat(FigMisc.vminMagSurfDiff_nT, nts)

        if FigMisc.nMagContours is not None:
            for iEval in range(nts):
                cLevels['comp'][iEval] = np.linspace(vmin['comp'][iEval], vmax['comp'][iEval], FigMisc.nMagContours)

    return cmap, vmin, vmax, cLevels


def PlotAsym(PlanetList, Params):
    """ Plot a contour map showing the asymmetry present in selected boundaries.
        This function is time-consuming, so it is only applied to the first body
        in PlanetList (the primary body)
    """

    if Params.CALC_ASYM:
        if np.size(PlanetList[0].Magnetic.iAsymBds) > 0:
            mainPlanet = PlanetList[0]
        else:
            mainPlanet = next(Planet for Planet in PlanetList if np.size(Planet.Magnetic.iAsymBds) > 0)

        if mainPlanet.lonMap_deg is None:
            mainPlanet.lonMap_deg = FigMisc.lonMap_deg
            mainPlanet.nLonMap = np.size(FigMisc.lonMap_deg)
        if mainPlanet.latMap_deg is None:
            mainPlanet.latMap_deg = FigMisc.latMap_deg
            mainPlanet.nLatMap = np.size(FigMisc.latMap_deg)

        for i, zMean_km in enumerate(mainPlanet.Magnetic.zMeanAsym_km):
            fig = plt.figure(figsize=FigSize.asym)
            grid = GridSpec(1, 1)
            ax = fig.add_subplot(grid[0, 0])
            SetMap(ax)
            if mainPlanet.Magnetic.asymDescrip is not None and i < np.size(mainPlanet.Magnetic.asymDescrip):
                title = f'{mainPlanet.Magnetic.asymDescrip[i]}{FigLbl.asymAfterTitle}'
                if mainPlanet.Magnetic.asymDescrip[i] in mainPlanet.Magnetic.asymContours_km.keys():
                    cLevelsAsym = mainPlanet.Magnetic.asymContours_km[mainPlanet.Magnetic.asymDescrip[i]]
                else:
                    cLevelsAsym = None
            else:
                if mainPlanet.Magnetic.asymPlotType == 'ice':
                    title = f'{FigLbl.asymIceTitle}\SI{{{np.abs(zMean_km):.1f}}}{{km}}'
                elif mainPlanet.Magnetic.asymPlotType == 'surf':
                    title = f'{FigLbl.asymSurfTitle}{mainPlanet.name[0]}}}$ = \SI{{{mainPlanet.Bulk.R_m/1e3:.1f}}}{{km}}'
                elif mainPlanet.Magnetic.asymPlotType == 'ionos':
                    if np.size(mainPlanet.Magnetic.zMeanAsym_km) == 1:
                        hMeanAsym_km = np.squeeze(mainPlanet.Magnetic.zMeanAsym_km)
                    else:
                        hMeanAsym_km = mainPlanet.Magnetic.zMeanAsym_km[-1]
                    title = f'{FigLbl.asymIonosTitle}\SI{{{-hMeanAsym_km:.1f}}}{{km}}'
                elif mainPlanet.Magnetic.asymPlotType == 'depth':
                    if np.size(mainPlanet.Magnetic.zMeanAsym_km) == 1:
                        zMeanAsym_km = np.squeeze(mainPlanet.Magnetic.zMeanAsym_km)
                    else:
                        zMeanAsym_km = mainPlanet.Magnetic.zMeanAsym_km[0]
                    title = f'{FigLbl.asymDepthTitle}\SI{{{zMeanAsym_km:.1f}}}{{km}}'
                else:
                    log.warning(f'asymPlotType "{mainPlanet.Magnetic.asymPlotType}" not recognized. Using fallback.')
                    title = f'{FigLbl.asymSurfTitle}{mainPlanet.name[0]}}}$ = \SI{{{mainPlanet.R_m/1e3:.1f}}}{{km}}'
                cLevelsAsym = None
            if Params.TITLES:
                # Remove latex styling from labels if Latex is not installed
                if not FigMisc.TEX_INSTALLED:
                    title = FigLbl.StripLatexFromString(title)
                ax.set_title(title, size=FigMisc.mapTitleSize)

            if cLevelsAsym is None and FigMisc.nAsymContours is not None:
                cLevelsAsym = np.unique(np.round(np.linspace(np.min(mainPlanet.Magnetic.asymDevs_km[i, ...]),
                                                             np.max(mainPlanet.Magnetic.asymDevs_km[i, ...]),
                                                             FigMisc.nAsymContours)))

            asymMap = ax.pcolormesh(FigMisc.lonMap_deg, FigMisc.latMap_deg,
                                    mainPlanet.Magnetic.asymDevs_km[i, ...],
                                    shading='auto', cmap=Color.cmap['asymDev'], rasterized=FigMisc.PT_RASTER)
            asymContours = ax.contour(FigMisc.lonMap_deg, FigMisc.latMap_deg,
                                      mainPlanet.Magnetic.asymDevs_km[i, ...],
                                      levels=cLevelsAsym, colors='black')
            ax.clabel(asymContours, fmt=ticker.FuncFormatter(FigMisc.Cformat),
                      fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)

            if not FigMisc.LARGE_ADJUST:
                fig.colorbar(asymMap, ax=ax, label=FigLbl.asymCbarLabel)

            ax.set_aspect(1)
            plt.tight_layout()
            fName = f'{Params.FigureFiles.asym}z{zMean_km:.1f}km{FigMisc.xtn}'
            fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
            log.debug(f'Asymmetric boundary surface map for z = {zMean_km:.1f} km saved to file: {fName}')
            plt.close()

    return


def PlotMagSpectrumInduced(PlanetList, Params):
    """ Plot Fourier spectra of magnetic excitations, complex response
        amplitude, and induced dipole components.
    """

    # Get first entry in PlanetList for which FT_LOADED is True
    if PlanetList[0].Magnetic.FT_LOADED:
        mainPlanet = PlanetList[0]
    else:
        mainPlanet = next(Planet for Planet in PlanetList if Planet.Magnetic.FT_LOADED)

    # Create figure
    fig = plt.figure(figsize=FigSize.MagFT)
    grid = GridSpec(3, 1)
    axes = np.array([fig.add_subplot(grid[i, 0]) for i in range(3)])
    if Style.GRIDS:
        [ax.grid() for ax in axes]
        [ax.set_axisbelow(True) for ax in axes]

    # Labels and titles
    if FigMisc.Tmin_hr is None:
        Tmin_hr = np.min(mainPlanet.Magnetic.TexcFT_hr)
    else:
        Tmin_hr = FigMisc.Tmin_hr
    Tmax_hr = np.minimum(np.max(mainPlanet.Magnetic.TexcFT_hr), mainPlanet.Magnetic.TmaxFT_hr)
    Texc_hr = np.array(list(Excitations.Texc_hr[mainPlanet.bodyname].values()), dtype=np.float_)
    Texc_hr = Texc_hr[np.isfinite(Texc_hr)]
    if FigMisc.MAG_SPECTRA_PERIODS:
        freqLabel = FigLbl.TexcLabel
        freqVar = mainPlanet.Magnetic.TexcFT_hr
        freqLims = [Tmin_hr, Tmax_hr]
        mainExc = np.sort(Texc_hr)
    else:
        freqLabel = FigLbl.fExcLabel
        freqVar = 1 / mainPlanet.Magnetic.TexcFT_hr / 3600
        freqLims = [1 / Tmax_hr / 3600, 1 / Tmin_hr / 3600]
        mainExc = np.sort(1 / Texc_hr / 3600)
    [ax.set_xlim(freqLims) for ax in axes]
    [ax.set_xlabel(freqLabel) for ax in axes]
    [ax.set_xscale('log') for ax in axes]
    axes[1].set_ylim([0, 1])
    axes[0].set_ylabel(FigLbl.BeFTlabel)
    axes[1].set_ylabel(FigLbl.Ae1FTlabel)
    axes[2].set_ylabel(FigLbl.BiFTlabel)
    axes[0].set_title(FigLbl.BeFTtitle)
    axes[1].set_title(FigLbl.Ae1FTtitle)
    axes[2].set_title(FigLbl.BiFTtitle)
    [ax.set_yscale('log') for ax in [axes[0], axes[2]]]
    if Params.TITLES:
        if mainPlanet.Magnetic.coordTypeFT is not None:
            FTtitle = f'{mainPlanet.name} {FigLbl.MagFTtitle}, {mainPlanet.Magnetic.coordTypeFT} coordinates'
        else:
            FTtitle = f'{mainPlanet.name} {FigLbl.MagFTtitle}'
        fig.suptitle(FTtitle)

    # Plot the data
    [axes[0].plot(freqVar, np.abs(mainPlanet.Magnetic.Be1xyzFT_nT[vComp]), label=f'$B^e_{vComp}$',
                  color=Color.BeiFT[vComp], linestyle=Style.LS_FT, linewidth=Style.LW_FT) for vComp in xyzComps]
    [axes[2].plot(freqVar, np.abs(mainPlanet.Magnetic.Bi1xyzFT_nT[vComp]), label=f'$B^i_{vComp}$',
                  color=Color.BeiFT[vComp], linestyle=Style.LS_FT, linewidth=Style.LW_FT) for vComp in xyzComps]
    if np.size(PlanetList) > 1 and Params.COMPARE:
        # Only plot comparisons on amplitude plot
        [axes[1].plot(freqVar, np.abs(Planet.Magnetic.Ae1FT), label=Planet.label,
                      color=Color.Ae1FT, linestyle=Style.LS_FT, linewidth=Style.LW_FT)
         for Planet in PlanetList if Planet.Magnetic.FT_LOADED and Planet.bodyname == mainPlanet.bodyname]
    else:
        axes[1].plot(freqVar, np.abs(mainPlanet.Magnetic.Ae1FT), label=r'$\mathcal{A}^e_1$',
                     color=Color.Ae1FT, linestyle=Style.LS_FT, linewidth=Style.LW_FT)
    if FigMisc.MARK_TEXC:
        axes[1].plot([mainExc[0], mainExc[0]], [0, 1], label=r'Dominant $T_\mathrm{exc}$',
                     color=Color.TexcFT, linestyle=Style.LS_TexcFT, linewidth=Style.LW_TexcFT)
        [axes[1].plot([Tmark, Tmark], [0, 1], color=Color.TexcFT, linestyle=Style.LS_TexcFT,
                      linewidth=Style.LW_TexcFT) for Tmark in mainExc[1:]]

    if Params.LEGEND:
        [ax.legend() for ax in axes]

    plt.tight_layout()
    fig.savefig(Params.FigureFiles.MagFT, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
    log.debug(f'Magnetic Fourier spectra plot saved to file: {Params.FigureFiles.MagFT}')
    plt.close()

    return


def PlotMagSpectrum(PlanetList, Params):
    """ Plot Fourier spectra of magnetic excitations with optional highest peak annotated.
    """

    # Get first entry in PlanetList for which FT_LOADED is True
    if PlanetList[0].Magnetic.FT_LOADED:
        mainPlanet = PlanetList[0]
    else:
        mainPlanet = next(Planet for Planet in PlanetList if Planet.Magnetic.FT_LOADED)

    # Create figure
    fig = plt.figure(figsize=FigSize.MagFTexc)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    if Style.GRIDS:
        ax.grid(color='#CCCCCC', linestyle='-', linewidth=0.40, which='major')
        ax.grid(color='#CCCCCC', linestyle=':', linewidth=0.25, which='minor')
        ax.set_axisbelow(True)

    # Labels and titles
    if FigMisc.Tmin_hr is None:
        Tmin_hr = np.min(mainPlanet.Magnetic.TexcFT_hr)
    else:
        Tmin_hr = FigMisc.Tmin_hr
    Tmax_hr = np.minimum(np.max(mainPlanet.Magnetic.TexcFT_hr), mainPlanet.Magnetic.TmaxFT_hr)
    Texc_hr = np.array(list(Excitations.Texc_hr[mainPlanet.bodyname].values()), dtype=np.float_)
    Texc_hr = Texc_hr[np.isfinite(Texc_hr)]
    if FigMisc.MAG_SPECTRA_PERIODS:
        freqLabel = FigLbl.TexcLabel
        freqVar = mainPlanet.Magnetic.TexcFT_hr
        freqUnits = FigLbl.TexcUnits
        freqLims = [Tmin_hr, Tmax_hr]
        mainExc = np.sort(Texc_hr)
    else:
        freqLabel = FigLbl.fExcLabel
        freqVar = 1 / mainPlanet.Magnetic.TexcFT_hr / 3600
        freqUnits = FigLbl.fExcUnits
        freqLims = [1 / Tmax_hr / 3600, 1 / Tmin_hr / 3600]
        mainExc = np.sort(1 / Texc_hr / 3600)
    ax.set_xlim(freqLims)
    ax.set_xlabel(freqLabel)
    ax.set_ylabel(FigLbl.BeFTexcLabel)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    ax.xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.0f'))
    if Params.TITLES:
        if mainPlanet.Magnetic.coordTypeFT is not None:
            FTtitle = f'{mainPlanet.name} {FigLbl.MagFTexcTitle}, ' + \
                      f'{mainPlanet.Magnetic.coordTypeFT} coordinates'
        else:
            FTtitle = f'{mainPlanet.name} {FigLbl.MagFTExcTitle}'
        fig.suptitle(FTtitle)

    # Plot the data
    Bexc_nT = {vComp: np.abs(mainPlanet.Magnetic.Be1xyzFT_nT[vComp]) for vComp in xyzComps}
    [ax.plot(freqVar, Bexc_nT[vComp], label=f'$B^e_{vComp}$', color=Color.BeiFT[vComp],
             linestyle=Style.LS_FT, linewidth=Style.LW_FT) for vComp in xyzComps]
    if FigMisc.MARK_BEXC_MAX:
        iPeak = np.argmax(np.max([Bexc_nT[vComp] for vComp in xyzComps], axis=0))
        xMax = Bexc_nT[xyzComps[0]][iPeak]
        yMax = Bexc_nT[xyzComps[1]][iPeak]
        zMax = Bexc_nT[xyzComps[2]][iPeak]
        compMax = np.argmax([xMax, yMax, zMax])
        peakMax = Bexc_nT[xyzComps[compMax]][iPeak]
        peakLabel = f'$T = \SI{{{freqVar[iPeak]:.2f}}}{{{freqUnits}}}$\n' + \
                    f'$|B_{xyzComps[compMax]}| = \SI{{{peakMax:.2f}}}{{nT}}$'
        hadj = FigMisc.peakLblSize * 8
        vadj = FigMisc.peakLblSize * 2
        ax.annotate(peakLabel, (freqVar[iPeak], peakMax), (-hadj, -vadj), textcoords='offset points',
                    fontsize=FigMisc.peakLblSize, arrowprops={'width': Style.LW_FT*2,
                    'shrink': 0.05, 'headwidth': Style.LW_FT*6, 'headlength': Style.LW_FT*6,
                    'fc': 'black'}, bbox={'boxstyle': 'square', 'fc': 'white', 'ec': 'black',
                    'linewidth': 0.5})

    if Params.LEGEND:
        ax.legend()

    plt.tight_layout()
    fig.savefig(Params.FigureFiles.MagFTexc, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
    log.debug(f'Magnetic excitation spectrum plot saved to file: {Params.FigureFiles.MagFTexc}')
    plt.close()

    return
