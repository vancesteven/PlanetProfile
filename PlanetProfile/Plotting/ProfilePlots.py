import os
import numpy as np
import logging
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Wedge, Rectangle
from scipy.interpolate import interp1d
from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc
from PlanetProfile.Thermodynamics.RefProfiles.RefProfiles import CalcRefProfiles, ReloadRefProfiles
from PlanetProfile.Thermodynamics.HydroEOS import GetPhaseIndices
from PlanetProfile.Thermodynamics.InnerEOS import GetInnerEOS
from PlanetProfile.Utilities.defineStructs import Constants

# Assign logger
log = logging.getLogger('PlanetProfile')

def GeneratePlots(PlanetList, Params):
    
    # Remove latex styling from legend labels if Latex is not installed
    if not FigMisc.TEX_INSTALLED:
        for Planet in PlanetList:
            Planet.label = FigLbl.StripLatexFromString(Planet.label)

    # Handle refprofiles first, so we can print log messages before silencing them
    if Params.PLOT_HYDROSPHERE and not Params.ALL_NO_H2O:
        if Params.CALC_NEW_REF:
            # Calculate reference profiles showing melting curves for
            # several salinities specified in configPP.py
            Params = CalcRefProfiles(PlanetList, Params)
        else:
            # Reload refprofiles for this composition
            Params = ReloadRefProfiles(PlanetList, Params)

    if Params.PLOT_GRAVITY:
        PlotGravPres(PlanetList, Params)
    if Params.PLOT_HYDROSPHERE and np.any([not Planet.Do.NO_H2O for Planet in PlanetList]):
        PlotHydrosphereProps(PlanetList, Params)
    if Params.PLOT_TRADEOFF:
        PlotSilTradeoff(PlanetList, Params)
        if np.any([Planet.Do.Fe_CORE for Planet in PlanetList]):
            PlotCoreTradeoff(PlanetList, Params)
    if Params.PLOT_POROSITY and np.any([Planet.Do.POROUS_ROCK or Planet.Do.POROUS_ICE for Planet in PlanetList]):
        PlotPorosity(PlanetList, Params)
    if Params.PLOT_SEISMIC and Params.CALC_SEISMIC:
        PlotSeismic(PlanetList, Params)
    if Params.PLOT_WEDGE:
        PlotWedge(PlanetList, Params)
    if Params.PLOT_PVT and not Params.SKIP_INNER:
        PlotPvT(PlanetList, Params)

    return


def PlotGravPres(PlanetList, Params):

    fig = plt.figure(figsize=FigSize.vgrav)
    grid = GridSpec(1, 2)
    axes = [fig.add_subplot(grid[0, i]) for i in range(2)]
    if Style.GRIDS:
        [ax.grid() for ax in axes]
        [ax.set_axisbelow(True) for ax in axes]

    axes[0].set_xlabel(FigLbl.gLabel)
    axes[1].set_xlabel(FigLbl.PlabelFull)
    [ax.set_ylabel(FigLbl.rLabel) for ax in axes]
    if Params.ALL_ONE_BODY:
        fig.suptitle(f'{PlanetList[0].name}{FigLbl.gravTitle}')
    else:
        fig.suptitle(FigLbl.gravCompareTitle)

    for Planet in PlanetList:
        legLbl = Planet.label
        if (not Params.ALL_ONE_BODY) and FigLbl.BODYNAME_IN_LABEL:
            legLbl = f'{Planet.name} {legLbl}'
        axes[0].plot(Planet.g_ms2, Planet.r_m[:-1]/1e3,
                     label=legLbl, linewidth=Style.LW_std)
        axes[1].plot(Planet.P_MPa*FigLbl.PmultFull, Planet.r_m[:-1]/1e3,
                     label=legLbl, linewidth=Style.LW_std)

    if Params.LEGEND:
        axes[1].legend()

    plt.tight_layout()
    fig.savefig(Params.FigureFiles.vgrav, format=FigMisc.figFormat, dpi=FigMisc.dpi)
    log.debug(f'Gravity and pressure plot saved to file: {Params.FigureFiles.vgrav}')
    plt.close()

    return


def PlotHydrosphereProps(PlanetList, Params):

    vRow = 1
    if Params.PLOT_SIGS and Params.CALC_CONDUCT:
        DO_SIGS = True
        vRow += 1
    else:
        DO_SIGS = False
        axsigz = None
    if Params.PLOT_SOUNDS and Params.CALC_SEISMIC:
        DO_SOUNDS = True
        vRow += 1
    else:
        DO_SOUNDS = False
        axv = None

    # Generate canvas and add labels
    fig = plt.figure(figsize=FigSize.vhydro)
    grid = GridSpec(vRow, 6)

    axPrho = fig.add_subplot(grid[:, :3])
    axTz = fig.add_subplot(grid[0, 3:])

    axPrho.set_xlabel(FigLbl.rhoLabel)
    axPrho.set_ylabel(FigLbl.PlabelHydro)
    axPrho.invert_yaxis()
    axTz.set_xlabel(FigLbl.Tlabel)
    axTz.set_ylabel(FigLbl.zLabel)
    axTz.invert_yaxis()
    zMax = np.max([Planet.z_m[Planet.Steps.nHydro-1]/1e3 for Planet in PlanetList if not Planet.Do.NO_H2O], initial=0) * 1.05
    axTz.set_ylim([zMax, 0])

    axes = [axPrho, axTz]
    if DO_SIGS:
        axsigz = fig.add_subplot(grid[-1, 3:])
        axsigz.set_xlabel(FigLbl.sigLabel)
        axsigz.set_ylabel(FigLbl.zLabel)
        axsigz.invert_yaxis()
        axsigz.set_ylim([zMax, 0])
        axes.append(axsigz)

    if DO_SOUNDS:
        axv = [fig.add_subplot(grid[1, i]) for i in range(3, 6)]
        axv[0].set_ylabel(FigLbl.zLabel)
        axv[0].set_xlabel(FigLbl.vPoceanLabel)
        axv[1].set_xlabel(FigLbl.vPiceLabel)
        axv[2].set_xlabel(FigLbl.vSiceLabel)
        [ax.invert_yaxis() for ax in axv]
        [ax.set_ylim([zMax, 0]) for ax in axv]
        axes = axes + axv

    if Style.GRIDS:
        [ax.grid() for ax in axes]
        [ax.set_axisbelow(True) for ax in axes]

    if Params.ALL_ONE_BODY:
        fig.suptitle(f'{PlanetList[0].name}{FigLbl.hydroTitle}')
    else:
        fig.suptitle(FigLbl.hydroCompareTitle)

    # Plot reference profiles first, so they plot on bottom of everything
    comps = np.unique([Planet.Ocean.comp for Planet in PlanetList])
    if Params.PLOT_REF:
        # Keep track of which reference profiles have been plotted so that we do each only once
        newRef = {comp:True for comp in comps}

        # Get max pressure among all profiles so we know how far out to plot refs
        Plist = np.concatenate([Planet.P_MPa[:Planet.Steps.nHydro] for Planet in PlanetList])
        Pmax_MPa = np.max(Plist)

        for Planet in PlanetList:
            if newRef[Planet.Ocean.comp] and Planet.Ocean.comp != 'none':
                # Get strings for referencing and labeling
                wList = f'$\\rho_\mathrm{{melt}}$ \ce{{{Planet.Ocean.comp}}} \\{{'
                wList += ', '.join([f'{w*FigLbl.wMult:.0f}' for w in Params.wRef_ppt[Planet.Ocean.comp]])
                wList += '\}\,$\si{' + FigLbl.wUnits + '}$'
                if not FigMisc.TEX_INSTALLED:
                    wList = FigLbl.StripLatexFromString(wList)
                # Take care to only plot the values consistent with layer solutions
                iPlot = Params.Pref_MPa[Planet.Ocean.comp] < Pmax_MPa
                # Plot all reference melting curve densities
                for i in range(Params.nRef[Planet.Ocean.comp]):
                    thisRef, = axPrho.plot(Params.rhoRef_kgm3[Planet.Ocean.comp][i,iPlot],
                                           Params.Pref_MPa[Planet.Ocean.comp][iPlot]*FigLbl.PmultHydro,
                                           color=Color.ref,
                                           lw=Style.LW_ref,
                                           ls=Style.LS_ref[Planet.Ocean.comp])
                    if FigMisc.REFS_IN_LEGEND and i == 0: thisRef.set_label(wList)
                newRef[Planet.Ocean.comp] = False

    wMinMax_ppt = {}
    TminMax_K = {}
    if FigMisc.SCALE_HYDRO_LW or FigMisc.MANUAL_HYDRO_COLORS:
        # Get min and max salinities and temps for each comp for scaling
        for comp in comps:
            if comp != 'none':
                wAll_ppt = [Planet.Ocean.wOcean_ppt for Planet in PlanetList if Planet.Ocean.comp == comp]
                wMinMax_ppt[comp] = [np.min(wAll_ppt), np.max(wAll_ppt)]
                Tall_K = [Planet.Bulk.Tb_K for Planet in PlanetList if Planet.Ocean.comp == comp]
                TminMax_K[comp] = [np.min(Tall_K), np.max(Tall_K)]
                # Reset to default if all models are the same or if desired
                if not FigMisc.RELATIVE_Tb_K or TminMax_K[comp][0] == TminMax_K[comp][1]:
                    TminMax_K[comp] = Color.Tbounds_K

    # Now plot all profiles together
    for Planet in PlanetList:
        # This is a hydrosphere-only plot, so skip waterless bodies
        if Planet.Ocean.comp != 'none':
            legLbl = Planet.label
            if (not Params.ALL_ONE_BODY) and FigLbl.BODYNAME_IN_LABEL:
                legLbl = f'{Planet.name} {legLbl}'

            # Set style options
            if FigMisc.MANUAL_HYDRO_COLORS:
                Color.Tbounds_K = TminMax_K[Planet.Ocean.comp]
                thisColor = Color.cmap[Planet.Ocean.comp](Color.GetNormT(Planet.Bulk.Tb_K))
            else:
                thisColor = None
            if FigMisc.SCALE_HYDRO_LW and wMinMax_ppt[Planet.Ocean.comp][0] != wMinMax_ppt[Planet.Ocean.comp][1]:
                thisLW = Style.GetLW(Planet.Ocean.wOcean_ppt, wMinMax_ppt[Planet.Ocean.comp])
            else:
                thisLW = Style.LW_std

            # Plot density vs. pressure curve for hydrosphere
            axPrho.plot(Planet.rho_kgm3[:Planet.Steps.nHydro],
                        Planet.P_MPa[:Planet.Steps.nHydro]*FigLbl.PmultHydro,
                        label=legLbl, color=thisColor, linewidth=thisLW,
                        linestyle=Style.LS[Planet.Ocean.comp])
            # Plot thermal profile vs. depth in hydrosphere
            therm = axTz.plot(Planet.T_K[:Planet.Steps.nHydro] - FigLbl.Tsub,
                              Planet.z_m[:Planet.Steps.nHydro]/1e3,
                              color=thisColor, linewidth=thisLW,
                              linestyle=Style.LS[Planet.Ocean.comp])
            # Make a dot at the end of the thermal profile
            axTz.scatter(np.max(Planet.T_K[:Planet.Steps.nHydro] - FigLbl.Tsub),
                            np.max(Planet.z_m[:Planet.Steps.nHydro]/1e3),
                            color=therm[-1].get_color(), edgecolors=therm[-1].get_color(),
                            marker=Style.MS_hydro, s=Style.MW_hydro**2*thisLW)

            if DO_SIGS or DO_SOUNDS:
                indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund, \
                indsClath, indsClathWet, _, indsSilLiq, _, _, _, _, _, _ = GetPhaseIndices(Planet.phase)

            if DO_SIGS:
                # Plot electrical conductivity vs. depth for hydrosphere
                axsigz.plot(Planet.sigma_Sm[indsLiq], Planet.z_m[indsLiq]/1e3,
                            color=thisColor, linewidth=thisLW,
                            linestyle=Style.LS[Planet.Ocean.comp])

            if DO_SOUNDS:
                indsIce = np.sort(np.concatenate((indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund,
                                                  indsV, indsVund, indsVI, indsVIund, indsClath, indsClathWet)))
                # Plot sound speeds in ocean and ices vs. depth in hydrosphere
                axv[0].plot(Planet.Seismic.VP_kms[indsLiq], Planet.z_m[indsLiq]/1e3,
                            color=thisColor, linewidth=Style.LW_sound,
                            linestyle=Style.LS[Planet.Ocean.comp])
                axv[1].plot(Planet.Seismic.VP_kms[indsIce], Planet.z_m[indsIce]/1e3,
                            color=thisColor, linewidth=Style.LW_sound,
                            linestyle=Style.LS[Planet.Ocean.comp])
                axv[2].plot(Planet.Seismic.VS_kms[indsIce], Planet.z_m[indsIce]/1e3,
                            color=thisColor, linewidth=Style.LW_sound,
                            linestyle=Style.LS[Planet.Ocean.comp])

    # Limit Tmin so the relevant plot can better show what's going on in the ocean
    axTz.set_xlim(left=FigMisc.TminHydro)

    if Params.LEGEND:
        handles, lbls = axPrho.get_legend_handles_labels()
        axPrho.legend(handles, lbls)
    plt.tight_layout()
    fig.savefig(Params.FigureFiles.vhydro, format=FigMisc.figFormat, dpi=FigMisc.dpi)
    log.debug(f'Hydrosphere plot saved to file: {Params.FigureFiles.vhydro}')
    plt.close()

    return


def PlotCoreTradeoff(PlanetList, Params):

    fig, ax = plt.subplots(1, 1, figsize=FigSize.vcore)
    if Style.GRIDS:
        ax.grid()
        ax.set_axisbelow(True)
    ax.set_xlabel(FigLbl.RsilLabel)
    ax.set_ylabel(FigLbl.RcoreLabel)
    ALL_SAME_CMR2 = np.all([Planet.Bulk.Cmeasured == PlanetList[0].Bulk.Cmeasured for Planet in PlanetList]) \
                and np.all([Planet.Bulk.Cuncertainty == PlanetList[0].Bulk.Cuncertainty for Planet in PlanetList])
    if ALL_SAME_CMR2:
        CMR2str = f', $C/MR^2 = {PlanetList[0].Bulk.Cmeasured}\pm{PlanetList[0].Bulk.Cuncertainty}$'
    else:
        CMR2str = ''
    if Params.ALL_ONE_BODY:
        title = f'{PlanetList[0].name}{FigLbl.coreTitle}{CMR2str}'
    else:
        title = FigLbl.coreCompareTitle + CMR2str
    fig.suptitle(title)


    for Planet in PlanetList:
        if Planet.Do.Fe_CORE:
            if ALL_SAME_CMR2:
                legLbl = Planet.label
            else:
                legLbl = Planet.tradeLabel
            if (not Params.ALL_ONE_BODY) and FigLbl.BODYNAME_IN_LABEL:
                legLbl = f'{Planet.name} {legLbl}'
            ax.plot(Planet.Sil.Rtrade_m/1e3, Planet.Core.Rtrade_m/1e3,
                    label=legLbl, linewidth=Style.LW_std)

    if Params.LEGEND:
        ax.legend()

    plt.tight_layout()
    fig.savefig(Params.FigureFiles.vcore, format=FigMisc.figFormat, dpi=FigMisc.dpi)
    log.debug(f'Core trade plot saved to file: {Params.FigureFiles.vcore}')
    plt.close()

    return


def PlotSilTradeoff(PlanetList, Params):

    fig, ax = plt.subplots(1, 1, figsize=FigSize.vmant)
    if Style.GRIDS:
        ax.grid()
        ax.set_axisbelow(True)

    ax.set_xlabel(FigLbl.RsilLabel)
    ax.set_ylabel(FigLbl.rhoSilLabel)
    ALL_SAME_CMR2 = np.all([Planet.Bulk.Cmeasured == PlanetList[0].Bulk.Cmeasured for Planet in PlanetList]) \
                    and np.all([Planet.Bulk.Cuncertainty == PlanetList[0].Bulk.Cuncertainty for Planet in PlanetList])
    if ALL_SAME_CMR2:
        CMR2str = f', $C/MR^2 = {PlanetList[0].Bulk.Cmeasured}\pm{PlanetList[0].Bulk.Cuncertainty}$'
    else:
        CMR2str = ''
    if Params.ALL_ONE_BODY:
        title = f'{PlanetList[0].name}{FigLbl.mantTitle}{CMR2str}'
    else:
        title = FigLbl.mantCompareTitle + CMR2str
    fig.suptitle(title)

    for Planet in PlanetList:
        if ALL_SAME_CMR2:
            legLbl = Planet.label
        else:
            legLbl = Planet.tradeLabel
        if (not Params.ALL_ONE_BODY) and FigLbl.BODYNAME_IN_LABEL:
            legLbl = f'{Planet.name} {legLbl}'
        ax.plot(Planet.Sil.Rtrade_m/1e3, Planet.Sil.rhoTrade_kgm3,
                label=legLbl, linewidth=Style.LW_std)

    if Params.LEGEND:
        ax.legend()

    plt.tight_layout()
    fig.savefig(Params.FigureFiles.vmant, format=FigMisc.figFormat, dpi=FigMisc.dpi)
    log.debug(f'Mantle trade plot saved to file: {Params.FigureFiles.vmant}')
    plt.close()

    return


def PlotPorosity(PlanetList, Params):

    # Plot dual-axis plot for first entry in PlanetList (usually a main profile), unless we are already doing comparison plots
    if os.path.dirname(Params.FigureFiles.vporeDbl) != 'Comparison':
        Planet = PlanetList[0]
        fig = plt.figure(figsize=FigSize.vpore)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)

        ax.set_xlabel(FigLbl.phiLabel)
        ax.set_ylabel(FigLbl.zLabel)
        ax.invert_yaxis()
        P_from_z = interp1d(Planet.z_m[:-1]/1e3, Planet.P_MPa*FigLbl.PmultFull, bounds_error=False, fill_value='extrapolate')
        z_from_P = interp1d(Planet.P_MPa*FigLbl.PmultFull, Planet.z_m[:-1]/1e3, bounds_error=False, fill_value='extrapolate')
        Pax = ax.secondary_yaxis('right', functions=(P_from_z, z_from_P))
        Pax.set_ylabel(FigLbl.PlabelFull)
        fig.suptitle(f'{Planet.name}{FigLbl.poreTitle}')

        ax.plot(Planet.phi_frac[Planet.phi_frac >= Planet.Sil.phiMin_frac]*FigLbl.phiMult,
                Planet.z_m[:-1][Planet.phi_frac >= Planet.Sil.phiMin_frac]/1e3,
                label=Planet.label, linewidth=Style.LW_std)

        if Params.LEGEND:
            ax.legend()

        fig.savefig(Params.FigureFiles.vporeDbl, format=FigMisc.figFormat, dpi=FigMisc.dpi)
        log.debug(f'Porosity plot (dual axis) saved to file: {Params.FigureFiles.vporeDbl}')
        plt.close()

    # Plot standard config with all passed Planet objects
    fig = plt.figure(figsize=FigSize.vpore)
    grid = GridSpec(1, 2)
    axes = [fig.add_subplot(grid[0, i]) for i in range(2)]
    if Style.GRIDS:
        [ax.grid() for ax in axes]
        [ax.set_axisbelow(True) for ax in axes]

    # Labels and titles
    [ax.set_xlabel(FigLbl.phiLabel) for ax in axes]
    axes[0].set_ylabel(FigLbl.zLabel)
    axes[1].set_ylabel(FigLbl.PlabelFull)
    [ax.invert_yaxis() for ax in axes]
    if Params.ALL_ONE_BODY:
        fig.suptitle(f'{PlanetList[0].name}{FigLbl.poreTitle}')
    else:
        fig.suptitle(FigLbl.poreCompareTitle)

    for Planet in PlanetList:
        if Planet.Do.POROUS_ROCK or Planet.Do.POROUS_ICE:
            legLbl = Planet.label
            if (not Params.ALL_ONE_BODY) and FigLbl.BODYNAME_IN_LABEL:
                legLbl = f'{Planet.name} {legLbl}'
            axes[0].plot(Planet.phi_frac[Planet.phi_frac >= Planet.Sil.phiMin_frac]*FigLbl.phiMult,
                         Planet.z_m[:-1][Planet.phi_frac >= Planet.Sil.phiMin_frac]/1e3,
                         label=legLbl, linewidth=Style.LW_std)
            axes[1].plot(Planet.phi_frac[Planet.phi_frac >= Planet.Sil.phiMin_frac]*FigLbl.phiMult,
                         Planet.P_MPa[Planet.phi_frac >= Planet.Sil.phiMin_frac]*FigLbl.PmultFull,
                         label=legLbl, linewidth=Style.LW_std)

    if Params.LEGEND:
        axes[1].legend()

    plt.tight_layout()
    fig.savefig(Params.FigureFiles.vpore, format=FigMisc.figFormat, dpi=FigMisc.dpi)
    log.debug(f'Porosity plot saved to file: {Params.FigureFiles.vpore}')
    plt.close()

    return


def PlotSeismic(PlanetList, Params):

    fig = plt.figure(figsize=FigSize.vseis)
    grid = GridSpec(2, 2)
    axes = np.array([[fig.add_subplot(grid[i, j]) for j in range(2)] for i in range(2)])
    axf = axes.flatten()
    if Style.GRIDS:
        [ax.grid() for ax in axf]
        [ax.set_axisbelow(True) for ax in axf]

    # Labels and titles
    [ax.set_ylabel(FigLbl.rLabel) for ax in axf]
    axes[0,0].set_xlabel(FigLbl.GSKSlabel)
    axes[0,1].set_xlabel(FigLbl.PTrhoLabel)
    axes[1,0].set_xlabel(FigLbl.vSoundLabel)
    axes[1,1].set_xlabel(FigLbl.QseisLabel)
    axes[1,1].set_xscale('log')
    if Params.ALL_ONE_BODY:
        fig.suptitle(f'{PlanetList[0].name}{FigLbl.seisTitle}')
    else:
        fig.suptitle(FigLbl.seisCompareTitle)

    for Planet in PlanetList:
        if Params.ALL_ONE_BODY:
            legLbl = ''
        else:
            legLbl = Planet.label
            if FigLbl.BODYNAME_IN_LABEL:
                legLbl = f'{Planet.name} {legLbl}'
        axes[0,0].plot(Planet.Seismic.KS_GPa, Planet.r_m[:-1]/1e3,
                       label=legLbl+r' $K_S$', linewidth=Style.LW_seis,
                       linestyle=Style.LS_seis['KS'])
        axes[0,0].plot(Planet.Seismic.GS_GPa, Planet.r_m[:-1]/1e3,
                       label=legLbl+r' $G_S$', linewidth=Style.LW_seis,
                       linestyle=Style.LS_seis['GS'])

        axes[0,1].plot(Planet.P_MPa, Planet.r_m[:-1]/1e3,
                       label=legLbl+r' $P$', linewidth=Style.LW_seis,
                       linestyle=Style.LS_seis['P'])
        axes[0,1].plot(Planet.T_K, Planet.r_m[:-1]/1e3,
                       label=legLbl+r' $T$', linewidth=Style.LW_seis,
                       linestyle=Style.LS_seis['T'])
        axes[0,1].plot(Planet.rho_kgm3, Planet.r_m[:-1]/1e3,
                       label=legLbl+r' $\rho$', linewidth=Style.LW_seis,
                       linestyle=Style.LS_seis['rho'])

        axes[1,0].plot(Planet.Seismic.VP_kms, Planet.r_m[:-1]/1e3,
                       label=legLbl + r' $V_P$', linewidth=Style.LW_seis,
                       linestyle=Style.LS_seis['VP'])
        axes[1,0].plot(Planet.Seismic.VS_kms, Planet.r_m[:-1]/1e3,
                       label=legLbl + r' $V_S$', linewidth=Style.LW_seis,
                       linestyle=Style.LS_seis['VS'])

        axes[1,1].plot(Planet.Seismic.QS, Planet.r_m[:-1]/1e3,
                       label=legLbl + f' ${FigLbl.QseisVar}$', linewidth=Style.LW_seis,
                       linestyle=Style.LS_seis['QS'])

    if Params.LEGEND and np.size(PlanetList) == 1:
        axes[0,0].legend()
        axes[0,1].legend()
        axes[1,0].legend()
        axes[1,1].legend()

    plt.tight_layout()
    fig.savefig(Params.FigureFiles.vseis, format=FigMisc.figFormat, dpi=FigMisc.dpi)
    log.debug(f'Seismic plot saved to file: {Params.FigureFiles.vseis}')
    plt.close()

    return


def PlotWedge(PlanetList, Params):
    """ Plot a wedge diagram showing a visual representation of the major layer structure.
    """

    nWedges = np.size(PlanetList)
    fig = plt.figure(figsize=(FigSize.vwedg[0]*nWedges - (nWedges-1), FigSize.vwedg[1]))
    grid = GridSpec(1, nWedges)
    axes = [fig.add_subplot(grid[0, i]) for i in range(nWedges)]

    # Set plot title based on possible comparison conditions
    if Params.ALL_ONE_BODY:
        title = f'{PlanetList[0].name} {FigLbl.wedgeTitle}'
        if nWedges > 1:
            title = f'{title}s'
            fig.suptitle(f'\\textbf{{{title}}}')

    ang1 = 90 - Style.wedgeAngle_deg
    ang2 = 90 + Style.wedgeAngle_deg

    # Get ionosBounds_km for all bodies without affecting what's in Planet.Magnetic
    if FigMisc.IONOSPHERE_IN_WEDGE:
        ionosUpper_km = np.array([np.max(Planet.Magnetic.ionosBounds_m) / 1e3 if Planet.Magnetic.ionosBounds_m is not None else 0 for Planet in PlanetList])
        ionosLower_km = np.array([np.min(Planet.Magnetic.ionosBounds_m) / 1e3 if np.size(Planet.Magnetic.ionosBounds_m) > 1 else 0 for Planet in PlanetList], dtype=np.float_)
    else:
        ionosUpper_km, ionosLower_km = (np.zeros(nWedges) for _ in range(2))

    # Get largest radius across all wedges
    rMax_km = np.max([ionoTop_km + Planet.Bulk.R_m/1e3 for ionoTop_km, Planet in zip(ionosUpper_km, PlanetList)])

    # Optional boundaries
    if FigMisc.DRAW_CONVECTION_BOUND:
        iceConvBd = Color.wedgeBd
        clathConvBd = Color.wedgeBd
        silConvBd = Color.wedgeBd
    else:
        iceConvBd = Color.iceIcond
        clathConvBd = Color.clathCond
        silConvBd = Color.silCondCmap(1.0)

    # Plot each significant layer for each model, from the outside inward
    for Planet, ax, ionoTop_km, ionoBot_km in zip(PlanetList, axes, ionosUpper_km, ionosLower_km):

        # Construct labels
        if Planet.Do.Fe_CORE:
            if Planet.Core.wS_ppt is None:
                if Planet.Core.wFe_ppt is None:
                    wFeCore_ppt = Constants.wFeDef_ppt
                else:
                    wFeCore_ppt = Planet.Core.wFe_ppt
                wScore_ppt = 1e3 - wFeCore_ppt
            else:
                wScore_ppt = Planet.Core.wS_ppt
            Planet.Core.xS_frac = wScore_ppt / 1e3
            if FigLbl.w_IN_WTPCT:
                xStr = f'{Planet.Core.xS_frac * 1e3 * FigLbl.wMult:.0f}'
            else:
                xStr = f'{Planet.Core.xS_frac * 1e3 * FigLbl.wMult:.2f}'
            coreLine = f'\ce{{Fe}} core with \SI{{{xStr}}}{{{FigLbl.wUnits}}}~\ce{{S}}'
        elif Planet.Sil.EOS is not None and 'undifferentiated' in Planet.Sil.EOS.comp:
            coreLine = 'undifferentiated'
        else:
            coreLine = ''

        if 'Comet' in Planet.Sil.mantleEOS:
            silLine = 'Comet 67P'
        else:
            silLine = f'{Planet.Sil.mantleEOS[:2]} chondrite'

        if Planet.Do.NO_H2O:
            wedgeLabel = f'{silLine}\n{coreLine}\n$q_\mathrm{{surf}}$~\SI{{{Planet.Bulk.qSurf_Wm2*1e3}}}{{{FigLbl.fluxUnits}}}'
        else:
            if Planet.Ocean.comp == 'PureH2O':
                compStr = r'Pure \ce{H2O} ocean'
            else:
                compStr = f'\SI{{{Planet.Ocean.wOcean_ppt:.1f}}}{{{FigLbl.wUnits}}}~\ce{{{Planet.Ocean.comp}}}'
            wedgeLabel = f'{silLine} mantle\n{coreLine}\n{compStr}, $z_b$~\SI{{{Planet.zb_km:.1f}}}{{km}}'

        if Planet.Do.POROUS_ROCK:
            wedgeLabel = f'Porous {wedgeLabel}'

        if Params.ALL_ONE_BODY and not nWedges == 1:
            indivTitle = wedgeLabel
        else:
            indivTitle = f'\\textbf{{{Planet.name}}}\n{wedgeLabel}'

        if not FigMisc.TEX_INSTALLED:
            indivTitle = FigLbl.StripLatexFromString(indivTitle)
        ax.set_title(indivTitle)
        R_km = Planet.Bulk.R_m / 1e3
        rTicks = []
        rTickRefs = []

        # @@@@@@@@@@
        # Ionosphere
        # @@@@@@@@@@
        if FigMisc.IONOSPHERE_IN_WEDGE and ionoTop_km != 0:
            # Ionosphere gradient layers
            dzIonos_km = ionoTop_km - ionoBot_km
            ionosGrad, dz = np.linspace(0, 1, Color.ionoN+1, retstep=True)
            # Outer boundary around ionosphere to use as clip path
            ionosOuter = ax.add_patch(Wedge((0.5,0), (R_km + ionoTop_km)/rMax_km, ang1, ang2,
                               width=dzIonos_km/rMax_km,
                               fc=Color.none, lw=Style.LW_wedgeMajor, ec=Color.none))
            for thisIonosFrac in ionosGrad[:-1]:
                ax.add_patch(Wedge((0.5, 0), (R_km + ionoTop_km)/rMax_km - thisIonosFrac*dzIonos_km/rMax_km, ang1, ang2,
                                   width=dz*dzIonos_km/rMax_km, clip_path=ionosOuter,
                                   fc=Color.ionoCmap(thisIonosFrac), ec=Color.ionoCmap(thisIonosFrac)))

            # Draw outer boundary for ionosphere
            if FigMisc.DRAW_IONOS_BOUND:
                ionosOuter.set_edgecolor(Color.wedgeBd)
                ax.add_patch(ionosOuter)

        if not Planet.Do.NO_H2O:
            # @@@@@@@@@
            # Ice shell
            # @@@@@@@@@
            rTicks.append(R_km)
            rTickRefs.append(R_km/rMax_km)

            # Starting with ice I or clathrates
            if Planet.Do.CLATHRATE:
                if Planet.Bulk.clathType == 'top' or Planet.Bulk.clathType == 'whole':
                    # Clathrates at the surface in this case
                    ax.add_patch(Wedge((0.5,0), R_km/rMax_km, ang1, ang2,
                                       width=Planet.eLid_m/1e3/rMax_km,
                                       fc=Color.clathCond, lw=Style.LW_wedge, ec=clathConvBd))
                    if Planet.Bulk.clathType == 'top':
                        # Ice I boundary line
                        ax.add_patch(Wedge((0.5,0), (R_km - Planet.zIceI_m/1e3)/rMax_km, ang1, ang2,
                                           width=Planet.dzIceI_km/rMax_km,
                                           fc=Color.none, lw=Style.LW_wedgeMajor, ec=iceConvBd))
                        # Conductive ice I underneath clathrates
                        if Planet.zIceI_m < Planet.eLid_m:
                            ax.add_patch(Wedge((0.5,0), (R_km - Planet.zIceI_m/1e3)/rMax_km, ang1, ang2,
                                               width=(Planet.dzIceI_km - (Planet.Dconv_m + Planet.deltaTBL_m)/1e3)/rMax_km,
                                               fc=Color.iceIcond, lw=Style.LW_wedge, ec=iceConvBd))
                        # Convective ice I underneath clathrates
                        if (Planet.Dconv_m + Planet.deltaTBL_m) > 0:
                            ax.add_patch(Wedge((0.5,0), (R_km - Planet.zIceI_m/1e3)/rMax_km, ang1, ang2,
                                               width=(Planet.Dconv_m + Planet.deltaTBL_m)/1e3/rMax_km,
                                               fc=Color.iceIconv, lw=Style.LW_wedge, ec=iceConvBd))
                    else:
                        # Convective clathrates
                        if (Planet.Dconv_m + Planet.deltaTBL_m) > 0:
                            ax.add_patch(Wedge((0.5,0), (R_km - Planet.eLid_m/1e3)/rMax_km, ang1, ang2,
                                               width=(Planet.Dconv_m + Planet.deltaTBL_m)/1e3/rMax_km,
                                               fc=Color.clathConv, lw=Style.LW_wedge, ec=clathConvBd))
                else:
                    # Clathrates in an underplate in this case, always conductive
                    # Conductive ice I at the surface
                    ax.add_patch(Wedge((0.5, 0), R_km/rMax_km, ang1, ang2,
                                       width=Planet.dzIceI_km/rMax_km,
                                       fc=Color.iceIcond, lw=Style.LW_wedge, ec=iceConvBd))
                    # Clathrate underplate
                    ax.add_patch(Wedge((0.5,0), (R_km - Planet.zClath_km)/rMax_km, ang1, ang2,
                                       width=Planet.dzClath_km/rMax_km,
                                       fc=Color.clathCond, lw=Style.LW_wedge, ec=clathConvBd))
                    
                # Outer boundary around clathrates
                ax.add_patch(Wedge((0.5,0), R_km/rMax_km, ang1, ang2,
                                   width=Planet.dzClath_km/rMax_km,
                                   fc=Color.none, lw=Style.LW_wedgeMajor, ec=Color.wedgeBd))
            else:
                # Ice Ih at the surface in this case
                # Conductive ice I
                ax.add_patch(Wedge((0.5,0), R_km/rMax_km, ang1, ang2,
                                   width=Planet.eLid_m/1e3/rMax_km,
                                   fc=Color.iceIcond, lw=Style.LW_wedge, ec=iceConvBd))
                # Convective ice I
                if (Planet.Dconv_m + Planet.deltaTBL_m) > 0:
                    ax.add_patch(Wedge((0.5,0), (R_km - Planet.eLid_m/1e3)/rMax_km, ang1, ang2,
                                       width=(Planet.Dconv_m + Planet.deltaTBL_m)/1e3/rMax_km,
                                       fc=Color.iceIconv, lw=Style.LW_wedge, ec=iceConvBd))
            
            # Outer boundary around ice I
            if Planet.dzIceI_km > 0:
                ax.add_patch(Wedge((0.5,0), (R_km - Planet.zIceI_m/1e3)/rMax_km, ang1, ang2,
                                   width=Planet.dzIceI_km/rMax_km,
                                   fc=Color.none, lw=Style.LW_wedgeMajor, ec=Color.wedgeBd))
                
            # Surface HP ices
            if Planet.dzIceIII_km > 0:
                ax.add_patch(Wedge((0.5,0), (R_km - Planet.zIceIII_m/1e3)/rMax_km, ang1, ang2,
                                   width=Planet.dzIceIII_km/rMax_km,
                                   fc=Color.iceIII, lw=Style.LW_wedgeMajor, ec=Color.wedgeBd))
            if Planet.dzIceVund_km > 0:
                ax.add_patch(Wedge((0.5,0), (R_km - Planet.zIceVund_m/1e3)/rMax_km, ang1, ang2,
                                   width=Planet.dzIceVund_km/rMax_km,
                                   fc=Color.iceV, lw=Style.LW_wedgeMajor, ec=Color.wedgeBd))
            
            # @@@@@@@@@@@
            # Ocean layer
            # @@@@@@@@@@@
            if Planet.D_km > 0:
                if FigMisc.WEDGE_ICE_TICKS or (Planet.zb_km/R_km >= FigMisc.minzbRratio_frac and Planet.D_km/R_km >= FigMisc.minzbRratio_frac):
                    rTicks.append(R_km - Planet.zb_km)
                    rTickRefs.append((R_km - Planet.zb_km)/rMax_km)

                # Ocean outer boundary to use as clip path
                oceanOuter = ax.add_patch(Wedge((0.5,0), (R_km - Planet.zb_km)/rMax_km, ang1, ang2,
                                                width=Planet.D_km/rMax_km,
                                                fc=Color.none, lw=Style.LW_wedgeMajor, ec=Color.none))

                # Ocean gradient layers
                oceanGrad, dz = np.linspace(0, 1, Color.oceanN+1, retstep=True)
                for thisOceanFrac in oceanGrad[:-1]:
                    ax.add_patch(Wedge((0.5, 0), ((R_km - Planet.zb_km) - thisOceanFrac*Planet.D_km)/rMax_km, ang1, ang2,
                                       width=dz*Planet.D_km/rMax_km, clip_path=oceanOuter,
                                       fc=Color.oceanCmap(thisOceanFrac), ec=Color.oceanCmap(thisOceanFrac)))
    
                # Draw outer boundary
                oceanOuter.set_edgecolor(Color.wedgeBd)
                ax.add_patch(oceanOuter)
                    
            # Undersea HP ices
            if Planet.dzIceV_km > 0:
                ax.add_patch(Wedge((0.5,0), (R_km - Planet.zIceV_m/1e3)/rMax_km, ang1, ang2,
                                   width=Planet.dzIceV_km/rMax_km,
                                   fc=Color.iceV, lw=Style.LW_wedgeMajor, ec=Color.wedgeBd))
            if Planet.dzIceVI_km > 0:
                ax.add_patch(Wedge((0.5,0), (R_km - Planet.zIceVI_m/1e3)/rMax_km, ang1, ang2,
                                   width=Planet.dzIceVI_km/rMax_km,
                                   fc=Color.iceVI, lw=Style.LW_wedgeMajor, ec=Color.wedgeBd))
             
        # @@@@@@@@@
        # Silicates
        # @@@@@@@@@
        rTicks.append(Planet.Sil.Rmean_m/1e3)
        rTickRefs.append(Planet.Sil.Rmean_m/1e3/rMax_km)

        # Outer boundary around silicate layer to use as clip path for gradients
        silOuter = ax.add_patch(Wedge((0.5,0), Planet.Sil.Rmean_m/1e3/rMax_km, ang1, ang2,
                                      width=(Planet.Sil.Rmean_m - Planet.Core.Rmean_m)/rMax_km,
                                      fc=Color.none, lw=Style.LW_wedgeMajor, ec=Color.none))

        if Planet.Do.POROUS_ROCK:
            # Outer boundary around porous silicate layer to use ase clip path
            silPorousOuter = ax.add_patch(Wedge((0.5,0), Planet.Sil.Rmean_m/1e3/rMax_km, ang1, ang2,
                                          width=Planet.dzSilPorous_km/rMax_km,
                                          fc=Color.none, lw=Style.LW_wedgeMajor, ec=Color.none))

            # Porous rock gradient
            porousGrad, dz = np.linspace(0, 1, Color.silPorousN+1, retstep=True)
            for thisPorousFrac in porousGrad[:-1]:
                ax.add_patch(Wedge((0.5,0), (Planet.Sil.Rmean_m/1e3 - thisPorousFrac*Planet.dzSilPorous_km)/rMax_km, ang1, ang2,
                                   width=dz*Planet.dzSilPorous_km/rMax_km, clip_path=silPorousOuter,
                                   fc=Color.silPorousCmap(thisPorousFrac), ec=Color.silPorousCmap(thisPorousFrac)))

            # Draw outer boundary around porous silicate layer
            if FigMisc.DRAW_POROUS_BOUND:
                silPorousOuter.set_edgecolor(Color.wedgeBd)
                ax.add_patch(silPorousOuter)

        # Only conductive silicates are currently modeled, no convective
        # Conductive silicate gradient
        silCondGrad, dz = np.linspace(0, 1, Color.silCondN+1, retstep=True)
        dzSilCond_km = (Planet.Sil.Rmean_m - Planet.Core.Rmean_m) / 1e3 - Planet.dzSilPorous_km
        for thisSilFrac in silCondGrad[:-1]:
            ax.add_patch(Wedge((0.5, 0), (Planet.Sil.Rmean_m/1e3 - Planet.dzSilPorous_km - thisSilFrac*dzSilCond_km)/rMax_km, ang1, ang2,
                               width=dz*dzSilCond_km/rMax_km, clip_path=silOuter,
                               fc=Color.silCondCmap(thisSilFrac), ec=Color.silCondCmap(thisSilFrac)))

        # Draw outer boundary
        silOuter.set_edgecolor(Color.wedgeBd)
        ax.add_patch(silOuter)

        # @@@@@@@@@
        # Iron core
        # @@@@@@@@@
        if Planet.Do.Fe_CORE:
            rTicks.append(Planet.Core.Rmean_m/1e3)
            rTickRefs.append(Planet.Core.Rmean_m/1e3/rMax_km)

            # FeS layer
            if FigMisc.DRAW_FeS_BOUND:
                FeSbd = Color.wedgeBd
            else:
                FeSbd = Color.FeS
            if Planet.dzFeS_km > 0:
                ax.add_patch(Wedge((0.5,0), Planet.Core.Rmean_m/1e3/rMax_km, ang1, ang2,
                               width=Planet.dzFeS_km/rMax_km,
                               fc=Color.FeS, lw=Style.LW_wedge, ec=FeSbd))

            # Remaining iron (pure or mixed)
            ax.add_patch(Wedge((0.5,0), Planet.Core.Rmean_m/1e3/rMax_km, ang1, ang2,
                               fc=Color.Fe, lw=Style.LW_wedge, ec=FeSbd))

        # Outer boundary around core layer
        ax.add_patch(Wedge((0.5,0), Planet.Core.Rmean_m/1e3/rMax_km, ang1, ang2,
                           fc=Color.none, lw=Style.LW_wedgeMajor, ec=Color.wedgeBd))

        # Adjust plots to look nice
        ax.set_yticks(rTickRefs)
        ax.set_yticklabels(np.array(rTicks, dtype=np.int_))
        [ax.spines[side].set_visible(False) for side in ['top', 'right', 'bottom']]
        ax.get_xaxis().set_visible(False)
        ax.set_ylabel(FigLbl.wedgeRadius)
        ax.set_aspect('equal')

    fig.tight_layout()
    fig.savefig(Params.FigureFiles.vwedg, format=FigMisc.figFormat, dpi=FigMisc.dpi)
    log.debug(f'Wedge plot saved to file: {Params.FigureFiles.vwedg}')
    plt.close()

    return


def PlotPvT(PlanetList, Params):

    # Unlike most other plotting routines, here we can't plot multiple bodies together.
    # Just make these plots for the first model, which is the primary.
    if os.path.dirname(Params.FigureFiles.vpvt) != 'Comparison':
        Planet = PlanetList[0]
        if Planet.Sil.EOS is None:
            Planet.Sil.EOS = GetInnerEOS(Planet.Sil.mantleEOS, EOSinterpMethod=Params.lookupInterpMethod,
                        kThermConst_WmK=Planet.Sil.kTherm_WmK, HtidalConst_Wm3=Planet.Sil.Htidal_Wm3,
                        porosType=Planet.Sil.porosType, phiTop_frac=Planet.Sil.phiRockMax_frac,
                        Pclosure_MPa=Planet.Sil.Pclosure_MPa, phiMin_frac=Planet.Sil.phiMin_frac,
                        EXTRAP=Params.EXTRAP_SIL)
        INCLUDING_CORE = FigMisc.PVT_INCLUDE_CORE and Planet.Do.Fe_CORE
        if INCLUDING_CORE and Planet.Core.EOS is None:
            Planet.Core.EOS = GetInnerEOS(Planet.Core.coreEOS, EOSinterpMethod=Params.lookupInterpMethod, Fe_EOS=True,
                        kThermConst_WmK=Planet.Core.kTherm_WmK, EXTRAP=Params.EXTRAP_Fe,
                        wFeCore_ppt=Planet.Core.wFe_ppt, wScore_ppt=Planet.Core.wS_ppt)

        # Check that it's worth converting to GPa if that setting has been selected -- reset labels if not
        if Planet.P_MPa[-1] < 100 and FigLbl.PFULL_IN_GPa:
            log.debug('FigLbl.PFULL_IN_GPa is True, but Pmax is less than 0.1 GPa. Pressures will be plotted in MPa.')
            FigLbl.PFULL_IN_GPa = False
            FigLbl.SetUnits()

        fig = plt.figure(figsize=FigSize.vpvt)
        grid = GridSpec(2, 4)
        axes = np.array([[fig.add_subplot(grid[i, j]) for j in range(4)] for i in range(2)])
        axes[0,1].set_axis_off()
        axf = axes[axes != axes[0,1]].flatten()
        if Style.GRIDS:
            [ax.grid() for ax in axf]
            [ax.set_axisbelow(False) for ax in axf]
        # Labels and titles
        [ax.set_xlabel(FigLbl.Tlabel) for ax in axes[1, :]]
        [ax.set_ylabel(FigLbl.PlabelFull) for ax in axes[:, 0]]
        [ax.invert_yaxis() for ax in axf]
        axes[0,0].set_title(FigLbl.rhoLabel)
        axes[1,0].set_title(FigLbl.CpLabel)
        axes[1,1].set_title(FigLbl.alphaLabel)
        # axes[1,1].set_title()  # This axes intentionally left blank.
        axes[0,2].set_title(FigLbl.VPlabel)
        axes[1,2].set_title(FigLbl.VSlabel)
        axes[0,3].set_title(FigLbl.KSlabel)
        axes[1,3].set_title(FigLbl.GSlabel)
        # Get pressure and temperature bounds and eval pts and set title

        iSil = np.logical_and(Planet.phase >= Constants.phaseSil,
                              Planet.phase < Constants.phaseSil + 10)
        if INCLUDING_CORE:
            fig.suptitle(f'{Planet.name}{FigLbl.PvTtitleCore}')
            iCore = Planet.phase >= Constants.phaseFe
            iInner = np.logical_or(iSil, iCore)
        else:
            fig.suptitle(f'{Planet.name}{FigLbl.PvTtitleSil}')
            iCore = np.zeros_like(Planet.phase).astype(bool)
            iInner = iSil
        Pgeo = Planet.P_MPa[iInner]
        Tgeo = Planet.T_K[iInner]
        Pmin = np.min(Pgeo)
        Pmax = np.max(Pgeo)
        Tmin = np.min(Tgeo)
        Tmax = np.max(Tgeo)
        if INCLUDING_CORE:
            nPsil = FigMisc.nPgeo - FigMisc.nPgeoCore
            nPcore = FigMisc.nPgeoCore
            Pcore = np.linspace(np.min(Planet.P_MPa[iCore]), Pmax, nPcore)
        else:
            nPsil = FigMisc.nPgeo
            Pcore = np.empty(0)
        Psil = np.linspace(Pmin, np.max(Planet.P_MPa[iSil]), nPsil)
        Tinner = np.linspace(Tmin, Tmax, FigMisc.nTgeo)

        # Get data to plot
        rhoSil = Planet.Sil.EOS.fn_rho_kgm3(Psil, Tinner, grid=True)
        CpSil = Planet.Sil.EOS.fn_Cp_JkgK(Psil, Tinner, grid=True)
        alphaSil = Planet.Sil.EOS.fn_alpha_pK(Psil, Tinner, grid=True)
        VPsil = Planet.Sil.EOS.fn_VP_kms(Psil, Tinner, grid=True)
        VSsil = Planet.Sil.EOS.fn_VS_kms(Psil, Tinner, grid=True)
        KSsil = Planet.Sil.EOS.fn_KS_GPa(Psil, Tinner, grid=True)
        GSsil = Planet.Sil.EOS.fn_GS_GPa(Psil, Tinner, grid=True)

        # Plot colormaps of Perple_X data
        PsilScaled = Psil * FigLbl.PmultFull
        rhoPlotSil =   axes[0,0].pcolormesh(Tinner, PsilScaled, rhoSil, cmap=Color.PvTsilCmap)
        CpPlotSil =    axes[1,0].pcolormesh(Tinner, PsilScaled, CpSil, cmap=Color.PvTsilCmap)
        alphaPlotSil = axes[1,1].pcolormesh(Tinner, PsilScaled, alphaSil, cmap=Color.PvTsilCmap)
        VPplotSil =    axes[0,2].pcolormesh(Tinner, PsilScaled, VPsil, cmap=Color.PvTsilCmap)
        VSplotSil =    axes[1,2].pcolormesh(Tinner, PsilScaled, VSsil, cmap=Color.PvTsilCmap)
        KSplotSil =    axes[0,3].pcolormesh(Tinner, PsilScaled, KSsil, cmap=Color.PvTsilCmap)
        GSplotSil =    axes[1,3].pcolormesh(Tinner, PsilScaled, GSsil, cmap=Color.PvTsilCmap)

        if INCLUDING_CORE:
            rhoCore = Planet.Core.EOS.fn_rho_kgm3(Pcore, Tinner, grid=True)
            CpCore = Planet.Core.EOS.fn_Cp_JkgK(Pcore, Tinner, grid=True)
            alphaCore = Planet.Core.EOS.fn_alpha_pK(Pcore, Tinner, grid=True)
            VPcore = Planet.Core.EOS.fn_VP_kms(Pcore, Tinner, grid=True)
            VScore = Planet.Core.EOS.fn_VS_kms(Pcore, Tinner, grid=True)
            KScore = Planet.Core.EOS.fn_KS_GPa(Pcore, Tinner, grid=True)
            GScore = Planet.Core.EOS.fn_GS_GPa(Pcore, Tinner, grid=True)

            # Plot colormaps of core data
            PcoreScaled = Pcore * FigLbl.PmultFull
            rhoPlotCore =   axes[0,0].pcolormesh(Tinner, PcoreScaled, rhoCore, cmap=Color.PvTcoreCmap)
            CpPlotCore =    axes[1,0].pcolormesh(Tinner, PcoreScaled, CpCore, cmap=Color.PvTcoreCmap)
            alphaPlotCore = axes[1,1].pcolormesh(Tinner, PcoreScaled, alphaCore, cmap=Color.PvTcoreCmap)
            VPplotCore =    axes[0,2].pcolormesh(Tinner, PcoreScaled, VPcore, cmap=Color.PvTcoreCmap)
            VSplotCore =    axes[1,2].pcolormesh(Tinner, PcoreScaled, VScore, cmap=Color.PvTcoreCmap)
            KSplotCore =    axes[0,3].pcolormesh(Tinner, PcoreScaled, KScore, cmap=Color.PvTcoreCmap)
            GSplotCore =    axes[1,3].pcolormesh(Tinner, PcoreScaled, GScore, cmap=Color.PvTcoreCmap)

            # Add core colorbars for each plot
            cbarsCore = [
                fig.colorbar(rhoPlotCore, ax=axes[0,0]),
                fig.colorbar(CpPlotCore, ax=axes[1,0]),
                fig.colorbar(alphaPlotCore, ax=axes[1,1]),
                fig.colorbar(VPplotCore, ax=axes[0,2]),
                fig.colorbar(VSplotCore, ax=axes[1,2]),
                fig.colorbar(KSplotCore, ax=axes[0,3]),
                fig.colorbar(GSplotCore, ax=axes[1,3])
            ]

            if FigLbl.PVT_CBAR_LABELS:
                [cbar.ax.set_title(FigLbl.core) for cbar in cbarsCore]

        # Add colorbars for each silicate plot (second, so that silicate bar appears closest to the plot)
        cbarsSil = [
            fig.colorbar(rhoPlotSil, ax=axes[0,0]),
            fig.colorbar(CpPlotSil, ax=axes[1,0]),
            fig.colorbar(alphaPlotSil, ax=axes[1,1]),
            fig.colorbar(VPplotSil, ax=axes[0,2]),
            fig.colorbar(VSplotSil, ax=axes[1,2]),
            fig.colorbar(KSplotSil, ax=axes[0,3]),
            fig.colorbar(GSplotSil, ax=axes[1,3])
        ]

        if FigLbl.PVT_CBAR_LABELS:
            [cbar.ax.set_title(FigLbl.sil) for cbar in cbarsSil]

        # Plot geotherm on top of colormaps
        [ax.plot(Tgeo, Pgeo * FigLbl.PmultFull, linewidth=Style.LW_geotherm, linestyle=Style.LS_geotherm,
                 color=Color.geotherm) for ax in axf]

        plt.tight_layout()
        fig.savefig(Params.FigureFiles.vpvt, format=FigMisc.figFormat, dpi=FigMisc.dpi)
        log.debug(f'Silicate/core PT properties plot saved to file: {Params.FigureFiles.vpvt}')
        plt.close()

    return


def PlotExploreOgram(ExplorationList, Params):
    """ For plotting points showing the various models used in making
        exploreogram plots.
    """

    FigLbl.SetExploration(ExplorationList[0].bodyname, ExplorationList[0].xName,
                          ExplorationList[0].yName, ExplorationList[0].zName)

    for Exploration in ExplorationList:
        fig, ax = plt.subplots(1, 1, figsize=FigSize.explore)
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)
    
        fig.suptitle(FigLbl.explorationTitle)
        ax.set_xlabel(FigLbl.xLabelExplore)
        ax.set_ylabel(FigLbl.yLabelExplore)
        ax.set_xscale(FigLbl.xScaleExplore)
        ax.set_yscale(FigLbl.yScaleExplore)
    
        x = Exploration.__getattribute__(Exploration.xName) * FigLbl.xMultExplore
        y = Exploration.__getattribute__(Exploration.yName) * FigLbl.yMultExplore
        z = Exploration.__getattribute__(Exploration.zName) * FigLbl.zMultExplore
        ax.set_xlim([np.min(x), np.max(x)])
        ax.set_ylim([np.min(y), np.max(y)])
        mesh = ax.pcolormesh(x, y, z, shading='auto', cmap=Color.cmap['default'])
        cont = ax.contour(x, y, z, colors='black')
        lbls = plt.clabel(cont, fmt='%1.0f')
        cbar = fig.colorbar(mesh, ax=ax)
        # Add the min and max values to the colorbar for reading convenience
        # We compare z values to z values to exclude nans from the max finding,
        # exploiting the fact that nan == nan is False.
        zValid = z[z == z]
        if np.size(zValid) > 0:
            new_ticks = np.insert(np.append(cbar.get_ticks(), np.max(zValid)), 0, np.min(zValid))
            cbar.set_ticks(np.unique(new_ticks))
        cbar.set_label(FigLbl.cbarLabelExplore, size=12)

        plt.tight_layout()
        fig.savefig(Params.FigureFiles.explore, format=FigMisc.figFormat, dpi=FigMisc.dpi)
        log.debug(f'Plot saved to file: {Params.FigureFiles.explore}')
        plt.close()

    # Plot combination
    if Params.COMPARE and np.size(ExplorationList) > 1:
        fig, ax = plt.subplots(1, 1, figsize=FigSize.explore)
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)

        fig.suptitle(FigLbl.explorationTitle)
        ax.set_xlabel(FigLbl.xLabelExplore)
        ax.set_ylabel(FigLbl.yLabelExplore)
        ax.set_xscale(FigLbl.xScaleExplore)
        ax.set_yscale(FigLbl.yScaleExplore)

        x = ExplorationList[0].__getattribute__(ExplorationList[0].xName) * FigLbl.xMult
        y = ExplorationList[0].__getattribute__(ExplorationList[0].yName) * FigLbl.yMult
        z = ExplorationList[0].__getattribute__(ExplorationList[0].zName) * FigLbl.zMult
        for Exploration in ExplorationList[1:]:
            x = np.append(x, Exploration.__getattribute__(Exploration.xName)) * FigLbl.xMult
            y = np.append(y, Exploration.__getattribute__(Exploration.yName)) * FigLbl.yMult
            z = np.append(z, Exploration.__getattribute__(Exploration.zName)) * FigLbl.zMult
        mesh = ax.pcolormesh(x, y, z, shading='auto', cmap=Color.cmap['default'])
        cont = ax.contour(x, y, z, colors='black')
        lbls = plt.clabel(cont, fmt='%1.0f')
        cbar = fig.colorbar(mesh, ax=ax)
        # Append the max value to the colorbar for reading convenience
        # We compare z values to z values to exclude nans from the max finding,
        # exploiting the fact that nan == nan is False.
        cbar.set_ticks(np.append(cbar.get_ticks(), np.max(z[z == z])))
        cbar.set_label(FigLbl.cbarLabelExplore, size=12)

        plt.tight_layout()
        fig.savefig(Params.FigureFiles.explore, format=FigMisc.figFormat, dpi=FigMisc.dpi)
        log.debug(f'Plot saved to file: {Params.FigureFiles.explore}')
        plt.close()

    return


