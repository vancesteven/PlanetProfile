import os
import numpy as np
import logging
from collections.abc import Iterable
import matplotlib.pyplot as plt
import matplotlib.colorbar as mcbar
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Wedge, Rectangle
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp1d
from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc, InductParams as IndParams
from PlanetProfile.Thermodynamics.RefProfiles.RefProfiles import CalcRefProfiles, ReloadRefProfiles
from PlanetProfile.Thermodynamics.HydroEOS import GetPhaseIndices
from PlanetProfile.Thermodynamics.InnerEOS import GetInnerEOS
from PlanetProfile.Utilities.SetupInit import SetupFilenames
from PlanetProfile.Utilities.defineStructs import Constants
from PlanetProfile.MagneticInduction.Moments import Excitations
from MoonMag.asymmetry_funcs import getMagSurf as GetMagSurf

# Assign logger
log = logging.getLogger('PlanetProfile')

def GeneratePlots(PlanetList, Params):

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


def GenerateMagPlots(PlanetList, Params):

    if Params.PLOT_BDIP:
        PlotComplexBdip(PlanetList, Params)
    if Params.PLOT_MAG_SPECTRUM and np.any([Planet.Magnetic.FT_LOADED for Planet in PlanetList]):
        PlotMagSpectrum(PlanetList, Params)
    if Params.PLOT_BSURF:
        PlotMagSurface(PlanetList, Params)
    if Params.PLOT_ASYM and Params.CALC_ASYM:
        PlotAsym(PlanetList, Params)

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
        fig, ax = plt.subplots(1, 1, figsize=FigSize.vpore)
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


def PlotInductOgramPhaseSpace(InductionList, Params):
    """ For plotting points showing the various models used in making
        inductogram plots.
    """

    if InductionList[0].SINGLE_COMP:
        FigLbl.singleComp(InductionList[0].comps[0])
    FigLbl.SetInduction(InductionList[0].bodyname, Params.Induct, InductionList[0].Texc_hr.values())

    sigma_Sm, D_km, ptColors = (np.empty_like(InductionList) for _ in range(3))
    for i, Induction in enumerate(InductionList):
        sigma_Sm[i] = Induction.sigmaMean_Sm.flatten()
        D_km[i] = Induction.D_km.flatten()
        if Params.Induct.inductOtype == 'sigma':
            # In this case, we likely don't have salinity and ocean temp information
            # so we need to set the colormap to use the info we do have
            sigmaNorm = sigma_Sm[i] / 10**Params.Induct.sigmaMax[Induction.bodyname]
            Dnorm = D_km[i] / np.max(D_km)
            ptColors[i] = Color.OceanCmap(Induction.compsList, sigmaNorm, Dnorm,
                                          DARKEN_SALINITIES=FigMisc.DARKEN_SALINITIES)
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
                ptColors[i] = Color.OceanCmap(Induction.compsList, w_normFrac, Tmean_normFrac,
                                              DARKEN_SALINITIES=FigMisc.DARKEN_SALINITIES)
            elif Params.Induct.colorType == 'zb':
                zb_km = Induction.zb_km.flatten()
                zb_normFrac = interp1d([np.min(zb_km), np.max(zb_km)], [0.0, 1.0])(zb_km)
                ptColors[i] = Color.OceanCmap(Induction.compsList, w_normFrac, zb_normFrac,
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
        cbarUnits = InductionList[0].zb_km.flatten()
        cbarLabel = FigLbl.iceThickLbl
    else:
        comps = np.unique(InductionList[0].comps)
        w_ppt = InductionList[0].x.flatten()
        yFlat = InductionList[0].y.flatten()
        if Params.Induct.colorType == 'Tmean':
            cbarUnits = InductionList[0].Tmean_K.flatten()
            cbarLabel = FigLbl.oceanTempLbl
        elif Params.Induct.colorType == 'zb':
            cbarUnits = InductionList[0].zb_km.flatten()
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
        axes[1].scatter(w_ppt, yFlat, s=Style.MW_Induction,
                        marker=Style.MS_Induction, c=ptColors[0])

    # Labels and titles
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
            thisComp = InductionList[0].compsList == comp
            pts[comp] = axes[0].scatter(sigma_Sm[0][thisComp], D_km[0][thisComp], s=Style.MW_Induction,
                            marker=Style.MS_Induction, c=cbarUnits[thisComp], cmap=Color.cmap[comp])
            cbar[comp] = fig.colorbar(pts[comp], cax=cbarAx)
            cbarAx.set_title(f'\ce{{{comp}}}', fontsize=FigMisc.cbarTitleSize)

    cbar[comps[-1]].set_label(cbarLabel)
    fig.savefig(Params.FigureFiles.phaseSpace, format=FigMisc.figFormat, dpi=FigMisc.dpi)
    log.debug(f'InductOgram phase space plot saved to file: {Params.FigureFiles.phaseSpace}')
    plt.close()

    # Plot combination
    if Params.COMPARE and np.size(InductionList) > 1 and Params.Induct.inductOtype != 'sigma':
        comps = np.unique(np.append([],[Induction.comps for Induction in InductionList]))
        nComps = np.size(comps)
        figWidth = FigSize.phaseSpaceSolo[0] + nComps * FigMisc.cbarSpace
        fig, ax = plt.subplots(1, 1, figsize=(figWidth, FigSize.phaseSpaceSolo[1]))
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)

        # Labels and titles
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
            cbar = mcbar.ColorbarBase(cbarAx, cmap=Color.cmap[comp], format=FigMisc.cbarFmt,
                                      values=np.linspace(np.min(comboCbarUnits[thisComp]), np.max(comboCbarUnits[thisComp]), FigMisc.nCbarPts))
            fig.add_axes(cbarAx)
            cbarAx.set_title(f'\ce{{{comp}}}', fontsize=FigMisc.cbarTitleSize)

        cbar.set_label(cbarLabel, size=12)
        plt.tight_layout()
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
    FigLbl.SetInduction(Induction.bodyname, Params.Induct, Induction.Texc_hr.values())
    iSort = np.argsort(list(Induction.Texc_hr.values()))

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

        fig.suptitle(FigLbl.inductionTitle)
        # Only label the bottom-left sides of axes
        [ax.set_xlabel(FigLbl.sigMeanLabel) for ax in (axes[1,0], axes[1,1])]
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
            ax.set_title(name)
            zContours = [ax.contour(Induction.sigmaMean_Sm, Induction.D_km, z[i, ...],
                           colors=Color.Induction[T], linestyles=Style.LS_Induction[T],
                           linewidths=Style.LW_Induction[T], levels=IndParams.GetClevels(fLabel, T))
                           for i, T in enumerate(Induction.Texc_hr.keys())]
            if Params.Induct.inductOtype == 'sigma':
                [ax.clabel(zContours[i], fmt=IndParams.GetCfmt(fLabel, T),
                           fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                           for i, T in enumerate(Induction.Texc_hr.keys())]

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
        fig.savefig(fNameSigma, format=FigMisc.figFormat, dpi=FigMisc.dpi)
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
                                        levels=IndParams.GetClevels(fLabel, T))
                             for i, T in enumerate(Induction.Texc_hr.keys())]
                [ax.clabel(zContours[i], fmt=IndParams.GetCfmt(fLabel, T),
                           fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
                           for i, T in enumerate(Induction.Texc_hr.keys())]

            if FigMisc.PLOT_V2021 and Induction.bodyname in ['Europa', 'Ganymede', 'Callisto']:
                AddV2021points(Params.Induct, Induction.bodyname, Params.Induct.inductOtype, allAxes)

            if Params.LEGEND:
                lines = np.array([contour.legend_elements()[0][0] for contour in zContours])
                axes[1,1].legend(lines[iSort], FigLbl.legendTexc[iSort], framealpha=FigMisc.cLegendOpacity)

            plt.tight_layout()
            fig.savefig(Params.FigureFiles.induct['Bcomps'], format=FigMisc.figFormat, dpi=FigMisc.dpi)
            log.debug(f'Plot saved to file: {Params.FigureFiles.induct["Bcomps"]}')
            plt.close()

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

            if FigMisc.PLOT_V2021 and Induction.bodyname in ['Europa', 'Ganymede', 'Callisto']:
                AddV2021points(Params.Induct, Induction.bodyname, Params.Induct.inductOtype, axes[0,:])
                AddV2021points(Params.Induct, Induction.bodyname, 'sigma', axes[1,:])

            if Params.LEGEND:
                lines = np.array([contour.legend_elements()[0][0] for contour in zContours])
                axes[1,1].legend(lines[iSort], FigLbl.legendTexc[iSort], framealpha=FigMisc.cLegendOpacity)

            plt.tight_layout()
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
        fig = plt.figure(figsize=FigSize.induct)
        grid = GridSpec(1, 2)
        axes = [fig.add_subplot(grid[0, j]) for j in range(2)]
        if Style.GRIDS:
            [ax.grid() for ax in axes]
            [ax.set_axisbelow(True) for ax in axes]

        # Labels and titles
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

        if FigMisc.PLOT_V2021 and Induction.bodyname in ['Europa', 'Ganymede', 'Callisto']:
            AddV2021points(Params.Induct, Induction.bodyname, 'sigma', axes)

        if Params.LEGEND:
            lines = np.array([contour.legend_elements()[0][0] for contour in phaseContours])
            axes[1].legend(lines[iSort], FigLbl.legendTexc[iSort], framealpha=FigMisc.cLegendOpacity)

        plt.tight_layout()
        fig.savefig(fNameSigma, format=FigMisc.figFormat, dpi=FigMisc.dpi)
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
            fig.suptitle(FigLbl.inductionTitle)
            axes[0].set_title(name)
            axes[1].set_title(FigLbl.phaseTitle)
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

            if FigMisc.PLOT_V2021 and Induction.bodyname in ['Europa', 'Ganymede', 'Callisto']:
                AddV2021points(Params.Induct, Induction.bodyname, Params.Induct.inductOtype, axes)

            if Params.LEGEND:
                lines = np.array([contour.legend_elements()[0][0] for contour in phaseContours])
                axes[1].legend(lines[iSort], FigLbl.legendTexc[iSort], framealpha=FigMisc.cLegendOpacity)

            plt.tight_layout()
            fig.savefig(Params.FigureFiles.induct[fLabel], format=FigMisc.figFormat, dpi=FigMisc.dpi)
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


def PlotComplexBdip(PlanetList, Params):

    if PlanetList[0].bodyname in Excitations.Texc_hr.keys():
        axComps = ['x', 'y', 'z']
        refs = {}
        nPeaksToPlot = np.sum([Tdo and Thave for Tdo, Thave, Texc in zip(Params.Induct.excSelectionPlot.values(),
                                                                         Params.Induct.excSelectionCalc.values(),
                                                                         Excitations.Texc_hr[PlanetList[0].bodyname].values()) if Texc is not None])
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
                axc = {vComp: row for vComp, row in zip(axComps, axes)}
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
                if Params.ALL_ONE_BODY:
                    fig.suptitle(f'{PlanetList[0].name}{titleToUse[0]}')
                else:
                    fig.suptitle(titleToUse[1])

                insetx, insety = ({vComp: 0 for vComp in axComps} for _ in range(2))
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

                        for iRow, axComp in enumerate(axComps):
                            xPlotted, yPlotted, absBiPlotted = (np.empty(0) for _ in range(3))
                            for iPeak, Tkey in enumerate(Planet.Magnetic.calcedExc):
                                if Params.Induct.excSelectionPlot[Tkey]:
                                    pts = [ax.scatter(np.real(Planet.Magnetic.Bi1xyz_nT[axComp][iPeak]), np.imag(Planet.Magnetic.Bi1xyz_nT[axComp][iPeak]),
                                                      marker=Style.MS_dip[Tkey], s=Style.MW_dip[Tkey]**2,
                                                      facecolor=thisFaceColor, edgecolor=thisEdgeColor) for ax in axc[axComp]][0].get_offsets()[0,:]

                                    if iPlanet == 0 and iRow == 0:
                                        refs[Tkey] = axes[0, -1].scatter(0, 0, label=Tkey.capitalize(), marker=Style.MS_dip[Tkey],
                                                                     s=Style.MW_dip[Tkey]**2, facecolor=Color.ref, edgecolor=Color.ref)

                                    xPlotted = np.append(xPlotted, pts[0])
                                    yPlotted = np.append(yPlotted, pts[1])
                                    absBiPlotted = np.append(absBiPlotted, np.abs(Planet.Magnetic.Bi1xyz_nT[axComp][iPeak]))

                            if DO_ZOOM:
                                iAllButLargest = np.argsort(absBiPlotted)[:-1]
                                secondLargestRe = np.max(xPlotted[iAllButLargest])
                                secondLargestIm = np.max(yPlotted[iAllButLargest])
                                insetx[axComp] = np.maximum(insetx[axComp], secondLargestRe * FigMisc.BdipZoomMult)
                                insety[axComp] = np.maximum(insety[axComp], secondLargestIm * FigMisc.BdipZoomMult)
                                axes[iRow, 0].set_xlim([0, insetx[axComp]])
                                axes[iRow, 0].set_ylim([0, insety[axComp]])

                if DO_ZOOM and FigMisc.SHOW_INSET:
                    refs[f'inset'] = axes[0, -1].plot([0,1], [1,0], color=Color.BdipInset, linewidth=Style.LW_BdipInset,
                                              linestyle=Style.LS_BdipInset, label='Inset region')[0]
                    for iRow, axComp in enumerate(axComps):
                        axes[iRow, -1].add_patch(Rectangle((0,0), insetx[axComp], insety[axComp], edgecolor=Color.BdipInset, zorder=-1,
                                           linewidth=Style.LW_BdipInset, linestyle=Style.LS_BdipInset, facecolor='None'))

                if Params.LEGEND:
                    axes[0, -1].legend()

                for refPt in refs.values():
                    refPt.remove()

                [axes[iRow, -1].set_xlim(left=0) for iRow in range(3)]
                [axes[iRow, -1].set_ylim(bottom=0) for iRow in range(3)]
                plt.tight_layout()
                fig.savefig(Params.FigureFiles.Bdip['all'], format=FigMisc.figFormat, dpi=FigMisc.dpi)
                log.debug(f'Induced dipole surface strength plot saved to file: {Params.FigureFiles.Bdip["all"]}')
                plt.close()

            else:
                for axComp in axComps:
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
                                    pts = [ax.scatter(np.real(Planet.Magnetic.Bi1xyz_nT[axComp][iPeak]), np.imag(Planet.Magnetic.Bi1xyz_nT[axComp][iPeak]),
                                                      marker=Style.MS_dip[Tkey], s=Style.MW_dip[Tkey]**2,
                                                      facecolor=thisFaceColor, edgecolor=thisEdgeColor) for ax in axes][0].get_offsets()[0,:]

                                    if iPlanet == 0:
                                        refs[Tkey] = axes[-1].scatter(0, 0, label=Tkey.capitalize(), marker=Style.MS_dip[Tkey],
                                                                      s=Style.MW_dip[Tkey]**2, facecolor=Color.ref, edgecolor=Color.ref)

                                    xPlotted = np.append(xPlotted, pts[0])
                                    yPlotted = np.append(yPlotted, pts[1])
                                    absBiPlotted = np.append(absBiPlotted, np.abs(Planet.Magnetic.Bi1xyz_nT[axComp][iPeak]))

                            if DO_ZOOM:
                                iAllButLargest = np.argsort(absBiPlotted)[:-1]
                                secondLargestRe = np.max(xPlotted[iAllButLargest])
                                secondLargestIm = np.max(yPlotted[iAllButLargest])
                                insetx = np.maximum(insetx, secondLargestRe * FigMisc.BdipZoomMult)
                                insety = np.maximum(insety, secondLargestIm * FigMisc.BdipZoomMult)
                                axes[0].set_xlim([0, insetx])
                                axes[0].set_ylim([0, insety])

                    if DO_ZOOM and FigMisc.SHOW_INSET:
                        axes[-1].add_patch(Rectangle((0,0), insetx, insety, edgecolor=Color.BdipInset, zorder=-1,
                                           linewidth=Style.LW_BdipInset, linestyle=Style.LS_BdipInset, facecolor='None'))
                        refs['inset'] = axes[-1].plot([0,1], [1,0], color=Color.BdipInset, linewidth=Style.LW_BdipInset,
                                                      linestyle=Style.LS_BdipInset, label='Inset region')[0]

                    if Params.LEGEND:
                        axes[-1].legend()

                    for refPt in refs.values():
                        refPt.remove()

                    axes[-1].set_xlim(left=0)
                    axes[-1].set_ylim(bottom=0)
                    plt.tight_layout()
                    fig.savefig(Params.FigureFiles.Bdip[axComp], format=FigMisc.figFormat, dpi=FigMisc.dpi)
                    log.debug(f'Induced dipole surface strength plot saved to file: {Params.FigureFiles.Bdip[axComp]}')
                    plt.close()

    return


def PlotMagSurface(PlanetList, Params):
    """ Plot magnetic field on a spherical surface of fixed radius using 
        calculations from MoonMag. 
    """
    PlanetListSubset = [Planet for Planet in PlanetList if Planet.bodyname in Excitations.Texc_hr.keys()]

    if not isinstance(FigMisc.tMagEval_s, Iterable):
        FigMisc.tMagEval_s = np.array([FigMisc.tMagEval_s])    
    for Planet in PlanetListSubset:
        if Planet.lonMap_deg is None:
            Planet.lonMap_deg = FigMisc.lonMap_deg
            Planet.nLonMap = np.size(Planet.lonMap_deg)
        if Planet.phiMap_rad is None:
            Planet.phiMap_rad = np.radians(Planet.lonMap_deg)
        if Planet.latMap_deg is None:
            Planet.latMap_deg = FigMisc.latMap_deg
            Planet.nLatMap = np.size(Planet.latMap_deg)
        if Planet.thetaMap_rad is None:
            Planet.thetaMap_rad = np.radians(90 - Planet.latMap_deg)
        rMagEvalLbl, rMagEvalPrint = FigLbl.rStr(FigMisc.rMagEval_Rp, Planet.bodyname)
        if Params.Sig.INCLUDE_ASYM:
            asymStr = ',asym'
            # Get separate spherically symmetric induced moments
            BinmSph_nT = np.zeros_like(Planet.Magnetic.Benm_nT)
            for n in range(1, Planet.Magnetic.nprmMax+1):
                for iExc in range(Planet.Magnetic.nExc):
                    BinmSph_nT[iExc,:,n,:] = n / (n+1) * Planet.Magnetic.Aen[iExc,n] * Planet.Magnetic.Benm_nT[iExc,:,n,:]
                    
        else:
            asymStr = ''
            BinmSph_nT = None
    
        for tEval_s in FigMisc.tMagEval_s:
            BvecRe_nT = {vComp: np.zeros((Planet.nLatMap, Planet.nLonMap)) for vComp in ['x', 'y', 'z', 'mag']}
            tMagEvalLbl, tMagEvalPrint, tFnameEnd = FigLbl.tStr(tEval_s)
            phaseNow = np.exp(-1j * Planet.Magnetic.omegaExc_radps * tEval_s)

            log.debug(f'Evaluating induced magnetic field at {rMagEvalPrint}, {tMagEvalPrint}.')
            for iExc in range(Planet.Magnetic.nExc):
                BinmNow_nT = Planet.Magnetic.Binm_nT[iExc,...] * phaseNow[iExc]
                Bx_nT, By_nT, Bz_nT = GetMagSurf(Planet.Magnetic.nLin, Planet.Magnetic.mLin, BinmNow_nT, 
                                                 FigMisc.rMagEval_Rp, Planet.thetaMap_rad, Planet.phiMap_rad,
                                                 do_parallel=Params.DO_PARALLEL)
                BvecRe_nT['x'] = BvecRe_nT['x'] + np.real(Bx_nT)
                BvecRe_nT['y'] = BvecRe_nT['y'] + np.real(By_nT)
                BvecRe_nT['z'] = BvecRe_nT['z'] + np.real(Bz_nT)
                
            if FigMisc.vCompMagSurf == 'mag':
                BvecRe_nT['mag'] = np.sqrt(BvecRe_nT['x']**2 + BvecRe_nT['y']**2 + BvecRe_nT['z']**2)
                
            if FigMisc.BASYM_WITH_SYM:
                fig = plt.figure(figsize=FigSize.MagSurfCombo)
                grid = GridSpec(1, 2)
            else:
                fig = plt.figure(figsize=FigSize.MagSurf)
                grid = GridSpec(1, 1)
            ax = fig.add_subplot(grid[0, -1])
            SetMap(ax)
            if FigMisc.LARGE_ADJUST:
                if FigMisc.vCompMagSurf == 'mag':
                    compLabel = f'$B_\mathrm{{{FigMisc.vCompMagSurf}{asymStr}}}$'
                else:
                    compLabel = f'$B_{{{FigMisc.vCompMagSurf}\mathrm{{{asymStr}}}}}$'
                title = f'{FigLbl.MagSurfShortTitle} {compLabel} ' + \
                        f'at {rMagEvalLbl} ($\si{{nT}}$)'
            else:
                if FigMisc.vCompMagSurf == 'mag':
                    compLabel = 'magnitude'
                else:
                    compLabel = f'${FigMisc.vCompMagSurf}$ component'
                title = f'{Planet.name} {FigLbl.MagSurfTitle} {compLabel}, ' + \
                        f'{rMagEvalLbl} ($\si{{nT}}$), {tMagEvalLbl}'
            ax.set_title(title, size=FigMisc.mapTitleSize)

            vmin = FigMisc.vminMagSurf_nT
            vmax = FigMisc.vmaxMagSurf_nT
            cLevels = None
            if np.all(BvecRe_nT[FigMisc.vCompMagSurf] >= 0):
                cmap = Color.cmap['BmapPos']
            elif np.all(BvecRe_nT[FigMisc.vCompMagSurf] < 0):
                cmap = Color.cmap['BmapNeg']
            else:
                cmap = Color.cmap['BmapDiv']
                if vmin is None and vmax is None:
                    vAbsMax = np.max(np.abs(BvecRe_nT[FigMisc.vCompMagSurf]))
                    vmin = -vAbsMax
                    vmax = vAbsMax
                    cLevels = np.linspace(vmin, vmax, FigMisc.nMagContours)

            Bmap = ax.pcolormesh(Planet.lonMap_deg, Planet.latMap_deg, BvecRe_nT[FigMisc.vCompMagSurf],
                          shading='auto', cmap=cmap, vmin=vmin, vmax=vmax)
            BmapContours = ax.contour(Planet.lonMap_deg, Planet.latMap_deg, BvecRe_nT[FigMisc.vCompMagSurf],
                              levels=cLevels, colors='black')
            ax.clabel(BmapContours, fmt=ticker.FuncFormatter(FigMisc.Cformat), 
                      fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)

            if not FigMisc.LARGE_ADJUST:
                cbar = fig.colorbar(Bmap, ax=ax)
                cbar.ax.set_title(FigLbl.MagSurfCbarTitle)

            ax.set_aspect(1)
            if not FigMisc.BASYM_WITH_SYM:
                plt.tight_layout()
                fName = f'{Params.FigureFiles.MagSurf[FigMisc.vCompMagSurf]}{tFnameEnd}{FigMisc.xtn}'
                fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat, dpi=FigMisc.dpi)
                log.debug(f'Induced field surface map saved to file: {fName}')
                plt.close()

            if Params.Sig.INCLUDE_ASYM:
                log.debug(f'Evaluating symmetric induced magnetic field at {rMagEvalPrint}, {tMagEvalPrint}.')
                BvecReSym_nT = {vComp: np.zeros((Planet.nLatMap, Planet.nLonMap)) for vComp in ['x', 'y', 'z', 'mag']}
                for iExc in range(Planet.Magnetic.nExc):
                    BinmSphNow_nT = BinmSph_nT[iExc,...] * phaseNow[iExc]
                    Bx_nT, By_nT, Bz_nT = GetMagSurf(Planet.Magnetic.nprmLin, Planet.Magnetic.mprmLin, BinmSphNow_nT, 
                                                     FigMisc.rMagEval_Rp, Planet.thetaMap_rad, Planet.phiMap_rad,
                                                     do_parallel=Params.DO_PARALLEL)
                    BvecReSym_nT['x'] = BvecReSym_nT['x'] + np.real(Bx_nT)
                    BvecReSym_nT['y'] = BvecReSym_nT['y'] + np.real(By_nT)
                    BvecReSym_nT['z'] = BvecReSym_nT['z'] + np.real(Bz_nT)
                
                if FigMisc.vCompMagSurf == 'mag':
                    BvecReSym_nT['mag'] = np.sqrt(BvecReSym_nT['x']**2 + BvecReSym_nT['y']**2 + BvecReSym_nT['z']**2)

                if not FigMisc.BASYM_WITH_SYM:
                    fig = plt.figure(figsize=FigSize.MagSurf)
                    grid = GridSpec(1, 1)

                ax = fig.add_subplot(grid[0, 0])
                SetMap(ax)
                if FigMisc.LARGE_ADJUST:
                    if FigMisc.vCompMagSurf == 'mag':
                        compLabel = f'$B_\mathrm{{{FigMisc.vCompMagSurf},sym}}$'
                    else:
                        compLabel = f'$B_{{{FigMisc.vCompMagSurf}\mathrm{{,sym}}}}$'
                    title = f'{FigLbl.MagSurfShortTitle} {compLabel} ' + \
                            f'at {rMagEvalLbl} ($\si{{nT}}$)'
                else:
                    # compLabel able to be reused from asymmetric in this case
                    title = f'{Planet.name} {FigLbl.MagSurfSymTitle} {compLabel}, ' + \
                            f'{rMagEvalLbl} ($\si{{nT}}$), {tMagEvalLbl}'
                ax.set_title(title, size=FigMisc.mapTitleSize)

                vmin = FigMisc.vminMagSurf_nT
                vmax = FigMisc.vmaxMagSurf_nT
                cLevels = None
                if np.all(BvecReSym_nT[FigMisc.vCompMagSurf] > 0):
                    cmap = Color.cmap['BmapPos']
                elif np.all(BvecReSym_nT[FigMisc.vCompMagSurf] < 0):
                    cmap = Color.cmap['BmapNeg']
                else:
                    cmap = Color.cmap['BmapDiv']
                    if vmin is None and vmax is None:
                        vAbsMax = np.max(np.abs(BvecReSym_nT[FigMisc.vCompMagSurf]))
                        vmin = -vAbsMax
                        vmax = vAbsMax
                        cLevels = np.linspace(vmin, vmax, FigMisc.nMagContours)

                Bmap = ax.pcolormesh(Planet.lonMap_deg, Planet.latMap_deg, BvecReSym_nT[FigMisc.vCompMagSurf],
                                     shading='auto', cmap=cmap, vmin=vmin, vmax=vmax)
                BmapContours = ax.contour(Planet.lonMap_deg, Planet.latMap_deg, BvecReSym_nT[FigMisc.vCompMagSurf],
                                          levels=cLevels, colors='black')
                ax.clabel(BmapContours, fmt=ticker.FuncFormatter(FigMisc.Cformat),
                          fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)

                if not FigMisc.LARGE_ADJUST:
                    cbar = fig.colorbar(Bmap, ax=ax)
                    cbar.ax.set_title(FigLbl.MagSurfCbarTitle)

                ax.set_aspect(1)
                plt.tight_layout()
                if FigMisc.BASYM_WITH_SYM:
                    fName = f'{Params.FigureFiles.MagSurfCombo[FigMisc.vCompMagSurf]}{tFnameEnd}{FigMisc.xtn}'
                    fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat, dpi=FigMisc.dpi)
                    log.debug(f'Symmetric/asymmetric induced field surface maps saved to file: {fName}')
                else:
                    fName = f'{Params.FigureFiles.MagSurfSym[FigMisc.vCompMagSurf]}{tFnameEnd}{FigMisc.xtn}'
                    fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat, dpi=FigMisc.dpi)
                    log.debug(f'Symmetric induced field surface map saved to file: {fName}')
                
                plt.close()
                
                # Now plot the difference between asymmetric and symmetric for this component
                Bdiff_nT = BvecRe_nT[FigMisc.vCompMagSurf] - BvecReSym_nT[FigMisc.vCompMagSurf]

                fig = plt.figure(figsize=FigSize.MagSurf)
                grid = GridSpec(1, 1)
                ax = fig.add_subplot(grid[0, 0])
                SetMap(ax)
                if FigMisc.LARGE_ADJUST:
                    if FigMisc.vCompMagSurf == 'mag':
                        compLabel = f'$\Delta B_\mathrm{{{FigMisc.vCompMagSurf},sym}}$'
                    else:
                        compLabel = f'$\Delta B_{{{FigMisc.vCompMagSurf}\mathrm{{,sym}}}}$'
                    title = f'{FigLbl.MagSurfShortTitle} {compLabel} ' + \
                            f'at {rMagEvalLbl} ($\si{{nT}}$)'
                else:
                    # compLabel able to be reused from above in this case
                    title = f'{Planet.name} {FigLbl.MagSurfTitle} {compLabel} {FigLbl.MagSurfDiffTitle}, ' + \
                            f'{rMagEvalLbl} ($\si{{nT}}$), {tMagEvalLbl}'
                ax.set_title(title, size=FigMisc.mapTitleSize)

                vmin = FigMisc.vminMagSurfDiff_nT
                vmax = FigMisc.vmaxMagSurfDiff_nT
                cLevels = None
                cmap = Color.cmap['BmapDiv']
                if vmin is None and vmax is None:
                    vAbsMax = np.max(np.abs(Bdiff_nT))
                    vmin = -vAbsMax
                    vmax = vAbsMax
                    cLevels = np.linspace(vmin, vmax, FigMisc.nMagContours)

                Bmap = ax.pcolormesh(Planet.lonMap_deg, Planet.latMap_deg, Bdiff_nT,
                                     shading='auto', cmap=cmap, vmin=vmin, vmax=vmax)
                BmapContours = ax.contour(Planet.lonMap_deg, Planet.latMap_deg, Bdiff_nT,
                                          levels=cLevels, colors='black')
                ax.clabel(BmapContours, fmt=ticker.FuncFormatter(FigMisc.Cformat),
                          fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)

                if not FigMisc.LARGE_ADJUST:
                    fig.colorbar(Bmap, ax=ax, label=FigLbl.MagSurfCbarDiffLabel)

                ax.set_aspect(1)
                plt.tight_layout()
                fName = f'{Params.FigureFiles.MagSurfDiff[FigMisc.vCompMagSurf]}{tFnameEnd}{FigMisc.xtn}'
                fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat, dpi=FigMisc.dpi)
                log.debug(f'Induced field difference surface map saved to file: {fName}')
                plt.close()

    return


def SetMap(ax):
    ax.set_xticks(FigMisc.lonMapTicks_deg)
    ax.set_yticks(FigMisc.latMapTicks_deg)
    ax.tick_params(axis='both', which='major', labelsize=FigMisc.latlonSize)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(FigMisc.LonMapFormatter))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(FigMisc.LatMapFormatter))

    return


def PlotAsym(PlanetList, Params):
    """ Plot a contour map showing the asymmetry present in selected boundaries.
        This function is time-consuming, so it is only applied to the first body
        in PlanetList (the primary body)
    """
    
    if Params.Sig.INCLUDE_ASYM:
        if np.size(PlanetList[0].Magnetic.iAsymBds) > 0:
            mainPlanet = PlanetList[0]
        else:
            mainPlanet = next(Planet for Planet in PlanetList if np.size(Planet.Magnetic.iAsymBds) > 0)

        if mainPlanet.lonMap_deg is None:
            mainPlanet.lonMap_deg = FigMisc.lonMap_deg
            mainPlanet.nLonMap = np.size(mainPlanet.lonMap_deg)
        if mainPlanet.latMap_deg is None:
            mainPlanet.latMap_deg = FigMisc.latMap_deg
            mainPlanet.nLatMap = np.size(mainPlanet.latMap_deg)

        for i, zMean_km in enumerate(mainPlanet.Magnetic.zMeanAsym_km):
            fig = plt.figure(figsize=FigSize.asym)
            grid = GridSpec(1, 1)
            ax = fig.add_subplot(grid[0, 0])
            SetMap(ax)
            if i == mainPlanet.Magnetic.nAsymBds-1:
                title = FigLbl.asymGravTitle
                cLevelsAsym = None
            else:
                if mainPlanet.Magnetic.asymDescrip is not None and i < np.size(mainPlanet.Magnetic.asymDescrip):
                    descrip = f'{mainPlanet.Magnetic.asymDescrip[i]}{FigLbl.asymAfterDescrip}'
                    if mainPlanet.Magnetic.asymDescrip[i] in mainPlanet.Magnetic.asymContours_km.keys():
                        cLevelsAsym = mainPlanet.Magnetic.asymContours_km[mainPlanet.Magnetic.asymDescrip[i]]
                    else:
                        cLevelsAsym = None
                else:
                    descrip = FigLbl.asymTitle
                    cLevelsAsym = None
                title = f'{descrip}{np.abs(zMean_km):.1f}'
            ax.set_title(title, size=FigMisc.mapTitleSize)

            if cLevelsAsym is None:
                cLevelsAsym = np.unique(np.round(np.linspace(np.min(mainPlanet.Magnetic.asymDevs_km[i,...]),
                                                             np.max(mainPlanet.Magnetic.asymDevs_km[i,...]),
                                                             FigMisc.nAsymContours)))

            asymMap = ax.pcolormesh(mainPlanet.lonMap_deg, mainPlanet.latMap_deg, mainPlanet.Magnetic.asymDevs_km[i,...],
                                 shading='auto', cmap=Color.cmap['asymDev'])
            asymContours = ax.contour(mainPlanet.lonMap_deg, mainPlanet.latMap_deg, mainPlanet.Magnetic.asymDevs_km[i,...],
                                      levels=cLevelsAsym, colors='black')
            ax.clabel(asymContours, fmt=ticker.FuncFormatter(FigMisc.Cformat),
                      fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)

            if not FigMisc.LARGE_ADJUST:
                fig.colorbar(asymMap, ax=ax, label=FigLbl.asymCbarLabel)

            ax.set_aspect(1)
            plt.tight_layout()
            fName = f'{Params.FigureFiles.asym}{zMean_km:.1f}{FigMisc.xtn}'
            fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat, dpi=FigMisc.dpi)
            log.debug(f'Asymmetric boundary surface map for z = {zMean_km:.1f} km saved to file: {fName}')
            plt.close()

    return


def PlotMagSpectrum(PlanetList, Params):
    """ Plot Fourier spectra of magnetic excitations, complex response
        amplitude, and induced dipole components.
    """

    axComps = ['x', 'y', 'z']
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
    if Params.MagSpectrum.Tmin_hr is None:
        Tmin_hr = np.min(mainPlanet.Magnetic.TexcFT_hr)
    else:
        Tmin_hr = Params.MagSpectrum.Tmin_hr
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
    if mainPlanet.Magnetic.coordTypeFT is not None:
        FTtitle = f'{mainPlanet.name} {FigLbl.MagFTtitle}, {mainPlanet.Magnetic.coordTypeFT} coordinates'
    else:
        FTtitle = f'{mainPlanet.name} {FigLbl.MagFTtitle}'
    fig.suptitle(FTtitle)

    # Plot the data
    [axes[0].plot(freqVar, np.abs(mainPlanet.Magnetic.Be1xyzFT_nT[vComp]), label=f'$B^e_{vComp}$',
                  color=Color.BeiFT[vComp], linestyle=Style.LS_FT, linewidth=Style.LW_FT) for vComp in axComps]
    [axes[2].plot(freqVar, np.abs(mainPlanet.Magnetic.Bi1xyzFT_nT[vComp]), label=f'$B^i_{vComp}$',
                  color=Color.BeiFT[vComp], linestyle=Style.LS_FT, linewidth=Style.LW_FT) for vComp in axComps]
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
    fig.savefig(Params.FigureFiles.MagFT, format=FigMisc.figFormat, dpi=FigMisc.dpi)
    log.debug(f'Magnetic Fourier spectra plot saved to file: {Params.FigureFiles.MagFT}')
    plt.close()

    return
