import os
import numpy as np
import logging
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap as DiscreteCmap, to_rgb, BoundaryNorm
from scipy.interpolate import interp1d
from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc
from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS, GetIceEOS
from PlanetProfile.Utilities.Indexing import PhaseConv, PhaseInv
from PlanetProfile.Thermodynamics.InnerEOS import GetInnerEOS
from PlanetProfile.Utilities.defineStructs import Constants

# Assign logger
log = logging.getLogger('PlanetProfile')

# Unlike most other plotting routines, for those below we can't plot multiple bodies together.
# Just make these plots for the first model, which is the primary.

def PlotHydroPhase(PlanetList, Params):

    if FigMisc.PminHydro_MPa is None:
        Pmin_MPa = np.min([Planet.P_MPa[0] for Planet in PlanetList])
    else:
        Pmin_MPa = FigMisc.PminHydro_MPa
    if FigMisc.PmaxHydro_MPa is None:
        Pmax_MPa = np.max([Planet.P_MPa[Planet.Steps.nHydro-1] for Planet in PlanetList])
    else:
        Pmax_MPa = FigMisc.PmaxHydro_MPa
    if FigMisc.TminHydro_K is None:
        if not np.any([Planet.Do.NO_OCEAN for Planet in PlanetList]):
            Tmin_K = np.min([np.min(Planet.T_K[:Planet.Steps.nHydro]) for Planet in PlanetList])
        else:
            Tmin_K = np.min([np.min(Planet.T_K) for Planet in PlanetList])
    else:
        Tmin_K = FigMisc.TminHydro_K
    if FigMisc.TmaxHydro_K is None:
        if np.all([Planet.Do.NO_OCEAN for Planet in PlanetList]):
            Tmax_K = np.max([Planet.Sil.THydroMax_K for Planet in PlanetList])
        else:
            Tmax_K = np.max([np.max(Planet.T_K[:Planet.Steps.nHydro])
                             for Planet in PlanetList if not Planet.Do.NO_OCEAN])
    else:
        Tmax_K = FigMisc.TmaxHydro_K

    if os.path.dirname(Params.FigureFiles.vpvtHydro) != 'Comparison':
        Planet = PlanetList[0]
        if FigMisc.nPphase is None:
            P_MPa = Planet.P_MPa[:Planet.Steps.nHydro]
            P_MPa = P_MPa[np.logical_and(P_MPa >= Pmin_MPa, P_MPa <= Pmax_MPa)]
        else:
            P_MPa = np.linspace(Pmin_MPa, Pmax_MPa, FigMisc.nPphase)
        if FigMisc.nTphase is None:
            T_K = Planet.T_K[:Planet.Steps.nHydro]
            T_K = T_K[np.logical_and(T_K >= Tmin_K, T_K <= Tmax_K)]
        else:
            T_K = np.linspace(Tmin_K, Tmax_K, FigMisc.nTphase)
        # Load EOS independently from model run, because we will query wider ranges of conditions
        oceanEOS = GetOceanEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, P_MPa, T_K,
                               Planet.Ocean.MgSO4elecType, rhoType=Planet.Ocean.MgSO4rhoType,
                               scalingType=Planet.Ocean.MgSO4scalingType, FORCE_NEW=Params.FORCE_EOS_RECALC,
                               phaseType=Planet.Ocean.phaseType, EXTRAP=Params.EXTRAP_OCEAN,
                               sigmaFixed_Sm=Planet.Ocean.sigmaFixed_Sm, LOOKUP_HIRES=Planet.Do.OCEAN_PHASE_HIRES)

        phases = oceanEOS.fn_phase(P_MPa, T_K, grid=True).astype(int)
        # Add clathrates to phase and property diagrams where it is stable (if modeled)
        if Planet.Do.CLATHRATE:
            clath = PhaseConv(Constants.phaseClath)
            clathEOS = GetIceEOS(P_MPa, T_K, clath,
                                        porosType=Planet.Ocean.porosType[clath],
                                        phiTop_frac=Planet.Ocean.phiMax_frac[clath],
                                        Pclosure_MPa=Planet.Ocean.Pclosure_MPa[clath],
                                        phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE[clath])
            clathStable = clathEOS.fn_phase(P_MPa, T_K, grid=True).astype(int)
            phases[np.where(np.logical_and(clathStable == Constants.phaseClath, 
                                           phases == 1))] = Constants.phaseClath
        
        fig = plt.figure(figsize=FigSize.vphase)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(False)
        # Labels and titles
        ax.set_xlabel(FigLbl.Tlabel)
        ax.set_ylabel(FigLbl.PlabelHydro)

        # Set overall figure title
        if Params.TITLES:
            fig.suptitle(f'{Planet.compStr}{FigLbl.hydroPhaseTitle}')

        # Plot as colormesh
        phaseColors = DiscreteCmap.from_list('icePhases',
                                             [Color.oceanCmap(1),
                                              Color.iceIcond,
                                              Color.iceII,
                                              Color.iceIII,
                                              Color.iceV,
                                              Color.iceVI,
                                              Color.clathCond], N=7)
        phaseBounds = np.array([0, 0.5, 1.5, 2.5, 4, 5.5, (Constants.phaseClath+6)/2, Constants.phaseClath])
        bNorm = BoundaryNorm(boundaries=phaseBounds, ncolors=7)
        ax.pcolormesh(T_K, P_MPa * FigLbl.PmultHydro, phases, norm=bNorm, cmap=phaseColors, rasterized=FigMisc.PT_RASTER)

        ices = np.unique(phases[phases != 0])
        P2D_MPa, T2D_K = np.meshgrid(P_MPa, T_K, indexing='ij')
        for ice in ices:
            theseP = P2D_MPa[np.where(phases == ice)]
            theseT = T2D_K[np.where(phases == ice)]
            P = (np.max(theseP) + np.min(theseP)) / 2
            T = (np.max(theseT) + np.min(theseT)) / 2
            ax.text(T, P, PhaseConv(ice), ha='center', va='center', fontsize=FigLbl.hydroPhaseSize)

        if np.any(phases == 0):
            Pliq = P2D_MPa[np.where(phases == 0)]
            Tliq = T2D_K[np.where(phases == 0)]
            P = (np.max(Pliq) + np.min(Pliq)) / 2
            T = (np.max(Tliq) + np.min(Tliq)) / 2
            ax.text(T, P, 'liquid', ha='center', va='center', fontsize=FigLbl.hydroPhaseSize)

        # Plot geotherm(s) on top of colormaps
        for eachPlanet in PlanetList:
            # Geotherm curve
            if np.size(PlanetList) > 1:
                thisColor = None
            else:
                thisColor = Color.geothermHydro
            if eachPlanet.Do.NO_DIFFERENTIATION or eachPlanet.Do.PARTIAL_DIFFERENTIATION:
                Pgeo = eachPlanet.P_MPa * FigLbl.PmultHydro
                Tgeo = eachPlanet.T_K
            else:
                Pgeo = eachPlanet.P_MPa[:eachPlanet.Steps.nHydro] * FigLbl.PmultHydro
                Tgeo = eachPlanet.T_K[:eachPlanet.Steps.nHydro]
            ax.plot(Tgeo, Pgeo, linewidth=Style.LW_geotherm, linestyle=Style.LS_geotherm,
                     color=thisColor, label=eachPlanet.label)

            if Params.LEGEND and np.size(PlanetList) > 1:
                handles, lbls = ax.get_legend_handles_labels()
                ax.legend(handles, lbls)

        ax.set_xlim([Tmin_K, Tmax_K])
        ax.set_ylim([Pmin_MPa, Pmax_MPa])
        ax.invert_yaxis()
        plt.tight_layout()
        fig.savefig(Params.FigureFiles.vphase, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Hydrosphere phase diagram saved to file: {Params.FigureFiles.vphase}')
        plt.close()
    
    return


def PlotPvThydro(PlanetList, Params):

    if os.path.dirname(Params.FigureFiles.vpvtHydro) != 'Comparison':
        Planet = PlanetList[0]
        if FigMisc.PminHydro_MPa is None:
            Pmin_MPa = np.min([Planet.P_MPa[0] for Planet in PlanetList])
        else:
            Pmin_MPa = FigMisc.PminHydro_MPa
        if FigMisc.PmaxHydro_MPa is None:
            Pmax_MPa = np.max([Planet.P_MPa[Planet.Steps.nHydro-1] for Planet in PlanetList])
        else:
            Pmax_MPa = FigMisc.PmaxHydro_MPa
        if FigMisc.TminHydro_K is None:
            if not np.any([Planet.Do.NO_OCEAN for Planet in PlanetList]):
                Tmin_K = np.min([np.min(Planet.T_K[:Planet.Steps.nHydro]) for Planet in PlanetList])
            else:
                Tmin_K = np.min([np.min(Planet.T_K) for Planet in PlanetList])
        else:
            Tmin_K = FigMisc.TminHydro_K
        if FigMisc.TmaxHydro_K is None:
            if np.all([Planet.Do.NO_OCEAN for Planet in PlanetList]):
                Tmax_K = np.max([Planet.Sil.THydroMax_K for Planet in PlanetList])
            else:
                Tmax_K = np.max([np.max(Planet.T_K[:Planet.Steps.nHydro])
                                 for Planet in PlanetList if not Planet.Do.NO_OCEAN])
        else:
            Tmax_K = FigMisc.TmaxHydro_K

        if FigMisc.nPhydro is None:
            P_MPa = Planet.P_MPa[:Planet.Steps.nHydro]
            P_MPa = P_MPa[np.logical_and(P_MPa >= Pmin_MPa, P_MPa <= Pmax_MPa)]
        else:
            P_MPa = np.linspace(Pmin_MPa, Pmax_MPa, FigMisc.nPhydro)
        if FigMisc.nThydro is None:
            T_K = Planet.T_K[:Planet.Steps.nHydro]
            T_K = T_K[np.logical_and(T_K >= Tmin_K, T_K <= Tmax_K)]
        else:
            T_K = np.linspace(Tmin_K, Tmax_K, FigMisc.nThydro)

        # Load EOS independently from model run, because we will query wider ranges of conditions
        oceanEOS = GetOceanEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, P_MPa, T_K,
                               Planet.Ocean.MgSO4elecType, rhoType=Planet.Ocean.MgSO4rhoType,
                               scalingType=Planet.Ocean.MgSO4scalingType, FORCE_NEW=Params.FORCE_EOS_RECALC,
                               phaseType=Planet.Ocean.phaseType, EXTRAP=Params.EXTRAP_OCEAN,
                               sigmaFixed_Sm=Planet.Ocean.sigmaFixed_Sm, LOOKUP_HIRES=Planet.Do.OCEAN_PHASE_HIRES)

        phases = oceanEOS.fn_phase(P_MPa, T_K, grid=True).astype(int)
        ices = [PhaseConv(ice) for ice in np.unique(np.append(phases[phases != 0], 1))]
        iceEOS = {PhaseInv(ice): GetIceEOS(P_MPa, T_K, ice,
                                 porosType=Planet.Ocean.porosType[ice],
                                 phiTop_frac=Planet.Ocean.phiMax_frac[ice],
                                 Pclosure_MPa=Planet.Ocean.Pclosure_MPa[ice],
                                 phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE[ice],
                                 ICEIh_DIFFERENT=Planet.Do.ICEIh_DIFFERENT)
                  for ice in ices}
        # Add clathrates to phase and property diagrams where it is stable (if modeled)
        if Planet.Do.CLATHRATE:
            clath = PhaseConv(Constants.phaseClath)
            ices.append(clath)
            iceEOS[Constants.phaseClath] = GetIceEOS(P_MPa, T_K, clath,
                                        porosType=Planet.Ocean.porosType[clath],
                                        phiTop_frac=Planet.Ocean.phiMax_frac[clath],
                                        Pclosure_MPa=Planet.Ocean.Pclosure_MPa[clath],
                                        phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE[clath])
            clathStable = iceEOS[Constants.phaseClath].fn_phase(P_MPa, T_K, grid=True)
            phases[clathStable == Constants.phaseClath] = Constants.phaseClath

        fig = plt.figure(figsize=FigSize.vpvt)
        grid = GridSpec(2, 4)
        axes = np.array([[fig.add_subplot(grid[i, j]) for j in range(4)] for i in range(2)])
        axf = axes.flatten()
        if Style.GRIDS:
            [ax.grid() for ax in axf]
            [ax.set_axisbelow(False) for ax in axf]
        # Labels and titles
        [ax.set_xlabel(FigLbl.Tlabel) for ax in axes[1, :]]
        [ax.set_ylabel(FigLbl.PlabelHydro) for ax in axes[:, 0]]
        [ax.set_xlim([Tmin_K, Tmax_K]) for ax in axf]
        [ax.set_ylim([Pmin_MPa, Pmax_MPa]) for ax in axf]
        [ax.invert_yaxis() for ax in axf]
        axes[0,0].set_title(FigLbl.rhoLabel)
        axes[1,0].set_title(FigLbl.CpLabel)
        axes[1,1].set_title(FigLbl.alphaLabel)
        axes[0,1].set_title(FigLbl.sigLabel)
        axes[0,2].set_title(FigLbl.VPlabel)
        axes[1,2].set_title(FigLbl.VSlabel)
        axes[0,3].set_title(FigLbl.KSlabel)
        axes[1,3].set_title(FigLbl.GSlabel)

        # Set overall figure title
        if Params.TITLES:
            fig.suptitle(f'{Planet.name}{FigLbl.PvTtitleHydro}')

        # Get data to plot -- ocean EOS properties first
        rho = oceanEOS.fn_rho_kgm3(P_MPa, T_K, grid=True)
        Cp = oceanEOS.fn_Cp_JkgK(P_MPa, T_K, grid=True)
        alpha = oceanEOS.fn_alpha_pK(P_MPa, T_K, grid=True)
        VP, KS = oceanEOS.fn_Seismic(P_MPa, T_K, grid=True)
        sig = oceanEOS.fn_sigma_Sm(P_MPa, T_K, grid=True)
        VS, GS = (np.empty_like(rho)*np.nan for _ in range(2))

        # Exclude obviously erroneous Cp values, which happen when extending beyond the knots
        # for the input EOS. This is mainly a problem with GSW and Seawater compositions.
        Cp[np.logical_or(Cp < 3200, Cp > 5200)] = np.nan

        # Now get data for all ice EOSs and replace in grid
        for iceStr in ices:
            ice = PhaseInv(iceStr)
            rho[np.where(phases == ice)] = iceEOS[ice].fn_rho_kgm3(P_MPa, T_K, grid=True)[np.where(phases == ice)]
            Cp[np.where(phases == ice)] = iceEOS[ice].fn_Cp_JkgK(P_MPa, T_K, grid=True)[np.where(phases == ice)]
            alpha[np.where(phases == ice)] = iceEOS[ice].fn_alpha_pK(P_MPa, T_K, grid=True)[np.where(phases == ice)]
            VPice, VSice, KSice, GSice = iceEOS[ice].fn_Seismic(P_MPa, T_K, grid=True)
            VP[np.where(phases == ice)] = VPice[np.where(phases == ice)]
            VS[np.where(phases == ice)] = VSice[np.where(phases == ice)]
            KS[np.where(phases == ice)] = KSice[np.where(phases == ice)]
            GS[np.where(phases == ice)] = GSice[np.where(phases == ice)]
            sig[np.where(phases == ice)] = Planet.Ocean.sigmaIce_Sm[iceStr]

        # Highlight places where alpha is negative with opposite side of diverging colormap, 0 pegged to middle
        minAlpha = np.minimum(0, np.min(alpha))
        alphaCmap = Color.ComboPvThydroCmap(minAlpha, np.max(alpha))

        # Plot colormaps of hydrosphere data
        Pscaled = P_MPa * FigLbl.PmultHydro
        rhoPlot =   axes[0,0].pcolormesh(T_K, Pscaled, rho, cmap=Color.PvThydroCmap, rasterized=FigMisc.PT_RASTER)
        CpPlot =    axes[1,0].pcolormesh(T_K, Pscaled, Cp, cmap=Color.PvThydroCmap, rasterized=FigMisc.PT_RASTER)
        alphaPlot = axes[1,1].pcolormesh(T_K, Pscaled, alpha, cmap=alphaCmap, rasterized=FigMisc.PT_RASTER)
        sigPlot =   axes[0,1].pcolormesh(T_K, Pscaled, sig, cmap=Color.PvThydroCmap, rasterized=FigMisc.PT_RASTER)
        VPplot =    axes[0,2].pcolormesh(T_K, Pscaled, VP, cmap=Color.PvThydroCmap, rasterized=FigMisc.PT_RASTER)
        VSplot =    axes[1,2].pcolormesh(T_K, Pscaled, VS, cmap=Color.PvThydroCmap, rasterized=FigMisc.PT_RASTER)
        KSplot =    axes[0,3].pcolormesh(T_K, Pscaled, KS, cmap=Color.PvThydroCmap, rasterized=FigMisc.PT_RASTER)
        GSplot =    axes[1,3].pcolormesh(T_K, Pscaled, GS, cmap=Color.PvThydroCmap, rasterized=FigMisc.PT_RASTER)

        # Add colorbars for each plot
        cbars = [
            fig.colorbar(rhoPlot, ax=axes[0,0]),
            fig.colorbar(CpPlot, ax=axes[1,0]),
            fig.colorbar(alphaPlot, ax=axes[1,1]),
            fig.colorbar(sigPlot, ax=axes[0,1]),
            fig.colorbar(VPplot, ax=axes[0,2]),
            fig.colorbar(VSplot, ax=axes[1,2]),
            fig.colorbar(KSplot, ax=axes[0,3]),
            fig.colorbar(GSplot, ax=axes[1,3])
        ]

        # Plot geotherm on top of colormaps
        for eachPlanet in PlanetList:
            # Geotherm curve
            if np.size(PlanetList) > 1:
                thisColor = None
            else:
                thisColor = Color.geothermHydro
            if Planet.Do.NO_DIFFERENTIATION or Planet.Do.PARTIAL_DIFFERENTIATION:
                Pgeo = eachPlanet.P_MPa * FigLbl.PmultHydro
                Tgeo = eachPlanet.T_K
            else:
                Pgeo = eachPlanet.P_MPa[:eachPlanet.Steps.nHydro] * FigLbl.PmultHydro
                Tgeo = eachPlanet.T_K[:eachPlanet.Steps.nHydro]
            [ax.plot(Tgeo, Pgeo, linewidth=Style.LW_geotherm, linestyle=Style.LS_geotherm,
                     color=thisColor, label=eachPlanet.label) for ax in axf]

            if Params.LEGEND and np.size(PlanetList) > 1:
                handles, lbls = axes[-1,0].get_legend_handles_labels()
                axes[0,-1].legend(handles, lbls)


        plt.tight_layout()
        fig.savefig(Params.FigureFiles.vpvtHydro, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Hydrosphere PT properties plot saved to file: {Params.FigureFiles.vpvtHydro}')
        plt.close()

    return


def PlotPvTPerpleX(PlanetList, Params):

    if os.path.dirname(Params.FigureFiles.vpvtPerpleX) != 'Comparison':
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
            if not FigMisc.TEX_INSTALLED:
                FigLbl.StripLatex()

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
        # axes[0,1].set_title()  # This axes intentionally left blank.
        axes[0,2].set_title(FigLbl.VPlabel)
        axes[1,2].set_title(FigLbl.VSlabel)
        axes[0,3].set_title(FigLbl.KSlabel)
        axes[1,3].set_title(FigLbl.GSlabel)

        # Get pressure and temperature bounds and eval pts and set title
        iSil = np.logical_and(Planet.phase >= Constants.phaseSil,
                              Planet.phase < Constants.phaseSil + 10)
        if INCLUDING_CORE:
            if Params.TITLES:
                fig.suptitle(f'{Planet.name}{FigLbl.PvTtitleCore}')
            iCore = Planet.phase >= Constants.phaseFe
            iInner = np.logical_or(iSil, iCore)
        else:
            if Params.TITLES:
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
        rhoPlotSil =   axes[0,0].pcolormesh(Tinner, PsilScaled, rhoSil, cmap=Color.PvTsilCmap, rasterized=FigMisc.PT_RASTER)
        CpPlotSil =    axes[1,0].pcolormesh(Tinner, PsilScaled, CpSil, cmap=Color.PvTsilCmap, rasterized=FigMisc.PT_RASTER)
        alphaPlotSil = axes[1,1].pcolormesh(Tinner, PsilScaled, alphaSil, cmap=Color.PvTsilCmap, rasterized=FigMisc.PT_RASTER)
        VPplotSil =    axes[0,2].pcolormesh(Tinner, PsilScaled, VPsil, cmap=Color.PvTsilCmap, rasterized=FigMisc.PT_RASTER)
        VSplotSil =    axes[1,2].pcolormesh(Tinner, PsilScaled, VSsil, cmap=Color.PvTsilCmap, rasterized=FigMisc.PT_RASTER)
        KSplotSil =    axes[0,3].pcolormesh(Tinner, PsilScaled, KSsil, cmap=Color.PvTsilCmap, rasterized=FigMisc.PT_RASTER)
        GSplotSil =    axes[1,3].pcolormesh(Tinner, PsilScaled, GSsil, cmap=Color.PvTsilCmap, rasterized=FigMisc.PT_RASTER)

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
            rhoPlotCore =   axes[0,0].pcolormesh(Tinner, PcoreScaled, rhoCore, cmap=Color.PvTcoreCmap, rasterized=FigMisc.PT_RASTER)
            CpPlotCore =    axes[1,0].pcolormesh(Tinner, PcoreScaled, CpCore, cmap=Color.PvTcoreCmap, rasterized=FigMisc.PT_RASTER)
            alphaPlotCore = axes[1,1].pcolormesh(Tinner, PcoreScaled, alphaCore, cmap=Color.PvTcoreCmap, rasterized=FigMisc.PT_RASTER)
            VPplotCore =    axes[0,2].pcolormesh(Tinner, PcoreScaled, VPcore, cmap=Color.PvTcoreCmap, rasterized=FigMisc.PT_RASTER)
            VSplotCore =    axes[1,2].pcolormesh(Tinner, PcoreScaled, VScore, cmap=Color.PvTcoreCmap, rasterized=FigMisc.PT_RASTER)
            KSplotCore =    axes[0,3].pcolormesh(Tinner, PcoreScaled, KScore, cmap=Color.PvTcoreCmap, rasterized=FigMisc.PT_RASTER)
            GSplotCore =    axes[1,3].pcolormesh(Tinner, PcoreScaled, GScore, cmap=Color.PvTcoreCmap, rasterized=FigMisc.PT_RASTER)

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
                 color=Color.geothermInner) for ax in axf]

        plt.tight_layout()
        fig.savefig(Params.FigureFiles.vpvtPerpleX, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Silicate/core PT properties plot saved to file: {Params.FigureFiles.vpvtPerpleX}')
        plt.close()

    return
