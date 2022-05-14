import os
import numpy as np
import logging as log
import matplotlib.pyplot as plt
import matplotlib.colorbar as mcbar
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Wedge
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp1d
from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc, InductParams as IndParams
from PlanetProfile.Thermodynamics.RefProfiles.RefProfiles import CalcRefProfiles, ReloadRefProfiles
from PlanetProfile.Utilities.SetupInit import SetupFilenames
from PlanetProfile.Utilities.defineStructs import Constants

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
                    if FigMisc.REFS_IN_LEGEND and i == 0: thisRef.set_label(wList)
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
        fig.legend(loc=FigMisc.hydroLegendPos, bbox_to_anchor=FigMisc.hydroLegendBox)
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
            Planet.Core.xS_ppt = (100 -int(Planet.Core.coreEOS[2:5])) * 10
            coreLine = f'\ce{{Fe}} core with \SI{{{Planet.Core.xS_ppt / FigLbl.xDiv}}}{{{FigLbl.xUnits}}}~\ce{{S}}'
        elif 'undifferentiated' in Planet.Sil.EOS.comp:
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
                if FigMisc.WEDGE_ICE_TICKS:
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
        cbar[comps[0]] = mcbar.ColorbarBase(cbarAx, cmap=Color.cmap[comps[0]],
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
            cbar[comp] = mcbar.ColorbarBase(cbarAx, cmap=Color.cmap[comp],
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
            cbar = mcbar.ColorbarBase(cbarAx, cmap=Color.cmap[comp],
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

    # Adjust phi values in case we're plotting void volume % instead of void fraction
    if Params.Induct.inductOtype == 'phi':
        Induction.y = Induction.y * FigLbl.phiMult

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
