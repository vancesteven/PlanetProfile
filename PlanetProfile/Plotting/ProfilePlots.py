import os
import numpy as np
import logging
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Wedge
from scipy.interpolate import interp1d
from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc
from PlanetProfile.Plotting.PTPlots import PlotHydroPhase, PlotPvThydro, PlotPvTPerpleX, PlotHydrosphereSpecies, PlotIsoThermalPvThydro
from PlanetProfile.Thermodynamics.RefProfiles.RefProfiles import CalcRefProfiles, ReloadRefProfiles
from PlanetProfile.Utilities.Indexing import GetPhaseIndices, PhaseConv
from PlanetProfile.Utilities.defineStructs import Constants
from PlanetProfile.Thermodynamics.Reaktoro.CustomSolution import SetupCustomSolutionPlotSettings

# Assign logger
log = logging.getLogger('PlanetProfile')

def GeneratePlots(PlanetList, Params):
    # Remove latex styling from legend labels if Latex is not installed
    if not FigMisc.TEX_INSTALLED:
        for Planet in PlanetList:
            Planet.label = FigLbl.StripLatexFromString(Planet.label)
    # Generate CustomSolution plot settings
    if any('CustomSolution' in Planet.Ocean.comp for Planet in PlanetList):
        Params = SetupCustomSolutionPlotSettings(np.array([Planet.Ocean.comp for Planet in PlanetList]), Params)

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
    if Params.PLOT_HYDROSPHERE and np.any([not Planet.Do.NO_OCEAN for Planet in PlanetList]):
        PlotHydrosphereProps(PlanetList, Params)
    if Params.PLOT_SPECIES_HYDROSPHERE and np.any([not Planet.Do.NO_OCEAN for Planet in PlanetList]):
        PlotHydrosphereSpecies(PlanetList, Params)
   # if Params.PLOT_CUSTOMSOLUTION_EOS_PROPERTIES_TABLE and np.any([not Planet.Do.NO_OCEAN for Planet in PlanetList]) and np.any(
      #      ["CustomSolution" in Planet.Ocean.comp for Planet in PlanetList]):
      #  PlotCustomSolutionProperties(PlanetList, Params)
    if Params.PLOT_TRADEOFF:
        PlotSilTradeoff(PlanetList, Params)
        if np.any([Planet.Do.Fe_CORE for Planet in PlanetList]):
            PlotCoreTradeoff(PlanetList, Params)
    if Params.PLOT_POROSITY and np.any([Planet.Do.POROUS_ROCK or Planet.Do.POROUS_ICE for Planet in PlanetList]):
        PlotPorosity(PlanetList, Params)
    if Params.PLOT_SEISMIC and Params.CALC_SEISMIC:
        PlotSeismic(PlanetList, Params)
    if Params.PLOT_VISCOSITY and Params.CALC_VISCOSITY:
        PlotViscosity(PlanetList, Params)
    if Params.PLOT_WEDGE:
        PlotWedge(PlanetList, Params)
    if Params.PLOT_HYDRO_PHASE and np.any([not Planet.Do.NO_H2O for Planet in PlanetList]):
        PlotHydroPhase(PlanetList, Params)
    if Params.PLOT_PVT_HYDRO and np.any([not Planet.Do.NO_H2O for Planet in PlanetList]):
        PlotPvThydro(PlanetList, Params)
    if Params.PLOT_PVT_ISOTHERMAL_HYDRO and np.any([not Planet.Do.NO_H2O for Planet in PlanetList]):
        PlotIsoThermalPvThydro(PlanetList, Params)
    if Params.PLOT_PVT_INNER and not Params.SKIP_INNER:
        PlotPvTPerpleX(PlanetList, Params)

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
    if Params.TITLES:
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

    if FigMisc.FORCE_0_EDGES:
        [ax.set_ylim(bottom=0) for ax in axes]
        [ax.set_xlim(left=0) for ax in axes]

    if Params.LEGEND:
        axes[1].legend()

    plt.tight_layout()
    fig.savefig(Params.FigureFiles.vgrav, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
    log.debug(f'Gravity and pressure plot saved to file: {Params.FigureFiles.vgrav}')
    plt.close()

    return


def PlotHydrosphereProps(PlanetList, Params):

    vRow = 1
    if Params.PLOT_SIGS and Params.CALC_CONDUCT:
        if FigMisc.lowSigCutoff_Sm is None:
            sigCutoff_Sm = 0
        else:
            sigCutoff_Sm = FigMisc.lowSigCutoff_Sm
        maxSig_Sm = np.max([np.max(Planet.sigma_Sm[:Planet.Steps.nHydro]) for Planet in PlanetList if not Planet.Do.NO_OCEAN])
        if maxSig_Sm > sigCutoff_Sm:
            DO_SIGS = True
            vRow += 1
        else:
            log.warning(f'Attempted to plot conductivities, but no profile had above the cutoff ' +
                        f'setting of {FigMisc.lowSigCutoff_Sm}. Excluding sigma plot.')
            DO_SIGS = False
            axsigz = None
    else:
        DO_SIGS = False
        axsigz = None
        sigCutoff_Sm = None
    if Params.PLOT_SOUNDS and Params.CALC_SEISMIC:
        DO_SOUNDS = True
        vRow += 1
    else:
        DO_SOUNDS = False
        axv = None
    
    # Check if we should plot seismic properties (separate from viscosity now) 
    if Params.PLOT_SEISMIC and Params.CALC_SEISMIC:
        DO_SEISMIC_PROPS = True
        vRow += 1
    else:
        DO_SEISMIC_PROPS = False
        axseismic = None

    # Check if we should plot viscosity (will be in conductivity row)
    DO_VISCOSITY = Params.PLOT_VISCOSITY and Params.CALC_VISCOSITY

    # Check if we should plot pressure vs depth
    DO_PRESSURE = Params.PLOT_PRESSURE_DEPTH

    # Generate canvas and add labels - Always 6×7 grid
    fig = plt.figure(figsize=FigSize.vhydro)
    grid = GridSpec(6, 7)

    # Fixed layout: 
    # - Density: top 4 rows (0-3), left 4 columns (0-3)
    # - Thermal: bottom 2 rows (4-5), positioning depends on pressure plot
    # - Pressure: bottom 2 rows (4-5), next 2 columns (2-3) - if enabled
    # - Properties: right 3 columns (4-6), distributed vertically
    
    # Count property plots to determine right column layout
    right_plots = []
    if DO_SOUNDS:
        right_plots.append('sounds')
    if DO_SEISMIC_PROPS:
        right_plots.append('seismic')
    if DO_SIGS:
        right_plots.append('sigs')
    
    num_right_plots = len(right_plots)
    
    if num_right_plots == 0:
        # No property plots - density spans top, thermal and pressure span bottom
        axPrho = fig.add_subplot(grid[0:4, :])      # Top 4 rows, all columns
        if DO_PRESSURE:
            axTz = fig.add_subplot(grid[4:6, :4])   # Bottom 2 rows, left 4 columns
            axPz = fig.add_subplot(grid[4:6, 4:])   # Bottom 2 rows, right 3 columns
        else:
            axTz = fig.add_subplot(grid[4:6, :])    # Bottom 2 rows, all columns
            axPz = None
    else:
        # Property plots present - left 4 columns for density/thermal/pressure, right 3 for properties
        axPrho = fig.add_subplot(grid[0:4, 0:4])    # Top 4 rows, left 4 columns
        if DO_PRESSURE:
            axTz = fig.add_subplot(grid[4:6, 0:2])  # Bottom 2 rows, first 2 columns
            axPz = fig.add_subplot(grid[4:6, 2:4])  # Bottom 2 rows, next 2 columns
        else:
            axTz = fig.add_subplot(grid[4:6, 0:4])  # Bottom 2 rows, all left 4 columns
            axPz = None

    axPrho.set_xlabel(FigLbl.rhoLabel)
    if FigMisc.PLOT_DENSITY_VERSUS_DEPTH:
        axPrho.set_ylabel(FigLbl.zLabel)
    else:
        axPrho.set_ylabel(FigLbl.PlabelHydro)
    axPrho.invert_yaxis()
    
    axTz.set_xlabel(FigLbl.Tlabel)
    axTz.set_ylabel(FigLbl.zLabel)
    axTz.invert_yaxis()
    
    zMax = np.max([Planet.z_m[Planet.Steps.nHydro-1]/1e3 for Planet in PlanetList if not Planet.Do.NO_H2O], initial=0) * 1.05
    axTz.set_ylim([zMax, 0])
    
    axes = [axPrho, axTz]
    
    if DO_PRESSURE:
        axPz.set_xlabel(FigLbl.PlabelHydro)
        axPz.invert_yaxis()
        axPz.set_ylim([zMax, 0])
        axes.append(axPz)

    # Add property plots to right 3 columns (4-6)
    if num_right_plots > 0:
        # Determine row distribution for property plots
        if num_right_plots == 1:
            # 1 property: spans all 6 rows
            row_ranges = [(0, 6)]
        elif num_right_plots == 2:
            # 2 properties: each gets 3 rows
            row_ranges = [(0, 3), (3, 6)]
        elif num_right_plots == 3:
            # 3 properties: each gets 2 rows
            row_ranges = [(0, 2), (2, 4), (4, 6)]
        
        for i, plot_type in enumerate(right_plots):
            start_row, end_row = row_ranges[i]
            
            if plot_type == 'sounds':
                axv = [fig.add_subplot(grid[start_row:end_row, j]) for j in range(4, 7)]
                axv[0].set_xlabel(FigLbl.vPoceanLabel)
                axv[1].set_xlabel(FigLbl.vPiceLabel)
                axv[2].set_xlabel(FigLbl.vSiceLabel)
                [ax.invert_yaxis() for ax in axv]
                [ax.set_ylim([zMax, 0]) for ax in axv]
                axes = axes + axv
                
            elif plot_type == 'seismic':
                # Split seismic into 3 plots: Ocean KS, Ice KS, Ice GS
                axseismic = [fig.add_subplot(grid[start_row:end_row, j]) for j in range(4, 7)]
                axseismic[0].set_xlabel(FigLbl.KSoceanLabel)
                axseismic[1].set_xlabel(FigLbl.KSiceLabel)
                axseismic[2].set_xlabel(FigLbl.GSiceLabel)
                [ax.invert_yaxis() for ax in axseismic]
                [ax.set_ylim([zMax, 0]) for ax in axseismic]
                axes = axes + axseismic
                
            elif plot_type == 'sigs':
                # Conductivity plot - check if viscosity should be added
                if DO_VISCOSITY:
                    # Split conductivity row: columns 4-5 for sigma, column 6 for viscosity
                    axsigz = fig.add_subplot(grid[start_row:end_row, 4:6])
                    axviscz = fig.add_subplot(grid[start_row:end_row, 6:7])
                    axviscz.set_xlabel(FigLbl.etaLabel)
                    axviscz.invert_yaxis()
                    axviscz.set_xscale('log')
                    axviscz.set_ylim([zMax, 0])
                    axes.append(axviscz)
                else:
                    # Standard conductivity plot spanning all 3 right columns
                    axsigz = fig.add_subplot(grid[start_row:end_row, 4:7])
                axsigz.set_xlabel(FigLbl.sigLabel)
                axsigz.invert_yaxis()
                if FigMisc.LOG_SIG:
                    axsigz.set_xscale('log')
                axsigz.set_ylim([zMax, 0])
                axes.append(axsigz)

    if Style.GRIDS:
        [ax.grid() for ax in axes]
        [ax.set_axisbelow(True) for ax in axes]

    if Params.TITLES:
        if Params.ALL_ONE_BODY:
            fig.suptitle(f'{PlanetList[0].name}{FigLbl.hydroTitle}', fontsize=FigLbl.hydroTitleSize)
        else:
            fig.suptitle(FigLbl.hydroCompareTitle, fontsize=FigLbl.hydroTitleSize)

    # Plot reference profiles first, so they plot on bottom of everything
    # Ensure that we only have one unique CustomSolution identifier in comps
    ###TODO FIX THIS CODE
    comps = []
    for Planet in PlanetList:
        if "CustomSolution" in Planet.Ocean.comp:
            if Planet.Ocean.comp.split('=')[0] not in comps:
                comps.append(Planet.Ocean.comp)
        else:
            comps.append(Planet.Ocean.comp)
    comps = np.unique(comps)
    if Params.PLOT_REF:
        # Keep track of which reference profiles have been plotted so that we do each only once
        newRef = {comp:True for comp in comps}

        # Get max pressure among all profiles so we know how far out to plot refs
        Plist = np.concatenate([Planet.P_MPa[:Planet.Steps.nHydro] for Planet in PlanetList])
        Pmax_MPa = np.max(Plist)

        for Planet in PlanetList:
            if newRef[Planet.Ocean.comp] and Planet.Ocean.comp != 'none':
                # Get strings for referencing and labeling
                # If using CustomSolution, then adjust label so compatible with Latex formating
                if "CustomSolution" in Planet.Ocean.comp:
                    wList = f"$\\rho_\mathrm{{melt}}$ \ce{{{Planet.Ocean.comp.split('=')[0].strip()}}} \\{{"
                else:
                    wList = f"$\\rho_\mathrm{{melt}}$ \ce{{{Planet.Ocean.comp}}} \\{{"
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
    rhodots_kgm3 = np.empty(np.size(PlanetList))
    conddots_Sm = np.empty(np.size(PlanetList))
    Tdots_K = np.empty(np.size(PlanetList))
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
    for i,Planet in enumerate(PlanetList):
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
            if FigMisc.PLOT_DENSITY_VERSUS_DEPTH:
                # Plot density vs. depth for hydrosphere
                density = axPrho.plot(Planet.rho_kgm3[:Planet.Steps.nHydro], Planet.z_m[:Planet.Steps.nHydro] / 1e3,
                            label=legLbl, color=thisColor, linewidth=thisLW,
                            linestyle=Style.LS[Planet.Ocean.comp])
                # Make a dot at the end of the thermal profile, if there's an ocean
                if Planet.Steps.nHydro > 0:
                    rhodots_kgm3[i] = np.max(Planet.rho_kgm3[:Planet.Steps.nHydro])
                    axPrho.scatter(rhodots_kgm3[i],
                                np.max(Planet.z_m[:Planet.Steps.nHydro]/1e3),
                                color=density[-1].get_color(), edgecolors=density[-1].get_color(),
                                marker=Style.MS_hydro, s=Style.MW_hydro**2*thisLW)
            else:
                # Plot density vs. pressure curve for hydrosphere
                axPrho.plot(Planet.rho_kgm3[:Planet.Steps.nHydro],
                            Planet.P_MPa[:Planet.Steps.nHydro] * FigLbl.PmultHydro, label=legLbl, color=thisColor,
                            linewidth=thisLW, linestyle=Style.LS[Planet.Ocean.comp])
            # Plot thermal profile vs. depth in hydrosphere
            therm = axTz.plot(Planet.T_K[:Planet.Steps.nHydro] - FigLbl.Tsub,
                              Planet.z_m[:Planet.Steps.nHydro]/1e3,
                              color=thisColor, linewidth=thisLW,
                              linestyle=Style.LS[Planet.Ocean.comp])
            # Make a dot at the end of the thermal profile, if there's an ocean
            if Planet.Steps.nHydro > 0:
                Tdots_K[i] = np.max(Planet.T_K[:Planet.Steps.nHydro] - FigLbl.Tsub)
                axTz.scatter(Tdots_K[i],
                             np.max(Planet.z_m[:Planet.Steps.nHydro]/1e3),
                             color=therm[-1].get_color(), edgecolors=therm[-1].get_color(),
                             marker=Style.MS_hydro, s=Style.MW_hydro**2*thisLW)

            # Plot pressure profile vs. depth in hydrosphere
            if DO_PRESSURE:
                press = axPz.plot(Planet.P_MPa[:Planet.Steps.nHydro] * FigLbl.PmultHydro,
                                  Planet.z_m[:Planet.Steps.nHydro]/1e3,
                                  color=thisColor, linewidth=thisLW,
                                  linestyle=Style.LS[Planet.Ocean.comp])
                # Make a dot at the end of the pressure profile, if there's an ocean
                if Planet.Steps.nHydro > 0:
                    axPz.scatter(np.max(Planet.P_MPa[:Planet.Steps.nHydro] * FigLbl.PmultHydro),
                                 np.max(Planet.z_m[:Planet.Steps.nHydro]/1e3),
                                 color=press[-1].get_color(), edgecolors=press[-1].get_color(),
                                 marker=Style.MS_hydro, s=Style.MW_hydro**2*thisLW)

            if DO_SIGS or DO_SOUNDS or DO_SEISMIC_PROPS:
                indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund, \
                indsClath, indsClathWet, indsMixedClathrateIh, indsMixedClathrateII, indsMixedClathrateIII, indsMixedClathrateV, indsMixedClathrateVI, \
                indsMixedClathrateIhwet, indsMixedClathrateIIund, indsMixedClathrateIIIund, indsMixedClathrateVund, indsMixedClathrateVIund, \
                indsSil, indsSilLiq, _, _, _, _, _, _ = GetPhaseIndices(Planet.phase)

                indsIce = np.sort(np.concatenate((indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund,
                                                  indsV, indsVund, indsVI, indsVIund, indsClath, indsClathWet, 
                                                  indsMixedClathrateIh, indsMixedClathrateII, indsMixedClathrateIII, indsMixedClathrateV, indsMixedClathrateVI, 
                                                  indsMixedClathrateIhwet, indsMixedClathrateIIund, indsMixedClathrateIIIund, indsMixedClathrateVund, indsMixedClathrateVIund)))

                if DO_SIGS:
                    # Plot electrical conductivity vs. depth for hydrosphere
                    sigma_Sm = Planet.sigma_Sm[:Planet.Steps.nHydro] + 0
                    z_km = Planet.z_m[:Planet.Steps.nHydro]/1e3
                    if not FigMisc.SHOW_ICE_CONDUCT:
                        sigma_Sm[indsIce] = np.nan
                    if Planet.Do.POROUS_ICE:
                        indsWet = np.sort(np.concatenate((indsIwet, indsII, indsIII, indsV, indsVI, indsClathWet, 
                                                          indsMixedClathrateIhwet, indsMixedClathrateII, indsMixedClathrateIII, indsMixedClathrateV, indsMixedClathrateVI)))
                        sigma_Sm[indsWet] = Planet.sigma_Sm[indsWet]
                    sigma_Sm[sigma_Sm < sigCutoff_Sm] = np.nan
                    conductivity = axsigz.plot(sigma_Sm, z_km,
                                color=thisColor, linewidth=thisLW,
                                linestyle=Style.LS[Planet.Ocean.comp])
                     # Make a dot at the end of the thermal profile, if there's an ocean
                    if Planet.Steps.nHydro > 0:
                        conddots_Sm[i] = np.nanmax(sigma_Sm)
                        axsigz.scatter(conddots_Sm[i],
                                    np.max(Planet.z_m[:Planet.Steps.nHydro]/1e3),
                                    color=conductivity[-1].get_color(), edgecolors=conductivity[-1].get_color(),
                                    marker=Style.MS_hydro, s=Style.MW_hydro**2*thisLW)

                if DO_SOUNDS:
                    # Plot sound speeds in ocean and ices vs. depth in hydrosphere
                    indsHydro = np.sort(np.concatenate((indsIce, indsLiq)))
                    VPice = Planet.Seismic.VP_kms[indsHydro]
                    VSice = Planet.Seismic.VS_kms[indsHydro]
                    VPliq = VPice + 0
                    # Set non-matching values to nan to avoid gap-spanning lines in plots
                    VPliq[indsIce] = np.nan
                    VPice[indsLiq] = np.nan
                    VSice[indsLiq] = np.nan
                    axv[0].plot(VPliq, Planet.z_m[indsHydro]/1e3,
                                color=thisColor, linewidth=Style.LW_sound,
                                linestyle=Style.LS[Planet.Ocean.comp])
                    axv[1].plot(VPice, Planet.z_m[indsHydro]/1e3,
                                color=thisColor, linewidth=Style.LW_sound,
                                linestyle=Style.LS[Planet.Ocean.comp])
                    axv[2].plot(VSice, Planet.z_m[indsHydro]/1e3,
                                color=thisColor, linewidth=Style.LW_sound,
                                linestyle=Style.LS[Planet.Ocean.comp])

                if DO_SEISMIC_PROPS:
                    # Plot seismic properties (GS, KS) and viscosity vs. depth in hydrosphere
                    # Safely concatenate indices, handling cases where one might be empty
                    if np.size(indsIce) > 0 and np.size(indsLiq) > 0:
                        indsHydro = np.sort(np.concatenate((indsIce, indsLiq)))
                    elif np.size(indsIce) > 0:
                        indsHydro = indsIce
                    elif np.size(indsLiq) > 0:
                        indsHydro = indsLiq
                    else:
                        # Use all hydrosphere indices as fallback
                        indsHydro = np.arange(Planet.Steps.nHydro)
                    
                    z_vals = Planet.z_m[indsHydro]/1e3
                    
                    # Plot bulk modulus KS and shear modulus GS if seismic calculations were done
                    if Params.CALC_SEISMIC:
                        GS_vals = Planet.Seismic.GS_GPa[indsHydro]
                        KS_vals = Planet.Seismic.KS_GPa[indsHydro]
                        
                        # Separate ocean and ice phases for different plots
                        ocean_KS = KS_vals.copy()
                        ice_KS = KS_vals.copy()
                        ice_GS = GS_vals.copy()
                        
                        # Set values to NaN for inappropriate phases
                        if np.size(indsLiq) > 0:
                            # Find which of indsHydro correspond to liquid phases
                            liquid_mask = np.isin(indsHydro, indsLiq)
                            # Ocean plot shows liquid phases only
                            ocean_KS[~liquid_mask] = np.nan
                            # Ice plots show solid phases only
                            ice_KS[liquid_mask] = np.nan
                            ice_GS[liquid_mask] = np.nan
                        
                        if np.size(indsIce) > 0:
                            # Find which of indsHydro correspond to ice phases
                            ice_mask = np.isin(indsHydro, indsIce)
                            # Ocean plot shows liquid phases only
                            ocean_KS[ice_mask] = np.nan
                        
                        # Plot to appropriate subplots
                        axseismic[0].plot(ocean_KS, z_vals,
                                          color=thisColor, linewidth=thisLW,
                                          linestyle=Style.LS[Planet.Ocean.comp])
                        axseismic[1].plot(ice_KS, z_vals,
                                          color=thisColor, linewidth=thisLW,
                                          linestyle=Style.LS[Planet.Ocean.comp])
                        axseismic[2].plot(ice_GS, z_vals,
                                          color=thisColor, linewidth=thisLW,
                                          linestyle=Style.LS[Planet.Ocean.comp])
                    
                    # Plot viscosity if viscosity calculations were done and viscosity plotting is enabled
                    if DO_VISCOSITY and Params.CALC_VISCOSITY:
                        eta_vals = Planet.eta_Pas[indsHydro]
                        # Handle potential NaN values and zero/negative values for log scale
                        eta_plot = eta_vals.copy()
                        eta_plot[eta_plot <= 0] = np.nan
                        
                        axviscz.plot(eta_plot, z_vals,
                                     color=thisColor, linewidth=thisLW,
                                     linestyle=Style.LS[Planet.Ocean.comp])


    if FigMisc.FORCE_0_EDGES:
        axPrho.set_ylim(top=0)
        if DO_PRESSURE:
            axPz.set_ylim(top=0)

    # Limit Tmin so the relevant plot can better show what's going on in the ocean
    Tmax = np.max(Tdots_K)
    Tlims = [FigMisc.TminHydro, FigMisc.TminHydro + 1.05*(Tmax - FigMisc.TminHydro)]
    axTz.set_xlim([np.min(Tlims), np.max(Tlims)])

    if FigMisc.PHASE_LABELS:
        # Label the phases found in the hydrosphere
        phases = np.concatenate([Planet.phase[:Planet.Steps.nHydro] for Planet in PlanetList])
        Pall_MPa = np.concatenate([Planet.P_MPa[:Planet.Steps.nHydro] for Planet in PlanetList])
        rhoAll_kgm3 = np.concatenate([Planet.rho_kgm3[:Planet.Steps.nHydro] for Planet in PlanetList])
        phaseList = np.unique(phases)
        rhoRange = np.diff(axPrho.get_xlim())[0]
        for phase in phaseList:
            if phase < 0 or (phase == 1 and np.any([Planet.Do.POROUS_ICE for Planet in PlanetList])) \
                    or (phase == 5 and np.all([Planet.THIN_OCEAN for Planet in PlanetList])):
                # Underplate layers. Set text on the left side
                adj = -0.06*rhoRange
            else:
                adj = 0.06*rhoRange
            Padj = 0.03*np.abs(np.diff(axPrho.get_ylim())[0])

            if phase == 0:
                # loc = np.where(Pall_MPa == np.min(Pall_MPa[phases == 0]))[0]
                # if np.size(loc) > 1: loc = loc[0]
                rhoLoc_kgm3 = np.mean(rhoAll_kgm3[phases == 0])
                Ploc_MPa = np.mean(Pall_MPa[phases == 0])
                # Shift label closer in the case of a thin ocean, where the curve will be short
                if np.all([Planet.THIN_OCEAN for Planet in PlanetList]):
                    adj2 = -adj
                else:
                    adj2 = 0
                axPrho.text(rhoLoc_kgm3 + adj2, Ploc_MPa*FigLbl.PmultHydro,
                            'liquid', ha='center', va='center', fontsize=FigLbl.TS_hydroLabels)
            else:
                if phase >= Constants.phaseClath and phase < Constants.phaseClath + 10:
                    adj2 = 2*adj
                else:
                    adj2 = 0
                rhomin = np.min(rhoAll_kgm3[phases==phase])
                rhomax = np.max(rhoAll_kgm3[phases==phase])
                Pmin = np.min(Pall_MPa[phases==phase])
                Pmax = np.max(Pall_MPa[phases==phase])
                axPrho.text((rhomin+rhomax)/2 + adj + adj2, (Pmin+Pmax)/2*FigLbl.PmultHydro,
                            PhaseConv(phase), ha='center', va='center', fontsize=FigLbl.TS_hydroLabels)

    if Params.LEGEND:
        handles, lbls = axPrho.get_legend_handles_labels()
        axPrho.legend(handles, lbls, loc='upper right')

    # Set y-limits for right column plots
    if DO_SIGS:
        axsigz.set_ylim(top=0)
        if FigMisc.COMMON_ZMAX_SIG:
            axsigz.set_ylim(bottom=zMax)

    if DO_SEISMIC_PROPS:
        [ax.set_ylim(top=0) for ax in axseismic]
        [ax.set_ylim(bottom=zMax) for ax in axseismic]

    if DO_VISCOSITY:
        axviscz.set_ylim(top=0)
        axviscz.set_ylim(bottom=zMax)

    plt.tight_layout()
    fig.savefig(Params.FigureFiles.vhydro, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
    log.debug(f'Hydrosphere plot saved to file: {Params.FigureFiles.vhydro}')
    plt.close()

    return


def PlotCoreTradeoff(PlanetList, Params):
    fig = plt.figure(figsize=FigSize.vcore)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    if Style.GRIDS:
        ax.grid()
        ax.set_axisbelow(True)
    ax.set_xlabel(FigLbl.RsilLabel)
    ax.set_ylabel(FigLbl.RcoreLabel)
    ALL_SAME_CMR2 = np.all([Planet.CMR2str == PlanetList[0].CMR2str for Planet in PlanetList])
    if ALL_SAME_CMR2:
        CMR2str = f', $C/MR^2 = {PlanetList[0].CMR2str}$'
    else:
        CMR2str = ''
    if Params.TITLES:
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
    fig.savefig(Params.FigureFiles.vcore, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
    log.debug(f'Core trade plot saved to file: {Params.FigureFiles.vcore}')
    plt.close()

    return

"""
def PlotCustomSolutionSpeciesTable(PlanetList, Params):
    # Step 1: Get all unique species across all planets
    species_set = set()
    planet_compositions = []
    CustomSolutionPlanets = [Planet for Planet in PlanetList if "CustomSolution" in Planet.Ocean.comp]

    for Planet in CustomSolutionPlanets:
        species = Planet.Ocean.speciation_ratio_per_kg
        planet_compositions.append(species)
        species_set.update(species.keys())

    # Step 2: Create a list of species in sorted order
    sorted_species = sorted(species_set)

    # Step 3: Prepare the data for the DataFrame (rows: species, columns: planets)
    data = []

    for species in sorted_species:
        row = [species]  # Start with the species name as the first column
        for i, planet_compositions_i in enumerate(planet_compositions):
            # If the species exists in this planet's composition, add the value, else add '---'
            row.append(planet_compositions_i.get(species, '---'))
        data.append(row)

    # Step 4: Create a list of column names
    planet_names = [Planet.label 
    columns = ["Species"] + planet_names

    # Step 5: Create a DataFrame from the data
    df = pd.DataFrame(data, columns=columns)
"""
def PlotSilTradeoff(PlanetList, Params):
    fig = plt.figure(figsize=FigSize.vmant)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    if Style.GRIDS:
        ax.grid()
        ax.set_axisbelow(True)

    ax.set_xlabel(FigLbl.RsilLabel)
    ax.set_ylabel(FigLbl.rhoSilLabel)
    ALL_SAME_CMR2 = np.all([Planet.CMR2str == PlanetList[0].CMR2str for Planet in PlanetList])
    if ALL_SAME_CMR2:
        CMR2str = f', $C/MR^2 = {PlanetList[0].CMR2str}$'
    else:
        CMR2str = ''
    if Params.TITLES:
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
    fig.savefig(Params.FigureFiles.vmant, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
    log.debug(f'Mantle trade plot saved to file: {Params.FigureFiles.vmant}')
    plt.close()

    return


def PlotPorosity(PlanetList, Params):

    for Planet in PlanetList:
        nonzeroPhi = Planet.phi_frac > 0
        Planet.phiPlot = Planet.phi_frac[nonzeroPhi]
        Planet.zPlot = Planet.z_m[:-1][nonzeroPhi]
        Planet.Pplot = Planet.P_MPa[nonzeroPhi]
        iPhaseChanges = np.where(abs(np.diff(Planet.phase[nonzeroPhi])) > 10)[0] + 1
        if np.size(iPhaseChanges) > 0:
            # Add a nan so we don't get a line from porous ice to rock if both are modeled
            Planet.phiPlot = np.insert(Planet.phiPlot, iPhaseChanges, np.nan)
            Planet.zPlot = np.insert(Planet.zPlot, iPhaseChanges, np.nan)
            Planet.Pplot = np.insert(Planet.Pplot, iPhaseChanges, np.nan)

    # Plot dual-axis plot for first entry in PlanetList (usually a main profile), unless we are already doing comparison plots
    if os.path.dirname(Params.FigureFiles.vporeDbl) != 'Comparison':
        Planet = PlanetList[0]
        fig = plt.figure(figsize=FigSize.vpore)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)

        # Check that it's worth converting to GPa if that setting has been selected -- reset labels if not
        if Planet.P_MPa[-1] < 100 and FigLbl.PFULL_IN_GPa:
            log.debug('FigLbl.PFULL_IN_GPa is True, but Pmax is less than 0.1 GPa. Pressures will be plotted in MPa.')
            FigLbl.PFULL_IN_GPa = False
            FigLbl.SetUnits()
            if not FigMisc.TEX_INSTALLED:
                FigLbl.StripLatex()

        ax.set_xlabel(FigLbl.phiLabel)
        ax.set_ylabel(FigLbl.zLabel)
        ax.invert_yaxis()
        P_from_z = interp1d(Planet.z_m[:-1]/1e3, Planet.P_MPa*FigLbl.PmultFull, bounds_error=False, fill_value='extrapolate')
        z_from_P = interp1d(Planet.P_MPa*FigLbl.PmultFull, Planet.z_m[:-1]/1e3, bounds_error=False, fill_value='extrapolate')
        Pax = ax.secondary_yaxis('right', functions=(P_from_z, z_from_P))
        Pax.set_ylabel(FigLbl.PlabelFull)
        if Params.TITLES:
            fig.suptitle(f'{Planet.name}{FigLbl.poreTitle}')

        ax.plot(Planet.phiPlot*FigLbl.phiMult, Planet.zPlot/1e3,
                label=Planet.label, linewidth=Style.LW_std)
        if FigMisc.FORCE_0_EDGES:
            ax.set_ylim(top=0)
            ax.set_xlim(left=0)

        # Prevent uniform porosity from overlapping the right border
        phiMax = np.max(Planet.phiPlot)
        if phiMax - np.min(Planet.phiPlot) < 0.05 and phiMax > 0.15:
            ax.set_xlim(right=np.minimum(phiMax*1.3, 1.0))

        if Params.LEGEND:
            ax.legend()

        fig.savefig(Params.FigureFiles.vporeDbl, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
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
    if Params.TITLES:
        if Params.ALL_ONE_BODY:
            fig.suptitle(f'{PlanetList[0].name}{FigLbl.poreTitle}')
        else:
            fig.suptitle(FigLbl.poreCompareTitle)

    for Planet in PlanetList:
        if Planet.Do.POROUS_ROCK or Planet.Do.POROUS_ICE:
            legLbl = Planet.label
            if (not Params.ALL_ONE_BODY) and FigLbl.BODYNAME_IN_LABEL:
                legLbl = f'{Planet.name} {legLbl}'
            axes[0].plot(Planet.phiPlot*FigLbl.phiMult, Planet.zPlot/1e3,
                         label=legLbl, linewidth=Style.LW_std)
            axes[1].plot(Planet.phiPlot*FigLbl.phiMult, Planet.Pplot*FigLbl.PmultFull,
                         label=legLbl, linewidth=Style.LW_std)

    if FigMisc.FORCE_0_EDGES:
        [ax.set_ylim(top=0) for ax in axes]
        [ax.set_xlim(left=0) for ax in axes]

    # Prevent uniform porosity from overlapping the right border
    for Planet in PlanetList:
        phiMax = np.max(Planet.phiPlot)
        if phiMax - np.min(Planet.phiPlot) < 0.05 and phiMax > 0.15:
            [ax.set_xlim(right=np.minimum(phiMax*1.3, 1.0)) for ax in axes]

    if Params.LEGEND:
        axes[1].legend()

    plt.tight_layout()
    fig.savefig(Params.FigureFiles.vpore, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
    log.debug(f'Porosity plot saved to file: {Params.FigureFiles.vpore}')
    plt.close()

    return

def PlotViscosity(PlanetList, Params):

    fig = plt.figure(figsize=FigSize.vvisc)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    if Style.GRIDS:
        ax.grid()
        ax.set_axisbelow(True)

    ax.set_xlabel(FigLbl.etaLabel)
    ax.set_ylabel(FigLbl.rLabel)
    ax.set_xscale('log')
    if Params.TITLES:
        if Params.ALL_ONE_BODY:
            fig.suptitle(f'{PlanetList[0].name}{FigLbl.viscTitle}')
        else:
            fig.suptitle(FigLbl.viscCompareTitle)

    for Planet in PlanetList:
        legLbl = Planet.label
        if (not Params.ALL_ONE_BODY) and FigLbl.BODYNAME_IN_LABEL:
            legLbl = f'{Planet.name} {legLbl}'
        ax.plot(Planet.eta_Pas, Planet.r_m[:-1]/1e3,
                label=legLbl, linewidth=Style.LW_std)

    if FigMisc.FORCE_0_EDGES:
        ax.set_ylim(bottom=0)

    if Params.LEGEND:
        ax.legend()

    plt.tight_layout()
    fig.savefig(Params.FigureFiles.vvisc, format=FigMisc.figFormat, dpi=FigMisc.dpi,
                metadata=FigLbl.meta)
    log.debug(f'Viscosity plot saved to file: {Params.FigureFiles.vvisc}')
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
    if Params.TITLES:
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

    if FigMisc.FORCE_0_EDGES:
        [ax.set_xlim(left=0) for ax in axf if ax.get_xscale() != 'log']
        [ax.set_ylim(bottom=0) for ax in axf]

    plt.tight_layout()
    fig.savefig(Params.FigureFiles.vseis, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
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
    if Params.ALL_ONE_BODY and Params.TITLES:
        title = f'{PlanetList[0].name} {FigLbl.wedgeTitle}'
        if nWedges > 1:
            title = f'{title}s'
            fig.suptitle(f'\\textbf{{{title}}}', fontsize=Style.TS_super)

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


    # Plot each significant layer for each model, from the outside inward
    for Planet, ax, ionoTop_km, ionoBot_km in zip(PlanetList, axes, ionosUpper_km, ionosLower_km):
        
        # Optional boundaries
        if FigMisc.DRAW_CONVECTION_BOUND:
            iceConvBd = Color.wedgeBd
            clathConvBd = Color.wedgeBd
            silConvBd = Color.wedgeBd
        else:
            iceConvBd = Color.iceIcond
            clathConvBd = Color.mixedClathconv if Planet.Do.MIXED_CLATHRATE_ICE else Color.clathConv
            silConvBd = Color.silCondCmap(1.0)

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
        elif Planet.Sil.EOS is not None and 'undifferentiated' in Planet.Sil.EOS.comp and not (
            Planet.Do.NO_DIFFERENTIATION or Planet.Do.PARTIAL_DIFFERENTIATION):
            coreLine = 'undifferentiated'
        else:
            coreLine = ''

        if FigMisc.LABEL_RADII:
            RionosLbl = f'{FigLbl.ionosTickLbl}: '
            RsurfLbl = f'{FigLbl.surfTickLbl}: '
            RclathLbl = f'{FigLbl.clathTickLbl}: '
            RconvLbl = f'{FigLbl.convIceTickLbl}: '
            RoceanLbl = f'{FigLbl.oceanTickLbl}: '
            RmantLbl = f'{FigLbl.mantTickLbl}: '
            RcoreLbl = f'{FigLbl.coreTickLbl}: '
        else:
            RionosLbl = ''
            RsurfLbl = ''
            RclathLbl = ''
            RconvLbl = ''
            RoceanLbl = ''
            RmantLbl = ''
            RcoreLbl = ''

        if 'Comet' in Planet.Sil.mantleEOS:
            silLine = 'Comet 67P'
        else:
            silLine = f'{Planet.Sil.mantleEOS[:2]} chondrite'

        if Planet.Do.NO_H2O:
            wedgeLabel = f'{silLine}\n{coreLine}\n$q_\mathrm{{surf}}$~\SI{{{Planet.Bulk.qSurf_Wm2*1e3}}}{{{FigLbl.fluxUnits}}}'
        else:
            if Planet.Ocean.comp == 'PureH2O':
                compStr = r'Pure \ce{H2O} ocean'
            elif 'CustomSolution' in Planet.Ocean.comp:
                solutionTitle = Planet.Ocean.comp.split('=')[0].strip()
                solutionTitle = solutionTitle.replace("CustomSolution", "")
                compStr = f'{solutionTitle}'
            else:
                compStr = f'\SI{{{Planet.Ocean.wOcean_ppt:.1f}}}{{{FigLbl.wUnits}}}~\ce{{{Planet.Ocean.comp}}}'

            if Planet.Do.NO_DIFFERENTIATION:
                wedgeLabel = f'{silLine}\n$q_\mathrm{{surf}}$' + \
                             f'~\SI{{{Planet.Bulk.qSurf_Wm2*1e3}}}{{{FigLbl.fluxUnits}}}\n' + \
                             f'{compStr}'
            elif Planet.Do.PARTIAL_DIFFERENTIATION:
                if Planet.Do.DIFFERENTIATE_VOLATILES:
                    wedgeLabel = f'Undifferentiated ice+{silLine}\n$q_\mathrm{{surf}}$' + \
                                 f'~\SI{{{Planet.Bulk.qSurf_Wm2*1e3}}}{{{FigLbl.fluxUnits}}}\n' + \
                                 f'{compStr}, $z_b$~\SI{{{Planet.zb_km:.1f}}}{{km}}'
                else:
                    wedgeLabel = f'Partially differentiated ice+{silLine}\n$q_\mathrm{{surf}}$' + \
                                 f'~\SI{{{Planet.Bulk.qSurf_Wm2*1e3}}}{{{FigLbl.fluxUnits}}}\n' + \
                                 f'{compStr}'
            else:
                wedgeLabel = f'{silLine} mantle\n{coreLine}\n{compStr}, $z_b$~\SI{{{Planet.zb_km:.1f}}}{{km}}'

        if Planet.Do.POROUS_ROCK and not (Planet.Do.NO_DIFFERENTIATION
                                          or Planet.Do.PARTIAL_DIFFERENTIATION):
            wedgeLabel = f'Porous {wedgeLabel}'

        if Params.ALL_ONE_BODY and not nWedges == 1:
            indivTitle = wedgeLabel
        else:
            indivTitle = f'\\textbf{{{Planet.name}}}\n{wedgeLabel}'

        if not FigMisc.TEX_INSTALLED:
            indivTitle = FigLbl.StripLatexFromString(indivTitle)
        ax.set_title(indivTitle, fontsize=Style.TS_desc)
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
            rTicks.append(f'{RsurfLbl}{R_km:.0f}')
            rTickRefs.append(R_km/rMax_km)

            # Starting with ice I or clathrates
            if Planet.Do.CLATHRATE:
                clathCond = Color.mixedClathcond if Planet.Do.MIXED_CLATHRATE_ICE else Color.clathCond
                clathConv = Color.mixedClathconv if Planet.Do.MIXED_CLATHRATE_ICE else Color.clathConv
                    
                if Planet.Bulk.clathType == 'top' or Planet.Bulk.clathType == 'whole':
                    # Clathrates at the surface in this case
                    ax.add_patch(Wedge((0.5,0), R_km/rMax_km, ang1, ang2,
                                       width=Planet.eLid_m/1e3/rMax_km,
                                       fc=clathCond, lw=Style.LW_wedge, ec=clathConvBd))
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
                                               fc=clathConv, lw=Style.LW_wedge, ec=clathConvBd))
                else:
                    # Clathrates in an underplate in this case, always conductive
                    # Conductive ice I at the surface
                    ax.add_patch(Wedge((0.5, 0), R_km/rMax_km, ang1, ang2,
                                       width=Planet.dzIceI_km/rMax_km,
                                       fc=Color.iceIcond, lw=Style.LW_wedge, ec=iceConvBd))
                    # Clathrate underplate
                    ax.add_patch(Wedge((0.5,0), (R_km - Planet.zClath_km)/rMax_km, ang1, ang2,
                                       width=Planet.dzClath_km/rMax_km,
                                       fc=clathCond, lw=Style.LW_wedge, ec=clathConvBd))
                # Outer boundary around clathrates
                ax.add_patch(Wedge((0.5,0), (R_km - Planet.zClath_km)/rMax_km, ang1, ang2,
                                   width=Planet.dzClath_km/rMax_km,
                                   fc=Color.none, lw=Style.LW_wedgeMajor, ec=Color.wedgeBd))
            else:
                # Ice Ih at the surface in this case
                # Conductive ice I
                if Planet.Bulk.asymIce is None:
                    ax.add_patch(Wedge((0.5,0), R_km/rMax_km, ang1, ang2,
                                       width=Planet.eLid_m/1e3/rMax_km,
                                       fc=Color.iceIcond, lw=Style.LW_wedge, ec=iceConvBd))
                    # Convective ice I
                    if (Planet.Dconv_m + Planet.deltaTBL_m) > 0:
                        ax.add_patch(Wedge((0.5,0), (R_km - Planet.eLid_m/1e3)/rMax_km, ang1, ang2,
                                           width=(Planet.Dconv_m + Planet.deltaTBL_m)/1e3/rMax_km,
                                           fc=Color.iceIconv, lw=Style.LW_wedge, ec=iceConvBd))
                else:
                    nWdg = np.size(Planet.Bulk.asymIce)
                    angWdg = 2*Style.wedgeAngle_deg/nWdg
                    for iWdg, thickDiff in enumerate(Planet.Bulk.asymIce):
                        ax.add_patch(Wedge((0.5, 0), R_km/rMax_km, ang1 + iWdg*angWdg, ang1 + (iWdg+1)*angWdg,
                                           width=(Planet.dzIceI_km + thickDiff)/rMax_km, zorder=100,
                                           fc=Color.iceIcond, lw=Style.LW_wedge, ec=iceConvBd))
            # Outer boundary around ice I
            if Planet.dzIceI_km > 0:
                if Planet.Bulk.asymIce is None:
                    ax.add_patch(Wedge((0.5,0), (R_km - Planet.zIceI_m/1e3)/rMax_km, ang1, ang2,
                                       width=Planet.dzIceI_km/rMax_km,
                                       fc=Color.none, lw=Style.LW_wedgeMajor, ec=Color.wedgeBd))
                else:
                    for iWdg, thickDiff in enumerate(Planet.Bulk.asymIce):
                        ax.add_patch(Wedge((0.5, 0), R_km/rMax_km, ang1 + iWdg*angWdg, ang1 + (iWdg+1)*angWdg,
                                           width=(Planet.dzIceI_km + thickDiff)/rMax_km, zorder=100,
                                           fc=Color.none, lw=Style.LW_wedgeMajor, ec=Color.wedgeBd))
            # Surface HP ices
            if Planet.dzIceIIIund_km > 0:
                ax.add_patch(Wedge((0.5,0), (R_km - Planet.zIceIIIund_m/1e3)/rMax_km, ang1, ang2,
                                   width=Planet.dzIceIIIund_km/rMax_km,
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
                    rTicks.append(f'{RoceanLbl}{R_km - Planet.zb_km:.0f}')
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
            if Planet.dzIceII_km > 0:
                ax.add_patch(Wedge((0.5,0), (R_km - Planet.zIceII_m/1e3)/rMax_km, ang1, ang2,
                                   width=Planet.dzIceII_km/rMax_km,
                                   fc=Color.iceIII, lw=Style.LW_wedgeMajor, ec=Color.wedgeBd))
            if Planet.dzIceIII_km > 0:
                ax.add_patch(Wedge((0.5,0), (R_km - Planet.zIceIII_m/1e3)/rMax_km, ang1, ang2,
                                   width=Planet.dzIceIII_km/rMax_km,
                                   fc=Color.iceIII, lw=Style.LW_wedgeMajor, ec=Color.wedgeBd))
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
        rTicks.append(f'{RmantLbl}{Planet.Sil.Rmean_m/1e3:.0f}')
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
        # Conductive silicate gradient (beneath porous layer)
        silCondGrad, dz = np.linspace(0, 1, Color.silCondN+1, retstep=True)
        dzSilCond_km = (Planet.Sil.Rmean_m - Planet.Core.Rmean_m) / 1e3 - Planet.dzSilPorous_km
        # Only plot conductive silicate gradient if layer thickness is nonzero
        if dzSilCond_km > 0:
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
            rTicks.append(f'{RcoreLbl}{Planet.Core.Rmean_m/1e3:.0f}')
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
        if FigMisc.MARK_RADII:
            for rTick in rTickRefs:
                ax.axhline(y=rTick, xmin=0, xmax=0.5, color=Color.wedgeMarkRadii,
                           ls=Style.LS_markRadii, lw=Style.LW_markRadii)
        ax.set_yticklabels(np.array(rTicks), fontsize=Style.TS_ticks)
        [ax.spines[side].set_visible(False) for side in ['top', 'right', 'bottom']]
        ax.get_xaxis().set_visible(False)
        ax.set_ylabel(FigLbl.wedgeRadius, fontsize=Style.TS_ticks)
        ax.set_aspect('equal')

    fig.tight_layout()
    fig.savefig(Params.FigureFiles.vwedg, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
    log.debug(f'Wedge plot saved to file: {Params.FigureFiles.vwedg}')
    plt.close()

    return


def PlotExploreOgram(ExplorationList, Params):
    """ For plotting points showing the various models used in making
        exploreogram plots.
    """

    FigLbl.SetExploration(ExplorationList[0].bodyname, ExplorationList[0].xName,
                          ExplorationList[0].yName, ExplorationList[0].zName)
    if not FigMisc.TEX_INSTALLED:
        FigLbl.StripLatex()

    for Exploration in ExplorationList:
        fig = plt.figure(figsize=FigSize.explore)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)

        if Exploration.zName == 'CMR2calc':
            FigLbl.SetExploreTitle(Exploration.bodyname, Exploration.zName, Exploration.CMR2str)

        if Params.TITLES:
            fig.suptitle(FigLbl.explorationTitle)
        ax.set_xlabel(FigLbl.xLabelExplore)
        ax.set_ylabel(FigLbl.yLabelExplore)
        ax.set_xscale(FigLbl.xScaleExplore)
        ax.set_yscale(FigLbl.yScaleExplore)


        x = Exploration.__getattribute__(Exploration.xName)
        if np.issubdtype(x.dtype, np.number):
            x = x * FigLbl.xMultExplore
        else:
            x = np.zeros(x.shape)
            # Loop over each row and assign the same increasing integer
            for i in range(x.shape[0]):
                x[i, :] = i
        y = Exploration.__getattribute__(Exploration.yName)
        if  np.issubdtype(y.dtype, np.number):
            y = y * FigLbl.yMultExplore
        else:
            y = np.zeros(y.shape)
            # Loop over each row and assign the same increasing integer
            for i in range(y.shape[1]):
                y[:, i] = i
        z = Exploration.__getattribute__(Exploration.zName) * FigLbl.zMultExplore

        ax.set_xlim([np.min(x), np.max(x)])
        ax.set_ylim([np.min(y), np.max(y)])
        # Only keep data points for which a valid model was determined
        zShape = np.shape(z)
        z = np.reshape(z, -1).astype(np.float64)
        INVALID = np.logical_not(np.reshape(Exploration.VALID, -1))
        z[INVALID] = np.nan
        # Return data to original organization
        z = np.reshape(z, zShape)
        
        # Calculate valid data range for colorbar
        zValid = z[z == z]  # Exclude NaNs
        if np.size(zValid) > 0:
            vmin = np.min(zValid)
            vmax = np.max(zValid)
        else:
            vmin = vmax = None
            
        mesh = ax.pcolormesh(x, y, z, shading='auto', cmap=Color.cmap['default'], 
                           vmin=vmin, vmax=vmax, rasterized=FigMisc.PT_RASTER)
        cont = ax.contour(x, y, z, colors='black')
        lbls = plt.clabel(cont, fmt=FigLbl.cfmt)
        cbar = fig.colorbar(mesh, ax=ax, format=FigLbl.cbarFmt)
        # Add the min and max values to the colorbar for reading convenience
        # We compare z values to z values to exclude nans from the max finding,
        # exploiting the fact that nan == nan is False.
        if np.size(zValid) > 0:
            mesh.set_clim(vmin=vmin, vmax=vmax)
            # Filter existing ticks to only include those within valid data range
            existing_ticks = cbar.get_ticks()
            valid_ticks = existing_ticks[(existing_ticks >= vmin) & (existing_ticks <= vmax)]
            # Add min and max values to the filtered ticks
            new_ticks = np.insert(np.append(valid_ticks, vmax), 0, vmin)
            cbar.set_ticks(np.unique(new_ticks))
        cbar.set_label(FigLbl.cbarLabelExplore, size=12)

        plt.tight_layout()
        fig.savefig(Params.FigureFiles.explore, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Plot saved to file: {Params.FigureFiles.explore}')
        plt.close()

    # Plot combination
    if Params.COMPARE and np.size(ExplorationList) > 1:
        Exploration = ExplorationList[0]
        fig = plt.figure(figsize=FigSize.explore)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)

        if Params.TITLES:
            fig.suptitle(FigLbl.explorationTitle)
        ax.set_xlabel(FigLbl.xLabelExplore)
        ax.set_ylabel(FigLbl.yLabelExplore)
        ax.set_xscale(FigLbl.xScaleExplore)
        ax.set_yscale(FigLbl.yScaleExplore)

        # Initialize with first exploration's data
        initial_x_values = Exploration.__getattribute__(Exploration.xName)
        if np.issubdtype(initial_x_values.dtype, np.number):
            x = initial_x_values * FigLbl.xMultExplore
        else:
            x = np.zeros(initial_x_values.shape)
            # Loop over each row and assign the same increasing integer
            for i in range(initial_x_values.shape[0]):
                x[i, :] = i
        initial_y_values = Exploration.__getattribute__(Exploration.yName)
        if np.issubdtype(initial_y_values.dtype, np.number):
            y = initial_y_values * FigLbl.yMultExplore
        else:
            y = np.zeros(initial_y_values.shape)
            # Loop over each row and assign the same increasing integer
            for i in range(initial_y_values.shape[1]):
                y[:, i] = i
        z = Exploration.__getattribute__(Exploration.zName) * FigLbl.zMultExplore
        # Only keep data points for which a valid model was determined
        zShape = np.shape(z)
        z = np.reshape(z, -1).astype(np.float64)
        INVALID = np.logical_not(np.reshape(Exploration.VALID, -1))
        z[INVALID] = np.nan
        # Return data to original organization
        z = np.reshape(z, zShape)
        
        for Exploration in ExplorationList[1:]:
            # Get new exploration data
            x_values = Exploration.__getattribute__(Exploration.xName)
            if np.issubdtype(x_values.dtype, np.number):
                new_x = x_values * FigLbl.xMultExplore
            else:
                new_x = np.zeros(x_values.shape)
                for i in range(x_values.shape[0]):
                    new_x[i, :] = i
                    
            y_values = Exploration.__getattribute__(Exploration.yName)
            if np.issubdtype(y_values.dtype, np.number):
                new_y = y_values * FigLbl.yMultExplore
            else:
                new_y = np.zeros(y_values.shape)
                for i in range(y_values.shape[1]):
                    new_y[:, i] = i
                    
            thisz = Exploration.__getattribute__(Exploration.zName) * FigLbl.zMultExplore
            # Only keep data points for which a valid model was determined
            zShape = np.shape(thisz)
            thisz = np.reshape(thisz, -1).astype(np.float64)
            INVALID = np.logical_not(np.reshape(Exploration.VALID, -1))
            thisz[INVALID] = np.nan
            # Return data to original organization
            thisz = np.reshape(thisz, zShape)
            
            # Determine concatenation axis based on which dimension matches
            if new_x.shape[0] == x.shape[0] and (
                (np.issubdtype(x_values.dtype, np.number) and np.allclose(x_values[:, 0], initial_x_values[:, 0], equal_nan=True)) or
                (not np.issubdtype(x_values.dtype, np.number) and np.array_equal(x_values[:, 0], initial_x_values[:, 0]))
            ):
                # x arrays match, concatenate along y-axis (axis=1)
                # Find insertion point to maintain monotonic order
                y_existing = y[0, :]
                y_new_start = new_y[0, 0]
                insert_idx = np.searchsorted(y_existing, y_new_start)
                # Insert new data at the appropriate position
                x = np.concatenate([x[:, :insert_idx], new_x, x[:, insert_idx:]], axis=1)
                y = np.concatenate([y[:, :insert_idx], new_y, y[:, insert_idx:]], axis=1)
                z = np.concatenate([z[:, :insert_idx], thisz, z[:, insert_idx:]], axis=1)
            elif new_y.shape[1] == y.shape[1] and (
                (np.issubdtype(y_values.dtype, np.number) and np.allclose(y_values[0, :], initial_y_values[0, :], equal_nan=True)) or
                (not np.issubdtype(y_values.dtype, np.number) and np.array_equal(y_values[0, :], initial_y_values[0, :]))
            ):
                # y arrays match, concatenate along x-axis (axis=0)
                # Find insertion point to maintain monotonic order
                x_existing = x[:, 0]
                x_new_start = new_x[0, 0]
                insert_idx = np.searchsorted(x_existing, x_new_start)
                # Insert new data at the appropriate position
                x = np.concatenate([x[:insert_idx, :], new_x, x[insert_idx:, :]], axis=0)
                y = np.concatenate([y[:insert_idx, :], new_y, y[insert_idx:, :]], axis=0)
                z = np.concatenate([z[:insert_idx, :], thisz, z[insert_idx:, :]], axis=0)
            else:
                log.warning('The exploreogram comparison plot cannot be made because the x or y axes are not the same. Skipping the comparison plot.')
                break
            
        # Calculate valid data range for colorbar
        zValid = z[z == z]  # Exclude NaNs
        if np.size(zValid) > 0:
            vmin = np.min(zValid)
            vmax = np.max(zValid)
        else:
            vmin = vmax = None
            
        mesh = ax.pcolormesh(x, y, z, shading='auto', cmap=Color.cmap['default'], 
                           vmin=vmin, vmax=vmax, rasterized=FigMisc.PT_RASTER)
        cont = ax.contour(x, y, z, colors='black')
        lbls = plt.clabel(cont, fmt=FigLbl.cfmt)
        cbar = fig.colorbar(mesh, ax=ax, format=FigLbl.cbarFmt)
        # Add the min and max values to the colorbar for reading convenience
        # We compare z values to z values to exclude nans from the max finding,
        # exploiting the fact that nan == nan is False.
        if np.size(zValid) > 0:
            mesh.set_clim(vmin=vmin, vmax=vmax)
            # Filter existing ticks to only include those within valid data range
            existing_ticks = cbar.get_ticks()
            valid_ticks = existing_ticks[(existing_ticks >= vmin) & (existing_ticks <= vmax)]
            # Add min and max values to the filtered ticks
            new_ticks = np.insert(np.append(valid_ticks, vmax), 0, vmin)
            cbar.set_ticks(np.unique(new_ticks))
        cbar.set_label(FigLbl.cbarLabelExplore, size=12)

        plt.tight_layout()
        fig.savefig(Params.FigureFiles.explore, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Plot saved to file: {Params.FigureFiles.explore}')
        plt.close()

    return


def PlotExploreOgramDsigma(ExplorationList, Params):
    """ Plot a scatter showing the evaluated ocean mean conductivity and layer thickness,
        for comparison against canonical D/sigma exploration plots.
    """

    ExplorationList[0].xName = 'D_km'
    ExplorationList[0].yName = 'sigmaMean_Sm'
    ExplorationList[0].zName = 'zb_km'
    FigLbl.SetExploration(ExplorationList[0].bodyname, ExplorationList[0].xName,
                          ExplorationList[0].yName, ExplorationList[0].zName)
    if not FigMisc.TEX_INSTALLED:
        FigLbl.StripLatex()

    for Exploration in (ex for ex in ExplorationList if not ex.NO_H2O):
        Exploration.xName = 'D_km'
        Exploration.yName = 'sigmaMean_Sm'
        Exploration.zName = 'zb_km'
        fig = plt.figure(figsize=FigSize.explore)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)

        if Params.TITLES:
            fig.suptitle(FigLbl.explorationDsigmaTitle)
        ax.set_xlabel(FigLbl.xLabelExplore)
        ax.set_ylabel(FigLbl.yLabelExplore)
        # Override standard settings for this type of plot
        ax.set_xscale('linear')
        ax.set_yscale('log')
        ax.set_ylim(FigMisc.DSIGMA_YLIMS)

        x = np.reshape(Exploration.__getattribute__(Exploration.xName) * FigLbl.xMultExplore, -1)
        y = np.reshape(Exploration.__getattribute__(Exploration.yName) * FigLbl.yMultExplore, -1)
        # Only keep data points for which a valid model was determined
        VALID = np.logical_not(np.logical_or(np.isnan(x), np.isnan(y)))
        x = x[VALID]
        y = y[VALID]
        if np.size(x) > 0:
            ax.set_xlim([np.min(x), np.max(x)])
        linzb = np.reshape(Exploration.zb_km, -1)[VALID]
        
        # Get ocean composition for composition line drawing
        ocean_comp = np.reshape(Exploration.oceanComp, -1)[VALID]
        
        # Draw composition lines if enabled
        if FigMisc.DRAW_COMPOSITION_LINE:
            # Group by unique ocean compositions and plot connecting lines
            unique_comps = np.unique(ocean_comp)
            for comp in unique_comps:
                comp_indices = np.where(ocean_comp == comp)[0]
                x_line, y_line = x[comp_indices], y[comp_indices]
                
                # Sort points by x for connected lines
                sorted_idx = np.argsort(x_line)
                x_line = x_line[sorted_idx]
                y_line = y_line[sorted_idx]
                
                # Set color based on composition
                if FigMisc.MANUAL_HYDRO_COLORS:
                    thisColor = Color.cmap[comp](Color.GetNormT(np.nanmax(linzb[comp_indices])))
                else:
                    thisColor = Color.cmap[comp](0.5)
                
                # Clean composition label
                if 'CustomSolution' in comp:
                    comp_label = comp.split('=')[0].replace('CustomSolution', '')
                else:
                    comp_label = comp
                
                # Plot line for this composition
                ax.plot(x_line, y_line, color=thisColor, linewidth=FigMisc.DSIGMA_COMP_LINE_WIDTH, 
                       alpha=FigMisc.DSIGMA_COMP_LINE_ALPHA, label=f'{comp_label}', zorder=2)
        
        # Handle ice thickness coloring or regular colorbar
        if FigMisc.SHOW_ICE_THICKNESS_DOTS:
            # Get ice shell thickness and normalize for coloring
            ice_thickness = np.reshape(Exploration.zb_km, -1)[VALID]
            if np.size(ice_thickness) > 0:
                # Set bounds for normalization (min and max of ice thickness)
                Tbound_lower = np.min(ice_thickness)
                Tbound_upper = np.max(ice_thickness)
                # Get normalized values using GetNormT
                norm_thickness = Color.GetNormT(ice_thickness, Tbound_lower, Tbound_upper)
            else:
                norm_thickness = None
            
            # Plot scatter with ice thickness coloring
            pts = ax.scatter(x, y, c=norm_thickness,
                            cmap=FigMisc.DSIGMA_ICE_THICKNESS_CMAP, marker=Style.MS_Induction, 
                            s=Style.MW_Induction**2, edgecolors=FigMisc.DSIGMA_DOT_EDGE_COLOR, 
                            linewidths=FigMisc.DSIGMA_DOT_EDGE_WIDTH, zorder=3)
            
            # Create legend showing ice thickness values instead of colorbar
            if np.size(ice_thickness) > 0:
                # Round thickness values to avoid floating-point duplicates
                rounded_thicknesses = np.round(ice_thickness)
                unique_thicknesses = np.unique(rounded_thicknesses)
                
                if len(unique_thicknesses) > 10:  # Limit number of legend entries
                    # Select representative values
                    indices = np.linspace(0, len(unique_thicknesses)-1, 10, dtype=int)
                    selected_thicknesses = unique_thicknesses[indices]
                else:
                    selected_thicknesses = unique_thicknesses
                
                # Create legend elements
                legend_elements = []
                for thickness in selected_thicknesses:
                    norm_val = Color.GetNormT(thickness, Tbound_lower, Tbound_upper)
                    color = plt.cm.Greys(norm_val)
                    legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                                    markerfacecolor=color, markeredgecolor='black',
                                                    markersize=8, label=f'{thickness:.0f} km'))
                
                # Add legend for ice thickness
                thickness_legend = ax.legend(handles=legend_elements, title="Ice Shell Thickness", 
                                           loc='upper right', bbox_to_anchor=(1.0, 1.0), 
                                           fontsize=8, title_fontsize=10)
                ax.add_artist(thickness_legend)  # Keep this legend when adding composition legend
        else:
            # Standard colorbar approach
            pts = ax.scatter(x, y, c=linzb,
                            cmap=Color.cmap[Exploration.oceanComp[0,0]],
                            marker=Style.MS_Induction, s=Style.MW_Induction**2, zorder=3)

            cbar = fig.colorbar(pts, ax=ax, format=FigLbl.cbarFmt)
            # Append the max value to the colorbar for reading convenience
            # We compare z values to z values to exclude nans from the max finding,
            # exploiting the fact that nan == nan is False.
            if np.size(linzb) > 0:
                new_ticks = np.insert(np.append(cbar.get_ticks(), np.max(linzb)), 0, np.min(linzb))
                cbar.set_ticks(np.unique(new_ticks))
            cbar.set_label(FigLbl.cbarLabelExplore, size=12)
        
        # Add legend for composition lines if enabled
        if FigMisc.DRAW_COMPOSITION_LINE and Params.LEGEND:
            if FigMisc.SHOW_ICE_THICKNESS_DOTS:
                # Position composition legend to avoid overlap with ice thickness legend
                ax.legend(title="Ocean Composition", fontsize=FigMisc.DSIGMA_COMP_LEGEND_FONT_SIZE, 
                         title_fontsize=FigMisc.DSIGMA_COMP_LEGEND_TITLE_SIZE,
                         loc='upper left', bbox_to_anchor=(0.0, 1.0))
            else:
                ax.legend(title="Ocean Composition", fontsize=FigMisc.DSIGMA_COMP_LEGEND_FONT_SIZE, 
                         title_fontsize=FigMisc.DSIGMA_COMP_LEGEND_TITLE_SIZE)

        plt.tight_layout()
        fig.savefig(Params.FigureFiles.exploreDsigma, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Plot saved to file: {Params.FigureFiles.exploreDsigma}')
        plt.close()

    # Plot combination
    if Params.COMPARE and np.size(ExplorationList) > 1:
        # Filter out explorations with no H2O
        ValidExplorations = [ex for ex in ExplorationList if not ex.NO_H2O]
        
        if len(ValidExplorations) > 1:
            # Set up exploration parameters for the first valid exploration
            Exploration = ValidExplorations[0]
            Exploration.xName = 'D_km'
            Exploration.yName = 'sigmaMean_Sm'
            Exploration.zName = 'zb_km'
            FigLbl.SetExploration(Exploration.bodyname, Exploration.xName, 
                                  Exploration.yName, Exploration.zName)
            if not FigMisc.TEX_INSTALLED:
                FigLbl.StripLatex()

            fig = plt.figure(figsize=FigSize.explore)
            grid = GridSpec(1, 1)
            ax = fig.add_subplot(grid[0, 0])
            if Style.GRIDS:
                ax.grid()
                ax.set_axisbelow(True)

            if Params.TITLES:
                fig.suptitle(FigLbl.explorationDsigmaTitle)
            ax.set_xlabel(FigLbl.xLabelExplore)
            ax.set_ylabel(FigLbl.yLabelExplore)
            # Override standard settings for this type of plot
            ax.set_xscale('linear')
            ax.set_yscale('log')
            ax.set_ylim(FigMisc.DSIGMA_YLIMS)

            # Collect data from all explorations
            all_x = []
            all_y = []
            all_zb = []
            all_ocean_comp = []
            
            for Exploration in ValidExplorations:
                x = np.reshape(Exploration.__getattribute__(Exploration.xName) * FigLbl.xMultExplore, -1)
                y = np.reshape(Exploration.__getattribute__(Exploration.yName) * FigLbl.yMultExplore, -1)
                # Only keep data points for which a valid model was determined
                VALID = np.logical_not(np.logical_or(np.isnan(x), np.isnan(y)))
                x, y = x[VALID], y[VALID]
                
                # Get zb_km and ocean composition data
                zb = np.reshape(Exploration.zb_km, -1)[VALID]
                ocean_comp = np.reshape(Exploration.oceanComp, -1)[VALID]
                
                all_x.extend(x)
                all_y.extend(y)
                all_zb.extend(zb)
                all_ocean_comp.extend(ocean_comp)

            # Convert to numpy arrays
            all_x = np.array(all_x)
            all_y = np.array(all_y)
            all_zb = np.array(all_zb)
            all_ocean_comp = np.array(all_ocean_comp)

            # Draw composition lines if enabled
            if FigMisc.DRAW_COMPOSITION_LINE:
                # Group by unique ocean compositions and plot connecting lines
                unique_comps = np.unique(all_ocean_comp)
                plotted_labels = set()  # Keep track of which labels have been plotted
                
                for comp in unique_comps:
                    comp_indices = np.where(all_ocean_comp == comp)[0]
                    x_line, y_line = all_x[comp_indices], all_y[comp_indices]
                    
                    # Sort points by x for connected lines
                    sorted_idx = np.argsort(x_line)
                    x_line = x_line[sorted_idx]
                    y_line = y_line[sorted_idx]
                    
                    # Set color based on composition
                    if FigMisc.MANUAL_HYDRO_COLORS:
                        thisColor = Color.cmap[comp](Color.GetNormT(np.nanmax(all_zb[comp_indices])))
                    else:
                        thisColor = Color.cmap[comp](0.5)
                    
                    # Clean composition label
                    if 'CustomSolution' in comp:
                        comp_label = comp.split('=')[0].replace('CustomSolution', '')
                    else:
                        comp_label = comp
                    
                    # Only add label if this composition hasn't been plotted yet
                    line_label = comp_label if comp_label not in plotted_labels else None
                    if line_label:
                        plotted_labels.add(comp_label)
                    
                    # Plot line for this composition
                    ax.plot(x_line, y_line, color=thisColor, linewidth=FigMisc.DSIGMA_COMP_LINE_WIDTH, 
                           alpha=FigMisc.DSIGMA_COMP_LINE_ALPHA, label=line_label, zorder=2)

            # Plot scatter points using the same marker for all explorations
            if np.size(all_x) > 0:
                ax.set_xlim([np.min(all_x), np.max(all_x)])
            
            # Handle ice thickness coloring or regular colorbar for combined plot
            if FigMisc.SHOW_ICE_THICKNESS_DOTS:
                # Get ice shell thickness and normalize for coloring
                if np.size(all_zb) > 0:
                    # Set bounds for normalization (min and max of ice thickness across all explorations)
                    Tbound_lower = np.min(all_zb)
                    Tbound_upper = np.max(all_zb)
                    # Get normalized values using GetNormT
                    norm_thickness = Color.GetNormT(all_zb, Tbound_lower, Tbound_upper)
                else:
                    norm_thickness = None
                
                # Plot scatter with ice thickness coloring
                pts = ax.scatter(all_x, all_y, c=norm_thickness,
                                cmap=FigMisc.DSIGMA_ICE_THICKNESS_CMAP, marker=Style.MS_Induction, 
                                s=Style.MW_Induction**2, edgecolors=FigMisc.DSIGMA_DOT_EDGE_COLOR, 
                                linewidths=FigMisc.DSIGMA_DOT_EDGE_WIDTH, zorder=3)
                
                # Create legend showing ice thickness values instead of colorbar
                if np.size(all_zb) > 0:
                    # Round thickness values to avoid floating-point duplicates
                    rounded_thicknesses = np.round(all_zb)
                    unique_thicknesses = np.unique(rounded_thicknesses)
                    
                    if len(unique_thicknesses) > FigMisc.DSIGMA_MAX_LEGEND_ENTRIES:  # Limit number of legend entries
                        # Select representative values
                        indices = np.linspace(0, len(unique_thicknesses)-1, FigMisc.DSIGMA_MAX_LEGEND_ENTRIES, dtype=int)
                        selected_thicknesses = unique_thicknesses[indices]
                else:
                    selected_thicknesses = unique_thicknesses
                
                # Create legend elements
                legend_elements = []
                for thickness in selected_thicknesses:
                    norm_val = Color.GetNormT(thickness, Tbound_lower, Tbound_upper)
                    color = getattr(plt.cm, FigMisc.DSIGMA_ICE_THICKNESS_CMAP)(norm_val)
                    legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                                    markerfacecolor=color, markeredgecolor=FigMisc.DSIGMA_DOT_EDGE_COLOR,
                                                    markersize=8, label=f'{thickness:.0f}'))
                
                # Add legend for ice thickness
                thickness_legend = ax.legend(handles=legend_elements, title="Ice Shell Thickness (km)", 
                                           loc='upper right', bbox_to_anchor=(1.0, 1.0), 
                                           fontsize=FigMisc.DSIGMA_ICE_LEGEND_FONT_SIZE, title_fontsize=FigMisc.DSIGMA_ICE_LEGEND_TITLE_SIZE)
                ax.add_artist(thickness_legend)  # Keep this legend when adding composition legend
            else:
                # Standard colorbar approach
                pts = ax.scatter(all_x, all_y, c=all_zb,
                                cmap=Color.cmap['default'], marker=Style.MS_Induction, 
                                s=Style.MW_Induction**2, zorder=3)

                cbar = fig.colorbar(pts, ax=ax, format=FigLbl.cbarFmt)
                # Append the max value to the colorbar for reading convenience
                if np.size(all_zb) > 0:
                    new_ticks = np.insert(np.append(cbar.get_ticks(), np.max(all_zb)), 0, np.min(all_zb))
                    cbar.set_ticks(np.unique(new_ticks))
                cbar.set_label(FigLbl.cbarLabelExplore, size=12)
            
            # Add legend for composition lines if enabled
            if FigMisc.DRAW_COMPOSITION_LINE and Params.LEGEND:
                if FigMisc.SHOW_ICE_THICKNESS_DOTS:
                    # Position composition legend to avoid overlap with ice thickness legend
                    ax.legend(title="Ocean Composition", fontsize=FigMisc.DSIGMA_COMP_LEGEND_FONT_SIZE, 
                             title_fontsize=FigMisc.DSIGMA_COMP_LEGEND_TITLE_SIZE,
                             bbox_to_anchor=(0.0, 1.0))
                else:
                    ax.legend(title="Ocean Composition", fontsize=FigMisc.DSIGMA_COMP_LEGEND_FONT_SIZE, 
                             title_fontsize=FigMisc.DSIGMA_COMP_LEGEND_TITLE_SIZE)

            plt.tight_layout()
            fig.savefig(Params.FigureFiles.exploreDsigma, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
            log.debug(f'Combined D/sigma plot saved to file: {Params.FigureFiles.exploreDsigma}')
            plt.close()

    return

def PlotExploreOgramZbD(ExplorationList, Params):
    """ Plot a scatter showing ice shell thickness vs ocean thickness with user-configurable z variable,
        with optional ocean composition lines.
    """

    FigLbl.SetExploration(ExplorationList[0].bodyname, 'zb_km', 'D_km', ExplorationList[0].zName)
    if not FigMisc.TEX_INSTALLED:
        FigLbl.StripLatex()

    for Exploration in (ex for ex in ExplorationList if not ex.NO_H2O):
        Exploration.xName = 'zb_km'
        Exploration.yName = 'D_km'
        # zName should already be set by user
        
        fig = plt.figure(figsize=FigSize.explore)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)

        if Params.TITLES:
            fig.suptitle(FigLbl.explorationTitle)
        ax.set_xlabel(FigLbl.xLabelExplore)
        ax.set_ylabel(FigLbl.yLabelExplore)
        ax.set_xscale('linear')
        ax.set_yscale('linear')

        x = np.reshape(Exploration.__getattribute__(Exploration.xName) * FigLbl.xMultExplore, -1)
        y = np.reshape(Exploration.__getattribute__(Exploration.yName) * FigLbl.yMultExplore, -1)
        z = np.reshape(Exploration.__getattribute__(Exploration.zName) * FigLbl.zMultExplore, -1)
        # Only keep data points for which a valid model was determined
        VALID = np.logical_not(np.logical_or(np.isnan(x), np.isnan(y)))
        x, y, z = x[VALID], y[VALID], z[VALID]
        
        if np.size(x) > 0:
            # Add padding to prevent marker cutoff at edges
            x_range = np.max(x) - np.min(x)
            y_range = np.max(y) - np.min(y)
            x_padding = x_range * FigMisc.ZBD_AXIS_PADDING
            y_padding = y_range * FigMisc.ZBD_AXIS_PADDING
            ax.set_xlim([np.min(x) - x_padding, np.max(x) + x_padding])
            ax.set_ylim([np.min(y) - y_padding, np.max(y) + y_padding])

        # Get ocean composition for composition line drawing
        ocean_comp = np.reshape(Exploration.oceanComp, -1)[VALID]
        
        # Draw composition lines if enabled
        if FigMisc.DRAW_COMPOSITION_LINE:
            # Group by unique ocean compositions and plot connecting lines
            unique_comps = np.unique(ocean_comp)
            for comp in unique_comps:
                comp_indices = np.where(ocean_comp == comp)[0]
                x_line, y_line = x[comp_indices], y[comp_indices]
                
                # Sort points by x for connected lines
                sorted_idx = np.argsort(x_line)
                x_line = x_line[sorted_idx]
                y_line = y_line[sorted_idx]
                
                # Set color based on composition
                if FigMisc.MANUAL_HYDRO_COLORS:
                    thisColor = Color.cmap[comp](Color.GetNormT(np.nanmax(z[comp_indices])))
                else:
                    thisColor = Color.cmap[comp](0.5)
                
                # Clean composition label
                if 'CustomSolution' in comp:
                    comp_label = comp.split('=')[0].replace('CustomSolution', '')
                else:
                    comp_label = comp
                
                # Plot line for this composition
                ax.plot(x_line, y_line, color=thisColor, linewidth=FigMisc.ZBD_COMP_LINE_WIDTH, 
                       alpha=FigMisc.ZBD_COMP_LINE_ALPHA, label=f'{comp_label}', zorder=2)

        # Separate NaN and non-NaN points for different plotting
        nan_mask = np.isnan(z)
        valid_mask = ~nan_mask
        
        # Plot non-NaN points colored by z variable
        if np.any(valid_mask):
            z_valid = z[valid_mask]
            pts = ax.scatter(x[valid_mask], y[valid_mask], c=z_valid,
                            cmap=FigMisc.ZBD_COLORMAP, marker=Style.MS_Induction, 
                            s=Style.MW_Induction**2, edgecolors=FigMisc.ZBD_DOT_EDGE_COLOR,
                            linewidths=FigMisc.ZBD_DOT_EDGE_WIDTH, zorder=3, 
                            vmin=np.nanmin(z_valid), vmax=np.nanmax(z_valid))
            
            # Add colorbar
            cbar = fig.colorbar(pts, ax=ax, format=FigLbl.cbarFmt)
            # Extend colorbar to cover the full range of default ticks
            default_ticks = cbar.get_ticks()
            if len(default_ticks) > 0:
                pts.set_clim(vmin=np.min(default_ticks), vmax=np.max(default_ticks))
            cbar.set_label(FigLbl.cbarLabelExplore, size=12)
        
        # Plot NaN points with distinct color
        if np.any(nan_mask):
            ax.scatter(x[nan_mask], y[nan_mask], c=FigMisc.ZBD_NAN_COLOR,
                      marker=FigMisc.ZBD_NAN_MARKER, s=Style.MW_Induction**2, 
                      edgecolors=FigMisc.ZBD_DOT_EDGE_COLOR, linewidths=FigMisc.ZBD_DOT_EDGE_WIDTH,
                      zorder=3, label='Invalid data')
        
        # Add legend for composition lines if enabled
        if FigMisc.DRAW_COMPOSITION_LINE and Params.LEGEND:
            ax.legend(title="Ocean Composition", fontsize=FigMisc.ZBD_COMP_LEGEND_FONT_SIZE, 
                     title_fontsize=FigMisc.ZBD_COMP_LEGEND_TITLE_SIZE)
        elif np.any(nan_mask) and Params.LEGEND:
            # Show legend for invalid data points if no composition lines are shown
            ax.legend(fontsize=FigMisc.ZBD_COMP_LEGEND_FONT_SIZE)

        plt.tight_layout()
        fig.savefig(Params.FigureFiles.exploreZbD, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Plot saved to file: {Params.FigureFiles.exploreZbD}')
        plt.close()

    # Plot combination
    if Params.COMPARE and np.size(ExplorationList) > 1:
        # Filter out explorations with no H2O
        ValidExplorations = [ex for ex in ExplorationList if not ex.NO_H2O]
        
        if len(ValidExplorations) > 1:
            # Set up exploration parameters for the first valid exploration
            Exploration = ValidExplorations[0]
            Exploration.xName = 'zb_km'
            Exploration.yName = 'D_km'
            # zName should already be set by user
            FigLbl.SetExploration(Exploration.bodyname, Exploration.xName, 
                                  Exploration.yName, Exploration.zName)
            if not FigMisc.TEX_INSTALLED:
                FigLbl.StripLatex()

            fig = plt.figure(figsize=FigSize.explore)
            grid = GridSpec(1, 1)
            ax = fig.add_subplot(grid[0, 0])
            if Style.GRIDS:
                ax.grid()
                ax.set_axisbelow(True)

            if Params.TITLES:
                fig.suptitle(FigLbl.explorationTitle)
            ax.set_xlabel(FigLbl.xLabelExplore)
            ax.set_ylabel(FigLbl.yLabelExplore)
            ax.set_xscale(FigLbl.xScaleExplore)
            ax.set_yscale(FigLbl.yScaleExplore)

            # Collect data from all explorations
            all_x = []
            all_y = []
            all_z = []
            all_ocean_comp = []
            
            for Exploration in ValidExplorations:
                x = np.reshape(Exploration.__getattribute__(Exploration.xName) * FigLbl.xMultExplore, -1)
                y = np.reshape(Exploration.__getattribute__(Exploration.yName) * FigLbl.yMultExplore, -1)
                z = np.reshape(Exploration.__getattribute__(Exploration.zName) * FigLbl.zMultExplore, -1)
                # Only keep data points for which a valid model was determined
                VALID = np.logical_not(np.logical_or(np.isnan(x), np.isnan(y)))
                x, y, z = x[VALID], y[VALID], z[VALID]
                
                if np.size(x) > 0:
                    ax.set_xlim([np.min(x), np.max(x)])
                    ax.set_ylim([np.min(y), np.max(y)])
                
                # Get ocean composition data
                ocean_comp = np.reshape(Exploration.oceanComp, -1)[VALID]
                
                all_x.extend(x)
                all_y.extend(y)
                all_z.extend(z)
                all_ocean_comp.extend(ocean_comp)

            # Convert to numpy arrays
            all_x = np.array(all_x)
            all_y = np.array(all_y)
            all_z = np.array(all_z)
            all_ocean_comp = np.array(all_ocean_comp)

            # Draw composition lines if enabled
            if FigMisc.DRAW_COMPOSITION_LINE:
                # Group by unique ocean compositions and plot connecting lines
                unique_comps = np.unique(all_ocean_comp)
                plotted_labels = set()  # Keep track of which labels have been plotted
                
                for comp in unique_comps:
                    comp_indices = np.where(all_ocean_comp == comp)[0]
                    x_line, y_line = all_x[comp_indices], all_y[comp_indices]
                    
                    # Sort points by x for connected lines
                    sorted_idx = np.argsort(x_line)
                    x_line = x_line[sorted_idx]
                    y_line = y_line[sorted_idx]
                    
                    # Set color based on composition
                    if FigMisc.MANUAL_HYDRO_COLORS:
                        thisColor = Color.cmap[comp](Color.GetNormT(np.nanmax(all_z[comp_indices])))
                    else:
                        thisColor = Color.cmap[comp](0.5)
                    
                    # Clean composition label
                    if 'CustomSolution' in comp:
                        comp_label = comp.split('=')[0].replace('CustomSolution', '')
                    else:
                        comp_label = comp
                    
                    # Only add label if this composition hasn't been plotted yet
                    line_label = comp_label if comp_label not in plotted_labels else None
                    if line_label:
                        plotted_labels.add(comp_label)
                    
                    # Plot line for this composition
                    ax.plot(x_line, y_line, color=thisColor, linewidth=FigMisc.ZBD_COMP_LINE_WIDTH, 
                           alpha=FigMisc.ZBD_COMP_LINE_ALPHA, label=line_label, zorder=2)

            # Plot scatter points using combined data
            if np.size(all_x) > 0:
                # Add padding to prevent marker cutoff at edges
                x_range = np.max(all_x) - np.min(all_x)
                y_range = np.max(all_y) - np.min(all_y)
                x_padding = x_range * FigMisc.ZBD_AXIS_PADDING
                y_padding = y_range * FigMisc.ZBD_AXIS_PADDING
                ax.set_xlim([np.min(all_x) - x_padding, np.max(all_x) + x_padding])
                ax.set_ylim([np.min(all_y) - y_padding, np.max(all_y) + y_padding])
            
            # Separate NaN and non-NaN points for different plotting
            nan_mask = np.isnan(all_z)
            valid_mask = ~nan_mask
            
            # Plot non-NaN points colored by z variable
            if np.any(valid_mask):
                z_valid = all_z[valid_mask]
                pts = ax.scatter(all_x[valid_mask], all_y[valid_mask], c=z_valid,
                                cmap=FigMisc.ZBD_COLORMAP, marker=Style.MS_Induction, 
                                s=Style.MW_Induction**2, edgecolors=FigMisc.ZBD_DOT_EDGE_COLOR,
                                linewidths=FigMisc.ZBD_DOT_EDGE_WIDTH, zorder=3, 
                                vmin=np.nanmin(z_valid), vmax=np.nanmax(z_valid))

                # Add colorbar
                cbar = fig.colorbar(pts, ax=ax, format=FigLbl.cbarFmt)
                # Extend colorbar to cover the full range of default ticks
                default_ticks = cbar.get_ticks()
                if len(default_ticks) > 0:
                    pts.set_clim(vmin=np.min(default_ticks), vmax=np.max(default_ticks))
                cbar.set_label(FigLbl.cbarLabelExplore, size=12)
            
            # Plot NaN points with distinct color
            if np.any(nan_mask):
                ax.scatter(all_x[nan_mask], all_y[nan_mask], c=FigMisc.ZBD_NAN_COLOR,
                          marker=FigMisc.ZBD_NAN_MARKER, s=Style.MW_Induction**2, 
                          edgecolors=FigMisc.ZBD_DOT_EDGE_COLOR, linewidths=FigMisc.ZBD_DOT_EDGE_WIDTH,
                          zorder=3, label='Invalid data')
            
            # Add legend for composition lines if enabled
            if FigMisc.DRAW_COMPOSITION_LINE and Params.LEGEND:
                ax.legend(title="Ocean Composition", fontsize=FigMisc.ZBD_COMP_LEGEND_FONT_SIZE, 
                         title_fontsize=FigMisc.ZBD_COMP_LEGEND_TITLE_SIZE)
            elif np.any(nan_mask) and Params.LEGEND:
                # Show legend for invalid data points if no composition lines are shown
                ax.legend(fontsize=FigMisc.ZBD_COMP_LEGEND_FONT_SIZE)

            plt.tight_layout()
            fig.savefig(Params.FigureFiles.exploreZbD, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
            log.debug(f'Combined ZbD plot saved to file: {Params.FigureFiles.exploreZbD}')
            plt.close()

    return


def PlotExploreOgramLoveComparison(ExplorationList, Params):
    """ Plot a scatter showing the evaluated tidal love number k2 versus delta (1+k2-h2),
        for comparison against canonical k2/delta exploration plots.
    """

    for Exploration in (ex for ex in ExplorationList if not ex.NO_H2O):
        Exploration.xName = 'delta_love_number_relation'
        Exploration.yName = 'k_love_number'
        Exploration.zName = 'oceanComp'
        FigLbl.SetExploration(Exploration.bodyname, Exploration.xName, Exploration.yName, Exploration.zName)
        if not FigMisc.TEX_INSTALLED:
            FigLbl.StripLatex()

        fig = plt.figure(figsize=FigSize.explore)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)

        if Params.TITLES:
            fig.suptitle(FigLbl.explorationLoveComparisonTitle)
        ax.set_xlabel(FigLbl.xLabelExplore)
        ax.set_ylabel(FigLbl.yLabelExplore)
        # Override standard settings for this type of plot
        ax.set_xscale('linear')
        ax.set_yscale('linear')

        x = np.reshape(Exploration.__getattribute__(Exploration.xName) * FigLbl.xMultExplore, -1)
        y = np.reshape(Exploration.__getattribute__(Exploration.yName) * FigLbl.yMultExplore, -1)
        z = np.reshape(Exploration.__getattribute__(Exploration.zName), -1)
        # Only keep data points for which a valid model was determined
        VALID = np.logical_not(np.logical_or(np.isnan(x), np.isnan(y)))
        x, y, z = x[VALID], y[VALID], z[VALID]

        # Draw composition lines if enabled
        if FigMisc.DRAW_COMPOSITION_LINE:
            # Group by unique z values and plot connecting lines
            unique_z = np.unique(z)
            for z_val in unique_z:
                TminMax_K = {}
                Tall_K = np.reshape(Exploration.__getattribute__('Tb_K'), -1)
                TminMax_K[z_val] = [np.nanmin(Tall_K), np.nanmax(Tall_K)]
                # Set style options
                if FigMisc.MANUAL_HYDRO_COLORS:
                    Color.Tbounds_K = TminMax_K[z_val]
                    thisColor = Color.cmap[z_val](Color.GetNormT(np.nanmax(Tall_K)))
                else:
                    thisColor = None
                indices = np.where(z == z_val)[0]
                x_line, y_line = x[indices], y[indices]
                # Sort points by x for connected lines
                sorted_idx = np.argsort(x_line)
                x_line = x_line[sorted_idx]
                y_line = y_line[sorted_idx]
                if 'CustomSolution' in z_val:
                    z_label = z_val.split('=')[0].replace('CustomSolution', '')
                else:
                    z_label = z_val
                # Plot line for each z group
                ax.plot(x_line, y_line, color=thisColor, label=f'{z_label}', 
                       linewidth=FigMisc.LOVE_COMP_LINE_WIDTH, alpha=FigMisc.LOVE_COMP_LINE_ALPHA, zorder=2)
                # Add error bars if enabled
                if FigMisc.SHOW_ERROR_BARS:
                    ax.errorbar(x_line, y_line, yerr=FigMisc.ERROR_BAR_MAGNITUDE, fmt='none', color=thisColor, capsize=3)

        # Handle ice thickness coloring or regular scatter
        if FigMisc.SHOW_ICE_THICKNESS_DOTS:
            # Get ice shell thickness and normalize for coloring
            ice_thickness = np.reshape(Exploration.__getattribute__('zb_approximate_km'), -1)
            ice_thickness = ice_thickness[VALID]  # Filter invalid data like x and y
            if np.size(ice_thickness) > 0:
                # Set bounds for normalization (min and max of ice thickness)
                Tbound_lower = np.min(ice_thickness)
                Tbound_upper = np.max(ice_thickness)
                # Get normalized values using GetNormT
                norm_thickness = Color.GetNormT(ice_thickness, Tbound_lower, Tbound_upper)
            else:
                norm_thickness = None
            
            # Apply colormap using normalized thickness
            pts = ax.scatter(x, y, c=norm_thickness,
                            cmap=FigMisc.LOVE_ICE_THICKNESS_CMAP, marker=Style.MS_Induction, 
                            s=Style.MW_Induction**2, edgecolors=FigMisc.LOVE_DOT_EDGE_COLOR, 
                            linewidths=FigMisc.LOVE_DOT_EDGE_WIDTH, zorder=3)
            
            # Create legend showing ice thickness values instead of colorbar
            if np.size(ice_thickness) > 0:
                # Round thickness values to avoid floating-point duplicates
                rounded_thicknesses = np.round(ice_thickness)
                unique_thicknesses = np.unique(rounded_thicknesses)
                
                if len(unique_thicknesses) > FigMisc.LOVE_MAX_LEGEND_ENTRIES:  # Limit number of legend entries
                    # Select representative values
                    indices = np.linspace(0, len(unique_thicknesses)-1, FigMisc.LOVE_MAX_LEGEND_ENTRIES, dtype=int)
                    selected_thicknesses = unique_thicknesses[indices]
                else:
                    selected_thicknesses = unique_thicknesses
                
                # Create legend elements
                legend_elements = []
                for thickness in selected_thicknesses:
                    norm_val = Color.GetNormT(thickness, Tbound_lower, Tbound_upper)
                    color = getattr(plt.cm, FigMisc.LOVE_ICE_THICKNESS_CMAP)(norm_val)
                    legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                                    markerfacecolor=color, markeredgecolor=FigMisc.LOVE_DOT_EDGE_COLOR,
                                                    markersize=8, label=f'{thickness:.0f}'))
                
                # Add legend for ice thickness
                if Params.LEGEND:
                    thickness_legend = ax.legend(handles=legend_elements, title="Ice Shell Thickness (km)", 
                                               loc='upper right', bbox_to_anchor=(1.0, 1.0), 
                                               fontsize=FigMisc.LOVE_ICE_LEGEND_FONT_SIZE, title_fontsize=FigMisc.LOVE_ICE_LEGEND_TITLE_SIZE)
                    ax.add_artist(thickness_legend)  # Keep this legend when adding composition legend
        else:
            # Standard scatter plot without ice thickness coloring
            # Get ice shell thickness for backward compatibility
            ice_thickness = np.reshape(Exploration.__getattribute__('zb_approximate_km'), -1)
            ice_thickness = ice_thickness[VALID]
            if np.size(ice_thickness) > 0:
                # Set bounds for normalization (min and max of ice thickness)
                Tbound_lower = np.min(ice_thickness)
                Tbound_upper = np.max(ice_thickness)
                # Get normalized values using GetNormT
                norm_thickness = Color.GetNormT(ice_thickness, Tbound_lower, Tbound_upper)
            else:
                norm_thickness = None
            
            # Apply colormap using normalized thickness (maintaining backward compatibility)
            pts = ax.scatter(x, y, c=norm_thickness,
                            cmap=FigMisc.LOVE_ICE_THICKNESS_CMAP, marker=Style.MS_Induction, 
                            s=Style.MW_Induction**2, edgecolors=FigMisc.LOVE_DOT_EDGE_COLOR, 
                            linewidths=FigMisc.LOVE_DOT_EDGE_WIDTH, zorder=3)

        # Set axis limits
        if np.size(x) > 0:
            ax.set_xlim([0, np.max(x) + 0.01])
        if np.size(y) > 0:
            # Account for error bars when setting y limits
            if FigMisc.SHOW_ERROR_BARS:
                error_magnitude = FigMisc.ERROR_BAR_MAGNITUDE
                y_min_with_error = np.min(y) - error_magnitude
                y_max_with_error = np.max(y) + error_magnitude
                ax.set_ylim([np.floor(y_min_with_error*100)/100, np.ceil(y_max_with_error*100)/100])
            else:
                ax.set_ylim([np.floor(np.min(y)*100)/100, np.ceil(np.max(y)*100)/100])

        # Add legend for composition lines if enabled
        if FigMisc.DRAW_COMPOSITION_LINE and Params.LEGEND:
            if FigMisc.SHOW_ICE_THICKNESS_DOTS:
                # Position composition legend to avoid overlap with ice thickness legend
                ax.legend(title="Ocean Composition", fontsize=FigMisc.LOVE_COMP_LEGEND_FONT_SIZE, 
                         title_fontsize=FigMisc.LOVE_COMP_LEGEND_TITLE_SIZE,
                         loc='upper left', bbox_to_anchor=(0.0, 1.0))
            else:
                ax.legend(title="Ocean Composition", fontsize=FigMisc.LOVE_COMP_LEGEND_FONT_SIZE, 
                         title_fontsize=FigMisc.LOVE_COMP_LEGEND_TITLE_SIZE)

        plt.tight_layout()
        fig.savefig(Params.FigureFiles.exploreLoveComparison, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Plot saved to file: {Params.FigureFiles.exploreLoveComparison}')
        plt.close()

    # Plot combination
    if Params.COMPARE and np.size(ExplorationList) > 1:
        # Filter out explorations with no H2O
        ValidExplorations = [ex for ex in ExplorationList if not ex.NO_H2O]
        
        if len(ValidExplorations) > 1:
            # Set up exploration parameters for the first valid exploration
            Exploration = ValidExplorations[0]
            Exploration.xName = 'delta_love_number_relation'
            Exploration.yName = 'k_love_number'
            Exploration.zName = 'oceanComp'
            FigLbl.SetExploration(Exploration.bodyname, Exploration.xName, Exploration.yName, Exploration.zName)
            if not FigMisc.TEX_INSTALLED:
                FigLbl.StripLatex()

            fig = plt.figure(figsize=FigSize.explore)
            grid = GridSpec(1, 1)
            ax = fig.add_subplot(grid[0, 0])
            if Style.GRIDS:
                ax.grid()
                ax.set_axisbelow(True)

            if Params.TITLES:
                fig.suptitle(FigLbl.explorationLoveComparisonTitle)
            ax.set_xlabel(FigLbl.xLabelExplore)
            ax.set_ylabel(FigLbl.yLabelExplore)
            # Override standard settings for this type of plot
            ax.set_xscale('linear')
            ax.set_yscale('linear')

            # Collect data from all explorations
            all_x = []
            all_y = []
            all_z = []
            all_ice_thickness = []
            exploration_labels = []
            
            for i, Exploration in enumerate(ValidExplorations):
                x = np.reshape(Exploration.__getattribute__(Exploration.xName) * FigLbl.xMultExplore, -1)
                y = np.reshape(Exploration.__getattribute__(Exploration.yName) * FigLbl.yMultExplore, -1)
                z = np.reshape(Exploration.__getattribute__(Exploration.zName), -1)
                # Only keep data points for which a valid model was determined
                VALID = np.logical_not(np.logical_or(np.isnan(x), np.isnan(y)))
                x, y, z = x[VALID], y[VALID], z[VALID]
                
                # Get ice thickness data
                ice_thickness = np.reshape(Exploration.__getattribute__('zb_approximate_km'), -1)
                ice_thickness = ice_thickness[VALID]
                
                # Add exploration identifier
                exploration_id = np.full(len(x), i)
                
                all_x.extend(x)
                all_y.extend(y)
                all_z.extend(z)
                all_ice_thickness.extend(ice_thickness)
                exploration_labels.extend(exploration_id)

            # Convert to numpy arrays
            all_x = np.array(all_x)
            all_y = np.array(all_y)
            all_z = np.array(all_z)
            all_ice_thickness = np.array(all_ice_thickness)
            exploration_labels = np.array(exploration_labels)

            # Draw composition lines if enabled
            if FigMisc.DRAW_COMPOSITION_LINE:
                # Group by unique z values and plot connecting lines
                unique_z = np.unique(all_z)
                plotted_labels = set()  # Keep track of which labels have been plotted
                
                for z_val in unique_z:
                    z_indices = np.where(all_z == z_val)[0]
                    
                    # Collect all points for this z value across all explorations
                    x_line, y_line = all_x[z_indices], all_y[z_indices]
                    
                    # Sort points by x for connected lines
                    sorted_idx = np.argsort(x_line)
                    x_line = x_line[sorted_idx]
                    y_line = y_line[sorted_idx]
                    
                    # Set style options (use first exploration's temperature bounds)
                    TminMax_K = {}
                    Tall_K = np.reshape(ValidExplorations[0].__getattribute__('Tb_K'), -1)
                    TminMax_K[z_val] = [np.nanmin(Tall_K), np.nanmax(Tall_K)]
                    
                    if FigMisc.MANUAL_HYDRO_COLORS:
                        Color.Tbounds_K = TminMax_K[z_val]
                        thisColor = Color.cmap[z_val](Color.GetNormT(np.nanmax(Tall_K)))
                    else:
                        thisColor = None
                    
                    if 'CustomSolution' in z_val:
                        z_label = z_val.split('=')[0].replace('CustomSolution', '')
                    else:
                        z_label = z_val
                    
                    # Only add label if this composition hasn't been plotted yet
                    line_label = z_label if z_label not in plotted_labels else None
                    if line_label:
                        plotted_labels.add(z_label)
                    
                    # Plot line for each z group combining all explorations
                    ax.plot(x_line, y_line, color=thisColor, label=line_label, 
                           linewidth=FigMisc.LOVE_COMP_LINE_WIDTH, alpha=FigMisc.LOVE_COMP_LINE_ALPHA, zorder=2)
                    
                    # Add error bars if enabled
                    if FigMisc.SHOW_ERROR_BARS:
                        ax.errorbar(x_line, y_line, yerr=FigMisc.ERROR_BAR_MAGNITUDE, fmt='none', color=thisColor, capsize=3)

            # Handle ice thickness coloring or regular scatter
            if FigMisc.SHOW_ICE_THICKNESS_DOTS and len(all_ice_thickness) > 0:
                # Set bounds for normalization (min and max of ice thickness across all explorations)
                Tbound_lower = np.min(all_ice_thickness)
                Tbound_upper = np.max(all_ice_thickness)
                
                # Get normalized values using GetNormT
                norm_thickness = Color.GetNormT(all_ice_thickness, Tbound_lower, Tbound_upper)
            else:
                norm_thickness = None

            # Plot scatter points using the same marker for all explorations
            ax.scatter(all_x, all_y, c=norm_thickness,
                      cmap=FigMisc.LOVE_ICE_THICKNESS_CMAP, marker=Style.MS_Induction, 
                      s=Style.MW_Induction**2, edgecolors=FigMisc.LOVE_DOT_EDGE_COLOR, 
                      linewidths=FigMisc.LOVE_DOT_EDGE_WIDTH, zorder=3)
                
            # Create legend showing ice thickness values instead of colorbar
            if FigMisc.SHOW_ICE_THICKNESS_DOTS and len(all_ice_thickness) > 0 and Params.LEGEND:
                # Round thickness values to avoid floating-point duplicates
                rounded_thicknesses = np.round(all_ice_thickness)
                unique_thicknesses = np.unique(rounded_thicknesses)
                
                if len(unique_thicknesses) > FigMisc.LOVE_MAX_LEGEND_ENTRIES:  # Limit number of legend entries
                    # Select representative values
                    indices = np.linspace(0, len(unique_thicknesses)-1, FigMisc.LOVE_MAX_LEGEND_ENTRIES, dtype=int)
                    selected_thicknesses = unique_thicknesses[indices]
                else:
                    selected_thicknesses = unique_thicknesses
                
                # Create legend elements
                legend_elements = []
                for thickness in selected_thicknesses:
                    norm_val = Color.GetNormT(thickness, Tbound_lower, Tbound_upper)
                    color = getattr(plt.cm, FigMisc.LOVE_ICE_THICKNESS_CMAP)(norm_val)
                    legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                                    markerfacecolor=color, markeredgecolor=FigMisc.LOVE_DOT_EDGE_COLOR,
                                                    markersize=8, label=f'{thickness:.0f}'))
                
                # Add legend for ice thickness
                thickness_legend = ax.legend(handles=legend_elements, title="Ice Shell Thickness (km)", 
                                           loc='upper right', bbox_to_anchor=(1.0, 1.0), 
                                           fontsize=FigMisc.LOVE_ICE_LEGEND_FONT_SIZE, title_fontsize=FigMisc.LOVE_ICE_LEGEND_TITLE_SIZE)
                ax.add_artist(thickness_legend)  # Keep this legend when adding composition legend

            # Set axis limits
            if np.size(all_x) > 0:
                ax.set_xlim([0, np.max(all_x) + 0.01])
            if np.size(all_y) > 0:
                # Account for error bars when setting y limits
                if FigMisc.SHOW_ERROR_BARS:
                    error_magnitude = FigMisc.ERROR_BAR_MAGNITUDE
                    y_min_with_error = np.min(all_y) - error_magnitude
                    y_max_with_error = np.max(all_y) + error_magnitude
                    ax.set_ylim([np.floor(y_min_with_error*100)/100, np.ceil(y_max_with_error*100)/100])
                else:
                    ax.set_ylim([np.floor(np.min(all_y)*100)/100, np.ceil(np.max(all_y)*100)/100])

            # Add legend for composition lines if enabled
            if FigMisc.DRAW_COMPOSITION_LINE and Params.LEGEND:
                if FigMisc.SHOW_ICE_THICKNESS_DOTS:
                    # Position composition legend to avoid overlap with ice thickness legend
                    ax.legend(title="Ocean Composition", fontsize=FigMisc.LOVE_COMP_LEGEND_FONT_SIZE, 
                             title_fontsize=FigMisc.LOVE_COMP_LEGEND_TITLE_SIZE,
                             loc='upper left', bbox_to_anchor=(0.0, 1.0))
                else:
                    ax.legend(title="Ocean Composition", fontsize=FigMisc.LOVE_COMP_LEGEND_FONT_SIZE, 
                             title_fontsize=FigMisc.LOVE_COMP_LEGEND_TITLE_SIZE)

            plt.tight_layout()
            fig.savefig(Params.FigureFiles.exploreLoveComparison, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
            log.debug(f'Combined Love comparison plot saved to file: {Params.FigureFiles.exploreLoveComparison}')
            plt.close()

    return
