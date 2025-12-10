import os
import numpy as np
import re
import logging
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap as DiscreteCmap, to_rgb, BoundaryNorm, ListedColormap, TwoSlopeNorm
from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc
from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS, GetIceEOS
from PlanetProfile.Utilities.Indexing import PhaseConv, PhaseInv, GetPhaseIndices
from PlanetProfile.Thermodynamics.InnerEOS import GetInnerEOS
from PlanetProfile.Utilities.defineStructs import Constants
from typing import Optional, List
import reaktoro as rkt
import itertools
import copy

# Assign logger
log = logging.getLogger('PlanetProfile')



# Unlike most other plotting routines, for those below we can't plot multiple bodies together.
# Just make these plots for the first model, which is the primary.

def PlotHydrosphereSpecies(PlanetList, Params):
    """
    Plot species in hydrosphere.
    """
    # Only make this plot once for a planet
    Planet = PlanetList[0]
    if not Planet.Do.NO_H2O and 'CustomSolution' in Planet.Ocean.comp:
        fig = plt.figure(figsize=FigSize.vhydroSpecies)
        grid = GridSpec(4, 2)
        # Get supcrt database to check phase type
        supcrt_database = rkt.SupcrtDatabase(Params.CustomSolution.SUPCRT_DATABASE)
        speciesPhase = [supcrt_database.species(species).aggregateState() for species in Planet.Ocean.aqueousSpecies]
        if np.any([speciesPhase[i] == rkt.AggregateState.Solid for i in range(len(speciesPhase))]):
            solidAx = True
            solidspeciesax = fig.add_subplot(grid[0:3, 0])
            aqueousSpeciesAx = fig.add_subplot(grid[0:3, 1])
        else:
            solidAx = False
            aqueousSpeciesAx = fig.add_subplot(grid[0:3, 0:2])
        # If we have a reaction with affinities to plot, then we should split the second axis into two columns
        plot_reaction_marker = Planet.Ocean.Reaction.reaction != "NaN"
        if plot_reaction_marker:
            pHax = fig.add_subplot(grid[3, 0])
            affinityax = fig.add_subplot(grid[3, 1])
            affinityax.set_xlabel(FigLbl.zLabel)
            affinityax.set_ylabel(FigLbl.rxnAffinityLabel)
        # If not, then just plot pH
        else:
            pHax = fig.add_subplot(grid[3, :])
        if solidAx:
            axs = [solidspeciesax, aqueousSpeciesAx]
            solidspeciesax.set_xlabel(FigLbl.solidSpeciesLabel)
            solidspeciesax.set_ylabel(FigLbl.zLabel)
        else:
            axs = [aqueousSpeciesAx]
            aqueousSpeciesAx.set_ylabel(FigLbl.zLabel)
        aqueousSpeciesAx.set_xlabel(FigLbl.aqueousSpeciesLabel)
        for ax in axs:
            if Style.GRIDS:
                ax.grid()
                ax.set_axisbelow(True)
        pHax.set_xlabel(FigLbl.zLabel)
        pHax.set_ylabel(FigLbl.pHLabel)

        # Get liquid indices
        indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund, \
            indsClath, indsClathWet, indsMixedClathrateIh, indsMixedClathrateII, indsMixedClathrateIII, indsMixedClathrateV, indsMixedClathrateVI, \
            indsMixedClathrateIhwet, indsMixedClathrateIIund, indsMixedClathrateIIIund, indsMixedClathrateVund, indsMixedClathrateVIund, \
            indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI, \
            indsFe = GetPhaseIndices(Planet.phase)
        # If we have liquid indices then let's plot hydrosphere species
        if np.size(indsLiq) != 0:
            aqueousSpecies = Planet.Ocean.aqueousSpecies.copy()
            allSpeciesData = Planet.Ocean.aqueousSpeciesAmount_mol.copy()
            
            # Set overall figure title
            if Params.TITLES:
                fig.suptitle(
                    f'{PlanetList[0].bodyname} {PlanetList[0].label},{FigLbl.hydroSpeciesTitle}')
            ocean_depth = Planet.z_m[indsLiq]
            
            # Go through each species and plot
            text_objects = []  # Store text objects for adjust_text if available
                
            
            for i, species in enumerate(aqueousSpecies):
                speciesAmountData = allSpeciesData[i]
                species_phase = ''
                if supcrt_database.species(species).aggregateState() == rkt.AggregateState.Aqueous:
                    species_phase = 'aqueous'
                elif supcrt_database.species(species).aggregateState() == rkt.AggregateState.Gas:
                    species_phase = 'gas'
                else:
                    species_phase = 'solid'
                if FigMisc.TEX_INSTALLED:
                    species_name = re.sub(r'(\w)(\+|\-)(\d+)', r'\1^{\3\2}', species)
                    species_label = rf"$\ce{{{species_name}}}$"
                else:
                    species_label = species
                # Plot species labels - But only if they are above FigMisc.minVolSolidThreshold_cm3
                if species_phase == 'solid':
                    indices_above_min_threshold = np.where(speciesAmountData > FigMisc.minVolSolidThreshold_cm3)[0]
                else:
                    # Plot species labels - But only if they are above FigMisc.minThreshold
                    indices_above_min_threshold = np.where(speciesAmountData > FigMisc.minAqueousThreshold)[0]
                
                if indices_above_min_threshold.size > 0 and species not in FigMisc.excludeSpeciesFromHydrospherePlot:
                    indexPosition = indices_above_min_threshold[(i*4) % len(indices_above_min_threshold)]
                    x_label_pos = speciesAmountData[indexPosition]
                    y_label_pos = ocean_depth[indexPosition] / 1e3
                    style = Style.LS_hydroSpecies[species_phase]
                    linewidth = Style.LW_hydroSpecies[species_phase]
                    cmap = Color.cmap['hydroSpecies']
                    if isinstance(cmap, ListedColormap):
                        nColors = len(cmap.colors)
                        color = Color.cmap['hydroSpecies'](i % nColors)
                    else:
                        color = Color.cmap['hydroSpecies'](i % len(aqueousSpecies)) / len(aqueousSpecies)
                    if species_phase == 'solid':
                        line, = solidspeciesax.plot(speciesAmountData, ocean_depth / 1e3, linestyle=style, color=color, linewidth = linewidth)
                        # Add text to all species axis
                        text_obj = solidspeciesax.text(x_label_pos, y_label_pos, species_label,
                                        color=line.get_color(),
                                        verticalalignment='bottom',
                                        horizontalalignment='center',
                                        fontsize=FigLbl.speciesSize, bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.5, edgecolor='none'), zorder=10)
                        text_objects.append(text_obj)
                    elif species_phase == 'aqueous':
                            line, = aqueousSpeciesAx.plot(speciesAmountData, ocean_depth / 1e3, linestyle=style, color=color, linewidth = linewidth)
                            if x_label_pos >= 0:
                                aqueousSpeciesAx.text(x_label_pos, y_label_pos, species_label,
                                                    color=line.get_color(),
                                                    verticalalignment='bottom',
                                                    horizontalalignment='center',
                                                    fontsize=FigLbl.speciesSize,
                                                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.5, edgecolor='none'), zorder=10)
            """if Planet.Ocean.Reaction.reaction != 'NaN':
                for species in Planet.Ocean.Reaction.disequilibriumConcentrations.keys():
                    if not Planet.Ocean.Reaction.useReferenceSpecies:
                        speciesDisequilibriumAmount = Planet.Ocean.Reaction.disequilibriumConcentrations[species]
                        speciesIndex = np.where(Planet.Ocean.aqueousSpecies == species)[0][0]
                        speciesEquilibriumAmount = Planet.Ocean.aqueousSpeciesAmount_mol[speciesIndex]
                        if speciesAmountData is None:
                            speciesAmountData = speciesEquilibriumAmount
                        else:
                            speciesAmountData = np.repeat(speciesDisequilibriumAmount, len(speciesAmountData))
                    else:
                        referenceSpecies = Planet.Ocean.Reaction.referenceSpecies
                        referenceSpeciesIndex = np.where(Planet.Ocean.aqueousSpecies == referenceSpecies)[0][0]
                        referenceSpeciesAmount = Planet.Ocean.aqueousSpeciesAmount_mol[referenceSpeciesIndex]
                        speciesRatioToReferenceSpecies = Planet.Ocean.Reaction.disequilibriumConcentrations[species]
                        if speciesRatioToReferenceSpecies is None:
                            speciesIndex = np.where(Planet.Ocean.aqueousSpecies == species)[0][0]
                            speciesAmountData = Planet.Ocean.aqueousSpeciesAmount_mol[speciesIndex]
                        else:   
                            speciesAmountData = speciesRatioToReferenceSpecies * referenceSpeciesAmount
                    style = Style.LS_hydroSpeciesDisequilibrium
                    linewidth = Style.LW_hydroSpeciesDisequilibrium
                    cmap = Color.cmap['hydroSpecies']
                    line = aqueousSpeciesAx.plot(speciesAmountData, ocean_depth / 1e3, linestyle=style, color=color, linewidth = linewidth)
                    """
                    
                
            for ax in axs:
                ax.set_xscale('log')
                current_xlim = ax.get_xlim()
                if solidAx and ax == solidspeciesax: 
                    new_xmin = max(current_xlim[0], FigMisc.minVolSolidThreshold_cm3)
                else:
                    new_xmin = max(current_xlim[0], FigMisc.minAqueousThreshold)
                new_xmax = 10 ** np.ceil(np.log10(current_xlim[1]))
                ax.set_xlim([new_xmin, new_xmax])
                ax.invert_yaxis()
            # Plot pH plot
            bulk_pH_not_nan = np.where(~np.isnan(Planet.Ocean.Bulk_pHs))[0]
            bulk_line, = pHax.plot(ocean_depth[bulk_pH_not_nan] / 1e3, Planet.Ocean.Bulk_pHs[bulk_pH_not_nan], linestyle ='-',
                              color = 'black', label = 'Bulk Ocean pH')
            # If we should plot reaction, then let's plot reaction pHs and affinity
            if plot_reaction_marker:
                if FigMisc.TEX_INSTALLED:
                    reaction_label = rf"$\ce{{{Planet.Ocean.Reaction.reaction}}}$"
                else:
                    reaction_label = Planet.Ocean.Reaction.reaction
                affinityax.plot(ocean_depth[bulk_pH_not_nan] / 1e3, Planet.Ocean.affinity_kJ[bulk_pH_not_nan], linestyle ='-',
                              color = 'black', label = reaction_label)
                if Params.LEGEND:
                    handles, lbls = affinityax.get_legend_handles_labels()
                    affinityax.legend(handles, lbls, fontsize = 5)

            plt.tight_layout()
            fig.savefig(Params.FigureFiles.hydroSpecies, format=FigMisc.figFormat, dpi=FigMisc.dpi,
                        metadata=FigLbl.meta)
            log.debug(f'Ocean aqueous species plot saved to file: {Params.FigureFiles.hydroSpecies}')
            plt.close()
        else:
            log.warning("There is no ocean, thus will not plot a hydrosphere species plot for this model.")


def PlotHydroPhase(PlanetList, Params):
    if os.path.dirname(Params.FigureFiles.vpvtHydro) != 'Comparison':
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
            Planet = PlanetList[0]
            P_MPa = Planet.P_MPa[:Planet.Steps.nHydro]
            P_MPa = P_MPa[np.logical_and(P_MPa >= Pmin_MPa, P_MPa <= Pmax_MPa)]
            for Planet in PlanetList:
                # Get highest fidelity pressure range
                P = Planet.P_MPa[:Planet.Steps.nHydro]
                P = P[np.logical_and(P >= Pmin_MPa, P <= Pmax_MPa)]
                if P.size > P_MPa.size:
                    P_MPa = P
        else:
            P_MPa = np.linspace(Pmin_MPa, Pmax_MPa, FigMisc.nPhydro)
        if FigMisc.nThydro is None:
            Planet = PlanetList[0]
            T_K = Planet.T_K[:Planet.Steps.nHydro]
            T_K = T_K[np.logical_and(T_K >= Tmin_K, T_K <= Tmax_K)]
            for Planet in PlanetList:
                # Get highest fidelity pressure range
                T = Planet.T_K[:Planet.Steps.nHydro]
                T = T[np.logical_and(T >= Tmin_K, T <= Tmax_K)]
                if T.size > T_K.size:
                    T_K = T
        else:
            T_K = np.linspace(Tmin_K, Tmax_K, FigMisc.nThydro)
        oceanListEOS = []
        phasesList = []
        ices = {}
        iceEOS = {}
        for Planet in PlanetList:
            # Load EOS independently from model run, because we will query wider ranges of conditions
            oceanEOS = GetOceanEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, P_MPa, T_K,
                                   Planet.Ocean.MgSO4elecType, rhoType=Planet.Ocean.MgSO4rhoType,
                                   scalingType=Planet.Ocean.MgSO4scalingType, FORCE_NEW=Params.FORCE_EOS_RECALC,
                                   phaseType=Planet.Ocean.phaseType, EXTRAP=Params.EXTRAP_OCEAN,
                                   sigmaFixed_Sm=Planet.Ocean.sigmaFixed_Sm, LOOKUP_HIRES=Planet.Do.OCEAN_PHASE_HIRES, kThermConst_WmK=Planet.Ocean.kThermWater_WmK,
                                   propsStepReductionFactor=Planet.Ocean.propsStepReductionFactor)

            phases = oceanEOS.fn_phase(P_MPa, T_K, grid=True).astype(int)
            new_ices = set([PhaseConv(ice) for ice in np.unique(np.append(phases[phases != 0], 1))])
            if not new_ices.issubset(ices):
                ices = new_ices
                iceEOS = {PhaseInv(ice): GetIceEOS(P_MPa, T_K, ice,
                                                   porosType=Planet.Ocean.porosType[ice],
                                                   phiTop_frac=Planet.Ocean.phiMax_frac[ice],
                                                   Pclosure_MPa=Planet.Ocean.Pclosure_MPa[ice],
                                                   phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE[ice],
                                                   ICEIh_DIFFERENT=Planet.Do.ICEIh_DIFFERENT, kThermConst_WmK=Planet.Ocean.kThermIce_WmK)
                          for ice in ices}
            # Add clathrates to phase and property diagrams where it is stable (if modeled)
            if Planet.Do.CLATHRATE:
                clath = PhaseConv(Constants.phaseClath)
                if Planet.Do.MIXED_CLATHRATE_ICE:
                    phaseIndex = Constants.phaseClath + 1
                    phaseStr = PhaseConv(phaseIndex)
                else:
                    phaseIndex = Constants.phaseClath
                    phaseStr = PhaseConv(phaseIndex)
                if phaseStr not in ices:
                    ices.add(phaseStr)
                    iceEOS[phaseStr] = GetIceEOS(P_MPa, T_K, phaseStr,
                                                             porosType=Planet.Ocean.porosType[clath],
                                                             phiTop_frac=Planet.Ocean.phiMax_frac[clath],
                                                             Pclosure_MPa=Planet.Ocean.Pclosure_MPa[clath],
                                                             phiMin_frac=Planet.Ocean.phiMin_frac,
                                                             EXTRAP=Params.EXTRAP_ICE[phaseStr], kThermConst_WmK=Planet.Ocean.kThermIce_WmK,
                                                             mixParameters={'mixFrac': Planet.Bulk.volumeFractionClathrate, 'JmixedRheologyConstant': Planet.Bulk.JmixedRheologyConstant})
                clathStable = iceEOS[phaseIndex].fn_phase(P_MPa, T_K, grid=True)
                phases[clathStable == phaseIndex] = phaseIndex
            oceanListEOS.append(oceanEOS)
            phasesList.append(phases)

        # Go through every two combination
        SinglePlanetPlot = (len(PlanetList) == 1)
        # If we are plotting for only one run, then compare it to a blank HydroEOS
        if SinglePlanetPlot:
            combinations = [(0, 0)]
        else:
            combinations = itertools.combinations(range(0, PlanetList.size), 2)
        for comparison in combinations:
            FirstPlanetIndex = comparison[0]
            SecondPlanetIndex = comparison[1]
            FirstPlanet = PlanetList[FirstPlanetIndex]
            SecondPlanet = PlanetList[SecondPlanetIndex]
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

                if SinglePlanetPlot:
                    fig.suptitle(f'{FirstPlanet.compStr}{FigLbl.hydroPhaseTitle}')
                else:
                    fig.suptitle(
                        f'Comparison of {FirstPlanet.compStr} and {SecondPlanet.compStr}{FigLbl.hydroPhaseTitle}')

            # Plot as colormesh
            if SinglePlanetPlot:
                phases = phasesList[0]
                phaseColors = DiscreteCmap.from_list('icePhases',
                                                     [Color.oceanCmap(1),
                                                      Color.iceIcond,
                                                      Color.iceII,
                                                      Color.iceIII,
                                                      Color.iceV,
                                                      Color.iceVI,
                                                      Color.clathCond], N=7)
                phaseBounds = np.array([0, 0.5, 1.5, 2.5, 4, 5.5, (Constants.phaseClath + 6) / 2, Constants.phaseClath])
                bNorm = BoundaryNorm(boundaries=phaseBounds, ncolors=7)
            else:
                # If we are doing a comparison, then we must find difference between ices
                phases = phasesList[0] - phasesList[1]
                phaseColors = DiscreteCmap.from_list('icePhases',
                                                     [Color.diffPhase,
                                                      Color.correctPhase, Color.diffPhase], N=3)
                phaseBounds = np.array([-1.5, -0.5, 0.5, 1.5])
                bNorm = BoundaryNorm(boundaries=phaseBounds, ncolors=3)

            c = ax.pcolormesh(T_K, P_MPa * FigLbl.PmultHydro, phases, norm=bNorm, cmap=phaseColors,
                          rasterized=FigMisc.PT_RASTER)

            ices = np.unique(phases[phases != 0])
            if SinglePlanetPlot:
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
                    if FigMisc.SHOW_GEOTHERM:
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

            else:
                if Params.LEGEND:
                    colorbar = plt.colorbar(c, ax=ax, boundaries=phaseBounds, ticks=[-1, 0,  1])
                    colorbar.set_ticklabels([f'Ice vs Water Difference', 'Same Phase Prediction', 'Water vs Ice Difference'])
            ax.set_xlim([Tmin_K, Tmax_K])
            ax.set_ylim([Pmin_MPa, Pmax_MPa])
            ax.invert_yaxis()
            plt.tight_layout()
            fig.savefig(Params.FigureFiles.vphase, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
            log.debug(f'Hydrosphere phase diagram saved to file: {Params.FigureFiles.vphase}')
            plt.close()
class HydrosphereProp:
    def __init__(self,
                 prop: str,
                 fn_name: str,
                 fig_location: tuple,
                 ax: Axes,
                 prop_label: str,
                 ocean_prop_index: Optional[int] = None,
                 ice_prop_index: Optional[int] = None):
        self.prop = prop
        self.fn_name = fn_name
        self.fig_location = fig_location
        self.ax = ax
        self.prop_label = prop_label
        self.ocean_prop_index = ocean_prop_index
        self.ice_prop_index = ice_prop_index
        # Initialize prop_data as an empty list
        self.prop_data = []


def PlotIsoThermalPvThydro(PlanetList, Params):
    props = {'rho': HydrosphereProp('rho', 'fn_rho_kgm3', None, None, FigLbl.rhoLabel),
             'Cp': HydrosphereProp('Cp', 'fn_Cp_JkgK', None, None, FigLbl.CpLabel),
             'alpha': HydrosphereProp('alpha', 'fn_alpha_pK', None, None, FigLbl.alphaLabel),
             'VP': HydrosphereProp('VP', 'fn_Seismic', None, None, FigLbl.VPlabel, 0, 0),
             'KS': HydrosphereProp('KS', 'fn_Seismic', None, None, FigLbl.KSlabel, 0, 2),
             'sig': HydrosphereProp('sig', 'fn_sigma_Sm', None, None, FigLbl.sigLabel, ),
             'VS': HydrosphereProp('VS', 'fn_Seismic', None, None, FigLbl.VSlabel, None, 1),
             'GS': HydrosphereProp('GS', 'fn_Seismic', None, None, FigLbl.GSlabel, None, 3)}
    # This will get the actual HydrosphereProp instances, not just their names
    props_to_plot = [props[prop] for prop in FigMisc.propsToPlot if prop in props]
    if len(props_to_plot) == 0:
        log.warning("No valid props provided to plot. Check that FigMisc.propsToPlot has valid names. Will plot all available properties.")
        props_to_plot = [props[prop] for prop in props]
    elif len(props_to_plot) <= 4:
        rows = 1
        cols = len(props_to_plot)
    else:
        rows = 2
        cols = (len(props_to_plot) + 1) // 2
    if FigMisc.PminHydro_MPa is None:
        Pmin_MPa = np.min([Planet.P_MPa[0] for Planet in PlanetList])
    else:
        Pmin_MPa = FigMisc.PminHydro_MPa
    if FigMisc.PmaxHydro_MPa is None:
        Pmax_MPa = np.max([Planet.P_MPa[Planet.Steps.nHydro - 1] for Planet in PlanetList])
    else:
        Pmax_MPa = FigMisc.PmaxHydro_MPa
    # Get highest Tmin
    if not np.any([Planet.Do.NO_OCEAN for Planet in PlanetList]):
        Tmin_K = np.max([np.min(Planet.T_K[:Planet.Steps.nHydro]) for Planet in PlanetList])
    else:
        Tmin_K = np.max([np.min(Planet.T_K) for Planet in PlanetList])
    # Get lowest Tmax
    if not np.any([Planet.Do.NO_OCEAN for Planet in PlanetList]):
        Tmax_K = np.min([Planet.Sil.THydroMax_K for Planet in PlanetList])
    else:
        Tmax_K = np.min([np.max(Planet.T_K[:Planet.Steps.nHydro]) for Planet in PlanetList if not Planet.Do.NO_OCEAN])

    if FigMisc.nPhydro is None:
        Planet = PlanetList[0]
        P_MPa = Planet.P_MPa[:Planet.Steps.nHydro]
        P_MPa = P_MPa[np.logical_and(P_MPa >= Pmin_MPa, P_MPa <= Pmax_MPa)]
        for Planet in PlanetList:
            # Get highest fidelity pressure range
            P = Planet.P_MPa[:Planet.Steps.nHydro]
            P = P[np.logical_and(P >= Pmin_MPa, P <= Pmax_MPa)]
            if P.size > P_MPa.size:
                P_MPa = P
    else:
        P_MPa = np.linspace(Pmin_MPa, Pmax_MPa, FigMisc.nPhydro)

    # Get temperatures to plot that are within bounds
    T_K = np.array([T for T in FigMisc.TtoPlot_K if Tmin_K <= T <= Tmax_K])
    # If no temperatures in T_K, then should just use Tmax_K and plot a single isotherm
    if T_K.size == 0:
        T_K = np.array([Tmax_K])


    oceanListEOS = []
    phasesList = []
    ices = {}
    iceEOS = {}
    for Planet in PlanetList:
        # Load EOS independently from model run, because we will query wider ranges of conditions
        oceanEOS = GetOceanEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, P_MPa, T_K, Planet.Ocean.MgSO4elecType,
                               rhoType=Planet.Ocean.MgSO4rhoType, scalingType=Planet.Ocean.MgSO4scalingType,
                               FORCE_NEW=Params.FORCE_EOS_RECALC, phaseType=Planet.Ocean.phaseType,
                               EXTRAP=Params.EXTRAP_OCEAN, sigmaFixed_Sm=Planet.Ocean.sigmaFixed_Sm,
                               LOOKUP_HIRES=Planet.Do.OCEAN_PHASE_HIRES, kThermConst_WmK=Planet.Ocean.kThermWater_WmK,
                               propsStepReductionFactor=Planet.Ocean.propsStepReductionFactor)

        phases = oceanEOS.fn_phase(P_MPa, T_K, grid=True).astype(int)
        new_ices = set([PhaseConv(ice) for ice in np.unique(np.append(phases[phases != 0], 1))])
        if not new_ices.issubset(ices):
            ices = new_ices
            iceEOS = {PhaseInv(ice): GetIceEOS(P_MPa, T_K, ice, porosType=Planet.Ocean.porosType[ice],
                                               phiTop_frac=Planet.Ocean.phiMax_frac[ice],
                                               Pclosure_MPa=Planet.Ocean.Pclosure_MPa[ice],
                                               phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE[ice],
                                               ICEIh_DIFFERENT=Planet.Do.ICEIh_DIFFERENT,
                                               mixParameters={'mixFrac': Planet.Bulk.volumeFractionClathrate, 'JmixedRheologyConstant': Planet.Bulk.JmixedRheologyConstant}) for ice in ices}
        # Add clathrates to phase and property diagrams where it is stable (if modeled)
        if Planet.Do.CLATHRATE:
            clath = PhaseConv(Constants.phaseClath)
            if Planet.Do.MIXED_CLATHRATE_ICE:
                phaseIndex = Constants.phaseClath + 1
                phaseStr = PhaseConv(phaseIndex)
            else:
                phaseIndex = Constants.phaseClath
                phaseStr = PhaseConv(phaseIndex)
            if phaseStr not in ices:
                ices.add(phaseStr)
                iceEOS[phaseStr] = GetIceEOS(P_MPa, T_K, phaseStr, porosType=Planet.Ocean.porosType[clath],
                                                         phiTop_frac=Planet.Ocean.phiMax_frac[clath],
                                                         Pclosure_MPa=Planet.Ocean.Pclosure_MPa[clath],
                                                         phiMin_frac=Planet.Ocean.phiMin_frac,
                                                         EXTRAP=Params.EXTRAP_ICE[phaseStr],
                                                         mixParameters={'mixFrac': Planet.Bulk.volumeFractionClathrate, 'JmixedRheologyConstant': Planet.Bulk.JmixedRheologyConstant})
            clathStable = iceEOS[phaseStr].fn_phase(P_MPa, T_K, grid=True)
            phases[clathStable == phaseIndex] = phaseIndex
        oceanListEOS.append(oceanEOS)
        phasesList.append(phases)
    # Obtain data for each planet
    for i, Planet in enumerate(PlanetList):
        oceanEOS = oceanListEOS[i]
        phases = phasesList[i]
        # Loop over each property in props_to_plot
        for prop in props_to_plot:
            if prop.prop in ['VS', 'GS']:
                # For VS and GS, these properties only apply to solid ice phases
                prop_data = np.empty((len(P_MPa), len(T_K))) * np.nan
            else:
                # Call the function dynamically with the required arguments
                prop_data = getattr(oceanEOS, prop.fn_name)(P_MPa, T_K, grid=True)
                # In case of VP or KS, we need to get the appropriate data since fn_seismic returns both VP and KS
                if prop.prop in ['VP', 'KS']:
                    prop_data = prop_data[prop.ocean_prop_index]
                # Exclude obviously erroneous Cp values, which happen when extending beyond the knots
                # for the input EOS. This is mainly a problem with GSW and Seawater compositions.
                if prop.prop == 'Cp':
                    prop_data[np.logical_or(prop_data < 3200, prop_data > 5200)] = np.nan
            # Now get data for all ice EOSs and replace in grid
            for iceStr in ices:
                ice = PhaseInv(iceStr)
                if np.any(phases == ice):
                    if prop.prop != 'sig':
                        ice_prop_data = getattr(iceEOS[ice], prop.fn_name)(P_MPa, T_K, grid=True)
                        if prop.prop in ['VP', 'KS', 'VS', 'GS']:
                            prop_data[np.where(phases == ice)] = ice_prop_data[prop.ice_prop_index][
                                np.where(phases == ice)]
                        else:
                            prop_data[np.where(phases == ice)] = ice_prop_data[np.where(phases == ice)]
                    else:
                        # In the case of sig, we need to set a different way
                        prop_data[np.where(phases == ice)] = Planet.Ocean.sigmaIce_Sm[iceStr]
            # Add the data for this planet to the prop.prop_data list
            prop.prop_data.append(prop_data)
    # Go through every two combination and plot
    SinglePlanetPlot = (len(PlanetList) == 1)
    # If we are plotting for only one run, then do not compare it to anything
    if SinglePlanetPlot:
        combinations = [(0, 0)]
    else:
        combinations = itertools.combinations(range(0, PlanetList.size), 2)

    for comparison in combinations:
        FirstPlanetIndex = comparison[0]
        SecondPlanetIndex = comparison[1]
        FirstPlanet = PlanetList[FirstPlanetIndex]
        SecondPlanet = PlanetList[SecondPlanetIndex]
        # Create the figure grid
        fig = plt.figure(figsize=FigSize.vpvt)
        grid = GridSpec(rows, cols)
        axes = np.array([[fig.add_subplot(grid[i, j]) for j in range(cols)] for i in range(rows)])
        axf = axes.flatten()
        # Assign fig_location dynamically to each property and update ax
        for idx, prop in enumerate(props_to_plot):
            row, col = divmod(idx, cols)
            prop.fig_location = (row, col)
            ax = axes[row, col]
            prop.ax = ax  # Update ax in HydrosphereProp
        # Hide extra subplots if there are fewer properties than subplots
        for idx in range(len(props_to_plot), len(axf)):
            axf[idx].set_visible(False)
        if Style.GRIDS:
            [ax.grid() for ax in axf]
            [ax.set_axisbelow(False) for ax in axf]
        # Set overall figure title
        if Params.TITLES:
            if SinglePlanetPlot:
                fig.suptitle(f'{FirstPlanet.compStr}{FigLbl.isoThermPvTtitleHydro}')
            else:
                fig.suptitle(
                    f'Weight Difference Comparison of {FirstPlanet.compStr} and {SecondPlanet.compStr}{FigLbl.isoThermPvTtitleHydro}')  # Now plot the data dynamically on the correct ax for each property
        for idx, prop in enumerate(props_to_plot):
            ax = prop.ax  # Get the axis for this property
            if SinglePlanetPlot:
                prop_data_to_plot = prop.prop_data[FirstPlanetIndex]
            else:
                firstPlanetData = prop.prop_data[FirstPlanetIndex]
                secondPlanetData = prop.prop_data[SecondPlanetIndex]
                prop_data_to_plot = np.abs(firstPlanetData - secondPlanetData) / (
                            (firstPlanetData + secondPlanetData) / 2) * 100
                # Handle phase differences
                phase_difference = (phasesList[FirstPlanetIndex] != phasesList[SecondPlanetIndex])
                phase_different_indices = np.where(phase_difference)
                prop_data_to_plot[phase_different_indices] = np.nan
            # Plot scatter points and isothermal lines
            for i, T in enumerate(T_K):
                # Set style options
                thisColor = Color.PvThydroCmap(Color.GetNormT(T, T_K[0], T_K[-1]))
                thisMarker = Style.MS_isoThermPvTHydro[i % len(Style.MS_isoThermPvTHydro)]

                # Plot the isothermal line
                ax.plot(P_MPa * FigLbl.PmultHydro, prop_data_to_plot[:, i], color=thisColor, linestyle='-')

                # Add a scatter point at a representative pressure on the line
                # Choose a index for visibility and add index to offset for each temperature (modulus by len(P_MPa) to ensure we are within bounds of P_MPa
                scatter_idx = (len(P_MPa) // 4 + i) % len(P_MPa)
                ax.scatter(P_MPa[scatter_idx] * FigLbl.PmultHydro, prop_data_to_plot[scatter_idx, i], color=thisColor,
                           marker=thisMarker, label=f'{T} K', s=100)
            # Set labels, title, etc.
            ax.set_xlabel(FigLbl.PlabelHydro)
            if SinglePlanetPlot:
                ax.set_ylabel(prop.prop_label)
            else:
                ax.set_ylabel('Weight Difference')
            ax.set_xlim([Pmin_MPa, Pmax_MPa])
            # In case of sig or alpha plots, we should increment ymax a bit to not cut off top line
            if prop.prop in ['sig', 'alpha']:
                ymin = np.nanmin(prop_data_to_plot) * 0.95
                ymax = np.nanmax(prop_data_to_plot) * 1.05
            else:
                ymin = np.nanmin(prop_data_to_plot)
                ymax = np.nanmax(prop_data_to_plot)
            if np.isnan(ymin) and np.isnan(ymin):
                ax.set_title(prop.prop_label)
                ax.text(0.5, 0.5, "No valid values\nfor this property.",
                        transform=ax.transAxes, ha='center', va='center')
            else:
                ax.set_ylim([ymin, ymax])
                ax.set_title(prop.prop_label)
        if Params.LEGEND:
            # Add a legend to the plot
            axf[0].legend(title="Temperature (K)")
        plt.tight_layout()
        if SinglePlanetPlot:
            saveFile = Params.FigureFiles.isoThermvpvtHydro
        else:
            saveFile = Params.FigureFiles.comparisonFileGenerator(FirstPlanet.saveLabel, SecondPlanet.saveLabel,
                                                                  'isoThermvpvtHydro')
        fig.savefig(saveFile, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'IsoTherm Hydrosphere PT properties plot saved to file: {saveFile}')
        plt.close()

    return


def PlotPvThydro(PlanetList, Params):
    props = {'rho': HydrosphereProp('rho', 'fn_rho_kgm3', None, None, FigLbl.rhoLabel),
            'Cp': HydrosphereProp('Cp', 'fn_Cp_JkgK', None, None, FigLbl.CpLabel),
            'alpha': HydrosphereProp('alpha', 'fn_alpha_pK', None, None, FigLbl.alphaLabel),
            'VP': HydrosphereProp('VP', 'fn_Seismic', None, None, FigLbl.VPlabel, 0, 0),
             'KS': HydrosphereProp('KS', 'fn_Seismic', None, None, FigLbl.KSlabel, 0, 2),
            'sig': HydrosphereProp('sig', 'fn_sigma_Sm', None, None, FigLbl.sigLabel,),
            'VS': HydrosphereProp('VS', 'fn_Seismic', None, None, FigLbl.VSlabel, None, 1),
             'GS': HydrosphereProp('GS', 'fn_Seismic', None, None, FigLbl.GSlabel, None, 3)}
    # This will get the actual HydrosphereProp instances, not just their names
    props_to_plot = [props[prop] for prop in props if prop in FigMisc.propsToPlot]
    if len(props_to_plot) == 0:
        log.warning("No valid props provided to plot. Check that FigMisc.propsToPlot has valid names. Will plot all available properties.")
        props_to_plot = [props[prop] for prop in props]
    elif len(props_to_plot) <= 4:
        rows = 1
        cols = len(props_to_plot)
    else:
        rows = 2
        cols = (len(props_to_plot) + 1) // 2
    if FigMisc.PminHydro_MPa is None:
        Pmin_MPa = np.min([Planet.P_MPa[0] for Planet in PlanetList])
    else:
        Pmin_MPa = FigMisc.PminHydro_MPa
    if FigMisc.PmaxHydro_MPa is None:
        Pmax_MPa = np.max([Planet.P_MPa[Planet.Steps.nHydro - 1] for Planet in PlanetList])
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
            Tmax_K = np.max(
                [np.max(Planet.T_K[:Planet.Steps.nHydro]) for Planet in PlanetList if not Planet.Do.NO_OCEAN])
    else:
        Tmax_K = FigMisc.TmaxHydro_K

    if FigMisc.nPhydro is None:
        Planet = PlanetList[0]
        P_MPa = Planet.P_MPa[:Planet.Steps.nHydro]
        P_MPa = P_MPa[np.logical_and(P_MPa >= Pmin_MPa, P_MPa <= Pmax_MPa)]
        for Planet in PlanetList:
            # Get highest fidelity pressure range
            P = Planet.P_MPa[:Planet.Steps.nHydro]
            P = P[np.logical_and(P >= Pmin_MPa, P <= Pmax_MPa)]
            if P.size > P_MPa.size:
                P_MPa = P
    else:
        P_MPa = np.linspace(Pmin_MPa, Pmax_MPa, FigMisc.nPhydro)
    if FigMisc.nThydro is None:
        Planet = PlanetList[0]
        T_K = Planet.T_K[:Planet.Steps.nHydro]
        T_K = T_K[np.logical_and(T_K >= Tmin_K, T_K <= Tmax_K)]
        for Planet in PlanetList:
            # Get highest fidelity temperature range
            T = Planet.T_K[:Planet.Steps.nHydro]
            T = T[np.logical_and(T >= Tmin_K, T <= Tmax_K)]
            if T.size > T_K.size:
                T_K = T
    else:
        T_K = np.linspace(Tmin_K, Tmax_K, FigMisc.nThydro)
    oceanListEOS = []
    phasesList = []
    ices = {}
    iceEOS = {}
    for Planet in PlanetList:
        # Load EOS independently from model run, because we will query wider ranges of conditions
        oceanEOS = GetOceanEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, P_MPa, T_K, Planet.Ocean.MgSO4elecType,
                               rhoType=Planet.Ocean.MgSO4rhoType, scalingType=Planet.Ocean.MgSO4scalingType,
                               FORCE_NEW=Params.FORCE_EOS_RECALC, phaseType=Planet.Ocean.phaseType,
                               EXTRAP=Params.EXTRAP_OCEAN, sigmaFixed_Sm=Planet.Ocean.sigmaFixed_Sm,
                               LOOKUP_HIRES=Planet.Do.OCEAN_PHASE_HIRES, kThermConst_WmK=Planet.Ocean.kThermWater_WmK,
                               propsStepReductionFactor=Planet.Ocean.propsStepReductionFactor)

        phases = oceanEOS.fn_phase(P_MPa, T_K, grid=True).astype(int)
        new_ices = set([PhaseConv(ice) for ice in np.unique(np.append(phases[phases != 0], 1))])
        if not new_ices.issubset(ices):
            ices = new_ices
            iceEOS = {PhaseInv(ice): GetIceEOS(P_MPa, T_K, ice, porosType=Planet.Ocean.porosType[ice],
                                               phiTop_frac=Planet.Ocean.phiMax_frac[ice],
                                               Pclosure_MPa=Planet.Ocean.Pclosure_MPa[ice],
                                               phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE[ice],
                                               kThermConst_WmK=Planet.Ocean.kThermIce_WmK,
                                               ICEIh_DIFFERENT=Planet.Do.ICEIh_DIFFERENT, mixParameters={'mixFrac': Planet.Bulk.volumeFractionClathrate, 'JmixedRheologyConstant': Planet.Bulk.JmixedRheologyConstant})
                      for ice in ices}
        # Add clathrates to phase and property diagrams where it is stable (if modeled)
        if Planet.Do.CLATHRATE:
            clath = PhaseConv(Constants.phaseClath)
            if Planet.Do.MIXED_CLATHRATE_ICE:
                phaseIndex = Constants.phaseClath + 1
                phaseStr = PhaseConv(phaseIndex)
            else:
                phaseIndex = Constants.phaseClath
                phaseStr = PhaseConv(phaseIndex)
            if phaseStr not in ices:
                ices.add(phaseStr)
                iceEOS[phaseStr] = GetIceEOS(P_MPa, T_K, phaseStr, porosType=Planet.Ocean.porosType[clath],
                                                         phiTop_frac=Planet.Ocean.phiMax_frac[clath],
                                                         Pclosure_MPa=Planet.Ocean.Pclosure_MPa[clath],
                                                         phiMin_frac=Planet.Ocean.phiMin_frac,
                                                         EXTRAP=Params.EXTRAP_ICE[phaseStr],
                                                         mixParameters={'mixFrac': Planet.Bulk.volumeFractionClathrate, 'JmixedRheologyConstant': Planet.Bulk.JmixedRheologyConstant})
            clathStable = iceEOS[phaseStr].fn_phase(P_MPa, T_K, grid=True)
            phases[clathStable == phaseIndex] = phaseIndex
        oceanListEOS.append(oceanEOS)
        phasesList.append(phases)
    # Obtain data for each planet
    for i, Planet in enumerate(PlanetList):
        oceanEOS = oceanListEOS[i]
        phases = phasesList[i]
        # Loop over each property in props_to_plot
        for prop in props_to_plot:
            if prop.prop in ['VS', 'GS']:
                # For VS and GS, these properties only apply to solid ice phases
                prop_data = np.empty((len(P_MPa), len(T_K))) * np.nan
            else:
                # Call the function dynamically with the required arguments
                prop_data = getattr(oceanEOS, prop.fn_name)(P_MPa, T_K, grid=True)
                # In case of VP or KS, we need to get the appropriate data since fn_seismic returns both VP and KS
                if prop.prop in ['VP', 'KS']:
                    prop_data = prop_data[prop.ocean_prop_index]
            # Now get data for all ice EOSs and replace in grid
            for iceStr in ices:
                ice = PhaseInv(iceStr)
                if np.any(phases == ice):
                    if prop.prop != 'sig':
                        ice_prop_data = getattr(iceEOS[ice], prop.fn_name)(P_MPa, T_K, grid=True)
                        if prop.prop in ['VP', 'KS', 'VS', 'GS']:
                            prop_data[np.where(phases == ice)] = ice_prop_data[prop.ice_prop_index][
                                np.where(phases == ice)]
                        else:
                            prop_data[np.where(phases == ice)] = ice_prop_data[
                                np.where(phases == ice)]
                    else:
                        # In the case of sig, we need to set a different way
                        prop_data[np.where(phases == ice)] = Planet.Ocean.sigmaIce_Sm[iceStr]
            # Add the data for this planet to the prop.prop_data list
            prop.prop_data.append(prop_data)
    # Go through every two combination and plot
    SinglePlanetPlot = (len(PlanetList) == 1)
    # If we are plotting for only one run, then do not compare it to anything
    if SinglePlanetPlot:
        combinations = [(0, 0)]
    else:
        combinations = itertools.combinations(range(0, PlanetList.size), 2)
    for comparison in combinations:
        FirstPlanetIndex = comparison[0]
        SecondPlanetIndex = comparison[1]
        FirstPlanet = PlanetList[FirstPlanetIndex]
        SecondPlanet = PlanetList[SecondPlanetIndex]
        # Create the figure grid
        fig = plt.figure(figsize=FigSize.vpvt)
        grid = GridSpec(rows, cols)
        axes = np.array([[fig.add_subplot(grid[i, j]) for j in range(cols)] for i in range(rows)])
        axf = axes.flatten()
        # Assign fig_location dynamically to each property and update ax
        for idx, prop in enumerate(props_to_plot):
            row, col = divmod(idx, cols)
            prop.fig_location = (row, col)
            ax = axes[row, col]
            prop.ax = ax  # Update ax in HydrosphereProp
        # Hide extra subplots if there are fewer properties than subplots
        for idx in range(len(props_to_plot), len(axf)):
            axf[idx].set_visible(False)
        if Style.GRIDS:
            [ax.grid() for ax in axf]
            [ax.set_axisbelow(False) for ax in axf]
        # Set overall figure title
        if Params.TITLES:
            if SinglePlanetPlot:
                fig.suptitle(f'{FirstPlanet.compStr}{FigLbl.PvTtitleHydro}')
            else:
                fig.suptitle(
                    f'Comparison of {FirstPlanet.compStr} and {SecondPlanet.compStr}{FigLbl.PvTtitleHydroComparison}')  # Now plot the data dynamically on the correct ax for each property
        for idx, prop in enumerate(props_to_plot):
            ax = prop.ax  # Get the axis for this property
            if SinglePlanetPlot:
                prop_data_to_plot = prop.prop_data[FirstPlanetIndex]
            else:
                firstPlanetData = prop.prop_data[FirstPlanetIndex]
                secondPlanetData = prop.prop_data[SecondPlanetIndex]
                prop_data_to_plot = (firstPlanetData - secondPlanetData) / ((firstPlanetData + secondPlanetData) / 2) * 100
                # Handle phase differences
                phase_difference = (phasesList[FirstPlanetIndex] != phasesList[SecondPlanetIndex])
                phase_different_indices = np.where(phase_difference)
                prop_data_to_plot[phase_different_indices] = np.nan
            # Plot property data
            cmap = Color.PvThydroCmap
            if SinglePlanetPlot:
                vmin = np.nanmin(prop_data_to_plot)
                vmax = np.nanmax(prop_data_to_plot)
                vmean = float(np.nanmean(prop_data_to_plot))
                norm = TwoSlopeNorm(vmin = vmin, vcenter = vmean, vmax = vmax)
            else:
                # Determine the limits for each property so we can peg 0 to middle of color plot
                abs_max = max(abs(np.nanmin(prop_data_to_plot)), abs(np.nanmax(prop_data_to_plot)))
                if abs_max == 0:
                    abs_max = 1e-14
                
                # Calculate a more reasonable range using 95th percentile to avoid extreme outliers
                valid_data = prop_data_to_plot[~np.isnan(prop_data_to_plot)]
                if len(valid_data) > 0:
                    # Use 95th percentile for both positive and negative values
                    vmin_95th = np.nanpercentile(valid_data, 5)   # 5th percentile for lower bound
                    vmax_95th = np.nanpercentile(valid_data, 95)  # 95th percentile for upper bound
                    vcenter = 0

                    # Make range symmetric around zero for better visualization
                    abs_display_max = max(abs(vmin_95th), abs(vmax_95th))
                    vmin_symmetric = -abs_display_max
                    vmax_symmetric = abs_display_max
                    
                    data_min = np.nanmin(prop_data_to_plot)
                    data_max = np.nanmax(prop_data_to_plot)
                else:
                    vmin_95th = -abs_max
                    vmax_95th = abs_max
                    vmin_symmetric = -abs_max
                    vmax_symmetric = abs_max
                    vcenter = 0
                    data_min = -abs_max
                    data_max = abs_max
                
                norm = TwoSlopeNorm(vmin=vmin_symmetric, vcenter=vcenter, vmax=vmax_symmetric)
            pcolormesh = ax.pcolormesh(T_K, P_MPa * FigLbl.PmultHydro, prop_data_to_plot, cmap=cmap,
                                       rasterized=FigMisc.PT_RASTER, norm=norm)
            # Add colorbar with extend to indicate values beyond visible range
            cbar = fig.colorbar(pcolormesh, ax=ax, extend='both')
            
            # Add annotations for 5th and 95th percentiles on the colorbar using lines and text
            if not SinglePlanetPlot and len(valid_data) > 0:
                # Normalize percentile values to [0, 1] range for the colorbar
                cbar_min, cbar_max = cbar.mappable.get_clim()
                p5_pos = (vmin_95th - cbar_min) / (cbar_max - cbar_min)
                p95_pos = (vmax_95th - cbar_min) / (cbar_max - cbar_min)

                # Threshold below which we offset to avoid overlap (in axes fraction)
                min_sep = 0.05  # adjust as needed

                # Initial positions
                p5_text_y = p5_pos
                p95_text_y = p95_pos

                # If too close, offset them slightly
                if abs(p95_text_y - p5_text_y) < min_sep:
                    offset = min_sep / 2
                    p5_text_y -= offset
                    p95_text_y += offset
                    # Clip to valid [0,1] range
                    p5_text_y = max(min(p5_text_y, 1.0), 0.0)
                    p95_text_y = max(min(p95_text_y, 1.0), 0.0)

                # Optional: background box to improve readability
                bbox_props = dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2', alpha=0.8)

                # Add text labels with safe positioning
                cbar.ax.text(1.25, p5_text_y, f'5th: {vmin_95th:.2g}',
                            transform=cbar.ax.transAxes, fontsize=8, color='red',
                            va='center', ha='left', bbox=bbox_props)

                cbar.ax.text(1.25, p95_text_y, f'95th: {vmax_95th:.2g}',
                            transform=cbar.ax.transAxes, fontsize=8, color='red',
                            va='center', ha='left', bbox=bbox_props)

                # Add horizontal lines (optional)
                cbar.ax.hlines(p5_pos, 0, 1, transform=cbar.ax.transAxes, colors='red', linestyles='--', linewidth=1)
                cbar.ax.hlines(p95_pos, 0, 1, transform=cbar.ax.transAxes, colors='red', linestyles='--', linewidth=1)
                # Display the actual data min and max in the colorbar label
                cbar.ax.text(1.05, 1.05, f'Max: {data_max:.2e}', transform=cbar.ax.transAxes, 
                            fontsize=8, verticalalignment='bottom')
                cbar.ax.text(1.05, -0.05, f'Min: {data_min:.2e}', transform=cbar.ax.transAxes, 
                            fontsize=8, verticalalignment='top')
            # Set labels, title, etc.
            ax.set_xlabel(FigLbl.Tlabel)
            ax.set_ylabel(FigLbl.PlabelHydro)
            ax.set_xlim([Tmin_K, Tmax_K])
            ax.set_ylim([Pmin_MPa, Pmax_MPa])
            ax.invert_yaxis()
            ax.set_title(prop.prop_label)

        # Plot geotherm on top of colormaps
        if FigMisc.SHOW_GEOTHERM:
            geothermList = [FirstPlanet, SecondPlanet] if not SinglePlanetPlot else [FirstPlanet]
            for eachPlanet in geothermList:
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
                [ax.plot(Tgeo, Pgeo, linewidth=Style.LW_geotherm, linestyle=Style.LS_geotherm, color=thisColor,
                         label=eachPlanet.label) for ax in axf]

                if Params.LEGEND and np.size(PlanetList) > 1:
                    handles, lbls = axes[-1, 0].get_legend_handles_labels()
                    axes[0, -1].legend(handles, lbls)

        plt.tight_layout()
        if SinglePlanetPlot:
            saveFile = Params.FigureFiles.vpvtHydro
        else:
            saveFile = Params.FigureFiles.comparisonFileGenerator(FirstPlanet.saveLabel, SecondPlanet.saveLabel, 'vpvtHydro')
        fig.savefig(saveFile, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Hydrosphere PT properties plot saved to file: {saveFile}')
        plt.close()

    return

def PlotCustomSolutionProperties(PlanetList, Params):
    Planet = PlanetList[0]
    # This is a hydrosphere-only plot for Reaktoro, so skip waterless bodies or bodies not utilizing Reaktoro
    if "CustomSolution" in Planet.Ocean.comp:
        fig = plt.figure(figsize=FigSize.vmant)
        grid = GridSpec(1, 1)
        allspeciesax = fig.add_subplot(grid[0, 0])
        aqueouspseciesax = fig.add_subplot(grid[0, 1])
        axs = [allspeciesax, aqueouspseciesax]
        if Style.GRIDS:
            allspeciesax.grid()
            allspeciesax.set_axisbelow(True)
        allspeciesax.set_xlabel("All species (moles)")
        allspeciesax.set_ylabel("Depth z (km)")
        aqueouspseciesax.set_xlabel("Aqueous species (moles per kg fluid)")

def PlotPvTPerpleX(PlanetList, Params):

    if os.path.dirname(Params.FigureFiles.vpvtPerpleX) != 'Comparison':
        Planet = PlanetList[0]
        if Planet.Sil.EOS is None:
            Planet.Sil.EOS = GetInnerEOS(Planet.Sil.mantleEOS, EOSinterpMethod=Params.lookupInterpMethod,
                        kThermConst_WmK=Planet.Sil.kTherm_WmK, HtidalConst_Wm3=Planet.Sil.Htidal_Wm3,
                        porosType=Planet.Sil.porosType, phiTop_frac=Planet.Sil.phiRockMax_frac,
                        Pclosure_MPa=Planet.Sil.Pclosure_MPa, phiMin_frac=Planet.Sil.phiMin_frac,
                        EXTRAP=Params.EXTRAP_SIL, etaSilFixed_Pas=Planet.Sil.etaRock_Pas, etaCoreFixed_Pas=[Planet.Core.etaFeSolid_Pas, Planet.Core.etaFeLiquid_Pas])
        INCLUDING_CORE = FigMisc.PVT_INCLUDE_CORE and Planet.Do.Fe_CORE
        if INCLUDING_CORE and Planet.Core.EOS is None:
            Planet.Core.EOS = GetInnerEOS(Planet.Core.coreEOS, EOSinterpMethod=Params.lookupInterpMethod, Fe_EOS=True,
                        kThermConst_WmK=Planet.Core.kTherm_WmK, EXTRAP=Params.EXTRAP_Fe,
                        wFeCore_ppt=Planet.Core.wFe_ppt, wScore_ppt=Planet.Core.wS_ppt, etaSilFixed_Pas=Planet.Sil.etaRock_Pas, etaCoreFixed_Pas=[Planet.Core.etaFeSolid_Pas, Planet.Core.etaFeLiquid_Pas])

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

def PlotMeltingCurves(PlanetList, Params):
    """
    Plot melting curves for multiple PlanetProfile models on the same figure.
    
    Args:
        PlanetList: List of Planet objects to plot melting curves for
        Params: Parameters object containing figure settings
    """
    if os.path.dirname(Params.FigureFiles.vmeltingCurves) != 'Comparison':
        # Determine pressure and temperature ranges similar to PlotHydroPhase
        if FigMisc.PminHydro_MPa is None:
            Pmin_MPa = np.min([Planet.P_MPa[0] for Planet in PlanetList])
        else:
            Pmin_MPa = FigMisc.PminHydro_MPa
        if FigMisc.PmaxHydro_MPa is None:
            Pmax_MPa = np.max([Planet.P_MPa[Planet.Steps.nHydro-1] for Planet in PlanetList])
        else:
            Pmax_MPa = FigMisc.PmaxHydro_MPa
        if FigMisc.TminMeltingCurve_K is None:
            if not np.any([Planet.Do.NO_OCEAN for Planet in PlanetList]):
                Tmin_K = np.min([np.min(Planet.T_K[:Planet.Steps.nHydro]) for Planet in PlanetList])
            else:
                Tmin_K = np.min([np.min(Planet.T_K) for Planet in PlanetList])
        else:
            Tmin_K = FigMisc.TminMeltingCurve_K
        if FigMisc.TmaxMeltingCurve_K is None:
            if np.all([Planet.Do.NO_OCEAN for Planet in PlanetList]):
                Tmax_K = np.max([Planet.Sil.THydroMax_K for Planet in PlanetList])
            else:
                Tmax_K = np.max([np.max(Planet.T_K[:Planet.Steps.nHydro])
                                 for Planet in PlanetList if not Planet.Do.NO_OCEAN])
            if Pmax_MPa <= Constants.PminHPices_MPa:
                Tmax_K = np.min([Tmax_K, 280]) # Set the cut off at 280K since melting curve won't be above this point for HP ices
        else:
            Tmax_K = FigMisc.TmaxMeltingCurve_K

        # Create pressure and temperature grids for melting curve calculation
        P_MPa = np.linspace(Pmin_MPa, Pmax_MPa, FigMisc.nPmeltingCurve)
        T_K = np.linspace(Tmin_K, Tmax_K, FigMisc.nTmeltingCurve)

        # Get unique ocean compositions and salinities
        comps = []
        for Planet in PlanetList:
            if "CustomSolution" in Planet.Ocean.comp:
                if Planet.Ocean.comp.split('=')[0] not in comps:
                    comps.append(Planet.Ocean.comp)
            else:
                comps.append(Planet.Ocean.comp)
        comps = np.unique(comps)

        # Create figure
        fig = plt.figure(figsize=FigSize.vmeltingCurves)
        ax = fig.add_subplot(1, 1, 1)
        
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)
        
        # Labels and titles
        ax.set_xlabel(FigLbl.Tlabel)
        ax.set_ylabel(FigLbl.PlabelHydro)
        
        # Set overall figure title
        if Params.TITLES:
            if Params.ALL_ONE_BODY:
                fig.suptitle(f'{PlanetList[0].name}{FigLbl.meltingCurvesTitle}', fontsize=FigLbl.hydroTitleSize)
            else:
                fig.suptitle(f'{FigLbl.meltingCurvesTitle}', fontsize=FigLbl.hydroTitleSize)

        # Get min and max salinities and temps for each comp for scaling (similar to PlotHydrosphereProps)
        wMinMax_ppt = {}
        TminMax_K = {}
        for comp in comps:
            if comp != 'none':
                wAll_ppt = [Planet.Ocean.wOcean_ppt for Planet in PlanetList if Planet.Ocean.comp == comp]
                wMinMax_ppt[comp] = [np.min(wAll_ppt), np.max(wAll_ppt)]
                Tall_K = [Planet.Bulk.Tb_K for Planet in PlanetList if Planet.Ocean.comp == comp]
                TminMax_K[comp] = [np.min(Tall_K), np.max(Tall_K)]
                # Reset to default if all models are the same or if desired
                if not FigMisc.RELATIVE_Tb_K or TminMax_K[comp][0] == TminMax_K[comp][1]:
                    TminMax_K[comp] = Color.Tbounds_K

        # Calculate and plot melting curves for each unique composition
        melting_curves = {}  # Store melting curves by composition
        for comp in comps:
            if comp == 'none':
                continue
                
            # Get the first planet with this composition to use for EOS
            comp_planets = [Planet for Planet in PlanetList if Planet.Ocean.comp == comp]
            if not comp_planets:
                continue
                
            ref_planet = comp_planets[0]
            
            # Load EOS for this composition
            oceanEOS = GetOceanEOS(ref_planet.Ocean.comp, ref_planet.Ocean.wOcean_ppt, P_MPa, T_K,
                                   ref_planet.Ocean.MgSO4elecType, rhoType=ref_planet.Ocean.MgSO4rhoType,
                                   scalingType=ref_planet.Ocean.MgSO4scalingType, FORCE_NEW=Params.FORCE_EOS_RECALC,
                                   phaseType=ref_planet.Ocean.phaseType, EXTRAP=Params.EXTRAP_OCEAN,
                                   sigmaFixed_Sm=ref_planet.Ocean.sigmaFixed_Sm, LOOKUP_HIRES=ref_planet.Do.OCEAN_PHASE_HIRES, kThermConst_WmK=ref_planet.Ocean.kThermWater_WmK,
                                   propsStepReductionFactor=ref_Planet.Ocean.propsStepReductionFactor)

            # Calculate phase diagram
            phases = oceanEOS.fn_phase(P_MPa, T_K, grid=True).astype(int)
            
            # Find melting curve (transition from ice to liquid)
            melting_points = []
            for i, P in enumerate(P_MPa):
                # Find where phase changes from ice (non-zero) to liquid (zero)
                phase_col = phases[i, :]
                transitions = np.where(np.diff(phase_col) == -1)[0]  # Ice to liquid transitions
                
                if len(transitions) > 0:
                    # Use the first transition (lowest temperature)
                    T_melt = T_K[transitions[0]]
                    melting_points.append((T_melt, P))
            
            if melting_points:
                melting_curves[comp] = np.array(melting_points)
            else:
                log.warning(f'No melting curve found for composition {comp}')

        # Plot melting curves with appropriate styling
        for comp in comps:
            if comp == 'none' or comp not in melting_curves:
                continue
                
            melting_data = melting_curves[comp]
            T_melt = melting_data[:, 0]
            P_melt = melting_data[:, 1]
            
            # Set style options (similar to PlotHydrosphereProps)
            if FigMisc.MANUAL_HYDRO_COLORS:
                Color.Tbounds_K = TminMax_K[comp]
                thisColor = Color.cmap[comp](Color.GetNormT(np.mean(T_melt)))
            else:
                thisColor = None
                
            if FigMisc.SCALE_HYDRO_LW and wMinMax_ppt[comp][0] != wMinMax_ppt[comp][1]:
                thisLW = Style.GetLW(np.mean([Planet.Ocean.wOcean_ppt for Planet in PlanetList if Planet.Ocean.comp == comp]), wMinMax_ppt[comp])
            else:
                thisLW = FigMisc.MELTING_CURVE_LINE_WIDTH
            if FigMisc.LS_SOLID_MELTING_CURVES:
                thisLS = 'solid'
            else:
                thisLS = Style.LS[comp]
            # Create label for this composition
            if "CustomSolution" in comp:
                comp_label = comp.split('=')[0].strip()
            else:
                comp_label = comp
                
            # Plot melting curve
            ax.plot(T_melt, P_melt * FigLbl.PmultHydro, 
                   label=comp_label, color=thisColor, linewidth=thisLW,
                   linestyle=thisLS)

        # Mark model melting points if requested
        if FigMisc.MARK_MODEL_POINTS:
            for Planet in PlanetList:
                if Planet.Ocean.comp == 'none':
                    continue
                    
                # Get model melting point (Tb_K, Pb_MPa)
                T_model = Planet.Bulk.Tb_K
                P_model = Planet.Pb_MPa
                
                # Set style for this planet
                if FigMisc.MANUAL_HYDRO_COLORS:
                    Color.Tbounds_K = TminMax_K[Planet.Ocean.comp]
                    thisColor = Color.cmap[Planet.Ocean.comp](Color.GetNormT(T_model))
                else:
                    thisColor = None
                    
                if FigMisc.SCALE_HYDRO_LW and wMinMax_ppt[Planet.Ocean.comp][0] != wMinMax_ppt[Planet.Ocean.comp][1]:
                    thisLW = Style.GetLW(Planet.Ocean.wOcean_ppt, wMinMax_ppt[Planet.Ocean.comp])
                else:
                    thisLW = FigMisc.MELTING_CURVE_LINE_WIDTH
                
                # Plot model point
                ax.scatter(T_model, P_model * FigLbl.PmultHydro,
                          color=thisColor, s=FigMisc.MODEL_POINT_SIZE,
                          marker='o', edgecolors='black', linewidth=0.5)

        # Plot geotherms if requested
        if FigMisc.SHOW_GEOTHERM:
            for Planet in PlanetList:
                if Planet.Ocean.comp == 'none':
                    continue
                    
                # Set style for this planet
                if FigMisc.MANUAL_HYDRO_COLORS:
                    Color.Tbounds_K = TminMax_K[Planet.Ocean.comp]
                    thisColor = Color.cmap[Planet.Ocean.comp](Color.GetNormT(Planet.Bulk.Tb_K))
                else:
                    thisColor = Color.geothermHydro
                    
                if FigMisc.SCALE_HYDRO_LW and wMinMax_ppt[Planet.Ocean.comp][0] != wMinMax_ppt[Planet.Ocean.comp][1]:
                    thisLW = Style.GetLW(Planet.Ocean.wOcean_ppt, wMinMax_ppt[Planet.Ocean.comp])
                else:
                    thisLW = Style.LW_geotherm
                
                # Plot geotherm
                if Planet.Do.NO_DIFFERENTIATION or Planet.Do.PARTIAL_DIFFERENTIATION:
                    Pgeo = Planet.P_MPa * FigLbl.PmultHydro
                    Tgeo = Planet.T_K
                else:
                    Pgeo = Planet.P_MPa[:Planet.Steps.nHydro] * FigLbl.PmultHydro
                    Tgeo = Planet.T_K[:Planet.Steps.nHydro]
                
                ax.plot(Tgeo, Pgeo, linewidth=thisLW, linestyle=Style.LS_geotherm,
                       color=thisColor, alpha=0.7, label=f'{Planet.label} geotherm')

        # Set axis limits
        ax.set_xlim([Tmin_K, Tmax_K])
        ax.set_ylim([Pmin_MPa * FigLbl.PmultHydro, Pmax_MPa * FigLbl.PmultHydro])
        ax.invert_yaxis()
        
        # Add legend if requested
        if Params.LEGEND:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels, loc='upper left')
        
        plt.tight_layout()
        fig.savefig(Params.FigureFiles.vmeltingCurves, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Melting curves plot saved to file: {Params.FigureFiles.vmeltingCurves}')
        plt.close()

    return
