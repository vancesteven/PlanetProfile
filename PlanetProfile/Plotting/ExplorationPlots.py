import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc
from PlanetProfile.Plotting.EssentialHelpers import *
from PlanetProfile.Thermodynamics.Reaktoro.CustomSolution import SetupCustomSolutionPlotSettings

import logging
log = logging.getLogger('PlanetProfile')

def GenerateExplorationPlots(ExplorationList, Params):
    """
    Generate all exploration plots for a given list of results.
    """
    PLOT_EXPLORATION = not (Params.Explore.nx == 1 or Params.Explore.ny == 1) # If either axis is of size 1, then we cannot plot 
    # Setup CustomSolution plot settings
    all_ocean_comps = []
    for Exploration in ExplorationList:
        all_ocean_comps.extend(np.array(Exploration.base.oceanComp).flatten())
    Params = SetupCustomSolutionPlotSettings(np.array(all_ocean_comps), Params)
    if PLOT_EXPLORATION:
        # Use multi-subplot function only for multiple z-variables
        if isinstance(Params.Explore.zName, list):
            PlotExploreOgramMultiSubplot(ExplorationList, Params)
        else:
            # Use original single-plot function for single z-variable
            for Exploration in ExplorationList:
                Exploration.zName = Params.Explore.zName
            PlotExploreOgram(ExplorationList, Params)
        # Now we plot the ZbD plots (must plot after exploreogram plots since we change the x and y variables)
        if Params.PLOT_Zb_D:
            if isinstance(Params.Explore.zName, list):
                figNames = Params.FigureFiles.exploreZbD + []
                for zName, figName in zip(Params.Explore.zName, figNames):
                    for Exploration in ExplorationList:
                        Exploration.zName = zName
                    Params.FigureFiles.exploreZbD = figName
                    PlotExploreOgramZbD(ExplorationList, Params)
            else:
                for Exploration in ExplorationList:
                    Exploration.zName = Params.Explore.zNameZbD
                PlotExploreOgramZbD(ExplorationList, Params)
        # Now we plot the ZbY plots (ice shell thickness vs specified property with ocean thickness as color)
        if Params.PLOT_Zb_Y:
            if isinstance(Params.Explore.zName, list):
                figNames = Params.FigureFiles.exploreZbY + []
                for zName, figName in zip(Params.Explore.zName, figNames):
                    for Exploration in ExplorationList:
                        Exploration.zName = zName
                    Params.FigureFiles.exploreZbY = figName
                    PlotExploreOgramXZ(ExplorationList, Params)
            else:
                for Exploration in ExplorationList:
                    Exploration.zName = Params.Explore.zNameZbY
                PlotExploreOgramXZ(ExplorationList, Params)
        PlotExploreOgramDsigma(ExplorationList, Params)
        PlotExploreOgramLoveComparison(ExplorationList, Params)

    
def PlotExploreOgramDsigma(results_list, Params):
    """
    TODO: Create comparison structure creation function that:
    - Combines data from all individual results into a single comparison structure
    - Sets appropriate comparison titles and labels
    - Appends to results_list so plotting treats it like any other result
    """
    
    # Step 1: Configure FigLbl using existing system (original pattern)
    results_list[0].xName = 'D_km'
    results_list[0].yName = 'sigmaMean_Sm' 
    results_list[0].zName = 'zb_km'
    
    FigLbl.SetExploration(results_list[0].base.bodyname, results_list[0].xName,
                          results_list[0].yName, results_list[0].zName)
    if not FigMisc.TEX_INSTALLED:
        FigLbl.StripLatex()

    # Step 2: Create plots for each result (individual + comparison if present)
    for i, result in enumerate(r for r in results_list if not r.base.NO_H2O):
        result.xName = 'D_km'
        result.yName = 'sigmaMean_Sm'
        result.zName = 'zb_km'
        
        # Determine if this is the comparison plot (last result when COMPARE is enabled)
        is_comparison_plot = (Params.COMPARE and i == len([r for r in results_list if not r.base.NO_H2O]) - 1 
                             and len(results_list) > 1)
        
        # Extract and validate data using helper
        plot_data = extract_and_validate_plot_data(result_obj = result, x_field = result.xName, y_field = result.yName, c_field = result.zName,
                                                  x_multiplier = FigLbl.xMultExplore, y_multiplier = FigLbl.yMultExplore, c_multiplier = 1.0, contour_multiplier = 1.0,
                                                  custom_x_axis = FigLbl.xCustomAxis, custom_y_axis = FigLbl.yCustomAxis)
        
        if len(plot_data['x']) == 0:
            log.warning(f"No valid data points for {result.base.bodyname}")
            continue
        
        # Create figure and axis (original pattern)
        fig = plt.figure(figsize=FigSize.explore)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)

        # Set up axis properties using FigLbl + D/sigma specific overrides
        if Params.TITLES:
            # Use comparison title for comparison plots, regular title for individual plots
            if is_comparison_plot:
                fig.suptitle(FigLbl.exploreCompareTitle)  # Comparison-specific title
            else:
                fig.suptitle(FigLbl.explorationDsigmaTitle)  # Individual plot title
        ax.set_xlabel(FigLbl.xLabelExplore)
        ax.set_ylabel(FigLbl.yLabelExplore)
        # D/sigma specific settings (override standard)
        ax.set_xscale('linear')
        ax.set_yscale('log')
        ax.set_ylim(FigMisc.DSIGMA_YLIMS)
        ax.set_xlim([np.min(plot_data['x']), np.max(plot_data['x'])])
        
        # Get ocean composition data for composition lines and legends
        ocean_comp = result.base.oceanComp.flatten()
        
        # Draw composition lines if enabled (using helper for 20+ line block)
        plotted_labels = set()
        if FigMisc.DRAW_COMPOSITION_LINE:
            plotted_labels = draw_ocean_composition_lines(
                ax, plot_data['x'], plot_data['y'], plot_data['c'], ocean_comp,
                use_manual_colors=FigMisc.MANUAL_HYDRO_COLORS,
                line_width=FigMisc.DSIGMA_COMP_LINE_WIDTH,
                line_alpha=FigMisc.DSIGMA_COMP_LINE_ALPHA
            )
        
        # Create scatter plot with clear if/else for different modes (original logic)
        if FigMisc.SHOW_ICE_THICKNESS_DOTS:
            # Ice thickness coloring mode
            pts = ax.scatter(plot_data['x'], plot_data['y'], c=plot_data['c'],
                            cmap=FigMisc.DSIGMA_ICE_THICKNESS_CMAP, 
                            marker=Style.MS_Induction, s=Style.MW_Induction**2, 
                            edgecolors=FigMisc.DSIGMA_DOT_EDGE_COLOR, 
                            linewidths=FigMisc.DSIGMA_DOT_EDGE_WIDTH, 
                            zorder=3)
            
            # Create ice thickness legend using helper (15+ line block)
            ice_legend = create_ice_thickness_colorbar(
                ax, plot_data['c'],
                cmap_name=FigMisc.DSIGMA_ICE_THICKNESS_CMAP,
            )
                
        else:
            # Standard colorbar mode
            # Get first ocean composition for colormap selection
            first_comp = ocean_comp[0] if len(ocean_comp) > 0 else 'default'
            cmap_to_use = Color.cmap.get(first_comp, 'default')
            
            pts = ax.scatter(plot_data['x'], plot_data['y'], c=plot_data['c'],
                            cmap=cmap_to_use, marker=Style.MS_Induction, 
                            s=Style.MW_Induction**2, zorder=3)

            # Add colorbar (original logic)
            cbar = fig.colorbar(pts, ax=ax, format=FigLbl.cbarFmt)
            if plot_data['c'] is not None and len(plot_data['c']) > 0:
                # Append min and max values to colorbar for reading convenience
                new_ticks = np.insert(np.append(cbar.get_ticks(), np.max(plot_data['c'])), 0, np.min(plot_data['c']))
                cbar.set_ticks(np.unique(new_ticks))
            cbar.set_label(FigLbl.cbarLabelExplore, size=12)
        
        # Add composition legend if enabled (using helper)
        if FigMisc.DRAW_COMPOSITION_LINE and Params.LEGEND:
            add_composition_legend(ax, ocean_comp,
                                 show_ice_thickness_legend=FigMisc.SHOW_ICE_THICKNESS_DOTS,
                                 font_size=FigMisc.DSIGMA_COMP_LEGEND_FONT_SIZE)

        # Save figure (original pattern)
        plt.tight_layout()
        fig.savefig(Params.FigureFiles.exploreDsigma, format=FigMisc.figFormat, 
                   dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Plot saved to file: {Params.FigureFiles.exploreDsigma}')
        plt.close()


def PlotExploreOgram(results_list, Params, ax=None):
    # Set up basic figure labels using first result
    first_result = results_list[0]
    FigLbl.SetExploration(first_result.base.bodyname, first_result.xName,
                          first_result.yName, first_result.zName, Params.Explore.contourName, first_result.excName)
    if not FigMisc.TEX_INSTALLED:
        FigLbl.StripLatex()
    
    # Plot each result (individual exploration results + optional comparison)
    last_index = len(results_list) - 1
    for i, exploration in enumerate(results_list):
        
        
        # Detect if this is a comparison plot (last result when COMPARE=True)
        is_comparison_plot = Params.COMPARE and i == last_index
        
        # Handle axis creation - support for multi-subplot usage
        create_new_figure = ax is None
        if create_new_figure:
            # Create new figure/axes and save it (original behavior)
            fig = plt.figure(figsize=FigSize.explore)
            grid = GridSpec(1, 1)
            ax = fig.add_subplot(grid[0, 0])
        elif isinstance(ax, list):
            ax = ax[i]
            fig = ax.figure  # Get figure from provided axis
        else:
            fig = ax.figure  # Get figure from provided axis
        
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)
        
        # Set title based on plot type - only set individual titles, not figure titles
        if exploration.zName == 'CMR2calc':
            FigLbl.SetExploreTitle(exploration.base.bodyname, exploration.zName, exploration.CMR2str)
        
        # Only set figure title if we created the figure (not in subplot mode)
        if Params.TITLES and create_new_figure:
            if is_comparison_plot:
                fig.suptitle(FigLbl.exploreCompareTitle)
            else:
                fig.suptitle(FigLbl.explorationTitle)
        elif Params.TITLES and not create_new_figure:
            # In subplot mode, set the subplot title instead of figure title
            ax.set_title(FigLbl.cbarLabelExplore)
        
        # Set up axis labels and scales
        ax.set_xlabel(FigLbl.xLabelExplore)
        ax.set_ylabel(FigLbl.yLabelExplore)
        ax.set_xscale(FigLbl.xScaleExplore)
        ax.set_yscale(FigLbl.yScaleExplore)
        
        # Extract x, y, z, and contour data using enhanced helper function
        contour_field = Params.Explore.contourName
            
        plot_data = extract_and_validate_plot_data(exploration, exploration.xName, exploration.yName, exploration.zName, contour_field,
                                  x_multiplier=FigLbl.xMultExplore, y_multiplier=FigLbl.yMultExplore, c_multiplier=FigLbl.zMultExplore, 
                                  contour_multiplier=1.0, custom_x_axis = FigLbl.xCustomAxis, custom_y_axis = FigLbl.yCustomAxis)
        original_shape = plot_data['original_shape']
        x = plot_data['x'].reshape(original_shape)
        y = plot_data['y'].reshape(original_shape)
        z = plot_data['c'].reshape(original_shape)
        
        # Set axis limits
        ax.set_xlim([np.min(x), np.max(x)])
        ax.set_ylim([np.min(y), np.max(y)])
        
        # Only keep data points for which a valid model was determined
        z_shape = np.shape(z)
        z = np.reshape(z, -1).astype(np.float64)
        INVALID = np.logical_not(np.reshape(exploration.base.VALID, -1))
        z[INVALID] = np.nan
        # Return data to original organization
        z = np.reshape(z, z_shape)
        

        
        contour_levels=[]
        # Calculate valid data range for colorbar
        z_valid = z[z == z]  # Exclude NaNs
        if np.size(z_valid) > 0:
                            
            data_min = np.min(z_valid)
            data_max = np.max(z_valid)
            # Check if we should pin colormap center to zero
            if (data_min is not None and data_max is not None and 
                exploration.zName in FigLbl.cMapZero and 
                data_min < 0 and data_max > 0):
                # Pin colormap center to zero for variables that can be positive/negative
                abs_max = max(abs(data_min), abs(data_max))
                vmin = -abs_max
                vmax = abs_max
            else:
                # Use actual data range for colormap
                vmin = data_min
                vmax = data_max
                        # Create the main plot
            mesh = ax.pcolormesh(x, y, z, shading='auto', cmap=Color.cmap['default'], 
                            vmin=vmin, vmax=vmax, rasterized=FigMisc.PT_RASTER)
             # Add colorbar
            cbar = fig.colorbar(mesh, ax=ax, format=FigLbl.cbarFmt)
            # Get colorbar tick values to use as contour levels
            cbar_ticks = cbar.get_ticks()
            if FigLbl.cTicksSpacingExplore is not None:
                cbar_ticks = np.arange(np.ceil(data_min / FigLbl.cTicksSpacingExplore) * FigLbl.cTicksSpacingExplore, np.floor(data_max / FigLbl.cTicksSpacingExplore) * FigLbl.cTicksSpacingExplore + 0.0001, FigLbl.cTicksSpacingExplore)
            
                if len(cbar_ticks) > 10:
                    tooManyTicks = True
                    cTicksSpacingExplore = FigLbl.cTicksSpacingExplore
                    while tooManyTicks:
                        cTicksSpacingExplore = cTicksSpacingExplore + FigLbl.cTicksSpacingExplore
                        cbar_ticks = np.arange(np.ceil(data_min / cTicksSpacingExplore) * cTicksSpacingExplore, np.floor(data_max / cTicksSpacingExplore) * cTicksSpacingExplore + 0.0001, cTicksSpacingExplore)
                        if len(cbar_ticks) <= 10:
                            tooManyTicks = False
                cbar.set_ticks(cbar_ticks)
            if vmin is not None and vmax is not None:
                contour_levels = cbar_ticks[(cbar_ticks > vmin) & (cbar_ticks < vmax)]
            else:
                contour_levels = cbar_ticks
            if len(contour_levels) > 0:
                contourMin = np.min(contour_levels)
                contourMax = np.max(contour_levels)
                if FigLbl.cSpacingExplore is not None:
                    contour_levels = np.arange(np.ceil(data_min / FigLbl.cSpacingExplore) * FigLbl.cSpacingExplore, np.floor(data_max / FigLbl.cSpacingExplore) * FigLbl.cSpacingExplore + 0.0001, FigLbl.cSpacingExplore)
            # Determine contour variable and levels
            contour_variable = z  # Default to using z variable for contours
            if plot_data['contour'] is not None:
                # Use the separate contour variable
                contour_variable = plot_data['contour'].reshape(original_shape)
                
                # Apply same validity mask to contour variable
                contour_shape = np.shape(contour_variable)
                contour_variable = np.reshape(contour_variable, -1).astype(np.float64)
                contour_variable[INVALID] = np.nan
                contour_variable = np.reshape(contour_variable, contour_shape)
                # Calculate valid data range for colorbar
                contour_valid = contour_variable[contour_variable == contour_variable]  # Exclude NaNs
                if np.size(contour_valid) > 0:
                    contourMin = np.min(contour_valid)
                    contourMax = np.max(contour_valid)
                else:
                    contourMin = contourMax = None
                if contourMin is not None and contourMax is not None:
                    contour_levels = np.linspace(contourMin, contourMax, 5)
                
        if FigLbl.cSpacingExplore is not None:
            if len(contour_levels) > 20:
                tooManyContours = True
                cSpacingExplore = FigLbl.cSpacingExplore
                while tooManyContours:
                    cSpacingExplore = 9
                    contour_levels = np.arange(np.ceil(contourMin / cSpacingExplore) * cSpacingExplore, np.floor(contourMax / cSpacingExplore) * cSpacingExplore + 0.0001, cSpacingExplore)
                    if len(contour_levels) <= 20:
                        tooManyContours = False
            
        # Create contours at colorbar tick levels
        if len(contour_levels) > 0:
            cont = ax.contour(x, y, contour_variable, levels=contour_levels, colors='black')
            lbls = plt.clabel(cont, fmt=FigLbl.cfmt, fontsize=20, colors='black')
                            # Get pixel positions of all labels
            positions = np.array([ax.transData.transform(txt.get_position()) for txt in lbls])
            n = len(lbls)
            
            close_to_another = np.zeros(n, dtype=bool)
            
            for i in range(n):
                for j in range(i+1, n):
                    dist = np.linalg.norm(positions[i] - positions[j])
                    if dist < 25:
                        close_to_another[i] = True
                        close_to_another[j] = True
            
            for i, txt in enumerate(lbls):
                if close_to_another[i]:
                    txt.set_fontsize(2)
        # Add the min and max values to the colorbar for reading convenience
        if np.size(z_valid) > 0:
            # Use the adjusted vmin/vmax that may have been zero-pinned
            mesh.set_clim(vmin=vmin, vmax=vmax)
            
            # Check if we should truncate the colorbar (only for zero-pinned colormaps)
            if (exploration.zName in FigLbl.cMapZero and 
                data_min is not None and data_max is not None and
                data_min < 0 and data_max > 0):
                # Truncate colorbar to only show the portion with actual data
                # Calculate the fraction of the colormap that corresponds to actual data
                colormap_range = vmax - vmin
                data_range = data_max - data_min
                
                # Calculate the position of data within the full colormap
                data_start_fraction = (data_min - vmin) / colormap_range
                data_end_fraction = (data_max - vmin) / colormap_range
                existing_ticks = cbar.get_ticks()
                # Truncate the colorbar
                cbar.ax.set_ylim(data_start_fraction, data_end_fraction)
                
                valid_ticks = existing_ticks[(existing_ticks >= data_min) & (existing_ticks <= data_max)]
                # Add min and max values to the filtered ticks (using actual data range)
                new_ticks = np.insert(np.append(valid_ticks, data_max), 0, data_min)
                cbar.set_ticks(np.unique(new_ticks))
            else:
                # Normal behavior - filter ticks but don't truncate colorbar
                existing_ticks = cbar.get_ticks()
                if len(existing_ticks) > 1:
                    tick_diff = existing_ticks[1] - existing_ticks[0]
                    valid_ticks = existing_ticks[(existing_ticks >= data_min) & (existing_ticks <= data_max)]
                    if valid_ticks[0] - data_min < tick_diff * 0.2:
                        valid_ticks = np.delete(valid_ticks, 0)
                    if data_max - valid_ticks[-1] < tick_diff * 0.2:
                        valid_ticks = np.delete(valid_ticks, -1)
                else:
                    valid_ticks = existing_ticks
                # Add min and max values to the filtered ticks (using actual data range)
                new_ticks = np.insert(np.append(valid_ticks, data_max), 0, data_min)
                cbar.set_ticks(np.unique(new_ticks))  
            
            cbar.ax.tick_params(labelsize=20)
            if create_new_figure:
                cbar.set_label(FigLbl.cbarLabelExplore, size=25)
        
    
    # Save the plot only if we created a new figure
    if create_new_figure:
        plt.tight_layout()
        fig.savefig(Params.FigureFiles.explore, format=FigMisc.figFormat,
                   dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Plot saved to file: {Params.FigureFiles.explore}')
        plt.close()
    
    # Return the axis for multi-subplot usage (only if ax was provided)
    if not create_new_figure:
        return ax
    

def PlotExploreOgramMultiSubplot(results_list, Params):
    """
    Create multiple exploration subplots with different z-variables in a single figure.
    
    This function should only be called when Params.Explore.zName is a list with 
    multiple z-variables. It arranges them in a square-ish grid layout, calls 
    PlotExploreOgram for each z-variable, and lets that function handle 
    individual subplot titles. Only adds an overall figure title at the top.
    
    Args:
        results_list: List of ExplorationResults objects
        Params: Configuration parameters (Params.Explore.zName must be a list)
    """
    
    # Params.Explore.zName should be a list when this function is called
    zNames = []
    zExcNames = []
    n_subplots = 0
    for i, z_name in enumerate(Params.Explore.zName):
        if z_name in FigLbl.zNamePlotRealImag:
            real_name = FigLbl.zNamePlotRealImag[z_name][0]
            imag_name = FigLbl.zNamePlotRealImag[z_name][1]
            nExc_subplot, excNames = count_plottable_excitations(results_list[0].induction.calcedExc, Params.Induct)
            for excName in excNames:
                zNames.append(real_name)
                zExcNames.append(excName)
                zNames.append(imag_name)
                zExcNames.append(excName)
                n_subplots += 2
        elif 'Induction' in z_name:
            nExc_subplot, excNames = count_plottable_excitations(results_list[0].induction.calcedExc, Params.Induct)
            n_subplots += nExc_subplot
            for excName in excNames:
                zNames.append(z_name)
                zExcNames.append(excName)
        else:
            n_subplots += 1
            zNames.append(z_name)
            zExcNames.append(None)
    
    n_cols = int(np.ceil(np.sqrt(n_subplots)))
    n_rows = int(np.ceil(n_subplots / n_cols))
    
    # Calculate figure size with scaling
    base_size = FigSize.explore
    scale_factor = 1
    fig_width = base_size[0] * n_cols * scale_factor
    fig_height = base_size[1] * n_rows * scale_factor
    
    # Create figure with subplots
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    
    # Create subplots and let PlotExploreOgram handle individual plots
    axes = []
    for i in range(n_subplots):
        row = i // n_cols
        col = i % n_cols
        ax = fig.add_subplot(n_rows, n_cols, i + 1)
        axes.append(ax)
        
        # Set z-variable for all results
        z_name = zNames[i]
        exc_name = zExcNames[i]
        for result in results_list:
            result.zName = z_name
            result.excName = exc_name
        
        # Call PlotExploreOgram with this axis - let it handle the subplot title
        PlotExploreOgram(results_list, Params, ax=ax)
        
        # Add subplot label (a, b, c, etc.) if enabled
        if FigMisc.SUBPLOT_LABELS:
            # Generate label: 'a)', 'b)', 'c)', etc.
            label = f"{chr(ord('a') + i)}"
            ax.text(FigMisc.SUBPLOT_LABEL_X, FigMisc.SUBPLOT_LABEL_Y, label, 
                   transform=ax.transAxes, fontsize=FigMisc.SUBPLOT_LABEL_FONTSIZE,
                   fontweight='bold', ha='left', va='top',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        
        # After plotting, hide axis labels selectively
        is_bottom_row = (row == n_rows - 1) or (i >= n_subplots - n_cols)
        is_left_column = (col == 0) 
        
        if not is_bottom_row:
            ax.set_xlabel('')
            ax.tick_params(axis='x', labelbottom=False)
        
        if not is_left_column:
            ax.set_ylabel('')
            ax.tick_params(axis='y', labelleft=False)
    
    # Hide unused subplots
    total_subplots = n_rows * n_cols
    for i in range(n_subplots, total_subplots):
        ax = fig.add_subplot(n_rows, n_cols, i + 1)
        ax.set_visible(False)
        # Set overall title if configured
    fig.suptitle(FigLbl.subplotExplorationTitle, fontsize=15)
    # Save the combined figure
    plt.tight_layout()

    
    fig.savefig(Params.FigureFiles.exploreMultiSubplot, format=FigMisc.figFormat,
               dpi=FigMisc.dpi, metadata=FigLbl.meta)
    log.debug(f'Multi-subplot exploration plot saved to file: {Params.FigureFiles.exploreMultiSubplot}')
    plt.close()


def PlotExploreOgramZbD(results_list, Params):
    # Set up basic figure labels using first result
    first_result = results_list[0]
    FigLbl.SetExploration(first_result.base.bodyname, 'zb_km', 'D_km', first_result.zName)
    if not FigMisc.TEX_INSTALLED:
        FigLbl.StripLatex()
    
    # Plot each result (individual exploration results + optional comparison)
    last_index = len(results_list) - 1
    for i, result in enumerate(results_list):
        
        # Skip results with no H2O
        if result.base.NO_H2O:
            continue
            
        # Set axis variable names for this plot type
        _, result.xName = getIceShellThickness(result)
        result.yName = 'D_km' 
        # zName should already be set by user
        
        # Detect if this is a comparison plot (last result when COMPARE=True)
        is_comparison_plot = Params.COMPARE and i == last_index
        
        # Extract and validate data using helper
        plot_data = extract_and_validate_plot_data(result_obj = result, x_field = result.xName, y_field = result.yName, c_field = result.zName,
                                                  x_multiplier = FigLbl.xMultExplore, y_multiplier = FigLbl.yMultExplore, c_multiplier = FigLbl.zMultExplore,
                                                  custom_x_axis = FigLbl.xCustomAxis, custom_y_axis = FigLbl.yCustomAxis)
        
        if len(plot_data['x']) == 0:
            log.warning(f"No valid data points for {result.base.bodyname}")
            continue
        
        # Create figure and axis
        fig = plt.figure(figsize=FigSize.explore)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)
        
        # Set title based on plot type
        if Params.TITLES:
            if is_comparison_plot:
                fig.suptitle(FigLbl.exploreCompareTitle)
            else:
                fig.suptitle(FigLbl.explorationTitle)
        
        # Set up axis labels and scales
        ax.set_xlabel(FigLbl.xLabelExplore)
        ax.set_ylabel(FigLbl.yLabelExplore)
        ax.set_xscale('linear')  # ZbD plots use linear scales
        ax.set_yscale('linear')
        
        # Extract plot data
        x, y, z = plot_data['x'], plot_data['y'], plot_data['c']
        
        # Set axis limits with padding
        if np.size(x) > 0:
            x_range = np.max(x) - np.min(x)
            y_range = np.max(y) - np.min(y)
            x_padding = x_range * FigMisc.ZBD_AXIS_PADDING
            y_padding = y_range * FigMisc.ZBD_AXIS_PADDING
            ax.set_xlim([np.min(x) - x_padding, np.max(x) + x_padding])
            ax.set_ylim([np.min(y) - y_padding, np.max(y) + y_padding])
        
        # Get ocean composition data for composition lines
        ocean_comp = result.base.oceanComp.flatten()
        
        
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
                    # Draw composition lines if enabled using helper
        if FigMisc.DRAW_COMPOSITION_LINE:
            draw_ocean_composition_lines(ax, x, y, z, ocean_comp, 
                                       FigMisc.MANUAL_HYDRO_COLORS,
                                       FigMisc.ZBD_COMP_LINE_WIDTH, FigMisc.ZBD_COMP_LINE_ALPHA)
            
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
                      linewidths=FigMisc.ZBD_DOT_EDGE_WIDTH,
                      zorder=3)
        
        # Add legends if enabled
        if FigMisc.DRAW_COMPOSITION_LINE and Params.LEGEND:
            ax.legend(fontsize=FigMisc.ZBD_COMP_LEGEND_FONT_SIZE)
        # Save the plot
        plt.tight_layout()
        fig.savefig(Params.FigureFiles.exploreZbD, format=FigMisc.figFormat,
                   dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Plot saved to file: {Params.FigureFiles.exploreZbD}')
        plt.close()


def PlotExploreOgramXZ(results_list, Params):
    """
    Plots scatter showing any x variable vs any z variable with ocean thickness as color,
    with optional ocean composition lines. The x variable is configurable via Params.XZPLOT_X_VARIABLE
    (defaults to 'D_km' for ocean thickness), and the y variable is taken from the result's zName.
    Ocean thickness (zb_km) is used as the color variable.
    
    Similar to PlotExploreOgramZbD but with configurable axes and different color variable.
    
    Usage:
        # To change the x variable, set Params.XZPLOT_X_VARIABLE before calling this function:
        # Params.XZPLOT_X_VARIABLE = 'Tb_K'  # Use ocean temperature as x-axis
        # Params.XZPLOT_X_VARIABLE = 'wOcean_ppt'  # Use ocean salinity as x-axis
        # Params.XZPLOT_X_VARIABLE = 'Pseafloor_MPa'  # Use seafloor pressure as x-axis
    Args:
        results_list: List of ExplorationResults objects (individual + optional comparison)
        Params: Configuration parameters
    """ 
    
    # Set up basic figure labels using first result
    first_result = results_list[0]
    first_result.yName = first_result.zName
    # Use configurable x variable instead of hardcoded 'D_km'
    x_variable = Params.XZPLOT_X_VARIABLE
    FigLbl.SetExploration(first_result.base.bodyname, x_variable, first_result.yName, 'zb_km')
    if not FigMisc.TEX_INSTALLED:
        FigLbl.StripLatex()
    
    # Plot each result (individual exploration results + optional comparison)
    last_index = len(results_list) - 1
    for i, result in enumerate(results_list):
        
        # Skip results with no H2O
        if result.base.NO_H2O:
            continue
        # Use configurable x variable instead of hardcoded 'D_km'
        result.xName = Params.XZPLOT_X_VARIABLE
        # Set axis variable names for this plot type
        _, result.zName = getIceShellThickness(result)
        
        # Detect if this is a comparison plot (last result when COMPARE=True)
        is_comparison_plot = Params.COMPARE and i == last_index
        
        # Extract and validate data using helper
        plot_data = extract_and_validate_plot_data(result_obj = result, x_field = result.xName, y_field = result.yName, c_field = result.zName,
                                                  x_multiplier = FigLbl.xMultExplore, y_multiplier = FigLbl.yMultExplore, c_multiplier = FigLbl.zMultExplore,
                                                  custom_x_axis = FigLbl.xCustomAxis, custom_y_axis = FigLbl.yCustomAxis)
        
        if len(plot_data['x']) == 0:
            log.warning(f"No valid data points for {result.base.bodyname}")
            continue
        
        # Create figure and axis
        fig = plt.figure(figsize=FigSize.explore)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)
        
        # Set title based on plot type
        if Params.TITLES:
            if is_comparison_plot:
                fig.suptitle(FigLbl.exploreCompareTitle)
            else:
                fig.suptitle(FigLbl.explorationYTitle)
        
        # Set up axis labels and scales
        ax.set_xlabel(FigLbl.xLabelExplore)
        ax.set_ylabel(FigLbl.yLabelExplore)
        ax.set_xscale('linear')  # ZbY plots use linear scales
        ax.set_yscale('linear')
        
        # Extract plot data
        x, y, z = plot_data['x'], plot_data['y'], plot_data['c']
        
        # Set axis limits with padding
        if np.size(x) > 0:
            x_range = np.nanmax(x) - np.nanmin(x)
            y_range = np.nanmax(y) - np.nanmin(y)
            x_padding = x_range * FigMisc.ZBD_AXIS_PADDING
            y_padding = y_range * FigMisc.ZBD_AXIS_PADDING
            ax.set_xlim([np.nanmin(x) - x_padding, np.nanmax(x) + x_padding])
            ax.set_ylim([np.nanmin(y) - y_padding, np.nanmax(y) + y_padding])

        # Get ocean composition data for composition lines
        ocean_comp = result.base.oceanComp.flatten()
        
        
        # Separate NaN and non-NaN points for different plotting
        nan_mask = np.isnan(y)  # Check for NaN in y variable (the specified property)
        valid_mask = ~nan_mask
        
        # Plot non-NaN points colored by ocean thickness (z variable)
        if np.any(valid_mask):
            z_valid = z[valid_mask]
            pts = ax.scatter(x[valid_mask], y[valid_mask], c=z_valid,
                            cmap=FigMisc.ZBD_COLORMAP, marker=Style.MS_Induction, 
                            s=Style.MW_Induction**2, edgecolors=FigMisc.ZBD_DOT_EDGE_COLOR,
                            linewidths=FigMisc.ZBD_DOT_EDGE_WIDTH, zorder=3, 
                            vmin=np.nanmin(z_valid), vmax=np.nanmax(z_valid))
                        # Draw composition lines if enabled using helper
            if FigMisc.DRAW_COMPOSITION_LINE:
                draw_ocean_composition_lines(ax, x[valid_mask], y[valid_mask], z[valid_mask], ocean_comp[valid_mask], 
                                        FigMisc.MANUAL_HYDRO_COLORS,
                                        FigMisc.ZBD_COMP_LINE_WIDTH, FigMisc.ZBD_COMP_LINE_ALPHA)
            
            # Add colorbar
            cbar = fig.colorbar(pts, ax=ax, format=FigLbl.cbarFmt)
            # Set colorbar ticks to match every unique ice shell thickness value
            if np.size(z_valid) > 0:
                # Get unique ice shell thickness values (x contains the ice shell thicknesses)
                unique_ice_thicknesses = np.unique(z_valid)
                # Set ticks to every unique ice shell thickness
                cbar.set_ticks(unique_ice_thicknesses)
                pts.set_clim(vmin=np.min(z_valid), vmax=np.max(z_valid))
            cbar.set_label(FigLbl.cbarLabelExplore, size=12)
        
   
        
        # Add legends if enabled
        if FigMisc.DRAW_COMPOSITION_LINE and Params.LEGEND:
            ax.legend(title="Ocean Composition", fontsize=FigMisc.ZBD_COMP_LEGEND_FONT_SIZE, 
                     title_fontsize=FigMisc.ZBD_COMP_LEGEND_TITLE_SIZE)
        elif np.any(nan_mask) and Params.LEGEND:
            # Show legend for invalid data points if no composition lines are shown
            ax.legend(fontsize=FigMisc.ZBD_COMP_LEGEND_FONT_SIZE)
        
        # Save the plot
        plt.tight_layout()
        fig.savefig(Params.FigureFiles.exploreZbY, format=FigMisc.figFormat,
                   dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Plot saved to file: {Params.FigureFiles.exploreZbY}')
        plt.close()


def PlotExploreOgramLoveComparison(results_list, Params):
    """
    
    Plots scatter showing tidal love number k2 versus delta (1+k2-h2) for comparison 
    against canonical k2/delta exploration plots. Features ice thickness coloring,
    composition lines, dual legend system, and optional error bars.
    
    Args:
        results_list: List of ExplorationResults objects (individual + optional comparison)
        Params: Configuration parameters
    """
    
    # Set up basic figure labels using first result
    first_result = results_list[0]
    FigLbl.SetExploration(first_result.base.bodyname, 'deltaLoveAmp', 
                          'kLoveAmp', 'oceanComp')
    if not FigMisc.TEX_INSTALLED:
        FigLbl.StripLatex()
    
    # Plot each result (individual exploration results + optional comparison)
    last_index = len(results_list) - 1
    for i, result in enumerate(results_list):
        
        # Skip results with no H2O
        if result.base.NO_H2O:
            continue
            
        # Set axis variable names for this plot type
        result.xName = 'deltaLoveAmp'
        result.yName = 'kLoveAmp'
        result.zName = 'oceanComp'
        
        # Detect if this is a comparison plot (last result when COMPARE=True)
        is_comparison_plot = Params.COMPARE and i == last_index
        
        # Extract and validate data using helper
        plot_data = extract_and_validate_plot_data(result_obj = result, x_field = result.xName, y_field = result.yName, c_field = result.zName,
                                                  x_multiplier = FigLbl.xMultExplore, y_multiplier = FigLbl.yMultExplore, c_multiplier = FigLbl.zMultExplore,
                                                  custom_x_axis = FigLbl.xCustomAxis, custom_y_axis = FigLbl.yCustomAxis)
        
        if len(plot_data['x']) == 0:
            log.warning(f"No valid data points for {result.base.bodyname}")
            continue
        
        # Create figure and axis
        fig = plt.figure(figsize=FigSize.explore)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)
        
        # Set title based on plot type
        if Params.TITLES:
            if is_comparison_plot:
                fig.suptitle(FigLbl.exploreCompareTitle)
            else:
                fig.suptitle(FigLbl.explorationLoveComparisonTitle)
        
        # Set up axis labels and scales
        ax.set_xlabel(FigLbl.xLabelExplore)
        ax.set_ylabel(FigLbl.yLabelExplore)
        ax.set_xscale('linear')  # Love plots use linear scales
        ax.set_yscale('linear')
        
        # Extract plot data
        x, y, z = plot_data['x'], plot_data['y'], plot_data['c']
        
        # Get ice thickness data for coloring
        ice_thickness_data = extract_and_validate_plot_data(result_obj = result, x_field = 'zb_approximate_km', y_field = result.yName,
                                                           c_field = None, contour_field = None,
                                                           x_multiplier = 1.0, y_multiplier = 1.0, c_multiplier = 1.0, contour_multiplier = 1.0,
                                                           custom_x_axis = FigLbl.xCustomAxis, custom_y_axis = FigLbl.yCustomAxis)
        ice_thickness = ice_thickness_data['x']  # x field contains the ice thickness
        
        # Get ocean composition data for composition lines
        ocean_comp = result.base.oceanComp.flatten()
        
        # Draw composition lines if enabled using helper
        if FigMisc.DRAW_COMPOSITION_LINE:
            # Get temperature data for color scaling
            temp_data = extract_and_validate_plot_data(result_obj = result, x_field = 'Tb_K', y_field = result.yName, c_field = None, contour_field = None,
                                                       x_multiplier = 1.0, y_multiplier = 1.0, c_multiplier = 1.0, contour_multiplier = 1.0,
                                                       custom_x_axis = FigLbl.xCustomAxis, custom_y_axis = FigLbl.yCustomAxis)
            temp_values = temp_data['x']  # Temperature values for color scaling
            
            # Set up temperature bounds for color mapping
            if len(temp_values) > 0:
                TminMax_K = [np.nanmin(temp_values), np.nanmax(temp_values)]
                Color.Tbounds_K = TminMax_K
            
            draw_ocean_composition_lines(ax, x, y, temp_values, ocean_comp, 
                                       FigMisc.MANUAL_HYDRO_COLORS,
                                       FigMisc.LOVE_COMP_LINE_WIDTH, FigMisc.LOVE_COMP_LINE_ALPHA)
            
            # Add error bars if enabled
            if FigMisc.SHOW_ERROR_BARS:
                ax.errorbar(x, y, yerr=FigMisc.ERROR_BAR_MAGNITUDE, fmt='none', 
                           color='gray', capsize=3, alpha=0.5)
        
        # Handle ice thickness coloring
        if len(ice_thickness) > 0:
            # Set bounds for normalization
            Tbound_lower = np.min(ice_thickness)
            Tbound_upper = np.max(ice_thickness)
            # Get normalized values using GetNormT
            norm_thickness = Color.GetNormT(ice_thickness, Tbound_lower, Tbound_upper)
        else:
            norm_thickness = np.ones(np.shape(x))
        
        # Plot scatter points with ice thickness coloring and optional convection-based markers
        if FigMisc.SHOW_CONVECTION_WITH_SHAPE:
            # Plot with different markers for convection vs non-convection
            # Convection points (Dconv_m > 0): use circle marker
            # Non-convection points (Dconv_m <= 0): use square marker
            convection_data = extract_and_validate_plot_data(result_obj = result, x_field = result.xName, y_field = result.yName,
                                                           c_field = 'Dconv_m', contour_field = None,
                                                           x_multiplier = 1.0, y_multiplier = 1.0, c_multiplier = 1.0, contour_multiplier = 1.0,
                                                           custom_x_axis = FigLbl.xCustomAxis, custom_y_axis = FigLbl.yCustomAxis)
            # Create boolean array: True for convection (Dconv_m > 0), False for no convection
            convection_markers = convection_data['c'] > 0
            conv_mask = convection_markers
            non_conv_mask = ~convection_markers
            # 2. Create a Normalize object to map your data range to [0, 1]
            norm = Normalize(vmin=norm_thickness.min(), vmax=norm_thickness.max())

            # 3. Apply the colormap to normalized data to get RGBA values
            colors = get_cmap(FigMisc.LOVE_ICE_THICKNESS_CMAP)(norm(norm_thickness)) 
            
        
            # Plot convection points
            if np.any(conv_mask):
                ax.scatter(x[conv_mask], y[conv_mask], c=colors[conv_mask], marker='v', 
                          s=Style.MW_Induction**2, edgecolors=FigMisc.LOVE_DOT_EDGE_COLOR, 
                          linewidths=FigMisc.LOVE_DOT_EDGE_WIDTH, zorder=3)
            
            # Plot non-convection points
            if np.any(non_conv_mask):
                ax.scatter(x[non_conv_mask], y[non_conv_mask], c=colors[non_conv_mask], marker='^', 
                          s=Style.MW_Induction**2, edgecolors=FigMisc.LOVE_DOT_EDGE_COLOR, 
                          linewidths=FigMisc.LOVE_DOT_EDGE_WIDTH, zorder=3)
        else:
            # Plot with default marker
            pts = ax.scatter(x, y, c=norm_thickness,
                            cmap=FigMisc.LOVE_ICE_THICKNESS_CMAP, marker=Style.MS_Induction, 
                            s=Style.MW_Induction**2, edgecolors=FigMisc.LOVE_DOT_EDGE_COLOR, 
                            linewidths=FigMisc.LOVE_DOT_EDGE_WIDTH, zorder=3)
        
        # Create ice thickness legend if enabled
        if FigMisc.SHOW_ICE_THICKNESS_DOTS and len(ice_thickness) > 0 and Params.LEGEND:
            ice_legend = create_ice_thickness_colorbar(
                ax, ice_thickness,
                cmap_name=FigMisc.LOVE_ICE_THICKNESS_CMAP,
            )
        
        # Set axis limits with special Love plot formatting
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
            legend1 = ax.legend(fontsize=FigMisc.LOVE_COMP_LEGEND_FONT_SIZE,
                               loc='upper right', bbox_to_anchor=(1.0, 1.0))
            ax.add_artist(legend1)
        
        # Add second legend for convection markers if enabled
        if FigMisc.SHOW_CONVECTION_WITH_SHAPE and Params.LEGEND:
            legend_elements = [
                Line2D([0], [0], marker='v', color='w', markerfacecolor='gray', 
                       markersize=8, label='Convection'),
                Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', 
                       markersize=8, label='No Convection')
            ]
            legend2 = ax.legend(handles=legend_elements, fontsize=FigMisc.LOVE_COMP_LEGEND_FONT_SIZE,
                               loc='lower left', bbox_to_anchor=(0.0, 0.0))
            ax.add_artist(legend2)
        # Save the plot
        plt.tight_layout()
        fig.savefig(Params.FigureFiles.exploreLoveComparison, format=FigMisc.figFormat,
                   dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.debug(f'Plot saved to file: {Params.FigureFiles.exploreLoveComparison}')
        plt.close()
    

def PlotImaginaryVersusReal(results_list, Params, ax=None):
    """
    Creates subplots for each excitation/peak with ice thickness-based zoom insets similar to 
    PlotMonteCarloScatter. Features composition-based coloring, marker shapes per excitation,
    and configurable ice thickness highlighting with tolerance-based selection.
    
    Args:
        results_list: List of result objects (MonteCarloResults, ExplorationResults, etc.)
        Params: Configuration parameters
        ax: Optional existing axis (if None, creates new figure with subplots)
        
    Configuration flags:
        FigMisc.HIGHLIGHT_ICE_THICKNESSES: Enable ice thickness highlighting
        FigMisc.DO_SCATTER_INSET: Enable zoom insets
        FigMisc.ICE_THICKNESSES_TO_SHOW: Target ice thicknesses to highlight
        FigMisc.ICE_THICKNESS_TOLERANCE: Tolerance for ice thickness matching
    """
    
    if not results_list:
        log.warning("No results provided for complex plotting")
        return
    
    
    # For now, only handle magnetic data types
    data_types_to_try = ['magnetic']
    edgeLineWidth = 1 #TODO Make variable
    for data_type in data_types_to_try:
        # Extract complex data from first result to determine data type and setup
        first_result = results_list[0]
        plot_data = extract_complex_plot_data(first_result, data_type=data_type, Params=Params)
        
        if not plot_data or not plot_data.get('complex_data'):
            continue  # Try next data type
        
        data_type = plot_data['data_type']
        componentsAvailable = plot_data['components'] 
        excitation_names = plot_data['excitation_names']
        n_peaks = len(excitation_names)
        
        log.info(f"Plotting {data_type} complex data with {len(componentsAvailable)} components and {n_peaks} peaks")

                    
        # Check if ice thickness highlighting is enabled
        highlight_ice_thickness = (FigMisc.HIGHLIGHT_ICE_THICKNESSES and 
                                   first_result.base.zb_km is not None)
                
        if FigMisc.EDGE_COLOR_K_IN_COMPLEX_PLOTS and data_type == 'magnetic':
            kData = first_result.base.kLoveAmp
            if np.all(np.isnan(kData)):
                EDGE_COLOR_WITH_K = False
            else:
                EDGE_COLOR_WITH_K = True
                cmap = plt.cm.viridis  # Use viridis colormap for edge colors
                k_norm, edgeColor, cBarTitle, cBarFmt, cTicksSpacings = normalizeDataForColor(kData, 'kLoveAmp', cmap)
                edgeColor = edgeColor.reshape(-1, 4)
        else:
            EDGE_COLOR_WITH_K = False
            
        if not EDGE_COLOR_WITH_K:
            # Create array of black colors matching Tb_K shape
            edgeColor = np.full(first_result.base.Tb_K.shape + (4,), [0, 0, 0, 1]).reshape(-1, 4)
            
        # Plot each component separately  
        for comp in componentsAvailable:
            n_cols = int(np.ceil(np.sqrt(n_peaks)))
            n_rows = int(np.ceil(n_peaks / n_cols))
            
            
            # Create figure with subplots
            if ax is None:
                fig, axes = plt.subplots(n_rows, n_cols, figsize=(FigSize.imaginaryRealSoloCombo[0] * n_cols, FigSize.imaginaryRealSoloCombo[1] * n_rows))
                # Ensure axes is always 2D array for consistent indexing
                if n_rows == 1 and n_cols == 1:
                    axes = np.array([[axes]])
                elif n_rows == 1:
                    axes = axes.reshape(1, -1)
                elif n_cols == 1:
                    axes = axes.reshape(-1, 1)
            else:
                fig = ax.figure
                axes = np.array([[ax]])
                n_rows, n_cols = 1, 1
            # Set main title
            if Params.TITLES:
                if len(set(getattr(result.base, 'bodyname', 'Unknown') for result in results_list)) == 1:
                    bodyname = getattr(results_list[0].base, 'bodyname', 'Unknown')
                    if highlight_ice_thickness:
                        fig.suptitle(f'{bodyname} {FigLbl.BdipTitle} - {comp} (Ice Thickness Highlighted)', fontsize=14)
                    else:
                        fig.suptitle(f'{bodyname} {FigLbl.BdipTitle} - {comp}', fontsize=14)
                else:
                    if highlight_ice_thickness:
                        fig.suptitle(f'{FigLbl.BdipCompareTitle} - {comp} (Ice Thickness Highlighted)', fontsize=14)
                    else:
                        fig.suptitle(f'{FigLbl.BdipCompareTitle} - {comp}', fontsize=14)
            
            # Create subplot for each excitation/peak
            for peak_idx, exc_name in enumerate(excitation_names):
                row, col = peak_idx // n_cols, peak_idx % n_cols
                ax_subplot = axes[row, col]
                
                if Style.GRIDS:
                    ax_subplot.grid()
                    ax_subplot.set_axisbelow(True)
                
                # Set axis labels
                if data_type == 'magnetic':
                    if col == 0:
                        ax_subplot.set_ylabel(FigLbl.BdipImLabel[comp])
                    if row == n_rows-1:
                        ax_subplot.set_xlabel(FigLbl.BdipReLabel[comp])
                
                # Set subplot title
                if data_type == 'magnetic':
                    ax_subplot.set_title(f'{exc_name.capitalize()}')
                else:
                    ax_subplot.set_title(f'{exc_name.capitalize()} - {comp}')
                
                # Collect all data for this peak across all results
                all_real_data = np.array([])
                all_imag_data = np.array([])
                all_ice_thickness = np.array([])
                all_ocean_comps = np.array([])
                all_marker_sizes = np.array([])
                legendElements = []
                
                for i_result, result in enumerate(results_list):
                    plot_data = extract_complex_plot_data(result, data_type=data_type, Params=Params)
                    complex_data = plot_data['complex_data']
                    
                    # Get ice shell thickness for this result
                    iceShellThickness, zb_name = getIceShellThickness(result)
                    
                    # Get composition data
                    oceanComps = result.base.oceanComp
                    
                    # Normalize ice shell thickness and convert to marker sizes
                    if np.max(iceShellThickness) - np.min(iceShellThickness) > FigMisc.ICE_THICKNESS_TOLERANCE:
                        iceShellThicknessNormalized = (iceShellThickness - np.min(iceShellThickness)) / (np.max(iceShellThickness) - np.min(iceShellThickness))
                        min_size = 10
                        max_size = 90
                        marker_sizes = min_size + iceShellThicknessNormalized * (max_size - min_size)
                    else:
                        marker_sizes = np.full_like(iceShellThickness, fill_value=50)
                    
                    # Extract data for this peak
                    excData = complex_data[comp][peak_idx]
                    real_val = np.real(excData)
                    imag_val = np.imag(excData)
                    
                    # Accumulate data
                    all_real_data = np.append(all_real_data, real_val)
                    all_imag_data = np.append(all_imag_data, imag_val)
                    all_ice_thickness = np.append(all_ice_thickness, iceShellThickness)
                    all_ocean_comps = np.append(all_ocean_comps, oceanComps)
                    all_marker_sizes = np.append(all_marker_sizes, marker_sizes)
                
                # Determine ice thickness highlighting mask if enabled
                highlight_mask = None
                if highlight_ice_thickness:
                    highlight_mask = np.zeros(len(all_real_data), dtype=bool)
                    for target_thickness in FigMisc.ICE_THICKNESSES_TO_SHOW:
                        thickness_mask = np.abs(all_ice_thickness - target_thickness) <= FigMisc.ICE_THICKNESS_TOLERANCE
                        highlight_mask |= thickness_mask
                    
                    if EDGE_COLOR_WITH_K:
                        k_norm, edgeColor, cBarTitle, cBarFmt, cTicksSpacings = normalizeDataForColor(kData, 'kLoveAmp', cmap, highlight_mask)
                # Group by unique ocean compositions (maintaining order) and plot connecting lines
                uniqueOceanComps, indices = np.unique(all_ocean_comps, return_index=True)
                uniqueOceanComps = uniqueOceanComps[np.argsort(indices)]
                
                # Track plotted dimmed coordinates to prevent overlapping brightness
                plotted_dimmed_coords = set()
                
                for comp_idx, ocean_comp in enumerate(uniqueOceanComps):
                    oceanCompLabel = format_composition_label(ocean_comp)
                    comp_mask = all_ocean_comps == ocean_comp
                    
                    if not np.any(comp_mask):
                        continue
                    
                    color = Color.cmap[ocean_comp](0.5)
                    markerShape = Style.MS_dip[exc_name] if data_type == 'magnetic' else 'o'
                    
                    # Get data for this composition
                    comp_real = all_real_data[comp_mask]
                    comp_imag = all_imag_data[comp_mask]
                    comp_sizes = all_marker_sizes[comp_mask]
                    comp_edge_color = edgeColor[comp_mask]
                    
                    
                    if highlight_ice_thickness:
                        # Plot highlighted points (full alpha)
                        highlight_comp_mask = highlight_mask[comp_mask]
                        highlight_edge_color = comp_edge_color[highlight_comp_mask]
                        if np.any(highlight_comp_mask):
                            ax_subplot.scatter(comp_real[highlight_comp_mask], comp_imag[highlight_comp_mask],
                                             facecolor=color, edgecolor=highlight_edge_color, linewidths=edgeLineWidth,
                                             s=comp_sizes[highlight_comp_mask], marker=markerShape, zorder=3,
                                             alpha=1.0, label=f'{oceanCompLabel}' if comp_idx == 0 else "")
                        
                        # Plot dimmed points (lower alpha) - but only if not already plotted
                        dimmed_comp_mask = ~highlight_comp_mask
                        if np.any(dimmed_comp_mask):
                            # Filter out coordinates that have already been plotted
                            dimmed_real = comp_real[dimmed_comp_mask]
                            dimmed_imag = comp_imag[dimmed_comp_mask]
                            dimmed_sizes = comp_sizes[dimmed_comp_mask]
                            
                            # Check which points haven't been plotted yet
                            # This is a hack to prevent overlapping points from being plotted and increasing the alpha of the dimmed points
                            unplotted_mask = np.array([
                                (round(x, 1), round(y, 1)) not in plotted_dimmed_coords 
                                for x, y in zip(dimmed_real, dimmed_imag)
                            ])
                            
                            if np.any(unplotted_mask):
                                # Plot only unplotted points
                                unplotted_real = dimmed_real[unplotted_mask]
                                unplotted_imag = dimmed_imag[unplotted_mask]
                                unplotted_sizes = dimmed_sizes[unplotted_mask]
                                
                                ax_subplot.scatter(unplotted_real, unplotted_imag,
                                                 facecolor=color, edgecolor='none',
                                                 s=unplotted_sizes, marker=markerShape, zorder=1,
                                                 alpha=0.2, label="")
                                
                                # Add these coordinates to the plotted set
                                for x, y in zip(unplotted_real, unplotted_imag):
                                    plotted_dimmed_coords.add((round(x, 1), round(y, 1)))

                    else:
                        # Standard plotting without highlighting
                        ax_subplot.scatter(comp_real, comp_imag,
                                         facecolor=color, edgecolor=comp_edge_color, linewidths=edgeLineWidth,
                                         s=comp_sizes, marker=markerShape,
                                         alpha=1.0, label=f'{oceanCompLabel}')
                    
                    # Add k2 value labels above points if enabled and few enough points
                    if EDGE_COLOR_WITH_K and len(np.unique(kData[~np.isnan(kData)])) <= 5:
                        for i, (x, y) in enumerate(zip(comp_real, comp_imag)):
                            # Find the k value for this point
                            original_idx = np.where((all_real_data == x) & (all_imag_data == y))[0]
                            if len(original_idx) > 0:
                                k_idx = original_idx[0]  # Take first match if multiple
                                if k_idx < len(kData) and not np.isnan(kData.flatten()[k_idx]):
                                    ax_subplot.annotate(f'{kData.flatten()[k_idx]:.3f}', 
                                                      (x, y), 
                                                      xytext=(0, 8), 
                                                      textcoords='offset points',
                                                      ha='center', va='bottom',
                                                      fontsize=8, 
                                                      alpha=0.8)
                    legendElements.append(plt.Line2D([0], [0], marker=markerShape, color='none', markerfacecolor=color, markeredgecolor='black', markeredgewidth=edgeLineWidth,
                                    markersize=np.sqrt(comp_sizes[0]), label=f'{oceanCompLabel}'))
                if not FigMisc.DO_SCATTER_INSET:
                    DO_INSET = False
                elif not highlight_ice_thickness or not np.any(highlight_mask):
                    DO_INSET = False
                else:
                    DO_INSET = True
                
                # Create zoom inset if highlighting enabled and there are highlighted points
                inset_ax = None
                if DO_INSET:
                    inset_ax, optimal_location = create_zoom_inset(ax_subplot, all_real_data, all_imag_data, 
                                                                  highlight_mask, add_visual_indicators=True)
                    
                    if inset_ax is not None:
                        # Plot only highlighted points in inset
                        for comp_idx, ocean_comp in enumerate(uniqueOceanComps):
                            comp_mask = all_ocean_comps == ocean_comp
                            highlight_comp_mask = highlight_mask & comp_mask
                            
                            if np.any(highlight_comp_mask):
                                color = Color.cmap[ocean_comp](0.5)
                                markerShape = Style.MS_dip[exc_name] if data_type == 'magnetic' else 'o'
                                
                                inset_ax.scatter(all_real_data[highlight_comp_mask], all_imag_data[highlight_comp_mask],
                                               facecolor=color, edgecolor=edgeColor[highlight_comp_mask],
                                               s=10, marker=markerShape,
                                               alpha=1.0)
                ax_subplot.set_xlim(left=0)
                ax_subplot.set_ylim(bottom=0)
                
                
                # Set up grid with configurable spacing
                axs_to_format = [ax_subplot]
                if inset_ax is not None:
                    axs_to_format.append(inset_ax)
                
                for ax_format in axs_to_format:
                    """points = ax_format.dataLim.get_points()
                    # returns [[xmin, ymin], [xmax, ymax]]
                    xmax, ymax = points[1]
                    xmin, _ = ax_format.get_xlim()
                    ymin, _ = ax_format.get_ylim()
                    ax_format.set(xlim=(xmin, xmax), ylim=(ymin, ymax))"""
                    minor_spacing = 3.0
                    
                    # Calculate major spacing as a multiple of minor spacing based on axis limits
                    x_range = ax_format.get_xlim()[1] - ax_format.get_xlim()[0]
                    y_range = ax_format.get_ylim()[1] - ax_format.get_ylim()[0]
                    
                    # Determine appropriate major spacing multiplier based on range
                    x_major_multiplier = max(1, int(x_range / (minor_spacing * 4)))
                    y_major_multiplier = max(1, int(y_range / (minor_spacing * 4)))
                    major_multiplier = min(x_major_multiplier, y_major_multiplier)
                    
                    x_major_spacing = minor_spacing * major_multiplier
                    y_major_spacing = minor_spacing * major_multiplier
                    # Round axis limits to the nearest ceiling major spacing
                    current_xlim = ax_format.get_xlim()
                    current_ylim = ax_format.get_ylim()
                    
                    # Calculate ceiling values rounded to major spacing
                    new_xmax = np.ceil(current_xlim[1] / minor_spacing) * minor_spacing
                    new_ymax = np.ceil(current_ylim[1] / minor_spacing) * minor_spacing
                    
                    new_xlim = (current_xlim[0], new_xmax)
                    new_ylim = (current_ylim[0], new_ymax)
                    
                    ax_format.set_xlim(new_xlim)
                    ax_format.set_ylim(new_ylim)
                    
                    ax_format.xaxis.set_major_locator(ticker.MultipleLocator(x_major_spacing))
                    ax_format.yaxis.set_major_locator(ticker.MultipleLocator(y_major_spacing))
                    ax_format.xaxis.set_minor_locator(ticker.MultipleLocator(minor_spacing))
                    ax_format.yaxis.set_minor_locator(ticker.MultipleLocator(minor_spacing))
                    ax_format.grid(which='minor', alpha=0.3)
                    ax_format.grid(which='major', alpha=0.5)
            
            
            # Hide unused subplots
            for i in range(n_peaks, n_rows * n_cols):
                row, col = i // n_cols, i % n_cols
                axes[row, col].set_visible(False)

            if EDGE_COLOR_WITH_K:
                fig.colorbar(plt.cm.ScalarMappable(norm=k_norm, cmap=cmap), ax=axes[0, n_cols-1], label=cBarTitle, format=cBarFmt, ticks=cTicksSpacings)
            
             # Add marker size legend to show ice shell thickness relationship
            if Params.LEGEND:
                add_composition_legend(axes[0, n_cols-1], uniqueOceanComps, show_ice_thickness_legend=False, overrideLegend = legendElements)
                # Get unique ice thickness and marker size pairs that actually exist in the data
                unique_pairs = []
                for thickness, size in zip(all_ice_thickness, all_marker_sizes):
                    pair = (round(thickness, 1), round(size, 1))
                    if pair not in unique_pairs:
                        unique_pairs.append(pair)
                if len(unique_pairs) > 1:
                    # Sort by ice thickness for consistent ordering
                    unique_pairs.sort(key=lambda x: x[0])
                    
                    # Limit to max 5 legend entries to avoid overcrowding
                    if len(unique_pairs) > 5:
                        # Select evenly spaced entries including first and last
                        indices = np.linspace(0, len(unique_pairs)-1, 5, dtype=int)
                        unique_pairs = [unique_pairs[i] for i in indices]
                    
                    legend_elements = []
                    for thickness, size in unique_pairs:
                        legend_elements.append(
                            plt.Line2D([0], [0], marker=markerShape, color='w', markerfacecolor='gray', markeredgecolor='black', markeredgewidth=np.sqrt(3),
                                        markersize=np.sqrt(size), label=f'{thickness:.1f} {FigLbl.distanceUnits}')
                        )
                    
                    axes[n_rows-1, n_cols-1].legend(handles=legend_elements, fontsize='small',
                                    title=FigLbl.zbLabel, loc= 'center left')
            # Save the plot if we created a new figure
            if ax is None:
                plt.tight_layout()
                if data_type == 'magnetic':
                    save_path = Params.FigureFiles.BdipExplore[comp]
                else:
                    # For Love numbers, create file path
                    save_path = getattr(Params.FigureFiles, 'Love', {}).get(comp, f'love_numbers_{comp}.png')
                
                fig.savefig(save_path, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
                log.debug(f'Complex {data_type} {comp} subplot plot saved to file: {save_path}')
                plt.close()
    
    return