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
from PlanetProfile.MagneticInduction.Moments import Excitations as Mag
from MoonMag.asymmetry_funcs import getMagSurf as GetMagSurf

log = logging.getLogger('PlanetProfile')

def PlotMonteCarloResults(results_list, Params):
    """  
    Plots histograms of output parameters from Monte Carlo or exploration analysis.
    Uses user-defined output parameters and works with both Monte Carlo and exploration structures.
    
    Args:
        results_list: List of MonteCarloResults or ExplorationResults objects 
        Params: Configuration parameters with optional outputParams attribute
    """
    
    # Define default output parameters - user can override via Params.outputParams
    default_output_params = ['CMR2calc', 'Mtot_kg', 'kLoveAmp', 'hLoveAmp', 'zb_km']
    output_params = getattr(Params, 'outputParams', default_output_params)
    
    # Process each result
    for result in results_list:
        
        # Get valid mask - assume it exists if needed
        valid_mask = result.base.VALID
        
        # Check if we should plot by ocean composition
        plot_by_ocean_comp = (result.base.oceanComp is not None and 
                              Params.MonteCarlo.plotOceanComps)
        
        # Set up subplot grid
        n_params = len(output_params)
        n_cols = 4
        n_rows = int(np.ceil(n_params / n_cols))
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4*n_rows))
        if n_rows == 1:
            axes = axes.reshape(1, -1)
        
        fig.suptitle(f'Output Parameter Distributions - {result.base.bodyname}', fontsize=16)
        
        # Ocean composition setup if needed
        if plot_by_ocean_comp:
            valid_ocean_comps = result.base.oceanComp[valid_mask]
            unique_comps = np.unique(valid_ocean_comps)
            
            # Set up color mapping
            ocean_comp_list = result.base.oceanComp
            if isinstance(ocean_comp_list, np.ndarray):
                ocean_comp_list = ocean_comp_list.tolist()
            
            cmap = plt.cm.coolwarm
        
        # Plot each output parameter
        for i, param_name in enumerate(output_params):
            row, col = i // n_cols, i % n_cols
            ax = axes[row, col]
            
            # Extract parameter values from base structure
            values = getattr(result.base, param_name)
            
            # Apply valid mask
            valid_values = values[valid_mask]
            
            if len(valid_values) > 0 and np.any(np.isfinite(valid_values)):
                if plot_by_ocean_comp:
                    # Plot as jittered dots grouped by ocean composition
                    all_valid_values = valid_values[~np.isnan(valid_values)]
                    _, bin_edges = np.histogram(all_valid_values, bins=20)
                    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                    
                    for j, comp in enumerate(unique_comps):
                        comp_mask = valid_ocean_comps == comp
                        comp_values = valid_values[comp_mask]
                        
                        if len(comp_values) > 0:
                            comp_values_clean = comp_values[~np.isnan(comp_values)]
                            hist, _ = np.histogram(comp_values_clean, bins=bin_edges)
                            
                            # Set color based on composition order
                            try:
                                comp_index = ocean_comp_list.index(comp)
                                norm_val = comp_index / (len(ocean_comp_list) - 1) if len(ocean_comp_list) > 1 else 0.5
                                this_color = cmap(norm_val)
                            except ValueError:
                                norm_val = j / (len(unique_comps) - 1) if len(unique_comps) > 1 else 0.5
                                this_color = cmap(norm_val)
                            
                            # Plot jittered dots for each bin
                            for k, (center, count) in enumerate(zip(bin_centers, hist)):
                                if count > 0:
                                    # Handle composition label formatting
                                    comp_label = format_composition_label(comp)
                                    
                                    # Create jittered scatter points
                                    x_jittered = center + np.random.normal(0, (bin_edges[1] - bin_edges[0]) * 0.05, count)
                                    y_jittered = np.ones(count) * count + np.random.normal(0, 0.15, count)
                                    
                                    ax.scatter(x_jittered, y_jittered, 
                                             color=this_color, alpha=0.7, s=25,
                                             label=f'{comp_label}' if k == 0 else "",
                                             marker='o')
                    
                    ax.set_ylabel('Count')
                    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
                    ax.set_title(f'{param_name} by Ocean Composition')
                    
                    # Add legend only to first subplot
                    if i == 0 and Params.LEGEND:
                        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
                
                else:
                    # Plot as regular histogram
                    ax.hist(valid_values, bins=30, alpha=0.7, density=True, edgecolor='black')
                    ax.set_ylabel('Density')
                    ax.set_title(f'{param_name} Distribution')
                
                ax.set_xlabel(param_name)
                ax.grid(True, alpha=0.3)
            else:
                # No valid data for this parameter
                ax.text(0.5, 0.5, f'No valid data\nfor {param_name}', 
                       transform=ax.transAxes, ha='center', va='center')
                ax.set_xlabel(param_name)
                ax.set_ylabel('Density' if not plot_by_ocean_comp else 'Count')
                ax.set_title(f'{param_name} Distribution')
        
        # Hide unused subplots
        for i in range(n_params, n_rows * n_cols):
            row, col = i // n_cols, i % n_cols
            axes[row, col].set_visible(False)
        
        # Save the plot
        plt.tight_layout()
        
        # Save to Monte Carlo results file - works for both result types
        fig.savefig(Params.FigureFiles.montecarloResults, format=FigMisc.figFormat, 
                   dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.info(f'Output parameter distributions plot saved to: {Params.FigureFiles.montecarloResults}')
        
        plt.close()


def PlotMonteCarloDistributions(results_list, Params):
    """
    Modern implementation of PlotMonteCarloDistributions using streamlined architecture.
    
    Plots histograms of parameter distributions from Monte Carlo analysis.
    Uses the same single-loop architecture as exploration functions.
    
    Args:
        results_list: List of MonteCarloResults objects 
        Params: Configuration parameters
    """
    
    # Process each Monte Carlo result
    for mc_result in results_list:
        
        # Skip if no successful models
        if mc_result.base.nSuccess == 0:
            log.warning(f'No successful models to plot for {mc_result.base.bodyname}')
            continue
        
        # Get parameter information
        n_params = len(mc_result.base.paramsToSearch)
        n_cols = 4
        n_rows = int(np.ceil(n_params / n_cols))
        
        # Create figure and subplots
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4*n_rows))
        if n_rows == 1:
            axes = axes.reshape(1, -1)
        
        fig.suptitle(f'Monte Carlo Parameter Distributions - {mc_result.base.bodyname}', fontsize=16)
        
        # Plot each parameter distribution
        for i, param in enumerate(mc_result.base.paramsToSearch):
            row, col = i // n_cols, i % n_cols
            ax = axes[row, col]
            
            # Extract parameter values using hierarchical structure
            values = mc_result.base.paramValues[param]
            valid_values = values[mc_result.base.VALID]
            
            # Plot all samples with histogram
            ax.hist(values, bins=30, alpha=0.5, label='All samples', density=True)
            
            # Add styling and labels
            ax.set_xlabel(param)
            ax.set_ylabel('Density')
            ax.set_title(f'{param} Distribution')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        # Hide unused subplots
        for i in range(n_params, n_rows * n_cols):
            row, col = i // n_cols, i % n_cols
            axes[row, col].set_visible(False)
        
        # Save the plot
        plt.tight_layout()
        fig.savefig(Params.FigureFiles.montecarloDistributions, format=FigMisc.figFormat, 
                   dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.info(f'Monte Carlo distributions plot saved to: {Params.FigureFiles.montecarloDistributions}')
        plt.close()


def PlotMonteCarloScatter(results_list, Params):
    """
    Modern implementation of PlotMonteCarloScatter using streamlined architecture.
    
    Plots scatter plots of Monte Carlo parameter pairs with optional ocean composition coloring
    and ice thickness highlighting. Features zoom inset when highlighting ice thickness.
    
    Args:
        results_list: List of MonteCarloResults objects 
        Params: Configuration parameters
    """
    # Process each Monte Carlo result
    for result in results_list:
        
        # Skip if no scatter parameters specified
        if not Params.MonteCarlo.scatterParams:
            log.warning('No scatter parameters specified in Params.MonteCarlo.scatterParams.')
            continue
        
        # Get valid mask
        if Params.ALLOW_BROKEN_MODELS:
            log.warning('Params.ALLOW_BROKEN_MODELS is True. Will plot all models, not just successful ones.')
            valid_mask = np.ones_like(result.base.CMR2calc, dtype=bool)
        else:
            valid_mask = result.base.VALID
        
        # Check if we should color by ocean composition
        color_by_ocean_comp = (result.base.oceanComp is not None and 
                             Params.MonteCarlo.plotOceanComps)
        
        # Get available excitations using the helper function
        available_excitations = get_available_excitations(result)
        
        # Process scatter parameter pairs
        scatter_pairs = Params.MonteCarlo.scatterParams
        if not isinstance(scatter_pairs[0], list):
            scatter_pairs = [scatter_pairs]
        
        # Build plot_info with proper excitation expansion
        plot_info = []
        for x_param, y_param in scatter_pairs:
            # Check if parameters are magnetic and need excitation expansion
            x_is_magnetic = x_param in ['Bix_nT', 'Biy_nT', 'Biz_nT', 'Bi_nT', 'phase']
            y_is_magnetic = y_param in ['Bix_nT', 'Biy_nT', 'Biz_nT', 'Bi_nT', 'phase']
            
            if (x_is_magnetic or y_is_magnetic) and available_excitations:
                # Create plots for each available excitation
                for exc_name in available_excitations:
                    x_exc = exc_name if x_is_magnetic else None
                    y_exc = exc_name if y_is_magnetic else None
                    plot_info.append((x_param, y_param, x_exc, y_exc))
            else:
                # Non-magnetic parameters or no excitations available
                plot_info.append((x_param, y_param, None, None))
        
        total_plots = len(plot_info)
        
        if total_plots == 0:
            log.warning('No valid parameter combinations found for scatter plots.')
            continue
        
        # Set up subplot grid
        n_cols = min(2, total_plots)
        n_rows = int(np.ceil(total_plots / n_cols))
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
        
        # Ensure axes is always 2D array for consistent indexing
        if n_rows == 1 and n_cols == 1:
            axes = np.array([[axes]])
        elif n_rows == 1:
            axes = axes.reshape(1, -1)
        elif n_cols == 1:
            axes = axes.reshape(-1, 1)
        
        # Set main title
        if color_by_ocean_comp:
            fig.suptitle(f'Monte Carlo Parameter Scatter Plots by Ocean Composition - {result.base.bodyname}', 
                        fontsize=14)
        else:
            fig.suptitle(f'Monte Carlo Parameter Scatter Plots - {result.base.bodyname}', fontsize=14)
        
        # Create scatter plots
        for plot_idx, (x_param, y_param, x_exc, y_exc) in enumerate(plot_info):
            row, col = plot_idx // n_cols, plot_idx % n_cols
            ax = axes[row, col]
            
            # Extract parameter values using helper (pass full result for hierarchical structure)
            x_values = extract_monte_carlo_parameter_values(result, x_param, x_exc, valid_mask)
            y_values = extract_monte_carlo_parameter_values(result, y_param, y_exc, valid_mask)
            
            # Apply valid mask and remove NaN values
            valid_data_mask = valid_mask & np.isfinite(x_values) & np.isfinite(y_values)
            x_plot = x_values[valid_data_mask]
            y_plot = y_values[valid_data_mask]
            
            if len(x_plot) == 0:
                ax.text(0.5, 0.5, 'No valid data', transform=ax.transAxes, ha='center', va='center')
                ax.set_xlabel(x_param)
                ax.set_ylabel(y_param)
                continue
            
            # Check if we should highlight specific ice thicknesses
            highlight_ice_thickness = (FigMisc.HIGHLIGHT_ICE_THICKNESSES and 
                                     result.base.zb_km is not None)
            
            highlight_mask = None
            if highlight_ice_thickness:
                # Get ice thickness data for valid points
                ice_thickness_plot = result.base.zb_km[valid_data_mask]
                
                # Determine which points to highlight
                highlight_mask = np.zeros(len(x_plot), dtype=bool)
                for target_thickness in FigMisc.ICE_THICKNESSES_TO_SHOW:
                    thickness_mask = np.abs(ice_thickness_plot - target_thickness) <= FigMisc.ICE_THICKNESS_TOLERANCE
                    highlight_mask |= thickness_mask
                
                dimmed_mask = ~highlight_mask
            
            # Plot data based on coloring mode
            if color_by_ocean_comp:
                # Group by unique ocean compositions (maintaining order) and plot connecting lines
                uniqueComps, indices = np.unique(all_ocean_comps, return_index=True)
                uniqueComps = uniqueComps[np.argsort(indices)]
                
                # Track plotted dimmed coordinates to prevent overlapping brightness
                plotted_dimmed_coords = set()
                
                for comp in uniqueComps:
                    comp_mask = result.base.oceanComp[valid_data_mask] == comp
                    if not np.any(comp_mask):
                        continue
                    
                    comp_label = format_composition_label(comp)
                    color = Color.cmap[comp](0.5)
                    
                    if highlight_ice_thickness:
                        # Plot highlighted points (higher alpha)
                        highlight_comp_mask = comp_mask & highlight_mask
                        if np.any(highlight_comp_mask):
                            ax.scatter(x_plot[highlight_comp_mask], y_plot[highlight_comp_mask], 
                                     color=color, alpha=0.7, s=FigMisc.SCATTER_DOT_SIZE, 
                                     label=comp_label, edgecolors='black', linewidth=0.5, zorder=3)
                        
                        # Plot dimmed points (lower alpha) - but only if not already plotted
                        dimmed_comp_mask = comp_mask & dimmed_mask
                        if np.any(dimmed_comp_mask):
                            # Filter out coordinates that have already been plotted
                            dimmed_x = x_plot[dimmed_comp_mask]
                            dimmed_y = y_plot[dimmed_comp_mask]
                            
                            # Check which points haven't been plotted yet
                            unplotted_mask = np.array([
                                (round(x, 6), round(y, 6)) not in plotted_dimmed_coords 
                                for x, y in zip(dimmed_x, dimmed_y)
                            ])
                            
                            if np.any(unplotted_mask):
                                # Plot only unplotted points
                                unplotted_x = dimmed_x[unplotted_mask]
                                unplotted_y = dimmed_y[unplotted_mask]
                                
                                ax.scatter(unplotted_x, unplotted_y, 
                                         color=color, alpha=FigMisc.DIMMED_ALPHA, s=FigMisc.SCATTER_DOT_SIZE, 
                                         edgecolors='black', linewidth=0.5, zorder=1)
                                
                                # Add these coordinates to the plotted set
                                for x, y in zip(unplotted_x, unplotted_y):
                                    plotted_dimmed_coords.add((round(x, 6), round(y, 6)))
                    else:
                        ax.scatter(x_plot[comp_mask], y_plot[comp_mask], 
                                 color=color, alpha=0.7, s=FigMisc.SCATTER_DOT_SIZE, 
                                 label=comp_label, edgecolors='black', linewidth=0.5)
                
                # Add legend only to first subplot
                if plot_idx == 0 and Params.LEGEND:
                    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
            
            else:
                # Standard scatter plot
                if highlight_ice_thickness:
                    # Plot highlighted points (higher alpha)
                    if np.any(highlight_mask):
                        ax.scatter(x_plot[highlight_mask], y_plot[highlight_mask], 
                                 alpha=0.7, s=FigMisc.SCATTER_DOT_SIZE, edgecolors='black', 
                                 linewidth=0.5, zorder=3,
                                 label=f'Ice thickness: {FigMisc.ICE_THICKNESSES_TO_SHOW} ± {FigMisc.ICE_THICKNESS_TOLERANCE} km' if plot_idx == 0 else None)
                    
                    # Plot dimmed points (lower alpha)
                    if np.any(dimmed_mask):
                        ax.scatter(x_plot[dimmed_mask], y_plot[dimmed_mask], 
                                 alpha=FigMisc.DIMMED_ALPHA, s=FigMisc.SCATTER_DOT_SIZE, 
                                 edgecolors='black', linewidth=0.5, zorder=1,
                                 label='Other ice thicknesses' if plot_idx == 0 else None)
                    
                    # Add legend for ice thickness highlighting
                    if plot_idx == 0 and Params.LEGEND:
                        ax.legend(fontsize='small')
                
                else:
                    ax.scatter(x_plot, y_plot, alpha=0.7, s=FigMisc.SCATTER_DOT_SIZE, 
                             edgecolors='black', linewidth=0.5)
            
            # Create zoom inset if highlighting ice thickness
            if highlight_ice_thickness and np.any(highlight_mask) and FigMisc.DO_SCATTER_INSET:
                inset_ax, optimal_location = create_zoom_inset(ax, x_plot, y_plot, highlight_mask)
                
                if inset_ax is not None:
                    # Plot the same data in the inset
                    if color_by_ocean_comp:
                        comp_plot = result.base.oceanComp[valid_data_mask]
                        for comp in uniqueComps:
                            comp_mask = comp_plot == comp
                            if np.any(comp_mask):
                                color = Color.cmap[comp](0.5)
                                
                                # Only plot highlighted points in zoom
                                highlight_comp_mask = comp_mask & highlight_mask
                                if np.any(highlight_comp_mask):
                                    inset_ax.scatter(x_plot[highlight_comp_mask], y_plot[highlight_comp_mask], 
                                                   color=color, alpha=0.7, s=FigMisc.SCATTER_DOT_SIZE*0.7, 
                                                   edgecolors='black', linewidth=0.3)
                    else:
                        inset_ax.scatter(x_plot[highlight_mask], y_plot[highlight_mask], 
                                       alpha=0.7, s=FigMisc.SCATTER_DOT_SIZE*0.7, 
                                       edgecolors='black', linewidth=0.3)
            else:
                inset_ax = None
            
            # Set axis labels
            if x_param in ['Bix_nT', 'Biy_nT', 'Biz_nT', 'Bi_nT', 'phase']:
                ax.set_xlabel(x_param + '_' + x_exc)
            else:
                ax.set_xlabel(x_param)
            if y_param in ['Bix_nT', 'Biy_nT', 'Biz_nT', 'Bi_nT', 'phase']:
                ax.set_ylabel(y_param + '_' + y_exc)
            else:
                ax.set_ylabel(y_param)


 
            
            axList = [ax] if inset_ax is None else [ax, inset_ax]
            for a in axList:
                # Set up grid with configurable spacing
                major_x_spacing = 0.12
                major_y_spacing = 12
                minor_x_spacing = 0.04  
                minor_y_spacing = 3
                if a.get_ylim()[1] - a.get_ylim()[0] < major_y_spacing:
                    major_y_spacing = minor_y_spacing * 2
                if a.get_xlim()[1] - a.get_xlim()[0] < major_x_spacing:
                    major_x_spacing = minor_x_spacing * 2
                a.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
                a.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
                a.xaxis.set_minor_locator(ticker.MultipleLocator(minor_x_spacing))
                a.yaxis.set_minor_locator(ticker.MultipleLocator(minor_y_spacing))
                a.grid(which = 'minor', alpha = 0.3)
                a.grid(which = 'major', alpha = 0.5)
        
        # Hide unused subplots
        if total_plots < n_rows * n_cols:
            for i in range(total_plots, n_rows * n_cols):
                row, col = i // n_cols, i % n_cols
                axes[row, col].set_visible(False)
        
        plt.tight_layout()
        
        # Save plot
        scatter_file_name = Params.FigureFiles.montecarloResults.replace('_results', '_scatter')
        fig.savefig(scatter_file_name, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
        log.info(f'Monte Carlo scatter plots saved to: {scatter_file_name}')
        plt.close()