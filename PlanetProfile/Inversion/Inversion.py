from PlanetProfile.Utilities.defineStructs import InversionParamsStruct
from PlanetProfile.Main import WriteProfile, ReloadProfile, ExploreOgram
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.gridspec import GridSpec
from PlanetProfile.Plotting.EssentialHelpers import *
import logging
logger = logging.getLogger('PlanetProfile')
def InvertBestPlanetList(BestPlanetList, Params, fNames):
    ExplorationList = []
    if len(fNames) == 1:
        Exploration, Params = ExploreOgram(BestPlanetList[0].bodyname, Params, fNameOverride=fNames[0])
        ExplorationList.append(Exploration)
    else:
        for fName in fNames:
            Exploration, Params = ExploreOgram(BestPlanetList[0].bodyname, Params, fNameOverride=fName)
            ExplorationList.append(Exploration)
    InvertBestPlanetMultiplot(BestPlanetList, ExplorationList, Params)

def InvertBestPlanet(BestPlanet, Params, fNames):
    """
    Invert for best-fit interior structure/constraints based on a set of input parameters.
    
    Now supports both single planet and list of planets:
    - If BestPlanet is a single planet object, creates single uncertainty plot
    - If BestPlanet is a list of planets, creates multiplot with all planets
    """
    Params.Inversion.assignRealPlanetModel(BestPlanet, BestPlanet.xBestFit, BestPlanet.yBestFit)
    ExplorationList = []
    if len(fNames) == 1:
        Exploration, Params = ExploreOgram(BestPlanet.bodyname, Params, fNameOverride=fNames[0])
        ExplorationList.append(Exploration)
    else:
        for fName in fNames:
            Exploration, Params = ExploreOgram(BestPlanet.bodyname, Params, fNameOverride=fName)
            ExplorationList.append(Exploration)
    ExplorationList = FitWithinUncertainty(ExplorationList, Params)
    
    PlotUncertainty(ExplorationList, Params)


def InvertBestPlanetMultiplot(BestPlanetList, Exploration, Params):
    """
    Create uncertainty plots for multiple best-fit planets in a single multiplot figure.
    
    Similar to PlotExploreOgramMultiSubplot, this function:
    1. Takes a list of best-fit planets instead of a single planet
    2. Creates exploration results for each planet
    3. Arranges them in a square-ish grid layout
    4. Uses shared legends and axes (only left column and bottom row labels)
    5. Saves only the combined multiplot
    
    Args:
        BestPlanetList: List of best-fit planet objects
        Params: Parameters object with inversion settings
        fName: Base filename for exploration results
    """
    
    n_planets = len(BestPlanetList)
    if n_planets == 0:
        log.warning("No planets provided for uncertainty multiplot")
        return
    
    # Calculate square-ish grid layout
    n_cols = 4
    n_rows = 2
    
    # Calculate figure size with scaling (same as PlotExploreOgramMultiSubplot)
    base_size = (6, 4)  # Same as PlotUncertainty default
    scale_factor = 1
    fig_width = base_size[0] * n_cols * scale_factor
    fig_height = base_size[1] * n_rows * scale_factor
    
    # Create figure with subplots
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    # Create subplots and plot uncertainty for each planet
    axes = []
    legend_elements = None  # Store legend elements from first plot
    
    for i, Planet in enumerate(BestPlanetList):
        # Fill by column: first column (0, 2, 4, 6), then second column (1, 3, 5, 7)
        col = i // n_rows
        row = i % n_rows
        subplot_index = row * n_cols + col + 1
        ax = fig.add_subplot(n_rows, n_cols, subplot_index)
        axes.append(ax)
        
        Params.Inversion.assignRealPlanetModel(Planet, Planet.xBestFit, Planet.yBestFit)
        
        # Call PlotUncertainty with this specific axis
        ax_result = PlotUncertaintySubplot(Exploration, Params, ax=ax, planet_name=Planet.bodyname)
        
        # Remove legend from all plots except the first one
        if i > 0 and ax.get_legend() is not None:
            ax.get_legend().remove()
        
        
        # Add subplot label (a, b, c, etc.) if enabled
        if hasattr(FigMisc, 'SUBPLOT_LABELS') and FigMisc.SUBPLOT_LABELS:
            label = f"{chr(ord('a') + i)}"
            if hasattr(FigMisc, 'SUBPLOT_LABEL_X'):
                label_x = FigMisc.SUBPLOT_LABEL_X
                label_y = FigMisc.SUBPLOT_LABEL_Y
                label_fontsize = FigMisc.SUBPLOT_LABEL_FONTSIZE
            else:
                label_x, label_y, label_fontsize = 0.02, 0.98, 12
            
            ax.text(label_x, label_y, label, 
                   transform=ax.transAxes, fontsize=label_fontsize,
                   fontweight='bold', ha='left', va='top',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        
        # Hide axis labels selectively (only show left column and bottom row)
        is_bottom_row = (row == n_rows - 1) or (i >= n_planets - n_rows)
        is_left_column = (col == 0)
        
        if not is_bottom_row:
            ax.set_xlabel('')
            ax.tick_params(axis='x', labelbottom=False)
        
        if not is_left_column:
            ax.set_ylabel('')
            ax.tick_params(axis='y', labelleft=False)
    # Hide unused subplots
    total_subplots = n_rows * n_cols
    for i in range(n_planets, total_subplots):
        # Calculate subplot position for unused plots
        col = i // n_rows
        row = i % n_rows
        subplot_index = row * n_cols + col + 1
        ax = fig.add_subplot(n_rows, n_cols, subplot_index)
        ax.set_visible(False)
    # Set overall title
    if hasattr(Params, 'TITLES') and Params.TITLES:
        fig.suptitle('Uncertainty Analysis - Multiple Best-Fit Planets', fontsize=15)
    
    # Save the combined figure
    plt.tight_layout()
    
    # Save to uncertainty multiplot file
    if hasattr(Params, 'FigureFiles') and hasattr(Params.FigureFiles, 'path'):
        fig_path = f"{Params.FigureFiles.path}/uncertainty_multiplot.pdf"
        fig.savefig(fig_path, format='pdf', dpi=300, bbox_inches='tight')
        log.debug(f'Uncertainty multiplot saved to file: {fig_path}')
    
    plt.close()


def PlotUncertaintySubplot(Exploration, Params, ax=None, planet_name=None):
    """
    Subplot version of PlotUncertainty - simply calls PlotUncertainty with ax parameter.
    
    Args:
        Exploration: Exploration results object
        Params: Parameters object with inversion settings
        ax: Matplotlib axis to plot on (required for subplot)
        planet_name: Optional planet name for custom subplot title
    """
    return PlotUncertainty(Exploration, Params, ax=ax, planet_name=planet_name)
def FitWithinUncertainty(ExplorationList, Params):
    """
    Fit the best-fit model within the uncertainty of the best-fit model.
    """
    for Exploration in ExplorationList:
        xInductionResponseUncertainty = CalcGridWithinUncertainty(Params.Inversion.Bi1xyz_nT['x'], Exploration.induction.Bi1x_nT, Params.Inversion.InductionResponseUncertainty_nT)
        yInductionResponseUncertainty = CalcGridWithinUncertainty(Params.Inversion.Bi1xyz_nT['y'], Exploration.induction.Bi1y_nT, Params.Inversion.InductionResponseUncertainty_nT)
        zInductionResponseUncertainty = CalcGridWithinUncertainty(Params.Inversion.Bi1xyz_nT['z'], Exploration.induction.Bi1z_nT, Params.Inversion.InductionResponseUncertainty_nT)
        Exploration.inversion.gridWithinInductionResponseUncertainty = np.all([xInductionResponseUncertainty, yInductionResponseUncertainty, zInductionResponseUncertainty], axis=0)
        Exploration.inversion.gridWithinkLoveAmpUncertainty = CalcGridWithinUncertainty(Params.Inversion.kLoveAmp, Exploration.base.kLoveAmp, Params.Inversion.kLoveAmpUncertainity)
        Exploration.inversion.gridWithinhLoveAmpUncertainty = CalcGridWithinUncertainty(Params.Inversion.hLoveAmp, Exploration.base.hLoveAmp, Params.Inversion.hLoveAmpUncertainity)
        Exploration.inversion.gridWithinAllUncertainty = np.all([Exploration.inversion.gridWithinInductionResponseUncertainty, Exploration.inversion.gridWithinkLoveAmpUncertainty, Exploration.inversion.gridWithinhLoveAmpUncertainty], axis=0)
    return ExplorationList
    
def CalcGridWithinUncertainty(bestPlanetData, GridData, UncertaintyData):
    """
    Calculate the grid of models within the uncertainty of the best-fit model.
    For complex data, uses magnitude of difference (circular uncertainty region).
    For real data, uses direct comparison (linear uncertainty bounds).
    
    Args:
        bestPlanetData: The best-fit model data value (real or complex)
        GridData: 2D array of grid data to compare against (real or complex)
        UncertaintyData: Uncertainty range for the best-fit model
        
    Returns:
        Boolean 2D array indicating which grid points are within uncertainty
    """
    # Check if data is complex
    if np.iscomplexobj(bestPlanetData) or np.iscomplexobj(GridData):
        # For complex data, check if magnitude of difference is within uncertainty
        # This creates a circular uncertainty region in the complex plane
        # GridData shape: (nexc, planetGridWidth, planetGridHeight)
        # bestPlanetData shape: (nexc,)
        # Need to broadcast bestPlanetData to match GridData dimensions
        bestPlanetData_expanded = bestPlanetData[:, np.newaxis, np.newaxis]
        realGridData = np.real(GridData)
        imagGridData = np.imag(GridData)
        realBestPlanetData = np.real(bestPlanetData_expanded)
        imagBestPlanetData = np.imag(bestPlanetData_expanded)
        realDifference = np.abs(realGridData - realBestPlanetData)
        imagDifference = np.abs(imagGridData - imagBestPlanetData)
        realWithinUncertainty = realDifference <= UncertaintyData
        imagWithinUncertainty = imagDifference <= UncertaintyData
        within_uncertainty = np.logical_and(realWithinUncertainty, imagWithinUncertainty)
        
        # If within_uncertainty is 3D (nexc, planetGridWidth, planetGridHeight),
        # check if any of the nexc rows are within uncertainty
        if within_uncertainty.ndim == 3:
            if within_uncertainty.shape[0] == 0:
                within_uncertainty = np.zeros((within_uncertainty.shape[1], within_uncertainty.shape[2]), dtype=bool)
            else:
                within_uncertainty = np.all(within_uncertainty, axis=0)
    else:
        # For real data, use direct comparison (linear bounds)
        within_uncertainty = _check_within_bounds(bestPlanetData, GridData, UncertaintyData)
    return within_uncertainty

def _compute_interpolated_uncertainty_regions(Exploration, Params, x_data, y_data, x_interp, y_interp):
    """
    Compute interpolation-aware uncertainty regions for combined contours.
    
    Args:
        Exploration: Exploration results object
        Params: Parameters object with inversion settings
        x_data, y_data: Original coordinate grids
        x_interp, y_interp: Interpolated coordinate grids
        
    Returns:
        Dictionary with boolean arrays for different parameter group combinations
    """
    regions = {}
    
    # Compute individual parameter uncertainty regions on interpolated grid
    k_love_region = None
    h_love_region = None
    induction_regions = []
    
    # k Love number region
    k_love_interp = _interpolate_z_data(x_data, y_data, Exploration.base.kLoveAmp, x_interp, y_interp)
    k_love_region = _check_within_bounds(
        Params.Inversion.kLoveAmp, k_love_interp, Params.Inversion.kLoveAmpUncertainity
    )
    
    # h Love number region  

    h_love_interp = _interpolate_z_data(x_data, y_data, Exploration.base.hLoveAmp, x_interp, y_interp)
    h_love_region = _check_within_bounds(
        Params.Inversion.hLoveAmp, h_love_interp, Params.Inversion.hLoveAmpUncertainity
    )
    
    # Induction component regions
    nExc, nExcNames = count_plottable_excitations(Exploration.induction.calcedExc, Params.Induct)
    inductionData = np.zeros((nExc, x_data.shape[0], x_data.shape[1]), dtype=np.complex_)
    bestInductionData = np.zeros((nExc), dtype=np.complex_)
    for iExc, nExcName in enumerate(nExcNames):
        Exploration.excName = nExcName
        data = extract_magnetic_field_data(Exploration, 'Bi1x_nT')
        inductionData[iExc, :, :] = data
        bestInductionData[iExc] = Params.Inversion.Bi1xyz_nT['x'][iExc]
    data_interp = _interpolate_z_data(x_data, y_data, inductionData, x_interp, y_interp)
    regions['induction'] = CalcGridWithinUncertainty(
        bestInductionData, data_interp, Params.Inversion.InductionResponseUncertainty_nT
    )
    
    # Combine regions using intersection (logical AND)
    regions['gravity'] = k_love_region & h_love_region
    
    regions['all_data'] = regions['gravity'] & regions['induction']

    return regions


def _plot_boolean_contour(ax, x_data, y_data, boolean_region, color, alpha):
    """
    Plot a contour from a boolean region mask.
    
    Args:
        ax: Matplotlib axis
        x_data, y_data: Coordinate arrays
        boolean_region: Boolean array indicating region to contour
        color: Contour color
        alpha: Transparency
        
    Returns:
        Matplotlib contour object or None
    """
    if boolean_region is None:
        return None
    
    # Convert boolean mask to numeric for contouring
    z_mask = np.where(boolean_region, 1.0, 0.0)
    
    # Plot filled contour
    contour = ax.contourf(x_data, y_data, z_mask,
                         levels=[0.5, 1.5],  # Contour at value 1.0
                         colors=[color], alpha=alpha, zorder=2)
    
    # Plot contour outline
    contour_outline = ax.contour(x_data, y_data, z_mask,
                               levels=[0.5],  # Contour at value 0.5 (edge of region)
                               colors=[color], alpha=0.8, linewidths=2)
    
    return contour


def _check_within_bounds(best_value, grid_values, uncertainty):
    """
    Helper function to check if grid values are within uncertainty bounds.
    
    Args:
        best_value: Reference value (scalar)
        grid_values: Array of values to check
        uncertainty: Uncertainty range
        
    Returns:
        Boolean array indicating which values are within bounds
    """
    upper_bound = best_value + uncertainty
    lower_bound = best_value - uncertainty
    return (grid_values >= lower_bound) & (grid_values <= upper_bound)
    
    
def PlotUncertainty(ExplorationList, Params, ax=None, planet_name=None):
    """
    Plot the uncertainty of the best-fit model.
    
    Creates a plot showing uncertainty regions using overlapping filled contours
    of physical quantities at best-fit ± uncertainty bounds. Shows separate
    contours for x, y, z induction components plus Love number uncertainties.
    Supports combined contours for parameter groups (gravity, induction, all data).
    
    Args:
        Exploration: Exploration results object
        Params: Parameters object with inversion settings
        ax: Optional matplotlib axis for subplot usage
        planet_name: Optional planet name for custom subplot title
    """
    from scipy.interpolate import griddata
    
    # Determine if we're creating a new figure or using existing axis
    create_new_figure = (ax is None)
    multiExploration = len(ExplorationList) > 1
    # Plot Configuration Dictionary
    # Controls which contours are displayed and their visual properties
    plot_config = {
        # Combined parameter group contours (intersection of individual uncertainties)
        'combined': {
            'gravity': False,         # Combined k + h Love number constraints
            'induction': False,       # Combined induction response constraints
            'all_data': True,        # All constraints combined (gravity + induction)
        },
        # Visual properties
        'alpha': {
            'individual': 0.2,       # Transparency for individual contours
            'combined': 0.4,         # Higher transparency for combined contours (more prominent)
        }
    }
    # Configurable variables (will make these parameters later)
    induction_x_color = 'orange'
    induction_y_color = 'orange'
    induction_z_color = 'darkred'
    klove_color = 'blue'  
    hlove_color = 'green'
    gravity_color = 'purple'         # Blended color for combined gravity constraints
    induction_color = 'darkorange'   # Blended color for combined induction constraints
    all_data_color = 'darkviolet'    # Color for all data combined
    true_model_color = 'black'
    true_model_marker = '*'
    true_model_size = 200
    contour_alpha = plot_config['alpha']['individual']
    combined_alpha = plot_config['alpha']['combined']
    figure_size = (3, 4)
    legend_font_size = 10
    title_font_size = 14
    interpolation_factor = 10 # Factor to increase grid resolution for smoother contours
    
    # Set up exploration data
    Exploration = ExplorationList[0]
    Exploration.zName = Exploration.xName
    FigLbl.SetExploration(Exploration.base.bodyname, Exploration.xName,
                          Exploration.yName, Exploration.zName)
    
    # Create figure and axis (only if not provided)
    if create_new_figure:
        fig = plt.figure(figsize=figure_size)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)
    else:
        # Using provided axis
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)
    
    # Set up basic plot title
    if create_new_figure:
        # For standalone plots, use full title with "Uncertainty Analysis"
        if hasattr(Params, 'TITLES') and Params.TITLES:
            ax.set_title(f'{Exploration.base.bodyname} Uncertainty Analysis', fontsize=title_font_size)
    else:
        # For subplots, use custom planet name or body name with x/y variables
        if planet_name:
            title = f'{planet_name}'
        else:
            title = f'{Exploration.base.bodyname}'
        
        # Add x and y variable info to subtitle
        if hasattr(FigLbl, 'xLabelExplore') and hasattr(FigLbl, 'yLabelExplore'):
            title += f'\n{FigLbl.xLabelExplore} vs {FigLbl.yLabelExplore}'
        
        ax.set_title(title, fontsize=title_font_size-2)  # Slightly smaller for subplots
    legend_elements = []
    for i, Exploration in enumerate(ExplorationList):
        # Extract x and y data from exploration
        data = extract_and_validate_plot_data(result_obj=Exploration, x_field=Exploration.xName, y_field=Exploration.yName,
                                            x_multiplier=FigLbl.xMultExplore, y_multiplier=FigLbl.yMultExplore,
                                            custom_x_axis=FigLbl.xCustomAxis, custom_y_axis=FigLbl.yCustomAxis)
        x_data = data['x'].reshape(data['original_shape'])
        y_data = data['y'].reshape(data['original_shape'])
        
        # Create higher resolution grid for smoother contours
        x_interp, y_interp = _create_interpolated_grid(x_data, y_data, interpolation_factor)
        
        ax.set_xlabel(FigLbl.xLabelExplore)
        ax.set_ylabel(FigLbl.yLabelExplore)
        ax.set_xscale(FigLbl.xScaleExplore)
        ax.set_yscale(FigLbl.yScaleExplore)
    
        
        # Compute interpolation-aware uncertainty regions for combined contours
        uncertainty_regions = _compute_interpolated_uncertainty_regions(
            Exploration, Params, x_data, y_data, x_interp, y_interp
        )
        
        # Define individual contour specifications
        individual_specs = [
            {
                'data': Exploration.base.kLoveAmp,
                'color': klove_color,
                'label': 'k Love Number',
                'best_value': Params.Inversion.kLoveAmp,
                'uncertainty': Params.Inversion.kLoveAmpUncertainity,
                'config_key': 'k_love',
                'group': 'gravity'
            },
            {
                'data': Exploration.base.hLoveAmp,
                'color': hlove_color,
                'label': 'h Love Number', 
                'best_value': Params.Inversion.hLoveAmp,
                'uncertainty': Params.Inversion.hLoveAmpUncertainity,
                'config_key': 'h_love',
                'group': 'gravity'
            },
            {
                'data': Exploration.induction.Bi1x_nT[0], 
                'color': induction_x_color,
                'label': 'Bi1x (Orbital)',
                'best_value': Params.Inversion.Bi1xyz_nT['x'][0],
                'uncertainty': Params.Inversion.InductionResponseUncertainty_nT,
                'config_key': 'induction_orbital',
                'group': 'induction'
            },
            {
                'data': Exploration.induction.Bi1x_nT[1], 
                'color': 'red',
                'label': 'Bi1x (Synodic)',
                'best_value': Params.Inversion.Bi1xyz_nT['x'][1],
                'uncertainty': Params.Inversion.InductionResponseUncertainty_nT,
                'config_key': 'induction_synodic',
                'group': 'induction'
            }]
            
            # Note: Additional induction components can be enabled by setting their config flags to True


        # Plot individual uncertainty contours (if enabled)
        if not np.any([plot_config['combined'][key] for key in plot_config['combined']]):
            for spec in individual_specs:
                if (spec['data'] is not None):
                
                    contour_fill = _plot_uncertainty_contour(
                        ax, x_data, y_data, spec['data'], 
                        spec['best_value'], spec['uncertainty'],
                        spec['color'], contour_alpha, x_interp, y_interp
                    )
                    if contour_fill is not None:
                        legend_elements.append(plt.Rectangle((0,0),1,1, fc=spec['color'], 
                                                        alpha=contour_alpha, label=spec['label']))
            
        # Plot combined uncertainty contours (if enabled)
        combined_specs = [
            {
                'region': uncertainty_regions['gravity'],
                'color': gravity_color,
                'label': 'Gravity Constraints (k+h Love)',
                'config_key': 'gravity'
            },
            {
                'region': uncertainty_regions['induction'],
                'color': induction_color,
                'label': 'Induction Constraints',
                'config_key': 'induction'
            },
            {
                'region': uncertainty_regions['all_data'],
                'color': all_data_color,
                'label': 'All Constraints Combined',
                'config_key': 'all_data'
            }
        ]
        
        for spec in combined_specs:
            if (spec['region'] is not None and 
                plot_config['combined'].get(spec['config_key'], False)):
                
                    if multiExploration:
                        label = Params.Explore.titleName[i]
                        color = Params.Explore.color[i]
                    else:
                        label = spec['label']
                        color = spec['color']
                    contour_fill = _plot_boolean_contour(
                        ax, x_interp, y_interp, spec['region'],
                        color, combined_alpha
                    )
                    legend_elements.append(plt.Rectangle((0,0),1,1, fc=color, ec=color,
                                                    alpha=combined_alpha, label=label))
    
        
        ax.set_xlim(np.min(x_data), np.max(x_data))
        ax.set_ylim(np.min(y_data), np.max(y_data))
        # Set minor grid lines with low alpha
        ax.grid(True, which='minor', alpha=0.1, linestyle='-', linewidth=0.5)
        ax.grid(True, which='major', alpha=0.5, linestyle='-', linewidth=0.5)
    # Add marker for true model location if available
    true_marker = ax.scatter(Params.Inversion.xBestFit, Params.Inversion.yBestFit, c=true_model_color, 
                            marker=true_model_marker, s=true_model_size, 
                            zorder=5, label='True Europa', 
                            edgecolors='white', linewidth=2)
    legend_elements.append(true_marker)
    
    # Add legend (only for standalone plots or first subplot)
    if legend_elements:
        ax.legend(handles=legend_elements, loc='best', fontsize=legend_font_size,
                framealpha=0.9)
    # Save and show only if we created a new figure
    if create_new_figure:
        # Tight layout and save (if file path is available)
        plt.tight_layout()
        
        # Save figure if Params has figure files configured
        if hasattr(Params, 'FigureFiles') and hasattr(Params.FigureFiles, 'path'):
            fig_path = f"{Params.FigureFiles.path}/uncertainty_plot.pdf"
            fig.savefig(fig_path, format='pdf', dpi=600)
            print(f"Uncertainty plot saved to: {fig_path}")
        
        plt.show()
        return fig, ax
    else:
        # Return the axis for multi-subplot usage
        return ax


def _create_interpolated_grid(x_data, y_data, interpolation_factor):
    """
    Create a higher resolution interpolated grid for smoother contours.
    
    Args:
        x_data: Original x coordinate array (2D)
        y_data: Original y coordinate array (2D)
        interpolation_factor: Factor to increase grid resolution
        
    Returns:
        x_interp, y_interp: Higher resolution coordinate arrays
    """
    # Get original grid dimensions
    nx, ny = x_data.shape
    
    # Create higher resolution grid
    new_nx = nx * interpolation_factor
    new_ny = ny * interpolation_factor
    
    # Create new coordinate arrays
    x_min, x_max = np.min(x_data), np.max(x_data)
    y_min, y_max = np.min(y_data), np.max(y_data)
    
    x_new = np.linspace(x_min, x_max, new_nx)
    y_new = np.linspace(y_min, y_max, new_ny)
    
    x_interp, y_interp = np.meshgrid(x_new, y_new, indexing='ij')
    
    return x_interp, y_interp


def _interpolate_z_data(x_data, y_data, z_data, x_interp, y_interp):
    """
    Interpolate z data onto higher resolution grid using linear interpolation.
    
    Args:
        x_data, y_data: Original coordinate arrays (2D)
        z_data: Original z data to interpolate
        x_interp, y_interp: Target interpolation grid
        
    Returns:
        z_interp: Interpolated z data on higher resolution grid
    """
    from scipy.interpolate import RectBivariateSpline
    xRow = x_data[:, 0]
    yRow = y_data[0, :]
    xInterpRow = x_interp[:, 0]
    yInterpRow = y_interp[0, :]
    
    # Handle different dimensionalities of z_data
    if z_data.ndim == 3:  # (nexc, ny, nx)
        nexc, ny, nx = z_data.shape
        # Initialize interpolated array with correct dtype
        if np.iscomplexobj(z_data):
            z_interp = np.zeros((nexc, x_interp.shape[0], x_interp.shape[1]), dtype=complex)
        else:
            z_interp = np.zeros((nexc, x_interp.shape[0], x_interp.shape[1]), dtype=float)

        for i in range(nexc):
            z_slice = z_data[i, :, :]
            
            # Handle complex data by interpolating real and imaginary parts separately
            if np.iscomplexobj(z_slice):
                splineReal = RectBivariateSpline(xRow, yRow, np.real(z_slice), kx=3, ky=3)
                splineImag = RectBivariateSpline(xRow, yRow, np.imag(z_slice), kx=3, ky=3)
                real_interp = splineReal(xInterpRow, yInterpRow)
                imag_interp = splineImag(xInterpRow, yInterpRow)
                
                z_interp_flat = real_interp + 1j * imag_interp
            else:
                spline = RectBivariateSpline(xRow, yRow, z_slice, kx=3, ky=3)
                z_interp_flat = spline(xInterpRow, yInterpRow)
            
            z_interp[i, :, :] = z_interp_flat
            
    else:  # 2D data (ny, nx)
        # Handle complex data by interpolating real and imaginary parts separately
        if np.iscomplexobj(z_data):
            zReal = np.real(z_data)
            zImag = np.imag(z_data)
            splineReal = RectBivariateSpline(xRow, yRow, zReal, kx=3, ky=3)
            splineImag = RectBivariateSpline(xRow, yRow, zImag, kx=3, ky=3)
            real_interp = splineReal(xInterpRow, yInterpRow)
            imag_interp = splineImag(xInterpRow, yInterpRow)
            z_interp = real_interp + 1j * imag_interp
        else:
            spline = RectBivariateSpline(xRow, yRow, z_data, kx=3, ky=3)
            z_interp = spline(xInterpRow, yInterpRow)
    
    return z_interp


def _extract_contour_data(grid_data, best_value, uncertainty):
    """
    Extract and process data for contour plotting, handling complex and 3D data.
    
    Args:
        grid_data: The exploration grid data 
        best_value: Best-fit value(s)
        uncertainty: Uncertainty range
        
    Returns:
        Processed 2D array suitable for contouring
    """
    if grid_data is None:
        return None
        
    # Handle complex data by taking magnitude
    if np.iscomplexobj(grid_data) or np.iscomplexobj(best_value):
        if grid_data.ndim == 3:  # (nexc, nx, ny)
            # For 3D complex data, compute RMS magnitude across excitations
            grid_magnitude = np.sqrt(np.mean(np.abs(grid_data)**2, axis=0))
        else:
            grid_magnitude = np.abs(grid_data)
        return grid_magnitude
    else:
        # For real data, return as-is
        if grid_data.ndim == 3:  # (nexc, nx, ny)
            # For 3D real data, compute RMS across excitations  
            grid_rms = np.sqrt(np.mean(grid_data**2, axis=0))
            return grid_rms
        else:
            return grid_data


def _plot_uncertainty_contour(ax, x_data, y_data, z_data, best_value, uncertainty, color, alpha, x_interp=None, y_interp=None):
    """
    Plot a single uncertainty contour with proper handling of complex/3D data.
    Uses interpolation for higher fidelity contours if interpolated grids are provided.
    
    Args:
        ax: Matplotlib axis
        x_data, y_data: Original coordinate arrays
        z_data: Data to contour
        best_value: Best-fit value (scalar or array)
        uncertainty: Uncertainty range
        color: Contour color
        alpha: Transparency
        x_interp, y_interp: Optional higher resolution interpolated grids
        
    Returns:
        Matplotlib contour object or None
    """
    # Use interpolated grids if provided, otherwise use original
    if x_interp is not None and y_interp is not None:
        x_plot = x_interp
        y_plot = y_interp
        z_plot = _interpolate_z_data(x_data, y_data, z_data, x_interp, y_interp)
    else:
        x_plot = x_data
        y_plot = y_data
        z_plot = z_data
    
    def _plot_single_contour(data_slice, best_val, alpha_val):
        """Helper function to plot a single contour for real data."""
        lower_level = best_val - uncertainty
        upper_level = best_val + uncertainty
        contour = ax.contourf(x_plot, y_plot, data_slice,
                          levels=[lower_level, upper_level],
                          colors=[color], alpha=alpha_val, zorder=2)
        contourOutline = ax.contour(x_plot, y_plot, data_slice,
                             levels=contour.levels,  # Contour at value 1.0
                             colors=[color], alpha=0.8, linewidths=2)
        return contour
    
    def _plot_complex_contours(data_slice, best_val_slice, alpha_val):
        """Helper function to plot contours for complex data (overlap of real and imaginary parts)."""
        # Extract real and imaginary parts
        real_data = np.real(data_slice)
        imag_data = np.imag(data_slice)
        best_real = np.real(best_val_slice)
        best_imag = np.imag(best_val_slice)
        
        # Create masks for uncertainty regions
        real_mask = np.abs(real_data - best_real) <= uncertainty
        imag_mask = np.abs(imag_data - best_imag) <= uncertainty
        
        # Find overlap region where both real and imaginary parts are within uncertainty
        overlap_mask = real_mask & imag_mask
        
        # Create a new z variable for the overlap region
        z_overlap = np.where(overlap_mask, 1.0, 0.0)
        
        # Plot contour of the overlap region
        contour = ax.contourf(x_plot, y_plot, z_overlap,
                             levels=[0.5, 1.5],  # Contour at value 1.0
                             colors=[color], alpha=alpha_val, zorder=2)
        contourOutline = ax.contour(x_plot, y_plot, z_overlap,
                             levels=contour.levels,  # Contour at value 1.0
                             colors=[color], alpha=0.8, linewidths=2)
        return [contour]
    
    contours = []
    
    # Check if data is 3D
    if z_plot.ndim == 3:  # (nexc, nx, ny)
        # Iterate over each row of 3D data
        for i in range(z_plot.shape[0]):
            z_slice = z_plot[i, :, :]
            best_val_slice = best_value[i] if hasattr(best_value, '__len__') else best_value
            
            if np.iscomplexobj(z_slice):
                contours.extend(_plot_complex_contours(z_slice, best_val_slice, alpha/2))
            else:
                contours.append(_plot_single_contour(z_slice, best_val_slice, alpha))
    else:
        # 2D data
        if np.iscomplexobj(z_plot):
            contours.extend(_plot_complex_contours(z_plot, best_value, alpha/2))
        else:
            best_mag = np.abs(best_value)
            contours.append(_plot_single_contour(z_plot, best_mag, alpha))
    return contours[0] if len(contours) == 1 else contours
