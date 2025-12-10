"""
Essential Plot Helpers - Minimal set for duplicated logic only

This contains only the essential helper functions for code that is clearly
duplicated across multiple plotting functions. Emphasis on readability
and simplicity over extensibility.
"""

import numpy as np
import matplotlib.pyplot as plt
import logging
from matplotlib.gridspec import GridSpec
from PlanetProfile.MagneticInduction.Moments import Excitations
from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc
from PlanetProfile.Utilities.defineStructs import xyzComps, vecComps
from matplotlib.patches import Rectangle, ConnectionPatch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
log = logging.getLogger('PlanetProfile')


def get_excitation_indices_and_names(induction_obj, params_induct):
    """
    Extract excitation selection logic from PlotInductOgram for reuse.
    
    This consolidates the complex excitation filtering logic that appears in:
    - PlotInductOgram (MagPlots.py lines ~340-348)
    - extract_monte_carlo_parameter_values (current function)  
    - Other magnetic plotting functions
    
    Args:
        induction_obj: Object with Texc_hr attribute (dict of excitation periods)
        params_induct: Parameters object with excSelectionPlot attribute
        
    Returns:
        tuple: (iTexc, iTexcAvail, TexcPlotNames, Texc_hr, iSort)
            - iTexc: Indices in full Texc_hr array for selected excitations
            - iTexcAvail: Indices in available (finite) Texc_hr array for selected excitations  
            - TexcPlotNames: Names of excitations to plot
            - Texc_hr: Periods of excitations to plot
            - iSort: Sort indices for Texc_hr (ascending order)
    """
    # Get indices for the oscillations that we can and want to plot
    excSelectionCalc = {key: Texc for key, Texc in zip(induction_obj.calcedExc, induction_obj.Texc_hr)}
    whichTexc = excSelectionCalc and params_induct.excSelectionPlot
    allTexc_hr = np.fromiter(excSelectionCalc.values(), dtype=np.float_)
    allAvailableTexc_hr = allTexc_hr[np.isfinite(allTexc_hr)]
    iTexc = [np.where(allTexc_hr == Texc)[0][0] for key, Texc in excSelectionCalc.items() 
             if whichTexc[key] and np.size(np.where(allTexc_hr == Texc)[0]) > 0]
    iTexcAvail = [np.where(allAvailableTexc_hr == Texc)[0][0] for key, Texc in excSelectionCalc.items() 
                  if whichTexc[key] and np.size(np.where(allAvailableTexc_hr == Texc)[0]) > 0]
    TexcPlotNames = np.fromiter(excSelectionCalc.keys(), dtype='<U20')[iTexc]
    Texc_hr = allTexc_hr[iTexc]
    iSort = np.argsort(Texc_hr)
    
    return iTexc, iTexcAvail, TexcPlotNames, Texc_hr, iSort


def get_available_excitations(result_obj, exc_selection_calc=None):
    """
    Get list of available excitation names for a result object.
    
    This assumes the result object has the hierarchical structure with induction substruct
    containing excSelectionCalc that matches the format:
    excSelectionCalc = {
        'orbital': True, 'synodic': True, 'orbital 2nd': False, 
        'synodic 2nd': False, 'true anomaly': False, ...
    }
    
    Args:
        result_obj: Result object with hierarchical structure (result.base.*, result.induction.*)
        exc_selection_calc: Optional excitation selection dict, uses result_obj.induction if None
        
    Returns:
        list: Available excitation names that are calculated and set to True
    """
    # Get bodyname from result object
    if hasattr(result_obj, 'base') and hasattr(result_obj.base, 'bodyname'):
        bodyname = result_obj.base.bodyname
    else:
        # Fallback for older structures
        bodyname = getattr(result_obj, 'bodyname', 'Unknown')
    
    # Get excSelectionCalc from induction substruct if not provided
    if exc_selection_calc is None:
        if hasattr(result_obj, 'induction') and hasattr(result_obj.induction, 'excSelectionCalc'):
            exc_selection_calc = result_obj.induction.excSelectionCalc
        else:
            # Fallback for older structures
            exc_selection_calc = getattr(result_obj, 'excSelectionCalc', None)
    
    # Return excitations that are set to True
    if exc_selection_calc is not None:
        return [key for key, calc in exc_selection_calc.items() if calc]
    
    log.warning(f'No excSelectionCalc found for {bodyname}')
    return []


def count_plottable_excitations(calcedExc, params_induct):
    """
    Count number of excitations that can be plotted for PlotComplexBdip logic.
    
    This consolidates the complex counting logic from PlotComplexBdip that
    determines if we have enough excitations to plot.
    
    Args:
        bodyname: Name of celestial body
        params_induct: Induction parameters object
        
    Returns:
        int: Number of excitations that can be plotted
    """
    plottable_keys = []
    for key in calcedExc:
        if key in params_induct.excSelectionPlot.keys() and params_induct.excSelectionPlot[key]:
            plottable_keys.append(key)
    
    return len(plottable_keys), plottable_keys
        

def draw_ocean_composition_lines(ax, x, y, c_values, ocean_comp, 
                                use_manual_colors=False, 
                                line_width=1.0, line_alpha=0.7):
    """
    Draw connecting lines between points of the same ocean composition.
    
    This logic block appears identically in multiple plotting functions:
    - PlotExploreOgramDsigma (appears 2x: individual + combined plots)  
    - PlotExploreOgramZbD (appears 2x: individual + combined plots)
    - Other exploration plotting functions
    
    Args:
        ax: Matplotlib axis to plot on
        x, y: Data arrays (already flattened and valid)
        c_values: Color values (ice thickness, etc.) for manual coloring
        ocean_comp: Ocean composition labels for each point
        use_manual_colors: Whether to use manual hydro colors
        Color: Color configuration object (from PlanetProfile.Utilities)
        line_width: Line width for composition lines
        line_alpha: Line transparency
        
    Returns:
        plotted_labels: Set of composition labels that were plotted (for legend)
    """
    
    if len(x) == 0 or len(ocean_comp) == 0:
        return set()
    use_manual_colors = False
    # Group by unique ocean compositions (maintaining order) and plot connecting lines
    unique_comps, indices = np.unique(ocean_comp, return_index=True)
    unique_comps = unique_comps[np.argsort(indices)]
    plotted_labels = set()
    
    for comp in unique_comps:
        comp_indices = np.where(ocean_comp == comp)[0]
        x_line, y_line = x[comp_indices], y[comp_indices]
        
        # Sort points by x for connected lines
        sorted_idx = np.argsort(x_line)
        x_line = x_line[sorted_idx]
        y_line = y_line[sorted_idx]
        
        
        thisColor = Color.cmap[comp](0.5)
        
        # Clean composition label
        if 'CustomSolution' in comp:
            comp_label = comp.split('=')[0].replace('CustomSolution', '')
        else:
            comp_label = comp
        
        # Plot line for this composition
        ax.plot(x_line, y_line, color=thisColor, linewidth=line_width, 
               alpha=line_alpha, label=comp_label, zorder=2)
        
        plotted_labels.add(comp_label)
    
    return plotted_labels


def create_ice_thickness_colorbar(ax, ice_thickness_values, cmap_name='Greys',
                                 label='Ice Shell Thickness (km)', 
                                 orientation='vertical', shrink=0.8, pad=0.05):
    """
    Create colorbar showing ice thickness values.
    
    Args:
        ax: Matplotlib axis to add colorbar to
        ice_thickness_values: Array of ice thickness values
        cmap_name: Colormap name for ice thickness coloring
        label: Label for the colorbar
        orientation: Colorbar orientation ('vertical' or 'horizontal')
        shrink: Fraction by which to multiply the size of the colorbar
        pad: Fraction of original axes between colorbar and new image axes
        
    Returns:
        colorbar: The created colorbar object
    """
    
    if len(ice_thickness_values) == 0:
        return None
    
    # Get thickness bounds for normalization
    t_min, t_max = np.min(ice_thickness_values), np.max(ice_thickness_values)
    
    # Create colormap and normalization
    cmap = getattr(plt.cm, cmap_name)
    norm = plt.Normalize(vmin=t_min, vmax=t_max)
    
    # Create a ScalarMappable for the colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Required for ScalarMappable
    
    # Create and return the colorbar
    colorbar = plt.colorbar(sm, ax=ax, label=label, orientation=orientation,
                           shrink=shrink, pad=pad)
    
    return colorbar

def getIceShellThickness(result_obj):
    """
    Get ice shell thickness from result object, depending on whether we have a zb_approximate_km field (and thus used Planet.Do.zb_approximate_km) or not.
    """
    if np.any(np.isnan(result_obj.base.zb_approximate_km)):
        return result_obj.base.zb_km, 'zb_km'
    else:
        return result_obj.base.zb_approximate_km, 'zb_approximate_km'

def add_composition_legend(ax, ocean_comps,
                          show_ice_thickness_legend=False,
                          font_size=8, title_font_size=10, overrideLegend=None):
    """
    Add ocean composition legend to plot.
    
    This logic appears in multiple plotting functions with slight variations.
    
    Args:
        ax: Matplotlib axis to add legend to
        ocean_comp: Ocean composition array (for determining if legend needed)
        show_ice_thickness_legend: Whether ice thickness legend is also shown
        font_size: Font size for legend text
        title_font_size: Font size for legend title
    """   
    if len(ocean_comps) == 0:
        return

    if overrideLegend is not None and len(overrideLegend) > 0:
        
        legend = ax.legend(handles=overrideLegend, fontsize=font_size, 
                     title_fontsize=title_font_size)
        
    else:
        # Position legend to avoid overlap with ice thickness legend
        if show_ice_thickness_legend:
            legend =ax.legend(title="Ocean Composition", fontsize=font_size, 
                    title_fontsize=title_font_size,
                    loc='upper left', bbox_to_anchor=(0.0, 1.0))
        else:
            legend = ax.legend(title="Ocean Composition", fontsize=font_size, 
                        title_fontsize=title_font_size, loc = 'upper left')
    ax.add_artist(legend)


def extract_and_validate_plot_data(result_obj, x_field, y_field, c_field=None, contour_field=None,
                                  x_multiplier=1.0, y_multiplier=1.0, c_multiplier=1.0, contour_multiplier=1.0,
                                  custom_x_axis = None, custom_y_axis = None):
    """
    Extract and validate basic plot data from a result object.
    
    This follows the original pattern:
    x = np.reshape(result.__getattribute__(x_field) * multiplier, -1)
    VALID = np.logical_not(np.logical_or(np.isnan(x), np.isnan(y)))
    
    Args:
        result_obj: Result object with hierarchical structure (result.base.field_name)
        x_field, y_field: Field names to extract from result.base
        c_field: Optional color field name (for colormap)
        contour_field: Optional contour field name (for contour lines)
        x_multiplier, y_multiplier, c_multiplier, contour_multiplier: Scaling factors
        custom_x_axis, custom_y_axis: Custom x and y axes (for example, if we want to plot a custom x axis)
        
    Returns:
        dict: {'x': x_valid, 'y': y_valid, 'c': c_valid, 'contour': contour_valid, 'original_shape': shape}
    """
    
    # Extract data using original pattern - check base structure
    if hasattr(result_obj.base, x_field):
        x_raw = result_obj.base.__dict__[x_field]
        # Handle string vs numeric data (like in PlotExploreOgramModern)
        if np.issubdtype(x_raw.dtype, np.number):
            x = np.reshape(x_raw * x_multiplier, -1)
        else:
            x = np.array([])
            for i in range(x_raw.shape[0]):
                row = np.repeat(int(i), x_raw.shape[1])
                x = np.concatenate((x, row))
            x = x.reshape(-1)
        if custom_x_axis is not None:
            if len(custom_x_axis) == len(x):
                x = custom_x_axis
            elif len(custom_x_axis) == x_raw.shape[0]:
                x = np.repeat(custom_x_axis, x_raw.shape[1])
            else:
                log.warning(f"Custom x axis length {len(custom_x_axis)} does not match x length {len(x)}")
        else:
            x = x
    else:
        raise ValueError(f"Field {x_field} not found in result.base")
        
    if hasattr(result_obj.base, y_field):
        y_raw = result_obj.base.__dict__[y_field]
        # Handle string vs numeric data (like in PlotExploreOgramModern)
        if np.issubdtype(y_raw.dtype, np.number):
            y = np.reshape(y_raw * y_multiplier, -1)
        else:
            y = np.array([])
            for i in range(y_raw.shape[0]):
                y = np.concatenate((y, np.arange(0, y_raw.shape[1], dtype=np.float_)))
            y = y.reshape(-1)
        if custom_y_axis is not None:
            if len(custom_y_axis) == len(y):
                y = custom_y_axis
            elif len(custom_y_axis) == y_raw.shape[1]:
                y = np.tile(custom_y_axis, y_raw.shape[0])
            else:
                log.warning(f"Custom y axis length {len(custom_y_axis)} does not match y length {len(y)}")
        else:
            y = y
    else:
        raise ValueError(f"Field {y_field} not found in result.base")

    
    x_valid = x
    y_valid  = y
    original_shape = x_raw.shape
    
    # Extract color field if provided
    c_valid = None
    if c_field is not None:
        if 'Induction' in c_field:
            c_raw = extract_magnetic_field_data(result_obj, c_field)
        elif hasattr(result_obj.base, c_field):
            c_raw = result_obj.base.__dict__[c_field]
        else:
            log.warning(f"Color field {c_field} not found in result.base or result.induction, using default coloring")
        # Handle string vs numeric data (like in PlotExploreOgramModern)
        if np.issubdtype(c_raw.dtype, np.number):
            c = np.reshape(c_raw * c_multiplier, -1)
        else:
            c = np.reshape(c_raw, -1)  # Don't multiply string data
        c_valid = c

    # Extract contour field if provided
    contour_valid = None
    if contour_field is not None:
        if hasattr(result_obj.base, contour_field):
            contour_raw = result_obj.base.__dict__[contour_field]
            # Handle string vs numeric data (like in PlotExploreOgramModern)
            if np.issubdtype(contour_raw.dtype, np.number):
                contour = np.reshape(contour_raw * contour_multiplier, -1)
            else:
                contour = np.reshape(contour_raw, -1)  # Don't multiply string data
            contour_valid = contour
        elif hasattr(result_obj, 'induction') and hasattr(result_obj.induction, contour_field):
            # Handle magnetic induction fields from 3D data
            contour_valid = extract_magnetic_field_data(result_obj, contour_field, contour_multiplier)
        else:
            log.warning(f"Contour field {contour_field} not found in result.base or result.induction, contours will use color field")
    elif contour_field is None and c_field is not None:
        contour_valid = c_valid
    return {
        'x': x_valid,
        'y': y_valid, 
        'c': c_valid,
        'contour': contour_valid,
        'original_shape': original_shape
    }


def extract_complex_plot_data(result_obj, data_type, Params):
    """
    Extract complex number data from result objects for real vs imaginary plotting.
    
    Handles both magnetic induction and tidal Love numbers automatically.
    
    Args:
        result_obj: Result object with hierarchical structure
        data_type: 'magnetic' for magnetic induction or 'love' for love numbers
        
    Returns:
        dict: {
            'data_type': 'magnetic' or 'love',
            'complex_data': {comp: complex_array for comp in components},
            'excitation_names': [list of excitation names] or ['l2'] for Love,
            'n_peaks': number of peaks/excitations
        }
    """
    # Auto-detect data type based on available data
    has_magnetic = (hasattr(result_obj, 'induction') and 
                    hasattr(result_obj.induction, 'Bi1x_nT') and
                    result_obj.induction.Bi1x_nT is not None)
    has_love = (hasattr(result_obj, 'base') and 
                hasattr(result_obj.base, 'kLoveComplex') and
                result_obj.base.kLoveComplex is not None)
    
    if data_type == 'magnetic':
        if not has_magnetic:
            raise ValueError("Magnetic induction data not found")
        
        complex_data = {}
        complex_data['x'] = np.array(result_obj.induction.Bi1x_nT)
        complex_data['y'] = np.array(result_obj.induction.Bi1y_nT)
        complex_data['z'] = np.array(result_obj.induction.Bi1z_nT)
        excitationNamesCalc = [key for key, calc in result_obj.induction.excSelectionCalc.items() 
                              if calc and key != 'none']
        excictationNamesToPlot = [key for key, calc in Params.Induct.excSelectionPlot.items() 
                                  if calc and key != 'none']
        excitation_names = [key for key in excitationNamesCalc if key in excictationNamesToPlot]
        
        
        n_peaks = len(excitation_names)
        
    elif data_type == 'love':
        if not has_love:
            raise ValueError("Love number data not found")
        
        complex_data = {}
        complex_data['k'] = np.array(result_obj.base.kLoveComplex)
        complex_data['h'] = np.array(result_obj.base.hLoveComplex)
        complex_data['l'] = np.array(result_obj.base.lLoveComplex)
        complex_data['delta'] = np.array(result_obj.base.deltaLoveComplex)
        
        # Love numbers are typically for l=2 harmonic only
        excitation_names = [2]
        n_peaks = 1
        
    else:
        raise ValueError(f"Unknown data_type: {data_type}")
    
    components = list(complex_data.keys())
    return {
        'data_type': data_type,
        'components': components,
        'complex_data': complex_data,
        'excitation_names': excitation_names,
        'n_peaks': n_peaks
    }


def normalizeDataForColor(data, dataName, cmap, highlight_mask=None):
    """
    Normalize data for color mapping.
    
    Args:
        data: data to normalize
        cmap: colormap to use
    """
    if dataName in FigLbl.axisLabelsExplore.keys():
        cBarTitle = FigLbl.axisLabelsExplore[dataName]
    else:
        cBarTitle = dataName
    if dataName in FigLbl.cbarfmtExplore.keys():
        cBarFmt = FigLbl.cbarfmtExplore[dataName]
    else:
        cBarFmt = None

    vmin, vmax = np.nanmin(data), np.nanmax(data)
    if vmax > vmin:
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        color = cmap(norm(data))
        color = color.reshape(-1, 4)
    else:
        color = cmap(0.5)
        color = np.full((len(data), 4), color)
    if highlight_mask is not None:
        color[~highlight_mask] = [0, 0, 0, 0]
    if dataName in FigLbl.cTicksSpacingsExplore.keys():
        cTicksSpacings = FigLbl.cTicksSpacingsExplore[dataName]
        cbar_ticks = np.arange(np.ceil(vmin / cTicksSpacings) * cTicksSpacings, np.floor(vmax / cTicksSpacings) * cTicksSpacings + 0.0001, cTicksSpacings)
    else:
        cbar_ticks = None
    return norm, color, cBarTitle, cBarFmt, cbar_ticks


def extract_monte_carlo_parameter_values(result_obj, param_name, excitation_name=None, valid_mask=None):
    """
    Extract parameter values from Monte Carlo results with excitation support.
    
    This consolidates the complex parameter extraction logic that appears
    in multiple Monte Carlo plotting functions. Now handles hierarchical structure
    with induction substruct for magnetic parameters.
    
    Args:
        result_obj: Monte Carlo results object with hierarchical structure (result.base.*, result.induction.*)
        param_name: Name of parameter to extract
        excitation_name: Optional excitation name for magnetic parameters (currently not used - field values are pre-calculated)
        valid_mask: Optional mask for valid runs
        
    Returns:
        np.array: Parameter values
    """
    if valid_mask is None:
        # For hierarchical structure, get nRuns from base structure
        n_runs = getattr(result_obj.base, 'nRuns', len(result_obj.base.VALID))
        valid_mask = np.ones(n_runs, dtype=bool)
    
    # Handle magnetic parameters - these are stored in the induction substruct with 3D structure
    if param_name in ['Amp', 'Bix_nT', 'Biy_nT', 'Biz_nT', 'Bi_nT', 'phase']:
        # Get the parameter values from induction substruct - expect 3D structure (nPeaks, rows, cols)
        values = getattr(result_obj.induction, param_name)
        values_array = np.array(values)
        
        if excitation_name is not None:
            # Get available excitations to find the correct index
            available_excitations = get_available_excitations(result_obj)
            exc_index = available_excitations.index(excitation_name)
            
            # Extract the 2D slice for this excitation: values[exc_index, :, :]
            excitation_slice = values_array[exc_index, :, :]
            return excitation_slice
        else:
            # No excitation specified - return first excitation as default
            return values_array[0, :, :]
    
    else:
        # Handle regular parameters - these are in the base substruct
        values = getattr(result_obj.base, param_name)
        return np.array(values)


def create_zoom_inset(ax, x_data, y_data, zoom_mask, padding_factor=0.2, 
                     add_visual_indicators=True):
    """
    Create an inset zoom plot showing highlighted data points with optional visual indicators.
    
    Args:
        ax: Main axis to add inset to
        x_data, y_data: Full data arrays
        zoom_mask: Boolean mask for points to zoom in on
        padding_factor: Fraction of data range to add as padding
        add_visual_indicators: Whether to add rectangle and connection lines
        
    Returns:
        tuple: (inset_ax, optimal_location) or (None, None) if no zoom data
    """
    
    if not np.any(zoom_mask):
        return None, None
    
    # Get zoom data
    x_zoom = x_data[zoom_mask]
    y_zoom = y_data[zoom_mask]
    
    # Calculate zoom limits with padding
    x_range = np.max(x_zoom) - np.min(x_zoom)
    y_range = np.max(y_zoom) - np.min(y_zoom)
    
    x_padding = x_range * padding_factor if x_range > 0 else 0.1
    y_padding = y_range * padding_factor if y_range > 0 else 0.1
    
    # Check if zoom region is already a large portion of the main axis
    # If so, don't create inset as it won't provide meaningful zoom
    ax_xlim = ax.get_xlim()
    ax_ylim = ax.get_ylim()
    
    ax_x_range = ax_xlim[1] - ax_xlim[0]
    ax_y_range = ax_ylim[1] - ax_ylim[0]
    
    # Calculate what fraction of the axis the zoom region covers
    zoom_x_fraction = (x_range + 2 * x_padding) / ax_x_range if ax_x_range > 0 else 0
    zoom_y_fraction = (y_range + 2 * y_padding) / ax_y_range if ax_y_range > 0 else 0
    
    # If zoom region covers more than 50% in both directions, don't create inset
    if zoom_x_fraction > 0.5 and zoom_y_fraction > 0.5:
        return None, None
    
    x_min = max(0, np.min(x_zoom) - x_padding)  # Ensure non-negative for complex plots
    x_max = np.max(x_zoom) + x_padding
    y_min = max(0, np.min(y_zoom) - y_padding)  # Ensure non-negative for complex plots
    y_max = np.max(y_zoom) + y_padding
    
    # Find optimal inset position based on point density
    optimal_location = find_optimal_inset_position(ax, x_data, y_data, inset_size_frac=0.4)
    
    # Create inset in optimal position
    inset_ax = inset_axes(ax, width="50%", height="50%", loc=optimal_location)
    inset_ax.set_xlim(x_min, x_max)
    inset_ax.set_ylim(y_min, y_max)
    
    # Dynamically set tick positions based on inset location
    # This prevents ticks from overlapping with the main plot data
    if 'lower' in optimal_location:
        # For lower positions, put x-ticks on top to avoid main plot data
        inset_ax.xaxis.tick_top()
        inset_ax.xaxis.set_label_position('top')
    else:
        # For upper positions, keep x-ticks on bottom (default)
        inset_ax.xaxis.tick_bottom()
        inset_ax.xaxis.set_label_position('bottom')
    
    if 'right' in optimal_location:
        # For right positions, put y-ticks on left to avoid main plot data
        inset_ax.yaxis.tick_left()
        inset_ax.yaxis.set_label_position('left')
    else:
        # For left positions, put y-ticks on right
        inset_ax.yaxis.tick_right()
        inset_ax.yaxis.set_label_position('right')
    
    # Add visual indicators if requested
    if add_visual_indicators:
        
        try:
            # Add rectangle around zoom region on main plot
            zoom_rect = Rectangle((x_min, y_min), 
                                x_max - x_min, 
                                y_max - y_min,
                                fill=False, edgecolor=Color.BdipInset, 
                                linewidth=Style.LW_BdipInset, 
                                linestyle=Style.LS_BdipInset, 
                                zorder=10)
            ax.add_patch(zoom_rect)
            
            # Add border around inset axis
            inset_ax.spines['top'].set_color(Color.BdipInset)
            inset_ax.spines['top'].set_linewidth(Style.LW_BdipInset * 1.5)
            inset_ax.spines['bottom'].set_color(Color.BdipInset)
            inset_ax.spines['bottom'].set_linewidth(Style.LW_BdipInset * 1.5)
            inset_ax.spines['left'].set_color(Color.BdipInset)
            inset_ax.spines['left'].set_linewidth(Style.LW_BdipInset * 1.5)
            inset_ax.spines['right'].set_color(Color.BdipInset)
            inset_ax.spines['right'].set_linewidth(Style.LW_BdipInset * 1.5)
            
            # Add smart connection lines based on inset location
            if 'upper' in optimal_location:
                # Upper insets: connect top of rectangle to bottom of inset
                con1 = ConnectionPatch((x_min, y_max), (0, 0), 
                                     "data", "axes fraction",
                                     axesA=ax, axesB=inset_ax,
                                     color=Color.BdipInset, 
                                     linewidth=Style.LW_BdipInset * 0.8,
                                     linestyle='--', alpha=0.7)
                ax.figure.add_artist(con1)
                
                con2 = ConnectionPatch((x_max, y_max), (1, 0), 
                                     "data", "axes fraction",
                                     axesA=ax, axesB=inset_ax,
                                     color=Color.BdipInset, 
                                     linewidth=Style.LW_BdipInset * 0.8,
                                     linestyle='--', alpha=0.7)
                ax.figure.add_artist(con2)
                
            elif 'lower right' in optimal_location:
                # Lower right inset: connect right of rectangle to left of inset
                con1 = ConnectionPatch((x_max, y_min), (0, 0), 
                                     "data", "axes fraction",
                                     axesA=ax, axesB=inset_ax,
                                     color=Color.BdipInset, 
                                     linewidth=Style.LW_BdipInset * 0.8,
                                     linestyle='--', alpha=0.7)
                ax.figure.add_artist(con1)
                
                con2 = ConnectionPatch((x_max, y_max), (0, 1), 
                                     "data", "axes fraction",
                                     axesA=ax, axesB=inset_ax,
                                     color=Color.BdipInset, 
                                     linewidth=Style.LW_BdipInset * 0.8,
                                     linestyle='--', alpha=0.7)
                ax.figure.add_artist(con2)
                
            elif 'lower left' in optimal_location:
                # Lower left inset: connect left of rectangle to right of inset
                con1 = ConnectionPatch((x_min, y_min), (1, 0), 
                                     "data", "axes fraction",
                                     axesA=ax, axesB=inset_ax,
                                     color=Color.BdipInset, 
                                     linewidth=Style.LW_BdipInset * 0.8,
                                     linestyle='--', alpha=0.7)
                ax.figure.add_artist(con1)
                
                con2 = ConnectionPatch((x_min, y_max), (1, 1), 
                                     "data", "axes fraction",
                                     axesA=ax, axesB=inset_ax,
                                     color=Color.BdipInset, 
                                     linewidth=Style.LW_BdipInset * 0.8,
                                     linestyle='--', alpha=0.7)
                ax.figure.add_artist(con2)
                
            else:
                # Fallback: connect based on position
                if 'right' in optimal_location:
                    # Right side: connect right edge to left of inset
                    con1 = ConnectionPatch((x_max, y_min), (0, 0), 
                                         "data", "axes fraction",
                                         axesA=ax, axesB=inset_ax,
                                         color=Color.BdipInset, 
                                         linewidth=Style.LW_BdipInset * 0.8,
                                         linestyle='--', alpha=0.7)
                    ax.figure.add_artist(con1)
                    
                    con2 = ConnectionPatch((x_max, y_max), (0, 1), 
                                         "data", "axes fraction",
                                         axesA=ax, axesB=inset_ax,
                                         color=Color.BdipInset, 
                                         linewidth=Style.LW_BdipInset * 0.8,
                                         linestyle='--', alpha=0.7)
                    ax.figure.add_artist(con2)
                else:
                    # Left side: connect left edge to right of inset
                    con1 = ConnectionPatch((x_min, y_min), (1, 0), 
                                         "data", "axes fraction",
                                         axesA=ax, axesB=inset_ax,
                                         color=Color.BdipInset, 
                                         linewidth=Style.LW_BdipInset * 0.8,
                                         linestyle='--', alpha=0.7)
                    ax.figure.add_artist(con1)
                    
                    con2 = ConnectionPatch((x_min, y_max), (1, 1), 
                                         "data", "axes fraction",
                                         axesA=ax, axesB=inset_ax,
                                         color=Color.BdipInset, 
                                         linewidth=Style.LW_BdipInset * 0.8,
                                         linestyle='--', alpha=0.7)
                    ax.figure.add_artist(con2)
                    
        except Exception as e:
            log.debug(f"Could not add visual indicators to zoom inset: {e}")
            # Continue without visual indicators if they fail
            pass
    
    return inset_ax, optimal_location


def find_optimal_inset_position(ax, x_data, y_data, inset_size_frac=0.4):
    """
    Find the optimal position for an inset plot based on point density.
    
    Args:
        ax: Main axis object
        x_data, y_data: Full data arrays
        inset_size_frac: Fraction of axis space the inset will occupy
        
    Returns:
        str: Best location string for inset_axes
    """
    # Get axis limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    # Define potential inset positions and their corresponding regions
    # Each position maps to a region where we'll count points
    positions = {
        'upper right': (xlim[1] - (xlim[1] - xlim[0]) * inset_size_frac, xlim[1],
                       ylim[1] - (ylim[1] - ylim[0]) * inset_size_frac, ylim[1]),
        'upper left': (xlim[0], xlim[0] + (xlim[1] - xlim[0]) * inset_size_frac,
                      ylim[1] - (ylim[1] - ylim[0]) * inset_size_frac, ylim[1]),
        'upper center': (xlim[0] + (xlim[1] - xlim[0]) * (0.5 - inset_size_frac/2), 
                        xlim[0] + (xlim[1] - xlim[0]) * (0.5 + inset_size_frac/2),
                        ylim[1] - (ylim[1] - ylim[0]) * inset_size_frac, ylim[1]),
        'lower right': (xlim[1] - (xlim[1] - xlim[0]) * inset_size_frac, xlim[1],
                       ylim[0], ylim[0] + (ylim[1] - ylim[0]) * inset_size_frac),
        'lower left': (xlim[0], xlim[0] + (xlim[1] - xlim[0]) * inset_size_frac,
                      ylim[0], ylim[0] + (ylim[1] - ylim[0]) * inset_size_frac),
        'lower center': (xlim[0] + (xlim[1] - xlim[0]) * (0.5 - inset_size_frac/2), 
                        xlim[0] + (xlim[1] - xlim[0]) * (0.5 + inset_size_frac/2),
                        ylim[0], ylim[0] + (ylim[1] - ylim[0]) * inset_size_frac),
        'center right': (xlim[1] - (xlim[1] - xlim[0]) * inset_size_frac, xlim[1],
                        ylim[0] + (ylim[1] - ylim[0]) * 0.3, ylim[0] + (ylim[1] - ylim[0]) * 0.7),
        'center left': (xlim[0], xlim[0] + (xlim[1] - xlim[0]) * inset_size_frac,
                       ylim[0] + (ylim[1] - ylim[0]) * 0.3, ylim[0] + (ylim[1] - ylim[0]) * 0.7),
    }
    
    # Count points in each potential inset region
    point_counts = {}
    for location, (x_min, x_max, y_min, y_max) in positions.items():
        # Count points that fall within this region
        in_region = ((x_data >= x_min) & (x_data <= x_max) & 
                    (y_data >= y_min) & (y_data <= y_max))
        point_counts[location] = np.sum(in_region)
    
    # Find location with minimum point density
    best_location = min(point_counts, key=point_counts.get)
    
    # If all positions have similar point counts, prefer corner positions for aesthetics
    min_count = min(point_counts.values())
    corner_positions = ['upper right', 'upper left', 'lower right', 'lower left']
    
    # If multiple positions have the same minimal count, prefer corners
    minimal_positions = [loc for loc, count in point_counts.items() if count == min_count]
    corner_minimal = [loc for loc in minimal_positions if loc in corner_positions]
    
    if corner_minimal:
        # Among corners with minimal points, prefer upper right as default
        if 'upper right' in corner_minimal:
            return 'upper right'
        else:
            return corner_minimal[0]
    else:
        return best_location


def formatOceanCompositionLabel(comp):
    """
    Format ocean composition labels for display.
    
    Args:
        comp: Composition string
        
    Returns:
        str: Formatted label
    """
    if 'CustomSolution' in comp:
        return comp.split('=')[0].replace('CustomSolution', '').strip()
    else:
        return comp


def extract_magnetic_field_data(result_obj, field_name):

    excitation_name = result_obj.excName
    inductionFieldName = field_name.replace('Induction', '')
    calcedExc = result_obj.induction.calcedExc
    iExcName = calcedExc.index(excitation_name)
    fullMagneticInductionData = getattr(result_obj.induction, inductionFieldName)
    excExcitationData = fullMagneticInductionData[iExcName, :, :]
    return excExcitationData

