"""
Enhanced Plotly plotting for exploreograms that matches matplotlib styling from PlanetProfile.
"""
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
import os

def get_matplotlib_colormap_colors(cmap_name='viridis', n_colors=256):
    """
    Get colors from matplotlib colormap for use in Plotly.
    Falls back to viridis if matplotlib not available.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        from matplotlib.colors import LinearSegmentedColormap

        cmap = cm.get_cmap(cmap_name)
        colors = [cmap(i / (n_colors - 1)) for i in range(n_colors)]
        # Convert RGBA to hex
        hex_colors = ['#{:02x}{:02x}{:02x}'.format(
            int(r*255), int(g*255), int(b*255)
        ) for r, g, b, a in colors]
        return hex_colors
    except:
        # Fallback to plotly viridis if matplotlib not available
        return None


def create_exploreogram_plotly(Exploration, Params, FigLbl=None, smoothing=False, smooth_factor=2):
    """
    Create Plotly version of exploreogram that matches matplotlib styling.

    Args:
        Exploration: ExplorationResultsStruct with results
        Params: ParamsStruct with configuration
        FigLbl: Figure labels struct (optional, for advanced styling)
        smoothing: Whether to apply smoothing
        smooth_factor: Smoothing factor if smoothing enabled

    Returns:
        Plotly figure object
    """

    # Extract data with safety checks
    if not hasattr(Exploration, 'xData') or Exploration.xData is None:
        raise ValueError("Exploration.xData is None or missing")
    if not hasattr(Exploration, 'yData') or Exploration.yData is None:
        raise ValueError("Exploration.yData is None or missing")
    if not hasattr(Exploration, 'zName') or Exploration.zName is None:
        raise ValueError("Exploration.zName is None or missing")

    xData = Exploration.xData
    yData = Exploration.yData
    zName = Exploration.zName

    # Get z-variable data
    if not hasattr(Exploration, 'base') or Exploration.base is None:
        raise ValueError("Exploration.base is None or missing")

    # Check if this is a magnetic induction variable
    induction_vars = ['Amp_nT', 'Bix_nT', 'Biy_nT', 'Biz_nT', 'phase_deg',
                      'Bi1x_nT', 'Bi1y_nT', 'Bi1z_nT', 'Bi1Tot_nT',
                      'rBi1x_nT', 'rBi1y_nT', 'rBi1z_nT', 'rBi1Tot_nT',
                      'iBi1x_nT', 'iBi1y_nT', 'iBi1z_nT', 'iBi1Tot_nT']

    if zName in induction_vars:
        # Extract from induction substructure
        if not hasattr(Exploration, 'induction') or Exploration.induction is None:
            raise ValueError(f"Magnetic induction data not found for variable '{zName}'. Induction calculations may not have been enabled.")

        # Map GUI variable names to induction data structure names
        induction_var_map = {
            'Amp_nT': 'Amp',
            'phase_deg': 'Phase',
            'Bix_nT': 'Bix_nT',
            'Biy_nT': 'Biy_nT',
            'Biz_nT': 'Biz_nT',
            'Bi1x_nT': 'Bi1x_nT',
            'Bi1y_nT': 'Bi1y_nT',
            'Bi1z_nT': 'Bi1z_nT',
            'Bi1Tot_nT': 'Bi1Tot_nT',
            'rBi1x_nT': 'rBi1x_nT',
            'rBi1y_nT': 'rBi1y_nT',
            'rBi1z_nT': 'rBi1z_nT',
            'rBi1Tot_nT': 'rBi1Tot_nT',
            'iBi1x_nT': 'iBi1x_nT',
            'iBi1y_nT': 'iBi1y_nT',
            'iBi1z_nT': 'iBi1z_nT',
            'iBi1Tot_nT': 'iBi1Tot_nT',
        }

        induction_attr_name = induction_var_map.get(zName, zName)

        if not hasattr(Exploration.induction, induction_attr_name):
            raise ValueError(f"Could not find '{induction_attr_name}' in induction results")

        induction_data = getattr(Exploration.induction, induction_attr_name)

        if induction_data is None:
            raise ValueError(f"Induction data '{induction_attr_name}' is None")

        # Induction data is 3D: (nPeaks, ny, nx)
        # Use the first frequency peak
        if len(induction_data.shape) == 3:
            zData = induction_data[0, :, :]  # First frequency peak
        elif len(induction_data.shape) == 2:
            zData = induction_data  # Already 2D
        else:
            raise ValueError(f"Induction data '{induction_attr_name}' has unexpected shape: {induction_data.shape}")

    elif hasattr(Exploration.base, zName):
        zData = getattr(Exploration.base, zName)
        if zData is None:
            raise ValueError(f"'{zName}' exists but is None")
    else:
        raise ValueError(f"Could not find '{zName}' in exploration results")

    # Ensure 2D shape
    if len(xData.shape) == 1:
        xData, yData = np.meshgrid(xData, yData)

    # Get validity mask
    VALID = Exploration.base.VALID
    if VALID is not None:
        zData = zData.copy()
        zData[~VALID] = np.nan

    # Apply smoothing if requested
    if smoothing and smooth_factor > 1:
        try:
            from PlanetProfile.Utilities.DataManip import smoothGrid
            xData, yData, [zData] = smoothGrid(xData, yData, [zData], smooth_factor)
        except:
            pass  # Skip smoothing if not available

    # Get matplotlib colormap if available
    cmap_name = 'viridis'  # Default
    if FigLbl is not None and hasattr(FigLbl, 'cmap'):
        cmap_name = getattr(FigLbl, 'cmap', 'viridis')

    colorscale = get_matplotlib_colormap_colors(cmap_name)
    if colorscale is None:
        colorscale = 'Viridis'  # Plotly default

    # Calculate valid data range
    z_valid = zData[~np.isnan(zData)]
    if len(z_valid) == 0:
        raise ValueError("No valid data points to plot")

    data_min = np.nanmin(zData)
    data_max = np.nanmax(zData)

    # Check if we should center colormap at zero
    center_at_zero = False
    if FigLbl is not None and hasattr(FigLbl, 'cMapZero'):
        if zName in FigLbl.cMapZero and data_min < 0 and data_max > 0:
            center_at_zero = True
            abs_max = max(abs(data_min), abs(data_max))
            zmin, zmax = -abs_max, abs_max
        else:
            zmin, zmax = data_min, data_max
    else:
        zmin, zmax = data_min, data_max

    # Create heatmap
    fig = go.Figure()

    heatmap = go.Heatmap(
        x=xData[0, :] if len(xData.shape) == 2 else xData,
        y=yData[:, 0] if len(yData.shape) == 2 else yData,
        z=zData,
        colorscale=colorscale,
        zmin=zmin,
        zmax=zmax,
        colorbar=dict(
            title=dict(
                text=zName,
                side='right'
            ),
            thickness=20,
            len=0.7,
        ),
        hovertemplate=(
            f'X: %{{x}}<br>' +
            f'Y: %{{y}}<br>' +
            f'{zName}: %{{z}}<br>' +
            '<extra></extra>'
        )
    )

    fig.add_trace(heatmap)

    # Get axis labels
    if FigLbl is not None:
        xLabel = getattr(FigLbl, 'xLabelExplore', Exploration.xName)
        yLabel = getattr(FigLbl, 'yLabelExplore', Exploration.yName)
        title = getattr(FigLbl, 'explorationTitle', f'{zName} Exploreogram')
    else:
        xLabel = Exploration.xName
        yLabel = Exploration.yName
        title = f'{zName} vs {Exploration.xName} and {Exploration.yName}'

    # Update layout to match matplotlib style
    fig.update_layout(
        title=dict(
            text=title,
            x=0.5,
            xanchor='center',
            font=dict(size=16)
        ),
        xaxis=dict(
            title=xLabel,
            showgrid=True,
            gridcolor='rgba(128, 128, 128, 0.2)',
            zeroline=False,
        ),
        yaxis=dict(
            title=yLabel,
            showgrid=True,
            gridcolor='rgba(128, 128, 128, 0.2)',
            zeroline=False,
        ),
        width=900,
        height=700,
        plot_bgcolor='white',
        font=dict(family='Arial, sans-serif', size=12),
    )

    # Add contour lines if contour data exists
    if hasattr(Exploration, 'contourName') and Exploration.contourName is not None:
        contour_name = Exploration.contourName
        if hasattr(Exploration.base, contour_name):
            contour_data = getattr(Exploration.base, contour_name)

            # Add contour lines
            fig.add_trace(go.Contour(
                x=xData[0, :] if len(xData.shape) == 2 else xData,
                y=yData[:, 0] if len(yData.shape) == 2 else yData,
                z=contour_data,
                showscale=False,
                contours=dict(
                    coloring='lines',
                    showlabels=True,
                    labelfont=dict(size=10, color='black')
                ),
                line=dict(color='black', width=1),
                hoverinfo='skip'
            ))

    return fig


def create_multi_exploreogram_plotly(Exploration, Params, zNames, FigLbl=None):
    """
    Create multi-subplot exploreogram with multiple z-variables.

    Args:
        Exploration: ExplorationResultsStruct with results
        Params: ParamsStruct with configuration
        zNames: List of z-variable names to plot
        FigLbl: Figure labels struct (optional)

    Returns:
        Plotly figure object with subplots
    """

    n_plots = len(zNames)
    rows = int(np.ceil(n_plots / 2))
    cols = 2 if n_plots > 1 else 1

    # Create subplot titles
    subplot_titles = [f'{zName}' for zName in zNames]

    fig = make_subplots(
        rows=rows, cols=cols,
        subplot_titles=subplot_titles,
        vertical_spacing=0.12,
        horizontal_spacing=0.15
    )

    # Get matplotlib colormap
    cmap_name = 'viridis'
    if FigLbl is not None and hasattr(FigLbl, 'cmap'):
        cmap_name = getattr(FigLbl, 'cmap', 'viridis')

    colorscale = get_matplotlib_colormap_colors(cmap_name)
    if colorscale is None:
        colorscale = 'Viridis'

    # Extract common data
    xData = Exploration.xData
    yData = Exploration.yData
    VALID = Exploration.base.VALID

    # Ensure 2D shape
    if len(xData.shape) == 1:
        xData, yData = np.meshgrid(xData, yData)

    # Add each z-variable as a subplot
    for idx, zName in enumerate(zNames):
        row = idx // cols + 1
        col = idx % cols + 1

        # Get z-data
        if hasattr(Exploration.base, zName):
            zData = getattr(Exploration.base, zName).copy()
        else:
            continue

        # Apply validity mask
        if VALID is not None:
            zData[~VALID] = np.nan

        # Calculate range
        data_min = np.nanmin(zData)
        data_max = np.nanmax(zData)

        # Add heatmap
        fig.add_trace(
            go.Heatmap(
                x=xData[0, :],
                y=yData[:, 0],
                z=zData,
                colorscale=colorscale,
                zmin=data_min,
                zmax=data_max,
                colorbar=dict(
                    title=zName,
                    thickness=15,
                    len=0.4,
                    x=1.02 if col == 2 else 0.48,
                    y=1.0 - (row - 0.5) / rows
                ),
                hovertemplate=f'{zName}: %{{z}}<br><extra></extra>'
            ),
            row=row, col=col
        )

        # Update axes
        fig.update_xaxes(title_text=Exploration.xName, row=row, col=col, showgrid=True)
        fig.update_yaxes(title_text=Exploration.yName, row=row, col=col, showgrid=True)

    # Update layout
    fig.update_layout(
        title=dict(
            text=f'Multi-Variable Exploreogram: {Exploration.xName} vs {Exploration.yName}',
            x=0.5,
            xanchor='center'
        ),
        height=400 * rows,
        width=1200,
        showlegend=False,
        plot_bgcolor='white'
    )

    return fig
