"""
Enhanced Plotly plotting for exploreograms that matches matplotlib styling from PlanetProfile.
Includes multi-frequency inductogram support with contour lines.
"""
import re
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
import os


def latex_to_plotly(text):
    """
    Convert matplotlib-style LaTeX labels to Plotly-compatible Unicode/HTML.

    Converts everything — including $...$ math content — to Unicode so that
    labels render correctly without MathJax (which is unreliable in Streamlit).

    Examples:
        'Salinity $w$ ($g\\,kg^{-1}$)' → 'Salinity w (g kg⁻¹)'
        '\\textbf{Europa}' → '<b>Europa</b>'
        '$\\mathrm{Re}\\{B^i\\}$ (nT)' → 'Re{Bⁱ} (nT)'
    """
    if text is None:
        return ''

    # \textbf{...} -> <b>...</b>
    text = re.sub(r'\\textbf\{([^}]*)\}', r'<b>\1</b>', text)
    # \textit{...} -> <i>...</i>
    text = re.sub(r'\\textit\{([^}]*)\}', r'<i>\1</i>', text)

    # \si{...} -> clean unit text (e.g. \si{g\,kg^{-1}} -> g kg⁻¹)
    text = re.sub(r'\\si\{([^}]*)\}', lambda m: _latex_unit_to_text(m.group(1)), text)

    # Strip $...$ delimiters and convert math content to Unicode
    text = re.sub(r'\$([^$]*)\$', lambda m: _latex_math_to_unicode(m.group(1)), text)

    # Clean remaining LaTeX outside math: \, \; \\ etc.
    text = text.replace('\\,', ' ')
    text = text.replace('\\;', ' ')
    text = text.replace('\\.', '')
    # Stray ^{...} and _{...} outside math
    text = re.sub(r'\^\{([^}]*)\}', lambda m: _to_superscript(m.group(1)), text)
    text = re.sub(r'_\{([^}]*)\}', lambda m: _to_subscript(m.group(1)), text)

    return text


def _latex_math_to_unicode(math_content):
    """Convert LaTeX math content (without $ delimiters) to Unicode."""
    s = math_content
    # \mathrm{...} -> {content} (keep braces so _\mathrm{sea} becomes _{sea} for subscript)
    s = re.sub(r'\\mathrm\{([^}]*)\}', r'{\1}', s)
    # \mathbf{...} -> {content}
    s = re.sub(r'\\mathbf\{([^}]*)\}', r'{\1}', s)
    # \text{...} -> {content}
    s = re.sub(r'\\text\{([^}]*)\}', r'{\1}', s)
    # \ce{...} -> {content} (chemistry notation)
    s = re.sub(r'\\ce\{([^}]*)\}', r'{\1}', s)
    # Escaped braces: \{ \} -> literal braces
    s = s.replace('\\{', '{').replace('\\}', '}')
    # Greek letters
    _greek = {
        '\\alpha': '\u03b1', '\\beta': '\u03b2', '\\gamma': '\u03b3',
        '\\delta': '\u03b4', '\\epsilon': '\u03b5', '\\sigma': '\u03c3',
        '\\omega': '\u03c9', '\\Omega': '\u03a9', '\\mu': '\u03bc',
        '\\rho': '\u03c1', '\\theta': '\u03b8', '\\phi': '\u03c6',
        '\\pi': '\u03c0', '\\lambda': '\u03bb', '\\tau': '\u03c4',
    }
    for cmd, char in _greek.items():
        s = s.replace(cmd, char)
    # ^{\circ} or ^\circ -> degree symbol (must come before generic \circ)
    s = s.replace('^{\\circ}', '°')
    s = s.replace('^\\circ', '°')
    # \circ -> degree symbol
    s = s.replace('\\circ', '°')
    # Spacing
    s = s.replace('\\,', ' ')
    s = s.replace('\\;', ' ')
    s = s.replace('\\:', ' ')
    s = s.replace('\\ ', ' ')
    # Superscripts and subscripts
    s = re.sub(r'\^\{([^}]*)\}', lambda m: _to_superscript(m.group(1)), s)
    s = re.sub(r'_\{([^}]*)\}', lambda m: _to_subscript(m.group(1)), s)
    # Single-char super/sub without braces: ^2 -> ², _i -> ᵢ
    s = re.sub(r'\^([0-9a-zA-Z])', lambda m: _to_superscript(m.group(1)), s)
    s = re.sub(r'_([0-9a-zA-Z])', lambda m: _to_subscript(m.group(1)), s)
    # Strip leftover braces from {content} that weren't part of sub/superscripts
    s = re.sub(r'\{([^}]*)\}', r'\1', s)
    return s


# Unicode superscript/subscript maps for common characters
_SUPERSCRIPTS = str.maketrans('0123456789+-=()niab', '⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁿⁱᵃᵇ')
_SUBSCRIPTS = str.maketrans('0123456789+-=()aehijklmnoprstuvx',
                            '₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓ')


def _to_superscript(s):
    return s.translate(_SUPERSCRIPTS)


def _to_subscript(s):
    return s.translate(_SUBSCRIPTS)


def _latex_unit_to_text(unit_str):
    """Convert LaTeX unit string to plain Unicode text.
    e.g. 'g\\,kg^{-1}' -> 'g kg⁻¹'
    """
    s = unit_str
    s = s.replace('\\,', ' ')
    s = s.replace('\\;', ' ')
    s = s.replace('\\cdot', '\u00b7')
    # ^{...} -> superscript unicode
    s = re.sub(r'\^\{([^}]*)\}', lambda m: _to_superscript(m.group(1)), s)
    # _{...} -> subscript unicode
    s = re.sub(r'_\{([^}]*)\}', lambda m: _to_subscript(m.group(1)), s)
    return s


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
    except Exception:
        # Fallback to plotly viridis if matplotlib not available
        return None


def salinity_to_conductivity_axis(salinity_ticks):
    """
    Map salinity tick values (ppt) to conductivity labels (S/m) for a secondary axis.

    Uses the linear approximation sigma ≈ 0.18 * w (S/m per ppt), which is a
    display hint only — not a scientific claim.  Values are rounded to 2
    significant figures so the secondary-axis labels stay compact.

    Args:
        salinity_ticks: array-like of salinity values in ppt

    Returns:
        List of conductivity label strings, or None if the mapping is
        trivial (all-zero salinity) or the input is empty/invalid.
    """
    try:
        ticks = np.asarray(salinity_ticks, dtype=float)
    except (TypeError, ValueError):
        return None

    if ticks.size == 0:
        return None

    # Drop NaN ticks before any arithmetic to avoid ValueError in int(rounded)
    ticks = ticks[~np.isnan(ticks)]
    if ticks.size == 0:
        return None

    sigma_vals = 0.18 * ticks

    # If all values are effectively zero there is nothing useful to show
    if np.all(sigma_vals == 0.0):
        return None

    labels = []
    for v in sigma_vals:
        if v == 0.0:
            labels.append('0')
        else:
            # Round to 2 significant figures
            magnitude = 10 ** np.floor(np.log10(abs(v)))
            rounded = round(v / magnitude, 1) * magnitude
            # Format: drop trailing zeros after decimal
            if rounded == int(rounded):
                labels.append(str(int(rounded)))
            else:
                labels.append(f'{rounded:.2g}')
    return labels


def _build_secondary_axis_layout(x_var, y_var, x1d, y1d):
    """
    Return a dict of extra layout kwargs to attach a secondary axis when one of
    the plot axes is ``wOcean_ppt`` or ``sigmaMean_Sm``.

    The secondary axis is purely cosmetic (tick labels only); it overlays the
    primary axis so zoom/pan keep both axes in sync.

    Args:
        x_var: string name of the x-axis exploration variable
        y_var: string name of the y-axis exploration variable
        x1d:   1-D array of unique x values (used to compute secondary ticks)
        y1d:   1-D array of unique y values (used to compute secondary ticks)

    Returns:
        dict of Plotly layout kwargs (may be empty if no secondary axis applies)
    """
    extra = {}

    # Choose at most ~6 evenly-spaced tick positions from the primary array
    def _sparse_ticks(arr, n_max=6):
        arr = np.asarray(arr, dtype=float)
        if arr.size <= n_max:
            return arr
        idx = np.round(np.linspace(0, arr.size - 1, n_max)).astype(int)
        return arr[idx]

    SALINITY_VAR = 'wOcean_ppt'
    CONDUCTIVITY_VAR = 'sigmaMean_Sm'

    if x_var == SALINITY_VAR:
        ticks = _sparse_ticks(x1d)
        labels = salinity_to_conductivity_axis(ticks)
        if labels is not None:
            extra['xaxis2'] = dict(
                title=dict(text='Mean Conductivity (S/m)', standoff=4),
                overlaying='x',
                side='top',
                matches='x',
                tickmode='array',
                tickvals=list(ticks),
                ticktext=labels,
                showgrid=False,
                zeroline=False,
            )

    elif x_var == CONDUCTIVITY_VAR:
        ticks = _sparse_ticks(x1d)
        # Inverse: w = sigma / 0.18; guard against divide-by-zero
        w_vals = np.where(ticks > 0, ticks / 0.18, 0.0)
        labels = [f'{v:.2g}' if v > 0 else '0' for v in w_vals]
        extra['xaxis2'] = dict(
            title=dict(text='Salinity (ppt)', standoff=4),
            overlaying='x',
            side='top',
            matches='x',
            tickmode='array',
            tickvals=list(ticks),
            ticktext=labels,
            showgrid=False,
            zeroline=False,
        )

    if y_var == SALINITY_VAR:
        ticks = _sparse_ticks(y1d)
        labels = salinity_to_conductivity_axis(ticks)
        if labels is not None:
            extra['yaxis2'] = dict(
                title=dict(text='Mean Conductivity (S/m)', standoff=4),
                overlaying='y',
                side='right',
                matches='y',
                tickmode='array',
                tickvals=list(ticks),
                ticktext=labels,
                showgrid=False,
                zeroline=False,
            )

    elif y_var == CONDUCTIVITY_VAR:
        ticks = _sparse_ticks(y1d)
        w_vals = np.where(ticks > 0, ticks / 0.18, 0.0)
        labels = [f'{v:.2g}' if v > 0 else '0' for v in w_vals]
        extra['yaxis2'] = dict(
            title=dict(text='Salinity (ppt)', standoff=4),
            overlaying='y',
            side='right',
            matches='y',
            tickmode='array',
            tickvals=list(ticks),
            ticktext=labels,
            showgrid=False,
            zeroline=False,
        )

    return extra


def create_exploreogram_plotly(Exploration, Params, FigLbl=None, smoothing=False, smooth_factor=2, use_contours=True):
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

    # Get axis labels and convert LaTeX to Plotly format (before creating traces)
    if FigLbl is not None:
        xLabel = latex_to_plotly(getattr(FigLbl, 'xLabelExplore', Exploration.xName))
        yLabel = latex_to_plotly(getattr(FigLbl, 'yLabelExplore', Exploration.yName))
        zLabel = latex_to_plotly(getattr(FigLbl, 'cbarLabelExplore', zName))
        title = latex_to_plotly(getattr(FigLbl, 'explorationTitle', f'{zName} Exploreogram'))
    else:
        xLabel = Exploration.xName
        yLabel = Exploration.yName
        zLabel = zName
        title = f'{zName} vs {Exploration.xName} and {Exploration.yName}'

    # Create heatmap
    fig = go.Figure()

    x1d = xData[:, 0] if len(xData.shape) == 2 else xData
    y1d = yData[0, :] if len(yData.shape) == 2 else yData

    heatmap = go.Heatmap(
        x=x1d,
        y=y1d,
        z=zData.T,
        colorscale=colorscale,
        zmin=zmin,
        zmax=zmax,
        colorbar=dict(
            title=dict(
                text=zLabel,
                side='right'
            ),
            thickness=20,
            len=0.7,
        ),
        hovertemplate=(
            f'X: %{{x}}<br>' +
            f'Y: %{{y}}<br>' +
            f'{zLabel}: %{{z}}<br>' +
            '<extra></extra>'
        )
    )

    fig.add_trace(heatmap)

    # Add contour lines overlay on the heatmap if enabled
    if use_contours and zmax != zmin:
        n_contours = 8
        fig.add_trace(go.Contour(
            x=x1d, y=y1d, z=zData.T,
            showscale=False,
            contours=dict(
                start=zmin, end=zmax,
                size=(zmax - zmin) / n_contours,
                coloring='lines',
                showlabels=True,
                labelfont=dict(size=10, color='white'),
            ),
            line=dict(color='white', width=1),
            hoverinfo='skip',
        ))

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

    # Add contour lines from a separate contour variable if specified
    if use_contours and hasattr(Exploration, 'contourName') and Exploration.contourName is not None:
        contour_name = Exploration.contourName
        if hasattr(Exploration.base, contour_name):
            contour_data = getattr(Exploration.base, contour_name)

            # Add contour lines
            fig.add_trace(go.Contour(
                x=xData[:, 0] if len(xData.shape) == 2 else xData,
                y=yData[0, :] if len(yData.shape) == 2 else yData,
                z=contour_data.T,
                showscale=False,
                contours=dict(
                    coloring='lines',
                    showlabels=True,
                    labelfont=dict(size=10, color='black')
                ),
                line=dict(color='black', width=1),
                hoverinfo='skip'
            ))

    # Attach secondary axis when one axis is wOcean_ppt or sigmaMean_Sm
    secondary_axis_kwargs = _build_secondary_axis_layout(
        Exploration.xName, Exploration.yName, x1d, y1d
    )
    if secondary_axis_kwargs:
        fig.update_layout(**secondary_axis_kwargs)

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

    # Create subplot titles — convert any LaTeX to Unicode
    subplot_titles = [latex_to_plotly(zName) for zName in zNames]

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
                x=xData[:, 0],
                y=yData[0, :],
                z=zData.T,
                colorscale=colorscale,
                zmin=data_min,
                zmax=data_max,
                colorbar=dict(
                    title=latex_to_plotly(zName),
                    thickness=15,
                    len=0.4,
                    x=1.02 if col == 2 else 0.48,
                    y=1.0 - (row - 0.5) / rows
                ),
                hovertemplate=f'{latex_to_plotly(zName)}: %{{z}}<br><extra></extra>'
            ),
            row=row, col=col
        )

        # Update axes — convert labels
        x_title = latex_to_plotly(Exploration.xName)
        y_title = latex_to_plotly(Exploration.yName)
        if FigLbl is not None:
            x_title = latex_to_plotly(getattr(FigLbl, 'xLabelExplore', Exploration.xName))
            y_title = latex_to_plotly(getattr(FigLbl, 'yLabelExplore', Exploration.yName))
        fig.update_xaxes(title_text=x_title, row=row, col=col, showgrid=True)
        fig.update_yaxes(title_text=y_title, row=row, col=col, showgrid=True)

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


def create_inductogram_plotly(Exploration, Params, FigLbl=None,
                              display_mode='real_imaginary', use_contours=True):
    """
    Create multi-frequency inductogram with contour lines following
    published literature conventions.

    Each selected frequency gets a row of panels showing either:
    - Real + Imaginary components (display_mode='real_imaginary')
    - Amplitude + Phase (display_mode='amplitude_phase')

    Args:
        Exploration: ExplorationResultsStruct with induction results
        Params: ParamsStruct with configuration
        FigLbl: Figure labels struct (optional)
        display_mode: 'real_imaginary' or 'amplitude_phase'
        use_contours: If True, use contour lines (paper convention).
                      If False, use heatmap coloring.

    Returns:
        Plotly figure object
    """
    induction = Exploration.induction
    if induction is None:
        raise ValueError("No induction data available in Exploration results")

    nPeaks = induction.nPeaks
    if nPeaks is None or nPeaks == 0:
        raise ValueError("No frequency peaks found in induction data")

    # Get frequency labels
    if hasattr(induction, 'calcedExc') and induction.calcedExc is not None:
        freq_names = list(induction.calcedExc)
    else:
        freq_names = [f'Peak {i}' for i in range(nPeaks)]

    if hasattr(induction, 'Texc_hr') and induction.Texc_hr is not None:
        periods = induction.Texc_hr
    else:
        periods = [None] * nPeaks

    # Extract common grid data
    xData = Exploration.xData
    yData = Exploration.yData
    VALID = Exploration.base.VALID

    if len(xData.shape) == 1:
        xData, yData = np.meshgrid(xData, yData)

    x1d = xData[:, 0]
    y1d = yData[0, :]

    # Get axis labels
    if FigLbl is not None:
        xLabel = latex_to_plotly(getattr(FigLbl, 'xLabelExplore', Exploration.xName))
        yLabel = latex_to_plotly(getattr(FigLbl, 'yLabelExplore', Exploration.yName))
    else:
        xLabel = Exploration.xName
        yLabel = Exploration.yName

    # Determine subplot layout: nPeaks rows x 2 columns
    rows = int(nPeaks)
    cols = 2

    # Build subplot titles
    subplot_titles = []
    for iPeak in range(nPeaks):
        freq_label = freq_names[iPeak] if iPeak < len(freq_names) else f'Peak {iPeak}'
        T_str = f' (T = {periods[iPeak]:.2f} hr)' if periods[iPeak] is not None else ''

        if display_mode == 'real_imaginary':
            subplot_titles.append(f'Re(Bi) — {freq_label}{T_str}')
            subplot_titles.append(f'Im(Bi) — {freq_label}{T_str}')
        else:
            subplot_titles.append(f'Amplitude — {freq_label}{T_str}')
            subplot_titles.append(f'Phase — {freq_label}{T_str}')

    # Wide horizontal gap for colorbars between columns; generous vertical for titles
    # Gap sized so left-column colorbar label + 3× ylabel cap-height clears right ylabel
    v_spacing = 0.15 / rows + 0.06 if rows > 1 else 0.15
    fig = make_subplots(
        rows=rows, cols=cols,
        subplot_titles=subplot_titles,
        vertical_spacing=v_spacing,
        horizontal_spacing=0.36,
        column_widths=[0.45, 0.45],
    )

    # Get colorscale
    colorscale_pos = get_matplotlib_colormap_colors('viridis')
    if colorscale_pos is None:
        colorscale_pos = 'Viridis'

    # For diverging data (real/imaginary that can be negative), use RdBu
    colorscale_div = 'RdBu_r'

    for iPeak in range(nPeaks):
        row = iPeak + 1

        if display_mode == 'real_imaginary':
            # Real component — use rBi1Tot_nT (total induced field, real part)
            re_data = _get_induction_slice(induction, 'rBi1Tot_nT', iPeak, VALID)
            im_data = _get_induction_slice(induction, 'iBi1Tot_nT', iPeak, VALID)
            left_label = 'Re{B\u2071} (nT)'
            right_label = 'Im{B\u2071} (nT)'
            left_data, right_data = re_data, im_data
            left_cscale = colorscale_div
            right_cscale = colorscale_div
        else:
            # Amplitude + Phase
            amp_data = _get_induction_slice(induction, 'Amp', iPeak, VALID)
            phase_data = _get_induction_slice(induction, 'Phase', iPeak, VALID)
            left_label = 'Amplitude (nT)'
            right_label = 'Phase (\u00b0)'
            left_data, right_data = amp_data, phase_data
            left_cscale = colorscale_pos
            right_cscale = colorscale_pos

        # Add left panel
        _add_inductogram_panel(
            fig, x1d, y1d, left_data, row, 1, left_label,
            colorscale=left_cscale, use_contours=use_contours,
            nRows=rows, iPeak=iPeak
        )

        # Add right panel
        _add_inductogram_panel(
            fig, x1d, y1d, right_data, row, 2, right_label,
            colorscale=right_cscale, use_contours=use_contours,
            nRows=rows, iPeak=iPeak
        )

        # Set axis labels — only bottom row gets x-labels, only left column gets y-label
        is_bottom = (row == rows)
        fig.update_xaxes(title_text=xLabel if is_bottom else '', row=row, col=1)
        fig.update_xaxes(title_text=xLabel if is_bottom else '', row=row, col=2)
        fig.update_yaxes(title_text=yLabel, row=row, col=1)
        fig.update_yaxes(title_text=yLabel, row=row, col=2)

    # Build title
    bodyname = getattr(Exploration, 'bodyname', '')
    if not bodyname and hasattr(Exploration, 'xName'):
        bodyname = ''
    mode_str = 'Real + Imaginary' if display_mode == 'real_imaginary' else 'Amplitude + Phase'
    title_text = f'{bodyname} Inductogram ({mode_str})' if bodyname else f'Inductogram ({mode_str})'

    # Reduce subplot title font sizes and shift upward to avoid overlapping figure border
    for annotation in fig['layout']['annotations']:
        annotation['font'] = dict(size=10)
        # Shift title upward by adding offset to y position
        if 'y' in annotation:
            annotation['y'] = annotation['y'] + 0.015

    # Push overall title up by 3× its cap-height so it clears subplot titles
    fig_height = max(480, 420 * rows)
    title_font_size = 12
    title_lift_px = 3 * 0.7 * title_font_size  # 3× cap-height
    title_y = 0.98 + title_lift_px / fig_height
    top_margin = 50 + int(title_lift_px)

    fig.update_layout(
        title=dict(text=title_text, x=0.5, xanchor='center',
                   font=dict(size=title_font_size), y=title_y, yanchor='top'),
        height=fig_height,
        width=1300,
        showlegend=False,
        plot_bgcolor='white',
        font=dict(family='Arial, sans-serif', size=10),
        margin=dict(l=60, r=100, t=top_margin, b=40),
    )

    return fig


def _get_induction_slice(induction, attr_name, iPeak, VALID):
    """Extract a 2D slice from a 3D induction array for a given peak index."""
    data = getattr(induction, attr_name, None)
    if data is None:
        return None

    # Handle complex arrays — take real part if complex
    if np.iscomplexobj(data):
        data = np.real(data)

    if len(data.shape) == 3:
        sliced = data[iPeak, :, :].copy()
    elif len(data.shape) == 2:
        sliced = data.copy()
    else:
        return None

    # Apply validity mask
    if VALID is not None:
        sliced[~VALID] = np.nan

    return sliced


def _add_inductogram_panel(fig, x1d, y1d, zData, row, col, label,
                           colorscale='Viridis', use_contours=True,
                           nRows=1, iPeak=0):
    """Add a single inductogram panel (contour or heatmap) to a subplot."""
    if zData is None:
        return

    # Ensure label is clean Unicode (no raw LaTeX)
    label = latex_to_plotly(label)

    z_valid = zData[~np.isnan(zData)]
    if len(z_valid) == 0:
        return

    zmin = np.nanmin(zData)
    zmax = np.nanmax(zData)

    # Center diverging colormaps at zero when data spans both signs
    if isinstance(colorscale, str) and 'RdBu' in colorscale:
        if zmin < 0 and zmax > 0:
            abs_max = max(abs(zmin), abs(zmax))
            zmin, zmax = -abs_max, abs_max

    # Position colorbar to avoid overlap in multi-row layout
    cbar_len = max(0.12, 0.65 / nRows)
    cbar_y = 1.0 - (iPeak + 0.5) / nRows
    # Left colorbars sit just right of left column; right colorbars past right edge
    cbar_x = 0.36 if col == 1 else 1.04
    cbar_dict = dict(
        title=dict(text=label, side='right', font=dict(size=8)),
        thickness=8,
        len=cbar_len,
        y=cbar_y,
        x=cbar_x,
        yanchor='middle',
    )

    if use_contours:
        # Contour lines — matches paper convention
        n_contours = 12
        fig.add_trace(
            go.Contour(
                x=x1d,
                y=y1d,
                z=zData.T,
                colorscale=colorscale,
                contours=dict(
                    start=zmin,
                    end=zmax,
                    size=(zmax - zmin) / n_contours if zmax != zmin else 1,
                    showlabels=True,
                    labelfont=dict(size=9, color='black'),
                ),
                line=dict(width=1.5),
                colorbar=cbar_dict,
                hovertemplate=(
                    f'X: %{{x:.3g}}<br>'
                    f'Y: %{{y:.3g}}<br>'
                    f'{label}: %{{z:.4g}}<br>'
                    '<extra></extra>'
                ),
            ),
            row=row, col=col
        )
    else:
        # Heatmap fallback
        fig.add_trace(
            go.Heatmap(
                x=x1d,
                y=y1d,
                z=zData.T,
                colorscale=colorscale,
                zmin=zmin,
                zmax=zmax,
                colorbar=cbar_dict,
                hovertemplate=(
                    f'X: %{{x:.3g}}<br>'
                    f'Y: %{{y:.3g}}<br>'
                    f'{label}: %{{z:.4g}}<br>'
                    '<extra></extra>'
                ),
            ),
            row=row, col=col
        )
