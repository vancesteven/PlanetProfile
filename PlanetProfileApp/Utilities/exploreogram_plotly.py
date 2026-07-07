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

# Authoritative list of magnetic induction z-variable names used in the
# explore-o-gram GUI.  Import INDUCTION_Z_VARS from this module; do not
# maintain a separate copy in Exploreogram.py.
INDUCTION_Z_VARS = [
    'Amp_nT', 'Bix_nT', 'Biy_nT', 'Biz_nT', 'phase_deg',
    'Bi1x_nT', 'Bi1y_nT', 'Bi1z_nT', 'Bi1Tot_nT',
    'rBi1x_nT', 'rBi1y_nT', 'rBi1z_nT', 'rBi1Tot_nT',
    'iBi1x_nT', 'iBi1y_nT', 'iBi1z_nT', 'iBi1Tot_nT',
]


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
    # \overline{x} and \bar{x}: Greek was already substituted above, so
    # m.group(1) may contain Unicode (e.g. 'σ'). Apply one combining overline.
    def _with_overline(inner):
        return ''.join(c + '\u0305' for c in inner)
    s = re.sub(r'\\overline\{([^}]*)\}', lambda m: _with_overline(m.group(1)), s)
    s = re.sub(r'\\bar\{([^}]*)\}', lambda m: _with_overline(m.group(1)), s)
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
        import matplotlib as mpl
        from matplotlib.colors import LinearSegmentedColormap

        cmap = mpl.colormaps[cmap_name]
        colors = [cmap(i / (n_colors - 1)) for i in range(n_colors)]
        # Convert RGBA to hex
        hex_colors = ['#{:02x}{:02x}{:02x}'.format(
            int(r*255), int(g*255), int(b*255)
        ) for r, g, b, a in colors]
        return hex_colors
    except Exception:
        # Fallback to plotly viridis if matplotlib not available
        return None


def _interp_secondary_labels(primary_ticks, primary_data, secondary_data):
    """Interpolate secondary-axis tick labels from PP-computed (primary, secondary) pairs.

    Used to attach a salinity↔conductivity secondary axis: when the primary
    x-axis is wOcean_ppt, ``primary_data`` is the 1-D salinity values and
    ``secondary_data`` is the corresponding mean ocean conductivities (averaged
    across the orthogonal grid axis to produce a 1-D curve).

    Returns None (disabling the secondary axis) when:
      - either array is empty after NaN removal,
      - the two arrays differ in length,
      - primary_data is non-monotone (e.g. MgSO4 σ↔w near eutectic — np.interp
        would produce wrong labels in that case).

    Args:
        primary_ticks:  array-like of primary-axis tick positions to label.
        primary_data:   1-D array of PP-computed primary-axis values.
        secondary_data: 1-D array of PP-computed secondary-axis values,
                        same length as primary_data.

    Returns:
        list[str] of formatted tick labels, or None to omit the secondary axis.
    """
    try:
        p_ticks = np.asarray(primary_ticks, dtype=float).ravel()
        p_data  = np.asarray(primary_data,  dtype=float).ravel()
        s_data  = np.asarray(secondary_data, dtype=float).ravel()
    except (TypeError, ValueError):
        return None

    if p_data.size == 0 or s_data.size == 0 or p_data.size != s_data.size:
        return None

    finite = np.isfinite(p_data) & np.isfinite(s_data)
    p_data = p_data[finite]
    s_data = s_data[finite]
    if p_data.size < 2:
        return None

    order = np.argsort(p_data)
    p_sorted = p_data[order]
    s_sorted = s_data[order]

    if not np.all(np.diff(p_sorted) > 0):
        return None

    diffs = np.diff(s_sorted)
    if not (np.all(diffs > 0) or np.all(diffs < 0)):
        return None  # non-monotone σ↔w (e.g. MgSO4 near eutectic) — disable

    p_ticks_clean = p_ticks[np.isfinite(p_ticks)]
    if p_ticks_clean.size == 0:
        return None

    s_at_ticks = np.interp(p_ticks_clean, p_sorted, s_sorted)

    labels = []
    for v in s_at_ticks:
        if v == 0.0:
            labels.append('0')
        elif abs(v) < 0.01 or abs(v) >= 1000:
            labels.append(f'{v:.2g}')
        else:
            labels.append(f'{v:.3g}')
    return labels


def _attach_secondary_axis(fig, Exploration, x1d, y1d):
    """Attach a synchronized salinity↔conductivity secondary axis to a single-panel figure.

    Activates when Exploration.xName or Exploration.yName is 'wOcean_ppt' or
    'sigmaMean_Sm'.  Secondary labels are derived from PP-computed
    (w, σ) pairs via np.nanmean across the orthogonal grid axis — this
    gives a physically defensible 1-D σ-vs-w curve without a hard-coded
    approximation.  Silently skips when the σ↔w relation is non-monotone
    (see _interp_secondary_labels).

    Args:
        fig:         Plotly Figure (single-panel exploreogram).
        Exploration: ExplorationResultsStruct (.xName, .yName, .base).
        x1d:         1-D primary x-axis values used in the heatmap.
        y1d:         1-D primary y-axis values used in the heatmap.
    """
    SAL = 'wOcean_ppt'
    SIG = 'sigmaMean_Sm'

    base = getattr(Exploration, 'base', None)
    if base is None:
        return
    sigma_2d = getattr(base, 'sigmaMean_Sm', None)
    w_2d     = getattr(base, 'wOcean_ppt',  None)
    if sigma_2d is None or w_2d is None:
        return

    sigma_2d = np.asarray(sigma_2d, dtype=float)
    w_2d     = np.asarray(w_2d,     dtype=float)

    def _sparse(arr, n_max=6):
        arr = np.asarray(arr, dtype=float)
        finite = arr[np.isfinite(arr)]
        if finite.size <= n_max:
            return finite
        idx = np.round(np.linspace(0, finite.size - 1, n_max)).astype(int)
        return finite[idx]

    x_var = getattr(Exploration, 'xName', None)
    y_var = getattr(Exploration, 'yName', None)

    layout_updates = {}

    if x_var in (SAL, SIG):
        # Mean across y (axis=1) → 1-D of length nx, matching x1d
        with np.errstate(invalid='ignore'):
            sigma_x = np.nanmean(sigma_2d, axis=1)
            w_x     = np.nanmean(w_2d,     axis=1)
        if x_var == SAL:
            primary_data, secondary_data = w_x, sigma_x
            sec_title = 'Mean Conductivity (S/m)'
        else:
            primary_data, secondary_data = sigma_x, w_x
            sec_title = 'Salinity (ppt)'
        ticks = _sparse(x1d)
        labels = _interp_secondary_labels(ticks, primary_data, secondary_data)
        if labels is not None:
            layout_updates['xaxis2'] = dict(
                title=dict(text=sec_title, standoff=4),
                overlaying='x', side='top', matches='x',
                tickmode='array',
                tickvals=list(ticks),
                ticktext=labels,
                showgrid=False, zeroline=False,
            )

    if y_var in (SAL, SIG):
        # Mean across x (axis=0) → 1-D of length ny, matching y1d
        with np.errstate(invalid='ignore'):
            sigma_y = np.nanmean(sigma_2d, axis=0)
            w_y     = np.nanmean(w_2d,     axis=0)
        if y_var == SAL:
            primary_data, secondary_data = w_y, sigma_y
            sec_title = 'Mean Conductivity (S/m)'
        else:
            primary_data, secondary_data = sigma_y, w_y
            sec_title = 'Salinity (ppt)'
        ticks = _sparse(y1d)
        labels = _interp_secondary_labels(ticks, primary_data, secondary_data)
        if labels is not None:
            layout_updates['yaxis2'] = dict(
                title=dict(text=sec_title, standoff=4),
                overlaying='y', side='right', matches='y',
                tickmode='array',
                tickvals=list(ticks),
                ticktext=labels,
                showgrid=False, zeroline=False,
            )

    if layout_updates:
        fig.update_layout(**layout_updates)


def create_exploreogram_plotly(Exploration, Params, FigLbl=None, smoothing=False, smooth_factor=2, use_contours=True,
                               x_log=False, y_log=False, use_nT_amplitudes=True,
                               show_salinity_axis=True, induct_component='Total',
                               x_derived=None, y_derived=None,
                               x_derived_label=None, y_derived_label=None):
    """
    Create Plotly version of exploreogram that matches matplotlib styling.

    Args:
        Exploration: ExplorationResultsStruct with results
        Params: ParamsStruct with configuration
        FigLbl: Figure labels struct (optional, for advanced styling)
        smoothing: Whether to apply smoothing
        smooth_factor: Smoothing factor if smoothing enabled
        use_nT_amplitudes: If True (default), show induced-field amplitudes in nT
            using induction.Bi1Tot_nT (modulus); if False, use the dimensionless
            Ae stored in induction.Amp.
        induct_component: Which induced-field component to display when
            zName == 'Amp_nT' and use_nT_amplitudes is True.  One of
            'Total', 'Bx', 'By', 'Bz'.  Defaults to 'Total'.

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
    if zName in INDUCTION_Z_VARS:
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

        # Amplitude special case: when use_nT_amplitudes is True (default), show
        # the chosen per-component amplitude in nT rather than the dimensionless
        # modulus stored in induction.Amp.  The user can select the component
        # (Total, Bx, By, Bz) via induct_component; 'Total' uses Bi1Tot_nT,
        # the others use the corresponding per-axis array.
        if zName == 'Amp_nT' and use_nT_amplitudes:
            _component_map = {
                'Total': 'Bi1Tot_nT',
                'Bx':    'Bi1x_nT',
                'By':    'Bi1y_nT',
                'Bz':    'Bi1z_nT',
            }
            _attr = _component_map.get(induct_component, 'Bi1Tot_nT')
            _arr = getattr(Exploration.induction, _attr, None)
            if _arr is None:
                raise ValueError(f"{_attr} not available; cannot render nT amplitude for component '{induct_component}'")
            induction_data = np.abs(_arr)
            induction_attr_name = f'|{_attr}| (nT)'  # informational only
        else:
            if not hasattr(Exploration.induction, induction_attr_name):
                raise ValueError(f"Could not find '{induction_attr_name}' in induction results")

            induction_data = getattr(Exploration.induction, induction_attr_name)

            if induction_data is None:
                raise ValueError(f"Induction data '{induction_attr_name}' is None")

        # Induction data is 3D: (nPeaks, nx, ny)
        # Use the first frequency peak
        if len(induction_data.shape) == 3:
            _nPeaks, _nx, _ny = induction_data.shape
            assert _nx == len(xData) and _ny == len(yData), (
                f"Induction data shape {induction_data.shape} inconsistent with "
                f"xData len={len(xData)}, yData len={len(yData)}"
            )
            zData = induction_data[0, :, :]  # First frequency peak, shape (nx, ny)
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
        # defensive fallback for 1D input; production path uses 2D (nx, ny) arrays from ResultsIO
        x1d_pre = xData.copy()
        y1d_pre = yData.copy()
        xData, yData = np.meshgrid(xData, yData)
        # meshgrid convention: shape (ny, nx) — row=y, col=x
        assert xData.shape[0] == len(y1d_pre) and xData.shape[1] == len(x1d_pre), \
            f"meshgrid shape {xData.shape} inconsistent with x1d={len(x1d_pre)}, y1d={len(y1d_pre)}"

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

    # Override colorbar label when amplitude is shown — match the actual data unit
    if zName == 'Amp_nT':
        if use_nT_amplitudes:
            _comp_label_map = {
                'Total': '|Bi| (nT)',
                'Bx':    'Bi_x (nT)',
                'By':    'Bi_y (nT)',
                'Bz':    'Bi_z (nT)',
            }
            zLabel = _comp_label_map.get(induct_component, '|Bi| (nT)')
        else:
            zLabel = 'Amplitude Ae'  # Ae (dimensionless)

    # Create heatmap
    fig = go.Figure()

    # xData/yData shape (nx, ny): col 0 has all unique x; row 0 has all unique y
    x1d = xData[:, 0] if len(xData.shape) == 2 else xData
    y1d = yData[0, :] if len(yData.shape) == 2 else yData

    # Derived-axis tick-label substitution.
    # If a derived 2D array is provided (e.g. sigmaMean_Sm mapped from wOcean_ppt),
    # collapse it to a 1D mean and use it for axis tick labels while keeping the
    # rectilinear driver grid geometry intact.
    x_tickvals  = None
    x_ticklabels = None
    y_tickvals  = None
    y_ticklabels = None

    # 1D derived curves stored for secondary-axis mapping below.
    _x_derived_1d = None  # conductivity values at each x grid point (length nx)
    _y_derived_1d = None  # conductivity values at each y grid point (length ny)

    def _near_monotone(arr):
        """True if arr is monotone or has ≤5% violations (handles low-salinity plateaus)."""
        diffs = np.diff(arr)
        finite = diffs[np.isfinite(diffs)]
        if len(finite) == 0:
            return False
        n = len(diffs)
        return (
            np.all(diffs >= 0) or np.all(diffs <= 0) or
            (np.sum(diffs < 0) / n < 0.05) or
            (np.sum(diffs > 0) / n < 0.05)
        )

    if x_derived is not None:
        # mean across y dimension (axis=1) → length nx, matching x1d
        _x_ticks = np.nanmean(x_derived, axis=1)
        if len(_x_ticks) == len(x1d):
            if _near_monotone(_x_ticks):
                # Use DERIVED coordinates as axis positions so the axis extent matches
                # the derived quantity range (e.g. 0–100 S/m), not the driver range.
                x_tickvals    = _x_ticks
                x_ticklabels  = [f'{v:.3g}' for v in _x_ticks]
                _x_derived_1d = _x_ticks  # save for secondary-axis mapping
            else:
                import warnings
                warnings.warn(
                    f"Derived x-axis '{x_derived_label}' is non-monotone after averaging; "
                    "using driver tick values instead."
                )

    if y_derived is not None:
        # mean across x dimension (axis=0) → length ny, matching y1d
        _y_ticks = np.nanmean(y_derived, axis=0)
        if len(_y_ticks) == len(y1d):
            if _near_monotone(_y_ticks):
                # Use DERIVED coordinates as axis positions.
                y_tickvals    = _y_ticks
                y_ticklabels  = [f'{v:.3g}' for v in _y_ticks]
                _y_derived_1d = _y_ticks  # save for secondary-axis mapping
            else:
                import warnings
                warnings.warn(
                    f"Derived y-axis '{y_derived_label}' is non-monotone after averaging; "
                    "using driver tick values instead."
                )

    # When a derived 1D axis is available, position heatmap cells at derived coordinates
    # so the visual extent matches the derived quantity range (e.g. conductivity 0–100 S/m).
    _hmap_x = _x_derived_1d if _x_derived_1d is not None else x1d
    _hmap_y = _y_derived_1d if _y_derived_1d is not None else y1d

    heatmap = go.Heatmap(
        x=_hmap_x,
        y=_hmap_y,
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
            x=_hmap_x, y=_hmap_y, z=zData.T,
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
            type='log' if x_log else 'linear',
            showgrid=True,
            gridcolor='rgba(128, 128, 128, 0.2)',
            zeroline=False,
        ),
        yaxis=dict(
            title=yLabel,
            type='log' if y_log else 'linear',
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
            # xData/yData shape (nx, ny): col 0 → unique x; row 0 → unique y
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

    # Apply derived-axis tick overrides and title substitution.
    # These replace driver tick values with derived quantity values (e.g. σ in S/m
    # instead of salinity in ppt) while keeping the rectilinear grid geometry.
    _DERIVED_TITLES = {
        'sigmaMean_Sm':    'Mean Conductivity (S/m)',
        'D_km':            'Ocean Thickness (km)',
        'rhoSilMean_kgm3': 'Mean Silicate Density (kg/m³)',
    }

    xaxis_update = {}
    if x_tickvals is not None:
        xaxis_update.update(dict(
            tickmode='array',
            tickvals=list(x_tickvals),
            ticktext=x_ticklabels,
            type='linear',
            range=[float(np.nanmin(x_tickvals)), float(np.nanmax(x_tickvals))],
        ))
    if x_derived_label is not None:
        xaxis_update['title'] = _DERIVED_TITLES.get(x_derived_label, x_derived_label)

    yaxis_update = {}
    if y_tickvals is not None:
        yaxis_update.update(dict(
            tickmode='array',
            tickvals=list(y_tickvals),
            ticktext=y_ticklabels,
            type='linear',
            range=[float(np.nanmin(y_tickvals)), float(np.nanmax(y_tickvals))],
        ))
    if y_derived_label is not None:
        yaxis_update['title'] = _DERIVED_TITLES.get(y_derived_label, y_derived_label)

    if xaxis_update:
        fig.update_layout(xaxis=xaxis_update)
    if yaxis_update:
        fig.update_layout(yaxis=yaxis_update)

    # Secondary salinity↔conductivity axis.
    # When the primary axis already shows conductivity (sigmaMean_Sm via derived-tick
    # labels), and show_salinity_axis is checked, the secondary axis should show the
    # corresponding salinity values — which are exactly the driver tick positions in
    # x1d/y1d (those ARE the wOcean_ppt grid values since xData is never overwritten).
    # For the reverse direction (primary=wOcean_ppt), _attach_secondary_axis handles it.
    _SIG = 'sigmaMean_Sm'

    def _sparse_idx(n, n_max=6):
        return np.round(np.linspace(0, n - 1, min(n_max, n))).astype(int)

    if show_salinity_axis:
        if x_derived_label == _SIG and _x_derived_1d is not None:
            # Primary axis is now in conductivity space (_x_derived_1d).
            # Secondary axis shows salinity (x1d) at the matching conductivity positions.
            _idx = _sparse_idx(len(_x_derived_1d))
            _sig_ticks = [float(_x_derived_1d[i]) for i in _idx]
            _sal_at_sig = list(np.interp(_sig_ticks, _x_derived_1d, x1d))
            fig.update_layout(xaxis2=dict(
                title=dict(text='Salinity (ppt)', standoff=4),
                overlaying='x', side='top', matches='x',
                type='linear',
                anchor='y',
                tickmode='array',
                tickvals=_sig_ticks,
                ticktext=[f'{v:.3g}' for v in _sal_at_sig],
                showgrid=False, zeroline=False,
            ))
            # Plotly only renders axes that have at least one assigned trace.
            fig.add_trace(go.Scatter(
                x=[_sig_ticks[0], _sig_ticks[-1]],
                y=[None, None],
                mode='markers',
                marker=dict(opacity=0, size=1),
                showlegend=False,
                xaxis='x2',
                yaxis='y',
                hoverinfo='skip',
            ))
        if y_derived_label == _SIG and _y_derived_1d is not None:
            _idx = _sparse_idx(len(_y_derived_1d))
            _sig_ticks = [float(_y_derived_1d[i]) for i in _idx]
            _sal_at_sig = list(np.interp(_sig_ticks, _y_derived_1d, y1d))
            fig.update_layout(yaxis2=dict(
                title=dict(text='Salinity (ppt)', standoff=4),
                overlaying='y', side='right', matches='y',
                type='linear',
                anchor='x',
                tickmode='array',
                tickvals=_sig_ticks,
                ticktext=[f'{v:.3g}' for v in _sal_at_sig],
                showgrid=False, zeroline=False,
            ))
            fig.add_trace(go.Scatter(
                x=[None, None],
                y=[_sig_ticks[0], _sig_ticks[-1]],
                mode='markers',
                marker=dict(opacity=0, size=1),
                showlegend=False,
                xaxis='x',
                yaxis='y2',
                hoverinfo='skip',
            ))
        # _attach_secondary_axis handles salinity↔conductivity when neither side
        # is a derived-axis variable.  When a derived label IS present the
        # xaxis2/yaxis2 dummy-scatter path above already drew that secondary.
        # But if only ONE side is derived, the other side may still need a
        # salinity↔conductivity secondary — so only skip when BOTH are derived.
        if not (x_derived_label is not None and y_derived_label is not None):
            _attach_secondary_axis(fig, Exploration, x1d, y1d)

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
        x1d_pre = xData.copy()
        y1d_pre = yData.copy()
        xData, yData = np.meshgrid(xData, yData)
        # meshgrid convention: shape (ny, nx) — row=y, col=x
        assert xData.shape[0] == len(y1d_pre) and xData.shape[1] == len(x1d_pre), \
            f"meshgrid shape {xData.shape} inconsistent with x1d={len(x1d_pre)}, y1d={len(y1d_pre)}"

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
        # xData/yData shape (nx, ny): col 0 → unique x; row 0 → unique y
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
                              display_mode='real_imaginary', use_contours=True,
                              use_nT_amplitudes=True,
                              x_log=False, y_log=False,
                              show_salinity_axis=True,
                              induct_component='Total'):
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
        use_nT_amplitudes: If True (default), amplitudes / Re / Im are shown in
            nT (using rBi1Tot_nT, iBi1Tot_nT, |Bi1Tot_nT|).  If False, the
            dimensionless Ae representation is used (Amp + phase, or scaled
            Re/Im built from Ae·exp(iφ)).
        x_log: If True, use logarithmic x-axis on each subplot.
        y_log: If True, use logarithmic y-axis on each subplot.
        show_salinity_axis: Accepted for signature parity with
            create_exploreogram_plotly. Secondary axis attachment is not yet
            implemented for multi-panel inductogram subplots; this parameter
            is currently ignored.

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
        # defensive fallback for 1D input; production path uses 2D (nx, ny) arrays from ResultsIO
        x1d_pre = xData.copy()
        y1d_pre = yData.copy()
        xData, yData = np.meshgrid(xData, yData)
        # meshgrid convention: shape (ny, nx) — row=y, col=x
        assert xData.shape[0] == len(y1d_pre) and xData.shape[1] == len(x1d_pre), \
            f"meshgrid shape {xData.shape} inconsistent with x1d={len(x1d_pre)}, y1d={len(y1d_pre)}"

    # xData/yData shape (nx, ny): col 0 → unique x; row 0 → unique y
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

        # Map induct_component selector to field suffixes
        _comp_suffix = {'Total': 'Tot', 'Bx': 'x', 'By': 'y', 'Bz': 'z'}.get(induct_component, 'Tot')
        _comp_label  = '' if induct_component == 'Total' else f'_{induct_component}'

        if display_mode == 'real_imaginary':
            re_data = _get_induction_slice(induction, f'rBi1{_comp_suffix}_nT', iPeak, VALID)
            im_data = _get_induction_slice(induction, f'iBi1{_comp_suffix}_nT', iPeak, VALID)
            left_label  = f'Re{{B\u2071{_comp_label}}} (nT)'
            right_label = f'Im{{B\u2071{_comp_label}}} (nT)'
            if not use_nT_amplitudes:
                # Override: dimensionless A_e * exp(i*phi)
                amp_data = _get_induction_slice(induction, 'Amp', iPeak, VALID)
                phase_data = _get_induction_slice(induction, 'Phase', iPeak, VALID)
                if amp_data is not None and phase_data is not None:
                    phase_rad = np.deg2rad(phase_data)
                    re_data = amp_data * np.cos(phase_rad)
                    im_data = amp_data * np.sin(phase_rad)
                else:
                    re_data, im_data = None, None
                left_label  = 'Re{A\u2091}'
                right_label = 'Im{A\u2091}'
            left_data, right_data = re_data, im_data
            left_cscale = colorscale_div
            right_cscale = colorscale_div
        else:
            # Amplitude + Phase
            if use_nT_amplitudes:
                bi1_arr = getattr(induction, f'Bi1{_comp_suffix}_nT', None)
                if bi1_arr is not None:
                    abs_bi = np.abs(bi1_arr)
                    if len(abs_bi.shape) == 3:
                        amp_data = abs_bi[iPeak, :, :].copy()
                    elif len(abs_bi.shape) == 2:
                        amp_data = abs_bi.copy()
                    else:
                        amp_data = None
                    if amp_data is not None and VALID is not None:
                        amp_data[~VALID] = np.nan
                else:
                    amp_data = None
                left_label = f'|B\u2071{_comp_label}| (nT)'
            else:
                amp_data = _get_induction_slice(induction, 'Amp', iPeak, VALID)
                left_label = 'Amplitude A\u2091'
            phase_data = _get_induction_slice(induction, 'Phase', iPeak, VALID)
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

    # Apply per-axis log scaling to every subplot
    if x_log:
        fig.update_xaxes(type='log')
    if y_log:
        fig.update_yaxes(type='log')

    return fig


def _get_induction_slice(induction, attr_name, iPeak, VALID):
    """Extract a 2D slice from a 3D induction array for a given peak index.

    Stored attributes on InductionData (see ResultsStructs.py and
    ResultsIO.py:ExtractInductionData):
      - real-valued (np.iscomplexobj returns False): Amp, Phase,
        Bix_nT, Biy_nT, Biz_nT, rBi1Tot_nT, iBi1Tot_nT,
        rBi1{x,y,z}_nT, iBi1{x,y,z}_nT
      - complex-valued: Bi1x_nT, Bi1y_nT, Bi1z_nT, Bi1Tot_nT

    For the complex-stored attributes the only meaningful scalar reduction
    here is the modulus (amplitude).  Sign-bearing real/imaginary panels
    must be requested by their explicit r*/i* attribute names so they keep
    their sign — callers must NOT pass Bi1Tot_nT when they want Re or Im.
    """
    data = getattr(induction, attr_name, None)
    if data is None:
        return None

    if np.iscomplexobj(data):
        # Only Bi1{x,y,z,Tot}_nT are complex. Return modulus for those.
        # Real-stored attributes (rBi1Tot_nT, iBi1Tot_nT, Phase, Amp, etc.)
        # fall through unchanged — do NOT apply np.abs to them as that would
        # discard the sign of real/imaginary-component panels.
        data = np.abs(data)

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
