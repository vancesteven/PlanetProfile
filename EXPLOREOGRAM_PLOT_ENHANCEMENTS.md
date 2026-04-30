# Exploreogram Plot Type Enhancements

## Date: 2026-03-10

## Summary

Added interactive plot type selection to the Exploreogram feature, allowing users to choose between interactive Plotly plots (default) and publication-quality matplotlib plots (matching CLI output).

## Motivation

Users requested the ability to:
1. Have interactive plots for data exploration (already partially implemented)
2. Generate matplotlib plots that match the default CLI exploreogram output
3. Toggle between plot types easily
4. Get publication-quality figures directly from the GUI

## Implementation

### 1. Enhanced Plotly Plotting Module

**File**: `PlanetProfileApp/Utilities/exploreogram_plotly.py` (NEW)

**Features**:
- `create_exploreogram_plotly()`: Creates Plotly plots styled to match matplotlib defaults
- Matches PlanetProfile matplotlib colormaps (viridis, etc.)
- Handles zero-centered colormaps for bipolar variables
- Supports contour line overlays
- Optional smoothing support
- `create_multi_exploreogram_plotly()`: Multi-subplot support for multiple z-variables

**Styling**:
- Uses matplotlib colormaps converted to Plotly format
- Matches axis labels from FigLbl system
- Respects exploration-specific formatting
- Handles invalid models (NaN masking)

### 2. Exploreogram Page Updates

**File**: `PlanetProfileApp/pages/Exploreogram.py`

**Changes**:

#### Added Interactive Plot Toggle
```python
use_interactive_plots = st.checkbox(
    "🎯 Interactive Plots",
    value=True,
    help="Use interactive Plotly plots (default) or matplotlib static plots"
)
```

#### Dual Plotting System

**Interactive Mode (Default)**:
- Uses enhanced `create_exploreogram_plotly()` function
- Loads FigLbl for proper styling
- Fallback to basic Plotly if enhanced plotting fails
- Fast rendering, immediate display

**Static Mode**:
- Calls `GenerateExplorationPlots()` from PlanetProfile
- Generates publication-quality matplotlib PDFs
- Identical to CLI exploreogram output
- Displays rendered PDF in GUI

#### Enhanced Save Options

**Interactive mode**:
- Save Interactive (HTML): Fully interactive browser-viewable plot
- Save Image (PNG): Static image (requires kaleido)
- Save Data (PKL): Full results for reanalysis

**Static mode**:
- PDF automatically saved by GenerateExplorationPlots
- Save Data (PKL): Full results for reanalysis

### 3. Warning Suppression

**Issue**: PyALMA parallel processing warning was confusing users

**Solution**: Added warning filters:
```python
warnings.filterwarnings('ignore', message='.*parallel setting should be boolean.*')
warnings.filterwarnings('ignore', message='.*Reverting to serial operation.*')
```

**Explanation**: This warning is expected and doesn't indicate a problem. Models run in parallel across the grid (fast), while individual model calculations run serially (prevents nested parallelization).

### 4. Documentation Updates

**File**: `EXPLOREOGRAM_FEATURE.md`

Added sections on:
- Plot Type Selection (Interactive vs Static)
- Interactive Plots features and use cases
- Static Plots features and use cases
- Save options for each mode
- Tips for choosing plot type

## User Workflow

### Interactive Plots (Default)

1. Run exploreogram as normal
2. Results display as interactive Plotly heatmap
3. **Use**: Hover to see values, zoom to specific regions, pan around
4. **Save**: HTML for sharing, PNG for presentations, PKL for data

### Static Plots

1. Uncheck "🎯 Interactive Plots" at top of page
2. Run exploreogram as normal
3. Results display as rendered matplotlib PDF
4. **Use**: Publication figures, matching CLI output exactly
5. **Save**: PDF for papers (automatically saved), PKL for data

## Technical Details

### Matplotlib → Plotly Colormap Conversion

```python
def get_matplotlib_colormap_colors(cmap_name='viridis', n_colors=256):
    """Convert matplotlib colormap to Plotly format"""
    cmap = cm.get_cmap(cmap_name)
    colors = [cmap(i / (n_colors - 1)) for i in range(n_colors)]
    hex_colors = ['#{:02x}{:02x}{:02x}'.format(
        int(r*255), int(g*255), int(b*255)
    ) for r, g, b, a in colors]
    return hex_colors
```

### FigLbl Integration

```python
from PlanetProfile.GetConfig import FigLbl
FigLbl.SetExploration(bodyname, xName, yName, zName)
# Now FigLbl has proper labels and formatting
```

### Matplotlib Plot Generation

```python
from PlanetProfile.Plotting.ExplorationPlots import GenerateExplorationPlots
GenerateExplorationPlots(ExplorationList, FigureFilesList, Params)
# Generates PDF matching CLI output
```

## Advantages

### Interactive Plots
- ✅ Fast rendering and display
- ✅ Zoom/pan/hover for exploration
- ✅ Easy sharing (HTML)
- ✅ Good for presentations
- ✅ No additional dependencies

### Static Plots
- ✅ Publication-quality output
- ✅ Identical to CLI output
- ✅ Vector graphics (PDF)
- ✅ Advanced features (contours, smoothing)
- ✅ LaTeX formatting
- ⚠️ Slower rendering (matplotlib PDF generation)
- ⚠️ Requires pdf2image for display

## Dependencies

### Required
- `plotly`: Interactive plotting (already required)
- `matplotlib`: Static plotting (already required)
- `pdf2image`: PDF rendering for display (already required)

### Optional
- `kaleido`: PNG export from Plotly (`pip install kaleido`)

## Future Enhancements

Potential improvements:
1. **Multi-variable interactive plots**: Plotly version of multi-subplot exploreograms
2. **Smoothing controls**: GUI slider for smoothing factor
3. **Colormap selection**: Dropdown to choose colormap
4. **Contour line toggle**: Checkbox to add/remove contour lines in interactive mode
5. **3D surface plots**: Alternative visualization for parameter space
6. **Animation**: Animated sweep through parameter space

## Files Modified

```
PlanetProfileApp/
├── Utilities/
│   └── exploreogram_plotly.py (NEW - 350 lines)
├── pages/
│   └── Exploreogram.py (MODIFIED - enhanced plotting section)
EXPLOREOGRAM_FEATURE.md (MODIFIED - documentation updates)
EXPLOREOGRAM_PLOT_ENHANCEMENTS.md (NEW - this file)
```

## Testing

### Manual Test Procedure

1. **Test Interactive Mode (Default)**:
   ```bash
   streamlit run PlanetProfileApp/PlanetProfileApp.py
   ```
   - Select Europa
   - Navigate to Exploreogram page
   - Verify "🎯 Interactive Plots" is checked
   - Configure: X=wOcean_ppt (0-50), Y=Tb_K (265-275), Grid=5×5
   - Click "Run Exploreogram"
   - **Verify**:
     - ✅ Plotly heatmap displays
     - ✅ Hover shows values
     - ✅ Zoom/pan work
     - ✅ Colormap matches viridis
     - ✅ Can save as HTML
     - ✅ Can save as PNG (if kaleido installed)

2. **Test Static Mode**:
   - Uncheck "🎯 Interactive Plots"
   - Click "Run Exploreogram" (or view previous results)
   - **Verify**:
     - ✅ Matplotlib PDF displays
     - ✅ Matches CLI exploreogram styling
     - ✅ PDF saved to output/exploreograms/
     - ✅ Contours visible (if applicable)
     - ✅ Publication-quality formatting

3. **Test Warning Suppression**:
   - Run exploreogram (any mode)
   - **Verify**:
     - ✅ No "parallel setting should be boolean" warning
     - ✅ Models run in parallel (check timing)

## Known Limitations

1. **Matplotlib PDF display**: Requires pdf2image library
2. **PNG export**: Requires kaleido for Plotly image export
3. **Static mode speed**: Slower than interactive due to matplotlib rendering
4. **Multi-variable plots**: Currently only single z-variable in Plotly mode

## Backward Compatibility

- ✅ Existing exploreogram functionality unchanged
- ✅ Default is interactive (better UX)
- ✅ Can switch to matplotlib anytime
- ✅ PKL files compatible with CLI tools
- ✅ Configuration records still saved

## Summary

Successfully implemented dual plotting system for exploreograms:
- **Interactive (default)**: Fast, explorable Plotly plots with matplotlib styling
- **Static (optional)**: Publication-quality matplotlib plots matching CLI output
- **User choice**: Toggle between modes with single checkbox
- **Backward compatible**: All existing functionality preserved

Users now have the best of both worlds: interactive exploration AND publication-quality figures.

---

**Status**: ✅ **IMPLEMENTED**

**Version**: PlanetProfile 3.1.0+exploreogram-enhanced-plots

**Date**: 2026-03-10
