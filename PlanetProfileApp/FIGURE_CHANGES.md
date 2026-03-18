# PlanetProfileApp Figure Handling Updates

## Overview

PlanetProfileApp has been updated to provide better figure management with the following improvements:

1. **No automatic figure file generation** during runs
2. **Interactive Plotly figures** generated on-demand from data
3. **Save figures only when requested** by the user
4. **Unique timestamped data files** to prevent overwriting

## What Changed

### During Model Runs

When you run a PlanetProfile simulation through the app:

- ✅ **Data files ARE saved** with timestamps (e.g., `Europa_20260310_143022.txt`)
- ❌ **Figure PDFs are NOT generated** automatically (saves time and disk space)
- ✅ **Calculations are performed** and results are stored

Configuration changes in `RunPlanetProfile.py`:
```python
Params.SKIP_PLOTS = True          # Skip PDF generation
Params.TIME_AND_DATE_LABEL = True # Unique filenames with timestamps
Params.NO_SAVEFILE = False        # Still save data files
```

### Viewing Outputs

On the **PlanetProfile Outputs** page, you can now choose:

#### Interactive Mode (Recommended)
- Generates figures **from saved data files** on-demand
- Fully interactive (zoom, pan, hover for values)
- Multiple figure types:
  - **Hydrosphere Properties**: T, P, ρ, conductivity, sound speeds vs depth
  - **Seismic Properties**: Moduli, sound speeds, quality factor
  - **Custom Plots**: Choose any x/y variables with log scales
- **Save figures** individually in multiple formats (HTML, PNG, PDF)

#### Legacy Mode
- Displays pre-generated PDF figures (if any exist)
- Useful for viewing old runs that were created with PDF generation

## Benefits

1. **Faster runs**: Skipping PDF generation saves significant time
2. **Less disk space**: Only data files are saved, not large figure PDFs
3. **More flexible**: Generate any plot combination from the data
4. **Interactive exploration**: Zoom, pan, and hover to explore data
5. **No overwrites**: Timestamped filenames prevent data loss

## Installation

### Required
```bash
pip install streamlit pandas
```

### Optional (for interactive figures)
```bash
pip install plotly kaleido
```

Note: `kaleido` is needed for saving Plotly figures as PNG/PDF. HTML export works without it.

## Usage

### Running a Simulation

1. Configure your planet on the settings pages
2. Click "Run PlanetProfile" on the Run page
3. Data files are saved with timestamps (no PDF figures generated)

### Viewing Results

1. Go to **PlanetProfile Outputs** page
2. Select **"🎯 Interactive (Plotly)"** mode (if Plotly is installed)
3. Choose figure type from tabs
4. Select which run to visualize (shows 5 most recent)
5. Explore interactively with zoom/pan
6. **Click "💾 Save Figure"** to export in your desired format

### Custom Plots

Use the **"Custom Plot"** tab to:
- Select any X-axis variable
- Select multiple Y-axis variables
- Enable log scales for either axis
- Create publication-ready figures from your data

## Technical Details

### Figure Generator Module

Location: `PlanetProfileApp/Utilities/figure_generator.py`

Functions:
- `load_profile_data(file_path)`: Load PlanetProfile .txt output
- `load_ocean_data(file_path)`: Load liquidOceanProps.txt
- `create_hydrosphere_plot(df, planet_name)`: Generate 6-panel hydrosphere figure
- `create_seismic_plot(df, planet_name)`: Generate 4-panel seismic figure
- `create_custom_plot(df, x_col, y_cols, ...)`: Generate user-defined plot
- `save_plotly_figure(fig, filename, format)`: Export figure to file

### Data File Formats

Profile files contain:
- Lines 1-81: Header metadata
- Line 82: Column headers (if present)
- Lines 83+: Numerical data (whitespace-separated)

Columns (23 total):
- P (MPa), T (K), r (m), phase ID, rho (kg/m3)
- Cp (J/kg/K), alpha (1/K), g (m/s2), phi, sigma (S/m)
- k (W/m/K), VP (km/s), VS (km/s), QS
- KS (GPa), GS (GPa), Ppore (MPa)
- rhoMatrix, rhoPore, MLayer, VLayer, Htidal, eta (Pa s)

## Troubleshooting

### "Plotly not available" warning
```bash
pip install plotly
```

### "Failed to save figure" (PNG/PDF)
```bash
pip install kaleido
```

### No data files found
Run a simulation first from the "Run PlanetProfile" page.

### Old PDF figures still showing
Switch to "Interactive (Plotly)" mode at the top of the Outputs page.

## Migration from Old Runs

Old runs that generated PDF figures:
- PDF figures remain accessible in "Legacy Mode"
- Data files can be reloaded and used to generate interactive figures
- No data is lost

## Future Enhancements

Potential additions:
- Export data tables to CSV
- Side-by-side comparison of multiple runs
- Animation of parameter sweeps
- 3D visualizations for asymmetric models
- Batch figure export
