# Planet Profile App

GUI for PlanetProfile, found here: git@github.com:vancesteven/PlanetProfile.git

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://blank-app-template.streamlit.app/)

## Features

- 🎯 **Interactive Figures**: Generate Plotly figures on-demand from data
- 💾 **Smart File Management**: Timestamped data files prevent overwrites
- 🚀 **Fast Runs**: Skip figure generation during simulations
- 📊 **Custom Plots**: Create any plot from saved data with zoom/pan
- 💿 **Flexible Export**: Save figures in HTML, PNG, or PDF format

## Recent Updates

**NEW**: The app now generates interactive figures from data files instead of static PDFs:
- Model runs no longer generate figure PDFs (saves time and space)
- Interactive Plotly figures are created on-demand in the Outputs page
- Save individual figures only when you want them
- Data files use unique timestamps to prevent overwrites

See [FIGURE_CHANGES.md](FIGURE_CHANGES.md) for details.

### How to run it on your own machine

1. Install the requirements

   ```bash
   # Required
   $ pip install streamlit pdf2image Pillow pandas
   $ conda install poppler  # if you’re on macOS
   # On Windows, you’ll need to download Poppler from: https://github.com/oschwartz10612/poppler-windows,
   # Then add its /bin folder to your PATH.

   # Optional (for interactive figures - recommended)
   $ pip install plotly
   $ pip install kaleido  # for saving PNG/PDF from Plotly
   ```

2. Run the app

   ```bash
   $ streamlit run PlanetProfileApp.py
   ```

## Usage

1. **Configure** your planet on the settings pages (Main Settings, Bulk Properties, Ocean, Core, Layers)
2. **Run** the simulation from the "Run PlanetProfile" page
3. **View** interactive figures and data on the "PlanetProfile Outputs" page
4. **Save** individual figures in your preferred format
