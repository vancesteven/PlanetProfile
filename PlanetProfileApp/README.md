# Planet Profile App

GUI for PlanetProfile, found here: git@github.com:vancesteven/PlanetProfile.git

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://blank-app-template.streamlit.app/)

## Features

- **Interactive outputs:** Generate Plotly figures on demand from saved run
  data, with zoom and pan instead of pre-generating every PDF during a model
  run.
- **Posterior globe and interior explorer:** The inference results page links a
  posterior sample picker to an interactive, textured 3D cutaway globe and the
  radial profiles, wedge, heating, geotherm, mineralogy, and data-table tabs.
  The globe splits packaged `D_icePhase_results` into proportional Ice III,
  Ice V, and Ice VI shells rather than collapsing the high-pressure ice stack.
- **No-WebGL fallback:** The exact control label
  **“Static render (no WebGL)”** switches the globe to a server-rendered
  Matplotlib view with the same cutaway, layer colors, degree-2 deformation,
  and principal axes.
- **k2 measurement view:** The posterior scatter overlays the current run's
  1-sigma and 2-sigma observational ellipses. **“Zoom to measurement ellipse (±4σ)”**
  is enabled by default so the measurement remains legible when the posterior
  is broad; disabling it restores the full autoscaled sample cloud.
- **Traceable figure exports:** Crisp SVG display plus PDF/PNG downloads embed
  the model slot, artifact filename, conditioning values, and app Git SHA in
  backend-supported metadata. Download names follow
  `<figure>_<slot-short>_<UTC-date>.<ext>`.
- **Visible validation scope:** Artifact slots render their `scope_note` and
  **“Gate status:”** beside the conditioning controls. Split-ratified results
  carry the slot's `sector_warning` to the top of the results page, so panels
  that include quarantined quantities cannot hide the warning.
- **Smart file management:** Timestamped data files prevent overwrites, and
  saved inference results retain the slot/provenance context used by exports.

## Recent Updates

The app generates interactive figures from data files instead of static PDFs:
- Model runs no longer generate figure PDFs (saves time and space)
- Interactive Plotly figures are created on-demand in the Outputs page
- Save individual figures only when you want them
- Data files use unique timestamps to prevent overwrites

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
