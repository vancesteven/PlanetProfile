import streamlit as st
import os
import sys
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pickle
import time

_DEBUG_EXPLOREOGRAM = os.environ.get('PP_EXPLOREOGRAM_DEBUG', '0') == '1'

# --- Setting up Page ---
st.set_page_config(
    page_title="Exploreogram",
    page_icon="🔬",
    layout="wide"
)

# --- File Path Management ---
this_file = os.path.abspath(__file__)
pages_directory = os.path.dirname(this_file)
app_directory = os.path.dirname(pages_directory)
parent_directory = os.path.dirname(app_directory)

# Add both app directory (for configPP) and parent directory (for PlanetProfile) to path
if app_directory not in sys.path:
    sys.path.insert(0, app_directory)
if parent_directory not in sys.path:
    sys.path.insert(0, parent_directory)

# Import the authoritative list of magnetic induction z-variable names.
# This is defined in exploreogram_plotly to avoid duplication.
from Utilities.exploreogram_plotly import INDUCTION_Z_VARS


def _strongest_induction_axis(induction):
    """Return 'Bx', 'By', 'Bz', or 'Total' — whichever has highest finite mean amplitude at peak 0."""
    if induction is None:
        return 'Total'
    candidates = [('Bx', 'Bi1x_nT'), ('By', 'Bi1y_nT'), ('Bz', 'Bi1z_nT')]
    means = []
    for label, attr in candidates:
        arr = getattr(induction, attr, None)
        if arr is None:
            continue
        slice0 = np.abs(arr[0]) if arr.ndim == 3 else np.abs(arr)
        m = np.nanmean(slice0)
        if np.isfinite(m):
            means.append((label, m))
    if not means:
        return 'Total'
    return max(means, key=lambda t: t[1])[0]


st.title("🔬 Exploreogram: Parameter Space Exploration")

# Interactive plot toggle
col_title, col_interactive = st.columns([3, 1])
with col_interactive:
    use_interactive_plots = st.checkbox(
        "🎯 Interactive Plots",
        value=False,
        help="Use interactive Plotly plots or matplotlib static plots (default, matches CLI output)"
    )

st.markdown("---")

st.markdown("""
Exploreograms run PlanetProfile across a 2D grid of parameters to explore how interior
properties vary across parameter space. This is useful for sensitivity studies,
parameter sweeps, and understanding model behavior.

**Note**: Exploreograms can be computationally expensive. A 10×10 grid runs 100 models!
""")

# Check if planet is selected
if "Planet" not in st.session_state or st.session_state["Planet"] is None:
    st.warning("⚠️ Please select a planet on the Main Settings page first!")
    st.stop()

# Available exploration parameters (from ExploreParamsStruct)
HYDRO_PARAMS = {
    'wOcean_ppt': {'label': 'Ocean Salinity (ppt)', 'default_range': [0, 100], 'desc': 'Ocean salt concentration'},
    'Tb_K': {'label': 'Ocean Bottom Temperature (K)', 'default_range': [250, 280], 'desc': 'Temperature at ice-ocean interface'},
    'icePhi_frac': {'label': 'Ice Porosity (fraction)', 'default_range': [0, 0.3], 'desc': 'Ice shell porosity'},
    'icePclosure_MPa': {'label': 'Ice Closure Pressure (MPa)', 'default_range': [50, 200], 'desc': 'Pressure where ice pores close'},
    'zb_approximate_km': {'label': 'Ice Shell Thickness (km)', 'default_range': [5, 30], 'desc': 'Approximate ice shell thickness'},
}

INNER_PARAMS = {
    'xFeS': {'label': 'Core FeS Fraction', 'default_range': [0, 0.2], 'desc': 'Iron sulfide fraction in core'},
    'rhoSilInput_kgm3': {'label': 'Silicate Density (kg/m³)', 'default_range': [3000, 3600], 'desc': 'Mantle silicate density'},
    'Rcore_km': {'label': 'Core Radius (km)', 'default_range': [500, 800], 'desc': 'Radius of iron core'},
    'silPhi_frac': {'label': 'Silicate Porosity (fraction)', 'default_range': [0, 0.4], 'desc': 'Mantle porosity'},
    'silPclosure_MPa': {'label': 'Silicate Closure Pressure (MPa)', 'default_range': [200, 600], 'desc': 'Pressure where mantle pores close'},
    'Htidal_Wm3': {'label': 'Tidal Heating (W/m³)', 'default_range': [1e-11, 1e-9], 'desc': 'Volumetric tidal heating rate'},
    'Qrad_Wkg': {'label': 'Radiogenic Heating (W/kg)', 'default_range': [1e-12, 1e-10], 'desc': 'Mass-specific radiogenic heating'},
    'qSurf_Wm2': {'label': 'Surface Heat Flux (W/m²)', 'default_range': [0.01, 0.1], 'desc': 'Surface heat flux'},
}

IONOS_PARAMS = {
    'ionosTop_km': {'label': 'Ionosphere Top (km)', 'default_range': [50, 200], 'desc': 'Top altitude of ionosphere'},
    'sigmaIonos_Sm': {'label': 'Ionosphere Conductivity (S/m)', 'default_range': [1e-4, 1e-2], 'desc': 'Ionosphere electrical conductivity'},
}

DERIVED_PARAMS = {
    'D_km': {'label': 'Ocean Thickness (km)', 'default_range': [10, 100], 'desc': 'Ocean layer thickness (derived — requires a driver parameter)',
             'default_driver': 'Tb_K'},
    'rhoSilMean_kgm3': {'label': 'Mean Silicate Density (kg/m³)', 'default_range': [2500, 3600], 'desc': 'Model-computed mean rock density (derived — requires a driver parameter)',
                         'default_driver': 'rhoSilInput_kgm3'},
    'sigmaMean_Sm': {'label': 'Mean Ocean Conductivity (S/m)', 'default_range': [0.1, 10], 'desc': 'Model-computed mean ocean electrical conductivity (derived — requires a driver parameter)',
                     'default_driver': 'wOcean_ppt'},
}

# Input parameters that can actually be varied during exploration
INPUT_PARAMS = {**HYDRO_PARAMS, **INNER_PARAMS, **IONOS_PARAMS}

ALL_PARAMS = {**INPUT_PARAMS, **DERIVED_PARAMS}

# Clear excitation selection when the planet changes between sessions.
# This prevents stale frequency selections from carrying over to a new body.
_current_planet_name = st.session_state["Planet"].name if "Planet" in st.session_state and st.session_state["Planet"] is not None else None
if '_last_planet_name' not in st.session_state:
    st.session_state['_last_planet_name'] = _current_planet_name
elif st.session_state['_last_planet_name'] != _current_planet_name:
    st.session_state['selected_excitations'] = []
    st.session_state['_last_planet_name'] = _current_planet_name

# Initialize session state for exploreogram settings
if 'explore_xName' not in st.session_state:
    st.session_state.explore_xName = 'wOcean_ppt'
if 'explore_yName' not in st.session_state:
    st.session_state.explore_yName = 'Tb_K'
if 'explore_nx' not in st.session_state:
    st.session_state.explore_nx = 10
if 'explore_ny' not in st.session_state:
    st.session_state.explore_ny = 10
if 'explore_zName' not in st.session_state:
    st.session_state.explore_zName = 'zb_km'
if 'explore_running' not in st.session_state:
    st.session_state.explore_running = False
if 'explore_results' not in st.session_state:
    st.session_state.explore_results = None
if 'explore_start_time' not in st.session_state:
    st.session_state.explore_start_time = None
if 'explore_error' not in st.session_state:
    st.session_state.explore_error = None
if 'explore_error_traceback' not in st.session_state:
    st.session_state.explore_error_traceback = None
if 'explore_cache_key' not in st.session_state:
    st.session_state.explore_cache_key = None
if 'explore_force_rerun' not in st.session_state:
    st.session_state.explore_force_rerun = False

# Main layout: settings and results
col1, col2 = st.columns([1, 2])

with col1:
    st.subheader("⚙️ Exploration Settings")

    # X-axis parameter
    st.markdown("#### X-Axis Parameter")
    x_param = st.selectbox(
        "Select X parameter",
        options=list(ALL_PARAMS.keys()),
        index=list(ALL_PARAMS.keys()).index(st.session_state.explore_xName),
        format_func=lambda x: ALL_PARAMS[x]['label'],
        key='x_param_select'
    )
    st.session_state.explore_xName = x_param
    st.caption(ALL_PARAMS[x_param]['desc'])

    # If derived parameter, show driver parameter selector
    if x_param in DERIVED_PARAMS:
        st.info(f"{ALL_PARAMS[x_param]['label']} is a model output. Select which input parameter to vary — the X-axis will show the resulting {ALL_PARAMS[x_param]['label']}.")
        x_driver_default = DERIVED_PARAMS[x_param].get('default_driver', list(INPUT_PARAMS.keys())[0])
        x_driver = st.selectbox(
            "X driver parameter (varied during exploration)",
            options=list(INPUT_PARAMS.keys()),
            index=list(INPUT_PARAMS.keys()).index(x_driver_default),
            format_func=lambda k: INPUT_PARAMS[k]['label'],
            key='x_driver_select'
        )
        st.session_state['explore_x_driver'] = x_driver
        x_range_default = INPUT_PARAMS[x_driver]['default_range']
        x_min = st.number_input(
            f"{INPUT_PARAMS[x_driver]['label']} - Min",
            value=float(x_range_default[0]),
            key='x_min'
        )
        x_max = st.number_input(
            f"{INPUT_PARAMS[x_driver]['label']} - Max",
            value=float(x_range_default[1]),
            key='x_max'
        )
    else:
        st.session_state['explore_x_driver'] = None
        x_range_default = ALL_PARAMS[x_param]['default_range']
        x_min = st.number_input(
            f"{ALL_PARAMS[x_param]['label']} - Min",
            value=float(x_range_default[0]),
            key='x_min'
        )
        x_max = st.number_input(
            f"{ALL_PARAMS[x_param]['label']} - Max",
            value=float(x_range_default[1]),
            key='x_max'
        )

    # Y-axis parameter
    st.markdown("#### Y-Axis Parameter")
    y_param = st.selectbox(
        "Select Y parameter",
        options=list(ALL_PARAMS.keys()),
        index=list(ALL_PARAMS.keys()).index(st.session_state.explore_yName),
        format_func=lambda x: ALL_PARAMS[x]['label'],
        key='y_param_select'
    )
    st.session_state.explore_yName = y_param
    st.caption(ALL_PARAMS[y_param]['desc'])

    # If derived parameter, show driver parameter selector
    if y_param in DERIVED_PARAMS:
        st.info(f"{ALL_PARAMS[y_param]['label']} is a model output. Select which input parameter to vary — the Y-axis will show the resulting {ALL_PARAMS[y_param]['label']}.")
        y_driver_default = DERIVED_PARAMS[y_param].get('default_driver', list(INPUT_PARAMS.keys())[0])
        y_driver = st.selectbox(
            "Y driver parameter (varied during exploration)",
            options=list(INPUT_PARAMS.keys()),
            index=list(INPUT_PARAMS.keys()).index(y_driver_default),
            format_func=lambda k: INPUT_PARAMS[k]['label'],
            key='y_driver_select'
        )
        st.session_state['explore_y_driver'] = y_driver
        y_range_default = INPUT_PARAMS[y_driver]['default_range']
        y_min = st.number_input(
            f"{INPUT_PARAMS[y_driver]['label']} - Min",
            value=float(y_range_default[0]),
            key='y_min'
        )
        y_max = st.number_input(
            f"{INPUT_PARAMS[y_driver]['label']} - Max",
            value=float(y_range_default[1]),
            key='y_max'
        )
    else:
        st.session_state['explore_y_driver'] = None
        y_range_default = ALL_PARAMS[y_param]['default_range']
        y_min = st.number_input(
            f"{ALL_PARAMS[y_param]['label']} - Min",
            value=float(y_range_default[0]),
            key='y_min'
        )
        y_max = st.number_input(
            f"{ALL_PARAMS[y_param]['label']} - Max",
            value=float(y_range_default[1]),
            key='y_max'
        )

    # Grid resolution
    st.markdown("#### Grid Resolution")
    col_nx, col_ny = st.columns(2)
    with col_nx:
        nx = st.number_input("X points", min_value=2, max_value=50, value=st.session_state.explore_nx, key='nx_input')
        st.session_state.explore_nx = nx
    with col_ny:
        ny = st.number_input("Y points", min_value=2, max_value=50, value=st.session_state.explore_ny, key='ny_input')
        st.session_state.explore_ny = ny

    total_runs = nx * ny
    st.info(f"Total runs: **{total_runs}** models")

    # Z-axis (color) parameter selection
    st.markdown("#### Color Variable (Z-axis)")

    Z_VARIABLES = {
        'zb_km': 'Ice Shell Thickness (km)',
        'D_km': 'Ocean Layer Thickness (km)',
        'Tmean_K': 'Mean Ocean Temperature (K)',
        'Pseafloor_MPa': 'Seafloor Pressure (MPa)',
        'sigmaMean_Sm': 'Mean Ocean Conductivity (S/m)',
        'Amp_nT': 'Induced Magnetic Field (nT)',
        'kLoveAmp': 'k₂ Love Number Amplitude',
        'Rcore_km': 'Core Radius (km)',
    }

    z_param = st.selectbox(
        "Select color variable",
        options=list(Z_VARIABLES.keys()),
        index=0,
        format_func=lambda x: Z_VARIABLES[x],
        key='z_param_select'
    )
    st.session_state.explore_zName = z_param

    # Update cached results zName instantly (no recomputation needed)
    if st.session_state.explore_results is not None:
        st.session_state.explore_results['zName'] = z_param
        if hasattr(st.session_state.explore_results.get('Planet', None), 'Exploration'):
            st.session_state.explore_results['Planet'].Exploration.zName = z_param

    # Show induction frequency selection if magnetic variable is chosen
    induction_z_vars = INDUCTION_Z_VARS
    if z_param in induction_z_vars:
        st.markdown("#### Magnetic Induction Configuration")
        st.info("Magnetic induction calculations enabled for selected output variable")

        # Get available excitations for this body
        try:
            from PlanetProfile.MagneticInduction.Moments import Excitations
            bodyname = st.session_state.get("Planet").name if "Planet" in st.session_state else "Europa"
            available_exc = Excitations(bodyname)

            if available_exc:
                st.write(f"**Available excitation frequencies for {bodyname}:**")

                # Define key excitations (most commonly used)
                key_excitations = ['orbital', 'synodic', 'true anomaly', 'synodic 2nd']

                # Initialize session state for excitation selection if not exists
                if 'selected_excitations' not in st.session_state:
                    # Default to key excitations that are available
                    st.session_state.selected_excitations = [
                        exc for exc in key_excitations if exc in available_exc
                    ]

                # Show checkboxes for each available excitation
                st.write("Select excitation frequencies to calculate:")

                selected = []
                for exc_name, period_hr in sorted(available_exc.items(), key=lambda x: x[1]):
                    # Mark key excitations
                    is_key = exc_name in key_excitations
                    label_suffix = " ⭐ (recommended)" if is_key else ""
                    default_checked = exc_name in st.session_state.selected_excitations

                    if st.checkbox(
                        f"{exc_name}: {period_hr:.2f} hr{label_suffix}",
                        value=default_checked,
                        key=f"exc_{exc_name}"
                    ):
                        selected.append(exc_name)

                st.session_state.selected_excitations = selected

                if len(selected) == 0:
                    st.warning("⚠️ No excitations selected! Select at least one frequency to calculate induction.")
                else:
                    st.success(f"✓ {len(selected)} excitation(s) selected")

            else:
                st.warning(f"No excitation data available for {bodyname}. Using default configuration.")
                if 'selected_excitations' not in st.session_state:
                    st.session_state.selected_excitations = ['orbital', 'synodic', 'true anomaly']

        except Exception as e:
            st.error(f"Error loading excitation data: {e}")
            if 'selected_excitations' not in st.session_state:
                st.session_state.selected_excitations = ['orbital', 'synodic', 'true anomaly']

        # Inductogram display mode
        st.markdown("#### Inductogram Display")
        induct_display_mode = st.radio(
            "Component display",
            options=['real_imaginary', 'amplitude_phase'],
            format_func=lambda x: 'Real + Imaginary' if x == 'real_imaginary' else 'Amplitude + Phase',
            index=0,
            key='induct_display_mode',
            help="Real+Imaginary is preferred for modern publications. Amplitude+Phase is the traditional format."
        )

    # Use-nT-amplitudes option — controls whether amplitude/Re/Im are shown
    # in nT (default) or dimensionless Ae (0–1).  Applies to amplitude, real,
    # and imaginary components for both Plotly and matplotlib paths.
    if 'exploreogram_use_nT_amplitudes' not in st.session_state:
        st.session_state['exploreogram_use_nT_amplitudes'] = True
    if z_param in INDUCTION_Z_VARS:
        use_nT_amplitudes = st.checkbox(
            "Use nT amplitudes",
            value=True,
            key='exploreogram_use_nT_amplitudes',
            help="Show amplitude / Re / Im induced-field components in nT (default). "
                 "Uncheck to use dimensionless Ae (0–1)."
        )
    else:
        use_nT_amplitudes = st.session_state['exploreogram_use_nT_amplitudes']

    # Component selector — only meaningful for Amp_nT
    if z_param == 'Amp_nT':
        _induct_component_options = ['Total', 'Bx', 'By', 'Bz']
        _induct_component_default = st.session_state.get('exploreogram_induct_component_default', 'Total')
        _default_idx = _induct_component_options.index(_induct_component_default) if _induct_component_default in _induct_component_options else 0
        induct_component = st.radio(
            'Induction component',
            _induct_component_options,
            index=_default_idx,
            key='exploreogram_induct_component',
            help='Which induced field axis to display. "Total" = |Bi1Tot|. Default auto-selects the strongest axis after a run.',
            horizontal=True,
        )
    else:
        induct_component = st.session_state.get('exploreogram_induct_component', 'Total')

    # Salinity↔conductivity secondary axis — auto-activates when one of the
    # plot axes is wOcean_ppt or sigmaMean_Sm.  Uses np.interp on PP-computed
    # (w, σ) pairs; silently disables for non-monotone compositions (e.g.
    # MgSO4 near eutectic).  Only applies to single-panel Plotly plots;
    # ignored for multi-panel inductograms.
    show_salinity_axis = st.checkbox(
        "Show salinity ↔ conductivity axis",
        value=True,
        key='exploreogram_show_salinity_axis',
        help="When the x- or y-axis is salinity (wOcean_ppt) or mean ocean "
             "conductivity (sigmaMean_Sm), show a synchronized secondary axis "
             "with the linked variable. Uses PP-computed values (not a linear "
             "approximation); disables automatically for non-monotone compositions."
    )

    # Contour option — applies to ALL exploreograms (not just induction)
    induct_use_contours = st.checkbox(
        "Use contour lines",
        value=True,
        key='induct_use_contours',
        help="Draw contour lines on the colormap. Uncheck for colormap only."
    )

    # Per-axis log-scale toggles — applies to both interactive Plotly and matplotlib paths
    x_log = st.checkbox(
        "Log x-axis",
        value=False,
        key="exploreogram_x_log",
        help="Use logarithmic scale on the x-axis. Requires positive values."
    )
    if x_log and x_min <= 0:
        st.caption("Log x-axis ignored: x min must be > 0")
    y_log = st.checkbox(
        "Log y-axis",
        value=False,
        key="exploreogram_y_log",
        help="Use logarithmic scale on the y-axis. Requires positive values."
    )
    if y_log and y_min <= 0:
        st.caption("Log y-axis ignored: y min must be > 0")

    st.markdown("---")

    # Show warning if parameters are the same
    if st.session_state.explore_xName == st.session_state.explore_yName:
        st.warning("⚠️ X and Y parameters are the same! Please select different parameters.")

    # Add helpful information box before run button
    with st.expander("📋 How Progress Monitoring Works", expanded=False):
        st.markdown("""
        When you click "Run Exploreogram":

        1. **This GUI window** will show:
           - Overall status ("Running...", "Complete!", etc.)
           - Estimated time
           - Configuration details
           - Final results when done

        2. **Your terminal/console window** will show:
           - Detailed model-by-model progress
           - Info messages for each model
           - Any warnings or errors

        **To see your terminal:**
        - **Mac**: Look for the Terminal app window where you typed `streamlit run`
        - **Windows**: Look for the Command Prompt or PowerShell window
        - **Linux**: Look for your terminal window

        **The GUI may appear frozen during computation - this is normal!**
        The models are running, and the GUI will update automatically when complete.
        """)

    # Show previous error if any
    if st.session_state.explore_error is not None:
        st.error(f"❌ Previous run failed: {st.session_state.explore_error}")

        # Show full traceback if available
        if st.session_state.explore_error_traceback is not None:
            with st.expander("📋 View Full Error Details and Diagnostics", expanded=True):
                st.code(st.session_state.explore_error_traceback)
                st.info("💡 The diagnostic messages above show exactly where the error occurred. Look for the last message before the error to identify which step failed.")

        if st.button("🔄 Clear Error and Try Again"):
            st.session_state.explore_error = None
            st.session_state.explore_error_traceback = None
            st.session_state.explore_running = False
            st.rerun()

    # Run / Force Re-run buttons
    col_run, col_force = st.columns([2, 1])
    with col_run:
        run_button = st.button("Run Exploreogram", type="primary", disabled=st.session_state.explore_running)
    with col_force:
        force_rerun = st.button("Force Re-run", disabled=st.session_state.explore_running,
                                help="Ignore cache and recompute all models")

# Cache utilities — computed once, at module level, from stable session-state
# widget values set above.  Must be outside both column blocks so col2 can
# reference current_cache_key without depending on col1's execution order.
from Utilities.explore_cache import generate_cache_key, get_cache_path, save_to_cache, load_from_cache

cache_dir = os.path.join(parent_directory, 'output', 'exploreograms', 'cache')
Planet = st.session_state["Planet"]

# Determine if induction is needed for current z-variable
needs_induction_for_z = z_param in induction_z_vars

# Generate cache key from current computational parameters
# Include excitation selection so changing frequencies triggers recompute
selected_exc_for_key = st.session_state.get('selected_excitations', []) if needs_induction_for_z else []
_x_driver_key = st.session_state.get('explore_x_driver', None)
_y_driver_key = st.session_state.get('explore_y_driver', None)
current_cache_key = generate_cache_key(
    bodyname=Planet.name,
    xName=st.session_state.explore_xName,
    xRange=[x_min, x_max],
    yName=st.session_state.explore_yName,
    yRange=[y_min, y_max],
    nx=nx, ny=ny,
    skip_induction=not needs_induction_for_z,
    exc_names=selected_exc_for_key,
)
# Fix F — append driver parameters so derived-axis combinations get distinct keys
if _x_driver_key or _y_driver_key:
    _driver_suffix = f"_xdrv{_x_driver_key or 'none'}_ydrv{_y_driver_key or 'none'}"
    current_cache_key += _driver_suffix.replace('/', '_').replace(' ', '_')

if (run_button or force_rerun) and not st.session_state.explore_running:
    if force_rerun:
        st.session_state.explore_force_rerun = True
        st.session_state.explore_results = None
        st.session_state.explore_cache_key = None
        # Fix D — delete the disk cache file so the force re-run truly recomputes
        from Utilities.explore_cache import get_cache_path as _get_cache_path_d
        _cache_file = _get_cache_path_d(cache_dir, current_cache_key)
        if os.path.isfile(_cache_file):
            os.remove(_cache_file)
    st.session_state.explore_running = True
    st.session_state.explore_start_time = None
    st.session_state.explore_error = None
    st.session_state.explore_error_traceback = None
    st.rerun()

with col1:
    # Check if we already have matching cached results (session or disk)
    has_session_cache = (
        st.session_state.explore_results is not None and
        st.session_state.explore_cache_key == current_cache_key
    )

    # Show cache status
    if has_session_cache:
        st.success("Cached results available — changing z-variable is instant")
    elif os.path.isfile(get_cache_path(cache_dir, current_cache_key)):
        st.info("Cached results found on disk — will load without recomputation")

with col2:
    st.subheader("📊 Exploration Results")

    if st.session_state.explore_running:
        # --- CACHE CHECK: try session cache, then disk cache, before computing ---
        cache_hit = False
        if not st.session_state.explore_force_rerun:
            # Session cache check
            if (st.session_state.explore_results is not None and
                    st.session_state.explore_cache_key == current_cache_key):
                cache_hit = True
                st.success("Using cached results (session) — no recomputation needed")
                Exploration = st.session_state.explore_results['Planet'].Exploration
                Params_cached = st.session_state.explore_results['Params']

            # Disk cache check
            if not cache_hit:
                cached_exploration = load_from_cache(cache_dir, current_cache_key)
                if cached_exploration is not None:
                    cache_hit = True
                    st.success("Loaded cached results from disk — no recomputation needed")
                    Exploration = cached_exploration
                    # Attach to Planet for downstream compatibility
                    Planet.Exploration = Exploration

                    # Build full Params for plotting (including Induct for matplotlib inductograms)
                    from configPP import configAssign
                    from PlanetProfile.MagneticInduction.defaultConfigInduct import inductAssign
                    Params_cached, ExploreParams = configAssign()
                    Params_cached.Explore = ExploreParams
                    Params_cached.Explore.xName = st.session_state.explore_xName
                    Params_cached.Explore.yName = st.session_state.explore_yName
                    Params_cached.Explore.zName = st.session_state.explore_zName
                    Params_cached.Explore.nx = int(nx)
                    Params_cached.Explore.ny = int(ny)

                    # Load induction params so matplotlib path can access excSelectionPlot
                    _, _, InductParams_cached, _ = inductAssign()
                    Params_cached.Induct = InductParams_cached

                    st.session_state.explore_results = {
                        'Planet': Planet,
                        'Exploration': Exploration,
                        'Params': Params_cached,
                        'xName': st.session_state.explore_xName,
                        'yName': st.session_state.explore_yName,
                        'zName': st.session_state.explore_zName,
                    }
                    st.session_state.explore_cache_key = current_cache_key

        if cache_hit:
            # Update z-variable for visualization
            Exploration.zName = st.session_state.explore_zName
            st.session_state.explore_running = False
            st.session_state.explore_force_rerun = False
            time.sleep(0.5)
            st.rerun()

        # Reset force flag
        st.session_state.explore_force_rerun = False

        # --- FULL COMPUTATION (cache miss) ---
        # Estimate time
        est_time_per_model = 8  # seconds (conservative estimate)
        est_total_minutes = (total_runs * est_time_per_model) / 60

        # Create prominent information box
        st.info(f"""
        ### Exploreogram Running: {total_runs} models ({nx} x {ny} grid)

        **Estimated time: ~{est_total_minutes:.1f} minutes** ({est_time_per_model} sec/model)

        **Progress Monitoring:**
        - This window shows overall status
        - **Look at your TERMINAL WINDOW** (where you ran `streamlit run`) for detailed model-by-model progress
        - The GUI will automatically update when complete

        The GUI may appear frozen during computation — this is normal.
        """)

        # Use st.status() for long-running operation
        with st.status(f"Computing {total_runs} models...", expanded=True) as status:
            if _DEBUG_EXPLOREOGRAM:
                st.write("⚙️ Using parallel processing on available CPU cores")
                st.write("")

            try:
                start_time = time.time()

                # Show configuration
                if _DEBUG_EXPLOREOGRAM:
                    st.write("**Configuration:**")
                config_info = st.container()
                with config_info:
                    col_a, col_b = st.columns(2)
                    with col_a:
                        if _DEBUG_EXPLOREOGRAM:
                            st.write(f"🌍 Planet: {st.session_state['Planet'].name if 'Planet' in st.session_state else 'Unknown'}")
                            st.write(f"📊 X: {ALL_PARAMS[st.session_state.explore_xName]['label']}")
                            st.write(f"   Range: [{x_min}, {x_max}]")
                    with col_b:
                        if _DEBUG_EXPLOREOGRAM:
                            st.write(f"📈 Y: {ALL_PARAMS[st.session_state.explore_yName]['label']}")
                            st.write(f"   Range: [{y_min}, {y_max}]")
                            st.write(f"🎨 Color: {Z_VARIABLES.get(st.session_state.explore_zName, st.session_state.explore_zName)}")

                if _DEBUG_EXPLOREOGRAM:
                    st.write("")
                    st.write("**Progress:**")

                # Import PlanetProfile modules
                if _DEBUG_EXPLOREOGRAM:
                    st.write("📦 Importing PlanetProfile modules...")
                from PlanetProfile.Main import ParPlanetExplore
                from PlanetProfile.Utilities.ResultsStructs import ExplorationResultsStruct
                from PlanetProfile.Utilities.ResultsIO import ExtractResults
                from PlanetProfile.Utilities.SummaryTables import GetLayerMeans
                from PlanetProfile.MagneticInduction.defaultConfigInduct import inductAssign
                from PlanetProfile.TrajecAnalysis.defaultConfigTrajec import trajecAssign
                from PlanetProfile.CustomSolution.defaultConfigCustomSolution import customSolutionAssign
                from PlanetProfile.Gravity.defaultConfigGravity import gravityAssign
                from PlanetProfile.MonteCarlo.defaultConfigMonteCarlo import montecarloAssign
                from PlanetProfile.Plotting.defaultConfigPlots import plotAssign
                from configPP import configAssign

                # Get current planet and params
                if _DEBUG_EXPLOREOGRAM:
                    st.write("⚙️ Configuring parameters...")
                Planet = st.session_state["Planet"]
                bodyname = Planet.name

                # Load default config (same way as RunPlanetProfile page)
                Params, ExploreParams = configAssign()
                Params.Explore = ExploreParams  # Assign exploration parameters to Params

                # Verify Explore attributes are properly initialized
                if not hasattr(Params.Explore, 'exploreType') or Params.Explore.exploreType is None:
                    st.error("❌ Error: exploreType not properly initialized in ExploreParams")
                    st.stop()
                if not hasattr(Params.Explore, 'exploreLogScale'):
                    Params.Explore.exploreLogScale = ['mixingRatioToH2O']  # Default value

                # Load all sub-configuration modules (required for full functionality)
                if _DEBUG_EXPLOREOGRAM:
                    st.write("⚙️ Loading configuration parameters...")
                SigParams, ExcSpecParams, InductParams, _ = inductAssign()
                TrajecParams = trajecAssign()
                CustomSolutionParams = customSolutionAssign()
                GravityParams = gravityAssign()
                MonteCarloParams = montecarloAssign()
                Color, Style, FigLbl, FigSize, FigMisc = plotAssign()

                # Assign all sub-configs to Params (as done in GetConfig.py)
                Params.Sig = SigParams
                Params.Induct = InductParams
                Params.MagSpectrum = ExcSpecParams
                Params.Trajec = TrajecParams
                Params.CustomSolution = CustomSolutionParams
                Params.Gravity = GravityParams
                Params.MonteCarlo = MonteCarloParams

                # Set up logging level for parallel processing (required by ParPlanetExplore)
                import logging
                import multiprocessing as mtp
                logging.PROFILE = logging.WARN + 5
                Params.logParallel = logging.PROFILE + 0
                if Params.QUIET:
                    Params.logParallel += 10

                # Suppress PyALMA warnings about parallel processing during exploreogram
                import warnings
                warnings.filterwarnings('ignore', message='.*parallel setting should be boolean.*')
                warnings.filterwarnings('ignore', message='.*Reverting to serial operation.*')

                # Set up parallel processing (required by ParPlanetExplore)
                if Params.DO_PARALLEL:
                    Params.maxCores = mtp.cpu_count()
                else:
                    Params.maxCores = 1

                # Initialize other required attributes (from GetConfig.py)
                Params.tStart_s = time.time()
                Params.fNameRef = {comp: f'{comp}Ref.txt' for comp in getattr(Params, 'wRef_ppt', {}).keys()}
                Params.Pref_MPa = {}
                Params.rhoRef_kgm3 = {}
                Params.nRef = {}
                Params.nRefPts = {}

                # Set exploration parameters
                Params.DO_EXPLOREOGRAM = True
                Params.CALC_NEW = True
                Params.EXPLOREOGRAM_IN_PROGRESS = True

                # For derived parameters, use the driver parameter as the actual
                # exploration variable — the derived value is read from results after the run
                x_driver = st.session_state.get('explore_x_driver', None)
                y_driver = st.session_state.get('explore_y_driver', None)
                x_display_name = st.session_state.explore_xName  # What user wants on the axis
                y_display_name = st.session_state.explore_yName

                Params.Explore.xName = x_driver if x_driver else x_display_name
                Params.Explore.yName = y_driver if y_driver else y_display_name
                Params.Explore.zName = st.session_state.explore_zName
                Params.Explore.xRange = [x_min, x_max]
                Params.Explore.yRange = [y_min, y_max]

                if _DEBUG_EXPLOREOGRAM:
                    if x_driver:
                        st.write(f"   Derived X-axis: {ALL_PARAMS[x_display_name]['label']} (varying {INPUT_PARAMS[x_driver]['label']})")
                    if y_driver:
                        st.write(f"   Derived Y-axis: {ALL_PARAMS[y_display_name]['label']} (varying {INPUT_PARAMS[y_driver]['label']})")

                # Configure induction calculations based on z-variable selection
                # Check if we need induction for the selected output variable
                needs_induction = st.session_state.explore_zName in induction_z_vars

                if needs_induction:
                    Params.SKIP_INDUCTION = False
                    Params.CALC_NEW_INDUCT = True

                    # Configure which excitations to calculate based on user selection
                    selected_exc = st.session_state.get('selected_excitations', ['orbital', 'synodic', 'true anomaly'])

                    # Set all to False first
                    for exc_name in Params.Induct.excSelectionCalc.keys():
                        Params.Induct.excSelectionCalc[exc_name] = False

                    # Enable only selected excitations
                    for exc_name in selected_exc:
                        if exc_name in Params.Induct.excSelectionCalc:
                            Params.Induct.excSelectionCalc[exc_name] = True

                    if _DEBUG_EXPLOREOGRAM:
                        st.write(f"   ✓ Magnetic induction enabled: {len(selected_exc)} excitation(s) selected")
                        st.write(f"     Frequencies: {', '.join(selected_exc)}")

                    # Preload excitation moments ONCE before exploreogram starts
                    # (During exploreogram, GetBexc is skipped in SetupInduction, so we need to load it here)
                    if _DEBUG_EXPLOREOGRAM:
                        st.write("   📡 Loading magnetic excitation moments...")
                    try:
                        from PlanetProfile.MagneticInduction.MagneticInduction import GetBexc

                        # Ensure pMax is initialized (defaults to 0 for spherically symmetric)
                        if Planet.Magnetic.pMax is None:
                            Planet.Magnetic.pMax = 0

                        # Load excitation moments using Planet's magnetic configuration
                        Planet.Magnetic.Texc_hr, Planet.Magnetic.omegaExc_radps, Planet.Magnetic.Benm_nT, \
                        Planet.Magnetic.B0_nT, Planet.Magnetic.Bexyz_nT = GetBexc(
                            Planet.name,
                            Planet.Magnetic.SCera,
                            Planet.Magnetic.extModel,
                            Params.Induct.excSelectionCalc,
                            nprmMax=Planet.Magnetic.nprmMax,
                            pMax=Planet.Magnetic.pMax
                        )

                        Planet.Magnetic.nExc = np.size(Planet.Magnetic.Texc_hr)

                        if Planet.Magnetic.Benm_nT is None:
                            st.error("❌ Failed to load excitation moments - file not found")
                            st.info(f"Expected file: Be1xyz_{Planet.name}_{Planet.Magnetic.SCera}_{Planet.Magnetic.extModel}noMP.txt")
                            st.info(f"in directory: {Planet.name}/inductionData/")
                            st.stop()
                        else:
                            # Initialize Binm_nT array with same shape as Benm_nT
                            if Planet.Magnetic.Binm_nT is None:
                                if isinstance(Planet.Magnetic.Benm_nT, dict):
                                    Planet.Magnetic.Binm_nT = {SCera: np.zeros_like(Benm_nT)
                                                               for SCera, Benm_nT in Planet.Magnetic.Benm_nT.items()}
                                else:
                                    Planet.Magnetic.Binm_nT = np.zeros_like(Planet.Magnetic.Benm_nT)

                            if _DEBUG_EXPLOREOGRAM:
                                st.write(f"   ✓ Excitation moments loaded: {Planet.Magnetic.nExc} frequency(ies)")
                    except Exception as e:
                        st.error(f"❌ Failed to load excitation moments: {e}")
                        import traceback
                        st.code(traceback.format_exc())
                        st.stop()

                else:
                    # Skip induction for non-magnetic variables to speed up computation
                    Params.SKIP_INDUCTION = True
                    Params.CALC_NEW_INDUCT = False
                    if _DEBUG_EXPLOREOGRAM:
                        st.write("   ℹ️ Skipping magnetic induction (not needed for selected output)")
                Params.Explore.nx = int(nx)
                Params.Explore.ny = int(ny)

                # GUI-specific settings to speed up runs
                Params.SKIP_PLOTS = True
                Params.NO_SAVEFILE = True  # Don't save individual runs
                Params.SKIP_INNER = False  # Keep inner calculations
                Params.ALLOW_BROKEN_MODELS = True  # Continue even if some models fail

                if _DEBUG_EXPLOREOGRAM:
                    st.write(f"✓ Configuration complete")
                    st.write("")

                # Generate exploration lists
                if _DEBUG_EXPLOREOGRAM:
                    st.write("📊 Generating parameter grid...")

                # Check if exploreLogScale exists and is not None
                if hasattr(Params.Explore, 'exploreLogScale') and Params.Explore.exploreLogScale is not None:
                    use_xlog = Params.Explore.xName in Params.Explore.exploreLogScale
                    use_ylog = Params.Explore.yName in Params.Explore.exploreLogScale
                else:
                    use_xlog = False
                    use_ylog = False

                if use_xlog:
                    xList = np.logspace(np.log10(x_min), np.log10(x_max), int(nx))
                else:
                    xList = np.linspace(x_min, x_max, int(nx))

                if use_ylog:
                    yList = np.logspace(np.log10(y_min), np.log10(y_max), int(ny))
                else:
                    yList = np.linspace(y_min, y_max, int(ny))

                Params.nModels = int(nx) * int(ny)

                # Create exploration results structure
                Exploration = ExplorationResultsStruct()
                Exploration.xName = Params.Explore.xName
                Exploration.yName = Params.Explore.yName
                Exploration.zName = Params.Explore.zName
                Exploration.xData = xList
                Exploration.yData = yList

                if _DEBUG_EXPLOREOGRAM:
                    st.write(f"✓ Grid created: {len(xList)}×{len(yList)} = {total_runs} models")
                    st.write("")

                # Run parallel exploration
                if _DEBUG_EXPLOREOGRAM:
                    st.write(f"🔬 **Running {total_runs} models in parallel...**")
                    st.write(f"   Expected completion: ~{est_total_minutes:.1f} minutes")
                    st.write("")

                # Create progress display
                progress_container = st.container()
                with progress_container:
                    progress_bar = st.progress(0, text="Starting exploration...")
                    status_placeholder = st.empty()

                    # Show terminal reminder
                    st.info(f"""
                    💡 **Check the Terminal/Console for Live Progress!**

                    Look at the terminal window where you launched Streamlit to see real-time model-by-model progress like:

                    `[PROFILE] Profile 1/{total_runs} complete...`
                    `[PROFILE] Profile 2/{total_runs} complete...`
                    `[PROFILE] Profile 3/{total_runs} complete...`

                    ⏱️ **Estimated time: ~{est_total_minutes:.1f} minutes** ({total_runs} models × ~{est_time_per_model:.1f} sec/model)

                    ⚠️ Note: This browser window will appear frozen during computation - this is normal!
                    The page will update automatically when all models complete.
                    """)

                # Start timer
                start_time_explore = time.time()
                progress_bar.progress(5, text=f"Computing {total_runs} models in parallel...")

                # Call the parallel exploration function directly (bypasses file loading)
                try:
                    # Debug: Verify critical attributes before calling ParPlanetExplore
                    if _DEBUG_EXPLOREOGRAM:
                        st.write("🔍 Verifying configuration...")

                    # Check exploreType
                    if not hasattr(Params.Explore, 'exploreType'):
                        st.error("❌ Params.Explore.exploreType attribute missing!")
                        raise AttributeError("exploreType not found in Params.Explore")

                    if Params.Explore.exploreType is None:
                        st.error("❌ Params.Explore.exploreType is None!")
                        raise ValueError("exploreType cannot be None")

                    # Check if xName and yName exist in exploreType
                    if Params.Explore.xName not in Params.Explore.exploreType:
                        st.error(f"❌ X parameter '{Params.Explore.xName}' not found in exploreType dict!")
                        st.error(f"Available keys: {list(Params.Explore.exploreType.keys())}")
                        raise KeyError(f"xName '{Params.Explore.xName}' not in exploreType")

                    if Params.Explore.yName not in Params.Explore.exploreType:
                        st.error(f"❌ Y parameter '{Params.Explore.yName}' not found in exploreType dict!")
                        st.error(f"Available keys: {list(Params.Explore.exploreType.keys())}")
                        raise KeyError(f"yName '{Params.Explore.yName}' not in exploreType")

                    if _DEBUG_EXPLOREOGRAM:
                        st.write(f"✓ Configuration verified: X={Params.Explore.xName} ({Params.Explore.exploreType[Params.Explore.xName]}), Y={Params.Explore.yName} ({Params.Explore.exploreType[Params.Explore.yName]})")

                    # Update progress before starting
                    status_placeholder.success("✅ Starting parallel computation...")
                    progress_bar.progress(10, text=f"Running {total_runs} models in parallel...")

                    # Run the exploration (this blocks until complete)
                    PlanetGrid = ParPlanetExplore(Planet, Params, xList, yList)

                    # Calculate elapsed time
                    elapsed_explore = time.time() - start_time_explore
                    progress_bar.progress(100, text=f"✅ All {total_runs} models complete!")
                    status_placeholder.success(f"⏱️ Computation complete in {elapsed_explore:.1f} seconds ({elapsed_explore/total_runs:.1f} sec/model)")

                except (TypeError, KeyError, AttributeError) as e:
                    st.error(f"❌ Error during exploration: {type(e).__name__}: {e}")
                    st.error("**Debug Information:**")
                    st.error(f"  - xName: {getattr(Params.Explore, 'xName', 'NOT SET')}")
                    st.error(f"  - yName: {getattr(Params.Explore, 'yName', 'NOT SET')}")
                    st.error(f"  - exploreType exists: {hasattr(Params.Explore, 'exploreType')}")
                    if hasattr(Params.Explore, 'exploreType'):
                        st.error(f"  - exploreType is None: {Params.Explore.exploreType is None}")
                        if Params.Explore.exploreType is not None:
                            st.error(f"  - exploreType keys: {list(Params.Explore.exploreType.keys())}")

                    import traceback
                    with st.expander("📋 Full Error Traceback"):
                        st.code(traceback.format_exc())
                    raise

                # Update progress for post-processing
                progress_bar.progress(60, text="📐 Calculating layer means...")
                status_placeholder.info("Processing completed models...")

                # Calculate additional parameters from profiles
                try:
                    if _DEBUG_EXPLOREOGRAM:
                        st.write(f"  - Grid shape: {np.shape(PlanetGrid)}")
                    gridShape = np.shape(PlanetGrid)
                    if _DEBUG_EXPLOREOGRAM:
                        st.write(f"  - Flattening grid...")
                    PlanetList = np.reshape(PlanetGrid, -1)
                    if _DEBUG_EXPLOREOGRAM:
                        st.write(f"  - Calling GetLayerMeans with {len(PlanetList)} planets...")
                    PlanetList, Params = GetLayerMeans(PlanetList, Params)
                    if _DEBUG_EXPLOREOGRAM:
                        st.write(f"  - Reshaping back to grid...")
                    PlanetGrid = np.reshape(PlanetList, gridShape)
                    if _DEBUG_EXPLOREOGRAM:
                        st.write(f"  ✓ Layer means calculated")
                except Exception as e:
                    st.error(f"❌ Error in GetLayerMeans: {type(e).__name__}: {e}")
                    st.error(f"  - Error type: {type(e)}")
                    st.error(f"  - PlanetGrid type: {type(PlanetGrid)}")
                    st.error(f"  - PlanetGrid shape: {np.shape(PlanetGrid)}")
                    if 'PlanetList' in locals():
                        st.error(f"  - PlanetList length: {len(PlanetList)}")
                    import traceback
                    with st.expander("📋 GetLayerMeans Full Error"):
                        st.code(traceback.format_exc())
                    raise

                # Extract results into structured format
                progress_bar.progress(80, text="📊 Extracting results...")
                try:
                    if _DEBUG_EXPLOREOGRAM:
                        st.write("  - Checking Exploration object...")
                        st.write(f"    - xData: {Exploration.xData is not None} (shape: {np.shape(Exploration.xData)})")
                        st.write(f"    - yData: {Exploration.yData is not None} (shape: {np.shape(Exploration.yData)})")
                        st.write(f"    - zName: {Exploration.zName}")
                        st.write(f"  - PlanetGrid shape: {np.shape(PlanetGrid)}")

                    # Check if magnetic induction was actually calculated
                    if needs_induction:
                        has_benm = (hasattr(PlanetGrid[0, 0], 'Magnetic') and
                                   hasattr(PlanetGrid[0, 0].Magnetic, 'Benm_nT') and
                                   PlanetGrid[0, 0].Magnetic.Benm_nT is not None)
                        if not has_benm:
                            st.warning("⚠️ Magnetic induction data not found in results. The selected z-variable may show zeros/NaN.")

                    if _DEBUG_EXPLOREOGRAM:
                        st.write(f"  - Calling ExtractResults...")

                    Exploration = ExtractResults(Exploration, PlanetGrid, Params)

                    # For derived parameters: do NOT overwrite Exploration.xData/yData.
                    # Keep the rectilinear driver arrays intact for the Plotly renderer.
                    # The derived 2D arrays are passed separately for tick-label substitution only.
                    if _DEBUG_EXPLOREOGRAM:
                        if x_driver and hasattr(Exploration.base, x_display_name):
                            derived_x_dbg = getattr(Exploration.base, x_display_name)
                            st.write(f"  ✓ Derived X available: {x_display_name} (range: {np.nanmin(derived_x_dbg):.2f} – {np.nanmax(derived_x_dbg):.2f})")
                        if y_driver and hasattr(Exploration.base, y_display_name):
                            derived_y_dbg = getattr(Exploration.base, y_display_name)
                            st.write(f"  ✓ Derived Y available: {y_display_name} (range: {np.nanmin(derived_y_dbg):.2f} – {np.nanmax(derived_y_dbg):.2f})")

                    if _DEBUG_EXPLOREOGRAM:
                        st.write(f"  ✓ Results extracted successfully")
                        st.write(f"  - Exploration.base exists: {hasattr(Exploration, 'base')}")
                        if hasattr(Exploration, 'base'):
                            st.write(f"  - Exploration.base.VALID exists: {hasattr(Exploration.base, 'VALID')}")
                            if hasattr(Exploration.base, 'VALID'):
                                st.write(f"  - VALID is not None: {Exploration.base.VALID is not None}")

                except Exception as e:
                    st.error(f"❌ Error in ExtractResults: {type(e).__name__}: {e}")
                    st.error("**Debug Info:**")
                    st.error(f"  - Exploration type: {type(Exploration)}")
                    st.error(f"  - Exploration.xData exists: {hasattr(Exploration, 'xData')}")
                    if hasattr(Exploration, 'xData'):
                        st.error(f"  - Exploration.xData is None: {Exploration.xData is None}")
                    st.error(f"  - Exploration.yData exists: {hasattr(Exploration, 'yData')}")
                    if hasattr(Exploration, 'yData'):
                        st.error(f"  - Exploration.yData is None: {Exploration.yData is None}")
                    st.error(f"  - Exploration.zName: {getattr(Exploration, 'zName', 'NOT SET')}")
                    st.error(f"  - PlanetGrid type: {type(PlanetGrid)}")
                    st.error(f"  - PlanetGrid shape: {np.shape(PlanetGrid)}")
                    import traceback
                    with st.expander("📋 ExtractResults Full Error"):
                        st.code(traceback.format_exc())
                    raise

                # Store exploration results in the first Planet object
                progress_bar.progress(90, text="💾 Storing results...")
                try:
                    if _DEBUG_EXPLOREOGRAM:
                        st.write("  - Attaching Exploration to Planet...")
                    Planet.Exploration = Exploration
                    Params.EXPLOREOGRAM_IN_PROGRESS = False

                    # Check for valid models
                    if _DEBUG_EXPLOREOGRAM:
                        st.write("  - Checking VALID array...")
                        st.write(f"    - Exploration.base exists: {hasattr(Exploration, 'base')}")
                        st.write(f"    - Exploration.base is not None: {Exploration.base is not None}")
                        if Exploration.base is not None:
                            st.write(f"    - VALID exists: {hasattr(Exploration.base, 'VALID')}")
                            if hasattr(Exploration.base, 'VALID'):
                                st.write(f"    - VALID is not None: {Exploration.base.VALID is not None}")
                                if Exploration.base.VALID is not None:
                                    st.write(f"    - VALID type: {type(Exploration.base.VALID)}")
                                    st.write(f"    - VALID shape: {np.shape(Exploration.base.VALID)}")
                    # These checks are always active (not debug-only)
                    if Exploration.base is not None and hasattr(Exploration.base, 'VALID'):
                        if Exploration.base.VALID is not None:
                            if not np.any(Exploration.base.VALID):
                                st.warning('⚠️ No valid models found for the given input settings.')
                        else:
                            st.error("❌ Exploration.base.VALID is None!")
                    elif Exploration.base is not None:
                        st.error("❌ Exploration.base does not have VALID attribute!")
                    else:
                        st.error("❌ Exploration.base is None!")

                    # Auto-select the strongest induction component for next render
                    st.session_state['exploreogram_induct_component_default'] = _strongest_induction_axis(
                        getattr(Exploration, 'induction', None)
                    )

                    if _DEBUG_EXPLOREOGRAM:
                        st.write("  - Storing to session state and disk cache...")
                    st.session_state.explore_results = {
                        'Planet': Planet,
                        'Exploration': Exploration,
                        'Params': Params,
                        'xName': st.session_state.explore_xName,
                        'yName': st.session_state.explore_yName,
                        'zName': st.session_state.explore_zName,
                    }
                    st.session_state.explore_cache_key = current_cache_key

                    # Save to disk cache for cross-session persistence
                    try:
                        cache_path = save_to_cache(Exploration, cache_dir, current_cache_key)
                        if _DEBUG_EXPLOREOGRAM:
                            st.write(f"  - Cached to disk: {os.path.basename(cache_path)}")
                    except Exception as e:
                        st.warning(f"  - Disk cache save failed (non-critical): {e}")

                    if _DEBUG_EXPLOREOGRAM:
                        st.write("  - Calculating statistics...")
                    elapsed_time = time.time() - start_time
                    if _DEBUG_EXPLOREOGRAM:
                        st.write(f"    - Elapsed: {elapsed_time:.1f} sec")

                    # Calculate valid count with safety check
                    if Exploration.base is not None and hasattr(Exploration.base, 'VALID') and Exploration.base.VALID is not None:
                        if _DEBUG_EXPLOREOGRAM:
                            st.write(f"    - Summing VALID array...")
                        valid_count = np.sum(Exploration.base.VALID)
                        if _DEBUG_EXPLOREOGRAM:
                            st.write(f"    - Valid count: {valid_count}/{total_runs}")
                    else:
                        st.error("❌ Cannot calculate valid count - VALID array is None or missing!")
                        valid_count = 0

                except Exception as e:
                    st.error(f"❌ Error during storage/statistics: {type(e).__name__}: {e}")
                    import traceback
                    with st.expander("📋 Storage Error Traceback"):
                        st.code(traceback.format_exc())
                    raise

                # Complete progress bar
                progress_bar.progress(100, text="✅ Exploreogram complete!")
                status_placeholder.success(f"✅ Completed in {elapsed_time/60:.1f} minutes - {valid_count}/{total_runs} valid models")

                if _DEBUG_EXPLOREOGRAM:
                    st.write("")
                st.success(f"**✅ Exploreogram Complete!**")
                col_time1, col_time2 = st.columns(2)
                with col_time1:
                    st.metric("Total Time", f"{elapsed_time/60:.1f} min")
                with col_time2:
                    st.metric("Time/Model", f"{elapsed_time/total_runs:.1f} sec")

                if _DEBUG_EXPLOREOGRAM:
                    st.write(f"**Valid models:** {valid_count}/{total_runs}")

                # Validity warnings
                if valid_count == 0:
                    st.error("❌ No valid models produced. Try adjusting parameter ranges.")
                elif valid_count < total_runs * 0.5:
                    st.warning(f"⚠️ Only {valid_count}/{total_runs} models were valid. Consider adjusting parameter ranges.")

                # Save configuration file as a record
                output_dir = os.path.join(parent_directory, 'output', 'exploreograms')
                os.makedirs(output_dir, exist_ok=True)

                timestamp = time.strftime('%Y%m%d_%H%M%S')
                config_filename = f"{bodyname}_explore_config_{timestamp}.txt"
                config_filepath = os.path.join(output_dir, config_filename)

                with open(config_filepath, 'w') as f:
                    f.write(f"# PlanetProfile Exploreogram Configuration\n")
                    f.write(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                    f.write(f"Planet: {bodyname}\n")
                    f.write(f"X Parameter: {Params.Explore.xName} ({ALL_PARAMS[Params.Explore.xName]['label']})\n")
                    f.write(f"X Range: [{x_min}, {x_max}]\n")
                    f.write(f"X Points: {nx}\n\n")
                    f.write(f"Y Parameter: {Params.Explore.yName} ({ALL_PARAMS[Params.Explore.yName]['label']})\n")
                    f.write(f"Y Range: [{y_min}, {y_max}]\n")
                    f.write(f"Y Points: {ny}\n\n")
                    f.write(f"Color Variable: {Params.Explore.zName} ({Z_VARIABLES.get(Params.Explore.zName, Params.Explore.zName)})\n\n")
                    f.write(f"Total Models: {total_runs}\n")
                    f.write(f"Valid Models: {valid_count}\n")
                    f.write(f"Runtime: {elapsed_time/60:.2f} minutes\n")
                    f.write(f"Time per Model: {elapsed_time/total_runs:.2f} seconds\n")

                st.info(f"📝 Configuration saved: `{config_filename}`")

                status.update(label="Exploreogram complete!", state="complete", expanded=False)

            except Exception as e:
                import traceback
                error_traceback = traceback.format_exc()
                st.session_state.explore_error = str(e)
                st.session_state.explore_error_traceback = error_traceback
                st.error(f"❌ Error: {type(e).__name__}: {str(e)}")
                with st.expander("📋 Full Error Traceback", expanded=True):
                    st.code(error_traceback)
                status.update(label="Error occurred", state="error", expanded=True)
            finally:
                st.session_state.explore_running = False
                time.sleep(1)  # Brief pause to show completion
                st.rerun()

    # Display results if available
    if st.session_state.explore_results is not None and not st.session_state.explore_running:
        results = st.session_state.explore_results

        if hasattr(results['Planet'], 'Exploration'):
            Exploration = results['Planet'].Exploration
            Params = results['Params']

            # Extract data for plotting
            try:
                xData = Exploration.xData
                yData = Exploration.yData
                zName = results['zName']

                # Verify data is not None
                if xData is None:
                    st.error("❌ xData is None - cannot plot results")
                    st.stop()
                if yData is None:
                    st.error("❌ yData is None - cannot plot results")
                    st.stop()
                if zName is None:
                    st.error("❌ zName is None - cannot plot results")
                    st.stop()

            except Exception as e:
                st.error(f"Error extracting plot data: {e}")
                import traceback
                st.code(traceback.format_exc())
                st.stop()

            # --- INDUCTION GAP DETECTION ---
            if zName in induction_z_vars:
                has_induction_data = (
                    hasattr(Exploration, 'induction') and
                    Exploration.induction is not None and
                    getattr(Exploration.induction, 'nPeaks', 0) > 0
                )
                if not has_induction_data:
                    st.warning("Magnetic induction data not available in cached results")
                    st.info(
                        "This exploreogram was run without induction calculations. "
                        "To view magnetic field data, re-run with induction enabled."
                    )
                    if st.button("Re-run with Induction Enabled"):
                        st.session_state.explore_results = None
                        st.session_state.explore_cache_key = None
                        st.session_state.explore_force_rerun = True
                        st.session_state.explore_running = True
                        st.rerun()
                    st.stop()

            # --- Derived-axis decoration (Plotly only) ---
            # Compute derived-axis arrays for tick-label substitution without
            # mutating Exploration.xData / Exploration.yData / Exploration.xName / Exploration.yName.
            # Use session-state keys so this works both during and after a run.
            _plot_x_driver = st.session_state.get('explore_x_driver', None)
            _plot_y_driver = st.session_state.get('explore_y_driver', None)
            _plot_x_param  = st.session_state.get('explore_xName', None)
            _plot_y_param  = st.session_state.get('explore_yName', None)

            x_derived_arr   = None
            x_derived_label = None
            y_derived_arr   = None
            y_derived_label = None

            if _plot_x_driver is not None and _plot_x_param is not None:
                _xd = getattr(Exploration.base, _plot_x_param, None)
                if _xd is not None and np.shape(_xd) == np.shape(Exploration.xData):
                    x_derived_arr   = _xd
                    x_derived_label = _plot_x_param

            if _plot_y_driver is not None and _plot_y_param is not None:
                _yd = getattr(Exploration.base, _plot_y_param, None)
                if _yd is not None and np.shape(_yd) == np.shape(Exploration.yData):
                    y_derived_arr   = _yd
                    y_derived_label = _plot_y_param

            # --- PLOTTING ---
            if use_interactive_plots:
                # Use enhanced Plotly plotting
                try:
                    from Utilities.exploreogram_plotly import (
                        create_exploreogram_plotly,
                        create_inductogram_plotly,
                        INDUCTION_Z_VARS,
                    )

                    # Try to load FigLbl for proper styling
                    try:
                        from PlanetProfile.GetConfig import FigLbl
                        FigLbl.SetExploration(results['Planet'].name, Exploration.xName,
                                            Exploration.yName, Exploration.zName)
                    except Exception:
                        FigLbl = None

                    # Check if this is a magnetic induction variable
                    is_induction = (zName in induction_z_vars and
                                   hasattr(Exploration, 'induction') and
                                   Exploration.induction is not None and
                                   getattr(Exploration.induction, 'nPeaks', 0) > 0)

                    if is_induction:
                        # Multi-frequency inductogram display
                        display_mode = st.session_state.get('induct_display_mode', 'real_imaginary')
                        use_contours = st.session_state.get('induct_use_contours', True)
                        use_nT = st.session_state.get('exploreogram_use_nT_amplitudes', True)

                        nPeaks = Exploration.induction.nPeaks
                        freq_names = list(Exploration.induction.calcedExc) if (hasattr(Exploration, 'induction') and Exploration.induction is not None and getattr(Exploration.induction, 'calcedExc', None) is not None) else []
                        st.info(f"Showing inductogram for {nPeaks} frequency peak(s): {', '.join(freq_names) if freq_names else 'unknown'}")

                        show_sal = st.session_state.get('exploreogram_show_salinity_axis', True)
                        _ic_ind = st.session_state.get('exploreogram_induct_component', 'Total')
                        fig = create_inductogram_plotly(
                            Exploration, Params, FigLbl=FigLbl,
                            display_mode=display_mode,
                            use_contours=use_contours,
                            use_nT_amplitudes=use_nT,
                            x_log=x_log,
                            y_log=y_log,
                            show_salinity_axis=show_sal,
                            induct_component=_ic_ind,
                        )
                        st.plotly_chart(fig, use_container_width=True,
                            key=f"inductogram_{zName}_{display_mode}_{_ic_ind}_nT{use_nT}_contours{use_contours}_xlog{x_log}_ylog{y_log}")
                    else:
                        # Non-induction variable — standard exploreogram heatmap
                        if zName in induction_z_vars:
                            st.warning("Magnetic induction data not found. Make sure induction calculations were enabled.")

                        use_contours = st.session_state.get('induct_use_contours', True)
                        use_nT = st.session_state.get('exploreogram_use_nT_amplitudes', True)
                        show_sal = st.session_state.get('exploreogram_show_salinity_axis', True)
                        _ic = st.session_state.get('exploreogram_induct_component', 'Total')
                        fig = create_exploreogram_plotly(
                            Exploration, Params, FigLbl=FigLbl,
                            smoothing=False, smooth_factor=2,
                            use_contours=use_contours,
                            x_log=x_log,
                            y_log=y_log,
                            use_nT_amplitudes=use_nT,
                            show_salinity_axis=show_sal,
                            induct_component=_ic,
                            x_derived=x_derived_arr,
                            y_derived=y_derived_arr,
                            x_derived_label=x_derived_label,
                            y_derived_label=y_derived_label,
                        )
                        st.plotly_chart(fig, use_container_width=True,
                            key=f"exploreogram_{zName}_{_ic}_nT{use_nT}_contours{use_contours}_xlog{x_log}_ylog{y_log}")

                except Exception as e:
                    st.error(f"Error creating interactive plot: {e}")
                    import traceback
                    with st.expander("Show error details"):
                        st.code(traceback.format_exc())
                    st.info("Falling back to basic Plotly plot...")

                    # Fallback to basic plot
                    if hasattr(Exploration.base, zName):
                        zData = getattr(Exploration.base, zName)
                    else:
                        zData = np.zeros_like(Exploration.base.VALID, dtype=float)

                    # Get axis labels safely
                    xLabel = ALL_PARAMS.get(results["xName"], {}).get("label", results["xName"])
                    yLabel = ALL_PARAMS.get(results["yName"], {}).get("label", results["yName"])
                    zLabel = Z_VARIABLES.get(zName, zName)

                    # Extract 1D axis values from 2D grids
                    x1d = xData[:, 0] if len(xData.shape) == 2 else xData
                    y1d = yData[0, :] if len(yData.shape) == 2 else yData

                    fig = go.Figure(data=go.Heatmap(
                        x=x1d,
                        y=y1d,
                        z=zData.T,
                        colorscale='Viridis',
                        colorbar=dict(title=zLabel),
                    ))

                    fig.update_layout(
                        title=f'{zLabel} Exploreogram',
                        xaxis_title=xLabel,
                        yaxis_title=yLabel,
                        width=800,
                        height=600,
                    )

                    st.plotly_chart(fig, use_container_width=True,
                        key=f"exploreogram_fallback_{zName}")

            else:
                # Use matplotlib plots via GenerateExplorationPlots
                if zName == 'Amp_nT' and st.session_state.get('exploreogram_induct_component', 'Total') != 'Total':
                    st.info("Component selection (Bx/By/Bz) is only supported in the interactive Plotly view. The matplotlib plot shows the total amplitude.")
                st.info("📊 Generating matplotlib plots... This may take a moment.")

                try:
                    import matplotlib
                    matplotlib.use('Agg')  # Non-interactive backend
                    import matplotlib.pyplot as plt
                    from PlanetProfile.Plotting.ExplorationPlots import GenerateExplorationPlots
                    from PlanetProfile.GetConfig import FigMisc
                    from pdf2image import convert_from_path

                    # Prepare for plotting
                    ExplorationList = [Exploration]

                    # Setup figure files
                    from PlanetProfile.Utilities.defineStructs import FigureFilesSubstruct
                    output_dir = os.path.join(parent_directory, 'output', 'exploreograms')
                    os.makedirs(output_dir, exist_ok=True)

                    fig_basename = f"{results['Planet'].name}_explore_{results['xName']}_vs_{results['yName']}"
                    _ic_mpl = st.session_state.get('exploreogram_induct_component', 'Total')
                    _pdf_append = results['zName'] if _ic_mpl == 'Total' else f"{results['zName']}_{_ic_mpl}"
                    FigureFiles = FigureFilesSubstruct(
                        figPath=output_dir,
                        figBase=fig_basename,
                        xtn='.pdf',
                        exploreAppend=_pdf_append
                    )

                    FigureFilesList = [FigureFiles]

                    # Assign FigureFiles to Params so GenerateExplorationPlots can access it
                    Params.FigureFiles = FigureFiles

                    # Ensure contourName is set
                    if not hasattr(Params.Explore, 'contourName'):
                        Params.Explore.contourName = None

                    # Always reset excName — stale induction excitation names
                    # from prior runs would otherwise append to non-induction labels
                    Exploration.excName = None
                    if not hasattr(Exploration, 'contourName'):
                        Exploration.contourName = Params.Explore.contourName

                    # Configure for induction multi-frequency display
                    is_induction_z = zName in induction_z_vars
                    if is_induction_z and hasattr(Exploration, 'induction') and Exploration.induction is not None:
                        # Read display mode from GUI — controls real+imag vs amp+phase
                        display_mode = st.session_state.get('induct_display_mode', 'real_imaginary')

                        # Map display_mode to appropriate z-variable list
                        if display_mode == 'amplitude_phase':
                            # Force amplitude + phase subplots regardless of selected zName
                            Params.Explore.zName = ['Amp_nT', 'phase_deg']
                        elif display_mode == 'real_imaginary':
                            # If user selected Amp or phase, substitute Bi1Tot for real+imag expansion
                            if zName in ('Amp_nT', 'phase_deg'):
                                Params.Explore.zName = ['Bi1Tot_nT']
                            else:
                                # zName should be in zNamePlotRealImag for automatic expansion
                                Params.Explore.zName = [zName]
                        else:
                            Params.Explore.zName = [zName]

                        Params.PLOT_COMBO_EXPLORATIONS = True

                        # Ensure Params.Induct is initialized (may be None if loaded from disk cache)
                        if Params.Induct is None:
                            try:
                                from PlanetProfile.MagneticInduction.defaultConfigInduct import inductAssign
                                _, _, InductParams_display, _ = inductAssign()
                                Params.Induct = InductParams_display
                            except Exception:
                                from PlanetProfile.Utilities.defineStructs import InductOgramParamsStruct
                                Params.Induct = InductOgramParamsStruct(inductOtype=None, cLevels=None, dftC=None, cfmt=None)

                        # Configure excSelectionPlot to match user-selected frequencies
                        selected_exc = st.session_state.get('selected_excitations', [])
                        # Ensure all selected excitations have entries in excSelectionPlot
                        for exc_name in selected_exc:
                            if exc_name not in Params.Induct.excSelectionPlot:
                                Params.Induct.excSelectionPlot[exc_name] = True
                        for exc_name in Params.Induct.excSelectionPlot.keys():
                            Params.Induct.excSelectionPlot[exc_name] = exc_name in selected_exc

                        # Rebuild FigureFiles with list-format exploreAppend
                        FigureFiles = FigureFilesSubstruct(
                            figPath=output_dir,
                            figBase=fig_basename,
                            xtn='.pdf',
                            exploreAppend=Params.Explore.zName
                        )
                        FigureFilesList = [FigureFiles]
                        Params.FigureFiles = FigureFiles
                    else:
                        # Non-induction: single variable mode
                        Params.Explore.zName = zName
                        Params.PLOT_COMBO_EXPLORATIONS = False

                    # Plot flags
                    if not hasattr(Params, 'PLOT_Zb_D'):
                        Params.PLOT_Zb_D = False
                    if not hasattr(Params, 'PLOT_D_SIGMA'):
                        Params.PLOT_D_SIGMA = False
                    if not hasattr(Params, 'PLOT_LOVE_COMPARISON'):
                        Params.PLOT_LOVE_COMPARISON = False

                    # Propagate contour toggle from GUI checkbox
                    Params.Explore.DRAW_CONTOURS = st.session_state.get('induct_use_contours', True)

                    # --- GUI log-scale toggles → matplotlib (Task 3) ---
                    # FigLbl.SetExploration() (called inside GenerateExplorationPlots) reads
                    # FigLbl.axisLogScalesExplore to decide xScaleExplore/yScaleExplore.
                    # Force the user's checkbox state by adding/removing xName/yName from
                    # that list, then restore in a finally-block below.
                    from PlanetProfile.GetConfig import FigLbl as _FigLbl_modlog
                    _orig_axisLogScales = list(_FigLbl_modlog.axisLogScalesExplore)
                    _xname_for_scale = Exploration.xName
                    _yname_for_scale = Exploration.yName
                    _new_log_list = [v for v in _orig_axisLogScales
                                     if v not in (_xname_for_scale, _yname_for_scale)]
                    if x_log:
                        _new_log_list.append(_xname_for_scale)
                    if y_log:
                        _new_log_list.append(_yname_for_scale)
                    _FigLbl_modlog.axisLogScalesExplore = _new_log_list

                    # --- GUI nT-amplitudes toggle → matplotlib (Task 2 / Bug A) ---
                    # The matplotlib pipeline reads:
                    #   * induction.Amp           — for Amp_nT (amplitude+phase mode)
                    #   * induction.rBi1Tot_nT    — for the real subplot in R+I mode
                    #   * induction.iBi1Tot_nT    — for the imaginary subplot in R+I mode
                    # Default (use_nT=True) renders in nT.  When unchecked we
                    # temporarily override these arrays with the dimensionless
                    # Aₑ representation:
                    #     Amp        ← Aₑ                   (already dimensionless)
                    #     rBi1Tot_nT ← Aₑ·cos(φ)            (Re component, dimensionless)
                    #     iBi1Tot_nT ← Aₑ·sin(φ)            (Im component, dimensionless)
                    # The originals are restored in a finally-block below.
                    use_nT_mpl = st.session_state.get('exploreogram_use_nT_amplitudes', True)
                    _amp_orig = None
                    _r_orig = None
                    _i_orig = None
                    _orig_amp_label = None
                    _orig_r_label = None
                    _orig_i_label = None
                    _amp_was_overridden = False
                    _ri_was_overridden = False
                    # Determine which induction arrays matplotlib will actually read,
                    # so we only override what is consumed.  Params.Explore.zName has
                    # already been resolved to the correct list for this display mode.
                    _plotted_z = list(Params.Explore.zName) if isinstance(Params.Explore.zName, (list, tuple)) else [Params.Explore.zName]
                    _will_plot_amp     = any(z in _plotted_z for z in ('Amp_nT', 'Bi1Tot_nT', 'Bi1x_nT', 'Bi1y_nT', 'Bi1z_nT'))
                    _will_plot_ri      = any(z in _plotted_z for z in ('Bi1Tot_nT', 'rBi1Tot_nT', 'iBi1Tot_nT',
                                                                        'rBi1x_nT', 'iBi1x_nT', 'rBi1y_nT', 'iBi1y_nT',
                                                                        'rBi1z_nT', 'iBi1z_nT'))

                    _ind = getattr(Exploration, 'induction', None)
                    if _ind is not None and getattr(_ind, 'Bi1Tot_nT', None) is not None:
                        if use_nT_mpl and _will_plot_amp:
                            # nT mode: replace dimensionless Aₑ modulus with |Bⁱ| in nT.
                            _amp_orig = _ind.Amp
                            _ind.Amp = np.abs(_ind.Bi1Tot_nT)
                            _amp_was_overridden = True
                        elif not use_nT_mpl and _will_plot_ri:
                            # Dimensionless Aₑ mode: override r/i components so matplotlib
                            # subplots show Re{Aₑ} and Im{Aₑ} instead of nT values.
                            _Amp_dimless = getattr(_ind, 'Amp', None)
                            _Phase_deg   = getattr(_ind, 'Phase', None)
                            if _Amp_dimless is not None and _Phase_deg is not None:
                                _phase_rad = np.deg2rad(_Phase_deg)
                                _r_orig = getattr(_ind, 'rBi1Tot_nT', None)
                                _i_orig = getattr(_ind, 'iBi1Tot_nT', None)
                                _ind.rBi1Tot_nT = _Amp_dimless * np.cos(_phase_rad)
                                _ind.iBi1Tot_nT = _Amp_dimless * np.sin(_phase_rad)
                                _ri_was_overridden = True

                        # Update axisLabelsExplore so colorbars/titles match the unit
                        try:
                            if (hasattr(_FigLbl_modlog, 'axisLabelsExplore') and
                                    _FigLbl_modlog.axisLabelsExplore is not None):
                                lbls = _FigLbl_modlog.axisLabelsExplore
                                if use_nT_mpl:
                                    if 'Amp_nT' in lbls:
                                        _orig_amp_label = lbls['Amp_nT']
                                        lbls['Amp_nT'] = r'$|B^i|$ (nT)'
                                else:
                                    if 'rBi1Tot_nT' in lbls:
                                        _orig_r_label = lbls['rBi1Tot_nT']
                                        lbls['rBi1Tot_nT'] = r'$\mathrm{Re}\{A_e\}$'
                                    if 'iBi1Tot_nT' in lbls:
                                        _orig_i_label = lbls['iBi1Tot_nT']
                                        lbls['iBi1Tot_nT'] = r'$\mathrm{Im}\{A_e\}$'
                        except Exception:
                            pass

                    # Monkey-patch Figure.savefig to unconditionally pin axis limits
                    # to the data range.  PlotExploreOgram's set_xlim([np.min(x), np.max(x)])
                    # silently no-ops when x contains NaN (derived axes like sigmaMean_Sm),
                    # and pcolormesh + sharex=True can then re-trigger autoscale.
                    # Restored in finally.
                    from matplotlib.figure import Figure as _MplFigure
                    _orig_savefig = _MplFigure.savefig

                    _xn = Exploration.xName
                    _yn = Exploration.yName

                    # Capture log-scale checkbox state in closure.
                    _force_x_log = bool(x_log)
                    _force_y_log = bool(y_log)

                    # Capture salinity-axis and derived-label state for secondary-axis logic.
                    _show_sal_axis   = st.session_state.get('exploreogram_show_salinity_axis', False)
                    _x_derived_label = x_derived_label   # 'sigmaMean_Sm' or None
                    _y_derived_label = y_derived_label   # 'sigmaMean_Sm' or None

                    # Substitute Exploration.xData/yData with derived arrays so that
                    # PlotExploreOgram (called below) positions pcolormesh cells in
                    # conductivity space instead of salinity space.  Without this,
                    # nonlinear σ(w) produces unequal column widths (jagged cells).
                    # Restored unconditionally in the finally block below.
                    _xdata_substituted = False
                    _orig_xData = None
                    _orig_xName = None
                    _ydata_substituted = False
                    _orig_yData = None
                    _orig_yName = None
                    _SIG_SUBST = 'sigmaMean_Sm'
                    if (_x_derived_label == _SIG_SUBST and
                            hasattr(Exploration, 'base') and Exploration.base is not None):
                        _derived_x_arr = getattr(Exploration.base, _SIG_SUBST, None)
                        if (_derived_x_arr is not None and
                                np.shape(_derived_x_arr) == np.shape(Exploration.xData)):
                            _orig_xData = Exploration.xData.copy()
                            _orig_xName = Exploration.xName
                            # Use column-mean conductivity (averaged over D) tiled to (nx, ny).
                            # The raw 2D sigmaMean_Sm varies with BOTH salinity and D; using it
                            # directly produces a non-rectilinear pcolormesh grid with staircase
                            # cell boundaries.  The mean gives a single representative σ per
                            # salinity column, consistent with the secondary-axis tick mapping.
                            _sig2d_raw = np.asarray(_derived_x_arr, dtype=float)
                            _sig1d_rep = np.nanmean(_sig2d_raw, axis=1)            # (nx,)
                            _sig2d_rect = np.tile(_sig1d_rep[:, np.newaxis],
                                                  (1, _sig2d_raw.shape[1]))        # (nx, ny)
                            Exploration.xData = _sig2d_rect
                            Exploration.xName = _SIG_SUBST
                            _xdata_substituted = True
                    if (_y_derived_label == _SIG_SUBST and
                            hasattr(Exploration, 'base') and Exploration.base is not None):
                        _derived_y_arr = getattr(Exploration.base, _SIG_SUBST, None)
                        if (_derived_y_arr is not None and
                                np.shape(_derived_y_arr) == np.shape(Exploration.yData)):
                            _orig_yData = Exploration.yData.copy()
                            _orig_yName = Exploration.yName
                            Exploration.yData = np.asarray(_derived_y_arr, dtype=float)
                            Exploration.yName = _SIG_SUBST
                            _ydata_substituted = True

                    # Pre-compute data-driven axis limits from PP-computed arrays.
                    # When a derived axis is active (e.g. sigmaMean_Sm replacing wOcean_ppt),
                    # use the DERIVED values for xlim so the axis spans conductivity space,
                    # not salinity space.
                    _xdata_arr = None   # valid x values for xlim
                    _ydata_arr = None   # valid y values for ylim
                    _SIG = 'sigmaMean_Sm'
                    if hasattr(Exploration, 'base') and Exploration.base is not None:
                        _valid_mask = None
                        if hasattr(Exploration.base, 'VALID') and Exploration.base.VALID is not None:
                            _vm = np.asarray(Exploration.base.VALID, dtype=bool).ravel()
                            _valid_mask = _vm

                        # For x: prefer derived quantity when label is set
                        _raw_x_name = _x_derived_label if _x_derived_label == _SIG else _xn
                        _raw_x = getattr(Exploration.base, _raw_x_name, None)
                        _raw_y_name = _y_derived_label if _y_derived_label == _SIG else _yn
                        _raw_y = getattr(Exploration.base, _raw_y_name, None)
                        if _raw_x is not None:
                            _arr = np.asarray(_raw_x, dtype=float).ravel()
                            _mask = np.isfinite(_arr)
                            if _valid_mask is not None and _valid_mask.size == _arr.size:
                                _mask = _mask & _valid_mask
                            if _force_x_log:
                                _mask = _mask & (_arr > 0)
                            _valid = _arr[_mask]
                            if _valid.size >= 2:
                                _xdata_arr = _valid
                        if _raw_y is not None:
                            _arr = np.asarray(_raw_y, dtype=float).ravel()
                            _mask = np.isfinite(_arr)
                            if _valid_mask is not None and _valid_mask.size == _arr.size:
                                _mask = _mask & _valid_mask
                            if _force_y_log:
                                _mask = _mask & (_arr > 0)
                            _valid = _arr[_mask]
                            if _valid.size >= 2:
                                _ydata_arr = _valid

                    def _patched_savefig(self_fig, *args, **kwargs):
                        try:
                            _data_axes = [a for a in self_fig.axes
                                          if getattr(a, '_label', '') != '<colorbar>']
                            for ax in _data_axes:
                                # Enforce x/y scale from checkbox state.
                                try:
                                    target_xscale = 'log' if _force_x_log else 'linear'
                                    target_yscale = 'log' if _force_y_log else 'linear'
                                    if ax.get_xscale() != target_xscale:
                                        ax.set_xscale(target_xscale)
                                    if ax.get_yscale() != target_yscale:
                                        ax.set_yscale(target_yscale)
                                except Exception:
                                    pass

                                # Unconditionally pin axis limits to the data range.
                                # PlotExploreOgram's set_xlim([np.min(x), np.max(x)]) silently
                                # no-ops when x has NaN entries (derived axes like sigmaMean_Sm),
                                # and pcolormesh + sharex=True can override it even when it
                                # succeeds.  Data-driven values from _xdata_arr are always correct.
                                try:
                                    if _xdata_arr is not None:
                                        xmin_d = float(_xdata_arr.min())
                                        xmax_d = float(_xdata_arr.max())
                                        if np.isfinite(xmin_d) and np.isfinite(xmax_d) and xmax_d > xmin_d:
                                            ax.set_autoscalex_on(False)
                                            ax.set_xlim([xmin_d, xmax_d])
                                except Exception:
                                    pass
                                try:
                                    if _ydata_arr is not None:
                                        ymin_d = float(_ydata_arr.min())
                                        ymax_d = float(_ydata_arr.max())
                                        if np.isfinite(ymin_d) and np.isfinite(ymax_d) and ymax_d > ymin_d:
                                            ax.set_autoscaley_on(False)
                                            ax.set_ylim([ymin_d, ymax_d])
                                except Exception:
                                    pass

                                # Derived-axis re-mapping: place axis in DERIVED coordinate
                                # space (conductivity S/m) so the range and tick positions
                                # are correct.  Then add a secondary axis for the driver
                                # quantity (salinity ppt).
                                if _show_sal_axis and hasattr(Exploration, 'base') and Exploration.base is not None:
                                    try:
                                        from matplotlib.ticker import MaxNLocator, LogLocator
                                        _base = Exploration.base
                                        if _x_derived_label == _SIG:
                                            _sig2d = getattr(_base, 'sigmaMean_Sm', None)
                                            _w2d   = getattr(_base, 'wOcean_ppt',   None)
                                            if _sig2d is not None and _w2d is not None:
                                                # 1-D σ(w): mean across y (axis=1) → length nx
                                                _sig1d = np.nanmean(np.asarray(_sig2d, dtype=float), axis=1)
                                                _w1d   = np.nanmean(np.asarray(_w2d,   dtype=float), axis=1)
                                                _sort  = np.argsort(_w1d)
                                                _ss = _sig1d[_sort]; _ws = _w1d[_sort]  # w ascending → σ monotone
                                                _sig_min = float(np.nanmin(_ss))
                                                _sig_max = float(np.nanmax(_ss))
                                                if np.isfinite(_sig_min) and np.isfinite(_sig_max) and _sig_max > _sig_min:
                                                    # 1. Set xlim to conductivity range
                                                    ax.set_autoscalex_on(False)
                                                    ax.set_xlim([_sig_min, _sig_max])
                                                    # 2. Place ticks in conductivity space
                                                    if _force_x_log:
                                                        _loc = LogLocator(base=10, numticks=6)
                                                        _sig_ticks = np.array([t for t in _loc.tick_values(_sig_min, _sig_max)
                                                                                if _sig_min <= t <= _sig_max])
                                                    else:
                                                        _loc = MaxNLocator(nbins=6, integer=False)
                                                        _sig_ticks = np.array([t for t in _loc.tick_values(_sig_min, _sig_max)
                                                                                if _sig_min <= t <= _sig_max])
                                                    if len(_sig_ticks) == 0:
                                                        _sig_ticks = np.linspace(_sig_min, _sig_max, 5)
                                                    ax.set_xticks(_sig_ticks)
                                                    ax.set_xticklabels([f'{t:.3g}' for t in _sig_ticks])
                                                    ax.set_xlabel('Mean Conductivity (S/m)')
                                                    # 3. Secondary axis: use sparse "nice" salinity values mapped
                                                    # to conductivity positions.  Driving from salinity → σ avoids
                                                    # crowded low-end labels that appear when conductivity decades
                                                    # (0.0001, 0.001, 0.01 …) map to very small, overlapping strings.
                                                    if _show_sal_axis and not getattr(ax, '_pp_secax_x_attached', False):
                                                        _w_pos = _ws[_ws > 0]
                                                        _w_min_s = float(np.nanmin(_w_pos)) if _w_pos.size > 0 else 1e-3
                                                        _w_max_s = float(np.nanmax(_ws))
                                                        _candidates = [0.001, 0.003, 0.01, 0.03, 0.1,
                                                                       0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
                                                        _w_nice = np.array([v for v in _candidates
                                                                            if _w_min_s <= v <= _w_max_s])
                                                        if len(_w_nice) < 2:
                                                            _w_nice = np.logspace(
                                                                np.log10(max(_w_min_s, 1e-4)),
                                                                np.log10(_w_max_s), 4)
                                                        _sig_at_sal = np.interp(_w_nice, _ws, _ss)
                                                        _in_rng = (_sig_at_sal >= _sig_min) & (_sig_at_sal <= _sig_max)
                                                        _w_nice = _w_nice[_in_rng]
                                                        _sig_at_sal = _sig_at_sal[_in_rng]
                                                        if len(_w_nice) >= 2:
                                                            _secax = ax.secondary_xaxis('top')
                                                            _secax.set_xticks(_sig_at_sal)
                                                            _secax.set_xticklabels([f'{v:.3g}' for v in _w_nice])
                                                            _secax.set_xlabel('Salinity (ppt)')
                                                            ax._pp_secax_x_attached = True
                                        if _y_derived_label == _SIG:
                                            _sig2d = getattr(_base, 'sigmaMean_Sm', None)
                                            _w2d   = getattr(_base, 'wOcean_ppt',   None)
                                            if _sig2d is not None and _w2d is not None:
                                                # 1-D σ(w): mean across x (axis=0) → length ny
                                                _sig1d = np.nanmean(np.asarray(_sig2d, dtype=float), axis=0)
                                                _w1d   = np.nanmean(np.asarray(_w2d,   dtype=float), axis=0)
                                                _sort  = np.argsort(_sig1d)
                                                _ss = _sig1d[_sort]; _ws = _w1d[_sort]
                                                _sig_min = float(np.nanmin(_ss))
                                                _sig_max = float(np.nanmax(_ss))
                                                if np.isfinite(_sig_min) and np.isfinite(_sig_max) and _sig_max > _sig_min:
                                                    ax.set_autoscaley_on(False)
                                                    ax.set_ylim([_sig_min, _sig_max])
                                                    if _force_y_log:
                                                        _loc = LogLocator(base=10, numticks=6)
                                                        _sig_ticks = np.array([t for t in _loc.tick_values(_sig_min, _sig_max)
                                                                                if _sig_min <= t <= _sig_max])
                                                    else:
                                                        _loc = MaxNLocator(nbins=6, integer=False)
                                                        _sig_ticks = np.array([t for t in _loc.tick_values(_sig_min, _sig_max)
                                                                                if _sig_min <= t <= _sig_max])
                                                    if len(_sig_ticks) == 0:
                                                        _sig_ticks = np.linspace(_sig_min, _sig_max, 5)
                                                    ax.set_yticks(_sig_ticks)
                                                    ax.set_yticklabels([f'{t:.3g}' for t in _sig_ticks])
                                                    ax.set_ylabel('Mean Conductivity (S/m)')
                                                    if _show_sal_axis and not getattr(ax, '_pp_secax_y_attached', False):
                                                        _sal_at_sig = np.interp(_sig_ticks, _ss, _ws)
                                                        _secay = ax.secondary_yaxis('right')
                                                        _secay.set_yticks(_sig_ticks)
                                                        _secay.set_yticklabels([f'{v:.3g}' for v in _sal_at_sig])
                                                        _secay.set_ylabel('Salinity (ppt)')
                                                        ax._pp_secax_y_attached = True
                                    except Exception:
                                        pass
                        except Exception:
                            pass
                        return _orig_savefig(self_fig, *args, **kwargs)

                    _MplFigure.savefig = _patched_savefig

                    # Temporarily enable plots
                    old_skip = Params.SKIP_PLOTS
                    Params.SKIP_PLOTS = False

                    # Generate matplotlib plots
                    try:
                        with st.spinner("Generating publication-quality matplotlib plots..."):
                            GenerateExplorationPlots(ExplorationList, FigureFilesList, Params)
                            plt.close('all')
                    finally:
                        # Restore monkey-patched savefig
                        _MplFigure.savefig = _orig_savefig
                        # Restore Exploration.xData/yData if substituted for derived-axis rendering
                        if _xdata_substituted and _orig_xData is not None:
                            Exploration.xData = _orig_xData
                            Exploration.xName = _orig_xName
                        if _ydata_substituted and _orig_yData is not None:
                            Exploration.yData = _orig_yData
                            Exploration.yName = _orig_yName
                        # Restore log-scale list
                        _FigLbl_modlog.axisLogScalesExplore = _orig_axisLogScales
                        # Restore amplitude/RI overrides
                        if _amp_was_overridden:
                            _ind.Amp = _amp_orig
                        if _ri_was_overridden:
                            _ind.rBi1Tot_nT = _r_orig
                            _ind.iBi1Tot_nT = _i_orig
                        # Restore labels
                        try:
                            if _orig_amp_label is not None:
                                _FigLbl_modlog.axisLabelsExplore['Amp_nT'] = _orig_amp_label
                            if _orig_r_label is not None:
                                _FigLbl_modlog.axisLabelsExplore['rBi1Tot_nT'] = _orig_r_label
                            if _orig_i_label is not None:
                                _FigLbl_modlog.axisLabelsExplore['iBi1Tot_nT'] = _orig_i_label
                        except Exception:
                            pass

                    Params.SKIP_PLOTS = old_skip

                    # Display generated plots
                    explore_path = FigureFiles.explore
                    if isinstance(explore_path, list):
                        # Multi-subplot generates one file
                        explore_path = FigureFiles.exploreMultiSubplot if hasattr(FigureFiles, 'exploreMultiSubplot') else explore_path[0]
                    if os.path.exists(explore_path):
                        st.success("Matplotlib plot generated!")
                        images = convert_from_path(explore_path)
                        for img in images:
                            st.image(img, use_container_width=True)
                    else:
                        # Try finding any generated PDF in output dir
                        import glob
                        pdfs = glob.glob(os.path.join(output_dir, f"{fig_basename}*.pdf"))
                        if pdfs:
                            for pdf_path in pdfs:
                                images = convert_from_path(pdf_path)
                                for img in images:
                                    st.image(img, use_container_width=True, caption=os.path.basename(pdf_path))
                        else:
                            st.warning(f"Plot file not found at: {explore_path}")
                            st.info("Try enabling 'Interactive Plots' checkbox for Plotly-based plots instead.")

                except Exception as e:
                    st.error(f"Error generating matplotlib plot: {e}")
                    st.info("💡 Try enabling 'Interactive Plots' instead")
                    import traceback
                    with st.expander("Show error details"):
                        st.code(traceback.format_exc())

            # Statistics (get zData for stats regardless of plot type)
            if Exploration.base.VALID is not None:
                induction_var_map = {
                    'Amp_nT': 'Amp', 'phase_deg': 'Phase',
                    'Bix_nT': 'Bix_nT', 'Biy_nT': 'Biy_nT', 'Biz_nT': 'Biz_nT',
                    'Bi1x_nT': 'Bi1x_nT', 'Bi1y_nT': 'Bi1y_nT', 'Bi1z_nT': 'Bi1z_nT',
                    'Bi1Tot_nT': 'Bi1Tot_nT',
                    'rBi1x_nT': 'rBi1x_nT', 'rBi1y_nT': 'rBi1y_nT', 'rBi1z_nT': 'rBi1z_nT',
                    'rBi1Tot_nT': 'rBi1Tot_nT',
                    'iBi1x_nT': 'iBi1x_nT', 'iBi1y_nT': 'iBi1y_nT', 'iBi1z_nT': 'iBi1z_nT',
                    'iBi1Tot_nT': 'iBi1Tot_nT',
                }
                if zName in induction_var_map and hasattr(Exploration, 'induction') and Exploration.induction is not None:
                    attr = induction_var_map[zName]
                    raw = getattr(Exploration.induction, attr, None)
                    if raw is not None and len(raw.shape) == 3:
                        zData_stats = np.real(raw[0, :, :])  # First peak for stats
                    else:
                        zData_stats = np.zeros_like(Exploration.base.VALID, dtype=float)
                elif hasattr(Exploration.base, zName):
                    zData_stats = getattr(Exploration.base, zName)
                else:
                    zData_stats = np.zeros_like(Exploration.base.VALID, dtype=float)

                st.markdown("#### Statistics")
                col_stat1, col_stat2, col_stat3, col_stat4 = st.columns(4)
                n_total = Exploration.base.VALID.size
                with col_stat1:
                    st.metric("Valid Models", f"{np.sum(Exploration.base.VALID)}/{n_total}")
                with col_stat2:
                    st.metric("Min", f"{np.nanmin(zData_stats):.4g}")
                with col_stat3:
                    st.metric("Max", f"{np.nanmax(zData_stats):.4g}")
                with col_stat4:
                    st.metric("Mean", f"{np.nanmean(zData_stats):.4g}")
            else:
                st.warning("VALID array not available — statistics unavailable.")

            # Download buttons
            st.markdown("---")
            st.markdown("#### Save Options")

            if use_interactive_plots:
                col_dl1, col_dl2, col_dl3 = st.columns(3)
                with col_dl1:
                    # Save interactive HTML
                    if st.button("💾 Save Interactive (HTML)"):
                        output_dir = os.path.join(parent_directory, 'output', 'exploreograms')
                        os.makedirs(output_dir, exist_ok=True)
                        filename = f"{results['Planet'].name}_explore_{results['xName']}_vs_{results['yName']}_interactive.html"
                        filepath = os.path.join(output_dir, filename)
                        if 'fig' in locals():
                            fig.write_html(filepath)
                            st.success(f"Saved to `{filepath}`")
                        else:
                            st.error("No figure available to save")
                with col_dl2:
                    # Save static image
                    if st.button("💾 Save Image (PNG)"):
                        output_dir = os.path.join(parent_directory, 'output', 'exploreograms')
                        os.makedirs(output_dir, exist_ok=True)
                        filename = f"{results['Planet'].name}_explore_{results['xName']}_vs_{results['yName']}.png"
                        filepath = os.path.join(output_dir, filename)
                        if 'fig' in locals():
                            try:
                                fig.write_image(filepath, width=1200, height=900)
                                st.success(f"Saved to `{filepath}`")
                            except Exception as e:
                                st.error(f"Error saving PNG: {e}")
                                st.info("Install kaleido for image export: pip install kaleido")
                        else:
                            st.error("No figure available to save")
                with col_dl3:
                    # Save results
                    if st.button("💾 Save Data (PKL)"):
                        output_dir = os.path.join(parent_directory, 'output', 'exploreograms')
                        os.makedirs(output_dir, exist_ok=True)
                        filename = f"{results['Planet'].name}_explore_{results['xName']}_vs_{results['yName']}.pkl"
                        filepath = os.path.join(output_dir, filename)
                        with open(filepath, 'wb') as f:
                            pickle.dump(results, f)
                        st.success(f"Saved to `{filepath}`")
            else:
                col_dl1, col_dl2 = st.columns(2)
                with col_dl1:
                    # PDF already saved by GenerateExplorationPlots
                    # For induction multi-subplot, use exploreMultiSubplot path
                    pdf_path = None
                    is_induction_z = zName in induction_z_vars
                    if hasattr(FigureFiles, 'exploreMultiSubplot') and is_induction_z:
                        candidate = FigureFiles.exploreMultiSubplot
                        if os.path.exists(candidate):
                            pdf_path = candidate
                    if pdf_path is None and hasattr(FigureFiles, 'explore'):
                        candidate = FigureFiles.explore
                        if isinstance(candidate, list):
                            for c in candidate:
                                if os.path.exists(c):
                                    pdf_path = c
                                    break
                        elif os.path.exists(candidate):
                            pdf_path = candidate
                    if pdf_path is None:
                        # Fallback: search output dir for any matching PDF
                        import glob as glob_mod
                        output_dir = os.path.join(parent_directory, 'output', 'exploreograms')
                        fig_basename = f"{results['Planet'].name}_explore_{results['xName']}_vs_{results['yName']}"
                        matches = glob_mod.glob(os.path.join(output_dir, f"{fig_basename}*.pdf"))
                        if matches:
                            pdf_path = matches[0]

                    if pdf_path is not None:
                        st.info(f"PDF saved at:\n`{pdf_path}`")
                    else:
                        st.warning("PDF not found — check output/exploreograms/ directory")

                with col_dl2:
                    # Save results
                    if st.button("💾 Save Data (PKL)"):
                        output_dir = os.path.join(parent_directory, 'output', 'exploreograms')
                        os.makedirs(output_dir, exist_ok=True)
                        filename = f"{results['Planet'].name}_explore_{results['xName']}_vs_{results['yName']}.pkl"
                        filepath = os.path.join(output_dir, filename)
                        with open(filepath, 'wb') as f:
                            pickle.dump(results, f)
                        st.success(f"Saved to `{filepath}`")

            # Export editable Python script for matplotlib plots
            if st.button("Export Editable Python Script"):
                output_dir = os.path.join(parent_directory, 'output', 'exploreograms')
                os.makedirs(output_dir, exist_ok=True)
                pkl_filename = f"{results['Planet'].name}_explore_{results['xName']}_vs_{results['yName']}.pkl"
                pkl_filepath = os.path.join(output_dir, pkl_filename)
                # Ensure Exploration is stored directly in results for script access
                if 'Exploration' not in results and hasattr(results.get('Planet'), 'Exploration'):
                    results['Exploration'] = results['Planet'].Exploration
                elif 'Exploration' not in results:
                    results['Exploration'] = Exploration
                # Save PKL (always overwrite to capture latest Exploration)
                with open(pkl_filepath, 'wb') as f:
                    pickle.dump(results, f)

                py_filename = f"plot_{results['Planet'].name}_{results['xName']}_vs_{results['yName']}.py"
                py_filepath = os.path.join(output_dir, py_filename)

                # Build z-variable info for the script
                is_induction_z_export = zName in induction_z_vars
                if is_induction_z_export:
                    display_mode_export = st.session_state.get('induct_display_mode', 'real_imaginary')
                    if display_mode_export == 'amplitude_phase':
                        zname_export = "['Amp_nT', 'phase_deg']"
                    elif zName in ('Amp_nT', 'phase_deg'):
                        zname_export = "['Bi1Tot_nT']"
                    else:
                        zname_export = f"['{zName}']"
                    combo_export = 'True'
                else:
                    zname_export = f"'{zName}'"
                    combo_export = 'False'

                use_contours_export = st.session_state.get('induct_use_contours', True)

                script_content = f'''#!/usr/bin/env python
"""
Editable matplotlib exploreogram plot script.
Generated by PlanetProfileApp — modify as needed for publication.

Data file: {pkl_filepath}
"""
import sys, os
import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# --- Add PlanetProfile to path and set working directory ---
# GetConfig.py looks for SPICE kernels relative to CWD, so we must
# change to the repo root before importing anything from PlanetProfile.
sys.path.insert(0, r'{parent_directory}')
os.chdir(r'{parent_directory}')

# --- Load data ---
PKL_PATH = r'{pkl_filepath}'
with open(PKL_PATH, 'rb') as f:
    results = pickle.load(f)

from PlanetProfile.GetConfig import FigLbl, FigMisc, Color, Style
from PlanetProfile.Plotting.ExplorationPlots import GenerateExplorationPlots
from PlanetProfile.Utilities.defineStructs import FigureFilesSubstruct, ParamsStruct, ExploreParamsStruct

# --- Illustrator-compatible fonts (AFTER GetConfig import, which sets STIX) ---
mpl.rcParams['text.usetex'] = False
mpl.rcParams['mathtext.fontset'] = 'dejavusans'
mpl.rcParams['pdf.fonttype'] = 42      # TrueType — Illustrator can edit
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']

# --- Reconstruct Params and Exploration ---
Exploration = results.get('Exploration') or results['Planet'].Exploration
if Exploration is None:
    raise RuntimeError("Exploration data not found in PKL. Re-export from PlanetProfileApp.")
Params = ParamsStruct()
Params.Explore = ExploreParamsStruct()
Params.Explore.xName = Exploration.xName
Params.Explore.yName = Exploration.yName
Params.Explore.zName = {zname_export}
Params.Explore.nx = Exploration.nx
Params.Explore.ny = Exploration.ny
Params.Explore.contourName = None
Params.Explore.DRAW_CONTOURS = {use_contours_export}
Params.PLOT_COMBO_EXPLORATIONS = {combo_export}
Params.TITLES = True
Params.SKIP_PLOTS = False
Params.PLOT_Zb_D = False
Params.PLOT_D_SIGMA = False
Params.PLOT_LOVE_COMPARISON = False

# --- Output paths (edit as needed) ---
OUTPUT_DIR = os.path.dirname(PKL_PATH)
FIG_BASENAME = '{results["Planet"].name}_explore_{results["xName"]}_vs_{results["yName"]}'
FigureFiles = FigureFilesSubstruct(
    figPath=OUTPUT_DIR,
    figBase=FIG_BASENAME,
    xtn='.pdf',
    exploreAppend=Params.Explore.zName
)
Params.FigureFiles = FigureFiles

'''
                # Add induction params if needed
                if is_induction_z_export:
                    selected_exc = st.session_state.get('selected_excitations', [])
                    exc_dict_str = ', '.join(f"'{e}': True" for e in selected_exc)
                    script_content += f'''# --- Induction configuration ---
from PlanetProfile.Utilities.defineStructs import InductOgramParamsStruct
Params.Induct = InductOgramParamsStruct(inductOtype=None, cLevels=None, dftC=None, cfmt=None)
Params.Induct.excSelectionPlot = {{{exc_dict_str}}}

'''

                script_content += f'''# --- Generate plots ---
ExplorationList = [Exploration]
FigureFilesList = [FigureFiles]

GenerateExplorationPlots(ExplorationList, FigureFilesList, Params)
plt.show()

print(f"Plots saved to {{OUTPUT_DIR}}")
'''

                with open(py_filepath, 'w') as f:
                    f.write(script_content)
                st.success(f"Python script saved to:\n`{py_filepath}`")
                st.info("Run it with: `python " + py_filename + "` from the output directory, or open and edit in your IDE.")

        else:
            st.error("No exploration results found in Planet object.")

    elif not st.session_state.explore_running:
        st.info("👈 Configure parameters and click 'Run Exploreogram' to start")

# Tips
with st.expander("💡 Tips for Using Exploreograms"):
    st.markdown("""
    - **Plot types**:
        - **Interactive (default)**: Plotly plots with zoom, pan, hover data. Best for exploration.
        - **Static (matplotlib)**: Publication-quality plots matching CLI output. Best for papers.
    - **Start small**: Try a 5×5 grid first to verify parameters before running larger grids
    - **Parameter types**:
        - **Hydro**: Affect ocean/ice shell (salinity, temperature, porosity)
        - **Inner**: Affect mantle/core (density, heating, core size)
        - **Ionos**: Affect ionosphere conductivity (for magnetic induction)
    - **Computational cost**: A 10×10 grid runs 100 models. Allow ~5-10 seconds per model.
    - **Mixed types**: Combining hydro+inner parameters requires more computation than hydro+hydro or inner+inner
    - **Results**: Use the color variable to visualize how different properties vary across parameter space
    - **Saving**: Interactive plots save as HTML (interactive) or PNG. Static plots save as PDF.
    """)

st.markdown("---")
st.caption("Exploreogram functionality powered by PlanetProfile's parallel computation engine")
