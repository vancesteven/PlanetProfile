import streamlit as st
import os
import sys
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pickle

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

ALL_PARAMS = {**HYDRO_PARAMS, **INNER_PARAMS, **IONOS_PARAMS}

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
        'Tmean_K': 'Mean Ocean Temperature (K)',
        'Pseafloor_MPa': 'Seafloor Pressure (MPa)',
        'sigmaOcean_Sm': 'Ocean Conductivity (S/m)',
        'Amp_nT': 'Induced Magnetic Field (nT)',
        'k2': 'k₂ Love Number',
        'Rcore_km': 'Core Radius (km)',
        'Rmantle_km': 'Mantle Radius (km)',
    }

    z_param = st.selectbox(
        "Select color variable",
        options=list(Z_VARIABLES.keys()),
        index=0,
        format_func=lambda x: Z_VARIABLES[x],
        key='z_param_select'
    )
    st.session_state.explore_zName = z_param

    # Show induction frequency selection if magnetic variable is chosen
    induction_vars = ['Amp_nT']  # Magnetic field amplitude - add more as needed
    if z_param in induction_vars:
        st.markdown("#### Magnetic Induction Configuration")
        st.info("🧲 Magnetic induction calculations enabled for selected output variable")

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

    # Run button
    run_button = st.button("🚀 Run Exploreogram", type="primary", disabled=st.session_state.explore_running)

    if run_button and not st.session_state.explore_running:
        st.session_state.explore_running = True
        st.session_state.explore_start_time = None
        st.session_state.explore_error = None
        st.session_state.explore_error_traceback = None
        st.rerun()

with col2:
    st.subheader("📊 Exploration Results")

    if st.session_state.explore_running:
        # Estimate time
        est_time_per_model = 8  # seconds (conservative estimate)
        est_total_minutes = (total_runs * est_time_per_model) / 60

        # Create prominent information box
        st.info(f"""
        ### 🔬 Exploreogram Running: {total_runs} models ({nx}×{ny} grid)

        **⏱️ Estimated time: ~{est_total_minutes:.1f} minutes** ({est_time_per_model} sec/model)

        **📊 Progress Monitoring:**
        - This window shows overall status
        - **Look at your TERMINAL WINDOW** (where you ran `streamlit run`) for detailed model-by-model progress
        - The terminal shows: "Running model 1/{total_runs}", "Running model 2/{total_runs}", etc.
        - The GUI will automatically update when complete

        **💡 The GUI may appear frozen** - this is normal! Models are running in the background.
        """)

        # Use st.status() for long-running operation
        with st.status(f"Computing {total_runs} models...", expanded=True) as status:
            st.write("⚙️ Using parallel processing on available CPU cores")
            st.write("")

            try:
                import time
                start_time = time.time()

                # Show configuration
                st.write("**Configuration:**")
                config_info = st.container()
                with config_info:
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(f"🌍 Planet: {st.session_state['Planet'].name if 'Planet' in st.session_state else 'Unknown'}")
                        st.write(f"📊 X: {ALL_PARAMS[st.session_state.explore_xName]['label']}")
                        st.write(f"   Range: [{x_min}, {x_max}]")
                    with col2:
                        st.write(f"📈 Y: {ALL_PARAMS[st.session_state.explore_yName]['label']}")
                        st.write(f"   Range: [{y_min}, {y_max}]")
                        st.write(f"🎨 Color: {Z_VARIABLES.get(st.session_state.explore_zName, st.session_state.explore_zName)}")

                st.write("")
                st.write("**Progress:**")

                # Import PlanetProfile modules
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
                Params.Explore.xName = st.session_state.explore_xName
                Params.Explore.yName = st.session_state.explore_yName
                Params.Explore.zName = st.session_state.explore_zName
                Params.Explore.xRange = [x_min, x_max]
                Params.Explore.yRange = [y_min, y_max]

                # Configure induction calculations based on z-variable selection
                # Check if we need induction for the selected output variable
                induction_vars = ['Amp_nT', 'Bix_nT', 'Biy_nT', 'Biz_nT', 'phase_deg']  # Add more as needed
                needs_induction = st.session_state.explore_zName in induction_vars

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

                    st.write(f"   ✓ Magnetic induction enabled: {len(selected_exc)} excitation(s) selected")
                    st.write(f"     Frequencies: {', '.join(selected_exc)}")

                    # Preload excitation moments ONCE before exploreogram starts
                    # (During exploreogram, GetBexc is skipped in SetupInduction, so we need to load it here)
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
                    st.write("   ℹ️ Skipping magnetic induction (not needed for selected output)")
                Params.Explore.nx = int(nx)
                Params.Explore.ny = int(ny)

                # GUI-specific settings to speed up runs
                Params.SKIP_PLOTS = True
                Params.NO_SAVEFILE = True  # Don't save individual runs
                Params.SKIP_INNER = False  # Keep inner calculations
                Params.ALLOW_BROKEN_MODELS = True  # Continue even if some models fail

                st.write(f"✓ Configuration complete")
                st.write("")

                # Generate exploration lists
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

                st.write(f"✓ Grid created: {len(xList)}×{len(yList)} = {total_runs} models")
                st.write("")

                # Run parallel exploration
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
                    st.write(f"  - Grid shape: {np.shape(PlanetGrid)}")
                    gridShape = np.shape(PlanetGrid)
                    st.write(f"  - Flattening grid...")
                    PlanetList = np.reshape(PlanetGrid, -1)
                    st.write(f"  - Calling GetLayerMeans with {len(PlanetList)} planets...")
                    PlanetList, Params = GetLayerMeans(PlanetList, Params)
                    st.write(f"  - Reshaping back to grid...")
                    PlanetGrid = np.reshape(PlanetList, gridShape)
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

                    st.write(f"  - Calling ExtractResults...")

                    Exploration = ExtractResults(Exploration, PlanetGrid, Params)

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
                    st.write("  - Attaching Exploration to Planet...")
                    Planet.Exploration = Exploration
                    Params.EXPLOREOGRAM_IN_PROGRESS = False

                    # Check for valid models
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

                                if not np.any(Exploration.base.VALID):
                                    st.warning('⚠️ No valid models found for the given input settings.')
                            else:
                                st.error("❌ Exploration.base.VALID is None!")
                        else:
                            st.error("❌ Exploration.base does not have VALID attribute!")
                    else:
                        st.error("❌ Exploration.base is None!")

                    st.write("  - Storing to session state...")
                    st.session_state.explore_results = {
                        'Planet': Planet,
                        'Params': Params,
                        'xName': st.session_state.explore_xName,
                        'yName': st.session_state.explore_yName,
                        'zName': st.session_state.explore_zName,
                        'xData': Exploration.xData,
                        'yData': Exploration.yData,
                    }

                    st.write("  - Calculating statistics...")
                    elapsed_time = time.time() - start_time
                    st.write(f"    - Elapsed: {elapsed_time:.1f} sec")

                    # Calculate valid count with safety check
                    if Exploration.base is not None and hasattr(Exploration.base, 'VALID') and Exploration.base.VALID is not None:
                        st.write(f"    - Summing VALID array...")
                        valid_count = np.sum(Exploration.base.VALID)
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

                st.write("")
                st.success(f"**✅ Exploreogram Complete!**")
                col_time1, col_time2 = st.columns(2)
                with col_time1:
                    st.metric("Total Time", f"{elapsed_time/60:.1f} min")
                with col_time2:
                    st.metric("Time/Model", f"{elapsed_time/total_runs:.1f} sec")

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
                st.error(f"❌ Error extracting plot data: {e}")
                import traceback
                st.code(traceback.format_exc())
                st.stop()

            # --- PLOTTING ---
            if use_interactive_plots:
                # Use enhanced Plotly plotting
                try:
                    from Utilities.exploreogram_plotly import create_exploreogram_plotly

                    # Get z-variable data
                    # Check if this is a magnetic induction variable
                    induction_vars = ['Amp_nT', 'Bix_nT', 'Biy_nT', 'Biz_nT', 'phase_deg',
                                      'Bi1x_nT', 'Bi1y_nT', 'Bi1z_nT', 'Bi1Tot_nT',
                                      'rBi1x_nT', 'rBi1y_nT', 'rBi1z_nT', 'rBi1Tot_nT',
                                      'iBi1x_nT', 'iBi1y_nT', 'iBi1z_nT', 'iBi1Tot_nT']

                    if zName in induction_vars:
                        # Extract from induction substructure
                        if hasattr(Exploration, 'induction') and Exploration.induction is not None:
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

                            if hasattr(Exploration.induction, induction_attr_name):
                                induction_data = getattr(Exploration.induction, induction_attr_name)

                                # Induction data is 3D: (nPeaks, ny, nx)
                                # For now, use the first frequency peak
                                if induction_data is not None and len(induction_data.shape) == 3:
                                    zData = induction_data[0, :, :]  # First frequency peak

                                    # Show which frequency we're plotting
                                    if hasattr(Exploration.induction, 'calcedExc') and Exploration.induction.calcedExc:
                                        freq_names = list(Exploration.induction.calcedExc)
                                        st.info(f"Showing {zName} for frequency: {freq_names[0] if freq_names else 'unknown'}")
                                else:
                                    st.error(f"Induction data '{induction_attr_name}' has unexpected shape or is None")
                                    zData = np.zeros((len(yData), len(xData)))
                            else:
                                st.error(f"Could not find '{induction_attr_name}' in induction results. Induction may not have been calculated.")
                                st.info("Available induction variables: " + str([attr for attr in dir(Exploration.induction) if not attr.startswith('_')]))
                                zData = np.zeros((len(yData), len(xData)))
                        else:
                            st.error(f"Magnetic induction data not found. Make sure induction calculations were enabled.")
                            zData = np.zeros((len(yData), len(xData)))
                    elif hasattr(Exploration.base, zName):
                        # Regular variable from base structure
                        zData = getattr(Exploration.base, zName)
                    else:
                        st.warning(f"Could not find '{zName}' in exploration results. Available variables:")
                        st.write([attr for attr in dir(Exploration.base) if not attr.startswith('_')])
                        zData = np.zeros((len(yData), len(xData)))

                    # Try to load FigLbl for proper styling
                    try:
                        from PlanetProfile.GetConfig import FigLbl
                        FigLbl.SetExploration(results['Planet'].name, Exploration.xName,
                                            Exploration.yName, Exploration.zName)
                    except:
                        FigLbl = None

                    # Create enhanced Plotly plot
                    fig = create_exploreogram_plotly(
                        Exploration, Params, FigLbl=FigLbl,
                        smoothing=False, smooth_factor=2
                    )

                    st.plotly_chart(fig, use_container_width=True)

                except Exception as e:
                    st.error(f"Error creating interactive plot: {e}")
                    st.info("Falling back to basic Plotly plot...")

                    # Fallback to basic plot
                    if hasattr(Exploration.base, zName):
                        zData = getattr(Exploration.base, zName)
                    else:
                        zData = np.zeros((len(yData), len(xData)))

                    # Get axis labels safely
                    xLabel = ALL_PARAMS.get(results["xName"], {}).get("label", results["xName"])
                    yLabel = ALL_PARAMS.get(results["yName"], {}).get("label", results["yName"])
                    zLabel = Z_VARIABLES.get(zName, zName)

                    fig = go.Figure(data=go.Heatmap(
                        x=xData,
                        y=yData,
                        z=zData,
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

                    st.plotly_chart(fig, use_container_width=True)

            else:
                # Use matplotlib plots via GenerateExplorationPlots
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
                    FigureFiles = FigureFilesSubstruct(
                        figPath=output_dir,
                        figBase=fig_basename,
                        xtn='pdf'
                    )

                    FigureFilesList = [FigureFiles]

                    # Temporarily enable plots
                    old_skip = Params.SKIP_PLOTS
                    Params.SKIP_PLOTS = False

                    # Generate matplotlib plots
                    with st.spinner("Generating publication-quality matplotlib plots..."):
                        GenerateExplorationPlots(ExplorationList, FigureFilesList, Params)
                        plt.close('all')  # Clean up

                    Params.SKIP_PLOTS = old_skip

                    # Display the generated PDF
                    if os.path.exists(FigureFiles.explore):
                        st.success("✅ Matplotlib plot generated!")
                        images = convert_from_path(FigureFiles.explore)
                        st.image(images[0], use_container_width=True, caption="Publication-quality matplotlib plot")
                    else:
                        st.warning("Plot file not found. Falling back to interactive plot.")
                        use_interactive_plots = True
                        st.rerun()

                except Exception as e:
                    st.error(f"Error generating matplotlib plot: {e}")
                    st.info("💡 Try enabling 'Interactive Plots' instead")
                    import traceback
                    with st.expander("Show error details"):
                        st.code(traceback.format_exc())

            # Statistics (get zData for stats regardless of plot type)
            if hasattr(Exploration.base, zName):
                zData_stats = getattr(Exploration.base, zName)
            else:
                zData_stats = np.zeros((len(yData), len(xData)))

            st.markdown("#### Statistics")
            col_stat1, col_stat2, col_stat3, col_stat4 = st.columns(4)
            with col_stat1:
                st.metric("Valid Models", f"{np.sum(Exploration.base.VALID)}/{len(xData)*len(yData)}")
            with col_stat2:
                st.metric("Min", f"{np.nanmin(zData_stats):.2f}")
            with col_stat3:
                st.metric("Max", f"{np.nanmax(zData_stats):.2f}")
            with col_stat4:
                st.metric("Mean", f"{np.nanmean(zData_stats):.2f}")

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
                            st.success(f"✅ Saved to {filename}")
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
                                st.success(f"✅ Saved to {filename}")
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
                        st.success(f"✅ Saved to {filename}")
            else:
                col_dl1, col_dl2 = st.columns(2)
                with col_dl1:
                    # PDF already saved by GenerateExplorationPlots
                    output_dir = os.path.join(parent_directory, 'output', 'exploreograms')
                    filename = f"{results['Planet'].name}_explore_{results['xName']}_vs_{results['yName']}.pdf"
                    filepath = os.path.join(output_dir, filename)
                    if os.path.exists(filepath):
                        st.info(f"📄 Matplotlib PDF saved at:\n{filename}")
                    else:
                        st.warning("PDF not found")

                with col_dl2:
                    # Save results
                    if st.button("💾 Save Data (PKL)"):
                        output_dir = os.path.join(parent_directory, 'output', 'exploreograms')
                        os.makedirs(output_dir, exist_ok=True)
                        filename = f"{results['Planet'].name}_explore_{results['xName']}_vs_{results['yName']}.pkl"
                        filepath = os.path.join(output_dir, filename)
                        with open(filepath, 'wb') as f:
                            pickle.dump(results, f)
                        st.success(f"✅ Saved to {filename}")
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
