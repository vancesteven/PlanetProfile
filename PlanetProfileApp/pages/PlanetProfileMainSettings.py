import streamlit as st
import os
import importlib
import sys
# reaktoro was imported here but never used; it is conda-forge-only and
# breaks pip-only deployments, so it is intentionally not imported.
from functools import partial

# --- Setting up Main Page ---
st.set_page_config(
    page_title="PlanetProfile Main",
    page_icon="🌙",
    layout="wide"
)

# Import utility modules
from Utilities.session_manager import show_session_manager_ui, show_recent_sessions
from Utilities.presets import get_presets_for_body, apply_preset, get_preset_summary
from Utilities.app_helpers import create_progress_indicator, show_tips_carousel, create_metric_card
from Utilities.help_system import create_help_button

st.title("PlanetProfile")
st.markdown("---")

# Show recent sessions and session management in sidebar
with st.sidebar:
    st.markdown("### 📁 Session Management")
    show_recent_sessions()
    with st.expander("💾 Save/Load Sessions"):
        show_session_manager_ui()

st.write("## Let's Start by Setting Up Your Planet!")

# Show helpful tips carousel
show_tips_carousel([
    "💡 Start with a planet template and customize settings on the following pages",
    "📊 Use presets for common scenarios (thin ice shell, thick ice shell, etc.)",
    "💾 Save your configuration to resume later",
    "❓ Click the help buttons (?) for detailed explanations"
])

st.write("Note on custom planets - start with an existing planet as a template, and then edit your settings on the following pages.")
st.markdown("---")

# --- File Path Management ---

# Get the path to the current script's directory
# /PlanetProfile/PlanetProfileApp/pages/PlanetProfileMainSettings.py
PlanetProfileMainSettings_directory = os.path.dirname(os.path.abspath(__file__))

# Get the app directory (/PlanetProfile/PlanetProfileApp)
app_directory = os.path.dirname(PlanetProfileMainSettings_directory)
if app_directory not in sys.path:
    sys.path.append(app_directory)

# Point to the root PlanetProfile directory (outermost one)
this_file = os.path.abspath(__file__)
PlanetProfileMainSettings_directory = os.path.dirname(this_file)
app_directory = os.path.dirname(PlanetProfileMainSettings_directory)
parent_directory = os.path.dirname(app_directory)

if parent_directory not in sys.path:
    sys.path.insert(0, parent_directory)
# Get the parent directory (/PlanetProfile)
#parent_directory  = os.path.dirname(app_directory)

# Add the parent directory to Python's search path.
#if parent_directory not in sys.path:
    #sys.path.append(parent_directory)


#--- Planet Selection - Custom or from dropdown ---
planet_list = ["-- Select a Planet --", "Ariel", "Callisto", "Dione", "Enceladus", "Europa", "Ganymede",
               "Iapetus", "Io", "Luna", "Mimas", "Miranda", "Oberon", "Pluto",
               "Rhea", "Tethys", "Titan", "Titania", "Triton", "Umbriel"]

# Initialzing all of the things to session state
if "planet_selectbox" not in st.session_state:
    st.session_state["ChosenPlanet"] = "-- Select a Planet --"

if "Planet" not in st.session_state:
    st.session_state["Planet"] = None

#if "run_custom_body" not in st.session_state:
    #st.session_state["run_custom_body"] = False


# ----- Removing for now - fully custom planet option. User will instead start from a planet template because that will make it more likely that they will have a successful model ------
# Option for user to select fully custom planet
#st.subheader("Run fully custom Planet?")
#st.write(
    #"Planet Profile has profiles of many moons ready for you to use. "
    #"If you want to create your own moon, check the box below. "
    #"Otherwise, choose from existing planetary bodies."
#)


# Checkbox manages if the user wants to run a custom planet or not
#run_custom_body = st.checkbox("Create fully custom Planet?", value=st.session_state["run_custom_body"], key="run_custom_body")

# Custom Planet path
#if run_custom_body:
    #from Utilities.CustomPlanetGenerator import generate_custom_pp_template
    #st.session_state["ChosenPlanet"] = "Custom"
    #st.session_state["Planet"] = None  # Custom settings handled elsewhere
    #st.success("Using Custom Planet. Configure settings in other tabs.")
    # Define folder and file paths (with capital 'Custom')
    #custom_dir = os.path.join(parent_directory, "PlanetProfile", "Default", "Custom")
    #custom_file = os.path.join(custom_dir, "PPCustom.py")

    # Create folder if it doesn't exist
    #if not os.path.exists(custom_dir):
        #os.makedirs(custom_dir)
        #st.info("Created 'Custom' folder for custom planet.")

    # Create PPCustom.py with default contents if it doesn't exist
    #if not os.path.isfile(custom_file):
        #generate_custom_pp_template(parent_directory)
        #st.info("Created PPCustom.py with default custom planet config.")


# This is how the app can find custom modules made by the user and the paths to them
def get_custom_module_paths(parent_dir, planet_name):
    """
    Returns a list of custom module file paths for a given planet,
    excluding the default module (PP{Planet}.py) that is a copy of the default.

    Example path scanned: PlanetProfile/Europa/PPEuropaLowTemp.py
    """
    custom_dir = os.path.join(parent_dir, planet_name)
    default_filename = f"PP{planet_name}.py"

    if not os.path.exists(custom_dir):
        return []

    all_files = os.listdir(custom_dir)

    # Get files like PPEuropaLowTemp.py but not PPEuropa.py
    custom_files = [
        f for f in all_files
        if f.startswith(f"PP{planet_name}") and f.endswith(".py") and f != default_filename
    ]

    return [os.path.join(custom_dir, f) for f in custom_files]

# Example usage

planet_name = "Europa"
get_custom_module_paths(parent_directory, planet_name)



st.subheader("Body Selection")

# Create two columns for body selection and presets
col1, col2 = st.columns([2, 1])

with col1:
    st.write("Please select a planetary body from the list of profiles below")
    selected_planet = st.selectbox(
        "Choose your Planetary Body:",
        planet_list,
        key="planet_selectbox",
        help="Select a planetary body to use as a template for your model"
    )
    if selected_planet != "-- Select a Planet --":
        st.session_state["ChosenPlanet"] = selected_planet

chosen_planet = st.session_state.get("ChosenPlanet", None)

# --- Checks and load planet ---
if chosen_planet == "-- Select a Planet --":
    st.warning("Please select a planetary body to continue.")
    st.stop()
else:
    # Updates the session state to the chosen planet
    st.sidebar.markdown(f"**Current Planet:**  {chosen_planet}")
    st.success(f"✅ Using Planet: {chosen_planet}")

    # Show presets for this body in col2
    with col2:
        st.write("**Quick Start Presets** ❓")
        create_help_button("presets", "Presets provide pre-configured scenarios for common modeling cases")

        presets = get_presets_for_body(chosen_planet)
        if presets:
            preset_names = ["-- No Preset --"] + list(presets.keys())
            selected_preset = st.selectbox(
                "Load a preset scenario:",
                preset_names,
                key="selected_preset"
            )

            if selected_preset != "-- No Preset --":
                if st.button("📋 Apply Preset", help="Apply this preset configuration to your model"):
                    preset_config = presets[selected_preset]
                    apply_preset(preset_config)
                    st.success(f"Applied preset: {selected_preset}")
                    st.rerun()

                # Show preset summary
                if selected_preset in presets:
                    with st.expander("Preview Preset Settings"):
                        st.info(get_preset_summary(presets[selected_preset]))
        else:
            st.info("No presets available for this body")

    # Progress indicator
    create_progress_indicator(
        current_step=1,
        total_steps=6,
        step_names=[
            "Planet Selection",
            "Bulk Properties",
            "Ocean Settings",
            "Core & Mantle",
            "Layer Steps",
            "Run & Review"
        ]
    )

# --- Load planet module ---
from Utilities.PlanetLoader import load_planet_module

# --- Load default + check for custom modules ---
# --- Get custom modules ---
parent_directory = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # adjust as needed
custom_modules = get_custom_module_paths(parent_directory, chosen_planet)

custom_modules = get_custom_module_paths(parent_directory, chosen_planet)


use_custom = False
custom_module_path = None

# --- If custom modules exist, offer choice ---
if custom_modules:
    st.info("Custom modules found for this planet.")

    choice = st.radio(
        "How would you like to load this planet?",
        options=["Use default module", "Load a custom module"],
        index=0,
        key="custom_module_choice"
    )

    use_custom = (choice == "Load a custom module")

    if use_custom:
        custom_names = [os.path.basename(path) for path in custom_modules]
        selected_custom_name = st.selectbox("Choose a custom module:", custom_names, key="custom_module_name")
        custom_module_path = next(
            (p for p in custom_modules if os.path.basename(p) == selected_custom_name),
            None
        )

# --- Load module (default or custom) ---
try:
    planet_module = load_planet_module(parent_directory, chosen_planet, custom_module_path=custom_module_path)
    st.success(f"{chosen_planet} data loaded{' (custom)' if use_custom else ''}.")
    st.session_state["Planet"] = planet_module.Planet
except Exception as e:
    st.error(f"Error loading planet module: {e}")
    st.stop()



Planet = st.session_state.get("Planet", None)



# --- Main Settings Options Begin here- Tb_K and  Ocean Compositions ---

st.markdown("---")

# Create header with help button
col_header, col_help = st.columns([4, 1])
with col_header:
    st.subheader("Ice Layer Thickness")
with col_help:
    create_help_button("ice_shell", "The ice shell can be defined by thickness or bottom temperature")

# Get the default Tb_K from Planet object
default_Tb_K = getattr(Planet.Bulk, 'Tb_K')

# Initializing everything to session_state
if "Planet.Bulk.Tb_K" not in st.session_state:
    st.session_state["Planet.Bulk.Tb_K"] = default_Tb_K

if 'thickness_or_Tb' not in st.session_state:
    st.session_state['thickness_or_Tb'] = "Input Bottom Temperature Tb_K"

if 'zb_approximate_km' not in st.session_state:
    st.session_state['zb_approximate_km'] = None

if 'ICEIh_THICKNESS' not in st.session_state:
    st.session_state['ICEIh_THICKNESS'] = False

# Ensure st.session_state is ready to track changes -- we are reusing the bulk settings flags since Tb is a bulk setting
if "changed_bulk_settings_flags" not in st.session_state:
    st.session_state["changed_bulk_settings_flags"] = {}
if "changed_bulk_settings" not in st.session_state:
    st.session_state["changed_bulk_settings"] = {}

# Re-using this function to track changes to bulk settings
def on_change_bulk_setting(bulk_setting_key):
    st.session_state["changed_bulk_settings_flags"][bulk_setting_key] = True
    st.session_state["changed_bulk_settings"][bulk_setting_key] = st.session_state[bulk_setting_key]

# UI with improved styling
st.info("ℹ️ Choose how to define your ice shell configuration")

st.session_state['thickness_or_Tb'] = st.selectbox(
    "Ice Shell Configuration Method:",
    ("Input Ice Shell thickness", "Input Bottom Temperature Tb_K"),
    index=("Input Ice Shell thickness", "Input Bottom Temperature Tb_K").index(st.session_state['thickness_or_Tb']),
    help="Choose whether to specify ice shell thickness directly or via bottom temperature"
)

# If user chooses they want to do Ice Shell Thickness
if st.session_state['thickness_or_Tb'] == "Input Ice Shell thickness":
    zb_key = "Planet.Bulk.zb_approximate_km"
    do_flag_key = "Planet.Do.ICEIh_THICKNESS"

    # Creating Planet.Bulk.zb_approximate_km attribute
    if not hasattr(Planet.Bulk, "zb_approximate_km"):
        Planet.Bulk.zb_approximate_km = None

    # Initialize session state
    if zb_key not in st.session_state:
        default_zb = Planet.Bulk.zb_approximate_km if Planet.Bulk.zb_approximate_km is not None else 30.0
        st.session_state[zb_key] = default_zb

    col_input, col_info = st.columns([2, 1])

    with col_input:
        st.number_input(
            "Ice Shell Thickness (km):",
            key=zb_key,
            min_value=0.1,
            max_value=200.0,
            step=1.0,
            on_change=on_change_bulk_setting,
            args=(zb_key,),
            help="Typical values: Europa ~20-30 km, Enceladus ~20-25 km"
        )

    with col_info:
        # Show estimated bottom temperature (if known from preset)
        if st.session_state[zb_key]:
            st.metric("Specified Thickness", f"{st.session_state[zb_key]:.1f} km")

    # Track ICEIh_THICKNESS flag for changed outputs later
    st.session_state[do_flag_key] = True
    st.session_state["changed_bulk_settings_flags"][do_flag_key] = True
    st.session_state["changed_bulk_settings"][do_flag_key] = True

    with st.expander("ℹ️ How This Works"):
        st.write(
            "Planet Profile will use the inputted ice layer thickness to generate the ice shell for your planet. "
            "Based on the ice shell thickness, the temperature at the bottom of the Ice shell will be calculated."
        )

# If user chooses they want to do Tb_K
elif st.session_state['thickness_or_Tb'] == "Input Bottom Temperature Tb_K":
    tb_key = "Planet.Bulk.Tb_K"

    col_input, col_info = st.columns([2, 1])

    with col_input:
        st.number_input(
            "Bottom Temperature (K):",
            key=tb_key,
            min_value=100.0,
            max_value=400.0,
            step=1.0,
            on_change=on_change_bulk_setting,
            args=(tb_key,),
            help="Typical values: 250-273 K (melting point of ice)"
        )

    with col_info:
        # Show the specified temperature
        if st.session_state[tb_key]:
            st.metric("Bottom Temperature", f"{st.session_state[tb_key]:.1f} K")

    st.session_state['ICEIh_THICKNESS'] = False

    with st.expander("ℹ️ How This Works"):
        st.write(
            "The temperature you select at the bottom of the ocean layer for your planet is used by Planet Profile "
            "to determine the thickness of the Ice Shell. Behind the scenes, this sets Planet.Bulk.Tb_K."
        )
