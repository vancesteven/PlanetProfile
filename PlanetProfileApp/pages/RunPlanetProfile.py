import streamlit as st
import os
import sys
from pdf2image import convert_from_path
from PIL import Image
import re
import subprocess
from copy import deepcopy
import pandas as pd
import re
from pathlib import Path
from dataclasses import asdict, is_dataclass
import pprint
import shutil
import time

# ----- Page setup stuff -----
st.set_page_config(
    page_title="Run PlanetProfile",
    page_icon="🚀",
    layout="wide"
)

from Utilities.planet_sidebar import show_planet_status
from Utilities.app_helpers import (
    estimate_runtime, create_progress_indicator,
    format_timestamp, create_metric_card
)
from Utilities.help_system import create_help_button

show_planet_status()

# Header with help
col_title, col_help = st.columns([5, 1])
with col_title:
    st.title("Run PlanetProfile")
with col_help:
    create_help_button("run_simulation", "Review your configuration and run the simulation")

# Progress indicator
create_progress_indicator(
    current_step=6,
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

st.subheader("Summary of Your Planet and Changes you have made from Defaults:")

# ----- File Path management -----
# Get the planet name from the session state
Planet = st.session_state.get("Planet", None)
if not Planet:
    st.error("Please Select a Planet on the Planet Profile Main Settings Page")
    st.stop()

chosen_planet = st.session_state.get("ChosenPlanet", None)

# Get the path to the current script's directory
# /PlanetProfile/PlanetProfileApp/RunPlanetProfile.py
RunPlanetProfile_directory = os.path.dirname(os.path.abspath(__file__))

# Get the app directory (/PlanetProfile/PlanetProfileApp)
app_directory = os.path.dirname(RunPlanetProfile_directory)

# Get the parent directory (/PlanetProfile)
parent_directory  = os.path.dirname(app_directory)

# Add the parent directory to Python's search path.
if parent_directory not in sys.path:
    sys.path.append(parent_directory)

figures_folder = os.path.join(parent_directory, chosen_planet, "figures")



# ----- Setting Up SemiCustomPlanet object and Module Saving  -----

# This initializes the changed settings dictionaries in the event that the user did not visit that particular page
for key in ["changed_bulk_settings_flags", "custom_ocean_flag", "changed_step_settings_flags", "changed_core_settings"]:
    if key not in st.session_state:
        # For the boolean custom_ocean_flag, initialize as False; for the others, use empty dict
        st.session_state[key] = False if key == "custom_ocean_flag" else {}


# If the user has changed any inputs, we will create a semi-custom Planet object here to push their changes to
# Check if any setting has been changed
any_changes_made = (
    st.session_state["custom_ocean_flag"]
    or any(st.session_state["changed_bulk_settings_flags"].values())
    or any(st.session_state["changed_step_settings_flags"].values())
    or any(st.session_state["changed_core_settings"].values())
)

#Flags for keeping track of if module has been created or saved
if "semi_custom_created" not in st.session_state:
    st.session_state["semi_custom_created"] = False
if "module_already_saved" not in st.session_state:
    st.session_state["module_already_saved"] = False

def sanitize_filename(name):
    #Sanitize the filename by replacing invalid characters with underscores.
    return re.sub(r'\W|^(?=\d)', '_', name)


# This is used to set attributes to the planet object when the settings aren't all the same subtype (step settings has both Planet.Steps and Planet.Ocean, core settings has Planet.Core and Planet.Sil)
def set_nested_attr(obj, attr_path, value):
    """
    Recursively sets a nested attribute given a dotted path.
    E.g. set_nested_attr(obj, "Steps.nIceI", 5) sets obj.Steps.nIceI = 5
    """
    parts = attr_path.split(".")
    for part in parts[:-1]:
        obj = getattr(obj, part)
    setattr(obj, parts[-1], value)



# ----- Changed Settings Summary and Updating of SemiCustomPlanet with new changes -----
changed_settings_for_SemiCustom = {}
default_values_for_comments = {}  # For optional inline comment

if any_changes_made: #and st.session_state.get('semi_custom_created', False):
    # TEMP OBJECT to avoid modifying real planet before confirmation
    SemiCustomPlanet = deepcopy(Planet)

    # ----- Changed Bulk Settings Summary -----
    st.markdown("## Custom Bulk Settings Applied")
    # If the user has changed any bulk planetary settings, they will be updated into the SemiCustomPlanet object here
    for key, changed in st.session_state["changed_bulk_settings_flags"].items():
        if changed:
            attr = key.split(".")[-1]
            val = st.session_state["changed_bulk_settings"][key]
            setattr(SemiCustomPlanet.Bulk, attr, val)

            full_key = f"Planet.Bulk.{attr}"
            changed_settings_for_SemiCustom[full_key] = val
            default_values_for_comments[full_key] = getattr(Planet.Bulk, attr, "N/A")

    # This prints out for the user any of the bulk planetary settings they have changed from the default here
    if any(st.session_state["changed_bulk_settings_flags"].values()):
        st.warning("You have changed the following settings from the defaults: ")
        for key, was_changed in st.session_state["changed_bulk_settings_flags"].items():
            if was_changed:
                bulk_setting_name = key.split(".")[-1]
                bulk_new_val = st.session_state["changed_bulk_settings"][key]
                bulk_default_val = getattr(Planet.Bulk, bulk_setting_name, "N/A")

                st.markdown(f"- **{bulk_setting_name}**: `{bulk_default_val}` → `{bulk_new_val}`")
    else:
        st.info("No bulk settings have been changed. Using default bulk planetary settings.")
    st.markdown("---")


    # ----- Changed Ocean Settings Summary -----
    # If the user has changed any ocean settings, they will be updated into the SemiCustomPlanet object here and printed for the user
    st.markdown("## Custom Ocean Settings Applied")
    running_custom_ocean = st.session_state.get("custom_ocean_flag")
    custom_ocean_type = st.session_state.get("custom_ocean_comp", None)
    custom_ocean_ppt = st.session_state.get("custom_ocean_concentration", None)

    default_ocean_type = getattr(Planet.Ocean, "comp", "N/A")
    default_ocean_ppt = getattr(Planet.Ocean, "wOcean_ppt", "N/A")

    if running_custom_ocean == True:
        st.warning("You are using an ocean different than the default ocean")

        if custom_ocean_type:
            st.markdown(f"- **Ocean Composition**: `{default_ocean_type}` → `{custom_ocean_type}`")
            SemiCustomPlanet.Ocean.comp = custom_ocean_type
            changed_settings_for_SemiCustom["Planet.Ocean.comp"] = custom_ocean_type
            default_values_for_comments["Planet.Ocean.comp"] = default_ocean_type
        else:
            SemiCustomPlanet.Ocean.comp = default_ocean_type  # reset

        if custom_ocean_ppt is not None:
            st.markdown(f"- **Ocean Salinity (ppt)**: `{default_ocean_ppt}` → `{custom_ocean_ppt}`")
            SemiCustomPlanet.Ocean.wOcean_ppt = custom_ocean_ppt
            changed_settings_for_SemiCustom["Planet.Ocean.wOcean_ppt"] = custom_ocean_ppt
            default_values_for_comments["Planet.Ocean.wOcean_ppt"] = default_ocean_ppt
        else:
            SemiCustomPlanet.Ocean.wOcean_ppt = default_ocean_ppt  # reset




    if running_custom_ocean == False:
        st.info("No ocean settings have been changed. Using default ocean.")
    st.markdown("---")

    # ----- Changed Layer Step Settings Summary -----
    st.markdown("## Custom Layer Step Settings Applied")
    # If the user has changed any layer step settings, they will be updated into the SemiCustomPlanet object here and printed for the user
    for key, new_val in st.session_state.get("changed_step_settings", {}).items():
        # Remove the "Planet." prefix to get the attribute path
        attr_path = key.replace("Planet.", "", 1)
        # Set the value in SemiCustomPlanet
        set_nested_attr(SemiCustomPlanet, attr_path, new_val)
        full_key = key  # e.g., Planet.Steps.nIceI
        changed_settings_for_SemiCustom[full_key] = new_val

        parts = attr_path.split(".")
        default_obj = Planet
        try:
            for part in parts[:-1]:
                default_obj = getattr(default_obj, part)
            default_values_for_comments[full_key] = getattr(default_obj, parts[-1], "N/A")
        except Exception:
            default_values_for_comments[full_key] = "N/A"

    changed_flags = st.session_state.get("changed_step_settings_flags", {})
    changed_settings = st.session_state.get("changed_step_settings", {})

    if any(changed_flags.values()):
        st.warning("You have changed the following settings from the defaults:")

        for key, changed in changed_flags.items():
            if changed:
                new_val = changed_settings.get(key, "N/A")
                default_val = st.session_state.get(key, "N/A")  # Default was stored here at init

                # Strip just the final part of the setting name
                setting_name = key.split(".")[-1]

                st.markdown(f"- `{setting_name}`: `{default_val}` → `{new_val}`")
    else:
        st.info("No step settings have been changed. All values are defaults.")


    st.markdown("---")

    # ----- Changed Core and Silicate Settings Summary -----
    st.markdown("## Custom Core and Silicate Settings Applied")

    changed_core_settings = st.session_state.get("changed_core_settings", {})

    if changed_core_settings:
        st.warning("You have changed the following settings from the defaults:")

        for key, was_changed in st.session_state["changed_core_settings_flags"].items():
            if was_changed:
                new_val = st.session_state["changed_core_settings"][key]
                full_key = key  # e.g., Planet.Core.wFe_ppt

                # Apply to runtime object
                attr_path = key.replace("Planet.", "", 1)
                set_nested_attr(SemiCustomPlanet, attr_path, new_val)

                # Store in changed settings
                changed_settings_for_SemiCustom[full_key] = new_val

                # Get default value for inline comment
                parts = attr_path.split(".")
                default_obj = Planet
                try:
                    for part in parts[:-1]:
                        default_obj = getattr(default_obj, part)
                    default_val = getattr(default_obj, parts[-1], "N/A")
                except Exception:
                    default_val = "N/A"

                default_values_for_comments[full_key] = default_val

                # For display
                setting_name = parts[-1]
                st.markdown(f"- **{setting_name}**: `{default_val}` → `{new_val}`")

    else:
        st.info("No core or silicate settings have been changed. All values are defaults.")
    st.markdown("---")

    # ----- Naming of SemiCustomModule
    st.subheader("📝 Name Your Custom Configuration")

    # Default name suggestion
    default_name = "SemiCustom" + chosen_planet

    # Text input for custom name
    custom_planet_name = st.text_input(
        "Configuration name (used for the module filename):",
        value=default_name,
        key="custom_planet_name",
        help="This name will be used to create the custom module file (e.g., PPSemiCustomEuropa.py)"
    )

    st.info("💡 Your custom module will be automatically saved when you click 'Run PlanetProfile' below")
    st.markdown("---")
else:
    st.info("No changes detected. Running with default Planet module.")
    SemiCustomPlanet = Planet
    st.markdown("---")

# ----- Printing All Figures -----
config_path = os.path.join(parent_directory, "configPP.py") #path to configPP.py
#st.write(config_path)

# Import config
from UserConfigs.configPP import configAssign  # qualified: bare configPP only resolved when CWD held UserConfigs/ (2026-07-12 fix)  # This brings in the current config state from the configPP.py file
# Call the function to get Params and ExploreParams
Params, ExploreParams = configAssign() #configAssign creates the ParamsStruct and ExploreParamsStruct

exclude_plot = "SKIP_PLOTS"
exclude_calcs = {"CALC_ASYM", "CALC_NEW_ASYM"}

# Enable LaTeX table output for publication use
Params.DISP_TABLE = True

# FOR GUI: Skip plotting during run - figures will be generated from data on-demand
Params.SKIP_PLOTS = True

# Enable timestamp in filenames to avoid overwriting previous runs
Params.TIME_AND_DATE_LABEL = True

# Still save data files (we need them for figure generation)
Params.NO_SAVEFILE = False

# Automatically set all CALC attributes to True (but skip plots)
for attr in dir(Params):
    if "CALC" in attr and not attr.startswith("__"):
        setattr(Params, attr, attr not in exclude_calcs)



# ----- Functions for reading PP terminal output LaTeX tabulars and configuring them into tables in the GUI ----

def strip_ansi_codes(s):
    ansi_escape = re.compile(r'\x1b\[[0-9;]*[a-zA-Z]')
    return ansi_escape.sub('', s)

def extract_section_titles(latex_str):
    return re.findall(r'\\section\*{([^}]*)}', latex_str)

def clean_latex_cell(cell):
    """Clean LaTeX cell content for plain text display"""
    cell = cell.strip()

    # Remove math mode markers
    cell = re.sub(r'\$', '', cell)

    # Handle common LaTeX commands - just extract content
    cell = re.sub(r'\\textbf{([^}]*)}', r'\1', cell)
    cell = re.sub(r'\\ce{([^}]*)}', r'\1', cell)
    cell = re.sub(r'\\num{([^}]*)}', r'\1', cell)
    cell = re.sub(r'\\si{([^}]*)}', r'(\1)', cell)

    # Handle substack as ±
    cell = re.sub(r'\\substack\{([^}]*)\}', lambda m: ' ± '.join(line.strip() for line in m.group(1).split(r'\\')), cell)

    # Replace special characters
    cell = cell.replace(r'\pm', '±')
    cell = cell.replace(r'\times', '×')
    cell = cell.replace('~', ' ')

    # Clean out formatting commands
    cell = re.sub(r'\\mathrm{([^}]*)}', r'\1', cell)
    cell = re.sub(r'\\overline{([^}]*)}', r'\1', cell)

    # Remove any remaining LaTeX commands
    cell = re.sub(r'\\[a-zA-Z]+', '', cell)
    cell = re.sub(r'[{}]', '', cell)

    return cell.strip()

def parse_all_latex_tables(latex_str):
    tables_raw = re.findall(r'\\begin{tabular}.*?\\hline(.*?)\\end{tabular}', latex_str, re.DOTALL)
    #st.write("Raw tabular blocks found:", tables_raw)  #for testing to see how the tables are getting parsed

    tables = []

    for block in tables_raw:
        block = block.replace('\n', ' ')
        rows = [line.strip() for line in block.strip().split('\\\\') if line.strip()]

        parsed_rows = []

        for row in rows:
            cells = [clean_latex_cell(cell) for cell in row.split('&')]
            parsed_rows.append(cells)

        # Remove short rows or empty ones
        valid_rows = [r for r in parsed_rows if len(r) == 2]

        if len(valid_rows) >= 2 and len(valid_rows) / len(parsed_rows) >= 0.8:
            df = pd.DataFrame(valid_rows, columns=["Property", "Value"])
            tables.append(df)
        elif len(parsed_rows) >= 2 and len(parsed_rows[0]) > 2:
            headers = parsed_rows[0]
            data_rows = parsed_rows[1:]
            try:
                df = pd.DataFrame(data_rows, columns=headers)
                tables.append(df)
            except Exception as e:
                st.warning(f"Failed to parse structured table: {e}")
        else:
            st.warning("Skipping unrecognized table format.")

    return tables

def remove_latex_tables(raw_text):
    # Remove everything between \begin{tabular} and \end{tabular}, including those lines
    return re.sub(r'\\begin{tabular}.*?\\end{tabular}', '', raw_text, flags=re.DOTALL).strip()

# ----- Run Planet Profile Button and Outputs -----
st.markdown("---")
st.subheader("🚀 Ready to Run")

# Determine which module to run
if any_changes_made:
    custom_name = st.session_state.get("custom_planet_name", f"SemiCustom{chosen_planet}").strip()
    module_to_run = f"{sanitize_filename(custom_name)}"
else:
    module_to_run = chosen_planet  # e.g., PPEuropa

# Estimate runtime
Planet = st.session_state.get("Planet", None)
if Planet:
    estimated_time = estimate_runtime(Planet)

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Configuration", module_to_run)
    with col2:
        st.metric("Estimated Runtime", estimated_time)
    with col3:
        changes_count = len([k for k, v in st.session_state.get("changed_bulk_settings_flags", {}).items() if v])
        changes_count += len([k for k, v in st.session_state.get("changed_step_settings_flags", {}).items() if v])
        st.metric("Custom Settings", changes_count)

if any_changes_made:
    st.info("💡 Your custom configuration will be automatically saved as a module file before running.")
else:
    st.info("💡 Running with default configuration. The simulation may take several minutes.")

# Output options
st.markdown("#### Output Options")
col_opt1, col_opt2 = st.columns(2)
with col_opt1:
    disp_latex_tables = st.checkbox(
        "Generate LaTeX tables in terminal output",
        value=True,
        key="disp_latex_tables",
        help="Enable this to include publication-ready LaTeX tables in the terminal output"
    )
with col_opt2:
    disp_layers = st.checkbox(
        "Display layer summary",
        value=True,
        key="disp_layers",
        help="Include detailed layer-by-layer summary in output"
    )

# Just the run button - no need for separate save button since saving is automatic
run_button = st.button("▶️ Run PlanetProfile", type="primary", use_container_width=True)

if run_button:
    start_time = time.time()

    # Progress bar setup
    progress_bar = st.progress(0)
    status_text = st.empty()

    status_text.text("🔧 Preparing simulation...")
    progress_bar.progress(10)

    # If user made changes, ensure the module file is saved before running
    if any_changes_made:
        # Get custom name from session state, or use default if not set
        name = st.session_state.get("custom_planet_name", f"SemiCustom{chosen_planet}").strip()
        if not name:  # If empty string
            name = f"SemiCustom{chosen_planet}"
        sanitized_name = sanitize_filename(name)
        module_filename = f"PP{sanitized_name}.py"
        planet_folder = os.path.join(parent_directory, chosen_planet)
        output_path = os.path.join(planet_folder, module_filename)
        original_path = os.path.join(planet_folder, f"PP{chosen_planet}.py")

        status_text.text("💾 Saving custom module...")

        # Load original file
        try:
            with open(original_path, "r") as f:
                original_lines = f.readlines()

            # Update lines with custom settings
            updated_lines = []
            for line in original_lines:
                updated = False
                for full_key, new_val in changed_settings_for_SemiCustom.items():
                    pattern = rf"^\s*{re.escape(full_key)}\s*="
                    if re.match(pattern, line):
                        if isinstance(new_val, str):
                            new_val_str = f'"{new_val}"'
                        elif isinstance(new_val, bool):
                            new_val_str = str(new_val)
                        else:
                            new_val_str = repr(new_val)
                        default_val = default_values_for_comments.get(full_key, "N/A")
                        line = f"{full_key} = {new_val_str}  # changed from default: {default_val}\n"
                        updated = True
                        break
                updated_lines.append(line)

            # Write updated module
            with open(output_path, "w") as f:
                f.writelines(updated_lines)

            st.info(f"✓ Saved custom module: {module_filename}")
        except Exception as e:
            st.error(f"Failed to save custom module: {e}")
            st.stop()
    else:
        name = chosen_planet
        sanitized_name = sanitize_filename(name)
        module_filename = f"PP{sanitized_name}.py"

    # Use the parent_directory that was already correctly calculated at the top of the file
    # This points to the PlanetProfile root directory
    status_text.text(f"📂 Working directory: {parent_directory}")
    progress_bar.progress(20)

    # Command to actually run PlanetProfile
    if not any_changes_made:
        command = ["python", "PlanetProfileCLI.py", str(module_to_run)]
    else:
        full_path = f"{chosen_planet}/{module_filename}"
        command = ["python", "PlanetProfileCLI.py", full_path]

    status_text.text(f"🚀 Running PlanetProfile: {module_to_run}")
    progress_bar.progress(30)

    # Create a container for real-time terminal output
    terminal_output = st.empty()
    output_lines = []

    # Set environment variables to control PlanetProfile output
    env = os.environ.copy()
    env['PP_DISP_TABLE'] = '1' if st.session_state.get('disp_latex_tables', True) else '0'
    env['PP_DISP_LAYERS'] = '1' if st.session_state.get('disp_layers', True) else '0'

    try:
        # Use Popen for streaming output
        process = subprocess.Popen(
            command,
            cwd=parent_directory,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
            env=env
        )

        # Stream output in real-time
        progress_value = 30
        for line in process.stdout:
            output_lines.append(line)
            # Update terminal display
            terminal_output.code('\n'.join(output_lines[-20:]), language='text')

            # Increment progress gradually
            if progress_value < 90:
                progress_value += 0.5
                progress_bar.progress(int(progress_value))

        process.wait(timeout=600)
        elapsed_time = time.time() - start_time

        output_str = ''.join(output_lines)
        result = type('obj', (object,), {'returncode': process.returncode, 'stdout': output_str, 'stderr': ''})()

        progress_bar.progress(100)
        status_text.empty()
        progress_bar.empty()
        terminal_output.empty()  # Clear the streaming output

        if result.returncode == 0:
            st.success(f"✅ PlanetProfile completed successfully in {elapsed_time:.1f} seconds!")

            # Show completion metrics
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Status", "Success ✅")
            with col2:
                st.metric("Runtime", f"{elapsed_time:.1f}s")
            with col3:
                estimated = estimate_runtime(Planet)
                st.metric("vs Estimated", estimated)

            st.markdown("---")
        else:
            st.error(f"❌ PlanetProfile encountered an error (exit code: {result.returncode})")

            # Show the command that was run
            st.code(f"Command: {' '.join(command)}\nWorking directory: {parent_directory}", language="bash")

            # Try to extract error message
            error_lines = [line for line in output_str.split('\n') if 'error' in line.lower() or 'exception' in line.lower() or 'traceback' in line.lower()]
            if error_lines:
                with st.expander("🔍 Error Details (filtered)", expanded=True):
                    st.code('\n'.join(error_lines[:20]))

            # Always show recent output
            recent_lines = output_str.split('\n')[-30:]  # Last 30 lines
            with st.expander("📋 Recent Output (last 30 lines)", expanded=True):
                st.code('\n'.join(recent_lines))

            st.markdown("---")

    except subprocess.TimeoutExpired:
        process.kill()
        progress_bar.empty()
        status_text.empty()
        terminal_output.empty()
        st.error("❌ Simulation timed out after 10 minutes. Check your configuration or increase timeout.")
        output_str = ''.join(output_lines) if output_lines else "Simulation timed out"
        result = type('obj', (object,), {'returncode': -1, 'stdout': output_str, 'stderr': 'Timeout'})()

    except Exception as e:
        if 'process' in locals():
            process.kill()
        progress_bar.empty()
        status_text.empty()
        if 'terminal_output' in locals():
            terminal_output.empty()
        st.error(f"❌ Unexpected error: {str(e)}")
        output_str = ''.join(output_lines) if output_lines else str(e)
        result = type('obj', (object,), {'returncode': -1, 'stdout': output_str, 'stderr': str(e)})()

    st.markdown("---")

    # Strip ANSI codes for display
    display_output = strip_ansi_codes(output_str)

    # --- Parse LaTeX tables from output ---
    tables = parse_all_latex_tables(display_output)
    titles = extract_section_titles(display_output)
    if tables:
        st.markdown(f"### PlanetProfile {chosen_planet}")
        for i, df in enumerate(tables):
            title = titles[i] if i < len(titles) else f"Table {i+1}"

            # First table is Inputs, second is Outputs
            if i == 0:
                section_type = "Inputs:"
                display_title = section_type  # Don't include planet name for Inputs
            elif i == 1:
                section_type = "Outputs:"
                display_title = section_type  # Don't include "Table 2" for Outputs
            else:
                section_type = f"Table {i+1}:"
                display_title = f"{section_type} {title}"

            st.subheader(display_title)
            st.dataframe(df, use_container_width=True)

    # Display raw terminal output - keep LaTeX intact for publication use
    st.markdown("---")
    st.markdown("### Raw Terminal Output")
    st.markdown(
        f"""
        <div style="height: 400px; overflow-y: scroll; background-color: #111; color: #eee; padding: 10px; font-family: monospace; white-space: pre-wrap; border: 1px solid #666;">
            {display_output}
        </div>
        """,
        unsafe_allow_html=True,
    )



else:
    st.warning("Choices have not yet been pushed to PlanetProfile")

st.markdown("---")
