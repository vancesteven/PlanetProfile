import streamlit as st
from pdf2image import convert_from_path
import os
import sys
import pandas as pd
import re
from io import StringIO
import altair as alt
import time
from collections import defaultdict
from datetime import datetime
import json

# Import figure generation utilities
from Utilities.figure_generator import (
    load_profile_data,
    load_ocean_data,
    create_hydrosphere_plot,
    create_seismic_plot,
    create_custom_plot,
    save_plotly_figure,
    PLOTLY_AVAILABLE
)

# ----- Page setup stuff -----
st.set_page_config(
    page_title="PlanetProfile Outputs",
    page_icon="📈",
    layout="wide"
)

from Utilities.planet_sidebar import show_planet_status
from Utilities.app_helpers import format_timestamp, create_metric_card, create_progress_indicator
from Utilities.help_system import create_help_button

show_planet_status()

# Header with help
col_title, col_help = st.columns([5, 1])
with col_title:
    st.title("PlanetProfile Outputs")
with col_help:
    create_help_button("outputs", "View and analyze simulation results")

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

st.subheader("Outputs Produced by PlanetProfile will Appear Below")

# Run naming section
st.markdown("---")
col1, col2 = st.columns([3, 1])
with col1:
    st.info("💡 Tip: Name your simulation runs to easily identify different configurations")
with col2:
    if st.button("🏷️ Name Current Run", use_container_width=True):
        st.session_state['show_name_dialog'] = True

# Show naming dialog if requested
if st.session_state.get('show_name_dialog', False):
    with st.form("name_run_form"):
        run_name = st.text_input(
            "Enter a descriptive name for this run:",
            placeholder="e.g., Thick Ice Shell Test, High Salinity Ocean, etc.",
            key="run_name_input"
        )
        col_submit, col_cancel = st.columns(2)
        with col_submit:
            submitted = st.form_submit_button("💾 Save Name", use_container_width=True)
        with col_cancel:
            cancelled = st.form_submit_button("Cancel", use_container_width=True)

        if submitted and run_name:
            # Save run name to metadata file
            metadata_file = os.path.join(planet_folder, ".run_metadata.json")
            metadata = {}
            if os.path.exists(metadata_file):
                with open(metadata_file, 'r') as f:
                    metadata = json.load(f)

            # Find most recent batch and associate name
            if pdf_file_paths:
                latest_time = max(t for _, t in pdf_file_paths)
                metadata[str(latest_time)] = {
                    'name': run_name,
                    'timestamp': format_timestamp(latest_time)
                }

            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)

            st.success(f"✅ Named this run: {run_name}")
            st.session_state['show_name_dialog'] = False
            st.rerun()

        if cancelled:
            st.session_state['show_name_dialog'] = False
            st.rerun()

st.markdown("---")

# ----- File Path management -----
# Get the planet name from the session state. Fresh sessions (e.g. a
# public visitor browsing the shipped demo outputs) default to Europa —
# same convention as the Exploreogram page — instead of stopping.
Planet = st.session_state.get("Planet", None)
chosen_planet = st.session_state.get("ChosenPlanet", None) or 'Europa'
if not Planet:
    st.info("No planet selected — showing the Europa demo outputs. "
            "Select a planet on the Main Settings page to browse your "
            "own runs.")

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
planet_folder = os.path.join(parent_directory, chosen_planet)
st.write(f"# {chosen_planet} Figures")


# ----- Batching and timestamping figures  -----
BATCH_TIME_THRESHOLD_SECONDS = 60  # Groups files created within 60 seconds (can change this threshold as needed)

# Loads and timestamps all PDFs
pdf_file_paths = []
for filename in os.listdir(figures_folder):
    if filename.endswith(".pdf"):
        path = os.path.join(figures_folder, filename)
        mod_time = os.path.getmtime(path)
        pdf_file_paths.append((path, mod_time))

# Sorts PDF files by modification time
pdf_file_paths.sort(key=lambda x: x[1])

# Empty dictionaries to add batches to
batches = []
current_batch = []
prev_time = None
# Grouping the PDF files into batches
for path, mod_time in pdf_file_paths:
    if prev_time is None or (mod_time - prev_time) <= BATCH_TIME_THRESHOLD_SECONDS:
        current_batch.append((path, mod_time))
    else:
        batches.append(current_batch)
        current_batch = [(path, mod_time)]
    prev_time = mod_time
# adds files to batches
if current_batch:
    batches.append(current_batch)

# Load run metadata if exists
metadata_file = os.path.join(planet_folder, ".run_metadata.json")
run_metadata = {}
if os.path.exists(metadata_file):
    try:
        with open(metadata_file, 'r') as f:
            run_metadata = json.load(f)
    except:
        run_metadata = {}

# Function to create tab titles - uses named runs if available, otherwise timestamps
def format_batch_label(batch, run_metadata):
    if not batch:
        return "Unknown Batch"
    batch_time = batch[0][1]
    timestamp = datetime.fromtimestamp(batch_time)

    # Check if this batch has a named run
    time_key = str(batch_time)
    if time_key in run_metadata and 'name' in run_metadata[time_key]:
        run_name = run_metadata[time_key]['name']
        return f"🏷️ {run_name} ({timestamp.strftime('%Y-%m-%d %H:%M')})"
    else:
        return f"📅 {timestamp.strftime('%Y-%m-%d %H:%M:%S')}"

batch_labels = [format_batch_label(batch, run_metadata) for batch in batches]

# Summary statistics
st.markdown("### 📊 Output Summary")
col1, col2, col3, col4 = st.columns(4)
with col1:
    st.metric("Total Runs", len(batches))
with col2:
    named_runs = len([k for k, v in run_metadata.items() if 'name' in v])
    st.metric("Named Runs", named_runs)
with col3:
    st.metric("Total Figures", len(pdf_file_paths))
with col4:
    txt_count = len([f for f in os.listdir(planet_folder) if f.endswith('.txt')])
    st.metric("Data Files", txt_count)

st.markdown("---")

# Add option for interactive vs static figures
if PLOTLY_AVAILABLE:
    figure_mode = st.radio(
        "Figure Mode:",
        ["🎯 Interactive (Plotly)", "📄 Static PDFs (Legacy)"],
        horizontal=True,
        help="Interactive figures are generated from data files and allow zooming/panning. Legacy mode shows pre-generated PDFs."
    )
    use_interactive = figure_mode == "🎯 Interactive (Plotly)"
else:
    st.warning("⚠️ Plotly not installed. Only static PDF view available. Install with: pip install plotly")
    use_interactive = False

# Add option to view in different modes
view_mode = st.radio(
    "View Mode:",
    ["📅 Chronological (Batches)", "📁 By Figure Type", "🔍 Search"],
    horizontal=True,
    help="Choose how to organize and view your outputs"
)

st.markdown("---")

# ----- Interactive Figures from Data -----
if use_interactive and PLOTLY_AVAILABLE:
    st.write(f"## Interactive {chosen_planet} Figures from Data")

    # Find most recent profile data file
    profile_files = [f for f in os.listdir(planet_folder)
                    if f.endswith('.txt') and 'liquidOceanProps' not in f and 'Core' not in f]

    if profile_files:
        # Sort by modification time to get most recent
        profile_files_with_time = [(f, os.path.getmtime(os.path.join(planet_folder, f)))
                                   for f in profile_files]
        profile_files_with_time.sort(key=lambda x: x[1], reverse=True)

        # Create tabs for different interactive figure types
        interactive_tabs = st.tabs(["Hydrosphere Properties", "Seismic Properties", "Custom Plot"])

        with interactive_tabs[0]:
            st.subheader("Hydrosphere Properties")

            # Let user select which run to visualize
            run_options = [f"{f} ({datetime.fromtimestamp(t).strftime('%Y-%m-%d %H:%M:%S')})"
                          for f, t in profile_files_with_time[:5]]  # Show 5 most recent
            if run_options:
                selected_run = st.selectbox("Select run to visualize:", run_options, key="hydro_run")
                selected_file = profile_files_with_time[run_options.index(selected_run)][0]

                try:
                    with st.spinner("Loading data and generating figure..."):
                        df = load_profile_data(os.path.join(planet_folder, selected_file))
                        fig = create_hydrosphere_plot(df, chosen_planet)

                        if fig:
                            st.plotly_chart(fig, use_container_width=True)

                            # Add save button
                            col1, col2, col3 = st.columns([1, 1, 2])
                            with col1:
                                save_format = st.selectbox("Format:", ["html", "png", "pdf"], key="hydro_format")
                            with col2:
                                if st.button("💾 Save Figure", key="save_hydro"):
                                    save_path = os.path.join(planet_folder, "figures",
                                                            f"{chosen_planet}_Hydrosphere_interactive.{save_format}")
                                    if save_plotly_figure(fig, save_path, save_format):
                                        st.success(f"✅ Saved to: {save_path}")
                except Exception as e:
                    st.error(f"❌ Error generating figure: {e}")
                    st.info("💡 Falling back to static PDF view below")

        with interactive_tabs[1]:
            st.subheader("Seismic Properties")

            run_options = [f"{f} ({datetime.fromtimestamp(t).strftime('%Y-%m-%d %H:%M:%S')})"
                          for f, t in profile_files_with_time[:5]]
            if run_options:
                selected_run = st.selectbox("Select run to visualize:", run_options, key="seismic_run")
                selected_file = profile_files_with_time[run_options.index(selected_run)][0]

                try:
                    with st.spinner("Loading data and generating figure..."):
                        df = load_profile_data(os.path.join(planet_folder, selected_file))
                        fig = create_seismic_plot(df, chosen_planet)

                        if fig:
                            st.plotly_chart(fig, use_container_width=True)

                            # Add save button
                            col1, col2, col3 = st.columns([1, 1, 2])
                            with col1:
                                save_format = st.selectbox("Format:", ["html", "png", "pdf"], key="seismic_format")
                            with col2:
                                if st.button("💾 Save Figure", key="save_seismic"):
                                    save_path = os.path.join(planet_folder, "figures",
                                                            f"{chosen_planet}_Seismic_interactive.{save_format}")
                                    if save_plotly_figure(fig, save_path, save_format):
                                        st.success(f"✅ Saved to: {save_path}")
                except Exception as e:
                    st.error(f"❌ Error generating figure: {e}")
                    st.info("💡 Falling back to static PDF view below")

        with interactive_tabs[2]:
            st.subheader("Custom Plot from Data")

            run_options = [f"{f} ({datetime.fromtimestamp(t).strftime('%Y-%m-%d %H:%M:%S')})"
                          for f, t in profile_files_with_time[:5]]
            if run_options:
                selected_run = st.selectbox("Select run to visualize:", run_options, key="custom_run")
                selected_file = profile_files_with_time[run_options.index(selected_run)][0]

                try:
                    df = load_profile_data(os.path.join(planet_folder, selected_file))

                    col1, col2 = st.columns(2)
                    with col1:
                        x_col = st.selectbox("X-axis:", df.columns, key="custom_x")
                    with col2:
                        y_cols = st.multiselect("Y-axis (can select multiple):",
                                               [c for c in df.columns if c != x_col],
                                               key="custom_y")

                    col3, col4 = st.columns(2)
                    with col3:
                        log_x = st.checkbox("Log scale X", key="custom_log_x")
                    with col4:
                        log_y = st.checkbox("Log scale Y", key="custom_log_y")

                    if y_cols:
                        with st.spinner("Generating custom plot..."):
                            fig = create_custom_plot(df, x_col, y_cols,
                                                    title=f"{chosen_planet} - {x_col} vs {', '.join(y_cols)}",
                                                    log_x=log_x, log_y=log_y)

                            if fig:
                                st.plotly_chart(fig, use_container_width=True)

                                # Add save button
                                col1, col2, col3 = st.columns([1, 1, 2])
                                with col1:
                                    save_format = st.selectbox("Format:", ["html", "png", "pdf"], key="custom_format")
                                with col2:
                                    if st.button("💾 Save Figure", key="save_custom"):
                                        save_path = os.path.join(planet_folder, "figures",
                                                                f"{chosen_planet}_Custom_interactive.{save_format}")
                                        if save_plotly_figure(fig, save_path, save_format):
                                            st.success(f"✅ Saved to: {save_path}")
                except Exception as e:
                    st.error(f"❌ Error generating custom plot: {e}")
    else:
        st.info("📊 No profile data files found. Run a simulation first to generate data.")

    st.markdown("---")

# ----- Static PDF Figures (Legacy or Fallback) -----
if not use_interactive or not PLOTLY_AVAILABLE:
    st.write(f"## {chosen_planet} Static Figures (PDFs)")

if view_mode == "📅 Chronological (Batches)":
    # Outer tabs - grouped by time
    st.write(f"## {chosen_planet} Figures by Run")
    outer_tabs = st.tabs(batch_labels)
elif view_mode == "📁 By Figure Type":
    # Group by figure type instead
    st.write(f"## {chosen_planet} Figures by Type")
    outer_tabs = None  # Will handle differently below
elif view_mode == "🔍 Search":
    st.write(f"## Search {chosen_planet} Figures")
    search_term = st.text_input("🔍 Search figures by name or type:")
    outer_tabs = None  # Will handle search results differently



# ----- Figure Printing and Setup of Captions and Figure Types -----
pdf_files = [f for f in os.listdir(figures_folder) if f.endswith(".pdf")]

if not pdf_files:
    st.warning(f"No figure PDFs found in: {figures_folder}")
    st.stop()

figure_dict = {} #empty so we can update it later

# This dictionary is used to parse the figure file names and link them to figure titles in the GUI
figure_types = {
    "Gravity" : "Gravity and Pressure",
    "Wedge" : "Interior Wedge Diagram",
    "Hydrosphere" : "Hydrosphere Properties",
    "CoreMantTrade" : "Silicate-Core Size Tradeoff",
    "MantleDens" : "Silicate Radius-Density Handoff",
    "Seismic" : "Seismic Properties",
    "Viscosity" : "Viscosity",
    "Porosity2axes": "Porosity with 2 axes",
    "Porosity" : "Porosity"}


# This does the actual parsing of the figure file names
for filename in pdf_files:
    if filename.endswith('.pdf'):
        name_only = filename[:-4]
        matched_keyword = None

        for keyword in figure_types:
            if keyword in name_only:
                matched_keyword = keyword
                break

        # Uses slightly more descriptive title from the figure_types dictionary above
        label = figure_types.get(matched_keyword, matched_keyword if matched_keyword else name_only)
        figure_dict[label] = os.path.join(figures_folder, filename)


#Below are descriptive captions for the figures, linked via dictionary to the titles for those figures
captions = {
    "Interior Wedge Diagram": "Shows an interior wedge diagram showing the calculations of the radii of each planet layer",
     "Gravity and Pressure": "Gravitational acceleration (g) and Pressure profiles as a function of radius",
    "Hydrosphere Properties": (
        "Interior Properties- \n\n"
        "Left: Phase diagram as a funciton of pressure and density. \n\n"
        "Right(Top): Temperature profile across different depths.\n\n"
        "Right(center): Longitudinal (p-wave) sound velocity Vp for each layer in km/s and shear (s-wave) sound velocity Vs for each layer in km/s as a funciton of depth, \n\n"
        "Right(Bottom): Electrical conductivity as a funciton of depth"),

    "Silicate-Core Size Tradeoff" : "Silicate–core size tradeoff - Based on given Moment of Inertia, calculates all profiles of silicate and core size pairs that fit within that gien MOI. The one determined by Planet Profile to most closely match the given MOI is marked on the figure",
    "Silicate Radius-Density Handoff" : "Silicate radius-density tradeoff - Based on given Moment of Inertia, calculates all profiles of silicate radii and densities that fit within that gien MOI. The one determined by Planet Profile to most closely match the given MOI is marked on the figure",
    "Porosity" : "Left: Porosity as a funciton of depth. Right: Porosity as a funciton of pressure",
    "Porosity with 2 axes": "Displays Porosity both as a function of depth and as a funciton of pressure on one figure" ,
    "Viscosity" : "Viscosity of the planet layers as a funciton of radius",
    "Seismic Properties" : (
        "Seismic Properties. \n\n"
        "Top Left - Bulk & shear moduli Ks & Gs as a function of radius. \n\n"
        "Top Right - Temperature, Pressure, and density as a funciton of radius. \n\n"
        "Bottom Left - Sound speeds Vp and Vs as a function of radius. \n\n"
        "Bottom Right - Seismic quality factor Qs as a funciton of raidus")
}


# Display figures based on view mode
if view_mode == "📅 Chronological (Batches)" and outer_tabs:
    # Original chronological view
    for tab, batch in zip(outer_tabs, batches):
        with tab:
            # --- Map figure types for this batch ---
            figure_dict = {}

            for filepath, _ in batch:
                filename = os.path.basename(filepath)
                name_only = filename[:-4]

                matched_keyword = None
                for keyword in figure_types:
                    if keyword in name_only:
                        matched_keyword = keyword
                        break

                label = figure_types.get(matched_keyword, matched_keyword if matched_keyword else name_only)
                figure_dict[label] = filepath

            # Sort figure labels alphabetically for consistent order
            figure_labels = sorted(figure_dict.keys())
            inner_tabs = st.tabs(figure_labels)

            # --- Display each figure in inner tab ---
            for fig_tab, label in zip(inner_tabs, figure_labels):
                with fig_tab:
                    pdf_path = figure_dict[label]
                    with st.spinner(f"Rendering: {label}..."):
                        try:
                            images = convert_from_path(pdf_path)
                            st.image(images[0], use_container_width=True)
                        except Exception as _pdf_exc:
                            st.warning(f"PDF preview unavailable ({_pdf_exc}); "
                                       "download instead:")
                            with open(pdf_path, 'rb') as _fpdf:
                                st.download_button("⬇️ Download PDF", _fpdf,
                                                   file_name=os.path.basename(pdf_path),
                                                   key=f"dl_{pdf_path}")
                    st.caption(f"**{captions.get(label, label)}**")

elif view_mode == "📁 By Figure Type":
    # Group all figures by type across all batches
    type_groups = defaultdict(list)

    for filepath, mod_time in pdf_file_paths:
        filename = os.path.basename(filepath)
        name_only = filename[:-4]

        matched_keyword = None
        for keyword in figure_types:
            if keyword in name_only:
                matched_keyword = keyword
                break

        label = figure_types.get(matched_keyword, matched_keyword if matched_keyword else name_only)
        type_groups[label].append((filepath, mod_time))

    # Create tabs for each figure type
    type_labels = sorted(type_groups.keys())
    type_tabs = st.tabs(type_labels)

    for type_tab, label in zip(type_tabs, type_labels):
        with type_tab:
            st.write(f"### {label}")
            st.caption(f"**{captions.get(label, label)}**")

            # Show all instances of this figure type
            for i, (filepath, mod_time) in enumerate(type_groups[label]):
                timestamp = datetime.fromtimestamp(mod_time).strftime('%Y-%m-%d %H:%M:%S')
                with st.expander(f"Run {i+1}: {timestamp}", expanded=(i==len(type_groups[label])-1)):
                    with st.spinner(f"Rendering {label}..."):
                        try:
                            images = convert_from_path(filepath)
                            st.image(images[0], use_container_width=True)
                        except Exception as _pdf_exc:
                            st.warning(f"PDF preview unavailable ({_pdf_exc}); "
                                       "download instead:")
                            with open(filepath, 'rb') as _fpdf:
                                st.download_button("⬇️ Download PDF", _fpdf,
                                                   file_name=os.path.basename(filepath),
                                                   key=f"dl_{filepath}_{mod_time}")

elif view_mode == "🔍 Search":
    # Search functionality
    if 'search_term' in locals() and search_term:
        search_lower = search_term.lower()
        matching_figures = []

        for filepath, mod_time in pdf_file_paths:
            filename = os.path.basename(filepath)
            name_only = filename[:-4]

            matched_keyword = None
            for keyword in figure_types:
                if keyword in name_only:
                    matched_keyword = keyword
                    break

            label = figure_types.get(matched_keyword, matched_keyword if matched_keyword else name_only)

            # Check if search term matches
            if (search_lower in label.lower() or
                search_lower in name_only.lower() or
                search_lower in captions.get(label, "").lower()):
                matching_figures.append((filepath, mod_time, label))

        if matching_figures:
            st.success(f"Found {len(matching_figures)} matching figure(s)")

            for filepath, mod_time, label in matching_figures:
                timestamp = datetime.fromtimestamp(mod_time).strftime('%Y-%m-%d %H:%M:%S')
                with st.expander(f"{label} - {timestamp}", expanded=True):
                    st.caption(f"**{captions.get(label, label)}**")
                    with st.spinner(f"Rendering {label}..."):
                        try:
                            images = convert_from_path(filepath)
                            st.image(images[0], use_container_width=True)
                        except Exception as _pdf_exc:
                            st.warning(f"PDF preview unavailable ({_pdf_exc}); "
                                       "download instead:")
                            with open(filepath, 'rb') as _fpdf:
                                st.download_button("⬇️ Download PDF", _fpdf,
                                                   file_name=os.path.basename(filepath),
                                                   key=f"dl_{filepath}_{mod_time}")
        else:
            st.warning(f"No figures found matching '{search_term}'")
    else:
        st.info("Enter a search term above to find specific figures")

st.markdown("---")




# ----- .txt. files loading for the GUI -----
st.write(f"# {chosen_planet} text files")
# Get list of .txt files and their modification times
# Step 1: Get list of .txt files and their modification times
txt_files = [
    (f, os.path.getmtime(os.path.join(planet_folder, f)))
    for f in os.listdir(planet_folder)
    if f.endswith(".txt")
]

# Sort by modification time
txt_files.sort(key=lambda x: x[1])

# Step 2: Group files into batches by timestamp (90 seconds apart)
threshold = 90  # seconds
batches = []
current_batch = []
last_ts = None

for file, ts in txt_files:
    if last_ts is None or (ts - last_ts) <= threshold:
        current_batch.append((file, ts))
    else:
        batches.append(current_batch)
        current_batch = [(file, ts)]
    last_ts = ts

if current_batch:
    batches.append(current_batch)

# Format batch labels
formatted_batches = []
for i, batch in enumerate(batches):
    timestamp = datetime.fromtimestamp(batch[0][1]).strftime('%Y-%m-%d %H:%M:%S')
    label = f"Batch {i + 1} ({timestamp})"
    file_list = [f for f, _ in batch]
    formatted_batches.append((label, file_list))

st.write(f"# {chosen_planet} text files grouped by batch and file type")

# Step 3: Outer tabs for batches
batch_labels = [label for label, files in formatted_batches]
batch_tabs = st.tabs(batch_labels)

for batch_idx, (batch_label, files) in enumerate(formatted_batches):
    with batch_tabs[batch_idx]:
        # Step 4: Inner tabs for file types
        inner_tabs = st.tabs(["Ocean text File", "Core text File", "Profile text File"])

        # Prepare dict to hold files by type
        files_by_type = {
            "Ocean": [f for f in files if f.endswith("liquidOceanProps.txt")],
            "Core": [f for f in files if f.endswith("Core.txt")],
            "Profile": [f for f in files if not (f.endswith("Core.txt") or f.endswith("liquidOceanProps.txt"))]
        }



# Read all .txt files in the folder
#txt_files = [f for f in os.listdir(planet_folder) if f.endswith(".txt")]



# ----- Printing of the liquidOceanProps.txt file table and scatter plot -----
     # ----- Ocean files -----
        with inner_tabs[0]:
            for ocean_file in files_by_type["Ocean"]:
                file_path = os.path.join(planet_folder, ocean_file)
                with open(file_path, "r", encoding="utf-8") as f:
                    lines = f.readlines()

                # The file declares its own header length ("nHeadLines = N"
                # on line 1); the Nth line is the column-name row. The old
                # hardcoded [1:4]/[4]/[5:] slices broke on files with more
                # header metadata (e.g. ocean speciation lines).
                n_head = 5
                if 'nHeadLines' in lines[0]:
                    try:
                        n_head = int(lines[0].split('=')[1])
                    except (IndexError, ValueError):
                        pass
                header_lines = "".join(lines[1:n_head - 1]).strip()
                header_line = lines[n_head - 1].strip()
                column_names = re.split(r'\s{2,}|\t+', header_line)
                data_lines = "".join(lines[n_head:])

                try:
                    # Try parsing with flexible whitespace separator
                    df_ocean = pd.read_csv(StringIO(data_lines), sep=r'\s+', header=None,
                                          engine='python', on_bad_lines='warn')
                except Exception as e:
                    st.error(f"❌ Failed to parse {ocean_file}: {str(e)}")
                    st.info("Attempting to show first few lines of data for debugging:")
                    st.code('\n'.join(lines[5:15]))
                    continue

                st.subheader(f"Ocean File Header Info: {ocean_file}")
                with st.expander("Show Ocean Header"):
                    st.code(header_lines)

                if len(column_names) != df_ocean.shape[1]:
                    st.warning(f"⚠️ Column count mismatch: Header has {len(column_names)} columns, data has {df_ocean.shape[1]}")
                    df_ocean.columns = [f"Col_{i+1}" for i in range(df_ocean.shape[1])]
                else:
                    df_ocean.columns = column_names

                st.dataframe(df_ocean, use_container_width=True)
                st.subheader("Scatter Plot")

                persist_plot = st.checkbox("Plot multiple Y-axes against same X-axis", key=f"persist_{ocean_file}")
                x_axis = st.selectbox("X-axis", df_ocean.columns, key=f"x_{ocean_file}")

                if persist_plot:
                    y_axes = st.multiselect(
                        "Y-axis (select one or more)",
                        df_ocean.columns.drop(x_axis),
                        default=[df_ocean.columns[1]] if x_axis != df_ocean.columns[1] else df_ocean.columns[2:3],
                        key=f"multi_y_{ocean_file}"
                    )
                else:
                    y_axes = [st.selectbox("Y-axis", df_ocean.columns.drop(x_axis), key=f"single_y_{ocean_file}")]

                valid_y_axes = [y for y in y_axes if y != x_axis]

                if valid_y_axes:
                    # Convert columns to numeric (coerce errors to NaN)
                    df_ocean[x_axis] = pd.to_numeric(df_ocean[x_axis], errors='coerce')
                    for col in valid_y_axes:
                        df_ocean[col] = pd.to_numeric(df_ocean[col], errors='coerce')

                    df_plot = df_ocean[[x_axis] + valid_y_axes].dropna()
                    df_plot = df_plot.set_index(x_axis)[valid_y_axes]

                    if not df_plot.empty:
                        st.scatter_chart(df_plot)
                    else:
                        st.info("No valid data to plot after cleaning.")
                else:
                    st.info("Please select a Y-axis different from the X-axis to plot.")

        # ----- Core files -----
        with inner_tabs[1]:
            core_column_names = ["RsilTrade (m)", "RcoreTrade (m)", "rhoSilTrade (kg/m³)"]
            for core_file in files_by_type["Core"]:
                file_path = os.path.join(planet_folder, core_file)
                with open(file_path, "r", encoding="utf-8") as f:
                    lines = f.readlines()

                data_only = "".join(lines[1:])
                df_core = pd.read_csv(StringIO(data_only), sep='\s+', header=None)
                df_core.columns = core_column_names

                st.subheader(f"Core Data Table: {core_file}")
                st.dataframe(df_core, use_container_width=True)

                st.subheader("Scatter Plot")

                persist_plot = st.checkbox("Plot multiple Y-axes against same X-axis", key=f"persist_core_{core_file}")
                x_axis = st.selectbox("X-axis", df_core.columns, key=f"x_core_{core_file}")

                if persist_plot:
                    y_axes = st.multiselect(
                        "Y-axis (select one or more)",
                        df_core.columns.drop(x_axis),
                        default=[df_core.columns[1]] if x_axis != df_core.columns[1] else df_core.columns[2:3],
                        key=f"multi_y_core_{core_file}"
                    )
                else:
                    y_axes = [st.selectbox("Y-axis", df_core.columns.drop(x_axis), key=f"single_y_core_{core_file}")]

                valid_y_axes = [y for y in y_axes if y != x_axis]

                if valid_y_axes:
                    df_plot = df_core.set_index(x_axis)[valid_y_axes]
                    st.scatter_chart(df_plot)
                else:
                    st.info("Please select a Y-axis different from the X-axis to plot.")

        # ----- Profile files -----
        with inner_tabs[2]:
            profile_column_names = [
                "P (MPa)", "T (K)", "r (m)", "phase ID", "rho (kg/m3)",
                "Cp (J/kg/K)", "alpha (1/K)", "g (m/s2)", "phi (void/solid frac)",
                "sigma (S/m)", "k (W/m/K)", "VP (km/s)", "VS (km/s)", "QS",
                "KS (GPa)", "GS (GPa)", "Ppore (MPa)", "rhoMatrix (kg/m3)",
                "rhoPore (kg/m3)", "MLayer (kg)", "VLayer (m3)", "Htidal (W/m3)", "eta (Pa s)"
            ]

            for profile_file in files_by_type["Profile"]:
                file_path = os.path.join(planet_folder, profile_file)
                with open(file_path, "r", encoding="utf-8") as f:
                    lines = f.readlines()

                # Header length varies with body/config — the file declares
                # it ("nHeadLines = N" near the top; the Nth line is the
                # column-name row, data follows). The old hardcoded 81/83
                # slice broke on any profile whose header differed.
                n_head = None
                for ln in lines[:10]:
                    if 'nHeadLines' in ln:
                        try:
                            n_head = int(ln.split('=')[1])
                        except (IndexError, ValueError):
                            pass
                        break
                if n_head is None or not (0 < n_head < len(lines)):
                    n_head = 83  # legacy fallback
                header_info = "".join(lines[:n_head]).strip()
                data_str = "".join(lines[n_head:])
                df_profile = pd.read_csv(StringIO(data_str), sep=r'\s+', header=None, na_values=['nan'])

                st.subheader(f"Profile File: {profile_file}")
                with st.expander("Show Extended Profile Header"):
                    st.code(header_info)

                if df_profile.shape[1] == len(profile_column_names):
                    df_profile.columns = profile_column_names
                else:
                    st.warning(
                        f"⚠️ Column mismatch: Expected {len(profile_column_names)} headers vs {df_profile.shape[1]} data columns"
                    )
                    df_profile.columns = [f"Col_{i+1}" for i in range(df_profile.shape[1])]

                st.dataframe(df_profile, use_container_width=True)
                st.subheader("Plot Profile Data")

                chart_placeholder = st.empty()

                persist_plot = st.checkbox(
                    "Plot multiple Y‑axes against same X‑axis", key=f"persist_profile_{profile_file}"
                )

                x_axis = st.selectbox("X-axis", df_profile.columns, key=f"x_profile_{profile_file}")

                available_y_columns = [col for col in df_profile.columns if col != x_axis]

                if persist_plot:
                    y_axes = st.multiselect(
                        "Y-axis (select one or more)",
                        options=available_y_columns,
                        default=available_y_columns[:1],
                        key=f"multi_y_profile_{profile_file}"
                    )
                else:
                    y_axes = [st.selectbox(
                        "Y-axis",
                        options=available_y_columns,
                        key=f"single_y_profile_{profile_file}"
                    )]

                if y_axes:
                    try:
                        df_profile[x_axis] = pd.to_numeric(df_profile[x_axis], errors='coerce')
                        for col in y_axes:
                            df_profile[col] = pd.to_numeric(df_profile[col], errors='coerce')

                        chart_data = df_profile[[x_axis] + y_axes].dropna()
                        chart_data = chart_data.set_index(x_axis)[y_axes]

                        if not chart_data.empty:
                            chart_placeholder.scatter_chart(chart_data)
                        else:
                            chart_placeholder.empty()
                            st.info("No valid data to plot after cleaning.")
                    except Exception as e:
                        chart_placeholder.empty()
                        st.warning(f"⚠️ Unable to plot: {e}")
                else:
                    chart_placeholder.empty()
                    st.info("Please select a Y-axis different from the X-axis to plot.")
