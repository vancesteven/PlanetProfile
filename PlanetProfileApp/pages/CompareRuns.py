import streamlit as st
from pdf2image import convert_from_path
import os
import sys
import pandas as pd
import json
from datetime import datetime
from collections import defaultdict

# ----- Page setup -----
st.set_page_config(
    page_title="Compare Runs",
    page_icon="⚖️",
    layout="wide"
)

# Import utilities
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Utilities.planet_sidebar import show_planet_status
from Utilities.app_helpers import create_metric_card, format_timestamp
from Utilities.help_system import create_help_button

show_planet_status()

# Header
col_title, col_help = st.columns([5, 1])
with col_title:
    st.title("⚖️ Compare Simulation Runs")
with col_help:
    create_help_button("compare", "Compare multiple simulation runs side-by-side")

st.info("💡 Select two or more runs to compare their configurations and outputs")
st.markdown("---")

# ----- Get Planet and folder paths -----
Planet = st.session_state.get("Planet", None)
if not Planet:
    st.error("Please Select a Planet on the Planet Profile Main Settings Page")
    st.stop()

chosen_planet = st.session_state.get("ChosenPlanet", None)

# Get directories
RunPlanetProfile_directory = os.path.dirname(os.path.abspath(__file__))
app_directory = os.path.dirname(RunPlanetProfile_directory)
parent_directory = os.path.dirname(app_directory)

if parent_directory not in sys.path:
    sys.path.append(parent_directory)

figures_folder = os.path.join(parent_directory, chosen_planet, "figures")
planet_folder = os.path.join(parent_directory, chosen_planet)

# ----- Load available runs -----
pdf_file_paths = []
for filename in os.listdir(figures_folder):
    if filename.endswith(".pdf"):
        path = os.path.join(figures_folder, filename)
        mod_time = os.path.getmtime(path)
        pdf_file_paths.append((path, mod_time))

pdf_file_paths.sort(key=lambda x: x[1])

# Group into batches (runs)
BATCH_TIME_THRESHOLD_SECONDS = 60
batches = []
current_batch = []
prev_time = None

for path, mod_time in pdf_file_paths:
    if prev_time is None or (mod_time - prev_time) <= BATCH_TIME_THRESHOLD_SECONDS:
        current_batch.append((path, mod_time))
    else:
        batches.append(current_batch)
        current_batch = [(path, mod_time)]
    prev_time = mod_time

if current_batch:
    batches.append(current_batch)

# Load run metadata
metadata_file = os.path.join(planet_folder, ".run_metadata.json")
run_metadata = {}
if os.path.exists(metadata_file):
    try:
        with open(metadata_file, 'r') as f:
            run_metadata = json.load(f)
    except:
        run_metadata = {}

# Create run labels
run_labels = []
for i, batch in enumerate(batches):
    batch_time = batch[0][1]
    timestamp = datetime.fromtimestamp(batch_time)
    time_key = str(batch_time)

    if time_key in run_metadata and 'name' in run_metadata[time_key]:
        label = f"{run_metadata[time_key]['name']} ({timestamp.strftime('%Y-%m-%d %H:%M')})"
    else:
        label = f"Run {i+1}: {timestamp.strftime('%Y-%m-%d %H:%M:%S')}"

    run_labels.append((label, i, batch))

# ----- Run Selection -----
st.subheader("Select Runs to Compare")

if len(batches) < 2:
    st.warning("⚠️ You need at least 2 simulation runs to compare. Please run PlanetProfile multiple times with different configurations.")
    st.stop()

selected_runs = st.multiselect(
    "Choose runs to compare (select 2-4 for best results):",
    options=[label for label, _, _ in run_labels],
    default=[run_labels[-2][0], run_labels[-1][0]] if len(run_labels) >= 2 else [],
    help="Select multiple runs to compare side-by-side"
)

if not selected_runs:
    st.info("👆 Select at least 2 runs above to begin comparison")
    st.stop()

# Get selected batch data
selected_batch_data = [
    (label, idx, batch) for label, idx, batch in run_labels if label in selected_runs
]

st.markdown("---")

# ----- Comparison View -----
st.subheader("📊 Comparison Results")

# Parse figure types from selected runs
figure_types = {
    "Gravity": "Gravity and Pressure",
    "Wedge": "Interior Wedge Diagram",
    "Hydrosphere": "Hydrosphere Properties",
    "CoreMantTrade": "Silicate-Core Size Tradeoff",
    "MantleDens": "Silicate Radius-Density Handoff",
    "Seismic": "Seismic Properties",
    "Viscosity": "Viscosity",
    "Porosity2axes": "Porosity with 2 axes",
    "Porosity": "Porosity"
}

# Organize figures by type for each run
run_figures = defaultdict(lambda: defaultdict(list))

for label, idx, batch in selected_batch_data:
    for filepath, mod_time in batch:
        filename = os.path.basename(filepath)
        name_only = filename[:-4]

        matched_keyword = None
        for keyword in figure_types:
            if keyword in name_only:
                matched_keyword = keyword
                break

        fig_label = figure_types.get(matched_keyword, matched_keyword if matched_keyword else name_only)
        run_figures[label][fig_label].append(filepath)

# Get all unique figure types across selected runs
all_fig_types = set()
for run_dict in run_figures.values():
    all_fig_types.update(run_dict.keys())

all_fig_types = sorted(all_fig_types)

# Comparison mode selector
comparison_mode = st.radio(
    "Comparison Layout:",
    ["📊 Side-by-Side", "🔄 Sequential", "📋 Summary Table"],
    horizontal=True,
    help="Choose how to display the comparison"
)

st.markdown("---")

if comparison_mode == "📊 Side-by-Side":
    # Show figures side-by-side
    for fig_type in all_fig_types:
        st.subheader(f"📈 {fig_type}")

        # Check which runs have this figure type
        available_runs = [label for label in selected_runs if fig_type in run_figures[label]]

        if len(available_runs) < 2:
            st.info(f"ℹ️ This figure type is only available in {len(available_runs)} of the selected runs")
            continue

        # Create columns for side-by-side comparison
        cols = st.columns(len(available_runs))

        for col, run_label in zip(cols, available_runs):
            with col:
                st.markdown(f"**{run_label.split(':')[0]}**")
                if run_figures[run_label][fig_type]:
                    filepath = run_figures[run_label][fig_type][0]
                    try:
                        with st.spinner("Rendering..."):
                            images = convert_from_path(filepath)
                            st.image(images[0], use_container_width=True)
                    except Exception as e:
                        st.error(f"Error loading figure: {e}")

        st.markdown("---")

elif comparison_mode == "🔄 Sequential":
    # Show figures one after another with tabs
    fig_tabs = st.tabs(all_fig_types)

    for fig_tab, fig_type in zip(fig_tabs, all_fig_types):
        with fig_tab:
            st.subheader(f"{fig_type}")

            for run_label in selected_runs:
                if fig_type in run_figures[run_label] and run_figures[run_label][fig_type]:
                    with st.expander(f"📊 {run_label}", expanded=True):
                        filepath = run_figures[run_label][fig_type][0]
                        try:
                            with st.spinner("Rendering..."):
                                images = convert_from_path(filepath)
                                st.image(images[0], use_container_width=True)
                        except Exception as e:
                            st.error(f"Error loading figure: {e}")

elif comparison_mode == "📋 Summary Table":
    # Show availability table
    st.subheader("Figure Availability Matrix")

    availability_data = []
    for fig_type in all_fig_types:
        row = {"Figure Type": fig_type}
        for run_label in selected_runs:
            has_figure = fig_type in run_figures[run_label] and len(run_figures[run_label][fig_type]) > 0
            row[run_label.split(':')[0]] = "✅" if has_figure else "❌"
        availability_data.append(row)

    df = pd.DataFrame(availability_data)
    st.dataframe(df, use_container_width=True)

    st.markdown("---")
    st.info("💡 Switch to 'Side-by-Side' or 'Sequential' view to see the actual figures")

# ----- Configuration Comparison -----
st.markdown("---")
st.subheader("⚙️ Configuration Differences")

with st.expander("Compare Configuration Settings", expanded=False):
    st.info("Configuration comparison coming soon! This will show differences in bulk properties, ocean settings, and layer steps between runs.")

    # Placeholder for future implementation
    st.write("**Planned features:**")
    st.markdown("""
    - Compare bulk planetary properties (mass, radius, temperature, etc.)
    - Compare ocean composition and concentration
    - Compare layer step settings
    - Highlight differences between runs
    - Export comparison report
    """)
