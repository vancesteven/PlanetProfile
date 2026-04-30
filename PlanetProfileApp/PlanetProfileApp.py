import streamlit as st

# This sets up all of the pages of the GUI
about = st.Page("pages/About.py", title = "Welcome to PlanetProfile", icon = "🌍")
main_settings = st.Page("pages/PlanetProfileMainSettings.py", title = "PlanetProfile Main Settings", icon = "🌙" )
bulk_planetary_settings = st.Page("pages/BulkPlanetarySettings.py", title = "Bulk Planetary Settings", icon = "🪙")
ocean_settings = st.Page("pages/OceanSettings.py", title = "Ocean Settings", icon = "🌊" )
core_settings = st.Page("pages/CoreSettings.py", title = "Core and Silicate Settings", icon = "🌋" )
layer_step_settings = st.Page("pages/LayerStepSettings.py", title = "Layer Step Settings", icon = "📶")
#figure_settings = st.Page("pages/FigureSettings.py", title = "Figure Settings", icon = "📈")
run_PlanetProfile = st.Page("pages/RunPlanetProfile.py", title = "Run PlanetProfile", icon = "🚀")
PlanetProfile_outputs = st.Page("pages/PlanetProfileOutputs.py", title = "PlanetProfile Outputs", icon = "📈")
exploreogram = st.Page("pages/Exploreogram.py", title = "Exploreogram", icon = "🔬")
inference = st.Page("pages/Inference.py", title = "Inference", icon = "🧮")
# Set up navigation
pg = st.navigation([about, main_settings, bulk_planetary_settings, ocean_settings, core_settings, layer_step_settings, run_PlanetProfile, PlanetProfile_outputs, exploreogram, inference])

# Run the selected page
pg.run()
