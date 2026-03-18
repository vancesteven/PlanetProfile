import streamlit as st

# Configure page - must be first Streamlit command
st.set_page_config(
    page_title="Welcome to PlanetProfile!",
    page_icon="🌍",
    layout="wide"
)


st.title("Welcome to PlanetProfile!")

st.subheader("Getting Started")
st.write("PlanetProfile generates simulations of the interiors of planetary bodies. Navigate the tabs to the left to set up your planet.")
st.markdown("---")

st.subheader("Page Navigation")
st.markdown(
" - Use the Main Settings Page to select your planet or design a custom body.\n\n"
" - Designate properties such as the mass and radius on the Bulk Planetary Settings Page. \n\n"
" - Use the Ocean Settings page to configure your ocean by using a pre-defined ocean or setting your own salt species and concentrations. \n\n"
" - Configure the core and mantle of your planet with the  Core and Silicate Settings Page. \n\n"
" - Designate the granularity of your simulation with the Layer Step Settings page. \n\n"
" - Run your simulation on the Run PlanetProfile page \n\n"
" - View figure outputs and tables on the PlanetProfile Outputs page \n\n"
"See below for more information on Planet Profile.")
st.markdown("---")

st.subheader("Key Functionality")
st.write("Radial models of planetary interiors are generated from bulk properties based on geophysical models, lab data, and minimal assumptions")
st.markdown("---")

st.subheader("Plain Language Summary")
st.write("The software package PlanetProfile was developed in order to connect measurable properties of planetary bodies to each other and determine how planetary interiors might be structured.")
st.markdown("---")

st.subheader("Abstract")
st.write("The open-source PlanetProfile framework was developed to investigate the interior structure of icy moons based on collectively matching their observed properties and comparative planetology. \
The software relates observed and measured properties, assumptions such as the type of materials present, and laboratory equation-of-state (EOS) data through geophysical and thermodynamic models \
to evaluate radial profiles of mechanical, thermodynamic, and electrical properties, as self-consistently as possible.")


st.write("$From: Styczinski, S. D. Vance, and M. Melwani Daswani (2023)$ \n\n"
         "PlanetProfile: Self-consistent interior structure modeling for ocean worlds and rocky dwarf planets in Python. \n"
        "Earth and Space Science, 10(8), 10.1029/2022ea002748")
