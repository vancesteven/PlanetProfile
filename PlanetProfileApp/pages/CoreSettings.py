import streamlit as st
import os
import sys
import re
from Utilities.planet_sidebar import show_planet_status
from functools import partial
from PIL import Image

# ----- Path Setup -----
# Get absolute path to app directory for loading images
app_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ----- Streamlit Page Setup -----
st.set_page_config(
    page_title="Core and Silicate Settings",
    page_icon="./PPlogo.ico",
    layout="wide"
)
show_planet_status()
st.title("Core and Silicate Settings")
st.markdown("#### Configure the makeup of your planet's core and mantle here")
st.markdown("---")

# Get the planet name from the session state
Planet = st.session_state.get("Planet", None)
if not Planet:
    st.error("Please Select a Planet on the Planet Profile Main Settings Page")
    st.stop()

# Dictionary of core and silicate defaults
core_and_sil_defaults = {
    "Planet.Do.Fe_CORE" : Planet.Do.Fe_CORE,
    "Planet.Core.rhoFe_kgm3" : Planet.Core.rhoFe_kgm3,
    "Planet.Core.rhoFeS_kgm3" : Planet.Core.rhoFeS_kgm3,
    "Planet.Core.coreEOS" : Planet.Core.coreEOS,
    "Planet.Core.QScore" : getattr(Planet.Core, "QScore", None), #not every planet has this, hence the get attr
    "Planet.Core.wFe_ppt" : Planet.Core.wFe_ppt,
    "Planet.Core.xFeSmeteoritic" : Planet.Core.xFeSmeteoritic,
    "Planet.Core.xFeS" : Planet.Core.xFeS,
    "Planet.Core.xFeCore" : getattr(Planet.Core, "xFeCore", None), #not every planet has this, hence the get attr
    "Planet.Core.xH2O" : Planet.Core.xH2O,
    "Planet.Sil.Qrad_Wkg" : Planet.Sil.Qrad_Wkg,
    "Planet.Sil.mantleEOS" : Planet.Sil.mantleEOS,
    "Planet.Do.POROUS_ROCK": Planet.Do.POROUS_ROCK,
    "Planet.Sil.phiRockMax_frac" : Planet.Sil.phiRockMax_frac,
    "Planet.Sil.Pclosure_MPa" : Planet.Sil.Pclosure_MPa,
    "Planet.Sil.porosType" : Planet.Sil.porosType,
    "Planet.Sil.rhoSilWithCore_kgm3" : Planet.Sil.rhoSilWithCore_kgm3
}

# Core Settings dictionary
# e.g. "Planet.Core.rhoFe_kgm3" : ("Density of Iron (Fe), (kg/m³)", 6.378e6),
# the keys are the same as in core_and_sil_defaults (core_and_sil), the 'value' is (label, value) with label
# being the string listed (e.g. "Density of Iron (Fe), (kg/m³)") and value being the default value from core_and_sil_defaults (e.g. planet_module.Planet.Core.rhoFe_kgm3)
core_settings = {
    key: (label, core_and_sil_defaults[key]) for key, label in {
    "Planet.Do.Fe_CORE" : "Include an Iron Core?",
    "Planet.Core.rhoFe_kgm3" : "Density of Iron (Fe), (kg/m³)",
    "Planet.Core.rhoFeS_kgm3" : "Density of Iron Sulfide (FeS), (kg/m³)",
    "Planet.Core.coreEOS" : "EOS file for the core",
    "Planet.Core.QScore" : "Quality Factor Qs; describes how efficiently the core absorbs or dissipates energy",
    "Planet.Core.wFe_ppt" : "Mass fraction of Iron in the core (ppt)",
    "Planet.Core.rhoFeS_kgm3" : "Density of Iron Sulfide (FeS), (kg/m³)",
    "Planet.Core.xFeSmeteoritic" : "Molar Fraction of Iron Sulfide (FeS), meteoritic",
    "Planet.Core.xFeS" : "Molar Fraction of Iron Sulfide (FeS) in the core",
    "Planet.Core.xFeCore" : "Molar Fraction of Iron (Fe)in the core",
    "Planet.Core.xH2O" : "Molar Fraction of Water in the core"}.items()
}

#Initializing the changed_inputs to keep track of what variables the user has changed
if "changed_core_settings_flags" not in st.session_state:
    st.session_state["changed_core_settings_flags"] = {}  #initializing blank list for changed inputs to go into later
if "changed_core_settings" not in st.session_state:
    st.session_state["changed_core_settings"] = {}  # key: value

# Initialize session state for each core and silicate attribute (do this once at start)
#this block runs every time the page loads
for key, default_val in core_and_sil_defaults.items(): #initializes variables into session_state
    if key not in st.session_state: #this means only not already created keys will be added to session state
        st.session_state[key] = default_val #now, all core settings are in the session state


#initializing the reset_core_flag in the session state as False
if "reset_core_flag" not in st.session_state:
    st.session_state["reset_core_flag"] = False


# This function is used to keep track of what settings the user has changed, so that the
# code can print out what settings have been changed
def on_change_core_setting(core_setting_key):
    st.session_state["changed_core_settings_flags"][core_setting_key] = True
    st.session_state["changed_core_settings"][core_setting_key] = st.session_state[core_setting_key]

# This block is only executed when the user clicks the “Reset” button.
if st.session_state["reset_core_flag"]: #if flag is true (if user presses reset button)
    for key, val in core_and_sil_defaults.items():
        st.session_state[key] = val  # reloads all of the core_and_sil_defaults into session_state
    st.session_state["changed_core_settings_flags"] = {} #clears the changed_core_settings_flags dictionary
    st.session_state["changed_core_settings"] = {} #clears the changed_core_settings dictionary
    st.session_state["reset_core_flag"] = False #reset_core_flag now is set to false
    st.rerun()  # ensures Streamlit restarts before widgets render




# ----- Do Fe_Core toggle -----
st.subheader("Include an Iron Core?")
st.write("If you choose not to include an iron core, your planet's core will be populated entirely by silicate mantle instead.")
# Actual choosing of core/no core

user_do_core_flag = st.toggle("Does your planet have a core?", key = "Planet.Do.Fe_CORE", on_change = partial(on_change_core_setting, "Planet.Do.Fe_CORE"))

# Wedge images of core/no core from test runs 6 and 7
st.write("The following wedge diagrams show the difference of including an iron core or not on the planet interior. If an iron core is diabled, then silicates populate the entire planet core.")
# Diagram to show different wedges with/without core
# Get absolute paths to image files
app_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
no_core_wedge = Image.open(os.path.join(app_directory, 'wedgenocore.png'))
with_core_wedge = Image.open(os.path.join(app_directory, 'wedgewithcore.png'))
col1, mid, col2 = st.columns([20,1,20])
with col1:
    st.image(no_core_wedge, caption='Example Wedge without a core')
with col2:
    st.image(with_core_wedge, caption='Example Wedge with a core')

st.markdown("---")




# ----- User selects core EOS or sets densities for the core and silicate materials -----
# Get absolute path to the EOS directory relative to this file
eos_folder = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "PlanetProfile", "Thermodynamics", "EOStables", "Perple_X")
)

# Get list of mantle EOS files
eos_files = [f for f in os.listdir(eos_folder) if f.endswith(".tab")]
eos_descriptions = [f for f in os.listdir(eos_folder) if f.endswith(".txt")]

# mapping of EOS files to their readme descriptions made by Mohit
eos_mapping = {
    "CI_undifferentiated_hhph_DEW17_nofluid_nomelt_685.tab" : "CI chondrite readme.txt",
    "Comet_67P_CG_v7_excluding_fluid_properties.tab" : "Comet 67PCG readme.txt",
    "CM_hydrous_differentiated_Ganymede_Core080Fe020S_excluding_fluid_properties.tab" : "README - CM_hydrous_differentiated_Ganymede.txt",
    "CM_undifferentiated_hhph_DEW17_nofluid_nomelt_685.tab" : "CM chondrite readme.txt",
    "CV_undifferentiated_v4_687_DEW17_nofluid_nomelt_v2.tab" : "CV chondrite readme.txt"
}

# This reads the readmes into the GUI for users to see
eos_info = {}
for tab_file, txt_file in eos_mapping.items():
    txt_path = os.path.join(eos_folder, txt_file)
    if os.path.exists(txt_path):
        with open(txt_path, "r") as f:
            eos_info[tab_file] = f.read()
    else:
        eos_info[tab_file] = "Description file not found."

#This keeps track if the planet has a core_EOS set. If no core EOS is set, then we will be using densities instead
default_core_eos = getattr(Planet.Core, "coreEOS", "")
default_mode = "Use Equation of State file" if default_core_eos else "Manually enter densities"

# This loads the chosen EOS or density choice from the session state or from the default for that planet
use_eos_mode = st.session_state.get("use_core_definition", default_mode) # this adds the seletion of using EOS or manually entering desnities into session state or updates it into session state

# Gets the index for the selectbox based on current value
eos_options = ["Use Equation of State file", "Manually enter densities"]
use_eos_default_index = eos_options.index(use_eos_mode)



st.markdown("### Use Equations of State or Densities")
st.write("You can either specify the equation of state (EOS) for your core and silicates from a list of available core EOS's or input the densities of your core and silicate species.")
# Toggle between using EOS or manual densities
default_core_eos = Planet.Core.coreEOS or ""
use_eos = st.selectbox(
        "How do you want to define your core and silicate composition?",
        options = eos_options ,
        index = use_eos_default_index,
        key="use_core_definition" if "use_core_definition" in st.session_state else None,
        on_change=partial(on_change_core_setting, "use_core_definition") #this will keep track if the user changes from wanting to use EOS or use densities instead
    )
st.markdown("---")


# ----- Core settings begin here -----

if user_do_core_flag:

    if use_eos == "Use Equation of State file":
        # Single core EOS
        core_eos_file = "Fe-S_3D_EOS.mat"
        core_eos_path = os.path.join(eos_folder, core_eos_file)
        st.info("Using Fe-S_3D_EOS.mat EOS for the core...")  #currently, we only have the one model for iron core stuff so it will efault to using this
        #track_core_change("Planet.Core.coreEOS", core_eos_file, default_core_eos)  if eventually we have mulitple core models this would keep track of user changes

    else:
        # User wants to enter manual densities
        Planet.Core.coreEOS = "" #these both clear Planet.Core.coreEOS and session state in the event the user wants to input manual densities instead
        st.session_state["Planet.Core.coreEOS"] = ""  # sync session state
        Planet.Do.CONSTANT_INNER_DENSITY = True #---> this needs to be set to true if using densities


        st.markdown("#### Manually Enter Core Densities")
        for var, label, default in [
            ("Planet.Core.rhoFe_kgm3", "Density of Iron (Fe), (kg/m³)", Planet.Core.rhoFe_kgm3), # pulls default of densities from the planet
            ("Planet.Core.rhoFeS_kgm3", "Density of Iron Sulfide (FeS), (kg/m³)", Planet.Core.rhoFeS_kgm3)
            #("Planet.Core.rhoPoFeFCC", "Density of PoFeFCC (kg/m³)", Planet.Core.rhoPoFeFCC)  --> Marshall started implementing this, but so far it is never actually used for anything according to Scott
        ]:
            val = st.number_input(
            label,
            key = var,
            on_change = partial(on_change_core_setting, var))


    st.markdown("---")
    st.markdown("### Other Core Settings")

    for var, label, default in [
        ("Planet.Core.QScore", "Quality Factor Qs; describes how efficiently the core absorbs or dissipates energy", Planet.Core.QScore),
        ("Planet.Core.wFe_ppt", "Mass fraction of Iron in the core (ppt)", Planet.Core.wFe_ppt),
    ]:
        val = st.number_input(label, key = var, on_change = partial(on_change_core_setting, var))
    st.markdown("---")



    st.markdown("### Core Molar Fractions")
    for var, label, default in [
        ("Planet.Core.xFeSmeteoritic", "Molar Fraction of Iron Sulfide (FeS), meteoritic", Planet.Core.xFeSmeteoritic),
        ("Planet.Core.xFeS", "Molar Fraction of Iron Sulfide (FeS) in the core", Planet.Core.xFeS),
        ("Planet.Core.xFeCore", "Molar Fraction of Iron (Fe)in the core", Planet.Core.xFeCore),
        ("Planet.Core.xH2O", "Molar Fraction of Water in the core", Planet.Core.xH2O),
    ]:
        val = st.number_input(label, format="%.4f", key = var, on_change = partial(on_change_core_setting, var))
    st.markdown("---")


#----- Silicate Settings begin here -----

st.markdown("### Silicate EOS and Densities")
st.write("This manages the equation of state for the planet mantle")
# Defaults from Planet object
default_silicate_eos = Planet.Sil.mantleEOS # Current silicate EOS default from Planet

if use_eos == "Use Equation of State file":
    st.info("You selected to use equations of state. Choose your silicate EOS below. For more information on each silicate EOS, click 'EOS File Description' below")

    current_eos = Planet.Sil.mantleEOS
    selected_silicate_eos = st.selectbox(
        "Select Silicate Mantle EOS File",
        eos_files,
        index = eos_files.index(default_silicate_eos)
        if default_silicate_eos in eos_files else 0,
        on_change = partial(on_change_core_setting, "Planet.Sil.mantleEOS")
        )

    # Display description if available of what each particular EOS file does
    description = eos_info.get(selected_silicate_eos, "No description available.")
    with st.expander("EOS File Description"):
        st.text(description)

else:
    st.info("You selected to manually enter densities. Enter your silicate density below.")

    silicate_density = st.number_input(
            "Silicate Density with Core (kg/m³)",
            key = "Planet.Sil.rhoSilWithCore_kgm3",
            on_change = partial(on_change_core_setting, "Planet.Sil.rhoSilWithCore_kgm3")
        )
    # Clear EOS if switching to manual density
    Planet.Sil.mantleEOS = "" #these both clear Planet.Core.coreEOS and session state in the event the user wants to input manual densities instead
    st.session_state["Planet.Sil.mantleEOS"] = ""  # sync session state
st.markdown("---")


st.markdown("### Other Silicate Settings")

# Radiogenic heating
Qrad = st.number_input("Average radiogenic heating rate for silicates in W/kg.",  format="%.15E", key = "Planet.Sil.Qrad_Wkg", on_change = partial(on_change_core_setting, "Planet.Sil.Qrad_Wkg"))


# Tidal heating  ---> Scott says this isn't actually being used yet, so we will ignore it for now
#default_Htidal = Planet.Sil.Htidal_Wm3
#Htidal = st.number_input("Tidal Heating (W/m³)", value=default_Htidal, format="%.1e", on_change = partial(on_change_core_setting,"Planet.Sil.Htidal_Wm3" )


# Diagrams of wedges with/without porosity enabled
# Using test wedges 1 and 11
st.write("The following wedge diagrams show the difference of including porosity or not on the planet interior. If porosity is enabled, then water will populate the available porous holes in the silicate layers beneath the ocean (visible as a gradient beneath the ocean layer).")
# Use absolute paths for images
no_porous_wedge = Image.open(os.path.join(app_directory, 'wedgenoporous.png'))
with_porous_wedge = Image.open(os.path.join(app_directory, 'wedgewithporous.png'))
cola, mida, colb = st.columns([20,1,20])
with cola:
    st.image(no_porous_wedge, caption='Example Wedge with porosity disabled')
with colb:
    st.image(with_porous_wedge, caption='Example Wedge with porosity enabled')

# Rock porosity toggle
porous_flag = st.toggle("Enable Porous Rock - If the rock is porous and you have an ocean, water from the ocean will be able to flow into the porous areas", key = "Planet.Do.POROUS_ROCK", on_change = partial(on_change_core_setting,"Planet.Do.POROUS_ROCK"))

# If user wants porosity, then porosity settings appear
if porous_flag:
    default_porosType = Planet.Sil.porosType
    poros_options = ['Han et al. 2014', 'Vitovtova et al. 2014', 'Chen et al. 2020'] #these are the current porosity models we have available. If we expect to have many more, we can add dynamic loading here.


    phiRockMax_frac = st.number_input("Porosity (void fraction) of the rocks in vacuum. This is the expected value for core-less bodies, and porosity is modeled for a range around here to find a matching MoI. For bodies with a core, this is a fixed value for rock porosity at P=0.",
        key = "Planet.Sil.phiRockMax_frac", on_change = partial(on_change_core_setting,"Planet.Sil.phiRockMax_frac"))

    Pclosure_MPa = st.number_input("Pressure threshold in MPa beyond which pores in silicates shut completely and porosity drops to zero, for use in Han et al. (2014) model.", key = "Planet.Sil.Pclosure_MPa", on_change = partial(on_change_core_setting,"Planet.Sil.Pclosure_MPa"))
    #st.selectbox requires using indexing for the default value that shows up in the selectbox, so this loop looks for the index of the default porosity model from the list and passes it to the selectbox
    try:
        default_index = poros_options.index(default_porosType)
    except ValueError:
        default_index = 0  # fallback to first option if not found

    porosType = st.selectbox(
    "Porosity model to apply for silicates.",
    options=poros_options,
    index=default_index,
    on_change = partial(on_change_core_setting,"Planet.Sil.porosType"))

# Actual reset button widget at the bottom of the page
st.markdown("---")
if st.button("🔄 Reset to default core  (double click)"):
    st.session_state["reset_core_flag"] = True
