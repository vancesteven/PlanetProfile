import streamlit as st
import os
import importlib
import sys
import numpy as np
from functools import partial
from Utilities.planet_sidebar import show_planet_status
from Utilities.app_helpers import (
    validate_radius, validate_mass, suggest_mass_from_radius,
    create_metric_card, show_parameter_info, create_progress_indicator,
    show_validation_message, format_scientific
)
from Utilities.help_system import create_help_button, get_help

# ----- Streamlit Page Setup -----
st.set_page_config(
    page_title="Bulk Planetary Settings",
    page_icon="🪙",
    layout="wide"
)
show_planet_status()

# Header with help
col_title, col_help = st.columns([5, 1])
with col_title:
    st.title("Bulk Planetary Settings")
with col_help:
    create_help_button("bulk_settings", "Configure fundamental properties of your planetary body")

st.write("Default values for your selected body are displayed below.")
st.write("You may also change the Bulk Planetary Settings to custom values")
st.info("💡 Tip: Hover over each input for helpful guidance on typical value ranges")

# Progress indicator
create_progress_indicator(
    current_step=2,
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

st.markdown("---")


# Get the planet object from the session state
Planet = st.session_state.get("Planet", None)

if Planet in (None, "-- Select a Planet --"):
    st.error("Please Select a Planet on the Planet Profile Main Settings Page")
    st.stop()


# Pulls default values into the app for each PPPlanet, if the user sets a new
# value then it is saved as a session state variable

# planet_defaults dictionary
#keeping track of defaults for the planet with this dicitonary {"key": value}
#to get to the values, call the key form the dicitonary (planet_defaults["Planet.Bulk.R_m"] returns the value of planet_module.Planet.Bulk.R_m
planet_bulk_defaults = {
    "Planet.Bulk.R_m": Planet.Bulk.R_m,
    "Planet.Bulk.M_kg": Planet.Bulk.M_kg,
    "Planet.Bulk.Tsurf_K": Planet.Bulk.Tsurf_K,
    "Planet.Bulk.Psurf_MPa": Planet.Bulk.Psurf_MPa,
    "Planet.Bulk.Cmeasured": Planet.Bulk.Cmeasured,
    "Planet.Bulk.Cuncertainty": Planet.Bulk.Cuncertainty,
    #"Planet.Bulk.Tb_K": planet_module.Planet.Bulk.Tb_K, THIS IS INSTEAD SET ON THE MAIN SETTINGS PAGE
}


# Bulk Settings dictionary
# e.g. "Planet.Bulk.R_m": ("Radius of the body (m)", 6.378e6),
# the keys are the same as in planet_bulk_defaults ("Planet.Bulk.R_m"), the 'value' is (label, value) with label
# being the string listed (e.g. "Radius of the body (m)") and value being the default value from planet_bulk_defaults (e.g. planet_module.Planet.Bulk.R_m)
bulk_settings = {
    key: (label, planet_bulk_defaults[key]) for key, label in {
        "Planet.Bulk.R_m": "Radius of the body (m)",
        "Planet.Bulk.M_kg": "Mass of the body (kg)",
        "Planet.Bulk.Tsurf_K": "Temperature at the surface ($^\circ K$)",
        "Planet.Bulk.Psurf_MPa": "Pressure at the surface (MPa)",
        "Planet.Bulk.Cmeasured": "Normalized Axial Moment of Inertia $C$",
        "Planet.Bulk.Cuncertainty": "Uncertainty in $C$",
        #"Planet.Bulk.Tb_K": "Temperature at the bottom ($^\circ K$) -
    }.items() #The .items() method is used on the dictionary, which returns an iterable of key-value pairs like: ("Planet.Bulk.R_m", "Radius of the body (m)"),
}

# Initialize session state for all bulk_settings -> the key is the first part in the dictionary, this initialization
# loop prevents having to initialize every variable individually as below

#Initializing the changed_inputs to keep track of what variables the user has changed
if "changed_bulk_settings_flags" not in st.session_state:
    st.session_state["changed_bulk_settings_flags"] = {}  #initializing blank list for changed inputs to go into later
if "changed_bulk_settings" not in st.session_state:
    st.session_state["changed_bulk_settings"] = {}  # key: value

#this block runs every time the page loads
for key, (_, default_val) in bulk_settings.items(): #initializes variables into session_state
    if key not in st.session_state: #this means only not already created keys will be added to session state
        st.session_state[key] = default_val #now, all bulk settings are in the session state


#initializing the reset_bulk_flag in the session state as False
if "reset_bulk_flag" not in st.session_state:
    st.session_state["reset_bulk_flag"] = False


# This block is only executed when the user clicks the “Reset” button.
if st.session_state["reset_bulk_flag"]: #if flag is true (if user presses reset button)
    for key, val in planet_bulk_defaults.items():
        st.session_state[key] = val  # reloads all of the planet_defaults into session_state
    st.session_state["changed_bulk_settings_flags"] = {} #clears the changed_bulk_settings_flags dictionary
    st.session_state["changed_bulk_settings"] = {} #clears the changed_bulk_settings dictionary
    st.session_state["reset_bulk_flag"] = False #reset_bulk_flag now is set to false
    st.rerun()  # ensures Streamlit restarts before widgets render



# This function is used to keep track of what settings the user has changed, so that the
# code can print out what settings have been changed
def on_change_bulk_setting(bulk_setting_key):
    st.session_state["changed_bulk_settings_flags"][bulk_setting_key] = True
    st.session_state["changed_bulk_settings"][bulk_setting_key] = st.session_state[bulk_setting_key]
    #if a user changes an input,
    # This looks up the "changed_bulk_settings" dictionary inside st.session_state and
    # adds a new entry to this dictionary with the key bulk_setting_key (the name of the setting that was changed),
    # and flags it as a setting that has been changed --> this flag will be used on the RunPlanetProfile page to actually update the values of the semi-custom planet
    # The second line then keeps track of the values that the user has selected for those settings



# Create tabs for organized input
input_tabs = st.tabs(["📏 Physical Properties", "🌡️ Boundary Conditions", "🎯 Moment of Inertia"])

with input_tabs[0]:
    st.subheader("Physical Properties")

    # Radius and Mass in columns with validation
    col1, col2 = st.columns(2)

    with col1:
        radius_key = "Planet.Bulk.R_m"
        label, _ = bulk_settings[radius_key]

        help_text = get_help("radius")
        st.number_input(
            f"{label}",
            key=radius_key,
            format="%.2e",
            help=help_text if help_text else "Radius of the planetary body in meters",
            on_change=partial(on_change_bulk_setting, radius_key)
        )

        # Validate radius
        radius_val = st.session_state[radius_key]
        is_valid, msg, severity = validate_radius(radius_val)
        show_validation_message(is_valid, msg, severity)

        # Show formatted value
        if radius_val:
            st.caption(f"= {format_scientific(radius_val)} m = {radius_val/1000:.1f} km")

    with col2:
        mass_key = "Planet.Bulk.M_kg"
        label, _ = bulk_settings[mass_key]

        help_text = get_help("mass")
        st.number_input(
            f"{label}",
            key=mass_key,
            format="%.2e",
            help=help_text if help_text else "Mass of the planetary body in kilograms",
            on_change=partial(on_change_bulk_setting, mass_key)
        )

        # Validate mass
        mass_val = st.session_state[mass_key]
        is_valid, msg, severity = validate_mass(mass_val)
        show_validation_message(is_valid, msg, severity)

        # Show formatted value
        if mass_val:
            st.caption(f"= {format_scientific(mass_val)} kg")

            # Show calculated mean density if radius is also provided
            if radius_val and radius_val > 0:
                volume = (4/3) * 3.14159 * radius_val**3
                density = mass_val / volume
                st.caption(f"Mean density: {density:.0f} kg/m³")

with input_tabs[1]:
    st.subheader("Boundary Conditions")

    col1, col2 = st.columns(2)

    with col1:
        # Surface Temperature
        tsurf_key = "Planet.Bulk.Tsurf_K"
        label, _ = bulk_settings[tsurf_key]
        help_text = get_help("surface_temperature")

        st.number_input(
            label,
            key=tsurf_key,
            min_value=0.0,
            max_value=500.0,
            step=1.0,
            help=help_text if help_text else "Temperature at the surface in Kelvin",
            on_change=partial(on_change_bulk_setting, tsurf_key)
        )

        # Show in Celsius too
        if st.session_state[tsurf_key]:
            celsius = st.session_state[tsurf_key] - 273.15
            st.caption(f"= {celsius:.1f} °C")

    with col2:
        # Surface Pressure
        psurf_key = "Planet.Bulk.Psurf_MPa"
        label, _ = bulk_settings[psurf_key]
        help_text = get_help("surface_pressure")

        st.number_input(
            label,
            key=psurf_key,
            min_value=0.0,
            format="%.6f",
            help=help_text if help_text else "Pressure at the surface in MPa",
            on_change=partial(on_change_bulk_setting, psurf_key)
        )

        # Show in bar/atm
        if st.session_state[psurf_key]:
            bar = st.session_state[psurf_key] * 10
            atm = st.session_state[psurf_key] / 0.101325
            st.caption(f"= {bar:.4f} bar = {atm:.4f} atm")

with input_tabs[2]:
    st.subheader("Moment of Inertia")

    col1, col2 = st.columns(2)

    with col1:
        # Cmeasured
        cmeas_key = "Planet.Bulk.Cmeasured"
        label, _ = bulk_settings[cmeas_key]
        help_text = get_help("moi")

        st.number_input(
            label,
            key=cmeas_key,
            min_value=0.0,
            max_value=0.5,
            format="%.6f",
            help=help_text if help_text else "Normalized axial moment of inertia",
            on_change=partial(on_change_bulk_setting, cmeas_key)
        )

        # Show interpretation
        c_val = st.session_state[cmeas_key]
        if c_val:
            if c_val < 0.33:
                interpretation = "Highly differentiated (dense core)"
            elif c_val < 0.37:
                interpretation = "Moderately differentiated"
            elif c_val < 0.4:
                interpretation = "Weakly differentiated"
            else:
                interpretation = "Nearly homogeneous"
            st.caption(f"Interpretation: {interpretation}")

    with col2:
        # Cuncertainty
        cunc_key = "Planet.Bulk.Cuncertainty"
        label, _ = bulk_settings[cunc_key]

        st.number_input(
            label,
            key=cunc_key,
            min_value=0.0,
            max_value=0.1,
            format="%.6f",
            help="Uncertainty in the moment of inertia measurement",
            on_change=partial(on_change_bulk_setting, cunc_key)
        )

        # Show as percentage
        if st.session_state[cunc_key] and st.session_state.get(cmeas_key):
            pct = st.session_state[cunc_key] / st.session_state[cmeas_key] * 100
            st.caption(f"= {pct:.2f}% of C value")


st.markdown("---")

# Summary of changes
st.subheader("Changes Summary")
changed_count = len([k for k, v in st.session_state.get("changed_bulk_settings_flags", {}).items() if v])

if changed_count > 0:
    st.warning(f"⚠️ You have modified {changed_count} setting(s) from defaults")

    with st.expander("View Changed Settings"):
        for key, was_changed in st.session_state.get("changed_bulk_settings_flags", {}).items():
            if was_changed:
                setting_name = key.split(".")[-1]
                new_val = st.session_state.get("changed_bulk_settings", {}).get(key, "N/A")
                default_val = planet_bulk_defaults.get(key, "N/A")
                st.markdown(f"- **{setting_name}**: `{default_val}` → `{new_val}`")
else:
    st.success("✅ All settings are at default values")

# Actual reset button widget at the bottom of the page
col1, col2, col3 = st.columns([2, 1, 2])
with col2:
    if st.button("🔄 Reset to Defaults", type="secondary", use_container_width=True):
        st.session_state["reset_bulk_flag"] = True
        # which triggers the if st.session_state["reset_bulk_flag"] function above
