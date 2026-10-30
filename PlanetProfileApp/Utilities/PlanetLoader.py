import os
import importlib.util
import streamlit as st


def load_planet_module(parent_dir, planet_name, custom_module_path=None):
    """
    Loads either the default or a custom planet module.
    Default: PlanetProfile/PlanetProfile/Default/<planet>/PP<planet>.py
    Custom: PlanetProfile/Planet/<planet>/PP<planet>Custom.py
    """

    if planet_name in (None, "-- Select a Planet --"):
        st.warning("No planet selected.")
        st.stop()

    if custom_module_path:
        module_path = custom_module_path
        module_name = os.path.splitext(os.path.basename(custom_module_path))[0]
    else:
        module_path = os.path.join(parent_dir, "PlanetProfile", "Default", planet_name, f"PP{planet_name}.py")
        module_name = f"PP{planet_name}"

    if not os.path.exists(module_path):
        st.error(f"Could not find file at: {module_path}")
        st.stop()

    # Avoid reload if already loaded
    if (
        "planet_data" in st.session_state and
        st.session_state.get("planet_loaded") == module_name
    ):
        return st.session_state["planet_data"]

    # Dynamically import the module
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    planet_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(planet_module)

    st.session_state["planet_data"] = planet_module
    st.session_state["planet_loaded"] = module_name

    return planet_module
