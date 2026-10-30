import streamlit as st

"""This is just to help the Planet object properly load into each page"""
def get_planet():
    planet = st.session_state.get("Planet", None)

    if not planet or planet == "-- Select a Planet --":
        st.error("No valid planet selected. Please go to the main page to choose one.")
        st.stop()

    return planet