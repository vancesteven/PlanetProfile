import streamlit
def show_planet_status():
    chosen = streamlit.session_state.get("ChosenPlanet", "-- Select a Planet --")
    if chosen and chosen != "-- Select a Planet --":
        streamlit.sidebar.markdown(f"**Current Planet:** {chosen}")
    else:
        streamlit.sidebar.markdown("⚠️ **Current Planet:** ❌ Not Set")
