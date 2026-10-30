"""
Temporary configuration override for PlanetProfile when run from Streamlit app
This ensures LaTeX table output is enabled
"""

def configure_for_streamlit(Params):
    """Enable settings needed for Streamlit display"""
    # Enable LaTeX table output for publication use
    Params.DISP_TABLE = True
    Params.DISP_LAYERS = True
    return Params
