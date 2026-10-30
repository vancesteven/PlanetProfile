"""Serving-mode detection shared by app pages.

Public serving (Streamlit Community Cloud) sets PP_PUBLIC_MODE=1 via
st.secrets or the environment; heavy compute (MCMC sampling, exploreogram
grids) is hidden and amortized settings are capped. Local runs never set
the flag, so behavior is unchanged off-cloud.
"""
import os

import streamlit as st


def public_mode() -> bool:
    if os.environ.get('PP_PUBLIC_MODE', '') == '1':
        return True
    try:
        return str(st.secrets.get('PP_PUBLIC_MODE', '')) == '1'
    except Exception:
        return False
