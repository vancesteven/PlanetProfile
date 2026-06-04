"""
Amortized Inference (SBI) Tab for PlanetProfileApp

Instant Bayesian inference using pre-trained Simulation-Based Inference (SBI)
models (Normalizing Flows).

Author: PlanetProfile Team
Date: 2026-06-04
"""
import streamlit as st
import os
import sys

# Path setup
this_file = os.path.abspath(__file__)
pages_directory = os.path.dirname(this_file)
app_directory = os.path.dirname(pages_directory)
parent_directory = os.path.dirname(app_directory)

if app_directory not in sys.path:
    sys.path.insert(0, app_directory)
if parent_directory not in sys.path:
    sys.path.insert(0, parent_directory)

from Utilities.get_planet import get_planet
from Utilities.planet_sidebar import show_planet_status

def main():
    show_planet_status()
    st.title("🚀 Amortized Inference (SBI)")
    
    st.markdown("""
    **Simulation-Based Inference (SBI)** allows for near-instant posterior estimation 
    by leveraging deep learning models (Normalizing Flows) pre-trained on millions 
    of interior structure simulations.
    """)
    
    st.info("💡 **Status:** This tab is currently under development. Pre-trained models for Titan and Europa are being generated using the MCMC pipeline.")
    
    st.subheader("⚙️ Model Selection")
    model_choice = st.selectbox(
        "Select pre-trained model:",
        options=["Titan (Andrade, No-Ocean)", "Titan (Maxwell, Ocean)", "Europa (Andrade, Thin Shell)"],
        index=0
    )
    
    st.subheader("🎯 Observables")
    st.write("Enter the observed data to get an instant posterior estimate.")
    
    col1, col2 = st.columns(2)
    with col1:
        st.number_input("Re(k₂) — Real part:", value=0.608, format="%.4f")
    with col2:
        st.number_input("Im(k₂) — Imaginary part:", value=0.135, format="%.4f")
        
    st.button("⚡ Generate Instant Posterior", disabled=True)
    
    st.divider()
    
    with st.expander("📚 About SBI"):
        st.markdown("""
        **What is SBI?**
        
        Traditional MCMC requires running the forward model thousands of times for 
        every new observation. SBI "amortizes" this cost by training a neural network 
        to learn the mapping from observables back to parameters.
        
        **Advantages:**
        - **Instant:** Posteriors are generated in milliseconds.
        - **Interactive:** Perfect for "what-if" scenarios in the GUI.
        - **Robust:** Can handle complex likelihoods that are difficult for MCMC.
        """)

if __name__ == "__main__":
    main()
