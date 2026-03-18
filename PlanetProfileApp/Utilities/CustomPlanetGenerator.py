import os
import streamlit as st

def generate_custom_pp_template(parent_directory):
    """
    Creates a PPCustom.py file in the Custom folder with default placeholder values.
    Compatible with PlanetLoader.
    """
    custom_dir = os.path.join(parent_directory, "PlanetProfile", "Default", "Custom")
    os.makedirs(custom_dir, exist_ok=True)

    custom_file = os.path.join(custom_dir, "PPCustom.py")

    if os.path.exists(custom_file):
        return  # Do not overwrite if it already exists

    content = '''\
"""
PPCustom
Auto-generated placeholder for a fully custom planetary body.
"""

import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Custom')

""" Bulk planetary settings """
Planet.Bulk.R_m = 0
Planet.Bulk.M_kg = 0
Planet.Bulk.Tsurf_K = 0
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0
Planet.Bulk.Cuncertainty = 0
Planet.Bulk.Tb_K = 0

""" Layer step settings """
Planet.Steps.nIceI = 100
Planet.Steps.nSilMax = 100
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Ocean assumptions """
Planet.Ocean.comp = 'Seawater'
Planet.Ocean.wOcean_ppt = Constants.stdSeawater_ppt
Planet.Ocean.deltaP = 1.0
Planet.Ocean.deltaT = 0.1
Planet.Ocean.PHydroMax_MPa = 100.0

""" Mantle assumptions """
Planet.Sil.Qrad_Wkg = 1e-12
Planet.Sil.Htidal_Wm3 = 1e-10
Planet.Do.Fe_CORE = False
Planet.Sil.mantleEOS = None

""" Core assumptions """
Planet.Core.rhoFe_kgm3 = 8000.0
Planet.Core.QScore = 1e4
Planet.Core.coreEOS = None
Planet.Core.wFe_ppt = 800

""" Seismic properties """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic field """
Planet.Bulk.J2 = 0.0
Planet.Bulk.C22 = 0.0
Planet.Magnetic.ionosBounds_m = 100e3
Planet.Magnetic.sigmaIonosPedersen_Sm = 0.1
'''

    with open(custom_file, 'w') as f:
        f.write(content)

    st.info("Generated PPCustom.py with default structure.")
