"""
Minimal Titan-like NaCl-ocean template for test_cache_comp_threading.py.

Configured to reliably build an ocean-bearing structure at Tb ~ 243-246 K
for NaCl compositions at moderate salinities:
  - NONHYDROSTATIC = False           (avoids C/MR^2 = None failure)
  - ConstantProps['Inner'] = True    (bypasses MoI EOS search; required when
                                      POROUS_ROCK=False and Fe_CORE=False)
  - Cuncertainty = 0.05              (wide window; allows Tb sweep)
  - PfreezeUpper_MPa = 300           (NaCl freeze curve needs >230 MPa at these Tb)
  - No NO_OCEAN_EXCEPT_INNER_ICES   (allows liquid ocean to form)

This file is in PlanetProfile/Inference/tests/ so it is importable as
"PlanetProfile.Inference.tests.PPTest_NaClOcean" when the project root is on
sys.path.
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test_NaClOcean')

# Titan-like body params (non-scientific test use only)
Planet.parent = 'Saturn'
Planet.Bulk.R_m = 2574.73e3
Planet.Bulk.M_kg = 1.3452e23
Planet.Bulk.Tsurf_K = 94
Planet.Bulk.Psurf_MPa = 0.15
Planet.Bulk.Cmeasured = 0.341
Planet.Bulk.Cuncertainty = 0.05           # wide window for Tb sweep
Planet.Bulk.eccentricity = 0.0288
Planet.Bulk.Torb_s = 15.945 * 24 * 3600  # mean motion derived from Torb_s

Planet.Do.NONHYDROSTATIC = False
Planet.Bulk.Tb_K = 244.5                  # default; overridden per node

# NaCl EOS requires PfreezeUpper_MPa >= ~300 MPa at Tb~244 K;
# the default 230 MPa misses the Ih-liquid transition.
Planet.PfreezeUpper_MPa = 300

""" Layer step settings """
Planet.Steps.nIceI = 100            # reduced for speed
Planet.Steps.nSilMax = 100
Planet.Steps.nCore = 5
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere """
Planet.Ocean.comp = 'NaCl'
Planet.Ocean.wOcean_ppt = 100
Planet.Ocean.deltaP = 8.0
Planet.Ocean.deltaT = 0.5
Planet.Ocean.phaseType = 'lookup'
Planet.Ocean.PHydroMax_MPa = 1800.0
Planet.Ocean.THydroMax_K = 350.0

# No HP-ice underplating flags set — let PP find the ocean naturally.

""" Silicate Mantle — ConstantProps path for fast MoI matching """
Planet.Do.POROUS_ROCK = False
Planet.Do.Fe_CORE = False
# ConstantProps['Inner'] = True bypasses the EOS-driven MoI solver.
Planet.Do.ConstantProps['Inner'] = True
Planet.Sil.rhoSilWithCore_kgm3 = 3300.0
Planet.Sil.Qrad_Wkg = 1.5e-12
Planet.Sil.Htidal_Wm3 = 1e-10

""" Core (unused — Fe_CORE=False) """
Planet.Core.rhoFe_kgm3 = 8000.0

""" Seismic """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Bulk.J2 = 33.089e-6
Planet.Bulk.C22 = 10.385e-6
Planet.Magnetic.ionosBounds_m = [100e3, 250e3]
Planet.Magnetic.sigmaIonosPedersen_Sm = [1e-16, 80 / 150e3]
