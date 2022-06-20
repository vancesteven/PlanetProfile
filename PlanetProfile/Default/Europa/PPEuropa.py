"""
PPEuropa
Contains all body-specific parameters and information for PlanetProfile models of this body.
Import as a module and access information assigned to the attributes of the Planet struct.
Note that this file expects to be imported from the directory above.
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Europa')

# Reduce search range for melting pressure to values realistic for Europa
Planet.PfreezeUpper_MPa = 150

""" Bulk planetary settings """
Planet.Bulk.R_m = 1561.0e3
Planet.Bulk.M_kg = 4.7991e22
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.346  # Value from Anderson et al. (1998): https://doi.org/10.1126/science.281.5385.2019
Planet.Bulk.Cuncertainty = 0.005
Planet.Bulk.Tb_K = 268.305  # 30 km ice with 1.0x Seawater

""" Layer step settings """
Planet.Steps.nIceI = 200
Planet.Steps.nSilMax = 300
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'Seawater'
Planet.Ocean.wOcean_ppt = Constants.stdSeawater_ppt
Planet.Ocean.deltaP = 1.0
Planet.Ocean.deltaT = 0.1
Planet.Ocean.PHydroMax_MPa = 250.0

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 5.33e-12  # Estimate from Hussmann and Spohn (2004): https://doi.org/10.1016/j.icarus.2004.05.020, calculated from total heating rate, silicate density, and layer radii
Planet.Sil.Htidal_Wm3 = 1e-10  # Approximate max. tidal heating in silicates as modeled by Tobie et al. (2003): https://doi.org/10.1029/2003JE002099
# Rock porosity
Planet.Do.POROUS_ROCK = False
# Mantle equation of state model
#Planet.Sil.rhoSilWithCore_kgm3 = 3539.0  # This is the 1-bar, 275 K value from CV3hy1wt_678_1.tab
#Planet.Do.CONSTANT_INNER_DENSITY = True
#Planet.Sil.rhoSilWithCore_kgm3 = 3300.0
# Planet.Sil.mantleEOS = 'Simple_CI_HS_green_PP.tab'  # CI chondrite material minus Fe core, computed with Green et al. 2016 (JMG) solution model and Lodders and Fegley 1998
# Planet.Sil.rhoSilWithCore_kgm3 = 2975
# Planet.Sil.mantleEOS = 'Simple_CM_HS_green_PP.tab'  # CM chondrite material minus Fe core, computed with Green et al. 2016 (JMG) solution model and Lodders and Fegley 1998
# Planet.Sil.rhoSilWithCore_kgm3 = 2975
# Planet.Sil.mantleEOS = 'Simple_CV_HS_green_PP.tab'  # CV chondrite material minus Fe core, computed with Green et al. 2016 (JMG) solution model and Lodders and Fegley 1998
# Planet.Sil.rhoSilWithCore_kgm3 = 2975
Planet.Do.Fe_CORE = True
if Planet.Do.Fe_CORE:
    Planet.Sil.mantleEOS = 'CM_hydrous_differentiated_Ganymede_Core080Fe020S_excluding_fluid_properties.tab'  # (2900 for Q= 100 GW, 3240 for Q= 220 GW)
else:
    Planet.Sil.mantleEOS = 'CM_undifferentiated_hhph_DEW17_nofluid_nomelt_685.tab'

""" Core assumptions """
Planet.Core.rhoFe_kgm3 = 8000.0
Planet.Core.rhoFeS_kgm3 = 5150.0
Planet.Core.rhoPoFeFCC = 5455.0
Planet.Core.QScore = 1e4
Planet.Core.coreEOS = 'Fe-S_3D_EOS.mat'
Planet.Core.wFe_ppt = 800

Planet.Core.xFeSmeteoritic = 0.0405
#Planet.Core.xFeS = 0.55
Planet.Core.xFeS = 0.882  # Matching Hamish's 5377 kg/m^3 for constant inner density
Planet.Core.xFeCore = 0.0279
Planet.Core.xH2O = 0.0035

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Bulk.J2 = 435.5e-6  # J2 and C22 values from Anderson et al. (1998): https://doi.org/10.1126/science.281.5385.2019
Planet.Bulk.C22 = 131.0e-6
Planet.Magnetic.ionosBounds_m = 100e3
Planet.Magnetic.sigmaIonosPedersen_Sm = 30/100e3

# The block below should be made into one single function that returns the FTdata struct if the file is found, and warns the user/downloads if not.
# try:
#     load(['FTdata' Planet.name],'FTdata')
# catch:
#     dlNow = FTdataMissingWarning()
#     if dlNow; load(['FTdata' Planet.name],'FTdata'); end

""" Interior constraints imposed in Vance et al. 2014 """
# Planet.Sil.mSi = 28.0855
# Planet.Sil.mS = 32.065
# Planet.Sil.mFe = 55.845
# Planet.Sil.mMg = 24.305
# Planet.Sil.xOl = 0.44  # Percentage of olivine - Javoy (1995) - Fe/Si = 0.909 Mg/Si = 0.531, Mg# = 0.369
# Planet.Sil.xSi = (Planet.Sil.xOl+2*(1-Planet.Sil.xOl))*Planet.Sil.mSi/(Planet.Sil.xOl*184.295+(1-Planet.Sil.xOl)*244.3805) # mass fraction of sulfur in silicates
# Planet.Sil.MEarth_kg = 5.97e24
# Planet.Sil.xSiEarth = 0.1923  # Javoy in kg/kg in Exoplanets paper20052006-xSiSolar only in mole
# Planet.Sil.xK = 1.0  # Enrichment in K
