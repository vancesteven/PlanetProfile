"""
PPEuropa
Contains all body-specific parameters and information for PlanetProfile models of this body.
Import as a module and access information assigned to the attributes of the Planet and Params structs.
Note that this file expects to be imported from the directory above.
"""
import numpy as np
from Utilities.dataStructs import PlanetStruct, ParamsStruct, Constants

from config import Params

Planet = PlanetStruct('Europa')

# Reduce search range for melting pressure to values realistic for Europa
Planet.PfreezeUpper_MPa = 150

""" Bulk planetary settings """
Planet.Bulk.R_m = 1561.0e3
Planet.Bulk.M_kg = 4.7991e22
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.346
Planet.Bulk.Cuncertainty = 0.005
Planet.Bulk.Tb_K = 269.8

""" Layer step settings """
Planet.Steps.nIceI = 200
Planet.Steps.nRefRho = 30
Planet.Steps.nSilMax = 300
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'Seawater'
Planet.Ocean.wOcean_ppt = 0.0
Planet.Ocean.deltaP = 1.0
Planet.Ocean.PHydroMax_MPa = 250

""" Silicate Mantle """
Planet.Sil.kTherm_WmK = 4.0
Planet.Sil.Qrad_Wkg = 5.38e-12  # Estimate from Hussmann and Spohn (2004): https://doi.org/10.1016/j.icarus.2004.05.020
Planet.Sil.Htidal_Wm3 = 1e-10  # Approximate max. tidal heating in silicates as modeled by Tobie et al. (2003): https://doi.org/10.1029/2003JE002099
Planet.Do.EQUIL_Q = False
# Rock porosity
Planet.Do.POROUS_ROCK = False
Planet.Do.P_EFFECTIVE = False
# Mantle equation of state model
Planet.Sil.mantleEOS = 'CV3hy1wt_678_1.tab'  # (2900 for Q= 100 GW, 3240 for Q= 220 GW)
Planet.Sil.rhoSilWithCore_kgm3 = 3539.0  # This is the 1-bar, 275 K value from CV3hy1wt_678_1.tab
# Planet.Sil.mantleEOS = 'Simple_CI_HS_green_PP.tab'  # CI chondrite material minus Fe core, computed with Green et al. 2016 (JMG) solution model and Lodders and Fegley 1998
# Planet.Sil.rhoSilWithCore_kgm3 = 2975
# Planet.Sil.mantleEOS = 'Simple_CM_HS_green_PP.tab'  # CM chondrite material minus Fe core, computed with Green et al. 2016 (JMG) solution model and Lodders and Fegley 1998
# Planet.Sil.rhoSilWithCore_kgm3 = 2975
# Planet.Sil.mantleEOS = 'Simple_CV_HS_green_PP.tab'  # CV chondrite material minus Fe core, computed with Green et al. 2016 (JMG) solution model and Lodders and Fegley 1998
# Planet.Sil.rhoSilWithCore_kgm3 = 2975

""" Core assumptions """
Planet.Do.Fe_CORE = True
Planet.Core.rhoFe_kgm3 = 8000.0
Planet.Core.rhoFeS_kgm3 = 5150.0
Planet.Core.rhoPoFeFCC = 5455.0
Planet.Core.QScore = 1e4
Planet.Core.coreEOS = 'sulfur_core_partition_SE15_1pctSulfur.tab'
Planet.Core.xFeSmeteoritic = 0.0405
Planet.Core.xFeS = 0.55
Planet.Core.xFeCore = 0.0279
Planet.Core.xH2O = 0.0035

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0
# Attenuation Parameters Based on those Described in Cammarano et al. 2006
# ice I
Planet.Seismic.BIceI = 0.56
Planet.Seismic.gammaIceI = 0.2
Planet.Seismic.gIceI = 22.0  # C2006
# ice II
Planet.Seismic.BIceII = 0.56
Planet.Seismic.gammaIceII = 0.2
Planet.Seismic.gIceII = 30.0
# ice III
Planet.Seismic.BIceIII = 0.56
Planet.Seismic.gammaIceIII = 0.2
Planet.Seismic.gIceIII = 25.0
# ice V
Planet.Seismic.BIceV = 0.56
Planet.Seismic.gammaIceI = 0.2
Planet.Seismic.gIceV = 27.0
# ice VI
Planet.Seismic.BIceVI = 0.56
Planet.Seismic.gammaIceVI = 0.2
Planet.Seismic.gIceVI = 28.0
# mantle
Planet.Seismic.BSil = 0.56
Planet.Seismic.gammaSil = 0.2
Planet.Seismic.gSil = 30.0  # C2006

""" Magnetic induction """
Planet.Magnetic.peaks_Hz = np.array([4.946e-5, 2.473e-5, 3.259e-6])
Planet.Magnetic.fOrb_radps = 2*np.pi/3.55/86400
Planet.Magnetic.ionosBounds_m = 100e3
Planet.Magnetic.sigmaIonosPedersen_Sm = 30/100e3

""" Other parameter options """
Params.PLOT_SIGS = True
Params.wLims = [np.log10(0.001), np.log10(1000)]
Params.LEGEND = True
Params.LegendPosition = 'North'
Params.yLim = [910, 1230]
Params.LineStyle = Params.LS_Mg
Params.wRefLine = '--'
Params.wRef = [0, 5, 10, 15]

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
