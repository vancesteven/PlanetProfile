"""
defineStructs: Create empty classes and subclasses for holding body-specific data
Values are typically set to None as defaults; body-specific values should be set in the PPBody.py files,
all of which override the settings below.
Optional SWITCHES are in all caps. These typically have a default value of False.

Example usage:
Planet = PlanetStruct('nameOfBody')
Planet.Bulk.R_m = 1560e3
Planet.Ocean.comp = 'MgSO4'
Planet.Sil.mantleEOS = 'CV3hy1wt_678_1.tab'
Planet.Do.Fe_CORE = False
"""

import numpy as np
import os, shutil
import cmasher
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import rgb_to_hsv, hsv_to_rgb
from scipy.interpolate import interp1d

# We have to define subclasses first in order to make them instanced to each Planet object
""" Run settings """
class BulkSubstruct():

    def __init__(self):
        self.Tb_K = None  # Temperature at the bottom of the ice I layer (ice-ocean interface when there are no ice III or V underplate layers). Ranges from 238.5 to 261.165 K for ice III transition and 251.165 to 273.16 for melting temp. This must remain set to None here as a default. Exactly two out of three of Bulk.Tb_K, Bulk.zb_km, and Ocean.wOcean_ppt must be set for every model with surface H2O.
        self.rho_kgm3 = None  # Bulk density in kg/m^3 -- note that this is intended to be derived and not set.
        self.R_m = None  # Mean body outer radius in m
        self.M_kg = None  # Total body mass in kg
        self.Tsurf_K = None  # Surface temperature in K
        self.Psurf_MPa = None  # Surface pressure in MPa
        self.Cmeasured = None  # Axial moment of inertia C/MR^2, dimensionless
        self.Cuncertainty = None  # Uncertainty (std dev) of C/MR^2 (used to constrain models via consistency within the uncertainty), dimensionless
        self.phiSurface_frac = None  # Scaling value for the ice porosity at the surface (void fraction): falls within a range of 0 and 1. 0 is completely non-porous and larger than 0.2 is rare. From Han et al. (2014)
        self.qSurf_Wm2 = None  # Heat flux leaving the planetary surface. Currently required only for clathType = 'bottom'.
        self.clathMaxThick_m = None  # (Approximate) fixed limit for thickness of clathrate layer in m. Treated as an assumed layer thickness when clathType = 'bottom' or Do.NO_ICE_CONVECTION is True, and as a maximum for 'top', where clathrates are only modeled for the conductive lid.
        self.clathType = None  # Type of model for sI methane clathrates in outer ice shell. Options are 'top', 'bottom', and 'whole', and indicate where clathrates are allowed to be and which type of model to use.
        self.TbIII_K = None  # Temperature at bottom of ice III underplate layer in K. Ranges from 248.85 to 256.164 K for transition to ice V and from 251.165 to 256.164 K for melting temp.
        self.TbV_K = None  # Temperature at bottom of ice V underplate layer in K. Ranges from 256.164 to 272.99 K for melting temp, and 218 to 272.99 for transition to ice VI.
        self.J2 = None  # Gravitational coefficient associated with oblateness, in Schmidt normalization
        self.C20 = None  # Negative of J2, only one of them needs to be set.
        self.C22 = None  # Gravitational coefficient associated with elongation, in Schmidt normalization
        self.C21 = None  # Additional gravitational coefficients that are usually set to zero.
        self.S21 = None
        self.S22 = None
        self.zbChangeTol_frac = 0.05  # Fractional change tolerance, which if exceeded, triggers IceConvect to run a second time for better self-consistency


""" Runtime flags """
class DoSubstruct:

    def __init__(self):
        self.Fe_CORE = False  # Whether to model an iron core for this body
        self.CONSTANT_INNER_DENSITY = False  # Whether to use a fixed density in silicates and core instead of using Perple_X EOS for each
        self.CLATHRATE = False  # Whether to model clathrates
        self.NO_H2O = False  # Whether to model waterless worlds (like Io)
        self.BOTTOM_ICEIII = False  # Whether to allow Ice III between ocean and ice I layer, when ocean temp is set very low- default is that this is off, can turn on as an error condition
        self.BOTTOM_ICEV = False  # Same as above but also including ice V. Takes precedence (forces both ice III and V to be present).
        self.NO_ICE_CONVECTION = False  # Whether to suppress convection in ice layers
        self.EQUIL_Q = True  # Whether to set heat flux from interior to be consistent with heat released through convective profile
        self.POROUS_ICE = False  # Whether to model porosity in ice
        self.POROUS_ROCK = False  # Whether to model porosity in silicates
        self.P_EFFECTIVE = True  # Whether to use effective pressure, modeled as lithostatic less hydrostatic pressure, to determine pore closure behavior (see Vitovtova et al., 2014)
        self.IONOS_ONLY = False  # Whether to ignore conducting layers within the body and model magnetic induction happening only in the ionosphere
        self.TAUP_SEISMIC = False  # Whether to make TauP model files and some basic plots using obspy.taup
        self.FIXED_POROSITY = False  # Whether to force tidal heating to vary instead of porosity to find a matching MoI for bodies with no iron core
        self.PORE_EOS_DIFFERENT = False  # Whether a salinity and/or composition has been set for pores that differs from the ocean


""" Layer step settings """
class StepsSubstruct:

    def __init__(self):
        self.nIceI = None  # Fixed number of steps in outermost ice shell
        self.nClath = None  # Fixed number of steps in clathrates
        self.nHydroMax = None  # Derived working length of hydrosphere layers, gets truncated after layer calcs
        self.nOceanMax = None  # Derived working length of ocean layers, also truncated after layer calcs
        self.nHydro = None  # Derived final number of steps in hydrosphere
        self.nTotal = None  # Total number of layers in profile
        self.nIbottom = None  # Derived number of clathrate + ice I layers
        self.nIIIbottom = 0  # Derived number of clathrate + ice I + ice III layers
        self.nSurfIce = 0  # Derived number of outer ice layers (above ocean) -- sum of nIceI, nClath, nIceIIILitho, nIceVLitho
        self.nRefRho = 30  # Number of values for plotting reference density curves (sets resolution)
        self.iSilStart = None  # Hydrosphere index at which to start silicate size search
        self.nSilMax = None  # Fixed max number of steps in silicate layers
        self.nSil = None  # Derived final number of steps in silicate layers
        self.nCore = None  # Fixed number of steps in core layers, if present
        self.nIceIIILitho = 100  # Fixed number of layers to use for ice III when either BOTTOM_ICEIII or BOTTOM_ICEV is True.
        self.nIceVLitho = 100  # Fixed number of layers to use for ice V when BOTTOM_ICEV is True.
        self.nPsHP = 150  # Number of interpolation steps to use for getting HP ice EOS (pressures)
        self.nTsHP = 100  # Number of interpolation steps to use for getting HP ice EOS (temperatures)
        self.nPoros = 10  # Number of steps in porosity to use in geometric series between phiMin and phiMax for porous rock when no core is present
        self.iCond = []  # Logical array to select indices corresponding to surface ice I
        self.iConv = []  # As above, for convecting ice I (and lower TBL)
        self.iCondIII = []  # As above, for conducting ice III
        self.iConvIII = []  # As above, for convecting ice III
        self.iCondV = []  # As above, for conducting ice V
        self.iConvV = []  # As above, for convecting ice V


""" Hydrosphere assumptions """
class OceanSubstruct:

    def __init__(self):
        self.comp = None  # Type of dominant dissolved salt in ocean. Options: 'Seawater', 'MgSO4', 'PureH2O', 'NH3', 'NaCl', 'none'
        self.wOcean_ppt = None  # (Absolute) salinity: Mass concentration of above composition in parts per thousand (ppt)
        self.sigmaMean_Sm = np.nan  # Mean conductivity across all ocean layers (linear average, ignoring spherical geometry effects)
        self.sigmaTop_Sm = np.nan  # Conductivity of shallowest ocean layer
        self.deltaP = None  # Increment of pressure between each layer in lower hydrosphere/ocean (sets profile resolution)
        self.deltaT = None  # Step size in K for temperature values used in generating ocean EOS functions. If set, overrides calculations that otherwise use the specified precision in Tb_K to determine this.
        self.Vtot_m3 = None  # Total volume of all ocean layers
        self.rhoMean_kgm3 = None  # Mean density for ocean layers
        self.Tmean_K = None  # Mean temperature of ocean layers based on total thermal energy
        self.rhoCondMean_kgm3 = {phase: np.nan for phase in ['Ih', 'II', 'III', 'V', 'VI', 'Clath']}  # Mean density for conducting ice layers
        self.rhoConvMean_kgm3 = {phase: np.nan for phase in ['Ih', 'II', 'III', 'V', 'VI', 'Clath']}  # Mean density for convecting ice layers
        self.sigmaCondMean_Sm = {phase: np.nan for phase in ['Ih', 'II', 'III', 'V', 'VI', 'Clath']}  # Mean conductivity for conducting ice layers
        self.sigmaConvMean_Sm = {phase: np.nan for phase in ['Ih', 'II', 'III', 'V', 'VI', 'Clath']}  # Mean conductivity for convecting ice layers
        self.sigmaCondMean_Sm = {phase: np.nan for phase in ['Ih', 'II', 'III', 'V', 'VI', 'Clath']}  # Mean conductivity for conducting ice layers
        self.sigmaConvMean_Sm = {phase: np.nan for phase in ['Ih', 'II', 'III', 'V', 'VI', 'Clath']}  # Mean conductivity for convecting ice layers
        self.GScondMean_GPa = {phase: np.nan for phase in ['Ih', 'II', 'III', 'V', 'VI', 'Clath']}  # Mean shear modulus for conducting ice layers
        self.GSconvMean_GPa = {phase: np.nan for phase in ['Ih', 'II', 'III', 'V', 'VI', 'Clath']}  # Mean shear modulus for convecting ice layers
        self.rhoMeanVwet_kgm3 = np.nan  # Mean density for in-ocean ice V layers
        self.rhoMeanVI_kgm3 = np.nan  # Mean density for in-ocean ice VI layers
        self.sigmaMeanVwet_Sm = np.nan  # Mean electrical conductivity for in-ocean ice V layers
        self.sigmaMeanVI_Sm = np.nan  # Mean electrical conductivity for in-ocean ice VI layers
        self.GSmeanVwet_GPa = np.nan  # Mean shear modulus for in-ocean ice V layers
        self.GSmeanVI_GPa = np.nan  # Mean shear modulus for in-ocean ice VI layers
        self.TfreezeOffset_K = 0.01  # Offset from the freezing temperature to avoid overshooting in HP ices
        self.koThermI_WmK = 2.21  # Thermal conductivity of ice I at melting temp. Default is from Eq. 6.4 of Melinder (2007), ISBN: 978-91-7178-707-1
        self.dkdTI_WmK2 = -0.012  # Temperature derivative of ice I relative to the melting temp. Default is from Melinder (2007).
        self.sigmaIce_Sm = {'Ih':1e-8, 'II':1e-8, 'III':1e-8, 'V':1e-8, 'VI':1e-8, 'Clath':5e-5}  # Assumed conductivity of solid ice phases (see Constants.sigmaClath_Sm below)
        self.THydroMax_K = 320  # Assumed maximum ocean temperature for generating ocean EOS functions. For large bodies like Ganymede, Callisto, and Titan, larger values are required.
        self.PHydroMax_MPa = 200  # Guessed maximum pressure of the hydrosphere in MPa. Must be greater than the actual pressure, but ideally not by much. Sets initial length of hydrosphere arrays, which get truncated after layer calculations are finished.
        self.MgSO4elecType = 'Vance2018'  # Type of electrical conductivity model to use for MgSO4. Options: 'Vance2018', 'Pan2020'
        self.MgSO4scalingType = 'Vance2018'  # Type of scaling to apply to Larionov and Kryukov model. Options: 'Vance2018', 'LK1984'
        self.MgSO4rhoType = 'Millero'  # Type of water density model to use in Larionov and Kryukov model. Options: 'Millero', 'SeaFreeze'
        self.phaseType = 'lookup'  # Type of phase calculation to use for MgSO4 and pure water. Currently, "lookup" runs a fast lookup table like the Perplex EOS functions, and anything else forces a (slow) individual calc for each P, T point.
        self.QfromMantle_W = None  # Heat flow from mantle into hydrosphere (calculated from ice thermal profile and applied to mantle)
        self.EOS = None  # Equation of state data to use for ocean layers
        self.meltEOS = None  # EOS just for finding ice I/liquid transition pressure Pb
        self.surfIceEOS = {}  # Equation of state data to use for surface ice layers
        self.iceEOS = {}  # Equation of state data to use for ice layers within the ocean
        """ Porosity parameters """
        self.phiMax_frac = {'Ih':0.2, 'II':0.2, 'III':0.2, 'V':0.2, 'VI':0.2, 'Clath':0.2}  # Porosity (void fraction) of ices, extrapolated down to 0 pressure. All arbitrary at 0.2, as placeholders.
        self.Pclosure_MPa = {'Ih':20, 'II':200, 'III':200, 'V':200, 'VI':300, 'Clath':120}  # Pressure threshold in MPa beyond which pores in ice shut completely and porosity drops rapidly to zero, for use in Han et al. (2014) model.
        self.porosType = {'Ih':'Han2014', 'II':'Han2014', 'III':'Han2014', 'V':'Han2014', 'VI':'Han2014', 'Clath':'Han2014'}  # Porosity model to apply for ices. Valid options currently only include 'Han2014', because other models use Earth-crust-based PREM model lookup.
        self.phiMin_frac = 1e-4  # Minimum porosity to model in ice, i.e. below this we just set to zero
        self.alphaPeff = {'Ih':0.95, 'II':0.95, 'III':0.95, 'V':0.95, 'VI':0.95, 'Clath':0.95}  # Scaling factor for Peff in pores. See Sil.alphaPeff
        # See descriptions of these quantities in SilSubstruct -- most will be the same here and there but seismic properties especially may differ.
        self.Jrho = 1
        self.JkTherm = 1
        self.Jalpha = 1
        self.JCp = 1
        self.Jsigma = 1
        self.JKS = 0.35
        self.JGS = 0.35
        self.JVP = 0.75
        self.JVS = 0.85


""" Silicate layers """
class SilSubstruct:

    def __init__(self):
        self.alphaPeff = 0.95  # Scaling factor by which the hydrostatic pressure of pore-filling fluids reduces pore closure pressure. Value taken from the high end of Vitovtova et al. (2014), based on Bernabe (1987). Vitovtova et al. also note that granites vary from 0.56-1.05 and sandstones from 0.35-0.75.
        self.sigmaSil_Sm = 1e-16  # Assumed conductivity of silicate rock
        self.Qrad_Wkg = 0  # Average radiogenic heating rate for silicates in W/kg.
        self.Htidal_Wm3 = 0  # Average tidal heating rate for silicates in W/m^3.
        self.HtidalMin_Wm3 = 1e-12  # Minimum average tidal heating rate in silicates in W/m^3 above zero to start with for finding MoI match with no core.
        self.HtidalMax_Wm3 = 1e-7  # Maximum average tidal heating to stop MoI search
        self.deltaHtidal_logUnits = 1/3  # Step size by which to increment Htidal_Wm3 for finding MoI match with no core.
        self.kTherm_WmK = None  # Constant thermal conductivity to set for a specific body (overrides Constants.kThermSil_WmK)
        """ Porosity parameters """
        self.phiRockMax_frac = None  # Porosity (void fraction) of the rocks in vacuum. This is the expected value for core-less bodies, and porosity is modeled for a range around here to find a matching MoI. For bodies with a core, this is a fixed value for rock porosity at P=0.
        self.phiRangeMult = 8  # Factor by which to divide the distance from user-defined phiRockMax_frac to 0 and 1 to obtain the search range for phi
        self.Pclosure_MPa = 350  # Pressure threshold in MPa beyond which pores in silicates shut completely and porosity drops to zero, for use in Han et al. (2014) model. See Saito et al. (2016) for evidence of values up to ~750 MPa: https://doi.org/10.1016/j.tecto.2016.03.044
        self.porosType = 'Han2014'  # Porosity model to apply for silicates. Options are 'Han2014', 'Vitovtova2014', 'Chen2020'.
        self.poreComp = None  # Solute to use for pore fluids. If None (not set), Ocean.comp is used.
        self.wPore_ppt = None  # Salinity of pore fluids in ppt. If None, Ocean.wOcean_ppt is used.
        self.PHydroMax_MPa = 2250  # Maximum pressure to evaluate for pore fluid EOS in MPa
        self.THydroMax_K = 500  # Maximum temperature to evaluate for pore fluid EOS in K
        self.poreEOS = None  # EOS to use for pore fluids.
        self.poreH2Orho_kgm3 = 1023  # Arbitrary ocean water density used to calculate filled-pore rock densities. This value is that used in Vance et al. (2018).
        self.poreConductPrefac = 3.3  # Prefactor by which to multiply pore fluid conductivity at porosities below some threshold, based on measurements from Wong et al. (1984).
        self.poreConductBelowExp = 2.3  # Exponent to use for below-threshold dependence of conductivity on porosity for pore fluids
        self.poreConductAboveExp = 1.5  # Exponent to use for above-threshold dependence of conductivity on porosity for pore fluids -- 3/2 is that expected for spherical grains of the matrix material; Wong et al. (1984) suggest it fits well to above-threshold measurements, and this is supported by Golden et al. (2007): https://doi.org/10.1029/2007GL030447
        self.poreConductThresh_frac = self.poreConductPrefac**(1/(self.poreConductAboveExp - self.poreConductBelowExp))  # Threshold in porosity at which the below-threshold Archie's law fit from above parameters is equal to above-threshold dependence
        self.phiMin_frac = 1e-4  # Minimum porosity to model in silicates, i.e. below this we just set to zero
        self.sigmaPoreMean_Sm = np.nan  # Mean conductivity of pore fluids across all layers with porosity above poreConductThresh (about 5%)
        self.sigmaPorousLayerMean_Sm = np.nan  # Mean conductivity of matrix + pore fluid combined, for all layers with porosity of poreConductThresh
        # J values for exponent of pore property / overlaying matrix combination as in Yu et al. (2016): http://dx.doi.org/10.1016/j.jrmge.2015.07.004
        # Values X combine as Xtot^J = Xmatrix^J * (1 - phi) + Xpore^J * phi, where phi is the porosity (volume fraction) of pore space filled with the secondary material.
        # Values of J typically range from -1 to 1, where 1 is a direct mean, e.g. Xtot = Xrock*(1-phi) + Xpore*phi and -1 is a geometric mean, e.g. Xtot = Xrock || Xpore = 1/((1-phi)/Xrock + phi/Xpore).
        # J can range outside -1 and 1 for VP and VS, according to Yu et al.
        self.Jrho = 1  # Mass density, directly summative
        self.JkTherm = 1  # Should add like conductors in parallel
        self.Jalpha = 1  # Volume increase with temperature should be summative
        self.JCp = 1  # Heat capacity is summative, but since porosity is volume fraction the inputs need to be scaled by mass density
        self.Jsigma = 1  # Pore fluid and matrix provide parallel electrically conducting pathways as well
        # For the following 4 properties, Yu et al. (2016) have values for rock combinations that vary around these marks, but they vary quite a lot. If the values are known
        # for the specific assumed composition of the silicate layer, those values should be used instead for the particular body being simulated.
        self.JKS = 0.35
        self.JGS = 0.35
        self.JVP = 0.75
        self.JVS = 0.85
        """ Mantle Equation of State (EOS) model """
        self.mantleEOS = None  # Equation of state data to use for silicates
        self.mantleEOSName = None  # Same as above but containing keywords like clathrates in filenames
        self.mantleEOSDry = None  # Name of mantle EOS to use assuming non-hydrated silicates
        self.EOS = None  # Interpolator functions for evaluating Perple_X EOS model
        self.rhoSilWithCore_kgm3 = 3300  # Assumed density of silicates when a core is present in kg/m^3
        # Derived quantities
        self.Rmean_m = None  # Mantle radius for mean compatible moment of inertia (MoI)
        self.Rrange_m = None  # Mantle radius range for compatible MoI
        self.Rtrade_m = None  # Array of mantle radii for compatible MoIs
        self.rhoMean_kgm3 = None  # Mean mantle density determined from MoI calculations
        self.GSmean_GPa = None  # Mean shear modulus in silicate layers
        self.HtidalMean_Wm3 = None  # Mean tidal heating in silicate layers
        self.sigmaMean_Sm = np.nan  # Mean conductivity across all silicate layers, ignoring spherical effects
        self.rhoTrade_kgm3 = None  # Array of mantle densities for compatible MoIs for core vs. mantle tradeoff plot
        self.mFluids = None  # WIP for tracking loss of fluids along the geotherm -- needs a better name.


""" Core layers """
class CoreSubstruct:

    def __init__(self):
        self.rhoFe_kgm3 = 8000  # Assumed density of pure iron in kg/m^3
        self.rhoFeS_kgm3 = 5150  # Assumed density of iron sulfide in kg/m^3
        self.rhoMin_kgm3 = 5150  # Assumed minimum possible density for the core in kg/m^3. Sets maximum core size.
        self.rhoPoFeFCC = None  # Density of pyrrhottite plus face-centered cubic iron
        self.sigmaCore_Sm = 1e6  # Fixed electrical conductivity to apply to core (typically low, to ignore core impacts on induction)
        self.coreEOS = 'Fe075_S025.tab'  # Default core EOS to use
        self.EOS = None  # Interpolator functions for evaluating Perple_X EOS model
        self.kTherm_WmK = None  # Constant thermal conductivity to set for a specific body (overrides Constants.kThermFe_WmK)
        # Derived quantities
        self.rho_kgm3 = None  # Core bulk density consistent with assumed mixing ratios of Fe, FeS, etc.
        self.rhoMean_kgm3 = None  # Core bulk density calculated from final MoI match using EOS properties
        self.rhoMeanFe_kgm3 = np.nan  # Pure iron layer bulk density calculated from final MoI match using EOS properties
        self.rhoMeanFeS_kgm3 = np.nan  # FeS layer bulk density calculated from final MoI match using EOS properties
        self.sigmaMean_Sm = np.nan  # Mean conductivity across all core layers, ignoring spherical effects
        self.GSmean_GPa = np.nan  # Mean shear modulus in iron core layers
        self.GSmeanFe_GPa = np.nan  # Mean shear modulus in pure Fe core layers
        self.GSmeanFeS_GPa = np.nan  # Mean shear modulus in FeS core layers
        self.Rmean_m = None  # Core radius for mean compatible moment of inertia (MoI)
        self.Rrange_m = None  # Core radius range for compatible MoI
        self.Rtrade_m = None  # Array of core radii for compatible MoIs
        self.Rset_m = None  # Value to set the core outer radius to, when we have already found it via e.g. using CONSTANT_INNER_DENSITY = True. Used to recycle SilicateLayers when we don't want to do an MoI search with the EOS functions.
        self.xS_ppt = None  # Mass fraction of sulfur in the core in ppt
        # 2021-12-30: Judging by usage of various different fractional variables in the literature and in
        # the Matlab code, these x variables should be molar fractions (# this species/total # molecules).
        self.xFeSmeteoritic = None  # CM2 mean from Jarosewich 1990
        self.xH2O = None  # Total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model


""" Seismic properties """
class SeismicSubstruct:
    """ Calculations based on Cammarano et al., 2006 (DOI: 10.1029/2006JE002710)
        in a parameterization for the shear anelasticity quality factor QS (Eqs. 4-6).
        g : Homologous temperature scaling constant- dimensionless constant described as a temperature scaling relative to the melting temperature
        B : Shear anelasticity/quality factor normalization constant, dimensionless - helps quantify the effects of anelastic attenuation on the seismic wavelet caused by fluid movement and grain boundary friction
        gamma : Exponent on the seismic frequency omega, dimensionless - helps describe the frequency dependence of attenuation
        """

    def __init__(self):
        self.lowQDiv = None  # Factor by which to divide the seismic attenuation Q to test out a low-Q value, dimensionless
        # Attenuation Parameters Based on those Described in Cammarano et al. 2006
        # Placeholder for clathrates, matching ice I
        self.BClath = 0.56
        self.gammaClath = 0.2
        self.gClath = 22.0
        # Ice I
        self.BIceI = 0.56
        self.gammaIceI = 0.2
        self.gIceI = 22.0
        # Ice II
        self.BIceII = 0.56
        self.gammaIceII = 0.2
        self.gIceII = 30.0
        # Ice III
        self.BIceIII = 0.56
        self.gammaIceIII = 0.2
        self.gIceIII = 25.0
        # Ice V
        self.BIceV = 0.56
        self.gammaIceV = 0.2
        self.gIceV = 27.0
        # Ice VI
        self.BIceVI = 0.56
        self.gammaIceVI = 0.2
        self.gIceVI = 28.0
        # Silcate mantle
        self.BSil = 0.56
        self.gammaSil = 0.2
        self.gSil = 30.0
        # Core
        self.QScore = None  # For assigning a QS value to the core, in lieu of a calculation.
        # Derived quantities
        self.VP_kms = None  # Longitudinal (p-wave) sound velocity for each layer in km/s
        self.VS_kms = None  # Shear (s-wave) sound velocity for each layer in km/s
        self.QS = None  # Anelastic shear quality factor Q_S of each layer, divided by omega^gamma to remove frequency dependence. Essentially the ratio of total seismic energy to that lost per cycle, see Stevenson (1983).
        self.KS_GPa = None  # Bulk modulus of each layer in GPa
        self.GS_GPa = None  # Shear modulus of each layer in GPa


""" Magnetic induction """
class MagneticSubstruct:

    def __init__(self):
        # Input settings
        self.SCera = None  # Spacecraft era to use for excitation moments. Read from Be1xyz file name.
        self.extModel = None  # External field model to use for excitation moments. Read from Be1xyz file name.
        self.inductMethod = 'Srivastava1966'  # Type of magnetic induction model to use. Options are "Srivastava1966" or "layer" for layer method and "Eckhardt1963" or "numeric" for numeric method (symmetric only).
        self.Texc_hr = None  # Periods in hr of peaks in magnetic excitation spectrum
        self.omegaExc_radps = None  # Angular frequency of peaks in magnetic excitation spectrum in rad/s
        self.ionosBounds_m = None  # Upper altitude cutoff for ionosphere layers in m. Omit the surface (don't include 0 in the list).
        self.sigmaIonosPedersen_Sm = 1e-4  # Pedersen conductivity for ionospheric layers in S/m. Length must match ionosBounds_m. The default value (set here) is set uniform when ionosBounds_m has size 1, and set uniform between entries 1 and 2 when it has size 2 (with zero conductivity between).
        self.rSigChange_m = None  # Radii of outer boundary of each conducting layer in m (i.e., radii where sigma changes)
        self.sigmaLayers_Sm = None  # Reduced set of conductivity values compatible with rSigChange_m that will work in induction calculations (i.e. all non-zero)
        self.nBds = None  # Number of radial boundaries between conductors specified
        self.nExc = None  # Number of excitation frequencies modeled
        self.Benm_nT = None  # Excitation moments (amplitude and phase of magnetic oscillations) applied to the body in nT
        self.nprmMax = 1  # Maximum n' to use in excitation moments (1 for uniform).
        self.pMax = None  # Maximum p to use in asymmetric shape (0 for spherically symmetric; 2 for up to and incuding gravity coefficients).
        self.asymShape_m = None  # Asymmetric shape to use in induction calculations. Only used when inductType = "Srivastava1966".
        self.gravShape_m = None  # Asymmetric shape in p = 2 coefficients to use in induction calculations. Only used when inductType = "Srivastava1966".
        self.Xid = None  # Mixing coefficient table used in asymmetric induction calculations
        # Output calculations
        self.Aen = None  # Complex response amplitude of magnetic excitation for dipole moment for each excitation frequency (unitless)
        self.Amp = None  # Amplitude (modulus) of magnetic excitation for spherically symmetric approximation for each excitation frequency (unitless)
        self.phase = None  # Phase delay of magnetic excitation for spherically symmetric approximation for each excitation frequency in degrees
        self.Binm_nT = None  # Induced magnetic moments relative to the body surface in nT
        self.nLin = None  # Linear list of n values for output Binm with shape matching mLin, so that n,m pairs corresponding to each BinmLin are read as (n[j], m[j]) -> BinmLin[i, j]
        self.mLin = None  # m values corresponding to each n in nLin
        self.BinmLin_nT = None  # Linear form of Binm_nT, with shape (nExc, (nPrmMax+pMax+1)**2 - 1), such that BinmLin[i, j] = Binm[i, int(m[j]<0), n[j], m[j]]


""" Main body profile info--settings and variables """
class PlanetStruct:

    # Require a body name as an argument for initialization; define instance attributes
    def __init__(self, name):
        self.name = name
        if self.name[:4] == 'Test':
            self.bodyname = 'Test'
        else:
            self.bodyname = self.name

        self.Bulk = BulkSubstruct()
        self.Do = DoSubstruct()
        self.Steps = StepsSubstruct()
        self.Ocean = OceanSubstruct()
        self.Sil = SilSubstruct()
        self.Core = CoreSubstruct()
        self.Seismic = SeismicSubstruct()
        self.Magnetic = MagneticSubstruct()

        self.saveLabel = None # Label for savefile
        # Settings for GetPfreeze start, stop, and step size.
        # Shrink closer to expected melting pressure to improve run times.
        self.PfreezeLower_MPa = 0.05  # Lower boundary for GetPfreeze to search for ice Ih phase transition
        self.PfreezeUpper_MPa = 230  # Upper boundary for GetPfreeze to search for ice Ih phase transition
        self.PfreezeRes_MPa = 0.05  # Step size in pressure for GetPfreeze to use in searching for phase transition

        """ Derived quantities (assigned during PlanetProfile runs) """
        # Layer arrays
        self.phase = None  # Phase of the layer input as an integer: ocean=0, ice I through VI are 1 through 6, clathrate=Constants.phaseClath, silicates=Constants.phaseSil, iron=Constants.phaseFe.
        self.r_m = None  # Distance from center of body to the outer bound of current layer in m
        self.z_m = None  # Distance from surface of body to the outer bound of current layer in m
        self.T_K = None  # Temperature of each layer in K
        self.P_MPa = None  # Pressure at top of each layer in MPa
        self.rho_kgm3 = None  # Overall mass density of each layer in kg/m^3
        self.g_ms2 = None  # Gravitational acceleration at top of each layer, m/s^2
        self.Cp_JkgK = None  # Heat capacity at constant pressure for each layer's material in J/kg/K
        self.alpha_pK = None  # Thermal expansivity of layer material in K^-1
        self.phi_frac = None  # Porosity of each layer's material as a fraction of void/solid
        self.sigma_Sm = None  # Electrical conductivity (sigma) in S/m of each conducting layer
        self.MLayer_kg = None  # Mass of each layer in kg
        self.VLayer_m3 = None  # Volume of each layer in m^3
        self.kTherm_WmK = None  # Thermal conductivity of each layer in W/(m K)
        self.Htidal_Wm3 = None  # Tidal heating rate of each layer in W/m^3
        self.Ppore_MPa = None  # Pressure of fluids assumed to occupy full pore space
        self.rhoMatrix_kgm3 = None  # Mass density of matrix material (rock or ice)
        self.rhoPore_kgm3 = None  # Mass density of pore material (typically ocean water)
        # Individual calculated quantities
        self.Mtot_kg = None  # Total calculated mass selected from MoI matching
        self.zb_km = None  # Thickness of outer ice shell/depth of ice-ocean interface in km
        self.qSurf_Wm2 = None  # Heat flux at the surface derived from layer profile and input parameters
        self.qCon_Wm2 = None  # Heat flux of conducting or convecting ice layers at the bottom of the region used in calculating the Rayleigh number
        self.zClath_m = None  # Thickness of clathrate layer in surface ice in m
        self.D_km = None  # Thickness of ocean layer in km
        self.Pb_MPa = None  # Pressure at ice-ocean interface in MPa
        self.PbI_MPa = None  # Pressure at bottom of surface ice I/clathrate layer in MPa
        self.PbIII_MPa = None  # Pressure at ice III/ice V transition in MPa, only used when BOTTOM_ICEIII or BOTTOM_ICEV is True
        self.PbV_MPa = None  # Pressure at bottom of ice V layer in MPa, only used when BOTTOM_ICEV is True
        self.PbClathMax_MPa = None  # Max pressure at the bottom of a clathrate lid. Actual PbClath may be reduced in convection calculations
        self.TclathTrans_K = None  # Temperature at the transition from clathrates to ice I
        self.CMR2mean = None  # Mean value of axial moment of inertia that is consistent with profile core/mantle trades
        self.CMR2less = None  # Neighboring value to CMR2mean in MoI matching, just below it
        self.CMR2more = None  # Neighboring value above CMR2mean
        self.Tconv_K = None  # Temperature of "well-mixed" convecting region in ice I layer in K
        self.TconvIII_K = None  # Same as above but for ice III underplate layers.
        self.TconvV_K = None  # Same as above but for ice V underplate layers.
        self.etaConv_Pas = None  # Viscosity of ice I at Tconv_K
        self.etaConvIII_Pas = None  # Same as above but for ice III underplate layers.
        self.etaConvV_Pas = None  # Same as above but for ice V underplate layers.
        self.eLid_m = None  # Thickness of conducting stagnant lid layer in m.
        self.eLidIII_m = None  # Same as above but for ice III underplate layers.
        self.eLidV_m = None  # Same as above but for ice V underplate layers.
        self.Dconv_m = None  # Thickness of convecting layer in m.
        self.DconvIII_m = None  # Same as above but for ice III underplate layers.
        self.DconvV_m = None  # Same as above but for ice V underplate layers.
        self.deltaTBL_m = None  # Thickness of lower thermal boundary layer in m when Htidal = 0 in the ice
        self.deltaTBLIII_m = None  # Same as above but for ice III underplate layers.
        self.deltaTBLV_m = None  # Same as above but for ice V underplate layers.
        self.RaConvect = None  # Rayleigh number of putative convective layer within the ice I layers. If this number is below Planet.RaCrit, convection does not occur.
        self.RaConvectIII = None  # Same as above but for ice III underplate layers.
        self.RaConvectV = None  # Same as above but for ice V underplate layers.
        self.RaCrit = None  # Critical Rayleigh number that determines whether or not we model convection.
        self.RaCritIII = None  # Same as above but for ice III underplate layers.
        self.RaCritV = None  # Same as above but for ice V underplate layers.
        self.MH2O_kg = None  # Total mass of water molecules contained in ice, liquid, and pore spaces
        self.Mrock_kg = None  # Total mass contained in silicate rock (just the matrix, when layers are porous)
        self.Mcore_kg = None  # Total mass contained in iron core material
        self.Mice_kg = None  # Total mass contained in all ice phases, including gas trapped in clathrates
        self.Msalt_kg = None  # Total mass contained in solutes
        self.MoceanSalt_kg = None  # Total mass contained in ocean solute in ocean layers
        self.MporeSalt_kg = None  # Total mass contained in ocean solute in pore spaces
        self.Mocean_kg = None  # Total mass contained in ocean fluids, including H2O and salts, excluding pores
        self.Mfluid_kg = None  # Sum of the masses in ocean and pore spaces
        self.MporeFluid_kg = None  # Total mass contained in pore fluids, including H2O and salts
        self.Mclath_kg = None  # Total mass of clathrate layers
        self.MclathGas_kg = None  # Total mass of non-water molecules trapped in clathrates
        self.index = None  # Numeric indicator to aid in progress info in multi-model runs
        # Layer thicknesses for table printout
        self.zIceI_m = np.nan
        self.zClath_km = np.nan  # Note this one breaks with the pattern because zClath_m is already in use.
        self.zIceIII_m = np.nan
        self.zIceVund_m = np.nan
        self.zIceV_m = np.nan
        self.zIceVI_m = np.nan
        self.dzIceI_km = np.nan
        self.dzClath_km = np.nan
        self.dzIceIII_km = np.nan
        self.dzIceVund_km = np.nan
        self.dzIceV_km = np.nan
        self.dzIceVI_km = np.nan
        self.dzIceVandVI_km = np.nan
        self.dzSilPorous_km = np.nan
        self.dzFeS_km = np.nan


""" Params substructs """
# Construct filenames for data, saving/reloading
class DataFilesSubstruct:
    def __init__(self, datPath, saveBase, comp, inductBase=None):
        if inductBase is None:
            inductBase = saveBase
        self.path = datPath
        self.inductPath = os.path.join(self.path, 'inductionData')
        if not self.path == '' and not os.path.isdir(self.path):
            os.makedirs(self.path)
        if not self.path == '' and not os.path.isdir(self.inductPath):
            os.makedirs(self.inductPath)

        self.fName = os.path.join(self.path, saveBase)
        self.saveFile = self.fName + '.txt'
        self.mantCoreFile = self.fName + '_mantleCore.txt'
        self.permFile = self.fName + '_mantlePerm.txt'
        self.fNameInduct = os.path.join(self.inductPath, saveBase)
        self.inductLayersFile = self.fNameInduct + '_inductLayers.txt'
        self.inducedMomentsFile = self.fNameInduct + '_inducedMoments.txt'
        self.fNameInductOgram = os.path.join(self.path, 'inductionData', inductBase)
        self.inductOgramFile = self.fNameInductOgram + f'{comp}_inductOgram.mat'
        self.inductOgramSigmaFile = self.fNameInductOgram + '_sigma_inductOgram.mat'


# Construct filenames for figures etc.
class FigureFilesSubstruct:
    def __init__(self, figPath, figBase, xtn, comp=None, inductBase=None):
        if inductBase is None:
            self.inductBase = figBase
        else:
            self.inductBase = inductBase
        if comp is None:
            self.comp = ''
        else:
            self.comp = comp
        self.path = figPath
        self.inductPath = os.path.join(self.path, 'induction')
        if not self.path == '' and not os.path.isdir(self.path):
            os.makedirs(self.path)
        if not self.path == '' and not os.path.isdir(self.inductPath):
            os.makedirs(self.inductPath)
        self.fName = os.path.join(self.path, figBase)
        self.fNameInductOgram = os.path.join(self.inductPath, self.inductBase + self.comp)

        # Figure filename strings
        vsP = 'Porosity_vs_P'
        vsR = 'Porosity_vs_R'
        vperm = 'Permeability'
        vgsks = 'Gs_Ks'
        vseis = 'Seismic'
        vhydro = 'Hydrosphere'
        vgrav = 'Gravity'
        vmant = 'MantleDens'
        vcore = 'CoreMantTrade'
        vpvt4 = 'PTx4'
        vpvt6 = 'PTx6'
        vwedg = 'Wedge'
        induct = 'inductOgram'
        sigma = 'inductOgramSigma'
        # Construct Figure Filenames
        self.vwedg = self.fName + vwedg + xtn
        self.vsP = self.fName + vsP + xtn
        self.vsR= self.fName + vsR + xtn
        self.vperm = self.fName + vperm + xtn
        self.vgsks = self.fName + vgsks + xtn
        self.vseis = self.fName + vseis + xtn
        self.vhydro = self.fName + vhydro + xtn
        self.vgrav = self.fName + vgrav + xtn
        self.vmant = self.fName + vmant + xtn
        self.vcore = self.fName + vcore + xtn
        self.vpvt4 = self.fName + vpvt4 + xtn
        self.vpvt6 = self.fName + vpvt6 + xtn
        self.phaseSpace =            f'{self.fNameInductOgram}_{induct}_phaseSpace{xtn}'
        self.phaseSpaceCombo =       f'{os.path.join(self.inductPath, self.inductBase)}Compare_{induct}_phaseSpace{xtn}'
        self.induct =        {zType: f'{self.fNameInductOgram}_{induct}_{zType}{xtn}' for zType in ['Amp', 'Bx', 'By', 'Bz', 'Bcomps']}
        self.inductCompare = {zType: f'{self.fNameInductOgram}Compare_{zType}{xtn}' for zType in ['Amp', 'Bx', 'By', 'Bz', 'Bcomps']}
        self.sigma =         {zType: f'{self.fNameInductOgram}_{sigma}_{zType}{xtn}' for zType in ['Amp', 'Bx', 'By', 'Bz', 'Bcomps']}
        self.sigmaOnly =     {zType: f'{self.fNameInductOgram}_{sigma}Only_{zType}{xtn}' for zType in ['Amp', 'Bx', 'By', 'Bz', 'Bcomps']}


""" General parameter options """
class ParamsStruct:
    # Do not set any values below. All values are assigned in PlanetProfile.GetConfig.
    def __init__(self):
        self.DataFiles = DataFilesSubstruct('', '', '')
        self.FigureFiles = FigureFilesSubstruct('', '', '')
        self.Sig = None  # General induction settings
        self.Induct = None  # Induction calculation settings
        self.MagSpectrum = None  # Excitation spectrum settings
        self.cLevels = None  # Contour level specifications
        self.cFmt = None  # Format of contour labels
        self.compareDir = 'Comparison'
        
        
""" Inductogram settings """
class InductOgramParamsStruct:
    # Do not set any values below (except V2021 values). All other values are assigned in PlanetProfile.GetConfig.
    def __init__(self, inductOtype, cLevels, dftC, cFmt):
        self.bodyname = None
        self.inductOtype = inductOtype
        self.cLevels = cLevels
        self.dftC = dftC
        self.cFmt = cFmt
        self.colorType = 'zb'  # What parameter to use for color of points in phase space plots. Options are "Tmean", "zb".
        self.SPECIFIC_CLEVELS = False  # Whether to use the specific cLevels listed below or default numbers
        self.excSelectionCalc = {'synodic': True, 'orbital': True, 'true anomaly': True,  'synodic harmonic': True}  # Which magnetic excitations to include in calculations
        self.excSelectionPlot = {'synodic': True, 'orbital': True, 'true anomaly': False, 'synodic harmonic': True}  # Which magnetic excitations to include in plotting
        self.nwPts = None  # Resolution for salinity values in ocean salinity vs. other plots
        self.wMin = None
        self.wMax = None
        self.nTbPts = None  # Resolution for Tb values in ocean salinity/Tb plots
        self.TbMin = None
        self.TbMax = None
        self.nphiPts = None  # Resolution for phiRockMax values in ocean salinity/phiMax plots
        self.phiMin = None
        self.phiMax = None
        self.nrhoPts = None  # Resolution for silicate density values in ocean salinity/rho plots
        self.rhoMin = None
        self.rhoMax = None
        self.nSigmaPts = None  # Resolution for conductivity values in ocean conductivity/thickness plots
        self.sigmaMin = None
        self.sigmaMax = None
        self.nDpts = None  # Resolution for ocean thickness as for conductivity
        self.Dmin = None
        self.Dmax = None
        self.zbFixed_km = None
        self.EckhardtSolveMethod = 'RK45'  # Numerical solution method for scipy.integrate.solve_ivp. See https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
        self.rMinODE = 1e3  # Minimum radius to use for numerical solution. Cannot be zero because of singularity at the origin.
        self.oceanInterpMethod = 'linear'  # Interpolation method for determining ocean conductivities when REDUCED_INDUCT is True.
        self.nIntL = 5  # Number of ocean layers to use when REDUCED_INDUCT = 1

        # Plot settings to match inductograms from Vance et al. (2021): https://doi.org/10.1029/2020JE006418
        self.V2021_D_km = [91, 117, 96, 124,
                           91, 117, 91, 119]
        self.V2021_sigma_Sm = [0.4132, 0.4533, 3.3661, 3.7646,
                               0.3651, 0.3855, 2.8862, 3.0760]
        self.V2021_faceColors = [None, None, 'b', 'm',
                                 None, None, 'c', '#b000ff']
        self.V2021_edgeColors = ['b', 'm', 'k', 'k',
                                 'c', '#b000ff', 'k', 'k']
        self.V2021_symbols = ['^', 'v', '^', 'v',
                              '^', 'v', '^', 'v']

    def GetClevels(self, zName, Tname):
        if self.bodyname is None:
            bodyname = 'Europa'
        else:
            bodyname = self.bodyname
        if self.SPECIFIC_CLEVELS:
            theClevels = self.cLevels[bodyname][Tname][zName]
        else:
            theClevels = self.dftC
        return theClevels

    def GetCfmt(self, zName, Tname):
        if self.bodyname is None:
            bodyname = 'Europa'
        else:
            bodyname = self.bodyname
        return self.cFmt[bodyname][Tname][zName]


""" General induction calculation settings """
class ConductLayerParamsStruct:
    # Do not set any values below. All values are assigned in PlanetProfile.GetConfig.
    def __init__(self):
        self.REDUCED_INDUCT = True  # Whether to limit number of ocean layers for faster computation of layered induction
        self.INCLUDE_ASYM = False  # Whether to include asymmetry in the induction conductivity profile based on J2 and C22 values
        self.CONCENTRIC_ASYM = False  # Whether to map a single asymmetric shape to all layers, concentrically, scaling by their radii.
        self.ALLOW_LOW_PMAX = False  # Whether to allow Magnetic.pMax to be set to an integer less than 2.
        self.asymFstring = 'Shape_4piNormDepth'


""" Excitation spectrum settings (WIP) """
class ExcitationSpectrumParamsStruct:
    # Do not set any values below. All values are assigned in PlanetProfile.GetConfig.
    def __init__(self):
        self.nOmegaPts = 100  # Resolution in log frequency space for magnetic excitation spectra
        self.nOmegaFine = 1000  # Fine-spacing resolution for log frequency spectrum
        

""" Figure color options """
class ColorStruct:
    # Do not set any values below. All values are assigned in PlanetProfile.GetConfig.
    def __init__(self):
        self.Induction = {'synodic': None, 'orbital': None, 'true anomaly': None, 'synodic harmonic': None}  # Colors for inductOgram plots
        self.ref = None
        self.PALE_SILICATES = False  # Whether to use a lighter color scheme for silicate layers, or a more "orangey" saturated one

        # Wedge diagram color options
        self.none = '#FFFFFF00'
        self.wedgeBd = 'black'
        self.ionoCmapName = None
        self.ionoTop = None
        self.ionoBot = None
        self.ionoN = None
        self.iceIcond = None
        self.iceIconv = None
        self.iceII = None
        self.iceIII = None
        self.iceV = None
        self.iceVI = None
        self.clathCond = None
        self.clathConv = None
        self.oceanCmapName = None
        self.oceanTop = None  # Fraction of ocean colormap to start at (from 0 to 1)
        self.oceanBot = None  # Fraction of ocean colormap to end at (from 0 to 1)
        self.oceanN = None
        self.silPorousCmapName = None
        self.silPorousTop = None
        self.silPorousBot = None
        self.paleSilPorousCmapName = None
        self.paleSilPorousTop = None
        self.paleSilPorousBot = None
        self.silPorousN = None
        self.silCondCmapName = None
        self.silCondTop = None
        self.silCondBot = None
        self.paleSilCondCmapName = None
        self.paleSilCondTop = None
        self.paleSilCondBot = None
        self.silCondN = None
        self.silConvCmapName = None
        self.silConvTop = None
        self.silConvBot = None
        self.silConvN = None
        self.FeS = None
        self.Fe = None
        self.innerCmapName = None

        self.cmapName = {}  # Colormaps for inductogram phase space plots, hydrosphere plots, etc
        self.cmapBounds = {}  # Select only a subset of the available colormap, if we choose to
        self.Tbounds_K = [245.0, 300.0]  # Set temperature bounds to use for colormap normalization
        self.saturation = {}  # Set upper bounds for max concentrations
        # Saturation & color brightness ("value" in HSV) values for salinity/conductivity axis bounds
        self.fresh = [0.5, 1.0]
        self.salty = [1.0, 0.5]


    def SetCmaps(self):
        """ Assign colormaps to make use of the above parameters
        """
        self.ionoCmap = cmasher.get_sub_cmap(self.ionoCmapName, self.ionoTop, self.ionoBot)
        self.oceanCmap = cmasher.get_sub_cmap(self.oceanCmapName, self.oceanTop, self.oceanBot)
        if self.PALE_SILICATES:
            self.silPorousCmap = cmasher.get_sub_cmap(self.paleSilPorousCmapName, 
                                                      self.paleSilPorousTop, self.paleSilPorousBot)
            self.silCondCmap = cmasher.get_sub_cmap(self.paleSilCondCmapName, self.paleSilCondTop, self.paleSilCondBot)
        else:
            self.silPorousCmap = cmasher.get_sub_cmap(self.silPorousCmapName, 
                                                      self.silPorousTop, self.silPorousBot)
            self.silCondCmap = cmasher.get_sub_cmap(self.silCondCmapName, self.silCondTop, self.silCondBot)
        self.silConvCmap = cmasher.get_sub_cmap(self.silConvCmapName, self.silConvTop, self.silConvBot)
        self.innerCmap = get_cmap(self.innerCmapName)
        # Use cmasher to return colormap objects that do the down-select for us
        self.cmap = {comp: cmasher.get_sub_cmap(cmap, self.cmapBounds[comp][0], self.cmapBounds[comp][1])
                     for comp, cmap in self.cmapName.items()}


    def GetNormT(self, T_K):
        """ Calculate normalized temperature to use with colormaps
        """
        return interp1d([self.Tbounds_K[0], self.Tbounds_K[1]], [0.0, 1.0],
                 bounds_error=False, fill_value='extrapolate')(T_K)
    
    def GetSat(self, w_ppt):
        """ Calculate color saturation value based on salinity relative to
            saturation concentration
        """
        return interp1d([0.0, 1.0], [self.fresh[0], self.salty[0]], 
                        bounds_error=False, fill_value=self.salty[0])(w_ppt)

    def GetVal(self, w_ppt):
        """ Calculate color value (light/dark) based on salinity relative to
            saturation concentration
        """
        return interp1d([0.0, 1.0], [self.fresh[1], self.salty[1]], 
                        bounds_error=False, fill_value=self.salty[1])(w_ppt)

    def OceanCmap(self, comps, w_normFrac, Tmean_normFrac, DARKEN_SALINITIES=True):
        """ Get colormap RGBA vectors for each salinity/T combination.

            Args:
                comps (string, shape N): Ocean composition string for each point.
                w_normFrac (float, shape N): Normalized ocean salinities as a fraction
                    of the saturation/maximum concentration used for the colormap.
                Tmean_normFrac (float, shape N): Normalized mean ocean temperatures
                    as a fraction of the range of values to use for the colormap, where
                    0 is at the bottom of the range and 1 is at the top.

            Returns:
                cList (float, shape N x 3): RGB vectors as rows for each combination.
        """

        if DARKEN_SALINITIES:
            # Get the hue for each point by getting the colormap entry for the normalized
            # temperature, then stripping off the alpha channel, then converting to HSV,
            # then taking only the first value in the [H,S,V] output.
            hueMap = np.array([rgb_to_hsv(self.cmap[comp](T)[:3])[0] for comp, T in zip(comps, Tmean_normFrac)])
            # Get the saturation and value lists from the min/max bounds for the
            # normalized salinity
            satMap = self.GetSat(w_normFrac)
            valMap = self.GetVal(w_normFrac)

            cList = hsv_to_rgb(np.column_stack((hueMap, satMap, valMap)))
        else:
            # Just use standard evaluation of the colorbar
            cList = np.row_stack([self.cmap[comp](T) for comp, T in zip(comps, Tmean_normFrac)])

        return cList


""" Figure (line)style settings """
class StyleStruct:
    # Do not set any values below. All values are assigned in PlanetProfile.GetConfig.
    def __init__(self):
        self.LS_dft = None  # Default line style to use on plots
        self.LS_Sw = None  # Linestyle for Seawater
        self.LS_Mg = None  # Linestyle for MgSO4
        self.LS_sp = None  # Linestyle for special consideration models
        self.LW_sal = None  # Linewidth for higher salinity
        self.LW_dil = None  # Linewidth for dilute salinity
        self.LW_std = None  # Linewidth for standard salinity
        self.LW_sound = None  # LineWidth for sound speed plots
        self.LW_seism = None  # LineWidth for seismic plots (Attenuation)
        self.LS_ref = {}  # Style for reference profiles
        self.LW_ref = None  # Linewidth for reference profiles
        self.LS_Induction = {}  # Style for inductOgram plots
        self.LW_Induction = {}  # Widths for inductOgram plots
        self.MW_Induction = None  # Marker size to use for induction scatter plots
        self.MS_Induction = None  # Marker style for induction scatter plots

        self.wedgeAngle_deg = None  # Angular size of wedge diagrams in degrees
        self.LW_wedge = None  # Linewidth in pt for minor boundaries in wedge diagrams
        self.LW_wedgeMajor = None  # Linewidth in pt for major layer boundaries in wedge diagrams


""" Figure label settings """
class FigLblStruct:
    # Do not set any toggles below. These values are assigned in PlanetProfile.GetConfig.
    # Unlike other structs, the labels set below ARE the ones used in plots, so the user
    # only needs to bother with the toggles in the config file.
    def __init__(self):
        # Label display toggles
        self.NEGATIVE_UNIT_POWERS = True  # Whether to use negative powers for units in latex tables, or instead a backslash.
        self.NAN_FOR_EMPTY = False  # Whether to use nan (or -) for empty layer parameters that were not calculated or not present.
        self.w_IN_WTPCT = False  # Whether to print salinities in wt% (or g/kg) in tables
        self.x_IN_WTPCT = True  # Whether to print silicate/core mass fractions in wt% (or g/kg) in tables
        self.qSURF_IN_mW = True  # Whether to print qSurf in mW/m^2 (or W/m^2)
        self.phi_IN_VOLPCT = False  # Whether to print porosity (phi) in vol% (or unitless volume fraction)

        # Wedge diagram labels
        self.wedgeTitle = 'interior structure diagram'
        self.wedgeRadius = r'Radius ($\mathrm{km}$)'

        # Inductogram labels
        self.plotTitles = ['Amplitude $A$', '$B_x$ component', '$B_y$ component', '$B_z$ component']
        self.fLabels = ['Amp', 'Bx', 'By', 'Bz']
        self.compEnd = ''
        self.phaseTitle = r'Phase delay $\upphi$ ($^\circ$)'
        self.Dlabel = r'Ocean thickness $D$ ($\si{km}$)'
        self.TbLabel = r'Ice bottom temp $T_b$ ($\si{K}$)'
        self.iceThickLbl = r'Ice shell thickness ($\si{km}$)'
        self.oceanTempLbl = r'Mean ocean temp ($\si{K}$)'
        self.wScale = 'log'
        self.sigScale = 'log'
        self.Dscale = 'log'

        # Induction parameter-dependent settings
        self.phaseSpaceTitle = None
        self.inductionTitle = None
        self.inductCompareTitle = None
        self.sigLims = None
        self.Dlims = None
        self.legendTexc = None
        self.yLabelInduct = None
        self.yScaleInduct = None

        # Unit-dependent labels set by SetUnits
        self.rhoUnits = None
        self.sigUnits = None
        self.wUnits = None
        self.xUnits = None
        self.fluxUnits = None
        self.wDiv = None
        self.xDiv = None
        self.phiMult = None
        self.qMult = None
        self.NA = None
        self.sigLabel = None
        self.wLabel = None
        self.rhoLabel = None
        self.phiLabel = None


    def SetUnits(self):
        """ Set labels that depend on different selections for units """

        # Plot and table units
        if self.NEGATIVE_UNIT_POWERS:
            self.rhoUnits = r'kg\,m^{-3}'
            self.sigUnits = r'S\,m^{-1}'
            self.wUnits = r'g\,kg^{-1}'
            self.xUnits = r'g\,kg^{-1}'
            self.fluxUnits = r'W\,m^{-2}'
        else:
            self.rhoUnits = r'kg/m^3'
            self.sigUnits = r'S/m'
            self.wUnits = r'g/kg'
            self.xUnits = r'g/kg'
            self.fluxUnits = r'W/m^2'
        if self.w_IN_WTPCT:
            self.wUnits = r'wt\%'
            self.wDiv = 10
        else:
            self.wDiv = 1
        if self.x_IN_WTPCT:
            self.xUnits = r'wt\%'
            self.xDiv = 10
        else:
            self.xDiv = 1
        if self.phi_IN_VOLPCT:
            self.phiUnits = r'~(\si{vol\%})'
            self.phiMult = 100
        else:
            self.phiUnits = ''
            self.phiMult = 1
        if self.qSURF_IN_mW:
            self.fluxUnits = 'm' + self.fluxUnits
            self.qMult = 1e3
        else:
            self.qMult = 1
        # What to put for NA or not calculated numbers
        if self.NAN_FOR_EMPTY:
            self.NA = r'\num{nan}'
        else:
            self.NA = '-'
        self.sigLabel = r'Mean conductivity $\overline{\sigma}$ ($\si{' + self.sigUnits + '}$)'
        self.wLabel = r'Salinity $w$ ($\si{' + self.wUnits + '}$)'
        self.rhoLabel = r'Silicate density $\rho_\mathrm{sil}$ ($\si{' + self.rhoUnits + '}$)'
        self.phiLabel = r'Seafloor porosity $\phi_\mathrm{sil}$ ' + self.phiUnits

    def singleComp(self, comp):
        # Set a tag to append to titles in the event all of what we're plotting
        # has a single composition, for additional clarity.
        self.compEnd = f', \ce{{{comp}}} ocean'

    def setInduction(self, bodyname, IndParams, Texc_h):
        # Set titles, labels, and axis settings pertaining to inductogram plots
        self.phaseSpaceTitle = f'\\textbf{{{bodyname} interior phase space}}'
        self.inductionTitle = f'\\textbf{{{bodyname} induction response{self.compEnd}}}'
        self.inductCompareTitle = f'\\textbf{{{bodyname} induction response on different axes{self.compEnd}}}'

        self.sigLims = [10**IndParams.sigmaMin[bodyname], 10**IndParams.sigmaMax[bodyname]]
        self.Dlims = [10**IndParams.Dmin[bodyname], 10**IndParams.Dmax[bodyname]]
        self.legendTexc = np.array([f'{T_h:.2f} h' for T_h in Texc_h])

        if IndParams.inductOtype != 'sigma':
            if IndParams.inductOtype == 'Tb':
                self.yLabelInduct = self.TbLabel
                self.yScaleInduct = 'linear'
            elif IndParams.inductOtype == 'rho':
                self.yLabelInduct = self.rhoLabel
                self.yScaleInduct = 'linear'
            elif IndParams.inductOtype == 'phi':
                self.yLabelInduct = self.phiLabel
                self.yScaleInduct = 'log'


""" Figure size settings """
class FigSizeStruct:
    # Do not set any values below. All values are assigned in PlanetProfile.GetConfig.
    def __init__(self):
        self.vsP = None
        self.vsR = None
        self.vperm = None
        self.vgsks = None
        self.vseis = None
        self.vhydro = None
        self.vgrav = None
        self.vmant = None
        self.vcore = None
        self.vpvt4 = None
        self.vpvt6 = None
        self.vwedg = None
        self.phaseSpaceSolo = None
        self.phaseSpaceCombo = None
        self.induct = None
        self.inductCombo = None


""" Miscellaneous figure options """
class FigMiscStruct:
    # Do not set any values below. All values are assigned in PlanetProfile.GetConfig.
    def __init__(self):
        # General figure options
        self.figFormat = 'pdf'
        self.dpi = 300  # Resolution in dots per inch for raster images (.png). Ignored for vector images (.pdf, .eps)
        self.xtn = '.' + self.figFormat  # Figure file extension. Good options are .eps, .pdf, and .png
        self.defaultFontName = 'STIXGeneral'  # Default font variables--STIX is what is used in Icarus journal submissions
        self.defaultFontCode = 'stix'  # Code name for default font needed in some function calls
        self.backupFont = 'Times New Roman'  # Backup font that looks similar to STIX that most users are likely to have

        # Wedge diagrams
        self.IONOSPHERE_IN_WEDGE = False  # Whether to include specified ionosphere in wedge diagram
        self.WEDGE_ICE_TICKS = False  # Whether to print ticks for ice shell, which usually overlap with the body outer radius
        self.DRAW_IONOS_BOUND = False  # Whether to draw a boundary line around the ionosphere
        self.DRAW_CONVECTION_BOUND = False  # Whether to draw a boundary line between convecting and conducting regions
        self.DRAW_POROUS_BOUND = False  # Whether to draw a boundary line between porous and non-porous materials
        self.DRAW_FeS_BOUND = False  # Whether to draw a boundary line between Fe and FeS in the core

        # Inductogram phase space plots
        self.DARKEN_SALINITIES = False  # Whether to match hues to the colorbar, but darken points based on salinity, or to just use the colorbar colors.
        self.NORMALIZED_SALINITIES = False  # Whether to normalize salinities to absolute concentrations relative to the saturation limit for each salt
        self.NORMALIZED_TEMPERATURES = False  # Whether to normalize ocean mean temperatures to specified maxima and minima for the colormap
        # Inductograms
        self.PLOT_INDUCT_SURF = False  # Plot surfaces or the default, contours
        self.PLOT_V2021 = True  # Mark the selected ocean/conductivity combos used in Vance et al. 2021
        # Excitation spectra
        self.DO_PER = True  # Convert frequency axes to periods for FFT plots

        # Legends
        self.LEGEND = True  # Whether to plot legends
        self.REFS_IN_LEGEND = True  # Hydrosphere plot: Whether to include reference profiles in legend
        self.hydroLegendBox = None  # Hydrosphere plot: Bounding box for where to place legends. Values are x, y, dx, dy in fractions of the figure size, where x and y are for the bottom-left corner of the box.
        self.hydroLegendPos = None  # Hydrosphere plot: Where to place legends within bounding box
        self.wedgeLegendPos = None  # Wedge diagram: Where in axes added at right to place legend

        # Table printout settings
        self.PRINT_BULK = True  # Whether to print bulk body properties, like mass and MoI
        self.ALWAYS_SHOW_HP = True  # Whether to force HP ices and clathrates to be shown in DISP_* outputs to the terminal, even when none are present.
        self.ALWAYS_SHOW_PHI = True  # Whether to force porosity printout in DISP_* outputs
        self.LATEX_VLINES = False  # Whether to include vertical lines at table edges and between entries. Some journals do not allow them.
        self.LATEX_HLINES = False  # Whether to print horizontal lines between table entries.
        self.HF_HLINES = True  # Whether to print horizontal lines at head and foot of latex tables

        self.cLabelSize = 10  # Font size in pt for contour labels
        self.cLabelPad = 5  # Padding in pt to set beside contour labels
        self.cLegendOpacity = 1.0  # Opacity of legend backgrounds in contour plots.

        self.cbarSpace = 0.5  # Amount of whitespace in inches to use for colorbars
        self.cbarSize = '5%'  # Description of the size of colorbar to use with make_axes_locatable
        self.cbarHeight = 0.6  # Fraction of total figure height to use for colorbar size
        self.cbarPad = 0.25  # Padding in pt to use for colorbars
        self.extraPad = self.cbarSpace * 0.8  # Amount of extra padding to apply to secondary colorbars
        self.cbarFmt = '%.1f'  # Format string to use for colorbar units
        self.nCbarPts = 80  # Number of points to use for drawing colorbar gradient
        self.cbarBottom = (1 - self.cbarHeight - self.cbarPad*2/72)/2  # Fraction of total figure height to use for bottom edge of colorbar

        # Latex settings
        self.TEX_INSTALLED = None
        self.fontFamily = None
        self.latexPackages = []
        self.latexPreamble = None

    def SetLatex(self):
        # Latex settings
        plt.rcParams['font.family'] = self.fontFamily  # Choose serif font for figures to best match math mode variables in body text
        plt.rcParams['font.serif'] = self.defaultFontName  # Set plots to use the default font
        # Check if Latex executable is on the path so we can use backup options if Latex is not installed
        if shutil.which('latex'):
            plt.rcParams['text.usetex'] = True  # Use Latex interpreter to render text on plots
            # Load in font and option packages in Latex
            plt.rcParams['text.latex.preamble'] = f'\\usepackage{{{self.defaultFontCode}}}' \
                                                  + ''.join(self.latexPackages)
            self.TEX_INSTALLED = True
        else:
            print('A LaTeX installation was not found. Some plots may have fallback options in labels.')
            plt.rcParams['font.serif'] += ', ' + self.backupFont  # Set plots to use the default font if installed, or a backup if not
            plt.rcParams['mathtext.fontset'] = self.defaultFontCode
            self.TEX_INSTALLED = False

        packageCmds = r'\n        '.join(self.latexPackages)
        self.latexPreamble = f"""For table setup, add this to your document preamble:
        {packageCmds}
        \sisetup{{group-separator={{\,}}, group-minimum-digits={{5}}, group-digits={{integer}}}}
        """


""" Global EOS list """
class EOSlistStruct:
    def __init__(self):
        pass
    loaded = {}  # Dict listing the loaded EOSs. Since we define this attribute outside of __init__, it will be common to all EOSlist structs when set.
    ranges = {}  # Dict listing the P, T ranges of the loaded EOSs.


""" Physical constants """
class ConstantsStruct:
    def __init__(self):
        """ General physical constants """
        self.G = 6.673e-11  # "Big G" gravitational constant, m^3/kg/s
        self.bar2GPa = 1.01325e-4  # Multiply by this to convert pressure from bars to GPa
        self.bar2MPa = 1.01325e-1  # Same but for MPa
        self.erg2J = 1e-7  # Multiply by this to convert from ergs to joules
        self.T0 = 273.15  # The Celsius zero point in K.
        self.P0 = 0.101325  # One standard atmosphere in MPa
        self.R = 8.314  # Ideal gas constant in J/mol/K
        self.mu0 = 4e-7*np.pi
        self.stdSeawater_ppt = 35.16504  # Standard Seawater salinity in g/kg (ppt by mass)
        self.sigmaH2O_Sm = 1e-5  # Assumed conductivity of pure water (only used when wOcean_ppt == 0)
        self.mMgSO4_gmol = 120.4  # Molecular mass of MgSO4 in g/mol
        self.mNaCl_gmol = 58.44  # Molecular mass of NaCl in g/mol
        self.mNH3_gmol = 17.03  # Molecular mass of NH3 in g/mol
        self.mH2O_gmol = 18.02  # Molecular mass of pure water in g/mol
        self.mCH4_gmol = 16.04  # Molecular mass of methane in g/mol
        self.mCO2_gmol = 44.01  # Molecular mass of CO2 in g/mol
        self.mClathGas_gmol = self.mCH4_gmol + 5.75 * self.mH2O_gmol  # Molecular mass of clathrate unit cell
        self.clathGasFrac_ppt = 1e3 * self.mCH4_gmol / self.mClathGas_gmol  # Mass fraction of gases trapped in clathrates in ppt
        self.QScore = 1e4  # Fixed QS value to use for core layers if not set in PPBody.py file
        self.kThermWater_WmK = 0.55  # Fixed thermal conductivity of liquid water in W/(m K)
        self.kThermSil_WmK = 4.0  # Fixed thermal conductivity of silicates in W/(m K)
        self.kThermFe_WmK = 33.3  # Fixed thermal conductivity of core material in W/(m K)
        self.phaseClath = 30  # Phase ID to use for (sI methane) clathrates. Must be larger than 6 to keep space for pure ice phases
        self.phaseSil = 50  # Phase ID for silicates (non-porous or liquid-filled). We add the phase ID for ice phases in pores, e.g. 56 is silicates with interstitial ice VI.
        self.phaseFe = 100  # Phase ID to use for liquid iron core material
        self.phaseFeSolid = 105  # Phase ID for solid iron core material
        self.phaseFeS = 110  # Phase ID for liquid FeS
        self.phaseFeSsolid = 115  # Phase ID for solid FeS
        self.sigmaClath_Sm = 5e-5  # Roughly fixed conductivity in the range 260-281 K, as reported by Stern et al. (2021): https://doi.org/10.1029/2021GL093475
        self.sigmaCO2Clath_Sm = 6.5e-4  # Also from Stern et al. (2021), at 273 K and 25% gas-filled porosity
        self.EactCO2Clath_kJmol = 46.5  # Also from Stern et al. (2021)
        # Initialize activation energies and melting point viscosities, for use in convection calculations
        self.Eact_kJmol, self.etaMelt_Pas = (np.ones(self.phaseClath+1) * np.nan for _ in range(2))
        self.Eact_kJmol[1:7] = np.array([59.4, 76.5, 127, np.nan, 136, 110])  # Activation energy for diffusion of ice phases Ih-VI in kJ/mol
        self.Eact_kJmol[self.phaseClath] = 90.0  # From Durham et al. (2003), at 50 and 100 MPa and 260-283 K: https://doi.org/10.1029/2002JB001872
        self.etaMelt_Pas[1:7] = np.array([1e14, 1e18, 5e12, np.nan, 5e14, 5e14])  # Viscosity at the melting temperature of ice phases Ih-VI in Pa*s. Ice Ih range of 5e13-1e16 is from Tobie et al. (2003), others unknown
        self.etaMelt_Pas[self.phaseClath] = self.etaMelt_Pas[1] * 20  # Estimate of clathrate viscosity 20x that of ice Ih at comparable conditions from Durham et al. (2003): https://doi.org/10.1029/2002JB001872
        self.PminHPices_MPa = 200.0  # Min plausible pressure of high-pressure ices for any ocean composition in MPa
        self.PmaxLiquid_MPa = 2250.0  # Maximum plausible pressure for liquid water oceans
        self.sigmaDef_Sm = 1e-8  # Default minimum conductivity to use for layers with NaN or 0 conductivity
        self.sigmaMin_Sm = 1e-8  # Threshold conductivity below which we set to the default to reduce computational overhead

Constants = ConstantsStruct()
EOSlist = EOSlistStruct()
