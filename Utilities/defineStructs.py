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
import os

# We have to define subclasses first in order to make them instanced to each Planet object
""" Run settings """
class BulkSubstruct():

    def __init__(self):
        self.Tb_K = None  # Temperature at the bottom of the ice I layer (ice-ocean interface when there are no ice III or V underplate layers). Ranges from 238.5 to 261.165 K for ice III transition and 251.165 to 273.16 for melting temp.
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
        self.J2 = None  # Gravitational oblateness shape parameter
        self.C22 = None  # Gravitational elongation shape parameter
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
        self.nRefRho = None  # Number of values for plotting reference density curves (sets resolution)
        self.iSilStart = None  # Hydrosphere index at which to start silicate size search
        self.nSilMax = None  # Fixed max number of steps in silicate layers
        self.nSil = None  # Derived final number of steps in silicate layers
        self.nCore = None  # Fixed number of steps in core layers, if present
        self.nIceIIILitho = 100  # Fixed number of layers to use for ice III when either BOTTOM_ICEIII or BOTTOM_ICEV is True.
        self.nIceVLitho = 100  # Fixed number of layers to use for ice V when BOTTOM_ICEV is True.
        self.nPsHP = 150  # Number of interpolation steps to use for getting HP ice EOS (pressures)
        self.nTsHP = 100  # Number of interpolation steps to use for getting HP ice EOS (temperatures)
        self.nPoros = 10  # Number of steps in porosity to use in geometric series between phiMin and phiMax for porous rock when no core is present


""" Hydrosphere assumptions """
class OceanSubstruct:

    def __init__(self):
        self.comp = None  # Type of dominant dissolved salt in ocean. Options: 'Seawater', 'MgSO4', 'PureH2O', 'NH3', 'NaCl', 'none'
        self.wOcean_ppt = None  # (Absolute) salinity: Mass concentration of above composition in parts per thousand (ppt)
        self.sigmaMean_Sm = None  # Mean conductivity across all ocean layers (linear average, ignoring spherical geometry effects)
        self.sigmaTop_Sm = None  # Conductivity of shallowest ocean layer
        self.deltaP = None  # Increment of pressure between each layer in lower hydrosphere/ocean (sets profile resolution)
        self.deltaT = None  # Step size in K for temperature values used in generating ocean EOS functions. If set, overrides calculations that otherwise use the specified precision in Tb_K to determine this.
        self.koThermI_WmK = 2.21  # Thermal conductivity of ice I at melting temp. Default is from Eq. 6.4 of Melinder (2007), ISBN: 978-91-7178-707-1
        self.dkdTI_WmK2 = -0.012  # Temperature derivative of ice I relative to the melting temp. Default is from Melinder (2007).
        self.sigmaIce_Sm = {'Ih':1e-8, 'II':1e-8, 'III':1e-8, 'V':1e-8, 'VI':1e-8, 'Clath':5e-5}  # Assumed conductivity of ice phases (see Constants.sigmaClath_Sm below)
        self.THydroMax_K = 320  # Assumed maximum ocean temperature for generating ocean EOS functions. For large bodies like Ganymede, Callisto, and Titan, larger values are required.
        self.PHydroMax_MPa = None  # Guessed maximum pressure of the hydrosphere in MPa. Must be greater than the actual pressure, but ideally not by much. Sets initial length of hydrosphere arrays, which get truncated after layer calculations are finished.
        self.MgSO4elecType = 'Vance2018'  # Type of electrical conductivity model to use for MgSO4. Options: 'Vance2018', 'Pan2020'
        self.MgSO4scalingType = 'Vance2018'  # Type of scaling to apply to Larionov and Kryukov model. Options: 'Vance2018', 'LK1984'
        self.MgSO4rhoType = 'Millero'  # Type of water density model to use in Larionov and Kryukov model. Options: 'Millero', 'SeaFreeze'
        self.MgSO4phaseType = 'lookup'  # Type of phase calculation to use for MgSO4. Currently, "Margules" forces a (slow) individual calc for each P, T point, and any other option uses a lookup table like the Perplex EOS functions.
        self.QfromMantle_W = None  # Heat flow from mantle into hydrosphere (calculated from ice thermal profile and applied to mantle)
        self.EOS = None  # Equation of state data to use for ocean layers
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
        self.sigmaPoreMean_Sm = None  # Mean conductivity of pore fluids across all layers with porosity above poreConductThresh (about 5%)
        self.sigmaPorousLayerMean_Sm = None  # Mean conductivity of matrix + pore fluid combined, for all layers with porosity of poreConductThresh
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
        self.rhoTrade_kgm3 = None  # Array of mantle densities for compatible MoIs for core vs. mantle tradeoff plot
        self.mFluids = None  # WIP for tracking loss of fluids along the geotherm -- needs a better name.
        # The below not necessary to be implemented until later (says Steve), these 5 are based on DPS presentation in 2017 – 5 diff models of permeability
        #turn off this plot feature until later- create flag, Use POROSITY flag to turn off these plots
        #perm1 = None
        #perm2 = None
        #perm3 = None
        #perm4 = None
        #perm5 = None


""" Core layers """
class CoreSubstruct:

    def __init__(self):
        self.rhoFe_kgm3 = 8000  # Assumed density of pure iron in kg/m^3
        self.rhoFeS_kgm3 = 5150  # Assumed density of iron sulfide in kg/m^3
        self.rhoMin_kgm3 = 5150  # Assumed minimum possible density for the core in kg/m^3. Sets maximum core size.
        self.rhoPoFeFCC = None  # Density of pyrrhottite plus face-centered cubic iron
        self.sigmaCore_Sm = 1e-16  # Fixed electrical conductivity to apply to core (typically low, to ignore core impacts on induction)
        self.coreEOS = 'sulfur_core_partition_SE15_1pctSulfur.tab'  # Default core EOS to use
        self.EOS = None  # Interpolator functions for evaluating Perple_X EOS model
        self.kTherm_WmK = None  # Constant thermal conductivity to set for a specific body (overrides Constants.kThermFe_WmK)
        # Derived quantities
        self.rho_kgm3 = None  # Core bulk density consistent with assumed mixing ratios of Fe, FeS, etc.
        self.rhoMean_kgm3 = None  # Core bulk density calculated from final MoI match using EOS properties
        self.Rmean_m = None  # Core radius for mean compatible moment of inertia (MoI)
        self.Rrange_m = None  # Core radius range for compatible MoI
        self.Rtrade_m = None  # Array of core radii for compatible MoIs
        #Re Steve- put all mass fraction stuff into a separate file until implemented later- remove from dataStructs.py
        #To implement: possible Meteoritics file/class?
        # 2021-12-30: Judging by usage of various different fractional variables in the literature and in
        # the Matlab code, these x variables should be molar fractions (# this species/total # molecules).
        self.xFeSmeteoritic = None  # CM2 mean from Jarosewich 1990
        self.xFeS = None  # Mass fraction of sulfur in the core
        self.xFeCore = None  # This is the total Fe in Fe and FeS
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
        self.inductMethod = 'Eckhardt1963'  # Type of magnetic induction model to use. Options are "Srivastava1966" or "layer" for layer method and "Eckhardt1963" or "numeric" for numeric method (symmetric only).
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
        self.pMax = 0  # Maximum p to use in asymmetric shape (0 for spherically symmetric).
        self.asymShape = None  # Asymmetric shape to use in induction calculations based on Bulk.J2 and Bulk.C22 values. Only used when inductType = "Srivastava1966".
        self.xInductOgram = None  # Value of variable along x axis of induct-o-gram to use for this model.
        self.yInductOgram = None  # Value of variable along y axis of induct-o-gram to use for this model.
        # Output calculations
        self.Aen = None  # Complex response amplitude of magnetic excitation for dipole moment for each excitation frequency (unitless)
        self.Amp = None  # Amplitude (modulus) of magnetic excitation for spherically symmetric approximation for each excitation frequency (unitless)
        self.phase = None  # Phase delay of magnetic excitation for spherically symmetric approximation for each excitation frequency in degrees
        self.Binm_nT = None  # Induced magnetic moments relative to the body surface in nT


""" Main body profile info--settings and variables """
class PlanetStruct:

    # Require a body name as an argument for initialization; define instance attributes
    def __init__(self, name):
        self.name = name

        self.Bulk = BulkSubstruct()
        self.Do = DoSubstruct()
        self.Steps = StepsSubstruct()
        self.Ocean = OceanSubstruct()
        self.Sil = SilSubstruct()
        self.Core = CoreSubstruct()
        self.Seismic = SeismicSubstruct()
        self.Magnetic = MagneticSubstruct()
        # Settings for GetPfreeze start, stop, and step size.
        # Shrink closer to expected melting pressure to improve run times.
        self.PfreezeLower_MPa = 5  # Lower boundary for GetPfreeze to search for ice Ih phase transition
        self.PfreezeUpper_MPa = 230  # Upper boundary for GetPfreeze to search for ice Ih phase transition
        self.PfreezeRes_MPa = 0.1  # Step size in pressure for GetPfreeze to use in searching for phase transition

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
        self.kTherm_WmK = None  # Thermal conductivity of each layer in W/(m K)
        self.Htidal_Wm3 = None  # Tidal heating rate of each layer in W/m^3
        self.Ppore_MPa = None  # Pressure of fluids assumed to occupy full pore space
        self.rhoMatrix_kgm3 = None  # Mass density of matrix material (rock or ice)
        self.rhoPore_kgm3 = None  # Mass density of pore material (typically ocean water)
        # Individual calculated quantities
        self.Mtot_kg = None  # Total calculated mass selected from MoI matching
        self.zb_km = None  # Thickness of outer ice shell/depth of ice-ocean interface in km
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
        self.eLid_m = None  # Thickness of conducting stagnant lid layer
        self.eLidIII_m = None  # Same as above but for ice III underplate layers.
        self.eLidV_m = None  # Same as above but for ice V underplate layers.
        self.deltaTBL_m = None  # Thickness of lower thermal boundary layer when Htidal = 0 in the ice
        self.deltaTBLIII_m = None  # Same as above but for ice III underplate layers.
        self.deltaTBLV_m = None  # Same as above but for ice V underplate layers.
        self.RaConvect = None  # Rayleigh number of putative convective layer within the ice I layers. If this number is below Constants.RaCrit, convection does not occur.
        self.RaConvectIII = None  # Same as above but for ice III underplate layers.
        self.RaConvectV = None  # Same as above but for ice V underplate layers.
        self.MH2O_kg = None  # Total mass of water molecules contained in ice, liquid, and pore spaces
        self.Mrock_kg = None  # Total mass contained in silicate rock (just the matrix, when layers are porous)
        self.Mcore_kg = None  # Total mass contained in iron core material
        self.Mice_kg = None  # Total mass contained in all ice phases, including gas trapped in clathrates
        self.Msalt_kg = None  # Total mass contained in ocean solute
        self.MporeSalt_kg = None  # Total mass contained in ocean solute in pore spaces
        self.Mocean_kg = None  # Total mass contained in ocean fluids, including H2O and salts, excluding pores
        self.Mfluid_kg = None  # Sum of the masses in ocean and pore spaces
        self.MporeFluid_kg = None  # Total mass contained in pore fluids, including H2O and salts
        self.Mclath_kg = None  # Total mass of clathrate layers
        self.MclathGas_kg = None  # Total mass of non-water molecules trapped in clathrates
        self.index = None  # Numeric indicator to aid in progress info in multi-model runs


""" Params substructs """
# Construct filenames for data, saving/reloading
class DataFilesSubstruct:
    def __init__(self, datPath, saveBase, inductBase=None):
        if inductBase is None:
            inductBase = saveBase
        self.path = datPath
        self.fName = os.path.join(self.path, saveBase)
        self.saveFile = self.fName + '.txt'
        self.mantCoreFile = self.fName + '_mantleCore.txt'
        self.permFile = self.fName + '_mantlePerm.txt'
        self.fNameInduct = os.path.join(self.path, 'induction', saveBase)
        self.inductLayersFile = self.fNameInduct + '_inductLayers.txt'
        self.inducedMomentsFile = self.fNameInduct + '_inducedMoments.txt'
        self.fNameInductOgram = os.path.join(self.path, 'induction', inductBase)
        self.inductOgramFile = self.fNameInductOgram + '_inductOgram.mat'
        self.inductOgramSigmaFile = self.fNameInductOgram + '_sigma_inductOgram.mat'


# Construct filenames for figures etc.
class FigureFilesSubstruct:
    def __init__(self, figPath, figBase, xtn):
        self.path = figPath
        self.fName = os.path.join(self.path, figBase)

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
        induct = 'InductOgram'
        sigma = 'InductOgramSigma'
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
        self.induct = self.fName + induct + xtn
        self.sigma = self.fName + sigma + xtn


""" Figure size """
class FigSizeStruct:
    def __init__(self):
        pass


""" Figure color options """
class ColorsStruct:
    def __init__(self):
        pass


""" General parameter options """
class ParamsStruct:
    def __init__(self):
        self.Colors = ColorsStruct()
        self.FigSize = FigSizeStruct()


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
