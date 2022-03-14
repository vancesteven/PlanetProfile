"""
dataStructs: Create empty classes and subclasses for holding body-specific data
Values are typically set to None as defaults; body-specific values should be set in the PPBody.py files.
Optional SWITCHES are in all caps. These typically have a default value of False.

Example usage:
Planet = PlanetStruct('nameOfBody')
Planet.R_m = 1560e3
Planet.Ocean.comp = 'MgSO4'
Planet.Sil.mantleEOS = 'CV3hy1wt_678_1.tab'
Planet.Do.Fe_CORE = False
"""
import numpy as np

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


""" Layer step settings """
class StepsSubstruct:

    def __init__(self):
        self.nIceI = None  # Fixed number of steps in outermost ice shell
        self.nClath = None  # Fixed number of steps in clathrates
        self.nHydroMax = None  # Derived working length of hydrosphere layers, gets truncated after layer calcs
        self.nOceanMax = None  # Derived working length of ocean layers, also truncated after layer calcs
        self.nHydro = None  # Derived final number of steps in hydrosphere
        self.nIbottom = None  # Derived number of clathrate + ice I layers
        self.nIIIbottom = None  # Derived number of clathrate + ice I + ice III layers
        self.nSurfIce = None  # Derived number of outer ice layers (above ocean) -- sum of nIceI, nClath, nIceIIILitho, nIceVLitho
        self.nStepsRefRho = None  # Number of values for plotting reference density curves (sets resolution)
        self.iSilStart = None  # Hydrosphere index at which to start silicate size search
        self.nSilMax = None  # Fixed max number of steps in silicate layers
        self.nSil = None  # Derived final number of steps in silicate layers
        self.nCore = None  # Fixed number of steps in core layers, if present
        self.nIceIIILitho = 100  # Fixed number of layers to use for ice III when either BOTTOM_ICEIII or BOTTOM_ICEV is True.
        self.nIceVLitho = 100  # Fixed number of layers to use for ice V when BOTTOM_ICEV is True.
        self.nPsHP = 150  # Number of interpolation steps to use for getting HP ice EOS (pressures)
        self.nTsHP = 100  # Number of interpolation steps to use for getting HP ice EOS (temperatures)
        self.nPoros = 10  # Number of steps in porosity to use in geometric series between phiMin and phiMax


""" Hydrosphere assumptions """
class OceanSubstruct:

    def __init__(self):
        self.comp = None  # Type of dominant dissolved salt in ocean. Options: 'Seawater', 'MgSO4', 'NH3', 'NaCl'
        self.wOcean_ppt = None  # (Absolute) salinity: Mass concentration of above composition in parts per thousand (ppt)
        self.deltaP = None  # Increment of pressure between each layer in lower hydrosphere/ocean (sets profile resolution)
        self.deltaT = None  # Step size in K for temperature values used in generating ocean EOS functions. If set, overrides calculations that otherwise use the specified precision in Tb_K to determine this.
        self.koThermI_WmK = 2.21  # Thermal conductivity of ice I at melting temp. Default is from Eq. 6.4 of Melinder (2007), ISBN: 978-91-7178-707-1
        self.dkdTI_WmK2 = -0.012  # Temperature derivative of ice I relative to the melting temp. Default is from Melinder (2007).
        self.sigmaIce_Sm = 1e-8  # Assumed conductivity of ice layers
        self.THydroMax_K = 320  # Assumed maximum ocean temperature for generating ocean EOS functions. For large bodies like Ganymede, Callisto, and Titan, larger values are required.
        self.PHydroMax_MPa = None  # Guessed maximum pressure of the hydrosphere in MPa. Must be greater than the actual pressure, but ideally not by much. Sets initial length of hydrosphere arrays, which get truncated after layer calculations are finished.
        self.MgSO4elecType = 'Vance2018'  # Type of electrical conductivity model to use for MgSO4. Options: 'Vance2018', 'Pan2020'
        self.MgSO4scalingType = 'Vance2018'  # Type of scaling to apply to Larionov and Kryukov model. Options: 'Vance2018', 'LK1984'
        self.MgSO4rhoType = 'Millero'  # Type of water density model to use in Larionov and Kryukov model. Options: 'Millero', 'SeaFreeze'
        self.QfromMantle_W = None  # Heat flow from mantle into hydrosphere (calculated from ice thermal profile and applied to mantle)
        self.EOS = None  # Equation of state data to use for ocean layers
        self.surfIceEOS = {}  # Equation of state data to use for surface ice layers
        self.iceEOS = {}  # Equation of state data to use for ice layers within the ocean


""" Silicate layers """
class SilSubstruct:

    def __init__(self):
        self.alphaPeff = 0.95  # Scaling factor by which the hydrostatic pressure of pore-filling fluids reduces pore closure pressure. Value taken from the high end of Vitovtova et al. (2014).
        self.sigmaSil_Sm = 1e-16  # Assumed conductivity of silicate rock
        self.Qrad_Wkg = 0  # Average radiogenic heating rate for silicates in W/kg.
        self.Htidal_Wm3 = 0  # Average tidal heating rate for silicates in W/m^3.
        self.HtidalMin_Wm3 = 1e-16  # Minimum average tidal heating rate in silicates in W/m^3 above zero to start with for finding MoI match with no core.
        self.HtidalMax_Wm3 = 1e-6  # Maximum average tidal heating to stop MoI search
        self.deltaHtidal_logUnits = 1/3  # Step size by which to increment Htidal_Wm3 for finding MoI match with no core.
        self.kTherm_WmK = None  # Constant thermal conductivity to set for a specific body (overrides Constants.kThermSil_WmK)
        """ Porosity parameters """
        self.phiRockMax_frac = None  # Porosity (void fraction) of the rocks in vacuum. This is the expected value for core-less bodies, and porosity is modeled for a range around here to find a matching MoI. For bodies with a core, this is a fixed value for rock porosity at P=0.
        self.phiRangeMult = 8  # Factor by which to divide the distance from user-defined phiRockMax_frac to 0 and 1 to obtain the search range for phi
        self.Pclosure_MPa = 350  # Pressure threshold in MPa beyond which pores in silicates shut completely and porosity drops to zero, for use in Han et al. (2014) model. See Saito et al. (2016) for evidence of values up to ~750 MPa: https://doi.org/10.1016/j.tecto.2016.03.044
        self.porosType = 'Han2014'  # Porosity model to apply for silicates. Options are 'Han2014', 'Vitovtova2014', 'Chen2020'.
        self.poreH2Orho_kgm3 = 1023  # Arbitrary ocean water density used to calculate filled-pore rock densities. This value is that used in Vance et al. (2018).
        # J values for exponent of pore property / overlaying matrix combination as in Yu et al. (2016): http://dx.doi.org/10.1016/j.jrmge.2015.07.004
        # Values X combine as Xtot^J = Xmatrix^J * (1 - phi) + Xpore^J * phi, where phi is the porosity (volume fraction) of pore space filled with the secondary material.
        # Values of J typically range from -1 to 1, where 1 is a direct mean, e.g. Xtot = Xrock*(1-phi) + Xpore*phi and -1 is a geometric mean, e.g. Xtot = Xrock || Xpore = 1/((1-phi)/Xrock + phi/Xpore).
        # J can range outside -1 and 1 for VP and VS, according to Yu et al.
        self.Jrho = 1  # Mass density, directly summative
        self.JkTherm = 1  # Should add like conductors in parallel
        # For the following 4 properties, Yu et al. (2016) have values for rock combinations that vary around these marks, but they vary quite a lot. If the values are known
        # for the specific assumed composition of the silicate layer, those values should be used instead for the particular body being simulated.
        self.JKS = 0.35
        self.JGS = 0.35
        self.JVP = 0.75
        self.JVS = 0.85
        self.Jalpha = 1  # Volume increase with temperature should be summative
        self.JCp = 1  # Heat capacity is summative, but since porosity is volume fraction the inputs need to be scaled by mass density
        self.Jsigma = 1  # Pore fluid and matrix provide parallel electrically conducting pathways as well
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
        self.peaks_Hz = None  # Frequencies in Hz of peaks in Fourier spectrum of magnetic excitations
        self.fOrb_radps = None  # Angular frequency of orbital motion of a moon around its parent planet in radians per second
        self.ionosBounds_m = None  # Upper altitude cutoff for ionosphere layers in m. Omit the surface (don't include 0 in the list).
        self.sigmaIonosPedersen_Sm = None  # Pedersen conductivity for ionospheric layers in S/m. Length must match ionosBounds_m.


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
        self.rSigChange_m = None  # Radii of outer boundary of each conducting layer in m (i.e., radii where sigma changes)
        self.MLayer_kg = None  # Mass of each layer in kg
        self.kTherm_WmK = None  # Thermal conductivity of each layer in W/(m K)
        self.Htidal_Wm3 = None  # Tidal heating rate of each layer in W/m^3
        self.Ppore_MPa = None  # Pressure of fluids assumed to occupy full pore space
        self.rhoMatrix_kgm3 = None  # Mass density of matrix material (rock or ice)
        self.rhoPore_kgm3 = None  # Mass density of pore material (ocean water)
        # Individual calculated quantities
        self.zb_km = None  # Thickness of outer ice shell/depth of ice-ocean interface in km
        self.zClath_m = None  # Thickness of clathrate layer in surface ice in m
        self.Pb_MPa = None  # Pressure at ice-ocean interface in MPa
        self.PbI_MPa = None  # Pressure at bottom of surface ice I/clathrate layer in MPa
        self.PbIII_MPa = None  # Pressure at ice III/ice V transition in MPa, only used when BOTTOM_ICEIII or BOTTOM_ICEV is True
        self.PbV_MPa = None  # Pressure at bottom of ice V layer in MPa, only used when BOTTOM_ICEV is True
        self.PbClathMax_MPa = None  # Max pressure at the bottom of a clathrate lid. Actual PbClath may be reduced in convection calculations
        self.TclathTrans_K = None  # Temperature at the transition from clathrates to ice I
        self.CMR2mean = None  # Mean value of axial moment of inertia that is consistent with profile core/mantle trades
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


""" Params substructs """
# Construct filenames for data, saving/reloading
class DataFilesSubstruct:
    def __init__(self, fName):
        self.saveFile = fName + '.txt'
        self.mantCoreFile = fName + '_mantleCore.txt'
        self.permFile = fName + '_mantlePerm.txt'


# Construct filenames for figures etc.
class FigureFilesSubstruct:
    def __init__(self, fName, xtn):
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
        # Construct Figure Filenames
        self.vwedg = fName + vwedg + xtn
        self.vsP = fName + vsP + xtn
        self.vsR= fName + vsR + xtn
        self.vperm = fName + vperm + xtn
        self.vgsks = fName + vgsks + xtn
        self.vseis = fName + vseis + xtn
        self.vhydro = fName + vhydro + xtn
        self.vgrav = fName + vgrav + xtn
        self.vmant = fName + vmant + xtn
        self.vcore = fName + vcore + xtn
        self.vpvt4 = fName + vpvt4 + xtn
        self.vpvt6 = fName + vpvt6 + xtn


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
        self.stdSeawater_ppt = 35.16504  # Standard Seawater salinity in g/kg (ppt by mass)
        self.sigmaH2O_Sm = 1e-5  # Assumed conductivity of pure water (only used when wOcean_ppt == 0)
        self.mMgSO4_gmol = 120.4  # Molecular mass of MgSO4 in g/mol
        self.mNaCl_gmol = 58.44  # Molecular mass of NaCl in g/mol
        self.mNH3_gmol = 17.03  # Molecular mass of NH3 in g/mol
        self.mH2O_gmol = 18.02  # Molecular mass of pure water in g/mol
        self.QScore = 1e4  # Fixed QS value to use for core layers if not set in PPBody.py file
        self.kThermWater_WmK = 0.55  # Fixed thermal conductivity of liquid water in W/(m K)
        self.kThermSil_WmK = 4.0  # Fixed thermal conductivity of silicates in W/(m K)
        self.kThermFe_WmK = 33.3  # Fixed thermal conductivity of core material in W/(m K)
        self.phaseClath = 30  # Phase ID to use for (sI methane) clathrates. Must be larger than 6 to keep space for pure ice phases
        self.phaseSil = 50  # Phase ID for silicates
        self.phaseFe = 100  # Phase ID to use for iron core material
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

Constants = ConstantsStruct()
