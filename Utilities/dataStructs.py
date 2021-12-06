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
class PlanetStruct:
    # Require a body name as an argument for initialization
    def __init__(self, name):
        self.name = name

    # Settings for GetPfreeze start, stop, and step size.
    # Shrink closer to expected melting pressure to improve run times.
    PfreezeLower_MPa = 0.1  # Lower boundary for GetPfreeze to search for phase transition
    PfreezeUpper_MPa = 300  # Upper boundary for GetPfreeze to search for phase transition
    PfreezeRes_MPa = 0.1  # Step size in pressure for GetPfreeze to use in searching for phase transition

    """ Derived quantities (assigned during PlanetProfile runs) """
    # Layer arrays
    phase = None  # Phase of the layer input as an integer: ocean=0, ice I through VI are 1 through 6, clathrate=30, silicates=50, iron=100.
    r_m = None  # Distance from center of body to the outer bound of current layer in m
    z_m = None  # Distance from surface of body to the outer bound of current layer in m
    T_K = None  # Temperature of each layer in K
    P_MPa = None  # Pressure at top of each layer in MPa
    rho_kgm3 = None  # Density of each layer in kg/m^3
    g_ms2 = None  # Gravitational acceleration at top of each layer, m/s^2
    Cp_JkgK = None  # Heat capacity at constant pressure for each layer's material in J/kg/K
    alpha_pK = None  # Thermal expansivity of layer material in hydrosphere in K^-1
    phi_frac = None  # Porosity of each layer's material as a fraction of void/solid
    sigma_Sm = None  # Electrical conductivity (sigma) in S/m of each conducting layer
    MLayer_kg = None  # Mass of each layer in kg
    rSigChange_m = None  # Radii of outer boundary of each conducting layer in m (i.e., radii where sigma changes)
    # Individual calculated quantities
    zb_km = None  # Thickness of outer ice shell/depth of ice-ocean interface in km in accordance with Vance et al. (2014)
    zClath_m = None  # Thickness of clathrate layer at body surface in m
    Pb_MPa = None  # Pressure at ice-ocean interface in MPa
    PbClath_MPa = None  # Pressure at bottom of clathrate layer in MPa
    PbI_MPa = None  # Pressure at bottom of ice I layer in MPa
    PbIII_MPa = None  # Pressure at ice III/ice V transition in MPa, only used when BOTTOM_ICEV is True
    CMR2mean = None  # Mean value of axial moment of inertia that is consistent with profile core/mantle trades
    # Q_Wm2 = None  # ??? WAIT UNTIL IMPLEMENT heat flux at ice shell

    """ Run settings """
    class Bulk:
        """ Bulk planetary settings """
        Tb_K = None  # Temperature at the bottom of the ice I layer (ice-ocean interface when there are no ice III or V underplate layers)
        rho_kgm3 = None  # Bulk density in kg/m^3
        R_m = None  # Mean body outer radius in m
        M_kg = None  # Total body mass in kg
        Tsurf_K = None  # Surface temperature in K
        Psurf_MPa = None  # Surface pressure in MPa
        Cmeasured = None  # Axial moment of inertia C/MR^2, dimensionless
        Cuncertainty = None  # Uncertainty (std dev) of C/MR^2 (used to constrain models via consistency within the uncertainty), dimensionless
        phiSurface = None  # Scaling value for the ice porosity at the surface (void fraction): falls within a range of 0 and 1. 0 is completely non-porous and larger than 0.2 is rare. From Han et al. (2014)
        clathMaxDepth = None  # Fixed limit for thickness of clathrate layer in m
        TbIII_K = None  # Temperature at bottom of ice III underplate layer in K
        TbV_K = None  # Temperature at bottom of ice V underplate layer in K

    """ Runtime flags """
    class Do:
        Fe_CORE = False  # Whether to model an iron core for this body
        POROUS_ICE = False  # Whether to model porosity in ice
        CLATHRATE = False  # Whether to model clathrates
        NO_H2O = False  # Whether to model waterless worlds (like Io)
        BOTTOM_ICEIII = False  # Whether to allow Ice III between ocean and ice I layer, when ocean temp is set very low- default is that this is off, can turn on as an error condition
        BOTTOM_ICEV = False  # Same as above but also including ice V. Takes precedence (forces both ice III and V to be present).
        NO_ICEI_CONVECTION = False  # Whether to suppress convection in the ice I layer - if True, checks Rayleigh number to see if convection conditions are met
        FORCE_ICEI_CONVECTION = False  # Whether to force convection in the ice I layer- if True, doesn’t check Rayleigh number, just assumes convection conditions are met
        ALLOW_NEG_ALPHA = False  # Whether to permit modeling of a Melosh et. al. layer with negative thermal expansivity
        MANTLE_HYDRO_THERMAL_EQ = False  # Whether to set thermal equilibrium between mantle and hydrosphere, where the hydrosphere is not gaining external heat via tidal forcing or radiation
        POROUS_ROCK = False  # Whether to model silicates as porous
        P_EFFECTIVE = False  # Effective pressure due to presence of water in pores (modeled as lithostatic minus hydrostatic pressure).
        IONOS_ONLY = False  # Whether to ignore conducting layers within the body and model magnetic induction happening only in the ionosphere
        TAUP_SEISMIC = False  # Whether to make TauP model files and some basic plots using obspy.taup
        IONOS_ONLY = False  # Whether to model induction happening only in the ionosphere, as if ocean were totally frozen out

    """ Layer step settings """
    class Steps:
        nIceI = None  # Fixed number of steps in outermost ice shell
        nClath = None  # Fixed number of steps in clathrates
        nHydroMax = None  # Derived working length of hydrosphere layers, gets truncated after layer calcs
        nHydro = None  # Derived final number of steps in hydrosphere
        nIbottom = None  # Derived number of clathrate + ice I layers
        nIIIbottom = None  # Derived number of clathrate + ice I + ice V layers
        nSurfIce = None  # Derived number of outer ice layers (above ocean) -- sum of nIceI, nClath, nIceIIILitho, nIceVLitho
        nStepsRefRho = None  # Number of values for plotting reference density curves (sets resolution)
        nSil = None  # Fixed number of steps in silicate layers
        nCore = None  # Fixed number of steps in core layers, if present
        nIceIIILitho = 5  # Fixed number of layers to use for ice III when either BOTTOM_ICEIII or BOTTOM_ICEV is True.
        nIceVLitho = 5  # Fixed number of layers to use for ice V when BOTTOM_ICEV is True.

    class Ocean:
        """ Hydrosphere assumptions """
        comp = None  # Type of dominant dissolved salt in ocean. Options: 'Seawater', 'MgSO4', 'NH3', 'NaCl'
        wOcean_ppt = None  # Salinity: Concentration of above salt in parts per thousand (ppt)
        deltaP = None  # Increment of pressure between each layer in lower hydrosphere/ocean (sets profile resolution)
        PHydroMax_MPa = None  # Guessed maximum pressure of the hydrosphere in MPa. Must be greater than the actual pressure, but ideally not by much. Sets initial length of hydrosphere arrays, which get truncated after layer calculations are finished.
        electrical = 'Vance2018'  # Type of electrical conductivity model to use. Options: 'Vance2018', 'Pan2020'
        QfromMantle_Wm2 = None  # Heat flow from mantle into hydrosphere
        sigmaIce_Sm = 1e-8  # Assumed conductivity of ice layers

    class Sil:
        """ Silicate layer settings """
        kSil_WmK = None  # Thermal conductivity (k) of silicates in W/(mK)
        phiRockMax = None  # Porosity (void fraction) of the rocks at the “seafloor”, where the hydrosphere first comes into contact with rock
        sigmaSil_Sm = 1e-16  # Assumed conductivity of silicate rock
        """ Mantle Equation of State (EOS) model """
        mantleEOS = None  # Equation of state data to use for silicates
        mantleEOSName = None  # Same as above but containing keywords like clathrates in filenames
        mantleEOSDry = None  # Name of mantle EOS to use assuming non-hydrated silicates
        rhoSilWithCore_kgm3 = 3300  # Assumed density of silicates when a core is present in kg/m^3
        # Derived quantities
        RsilMean_m = None  # Mantle radius for mean compatible moment of inertia (MoI)
        RsilRange_m = None  # Mantle radius range for compatible MoI
        RsilTrade_m = None  # Array of mantle radii for compatible MoIs
        rhoSil_kgm3 = None  # Mantle density consistent with bulk density and MoI when no core is present
        rhoSilTrade_kgm3 = None  # Array of mantle densities for compatible MoIs for core vs. mantle tradeoff plot
        mFluids = None  # WIP for tracking loss of fluids along the geotherm -- needs a better name.
        # The below not necessary to be implemented until later (says Steve), these 5 are based on DPS presentation in 2017 – 5 diff models of permeability
        #turn off this plot feature until later- create flag, Use POROSITY flag to turn off these plots
        #perm1 = None
        #perm2 = None
        #perm3 = None
        #perm4 = None
        #perm5 = None

    class Core:
        """ Core settings """
        rhoFe_kgm3 = 8000  # Assumed density of pure iron in kg/m^3
        rhoFeS_kgm3 = 5150  # Assumed density of iron sulfide in kg/m^3
        rhoPoFeFCC = None  # Density of pyrrhottite plus face-centered cubic iron
        sigmaCore_Sm = 1e-16  # Fixed electrical conductivity to apply to core (typically low, to ignore core impacts on induction)
        coreEOS = 'sulfur_core_partition_SE15_1pctSulfur.tab'  # Default core EOS to use
        # Derived quantities
        RFeMean_m = None  # Core radius for mean compatible moment of inertia (MOI)
        RFeRange_m = None  # Core radius range for compatible MoI
        RFeTrade_m = None  # Array of core radii for compatible MoIs
        #Re Steve- put all mass fraction stuff into a separate file until implemented later- remove from dataStructs.py
        #To implement: possible Meteoritics file/class?
        xFeS_meteoritic = None  # CM2 mean from Jarosewich 1990
        xFeS = None  # mass fraction of sulfur in the core
        xFe_core = None  # this is the total Fe in Fe and FeS
        xH2O = None  # total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model

    class Seismic:
        """ Seismic properties of solids,
            based on Cammarano et al., 2006 (DOI: 10.1029/2006JE002710)
            in a parameterization for the shear anelasticity quality factor QS (Eqs. 4-6).
            g : Homologous temperature scaling constant- dimensionless constant described as a temperature scaling relative to the melting temperature
            B : Shear anelasticity/quality factor normalization constant, dimensionless - helps quantify the effects of anelastic attenuation on the seismic wavelet caused by fluid movement and grain boundary friction
            gamma : Exponent on the seismic frequency omega, dimensionless - helps describe the frequency dependence of attenuation
            """
        lowQDiv = None  # Factor by which to divide the seismic attenuation Q to test out a low-Q value, dimensionless
        # Ice I
        BIceI = None
        gammaIceI = None
        gIceI = None
        # Ice II
        BIceII = None
        gammaIceII = None
        gIceII = None
        # Ice III
        BIceIII = None
        gammaIceIII = None
        gIceIII = None
        # Ice V
        BIceV = None
        gammaIceV = None
        gIceV = None
        # Ice VI
        BIceVI = None
        gammaIceVI = None
        gIceVI = None
        # Silicates
        BSil = None
        gammaSil = None
        gSil = None
        # Derived quantities
        VP_kms = None  # Longitudinal (p-wave) sound velocity for each layer in km/s
        VS_kms = None  # Shear (s-wave) sound velocity for each layer in km/s
        QS = None  # Anelastic shear quality factor Q_S of each layer, divided by omega^gamma to remove frequency dependence. Essentially the fraction of seismic energy lost per cycle.
        KS_GPa = None  # Bulk modulus of each layer in GPa
        GS_GPs = None  # Shear modulus of each layer in GPa

    class Magnetic:
        """ Magnetic induction """
        peaks_Hz = None  # Frequencies in Hz of peaks in Fourier spectrum of magnetic excitations
        fOrb_radps = None  # Angular frequency of orbital motion of a moon around its parent planet in radians per second
        ionosBounds_m = None  # Upper altitude cutoff for ionosphere layers in m. Omit the surface (don't include 0 in the list).
        sigmaIonosPedersen_Sm = None  # Pedersen conductivity for ionospheric layers in S/m. Length must match ionosBounds_m.


class ParamsStruct:
    """ General parameter options """
    PLOT_SIGS = False  # Make a plot of conductivity as a function of radius
    wlims = None  # Minimum and maximum to use for frequency spectrum plots (magnetic induction)
    LEGEND = False  # Whether to force legends to appear
    LegendPosition = None  # Where to place legends when forced
    yLim = None  # y axis limits of hydrosphere density in "Conductivity with interior properties" plot
    LineStyle = None  # Default line style to use on plots
    wRefLine_temporary = None  # Style of lines showing reference melting curves of hydrosphere density plot-should be done in config.py instead, delete this once implemented there
    wRef = None  # Salinities in ppt of reference melting curves

    class lbls: # Not sure we need this in Params for the Python implementation.
        pass


class ConstantsStruct:
    """ General physical constants """
    G = 6.673e-11  # "Big G" gravitational constant, m^3/kg/s
    bar2GPa = 1.01325e-4  # Multiply by this to convert pressure from bars to GPa
    PbI_MPa = 210  # ~fixed transition pressure between ice Ih and ice III or V
    T0 = 273.15  # The Celsius zero point in K.
    P0 = 101325  # One standard atmosphere in Pa
    DThermalConductIceI_Wm = 632  # Thermal conductivity of ice Ih in W/m from Andersson et al. (2005)
    # Core modeling--Modifying thermal profile of mantle
    rhoMantleMean = 3000  # Density of silicate layer in mantle, roughly, in kg/m3
    alphaMantleMean = 0.2e-4  # Thermal expansivity of silicates, roughly, in 1/K
    CpMantleMean = 2e6  # Heat capacity of silicates, roughly, in J/(kgK)
    KappaMantleMean = 1e-6  # ???
    nu_mantle = 1e21  # Mantle viscosity in Pa*s, a common number for Earth's mantle
    DeltaT = 800  # Temperature differential in K between core and mantle (???)

Constants = ConstantsStruct()
