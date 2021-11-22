"""
dataStructs: Create empty classes and subclasses for holding body-specific data
Values are typically set to None as defaults; body-specific values should be set in the PPBody.py files.
Optional SWITCHES are in all caps. These typically have a default value of False.

Example usage:
Planet = PlanetStruct("nameOfBody")
Planet.R_m = 1560e3
Planet.Ocean.comp = "MgSO4"
Planet.Silicate.mantleEOS = "CV3hy1wt_678_1.tab"
Planet.Core.Fe_core = False
"""
class PlanetStruct:
    def __init__(self, name):
        self.name = name

    """ Derived quantities (assigned during PlanetProfile runs) """
    # Layer arrays
    phase = None #phase of the layer input as an integer- ocean, Ice I, etc.,
    r_m = None #distance from center of body to the layer step in meters
    z_m = None #distance from surface of body to layer step in meters
    T_K = None #temperature of each layer in Kelvin
    P_MPa = None #pressure at each step in MPa
    rho_kgm3 = None #density of each layer in kg/m^3
    g_ms2 = None #gravitational acceleration of each layer, m/s^2
    Cp = None  # Heat capacity at constant pressure for each layer's material in W/kg/K
    vfluid = None  # Sound speed in fluid layers in km/s
    Mabove_kg = None #total mass above each layer in kg
    Mbelow_kg = None #total mass below each layer in kg
    C2mean = None  # Mean value of axial moment of inertia consistent with profile core/mantle trades
    rSigChange_m  = None  # Radii of outer boundary of each conducting layer in m
    sigma_Sm = None  # Electrical conductivity (sigma) in S/m of each conducting layer
    # Individual calculated quantities
    zb_km = None  # Thickness of outer ice shell/depth of ice-ocean interface in km in accordance with 2014 paper
    zClath_m = None  # Thickness of clathrate layer at surface in m
    Pb_MPa = None  # Pressure at ice-ocean interface in MPa
    PbI_MPa = None  # Pressure at bottom of ice I layer in MPa, mainly used when BOTTOM_ICEIII is True
    deltaP = None  # Increment of pressure between each layer
    alpha_perK = None # thermal expansivity at each pressure step in hydrosphere [1/K]
    Qmantle_Wm2 = None  # Heat flow from mantle into hydrosphere
    # Q_Wm2 = None  # ??? WAIT UNTIL IMPLEMENT heat flux at ice shell

    """ Run settings """
    class Bulk:
        """ Bulk planetary settings """
        Tb_K = None  # Temperature at the ice-ocean interface
        rho_kgm3 = None  # bulk density in kg/m3
        R_m = None  # Mean body outer radius in m
        M_kg = None  # Total body mass in kg
        Tsurf_K = None  # Surface temperature in kelvin
        Psurf_MPa = None  # Surface pressure in MPa
        Cmeasured = None  # Axial moment of inertia C/MR^, dimensionless
        Cuncertainty = None  # Uncertainty in C/MR^2 (used to constrain models via consistency within the uncertainty), dimensionless
        phiSurface = None  # Scaling value for the ice porosity at the surface (void fraction[]): falls within a range of 0 and 1, larger than 0.2 is rare, from Han et.al. 2014
        clathrateSetDepth = None  # Fixed depth for limiting clathrate layer

    """ Runtime flags """
    class Do:
        Fe_CORE = None  # Whether to model an iron core for this body
        POROUS_ICE = False  # Whether to model porosity in ice
        CLATHRATE = False  # Whether to model clathrates
        NO_H2O = False  # Whether to model waterless worlds (like Io)
        BOTTOM_ICEIII = False  # Whether to allow Ice III between ocean and ice I layer, when ocean temp is set very low- default is that this is off, can turn on as an error condition
        BOTTOM_ICEV = False  # Same as above but also including ice V. Takes precedence (forces both ice III and V to be present).
        NO_ICEI_CONVECTION = False  # Whether to suppress convection in the ice I layer - if True, checks Rayleigh number to see if convection conditions are met
        FORCE_ICEI_CONVECTION = False  # Whether to force convection in the ice I layer- if True, doesn’t check Rayleigh number, just assumes convection conditions are met
        ALLOW_NEG_ALPHA = False  # Whether to permit modeling of a Melosh et. al. layer with negative thermal expansivity
        MANTLE_HYDRO_THERMAL_EQ = False  # option to have thermal equilibrium between mantle and hydrosphere where the hydrosphere is not gaining external heat via tidal forcing or radiation
        POROUS_ROCK = False  # Whether to model the rock as porous
        P_EFFECTIVE = None  # effective pressure due to presence of water in pores: lithostatic minus hydrostatic pressure
        IONOS_ONLY = False  # Whether to ignore conducting layers within the body and model magnetic induction happening only in the ionosphere

    """ Layer step settings """
    class Steps:
        nIceI = None  # Fixed number of steps in outermost ice shell
        nClath = None  # Fixed number of steps in clathrates
        nOceanMax = None  # Maximum number of steps in hydrosphere below ice I outer shell (actual number will be fewer)
        firstSilLayerIndex = None  # First index of silicates, i.e. derived actual number of layers in hydrosphere (based on Moment of Inetia), plus 1
        nStepsRefRho = None  # number of different reference densities will be calculated in the run
        nSil = None  # Fixed number of steps in silicate layer
        nCore = None  # Fixed number of steps in core layer, if present
        nIceIIILitho = None  # Fixed number of layers to use for ice III when either BOTTOM_ICEIII or BOTTOM_ICEV is True- empirical number- should be user-specified
        nIceVLitho = None  # Fixed number of layers to use for ice V when BOTTOM_ICEV is True- empirical number- should be user-specified

    class Ocean:
        """ Hydrosphere assumptions """
        comp = None  # Type of dominant dissolved salt in ocean
        wOcean_ppt = None  # salinity: Concentration of above salt in parts per thousand (ppt) for Seawater or weight percent (wtPct) for MgSO4
        electrical = None  # Type of electrical conductivity model to use, passed as a string
        PMaxHydrosphere_MPa = mantleEOSname # Maximum pressure (upper bound of pressure range) of the hydrosphere
        deltaPOcean = (PMaxHydrosphere_MPa – Pb_MPa) / nOceanMax #difference in pressure in each step through ocean in MPa

        def fnTfreeze_K(self, PPg, wwg, TT):
            # Somehow make an interpolator a la:
            # load L_Ice_MgSO4.mat
            # fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');
            pass

    class Sil:  # Quantities from the old "Seismic" struct should go in here.
        """ Silicate layer settings """
        # Heat flow
        kSil_WmK = None  # Thermal conductivity (k) of silicates in W/(mK)
        mFluids = None  # WIP for tracking loss of fluids along the geotherm
        # TauP
        tauP = None  # ????
        RsilMean_m = None  # Mantle radius for mean compatible Moment of Inertia (MoI)
        RsilRange_m = None  # Mantle radius range for compatible MoI
        RsilTrade_m = None  # Array of mantle radii for compatible MoIs
        rho_sil_trade_kgm3 = None  # Array of mantle densities for compatible MoIs for core vs. mantle tradeoff plot
        phiRockMax = None  # Porosity (void fraction) of the rocks at the “seafloor”, where the ocean first comes into contact with rock
        # Not necessary to be implemented until later re Steve, these 5 are based on DPS presentation in 2017 – 5 diff models of permeability
        #turn off this plot feature until later- create flag, Use POROSITY flag to turn off these plots
        #perm1 = None  # No idea what what these numbers are or why there are 5 of them -MJS
        #perm2 = None
        #perm3 = None
        #perm4 = None
        #perm5 = None

        """ Mantle Equation of State (EOS) model """
        mantleEOS = None  # Equation of State data to use for silicates
        mantleEOSName = None  # Same as above but containing keywords like clathrates in filenames
        mantleEOS_dry = None  # Name of mantle EOS to use assuming non-hydrated silicates
        rhoSilWithCore_kgm3 = None # Assumed density of silicates when a core is present in kg/m^3

    class Core:
        """ Core assumptions """
        rhoFe_kgm3 = None  # Assumed density of pure iron in kg/m^3
        rhoFeS_kgm3 = None  # Assumed density of iron sulfide in kg/m^3
        rhoPoFeFCC = None  # ±40. Density of pyrrhottite plus face-centered cubic iron
        coreEOS = 'sulfur_core_partition_SE15_1pctSulfur.tab'  # Default core EOS to use

        #Re Steve- put all mass fraction stuff into a separate file until implemented later- remove from dataStructs.py
        #To implement: possible Meteoritics file/class?
        xFeS_meteoritic = None  # CM2 mean from Jarosewich 1990
        xFeS = None  # mass fraction of sulfur in the core
        xFe_core = None  # this is the total Fe in Fe and FeS
        xH2O = None  # total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model

        sigmaCore_Sm = None  # Fixed electrical conductivity to apply to core (typically low, to ignore core impacts on induction)
        # Derived quantities
        RFeMean_m = None  # Core radius for mean compatible Moment of Inertia (MOI)
        RFeRange_m = None  # Core radius range for compatible MoI
        RFeTrade_m = None  # Array of core radii for compatible MoIs

    class Seismic:
        """ Seismic properties of solids """
        low_ice_Q = None  # Divide Ice Q value by this number (??)
        # ice I
        B_aniso_iceI = None
        gamma_aniso_iceI = None
        g_aniso_iceI = None
        # ice II
        B_aniso_iceII = None
        gamma_aniso_iceII = None
        g_aniso_iceII = None
        # ice III
        B_aniso_iceIII = None
        gamma_aniso_iceIII = None
        g_aniso_iceIII = None
        # ice V
        B_aniso_iceV = None
        gamma_aniso_iceV = None
        g_aniso_iceV = None
        # ice VI
        B_aniso_iceVI = None
        gamma_aniso_iceVI = None
        g_aniso_iceVI = None
        # mantle
        B_aniso_mantle = None
        gamma_aniso_mantle = None
        g_aniso_mantle = None

    class Magnetic:
        """ Magnetic induction """
        peaks_Hz = None  # Frequencies in Hz of peaks in Fourier spectrum of magnetic excitations
        fOrb = None  # orbital period of a moon around its parent planet in radians per second
        ionosBounds_m = None  # Upper altitude cutoff for ionosphere conduction in meters- zero is the body surface and the bounds increase as you move into space, can take arrays if ionosphere has multiple layers
        sigmaIonosPedersen_Sm = None  # Pedersen conductivity for ionospheric layer in S/m
        IONOS_ONLY = False  # Set to ionosphere bottom altitude in the case of no ocean induction, otherwise defaults to False
        ADD_TRANSITION_BOUNDS = False  # Whether to insert another layer entry for changing conductivity at the planetary surface, in the case of a nonconducting atmosphere at the surface.


class ParamsStruct:
    """ General parameter options """
    PLOT_SIGS = False  # Make a plot of conductivity as a function of radius
    wlims = None  # Minimum and maximum to use for frequency spectrum plots (magnetic induction)
    FOURSUBPLOTS = False  # Whether or not to create "Conductivity with interior properties" plot
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
    T0 = 273.15  # The Celcius zero point; 273.15 K.
    P0 = 101325  # one standard atmosphere in Pa, 1atm = 101325 Pa
    DThermalConductIceI_Wm = 632  # Thermal conductivity of ice Ih in W/m from Andersson et al. (2005)
    # Core modeling- Modifying thermal profile of mantle: Pre-allocate with empty values, then populate once needed calcs are done
    rhoMantleMean = 3000  # Density of silicate layer in mantle, roughly, in kg/m3
    alphaMantleMean = 0.2e-4  # Thermal expansivity of silicates, roughly, in 1/K
    CpMantleMean = 2e6  # Heat capacity of silicates, roughly, in J/(kgK)
    KappaMantleMean = 1e-6  # ???
    nu_mantle = 1e21  # mantle viscosity in Pa*S, a common number for Earth's mantle
    DeltaT = 800  # Temperature differential in K between core and mantle (???)

Constants = ConstantsStruct()
