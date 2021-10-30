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

    """ Bulk planetary settings """
    rho_kgm3 = None  # ±46 (Schubert et al. 2004)
    R_m = None  # ±8.0 km. Radius in m
    M_kg = None  # Mass in kg
    Tsurf_K = None  # Surface temperature in kelvin
    Psurf_MPa = None  # Surface pressure in MPa
    CMeasured = None  # Axial moment of inertia C/MR^2
    CUncertainty = None  # Uncertainty in C/MR^2 (used to constrain models via consistency within the uncertainty)
    Tb_K = None  # Temperature at the ice-ocean interface
    POROUS_ICE = False  # Whether to model porosity in ice
    phiSurface = None  # Ice porosity at the surface
    CLATHRATE = False  # Whether to model clathrates
    clathrateSetDepth = None  # Fixed depth for limiting clathrate layer
    NO_H2O = None  # Whether to model waterless worlds (like Io)
    BOTTOM_ICEIII = False  # Whether to allow Ice III between ocean and ice I layer, when ocean temp is set very low
    BOTTOM_ICEV = False  # Same as above but also including ice V. Takes precedence (forces both ice III and V to be present).
    NO_ICEI_CONVECTION = False  # Whether to suppress convection in the ice I layer
    FORCE_ICEI_CONVECTION = False  # Whether to force convection in the ice I layer
    ALLOW_NEGALPHA = False  # Whether to permit modeling of a Melosh layer with negative thermal expansivity

    """ Derived quantities (assigned during PlanetProfile runs) """
    Pb_MPa = None  # Pressure at ice-ocean interface in MPa
    PbI_MPa = None  # Pressure at bottom of ice I layer in MPa, mainly used when BOTTOM_ICEIII is True
    alpha_o = None  # Thermal expansivity of ice at ice-ocean interface
    boundaries = None  # Radii of outer boundary of each conducting layer in m
    sig_S_m = None  # Electrical conductivity (sigma) in S/m of each conducting layer
    Zb_km = None  # Thickness of outer ice shell/depth of ice-ocean interface in km
    zClath_m = None  # Thickness of clathrate layer at surface in m

    """ Layer step settings """
    Pseafloor_MPa = None  # Maximum pressure at the top of the silicate layer
    nStepsIceI = None  # Fixed number of steps in outermost ice shell
    nStepsClath = None  # Fixed number of steps in clathrates
    nStepsOcean = None  # Fixed number of steps in hydrosphere below ice I outer shell
    nStepsRefRho = None  # ?????
    nStepsMantle = None  # Maximum number of steps in silicate layer (actual number will be fewer)
    nStepsCore = None  # Fixed number of steps in core layer, if present
    nIceIIILitho = 5  # Fixed number of layers to use for ice III when either BOTTOM_ICEIII or BOTTOM_ICEV is True
    nIceVLitho = 5  # Fixed number of layers to use for ive V when BOTTOM_ICEV is True

    class Ocean:
        """ Hydrosphere assumptions """
        comp = None  # Type of dominant dissolved salt in ocean
        wtOcean_ppt = None  # Concentration of above salt in parts per thousand (ppt)
        electrical = None  # Type of electrical conductivity model to use

        def fnTfreeze_K(self, PPg, wwg, TT):
            # Somehow make an interpolator a la:
            # load L_Ice_MgSO4.mat
            # fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');
            pass

    class Silicate:  # Quantities from the old "Seismic" struct should go in here.
        """ Silicate layer settings """
        # Heat flow
        krMantle_WmK = None  # Thermal conductivity of silicates in W/m/K
        QHMantle = None
        EQUIL_Q = None
        SMOOTH_ROCK = None  # ?????
        SMOOTH_VROCK = None  # ?????
        mFluids = None  # WIP for tracking loss of fluids along the geotherm
        # Porosity
        POROUS_ROCK = False  # Whether to model the rock as porous
        PEFF = None
        # Derived quantities
        QMantle_Wm2 = None  # Heat flow from mantle into hydrosphere

        """ Mantle Equation of State (EOS) model """
        mantleEOS = None  # Equation of State data to use for silicates
        mantleEOSName = None  # Same as above but containing keywords like clathrates in filenames
        mantleEOS_dry = None  # Name of mantle EOS to use assuming non-hydrated silicates
        rhoSilWithCore_kgm3 = None

    class Core:
        """ Core assumptions """
        FeCORE = None  # Whether to model an iron core for this body
        rhoFe = None  # Assumed density of pure iron in kg/m^3
        rhoFeS = None  # Assumed density of iron sulfide in kg/m^3
        rhoPoFeFCC = None  # ±40. Density of pyrrhottite plus face-centered cubic iron
        QSCore = None  # (??)
        coreEOS = 'sulfur_core_partition_SE15_1pctSulfur.tab'  # Default core EOS to use
        xFeS_meteoritic = None  # CM2 mean from Jarosewich 1990
        xFeS = None  # mass fraction of sulfur in the core
        xFe_core = None  # this is the total Fe in Fe and FeS
        XH2O = None  # total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
        coreSig = None  # Fixed electrical conductivity to apply to core (typically low, to ignore core impacts on induction)

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
        fOrb = None  # Radians per second
        ionosBounds = None  # Upper altitude cutoff for ionosphere conduction
        ionosPedersenSig = None  # Pedersen conductivity for ionospheric layer
        ionosOnly = None  # Set to ionosphere bottom altitude in the case of no ocean induction, otherwise set to None
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
    wRefLine = None  # Style of lines showing reference melting curves of hydrosphere density plot
    wRef = None  # Salinities in ppt of reference melting curves

    tauP = None

    class lbls: # Not sure we need this in Params for the Python implementation.
        pass


class ConstantsStruct:
    """ General physical constants """
    G = 6.673e-11  # "Big G" gravitational constant, m^3/kg/s
    bar2GPa = 1.01325e-4  # Multiply by this to convert pressure from bars to GPa
    PbI_MPa = 210  # ~fixed transition pressure between ice Ih and ice III or V
    T0 = 273.15  # The Celcius zero point; 273.15 K.
    P0 = 101325  # one standard atmosphere in Pa, 1atm = 101325 Pa
    DThermalConductIceI_Wm = 632 # Thermal conductivity of ice Ih in W/m from Andersson et al. (2005)
    # Core modeling
    rhom_rough = 3000  # Density of silicate layer, roughly
    alpha_rough = 0.2e-4  # Thermal expansivity of silicates, roughly
    Cp_rough = 2e6  # Heat capacity of silicates, roughly
    Kappa_rough = 1e-6  # ???
    nu_mantle = 1e21  # mantle viscosity in Pa*S, a common number for Earth's mantle
    DeltaT = 800  # Temperature differential in K between core and mantle (???)
