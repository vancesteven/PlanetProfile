"""
dataStructs: Create empty classes and subclasses for holding body-specific data
Values are typically set to None as defaults; body-specific values should be set in the PPBody.py files.

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
    gsurf_ms2 = None  # Surface gravity in m/s^2
    Tsurf_K = None  # Surface temperature in kelvin
    Psurf_MPa = None  # Surface pressure in MPa
    Cmeasured = None  # Axial moment of inertia C/MR^2
    Cuncertainty = None  # Uncertainty in C/MR^2 (used to constrain models via consistency within the uncertainty)
    Tb_K = None  # Temperature at the ice-ocean interface
    POROUS_ICE = None  # Whether to model porosity in ice
    phi_surface = None  # Ice porosity at the surface
    CLATHRATE = None # Whether to model clathrates
    NO_H2O = None # Whether to model waterless worlds (like Io)

    """ Layer step settings """
    Pseafloor_MPa = None  # Maximum pressure at the top of the silicate layer
    nsteps_iceI = None # Fixed number of steps in outermost ice shell
    nsteps_ocean = None # Fixed number of steps in hydrosphere below ice I outer shell
    nsteps_ref_rho = None # ?????
    nsteps_mantle = None # Maximum number of steps in silicate layer (actual number will be fewer)
    nsteps_core = None # Fixed number of steps in core layer, if present

    class Ocean:
        """ Hydrosphere assumptions """
        comp = None # Type of dominant dissolved salt in ocean
        w_ocean_ppt = None # Concentration of above salt in parts per thousand (ppt)
        electrical = None # Type of electrical conductivity model to use

        def fnTfreeze_K(self, PPg, wwg, TT):
            # Somehow make an interpolator a la:
            # load L_Ice_MgSO4.mat
            # fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');
            pass

    class Silicate: # Quantities from the old "Seismic" struct should go in here.
        """ Silicate layer settings """
        # Heat flow
        kr_mantle = None # Thermal conductivity of silicates in W/m/K
        Qmantle_Wm2 = None # Heat flow from mantle into hydrosphere
        QHmantle = None
        EQUIL_Q = None
        SMOOTH_ROCK = None # ?????
        SMOOTH_VROCK = None # ?????
        mfluids = None # WIP for tracking loss of fluids along the geotherm
        # Porosity
        POROUS_ROCK = None # Whether to model the rock as porous
        PEFF = None

        """ Mantle Equation of State (EOS) model """
        mantleEOS = None # Equation of State data to use for silicates
        mantleEOSname = None # Same as above but containing keywords like clathrates in filenames
        mantleEOS_dry = None # Name of mantle EOS to use assuming non-hydrated silicates
        rho_sil_withcore_kgm3 = None

    class Core:
        """ Core assumptions """
        FeCore = None # Whether to model an iron core for this body
        rhoFe = None # Assumed density of pure iron in kg/m^3
        rhoFeS = None # Assumed density of iron sulfide in kg/m^3
        rhoPoFeFCC = None # ±40. Density of pyrrhottite plus face-centered cubic iron
        QScore = None # (??)
        coreEOS = 'sulfur_core_partition_SE15_1pctSulfur.tab' # Default core EOS to use
        xFeS_meteoritic = None # CM2 mean from Jarosewich 1990
        xFeS = None # mass fraction of sulfur in the core
        xFe_core = None # this is the total Fe in Fe and FeS
        XH2O = None # total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
        coreSig = None # Fixed electrical conductivity to apply to core (typically low, to ignore core impacts on induction)

    class Seismic:
        """ Seismic properties of solids """
        low_ice_Q = None # Divide Ice Q value by this number (??)
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
        Magnetic.peaks_Hz = None # Frequencies in Hz of peaks in Fourier spectrum of magnetic excitations
        Magnetic.f_orb = None # Radians per second
        Magnetic.ionos_bounds = None # Upper altitude cutoff for ionosphere conduction
        Magnetic.ionosPedersen_sig = None # Pedersen conductivity for ionospheric layer
        Magnetic.ionos_only = None # Set to ionosphere bottom altitude in the case of no ocean induction, otherwise set to None
        Magnetic.ADD_TRANSITION_BOUNDS = None # Whether to insert another layer entry for changing conductivity at the planetary surface, in the case of a nonconducting atmosphere at the surface.


class ParamsStruct:
    """ General parameter options """
    PLOT_SIGS = None # Make a plot of conductivity as a function of radius
    wlims = None # Minimum and maximum to use for frequency spectrum plots (magnetic induction)
    FOURSUBPLOTS = None # Whether or not to create "Conductivity with interior properties" plot
    LEGEND = None # Whether to force legends to appear
    LegendPosition = None # Where to place legends when forced
    ylim = None # y axis limits of hydrosphere density in "Conductivity with interior properties" plot
    LineStyle = None # Default line style to use on plots
    wrefLine = None # Style of lines showing reference melting curves of hydrosphere density plot
    wref = None # Salinities in ppt of reference melting curves

    tauP = None
    clathrateSetDepth = None # Fixed depth to which clathrates are present until

    BOTTOM_ICEIII = None # Whether to allow Ice III between ocean and ice I layer, when ocean temp is set very low
    BOTTOM_ICEV = None # Same as above but also including ice V
    NO_ICEI_CONVECTION = None # Whether to suppress convection in the ice I layer
    FORCE_ICEI_CONVECTION = None # Whether to force convection in the ice I layer

    class lbls:
        pass


class Constants:
    """ General physical constants """
    G = 6.673e-11 # "Big G" gravitational constant, m^3/kg/s
    bar2GPa = 1e-4 # Multiply by this to convert pressure from bars to GPa
    PbI_MPa = 210 # ~fixed transition pressure between ice Ih and ice III or V