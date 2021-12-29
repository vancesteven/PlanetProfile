import numpy as np
import scipy.interpolate as spi
import scipy.optimize as spo
from seafreeze import seafreeze as SeaFreeze
from seafreeze import whichphase as WhichPhase
from Utilities.dataStructs import Constants

class OceanEOSStruct:
    def __init__(self, compstr, wOcean_ppt, P_MPa, T_K):
        self.comp = compstr
        self.w_ppt = wOcean_ppt

        # Get tabular data from the appropriate source for this composition
        if wOcean_ppt == 0:
            self.type = 'SeaFreeze'
            PTgrid = np.array([P_MPa, T_K], dtype=object)
            seaOut = SeaFreeze(PTgrid, 'water1')
            self.rho_kgm3 = seaOut.rho
            self.Cp_JkgK = seaOut.Cp
            self.alpha_pK = seaOut.alpha
            self.kTherm_WmK = np.zeros_like(self.alpha_pK) + Constants.kThermWater_WmK  # Placeholder until we implement a self-consistent calculation
            self.phase = WhichPhase(PTgrid)
        elif compstr == 'Seawater':
            self.type = 'GSW'
            raise ValueError('Unable to load ocean EOS. Seawater is not implemented yet.')
        elif compstr == 'NH3':
            self.type = 'PlanetProfile'
            raise ValueError('Unable to load ocean EOS. NH3 is not implemented yet.')
        elif compstr == 'MgSO4':
            self.type = 'LBF'
            raise ValueError('Unable to load ocean EOS. MgSO4 is not implemented yet.')
        elif compstr == 'NaCl':
            self.type = 'PlanetProfile'
            raise ValueError('Unable to load ocean EOS. NaCl is not implemented yet.')
        else:
            raise ValueError('Unable to load ocean EOS. compstr="'+compstr+'" but options are Seawater, NH3, MgSO4, and NaCl.')

        self.fn_rho_kgm3 = spi.RectBivariateSpline(P_MPa, T_K, self.rho_kgm3)
        self.fn_Cp_JkgK = spi.RectBivariateSpline(P_MPa, T_K, self.Cp_JkgK)
        self.fn_alpha_pK = spi.RectBivariateSpline(P_MPa, T_K, self.alpha_pK)
        self.fn_kTherm_WmK = spi.RectBivariateSpline(P_MPa, T_K, self.kTherm_WmK)

        # Create phase finder -- note that the results from this function must be cast to int after retrieval
        Plin_MPa = np.array([P for P in P_MPa for _ in T_K])
        Tlin_K = np.array([T for _ in P_MPa for T in T_K])
        PTpairs = list(zip(Plin_MPa, Tlin_K))
        phase1D = np.reshape(self.phase, (-1))
        # Create phase finder -- note that the results from this function must be cast to int after retrieval
        self.fn_phase = spi.NearestNDInterpolator(PTpairs, phase1D)


class IceEOSStruct:
    def __init__(self, P_MPa, T_K, phaseStr):
        # Make sure arrays are long enough to interpolate
        nPs = np.size(P_MPa)
        nTs = np.size(T_K)
        if(nPs <= 3):
            P_MPa = np.linspace(P_MPa[0], P_MPa[-1], nPs*3)
        if(nTs <= 3):
            T_K = np.linspace(T_K[0], T_K[-1], nTs*3)

        # If input arrays are equal length, repeat final T value due to a quirk of numpy arrays
        # combined with SeaFreeze's particular implementation that requires gridded P,T values
        # to have different array lengths
        if(nPs == nTs):
            T_K = np.append(T_K, T_K[-1]*1.00001)

        if phaseStr == 'Clath':
            # Special functions for clathrate properties
            self.rho_kgm3, self.Cp_JkgK, self.alpha_pK, self.kTherm_WmK \
                = ClathProps(P_MPa, T_K)
            self.phase = ClathStableSloan1998(P_MPa, T_K)

            Plin_MPa = np.array([P for P in P_MPa for _ in T_K])
            Tlin_K = np.array([T for _ in P_MPa for T in T_K])
            PTpairs = list(zip(Plin_MPa, Tlin_K))
            phase1D = np.reshape(self.phase, (-1))
            # Create phase finder -- note that the results from this function must be cast to int after retrieval
            # Returns either Constants.phaseClath (stable) or 0 (not stable), making it compatible with GetTfreeze
            self.fn_phase = spi.NearestNDInterpolator(PTpairs, phase1D)
        else:
            # Get tabular data from SeaFreeze for all other ice phases
            PTgrid = np.array([P_MPa, T_K], dtype=object)
            iceOut = SeaFreeze(PTgrid, phaseStr)
            self.rho_kgm3 = iceOut.rho
            self.Cp_JkgK = iceOut.Cp
            self.alpha_pK = iceOut.alpha
            self.kTherm_WmK = np.array([kThermIsobaricAnderssonIbari2005(T_K, PhaseInv(phaseStr)) for _ in P_MPa])

        # Interpolate functions for this ice phase that can be queried for properties
        self.fn_rho_kgm3 = spi.RectBivariateSpline(P_MPa, T_K, self.rho_kgm3)
        self.fn_Cp_JkgK = spi.RectBivariateSpline(P_MPa, T_K, self.Cp_JkgK)
        self.fn_alpha_pK = spi.RectBivariateSpline(P_MPa, T_K, self.alpha_pK)
        self.fn_kTherm_WmK = spi.RectBivariateSpline(P_MPa, T_K, self.kTherm_WmK)
        # Assign phase ID and string for convenience in functions where iceEOS is passed
        self.phaseStr = phaseStr
        self.phaseID = PhaseInv(phaseStr)


def GetPfreeze(oceanEOS, Tb_K, PLower_MPa=20, PUpper_MPa=300, PRes_MPa=0.1,
               Pguess=None, guessRange=5):
    """ Returns the pressure at which surface ice changes phase based on temperature, salinity, and composition

        Args:
            oceanEOS (OceanEOSStruct): Interpolator functions for evaluating the ocean EOS
            Tb_K (float): Temperature of the phase transition in K
        Returns:
            Pfreeze_MPa (float): Pressure at the phase change interface consistent with Tb_K
    """
    Psearch = np.arange(PLower_MPa, PUpper_MPa, PRes_MPa)
    # Suggest a smaller range around a guessed value
    if Pguess is not None:
        DO_GUESS = True
        PguessRange = np.arange(Pguess-guessRange/2, Pguess+guessRange/2, PRes_MPa)
    else:
        DO_GUESS = False
        GUESS_FAILED = True

    # Get phase of each P for the Tb_K value using SeaFreeze
    if DO_GUESS:
        searchPhases = oceanEOS.fn_phase(PguessRange, Tb_K).astype(np.int_)
        # Check if we failed to encounter a phase transition
        if np.all(searchPhases==searchPhases[0]): GUESS_FAILED = True
    if not DO_GUESS or GUESS_FAILED:
        searchPhases = oceanEOS.fn_phase(Psearch, Tb_K)
    # Find the first index for a phase that's not ice Ih or clathrates
    # (note that clathrates, with phase Constants.phaseClath, are not yet implemented in SeaFreeze)
    try:
        indMelt = next((i[0] for i, val in np.ndenumerate(searchPhases) if val!=1 and val!=Constants.phaseClath))
    except StopIteration:
        raise ValueError('No transition pressure was found below '+str(PUpper_MPa)+' MPa '+
                         'for ice Ih/clathrates. Increase PUpper_MPa until one is found.')
    # Get the pressure of the first non-Ih layer
    Pfreeze_MPa = Psearch[indMelt]

    return Pfreeze_MPa


def GetPfreezeHP(oceanEOS, Pmin_MPa, TbHP_K, phase, PLower_MPa=180, PUpper_MPa=900, PRes_MPa=0.5):
    """ Returns the pressure at which a high-pressure ice changes phase based on
        temperature, salinity, and composition

        Args:
            oceanEOS (OceanEOSStruct): Interpolator functions for evaluating the ocean EOS
            Pmin_MPa (float): Minimum pressure to accept for a result
            TbHP_K (float): Temperature of the phase transition in K
        Returns:
            PfreezeHP_MPa (float): Pressure at phase change interface
    """
    # Narrow the search range for narrow-stability HP ices
    if phase==3:
        PUpper_MPa = np.min([360, PUpper_MPa])
        PLower_MPa = np.max([209, PLower_MPa, Pmin_MPa])
    elif phase==5:
        PUpper_MPa = np.min([640, PUpper_MPa])
        PLower_MPa = np.max([343, PLower_MPa, Pmin_MPa])

    Psearch = np.arange(PLower_MPa, PUpper_MPa, PRes_MPa)

    # Get phase of each P from the ocean EOS
    searchPhases = oceanEOS.fn_phase(Psearch, TbHP_K).astype(np.int_)
    # Find the first index of desired HP ice
    try:
        indIceHP = next((i[0] for i, val in np.ndenumerate(searchPhases) if val==phase))
    except StopIteration:
        raise ValueError('No ice '+str(PhaseConv(phase))+' was found within the range '+str(PLower_MPa)+
              ' < P < '+str(PUpper_MPa)+' for Tb = '+str(TbHP_K)+'.')
    # Find the first index for a phase that's not the desired one after we encounter the desired phase
    try:
        indNotThisHP = next((i[0] for i, val in np.ndenumerate(searchPhases) if i[0]>indIceHP and val!=phase))
    except StopIteration:
        raise ValueError('No transition pressure was found below '+str(PUpper_MPa)+' MPa '+
                         'for ice '+PhaseConv(phase)+' at T = '+str(round(TbHP_K,3))+' K.'+
                         ' Increase PUpper_MPa until one is found.')
    # Get the pressure of the first lower layer with a non-matching phase
    PfreezeHP_MPa = Psearch[indNotThisHP]

    return PfreezeHP_MPa


def GetPmeltHP(oceanEOS, Pmin_MPa, TbHP_K, phase, PLower_MPa=180, PUpper_MPa=900, PRes_MPa=0.5):
    """ Returns the pressure at which a high-pressure ice changes phase to liquid based on
        temperature, salinity, and composition

        Args:
            oceanEOS (OceanEOSStruct): Interpolator functions for evaluating the ocean EOS
            Pmin_MPa (float): Minimum pressure to accept for a result
            TbHP_K (float): Temperature of the phase transition in K
        Returns:
            PmeltHP_MPa (float): Pressure at phase change interface
    """
    # Narrow the search range for narrow-stability HP ices
    if phase==3:
        PUpper_MPa = np.min([360, PUpper_MPa])
        PLower_MPa = np.max([209, PLower_MPa, Pmin_MPa])
    elif phase==5:
        PUpper_MPa = np.min([640, PUpper_MPa])
        PLower_MPa = np.max([343, PLower_MPa, Pmin_MPa])

    Psearch = np.arange(PLower_MPa, PUpper_MPa, PRes_MPa)

    # Get phase of each P from the ocean EOS
    searchPhases = oceanEOS.fn_phase(Psearch, TbHP_K).astype(np.int_)
    # Find the first index of desired HP ice
    try:
        indIceHP = next((i[0] for i, val in np.ndenumerate(searchPhases) if val==phase))
    except StopIteration:
        raise ValueError('No ice '+str(PhaseConv(phase))+' was found within the range '+str(round(PLower_MPa,3))+
              ' < P < '+str(round(PUpper_MPa,3))+' for Tb = '+str(round(TbHP_K,3))+' with this ocean composition.')
    # Check if the next index below is liquid
    if(searchPhases[indIceHP-1] == 0):
        PmeltHP_MPa = Psearch[indIceHP-1]
    else:
        raise ValueError('Pmelt found a phase transition that was from '+str(PhaseConv(Psearch[indIceHP]))+
                         ' to '+PhaseConv(Psearch[indIceHP-1])+' instead of liquid. Try adjusting Tb_K values.')

    return PmeltHP_MPa


def GetTfreeze(oceanEOS, P_MPa, T_K, TfreezeRange_K=50, TfreezeRes_K=0.05):
    """ Returns the temperature at which a solid layer melts based on temperature, salinity, and composition

        Args:
            oceanEOS (OceanEOSStruct): Interpolator functions for evaluating the ocean EOS
            P_MPa (float): Pressure of the fluid in MPa
            T_K (float): Temperature of the fluid in K
        Returns:
            Tfreeze_K (float): Temperature of nearest higher-temperature phase transition between
                liquid and ice at this pressure
    """
    Tsearch = np.arange(T_K, T_K + TfreezeRange_K, TfreezeRes_K)

    searchPhases = oceanEOS.fn_phase(P_MPa, Tsearch).astype(np.int_)
    # Find the first index for liquid
    try:
        indLiquid = next((i[0] for i, val in np.ndenumerate(searchPhases) if val==0))
    except StopIteration:
        raise ValueError('No melting temperature was found above ' + str(round(T_K,3)) + ' K ' +
                         'for ice ' + PhaseConv(thisPhase) + ' at pressure ' + str(round(P_MPa,3)) +
                         ' MPa. Increase TfreezeRange_K until one is found.')
    # Get the temperature of the first liquid index
    Tfreeze_K = Tsearch[indLiquid]

    return Tfreeze_K


def PhaseConv(phase):
    """ Convert phase integers into strings compatible with SeaFreeze

        Arguments:
            phase (int): ID of phase for each layer
        Returns:
            phaseStr (string): Corresponding string for each phase ID
    """
    if phase == 0:
        phaseStr = 'water1'
    elif phase == 1:
        phaseStr = 'Ih'
    elif phase == 2:
        phaseStr = 'II'
    elif phase == 3:
        phaseStr = 'III'
    elif phase == 5:
        phaseStr = 'V'
    elif phase == 6:
        phaseStr = 'VI'
    elif phase == Constants.phaseClath:
        phaseStr = 'Clath'
    elif phase == Constants.phaseSil:
        phaseStr = 'Sil'
    elif phase == Constants.phaseFe:
        phaseStr = 'Fe'
    else:
        raise ValueError('PhaseConv does not have a definition for phase ID '+str(phase)+'.')

    return phaseStr


def PhaseInv(phaseStr):
    """ Convert phase strings compatible with SeaFreeze into integers

        Arguments:
            phaseStr (string): String for each phase ID
        Returns:
            phase (int): Corresponding ID of phase for each layer
    """
    if phaseStr == 'water1':
        phase = 0
    elif phaseStr == 'Ih':
        phase = 1
    elif phaseStr == 'II':
        phase = 2
    elif phaseStr == 'III':
        phase = 3
    elif phaseStr == 'V':
        phase = 5
    elif phaseStr == 'VI':
        phase = 6
    elif phaseStr == 'Clath':
        phase = Constants.phaseClath
    elif phaseStr == 'Sil':
        phase = Constants.phaseSil
    elif phaseStr == 'Fe':
        phase = Constants.phaseFe
    else:
        raise ValueError('PhaseInv does not have a definition for phase string "'+phaseStr+'".')

    return phase


def GetPhaseIndices(phase):
    """ Get indices for each phase of ice/liquid

        Args:
            phase (int, shape N)
        Returns:
            indsLiquid, indsIceI, ... indsFe (int, shape 0-M): lists of indices corresponding to each phase.
                Variable length.
    """
    # Avoid an annoying problem where np.where returns an empty array for a length-1 list,
    # by making sure the input value(s) are a numpy array:
    phase = np.array(phase)

    indsLiquid = np.where(phase==0)
    indsIceI = np.where(phase==1)
    indsIceII = np.where(phase==2)
    indsIceIII = np.where(phase==3)
    indsIceV = np.where(phase==5)
    indsIceVI = np.where(phase==6)
    indsClath = np.where(phase==Constants.phaseClath)
    indsSil = np.where(phase==Constants.phaseSil)
    indsFe = np.where(phase==Constants.phaseFe)

    return indsLiquid, indsIceI, indsIceII, indsIceIII, indsIceV, indsIceVI, indsClath, indsSil, indsFe


def kThermIsobaricAnderssonIbari2005(T_K, phase):
    """ Calculate thermal conductivity of ice at a fixed pressure according to
        Andersson and Ibari (2005) as a function of temperature.
        Range of validity is as follows:
        Phase:  P (MPa):    T range (K):
        Ih      0.1         40-180*
        II      240         120-240
        III     240         180-250
        V       530         240-270
        VI      1000        135-250
        *Andersson and Ibari give an alternate equation that accounts for the range 180-273 K
        for ice Ih at 0.1 MPa, but as this was not included in the Matlab version, it's
        skipped here too. This implementation does not apply at the relevant T and P values
        for icy moon shells except at specific points, so a more versatile and accurate
        model should be found and used to replace this.

        Args:
            T_K (float, shape N): Temperatures to evaluate in K
            phase (int): Phase ID
        Returns:
            kTherm_WmK (float, shape N): Thermal conductivity of desired phase at specified temperatures
                in W/(mK)
    """
    D = np.array([np.nan, 630, 695, 93.2, np.nan, 38.0, 50.9])
    X = np.array([np.nan, 0.995, 1.097, 0.822, np.nan, 0.612, 0.612])

    kTherm_WmK = D[phase] * T_K**(-X[phase])

    return kTherm_WmK


def kThermIsothermalAnderssonIbari2005(P_MPa, phase):
    """ Calculate thermal conductivity of ice at a fixed temperature according to
        Andersson and Ibari (2005) as a function of pressure.
        Range of validity is as follows:
        Phase:  P range (GPa):  T (K):
        Ih      0-0.5           130
        II      0-0.24          120
        III     0.2-0.35        240
        V       0.35-0.6        246
        VI      0.7-2.0         246
        This implementation does not apply at the relevant T and P values for icy moon
        shells except at specific points, so a more versatile and accurate model should
        be found and used to replace this.

        Args:
            P_MPa (float, shape N): Pressure to evaluate in MPa
            phase (int, shape N): Phase index
        Returns:
            kTherm_WmK (float, shape N): Thermal conductivity of desired phase at specified pressures
                in W/(mK)
    """
    E = np.array([np.nan, 1.60, 1.25, -0.02, np.nan, 0.16, 0.37])
    F = np.array([np.nan, -0.44, 0.2, 0.2, np.nan, 0.2, 0.16])

    # Note the 1e-3 factor because F has units of 1/GPa
    kTherm_WmK = np.exp(E[phase] + F[phase] * P_MPa * 1e-3)

    return kTherm_WmK


def kThermMelinder2007(T_K, Tmelt_K, ko_WmK=2.21, dkdT_WmK2=-0.012):
    """ Calculate thermal conductivity of ice Ih according to Melinder (2007).

        Args:
            T_K (float, shape N): Temperature in K
            Tmelt_K (float, shape N): Melting temperature at the evaluated pressure in K
            ko_WmK = 2.21 (float): Thermal conductivity at the melting temperature in W/(mK)
            dkdT_WmK2 = -0.012 (float): Constant temperature derivative of k in W/(mK^2)
        Returns:
            kTherm_WmK (float, shape N): Thermal conductivity of ice Ih at specified temperature
                in W/(mK)
    """

    kTherm_WmK = ko_WmK + dkdT_WmK2 * (T_K - Tmelt_K)
    return kTherm_WmK


def kThermHobbs1974(T_K):
    """ Calculate thermal conductivity of ice Ih according to Hobbs (1974), as
        reported by Ojakangas and Stevenson (1989).

        Args:
            T_K (float, shape N): Temperature value(s) in K
        Returns:
            kTherm_WmK (float, shape N): Thermal conductivities in W/(mK)
    """
    a0 = 4.68e4  # Units of ergs/(K cm s)
    a1 = 4.88e7  # Units of ergs/(cm s)
    a0_SI = a0 * Constants.erg2J * 1e2
    a1_SI = a1 * Constants.erg2J * 1e2
    kTherm_WmK = a1_SI/T_K + a0_SI

    return kTherm_WmK


def ClathProps(Plin_MPa, Tlin_K):
    """ Evaluate methane clathrate physical properties using Helgerud et al. (2009): https://doi.org/10.1029/2009JB006451
        for density, Ning et al. (2015): https://doi.org/10.1039/C4CP04212C for thermal expansivity and heat capacity,
        and Waite et al. (2005): https://www.researchgate.net/profile/W-Waite/publication/252708287_Thermal_Property_Measurements_in_Tetrahydrofuran_THF_Hydrate_Between_-25_and_4deg_C_and_Their_Application_to_Methane_Hydrate/links/57b900ae08aedfe0ec94abd7/Thermal-Property-Measurements-in-Tetrahydrofuran-THF-Hydrate-Between-25-and-4deg-C-and-Their-Application-to-Methane-Hydrate.pdf
        for thermal conductivity.
        Range of validity:
            rho_kgm3: P from 30.5 to 97.7 MPa, T from -20 to 15 C
            Cp_JkgK: P at 20 MPa, T from 5 to 292 K
            alpha_pK: P at 0.1 MPa, T from 5 to 268 K (note that Ning et al. also give a parameterization for
                alpha_pK at 20 MPa that differs slightly. Since clathrates only appear at the surface of icy moons,
                we use just the 1 bar value for simplicity.
            kTherm_WmK: P from 13.8 to 24.8 MPa, T from -30 to 20 C


        Args:
            Plin_MPa (float, shape M): Pressures to evaluate in MPa
            Tlin_K (float, shape N): Temperatures to evaluate in K
        Returns:
            rho_kgm3 (float, shape MxN): Mass density in kg/m^3 for each P and T. rho_gcm3 = aT_C + bP_MPa + c
            Cp_JkgK (float, shape MxN): Heat capacity in J/(kg K) for each P and T. Cp_JkgK = aT_K + b
            alpha_pK (float, shape MxN): Thermal expansivity in 1/K for each P and T. alpha_pK = (2aT_K + b)/(aT_K^2 + bT_K + c)
            kTherm_WmK (float, shape MxN): Thermal conductivity in W/(mK) for each P and T. kTherm_WmK = c, a constant, over the specified range.
    """
    P_MPa, T_K = np.meshgrid(Plin_MPa, Tlin_K, indexing='ij')

    T_C = T_K - Constants.T0

    rho_kgm3 = (-2.3815e-4*T_C + 1.1843e-4*P_MPa + 0.92435) * 1e3
    Cp_JkgK = 3.19*T_K + 2150
    alpha_pK = (3.5697e-4*T_K + 0.2558)/(3.5697e-4*T_K**2 + 0.2558*T_K + 1612.8597)
    kTherm_WmK = np.zeros_like(P_MPa) + 0.5

    return rho_kgm3, Cp_JkgK, alpha_pK, kTherm_WmK


""" Dissociation temperatures for clathrates based on a degree-2 fit to
    the dissociation curve at low pressure and a logarithmic fit at higher pressures,
    for curves from Sloan (1998) as reported in Choukron et al. (2010):
    https://doi.org/10.1016/j.icarus.2009.08.011
"""
# For use when P_MPa < 2.567 and T < 273
TclathDissocLower_K = lambda P_MPa: 212.33820985 + 43.37319252*P_MPa - 7.83348412*P_MPa**2
# For use when P_MPa >= 2.567 and T >= 273
TclathDissocUpper_K = lambda P_MPa: -20.3058036 + 8.09637199*np.log(P_MPa/4.56717945e-16)

def ClathStableSloan1998(P_MPa, T_K):
    """ Returns a grid of boolean values to say whether fully occupied methane
        clathrates are stable at the given P,T conditions, based on the dissocation
        curves above.

        Args:
            P_MPa (float, shape N): Pressures in MPa
            T_K (float, shape M): Temperatures in K
        Returns:
            stable (int, shape NxM): Set to Constants.phaseClath (phase ID) if clathrates are stable
                at this P,T and 0 if they are not (if not stable, Ocean.EOS.fn_phase should be
                queried to determine the phase)
    """

    # Avoid log of 0 for surface pressure
    P_MPa[P_MPa==0] = 1e-12
    # Get (P,T) pairs relevant for each portion of the dissociation curve
    Plow_MPa, Tlow_K = np.meshgrid(P_MPa[P_MPa < 2.567], T_K, indexing='ij')
    Pupp_MPa, Tupp_K = np.meshgrid(P_MPa[P_MPa >= 2.567], T_K, indexing='ij')
    # Get evaluation points for lower and upper dissociation curves
    TdissocLower_K = TclathDissocLower_K(Plow_MPa)
    TdissocUpper_K = TclathDissocUpper_K(Pupp_MPa)
    # Assign clathrate phase ID to (P,T) points below the dissociation curves and 0 otherwise
    stableLow = np.zeros_like(Tlow_K).astype(np.int_)
    stableUpp = np.zeros_like(Tupp_K).astype(np.int_)
    stableLow[Tlow_K < TdissocLower_K] = Constants.phaseClath
    stableUpp[Tupp_K < TdissocUpper_K] = Constants.phaseClath
    stable = np.concatenate((stableLow, stableUpp), axis=0)

    return stable


def GetPbClath(Tb_K):
    """ Calculate the pressure consistent with Tb_K when clathrates are assumed
        to be in contact with the ocean, i.e. for Bulk.clathType = 'bottom' or 'whole'.

        Args:
            Tb_K (float): Clathrate layer bottom temperature in K
        Returns:
            PbClath_MPa (float): Bottom temperature consistent with dissociation curve
                pressure at Tb_K
    """
    if Tb_K < 273:
        TbZero_K = lambda P_MPa: Tb_K - TclathDissocLower_K(P_MPa)
        Pends_MPa = [0.0, 2.567]
    else:
        TbZero_K = lambda P_MPa: Tb_K - TclathDissocUpper_K(P_MPa)
        Pends_MPa = [2.567, Constants.PmaxLiquid_MPa]

    PbClath_MPa = spo.root_scalar(TbZero_K, bracket=Pends_MPa).root

    return PbClath_MPa