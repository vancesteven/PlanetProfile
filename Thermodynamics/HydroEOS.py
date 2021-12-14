import numpy as np
from seafreeze import seafreeze as SeaFreeze
from seafreeze import whichphase as WhichPhase

def FluidEOS(compstr, w_ppt, P_MPa, T_K):
    """ Returns mass density, heat capacity, and thermal expansivity based on thermodynamics
     from SeaFreeze and input pressure, temperature, salinity, and composition

        Args:
            compstr (string): Composition of dissolved salt
            w_ppt (float): Salinity of fluid in ppt
            P_MPa (float, shape N): Pressure of the fluid in MPa
            T_K (float, shape N): Temperature of the fluid in K
        Returns:
            rho_kgm3 (float, shape N): Density of fluid in kg/m^3
            Cp_JkgK (float, shape N): Heat capacity at constant pressure of fluid in J/kg/K
            alpha_pK (float, shape N): Thermal expansivity of fluid in K^-1
            VP_kms (float, shape N): P-wave seismic velocity in km/s
    """
    # Arrange input data into (P,T) value pair tuples compatible with SeaFreeze
    PTpairs = np.array([(P_MPa[i], T_K[i]) for i in range(np.size(P_MPa))], dtype='f,f').astype(object)

    if w_ppt == 0:
        seaOut = SeaFreeze(PTpairs, 'water1')
        rho_kgm3 = seaOut.rho
        Cp_JkgK = seaOut.Cp
        alpha_pK = seaOut.alpha
        VP_kms = seaOut.vel
    elif compstr == 'Seawater':
        raise ValueError('Unable to set FluidEOS. Seawater is not implemented yet.')
    elif compstr == 'NH3':
        raise ValueError('Unable to set FluidEOS. NH3 is not implemented yet.')
    elif compstr == 'MgSO4':
        raise ValueError('Unable to set FluidEOS. MgSO4 is not implemented yet.')
    elif compstr == 'NaCl':
        raise ValueError('Unable to set FluidEOS. NaCl is not implemented yet.')
    else:
        raise ValueError('Unable to set FluidEOS. compstr="'+compstr+'" but options are Seawater, NH3, MgSO4, and NaCl.')

    return rho_kgm3, Cp_JkgK, alpha_pK, VP_kms


def GetPhase(compstr, w_ppt, P_MPa, T_K):
    """ Get phase for single (scalar) (P,T) pair

        Args:
            compstr (string): Composition of dissolved salt
            w_ppt (float): Salinity of fluid in ppt
            P_MPa (float): Pressure of the fluid in MPa
            T_K (float): Temperature of the fluid in K
        Returns:
            phase (int): Ice/liquid phase ID for this P,T combo
    """
    PT = np.array([(P_MPa, T_K)], dtype='f,f').astype(object)
    if w_ppt == 0:
        phase = WhichPhase(PT)
    elif compstr == 'Seawater':
        raise ValueError('Unable to set FluidEOS. Seawater is not implemented yet.')
    elif compstr == 'NH3':
        raise ValueError('Unable to set FluidEOS. NH3 is not implemented yet.')
    elif compstr == 'MgSO4':
        raise ValueError('Unable to set FluidEOS. MgSO4 is not implemented yet.')
    elif compstr == 'NaCl':
        raise ValueError('Unable to set FluidEOS. NaCl is not implemented yet.')
    else:
        raise ValueError(
            'Unable to set FluidEOS. compstr="' + compstr + '" but options are Seawater, NH3, MgSO4, and NaCl.')

    return phase


def GetPfreeze(compstr, w_ppt, Tb_K, PfreezeLower_MPa=30, PfreezeUpper_MPa=300, PfreezeRes_MPa=0.1,
               Pguess=None, guessRange=5):
    """ Returns the pressure at which ice Ih melts based on temperature, salinity, and composition

        Args:
            compstr (string): Composition of dissolved salt
            w_ppt (float): Salinity of fluid in ppt
            Tb_K (float): Temperature of the phase transition in K
        Returns:
            Pfreeze_MPa (float): Pressure at the melting interface consistent with Tb_K
    """
    Psearch = np.arange(PfreezeLower_MPa, PfreezeUpper_MPa, PfreezeRes_MPa)
    # Arrange input data into (P,T) value pair tuples compatible with SeaFreeze
    PTsearchPairs = np.array([(P, Tb_K) for P in Psearch], dtype='f,f').astype(object)
    # Suggest a smaller range around a guessed value
    if Pguess is not None:
        DO_GUESS = True
        PguessRange = np.arange(Pguess-guessRange/2, Pguess+guessRange/2, PfreezeRes_MPa)
        PTguessPairs = np.array([(P, Tb_K) for P in PguessRange], dtype='f,f').astype(object)
    else:
        DO_GUESS = False
        GUESS_FAILED = True

    if w_ppt == 0:
        # Get phase of each P for the Tb_K value using SeaFreeze
        if DO_GUESS:
            searchPhases = WhichPhase(PTguessPairs)
            # Check if we failed to encounter a phase transition
            if np.all(searchPhases==searchPhases[0]): GUESS_FAILED = True
        if not DO_GUESS or GUESS_FAILED:
            searchPhases = WhichPhase(PTsearchPairs)
        # Find the first index for a phase that's not ice Ih or clathrates
        # (note that clathrates, with phase 30, are not yet implemented in SeaFreeze)
        try:
            indMelt = next((i[0] for i, val in np.ndenumerate(searchPhases) if val!=1 and val!=30))
        except StopIteration:
            raise ValueError('No melting pressure was found below '+str(PfreezeUpper_MPa)+' MPa'+
                             'for ice Ih/clathrates. Increase PfreezeUpper_MPa until one is found.')
        # Get the pressure of the ice Ih adjacent to the first non-Ih layer
        Pfreeze_MPa = Psearch[indMelt-1]
    elif compstr == 'Seawater':
        raise ValueError('Unable to GetPfreeze. Seawater is not implemented yet.')
    elif compstr == 'NH3':
        raise ValueError('Unable to GetPfreeze. NH3 is not implemented yet.')
    elif compstr == 'MgSO4':
        raise ValueError('Unable to GetPfreeze. MgSO4 is not implemented yet.')
    elif compstr == 'NaCl':
        raise ValueError('Unable to GetPfreeze. NaCl is not implemented yet.')
    else:
        raise ValueError('Unable to GetPfreeze. compstr="'+compstr+'" but options are Seawater, NH3, MgSO4, and NaCl.')

    return Pfreeze_MPa


def GetPfreezeHP(compstr, w_ppt, TbHP_K, phase, PfreezeHPLower_MPa=180, PfreezeHPUpper_MPa=900, PfreezeHPRes_MPa=0.5):
    """ Returns the pressure at which a high-pressure ice changes phase based on
     temperature, salinity, and composition

        Args:
            compstr (string): Composition of dissolved salt
            w_ppt (float): Salinity of fluid in ppt
            TbHP_K (float): Temperature of the phase transition in K
        Returns:
            PfreezeHP_MPa (float): Pressure at phase change interface
    """
    # Narrow the search range for narrow-stability HP ices
    if phase==3:
        PfreezeHPUpper_MPa = 360
    elif phase==5:
        PfreezeHPLower_MPa = 330

    Psearch = np.arange(PfreezeHPLower_MPa, PfreezeHPUpper_MPa, PfreezeHPRes_MPa)
    # Arrange input data into (P,T) value pair tuples compatible with SeaFreeze
    PTsearchPairs = np.array([(P, TbHP_K) for P in Psearch], dtype='f,f').astype(object)

    if w_ppt == 0:
        # Get phase of each P for the TbHP_K value using SeaFreeze
        searchPhases = WhichPhase(PTsearchPairs)
        # Find the first index of desired HP ice
        try:
            indIceHP = next((i[0] for i, val in np.ndenumerate(searchPhases) if val==phase))
        except StopIteration:
            raise ValueError('No ice '+str(PhaseConv(phase))+' was found within the range '+str(PfreezeHPLower_MPa)+
                  ' < P < '+str(PfreezeHPUpper_MPa)+' for Tb = '+str(TbHP_K)+'.')
        # Find the first index for a phase that's not the desired one after we encounter the desired phase
        try:
            indNotHP = next((i[0] for i, val in np.ndenumerate(searchPhases) if i[0]>indIceHP and val!=phase))
        except StopIteration:
            raise ValueError('No melting pressure was found below '+str(PfreezeHPUpper_MPa)+' MPa'+
                             'for ice'+PhaseConv(phase)+'. Increase PfreezeHPUpper_MPa until one is found.')
        # Get the pressure of the HP ice adjacent to the first layer with a non-matching phase
        PfreezeHP_MPa = Psearch[indNotHP-1]
    elif compstr == 'Seawater':
        raise ValueError('Unable to GetPfreezeHP. Seawater is not implemented yet.')
    elif compstr == 'NH3':
        raise ValueError('Unable to GetPfreezeHP. NH3 is not implemented yet.')
    elif compstr == 'MgSO4':
        raise ValueError('Unable to GetPfreezeHP. MgSO4 is not implemented yet.')
    elif compstr == 'NaCl':
        raise ValueError('Unable to GetPfreezeHP. NaCl is not implemented yet.')
    else:
        raise ValueError('Unable to GetPfreezeHP. compstr="'+compstr+'" but options are Seawater, NH3, MgSO4, and NaCl.')

    return PfreezeHP_MPa


def GetTfreeze(compstr, w_ppt, P_MPa, T_K, TfreezeRange_K=50, TfreezeRes_K=0.05):
    """ Returns the pressure at which the fluid freezes based on temperature, salinity, and composition

        Args:
            compstr (string): Composition of dissolved salt
            w_ppt (float): Salinity of fluid in ppt
            P_MPa (float): Pressure of the fluid in MPa
            T_K (float): Temperature of the fluid in K
        Returns:
            Tfreeze_K (float): Temperature of nearest higher-temperature phase transition between
                liquid and ice at this pressure
    """
    Tsearch = np.arange(T_K, T_K + TfreezeRange_K, TfreezeRes_K)
    # Arrange input data into (P,T) value pair tuples compatible with SeaFreeze
    PTsearchPairs = np.array([(P_MPa, T) for T in Tsearch], dtype='f,f').astype(object)

    if w_ppt == 0:
        # Get phase of each T for this pressure using SeaFreeze
        searchPhases = WhichPhase(PTsearchPairs)
        thisPhase = searchPhases[0]
        # Find the first index for a phase that doesn't match the input value
        try:
            indNotThis = next((i[0] for i, val in np.ndenumerate(searchPhases) if val!=thisPhase))
        except StopIteration:
            raise ValueError('No melting temperature was found above ' + str(T_K) + ' K' +
                             'for ice' + PhaseConv(thisPhase) + '. Increase TfreezeRange_K until one is found.')
        # Get the pressure of the HP ice adjacent to the first layer with a non-matching phase
        Tfreeze_K = Tsearch[indNotThis - 1]
    elif compstr == 'Seawater':
        raise ValueError('Unable to GetTfreeze. Seawater is not implemented yet.')
    elif compstr == 'NH3':
        raise ValueError('Unable to GetTfreeze. NH3 is not implemented yet.')
    elif compstr == 'MgSO4':
        raise ValueError('Unable to GetTfreeze. MgSO4 is not implemented yet.')
    elif compstr == 'NaCl':
        raise ValueError('Unable to GetTfreeze. NaCl is not implemented yet.')
    else:
        raise ValueError('Unable to GetTfreeze. compstr="'+compstr+'" but options are Seawater, NH3, MgSO4, and NaCl.')

    return Tfreeze_K


def GetIceThermo(P_MPa, T_K, phase):
    """ Thermodynamic calculation of ice density, heat capacity, and
        expansivity using SeaFreeze.

        Arguments:
            P_MPa (float, shape N): Pressure values of ice layers
            T_K (float, shape N): Temperature values of ice layers
            phase (int, shape N): ID for ice phase of each layer
        Returns:
            rho_kgm3 (float, shape N): Mass density for each layer
            Cp_JkgK (float, shape N): Heat capacity at constant pressure for each layer
            alpha_pK (float, shape N): Thermal expansivity for each layer
    """
    # Initialize outputs
    rho_kgm3, Cp_JkgK, alpha_pK = (np.zeros_like(P_MPa) for _ in range(3))
    # Identify which indices correspond to which phases
    indsLiquid, indsI, indsII, indsIII, indsV, indsVI, indsClath, _, _ = GetPhaseIndices(phase)
    # Organize (P,T) value pairs into a list of tuples compatible with SeaFreeze
    PTlist = np.array([(P_MPa[i], T_K[i]) for i in range(np.size(phase))], dtype='f,f').astype(object)

    # Call SeaFreeze for each phase type in one go
    if len(indsI) != 0:
        seaOut = SeaFreeze(PTlist[indsI], 'Ih')
        rho_kgm3[indsI] = seaOut.rho
        Cp_JkgK[indsI] = seaOut.Cp
        alpha_pK[indsI] = seaOut.alpha
    if len(indsLiquid) != 0:
        seaOut = SeaFreeze(PTlist[indsLiquid], 'water1')
        rho_kgm3[indsLiquid] = seaOut.rho
        Cp_JkgK[indsLiquid] = seaOut.Cp
        alpha_pK[indsLiquid] = seaOut.alpha
        print('WARNING: Unexpected liquids are present within ice layers. Check Steps.n settings.')
    if len(indsII) != 0:
        seaOut = SeaFreeze(PTlist[indsII], 'II')
        rho_kgm3[indsII] = seaOut.rho
        Cp_JkgK[indsII] = seaOut.Cp
        alpha_pK[indsII] = seaOut.alpha
    if len(indsIII) != 0:
        seaOut = SeaFreeze(PTlist[indsIII], 'III')
        rho_kgm3[indsIII] = seaOut.rho
        Cp_JkgK[indsIII] = seaOut.Cp
        alpha_pK[indsIII] = seaOut.alpha
    if len(indsV) != 0:
        seaOut = SeaFreeze(PTlist[indsV], 'V')
        rho_kgm3[indsV] = seaOut.rho
        Cp_JkgK[indsV] = seaOut.Cp
        alpha_pK[indsV] = seaOut.alpha
    if len(indsVI) != 0:
        seaOut = SeaFreeze(PTlist[indsVI], 'VI')
        rho_kgm3[indsVI] = seaOut.rho
        Cp_JkgK[indsVI] = seaOut.Cp
        alpha_pK[indsVI] = seaOut.alpha
    if len(indsClath) != 0:
        seaOut = SeaFreeze(PTlist[indsClath], 'Clath')
        rho_kgm3[indsClath] = seaOut.rho
        Cp_JkgK[indsClath] = seaOut.Cp
        alpha_pK[indsClath] = seaOut.alpha

    return rho_kgm3, Cp_JkgK, alpha_pK


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
    elif phase == 30:
        phaseStr = 'Clath'
    else:
        raise ValueError('PhaseConv does not have a definition for phase ID '+str(phase)+'.')

    return phaseStr


def GetPhaseIndices(phase):
    """ Get indices for each phase of ice/liquid

        Args:
            phase (int, shape N)
        Returns:
            indsLiquid, indsI, ... indsFe (int, shape 0-M): lists of indices corresponding to each phase.
                Variable length.
    """
    indsLiquid = [i for i, val in enumerate(phase) if val == 0]
    indsI = [i for i, val in enumerate(phase) if val == 1]
    indsII = [i for i, val in enumerate(phase) if val == 2]
    indsIII = [i for i, val in enumerate(phase) if val == 3]
    indsV = [i for i, val in enumerate(phase) if val == 5]
    indsVI = [i for i, val in enumerate(phase) if val == 6]
    indsClath = [i for i, val in enumerate(phase) if val == 30]
    indsSil = [i for i, val in enumerate(phase) if val == 50]
    indsFe = [i for i, val in enumerate(phase) if val == 100]

    return indsLiquid, indsI, indsII, indsIII, indsV, indsVI, indsClath, indsSil, indsFe