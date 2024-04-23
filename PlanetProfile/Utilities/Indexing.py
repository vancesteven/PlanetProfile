"""
Indexing: Functions for converting between various means of indexing material phases
"""

import numpy as np
import logging
from PlanetProfile.Utilities.defineStructs import Constants

# Assign logger
log = logging.getLogger('PlanetProfile')


def PhaseConv(phase, PORE=False, liq='water1'):
    """ Convert phase integers into strings compatible with SeaFreeze

        Arguments:
            phase (int): ID of phase for each layer
            PORE = False (bool): Whether to return phase of pore material instead of matrix
            liq = 'water1' (string): SeaFreeze-compatible composition for liquid phase
        Returns:
            phaseStr (string): Corresponding string for each phase ID
    """
    if phase == 0:
        phaseStr = liq
    elif abs(phase) == 1:
        if PORE and phase < 0:
            phaseStr = liq
        else:
            phaseStr = 'Ih'
    elif abs(phase) == 2:
        phaseStr = 'II'
    elif abs(phase) == 3:
        phaseStr = 'III'
    elif abs(phase) == 5:
        phaseStr = 'V'
    elif abs(phase) == 6:
        phaseStr = 'VI'
    elif abs(phase) == Constants.phaseClath:
        if PORE and phase < 0:
            phaseStr = liq
        else:
            phaseStr = 'Clath'
    elif phase >= Constants.phaseSil and phase < Constants.phaseSil+10:
        if PORE and phase != Constants.phaseSil:
            phaseStr = PhaseConv(phase % 10, PORE=False, liq=liq)
        else:
            phaseStr = 'Sil'
    elif phase >= Constants.phaseFe:
        phaseStr = 'Fe'
    else:
        raise ValueError(f'PhaseConv does not have a definition for phase ID {phase:d}.')

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
        raise ValueError(f'PhaseInv does not have a definition for phase string "{phaseStr}".')

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

    indsLiquid = np.where(phase==0)[0]
    indsIceI = np.where(phase==1)[0]
    indsIceIwet = np.where(phase==-1)[0]
    indsIceII = np.where(phase==2)[0]
    indsIceIIund = np.where(phase==-2)[0]
    indsIceIII = np.where(phase==3)[0]
    indsIceIIIund = np.where(phase==-3)[0]
    indsIceV = np.where(phase==5)[0]
    indsIceVund = np.where(phase==-5)[0]
    indsIceVI = np.where(phase==6)[0]
    indsIceVIund = np.where(phase==-6)[0]
    indsClath = np.where(phase==Constants.phaseClath)[0]
    indsClathWet = np.where(phase==-Constants.phaseClath)[0]
    indsSil = np.where(np.logical_and(phase>=Constants.phaseSil, phase<Constants.phaseSil+10))[0]
    indsSilLiq = np.where(phase==Constants.phaseSil)[0]
    indsSilI = np.where(phase==Constants.phaseSil+1)[0]
    indsSilII = np.where(phase==Constants.phaseSil+2)[0]
    indsSilIII = np.where(phase==Constants.phaseSil+3)[0]
    indsSilV = np.where(phase==Constants.phaseSil+5)[0]
    indsSilVI = np.where(phase==Constants.phaseSil+6)[0]
    indsFe = np.where(phase>=Constants.phaseFe)[0]

    return indsLiquid, indsIceI, indsIceIwet, indsIceII, indsIceIIund, indsIceIII, indsIceIIIund, indsIceV, indsIceVund, \
               indsIceVI, indsIceVIund, indsClath, indsClathWet, indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, \
               indsSilV, indsSilVI, indsFe
