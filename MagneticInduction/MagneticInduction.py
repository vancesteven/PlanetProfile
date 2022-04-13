import os
import numpy as np
import logging as log
import scipy.interpolate as spi
from scipy.integrate import solve_ivp as ODEsolve
from collections.abc import Iterable
from Utilities.defineStructs import Constants, EOSlist
from config import Excitations
from moonInduction.asymmetry_funcs import read_Benm as GetBenm

def MagneticInduction(Planet, Params):
    """ Calculate induced magnetic moments for the body and prints them to disk.

        Requires Planet attributes:
            Magnetic.inductType
        Sets Planet attributes:

    """

    # Set Magnetic struct layer arrays as we need for induction calculations
    Planet = SetupInduction(Planet, Params)

    # Calculate induced magnetic moments
    Planet = CalcInducedMoments(Planet, Params)

    return Planet


def CalcInducedMoments(Planet, Params):
    """ Calculate induced magnetic moments based on conductivity profile,
        possible asymmetric shape, and excitation moments for this body.

        Sets Planet attributes:
            Magnetic.Ae, Magnetic.Amp, Magnetic.phase, Magnetic.Binm_nT
    """
    # Get wavenumbers for each layer for each frequency
    k_pm = np.array([np.sqrt(1j * Constants.mu0 * omega * Planet.Magnetic.sigmaLayers_Sm)
                     for omega in Planet.Magnetic.omegaExc_radps])

    if Planet.Magnetic.inductMethod == 'Eckhardt1963' or Planet.Magnetic.inductMethod == 'numeric':
        if Params.INCLUDE_ASYM:
            log.warning('Asymmetry can only be modeled with Srivastava1966 (layer) method, but ' +
                        'Eckhardt1963 (numeric) method is selected. Asymmetry will be ignored.')

        Planet.Magnetic.Aen = np.zeros((Planet.Magnetic.nExc, Planet.Magnetic.nprmMax+1), dtype=np.complex_)
        for nprm in range(1, Planet.Magnetic.nprmMax+1):
            Q = np.array([SolveForQ(nprm, k_pm[i,:], Planet.Magnetic.rChange_m, Planet.Bulk.R_m,
                                    Params.EckhardtSolveMethod, rMin=Params.rMinODE)
                                    for i,omega in enumerate(Planet.Magnetic.omegaExc_radps)])
            Planet.Magnetic.Aen[:,nprm] = Q * (nprm + 1) / nprm
            Planet.Magnetic.Binm_nT[:,:,nprm,:] = [Planet.Magnetic.Benm_nT[i,:,nprm,:] * Planet.Magnetic.Aen[i,nprm]
                                                   for i in range(Planet.Magnetic.nExc)]

    elif Planet.Magnetic.inductMethod == 'Srivastava1966' or Planet.Magnetic.inductMethod == 'layer':
        pass
    else:
        raise ValueError(f'Induction method "{Planet.Magnetic.inductMethod}" not defined.')

    Planet.Magnetic.Amp =    np.abs(  Planet.Magnetic.Aen[:,1])
    Planet.Magnetic.phase = -np.angle(Planet.Magnetic.Aen[:,1], deg=True)

    return Planet


def SolveForQ(n, kBlw_pm, rBds_m, R_m, solveMethod, rMin=1e3):
    dQdr = fn_dQdr(n, kBlw_pm, rBds_m)
    Q = ODEsolve(dQdr, (rMin, rBds_m[-1]), [0j], method=solveMethod)
    Q = Q['y'][0,-1] * (rBds_m[-1]/R_m)**(n+2)
    return Q


class fn_dQdr:
    def __init__(self, n, kBlw_pm, rBds_m):
        self.n = n
        self.kBlw_pm = kBlw_pm
        self.rBds_m = rBds_m

    def fn_k(self, r_m):
        if not isinstance(r_m, Iterable):
            r_m = np.array([r_m])
        iNextAbove = [next((i for i, rBd_m in enumerate(self.rBds_m) if rBd_m >= ri_m)) for ri_m in r_m]
        return self.kBlw_pm[iNextAbove]

    def __call__(self, r, Q):
        k = self.fn_k(r)
        return -k**2 * r * (self.n+1) / (2*self.n + 1) / self.n * (Q - self.n/(self.n + 1))**2 - (2*self.n + 1) / r * Q


def SetupInduction(Planet, Params):
    """ Reconfigure layer boundaries and conductivities into a format
        usable by magnetic induction calculation functions.
        Optionally also identify asymmetric shape information from gravity.

        Requires Planet attributes:
        Sets Planet attributes:
            Magnetic.rSigChange_m, Magnetic.sigmaLayers_Sm, Magnetic.asymShape
    """

    # Get excitation spectrum
    Planet.Magnetic.Texc_hr, Planet.Magnetic.omegaExc_radps, Planet.Magnetic.Benm_nT, Planet.Magnetic.B0_nT \
        = GetBexc(Planet.name, Planet.Magnetic.SCera, Planet.Magnetic.extModel, Params.excSelectionCalc,
                  nprmMax=Planet.Magnetic.nprmMax, pMax=Planet.Magnetic.pMax)

    # Initialize Binm array to have the same shape and data type as Benm
    Planet.Magnetic.Binm_nT = np.zeros_like(Planet.Magnetic.Benm_nT)

    # Reconfigure layer conducting boundaries as needed.
    # For inductOtype == 'sigma', we have already set these arrays.
    if not Params.inductOtype == 'sigma':
        # Append optional ionosphere
        if Planet.Magnetic.ionosBounds_m is None:
            zIonos_m = []
            sigmaIonos_Sm = []
        else:
            zIonos_m = Planet.Magnetic.ionosBounds_m
            sigmaIonos_Sm = Planet.Magnetic.sigmaIonosPedersen_Sm
            # Allow special case for specifying an ionosphere with 1 conductivity and 2 bounds, when
            # there is a substantial neutral atmosphere and the conducting region is at altitude
            # (e.g. for Triton)
            if np.size(Planet.Magnetic.sigmaIonosPedersen_Sm) == 1 and np.size(Planet.Magnetic.ionosBounds_m) == 2:
                sigmaIonos_Sm = np.append(0, sigmaIonos_Sm)
        # Flip arrays to be in radial ascending order as needed in induction calculations, then add ionosphere
        rLayers_m = np.append(np.flip(Planet.r_m), zIonos_m)
        sigmaInduct_Sm = np.append(np.flip(Planet.sigma_Sm), sigmaIonos_Sm)

        # Eliminate NaN values and 0 values, assigning them to a default minimum
        sigmaInduct_Sm[np.logical_or(np.isnan(sigmaInduct_Sm), sigmaInduct_Sm == 0)] = Constants.sigmaDef_Sm
        # Set low conductivities to all be the same default value so we can shrink them down to single layers
        sigmaInduct_Sm[sigmaInduct_Sm < Constants.sigmaMin_Sm] = Constants.sigmaDef_Sm

        # Optionally, further reduce computational overhead by shrinking the number of ocean layers modeled
        if Params.REDUCED_INDUCT:
            indsLiq = np.where(np.flip(Planet.phase) == 0)[0]
            if np.size(indsLiq) > 0 and not np.all(np.diff(indsLiq) == 1):
                log.warning('HP ices found in ocean while REDUCED_INDUCT is True. They will be ignored ' +
                            'in the interpolation.')
            if np.size(indsLiq) <= Params.nIntL:
                log.warning(f'Only {np.size(indsLiq)} layers found in ocean, but number of layers to ' +
                            f'interpolate over is {Params.nIntL}. This profile will be not be reduced.')
            else:
                # Get radius values from D/nIntL above the seafloor to the ice shell
                rBot_m = rLayers_m[indsLiq[0]-1]
                rTop_m = rLayers_m[indsLiq[-1]]
                rOcean_m = np.linspace(rBot_m, rTop_m, Params.nIntL+1)[1:]
                # Interpolate the conductivities corresponding to those radii
                sigmaOcean_Sm = spi.interp1d(rLayers_m[indsLiq], sigmaInduct_Sm[indsLiq], kind=Params.oceanInterpMethod)(rOcean_m)
                # Stitch together the r and σ arrays with the new ocean values
                rLayers_m = np.concatenate((rLayers_m[0:indsLiq[0]], rOcean_m, rLayers_m[indsLiq[-1]+1:]))
                sigmaInduct_Sm = np.concatenate((sigmaInduct_Sm[0:indsLiq[0]], sigmaOcean_Sm, sigmaInduct_Sm[indsLiq[-1]+1:]))

        # Get the indices of layers just below where changes happen
        iChange = [i for i,sig in enumerate(sigmaInduct_Sm) if sig != np.append(sigmaInduct_Sm, np.nan)[i+1]]
        Planet.Magnetic.sigmaLayers_Sm = sigmaInduct_Sm[iChange]
        Planet.Magnetic.rChange_m = rLayers_m[iChange]

    Planet.Magnetic.nBds = np.size(Planet.Magnetic.rChange_m)
    Planet.Magnetic.nExc = np.size(Planet.Magnetic.Texc_hr)

    # Set asymmetric shape if applicable
    if Params.INCLUDE_ASYM:
        Planet.Magnetic.pMax = 2
        # The next line is a placeholder for now.
        Planet.Magnetic.asymShape = np.zeros((Planet.Magnetic.nBds, 2, Planet.Magnetic.pMax+1, Planet.Magnetic.pMax+1))
    else:
        Planet.Magnetic.pMax = 0
        Planet.Magnetic.asymShape = np.zeros((Planet.Magnetic.nBds, 2, Planet.Magnetic.pMax+1, Planet.Magnetic.pMax+1))

    return Planet


def GetBexc(bodyname, era, model, excSelection, nprmMax=1, pMax=0):
    """ Read in magnetic excitation information, including oscillation
        frequencies/periods and complex amplitudes and phases (moments).

        Args:
            bodyname (str): Body name.
            era (str): Spacecraft era over which excitation spectrum was evaluated.
                Read from Be1xyz filename after the body name.
            model (str): External field model applied to the body that was used
                to evaluate the excitation moments. Read from the filename after the
                SC era.
            excSelection (dict): Boolean flags for the major excitations identifying which
                ones should be included in calculations. Keys must match those in
                ExcitationsList in config.
            nprmMax = 1 (int): Maximum n' value to use for Benm. 1 is uniform field
                applied by the parent planet (uniform across the body).
            pMax = 0 (int): Maximum p value to use for asymmetric shape chipq. 0 is
                spherically symmetric, 2 includes gravitational perturbations.
        Returns:
            Texc_hr (float, shape N): Periods of oscillation for magnetic excitation
                spectrum to model for this body in hr.
            omegaExc_radps (float, shape N): Angular frequencies of oscillation for
                the excitation spectrum to model in rad/s.
            Benm_nT (complex, shape N x 2 x nprmMax+pMax+1 x nprmMax+pMax+1): Complex
                excitation moments applied to the body relative to J2000 epoch. First index
                is over osc. period, second is sign of m' (0 = positive, 1 = negative),
                third index is n', fourth index is abs(m').
            B0_nT (float, shape 3): Constant background field applied to the body by the
                parent planet.
    """
    BeLabel = f'{bodyname}{era}{model}{excSelection}'

    if BeLabel in EOSlist.loaded.keys():
        log.debug(f'{bodyname} excitation spectrum for {model} model and {era} era already loaded. Reusing existing.')
        Texc_hr, omegaExc_radps, Benm_nT, B0_nT = EOSlist.loaded[BeLabel]
    else:
        log.debug(f'Loading {bodyname} excitation spectrum for {model} model and {era} era.')
        if era is None and model is None:
            ID = ''
        else:
            ID = '_'.join([f'{era}', f'{model}']).replace('_None', '').replace('None_', '')
        fPath = os.path.join(bodyname, 'induction')
        inpTexc_hr, inpBenm_nT, B0_nT = GetBenm(nprmMax, pMax, bodyname=bodyname, fpath=fPath, model=ID)

        nPeaks = sum(excSelection.values())
        Texc_hr = np.zeros(nPeaks)
        Benm_nT = np.zeros((nPeaks, 2, nprmMax+pMax+1, nprmMax+pMax+1), dtype=np.complex_)
        BeList = Excitations(bodyname)
        # Include in the excitation spectrum only the periods specified in config.py
        iPeak = 0
        for oscillation, included in excSelection.items():
            if included and BeList[oscillation] is not None:
                Texc_hr[iPeak] = inpTexc_hr[np.round(inpTexc_hr,2) == round(BeList[oscillation],2)]
                Benm_nT[iPeak, ...] = inpBenm_nT[inpTexc_hr == Texc_hr[iPeak], ...]
                iPeak += 1

        omegaExc_radps = 2*np.pi / Texc_hr / 3600

        EOSlist.loaded[BeLabel] = (Texc_hr, omegaExc_radps, Benm_nT, B0_nT)
        EOSlist.ranges[BeLabel] = Texc_hr

    return Texc_hr, omegaExc_radps, Benm_nT, B0_nT


def Benm2absBexyz(Benm):
    A1 = np.sqrt(2*np.pi/3)
    Bex = np.abs(-1 /2/A1 * (Benm[:,1,1,1] - Benm[:,0,1,1]))
    Bey = np.abs( 1j/2/A1 * (Benm[:,1,1,1] + Benm[:,0,1,1]))
    Bez = np.abs(-1 /np.sqrt(2)/A1 * Benm[:,0,1,0])
    return Bex, Bey, Bez


def ReloadInduction(Planet, Params):
    """ Reload induced magnetic moments that have been printed to disk.

        Sets Planet attributes:
            Magnetic.Binm_nT
    """
    Planet.Magnetic.Binm_nT = np.loadtxt(Params.inducedMomentsFile, skiprows=1, unpack=False)

    return Planet