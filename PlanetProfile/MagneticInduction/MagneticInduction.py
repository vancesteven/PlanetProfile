import os
import numpy as np
import logging
import scipy.interpolate as spi
import scipy.special as sps
import spiceypy as spice
from scipy.integrate import solve_ivp as ODEsolve
from collections.abc import Iterable
from glob import glob as FilesMatchingPattern
from scipy.io import savemat, loadmat
from PlanetProfile import _Test, _Defaults
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist
from PlanetProfile.MagneticInduction.Moments import Excitations
from PlanetProfile.GetConfig import FigMisc
from MoonMag.asymmetry_funcs import read_Benm as GetBenm, BiList as BiAsym, get_chipq_from_CSpq as GeodesyNorm2chipq, \
    get_all_Xid as LoadXid, get_rsurf as GetrSurf
from MoonMag.symmetry_funcs import InducedAeList as AeList

# Assign logger
log = logging.getLogger('PlanetProfile')

def MagneticInduction(Planet, Params, fNameOverride=None):
    """ Calculate induced magnetic moments for the body and prints them to disk.
    """

    SKIP = False
    if Planet.Do.VALID and Params.CALC_NEW_INDUCT:
        # Set Magnetic struct layer arrays as we need for induction calculations
        Planet, Params = SetupInduction(Planet, Params)

        # Skip remaining calcs if we couldn't load excitation moments
        if Planet.Magnetic.Benm_nT is not None:
            # Calculate induced magnetic moments
            Planet = CalcInducedMoments(Planet, Params)

            # Report progress for sigma inductograms
            if Params.INDUCTOGRAM_IN_PROGRESS and Params.Induct.inductOtype == 'sigma':
                if Planet.index % 10 == 0:
                    log.profile(f'Point {Planet.index}/{Params.nModels} complete.')

    elif Planet.Do.VALID:
        # Reload induced moments from disk
        if fNameOverride is None:
            momentsFile = Params.DataFiles.inducedMomentsFile
        else:
            momentsFile = fNameOverride

        if os.path.isfile(momentsFile):
            Planet = ReloadMoments(Planet, momentsFile)
        else:
            log.warning(f'CALC_NEW_INDUCT is False, but {momentsFile} was not found. ' +
                        f'Skipping induction calculations.')
            SKIP = True
            
    else:
        SKIP = True

    if SKIP:
        if Planet.Magnetic.ionosBounds_m is None:
            Planet.Magnetic.ionosBounds_m = [np.nan]
        if Planet.Magnetic.sigmaIonosPedersen_Sm is None:
            Planet.Magnetic.sigmaIonosPedersen_Sm = [np.nan]
        Planet.Magnetic.calcedExc = []
    else:
        if not (Params.DO_INDUCTOGRAM or Params.DO_EXPLOREOGRAM):
            if Params.PLOT_MAG_SPECTRUM:
                Planet = FourierSpectrum(Planet, Params)

            if Params.CALC_ASYM:
                if Params.CALC_NEW_ASYM:
                    if not Params.CALC_NEW_INDUCT:
                        Planet, Params = SetupInduction(Planet, Params)
                    Planet = CalcAsymContours(Planet, Params)
                else:
                    Planet = ReloadAsym(Planet, Params)
    
    # Must return both Planet and Params in order to use common infrastructure
    # for unpacking parallel runs
    return Planet, Params


def CalcInducedMoments(Planet, Params):
    """ Calculate induced magnetic moments based on conductivity profile,
        possible asymmetric shape, and excitation moments for this body.

        Sets Planet attributes:
            Magnetic.Ae, Magnetic.Amp, Magnetic.phase, Magnetic.Binm_nT
    """
    # Get lists of n and m values for linearizing Binm after we calculate. Also needed for
    # asymmetric layer calculations, so we do this first.
    nMax = Planet.Magnetic.nprmMax + Planet.Magnetic.pMax
    Planet.Magnetic.nLin = [n for n in range(1, nMax+1) for _ in range(-n, n+1)]
    Planet.Magnetic.mLin = [m for n in range(1, nMax+1) for m in range(-n, n+1)]
    Planet.Magnetic.nprmLin = [n for n in range(1, Planet.Magnetic.nprmMax+1) for _ in range(-n, n+1)]
    Planet.Magnetic.mprmLin = [m for n in range(1, Planet.Magnetic.nprmMax+1) for m in range(-n, n+1)]
    Nnm = np.size(Planet.Magnetic.nLin)

    Planet.Magnetic.Aen = np.zeros((Planet.Magnetic.nExc, Planet.Magnetic.nprmMax+1), dtype=np.complex_)
    Planet.Magnetic.BinmLin_nT = np.zeros((Planet.Magnetic.nExc, Nnm), dtype=np.complex_)
    if Planet.Magnetic.inductMethod == 'Eckhardt1963' or Planet.Magnetic.inductMethod == 'numeric':
        if Params.Sig.INCLUDE_ASYM:
            log.warning('Asymmetry can only be modeled with Srivastava1966 (layer) method, but ' +
                        'Eckhardt1963 (numeric) method is selected. Asymmetry will be ignored.')

        # Get wavenumbers for each layer for each frequency
        k_pm = np.array([np.sqrt(1j * Constants.mu0 * omega * Planet.Magnetic.sigmaLayers_Sm)
                     for omega in Planet.Magnetic.omegaExc_radps])

        for nprm in range(1, Planet.Magnetic.nprmMax+1):
            Q = np.array([SolveForQ(nprm, k_pm[i,:], Planet.Magnetic.rSigChange_m, Planet.Bulk.R_m,
                                    Params.Induct.EckhardtSolveMethod, rMin=Params.Induct.rMinODE)
                                    for i,omega in enumerate(Planet.Magnetic.omegaExc_radps)])
            Planet.Magnetic.Aen[:,nprm] = Q * (nprm+1) / nprm
            Planet.Magnetic.Binm_nT[:,:,nprm,:] = [Planet.Magnetic.Benm_nT[i,:,nprm,:] * Planet.Magnetic.Aen[i,nprm]
                                                   for i in range(Planet.Magnetic.nExc)]

        Planet.Magnetic.Amp = np.abs(Planet.Magnetic.Aen[:, 1])
        Planet.Magnetic.phase = -np.angle(Planet.Magnetic.Aen[:, 1], deg=True)

    elif Planet.Magnetic.inductMethod == 'Srivastava1966' or Planet.Magnetic.inductMethod == 'layer':
        # Evaluate 1st-order complex response amplitudes
        Planet.Magnetic.Aen[:,1], Planet.Magnetic.Amp, AeArg \
            = AeList(Planet.Magnetic.rSigChange_m, Planet.Magnetic.sigmaLayers_Sm,
                     Planet.Magnetic.omegaExc_radps, 1/Planet.Bulk.R_m, nn=1,
                     writeout=False, do_parallel=not (Params.INDUCTOGRAM_IN_PROGRESS or Params.DO_EXPLOREOGRAM))
        Planet.Magnetic.phase = -np.degrees(AeArg)
        for n in range(2, Planet.Magnetic.nprmMax):
            Planet.Magnetic.Aen[:,n], _, _ \
                = AeList(Planet.Magnetic.rSigChange_m, Planet.Magnetic.sigmaLayers_Sm,
                         Planet.Magnetic.omegaExc_radps, 1/Planet.Bulk.R_m, nn=n,
                         writeout=False, do_parallel=not (Params.INDUCTOGRAM_IN_PROGRESS or Params.DO_EXPLOREOGRAM))

        if Params.Sig.INCLUDE_ASYM:
            # Use a separate function for evaluating asymmetric induced moments, as Binm is not as simple as
            # Aen * Benm for this case.
            Planet.Magnetic.Binm_nT = BiAsym(Planet.Magnetic.rSigChange_m, Planet.Magnetic.sigmaLayers_Sm,
                                             Planet.Magnetic.omegaExc_radps, Planet.Magnetic.asymShape_m,
                                             Planet.Magnetic.gravShape_m, Planet.Magnetic.Benm_nT, 1/Planet.Bulk.R_m,
                                             Planet.Magnetic.nLin, Planet.Magnetic.mLin, Planet.Magnetic.pMax,
                                             nprm_max=Planet.Magnetic.nprmMax, writeout=False,
                                             do_parallel=not (Params.INDUCTOGRAM_IN_PROGRESS or Params.DO_EXPLOREOGRAM),
                                             Xid=Planet.Magnetic.Xid)
        else:
            # Multiply complex response by Benm to get Binm for spherically symmetric case
            for iPeak in range(Planet.Magnetic.nExc):
                for n in Planet.Magnetic.nprmLin:
                    for m in Planet.Magnetic.mprmLin:
                        Planet.Magnetic.Binm_nT[iPeak, int(m<0), n, m] \
                            = n/(n+1) * Planet.Magnetic.Benm_nT[iPeak, int(m<0), n, m] * Planet.Magnetic.Aen[iPeak, n]
    else:
        raise ValueError(f'Induction method "{Planet.Magnetic.inductMethod}" not defined.')

    # Get linear lists of Binm for more convenient post-processing
    Planet.Magnetic.BinmLin_nT[:,:Nnm] \
        = np.array([[Planet.Magnetic.Binm_nT[iPeak,int(m<0),n,m]
                   for n, m in zip(Planet.Magnetic.nLin, Planet.Magnetic.mLin)]
                   for iPeak in range(Planet.Magnetic.nExc)])

    # Get surface strength in IAU components for plotting, with conjugate phase to match
    # Zimmer et al. (2000) phase convention
    Bex, Bey, Bez = Benm2absBexyz(Planet.Magnetic.Benm_nT)
    Planet.Magnetic.Bi1xyz_nT = {
        'x': Bex * np.conj(Planet.Magnetic.Aen[:,1]),
        'y': Bey * np.conj(Planet.Magnetic.Aen[:,1]),
        'z': Bez * np.conj(Planet.Magnetic.Aen[:,1])
    }

    Planet.Magnetic.calcedExc = [key for key, CALCED in Params.Induct.excSelectionCalc.items()
                                 if CALCED and Excitations.Texc_hr[Planet.bodyname][key] is not None]
    # Save calculated magnetic moments to disk
    if not Params.NO_SAVEFILE:
        saveDict = {
            'Benm_nT': Planet.Magnetic.Benm_nT,
            'Binm_nT': Planet.Magnetic.Binm_nT,
            'omegaExc_radps': Planet.Magnetic.omegaExc_radps,
            'BinmLin_nT': Planet.Magnetic.BinmLin_nT,
            'nLin': Planet.Magnetic.nLin,
            'mLin': Planet.Magnetic.mLin,
            'nprmLin': Planet.Magnetic.nprmLin,
            'mprmLin': Planet.Magnetic.mprmLin,
            'Bi1x_nT': Planet.Magnetic.Bi1xyz_nT['x'],
            'Bi1y_nT': Planet.Magnetic.Bi1xyz_nT['y'],
            'Bi1z_nT': Planet.Magnetic.Bi1xyz_nT['z'],
            'Aen': Planet.Magnetic.Aen,
            'asymShape_m': Planet.Magnetic.asymShape_m,
            'ionosBounds_m': Planet.Magnetic.ionosBounds_m,
            'calcedExc': Planet.Magnetic.calcedExc
        }
        savemat(Params.DataFiles.inducedMomentsFile, saveDict)
        log.debug(f'Saved induced moments to file: {Params.DataFiles.inducedMomentsFile}')

    return Planet


def ReloadMoments(Planet, momentsFile):
    """ Reload induced moments from disk """
    
    reload = loadmat(momentsFile)
    Planet.Magnetic.Benm_nT = reload['Benm_nT']
    Planet.Magnetic.Binm_nT = reload['Binm_nT']
    Planet.Magnetic.omegaExc_radps = reload['omegaExc_radps'][0]
    Planet.Magnetic.BinmLin_nT = reload['BinmLin_nT']
    Planet.Magnetic.nLin = reload['nLin'][0]
    Planet.Magnetic.mLin = reload['mLin'][0]
    Planet.Magnetic.nprmLin = reload['nprmLin'][0]
    Planet.Magnetic.mprmLin = reload['mprmLin'][0]
    Planet.Magnetic.Bi1xyz_nT = {
        'x': reload['Bi1x_nT'][0],
        'y': reload['Bi1y_nT'][0],
        'z': reload['Bi1z_nT'][0]
    }
    Planet.Magnetic.Aen = reload['Aen']
    Planet.Magnetic.asymShape_m = reload['asymShape_m']
    Planet.Magnetic.ionosBounds_m = reload['ionosBounds_m'][0]
    Planet.Magnetic.calcedExc = np.char.strip(reload['calcedExc'])
    Planet.Magnetic.nExc = np.size(Planet.Magnetic.omegaExc_radps)
    
    return Planet


def CalcAsymContours(Planet, Params):

    if Params.Sig.INCLUDE_ASYM:
        # Make a calculation for each loaded asymmetry file and the gravity shape
        Planet.Magnetic.asymDevs_km = np.zeros((Planet.Magnetic.nAsymBds, FigMisc.nLatMap, FigMisc.nLonMap))
        for i, iLayer, zMean_km in zip(range(Planet.Magnetic.nAsymBds), Planet.Magnetic.iAsymBds, Planet.Magnetic.zMeanAsym_km):
            log.debug(f'Calculating topographic data for zMean = {zMean_km:.1f} km depth with {360/FigMisc.nLonMap:.1f}° resolution. This may take some time.')
            if i == Planet.Magnetic.nAsymBds - 1:
                rSurf_m = np.real(GetrSurf(Planet.Magnetic.pLin, Planet.Magnetic.qLin, Planet.Magnetic.gravShape_m[iLayer, ...],
                                           Planet.Bulk.R_m, FigMisc.thetaMap_rad, FigMisc.phiMap_rad,
                                           do_parallel=Params.DO_PARALLEL and not (Params.INDUCTOGRAM_IN_PROGRESS or Params.DO_EXPLOREOGRAM)))
            else:
                rSurf_m = np.real(GetrSurf(Planet.Magnetic.pLin, Planet.Magnetic.qLin, Planet.Magnetic.asymShape_m[iLayer, ...],
                                           Planet.Bulk.R_m - zMean_km*1e3, FigMisc.thetaMap_rad, FigMisc.phiMap_rad,
                                           do_parallel=Params.DO_PARALLEL and not (Params.INDUCTOGRAM_IN_PROGRESS or Params.DO_EXPLOREOGRAM)))

            # Plot as depth
            Planet.Magnetic.asymDevs_km[i, ...] = (Planet.Bulk.R_m - rSurf_m) / 1e3
            if zMean_km <= 0:
                # For asymmetry at the surface and above, plot as elevation/altitudes and not depth
                Planet.Magnetic.asymDevs_km[i, ...] = -1 * Planet.Magnetic.asymDevs_km[i, ...]

        # Save calculated boundary deviations to disk
        if not Params.NO_SAVEFILE:
            saveDict = {
                'nAsymBds': Planet.Magnetic.nAsymBds,
                'iAsymBds': Planet.Magnetic.iAsymBds,
                'zMeanAsym_km': Planet.Magnetic.zMeanAsym_km,
                'asymShape_m': Planet.Magnetic.asymShape_m,
                'gravShape_m': Planet.Magnetic.gravShape_m,
                'asymDevs_km': Planet.Magnetic.asymDevs_km,
                'latMap_deg': FigMisc.latMap_deg,
                'lonMap_deg': FigMisc.lonMap_deg
            }
            savemat(Params.DataFiles.asymFile, saveDict)
            log.debug(f'Saved asymmetric boundary deviations to file: {Params.DataFiles.asymFile}')

    return Planet


def ReloadAsym(Planet, Params, fNameOverride=None):
    # Reload calculated boundary deviations from disk

    if fNameOverride is not None:
        asymFile = fNameOverride
    else:
        asymFile = Params.DataFiles.asymFile

    if not os.path.isfile(asymFile):
        log.error(f'Asymmetric shape file {asymFile} not found. Asymmetric shape ' +
                  'plotting will be skipped.')
        Params.CALC_ASYM = False
    else:
        reload = loadmat(asymFile)
        Planet.Magnetic.nAsymBds = reload['nAsymBds'][0]
        Planet.Magnetic.iAsymBds = reload['iAsymBds'][0]
        Planet.Magnetic.zMeanAsym_km = reload['zMeanAsym_km'][0]
        Planet.Magnetic.asymShape_m = reload['asymShape_m']
        Planet.Magnetic.gravShape_m = reload['gravShape_m']
        Planet.Magnetic.asymDevs_km = reload['asymDevs_km']
        if FigMisc.nLatMap != np.size(reload['latMap_deg']):
            log.warning('Reloaded latitudes don\'t match those in FigMisc. nLatMap will be overriden.')
        Planet.latMap_deg = reload['latMap_deg'][0]
        Planet.nLatMap = np.size(Planet.latMap_deg)
        if FigMisc.nLonMap != np.size(reload['lonMap_deg']):
            log.warning('Reloaded longitudes don\'t match those in FigMisc. nLonMap will be overridden.')
        Planet.lonMap_deg = reload['lonMap_deg'][0]
        Planet.nLonMap = np.size(Planet.lonMap_deg)

        # Adjust longitudes in case we change DO_360 between calculations
        if (FigMisc.DO_360 and Planet.lonMap_deg[0] < 0) or \
           (not FigMisc.DO_360 and Planet.lonMap_deg[-1] > 180):
            lonAdj_deg = Planet.lonMap_deg + 0
            if FigMisc.DO_360:
                lonAdj_deg[lonAdj_deg == -180] = 180.1
                lonAdj_deg[lonAdj_deg < 0] = lonAdj_deg[lonAdj_deg < 0] + 360
                iSort = np.argsort(lonAdj_deg)
            else:
                lonAdj_deg[lonAdj_deg == 360] = -0.1
                lonAdj_deg[lonAdj_deg >= 180] = lonAdj_deg[lonAdj_deg >= 180] - 360
                iSort = np.argsort(lonAdj_deg)

            lonAdj_deg = lonAdj_deg[iSort]
            rDevsAdj_km = Planet.Magnetic.asymDevs_km[:, :, iSort]

            # Stitch endpoints together if we adjusted the plotting coords
            if not (np.max(lonAdj_deg) == 180 or np.max(lonAdj_deg) == 360):
                if FigMisc.DO_360:
                    lonAdj_deg = np.append(lonAdj_deg, 360)
                else:
                    lonAdj_deg = np.append(lonAdj_deg, 180)
                rDevsAdj_km = np.append(rDevsAdj_km, rDevsAdj_km[:, :, 0:1], axis=2)

            Planet.lonMap_deg = lonAdj_deg
            Planet.Magnetic.asymDevs_km = rDevsAdj_km
            Planet.thetaMap_rad = np.radians(90 - Planet.latMap_deg)
            Planet.phiMap_rad = np.radians(Planet.lonMap_deg)

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
        return -k**2 * r * (self.n+1) / (2*self.n+1) / self.n * (Q - self.n/(self.n+1))**2 - (2*self.n+1) / r * Q


def SetupInduction(Planet, Params):
    """ Reconfigure layer boundaries and conductivities into a format
        usable by magnetic induction calculation functions.
        Optionally also identify asymmetric shape information from gravity.

        Requires Planet attributes:
        Sets Planet attributes:
            Magnetic.rSigChange_m, Magnetic.sigmaLayers_Sm, Magnetic.asymShape_m
    """

    # Lots of errors can happen if we haven't calculated the electrical conductivity,
    # so we make this contingent on having it.
    if Params.CALC_CONDUCT and Planet.Do.VALID:
        # Reconfigure layer conducting boundaries as needed.
        # For inductOtype == 'sigma', we have already set these arrays.
        if not Params.Induct.inductOtype == 'sigma' or not Params.DO_INDUCTOGRAM:
            # Append optional ionosphere
            # We first check if these are unset here, then assign them to what they should be if unset
            if Planet.Magnetic.ionosBounds_m is None or Planet.Magnetic.sigmaIonosPedersen_Sm is None:
                Planet.Magnetic.ionosBounds_m = [np.nan]
                Planet.Magnetic.sigmaIonosPedersen_Sm = [np.nan]
            # Now, we handle unset ionospheres
            if np.all(np.isnan(Planet.Magnetic.ionosBounds_m)) or np.all(np.isnan(Planet.Magnetic.sigmaIonosPedersen_Sm)):
                # Make sure the arrays are both just length 1 of nan
                Planet.Magnetic.ionosBounds_m = [np.nan]
                Planet.Magnetic.sigmaIonosPedersen_Sm = [np.nan]
                zIonos_m = []
                sigmaIonos_Sm = []
            else:
                zIonos_m = Planet.Bulk.R_m + np.array(Planet.Magnetic.ionosBounds_m)
                sigmaIonos_Sm = np.array(Planet.Magnetic.sigmaIonosPedersen_Sm)
                # Allow special case for specifying an ionosphere with 1 conductivity and 2 bounds, when
                # there is a substantial neutral atmosphere and the conducting region is at altitude
                # (e.g. for Triton)
                if np.size(Planet.Magnetic.sigmaIonosPedersen_Sm) == 1 and np.size(Planet.Magnetic.ionosBounds_m) == 2:
                    sigmaIonos_Sm = np.append(0, sigmaIonos_Sm)
            # Flip arrays to be in radial ascending order as needed in induction calculations, then add ionosphere
            rLayers_m = np.append(np.flip(Planet.r_m[:-1]), zIonos_m)
            sigmaInduct_Sm = np.append(np.flip(Planet.sigma_Sm), sigmaIonos_Sm)

            # Eliminate NaN values and 0 values, assigning them to a default minimum
            sigmaInduct_Sm[np.logical_or(np.isnan(sigmaInduct_Sm), sigmaInduct_Sm == 0)] = Constants.sigmaDef_Sm
            # Set low conductivities to all be the same default value so we can shrink them down to single layers
            sigmaInduct_Sm[sigmaInduct_Sm < Constants.sigmaMin_Sm] = Constants.sigmaDef_Sm

            # Optionally, further reduce computational overhead by shrinking the number of ocean layers modeled
            if Params.Sig.REDUCED_INDUCT:
                indsLiq = np.where(np.flip(Planet.phase) == 0)[0]
                if np.size(indsLiq) > 0 and not np.all(np.diff(indsLiq) == 1):
                    log.warning('HP ices found in ocean while REDUCED_INDUCT is True. They will be ignored ' +
                                'in the interpolation.')
                if np.size(indsLiq) <= Params.Induct.nIntL:
                    log.warning(f'Only {np.size(indsLiq)} layers found in ocean, but number of layers to ' +
                                f'interpolate over is {Params.Induct.nIntL}. This profile will be not be reduced.')
                else:
                    # Get radius values from D/nIntL above the seafloor to the ice shell
                    rBot_m = Planet.Bulk.R_m - (Planet.zb_km + Planet.D_km) * 1e3
                    rTop_m = rLayers_m[indsLiq[-1]]
                    rOcean_m = np.linspace(rBot_m, rTop_m, Params.Induct.nIntL+1)[1:]
                    # Interpolate the conductivities corresponding to those radii
                    sigmaOcean_Sm = spi.interp1d(rLayers_m[indsLiq], sigmaInduct_Sm[indsLiq], kind=Params.Induct.oceanInterpMethod)(rOcean_m)
                    # Stitch together the r and σ arrays with the new ocean values
                    rLayers_m = np.concatenate((rLayers_m[:indsLiq[0]], rOcean_m, rLayers_m[indsLiq[-1]+1:]))
                    sigmaInduct_Sm = np.concatenate((sigmaInduct_Sm[:indsLiq[0]], sigmaOcean_Sm, sigmaInduct_Sm[indsLiq[-1]+1:]))

            # Get the indices of layers just below where changes happen
            iChange = [i for i,sig in enumerate(sigmaInduct_Sm) if sig != np.append(sigmaInduct_Sm, np.nan)[i+1]]
            Planet.Magnetic.sigmaLayers_Sm = sigmaInduct_Sm[iChange]
            Planet.Magnetic.rSigChange_m = rLayers_m[iChange]

        Planet.Magnetic.nBds = np.size(Planet.Magnetic.rSigChange_m)

        # Set asymmetric shape if applicable
        if Params.Sig.INCLUDE_ASYM:
            if Planet.Magnetic.pMax == 0:
                Params.Sig.INCLUDE_ASYM = False
                log.warning('Magnetic.pMax is 0, asymmetry calculations will be skipped.')
                Planet.Magnetic.asymShape_m = np.zeros((Planet.Magnetic.nBds, 2, 1, 1))
            else:
                if Planet.Magnetic.pMax is not None and (Planet.Magnetic.pMax < 2 and not SigParams.ALLOW_LOW_PMAX):
                    log.warning('SigParams.INCLUDE_ASYM is True, but Magnetic.pMax is less than 2. ' +
                                'This would fail to include gravity contributions, which are typically ' +
                                'among the largest contributors to asymmetric induction. Magnetic.pMax has been ' +
                                'increased to 2. Toggle this check with SigParams.ALLOW_LOW_PMAX in configInduct.')
                    Planet.Magnetic.pMax = 2
                Planet = SetAsymShape(Planet, Params)
                # Initialize the gravity shape array
                Planet.Magnetic.gravShape_m = np.zeros_like(Planet.Magnetic.asymShape_m)
                if Planet.Magnetic.pMax >= 2:
                    Planet = SetGravShape(Planet, Params)
                Planet.Magnetic.nAsymBds = np.size(Planet.Magnetic.zMeanAsym_km)

                # Fetch Xid array
                XidLabel = f''
                if XidLabel not in EOSlist.loaded.keys():
                    nMax = Planet.Magnetic.nprmMax + Planet.Magnetic.pMax
                    Planet.Magnetic.Xid = LoadXid(Planet.Magnetic.nprmMax, Planet.Magnetic.pMax, nMax,
                                  Planet.Magnetic.nLin, Planet.Magnetic.mLin, reload=True, do_parallel=False)
                    EOSlist.loaded[XidLabel] = Planet.Magnetic.Xid
                    EOSlist.ranges[XidLabel] = f'{Planet.Magnetic.nprmMax}x{Planet.Magnetic.pMax}x{nMax}'
                else:
                    Planet.Magnetic.Xid = EOSlist.loaded[XidLabel]
        else:
            Planet.Magnetic.pMax = 0
            Planet.Magnetic.asymShape_m = np.zeros((Planet.Magnetic.nBds, 2, 1, 1))
        
        # Get excitation spectrum
        Planet.Magnetic.Texc_hr, Planet.Magnetic.omegaExc_radps, Planet.Magnetic.Benm_nT, Planet.Magnetic.B0_nT \
            = GetBexc(Planet.name, Planet.Magnetic.SCera, Planet.Magnetic.extModel, Params.Induct.excSelectionCalc,
                      nprmMax=Planet.Magnetic.nprmMax, pMax=Planet.Magnetic.pMax)
        Planet.Magnetic.nExc = np.size(Planet.Magnetic.Texc_hr)

        # Initialize Binm array to have the same shape and data type as Benm
        if Planet.Magnetic.Binm_nT is None: 
            Planet.Magnetic.Binm_nT = np.zeros_like(Planet.Magnetic.Benm_nT)

    else:
        # Make sure explore-o-grams play nice when ionosphere properties are not set
        if Planet.Magnetic.ionosBounds_m is None or np.any(np.isnan(Planet.Magnetic.ionosBounds_m)):
            Planet.Magnetic.ionosBounds_m = [0]
        if Planet.Magnetic.sigmaIonosPedersen_Sm is None or np.any(np.isnan(Planet.Magnetic.sigmaIonosPedersen_Sm)):
            Planet.Magnetic.sigmaIonosPedersen_Sm = [0]

    return Planet, Params


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
    BeLabel = f'{bodyname}Be{nprmMax}{era}{model}{excSelection}'

    if BeLabel in EOSlist.loaded.keys():
        log.debug(f'{bodyname} excitation spectrum for {model} model and {era} era already loaded. Reusing existing.')
        Texc_hr, omegaExc_radps, Benm_nT, B0_nT = EOSlist.loaded[BeLabel]
    else:
        if bodyname[:4] == 'Test':
            fNames = [f'Be{npi}xyz_Test' for npi in range(1, nprmMax+1)]
            fPath = os.path.join(_Test, 'inductionData')
        else:
            fNames = [f'Be{npi}xyz_{bodyname}' for npi in range(1, nprmMax+1)]
            fPath = os.path.join(bodyname, 'inductionData')
        # Append era and model info
        if era is not None:
            fNames = [f'{fNamenp}_{era}' for fNamenp in fNames]
        if model is not None:
            fNames = [f'{fNamenp}_{model}' for fNamenp in fNames]
        log.debug(f'Loading {bodyname} excitation spectrum for {model} model and {era} era.')
        if nprmMax > 1:
            log.warning('n\'_max greater than 1 is not yet supported. Be only up to n=1 will be loaded.')
            nprmMax = 1

        if os.path.isfile(os.path.join(fPath, f'{fNames[0]}.txt')):
            inpTexc_hr, inpBenm_nT, B0_nT = GetBenm(nprmMax, pMax, fpath=fPath, fName=fNames[0])
            BeList = Excitations(bodyname)
            eachT = np.logical_and([excSelection[key] for key in BeList.keys()], [BeList[key] is not None for key in BeList.keys()])
            nPeaks = sum(eachT)
            Texc_hr = np.zeros(nPeaks)
            Benm_nT = np.zeros((nPeaks, 2, nprmMax+pMax+1, nprmMax+pMax+1), dtype=np.complex_)
            # Include in the excitation spectrum only the periods specified in configPPinduct.py
            iPeak = 0
            for oscillation, included in excSelection.items():
                if included and BeList[oscillation] is not None:
                    iClosest = np.argmin(abs(inpTexc_hr - BeList[oscillation]))
                    Texc_hr[iPeak] = inpTexc_hr[iClosest]
                    Benm_nT[iPeak, ...] = inpBenm_nT[inpTexc_hr == Texc_hr[iPeak], ...]
                    iPeak += 1

            omegaExc_radps = 2*np.pi / Texc_hr / 3600

            EOSlist.loaded[BeLabel] = (Texc_hr, omegaExc_radps, Benm_nT, B0_nT)
            EOSlist.ranges[BeLabel] = Texc_hr

        else:
            log.warning(f'Excitation moments file(s) not found in {fPath}. Induction calculations will be skipped.')
            Texc_hr, omegaExc_radps, Benm_nT, B0_nT = (None for _ in range(4))

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
    if os.path.isfile(Params.inducedMomentsFile):
        Planet.Magnetic.Binm_nT = np.loadtxt(Params.inducedMomentsFile, skiprows=1, unpack=False)
    else:
        log.warning(f'Params.CALC_NEW_INDUCT is True, but {Params.inducedMomentsFile} was not found. ' +
                    'Recalculating.')
        Planet, Params = MagneticInduction(Planet, Params)

    return Planet


def SetGravShape(Planet, Params):
    """ Evaluate the gravity chi_pq coefficients based on published J2 and C22 coefficients.
        Note that J2 and C22 coefficients are expected to be UNNORMALIZED spherical
        harmonic coefficients, as they are most commonly reported in the literature, e.g.
        in Anderson et al. (1998) for Europa: https://doi.org/10.1126/science.281.5385.2019 
    """
    # Check that J2 and C22 are set, and make use of other parameters if they are set
    if Planet.Bulk.J2 is None and Planet.Bulk.C20 is None:
        log.warning(f'Params.Sig.INCLUDE_ASYM is True but Bulk.J2 is not set. It will be treated as 0.')
        Planet.Bulk.J2 = 0.0
        Planet.Bulk.C20 = 0.0
    elif Planet.Bulk.C20 is None:
        Planet.Bulk.C20 = -Planet.Bulk.J2
    if Planet.Bulk.C22 is None:
        log.warning(f'Params.Sig.INCLUDE_ASYM is True but Bulk.C22 is not set. It will be treated as 0.')
        Planet.Bulk.C22 = 0.0
    if Planet.Bulk.C21 is None:
        Planet.Bulk.C21 = 0.0
    if Planet.Bulk.S21 is None:
        Planet.Bulk.S21 = 0.0
    if Planet.Bulk.S22 is None:
        Planet.Bulk.S22 = 0.0

    # Organize shape coefficients into cos and sin terms for this p.
    # p = 2 is set manually here, but we can in principle also include
    # further p = 4 gravity terms.
    p = 2
    C2q = np.array([Planet.Bulk.C20, Planet.Bulk.C21, Planet.Bulk.C22])
    S2q = np.array([            0.0, Planet.Bulk.S21, Planet.Bulk.S22])

    # Next, get the secular Love number k_f in the hydrostatic approximation. Note that k_f has no units.
    # To do this, we use equation S94 from Styczinski et al. (2021): https://doi.org/10.1016/j.icarus.2021.114840
    u = (5/2 * (1 - 3/2 * Planet.Bulk.Cmeasured))**2
    kf = (4 - u) / (1 + u)
    # The fluid Love number is just kf + 1:
    hf = kf + 1
    # Now calculate the tidal deformation terms Hpq from Eq. S93 from Styczinski et al. (2021):
    H2c_m = hf * C2q * Planet.Bulk.R_m
    H2s_m = hf * S2q * Planet.Bulk.R_m

    # These are the UNNORMALIZED deformation terms. They are incorrectly labeled as Schmidt semi-
    # normalized coefficients in Styczinski et al. (2021). To get the 4π-normalized terms we need to use
    # calculations from MoonMag from the unnormalized ones, we need to divide by the 4π-normalization factor:
    H2c_4pi_m = [H2c_m[q] / normFactor_4pi(p, q) for q in range(p+1)]
    H2s_4pi_m = [H2s_m[q] / normFactor_4pi(p, q) for q in range(p+1)]
    # Convert to fully normalized, complex coefficients with Condon-Shortley phase
    g_chipq_m = GeodesyNorm2chipq(p, H2c_4pi_m, H2s_4pi_m)
    # Construct the gravity shape array by scaling the shape to the surface radius
    radialScale = Planet.Magnetic.rSigChange_m / Planet.Bulk.R_m
    for iLayer, rScale in enumerate(radialScale):
        Planet.Magnetic.gravShape_m[iLayer, :, p, :p+1] = g_chipq_m * rScale

    Planet.Magnetic.zMeanAsym_km = np.append(Planet.Magnetic.zMeanAsym_km, 0)
    Planet.Magnetic.iAsymBds = np.append(Planet.Magnetic.iAsymBds, np.argmin(np.abs(radialScale - 1.0)))

    return Planet


def SetAsymShape(Planet, Params):
    """ Read one or more files from disk that describe asymmetric layers.
        The Bodyname/inductionData folder is searched for files starting with
        a search string (f"{Planet.name}Shape_4piNormDepth" by default), and the values are
        added to the layer with the nearest depth.
        
        If Params.Sig.CONCENTRIC_ASYM = True, only "BodynameShape_4piNormDepth.txt" is
        searched for, and this shape is scaled to all layer boundaries. 
    """
    SKIP_ASYM = False
    reason = ''
    shapeFname = f'{Planet.name}{Params.Sig.asymFstring}'
    searchPath = os.path.join(Params.DataFiles.inductPath, f'{shapeFname}*.txt')
    backupShapeFname = f'{Planet.bodyname}{Params.Sig.asymFstring}'
    backupSearchPath = os.path.join(Params.DataFiles.inductPath, f'{backupShapeFname}*.txt')
    shapeFiles = FilesMatchingPattern(searchPath)
    if Planet.name == Planet.bodyname:
        # Note this in case we need it, it's not necessarily true yet.
        reason = f'No shape file was found for {searchPath}.'
    elif np.size(shapeFiles) == 0:
        # Try backup path with just the bodyname and not the Planet.name set in the PPBody.py file
        shapeFiles = FilesMatchingPattern(backupSearchPath)
        reason = f'No shape file was found for either\n  {searchPath}\n  or\n  {backupSearchPath}. '
    if np.size(shapeFiles) == 0:
        SKIP_ASYM = True

    if Params.Sig.CONCENTRIC_ASYM:
        # Here, we load a single file and scale it to all boundaries.
        fName = os.path.join(Params.DataFiles.inductPath, f'{shapeFname}.txt')
        if fName in shapeFiles:
            log.debug(f'Using asymmetric shape file, concentrically: {fName}')
            if fName in EOSlist.loaded.keys():
                log.debug(f'Already loaded {fName}, reusing existing.')
                Planet.Magnetic.pLin, Planet.Magnetic.qLin, Cpq_km, Spq_km = EOSlist.loaded[fName]
            else:
                Planet.Magnetic.pLin, Planet.Magnetic.qLin, Cpq_km, Spq_km, _, _ \
                    = np.loadtxt(fName, skiprows=1, unpack=True, delimiter=',')
                EOSlist.loaded[fName] = (Planet.Magnetic.pLin, Planet.Magnetic.qLin, Cpq_km, Spq_km)
                EOSlist.ranges[fName] = ''

            if not (Planet.Magnetic.pLin[0] == 0 and Planet.Magnetic.qLin[0] == 0):
                raise RuntimeError(f'Shape file {fName} does not start at p,q = [0,0].')
            rAsym_m = Planet.Bulk.R_m - Cpq_km[0] * 1e3
            radialScale = Planet.Magnetic.rSigChange_m / rAsym_m
            
            # Limit Magnetic.pMax to the values we now have
            if Planet.Magnetic.pMax is None or Planet.Magnetic.pMax > np.max(Planet.Magnetic.pLin):
                Planet.Magnetic.pMax = int(np.max(Planet.Magnetic.pLin))
            # Initialize the asymShape array
            Planet.Magnetic.asymShape_m = np.zeros((Planet.Magnetic.nBds, 2, Planet.Magnetic.pMax+1, Planet.Magnetic.pMax+1),
                                                 dtype=np.complex_)
            
            # Convert 4π-normalized depth coefficients in km to fully normalized radial and in m
            for p in range(1, Planet.Magnetic.pMax+1):
                iMin = int(p * (p+1) / 2)
                iMax = iMin + p+1
                # Negate to convert from deviations of depth to radius
                Cpq_m = -Cpq_km[iMin:iMax] * 1e3
                Spq_m = -Spq_km[iMin:iMax] * 1e3
                chipq_m = GeodesyNorm2chipq(p, Cpq_m, Spq_m)
                # Scale according to radius and fill into asymShape
                for iLayer, rScale in enumerate(radialScale):
                    Planet.Magnetic.asymShape_m[iLayer, 0, 0, 0] = Planet.Magnetic.rSigChange_m[iLayer]
                    Planet.Magnetic.asymShape_m[iLayer, :, p, :p+1] = chipq_m * rScale

            Planet.Magnetic.zMeanAsym_km = np.append(Planet.Magnetic.zMeanAsym_km, Cpq_km[0])
            Planet.Magnetic.iAsymBds = np.append(Planet.Magnetic.iAsymBds, 
                                                 np.argmin(np.abs(Planet.Magnetic.rSigChange_m - rAsym_m)))
        else:
            reason = f'A shape file was not found at {fName}, but Params.Sig.CONCENTRIC_ASYM is True. '
            SKIP_ASYM = True
    else:
        # Here, we load as many files as are present and add them to the nearest-depth boundary.
        log.debug(f'Using asymmetric shape file(s):\n  ' + ',\n  '.join(shapeFiles))
        pLin, qLin, Cpq_km, Spq_km = (np.empty(np.size(shapeFiles), dtype=object) for _ in range(4))
        for i, fName in enumerate(shapeFiles):
            if fName in EOSlist.loaded.keys():
                log.debug(f'Already loaded {fName}, reusing existing.')
                pLin[i], qLin[i], Cpq_km[i], Spq_km[i] = EOSlist.loaded[fName]
            else:
                pLin[i], qLin[i], Cpq_km[i], Spq_km[i], _, _ \
                    = np.loadtxt(fName, skiprows=1, unpack=True, delimiter=',')
                EOSlist.loaded[fName] = (pLin[i], qLin[i], Cpq_km[i], Spq_km[i])
                EOSlist.ranges[fName] = ''
        
        # Limit Magnetic.pMax to the values we now have
        pMaxNonzero = int(np.max([np.max(pL[np.logical_or(Cpq_km[i] != 0, Spq_km[i] != 0)]) for i, pL in enumerate(pLin)]))
        if Planet.Magnetic.pMax is None or Planet.Magnetic.pMax > pMaxNonZero:
            Planet.Magnetic.pMax = pMaxNonzero
        # Initialize the shape array
        Planet.Magnetic.asymShape_m = np.zeros((Planet.Magnetic.nBds, 2, Planet.Magnetic.pMax+1, Planet.Magnetic.pMax+1),
                                                 dtype=np.complex_)

        for i, fName in enumerate(shapeFiles):
            pMaxNonzero = np.max(pLin[i][np.logical_or(Cpq_km[i] != 0, Spq_km[i] != 0)])
            pMax = int(np.minimum(pMaxNonzero, Planet.Magnetic.pMax))
            if not (pLin[i][0] == 0 and qLin[i][0] == 0):
                raise RuntimeError(f'Shape file {fName} does not start at p,q = [0,0].')
            rAsym_m = Planet.Bulk.R_m - Cpq_km[i][0] * 1e3
            radDiff = Planet.Magnetic.rSigChange_m - rAsym_m
            iLayer = np.argmin(abs(radDiff))
            # Convert 4π-normalized depth coefficients in km to fully normalized radial and in m
            Planet.Magnetic.asymShape_m[iLayer, 0, 0, 0] = Planet.Magnetic.rSigChange_m[iLayer]
            for p in range(1, pMax+1):
                iMin = int(p * (p+1) / 2)
                iMax = iMin + p+1
                # Negate to convert from deviations of depth to radius
                Cpq_m = -Cpq_km[i][iMin:iMax] * 1e3
                Spq_m = -Spq_km[i][iMin:iMax] * 1e3
                chipq_m = GeodesyNorm2chipq(p, Cpq_m, Spq_m)
                # Scale according to radius and fill into asymShape
                Planet.Magnetic.asymShape_m[iLayer, :, p, :p+1] = Planet.Magnetic.asymShape_m[iLayer, :, p, :p+1] \
                        + chipq_m * rAsym_m / Planet.Magnetic.rSigChange_m[iLayer]

                Planet.Magnetic.zMeanAsym_km = np.append(Planet.Magnetic.zMeanAsym_km, np.real(Planet.Bulk.R_m - Planet.Magnetic.asymShape_m[iLayer, 0, 0, 0])/1e3)
                Planet.Magnetic.iAsymBds = np.append(Planet.Magnetic.iAsymBds, iLayer)

    if SKIP_ASYM:
        log.warning(reason + ' Asymmetry will be modeled only for the gravity coefficients.')
        Planet.Magnetic.pMax = 2
        Planet.Magnetic.asymShape_m = np.zeros((Planet.Magnetic.nBds, 2, Planet.Magnetic.pMax+1, Planet.Magnetic.pMax+1),
                                             dtype=np.complex_)
    else:
        Planet.Magnetic.zMeanAsym_km = np.unique(Planet.Magnetic.zMeanAsym_km)
        Planet.Magnetic.iAsymBds = np.unique(Planet.Magnetic.iAsymBds)

    if Planet.Magnetic.pLin is None:
        Planet.Magnetic.pLin = [p for p in range(1, Planet.Magnetic.pMax+1) for _ in range(-p, p+1)]
        Planet.Magnetic.qLin = [q for p in range(1, Planet.Magnetic.pMax+1) for q in range(-p, p+1)]

    return Planet


def normFactor_4pi(n, m):
    """ Calculate the normalization factor for 4π-normalized spherical harmonics,
        without the Condon-Shortley phase, as needed for shape calculations that make
        use of infrastructure from MoonMag.
    """
    return np.sqrt((2*n+1) * sps.factorial(n - abs(m)) / sps.factorial(n + abs(m)))


def FourierSpectrum(Planet, Params):
    """ Load a Fourier spectrum of magnetic excitations applied to the body,
        calculate the spherically symmetric induction amplitude across this
        spectrum, and calculate the surface field amplitude for each component.
    """

    # Load Fourier spectrum data from disk
    Planet = LoadFTdata(Planet, Params)
    
    if Planet.Magnetic.FT_LOADED:
        if Params.CALC_NEW_INDUCT:
            log.debug('Calculating magnetic Fourier spectrum.')
            
            # Define a smaller subset of periods over which to evaluate complex response
            # amplitudes to save on computation time (behavior will be smooth, we'll interpolate)
            TexcReduced_hr = np.geomspace(np.min(Planet.Magnetic.TexcFT_hr), np.max(Planet.Magnetic.TexcFT_hr),
                                          Params.MagSpectrum.nOmegaPts)
            omegaReduced_radps = 2 * np.pi / TexcReduced_hr / 3600
            # Evaluate complex amplitudes
            Ae1FTreduced, _, _ \
                = AeList(Planet.Magnetic.rSigChange_m, Planet.Magnetic.sigmaLayers_Sm,
                         omegaReduced_radps, 1/Planet.Bulk.R_m, nn=1,
                         writeout=False, do_parallel=Params.DO_PARALLEL)

            # Interpolate to full excitation spectrum
            Planet.Magnetic.Ae1FT = spi.interp1d(TexcReduced_hr, Ae1FTreduced, kind=Params.MagSpectrum.interpMethod
                                                 )(Planet.Magnetic.TexcFT_hr)

            # Get complex induced dipole components
            Planet.Magnetic.Bi1xyzFT_nT = {vComp: Planet.Magnetic.Be1xyzFT_nT[vComp] * Planet.Magnetic.Ae1FT
                                           for vComp in ['x', 'y', 'z']}
            
            if not Params.NO_SAVEFILE:
                # Save to disk
                saveDict = {
                    'Ae1FT': Planet.Magnetic.Ae1FT,
                    'Bi1x_nT': Planet.Magnetic.Bi1xyzFT_nT['x'],
                    'Bi1y_nT': Planet.Magnetic.Bi1xyzFT_nT['y'],
                    'Bi1z_nT': Planet.Magnetic.Bi1xyzFT_nT['z']
                }
                savemat(Params.DataFiles.FTdata, saveDict)
                log.debug(f'Saved magnetic Fourier spectrum to file: {Params.DataFiles.FTdata}')
        else:
            ReloadSpectrum(Planet, Params)

    return Planet


def LoadFTdata(Planet, Params):
    """ Load Fourier spectrum data for this body """

    # Get path to excitation spectrum file
    if Planet.bodyname == 'Test':
        BeFTdataPath = os.path.join(_Test, Params.DataFiles.BeFTdata)
    else:
        BeFTdataPath = os.path.join(_Defaults, Planet.bodyname, Params.DataFiles.BeFTdata)

    if os.path.isfile(BeFTdataPath):
        BeFTdata = loadmat(BeFTdataPath)
        Planet.Magnetic.TexcFT_hr = BeFTdata['T_h'][0]
        Planet.Magnetic.TmaxFT_hr = BeFTdata['Tmax'][0,0]
        Planet.Magnetic.extModelFT = BeFTdata['magModelDescrip'][0]
        Planet.Magnetic.coordTypeFT = BeFTdata['coordType'][0]
        Planet.Magnetic.Be1xyzFT_nT = {
            'x': BeFTdata['B1x_nT'][0],
            'y': BeFTdata['B1y_nT'][0],
            'z': BeFTdata['B1z_nT'][0],
        }
        log.debug(f'Loaded magnetic excitation spectrum from file: {BeFTdataPath}')
        Planet.Magnetic.FT_LOADED = True
        
    else:
        log.warning(f'PLOT_MAG_SPECTRUM is True but {BeFTdataPath} was not found. Skipping.')
        Planet.Magnetic.FT_LOADED = False

    return Planet


def ReloadSpectrum(Planet, Params):
    """ Reload Fourier spectrum calculations from disk """

    if os.path.isfile(Params.DataFiles.FTdata):
        FTdata = loadmat(Params.DataFiles.FTdata)
        Planet.Magnetic.Ae1FT = FTdata['Ae1FT'][0]
        Planet.Magnetic.Bi1xyzFT_nT = {
            'x': FTdata['Bi1x_nT'][0],
            'y': FTdata['Bi1y_nT'][0],
            'z': FTdata['Bi1z_nT'][0],
        }
        log.debug(f'Loaded magnetic Fourier spectrum from file: {Params.DataFiles.FTdata}')
        Planet.Magnetic.FT_LOADED = True

    else:
        log.warning(f'Magnetic Fourier spectrum was intended to be reloaded, but {Params.DataFiles.FTdata} ' +
                    f'was not found. Re-run with CALC_NEW_INDUCT = True in configPP.py.')
        Planet.Magnetic.FT_LOADED = False

    return Planet


def FindBindCA(Planet, Params, scName, tStrList):
    """ Get induced magnetic field vector at a list of closest approach times
        for the named spacecraft.
    """
    etCAapprox = spice.str2et(tStrList)
    etCA = np.zeros_like(etCAapprox)
    for i, etCAnamed in enumerate(etCAapprox):
        etSearch = np.arange(etCAnamed-Params.tRangeCA_s, etCAnamed+Params.tRangeCA_s, Params.tSearchRes_s)
        _, _, _, rRange_km = BodyDist_km(scName, Planet.bodyname, etSearch)
        etCA[i] = etSearch[np.argmin(rRange_km)]

    return Planet
