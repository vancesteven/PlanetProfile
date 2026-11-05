import os
import numpy as np
import logging
from scipy.io import loadmat
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator, interp1d as Interp1D, griddata as GridData
from PlanetProfile import _ROOT, _PERPLEXCACHE
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist, ConstantPropsStruct
from PlanetProfile.Utilities.DataManip import (ResetNearestExtrap, ReturnZeros, EOSwrapper,
                                              EnforceEOSDomain, EOSDomainExceedance, EOSDomainError)
from PlanetProfile.Thermodynamics.ConstantEOS import GetConstantEOS

# Assign logger
log = logging.getLogger('PlanetProfile')

# Introspection bypass for EOS-domain enforcement.
#
# A few places probe an EOS with sentinel arguments (P=0, T=0) purely to
# discover its TYPE or its zero-porosity value -- notably GetphiCalc, which
# calls fn_phiEOS_frac(0, 0) to distinguish a real porosity EOS from a
# ReturnZeros stand-in. Those sentinels are outside every real table domain by
# construction and are not scientific evaluations, so they must not trip a
# strict EOS-domain policy.
#
# This bypass is deliberately narrow: it is enabled only inside the
# EOSIntrospection() context manager, never for a property evaluation on a
# physical (P, T) node.
_EOS_INTROSPECTION_DEPTH = 0


class EOSIntrospection:
    """Temporarily suspend EOS-domain enforcement for sentinel/type probes."""

    def __enter__(self):
        global _EOS_INTROSPECTION_DEPTH
        _EOS_INTROSPECTION_DEPTH += 1
        return self

    def __exit__(self, *exc):
        global _EOS_INTROSPECTION_DEPTH
        _EOS_INTROSPECTION_DEPTH -= 1
        return False


def EOSIntrospectionActive():
    return _EOS_INTROSPECTION_DEPTH > 0


# Candidate-propagation phase.
#
# PlanetProfile propagates every candidate silicate profile from its start
# radius down toward r = 0 as scratch storage for the core-radius search; only
# indices [0, nSilFinal[i]) of candidate i are ever part of a physical
# structure (see IronCore.IronCoreLayers and the Csil_kgm2 sum in
# FindInnerWithMoIAndEOS). Those scratch nodes routinely exceed the EOS
# pressure ceiling and any finite prescribed thermal range, so enforcing a
# strict policy during propagation aborts the search on storage that is never
# used scientifically.
#
# Inside EOSCandidatePropagation(), a strict 'error' policy is therefore
# downgraded to the historical clamp. Because the recursion is strictly
# downward (node j depends only on j-1) and active nodes are always shallower
# than scratch nodes, clamped scratch values can never propagate back into a
# physical node. Strictness is then enforced on the physical nodes themselves
# by ValidateCandidateEOSDomains once each candidate's extent is known.
#
# 'warn' is left untouched so the diagnostics remain visible.
_EOS_CANDIDATE_PROPAGATION_DEPTH = 0


class EOSCandidatePropagation:
    """Downgrade strict EOS-domain failure to clamping during candidate search."""

    def __enter__(self):
        global _EOS_CANDIDATE_PROPAGATION_DEPTH
        _EOS_CANDIDATE_PROPAGATION_DEPTH += 1
        return self

    def __exit__(self, *exc):
        global _EOS_CANDIDATE_PROPAGATION_DEPTH
        _EOS_CANDIDATE_PROPAGATION_DEPTH -= 1
        return False


def EOSCandidatePropagationActive():
    return _EOS_CANDIDATE_PROPAGATION_DEPTH > 0


def ValidateCandidateEOSDomains(Planet, Psil_MPa, Tsil_K, nSilFinal, indsSilValid,
                                Pcore_MPa=None, Tcore_K=None, policy=None):
    """Enforce the EOS-domain policy on PHYSICALLY ACTIVE nodes only.

    For candidate i the physically active silicate extent is indices
    [0, nSilFinal[i]) -- the same extent used for that candidate's mass and
    moment-of-inertia sums, so no second notion of "active node" is introduced
    here. Anything at or beyond nSilFinal[i] is scratch storage replaced by
    core.

    Within that extent, only nodes whose MATERIAL is mantle are evaluated by
    the silicate EOS: an explicit crust supplies prescribed constants instead,
    and its shallow nodes legitimately sit below the mantle table's P and T
    floors. The existing Crust.material_mask is reused to make that
    distinction rather than re-deriving it. All core nodes are physical.

    Returns a diagnostics dict. Under policy 'error', raises EOSDomainError if
    any physically active node of any viable candidate lies outside its EOS
    domain. Under 'warn' it logs; under 'legacy_clamp' it only reports.
    """
    if policy is None:
        policy = ResolveEOSdomainPolicy(Planet, type('P', (), {})())
    silEOS = EOSlist.loaded.get(getattr(Planet.Sil.EOS, 'key', None), None)
    coreEOS = EOSlist.loaded.get(getattr(Planet.Core.EOS, 'key', None), None)

    report = {'nCandidates': 0, 'nBadCandidates': 0, 'nBadNodes': 0, 'offenders': []}

    def _scan(P, T, eos, label, idx):
        if eos is None or np.size(P) == 0:
            return 0
        rep = EOSDomainExceedance(P, T, eos.Pmin, eos.Pmax, eos.Tmin, eos.Tmax)
        if rep['nOutside']:
            report['offenders'].append(
                f'candidate {idx} {label}: {rep["nOutside"]}/{rep["nTotal"]} active node(s) outside '
                f'(P {rep["min1Requested"]:.2f}-{rep["max1Requested"]:.2f} MPa vs [{eos.Pmin:.0f},{eos.Pmax:.0f}]; '
                f'T {rep["min2Requested"]:.2f}-{rep["max2Requested"]:.2f} K vs [{eos.Tmin:.0f},{eos.Tmax:.0f}])')
        return rep['nOutside']

    nSilFinal = np.atleast_1d(nSilFinal)
    matMask = getattr(Planet.Crust, 'material_mask', None) if Planet.Crust.enabled else None
    for i in np.atleast_1d(indsSilValid):
        i = int(i)
        if i >= np.size(nSilFinal):
            continue
        nAct = int(nSilFinal[i])
        if nAct <= 0:
            continue
        report['nCandidates'] += 1
        Pa = np.asarray(Psil_MPa)[i, :nAct]
        Ta = np.asarray(Tsil_K)[i, :nAct]
        if matMask is not None and np.asarray(matMask).ndim == 2 and i < np.shape(matMask)[0]:
            isMantle = np.asarray(matMask)[i, :nAct] == Constants.materialMantle
            Pa, Ta = Pa[isMantle], Ta[isMantle]
            report['nCrustExcluded'] = report.get('nCrustExcluded', 0) + int(np.count_nonzero(~isMantle))
        bad = _scan(Pa, Ta, silEOS, 'mantle', i)
        report['nBadNodes'] += bad
        if bad:
            report['nBadCandidates'] += 1

    if Pcore_MPa is not None and coreEOS is not None and np.size(Pcore_MPa):
        Pc = np.asarray(Pcore_MPa); Tc = np.asarray(Tcore_K)
        bad = _scan(Pc.ravel(), Tc.ravel(), coreEOS, 'core', -1)
        report['nBadNodes'] += bad

    if report['nBadNodes']:
        msg = ('EOS domain violated on PHYSICALLY ACTIVE nodes (scratch storage excluded):\n  '
               + '\n  '.join(report['offenders'][:10]))
        if len(report['offenders']) > 10:
            msg += f"\n  ... and {len(report['offenders'])-10} more"
        if policy == 'error':
            raise EOSDomainError(msg)
        elif policy == 'warn':
            log.warning(msg)
    return report


def ResolveEOSdomainPolicy(Planet, Params):
    """Per-body EOS-domain policy, falling back to the global default.

    A body configuration may demand a stricter policy than the run-wide
    default (e.g. a lunar model that must fail rather than silently saturate
    at the edge of its Perple_X table). Planet.Do.EOS_DOMAIN_POLICY = None
    means inherit Params.EOS_DOMAIN_POLICY.
    """
    bodyPolicy = getattr(getattr(Planet, 'Do', None), 'EOS_DOMAIN_POLICY', None)
    if bodyPolicy is not None:
        return bodyPolicy
    return getattr(Params, 'EOS_DOMAIN_POLICY', 'legacy_clamp')


def GetPlanetInnerEOS(Planet, Params, Fe_EOS = False):
    if Fe_EOS:
        if Planet.Do.ConstantProps['Inner']:
            innerEOSWrapper = GetConstantEOS(Planet.Core.ConstantProps, EOStype='inner', EOScomp='Fe')
        else:
            innerEOSWrapper = GetInnerEOS(Planet.Core.coreEOS, EOSinterpMethod=Params.lookupInterpMethod, Fe_EOS=True,
                                      kThermConst_WmK=Planet.Core.kTherm_WmK, EXTRAP=Params.EXTRAP_Fe,
                                      wFeCore_ppt=Planet.Core.wFe_ppt, wScore_ppt=Planet.Core.wS_ppt, etaSilFixed_Pas=Planet.Sil.etaRock_Pas, etaCoreFixed_Pas=[Planet.Core.etaFeSolid_Pas, Planet.Core.etaFeLiquid_Pas],
                                      TviscTrans_K=Planet.Core.TviscTrans_K,
                                      EOSdomainPolicy=ResolveEOSdomainPolicy(Planet, Params))
    else:
        if Planet.Do.ConstantProps['Inner']:
            innerEOSWrapper = GetConstantEOS(Planet.Sil.ConstantProps, EOStype='inner', EOScomp='Sil')
        else:
            innerEOSWrapper = GetInnerEOS(Planet.Sil.mantleEOS, EOSinterpMethod=Params.lookupInterpMethod,
                        kThermConst_WmK=Planet.Sil.kTherm_WmK, HtidalConst_Wm3=Planet.Sil.Htidal_Wm3,
                        porosType=Planet.Sil.porosType, phiTop_frac=Planet.Sil.phiRockMax_frac,
                        Pclosure_MPa=Planet.Sil.Pclosure_MPa, phiMin_frac=Planet.Sil.phiMin_frac,
                        EXTRAP=Params.EXTRAP_SIL, etaSilFixed_Pas=Planet.Sil.etaRock_Pas,
                        etaCoreFixed_Pas=[Planet.Core.etaFeSolid_Pas, Planet.Core.etaFeLiquid_Pas],
                        TviscTrans_K=Planet.Sil.TviscTrans_K,
                        EOSdomainPolicy=ResolveEOSdomainPolicy(Planet, Params))
    return innerEOSWrapper


def GetInnerEOS(EOSfname, EOSinterpMethod='nearest', nHeaders=13, Fe_EOS=False, kThermConst_WmK=None,
                HtidalConst_Wm3=0, porosType=None, phiTop_frac=0, Pclosure_MPa=350, phiMin_frac=None,
                EXTRAP=False, wFeCore_ppt=None, wScore_ppt=None, etaSilFixed_Pas=None, etaCoreFixed_Pas=None, TviscTrans_K=None,
                EOSdomainPolicy='legacy_clamp'):
    innerEOS = PerplexEOSStruct(EOSfname, EOSinterpMethod=EOSinterpMethod, nHeaders=nHeaders,
                                Fe_EOS=Fe_EOS, kThermConst_WmK=kThermConst_WmK,
                                HtidalConst_Wm3=HtidalConst_Wm3, porosType=porosType,
                                phiTop_frac=phiTop_frac, Pclosure_MPa=Pclosure_MPa,
                                phiMin_frac=phiMin_frac, EXTRAP=EXTRAP,
                                wFeCore_ppt=wFeCore_ppt, wScore_ppt=wScore_ppt,
                                etaSilFixed_Pas=etaSilFixed_Pas, etaCoreFixed_Pas=etaCoreFixed_Pas, TviscTrans_K=TviscTrans_K,
                                EOSdomainPolicy=EOSdomainPolicy)
    if innerEOS.ALREADY_LOADED:
        log.debug(f'{innerEOS.comp} EOS already loaded. Reusing existing EOS.')
        innerEOS = EOSlist.loaded[innerEOS.EOSlabel]
        # The domain policy is a runtime behavior, not a property of the
        # tabulated data, so it is deliberately excluded from EOSlabel. Re-apply
        # it here so a cached EOS never carries a stale policy.
        innerEOS.EOSdomainPolicy = EOSdomainPolicy

    # Ensure each EOSlabel is included in EOSlist, in case we have reused EOSs with
    # e.g. a smaller range that can reuse the larger-range already-loaded EOS.
    if innerEOS.EOSlabel not in EOSlist.loaded.keys():
        EOSlist.loaded[innerEOS.EOSlabel] = innerEOS

    innerEOSwrapper = EOSwrapper(innerEOS.EOSlabel)

    return innerEOSwrapper

class PerplexEOSStruct:
    """ Loads in Perple_X table and creates interpolators that can be called to
        obtain silicate/core properties as functions of P and T.
    """
    def __init__(self, EOSfname, EOSinterpMethod='nearest', nHeaders=13, Fe_EOS=False, kThermConst_WmK=None,
                 HtidalConst_Wm3=0, porosType=None, phiTop_frac=0, Pclosure_MPa=350, phiMin_frac=None,
                 EXTRAP=False, wFeCore_ppt=None, wScore_ppt=None, etaSilFixed_Pas=None, etaCoreFixed_Pas=None, TviscTrans_K=None, EOSdomainPolicy='legacy_clamp'):
        self.comp = EOSfname[:-4]
        # EOS-domain policy: 'legacy_clamp' (historical, silent), 'warn', or 'error'.
        # Stored on the instance so it is available to every fn_* property call.
        self.EOSdomainPolicy = EOSdomainPolicy
        self.EOSlabel = f'comp{self.comp}interp{EOSinterpMethod}kTherm{kThermConst_WmK}' + \
                f'Htidal{HtidalConst_Wm3}poros{porosType}phi{phiTop_frac}' + \
                f'Pclose{Pclosure_MPa}phiMin{phiMin_frac}extrap{EXTRAP}' + \
                f'wFeppt{wFeCore_ppt}wSppt{wScore_ppt}etaFixed{etaSilFixed_Pas}' + \
                f'TviscTrans{TviscTrans_K}'
        if self.EOSlabel in EOSlist.loaded.keys():
            self.ALREADY_LOADED = True
            # A direct constructor call on a cached table used to return an object
            # with none of the loaded state (no domain limits, no fn_*). Copy the
            # cached instance's state so the object is usable on its own; the
            # domain policy of this call is kept, since it is runtime behavior.
            cached = EOSlist.loaded[self.EOSlabel]
            if cached is not self:
                self.__dict__.update({k: v for k, v in vars(cached).items()
                                      if k not in ('EOSdomainPolicy', 'ALREADY_LOADED')})
        else:
            self.ALREADY_LOADED = False

        if not self.ALREADY_LOADED:
            if Fe_EOS:
                descrip = ' core composition.'
            else:
                descrip = f' silicate composition with Htidal_Wm3 = {HtidalConst_Wm3}, porosType ' + \
                          f'{porosType}, phiTop_frac = {phiTop_frac}, Pclosure_MPa = {Pclosure_MPa}, ' + \
                          f'and EXTRAP = {EXTRAP}.'
            log.debug(f'Loading Perplex EOS for {self.comp}{descrip}')
            cache_fpath = os.path.join(_PERPLEXCACHE, EOSfname)
            self.dir = _PERPLEXCACHE
            self.fpath = cache_fpath
            self.Fe_EOS = Fe_EOS
            self.EXTRAP = EXTRAP
            self.EOStype = 'inner'
            self.porosType = porosType
            self.phiTop_frac = phiTop_frac
            self.Pclosure_MPa = Pclosure_MPa
            self.phiMin_frac = phiMin_frac
            self.wFeCore_ppt = wFeCore_ppt
            self.wScore_ppt = wScore_ppt

            if '3D_EOS' in self.fpath:
                if wFeCore_ppt is None and wScore_ppt is None:
                    log.warning(f'3D EOS {self.fpath} is being loaded, but wFeCore_ppt and wScore_ppt ' +
                                f'are both None -- at least one must be specified to select the table. ' +
                                f'Default value of {Constants.wFeDef_ppt} will be used.')
                    wFeCore_ppt = Constants.wFeDef_ppt
                elif wFeCore_ppt is not None and wScore_ppt is not None:
                    log.warning('Both wFeCore_ppt and wScore_ppt were specified. wScore_ppt will be ignored.')
                if wFeCore_ppt is None:
                    wFeCore_ppt = 1e3 - wScore_ppt
                tableKey = f'{self.fpath}_wFe{wFeCore_ppt}'
            else:
                tableKey = self.fpath

            if tableKey in EOSlist.loaded.keys():
                log.debug('Reloading previously imported Perplex table.')
                P_MPa, T_K, self.ufn_rho_kgm3, self.ufn_VP_kms, self.ufn_VS_kms, \
                self.ufn_KS_GPa, self.ufn_GS_GPa, self.ufn_Cp_JkgK, self.ufn_alpha_pK \
                    = EOSlist.loaded[tableKey]
                self.Pmin, self.Pmax, self.Tmin, self.Tmax, \
                self.deltaP, self.deltaT, self.EOSdeltaP, self.EOSdeltaT \
                    = EOSlist.ranges[tableKey]
            else:
                if '3D_EOS' in self.fpath:
                    EOS3D = loadmat(self.fpath)
                    wFe_ppt = EOS3D['wFe_ppt'][0]
                    P1D_MPa = EOS3D['P_MPa'][0]
                    T1D_K = EOS3D['T_K'][0]
                    intPts = (wFe_ppt, P1D_MPa, T1D_K)
                    Peval_MPa, Teval_K = np.meshgrid(P1D_MPa, T1D_K, indexing='ij')
                    Peval_MPa = Peval_MPa.flatten()
                    Teval_K = Teval_K.flatten()
                    evalPts = np.array([wFeCore_ppt * np.ones_like(Peval_MPa), Peval_MPa, Teval_K]).T
                    rho_kgm3 = RegularGridInterpolator(intPts, EOS3D['rho_kgm3'], method='linear')(evalPts)
                    VP_kms = RegularGridInterpolator(intPts, EOS3D['VP_kms'], method='linear')(evalPts)
                    VS_kms = RegularGridInterpolator(intPts, EOS3D['VS_kms'], method='linear')(evalPts)
                    Cp_JkgK = RegularGridInterpolator(intPts, EOS3D['Cp_JkgK'], method='linear')(evalPts)
                    alpha_pK = RegularGridInterpolator(intPts, EOS3D['alpha_pK'], method='linear')(evalPts)
                    KS_GPa = RegularGridInterpolator(intPts, EOS3D['KS_GPa'], method='linear')(evalPts)
                    GS_GPa = RegularGridInterpolator(intPts, EOS3D['GS_GPa'], method='linear')(evalPts)

                    tableShape = (np.size(P1D_MPa), -1)
                    rho_kgm3 = np.reshape(rho_kgm3, tableShape)
                    VP_kms = np.reshape(VP_kms, tableShape)
                    VS_kms = np.reshape(VS_kms, tableShape)
                    Cp_JkgK = np.reshape(Cp_JkgK, tableShape)
                    alpha_pK = np.reshape(alpha_pK, tableShape)
                    KS_GPa = np.reshape(KS_GPa, tableShape)
                    GS_GPa = np.reshape(GS_GPa, tableShape)

                    self.Pmin = np.min(P1D_MPa)
                    self.Pmax = np.max(P1D_MPa)
                    self.Tmin = np.min(T1D_K)
                    self.Tmax = np.max(T1D_K)
                    self.deltaP = P1D_MPa[1] - P1D_MPa[0]
                    self.deltaT = T1D_K[1] - T1D_K[0]
                    self.EOSdeltaP = self.deltaP
                    self.EOSdeltaT = self.deltaT
                else:
                    # Load in Perple_X data. Note that all P, KS, and GS are stored as bar
                    firstPT, secondPT, rho_kgm3, VP_kms, VS_kms, Cp_Jm3K, alpha_pK, KS_bar, GS_bar \
                        = np.loadtxt(self.fpath, skiprows=nHeaders, unpack=True)
                    # We don't know yet whether P or T is in the first column, which increments the fastest.
                    # The second column increments once for each time the first column runs through the whole range.
                    dim2 = next(i[0] for i, val in np.ndenumerate(secondPT) if val != secondPT[0])
                    # Check the column header to see if P or T is printed first
                    with open(self.fpath) as f:
                        [f.readline() for _ in range(nHeaders-1)]
                        colHeaderLine = f.readline().strip()
                    # Get necessary values to generate P, T arrays
                    if colHeaderLine[0] == 'P':
                        P_FIRST = True
                        Plin_MPa = firstPT * Constants.bar2MPa  # Pressures saved as bar
                        Tlin_K = secondPT
                        lenP = dim2
                    elif colHeaderLine[0] == 'T':
                        P_FIRST = False
                        Plin_MPa = secondPT * Constants.bar2MPa
                        Tlin_K = firstPT
                        lenP = int(len(Plin_MPa)/dim2)
                    else:
                        raise ValueError(f'Perple_X table {EOSfname} does not have T or P in the first column.')

                    # Make 1D arrays of P and T that just span the axes (unlike the 2D meshes of the dependent variables)
                    self.Pmin = Plin_MPa[0]
                    self.Pmax = Plin_MPa[-1]
                    self.Tmin = Tlin_K[0]
                    self.Tmax = Tlin_K[-1]
                    lenT = int(len(Tlin_K) / lenP)
                    P1D_MPa, self.deltaP = np.linspace(self.Pmin, self.Pmax, lenP, retstep=True)
                    T1D_K, self.deltaT = np.linspace(self.Tmin, self.Tmax, lenT, retstep=True)
                    self.EOSdeltaP = self.deltaP
                    self.EOSdeltaT = self.deltaT

                    # Set unphysical values to NaN so they will be caught by the next step (for all but alpha, which can cross zero)
                    rho_kgm3[rho_kgm3 <= 0] = np.nan
                    VP_kms[VP_kms <= 0] = np.nan
                    VS_kms[VS_kms <= 0] = np.nan
                    Cp_Jm3K[Cp_Jm3K <= 0] = np.nan
                    KS_bar[KS_bar <= 0] = np.nan
                    GS_bar[GS_bar <= 0] = np.nan

                    # Interpolate dependent variables where the values are NaN
                    PTpts = (Plin_MPa, Tlin_K)
                    thisVarValid = np.isfinite(rho_kgm3)
                    if not np.all(thisVarValid):
                        thisValidPs = Plin_MPa[np.where(thisVarValid)]
                        thisValidTs = Tlin_K[np.where(thisVarValid)]
                        rho_kgm3 = GridData((thisValidPs, thisValidTs), rho_kgm3[thisVarValid], PTpts, method=EOSinterpMethod)
                    thisVarValid = np.isfinite(VP_kms)
                    if not np.all(thisVarValid):
                        thisValidPs = Plin_MPa[np.where(thisVarValid)]
                        thisValidTs = Tlin_K[np.where(thisVarValid)]
                        VP_kms = GridData((thisValidPs, thisValidTs), VP_kms[thisVarValid], PTpts, method=EOSinterpMethod)
                    thisVarValid = np.isfinite(VS_kms)
                    if not np.all(thisVarValid):
                        thisValidPs = Plin_MPa[np.where(thisVarValid)]
                        thisValidTs = Tlin_K[np.where(thisVarValid)]
                        VS_kms = GridData((thisValidPs, thisValidTs), VS_kms[thisVarValid], PTpts, method=EOSinterpMethod)
                    thisVarValid = np.isfinite(Cp_Jm3K)
                    if not np.all(thisVarValid):
                        thisValidPs = Plin_MPa[np.where(thisVarValid)]
                        thisValidTs = Tlin_K[np.where(thisVarValid)]
                        # Note that Perple_X tables store heat capacities as J/m^3/K, so we convert to J/kg/K by dividing by rho in a moment
                        Cp_Jm3K = GridData((thisValidPs, thisValidTs), Cp_Jm3K[thisVarValid], PTpts, method=EOSinterpMethod)
                    thisVarValid = np.isfinite(alpha_pK)
                    if not np.all(thisVarValid):
                        thisValidPs = Plin_MPa[np.where(thisVarValid)]
                        thisValidTs = Tlin_K[np.where(thisVarValid)]
                        alpha_pK = GridData((thisValidPs, thisValidTs), alpha_pK[thisVarValid], PTpts, method=EOSinterpMethod)
                    thisVarValid = np.isfinite(KS_bar)
                    if not np.all(thisVarValid):
                        thisValidPs = Plin_MPa[np.where(thisVarValid)]
                        thisValidTs = Tlin_K[np.where(thisVarValid)]
                        KS_bar = GridData((thisValidPs, thisValidTs), KS_bar[thisVarValid], PTpts, method=EOSinterpMethod)
                    thisVarValid = np.isfinite(GS_bar)
                    if not np.all(thisVarValid):
                        thisValidPs = Plin_MPa[np.where(thisVarValid)]
                        thisValidTs = Tlin_K[np.where(thisVarValid)]
                        GS_bar = GridData((thisValidPs, thisValidTs), GS_bar[thisVarValid], PTpts, method=EOSinterpMethod)

                    # Check that NaN removal worked correctly
                    errNaNstart = 'Failed to interpolate over NaNs in PerplexEOS '
                    errNaNend = ' values. The NaN gap may be too large to use this Perple_X output.'
                    if np.any(np.isnan(rho_kgm3)): raise RuntimeError(errNaNstart + 'rho' + errNaNend)
                    if np.any(np.isnan(VP_kms)): raise RuntimeError(errNaNstart + 'VP' + errNaNend)
                    if np.any(np.isnan(VS_kms)): raise RuntimeError(errNaNstart + 'VS' + errNaNend)
                    if np.any(np.isnan(Cp_Jm3K)): raise RuntimeError(errNaNstart + 'Cp' + errNaNend)
                    if np.any(np.isnan(alpha_pK)): raise RuntimeError(errNaNstart + 'alpha' + errNaNend)
                    if np.any(np.isnan(KS_bar)): raise RuntimeError(errNaNstart + 'KS' + errNaNend)
                    if np.any(np.isnan(GS_bar)): raise RuntimeError(errNaNstart + 'GS' + errNaNend)

                    # Now make 2D grids of values.
                    rho_kgm3 = np.reshape(rho_kgm3, (-1,dim2))
                    VP_kms = np.reshape(VP_kms, (-1,dim2))
                    VS_kms = np.reshape(VS_kms, (-1,dim2))
                    Cp_JkgK = np.reshape(Cp_Jm3K, (-1,dim2)) / rho_kgm3
                    alpha_pK = np.reshape(alpha_pK, (-1,dim2))
                    KS_GPa = np.reshape(KS_bar, (-1,dim2)) * Constants.bar2GPa
                    GS_GPa = np.reshape(GS_bar, (-1,dim2)) * Constants.bar2GPa

                    if P_FIRST:
                        # Transpose 2D meshes if P is the first column.
                        rho_kgm3 = rho_kgm3.T
                        VP_kms = VP_kms.T
                        VS_kms = VS_kms.T
                        Cp_JkgK = Cp_JkgK.T
                        alpha_pK = alpha_pK.T
                        KS_GPa = KS_GPa.T
                        GS_GPa = GS_GPa.T

                P_MPa = P1D_MPa
                T_K = T1D_K

                # Assign temporary functions we will wrap with porosity if modeled
                self.ufn_rho_kgm3 = RectBivariateSpline(P1D_MPa, T1D_K, rho_kgm3)
                self.ufn_VP_kms = RectBivariateSpline(P1D_MPa, T1D_K, VP_kms)
                self.ufn_VS_kms = RectBivariateSpline(P1D_MPa, T1D_K, VS_kms)
                self.ufn_KS_GPa = RectBivariateSpline(P1D_MPa, T1D_K, KS_GPa)
                self.ufn_GS_GPa = RectBivariateSpline(P1D_MPa, T1D_K, GS_GPa)
                self.ufn_Cp_JkgK = RectBivariateSpline(P1D_MPa, T1D_K, Cp_JkgK)
                self.ufn_alpha_pK = RectBivariateSpline(P1D_MPa, T1D_K, alpha_pK)

                EOSlist.loaded[tableKey] = (P_MPa, T_K, self.ufn_rho_kgm3, self.ufn_VP_kms, self.ufn_VS_kms,
                                              self.ufn_KS_GPa, self.ufn_GS_GPa, self.ufn_Cp_JkgK, self.ufn_alpha_pK)
                EOSlist.ranges[tableKey] = (self.Pmin, self.Pmax, self.Tmin, self.Tmax, 
                                            self.deltaP, self.deltaT, self.EOSdeltaP, self.EOSdeltaT)
            self.EOSPmin = self.Pmin
            self.EOSPmax = self.Pmax
            self.EOSTmin = self.Tmin
            self.EOSTmax = self.Tmax
            self.rangeLabel = f'{self.Pmin},{self.Pmax},{self.deltaP},' + \
                              f'{self.Tmin},{self.Tmax},{self.deltaT}'

            # Assign thermal conductivity
            # (currently a placeholder)
            if kThermConst_WmK is None:
                if Fe_EOS:
                    kThermConst_WmK = Constants.kThermFe_WmK
                else:
                    kThermConst_WmK = Constants.kThermSil_WmK
            kTherm_WmK = np.zeros((np.size(P_MPa), np.size(T_K))) + kThermConst_WmK  # Placeholder until a self-consistent determination is implemented
            self.ufn_kTherm_WmK = RectBivariateSpline(P_MPa, T_K, kTherm_WmK)

            # Assign tidal heating function
            # (currently a placeholder)
            self.fn_Htidal_Wm3 = GetHtidalFunc(HtidalConst_Wm3)

            # Assign viscosity and porosity model functions, if applicable
            if Fe_EOS:
                # No porosity modeled
                self.ufn_phi_frac = ReturnZeros(1)

                # Assign viscosity function
                self.ufn_eta_Pas = ViscCoreUniform_Pas(etaSet_Pas=etaCoreFixed_Pas,
                                                       TviscTrans_K=TviscTrans_K)
            else:

                # Assign viscosity function
                self.ufn_eta_Pas = ViscRockUniform_Pas(etaSet_Pas=etaSilFixed_Pas,
                                                       TviscTrans_K=TviscTrans_K)

                if porosType is None or porosType == 'none':
                    # No porosity modeled, but need a dummy field for cross-compatibility
                    self.ufn_phi_frac = ReturnZeros(1)

                else:
                    if porosType == 'Vitovtova2014' or porosType == 'Chen2020':
                        # Pressure-depth lookup using the Preliminary Earth Reference Model,
                        # Dziewonski and Anderson (1981): https://doi.org/10.1016/0031-9201(81)90046-7
                        PREMzPfile = os.path.join(_ROOT, 'Thermodynamics', 'EOStables', 'PREMtable.txt')
                        zPREM_km, PPREM_kbar = np.loadtxt(PREMzPfile, skiprows=2, unpack=True, delimiter=',')
                        PPREM_MPa = PPREM_kbar * 1e3 * Constants.bar2MPa
                        self.PREMlookup = Interp1D(PPREM_MPa, zPREM_km)
                    else:
                        self.PREMlookup = None
                    self.ufn_phi_frac = GetphiFunc(self.porosType, self.phiTop_frac, self.Pclosure_MPa,
                                                  self.PREMlookup, P_MPa, T_K)

            # Store complete EOSStruct in global list of loaded EOSs
            EOSlist.loaded[self.EOSlabel] = self
            EOSlist.ranges[self.EOSlabel] = self.rangeLabel

    def _EnforceDomain(self, P_MPa, T_K):
        """Apply this EOS's domain policy before any property evaluation.

        Replaces the bare ResetNearestExtrap call so that an out-of-domain
        query can warn or raise instead of being silently saturated to the
        table edge. Clamping behavior itself is unchanged.
        """
        policy = getattr(self, 'EOSdomainPolicy', 'legacy_clamp')
        if policy == 'error' and EOSCandidatePropagationActive():
            # Scratch storage below a candidate's CMB is not scientifically
            # evaluated; strictness is applied to physical nodes afterwards.
            policy = 'legacy_clamp'
        if EOSIntrospectionActive():
            # Sentinel/type probe, not a physical evaluation -- fall back to
            # the historical clamp so introspection cannot trip a strict policy.
            policy = 'legacy_clamp'
        return EnforceEOSDomain(
            P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax,
            policy=policy,
            CLAMP=not self.EXTRAP,
            EOSname=f'{"core" if getattr(self, "Fe_EOS", False) else "silicate"} EOS {self.comp}')

    def fn_porosCorrect(self, propBulk, propPore, phi, J):
        # Combine pore fluid properties with matrix properties in accordance with
        # Yu et al. (2016): http://dx.doi.org/10.1016/j.jrmge.2015.07.004
        return (propBulk**J * (1 - phi) + propPore**J * phi) ** (1/J)

    def fn_rho_kgm3(self, P_MPa, T_K, grid=False):
        # Limit extrapolation to use nearest value from evaluated fit if desired
        P_MPa, T_K = self._EnforceDomain(P_MPa, T_K)
        return self.ufn_rho_kgm3(P_MPa, T_K, grid=grid)
    def fn_Cp_JkgK(self, P_MPa, T_K, grid=False):
        P_MPa, T_K = self._EnforceDomain(P_MPa, T_K)
        return self.ufn_Cp_JkgK(P_MPa, T_K, grid=grid)
    def fn_alpha_pK(self, P_MPa, T_K, grid=False):
        P_MPa, T_K = self._EnforceDomain(P_MPa, T_K)
        return self.ufn_alpha_pK(P_MPa, T_K, grid=grid)
    def fn_kTherm_WmK(self, P_MPa, T_K, grid=False):
        P_MPa, T_K = self._EnforceDomain(P_MPa, T_K)
        return self.ufn_kTherm_WmK(P_MPa, T_K, grid=grid)
    def fn_Seismic(self, P_MPa, T_K, grid=False):
        P_MPa, T_K = self._EnforceDomain(P_MPa, T_K)
        return (self.ufn_VP_kms(P_MPa, T_K, grid=grid),
                self.ufn_VS_kms(P_MPa, T_K, grid=grid),
                self.ufn_KS_GPa(P_MPa, T_K, grid=grid),
                self.ufn_GS_GPa(P_MPa, T_K, grid=grid))
    def fn_phi_frac(self, P_MPa, T_K, grid=False):
        P_MPa, T_K = self._EnforceDomain(P_MPa, T_K)
        return self.ufn_phi_frac(P_MPa, T_K, grid=grid)
    def fn_eta_Pas(self, P_MPa, T_K, grid=False):
        P_MPa, T_K = self._EnforceDomain(P_MPa, T_K)
        return self.ufn_eta_Pas(P_MPa, T_K, grid=grid)


def TsolidusHirschmann2000(P_MPa):
    """ Silicate melting temperature parameterization based on
        Hirschmann (2000): https://doi.org/10.1029/2000GC000070 .

        Args:
            P_MPa (float, shape N): Pressure values in MPa to evaluate
        Returns:
            Tsolidus_K (float, shape N): Solidus temperature for each P value
    """
    P_GPa = P_MPa * 1e-3
    a = -5.104
    b = 132.899
    c = 1120.661
    Tsolidus_K = a*P_GPa**2 + b*P_GPa + c + Constants.T0
    return Tsolidus_K


class GetHtidalFunc:
    def __init__(self, HtidalConst_Wm3):
        # Tidal heating as a function of density, gravity, and bulk and shear moduli
        # Currently a placeholder returning a constant value, awaiting a self-
        # consistent calculation.
        self.HtidalConst_Wm3 = HtidalConst_Wm3

    def __call__(self, rho_kgm3, g_ms2, KS_GPa, GS_GPa):
        return self.HtidalConst_Wm3


def GetphiFunc(porosType, phiTop_frac, Pclosure_MPa, PREMlookup, P1D_MPa, T1D_K):
    # Porosity models as functions of P and T, set by variable vacuum maximum
    # porosity, pore closure pressure, or model coefficients.

    # Get 2D grids of P, T values for constructing T-independent functions
    # of porosity
    P2D_MPa, T2D_K = np.meshgrid(P1D_MPa, T1D_K, indexing='ij')

    if porosType == 'Han2014':
        # Han et al. (2014) porosity model: https://doi.org/10.1002/2014GL059378
        # This is an exponential model based on a surface value and a pore closure
        # pressure. A more robust model might have a closure pressure that is
        # self-consistently determined as a function of temperature, instead of
        # using a value prescribed by the user.
        if phiTop_frac == 0 or Pclosure_MPa == 0:
            phi_frac = np.zeros_like(P2D_MPa)
        else:
            c = 6.15
            phi_frac = phiTop_frac * np.exp(-c * P2D_MPa / Pclosure_MPa)
    elif porosType == 'Vitovtova2014' or porosType == 'Chen2020':
        # Earth rock porosity parameterizations referenced to ocean worlds
        # via PREM lookup
        # Create an array of repeated (T-independent) arrays of PREM-equivalent depths
        zPREMequiv_km = np.tile(PREMlookup(P1D_MPa), (np.size(T1D_K), 1)).T
        if porosType == 'Vitovtova2014':
            # Vitovtova et al. (2014) porosity model: https://doi.org/10.1134/S1069351314040181
            # This is "generalized depth trend" model from Vitovtova et al., 2014
            # for porosity of Earth crust. Note that the equation reported in the text
            # does not match that in the figure, and the equation from the figure
            # matches the plot, so we assume the equation from the figure is
            # the correct one.
            phi_frac = 10 ** (-0.65 - 0.16 * zPREMequiv_km + 0.0019 * zPREMequiv_km ** 2)
            # Because the second-order term has a positive coefficient, porosity
            # will increase above about ~42 km PREM depth, but the porosity should
            # continue to decrease. For physical consistency, we set porosity to
            # zero beyond the minimum point.
            zPhiMin_km = 0.16 / (2 * 0.0019)
            phi_frac[zPREMequiv_km > zPhiMin_km] = 0.0
        else:
            # Chen et al. (2020) porosity model: https://doi.org/10.1007/s10040-020-02214-x
            # 1/(1 + z)^n model based on a surface value and an empirical parameterization
            # relevant to Earth. Here we use the "oceanic crust" fitting parameters.
            m = 0.008
            n = 89.53
            phi0 = 0.678
            phi_frac = phi0 / (1 + m * zPREMequiv_km) ** n
    else:
        # Invalid porosity type
        raise ValueError(f'Porosity type "{porosType}" is not supported.')

    # Create unchanging function for the expected porosity
    fn_phi_frac = RectBivariateSpline(P1D_MPa, T1D_K, phi_frac)
    return fn_phi_frac


class GetphiCalc:
    def __init__(self, phiMax_frac, fn_phiEOS_frac, phiMin_frac):
        self.phiMin_frac = phiMin_frac
        self.phiMax_frac = phiMax_frac
        # These (0, 0) calls are sentinel probes, not evaluations at a physical
        # state; suspend EOS-domain enforcement around them.
        with EOSIntrospection():
            isZeros = type(fn_phiEOS_frac(0, 0)) == ReturnZeros
        if isZeros:
            self.multFactor = 1.0
            # For some reason, binding a ReturnZeros func as a method requires an extra dummy call
            with EOSIntrospection():
                self.fn_phiEOS_frac = fn_phiEOS_frac(0, 0)
            self.EOSphiMax_frac = 0
            self.ZEROS = True
        else:
            with EOSIntrospection():
                self.EOSphiMax_frac = np.squeeze(fn_phiEOS_frac(0, 0))
            if self.EOSphiMax_frac == 0:
                self.multFactor = 1
            else:
                self.multFactor = phiMax_frac / self.EOSphiMax_frac
            self.fn_phiEOS_frac = fn_phiEOS_frac
            self.ZEROS = False

    def update(self, newPhiMax_frac):
        if not self.ZEROS:
            self.multFactor = newPhiMax_frac / self.EOSphiMax_frac

    def __call__(self, P_MPa, T_K, grid=False):
        if self.ZEROS:
            phi_frac = self.fn_phiEOS_frac(P_MPa, T_K, grid=grid)
        else:
            phi_frac = self.fn_phiEOS_frac(P_MPa, T_K, grid=grid) * self.multFactor
            if np.size(P_MPa) == 1:
                if phi_frac < self.phiMin_frac:
                    phi_frac = 0
            else:
                phi_frac[phi_frac < self.phiMin_frac] = 0

        return phi_frac


def nuPoisson(VP_kms, VS_kms):
    """
    Calculate Poisson ratio based on sound speeds after Eq. 1 of Ji et al. (2018):
    https://doi.org/10.1002/2017JB014606
    :param VP_kms: float, shape N, NxM, etc. P-wave sound speed in km/s.
    :param VS_kms: float, shape matching or broadcastable to VP_kms. S-wave sound speed
        in km/s.
    :return: nu: Dimensionless Poisson ratio, typically around 0.3 for ices and rocks.
    """
    return 0.5 * (1 - 1/((VP_kms/VS_kms)**2 - 1))


class ViscRockUniform_Pas:
    def __init__(self, etaSet_Pas=None, TviscTrans_K=None):
        if etaSet_Pas is None or any(eta is None for eta in etaSet_Pas):
            self.eta_Pas = Constants.etaRock_Pas
        else:
            self.eta_Pas = etaSet_Pas

        if TviscTrans_K is None:
            self.TviscTrans_K = Constants.TviscRock_K
        else:
            self.TviscTrans_K = TviscTrans_K

    def __call__(self, P_MPa, T_K, grid=False):
        Ttrans_K = np.insert([0.0, np.inf], 1, self.TviscTrans_K)
        if grid:
            eta_Pas = np.zeros((np.size(P_MPa), np.size(T_K)))
            for Tlow_K, Tupp_K, etaConst_Pas in zip(Ttrans_K[:-1], Ttrans_K[1:], self.eta_Pas):
                eta_Pas[:, np.logical_and(T_K >= Tlow_K, T_K < Tupp_K)] = etaConst_Pas
        else:
            eta_Pas = np.zeros_like(P_MPa)
            for Tlow_K, Tupp_K, etaConst_Pas in zip(Ttrans_K[:-1], Ttrans_K[1:], self.eta_Pas):
                eta_Pas[np.logical_and(T_K >= Tlow_K, T_K < Tupp_K)] = etaConst_Pas

        return eta_Pas


class ViscCoreUniform_Pas:
    def __init__(self, etaSet_Pas=None, TviscTrans_K=None):
        if etaSet_Pas is None or any(eta is None for eta in etaSet_Pas):
            self.eta_Pas = [Constants.etaFeSolid_Pas, Constants.etaFeLiquid_Pas]
        else:
            self.eta_Pas = etaSet_Pas

        if TviscTrans_K is None:
            self.TviscTrans_K = Constants.TviscFe_K
        else:
            self.TviscTrans_K = TviscTrans_K

    def __call__(self, P_MPa, T_K, grid=False):
        Ttrans_K = np.insert([0.0, np.inf], 1, self.TviscTrans_K)
        if grid:
            eta_Pas = np.zeros((np.size(P_MPa), np.size(T_K)))
            for Tlow_K, Tupp_K, etaConst_Pas in zip(Ttrans_K[:-1], Ttrans_K[1:], self.eta_Pas):
                eta_Pas[:, np.logical_and(T_K >= Tlow_K, T_K < Tupp_K)] = etaConst_Pas
        else:
            eta_Pas = np.zeros_like(P_MPa)
            for Tlow_K, Tupp_K, etaConst_Pas in zip(Ttrans_K[:-1], Ttrans_K[1:], self.eta_Pas):
                eta_Pas[np.logical_and(T_K >= Tlow_K, T_K < Tupp_K)] = etaConst_Pas

        return eta_Pas