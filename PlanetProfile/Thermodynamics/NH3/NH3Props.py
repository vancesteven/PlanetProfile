"""Aqueous-ammonia ocean properties via CoolProp (HEOS mixture).

Reinstates the NH3-ocean capability of the original (2018, MATLAB)
PlanetProfile, which used the REFPROP ammonia-water implementation.
CoolProp 8.x ships accurate PURE Helmholtz EOS for water (IAPWS-95) and
ammonia (Tillner-Roth), but no fitted binary-pair data for the mixture
— we enable the 'linear' reducing-parameter estimation scheme
(CP.apply_simple_mixing_rule) and validate it against the Melinder
(2010) aqueous-ammonia polynomials (CoolProp INCOMP::MAM):
density within ~1%, Cp within 2-5% for X <= 0.2 at low P
(tests/coolprop_nh3_test.py). The Lorentz-Berthelot scheme agrees with
'linear' to ~0.1% throughout the ocean domain, so the reducing-rule
choice is not the dominant uncertainty; the unfitted mixture DEPARTURE
function is, and results should be quoted with the few-percent Cp
caveat until a Tillner-Roth & Friend (1998)-grade binary model lands
in CoolProp.

Freezing (liquidus): the ammonia freezing-point depression of
Leliwa-Kopystynski et al. (2002, Icarus 159, 518), whose fitted
liquidus subtracts the SAME composition term from every ice-phase
melting branch:

    dT(X) = 53.8*X + 650*X^3   [K], X = NH3 mass fraction

valid for X <= ~0.23 (their fit range; the cubic underestimates the
~100 K depression approaching the eutectic X ~ 0.32). PlanetProfile
applies this as a shift of the SeaFreeze pure-water phase diagram:
phase(P, T; X) = phase_H2O(P, T + dT(X)), i.e. ices are pure H2O and
NH3 is rejected into the residual liquid (standard assumption).

Validated CoolProp flash coverage (probed 2026-07-25): 96-99% of a
230-330 K x 0.5-1500 MPa x X<=0.25 grid; isolated flash failures are
filled by 1-D interpolation along isotherms and flagged in the log.
"""
from __future__ import annotations

import logging
import numpy as np

log = logging.getLogger('PlanetProfile')

M_NH3_gmol = 17.031
M_H2O_gmol = 18.015

# Validated domain (see module docstring); HydroEOS clamps to these.
NH3_PMAX_MPA = 1500.0
NH3_TMIN_K = 230.0
NH3_TMAX_K = 330.0
NH3_WMAX_PPT = 250.0  # L-K liquidus fit bound (X ~ 0.23) with margin

_MIX_READY = False

# Liquidus mode: 'mu' (default) solves the chemical-potential equality
# per pressure/ice phase (muLiquidusCurve_K); 'LK2002' applies the
# scalar Leliwa-Kopystynski shift (see its accuracy caution).
NH3_MELT_MODE = 'mu'


def NH3liquidusShift_K(w_ppt):
    """Leliwa-Kopystynski et al. (2002) freezing-point depression (K)
    for NH3 mass fraction X = w_ppt/1000: dT = 53.8 X + 650 X^3.

    CAUTION (2026-07-25 cross-check): this polynomial UNDERESTIMATES the
    depression by ~2x at low-to-mid X relative to both the ideal
    colligative limit and the chemical-potential route below (e.g.
    5 wt%: 2.8 K here vs ~5.4 K mu-based vs ~5.7 K ideal-colligative);
    the coefficients/variable convention should be re-verified against
    the original paper before quantitative use. Retained as the
    'LK2002' fallback melt mode; the DEFAULT is the mu-based liquidus
    (NH3MuLiquidusShift_K)."""
    X = float(w_ppt) / 1000.0
    return 53.8 * X + 650.0 * X ** 3


def _sfG_scatter(mat, P_MPa, T_K):
    """SeaFreeze Gibbs energy (J/kg) at paired (P, T) points."""
    import seafreeze.seafreeze as sf
    pts = np.array(list(zip(np.atleast_1d(P_MPa).astype(float),
                            np.atleast_1d(T_K).astype(float))),
                   dtype='f,f').astype(object)
    return np.ravel(sf.getProp(pts, mat).G).astype(float)


def pureMeltingCurve_K(P_MPa):
    """Pure-water melting temperature and subjacent ice phase code at
    each pressure, from the SeaFreeze phase diagram (scan + refine)."""
    import seafreeze.seafreeze as sf
    P_MPa = np.atleast_1d(np.asarray(P_MPa, float))
    Tm = np.full(P_MPa.size, np.nan)
    ice = np.ones(P_MPa.size, dtype=int)
    Tg = np.linspace(239.5, 380.0, 300)
    for i, P in enumerate(P_MPa):
        grid = np.array([(P, t) for t in Tg], dtype='f,f').astype(object)
        ph = np.ravel(sf.whichphase(grid)).astype(int)
        liq = np.where(ph == 0)[0]
        if not liq.size:
            continue
        i0 = liq[0]
        if i0 == 0:
            Tm[i] = Tg[0]
            continue
        lo, hi = Tg[i0 - 1], Tg[i0]
        for _ in range(12):  # bisect the boundary to ~0.03 K
            mid = 0.5 * (lo + hi)
            g = np.array([(P, mid)], dtype='f,f').astype(object)
            if int(np.ravel(sf.whichphase(g))[0]) == 0:
                hi = mid
            else:
                lo = mid
        Tm[i] = 0.5 * (lo + hi)
        ice[i] = int(ph[i0 - 1])
    return Tm, ice


_SF_ICE_NAME = {1: 'Ih', 2: 'II', 3: 'III', 5: 'V', 6: 'VI'}
R_JMOLK = 8.31446


def _ln_aw(P_MPa, T_K, z_nh3):
    """ln water activity in the NH3 solution from CoolProp: chemical
    potential of water in the mixture minus pure water, evaluated in
    the SAME (HEOS) framework so reference states cancel."""
    import CoolProp.CoolProp as CP
    AS = CP.AbstractState('HEOS', 'Ammonia&Water')
    AS.set_mole_fractions([z_nh3, 1.0 - z_nh3])
    AS.specify_phase(CP.iphase_liquid)
    AS.update(CP.PT_INPUTS, P_MPa * 1e6, T_K)
    mu_mix = AS.chemical_potential(1)
    ASp = CP.AbstractState('HEOS', 'Water')
    ASp.specify_phase(CP.iphase_liquid)
    ASp.update(CP.PT_INPUTS, P_MPa * 1e6, T_K)
    return (mu_mix - ASp.chemical_potential(0)) / (R_JMOLK * T_K)


def muLiquidusCurve_K(P_MPa, w_ppt):
    """NH3-H2O liquidus T at each pressure by chemical-potential
    equality (user-preferred alternative to the L-K polynomial):

        [G_ice(P,T) - G_water(P,T)]_SeaFreeze * M_w = R T ln a_w(P,T,x)

    Left side: SeaFreeze ice and (metastable) liquid water Gibbs
    energies — same framework, references cancel. Right side: water
    activity from the CoolProp mixture (mixture-minus-pure difference,
    same framework, references cancel). Reduces to the pure melting
    curve as x -> 0; reproduces the ideal colligative dilute limit
    (~5.4 K at 5 wt%, 1 bar) and applies self-consistently under every
    ice phase (Ih through VI).

    Solved by sign-change scan + bisection on a T grid bounded below by
    SeaFreeze's water1 validity (239.5 K); pressures where no bracket
    is found (deep depression beyond validity, or CoolProp flash
    failure) return the L-K value as fallback with a log notice.
    """
    _ensure_mixture()
    z = _molefrac_NH3(w_ppt)
    P_MPa = np.atleast_1d(np.asarray(P_MPa, float))
    Tm_pure, ice_codes = pureMeltingCurve_K(P_MPa)
    Tliq = np.full(P_MPa.size, np.nan)
    n_fallback = 0
    for i, (P, T0, code) in enumerate(zip(P_MPa, Tm_pure, ice_codes)):
        if not np.isfinite(T0):
            continue
        ice = _SF_ICE_NAME.get(int(code), 'Ih')
        Tlo = max(239.6, T0 - 45.0)
        Tscan = np.linspace(Tlo, min(T0 + 0.5, 380.0), 24)
        try:
            dG = (_sfG_scatter(ice, np.full(Tscan.size, P), Tscan)
                  - _sfG_scatter('water1', np.full(Tscan.size, P),
                                 Tscan)) * (M_H2O_gmol * 1e-3)
            fvals = np.full(Tscan.size, np.nan)
            for k, T in enumerate(Tscan):
                try:
                    fvals[k] = dG[k] - R_JMOLK * T * _ln_aw(P, T, z)
                except Exception:
                    pass
            ok = np.isfinite(fvals)
            sign = np.sign(fvals)
            idx = None
            oki = np.where(ok)[0]
            for a, b in zip(oki[:-1], oki[1:]):
                if sign[a] != sign[b]:
                    idx = (a, b)
                    break
            if idx is None:
                raise RuntimeError('no bracket')
            lo, hi = Tscan[idx[0]], Tscan[idx[1]]
            flo = fvals[idx[0]]
            for _ in range(30):
                mid = 0.5 * (lo + hi)
                dGm = (_sfG_scatter(ice, [P], [mid])[0]
                       - _sfG_scatter('water1', [P], [mid])[0]) \
                    * (M_H2O_gmol * 1e-3)
                fm = dGm - R_JMOLK * mid * _ln_aw(P, mid, z)
                if np.sign(fm) == np.sign(flo):
                    lo, flo = mid, fm
                else:
                    hi = mid
                if hi - lo < 5e-3:
                    break
            Tliq[i] = 0.5 * (lo + hi)
        except Exception:
            Tliq[i] = T0 - NH3liquidusShift_K(w_ppt)
            n_fallback += 1
    if n_fallback:
        log.info(f'mu-based NH3 liquidus: L-K fallback at {n_fallback}/'
                 f'{P_MPa.size} pressures (bracket outside SeaFreeze '
                 'validity or CoolProp flash failure).')
    return Tliq, Tm_pure


def _ensure_mixture():
    """Register the estimated NH3-H2O binary pair once per process."""
    global _MIX_READY
    if _MIX_READY:
        return
    import CoolProp.CoolProp as CP
    CP.set_config_bool(CP.OVERWRITE_BINARY_INTERACTION, True)
    CP.apply_simple_mixing_rule('Ammonia', 'Water', 'linear')
    _MIX_READY = True


def _molefrac_NH3(w_ppt):
    X = float(w_ppt) / 1000.0
    return (X / M_NH3_gmol) / (X / M_NH3_gmol + (1.0 - X) / M_H2O_gmol)


def _fill_nan_1d(y, x):
    """Linear-interpolate over nan gaps along one axis; edge nans take
    the nearest finite value. Returns the fill count."""
    bad = ~np.isfinite(y)
    if not bad.any() or bad.all():
        return 0 if not bad.any() else -1
    y[bad] = np.interp(x[bad], x[~bad], y[~bad])
    return int(bad.sum())


def NH3Props(P_MPa, T_K, w_ppt):
    """Liquid ammonia-water property grids on (P_MPa x T_K).

    Returns (rho_kgm3, Cp_JkgK, alpha_pK, VP_kms, KS_GPa), each of
    shape (nP, nT). Isolated CoolProp flash failures are interpolated
    over along P at fixed T (then along T if a full isotherm failed)
    and reported via the log.
    """
    import CoolProp.CoolProp as CP
    _ensure_mixture()
    z = _molefrac_NH3(w_ppt)
    P_MPa = np.asarray(P_MPa, float)
    T_K = np.asarray(T_K, float)
    nP, nT = P_MPa.size, T_K.size
    rho = np.full((nP, nT), np.nan)
    Cp = np.full((nP, nT), np.nan)
    alpha = np.full((nP, nT), np.nan)
    vP = np.full((nP, nT), np.nan)

    AS = CP.AbstractState('HEOS', 'Ammonia&Water')
    AS.set_mole_fractions([z, 1.0 - z])
    AS.specify_phase(CP.iphase_liquid)
    for j, T in enumerate(T_K):
        for i, P in enumerate(P_MPa):
            try:
                AS.update(CP.PT_INPUTS, P * 1e6, T)
                r = AS.rhomass()
                if not (500.0 < r < 2000.0):
                    continue  # unconverged flash masquerading as liquid
                rho[i, j] = r
                Cp[i, j] = AS.cpmass()
                alpha[i, j] = AS.isobaric_expansion_coefficient()
                vP[i, j] = AS.speed_sound()
            except Exception:
                pass  # nan-filled below

    # Invalidate spurious liquid-branch roots: in the deep sub-liquidus
    # (solid) region CoolProp's forced-liquid flash can return a bogus
    # low-density root (e.g. 884 kg/m3 at 260 K / 1000 MPa) that would
    # poison the property splines. Liquid rho must be non-decreasing
    # with P along an isotherm — cells that break this are discarded and
    # refilled below (smooth metastable placeholders; the phase function
    # marks the region solid, so ocean profiles never sample it).
    for j in range(nT):
        run_max = -np.inf
        for i in range(nP):
            if not np.isfinite(rho[i, j]):
                continue
            if rho[i, j] < run_max - 1e-6:
                rho[i, j] = Cp[i, j] = alpha[i, j] = vP[i, j] = np.nan
            else:
                run_max = rho[i, j]

    n_fill = 0
    for arr in (rho, Cp, alpha, vP):
        for i in range(nP):   # fill along T first (cold+deep failures
            f = _fill_nan_1d(arr[i, :], T_K)  # borrow warmer isotherms)
            n_fill += max(f, 0)
        for j in range(nT):
            f = _fill_nan_1d(arr[:, j], P_MPa)
            n_fill += max(f, 0)
    if n_fill:
        log.debug(f'CoolProp NH3-H2O grid: interpolated over {n_fill} '
                  f'failed flash points (of {4 * nP * nT}).')
    if not np.all(np.isfinite(rho)):
        raise RuntimeError(
            'CoolProp NH3-H2O EOS failed across the requested grid '
            f'(P [{P_MPa.min()}, {P_MPa.max()}] MPa, '
            f'T [{T_K.min()}, {T_K.max()}] K, w = {w_ppt} ppt).')

    KS_GPa = rho * vP ** 2 * 1e-9  # adiabatic bulk modulus
    VP_kms = vP * 1e-3
    return rho, Cp, alpha, VP_kms, KS_GPa


class NH3Seismic:
    """(P_MPa, T_K, grid=) -> (VP_kms, KS_GPa), mirroring SFSeismic."""
    def __init__(self, P_MPa, T_K, VP_kms, KS_GPa, EXTRAP):
        from scipy.interpolate import RectBivariateSpline
        self.EXTRAP = EXTRAP
        self.ufn_VP_kms = RectBivariateSpline(P_MPa, T_K, VP_kms)
        self.ufn_KS_GPa = RectBivariateSpline(P_MPa, T_K, KS_GPa)

    def __call__(self, P_MPa, T_K, grid=False):
        return (self.ufn_VP_kms(P_MPa, T_K, grid=grid),
                self.ufn_KS_GPa(P_MPa, T_K, grid=grid))
