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
# Composition bound: Melinder (2010) anchor validity ends at X = 0.30
# and the NH3.2H2O eutectic/peritectic sits at X ~ 0.326 (176.2 K);
# 250 ppt keeps a margin inside both. (An earlier note attributed this
# bound to a Leliwa-Kopystynski 2002 polynomial; that polynomial's
# provenance could not be verified and its dilute-limit slope is
# provably ~2x too shallow — it has been removed. See the
# scientific-review record in plans/HANDOFF-2026-07-26-nh3-liquidus-
# defect.md, Machine A addendum 2026-07-28.)
NH3_WMAX_PPT = 250.0

_MIX_READY = False

# Activity model for the mu-equality liquidus:
#   'melinder-CG' (default): Redlich-Kister excess for the water
#       activity coefficient, AMPLITUDE anchored to the Melinder (2010)
#       experimental freezing-point curve (CoolProp INCOMP::MAM) at
#       1 bar, P,T SHAPE from the Choukroun & Grasset (2010) Margules
#       factor used by the original 2018 PlanetProfile
#       (Thermodynamics/NH3aq/getIcePhaseNH3.m). Analytic — no mixture
#       flash — so it evaluates at every (P, T) under all ice phases.
#       1-bar liquidus reproduces Melinder to < 0.05 K over
#       X in [0.02, 0.175].
#   'coolprop': the CoolProp HEOS estimated-mixture activity (legacy;
#       retained for regression comparison only). Its excess term has
#       the WRONG SIGN (gamma_w > 1; NH3-H2O requires gamma_w < 1), so
#       it under-depresses the liquidus by 9-37% over 30-150 ppt, and
#       its flash fails across much of the 200-400 MPa HP-ice band.
NH3_ACTIVITY_MODE = 'melinder-CG'

# Redlich-Kister coefficients for ln gamma_w = u^2 (a0 + a1 u + a2 u^2)
# at the 1-bar anchor (u = NH3 mole fraction), fitted 2026-07-28 to the
# Melinder anchor table below; independently reproduced on this machine
# (matches the scientific reviewer's g_req values to < 1e-5; max
# reconstruction residual 4.1e-5 in ln a_w ~ 0.005 K):
#   X (mass)  0.02      0.03      0.05      0.075     0.10
#   g_req    -0.00075  -0.00163  -0.00460  -0.01094  -0.02083
#   X (mass)  0.125     0.15      0.175
#   g_req    -0.03517  -0.05510  -0.08204
# Anchors are self-referenced to Melinder's own zero (T_freeze(0) =
# 273.1218 K) and applied at SeaFreeze's melting point (273.1527 K) —
# anchoring against 273.15 directly would inject a spurious excess
# that does not vanish as X -> 0.
_RK_A0, _RK_A1, _RK_A2 = -1.552208, -0.841384, -21.927892
_ANCHOR_P0_MPA = 0.1


def _cg_shape(P_MPa, T_K):
    """Choukroun & Grasset (2010) Margules P,T factor (dimensionless,
    used only as a NORMALIZED shape): (1 + 0.4 tanh(6e-3 P)) *
    (1 + 4.62e-4 T - 113/T). Saturates at 1.4x by ~500 MPa; monotone
    and bounded through 1 GPa; T factor has no zero above 108 K."""
    return ((1.0 + 0.4 * np.tanh(6.0e-3 * np.asarray(P_MPa, float)))
            * (1.0 + 4.62e-4 * np.asarray(T_K, float)
               - 113.0 / np.asarray(T_K, float)))


def _ln_aw_anchored(P_MPa, T_K, u, Tref_K):
    """ln water activity, Melinder-anchored model: ln gamma_w scales
    the 1-bar Redlich-Kister amplitude by the CG2010 shape normalized
    at (P0, Tref), plus the ideal ln x_w. Exact colligative limit by
    construction (ln gamma_w -> 0 as u^2)."""
    g_anchor = u * u * (_RK_A0 + _RK_A1 * u + _RK_A2 * u * u)
    ln_gam = g_anchor * _cg_shape(P_MPa, T_K) \
        / _cg_shape(_ANCHOR_P0_MPA, Tref_K)
    return ln_gam + np.log(1.0 - u)


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


_AS_CACHE = {}


def _ln_aw(P_MPa, T_K, z_nh3):
    """ln water activity in the NH3 solution from CoolProp: chemical
    potential of water in the mixture minus pure water, evaluated in
    the SAME (HEOS) framework so reference states cancel.

    GUARDED (Machine B defect report 2026-07-26): the forced-liquid
    flash of the estimated binary occasionally returns garbage at
    isolated (P, T) points, which used to seed spurious liquidus roots.
    A flash is rejected (ValueError) unless the mixture density lands in
    the liquid window AND ln a_w is physically plausible for x <= 0.25
    (ideal value ~ -x; window [-1.0, +0.02])."""
    import CoolProp.CoolProp as CP
    key = round(float(z_nh3), 10)
    if key not in _AS_CACHE:
        AS = CP.AbstractState('HEOS', 'Ammonia&Water')
        AS.set_mole_fractions([z_nh3, 1.0 - z_nh3])
        AS.specify_phase(CP.iphase_liquid)
        ASp = CP.AbstractState('HEOS', 'Water')
        ASp.specify_phase(CP.iphase_liquid)
        _AS_CACHE[key] = (AS, ASp)
    AS, ASp = _AS_CACHE[key]
    ASp.update(CP.PT_INPUTS, P_MPa * 1e6, T_K)
    rhop = ASp.rhomass()
    if not (500.0 < rhop < 2000.0):
        raise ValueError(f'bad pure flash (rho={rhop:.0f})')
    try:
        AS.update(CP.PT_INPUTS, P_MPa * 1e6, T_K)
        rho = AS.rhomass()
        if not (500.0 < rho < 2000.0):
            raise ValueError(f'bad mixture flash (rho={rho:.0f})')
    except Exception:
        # Known CoolProp mixture flash-failure pockets (e.g. ~250 K /
        # 200 MPa at low x): retry with a density guess seeded from the
        # pure-water liquid root at the same (P, T).
        g = CP.PyGuessesStructure()
        g.rhomolar = rhop / ASp.molar_mass() * 0.99
        AS.update_with_guesses(CP.PT_INPUTS, P_MPa * 1e6, T_K, g)
        rho = AS.rhomass()
        if not (500.0 < rho < 2000.0):
            raise ValueError(f'bad guessed mixture flash (rho={rho:.0f})')
    mu_mix = AS.chemical_potential(1)
    val = (mu_mix - ASp.chemical_potential(0)) / (R_JMOLK * T_K)
    if not (-1.0 < val < 0.02):
        raise ValueError(f'implausible ln a_w = {val:.3f}')
    return val


def muLiquidusCurve_K(P_MPa, w_ppt):
    """NH3-H2O liquidus T at each pressure by chemical-potential
    equality:

        [G_ice(P,T) - G_water(P,T)]_SeaFreeze * M_w = R T ln a_w(P,T,x)

    Dispatches on NH3_ACTIVITY_MODE (see its docstring): the default
    'melinder-CG' activity is analytic, so the solve is a plain
    monotone bisection per pressure — no flash guards, no repairs.
    Returns (Tliq_K, TmPure_K) arrays matching P_MPa."""
    if NH3_ACTIVITY_MODE == 'melinder-CG':
        return _muLiquidusAnchored_K(P_MPa, w_ppt)
    return _muLiquidusCoolProp_K(P_MPa, w_ppt)


def _anchoredTref_K(u):
    """1-bar anchored liquidus for NH3 mole fraction u: the reference
    temperature normalizing the CG2010 shape. Fixed-point iteration
    (the shape varies slowly in T, 3 iterations converge < 1e-3 K)."""
    Tm0 = float(pureMeltingCurve_K(np.array([_ANCHOR_P0_MPA]))[0][0])
    Tref = Tm0
    for _ in range(3):
        Tref = _bisectAnchored_K(_ANCHOR_P0_MPA, Tm0, 'Ih', u, Tref)
        if not np.isfinite(Tref):
            raise RuntimeError(
                f'anchored NH3 liquidus: no 1-bar reference solution '
                f'for x_NH3 = {u:.4f}')
    return Tref


def _bisectAnchored_K(P, Tm, ice, u, Tref, Tfloor=239.6):
    """Warmest root of f(T) = dG_fus*M_w - R T ln a_w (anchored model)
    below the pure melting point Tm. f(Tm-) > 0 (liquid stable, since
    ln a_w < 0); f decreases monotonically as T falls. Returns
    Tfloor + 0.05 if still liquid at the SeaFreeze water1 validity
    floor (true liquidus below the representable domain), nan only if
    the bracket itself fails (raised on by the caller)."""
    def f(T):
        dG = (_sfG_scatter(ice, [P], [T])[0]
              - _sfG_scatter('water1', [P], [T])[0]) * (M_H2O_gmol * 1e-3)
        return dG - R_JMOLK * T * _ln_aw_anchored(P, T, u, Tref)
    hi = Tm - 1e-3
    lo = max(Tfloor + 0.05, Tm - 60.0)
    if f(lo) > 0:
        return lo if lo == Tfloor + 0.05 else np.nan
    fhi = f(hi)
    if fhi < 0:
        return np.nan
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < 2e-3:
            break
    return 0.5 * (lo + hi)


def _muLiquidusAnchored_K(P_MPa, w_ppt):
    """Analytic (Melinder-anchored) liquidus: per-pressure monotone
    bisection against the SeaFreeze Gibbs difference. Raises on any
    unsolved pressure — an analytic activity leaves nothing to repair,
    so a failure means the domain itself is wrong."""
    u = _molefrac_NH3(w_ppt)
    P_MPa = np.atleast_1d(np.asarray(P_MPa, float))
    Tm_pure, ice_codes = pureMeltingCurve_K(P_MPa)
    if w_ppt <= 0:
        return Tm_pure.copy(), Tm_pure
    Tref = _anchoredTref_K(u)
    Tliq = np.full(P_MPa.size, np.nan)
    for i in range(P_MPa.size):
        if not np.isfinite(Tm_pure[i]):
            continue
        ice = _SF_ICE_NAME.get(int(ice_codes[i]), 'Ih')
        Tliq[i] = _bisectAnchored_K(P_MPa[i], Tm_pure[i], ice, u, Tref)
    bad = np.isfinite(Tm_pure) & ~np.isfinite(Tliq)
    if np.any(bad):
        raise RuntimeError(
            'anchored NH3 liquidus failed at P (MPa): '
            + ', '.join(f'{p:.0f}' for p in P_MPa[bad][:8])
            + f' (w = {w_ppt} ppt) — refusing to return a defective '
            'liquidus.')
    return Tliq, Tm_pure


def _muLiquidusCoolProp_K(P_MPa, w_ppt):
    """LEGACY CoolProp-activity liquidus (NH3_ACTIVITY_MODE =
    'coolprop'), retained for regression comparison only — its excess
    term has the wrong sign (see NH3_ACTIVITY_MODE docstring).

        [G_ice(P,T) - G_water(P,T)]_SeaFreeze * M_w = R T ln a_w(P,T,x)

    Left side: SeaFreeze ice and (metastable) liquid water Gibbs
    energies — same framework, references cancel. Right side: water
    activity from the CoolProp mixture (mixture-minus-pure difference,
    same framework, references cancel). Reduces to the pure melting
    curve as x -> 0; reproduces the ideal colligative dilute limit
    (~5.4 K at 5 wt%, 1 bar) and applies self-consistently under every
    ice phase (Ih through VI).

    Solver (rebuilt after the Machine B defect report, plans/HANDOFF-
    2026-07-26-nh3-liquidus-defect.md): guarded activity evaluations
    (bad CoolProp flashes rejected, brackets never built across them),
    warm-to-cold scan so the WARMEST sign change — the physical
    liquidus — is taken, continuation in P seeded from the previous
    pressure's root, and a per-ice-segment smoothness filter that
    repairs outliers (> 0.3 K off the local median) by interpolation
    across good pressures. The anchored analytic model provides the
    very first search seed; a segment with < 40% cleanly solved
    pressures raises instead of returning a defective curve. Bounded
    below by SeaFreeze water1 validity (239.6 K).
    """
    _ensure_mixture()
    z = _molefrac_NH3(w_ppt)
    P_MPa = np.atleast_1d(np.asarray(P_MPa, float))
    Tm_pure, ice_codes = pureMeltingCurve_K(P_MPa)

    def _law_robust(P, T, Tcap):
        """ln a_w with a nearest-valid-T substitution: CoolProp's
        mixture flash has dead zones (e.g. ~245-249 K at ~210 MPa) that
        can blanket the root, but the water activity is nearly
        T-independent (ideal-solution ln x_w exactly so; observed
        variation < 1e-3 over several K), so borrowing the value from
        the nearest temperature where the flash converges introduces
        far less error than the 0.3 K smoothness tolerance."""
        try:
            return _ln_aw(P, T, z)
        except Exception:
            pass
        for dTs in np.arange(0.5, 10.1, 0.5):
            for Tt in (T + dTs, T - dTs):
                if Tt > Tcap + 0.5 or Tt < 239.6:
                    continue
                try:
                    return _ln_aw(P, Tt, z)
                except Exception:
                    continue
        raise ValueError(f'no valid activity near P={P}, T={T}')

    def _f_at(P, T, ice, Tcap=380.0):
        """f(T) = dG_fus*M - R T ln a_w; raises on a bad activity."""
        dG = (_sfG_scatter(ice, [P], [T])[0]
              - _sfG_scatter('water1', [P], [T])[0]) * (M_H2O_gmol * 1e-3)
        return dG - R_JMOLK * T * _law_robust(P, T, Tcap)

    def _solve_at(P, T0, ice, Tguess):
        """Highest-T sign change near Tguess (widening on demand), then
        guarded bisection. Returns nan when no valid bracket exists.
        The liquidus is the WARMEST root — noise-seeded colder
        sign-change pairs are ignored by scanning from the warm end."""
        Tcap = min(T0 + 0.5, 380.0)
        Tfloor = max(239.6, T0 - 45.0)
        for half, step in ((6.0, 0.5), (20.0, 0.5), (45.0, 1.0)):
            lo_edge = max(Tfloor, min(Tguess - half, Tcap - 2.0))
            hi_edge = min(Tcap, Tguess + half)
            if hi_edge <= lo_edge:
                continue
            Ts = np.arange(hi_edge, lo_edge - step / 2, -step)  # warm→cold
            fprev = Tprev = None
            for T in Ts:
                try:
                    fT = _f_at(P, T, ice)
                except Exception:
                    fprev = Tprev = None  # break bracket across bad pts
                    continue
                if fprev is not None and np.sign(fT) != np.sign(fprev):  # noqa: E501  (bracket found)
                    lo, hi = T, Tprev            # lo colder than hi
                    flo = fT
                    for _ in range(40):
                        mid = 0.5 * (lo + hi)
                        try:
                            fm = _f_at(P, mid, ice)
                        except Exception:
                            mid += 0.07 * (hi - lo)  # nudge off bad pt
                            try:
                                fm = _f_at(P, mid, ice)
                            except Exception:
                                break
                        if np.sign(fm) == np.sign(flo):
                            lo, flo = mid, fm
                        else:
                            hi = mid
                        if hi - lo < 5e-3:
                            break
                    return 0.5 * (lo + hi)
                fprev, Tprev = fT, T
        # No sign change found: if the solution is still LIQUID at the
        # SeaFreeze water1 validity floor (f > 0 there), the true
        # liquidus lies below the representable domain — clamp to the
        # floor so the phase grid marks everything above it liquid.
        # (Understates the depression only for T below the floor, which
        # no EOS grid samples.)
        try:
            if _f_at(P, Tfloor + 0.05, ice, Tcap) > 0:
                return Tfloor + 0.05
        except Exception:
            pass
        return np.nan

    Tliq = np.full(P_MPa.size, np.nan)
    prev_dT = None
    prev_code = None
    for i in np.argsort(P_MPa):
        P, T0, code = P_MPa[i], Tm_pure[i], int(ice_codes[i])
        if not np.isfinite(T0):
            continue
        ice = _SF_ICE_NAME.get(code, 'Ih')
        # Continuation: the liquidus depression is smooth in P within an
        # ice-phase segment — seed the search from the previous root.
        g = prev_dT if (prev_dT is not None and code == prev_code) \
            else (prev_dT if prev_dT is not None
                  else float(Tm_pure[i]
                             - _bisectAnchored_K(
                                 P, Tm_pure[i],
                                 _SF_ICE_NAME.get(code, 'Ih'),
                                 z, Tm_pure[i])))
        Tliq[i] = _solve_at(P, T0, ice, T0 - g)
        if np.isfinite(Tliq[i]):
            prev_dT = T0 - Tliq[i]
            prev_code = code
    # Per-ice-segment outlier repair (Machine B defect report
    # 2026-07-26): dTliq(P) must be smooth within a segment. Points
    # deviating > 0.3 K from the local (5-point) median — or unsolved —
    # are replaced by interpolation across the segment's good points,
    # NEVER by the L-K polynomial. Per-segment counts are logged; a
    # segment with < 40% good points raises (a defective liquidus must
    # not silently reach cache builds).
    dT = Tm_pure - Tliq
    n_repair = 0
    for code in np.unique(ice_codes[np.isfinite(Tm_pure)]):
        seg = np.where((ice_codes == code) & np.isfinite(Tm_pure))[0]
        seg = seg[np.argsort(P_MPa[seg])]
        if seg.size < 2:
            continue
        d = dT[seg]
        good = np.isfinite(d)
        if good.sum() >= 3:
            med = np.array([
                np.nanmedian(d[max(0, k - 2):k + 3]) for k in range(d.size)])
            good &= np.abs(d - med) <= 0.3
        if good.sum() < max(2, int(0.4 * seg.size)):
            raise RuntimeError(
                f'mu-based NH3 liquidus: only {int(good.sum())}/{seg.size} '
                f'pressures solved cleanly under ice {code} at '
                f'w = {w_ppt} ppt — refusing to return a defective '
                'liquidus (see plans/HANDOFF-2026-07-26-nh3-liquidus-'
                'defect.md).')
        if not good.all():
            bad = seg[~good]
            log.debug('mu-based NH3 liquidus: repaired pressures (MPa): '
                      + ', '.join(f'{P_MPa[b]:.0f}' for b in bad))
            d_fix = np.interp(P_MPa[seg][~good], P_MPa[seg][good],
                              d[good])
            dT[bad] = d_fix
            Tliq[bad] = Tm_pure[bad] - d_fix
            n_repair += bad.size
    if n_repair:
        frac = n_repair / max(1, int(np.isfinite(Tm_pure).sum()))
        msg = (f'mu-based NH3 liquidus: {n_repair} of '
               f'{int(np.isfinite(Tm_pure).sum())} pressures repaired by '
               'segment interpolation (guarded-flash/bracket failures).')
        (log.warning if frac > 0.05 else log.info)(msg)
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
