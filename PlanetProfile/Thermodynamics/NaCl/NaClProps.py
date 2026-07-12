"""
Aqueous NaCl electrical conductivity for icy-moon subsurface oceans.

Implements the numerical model of Pan et al. (2021), "Electrical Conductivity
of Aqueous NaCl at High Pressure and Low Temperature: Application to Deep
Subsurface Oceans of Icy Moons", Geophys. Res. Lett. 48, e2021GL094020,
https://doi.org/10.1029/2021GL094020 (experiments: 212-1713 MPa, 233-295 K,
1/5/10 wt% NaCl(aq); model R^2 = 0.993).

Model provenance: the published Equation 3 typography does not survive PDF
text extraction, so the implementation below reproduces the AUTHORS' OWN
regression spreadsheet (Mendeley Data, https://doi.org/10.17632/g43xkvm3gx.6,
sheet 'regression'), which evaluates to the published sigma_calc values:

    log10(sigma) = A + B/T + C*log10(chat) + D*log10(rho)
                   + 2*log10(Lambda0) - log10(10*Lambda0 + 5584.2*chat_k)

    Lambda0(T, rho) = 3101.85 - 100.45*rho - 1481175.59/T + 185698255.2/T^2
    chat   = c * rho0 / (1 + M_NaCl*c)          (reduced concentration)
    chat_k = c * rho0 / sqrt(1 + M_NaCl*c)      (Kohlrausch-term concentration)

    A = -2.4857161222849897    (intercept)
    B = +902.6050662042346     (1/T coefficient; note the PDF's abstract
                                renders the sign ambiguously - the spreadsheet
                                regression table is ground truth)
    C = +1.7644668238775978    (log10 chat coefficient)
    D = +1.5569422624900764    (log10 rho coefficient)
    rho0 = 1.004205 g/cm^3, M_NaCl = 0.05844 kg/mol

with sigma in S/m, T in K, c the NaCl molality in mol/kg, rho the WATER
density in g/cm^3 (the paper deliberately uses pure-water density from
SeaFreeze rather than solution density), Lambda0 the limiting equivalent
conductance in S cm^2 equiv^-1.

Validated against the authors' 69-row regression dataset in
tests/nacl_conductivity_test.py.

Validity: P in [212, 1713] MPa, T in [233, 295] K, molality <= 1.90 mol/kg
(10 wt%). Below ~212 MPa the paper's MD simulations show a conductivity
maximum the regression does not capture; evaluation there is a smooth
extrapolation and a warning is logged once. Callisto's 100 ppt NaCl ocean
(1.90 mol/kg) sits exactly at the fitted upper concentration.
"""
import logging

import numpy as np
import seafreeze.seafreeze as sfz

from PlanetProfile.Utilities.defineStructs import Constants

log = logging.getLogger('PlanetProfile')

# Regression coefficients (authors' spreadsheet, sheet 'regression').
_A_INTERCEPT = -2.4857161222849897
_B_INV_T = 902.6050662042346
_C_LOG_CHAT = 1.7644668238775978
_D_LOG_RHO = 1.5569422624900764
_RHO0_GCM3 = 1.004205
_M_NACL_KGMOL = 0.05844
_KOHLRAUSCH = 5584.2

# Experimental validity envelope.
PAN2021_P_MPA = (212.0, 1713.0)
PAN2021_T_K = (233.0, 295.0)
PAN2021_MOLAL_MAX = 1.902


def Pan2021_sigma_Sm(molality_molkg, T_K, rho_gcm3):
    """Pan et al. (2021) NaCl(aq) conductivity from molality, T, and water
    density. Vectorized over T/rho; molality scalar or broadcastable.

    Returns sigma in S/m.
    """
    c = np.asarray(molality_molkg, dtype=np.float64)
    T = np.asarray(T_K, dtype=np.float64)
    rho = np.asarray(rho_gcm3, dtype=np.float64)

    denom = 1.0 + _M_NACL_KGMOL * c
    chat = c * _RHO0_GCM3 / denom
    chat_k = c * _RHO0_GCM3 / np.sqrt(denom)
    Lambda0 = (3101.85 - 100.45 * rho - 1481175.59 / T
               + 185698255.2 / T ** 2)
    log_sigma = (_A_INTERCEPT + _B_INV_T / T
                 + _C_LOG_CHAT * np.log10(chat)
                 + _D_LOG_RHO * np.log10(rho)
                 + 2.0 * np.log10(Lambda0)
                 - np.log10(10.0 * Lambda0 + _KOHLRAUSCH * chat_k))
    return 10.0 ** log_sigma


class NaClConductPan2021:
    """Callable conductivity model for a NaCl ocean of fixed salinity,
    following the SwConduct/MgSO4Conduct interface: construct with
    ``wOcean_ppt``, call with ``(P_MPa, T_K, grid=False)``, get S/m.

    Water density rho(P, T) comes from SeaFreeze's 'water1' EOS (the
    paper's own density source).
    """

    def __init__(self, wOcean_ppt):
        self.w_ppt = float(wOcean_ppt)
        m_kgmol = Constants.m_gmol['NaCl'] / 1e3
        w_frac = self.w_ppt / 1e3
        self.molality_molkg = w_frac / (1.0 - w_frac) / m_kgmol
        self._warned = False
        if self.molality_molkg > PAN2021_MOLAL_MAX + 1e-6:
            log.warning(
                f'NaCl molality {self.molality_molkg:.3f} mol/kg exceeds the '
                f'Pan et al. (2021) fitted maximum {PAN2021_MOLAL_MAX} '
                f'(10 wt%); conductivity is an extrapolation.')

    def _warn_once(self, P_MPa, T_K):
        if self._warned:
            return
        P = np.atleast_1d(P_MPa)
        T = np.atleast_1d(T_K)
        finite = np.isfinite(P) & np.isfinite(T)
        if not np.any(finite):
            return
        if (np.any(P[finite] < PAN2021_P_MPA[0]) or np.any(P[finite] > PAN2021_P_MPA[1])
                or np.any(T[finite] < PAN2021_T_K[0]) or np.any(T[finite] > PAN2021_T_K[1])):
            log.warning(
                f'Pan et al. (2021) NaCl conductivity evaluated outside the '
                f'experimental envelope P {PAN2021_P_MPA} MPa, T {PAN2021_T_K} K '
                f'(requested P [{np.nanmin(P[finite]):.0f}, {np.nanmax(P[finite]):.0f}] MPa, '
                f'T [{np.nanmin(T[finite]):.1f}, {np.nanmax(T[finite]):.1f}] K). '
                f'Values are smooth extrapolations; below ~212 MPa the model '
                f'misses the MD-predicted conductivity maximum.')
            self._warned = True

    def _rho_water_gcm3(self, P_MPa, T_K, grid):
        if grid:
            PT = np.array([np.atleast_1d(P_MPa), np.atleast_1d(T_K)],
                          dtype=object)
            rho = sfz.getProp(PT, 'water1').rho  # kg/m^3, (nP, nT)
        else:
            P = np.atleast_1d(np.asarray(P_MPa, dtype=np.float64))
            T = np.atleast_1d(np.asarray(T_K, dtype=np.float64))
            PT = np.empty(P.size, dtype=object)
            for i in range(P.size):
                PT[i] = (P[i], T[i])
            rho = sfz.getProp(PT, 'water1').rho
        return np.asarray(rho, dtype=np.float64) / 1e3

    def __call__(self, P_MPa, T_K, grid=False):
        self._warn_once(P_MPa, T_K)
        rho_gcm3 = self._rho_water_gcm3(P_MPa, T_K, grid)
        if grid:
            T_eval = np.broadcast_to(
                np.atleast_1d(T_K)[np.newaxis, :], rho_gcm3.shape)
        else:
            T_eval = np.atleast_1d(np.asarray(T_K, dtype=np.float64))
        sigma = Pan2021_sigma_Sm(self.molality_molkg, T_eval, rho_gcm3)
        if not grid and np.isscalar(P_MPa) and np.isscalar(T_K):
            return float(sigma[0])
        return sigma
