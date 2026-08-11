"""Unit tests for PlanetProfile.Gravity.isostasy (Hemingway & Mittal 2019
forward model; spec plans/active/enceladus-isostasy-module-spec.md)."""
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Gravity.isostasy import (  # noqa: E402
    airy_root, eccentricities_to_axes, finite_amplitude_coeffs,
    interface_gravity_coeff, isostatic_gravity, project_powers,
    rescale_coeff, triaxial_to_H2m, gravity_g, G)

# Enceladus constants (Archinal mean radius; mass per program table)
R_ENC = 252.1e3
M_ENC = 1.08022e20
RHO_ICE = 925.0
RHO_OCEAN = 1020.0


def test_triaxial_to_H2m_enceladus_class():
    # Nimmo et al. 2011-class axes: a=256.6, b=251.4, c=248.3 km
    H20, H22 = triaxial_to_H2m(256.6e3, 251.4e3, 248.3e3)
    # published H&M Table 1: H20 = -3846 m, H22 = 917 m (Nimmo)
    assert H20 == pytest.approx(-3766.7, rel=0.05)   # (2/3)(c-(a+b)/2)
    assert H22 == pytest.approx(866.7, rel=0.01)
    # class agreement with the published values (conventions match to ~6%)
    assert abs(H20 - (-3846)) < 250
    assert abs(H22 - 917) < 80


def test_eccentricities_axes_roundtrip():
    a, b, c = eccentricities_to_axes(
        R_ENC, ep=np.sqrt(1 - (248.3 / 256.6) ** 2),
        eq=np.sqrt(1 - (251.4 / 256.6) ** 2))
    assert (a * b * c) ** (1 / 3) == pytest.approx(R_ENC, rel=1e-12)
    assert a > b > c
    assert (a - b) / (256.6e3 - 251.4e3) == pytest.approx(1.0, rel=0.02)


def test_projection_roundtrip_and_powers():
    H = {(2, 0): -3510.0, (2, 2): 857.0, (3, 0): 420.0}
    powers = project_powers(H, 2)
    for lm, v in H.items():
        assert powers[1][lm] == pytest.approx(v, rel=1e-10), lm
    # H^2 of a pure Y22 field has no Y22 power (cos^2 has no cos term)
    p2 = project_powers({(2, 2): 1000.0}, 2)
    assert abs(p2[2][(2, 2)]) < 1e-6


def test_finite_amplitude_small_amplitude_limit():
    # tiny topography: adjusted coefficients -> raw coefficients
    H = {(2, 0): -3.5, (2, 2): 0.9, (3, 0): 0.4}   # metres on R=252 km
    adj = finite_amplitude_coeffs(H, R_ENC)
    for lm, v in H.items():
        assert adj[lm] == pytest.approx(v, rel=1e-4), lm


def test_finite_amplitude_percent_level_for_root_scale():
    # basal-root-scale topography (km-class on a 227 km interface) gets a
    # correction of a few percent (H&M: ~5% on J2 when omitted end-to-end)
    H = {(2, 0): -8400.0, (2, 2): 2100.0}
    adj = finite_amplitude_coeffs(H, 227.1e3)
    rel20 = abs(adj[(2, 0)] / H[(2, 0)] - 1.0)
    assert 0.002 < rel20 < 0.12, rel20


def test_airy_bracket_matches_reviewer_values():
    # Reviewer-computed cancellation bracket f(D) = R_t^4 - (g_t/g_b) R_b^4
    # at l=2, R_t = 252.1 km (units km^4): 20 km -> 1.2245e9;
    # 25 km -> 1.4859e9; 30 km -> 1.7314e9, with self-consistent g ratio.
    for D_km, f_expect, g_ratio_expect in ((20.0, 1.2245e9, 0.9699),
                                           (25.0, 1.4859e9, 0.9599),
                                           (30.0, 1.7314e9, 0.9484)):
        R_t = R_ENC
        R_b = R_ENC - D_km * 1e3
        M_below = M_ENC - RHO_ICE * (4 * np.pi / 3) * (R_t**3 - R_b**3)
        g_t = gravity_g(M_ENC, R_t)
        g_b = gravity_g(M_below, R_b)
        assert g_t / g_b == pytest.approx(g_ratio_expect, rel=2e-3)
        f = (R_t / 1e3) ** 4 - (g_t / g_b) * (R_b / 1e3) ** 4
        assert f == pytest.approx(f_expect, rel=5e-3), D_km


def test_airy_root_sign_and_amplification():
    H_t = {(2, 0): -600.0}
    root = airy_root(H_t, RHO_ICE, RHO_OCEAN, g_t=0.1134, g_b=0.109,
                     C2=1.0)
    # opposite sign, amplified by (rho_ice/drho)*(g_t/g_b) ~ 9.7*1.04
    assert root[(2, 0)] > 0
    amp = abs(root[(2, 0)] / H_t[(2, 0)])
    assert amp == pytest.approx((925.0 / 95.0) * (0.1134 / 0.109),
                                rel=1e-12)


def test_isostatic_gravity_cancellation_structure():
    # Full Airy (C2=1) must strongly reduce |C20_nh| vs rigid (C2=0 /
    # frozen), by the bracket ratio f(D)/R_t^4 ~ 0.368 at D=25 km.
    H_obs = {(2, 0): -3510.0, (2, 2): 857.0, (3, 0): 420.0}
    H_hyd = {(2, 0): -2900.0, (2, 2): 800.0}
    kw = dict(H_obs_lm=H_obs, H_hyd_t_lm=H_hyd, R_t=R_ENC,
              rho_ice=RHO_ICE, M_body=M_ENC, R_ref=252.22e3,
              finite_amplitude=False)
    rigid = isostatic_gravity(R_b=None, rho_ocean=None, **kw)
    airy = isostatic_gravity(R_b=R_ENC - 25e3, rho_ocean=RHO_OCEAN,
                             C2=1.0, **kw)
    ratio = airy[(2, 0)] / rigid[(2, 0)]
    assert ratio == pytest.approx(0.368, abs=0.03)
    # same sign as the (negative) non-hydrostatic surface topography
    assert rigid[(2, 0)] < 0 and airy[(2, 0)] < 0
    # C2=0 through the ocean branch == frozen branch exactly
    airy_c0 = isostatic_gravity(R_b=R_ENC - 25e3, rho_ocean=RHO_OCEAN,
                                C2=0.0, **kw)
    for lm in rigid:
        assert airy_c0[lm] == pytest.approx(rigid[lm], rel=1e-12)


def test_isostatic_gravity_magnitude_reviewer_ballpark():
    # Reviewer ballpark: H20_nh ~ -866 m at D=25 km -> C20_nh ~ -4.4e-4
    H_obs = {(2, 0): -3846.0}
    H_hyd = {(2, 0): -2980.0}
    out = isostatic_gravity(
        H_obs_lm=H_obs, H_hyd_t_lm=H_hyd, R_t=R_ENC,
        R_b=R_ENC - 25e3, rho_ice=RHO_ICE, rho_ocean=RHO_OCEAN,
        M_body=M_ENC, R_ref=252.1e3, C2=1.0, finite_amplitude=False)
    assert out[(2, 0)] == pytest.approx(-4.4e-4, rel=0.25)


def test_c30_is_pure_surface_signal_when_uncompensated():
    H_obs = {(3, 0): 420.0}
    out = isostatic_gravity(
        H_obs_lm=H_obs, H_hyd_t_lm={}, R_t=R_ENC, R_b=None,
        rho_ice=RHO_ICE, rho_ocean=None, M_body=M_ENC,
        R_ref=252.22e3, finite_amplitude=False)
    expect = interface_gravity_coeff(420.0, RHO_ICE, R_ENC, 3, M_ENC,
                                     252.22e3)
    assert out[(3, 0)] == pytest.approx(expect, rel=1e-12)
    assert out[(3, 0)] > 0


def test_rescale_coeff_b14():
    # l=2 factor (256.6/252.22)^2 = 1.0351 (3.5% on C22)
    c = rescale_coeff(1.0, 2, 256.6, 252.22)
    assert c == pytest.approx(1.0351, abs=2e-4)
    # round trip
    assert rescale_coeff(c, 2, 252.22, 256.6) == pytest.approx(1.0)


def test_airy_requires_density_contrast():
    with pytest.raises(ValueError, match='exceed'):
        airy_root({(2, 0): 1.0}, 1020.0, 925.0, 0.11, 0.11)
