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
    interface_gravity_coeff, isostatic_gravity, mass_neutral_shell_density,
    project_powers, rescale_coeff, triaxial_to_H2m, gravity_g, G)

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


# --- mass_neutral_shell_density (rho_ice mass-neutrality fix) -----------
#
# See metadata.blockers_open.NEW_MAJOR_1_rho_ice_mass_neutrality in
# PlanetProfile/Inference/configs/enceladus_cassini_isostasy_7D.json for the
# finding this helper closes: an un-compensated rho_ice_kgm3 EOS override
# breaks the campaign's own per-node mass invariant by 3-7x.

R_T_ENC = 252.22e3
M_ENC_FIDUCIAL = 1.08022e20
RHO_OCEAN_ENC = 1005.0
RHO_ICE_ENC = 925.0


def _ocean_stack(zb_km, D_ocean_km, rho_ice=RHO_ICE_ENC,
                  rho_ocean=RHO_OCEAN_ENC, R_T=R_T_ENC, M=M_ENC_FIDUCIAL):
    """Mass-conserved 3-layer (core, ocean, shell) Enceladus-class stack."""
    R_b = R_T - zb_km * 1e3
    R_c = R_b - D_ocean_km * 1e3
    V = 4.0 * np.pi / 3.0
    m_shell = rho_ice * V * (R_T ** 3 - R_b ** 3)
    m_ocean = rho_ocean * V * (R_b ** 3 - R_c ** 3)
    rho_core = (M - m_shell - m_ocean) / (V * R_c ** 3)
    return np.array([R_c, R_b, R_T]), np.array([rho_core, rho_ocean, rho_ice])


def _frozen_stack(R_c_km, R_T=R_T_ENC, rho_shell=RHO_ICE_ENC,
                  M=M_ENC_FIDUCIAL):
    """Mass-conserved 2-layer (core, shell) frozen stack."""
    R_c = R_c_km * 1e3
    V = 4.0 * np.pi / 3.0
    rho_core = (M - rho_shell * V * (R_T ** 3 - R_c ** 3)) / (V * R_c ** 3)
    return np.array([R_c, R_T]), np.array([rho_core, rho_shell])


def _mass_of(radial, rho):
    inner = np.concatenate(([0.0], radial[:-1]))
    vol = 4.0 / 3.0 * np.pi * (radial ** 3 - inner ** 3)
    return float(np.sum(rho * vol))


@pytest.mark.parametrize('zb_km,D_ocean_km,rho_shell_new', [
    (25.0, 36.0, 915.0),
    (25.0, 36.0, 935.0),
    (30.0, 40.0, 918.0),
    (20.0, 30.0, 933.0),
])
def test_mass_neutral_shell_density_conserves_mass_ocean_branch(
        zb_km, D_ocean_km, rho_shell_new):
    radial, rho = _ocean_stack(zb_km, D_ocean_km)
    rho_out, resid = mass_neutral_shell_density(
        radial, rho, rho_shell_new, M_ENC_FIDUCIAL, shell_idx=2,
        interior_idx=0)
    assert rho_out[2] == pytest.approx(rho_shell_new, rel=1e-12)
    m_new = _mass_of(radial, rho_out)
    assert abs(m_new / M_ENC_FIDUCIAL - 1.0) <= 1e-12
    assert abs(resid) <= 1e-12


@pytest.mark.parametrize('R_c_km,rho_shell_new', [
    (200.0, 915.0),
    (200.0, 935.0),
    (150.0, 918.0),
    (220.0, 933.0),
])
def test_mass_neutral_shell_density_conserves_mass_frozen_branch(
        R_c_km, rho_shell_new):
    radial, rho = _frozen_stack(R_c_km)
    rho_out, resid = mass_neutral_shell_density(
        radial, rho, rho_shell_new, M_ENC_FIDUCIAL, shell_idx=1,
        interior_idx=0)
    assert rho_out[1] == pytest.approx(rho_shell_new, rel=1e-12)
    m_new = _mass_of(radial, rho_out)
    assert abs(m_new / M_ENC_FIDUCIAL - 1.0) <= 1e-12
    assert abs(resid) <= 1e-12


def test_mass_neutral_shell_density_identity_is_noop():
    radial, rho = _ocean_stack(25.0, 36.0)
    rho_out, resid = mass_neutral_shell_density(
        radial, rho, rho[2], M_ENC_FIDUCIAL, shell_idx=2, interior_idx=0)
    assert np.allclose(rho_out, rho, rtol=1e-12)
    assert abs(resid) <= 1e-12


def test_mass_neutral_shell_density_does_not_mutate_input():
    radial, rho = _ocean_stack(25.0, 36.0)
    rho_copy = rho.copy()
    mass_neutral_shell_density(radial, rho, 915.0, M_ENC_FIDUCIAL,
                               shell_idx=2, interior_idx=0)
    assert np.array_equal(rho, rho_copy)


def test_mass_neutral_shell_density_unphysical_raises():
    radial, rho = _frozen_stack(200.0)
    # A shell density this large cannot be compensated by any positive
    # interior density: it would require rho_interior < 0.
    with pytest.raises(ValueError):
        mass_neutral_shell_density(radial, rho, 1.0e7, M_ENC_FIDUCIAL,
                                   shell_idx=1, interior_idx=0)
