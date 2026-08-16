"""Standing density-contrast invariant for rigid ocean libration torque.

The exact telescoping form is copied from the user-ratified B2' ruling in
``validation_reports/enceladus_isostasy/
b2prime_ADJUDICATED_drho_weighting.json``.  It guards against perturbing a
moment difference instead of the density jump carried by an interface.

This invariant uses the exact interface pressure term ``rho_ocean*f_base``.
It intentionally does not test the separately documented linearized
``Bsp_Asp`` coefficient in the production libration branch.


Scope note (reviewer sign-off 2026-08-15): this invariant is
verified on the EXACT per-interface form. The production libration
path meets it only to ~80% because Bsp_Asp is the linearized
8pi/15 figure moment (1.8% low on rho_ocean*f_base = ~23% of
Delta_rho). Do NOT read a pass here as certifying the production
coefficients; that residual is a recorded freeze-budget term
(config B2prime_signoff_2026_08_15.condition_1).
"""

import numpy as np

from PlanetProfile.Gravity.Librations import (
    compute_eccentricities,
    moi_ellps,
    radial2abc,
)


VOLUME_FACTOR = 4.0 * np.pi / 3.0


def _mass_conserved_stack(
    radius_m,
    mass_kg,
    shell_thickness_m,
    ocean_thickness_m,
    rho_ice_kgm3,
    rho_ocean_kgm3,
):
    """Return a three-layer stack whose core density closes total mass."""
    shell_base_m = radius_m - shell_thickness_m
    core_radius_m = shell_base_m - ocean_thickness_m
    shell_mass_kg = (
        rho_ice_kgm3
        * VOLUME_FACTOR
        * (radius_m**3 - shell_base_m**3)
    )
    ocean_mass_kg = (
        rho_ocean_kgm3
        * VOLUME_FACTOR
        * (shell_base_m**3 - core_radius_m**3)
    )
    rho_core_kgm3 = (
        mass_kg - shell_mass_kg - ocean_mass_kg
    ) / (VOLUME_FACTOR * core_radius_m**3)

    radial_m = np.array([core_radius_m, shell_base_m, radius_m])
    rho_kgm3 = np.array([rho_core_kgm3, rho_ocean_kgm3, rho_ice_kgm3])

    inner_radius_m = np.concatenate(([0.0], radial_m[:-1]))
    closed_mass_kg = np.sum(
        rho_kgm3
        * VOLUME_FACTOR
        * (radial_m**3 - inner_radius_m**3)
    )
    np.testing.assert_allclose(closed_mass_kg, mass_kg, rtol=2e-15)
    return radial_m, rho_kgm3


def _assert_density_contrast_identity(radial_m, rho_kgm3, omega_radps):
    """Assert the B2' exact shell-torque decomposition."""
    polar_ecc, equatorial_ecc = compute_eccentricities(
        radial_m, rho_kgm3, omega_radps
    )
    semi_a_m, semi_b_m, semi_c_m = radial2abc(
        radial_m, polar_ecc, equatorial_ecc
    )
    moment_a, moment_b, _ = moi_ellps(
        semi_a_m, semi_b_m, semi_c_m, rho_kgm3
    )

    def interface_factor(index):
        return (
            4.0
            * np.pi
            / 15.0
            * semi_a_m[index]
            * semi_b_m[index]
            * semi_c_m[index]
            * (semi_a_m[index] ** 2 - semi_b_m[index] ** 2)
        )

    f_top = interface_factor(2)
    f_base = interface_factor(1)
    rho_ocean = rho_kgm3[1]
    rho_ice = rho_kgm3[2]

    shell_moment_difference = moment_b[2] - moment_a[2]
    ks_over_3omega2 = shell_moment_difference + rho_ocean * f_base
    density_contrast_form = (
        rho_ice * f_top + (rho_ocean - rho_ice) * f_base
    )
    np.testing.assert_allclose(
        ks_over_3omega2,
        density_contrast_form,
        rtol=1e-12,
        atol=0.0,
    )


def test_density_contrast_identity_across_enceladus_sweep():
    """The invariant holds for 27 mass-conserved Enceladus stacks."""
    for shell_thickness_km in np.arange(5.0, 46.0, 5.0):
        for ocean_thickness_km in (20.0, 36.0, 50.0):
            radial_m, rho_kgm3 = _mass_conserved_stack(
                radius_m=252.22e3,
                mass_kg=1.08022e20,
                shell_thickness_m=shell_thickness_km * 1e3,
                ocean_thickness_m=ocean_thickness_km * 1e3,
                rho_ice_kgm3=925.0,
                rho_ocean_kgm3=1005.0,
            )
            _assert_density_contrast_identity(
                radial_m, rho_kgm3, omega_radps=5.307e-5
            )


def test_density_contrast_identity_for_titan_class_stack():
    """The invariant is not accidentally Enceladus-scale-specific."""
    radial_m, rho_kgm3 = _mass_conserved_stack(
        radius_m=2574.73e3,
        mass_kg=1.3452e23,
        shell_thickness_m=100.0e3,
        ocean_thickness_m=200.0e3,
        rho_ice_kgm3=930.0,
        rho_ocean_kgm3=1100.0,
    )
    titan_omega_radps = 2.0 * np.pi / (15.945 * 86400.0)
    _assert_density_contrast_identity(
        radial_m, rho_kgm3, omega_radps=titan_omega_radps
    )
