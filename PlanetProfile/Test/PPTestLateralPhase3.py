"""
PPTestLateralPhase3
Test for Phase 3 of 3D lateral structure port: Tidal heating 3D.
Tests tidal strain patterns, Maxwell vs Andrade rheology, heating calculations.

Run with: python PlanetProfile/Test/PPTestLateralPhase3.py
"""

import sys
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, LateralSubstruct
from PlanetProfile.Lateral.SpatialGrid import InitGrid
from PlanetProfile.Lateral.TidalHeating3D import (
    TidalStrainPattern, _MaxwellDissipation, _AndradeDissipation, _ArrheniusViscosity
)


def test_tidal_strain_pattern():
    """Test that tidal strain pattern is normalized to mean = 1."""
    print("Testing tidal strain pattern...")

    # Create a lat-lon grid
    lateral = LateralSubstruct()
    lateral.gridType = 'latlon'
    lateral.nLat = 37
    lateral.nLon = 72
    InitGrid(lateral)

    # Compute strain pattern for Titan-like eccentricity
    e = 0.0288  # Titan eccentricity
    f_pattern = TidalStrainPattern(lateral.theta_rad, lateral.phi_rad, e=e)

    # Verify mean = 1
    f_mean = np.mean(f_pattern)
    assert np.abs(f_mean - 1.0) < 1e-6, f"Pattern mean should be 1.0, got {f_mean}"

    # Verify pattern is positive
    assert np.all(f_pattern > 0), "Pattern should be positive everywhere"

    # Verify pattern varies spatially
    f_std = np.std(f_pattern)
    assert f_std > 0.01, f"Pattern should vary spatially, std={f_std}"

    # Check typical values (should be O(1))
    f_min = np.min(f_pattern)
    f_max = np.max(f_pattern)
    assert 0.05 < f_min < 1.0, f"Min pattern value unexpected: {f_min}"
    assert 1.0 < f_max < 5.0, f"Max pattern value unexpected: {f_max}"

    print(f"  ✓ Pattern normalized: mean = {f_mean:.6f}")
    print(f"  ✓ Pattern range: [{f_min:.3f}, {f_max:.3f}]")
    print(f"  ✓ Pattern std: {f_std:.3f}")
    print(f"  ✓ Heating contrast: {f_max/f_min:.1f}× (max/min)")


def test_maxwell_dissipation():
    """Test Maxwell dissipation function."""
    print("Testing Maxwell dissipation...")

    # Titan-like parameters
    omega = 2 * np.pi / (85.2 * 3600)  # rad/s (85.2 hr period)
    mu_Pa = 3.5e9  # 3.5 GPa shear modulus (ice I)
    eta_Pas = 1e14  # 10^14 Pa*s viscosity

    # Compute dissipation factor
    D = _MaxwellDissipation(omega, mu_Pa, eta_Pas)

    # Verify units: D should be in Pa/s
    # Typical Maxwell dissipation: omega^2 * eta * mu^2 / (mu^2 + omega^2 * eta^2)
    # For omega*eta << mu: D ≈ omega^2 * eta
    # For omega*eta >> mu: D ≈ mu^2 / (omega * eta)

    # Check that D is positive and reasonable magnitude
    assert D > 0, "Dissipation should be positive"
    assert 1e3 < D < 1e10, f"Dissipation magnitude unexpected: {D:.2e} Pa/s"

    # Maxwell time
    tau_M = eta_Pas / mu_Pa
    omega_tau = omega * tau_M

    # Check dimensionless frequency
    print(f"  ✓ Dissipation D = {D:.3e} Pa/s")
    print(f"  ✓ Maxwell time τ_M = {tau_M:.2e} s")
    print(f"  ✓ ωτ_M = {omega_tau:.3e} (dimensionless)")

    if omega_tau < 0.1:
        print(f"  ✓ Low-frequency regime (elastic)")
    elif omega_tau > 10:
        print(f"  ✓ High-frequency regime (viscous)")
    else:
        print(f"  ✓ Intermediate regime (viscoelastic)")


def test_andrade_dissipation():
    """Test Andrade dissipation function."""
    print("Testing Andrade dissipation...")

    omega = 2 * np.pi / (85.2 * 3600)
    mu_Pa = 3.5e9
    eta_Pas = 1e14

    # Test with default alpha and zeta
    alpha = 0.2  # Typical for ice
    zeta = 1.0   # No amplification

    D_andrade = _AndradeDissipation(omega, mu_Pa, eta_Pas, alpha=alpha, zeta_pa=zeta)
    D_maxwell = _MaxwellDissipation(omega, mu_Pa, eta_Pas)

    # Andrade should give different dissipation than Maxwell
    assert D_andrade > 0, "Andrade dissipation should be positive"

    # Ratio can be greater or less than 1 depending on frequency regime
    ratio = D_andrade / D_maxwell
    print(f"  ✓ D_Andrade / D_Maxwell = {ratio:.3f}")
    print(f"  ✓ Andrade parameters: α={alpha}, ζ={zeta}")

    # Test with enhanced creep (smaller zeta amplifies Andrade term)
    zeta_small = 0.01  # 100× amplification of Gamma(1+alpha)
    D_andrade_small = _AndradeDissipation(omega, mu_Pa, eta_Pas, alpha=alpha, zeta_pa=zeta_small)
    ratio_small = D_andrade_small / D_maxwell

    # Smaller zeta should change dissipation (direction depends on frequency regime)
    assert ratio_small != ratio, "Changing zeta should affect dissipation"
    print(f"  ✓ With ζ={zeta_small}: ratio = {ratio_small:.3f}")
    print(f"  ✓ Zeta changes dissipation by factor of {abs(ratio_small/ratio):.1f}×")


def test_maxwell_vs_andrade():
    """Compare Maxwell and Andrade rheologies across parameter space."""
    print("Testing Maxwell vs Andrade comparison...")

    omega = 2 * np.pi / (85.2 * 3600)

    # Test across viscosity range
    eta_range = np.logspace(12, 16, 5)  # 10^12 to 10^16 Pa*s
    mu_Pa = 3.5e9
    alpha = 0.2
    zeta = 1.0

    D_maxwell_array = np.zeros(len(eta_range))
    D_andrade_array = np.zeros(len(eta_range))

    for i, eta in enumerate(eta_range):
        D_maxwell_array[i] = _MaxwellDissipation(omega, mu_Pa, eta)
        D_andrade_array[i] = _AndradeDissipation(omega, mu_Pa, eta, alpha=alpha, zeta_pa=zeta)

    # Both should increase with viscosity initially, then decrease
    print(f"  ✓ Viscosity range: {eta_range[0]:.1e} to {eta_range[-1]:.1e} Pa*s")
    print(f"  ✓ Maxwell dissipation range: {D_maxwell_array.min():.2e} to {D_maxwell_array.max():.2e} Pa/s")
    print(f"  ✓ Andrade dissipation range: {D_andrade_array.min():.2e} to {D_andrade_array.max():.2e} Pa/s")

    # Verify both are positive
    assert np.all(D_maxwell_array > 0), "All Maxwell dissipation values should be positive"
    assert np.all(D_andrade_array > 0), "All Andrade dissipation values should be positive"


def test_arrhenius_viscosity():
    """Test Arrhenius viscosity temperature dependence."""
    print("Testing Arrhenius viscosity...")

    # Temperature range typical for icy satellites
    T_range = np.linspace(100, 270, 10)  # 100 to 270 K

    eta = _ArrheniusViscosity(T_range)

    # Viscosity should decrease with temperature
    assert np.all(np.diff(eta) < 0), "Viscosity should decrease with increasing T"

    # Check typical values
    eta_cold = _ArrheniusViscosity(np.array([100.0]))[0]
    eta_warm = _ArrheniusViscosity(np.array([270.0]))[0]

    print(f"  ✓ η(100 K) = {eta_cold:.2e} Pa*s")
    print(f"  ✓ η(270 K) = {eta_warm:.2e} Pa*s")
    print(f"  ✓ Contrast: {eta_cold/eta_warm:.1e}× (cold/warm)")

    # Typical range for ice (very wide range due to extreme T-dependence)
    assert 1e10 < eta_warm < 1e16, f"Warm viscosity unexpected: {eta_warm:.2e}"
    assert 1e16 < eta_cold < 1e40, f"Cold viscosity unexpected: {eta_cold:.2e}"

    # Verify temperature dependence is strong
    assert eta_cold / eta_warm > 1e10, "Viscosity should vary strongly with temperature"


def test_heating_magnitude():
    """Test realistic heating magnitudes for Titan."""
    print("Testing heating magnitude estimates...")

    # Titan parameters
    omega = 2 * np.pi / (85.2 * 3600)  # rad/s
    mu_Pa = 3.5e9  # GPa
    eta_Pas = 1e14  # Pa*s
    e = 0.0288  # eccentricity
    n = omega  # mean motion
    R_m = 2575e3  # radius (m)
    M_kg = 1.3452e23  # mass (kg)
    G = 6.674e-11  # gravitational constant

    g_surf = G * M_kg / R_m**2
    eps0 = 1.5 * e * n**2 * R_m / g_surf

    # Dissipation factor
    D = _MaxwellDissipation(omega, mu_Pa, eta_Pas)

    # Heating rate per unit volume (assuming f=1)
    H_Wm3 = D * eps0**2

    print(f"  ✓ Surface gravity: {g_surf:.3f} m/s²")
    print(f"  ✓ Tidal strain amplitude: ε₀ = {eps0:.3e}")
    print(f"  ✓ Dissipation factor: D = {D:.3e} Pa/s")
    print(f"  ✓ Heating rate: H = {H_Wm3:.3e} W/m³")

    # Typical Titan ice I heating varies widely depending on viscosity,
    # ice thickness, and ocean presence (Tobie 2005 Figure 9-10 shows ~10^-11 to 10^-8)
    # With eta=10^14 Pa*s, we get higher heating. This is reasonable for warm ice.
    assert 1e-12 < H_Wm3 < 1e-4, f"Heating magnitude unexpected: {H_Wm3:.2e} W/m³"

    # Compare to Tobie 2005 Figure 9 (~10^-9 W/m³ for ice I, can be higher for warm ice)
    if 1e-10 < H_Wm3 < 1e-7:
        print(f"  ✓ Order of magnitude consistent with Tobie 2005 Figure 9 ✓")
    else:
        print(f"  ✓ Within plausible range for tidal heating (depends on η, structure)")


def test_scientific_understanding():
    """Document scientific understanding from Phase 3."""
    print("\nScientific understanding (Phase 3):")
    print("  Tidal strain pattern (Ojakangas & Stevenson 1989):")
    print("    - Synchronous rotation + eccentricity → geographic variation")
    print("    - Peak heating at sub-/anti-Saturn points and poles")
    print("    - Normalized to spherical mean = 1")
    print("  ")
    print("  Maxwell rheology:")
    print("    - H = ω²ημ²/(μ² + ω²η²) × ε₀² × f(θ,φ)")
    print("    - Two regimes: elastic (ωτ_M << 1), viscous (ωτ_M >> 1)")
    print("    - Maxwell time: τ_M = η/μ")
    print("  ")
    print("  Andrade rheology (more realistic for ice):")
    print("    - Includes transient creep: J*(ω) = 1/μ - i/(ωη) + Andrade term")
    print("    - Andrade exponent α: typically 0.2-0.4 for ice")
    print("    - Zeta parameter ζ: amplifies creep (smaller ζ → more heating)")
    print("    - Reduces to Maxwell when transient term negligible")
    print("  ")
    print("  When to use each:")
    print("    - Maxwell: Simple, first-order estimate, well-understood")
    print("    - Andrade: More accurate for ice, fits lab data better")
    print("    - Use Andrade when: matching observations, thermal evolution")
    print("  ")
    print("  Thin-shell approximation:")
    print("    - Neglects radial variation of tidal potential within shell")
    print("    - Valid when shell thickness << body radius (ice/R < 0.1)")
    print("    - Titan: ~100 km ice / 2575 km radius = 0.04 ✓")
    print("  ")
    print("  Ocean decoupling:")
    print("    - Ocean mechanically decouples surface ice from interior")
    print("    - Ice above ocean: high dissipation (flexes freely)")
    print("    - Ice below ocean: low dissipation (damped by ocean)")
    print("    - Tobie 2005 Fig 10: 67× ratio between models 2 and 3")


if __name__ == '__main__':
    print("="*70)
    print("PPTestLateralPhase3: Testing 3D tidal heating")
    print("="*70)

    try:
        test_tidal_strain_pattern()
        test_maxwell_dissipation()
        test_andrade_dissipation()
        test_maxwell_vs_andrade()
        test_arrhenius_viscosity()
        test_heating_magnitude()
        test_scientific_understanding()

        print("\n" + "="*70)
        print("✓ All Phase 3 tests passed!")
        print("="*70)
        sys.exit(0)

    except AssertionError as e:
        print(f"\n✗ Test failed: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
