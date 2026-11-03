"""
Unit tests for ``PlanetProfile.Inference.structure_derivation``.

Covers Phase C1 Stage 2 mass-conservation helpers:
    - compute_cmr2: exact analytic CMR² for piecewise-constant shells
    - derive_silicate_density: mass-conservation solve with bounds rejection
    - extract_hydrosphere_layers: PP cache → (layers, R_oceanbot, M_hydro)

Tests use hand-computed expected values where possible to keep them robust
against any future formula refactor.
"""
import os
import sys

import numpy as np
import pytest

# Make sure the repo root is importable (consistent with sibling test files)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

from PlanetProfile.Inference.structure_derivation import (
    compute_cmr2,
    derive_silicate_density,
    extract_hydrosphere_layers,
)


# ---------------------------------------------------------------------------
# compute_cmr2
# ---------------------------------------------------------------------------


def test_uniform_sphere_cmr2_is_0p4_exact():
    """C/MR² = 2/5 for a uniform-density sphere — independent of R and ρ."""
    R = 1.0e6
    layers = [(0.0, R, 3000.0)]
    np.testing.assert_allclose(compute_cmr2(layers), 0.4, rtol=1e-12)

    # Same answer with different R, different ρ.
    np.testing.assert_allclose(
        compute_cmr2([(0.0, 7.0e6, 5500.0)]), 0.4, rtol=1e-12
    )


def test_two_shell_cmr2_hand_computed():
    """Two-shell body with hand-computed expected CMR².

    R = 1e6 m, R_core = 0.5 R, ρ_core = 8000, ρ_mantle = 3000.

    Working in normalized units (R = 1, drop the 4π/3 factor that cancels):
        M_core   ∝ (0.5)³ × 8000 = 0.125 × 8000 = 1000
        M_mantle ∝ (1   - 0.5³) × 3000 = 0.875 × 3000 = 2625
        M_total  ∝ 3625

    For inertia (drop the 8π/15 factor that cancels in the ratio with M·R²,
    where M's prefactor is 4π/3):
        I_core   ∝ (0.5)⁵ × 8000 = 0.03125 × 8000 = 250
        I_mantle ∝ (1 - 0.5⁵) × 3000 = 0.96875 × 3000 = 2906.25
        I_total  ∝ 3156.25

    CMR² = I/(M R²) = (8π/15 · I_norm) / ((4π/3 · M_norm) · 1²)
         = (8/15) · I_norm / ((4/3) · M_norm)
         = (8 · 3 · 3156.25) / (15 · 4 · 3625)
    """
    R = 1.0e6
    R_core = 0.5 * R
    layers = [(0.0, R_core, 8000.0), (R_core, R, 3000.0)]
    cmr2 = compute_cmr2(layers)
    expected = (8.0 * 3.0 * 3156.25) / (15.0 * 4.0 * 3625.0)
    np.testing.assert_allclose(cmr2, expected, rtol=1e-12)


def test_compute_cmr2_explicit_M_overrides_integrated():
    """Passing explicit M_body_kg should override the integrated mass.

    Halving M doubles CMR² (linear in 1/M).
    """
    layers = [(0.0, 1.0e6, 3000.0)]
    M_int = (4.0 / 3.0) * np.pi * (1.0e6) ** 3 * 3000.0
    cmr2_default = compute_cmr2(layers)
    cmr2_half_mass = compute_cmr2(layers, M_body_kg=M_int / 2.0)
    np.testing.assert_allclose(cmr2_half_mass, 2.0 * cmr2_default, rtol=1e-12)


def test_compute_cmr2_explicit_R_overrides_max():
    """Passing R_body_m larger than the outer layer should scale CMR² by
    (R_max / R_body)²."""
    layers = [(0.0, 1.0e6, 3000.0)]
    cmr2_default = compute_cmr2(layers)
    cmr2_double_R = compute_cmr2(layers, R_body_m=2.0e6)
    np.testing.assert_allclose(
        cmr2_double_R, cmr2_default / 4.0, rtol=1e-12
    )


def test_compute_cmr2_rejects_bad_layer():
    """Inner ≥ outer should raise."""
    with pytest.raises(ValueError):
        compute_cmr2([(1.0e6, 5.0e5, 3000.0)])
    with pytest.raises(ValueError):
        compute_cmr2([(-1.0, 5.0e5, 3000.0)])


# ---------------------------------------------------------------------------
# derive_silicate_density
# ---------------------------------------------------------------------------


def test_derive_silicate_density_roundtrip():
    """Build a body with known ρ_sil_true; decompose; check ρ_sil recovered."""
    R = 2.0e6
    R_oceanbot = 0.85 * R
    R_core = 0.45 * R
    rho_core = 7500.0
    rho_sil_true = 3100.0
    rho_hydro = 1050.0

    M_hydro = (4.0 / 3.0) * np.pi * (R ** 3 - R_oceanbot ** 3) * rho_hydro
    M_core = (4.0 / 3.0) * np.pi * R_core ** 3 * rho_core
    V_sil = (4.0 / 3.0) * np.pi * (R_oceanbot ** 3 - R_core ** 3)
    M_sil = V_sil * rho_sil_true
    M_total = M_hydro + M_core + M_sil

    rho_sil, accepted = derive_silicate_density(
        M_total_kg=M_total,
        M_hydrosphere_kg=M_hydro,
        R_oceanbot_m=R_oceanbot,
        R_core_m=R_core,
        rho_core_kgm3=rho_core,
        bounds=(2200.0, 3500.0),
    )
    np.testing.assert_allclose(rho_sil, rho_sil_true, rtol=1e-12)
    assert accepted is True


def test_derive_silicate_density_rejects_below_lower_bound():
    """Body too light for plausible silicate → ρ_sil < 2200 → rejected.

    With R=1e6, R_oceanbot=0.9 R, R_core=0.4 R, ρ_core=8000:
        V_sil  ≈ (4π/3)(0.9³ − 0.4³)R³ ≈ 2.785e18 m³
        M_core ≈ (4π/3)(0.4 R)³ × 8000 ≈ 2.14e21 kg
    For ρ_sil to fall below 2200 we need
        M_total − M_hydro − M_core < 2200 × V_sil ≈ 6.13e21 kg.
    Setting M_total = 1.0e22, M_hydro = 5.0e21 gives M_sil ≈ 2.86e21 →
    ρ_sil ≈ 1027 kg/m³, well below the lower bound.
    """
    R = 1.0e6
    R_oceanbot = 0.9 * R
    R_core = 0.4 * R
    rho_sil, accepted = derive_silicate_density(
        M_total_kg=1.0e22,
        M_hydrosphere_kg=5.0e21,
        R_oceanbot_m=R_oceanbot,
        R_core_m=R_core,
        rho_core_kgm3=8000.0,
        bounds=(2200.0, 3500.0),
    )
    assert accepted is False
    assert rho_sil < 2200.0


def test_derive_silicate_density_rejects_above_upper_bound():
    """Body too dense → ρ_sil > 3500 → rejected."""
    R = 1.0e6
    R_oceanbot = 0.9 * R
    R_core = 0.1 * R
    rho_sil, accepted = derive_silicate_density(
        M_total_kg=1.0e23,        # tuned high
        M_hydrosphere_kg=1.0e21,
        R_oceanbot_m=R_oceanbot,
        R_core_m=R_core,
        rho_core_kgm3=8000.0,
        bounds=(2200.0, 3500.0),
    )
    assert accepted is False
    assert rho_sil > 3500.0


def test_derive_silicate_density_no_core_limit():
    """R_core = 0: ρ_sil = (M_total - M_hydro) / V_oceanbot."""
    R = 1.0e6
    R_oceanbot = 0.9 * R
    M_total = 5.0e22
    M_hydro = 1.0e22
    rho_sil, _ = derive_silicate_density(
        M_total_kg=M_total,
        M_hydrosphere_kg=M_hydro,
        R_oceanbot_m=R_oceanbot,
        R_core_m=0.0,
        rho_core_kgm3=7000.0,        # value irrelevant when R_core = 0
        bounds=(2200.0, 3500.0),
    )
    expected = (M_total - M_hydro) / ((4.0 / 3.0) * np.pi * R_oceanbot ** 3)
    np.testing.assert_allclose(rho_sil, expected, rtol=1e-12)


def test_derive_silicate_density_zero_volume_rejects():
    """R_core = R_oceanbot ⇒ V_sil = 0 ⇒ NaN, rejected."""
    R = 1.0e6
    rho_sil, accepted = derive_silicate_density(
        M_total_kg=1.0e23,
        M_hydrosphere_kg=1.0e22,
        R_oceanbot_m=0.5 * R,
        R_core_m=0.5 * R,
        rho_core_kgm3=7000.0,
        bounds=(2200.0, 3500.0),
    )
    assert accepted is False
    assert np.isnan(rho_sil)


def test_derive_silicate_density_inverted_geometry_rejects():
    """R_core > R_oceanbot is non-physical → reject (V_sil < 0)."""
    R = 1.0e6
    rho_sil, accepted = derive_silicate_density(
        M_total_kg=1.0e23,
        M_hydrosphere_kg=1.0e22,
        R_oceanbot_m=0.4 * R,
        R_core_m=0.6 * R,
        rho_core_kgm3=7000.0,
        bounds=(2200.0, 3500.0),
    )
    assert accepted is False


# ---------------------------------------------------------------------------
# extract_hydrosphere_layers
# ---------------------------------------------------------------------------


def test_extract_hydrosphere_simple_ascending_cache():
    """5-cell mock cache: 2 silicate + 3 hydrosphere; ascending r_m."""
    structure = {
        "r_m": np.array([1.0e5, 5.0e5, 8.0e5, 9.0e5, 1.0e6]),
        "rho": np.array([3300.0, 3300.0, 1000.0, 920.0, 920.0]),
        "phases": np.array([50, 50, 0, 1, 1]),
    }
    layers, R_oceanbot, M_hydro = extract_hydrosphere_layers(structure)

    assert len(layers) == 3
    np.testing.assert_allclose(R_oceanbot, 5.0e5)
    # Hand-computed mass:
    expected_M = (4.0 / 3.0) * np.pi * (
        (8.0e5 ** 3 - 5.0e5 ** 3) * 1000.0
        + (9.0e5 ** 3 - 8.0e5 ** 3) * 920.0
        + (1.0e6 ** 3 - 9.0e5 ** 3) * 920.0
    )
    np.testing.assert_allclose(M_hydro, expected_M, rtol=1e-12)
    # Layer geometry:
    assert layers[0] == (5.0e5, 8.0e5, 1000.0)
    assert layers[1] == (8.0e5, 9.0e5, 920.0)
    assert layers[2] == (9.0e5, 1.0e6, 920.0)


def test_extract_hydrosphere_descending_cache_normalized():
    """Cache stored surface-to-center (descending r); helper sorts ascending
    internally and returns the same result as ascending input."""
    asc = {
        "r_m": np.array([1.0e5, 5.0e5, 8.0e5, 9.0e5, 1.0e6]),
        "rho": np.array([3300.0, 3300.0, 1000.0, 920.0, 920.0]),
        "phases": np.array([50, 50, 0, 1, 1]),
    }
    desc = {
        "r_m": asc["r_m"][::-1].copy(),
        "rho": asc["rho"][::-1].copy(),
        "phases": asc["phases"][::-1].copy(),
    }
    asc_out = extract_hydrosphere_layers(asc)
    desc_out = extract_hydrosphere_layers(desc)
    assert asc_out[0] == desc_out[0]
    np.testing.assert_allclose(asc_out[1], desc_out[1])
    np.testing.assert_allclose(asc_out[2], desc_out[2])


def test_extract_hydrosphere_includes_clathrate_and_hp_ices():
    """Phase codes 30 (clathrate) and 6 (Ice VI) are still hydrosphere."""
    structure = {
        "r_m": np.array([2.0e5, 5.0e5, 7.0e5, 9.0e5, 1.0e6]),
        "rho": np.array([3500.0, 3500.0, 1500.0, 1100.0, 920.0]),
        "phases": np.array([50, 50, 6, 30, 1]),
    }
    layers, R_oceanbot, _ = extract_hydrosphere_layers(structure)
    assert len(layers) == 3
    np.testing.assert_allclose(R_oceanbot, 5.0e5)


def test_extract_hydrosphere_excludes_iron_core_phases():
    """Phase 100 (liquid Fe) and 105 (solid Fe) are interior, not hydrosphere."""
    structure = {
        "r_m": np.array([1.0e5, 3.0e5, 5.0e5, 7.0e5, 1.0e6]),
        "rho": np.array([10000.0, 8000.0, 3300.0, 1000.0, 920.0]),
        "phases": np.array([100, 105, 50, 0, 1]),
    }
    layers, R_oceanbot, _ = extract_hydrosphere_layers(structure)
    assert len(layers) == 2
    np.testing.assert_allclose(R_oceanbot, 5.0e5)


def test_extract_hydrosphere_no_water_returns_empty():
    """All-rocky cache → empty layers, zeros, no exception."""
    structure = {
        "r_m": np.array([5.0e5, 1.0e6]),
        "rho": np.array([3500.0, 3500.0]),
        "phases": np.array([50, 50]),
    }
    layers, R_oceanbot, M_hydro = extract_hydrosphere_layers(structure)
    assert layers == []
    assert R_oceanbot == 0.0
    assert M_hydro == 0.0


def test_extract_hydrosphere_rejects_non_contiguous_water():
    """A rocky cell sandwiched between water → ValueError (defensive)."""
    structure = {
        "r_m": np.array([1.0e5, 3.0e5, 5.0e5, 7.0e5, 1.0e6]),
        "rho": np.array([3300.0, 1000.0, 3300.0, 1000.0, 920.0]),
        "phases": np.array([50, 0, 50, 0, 1]),
    }
    with pytest.raises(ValueError, match="not contiguous"):
        extract_hydrosphere_layers(structure)


# ---------------------------------------------------------------------------
# Integration: compute_cmr2 of an assembled v2 body
# ---------------------------------------------------------------------------


def test_assembled_v2_body_cmr2_self_consistent():
    """Build a body, decompose, derive ρ_sil, recompose, recompute CMR² —
    should match a direct CMR² calculation on the original layered profile.

    This is the end-to-end consistency check the runner will perform per
    sample.
    """
    R = 2.5e6  # Ganymede-ish
    R_oceanbot = 0.7 * R
    R_core = 0.4 * R
    rho_core = 7500.0
    rho_sil_true = 3300.0
    # Two-layer hydrosphere: deep HP ice + shallow Ih
    R_hp_top = R_oceanbot + 0.15 * R
    rho_hp = 1280.0
    rho_ih = 920.0

    direct_layers = [
        (0.0, R_core, rho_core),
        (R_core, R_oceanbot, rho_sil_true),
        (R_oceanbot, R_hp_top, rho_hp),
        (R_hp_top, R, rho_ih),
    ]
    # M_total / R from the direct geometry (this is what would be observed)
    M_total = sum(
        (4.0 / 3.0) * np.pi * (ro ** 3 - ri ** 3) * rh
        for ri, ro, rh in direct_layers
    )
    cmr2_direct = compute_cmr2(direct_layers, R_body_m=R, M_body_kg=M_total)

    # Now simulate what the runner does: derive ρ_sil from M_total + cached
    # hydrosphere + sampled core, then assemble layers and recompute CMR².
    M_hydro = sum(
        (4.0 / 3.0) * np.pi * (ro ** 3 - ri ** 3) * rh
        for ri, ro, rh in direct_layers[2:]    # only hydrosphere shells
    )
    rho_sil, accepted = derive_silicate_density(
        M_total_kg=M_total,
        M_hydrosphere_kg=M_hydro,
        R_oceanbot_m=R_oceanbot,
        R_core_m=R_core,
        rho_core_kgm3=rho_core,
        bounds=(2200.0, 3500.0),
    )
    assert accepted is True
    np.testing.assert_allclose(rho_sil, rho_sil_true, rtol=1e-12)

    derived_layers = [
        (0.0, R_core, rho_core),
        (R_core, R_oceanbot, rho_sil),
        (R_oceanbot, R_hp_top, rho_hp),
        (R_hp_top, R, rho_ih),
    ]
    cmr2_derived = compute_cmr2(derived_layers, R_body_m=R, M_body_kg=M_total)
    np.testing.assert_allclose(cmr2_derived, cmr2_direct, rtol=1e-12)
