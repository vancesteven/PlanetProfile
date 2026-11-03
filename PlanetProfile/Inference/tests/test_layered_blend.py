"""
Tests for the B-layered blend in forward_models.py.

Covers:
  1. Boundary radii blend linearly in w.
  2. Variable cell counts between s0 and s1 — no shape mismatch.
  3. Across-transition bracket — runner returns nearest neighbour, never blends.
  4. Cell count of blended structure equals s0's cell count.
  5. Continuous fields interpolate within each blended layer (no smearing
     across moving boundaries).
  6. Per-Tb scalar fields linear-blend; Tb-invariant fields copy from s0.
  7. _regions_match is the gate: differing region_phases → not matched.
"""
from __future__ import annotations

import numpy as np
import pytest

from PlanetProfile.Inference.forward_models import (
    _blend_b_layered, _regions_match, _copy_struct,
    apply_bottom_temperature,
)


# -------------------- Synthetic structure helpers --------------------

def _make_two_layer_struct(
    R_body_m: float = 2_400_000.0,
    R_sil_m: float = 2_000_000.0,
    rho_sil: float = 3300.0,
    rho_ice: float = 920.0,
    T_top: float = 100.0,
    T_bot: float = 270.0,
    n_sil: int = 50,
    n_ice: int = 30,
    Tb_K: float = 260.0,
    region_phases=("Sil", "Ih"),
) -> dict:
    """Build a synthetic 2-layer structure (silicate core + ice shell).

    All boundaries are sharp; ρ is layer-constant; T is linear from
    T_top (surface) down to T_bot (silicate top — i.e., the inner side of
    the ice layer). The ice layer occupies r ∈ [R_sil_m, R_body_m].
    """
    r_sil = np.linspace(0.0, R_sil_m, n_sil)
    r_ice = np.linspace(R_sil_m, R_body_m, n_ice + 1)[1:]  # exclude duplicate at R_sil
    r_m = np.concatenate([r_sil, r_ice])
    rho = np.concatenate([np.full(n_sil, rho_sil), np.full(n_ice, rho_ice)])
    phases = np.concatenate([np.full(n_sil, 50, dtype=int), np.full(n_ice, 1, dtype=int)])
    # Linear T profile across each layer (just to give us something to interpolate).
    T_K = np.empty_like(r_m)
    T_K[:n_sil] = np.linspace(T_bot + 1500, T_bot, n_sil)  # silicate hot at center
    T_K[n_sil:] = np.linspace(T_bot, T_top, n_ice)         # ice cool toward surface
    P_MPa = np.linspace(0.0, 4000.0, len(r_m))
    K_Pa = np.where(phases == 1, 1.0e10, 2.5e11)
    mu_Pa = np.where(phases == 1, 3.0e9, 1.0e11)
    eta_Pa_base = np.where(phases == 1, 1.0e14, 1.0e21)
    bulk_visc = np.zeros_like(r_m)

    changeIndices = np.array([0, n_sil, n_sil + n_ice], dtype=int)
    return {
        "r_m": r_m, "rho": rho, "T_K": T_K, "P_MPa": P_MPa,
        "K_Pa": K_Pa, "mu_Pa": mu_Pa, "eta_Pa_base": eta_Pa_base,
        "bulk_visc": bulk_visc, "phases": phases,
        "changeIndices": changeIndices,
        "n_layers": 2,
        "layer_types": ("solid", "solid"),
        "region_phases": list(region_phases),
        "layer_upper_radii": (R_sil_m, R_body_m),
        "Tb_K": Tb_K,
        "CMR2": 0.34,
        "rhoSil_kgm3": rho_sil,
        "D_hsphere_km": (R_body_m - R_sil_m) / 1e3,
        "D_iceIh_km": (R_body_m - R_sil_m) / 1e3,
        "D_iceIII_km": 0.0, "D_iceV_km": 0.0, "D_iceVI_km": 0.0,
        "R_body_m": R_body_m, "Mtot_kg": 1.0e23,
        "omega": 1e-5, "eccentricity": 0.01, "host_mass": 1.9e27, "a_m": 1.9e9,
    }


# -------------------- Tests --------------------

def test_regions_match_basic():
    s0 = _make_two_layer_struct(region_phases=("Sil", "Ih"))
    s1 = _make_two_layer_struct(region_phases=("Sil", "Ih"))
    assert _regions_match(s0, s1)
    s2 = _make_two_layer_struct(region_phases=("Sil", "0", "Ih"))
    s2["n_layers"] = 3  # mismatched
    assert not _regions_match(s0, s2)
    s3 = _make_two_layer_struct(region_phases=("Fe", "Ih"))
    assert not _regions_match(s0, s3)


def test_blend_boundary_linear_in_w():
    """layer_upper_radii (the physical phase interfaces) blend linearly in w."""
    s0 = _make_two_layer_struct(R_sil_m=2_000_000.0, Tb_K=256.0)
    s1 = _make_two_layer_struct(R_sil_m=2_010_000.0, Tb_K=257.0)
    for w in (0.0, 0.25, 0.5, 0.75, 1.0):
        out = _blend_b_layered(s0, s1, w)
        expected_sil_top = (1.0 - w) * 2_000_000.0 + w * 2_010_000.0
        sil_top = float(out["layer_upper_radii"][0])
        assert abs(sil_top - expected_sil_top) < 1.0, (
            f"w={w}: layer_upper_radii[0]={sil_top}, expected={expected_sil_top}"
        )
        # Top of ice layer (= surface) is body-fixed: should be R_body_m for both
        assert abs(float(out["layer_upper_radii"][1]) - 2_400_000.0) < 1.0
        # And the last cell of layer k must equal layer_upper_radii[k]
        ci = np.asarray(out["changeIndices"], dtype=int)
        for k in range(int(out["n_layers"])):
            last_cell_r = float(out["r_m"][ci[k + 1] - 1])
            ur = float(out["layer_upper_radii"][k])
            assert abs(last_cell_r - ur) < 1.0, (
                f"w={w}: layer {k} last-cell r={last_cell_r}, layer_upper={ur}"
            )


def test_blend_handles_variable_cell_counts():
    """The bug we hit: s0 and s1 have different len(r_m). Must not crash."""
    s0 = _make_two_layer_struct(n_sil=50, n_ice=30)   # 80 cells
    s1 = _make_two_layer_struct(n_sil=55, n_ice=35)   # 90 cells
    out = _blend_b_layered(s0, s1, 0.5)
    # Output cell count must equal s0's.
    assert len(out["r_m"]) == 80
    for f in ("rho", "T_K", "P_MPa", "K_Pa", "mu_Pa", "eta_Pa_base"):
        assert len(out[f]) == 80, f"field {f} has wrong length"


def test_blend_no_mush_cells_at_boundary():
    """Within each layer, ρ stays in the layer's expected range — no
    smeared cells at the moving boundary."""
    s0 = _make_two_layer_struct(R_sil_m=2_000_000.0)
    s1 = _make_two_layer_struct(R_sil_m=2_010_000.0)
    out = _blend_b_layered(s0, s1, 0.5)
    rho = out["rho"]
    phases = out["phases"]
    # Silicate cells: ρ ≈ 3300 (no mush near the interface)
    sil_mask = phases == 50
    assert np.all(np.abs(rho[sil_mask] - 3300.0) < 1e-6)
    # Ice cells: ρ ≈ 920
    ice_mask = phases == 1
    assert np.all(np.abs(rho[ice_mask] - 920.0) < 1e-6)


def test_blend_T_continuous_within_layer():
    """T_K interpolates within each layer (continuous), not across."""
    s0 = _make_two_layer_struct(T_top=100.0, T_bot=270.0)
    s1 = _make_two_layer_struct(T_top=110.0, T_bot=275.0)
    out = _blend_b_layered(s0, s1, 0.5)
    T = out["T_K"]
    phases = out["phases"]
    # Ice top (surface) cell: should be ≈ (100 + 110)/2 = 105 K
    surface_T = T[phases == 1][-1]
    assert abs(surface_T - 105.0) < 0.5
    # Ice bottom cell: should be ≈ (270 + 275)/2 = 272.5 K
    ice_bot_T = T[phases == 1][0]
    assert abs(ice_bot_T - 272.5) < 0.5


def test_blend_scalars_linear():
    s0 = _make_two_layer_struct(Tb_K=256.0)
    s0["CMR2"] = 0.340
    s1 = _make_two_layer_struct(Tb_K=257.0)
    s1["CMR2"] = 0.346
    out = _blend_b_layered(s0, s1, 0.5)
    assert abs(out["Tb_K"] - 256.5) < 1e-9
    assert abs(out["CMR2"] - 0.343) < 1e-9


def test_blend_meta_taken_from_s0():
    s0 = _make_two_layer_struct(); s0["R_body_m"] = 2_400_000.0
    s1 = _make_two_layer_struct(); s1["R_body_m"] = 2_400_000.0
    out = _blend_b_layered(s0, s1, 0.5)
    assert out["R_body_m"] == s0["R_body_m"]
    assert out["region_phases"] == s0["region_phases"]
    assert out["n_layers"] == 2


def test_apply_bottom_temperature_blends_when_regions_match():
    """End-to-end: apply_bottom_temperature picks B-layered when regions match."""
    s0 = _make_two_layer_struct(R_sil_m=2_000_000.0, Tb_K=256.0)
    s1 = _make_two_layer_struct(R_sil_m=2_010_000.0, Tb_K=257.0)
    structure_data = {
        "Tb_K_grid": np.array([256.0, 257.0]),
        "structures": [s0, s1],
    }
    out = apply_bottom_temperature({"Tb_K": 256.5}, structure_data)
    # Cell count = s0's
    assert len(out["r_m"]) == len(s0["r_m"])
    # Phase-boundary radius (silicate top) = midpoint of s0 / s1 R_sil_m
    assert abs(float(out["layer_upper_radii"][0]) - 2_005_000.0) < 1.0


def test_apply_bottom_temperature_nearest_across_transition():
    """When bracketing grid points have different region_phases, the
    runner returns the nearest neighbour (no blend)."""
    s0 = _make_two_layer_struct(region_phases=("Sil", "Ih"))
    s0["Tb_K"] = 256.0
    # s1 is qualitatively different: 3 layers (Sil, ocean, Ih).
    s1 = _make_two_layer_struct(region_phases=("Sil", "0", "Ih"))
    s1["Tb_K"] = 257.0
    s1["n_layers"] = 3
    structure_data = {
        "Tb_K_grid": np.array([256.0, 257.0]),
        "structures": [s0, s1],
    }
    # Sample closer to s0
    out = apply_bottom_temperature({"Tb_K": 256.4}, structure_data)
    assert out["region_phases"] == ["Sil", "Ih"]      # took s0
    assert out["n_layers"] == 2
    # Sample closer to s1
    out = apply_bottom_temperature({"Tb_K": 256.7}, structure_data)
    assert out["region_phases"] == ["Sil", "0", "Ih"]
    assert out["n_layers"] == 3


def test_apply_bottom_temperature_endpoints():
    s0 = _make_two_layer_struct(Tb_K=256.0)
    s1 = _make_two_layer_struct(Tb_K=257.0)
    structure_data = {
        "Tb_K_grid": np.array([256.0, 257.0]),
        "structures": [s0, s1],
    }
    # Below low end → s0
    out = apply_bottom_temperature({"Tb_K": 250.0}, structure_data)
    assert abs(out["Tb_K"] - 256.0) < 1e-9
    # Above high end → s1
    out = apply_bottom_temperature({"Tb_K": 260.0}, structure_data)
    assert abs(out["Tb_K"] - 257.0) < 1e-9


def test_blend_w_zero_returns_s0_values():
    """w=0 → blended values match s0 within float tolerance."""
    s0 = _make_two_layer_struct(R_sil_m=2_000_000.0, T_top=100.0, T_bot=270.0)
    s1 = _make_two_layer_struct(R_sil_m=2_010_000.0, T_top=110.0, T_bot=275.0)
    out = _blend_b_layered(s0, s1, 0.0)
    np.testing.assert_allclose(out["r_m"], s0["r_m"], atol=1e-6)
    np.testing.assert_allclose(out["rho"], s0["rho"])
    np.testing.assert_allclose(out["T_K"], s0["T_K"])


def test_blend_w_one_returns_s1_values_resampled_to_s0_grid():
    """w=1 → continuous fields equal s1's values resampled to s0's grid."""
    s0 = _make_two_layer_struct(n_sil=50, n_ice=30, T_top=100.0, T_bot=270.0)
    s1 = _make_two_layer_struct(n_sil=55, n_ice=35, T_top=110.0, T_bot=275.0)
    out = _blend_b_layered(s0, s1, 1.0)
    # T at surface (last ice cell) should match s1's surface T (resampled is exact at endpoints)
    phases = out["phases"]
    assert abs(out["T_K"][phases == 1][-1] - 110.0) < 1e-6
    assert abs(out["T_K"][phases == 1][0]  - 275.0) < 1e-6
    # r_m should match s1's boundaries at endpoints (Tb-invariant since R_body fixed)
    assert abs(out["r_m"][-1] - s0["R_body_m"]) < 1e-6


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
