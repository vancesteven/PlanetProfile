"""Regression tests for the explicit-crust Moho/material-boundary handling.

These guard the properties the lunar structural work depends on: the Moho must
land exactly on a structural node, every silicate shell must belong to exactly
one material, and a grid that straddles the Moho must fail loudly rather than
silently producing a mixed-material shell.

Unit-level by design -- they exercise the boundary functions directly and need
no Perple_X table or full PlanetProfile run.
"""
import numpy as np
import pytest

from PlanetProfile.Thermodynamics.Geophysical import (
    BuildMohoAwareSilicateGrid,
    GetSilicateMaterialMask,
    MohoBoundaryTolerance_m,
)
from PlanetProfile.Utilities.defineStructs import Constants


# Lunar-like geometry: 1737.4 km surface, 40 km crust -> Moho at 1697.4 km.
R_SURF_M = 1737.4e3
R_MOHO_M = R_SURF_M - 40.0e3
R_INNER_M = 350.0e3


def _grid(nShells=20, rMoho_m=R_MOHO_M):
    return BuildMohoAwareSilicateGrid(R_SURF_M, R_INNER_M, nShells, rMoho_m)


# --- 1. Moho lies exactly on a structural node ---------------------------

def test_moho_lands_exactly_on_a_node():
    grid = _grid()
    tol = MohoBoundaryTolerance_m(R_SURF_M, R_INNER_M, R_MOHO_M)
    assert np.min(np.abs(grid - R_MOHO_M)) <= tol, 'Moho is not on a grid node'


def test_moho_node_is_unique():
    grid = _grid()
    tol = MohoBoundaryTolerance_m(R_SURF_M, R_INNER_M, R_MOHO_M)
    assert np.count_nonzero(np.abs(grid - R_MOHO_M) <= tol) == 1


def test_insertion_preserves_legacy_boundaries_and_adds_exactly_one():
    nShells = 20
    legacy = np.linspace(R_SURF_M, R_INNER_M, nShells + 1)
    grid = _grid(nShells)
    assert grid.size == legacy.size + 1, 'expected exactly one inserted boundary'
    # Every legacy boundary must survive bit-for-bit.
    for value in legacy:
        assert np.any(grid == value), 'a legacy boundary was moved by insertion'


def test_already_aligned_moho_leaves_grid_unchanged():
    nShells = 20
    legacy = np.linspace(R_SURF_M, R_INNER_M, nShells + 1)
    aligned_moho = legacy[5]
    grid = BuildMohoAwareSilicateGrid(R_SURF_M, R_INNER_M, nShells, aligned_moho)
    assert np.array_equal(grid, legacy), 'aligned Moho should not alter the grid'


# --- 2 & 3. Unique crust/mantle assignment, no mixed shell ---------------

def test_every_shell_is_exactly_one_material():
    grid = _grid()
    mask = GetSilicateMaterialMask(grid[None, :], R_MOHO_M)
    assert mask.shape[1] == grid.size - 1
    valid = {Constants.materialCrust, Constants.materialMantle}
    assert set(np.unique(mask)).issubset(valid)
    # Both materials must be present -- a crust model with no crust (or no
    # mantle) is not a crust model.
    assert Constants.materialCrust in mask
    assert Constants.materialMantle in mask


def test_material_split_matches_moho_radius():
    grid = _grid()
    mask = GetSilicateMaterialMask(grid[None, :], R_MOHO_M)[0]
    rInner = grid[1:]
    # Every crust shell sits above the Moho; every mantle shell below it.
    assert np.all(rInner[mask == Constants.materialCrust] >= R_MOHO_M - 1.0)
    assert np.all(grid[:-1][mask == Constants.materialMantle] <= R_MOHO_M + 1.0)


# --- 4. A straddling grid must fail loudly ------------------------------

def test_straddling_shell_raises():
    # Hand-built grid whose middle shell spans the Moho.
    bad = np.array([[R_SURF_M, R_MOHO_M + 5.0e3, R_MOHO_M - 5.0e3, R_INNER_M]])
    with pytest.raises(ValueError, match='straddles the Moho'):
        GetSilicateMaterialMask(bad, R_MOHO_M)


def test_moho_outside_silicate_range_raises():
    with pytest.raises(ValueError):
        BuildMohoAwareSilicateGrid(R_SURF_M, R_INNER_M, 20, R_SURF_M + 1.0e3)


def test_all_crust_or_all_mantle_raises():
    # Moho below the innermost boundary -> no mantle shell at all.
    grid = np.linspace(R_SURF_M, R_INNER_M, 6)[None, :]
    with pytest.raises(ValueError):
        GetSilicateMaterialMask(grid, R_INNER_M - 10.0e3)


def test_non_monotonic_grid_raises():
    bad = np.array([[R_SURF_M, R_MOHO_M - 5.0e3, R_MOHO_M + 5.0e3, R_INNER_M]])
    with pytest.raises(ValueError):
        GetSilicateMaterialMask(bad, R_MOHO_M)


# --- 5. Crust contributes correctly to shell mass -----------------------

def test_crust_shells_contribute_expected_mass():
    """Shell masses built from the material mask must reproduce the analytic
    spherical-shell mass, and the crust total must match a 40 km shell at the
    crust density."""
    rho_crust = 2890.0   # Garcia M3 @ 40 km, kg/m^3
    rho_mantle = 3300.0
    grid = _grid(nShells=40)
    mask = GetSilicateMaterialMask(grid[None, :], R_MOHO_M)[0]

    rho = np.where(mask == Constants.materialCrust, rho_crust, rho_mantle)
    mShell = rho * 4.0 / 3.0 * np.pi * (grid[:-1] ** 3 - grid[1:] ** 3)

    mCrust = np.sum(mShell[mask == Constants.materialCrust])
    mCrustExact = rho_crust * 4.0 / 3.0 * np.pi * (R_SURF_M ** 3 - R_MOHO_M ** 3)
    assert mCrust == pytest.approx(mCrustExact, rel=1e-12), \
        'crust mass does not match the exact 40 km shell'

    mTot = np.sum(mShell)
    mTotExact = mCrustExact + rho_mantle * 4.0 / 3.0 * np.pi * (R_MOHO_M ** 3 - R_INNER_M ** 3)
    assert mTot == pytest.approx(mTotExact, rel=1e-12), 'total silicate mass not conserved'


# --- 6. Mantle EOS is never queried for a crust shell -------------------

def test_mantle_eos_not_queried_for_crust_shells():
    """Mirror the SilRecursionSolid property-assignment branch and assert the
    mantle EOS is called only with mantle-shell indices."""
    grid = _grid(nShells=30)
    mask = GetSilicateMaterialMask(grid[None, :], R_MOHO_M)

    queried = []

    def fake_mantle_eos(P, T):
        P = np.atleast_1d(P)
        queried.append(np.size(P))
        return np.full(np.size(P), 3300.0)

    nProfiles, nShells = mask.shape
    rho = np.zeros((nProfiles, nShells))
    for j in range(nShells):
        isCrust = mask[:, j] == Constants.materialCrust
        isMantle = np.logical_not(isCrust)
        rho[isCrust, j] = 2890.0
        if np.any(isMantle):
            rho[isMantle, j] = fake_mantle_eos(np.zeros(np.count_nonzero(isMantle)),
                                               np.zeros(np.count_nonzero(isMantle)))

    nCrustCells = int(np.count_nonzero(mask == Constants.materialCrust))
    nMantleCells = int(np.count_nonzero(mask == Constants.materialMantle))
    assert sum(queried) == nMantleCells, 'mantle EOS query count != mantle cell count'
    assert nCrustCells > 0, 'test is vacuous without crust cells'
    # No crust cell may carry the mantle density.
    assert np.all(rho[mask == Constants.materialCrust] == 2890.0)
