"""
Physics regression test for Europa ice shell structure.

Validates computed ice shell thickness against Roberts & Nimmo (2008) and
Anderson et al. (1998) observational constraints. This ensures the hydrosphere
calculations remain physically consistent with published Europa models.

References:
    Roberts, J. H., & Nimmo, F. (2008). Tidal heating and the long-term
    stability of a subsurface ocean on Enceladus. Icarus, 194(2), 675-689.

    Anderson, J. D., et al. (1998). Europa's differentiated internal structure:
    Inferences from four Galileo encounters. Science, 281(5385), 2019-2022.
"""
import numpy as np


def test_europa_ice_shell_thickness():
    """ Test that Europa ice shell thickness falls within observational constraints.

        Anderson et al. (1998) constrain the ice shell + ocean thickness to ~100-150 km
        based on Galileo gravity measurements. Most thermal models (Roberts & Nimmo 2008,
        Hussmann & Spohn 2004) predict ice shells of 15-40 km for reasonable heat fluxes.

        This test validates that PlanetProfile produces physically reasonable values.
    """
    from PlanetProfile.Utilities.defineStructs import PlanetStruct, ParamsStruct
    from PlanetProfile.Main import PlanetProfile
    from PlanetProfile import GetConfig

    # Minimal Europa configuration
    Planet = PlanetStruct('Europa')

    # Standard Europa parameters (Anderson et al. 1998)
    Planet.Bulk.R_m = 1560.8e3  # Mean radius
    Planet.Bulk.M_kg = 4.7998e22  # Mass
    Planet.Bulk.Tsurf_K = 110.0  # Surface temperature
    Planet.Bulk.Psurf_MPa = 0.0  # No atmosphere

    # Ocean parameters (Khurana et al. 1998 inferred salinity)
    Planet.Ocean.comp = 'Seawater'
    Planet.Ocean.wOcean_ppt = 35.0  # Standard seawater
    Planet.Bulk.Tb_K = 270.0  # Ice-ocean boundary temperature

    # Get default configuration
    Params = GetConfig(Planet, '.')
    Params.CALC_NEW = 1
    Params.CALC_NEW_REF = 1

    # Run 1D model
    Planet, Params = PlanetProfile(Planet, Params)

    # Extract ice shell thickness
    ice_thickness_km = Planet.zb_km

    # Observational constraints
    ICE_MIN_KM = 10.0  # Minimum from tidal heating models
    ICE_MAX_KM = 50.0  # Maximum from gravity + thermal models

    assert ice_thickness_km >= ICE_MIN_KM, (
        f'Europa ice shell too thin: {ice_thickness_km:.1f} km < {ICE_MIN_KM} km. '
        f'Unrealistically high heat flux or incorrect thermodynamics.'
    )

    assert ice_thickness_km <= ICE_MAX_KM, (
        f'Europa ice shell too thick: {ice_thickness_km:.1f} km > {ICE_MAX_KM} km. '
        f'Unrealistically low heat flux or missing tidal heating.'
    )

    print(f'Europa ice shell thickness: {ice_thickness_km:.2f} km')
    print(f'Constraints: [{ICE_MIN_KM}, {ICE_MAX_KM}] km')
    print(f'PASS: Within observational bounds')


def test_europa_ocean_thickness():
    """ Test that Europa ocean thickness is consistent with gravity constraints.

        Anderson et al. (1998) constrain the total H2O layer (ice + ocean) to
        ~80-170 km based on MoI = 0.346. This test validates that the combined
        ice + ocean thickness falls in this range.
    """
    from PlanetProfile.Utilities.defineStructs import PlanetStruct, ParamsStruct
    from PlanetProfile.Main import PlanetProfile
    from PlanetProfile import GetConfig

    Planet = PlanetStruct('Europa')
    Planet.Bulk.R_m = 1560.8e3
    Planet.Bulk.M_kg = 4.7998e22
    Planet.Bulk.Tsurf_K = 110.0
    Planet.Ocean.comp = 'Seawater'
    Planet.Ocean.wOcean_ppt = 35.0
    Planet.Bulk.Tb_K = 270.0

    Params = GetConfig(Planet, '.')
    Params.CALC_NEW = 1
    Params.CALC_NEW_REF = 1

    Planet, Params = PlanetProfile(Planet, Params)

    # Extract layer thicknesses
    ice_km = Planet.zb_km
    ocean_km = Planet.Bulk.R_m / 1e3 - Planet.zb_km - Planet.Steps.Rmantle_m / 1e3
    h2o_total_km = ice_km + ocean_km

    # Anderson et al. (1998) constraints
    H2O_MIN_KM = 80.0
    H2O_MAX_KM = 170.0

    assert h2o_total_km >= H2O_MIN_KM, (
        f'Total H2O layer too thin: {h2o_total_km:.1f} km < {H2O_MIN_KM} km. '
        f'Inconsistent with gravity data (MoI = 0.346 ± 0.005).'
    )

    assert h2o_total_km <= H2O_MAX_KM, (
        f'Total H2O layer too thick: {h2o_total_km:.1f} km > {H2O_MAX_KM} km. '
        f'Inconsistent with gravity data (MoI = 0.346 ± 0.005).'
    )

    print(f'Europa H2O layer: ice={ice_km:.1f} km + ocean={ocean_km:.1f} km = {h2o_total_km:.1f} km')
    print(f'Constraints: [{H2O_MIN_KM}, {H2O_MAX_KM}] km')
    print(f'PASS: Within gravity constraints')


if __name__ == '__main__':
    print('=== Europa Structure Regression Tests ===\n')

    print('Test 1: Ice shell thickness')
    test_europa_ice_shell_thickness()
    print()

    print('Test 2: Ocean thickness')
    test_europa_ocean_thickness()
    print()

    print('All tests PASSED!')
