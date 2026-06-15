"""
Physics regression test for 3D tidal heating.

Validates computed tidal heating against Tobie et al. (2005) Fig. 10 ice I
heating values for Titan (~3e-8 W/m^3 for eccentricity tides). This is an
order-of-magnitude sanity check, not a precise reproduction of that specific
model, since the exact viscosity, shell thickness, and rheology may differ.

Reference:
    Tobie, G., Mocquet, A., & Sotin, C. (2005). Tidal dissipation within large
    icy satellites: Applications to Europa and Titan. Icarus, 177(2), 534-549.
"""
import numpy as np


def test_tobie2005_mean_heating():
    """ Test that 3D tidal heating produces order-of-magnitude agreement with
        Tobie et al. (2005) reference values for Titan ice I dissipation.
    """
    from PlanetProfile.Utilities.defineStructs import PlanetStruct, ParamsStruct
    from PlanetProfile.Lateral.TidalHeating3D import ComputeTidalHeating3D

    # Minimal Titan-like configuration
    Planet = PlanetStruct('Titan')
    Planet.Bulk.R_m = 2574.7e3
    Planet.Bulk.M_kg = 1.3452e23
    Planet.Bulk.Tb_K = 255.0
    Planet.Bulk.eccentricity = 0.0288
    Planet.Bulk.meanMotion_radps = 2 * np.pi / (15.945 * 86400)  # Orbital period 15.945 days

    Planet.Ocean.comp = 'MgSO4'
    Planet.Ocean.wOcean_ppt = 100.0

    # HEALPix grid with 192 pixels
    from PlanetProfile.Lateral.SpatialGrid import InitGrid
    Planet.Lateral.gridType = 'healpix'
    Planet.Lateral.nSide = 4
    InitGrid(Planet.Lateral)
    nPix = Planet.Lateral.nPix
    assert nPix == 192, f'Expected 192 pixels for nSide=4, got {nPix}'

    # Uniform 50 km ice shell
    Planet.Lateral.dIce_m = np.full(nPix, 50e3)
    Planet.Lateral.DO_TIDAL_3D = True

    # Excitation period (Titan orbital period)
    Planet.Magnetic.Texc_hr = np.array([15.945 * 24.0])

    # Surface ice layer count (mock for minimal setup)
    Planet.Steps.nSurfIce = 10

    # Compute 3D tidal heating
    Params = ParamsStruct()
    Planet = ComputeTidalHeating3D(Planet, Params, columnPlanets=None)

    # Validate against Tobie 2005 reference
    TOBIE_2005_REF_Wm3 = 3e-8
    TOLERANCE = 0.5  # 50% — loose tolerance for order-of-magnitude validation

    mean_H = np.nanmean(Planet.Lateral.HtidalIce_Wm3)
    residual = abs(mean_H - TOBIE_2005_REF_Wm3) / TOBIE_2005_REF_Wm3

    assert residual < TOLERANCE, (
        f'Tobie 2005 regression: computed {mean_H:.3e} W/m³, '
        f'reference {TOBIE_2005_REF_Wm3:.3e} W/m³, '
        f'residual {residual*100:.1f}% > {TOLERANCE*100:.0f}%'
    )

    print(f'Tobie 2005 regression PASS: computed {mean_H:.3e} W/m³, '
          f'reference {TOBIE_2005_REF_Wm3:.3e} W/m³, residual {residual*100:.1f}%')


if __name__ == '__main__':
    test_tobie2005_mean_heating()
    print('PASS')
