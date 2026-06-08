

# Constants
=============

## T0 and P0

AUTHOR:
Trevor McDougall and Paul Barker                    [ help@teos-10.org ]

VERSION NUMBER: 3.05 (27th January 2015)

REFERENCES:
IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
seawater - 2010: Calculation and use of thermodynamic properties.
Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
See section 2.2, appendix A.2 and Table D.1 of this TEOS-10 Manual.

The software is available from http://www.TEOS-10.org

# Planet.Lateral Fields
=============

3D lateral structure fields for spatially-varying interior models. All fields are None by default unless DO_3D=True.

## Control Flags

- **DO_3D** (bool): Whether to compute 3D laterally-varying structure (default: False)
- **DO_CLATH_LATERAL** (bool): Whether to include lateral clathrate variation (default: False)
- **DO_TIDAL_3D** (bool): Whether to compute 3D tidal heating H(r,θ,φ) (default: False)
- **DO_MASS_CONSERVE** (bool): Whether to enforce mass conservation in 3D calculations (default: True)

## Grid Configuration

- **gridType** (str): Grid type, options: 'healpix' (equal-area) or 'latlon' (regular grid) (default: 'healpix')
- **nSide** (int): HEALPix NSIDE parameter where nPix = 12*nSide² (default: 8, → 768 pixels)
- **nLat** (int): Number of latitude points for latlon grid [dimensionless]
- **nLon** (int): Number of longitude points for latlon grid [dimensionless]
- **theta_rad** (array): Colatitude of each grid point [rad], shape (nPix,)
- **phi_rad** (array): Longitude of each grid point [rad], shape (nPix,)
- **nPix** (int): Total number of grid points [dimensionless]
- **pixArea_sr** (array): Area of each pixel [sr] (steradians), shape (nPix,)

## Ice Thickness Field

- **dIce_m** (array): Ice shell thickness at each grid point [m], shape (nPix,)
- **dIce_Cpq_km** (array): Cosine spherical harmonic coefficients for ice thickness [km], shape (pMax+1, pMax+1)
- **dIce_Spq_km** (array): Sine spherical harmonic coefficients for ice thickness [km], shape (pMax+1, pMax+1)
- **dIce_pMax** (int): Maximum spherical harmonic degree for ice thickness [dimensionless]
- **dIce_func** (callable): Optional function f(theta_rad) returning ice thickness [m]

## Clathrate Fraction Field

- **fClath** (array): Clathrate volume fraction at each grid point [dimensionless, 0-1], shape (nPix,)
- **fClath_Cpq** (array): Cosine spherical harmonic coefficients for clathrate fraction [dimensionless], shape (pMax+1, pMax+1)
- **fClath_Spq** (array): Sine spherical harmonic coefficients for clathrate fraction [dimensionless], shape (pMax+1, pMax+1)
- **fClath_pMax** (int): Maximum spherical harmonic degree for clathrate fraction [dimensionless]

## Tidal Heating

- **Htidal_Wm3** (array): Volumetric tidal heating rate [W/m³], shape (nPix, nRadial)
- **HtidalIce_Wm3** (array): Column-integrated ice tidal heating [W/m³], shape (nPix,)
- **HtidalIceI_top_Wm3** (array): Dissipation at top of ice I layer [W/m³], shape (nPix,)
- **HtidalIceI_bot_Wm3** (array): Dissipation at bottom of ice I layer [W/m³], shape (nPix,)
- **HtidalHP_top_Wm3** (array): Dissipation at top of HP ice below ocean [W/m³], shape (nPix,)
- **HtidalHP_bot_Wm3** (array): Dissipation at bottom of HP ice below ocean [W/m³], shape (nPix,)
- **HtidalIceI_profile_Wm3** (list): List of H(r) arrays in surface ice for each column [W/m³], length nPix
- **HtidalHP_profile_Wm3** (list): List of H(r) arrays in HP ice for each column [W/m³], length nPix
- **rIceI_profile_m** (list): List of radii arrays corresponding to ice I profiles [m], length nPix
- **rHP_profile_m** (list): List of radii arrays corresponding to HP ice profiles [m], length nPix

## Column Summary Fields

- **Tb_K** (array): Basal ice temperature at ice-ocean boundary [K], shape (nPix,)
- **qSurf_Wm2** (array): Surface heat flux [W/m²], shape (nPix,)
- **qBase_Wm2** (array): Basal heat flux at ice-ocean boundary [W/m²], shape (nPix,)
- **kThermEff_WmK** (array): Effective thermal conductivity (mixed clath/ice) [W/m/K], shape (nPix,)
- **sigma_mean_Sm** (array): Mean ocean electrical conductivity [S/m], shape (nPix,)

## Mass Conservation

- **Mtarget_kg** (float): Target total body mass [kg]
- **Mactual_kg** (float): Actual total mass from 3D model [kg]
- **massResidual_frac** (float): Fractional mass residual (Mactual - Mtarget)/Mtarget [dimensionless]

## Spherical Harmonic Convention

All spherical harmonic coefficients use 4π normalization (consistent with MoonMag):
- Y_pq(θ,φ) = P_pq(cos θ) × [C_pq cos(qφ) + S_pq sin(qφ)]
- ∫ Y_pq Y_p'q' dΩ = 4π δ_pp' δ_qq'
- p = degree (angular scale), q = order (longitudinal variation)
- Mean value = C_00 coefficient

REFERENCES:
- Górski et al. (2005): HEALPix equal-area pixelization scheme
- Styczinski et al. (2021): MoonMag 4π SH normalization convention (DOI: 10.1016/j.icarus.2021.114840)
- Tobie et al. (2005): Geographic tidal dissipation applications (DOI: 10.1016/j.icarus.2005.04.006)
- Ojakangas & Stevenson (1989): Tidal heating with eccentricity forcing (DOI: 10.1016/0019-1035(89)90052-3)
