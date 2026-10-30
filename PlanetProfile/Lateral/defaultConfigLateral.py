""" Default configuration for 3D lateral structure calculations. """

configLateralVersion = 1

def lateralAssign():
    """ Return default lateral structure parameters as a dict.

        These are applied to Planet.Lateral if DO_3D is True.
        Users can override individual settings in their PPBody.py config files.
    """
    return {
        'gridType': 'healpix',   # 'healpix' or 'latlon'
        'nSide': 8,              # HEALPix NSIDE (nPix = 12*nSide^2 = 768 for nSide=8)
        'nLat': 37,              # Latitude points for latlon grid (default 5-deg spacing)
        'nLon': 72,              # Longitude points for latlon grid (default 5-deg spacing)
        'DO_MASS_CONSERVE': True, # Enforce mass conservation after 3D column computation
        'DO_CLATH_LATERAL': False, # Include lateral clathrate variation
        'DO_TIDAL_3D': False,     # Compute 3D tidal heating
    }
