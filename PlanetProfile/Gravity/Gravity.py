import numpy as np
import logging
import ast
import os
from PlanetProfile.Gravity import pyALMA3Updated
from PlanetProfile.Gravity.Librations import librations
from PlanetProfile.Utilities.Indexing import PhaseConv
from PlanetProfile.Utilities.defineStructs import Constants, Timing
from PlanetProfile.Utilities.ResultsIO import ensure_parent_dir
import time
import mpmath as mp

# Assign logger
log = logging.getLogger('PlanetProfile')

ALMA_TIME_UNIT_S = 1000.0 * 365.25 * 24.0 * 3600.0

def GravityParameters(Planet, Params):
    """ Calculate induced gravity responses for the body and prints them to disk."""
    Timing.setFunctionTime(time.time())
    if (Planet.Do.VALID or (Params.ALLOW_BROKEN_MODELS and Planet.Do.STILL_CALCULATE_BROKEN_PROPERTIES)) and Params.CALC_NEW_GRAVITY and Params.CALC_VISCOSITY and Params.CALC_SEISMIC and not Params.SKIP_INNER:
        if not ValidateGravityOrbit(Planet):
            Timing.printFunctionTimeDifference('GravityParameters()', time.time())
            return Planet, Params
        # Check if there are any phase transitions in the model, or if this is a 1-layer model
        if not np.any(Planet.Reduced.phase[:-1] != Planet.Reduced.phase[1:]):
            # If we have only one phase in our model, ALMA cannot calculate love numbers so we must pass
            pass
        else:
            # Set Magnetic struct layer arrays as we need for induction calculations
            Planet, Params = SetupGravity(Planet, Params)
            Planet, Params = ComputeGravityObservations(Planet, Params)
            if (not Params.NO_SAVEFILE) and (not Params.INVERSION_IN_PROGRESS) and (not Params.DO_EXPLOREOGRAM):
                Planet, Params = WriteGravityParameters(Planet, Params)
    elif Planet.Do.VALID:
        if os.path.isfile(Params.DataFiles.gravityParametersFile):
            # Reload gravity parameters from disk
            Planet, Params = ReloadGravityParameters(Planet, Params)
        else:
            log.warning(
                f'CALC_NEW_GRAVITY is False, but {Params.DataFiles.gravityParametersFile} was not found. ' + f'Skipping gravity parameter calculations.')
    Timing.printFunctionTimeDifference('GravityParameters()', time.time())
    return Planet, Params




def SetupGravity(Planet, Params):
    """ Reconfigure layer boundaries and gravity model into a format
        usable by the DP ALMA gravity response backend.

        NOTE: Function follows the legacy read_model_pp layout so the model
        remains compatible with ALMA source functions.

        Requires Planet attributes:
    """
    """Combine data into model format required by the ALMA backend."""
    # Note we have to use r_m[:-1] since r_m has one extra value than other arrays
    Planet.Gravity.model = np.vstack(
        [Planet.Reduced.r_m, Planet.Reduced.phase, Planet.Reduced.rho_kgm3, Planet.Reduced.Seismic.VP_kms, Planet.Reduced.Seismic.VS_kms, Planet.Reduced.Seismic.GS_GPa, Planet.Reduced.eta_Pas]).T

    # Convert parameter units to Pa and meters
    for index, (header, unit) in enumerate(zip(Planet.Gravity.columns, Planet.Gravity.units_PyALMA3)):
        if header in Planet.Gravity.parameters_to_convert:
            Planet.Gravity.model[:, index] = Planet.Gravity.model[:, index] * Planet.Gravity.parameters_to_convert[header]
    # Get indices of used properties
    rIndex  = Planet.Gravity.columns.index('r')
    pIndex = Planet.Gravity.columns.index('phase')
    rhoIndex = Planet.Gravity.columns.index('rho')
    VPIndex = Planet.Gravity.columns.index('VP')
    VSIndex = Planet.Gravity.columns.index('VS')
    GSIndex = Planet.Gravity.columns.index('GS')
    etaIndex = Planet.Gravity.columns.index('eta')
    # Keep the legacy model hook for radius handling.
    Planet.Gravity.model[:, rIndex] = Planet.Gravity.model[:, rIndex]

    # Flip model: core at top and surface at bottom
    Planet.Gravity.model = np.flipud(Planet.Gravity.model)

    ## Calculate elastic parameters
    # LAMBDA = rho (Vp^2 - 2 * Vs^2)
    # 1st Lame parameter
    Planet.Gravity.LAMBDA_Pa = Planet.Gravity.model[:, rhoIndex] * (
            np.power(Planet.Gravity.model[:, VPIndex], 2) -
            2. * np.power(Planet.Gravity.model[:, VSIndex], 2))

    # shear modulus G or MU = rho Vs^2
    Planet.Gravity.MU_Pa = Planet.Gravity.model[:, GSIndex]

    # Poissons ratio sigma = lambda / 2*(lambda + mu)
    Planet.Gravity.SIGMA = Planet.Gravity.LAMBDA_Pa / (2 * Planet.Gravity.LAMBDA_Pa + 2 * Planet.Gravity.MU_Pa)

    # Youngs modulus Y = 2 * MU * (1 + sigma)
    Planet.Gravity.Y_Pa = 2. * Planet.Gravity.MU_Pa * (1 + Planet.Gravity.SIGMA)

    Planet.Gravity.VISCOSITY_kgms = Planet.Gravity.model[:, etaIndex]

    Planet.Gravity.ALMAModel = {'columns': Planet.Gravity.columns,
            'units': Planet.Gravity.units_PyALMA3,
            'model': Planet.Gravity.model,
            'lambda': Planet.Gravity.LAMBDA_Pa,
            'mu': Planet.Gravity.MU_Pa,
            'sigma': Planet.Gravity.SIGMA,
            'y': Planet.Gravity.Y_Pa,
            'vis': Planet.Gravity.VISCOSITY_kgms}

    # Set Planet time scale and harmonic degrees from Params
    Planet.Gravity.time_log_kyrs = Params.Gravity.time_log_kyrs
    Planet.Gravity.harmonic_degrees = Params.Gravity.harmonic_degrees

    # Finally, we must setup the rheology structure, from core to surface
    layers = [] # List of layer indices where layer change occurs (index right before the change)
    # Find where phase changes occur
    phases = Planet.Gravity.model[:, pIndex]
    # Flip changeIndices so that the largest index becomes 0 and 0 becomes the largest
    # This reverses the ordering while maintaining the relative spacing between indices
    changeIndices = np.max(Planet.Reduced.changeIndices) - np.flipud(Planet.Reduced.changeIndices)
    rheology_structure = []
    for start, end in zip(changeIndices[:-1], changeIndices[1:]):
        if end != changeIndices[-1]:  # Exclude the last index, which is the end of the model
            layers.append(end - 1)
        phase = phases[start]
        if phase >= Constants.phaseClath and phase < Constants.phaseClath + 10:
            phase = Constants.phaseClath # Reset phase to phase Clathrate so we can use that rheology model in the config - namely, we treat rheology of mixed layers as same as we treat clathrate
        convection = np.flipud(Planet.Reduced.iConv)[start]
        # Convert numerical phase to string representation for dictionary lookup
        phase_str = PhaseConv(phase, liq='0')
        if convection:
            phase_str += '_conv'
        if phase_str not in Params.Gravity.rheology_models:
            raise ValueError(f"Phase {phase_str} not found in rheology models.")
        else:
            rheology_model = Params.Gravity.rheology_models[phase_str]
        rheology_structure.append(rheology_model)
    
    # Store the compiled structures in Planet.Gravity
    Params.Gravity.rheology_structure = rheology_structure
    
    # Verify we have the right number of structural regions
    if len(rheology_structure) != len(layers) + 1:
        raise ValueError(f"Number of rheology structures ({len(rheology_structure)}) does not match "
                        f"number of layer regions ({len(layers) + 1})")
    
    # Create rheology array for each model layer
    rheo = []
    for layer_idx in range(len(rheology_structure)):
        if layer_idx == 0:
            # First region: from start to first transition
            end_idx = layers[layer_idx] + 1 if layers else len(phases)
            rheo.extend([rheology_structure[layer_idx] for _ in range(end_idx)])
            
        elif layer_idx < len(rheology_structure) - 1:
            # Intermediate regions: from previous transition to next transition
            start_idx = layers[layer_idx - 1] + 1
            end_idx = layers[layer_idx] + 1
            rheo.extend([rheology_structure[layer_idx] for _ in range(end_idx - start_idx)])
            
        else:
            # Last region: from last transition to end
            start_idx = layers[-1] + 1
            rheo.extend([rheology_structure[layer_idx] for _ in range(len(phases) - start_idx)])
    
    # Create parameters array (can be customized for Andrade/Burgers layers)
    params = np.zeros((len(rheo), 2))
    
    # Finally, create our own ALMA parameters structure for Andrade and Burgers layers
    # We do this because the infer_rheology_pp function doesn't let us customize
    # Andrade alpha/gamma or Burgers parameters directly.
    for i, rheol in enumerate(rheo):
        if rheol == 'andrade' or rheol == 6:
            if Planet.Gravity.andradAlpha is None:
                raise ValueError("Andrade rheology requires Planet.Gravity.andradAlpha.")
            if Planet.Gravity.andradGamma is None:
                Planet.Gravity.andradGamma = np.float64(mp.gamma(Planet.Gravity.andradAlpha + 1))
            params[i, 0] = Planet.Gravity.andradAlpha
            params[i, 1] = Planet.Gravity.andradGamma
        elif rheol == 'burgers' or rheol == 5:
            # Set Burgers parameters - can be customized
            params[i, 0] = Planet.Gravity.BurgerFirstParameter  # mu2/mu1 ratio
            params[i, 1] = Planet.Gravity.BurgerSecondParameter  # eta2/eta1 ratio
    
    # Store rheology and params for use in build_model
    Planet.Gravity.rheology = rheo
    Planet.Gravity.pyAlmaParams = params

    # Return Planet and Params
    return Planet, Params


def ComputeGravityObservations(Planet, Params):
    """Compute degree-2 complex Love numbers and libration amplitude."""
    if not ValidateGravityOrbit(Planet):
        return Planet, Params

    Torb_s = float(Planet.Bulk.Torb_s)
    Torb_kyr = Torb_s / ALMA_TIME_UNIT_S
    eccentricity = float(Planet.Bulk.eccentricity)
    r, phase, rho, mu, vis = GravityModelArrays(Planet)

    model_params = pyALMA3Updated.build_model(r, rho, mu, vis, Planet.Gravity.rheology,
                                       Planet.Gravity.pyAlmaParams, norm=False)

    Planet.Gravity.harmonic_degrees = [2]
    Planet.Gravity.time_log_kyrs = [Torb_kyr]
    Planet.Gravity.Torb_s = Torb_s
    Planet.Gravity.eccentricity = eccentricity

    # ALMA's time argument is in kyr, while Planet.Bulk.Torb_s is stored in seconds.
    h, l, k, y = pyALMA3Updated.love_numbers(model_params, [2], [Torb_kyr], 'tidal',
                                      'periodic', Params.Gravity.gorder,
                                      verbose=(not Params.QUIET_ALMA and Params.Gravity.verbose),
                                      internal=True)

    Planet.Gravity.h = np.complex128(h[0, 0, -1])
    Planet.Gravity.l = np.complex128(l[0, 0, -1])
    Planet.Gravity.k = np.complex128(k[0, 0, -1])
    Planet.Gravity.delta = 1 + Planet.Gravity.k - Planet.Gravity.h

    Planet.Gravity.y = y
    g0 = Constants.G * Planet.Bulk.M_kg / Planet.Bulk.R_m**2
    Planet.Gravity.y1 = np.real(y[:, 0]) * Planet.Bulk.M_kg / Planet.Bulk.R_m / g0

    omega = 2 * np.pi / Torb_s
    libration_r, libration_rho, ocean_idx, libration_y1 = LibrationModelInputs(
        r, phase, rho, Planet.Gravity.y1)
    Planet.Gravity.libration_m = librations(libration_r, libration_rho, omega,
                                            eccentricity, rigid=False,
                                            ocean=(ocean_idx is not None),
                                            ocean_idx=ocean_idx if ocean_idx is not None else 1,
                                            y1=libration_y1,
                                            rek2=np.real(Planet.Gravity.k))

    Planet = SetGravityDerivedQuantities(Planet)

    return Planet, Params


def ValidateGravityOrbit(Planet):
    """Return whether gravity/libration orbital inputs are available and valid."""
    Torb_s = getattr(Planet.Bulk, 'Torb_s', None)
    eccentricity = getattr(Planet.Bulk, 'eccentricity', None)
    body_name = getattr(Planet, 'name', 'Planet')

    if Torb_s is None or eccentricity is None:
        missing = []
        if Torb_s is None:
            missing.append('Planet.Bulk.Torb_s')
        if eccentricity is None:
            missing.append('Planet.Bulk.eccentricity')
        verb = 'are' if len(missing) > 1 else 'is'
        log.warning(
            f"{body_name} gravity/libration calculations skipped because "
            f"{' and '.join(missing)} {verb} not set in the planet defaults file."
        )
        return False

    try:
        Torb_s = float(Torb_s)
        eccentricity = float(eccentricity)
    except (TypeError, ValueError):
        log.warning(
            f"{body_name} gravity/libration calculations skipped because "
            "Planet.Bulk.Torb_s and Planet.Bulk.eccentricity must be numeric."
        )
        return False

    if not np.isfinite(Torb_s) or Torb_s <= 0:
        log.warning(
            f"{body_name} gravity/libration calculations skipped because "
            "Planet.Bulk.Torb_s must be a finite positive value in seconds."
        )
        return False
    if not np.isfinite(eccentricity) or eccentricity < 0:
        log.warning(
            f"{body_name} gravity/libration calculations skipped because "
            "Planet.Bulk.eccentricity must be a finite non-negative value."
        )
        return False

    return True


def GravityModelArrays(Planet):
    """Return core-to-surface arrays used by ALMA and libration calculations."""
    model = Planet.Gravity.ALMAModel['model']
    columns = Planet.Gravity.ALMAModel['columns']
    return (
        model[:, columns.index('r')],
        model[:, columns.index('phase')],
        model[:, columns.index('rho')],
        Planet.Gravity.ALMAModel['mu'],
        Planet.Gravity.ALMAModel['vis'],
    )


def LibrationModelInputs(r, phase, rho, y1):
    """Reduce PlanetProfile layers to the libration model expected by the helper."""
    ocean = np.where(phase == 0)[0]
    if ocean.size == 0:
        return r, rho, None, y1

    first_ocean = int(ocean[0])
    last_ocean = int(ocean[-1])
    if not np.all(phase[first_ocean:last_ocean + 1] == 0):
        raise ValueError("Libration calculations require contiguous ocean layers.")
    if first_ocean == 0 or last_ocean == len(r) - 1:
        raise ValueError(
            "Ocean libration calculations require interior, ocean, and outer shell regions."
        )

    reduced_r = np.array([r[first_ocean - 1], r[last_ocean], r[-1]], dtype=float)
    reduced_rho = np.array([
        VolumeWeightedDensity(r, rho, 0, first_ocean - 1),
        VolumeWeightedDensity(r, rho, first_ocean, last_ocean),
        VolumeWeightedDensity(r, rho, last_ocean + 1, len(r) - 1),
    ], dtype=float)
    reduced_y1 = np.array([y1[first_ocean - 1], y1[last_ocean], y1[-1]], dtype=float)

    return reduced_r, reduced_rho, 1, reduced_y1


def VolumeWeightedDensity(r, rho, start, end):
    """Compute volume-weighted density for inclusive core-to-surface layer bounds."""
    inner = np.concatenate(([0.0], r[:-1]))
    shell_volumes = r[start:end + 1]**3 - inner[start:end + 1]**3
    volume = np.sum(shell_volumes)
    if volume <= 0:
        raise ValueError("Cannot compute libration density for zero-volume layer region.")
    return np.sum(rho[start:end + 1] * shell_volumes) / volume


def SetGravityDerivedQuantities(Planet):
    """Update amplitudes and phases derived from complex Love numbers."""
    Planet.Gravity.hAmp = np.abs(Planet.Gravity.h)
    Planet.Gravity.hPhase = -np.degrees(np.angle(Planet.Gravity.h))
    Planet.Gravity.lAmp = np.abs(Planet.Gravity.l)
    Planet.Gravity.lPhase = -np.degrees(np.angle(Planet.Gravity.l))
    Planet.Gravity.kAmp = np.abs(Planet.Gravity.k)
    Planet.Gravity.kPhase = -np.degrees(np.angle(Planet.Gravity.k))
    Planet.Gravity.deltaAmp = np.abs(Planet.Gravity.delta)
    Planet.Gravity.deltaPhase = -np.degrees(np.angle(Planet.Gravity.delta))

    return Planet


def WriteGravityParameters(Planet, Params):
    """Write the gravity parameters to disk"""
    headerLines = [
        f'Orbital period (s) = {Planet.Gravity.Torb_s}',
        f'Eccentricity = {Planet.Gravity.eccentricity}',
        f'harmonic degrees = {Planet.Gravity.harmonic_degrees}',
        f'h love number = {Planet.Gravity.h}',
        f'l love number = {Planet.Gravity.l}',
        f'k love number = {Planet.Gravity.k}',
        f'delta love number = {Planet.Gravity.delta}',
        f'libration amplitude (m) = {Planet.Gravity.libration_m}'
        ]
    # Write out data from core/mantle trade
    ensure_parent_dir(Params.DataFiles.gravityParametersFile)
    with open(Params.DataFiles.gravityParametersFile, 'w') as f:
        f.write('\n  '.join(headerLines) + '\n')
    log.debug(f'Saved induced moments to file: {Params.DataFiles.gravityParametersFile}')

    return Planet, Params


def ReloadGravityParameters(Planet, Params):
    """Reload gravity paramters from disk"""
    with open(Params.DataFiles.gravityParametersFile) as f:
        reload = {}
        for line in f:
            if '=' in line:
                key, value = line.strip().split('=', 1)
                reload[key.strip()] = value.strip()

        Planet.Gravity.Torb_s = ParseOptionalFloat(reload.get('Orbital period (s)'))
        Planet.Gravity.eccentricity = ParseOptionalFloat(reload.get('Eccentricity'))
        if Planet.Gravity.Torb_s is not None and np.isfinite(Planet.Gravity.Torb_s):
            Planet.Gravity.time_log_kyrs = [Planet.Gravity.Torb_s / ALMA_TIME_UNIT_S]
        else:
            Planet.Gravity.time_log_kyrs = ParseOptionalLiteral(reload.get('Time range (log kyrs)'))
        Planet.Gravity.harmonic_degrees = ParseOptionalLiteral(reload.get('harmonic degrees')) or [2]
        Planet.Gravity.h = ParseComplex(reload['h love number'])
        Planet.Gravity.l = ParseComplex(reload['l love number'])
        Planet.Gravity.k = ParseComplex(reload['k love number'])
        if 'delta love number' in reload:
            Planet.Gravity.delta = ParseComplex(reload['delta love number'])
        else:
            Planet.Gravity.delta = 1 + Planet.Gravity.k - Planet.Gravity.h
        Planet.Gravity.libration_m = ParseOptionalFloat(reload.get('libration amplitude (m)'), default=np.nan)
        Planet = SetGravityDerivedQuantities(Planet)
    return Planet, Params


def ParseComplex(value):
    """Parse a scalar complex value from gravity save files."""
    try:
        return np.complex128(value)
    except (TypeError, ValueError):
        parsed = np.asarray(ast.literal_eval(value), dtype=np.complex128)
        if parsed.size == 1:
            return np.complex128(parsed.reshape(-1)[0])
        return parsed


def ParseOptionalFloat(value, default=None):
    """Parse optional scalar float values from gravity save files."""
    if value is None:
        return default
    try:
        return float(value)
    except ValueError:
        return default


def ParseOptionalLiteral(value):
    """Parse optional Python literal values from gravity save files."""
    if value is None:
        return None
    return ast.literal_eval(value)
