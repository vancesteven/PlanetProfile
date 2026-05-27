import numpy as np
from alma import build_model, love_numbers
import alma as _alma
import logging
import ast
import os
from PlanetProfile.Utilities.Indexing import PhaseConv
from PlanetProfile.Utilities.defineStructs import Constants, Timing
import time

# Assign logger
log = logging.getLogger('PlanetProfile')

# TidalPy conditional import — redirect config directory before TidalPy init creates ~/Documents/TidalPy
_TIDALPY_AVAILABLE = False
try:
    import platformdirs as _platformdirs
    _pp_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    _platformdirs.user_documents_dir = lambda: _pp_root
    from TidalPy.RadialSolver import radial_solver as _tidalpy_radial_solver
    from TidalPy.RadialSolver import build_rs_input_from_data as _tidalpy_build_rs_input
    from TidalPy.tides.multilayer.heating import calc_radial_volumetric_tidal_heating_from_rs_solution as _tidalpy_calc_heating
    from TidalPy.rheology import Andrade as _Andrade, Maxwell as _Maxwell, Elastic as _Elastic, Newton as _Newton
    _TIDALPY_AVAILABLE = True
except ImportError:
    pass

def GravityParameters(Planet, Params):
    """ Calculate induced gravity responses for the body and prints them to disk."""
    Timing.setFunctionTime(time.time())
    if (Planet.Do.VALID or (Params.ALLOW_BROKEN_MODELS and Planet.Do.STILL_CALCULATE_BROKEN_PROPERTIES)) and Params.CALC_NEW_GRAVITY and Params.CALC_VISCOSITY and Params.CALC_SEISMIC and not Params.SKIP_INNER:
        # Check if there are any phase transitions in the model, or if this is a 1-layer model
        if not np.any(Planet.Reduced.phase[:-1] != Planet.Reduced.phase[1:]):
            # If we have only one phase in our model, pyalma is not able to calculate love numbers so we must pass
            pass
        else:
            # Set Magnetic struct layer arrays as we need for induction calculations
            Planet, Params = SetupGravity(Planet, Params)

            # Dispatch to selected backend
            if Params.Gravity.backend == 'tidalpy':
                if not _TIDALPY_AVAILABLE:
                    raise ImportError('TidalPy backend selected but not installed. Install: pip install TidalPy')
                Planet = _run_tidalpy_backend(Planet, Params)
            else:
                # --- PyALMA3 backend ---
                # Apply per-point Andrade zeta by monkeypatching get_rheology.
                # PyALMA3 auto-sets params[i,1] = Gamma(alpha+1); dividing by zeta amplifies transient creep.
                # andrade_zeta_per_point is computed in SetupGravity from float or dict specification.
                _zeta_arr = Planet.Gravity.andrade_zeta_per_point
                _need_patch = np.any(_zeta_arr != 1.0)
                if _need_patch:
                    _orig_get_rheology = _alma.get_rheology
                    def _patched_get_rheology(rheology, nla, params, _zarr=_zeta_arr, _orig=_orig_get_rheology):
                        rheol, rpar = _orig(rheology, nla, params)
                        for i in range(nla):
                            if rheol[i] == 6 and _zarr[i] != 1.0:
                                rpar[i, 1] = rpar[i, 1] / _zarr[i]
                        return rheol, rpar
                    _alma.get_rheology = _patched_get_rheology
                model_params = build_model(Planet.Gravity.ALMAModel['model'][:, Planet.Gravity.ALMAModel['columns'].index('r')],
                                        Planet.Gravity.ALMAModel['model'][:, Planet.Gravity.ALMAModel['columns'].index('rho')],
                                        Planet.Gravity.ALMAModel['mu'],
                                        Planet.Gravity.ALMAModel['vis'],
                                        Planet.Gravity.rheology, Planet.Gravity.pyAlmaParams, ndigits=Params.Gravity.num_digits, verbose=(not Params.QUIET_ALMA and Params.Gravity.verbose), parallel=False) # False parallel for now since it is throwing a bad file desciptor
                # Restore original get_rheology
                if _need_patch:
                    _alma.get_rheology = _orig_get_rheology
                # Compute love numbers (use Planet.Gravity.time_log_kyrs which may be overridden by orbital period)
                Planet.Gravity.h, Planet.Gravity.l, Planet.Gravity.k = love_numbers(Params.Gravity.harmonic_degrees, Planet.Gravity.time_log_kyrs,
                                                                        Params.Gravity.loading_type, Params.Gravity.time_history_function, Params.Gravity.tau, model_params,
                                                                        Params.Gravity.output_type, Params.Gravity.gorder, verbose=(not Params.QUIET_ALMA and Params.Gravity.verbose),
                                                                        parallel=Params.Gravity.parallel and not (Params.INDUCTOGRAM_IN_PROGRESS or Params.DO_EXPLOREOGRAM))

            # Shared post-processing for both backends
            # Compute delta relation
            Planet.Gravity.delta = 1 + Planet.Gravity.k - Planet.Gravity.h

            # If our love numbers are 1x1 numpy array, let's convert to float - important for plotting and output
            if hasattr(Planet.Gravity.time_log_kyrs, '__len__') and hasattr(Planet.Gravity.harmonic_degrees, '__len__') \
                    and len(Planet.Gravity.time_log_kyrs) == 1 and len(Planet.Gravity.harmonic_degrees) == 1 \
                    and isinstance(Planet.Gravity.h, np.ndarray):
                Planet.Gravity.h = Planet.Gravity.h[0, 0]
                Planet.Gravity.l = Planet.Gravity.l[0, 0]
                Planet.Gravity.k = Planet.Gravity.k[0, 0]
                Planet.Gravity.delta = Planet.Gravity.delta[0, 0]
            # Convert love numbers from complex to magnitude and phase delay
            Planet.Gravity.hAmp = np.abs(Planet.Gravity.h)
            Planet.Gravity.hPhase = -np.degrees(np.angle(Planet.Gravity.h))
            Planet.Gravity.lAmp = np.abs(Planet.Gravity.l)
            Planet.Gravity.lPhase = -np.degrees(np.angle(Planet.Gravity.l))
            Planet.Gravity.kAmp = np.abs(Planet.Gravity.k)
            Planet.Gravity.kPhase = -np.degrees(np.angle(Planet.Gravity.k))
            Planet.Gravity.deltaAmp = np.abs(Planet.Gravity.delta)
            Planet.Gravity.deltaPhase = -np.degrees(np.angle(Planet.Gravity.delta))
            if (not Params.NO_SAVEFILE) and (not Params.INVERSION_IN_PROGRESS) and (not Params.DO_EXPLOREOGRAM):
                Planet, Params = WriteGravityParameters(Planet, Params)
    elif Planet.Do.VALID:
        if os.path.isfile(Params.DataFiles.gravityParametersFile):
            # Reload gravity parameters from disk
            Planet, Params = ReloadGravityParameters(Planet, Params)
        else:
            log.warning(
                f'CALC_NEW_GRAVITY is False, but {Params.DataFiles.gravityParametersFile} was not found. ' + f'Skipping gravity parameter calculations.')
    # Post-hoc check: compare user-specified ice tidal heating with value implied by computed k2
    Planet = CheckIceTidalHeatingConsistency(Planet)

    Timing.printFunctionTimeDifference('GravityParameters()', time.time())
    return Planet, Params


def CheckIceTidalHeatingConsistency(Planet):
    """ Compare user-specified Ocean.HtidalIce_Wm3 with value implied by computed Im(k2).

        Uses the Petricca et al. (2025) formula:
        H_tidal = (21/2) * n^5 * R^5 * e^2 * Im(k2) / (G * V_ice)

        Only runs if orbital parameters and Love numbers are available.
        Stores computed value in Planet.Ocean.HtidalIce_Wm3_computed for diagnostics.
    """
    # Check prerequisites
    if Planet.Bulk.eccentricity is None or Planet.Bulk.meanMotion_radps is None:
        return Planet
    if not hasattr(Planet.Gravity, 'k') or Planet.Gravity.k is None:
        return Planet

    k2 = Planet.Gravity.k
    if not np.isfinite(k2):
        return Planet

    Im_k2 = np.imag(k2) if np.iscomplex(k2) else 0.0
    if Im_k2 == 0:
        return Planet

    n = Planet.Bulk.meanMotion_radps
    R = Planet.Bulk.R_m
    e = Planet.Bulk.eccentricity
    G = Constants.G

    # Compute total tidal dissipation power: E_dot = (21/2) * n^5 * R^5 * e^2 * Im(k2) / G
    E_dot_W = (21.0 / 2.0) * n**5 * R**5 * e**2 * abs(Im_k2) / G

    # Estimate ice volume from layer structure
    V_ice_m3 = 0.0
    if hasattr(Planet, 'phase') and Planet.phase is not None:
        ice_phases = [1, 2, 3, 5, 6]  # All ice phases
        for i in range(len(Planet.phase) - 1):
            if Planet.phase[i] in ice_phases:
                # Shell volume between radii r[i] and r[i+1]
                r_outer = Planet.r_m[i]
                r_inner = Planet.r_m[i + 1]
                V_ice_m3 += (4.0 / 3.0) * np.pi * (r_outer**3 - r_inner**3)

    if V_ice_m3 <= 0:
        return Planet

    # Compute implied volumetric heating rate
    H_computed = E_dot_W / V_ice_m3

    # Store for diagnostics
    if not hasattr(Planet.Ocean, 'HtidalIce_Wm3_computed'):
        Planet.Ocean.HtidalIce_Wm3_computed = H_computed
    else:
        Planet.Ocean.HtidalIce_Wm3_computed = H_computed

    # Apply self-consistent override if enabled
    H_user = Planet.Ocean.HtidalIce_Wm3
    tidalpy_perPhase_W = getattr(Planet.Gravity, 'tidalpy_Htidal_perPhase_W', None)
    tidalpy_perPhase_Wm3 = getattr(Planet.Gravity, 'tidalpy_Htidal_perPhase_Wm3', None)

    if Planet.Do.DO_SELF_CONSISTENT_HTIDAL and tidalpy_perPhase_W is not None:
        # --- Per-phase heating from TidalPy (replaces uniform E_dot/V_ice) ---
        # When TidalPy per-phase data is available, use actual per-phase dissipation
        # to set both ice and silicate heating rates. This allows the heating
        # distribution to shift between layers depending on rheology (e.g., Maxwell
        # places most dissipation in silicates when eta_sil is near the Maxwell peak).
        ice_phase_names = ['Ih', 'II', 'III', 'V', 'VI']
        ice_total_W = sum(tidalpy_perPhase_W.get(p, 0) for p in ice_phase_names)
        sil_W = tidalpy_perPhase_W.get('Sil', 0)
        sil_Wm3 = tidalpy_perPhase_Wm3.get('Sil', 0)

        # Ice heating: volume-weighted average of all ice phases
        if V_ice_m3 > 0 and ice_total_W > 0:
            H_ice_from_tidalpy = ice_total_W / V_ice_m3
        else:
            H_ice_from_tidalpy = H_computed  # Fall back to uniform E_dot/V_ice

        Planet.Ocean.HtidalIce_Wm3 = H_ice_from_tidalpy
        Planet.Ocean.HtidalIce_Wm3_computed = H_ice_from_tidalpy

        # Silicate heating: use TidalPy prediction directly
        H_sil_user = getattr(Planet.Sil, 'Htidal_Wm3', 0) or 0
        if sil_Wm3 > 0:
            Planet.Sil.Htidal_Wm3 = sil_Wm3

        log.info(f'Self-consistent tidal heating from TidalPy per-phase data '
                 f'(E_dot = {E_dot_W/1e9:.1f} GW, Im(k2) = {Im_k2:.4f}):')
        log.info(f'  Ice:      {ice_total_W/1e9:.1f} GW → HtidalIce = {H_ice_from_tidalpy:.2e} W/m^3 '
                 f'(was {H_user:.2e})')
        log.info(f'  Silicate: {sil_W/1e9:.1f} GW → Htidal_sil = {sil_Wm3:.2e} W/m^3 '
                 f'(was {H_sil_user:.2e})')
        for p in ice_phase_names:
            if p in tidalpy_perPhase_W and tidalpy_perPhase_W[p] > 0:
                log.info(f'    {p:6s}: {tidalpy_perPhase_W[p]/1e9:8.1f} GW  '
                         f'({tidalpy_perPhase_Wm3.get(p, 0):.2e} W/m^3)')
        log.info(f'  Note: re-run for full convergence with updated heating.')
        _LogLayerHeatFluxes(Planet)
        return Planet

    elif Planet.Do.DO_SELF_CONSISTENT_HTIDAL and H_computed > 0:
        # Fallback: no TidalPy per-phase data — use uniform E_dot / V_ice
        Planet.Ocean.HtidalIce_Wm3 = H_computed
        log.info(f'Self-consistent tidal heating: overriding user value {H_user:.2e} W/m^3 '
                 f'with k2-implied {H_computed:.2e} W/m^3 (Im(k2) = {Im_k2:.4f}). '
                 f'Note: convection was computed with the original value; '
                 f're-run for full convergence.')
        _LogLayerHeatFluxes(Planet)
        return Planet
    elif H_user > 0 and H_computed > 0:
        ratio = H_computed / H_user
        if ratio > 10 or ratio < 0.1:
            log.warning(f'Ice tidal heating inconsistency: user-specified {H_user:.2e} W/m^3, '
                        f'k2-implied {H_computed:.2e} W/m^3 (ratio {ratio:.1f}x). '
                        f'Im(k2) = {Im_k2:.4f}, V_ice = {V_ice_m3:.3e} m^3')
        else:
            log.info(f'Ice tidal heating consistency check: user {H_user:.2e}, '
                     f'k2-implied {H_computed:.2e} W/m^3 (ratio {ratio:.1f}x)')
        _LogLayerHeatFluxes(Planet)
    elif H_computed > 0:
        log.info(f'k2-implied ice tidal heating: {H_computed:.2e} W/m^3 '
                 f'(Im(k2) = {Im_k2:.4f}). Set Ocean.HtidalIce_Wm3 to use in convection.')

    return Planet


def _LogLayerHeatFluxes(Planet):
    """ Log heat flux at each major layer boundary with per-layer HP ice breakdown.

        Approach:
        1. Compute silicate heat (radiogenic + silicate tidal) directly from
           layer densities, volumes, Qrad_Wkg, and Htidal_Wm3.
        2. For each HP ice phase (VI, V, III) independently, compute tidal
           heating power from shell volume * HtidalIce_Wm3.
        3. Track cumulative heat bottom-to-top: start at Q_sil, add each
           HP ice layer's contribution at its upper boundary.
        4. Compare the cumulative sum to the surface total (QfromMantle_W,
           set by the ice I thermal model) — any residual is the ice I
           contribution or a model consistency gap.

        No equilibrium assumption: each layer's heating is reported
        independently from its geometry and applied heating rate.
    """
    Q_total_W = Planet.Ocean.QfromMantle_W
    if Q_total_W is None or not np.isfinite(Q_total_W) or Q_total_W == 0:
        return

    R_m = Planet.Bulk.R_m

    # --- Silicate heat: radiogenic and tidal computed separately ---
    Q_rad_W = 0.0
    Q_sil_tidal_W = 0.0
    M_sil_kg = 0.0
    Qrad_Wkg = getattr(Planet.Sil, 'Qrad_Wkg', 0.0) or 0.0
    Htidal_sil = getattr(Planet.Sil, 'Htidal_Wm3', 0.0) or 0.0
    if hasattr(Planet.Sil, 'Rmean_m') and Planet.Sil.Rmean_m is not None:
        if hasattr(Planet, 'phase') and Planet.phase is not None:
            silPhases = [20, 21, 22, 23, 25, 26]
            for i in range(len(Planet.phase) - 1):
                if Planet.phase[i] in silPhases or (50 <= Planet.phase[i] <= 59):
                    r_o, r_i = Planet.r_m[i], Planet.r_m[i + 1]
                    dV = (4.0 / 3.0) * np.pi * (r_o**3 - r_i**3)
                    rho = Planet.rho_kgm3[i] if hasattr(Planet, 'rho_kgm3') else 3500.0
                    M_sil_kg += rho * dV
                    Q_rad_W += Qrad_Wkg * rho * dV
                    Q_sil_tidal_W += Htidal_sil * dV
    Q_sil_W = Q_rad_W + Q_sil_tidal_W

    # --- Per-layer HP ice tidal heating from geometry ---
    # When TidalPy per-phase data is available, use actual per-phase heating rates.
    # Otherwise, fall back to uniform HtidalIce_Wm3 * volume.
    tidalpy_perPhase_Wm3 = getattr(Planet.Gravity, 'tidalpy_Htidal_perPhase_Wm3', None)
    tidalpy_perPhase_W = getattr(Planet.Gravity, 'tidalpy_Htidal_perPhase_W', None)
    Htidal_ice = getattr(Planet.Ocean, 'HtidalIce_Wm3', 0.0) or 0.0
    hp_layers = []  # list of (phaseID, phaseName, r_top, r_bot, V, Q_tidal)
    if hasattr(Planet.Steps, 'nSurfIce') and hasattr(Planet, 'phase'):
        iOcStart = Planet.Steps.nSurfIce
        iOcEnd = iOcStart + getattr(Planet.Steps, 'nOceanMax', 0)
        if iOcEnd > iOcStart:
            phaseNames = {3: 'III', 5: 'V', 6: 'VI'}
            for phaseID in [6, 5, 3]:  # bottom to top
                inds = np.where(Planet.phase[iOcStart:iOcEnd] == phaseID)[0] + iOcStart
                if len(inds) == 0:
                    continue
                r_top = Planet.r_m[inds[0]]
                r_bot = Planet.r_m[inds[-1] + 1] if inds[-1] + 1 < len(Planet.r_m) else Planet.r_m[inds[-1]]
                V = (4.0 / 3.0) * np.pi * (r_top**3 - r_bot**3)
                pName = phaseNames[phaseID]
                # Use TidalPy per-phase data if available
                if tidalpy_perPhase_W is not None and pName in tidalpy_perPhase_W:
                    Q_tidal = tidalpy_perPhase_W[pName]
                else:
                    Q_tidal = Htidal_ice * V
                hp_layers.append((phaseID, pName, r_top, r_bot, V, Q_tidal))

    Q_hp_total_W = sum(layer[5] for layer in hp_layers)

    # --- Build display ---
    lines = []

    # Surface
    qSurf = 1e3 * Q_total_W / (4 * np.pi * R_m**2)
    lines.append(f'  Surface:            {qSurf:.2f} mW/m^2  ({Q_total_W/1e9:.1f} GW)')

    # Base of ice I
    if hasattr(Planet.Steps, 'nIbottom') and Planet.Steps.nIbottom is not None:
        r = Planet.r_m[Planet.Steps.nIbottom]
        if r > 0:
            # Cumulative heat at base of ice I = silicate + all HP ice tidal
            Q_base_iceI = Q_sil_W + Q_hp_total_W
            q = 1e3 * Q_base_iceI / (4 * np.pi * r**2)
            lines.append(f'  Base of ice I:      {q:.2f} mW/m^2  ({Q_base_iceI/1e9:.1f} GW, r = {r/1e3:.1f} km)')

    # Per-layer HP ice breakdown (top to bottom, matching spatial order from surface down)
    if hp_layers:
        # Track cumulative heat bottom-to-top for boundary values
        Q_cumulative = Q_sil_W
        boundary_heat = [(Q_cumulative, 'below ice VI')]  # at base of HP stack
        for phaseID, phaseName, r_top, r_bot, V, Q_tidal in hp_layers:
            Q_cumulative += Q_tidal
            boundary_heat.append((Q_cumulative, f'above ice {phaseName}'))

        # Display top-to-bottom (reverse of computation order)
        source = 'TidalPy per-layer' if tidalpy_perPhase_W is not None else f'uniform H_tidal = {Htidal_ice:.2e} W/m^3'
        lines.append(f'  --- HP ice layer breakdown ({source}) ---')
        for phaseID, phaseName, r_top, r_bot, V, Q_tidal in reversed(hp_layers):
            thick_km = (r_top - r_bot) / 1e3
            q_top = 1e3 * Q_tidal / (4 * np.pi * r_top**2) if r_top > 0 else 0
            lines.append(f'    Ice {phaseName:3s}: {thick_km:6.1f} km thick, '
                         f'V = {V:.2e} m^3, '
                         f'Q_tidal = {Q_tidal/1e9:.1f} GW '
                         f'({q_top:.2f} mW/m^2 equivalent)')

    # Top of silicates
    if hasattr(Planet.Sil, 'Rmean_m') and Planet.Sil.Rmean_m is not None and Planet.Sil.Rmean_m > 0:
        r = Planet.Sil.Rmean_m
        q = 1e3 * Q_sil_W / (4 * np.pi * r**2)
        lines.append(f'  Top of silicates:   {q:.2f} mW/m^2  ({Q_sil_W/1e9:.1f} GW, r = {r/1e3:.1f} km)')

    # Heat budget summary with consistency check
    Q_iceI_residual = Q_total_W - Q_sil_W - Q_hp_total_W
    lines.append(f'  --- Heat budget ---')
    lines.append(f'    Sil. radiogenic:       {Q_rad_W/1e9:8.1f} GW  (Qrad={Qrad_Wkg:.2e} W/kg, M_sil={M_sil_kg:.3e} kg)')
    lines.append(f'    Sil. tidal:            {Q_sil_tidal_W/1e9:8.1f} GW  (Htidal={Htidal_sil:.2e} W/m^3)')
    for phaseID, phaseName, r_top, r_bot, V, Q_tidal in hp_layers:
        lines.append(f'    Ice {phaseName:3s} tidal:        {Q_tidal/1e9:8.1f} GW')
    if Q_hp_total_W > 0:
        lines.append(f'    HP ice total:          {Q_hp_total_W/1e9:8.1f} GW  (HtidalIce={Htidal_ice:.2e} W/m^3)')
    lines.append(f'    Ice I residual:        {Q_iceI_residual/1e9:8.1f} GW')
    lines.append(f'    Surface total:         {Q_total_W/1e9:8.1f} GW')

    log.info('Heat flux at layer boundaries:\n' + '\n'.join(lines))



def _pyalma_zeta_to_tidalpy(zeta_pa, alpha):
    """ Convert PyALMA3-convention Andrade zeta to TidalPy-convention.

        PyALMA3 treats zeta as a linear inverse amplitude on the Andrade creep term:
            J_A = [Gamma(alpha+1) / zeta] * (1/mu) * (s*eta/mu)^(-alpha)

        TidalPy treats zeta as a timescale ratio entering the power law:
            J_A = Gamma(alpha+1) * (1/mu) * (omega*tau_M*zeta)^(-alpha)

        Equating: zeta_tidalpy = zeta_pyalma^(1/alpha)

        For zeta_pa = 1.0, the two are identical (zeta_tp = 1.0).
    """
    if zeta_pa == 1.0:
        return 1.0
    return zeta_pa ** (1.0 / alpha)


def _map_rheology_to_tidalpy(rheology_structure, region_phases, alpha, zeta_spec):
    """ Map PlanetProfile per-layer rheology strings to TidalPy rheology objects.

        Args:
            rheology_structure: list of str, one per region ('andrade', 'maxwell', 'elastic', 'newton')
            region_phases: list of str, phase string per region ('Ih', 'III', 'Sil', etc.)
            alpha: float, Andrade exponent (same for all Andrade layers)
            zeta_spec: float or dict, Andrade zeta per phase (PyALMA3 convention)

        Returns:
            shear_rheology_tuple: tuple of TidalPy rheology instances, one per layer
            bulk_rheology_tuple: tuple of TidalPy Elastic instances (bulk is always elastic)
    """
    shear_models = []
    for rheol_str, phase_str in zip(rheology_structure, region_phases):
        base_phase = phase_str.replace('_conv', '')
        if rheol_str == 'andrade':
            if isinstance(zeta_spec, dict):
                zeta_pa = zeta_spec.get(base_phase, zeta_spec.get('default', 1.0))
            else:
                zeta_pa = zeta_spec
            zeta_tp = _pyalma_zeta_to_tidalpy(zeta_pa, alpha)
            log.debug(f'TidalPy Andrade: phase={base_phase}, zeta_pyalma={zeta_pa:.4g} -> '
                      f'zeta_tidalpy={zeta_tp:.4e} (alpha={alpha})')
            shear_models.append(_Andrade(args=(alpha, zeta_tp)))
        elif rheol_str == 'maxwell':
            shear_models.append(_Maxwell())
        elif rheol_str == 'newton':
            shear_models.append(_Newton())
        elif rheol_str == 'elastic':
            shear_models.append(_Elastic())
        else:
            log.warning(f'Unknown rheology "{rheol_str}" for phase {phase_str}, defaulting to Maxwell')
            shear_models.append(_Maxwell())
    bulk_models = [_Elastic() for _ in rheology_structure]
    return tuple(shear_models), tuple(bulk_models)


def _run_tidalpy_backend(Planet, Params):
    """ Run TidalPy radial solver as alternative to PyALMA3.

        Uses the same SetupGravity model arrays (already called before dispatch).
        Produces Love numbers (h, k, l) and per-layer volumetric tidal heating.

        The model array after SetupGravity is in core-to-surface order (ascending r)
        with columns: r, phase, rho, VP(m/s), VS(m/s), GS(Pa), eta(Pa·s).
    """
    model = Planet.Gravity.ALMAModel['model']
    rIndex = Planet.Gravity.columns.index('r')
    rhoIndex = Planet.Gravity.columns.index('rho')
    VPIndex = Planet.Gravity.columns.index('VP')
    GSIndex = Planet.Gravity.columns.index('GS')
    etaIndex = Planet.Gravity.columns.index('eta')
    pIndex = Planet.Gravity.columns.index('phase')

    # Extract arrays (already in ascending r after flipud in SetupGravity)
    r_m = model[:, rIndex].astype(np.float64)
    rho = model[:, rhoIndex].astype(np.float64)
    mu_Pa = model[:, GSIndex].astype(np.float64)  # Shear modulus
    VP_ms = model[:, VPIndex].astype(np.float64)
    eta_Pas = model[:, etaIndex].astype(np.float64)
    phases = model[:, pIndex]

    # Bulk modulus: K = rho * VP^2 - (4/3) * mu
    # When VP is NaN (e.g. silicate layers with constant assumed properties),
    # estimate K from mu using typical Poisson's ratio: K = 2*mu*(1+nu)/(3*(1-2*nu))
    K_Pa = rho * VP_ms**2 - (4.0 / 3.0) * mu_Pa
    nan_mask = ~np.isfinite(K_Pa) | (K_Pa <= 0)
    if np.any(nan_mask):
        # Default Poisson's ratio by phase for K estimation
        for i in np.where(nan_mask)[0]:
            ph = int(phases[i])
            if ph >= 50 and ph < 100:  # silicate
                nu = 0.25
            elif ph >= 100:  # iron
                nu = 0.29
            else:  # ice/clathrate
                nu = 0.33
            K_Pa[i] = 2.0 * mu_Pa[i] * (1.0 + nu) / (3.0 * (1.0 - 2.0 * nu))
        log.info(f'TidalPy: estimated bulk modulus for {np.sum(nan_mask)} points with NaN VP')
    K_Pa = np.maximum(K_Pa, 1e6)  # Floor to avoid negative/zero

    # Bulk viscosity: set to 0 (elastic bulk response)
    bulk_visc = np.zeros_like(eta_Pas)

    # Forcing frequency
    omega = Planet.Bulk.meanMotion_radps
    if omega is None:
        raise ValueError('TidalPy backend requires Planet.Bulk.meanMotion_radps to be set')

    # Layer boundaries from changeIndices (flipped to match ascending-r model)
    changeIndices = np.max(Planet.Reduced.changeIndices) - np.flipud(Planet.Reduced.changeIndices)
    n_layers = len(changeIndices) - 1

    # Per-layer metadata
    layer_upper_radii = []
    layer_types = []
    layer_is_static = []
    layer_is_incompressible = []

    for i_layer in range(n_layers):
        start = changeIndices[i_layer]
        end = changeIndices[i_layer + 1]
        upper_r = r_m[end - 1]  # Last point in this layer
        layer_upper_radii.append(upper_r)

        # Determine layer type from phase: 0 = liquid, everything else = solid
        layer_phase = phases[start]
        if layer_phase == 0:
            layer_types.append('liquid')
        else:
            layer_types.append('solid')

        # All layers dynamic, compressible
        layer_is_static.append(False)
        layer_is_incompressible.append(False)

    layer_upper_radii_arr = np.array(layer_upper_radii, dtype=np.float64)

    # Ensure every layer has at least 5 radial points (TidalPy minimum).
    # Thin layers (e.g., lithospheric slivers at convection boundaries) may have
    # fewer points. Interpolate additional points within such layers.
    MIN_POINTS = 5
    needs_padding = False
    for i_layer in range(n_layers):
        n_pts = changeIndices[i_layer + 1] - changeIndices[i_layer]
        if n_pts < MIN_POINTS:
            needs_padding = True
            break

    if needs_padding:
        # Compute region_phases from ORIGINAL indices before padding changes them
        _orig_iConv = np.flipud(Planet.Reduced.iConv)
        region_phases = []
        for i_layer in range(n_layers):
            start = changeIndices[i_layer]
            phase = phases[start]
            if phase >= Constants.phaseClath and phase < Constants.phaseClath + 10:
                phase = Constants.phaseClath
            convection = _orig_iConv[start]
            phase_str = PhaseConv(phase, liq='0')
            if convection:
                phase_str += '_conv'
            region_phases.append(phase_str)

        new_r, new_rho, new_K, new_mu, new_eta, new_phases, new_bulk_visc = \
            [], [], [], [], [], [], []
        new_changeIndices = [0]

        for i_layer in range(n_layers):
            s, e = changeIndices[i_layer], changeIndices[i_layer + 1]
            n_pts = e - s

            if n_pts < MIN_POINTS and n_pts >= 2:
                # Interpolate to MIN_POINTS within this layer
                r_layer = r_m[s:e]
                r_interp = np.linspace(r_layer[0], r_layer[-1], MIN_POINTS)
                new_r.append(r_interp)
                new_rho.append(np.interp(r_interp, r_layer, rho[s:e]))
                new_K.append(np.interp(r_interp, r_layer, K_Pa[s:e]))
                new_mu.append(np.interp(r_interp, r_layer, mu_Pa[s:e]))
                new_eta.append(np.interp(r_interp, r_layer, eta_Pas[s:e]))
                new_phases.append(np.full(MIN_POINTS, phases[s]))
                new_bulk_visc.append(np.zeros(MIN_POINTS))
                new_changeIndices.append(new_changeIndices[-1] + MIN_POINTS)
                log.info(f'TidalPy: interpolated layer {i_layer} from {n_pts} to {MIN_POINTS} points '
                         f'(r = {r_layer[0]/1e3:.0f}-{r_layer[-1]/1e3:.0f} km)')
            else:
                new_r.append(r_m[s:e])
                new_rho.append(rho[s:e])
                new_K.append(K_Pa[s:e])
                new_mu.append(mu_Pa[s:e])
                new_eta.append(eta_Pas[s:e])
                new_phases.append(phases[s:e])
                new_bulk_visc.append(bulk_visc[s:e])
                new_changeIndices.append(new_changeIndices[-1] + n_pts)

        r_m = np.concatenate(new_r)
        rho = np.concatenate(new_rho)
        K_Pa = np.concatenate(new_K)
        mu_Pa = np.concatenate(new_mu)
        eta_Pas = np.concatenate(new_eta)
        phases = np.concatenate(new_phases)
        bulk_visc = np.concatenate(new_bulk_visc)
        changeIndices = np.array(new_changeIndices)

        # Update layer upper radii from new arrays
        for i_layer in range(n_layers):
            end = changeIndices[i_layer + 1]
            layer_upper_radii[i_layer] = r_m[end - 1]

    # Build TidalPy rheology objects per layer
    rheology_structure = Params.Gravity.rheology_structure
    # region_phases was already computed before padding (uses original indices)
    # If not padded, compute now from current arrays
    if not needs_padding:
        _orig_iConv = np.flipud(Planet.Reduced.iConv)
        _orig_changeIndices = np.max(Planet.Reduced.changeIndices) - np.flipud(Planet.Reduced.changeIndices)
        region_phases = []
        for i_layer in range(n_layers):
            start = _orig_changeIndices[i_layer]
            phase = phases[start]
            if phase >= Constants.phaseClath and phase < Constants.phaseClath + 10:
                phase = Constants.phaseClath
            convection = _orig_iConv[start]
            phase_str = PhaseConv(phase, liq='0')
            if convection:
                phase_str += '_conv'
            region_phases.append(phase_str)

    shear_rheology, bulk_rheology = _map_rheology_to_tidalpy(
        rheology_structure, region_phases,
        Planet.Gravity.andradExponent,
        Planet.Gravity.andrade_zeta
    )

    log.info(f'TidalPy backend: {n_layers} layers, omega = {omega:.4e} rad/s, '
             f'R = {r_m[-1]/1e3:.0f} km, {len(r_m)} radial points')
    for i, (rt, rp) in enumerate(zip(layer_types, region_phases)):
        log.debug(f'  Layer {i}: {rp}, type={rt}, upper_r={layer_upper_radii[i]/1e3:.0f} km')

    # Build input data (applies rheology to compute complex moduli)
    build_data = _tidalpy_build_rs_input(
        omega,
        np.ascontiguousarray(r_m),
        np.ascontiguousarray(rho),
        np.ascontiguousarray(K_Pa),
        np.ascontiguousarray(mu_Pa),
        np.ascontiguousarray(bulk_visc),
        np.ascontiguousarray(eta_Pas),
        tuple(layer_upper_radii),
        tuple(layer_types),
        tuple(layer_is_static),
        tuple(layer_is_incompressible),
        shear_rheology,
        bulk_rheology,
        perform_checks=True,
        warnings=True
    )

    # Run radial solver
    result = _tidalpy_radial_solver(
        build_data.radius_array,
        build_data.density_array,
        build_data.complex_bulk_modulus_array,
        build_data.complex_shear_modulus_array,
        build_data.frequency,
        build_data.planet_bulk_density,
        build_data.layer_types,
        build_data.is_static_bylayer,
        build_data.is_incompressible_bylayer,
        build_data.upper_radius_bylayer_array,
        degree_l=2,
        solve_for=('tidal',),
        verbose=not Params.QUIET_ALMA and Params.Gravity.verbose,
        raise_on_fail=False
    )

    if not result.success:
        log.error(f'TidalPy radial solver failed: {result.message} (error code {result.error_code})')
        return Planet

    # Extract Love numbers
    Planet.Gravity.k = complex(result.k)
    Planet.Gravity.h = complex(result.h)
    Planet.Gravity.l = complex(result.l)

    log.info(f'TidalPy Love numbers: k2 = {Planet.Gravity.k:.4f}, '
             f'h2 = {Planet.Gravity.h:.4f}, l2 = {Planet.Gravity.l:.4f}')

    # Per-layer tidal heating
    if Planet.Bulk.eccentricity is not None and Planet.Bulk.eccentricity > 0:
        parent = Planet.parent
        if parent in Constants.parentMass_kg:
            host_mass = Constants.parentMass_kg[parent]
            # Semi-major axis from Kepler III: a = (G*M_host / n^2)^(1/3)
            a_m = (Constants.G * host_mass / omega**2) ** (1.0 / 3.0)

            heating_profile = _tidalpy_calc_heating(
                Planet.Bulk.eccentricity,
                omega,
                a_m,
                host_mass,
                result,
                perform_checks=True
            )

            Planet.Gravity.tidalpy_Htidal_Wm3 = heating_profile

            # Aggregate per-phase heating using proper radial integral:
            # Total power = integral of H(r) * 4*pi*r^2 dr (not mean(H) * V)
            # The mean(H)*V approximation overestimates when H varies within
            # a layer because it weights all radial points equally instead of
            # by their shell volume element 4*pi*r^2*dr.
            result_radii = np.asarray(result.radius_array)
            phase_map = {0: '0', 1: 'Ih', 2: 'II', 3: 'III', 5: 'V', 6: 'VI'}
            perPhase_Wm3 = {}
            perPhase_W = {}
            perPhase_V = {}  # Track volume for computing volume-averaged Wm3

            for i_layer in range(n_layers):
                start_idx = changeIndices[i_layer]
                end_idx = changeIndices[i_layer + 1]
                layer_phase = int(phases[start_idx])
                phase_str = phase_map.get(layer_phase, PhaseConv(layer_phase, liq='0'))

                # Find radial points belonging to this layer in the result arrays
                r_lo = r_m[start_idx]
                r_hi = r_m[end_idx - 1]
                mask = (result_radii >= r_lo - 1.0) & (result_radii <= r_hi + 1.0)

                if np.any(mask):
                    layer_heating = heating_profile[mask]
                    layer_radii = result_radii[mask]

                    # Proper integral: total_power = integral H(r) * 4*pi*r^2 dr
                    if len(layer_radii) > 1:
                        total_power = np.trapezoid(layer_heating * 4.0 * np.pi * layer_radii**2, layer_radii)
                    else:
                        # Single point — fall back to mean*V
                        V_layer = (4.0 / 3.0) * np.pi * (r_hi**3 - r_lo**3)
                        total_power = layer_heating[0] * V_layer

                    V_layer = (4.0 / 3.0) * np.pi * (r_hi**3 - r_lo**3)

                    if phase_str in perPhase_W:
                        perPhase_W[phase_str] += total_power
                        perPhase_V[phase_str] += V_layer
                    else:
                        perPhase_W[phase_str] = total_power
                        perPhase_V[phase_str] = V_layer

            # Compute volume-averaged heating rate per phase
            for phase_str in perPhase_W:
                V = perPhase_V[phase_str]
                perPhase_Wm3[phase_str] = perPhase_W[phase_str] / V if V > 0 else 0.0

            Planet.Gravity.tidalpy_Htidal_perPhase_Wm3 = perPhase_Wm3
            Planet.Gravity.tidalpy_Htidal_perPhase_W = perPhase_W

            # Log per-phase breakdown
            total_heating = sum(perPhase_W.values())
            log.info(f'TidalPy per-layer tidal heating (total = {total_heating/1e9:.2f} GW):')
            for ps, power in sorted(perPhase_W.items()):
                log.info(f'  {ps:6s}: {power/1e9:8.2f} GW  ({perPhase_Wm3[ps]:.3e} W/m^3)')
        else:
            log.warning(f'Cannot compute per-layer heating: no parent mass for "{parent}"')

    return Planet


def SetupGravity(Planet, Params):
    """ Reconfigure layer boundaries and gravity model into a format
        usable by gravity response calculation functions of PyALMA3.

        NOTE: Function models the read_model_pp of PyALMA3 so that we can ensure compatibility with the package's functions.

        Requires Planet attributes:
    """
    """Combine data into model format that is required by PyALMA3"""
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
    # Round radius - In similar PyALMA3 function, it rounds radius so will do so here as well
    Planet.Gravity.model[:, rIndex] = np.round(Planet.Gravity.model[:, rIndex], 0)

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
    # Override time_log_kyrs with tidal forcing period if orbital parameters are available
    if Planet.Bulk.meanMotion_radps is not None:
        tidalPeriod_s = 2 * np.pi / Planet.Bulk.meanMotion_radps
        tidalPeriod_kyrs = tidalPeriod_s / (1000 * 365.25 * 24 * 3600)
        Planet.Gravity.time_log_kyrs = [tidalPeriod_kyrs]
        log.info(f'Using tidal forcing period from mean motion: {tidalPeriod_s/86400:.3f} days = {tidalPeriod_kyrs:.4e} kyrs')
    else:
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
    region_phases = []  # Track phase string per region for per-phase zeta
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
        region_phases.append(phase_str)
    
    # Store the compiled structures in Planet.Gravity
    Params.Gravity.rheology_structure = rheology_structure
    
    # Verify we have the right number of structural regions
    if len(rheology_structure) != len(layers) + 1:
        raise ValueError(f"Number of rheology structures ({len(rheology_structure)}) does not match "
                        f"number of layer regions ({len(layers) + 1})")
    
    # Create rheology array and per-point phase string array for each model layer
    rheo = []
    phase_per_point = []
    for layer_idx in range(len(rheology_structure)):
        if layer_idx == 0:
            # First region: from start to first transition
            end_idx = layers[layer_idx] + 1 if layers else len(phases)
            n = end_idx
        elif layer_idx < len(rheology_structure) - 1:
            # Intermediate regions: from previous transition to next transition
            start_idx = layers[layer_idx - 1] + 1
            end_idx = layers[layer_idx] + 1
            n = end_idx - start_idx
        else:
            # Last region: from last transition to end
            start_idx = layers[-1] + 1
            n = len(phases) - start_idx
        rheo.extend([rheology_structure[layer_idx]] * n)
        phase_per_point.extend([region_phases[layer_idx]] * n)
    
    # Create parameters array (can be customized for Andrade/Burgers layers)
    params = np.zeros((len(rheo), 2))
    
    # Finally, let's create our own PyALMA3 parameters structure for Andrade and Burgers layers
    # We do this because the infer_rheology_pp function doesn't let you customize the Andrade exponent
    for i, rheol in enumerate(rheo):
        if rheol == 'andrade' or rheol == 6:
            # Set Andrade exponent (alpha) - default to 0.2, can be customized
            params[i, 0] = Planet.Gravity.andradExponent
            # params[i, 1] will be auto-calculated as gamma(alpha + 1)
        elif rheol == 'burgers' or rheol == 5:
            # Set Burgers parameters - can be customized
            params[i, 0] = Planet.Gravity.BurgerFirstParameter  # mu2/mu1 ratio
            params[i, 1] = Planet.Gravity.BurgerSecondParameter  # eta2/eta1 ratio
    
    # Store rheology and params for use in build_model
    Planet.Gravity.rheology = rheo
    Planet.Gravity.pyAlmaParams = params

    # Compute per-point Andrade zeta array for use in monkeypatch.
    # andrade_zeta can be a float (uniform) or dict mapping phase strings to zeta values.
    # Phase strings follow PhaseConv convention: 'Ih', 'III', 'V', 'VI', 'Sil', 'Fe', 'Clath', etc.
    # Suffix '_conv' on phase strings is stripped for dict lookup; use base phase name.
    _zeta_spec = Planet.Gravity.andrade_zeta
    zeta_per_point = np.ones(len(rheo))
    if isinstance(_zeta_spec, dict):
        for i, ps in enumerate(phase_per_point):
            base_ps = ps.replace('_conv', '')
            if base_ps in _zeta_spec:
                zeta_per_point[i] = _zeta_spec[base_ps]
            elif 'default' in _zeta_spec:
                zeta_per_point[i] = _zeta_spec['default']
    elif _zeta_spec != 1.0:
        zeta_per_point[:] = _zeta_spec
    Planet.Gravity.andrade_zeta_per_point = zeta_per_point

    # Return Planet and Params
    return Planet, Params



def WriteGravityParameters(Planet, Params):
    """Write the gravity parameters to disk"""
    headerLines = [
        f'Time range (log kyrs) = {Planet.Gravity.time_log_kyrs}',
        f'harmonic degrees = {Planet.Gravity.harmonic_degrees}',
        f'h love number = {Planet.Gravity.h}',
        f'l love number = {Planet.Gravity.l}',
        f'k love number = {Planet.Gravity.k}'
        ]
    # Write out data from core/mantle trade
    with open(Params.DataFiles.gravityParametersFile, 'w') as f:
        f.write('\n  '.join(headerLines) + '\n')
    log.debug(f'Saved induced moments to file: {Params.DataFiles.gravityParametersFile}')

    return Planet, Params


def ReloadGravityParameters(Planet, Params):
    """Reload gravity paramters from disk"""
    with open(Params.DataFiles.gravityParametersFile) as f:
        Planet.Gravity.time_log_kyrs = ast.literal_eval(f.readline().split('=')[-1])
        Planet.Gravity.harmonic_degrees = ast.literal_eval(f.readline().split('=')[-1])
        if len(Planet.Gravity.time_log_kyrs) == 1 and len(Planet.Gravity.harmonic_degrees) == 1:
            # Make these variables floats
            Planet.Gravity.h = np.complex_((f.readline().split('=')[-1]))
            Planet.Gravity.l = np.complex_((f.readline().split('=')[-1]))
            Planet.Gravity.k = np.complex_((f.readline().split('=')[-1]))
        else:
            Planet.Gravity.h = np.array(ast.literal_eval(f.readline().split('=')[-1]), dtype=np.complex_)
            Planet.Gravity.l = np.array(ast.literal_eval(f.readline().split('=')[-1]), dtype=np.complex_)
            Planet.Gravity.k = np.array(ast.literal_eval(f.readline().split('=')[-1]), dtype=np.complex_)
        Planet.Gravity.delta = 1 + Planet.Gravity.k - Planet.Gravity.h
        # Convert love numbers from complex to magnitude and phase delay
        Planet.Gravity.hAmp = np.abs(Planet.Gravity.h)
        Planet.Gravity.hPhase = -np.degrees(np.angle(Planet.Gravity.h))
        Planet.Gravity.lAmp = np.abs(Planet.Gravity.l)
        Planet.Gravity.lPhase = -np.degrees(np.angle(Planet.Gravity.l))
        Planet.Gravity.kAmp = np.abs(Planet.Gravity.k)
        Planet.Gravity.kPhase = -np.degrees(np.angle(Planet.Gravity.k))
        Planet.Gravity.deltaAmp = np.abs(Planet.Gravity.delta)
        Planet.Gravity.deltaPhase = -np.degrees(np.angle(Planet.Gravity.delta))
    return Planet, Params
