"""
Body-agnostic structure-cache builder for the MCMC inference toolkit.

Phase C1 Stage 3 deliverable. Mirrors the structure-build logic in
``PlanetProfile.Test.Test50_mcmc_andrade_noocean_yao2014._build_single_structure``
with the two body-specific anchors (PP test-module name and surface radius)
parameterised. Test50 is left unchanged in this commit; a future cleanup can
delete its private copy and import from here.

Two public entry points:

    build_single_structure(planet_template_module, Tb_K)
        Run PlanetProfile once at one ``Tb_K`` and return the per-cell
        cache dict consumed by ``MCMCRunner._make_flexible_log_likelihood``.

    build_tb_grid_cache(planet_template_module, Tb_grid, output_path)
        Loop ``build_single_structure`` over a ``Tb_K`` grid, save as
        ``{'Tb_K_grid': [...], 'structures': [{...}, ...]}`` pickle.

The output dict shape matches ``Test50``'s exactly so the existing runner
paths (single-structure cache, grid-interpolated cache, v1 cached-CMR² and
v2 mass-conservation paths) all consume it without changes.
"""
from __future__ import annotations

import importlib
import logging
import os
import pickle
import sys
import tempfile
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

# Heavy PP imports are deferred to call-time so importing this module does
# not initialise EOS tables. Tests of the cache_builder API should avoid
# triggering ``build_single_structure`` unless a real PP environment is
# available.

log = logging.getLogger("PlanetProfile.Inference.cache_builder")


# Body orbital parameters used by the k2 forward model. PP's default body
# configs do not set ``Bulk.meanMotion_radps`` or ``Bulk.eccentricity``;
# Test50 injects them inline for Titan. We do the same here for Phase C1
# bodies via a lookup keyed on ``Planet.bodyname`` so the cache build
# works against the unmodified body defaults.
#
# Sources:
#   - Orbital periods from JPL satellite ephemerides (sat427/sat441).
#   - Eccentricities from Murray & Dermott (1999) Table A.1 / NASA fact sheets.
# Recompute meanMotion_radps as 2*pi / (T_days * 86400) when updating.
BODY_ORBITAL_PARAMS: Dict[str, Dict[str, float]] = {
    "Europa":   {"meanMotion_radps": 2.0479e-5, "eccentricity": 0.0094},
    "Ganymede": {"meanMotion_radps": 1.0163e-5, "eccentricity": 0.0013},
    "Callisto": {"meanMotion_radps": 4.358e-6,  "eccentricity": 0.0074},
    "Titan":    {"meanMotion_radps": 4.560e-6,  "eccentricity": 0.0288},
    # Enceladus: Torb 32.9 h (PPEnceladus/Archinal); ecc Jacobson 2006.
    "Enceladus": {"meanMotion_radps": 5.307e-5, "eccentricity": 0.0047},
}


def _save_cache_atomic(data: Any, filepath: str) -> None:
    """Atomic pickle write. Mirrors ``Test50._save_cache``."""
    filepath = str(filepath)
    parent = os.path.dirname(filepath) or "."
    os.makedirs(parent, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=parent, suffix=".pkl.tmp")
    try:
        with os.fdopen(fd, "wb") as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp_path, filepath)
    except BaseException:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise


def build_single_structure(
    planet_template_module: str,
    Tb_K: float,
    ocean_overrides: Dict[str, Any] | None = None,
    bulk_overrides: Dict[str, Any] | None = None,
    planet_overrides: Dict[str, Any] | None = None,
    do_overrides: Dict[str, Any] | None = None,
    extrap_ocean: bool = False,
) -> Dict[str, Any]:
    """Run PlanetProfile once and return the cache dict for one ``Tb_K``.

    Parameters
    ----------
    planet_template_module
        Dotted import path of a PP body config module that exposes a
        ``Planet`` object at module top-level (e.g. ``PlanetProfile.Default.
        Europa.PPEuropa``, ``PlanetProfile.Test.PPTest50``). The function
        deep-copies ``Planet`` and overrides ``Planet.Bulk.Tb_K``; the source
        module is not mutated.
    Tb_K
        Basal-temperature override applied to the deep-copied Planet.
    ocean_overrides
        Optional dict of ``Planet.Ocean.<key>`` overrides applied after
        the deep-copy. Use to switch composition (``{'comp': 'NaCl',
        'wOcean_ppt': 100.0}``) or sweep concentration without
        modifying the body's default template module. Each key/value
        is applied via ``setattr(Planet.Ocean, key, value)``.
    extrap_ocean
        Per-build override for ``Params.EXTRAP_OCEAN`` (ocean/HP-ice EOS
        density extrapolation above the composition's table cap). Default
        ``False`` matches ``defaultConfig``. Threaded explicitly and
        restored in a ``finally`` so a True value for one composition
        (e.g. MgSO4, whose Choukroun-Grasset density table caps at 800 MPa
        while Titan's deep ocean reaches ~1100 MPa) does NOT leak into other
        compositions built in the same process — ``configParams`` is a
        module-level singleton, so a caller-side global mutation would
        contaminate every later build. Set from the config, not the global.

    Returns
    -------
    A dict with keys: ``r_m``, ``rho``, ``K_Pa``, ``mu_Pa``, ``eta_Pa_base``,
    ``phases``, ``T_K``, ``P_MPa``, ``bulk_visc``, ``changeIndices``,
    ``n_layers``, ``layer_upper_radii``, ``layer_types``, ``region_phases``,
    ``omega``, ``eccentricity``, ``host_mass``, ``a_m``, ``R_body_m``,
    ``Mtot_kg``, ``Mtot_achieved_kg``, ``mass_residual_frac``, ``CMR2``,
    ``Tb_K``, ``rhoSil_kgm3``, ``D_hsphere_km``, ``D_iceIh_km``,
    ``D_iceIII_km``, ``D_iceV_km``, ``D_iceVI_km``.

    Mass fields (A4), three distinct quantities that must not be conflated:

    - ``Mtot_kg`` -- the TARGET, ``Planet.Bulk.M_kg`` (the measured body
      mass PP conserves against). Semantics UNCHANGED; it is consumed as
      the conservation target in six places including two ``Test/`` files.
    - ``Mtot_achieved_kg`` -- the mass actually carried by this structure's
      ``(r_m, rho)`` arrays, integrated in PlanetProfile's own OUTER-EDGE
      layer convention (``rho[i]`` fills ``[r[i-1], r[i]]``, implicit inner
      edge 0).
    - ``mass_residual_frac`` -- ``Mtot_achieved_kg / Mtot_kg - 1``. The
      solver's residual is one-signed BY DESIGN (``Silicates.py:118`` takes
      the first profile with ``Mtot_kg <= M_kg``), so a POSITIVE value is a
      distinct bug class, not noise.
    """
    from copy import deepcopy

    from PlanetProfile.GetConfig import Params as configParams
    from PlanetProfile.Gravity.Gravity import SetupGravity
    from PlanetProfile.Main import PlanetProfile as RunPP
    from PlanetProfile.Utilities.defineStructs import Constants, EOSlist
    from PlanetProfile.Utilities.Indexing import PhaseConv

    EOSlist.loaded.clear()

    # Re-import the template module to pick up any in-process edits, then
    # deep-copy its Planet object so this call doesn't pollute the source.
    if planet_template_module in sys.modules:
        importlib.reload(sys.modules[planet_template_module])
    else:
        importlib.import_module(planet_template_module)
    mod = sys.modules[planet_template_module]
    Planet = deepcopy(mod.Planet)

    Planet.Bulk.Tb_K = Tb_K

    # Apply ocean overrides (e.g. composition switch, concentration sweep).
    # Applied to the deep-copied Planet only — source module is unaffected.
    if ocean_overrides:
        for key, value in ocean_overrides.items():
            if not hasattr(Planet.Ocean, key):
                log.warning(
                    f"ocean_overrides: Planet.Ocean has no attribute {key!r}; "
                    "setting it anyway."
                )
            setattr(Planet.Ocean, key, value)
        log.info(f"Applied ocean_overrides: {ocean_overrides}")

    # Bulk overrides (e.g. {'Cuncertainty': 0.06} to widen the MoI matching
    # window for cache builds: small bodies like Enceladus swing MoI hugely
    # per 0.1 K of Tb, so the template's tight measurement window rejects
    # every non-default Tb node — the inference constrains MoI through its
    # CMR2/C20/C22 observables instead, so the cache must ACCEPT a range).
    if bulk_overrides:
        for key, value in bulk_overrides.items():
            if not hasattr(Planet.Bulk, key):
                log.warning(
                    f"bulk_overrides: Planet.Bulk has no attribute {key!r}; "
                    "setting it anyway."
                )
            setattr(Planet.Bulk, key, value)
        log.info(f"Applied bulk_overrides: {bulk_overrides}")

    # Planet-level overrides (top-level Planet attributes, NOT Planet.Ocean or
    # Planet.Bulk). Needed for e.g. PfreezeUpper_MPa (the GetPfreeze ice-Ih
    # transition search ceiling, defined on Planet directly in defineStructs;
    # default 230 MPa), which the cold/fresh NH3 grid corners exceed. Applied
    # to the deep-copied Planet only — source module is unaffected.
    if planet_overrides:
        for key, value in planet_overrides.items():
            if not hasattr(Planet, key):
                log.warning(
                    f"planet_overrides: Planet has no attribute {key!r}; "
                    "setting it anyway."
                )
            setattr(Planet, key, value)
        log.info(f"Applied planet_overrides: {planet_overrides}")

    # Do overrides (Planet.Do switch flags, NOT Planet / Planet.Ocean /
    # Planet.Bulk). Used by the 2D cache builder's frozen-node retry to set
    # ``NO_OCEAN_EXCEPT_INNER_ICES=True`` so a genuinely-frozen (Tb,w) node
    # builds as a real no-ocean structure (HP-ice column, own k2/gravity)
    # instead of raising — the mechanism behind the joint no-ocean+ocean Titan
    # NH3 posterior. Applied to the deep-copied Planet.Do only; the source
    # template module is unaffected.
    if do_overrides:
        for key, value in do_overrides.items():
            if not hasattr(Planet.Do, key):
                log.warning(
                    f"do_overrides: Planet.Do has no attribute {key!r}; "
                    "setting it anyway."
                )
            setattr(Planet.Do, key, value)
        log.info(f"Applied do_overrides: {do_overrides}")

    # Inject orbital params if the PP default config didn't set them.
    # Required by the k2 forward model (rheology path); not used by the v2
    # CMR² mass-conservation path itself.
    if getattr(Planet.Bulk, "meanMotion_radps", None) is None:
        # Synchronous rotators: derive mean motion from the template's
        # orbital period when set (e.g. PPEnceladus defines Torb_s but
        # not meanMotion_radps) before falling back to the table.
        Torb_s = getattr(Planet.Bulk, "Torb_s", None)
        if Torb_s is not None and np.isfinite(Torb_s) and Torb_s > 0:
            Planet.Bulk.meanMotion_radps = 2.0 * np.pi / float(Torb_s)
        body_orb = BODY_ORBITAL_PARAMS.get(getattr(Planet, "bodyname", ""))
        if getattr(Planet.Bulk, "meanMotion_radps", None) is not None:
            pass
        elif body_orb is not None:
            Planet.Bulk.meanMotion_radps = body_orb["meanMotion_radps"]
            if getattr(Planet.Bulk, "eccentricity", None) is None:
                Planet.Bulk.eccentricity = body_orb["eccentricity"]
        else:
            log.warning(
                f"Planet.Bulk.meanMotion_radps is None for body "
                f"{getattr(Planet, 'bodyname', '?')!r} and not in "
                "BODY_ORBITAL_PARAMS. The k2 forward model will fail "
                "for this cache; CMR² inference still works."
            )

    configParams.Gravity.backend = "tidalpy"
    configParams.CALC_NEW = True
    configParams.CALC_NEW_GRAVITY = False
    configParams.NO_SAVEFILE = True
    configParams.SKIP_PLOTS = True

    # EXTRAP_OCEAN is read only during EOS construction inside RunPP/SetupInit,
    # so set it from the per-build parameter and restore the module-level
    # singleton immediately after SetupGravity. try/finally guarantees the
    # restore even if PP raises — otherwise a True value would leak into every
    # subsequent build in this process (reviewer acc725ea R1, 2026-07-27).
    _prior_extrap_ocean = configParams.EXTRAP_OCEAN
    try:
        configParams.EXTRAP_OCEAN = bool(extrap_ocean)
        Planet, Params = RunPP(Planet, configParams)
        Params.CALC_NEW_GRAVITY = True
        Planet, Params = SetupGravity(Planet, Params)
    finally:
        configParams.EXTRAP_OCEAN = _prior_extrap_ocean

    # Surface radius — body-agnostic, pulled from the Planet object that PP
    # populated from the template module.
    R_body_m = float(Planet.Bulk.R_m)

    # ------------------------------------------------------------------
    # The remainder of this function is a verbatim copy of Test50's
    # _build_single_structure body, with TITAN_R_M replaced by R_body_m.
    # ------------------------------------------------------------------
    model = Planet.Gravity.ALMAModel["model"]
    cols = Planet.Gravity.columns
    rIndex = cols.index("r")
    rhoIndex = cols.index("rho")
    VPIndex = cols.index("VP")
    GSIndex = cols.index("GS")
    etaIndex = cols.index("eta")
    pIndex = cols.index("phase")

    r_m = model[:, rIndex].astype(np.float64)
    rho = model[:, rhoIndex].astype(np.float64)
    mu_Pa = model[:, GSIndex].astype(np.float64)
    VP_ms = model[:, VPIndex].astype(np.float64)
    eta_Pa_base = model[:, etaIndex].astype(np.float64)
    phases = model[:, pIndex]

    try:
        T_K_primary = np.asarray(Planet.T_K, dtype=np.float64)
        r_m_primary = np.asarray(Planet.r_m[: T_K_primary.size], dtype=np.float64)
        sort_idx = np.argsort(r_m_primary)
        T_K = np.interp(r_m, r_m_primary[sort_idx], T_K_primary[sort_idx])
    except (AttributeError, ValueError) as _exc:
        T_K = np.full(r_m.size, np.nan)
        log.warning(
            f"Planet.T_K extraction failed ({_exc}) — Arrhenius Ih viscosity will be skipped."
        )

    try:
        P_MPa_primary = np.asarray(Planet.P_MPa[: T_K_primary.size], dtype=np.float64)
        P_MPa = np.interp(r_m, r_m_primary[sort_idx], P_MPa_primary[sort_idx])
    except (AttributeError, ValueError, NameError) as _exc:
        P_MPa = np.full(r_m.size, np.nan)
        log.warning(
            f"Planet.P_MPa extraction failed ({_exc}) — no-ocean safeguard will be disabled."
        )

    K_Pa = rho * VP_ms ** 2 - (4.0 / 3.0) * mu_Pa
    nan_mask = ~np.isfinite(K_Pa) | (K_Pa <= 0)
    if np.any(nan_mask):
        for i in np.where(nan_mask)[0]:
            ph = int(phases[i])
            if 50 <= ph < 100:
                nu = 0.25
            elif ph >= 100:
                nu = 0.29
            else:
                nu = 0.33
            K_Pa[i] = 2.0 * mu_Pa[i] * (1.0 + nu) / (3.0 * (1.0 - 2.0 * nu))
    K_Pa = np.maximum(K_Pa, 1e6)

    changeIndices = np.max(Planet.Reduced.changeIndices) - np.flipud(
        Planet.Reduced.changeIndices
    )
    n_layers = len(changeIndices) - 1

    MIN_POINTS = 5
    needs_padding = any(
        changeIndices[i + 1] - changeIndices[i] < MIN_POINTS
        for i in range(n_layers)
    )

    _orig_iConv = np.flipud(Planet.Reduced.iConv)
    region_phases = []
    for i_layer in range(n_layers):
        start = changeIndices[i_layer]
        phase = phases[start]
        if phase >= Constants.phaseClath and phase < Constants.phaseClath + 10:
            phase = Constants.phaseClath
        convection = _orig_iConv[start]
        phase_str = PhaseConv(phase, liq="0")
        if convection:
            phase_str += "_conv"
        region_phases.append(phase_str)

    bulk_visc = np.zeros_like(eta_Pa_base)

    if needs_padding:
        new_r, new_rho, new_K, new_mu, new_eta, new_phases, new_bv, new_T, new_P = (
            [], [], [], [], [], [], [], [], []
        )
        new_ci = [0]
        for i_layer in range(n_layers):
            s, e = changeIndices[i_layer], changeIndices[i_layer + 1]
            n_pts = e - s
            if n_pts < MIN_POINTS and n_pts >= 2:
                r_layer = r_m[s:e]
                r_interp = np.linspace(r_layer[0], r_layer[-1], MIN_POINTS)
                new_r.append(r_interp)
                new_rho.append(np.interp(r_interp, r_layer, rho[s:e]))
                new_K.append(np.interp(r_interp, r_layer, K_Pa[s:e]))
                new_mu.append(np.interp(r_interp, r_layer, mu_Pa[s:e]))
                new_eta.append(np.interp(r_interp, r_layer, eta_Pa_base[s:e]))
                new_phases.append(np.full(MIN_POINTS, phases[s]))
                new_bv.append(np.zeros(MIN_POINTS))
                new_T.append(np.interp(r_interp, r_layer, T_K[s:e]))
                new_P.append(np.interp(r_interp, r_layer, P_MPa[s:e]))
                new_ci.append(new_ci[-1] + MIN_POINTS)
            elif n_pts == 1:
                r_top = r_m[s]
                r_bot = r_m[s - 1] if s > 0 else r_top - 1.0
                if r_bot >= r_top:
                    r_bot = r_top - 1.0
                r_interp = np.linspace(r_bot, r_top, MIN_POINTS)
                new_r.append(r_interp)
                new_rho.append(np.full(MIN_POINTS, rho[s]))
                new_K.append(np.full(MIN_POINTS, K_Pa[s]))
                new_mu.append(np.full(MIN_POINTS, mu_Pa[s]))
                new_eta.append(np.full(MIN_POINTS, eta_Pa_base[s]))
                new_phases.append(np.full(MIN_POINTS, phases[s]))
                new_bv.append(np.zeros(MIN_POINTS))
                T_lo = T_K[s - 1] if s > 0 and np.isfinite(T_K[s - 1]) else T_K[s]
                T_hi = T_K[s] if np.isfinite(T_K[s]) else T_lo
                new_T.append(np.linspace(T_lo, T_hi, MIN_POINTS))
                P_lo = P_MPa[s - 1] if s > 0 and np.isfinite(P_MPa[s - 1]) else P_MPa[s]
                P_hi = P_MPa[s] if np.isfinite(P_MPa[s]) else P_lo
                new_P.append(np.linspace(P_lo, P_hi, MIN_POINTS))
                new_ci.append(new_ci[-1] + MIN_POINTS)
            else:
                new_r.append(r_m[s:e])
                new_rho.append(rho[s:e])
                new_K.append(K_Pa[s:e])
                new_mu.append(mu_Pa[s:e])
                new_eta.append(eta_Pa_base[s:e])
                new_phases.append(phases[s:e])
                new_bv.append(bulk_visc[s:e])
                new_T.append(T_K[s:e])
                new_P.append(P_MPa[s:e])
                new_ci.append(new_ci[-1] + (e - s))

        r_m = np.concatenate(new_r)
        rho = np.concatenate(new_rho)
        K_Pa = np.concatenate(new_K)
        mu_Pa = np.concatenate(new_mu)
        eta_Pa_base = np.concatenate(new_eta)
        phases = np.concatenate(new_phases)
        bulk_visc = np.concatenate(new_bv)
        P_MPa = np.concatenate(new_P)
        T_K = np.concatenate(new_T)
        changeIndices = np.array(new_ci)

    layer_upper_radii: List[float] = []
    layer_types: List[str] = []
    for i_layer in range(n_layers):
        end = changeIndices[i_layer + 1]
        layer_upper_radii.append(r_m[end - 1])
        layer_types.append(
            "liquid" if phases[changeIndices[i_layer]] == 0 else "solid"
        )

    omega = Planet.Bulk.meanMotion_radps
    ecc = Planet.Bulk.eccentricity
    host_mass = Constants.parentMass_kg[Planet.parent]
    a_m = (Constants.G * host_mass / omega ** 2) ** (1.0 / 3.0)

    D_iceIh_km = D_iceIII_km = D_iceV_km = D_iceVI_km = 0.0
    D_ocean_km_direct = D_clath_km = 0.0
    for i_layer in range(n_layers):
        s = changeIndices[i_layer]
        e = changeIndices[i_layer + 1]
        ph = int(phases[s])
        thick_km = (r_m[e - 1] - r_m[s]) / 1e3
        if ph == 1:
            D_iceIh_km += thick_km
        elif ph == 3:
            D_iceIII_km += thick_km
        elif ph == 5:
            D_iceV_km += thick_km
        elif ph == 6:
            D_iceVI_km += thick_km
        elif ph == 0:
            D_ocean_km_direct += thick_km
        elif ph == Constants.phaseClath:
            D_clath_km += thick_km

    CMR2_pp = Planet.Bulk.Cmeasured if hasattr(Planet.Bulk, "Cmeasured") else np.nan
    try:
        CMR2_pp = float(Planet.CMR2mean)
    except (AttributeError, TypeError):
        pass

    R_sil = getattr(Planet.Sil, "Rmean_m", r_m[0])
    D_hsphere_km = (R_body_m - R_sil) / 1e3

    # C2-A: gravity-coefficient observables. Planet.J2 / Planet.C22 are set
    # by SetupGravity (or the layer-build) and may be None if the body
    # config didn't request a gravity calc. NaN-default keeps the cache
    # writeable; the runner rejects samples with non-finite predictions.
    J2_pred = getattr(Planet, "J2", None)
    if J2_pred is None or not np.isfinite(J2_pred):
        J2_pred = np.nan
    C22_pred = getattr(Planet, "C22", None)
    if C22_pred is None or not np.isfinite(C22_pred):
        C22_pred = np.nan

    # C2-B: electrical-conductivity profile for magnetic induction.
    # Planet.sigma_Sm (S/m) is per-cell on the original Reduced.r grid.
    # Resample it onto the cache's r_m grid (which has the layer-padding
    # applied) by linear interpolation. NaN cells in input are left as NaN
    # in output — induction code skips them.
    try:
        sigma_full = np.asarray(Planet.sigma_Sm, dtype=np.float64)
        if sigma_full.ndim != 1:
            raise ValueError("Planet.sigma_Sm is not 1-D")
        # Reduced.r is the same grid Reduced.changeIndices indexes into;
        # pre-padding it matches r_m_primary (which we built earlier from
        # Planet.r_m).  We interpolate against r_m_primary if it's available
        # in scope; otherwise use a uniform fallback.
        try:
            sigma_primary = sigma_full[: T_K_primary.size]
            sigma_Sm = np.interp(r_m, r_m_primary[sort_idx], sigma_primary[sort_idx])
        except (NameError, ValueError):
            sigma_Sm = np.full(r_m.size, np.nan)
    except (AttributeError, ValueError, TypeError) as _exc:
        log.debug(f"sigma_Sm unavailable ({_exc}); cache field set to NaN.")
        sigma_Sm = np.full(r_m.size, np.nan)

    # C2-B: layered representation for the induction forward model.
    # Planet.Magnetic.{rSigChange_m, sigmaLayers_Sm} are short arrays
    # (~5–10 entries, one per σ region). Planet.Magnetic.Texc_hr is a
    # numpy array of *periods only* — the canonical labels live in the
    # body-keyed lookup table in MagneticInduction.Moments.ExcitationsList.
    # We zip the two together below so the cache exposes a label→period
    # dict that the induction forward model can consume directly.
    # Note: PP skips MagneticInduction entirely when Ocean.comp == 'PureH2O'
    # (see PlanetProfile/Main.py:352), so for pure-water bodies these
    # fields will all be None — that is expected, not a bug.
    rSigChange_m = None
    sigmaLayers_Sm = None
    Texc_hr = None
    try:
        mag = Planet.Magnetic
        if mag.rSigChange_m is not None:
            rSigChange_m = np.asarray(mag.rSigChange_m, dtype=np.float64)
        if mag.sigmaLayers_Sm is not None:
            sigmaLayers_Sm = np.asarray(mag.sigmaLayers_Sm, dtype=np.float64)
        if mag.Texc_hr is not None:
            periods = np.asarray(mag.Texc_hr, dtype=np.float64).ravel()
            try:
                from PlanetProfile.MagneticInduction.Moments import (
                    ExcitationsList,
                )
                body_table = ExcitationsList().Texc_hr.get(
                    getattr(Planet, "bodyname", None), {}
                )
                # Match each cached period to its closest labelled period.
                Texc_hr = {}
                for p in periods:
                    if not np.isfinite(p):
                        continue
                    best_label = None
                    best_diff = np.inf
                    for label, ref_p in body_table.items():
                        if ref_p is None:
                            continue
                        d = abs(float(ref_p) - float(p))
                        if d < best_diff:
                            best_diff = d
                            best_label = label
                    key = best_label if best_label is not None else f"period_{p:.4f}hr"
                    Texc_hr[str(key)] = float(p)
            except (ImportError, AttributeError, TypeError) as _label_exc:
                log.debug(
                    f"Could not label Texc_hr periods ({_label_exc}); "
                    "falling back to numeric keys."
                )
                Texc_hr = {f"period_{p:.4f}hr": float(p)
                           for p in periods if np.isfinite(p)}
    except (AttributeError, TypeError) as _exc:
        log.debug(f"Magnetic substruct unavailable ({_exc}); induction "
                  "fields set to None.")

    # Phase-code → label mapping required by compute_heating (forward_models.py).
    # Same hardcoded mapping as PlanetProfileApp's extract_structure_from_planet
    # (structure_cache.py). PP phase codes: 0=ocean, 1=Ih, 2=II, 3=III, 5=V,
    # 6=VI; PhaseConv() in Electrical.py is the source of truth, but we keep
    # this dict here because compute_heating accesses it positionally and a
    # round-trip through PhaseConv would slow the per-sample heating loop.
    phase_map = {0: '0', 1: 'Ih', 2: 'II', 3: 'III', 5: 'V', 6: 'VI'}

    # ------------------------------------------------------------------
    # ACHIEVED mass (A4). ``Mtot_kg`` below is the TARGET (Planet.Bulk.M_kg,
    # the measured body mass PP conserves against) and its semantics are
    # DELIBERATELY UNCHANGED -- it is consumed as the conservation target in
    # six places including two Test/ files, so it must not be repurposed
    # (RESUME_NOTE "Not blocking the build, but owed"). The achieved mass is
    # published alongside it under distinct names so the mass-conservation
    # node invariant can be evaluated without re-integrating the profile at
    # every consumer, and so a violation can never again be MASKED by the
    # stored target (the defect behind MAJOR-1: the smoke cache's frozen
    # nodes stored the template Mtot_kg while being -22.16% / -9.21% off).
    #
    # Quadrature is PlanetProfile's OWN layer convention -- OUTER-EDGE:
    # rho[i] fills the shell [r[i-1], r[i]] with an implicit inner edge of 0
    # (reducedPlanetModel / Geophysical.py:170). Named here because the
    # convention swing across quadratures (~2% on one Titan node) is larger
    # than every effect being measured, and because "an audit of existing
    # data must adopt the convention the data was built with"
    # (mass_conservation_revision_history v3).
    _m_order = np.argsort(r_m)
    _r_sorted = r_m[_m_order]
    _rho_sorted = rho[_m_order]
    _r_inner = np.concatenate(([0.0], _r_sorted[:-1]))
    Mtot_achieved_kg = float(np.sum(
        _rho_sorted * (4.0 / 3.0) * np.pi
        * (_r_sorted ** 3 - _r_inner ** 3)))
    _M_target = float(Planet.Bulk.M_kg)
    mass_residual_frac = (
        float(Mtot_achieved_kg / _M_target - 1.0)
        if np.isfinite(_M_target) and _M_target > 0 else float("nan"))

    return {
        "r_m": np.ascontiguousarray(r_m),
        "rho": np.ascontiguousarray(rho),
        "K_Pa": np.ascontiguousarray(K_Pa),
        "mu_Pa": np.ascontiguousarray(mu_Pa),
        "eta_Pa_base": eta_Pa_base,
        "phases": phases,
        "T_K": np.ascontiguousarray(T_K),
        "P_MPa": np.ascontiguousarray(P_MPa),
        "bulk_visc": np.ascontiguousarray(bulk_visc),
        "sigma_Sm": np.ascontiguousarray(sigma_Sm),
        "changeIndices": changeIndices,
        "n_layers": n_layers,
        "layer_upper_radii": tuple(layer_upper_radii),
        "layer_types": tuple(layer_types),
        "region_phases": region_phases,
        "phase_map": phase_map,
        "omega": omega,
        "eccentricity": ecc,
        "host_mass": host_mass,
        "a_m": a_m,
        "R_body_m": R_body_m,
        # TARGET mass (unchanged semantics -- see the A4 note above).
        "Mtot_kg": Planet.Bulk.M_kg,
        # ACHIEVED mass by outer-edge integration of this structure's own
        # (r_m, rho), and its one-sided fractional residual against the
        # target. New in A4; never overwrite Mtot_kg with these.
        "Mtot_achieved_kg": Mtot_achieved_kg,
        "mass_residual_frac": mass_residual_frac,
        "CMR2": CMR2_pp,
        "J2": float(J2_pred),
        "C22": float(C22_pred),
        # Read back from Planet.Bulk (post-run) rather than echoing the
        # caller's input ``Tb_K`` unconditionally. Byte-identical to the old
        # behaviour for every Tb-driven build (Do.ICEIh_THICKNESS=False):
        # Planet.Bulk.Tb_K is never mutated on that path (verified against
        # SetupInit.py, which only branches on Tb_K when ICEIh_THICKNESS is
        # False, and LayerPropagators.py, which never reassigns it there).
        # For a zb-driven build (Do.ICEIh_THICKNESS=True, Bulk.
        # zb_approximate_km set -- B4 cache mode), GetIceShellTFreeze SOLVES
        # for Tb_K and reassigns Planet.Bulk.Tb_K to the root-found value
        # (LayerPropagators.py:53, ``Planet = GetIceShellTFreeze(Planet,
        # Params)``); the caller's input Tb_K is then just an unused
        # placeholder, and this field must read the solved value back.
        "Tb_K": float(getattr(Planet.Bulk, "Tb_K", Tb_K)),
        # Ocean salinity actually realized in this structure (g/kg). For v3
        # (Tb × w) caches this is the second grid axis; for 1D v2 caches it
        # simply records the fixed seawater value baked in at build time.
        # Read back from Planet.Ocean so an ocean_overrides={'wOcean_ppt': w}
        # is reflected faithfully (defineStructs sets it via 2-of-3 pairing).
        "wOcean_ppt": float(getattr(Planet.Ocean, "wOcean_ppt", np.nan)),
        "rhoSil_kgm3": getattr(Planet.Sil, "rhoMean_kgm3", np.nan),
        "D_hsphere_km": D_hsphere_km,
        "D_iceIh_km": D_iceIh_km,
        "D_iceIII_km": D_iceIII_km,
        "D_iceV_km": D_iceV_km,
        "D_iceVI_km": D_iceVI_km,
        # Clathrate shell thickness (phase code Constants.phaseClath). Summed
        # directly from clathrate layers; 0.0 when the model has none.
        "D_clath_km": D_clath_km,
        # Ocean thickness = SUM of actual liquid (phase-0) layer thicknesses,
        # matching find_tb_bounds.py:156-157 and hybrid_structure_cache.py.
        # Previously this was D_hsphere − Σ(4 ice phases), which (a) omitted the
        # clathrate shell and (b) dropped the ~few-km radial gaps between
        # adjacent layers' boundary nodes — so a FULLY FROZEN hydrosphere (no
        # phase-0 layer) reported a spurious ~16 km "ocean" (scientific-reviewer
        # 2026-07-24: verified the residual = clathrate 1.9 km + 14.8 km node
        # gaps, NOT liquid). Summing liquid layers gives a true 0 for no-ocean
        # models and the correct thickness for ocean-bearing ones. This is a
        # real computed 0, NOT the pt.get('D_ocean_km', 0.0) missing-key
        # fallback that caused the Europa "no-ocean 100%" panel — the key is
        # still present, so mcmc_plots (855/1069) reads a valid value.
        "D_ocean_km": D_ocean_km_direct,
        # C2-B induction layered representation
        "rSigChange_m": rSigChange_m,
        "sigmaLayers_Sm": sigmaLayers_Sm,
        "Texc_hr": Texc_hr,
    }


def _regions_match(s0: Dict[str, Any], s1: Dict[str, Any]) -> bool:
    """True if two structures share the same layer set (region_phases + n_layers)."""
    return (
        s0.get("region_phases") == s1.get("region_phases")
        and s0.get("n_layers") == s1.get("n_layers")
    )


def _bisect_transition(
    planet_template_module: str,
    Tb_lo: float,
    Tb_hi: float,
    regions_lo: Any,
    regions_hi: Any,
    eps_T: float = 0.01,
    max_iter: int = 20,
    ocean_overrides: Dict[str, Any] | None = None,
    bulk_overrides: Dict[str, Any] | None = None,
    planet_overrides: Dict[str, Any] | None = None,
    do_overrides: Dict[str, Any] | None = None,
    extrap_ocean: bool = False,
) -> List[Tuple[float, Dict[str, Any]]]:
    """Bisect ``[Tb_lo, Tb_hi]`` to localise a layer-set transition.

    The interval is shrunk until ``Tb_hi - Tb_lo < eps_T``. Each bisection
    midpoint runs PP (via :func:`build_single_structure`) and is appended
    to the returned list of ``(Tb, struct)`` pairs — the caller inserts
    these into the overall grid, rather than throwing them away. The
    converged anchor pair (final Tb_lo, final Tb_hi) bracket the actual
    transition by less than ``eps_T``.

    If a third layer set is discovered at a midpoint (multiple transitions
    in the original interval), the function recurses into both
    sub-intervals; the third regime's structure is also added to the
    returned list.

    Failures of the underlying PP build at a midpoint terminate refinement
    early — the partially-narrowed interval is returned (with whatever
    midpoints succeeded), and the caller still has the original endpoints
    so coverage is preserved.
    """
    discovered: List[Tuple[float, Dict[str, Any]]] = []
    iters = 0

    while (Tb_hi - Tb_lo) > eps_T and iters < max_iter:
        Tb_mid = 0.5 * (Tb_lo + Tb_hi)
        try:
            s_mid = build_single_structure(
                planet_template_module, float(Tb_mid),
                ocean_overrides=ocean_overrides,
                bulk_overrides=bulk_overrides,
                planet_overrides=planet_overrides,
                do_overrides=do_overrides,
                extrap_ocean=extrap_ocean,
            )
        except Exception as exc:
            log.warning(
                f"    × bisection midpoint Tb_K={Tb_mid:.4f} failed "
                f"({type(exc).__name__}: {exc}); stopping refinement of this interval."
            )
            break

        regions_mid = s_mid.get("region_phases")
        if regions_mid == regions_lo:
            Tb_lo = Tb_mid
            discovered.append((Tb_mid, s_mid))
        elif regions_mid == regions_hi:
            Tb_hi = Tb_mid
            discovered.append((Tb_mid, s_mid))
        else:
            log.info(
                f"    third regime at Tb_K={Tb_mid:.4f}: {regions_mid}; "
                "splitting bisection."
            )
            below = _bisect_transition(
                planet_template_module, Tb_lo, Tb_mid,
                regions_lo, regions_mid, eps_T, max_iter,
                ocean_overrides=ocean_overrides,
                bulk_overrides=bulk_overrides,
                planet_overrides=planet_overrides,
                do_overrides=do_overrides,
                extrap_ocean=extrap_ocean,
            )
            above = _bisect_transition(
                planet_template_module, Tb_mid, Tb_hi,
                regions_mid, regions_hi, eps_T, max_iter,
                ocean_overrides=ocean_overrides,
                bulk_overrides=bulk_overrides,
                planet_overrides=planet_overrides,
                do_overrides=do_overrides,
                extrap_ocean=extrap_ocean,
            )
            discovered.append((Tb_mid, s_mid))
            discovered.extend(below)
            discovered.extend(above)
            break
        iters += 1

    return discovered


def build_tb_grid_cache(
    planet_template_module: str,
    Tb_grid: Iterable[float],
    output_path: str,
    progress: bool = True,
    detect_transitions: bool = True,
    eps_T: float = 0.01,
    ocean_overrides: Dict[str, Any] | None = None,
    bulk_overrides: Dict[str, Any] | None = None,
    planet_overrides: Dict[str, Any] | None = None,
    do_overrides: Dict[str, Any] | None = None,
    extrap_ocean: bool = False,
) -> Dict[str, Any]:
    """Build a Tb-grid cache by repeatedly calling ``build_single_structure``.

    Two-phase build:

    1. **Regular grid.** Run PP at each ``Tb`` in ``Tb_grid`` (sorted).
       Skip-on-failure preserves coverage.
    2. **Transition refinement.** When ``detect_transitions`` is True
       (default), walk adjacent grid pairs and detect changes in
       ``region_phases``. For each transition, bisect to within
       ``eps_T`` Kelvin (default 0.01). The bisection midpoints are
       added to the grid as extra structures, and the converged
       endpoints flank the transition by < ``eps_T``. The returned
       cache also carries a ``transitions`` metadata list so the
       MCMC runner can dispatch nearest-neighbour across them
       instead of attempting a layered blend.

    Parameters
    ----------
    planet_template_module
        See :func:`build_single_structure`.
    Tb_grid
        Iterable of Tb_K values (K). Cast to a 1-D numpy float array
        and sorted ascending.
    output_path
        Pickle path. Atomic write via tempfile + os.replace.
    progress
        Per-grid-point INFO log lines (default True).
    detect_transitions
        If True (default), refine the grid near layer-set transitions.
        Set to False to reproduce the unrefined v2.0 cache shape.
    eps_T
        Bisection tolerance in Kelvin (default 0.01). The anchor pair
        flanking each transition is separated by < ``eps_T``.

    Returns
    -------
    The saved-to-disk dict::

        {'Tb_K_grid': np.ndarray,
         'structures': [structure_dict, ...],
         'transitions': [{'Tb_lo': float, 'Tb_hi': float,
                          'regions_lo': [...], 'regions_hi': [...]},
                         ...],
         'schema_version': 'v2.1'}

    The ``structures`` list is parallel to ``Tb_K_grid``; index ``i`` of
    one matches index ``i`` of the other. ``transitions`` is empty when
    no transitions were detected (or when ``detect_transitions`` is
    False).
    """
    Tb_arr = np.asarray(list(Tb_grid), dtype=np.float64)
    if Tb_arr.ndim != 1 or Tb_arr.size < 2:
        raise ValueError(
            f"Tb_grid must be a 1-D iterable with >= 2 points, got {Tb_arr!r}"
        )
    order = np.argsort(Tb_arr)
    Tb_arr = Tb_arr[order]

    # ----- Phase 1: regular grid -----
    structures: List[Dict[str, Any]] = []
    accepted_Tb: List[float] = []
    skipped: List[Tuple[float, str]] = []
    for i, Tb_K in enumerate(Tb_arr):
        if progress:
            log.info(
                f"[{i + 1}/{Tb_arr.size}] Building structure at Tb_K={Tb_K:.4f}"
            )
        try:
            struct = build_single_structure(
                planet_template_module, float(Tb_K),
                ocean_overrides=ocean_overrides,
                bulk_overrides=bulk_overrides,
                planet_overrides=planet_overrides,
                do_overrides=do_overrides,
                extrap_ocean=extrap_ocean,
            )
        except Exception as exc:
            log.warning(
                f"    × Tb_K={Tb_K:.4f} skipped — {type(exc).__name__}: {exc}"
            )
            skipped.append((float(Tb_K), f"{type(exc).__name__}: {exc}"))
            continue
        structures.append(struct)
        accepted_Tb.append(float(Tb_K))
        if progress:
            log.info(
                f"    → CMR²={struct.get('CMR2', float('nan')):.4f}, "
                f"D_iceIh={struct.get('D_iceIh_km', 0.0):.1f} km, "
                f"D_hsphere={struct.get('D_hsphere_km', 0.0):.1f} km"
            )

    if not structures:
        raise RuntimeError(
            f"All {Tb_arr.size} Tb points failed; no cache to write. "
            f"Skipped: {skipped!r}"
        )
    if skipped and progress:
        log.warning(
            f"Regular grid built with gaps: {len(skipped)}/{Tb_arr.size} "
            f"Tb points skipped — {[f'{t:.2f}' for t, _ in skipped]}"
        )

    # ----- Phase 2: transition refinement -----
    extra_Tb: List[float] = []
    extra_structs: List[Dict[str, Any]] = []
    if detect_transitions:
        for i in range(len(structures) - 1):
            s0, s1 = structures[i], structures[i + 1]
            if _regions_match(s0, s1):
                continue
            Tb_a, Tb_b = accepted_Tb[i], accepted_Tb[i + 1]
            if progress:
                log.info(
                    f"  Transition in [{Tb_a:.4f}, {Tb_b:.4f}]: "
                    f"{s0.get('region_phases')} → {s1.get('region_phases')}; "
                    f"refining to ε_T={eps_T:.4f} K"
                )
            new_points = _bisect_transition(
                planet_template_module, Tb_a, Tb_b,
                s0.get("region_phases"), s1.get("region_phases"),
                eps_T=eps_T,
                ocean_overrides=ocean_overrides,
                bulk_overrides=bulk_overrides,
                planet_overrides=planet_overrides,
                do_overrides=do_overrides,
                extrap_ocean=extrap_ocean,
            )
            for Tb_new, s_new in new_points:
                extra_Tb.append(float(Tb_new))
                extra_structs.append(s_new)

    # ----- Merge + sort -----
    all_Tb = accepted_Tb + extra_Tb
    all_structs = structures + extra_structs
    order2 = np.argsort(all_Tb)
    final_Tb = np.asarray([all_Tb[i] for i in order2], dtype=np.float64)
    final_structs = [all_structs[i] for i in order2]

    # ----- Build transitions metadata from the final, refined grid -----
    transitions: List[Dict[str, Any]] = []
    for i in range(len(final_structs) - 1):
        s0, s1 = final_structs[i], final_structs[i + 1]
        if not _regions_match(s0, s1):
            transitions.append({
                "Tb_lo": float(final_Tb[i]),
                "Tb_hi": float(final_Tb[i + 1]),
                "regions_lo": list(s0.get("region_phases", [])),
                "regions_hi": list(s1.get("region_phases", [])),
            })

    cache = {
        "Tb_K_grid": final_Tb,
        "structures": final_structs,
        "transitions": transitions,
        "schema_version": "v2.1",
    }
    _save_cache_atomic(cache, output_path)
    if progress:
        n_extra = len(extra_Tb)
        log.info(
            f"Cache written → {output_path} "
            f"({len(final_structs)} grid points: "
            f"{Tb_arr.size} regular + {n_extra} from transition refinement; "
            f"{len(transitions)} transition(s))"
        )
    return cache


def build_tbw_grid_cache(
    planet_template_module: str,
    Tb_grid: Iterable[float],
    wOcean_ppt_grid: Iterable[float],
    output_path: str,
    progress: bool = True,
    ocean_overrides: Dict[str, Any] | None = None,
    bulk_overrides: Dict[str, Any] | None = None,
    planet_overrides: Dict[str, Any] | None = None,
    do_overrides: Dict[str, Any] | None = None,
    extrap_ocean: bool = False,
    retry_frozen_as_no_ocean: bool = False,
) -> Dict[str, Any]:
    """Build a 2D (Tb × wOcean_ppt) structure-grid cache (schema ``v3.0``).

    Europa Clipper v3 samples ocean salinity ``log10_wOcean_ppt`` as an 8th
    parameter, so the 1D-in-Tb v2.1 cache becomes a 2D grid. Salinity is baked
    into each cached structure at build time (via ``ocean_overrides=
    {'wOcean_ppt': w}`` → ``Planet.Ocean.wOcean_ppt`` → the GSW/Seawater EOS →
    the per-cell ``sigma_Sm`` conductivity profile that drives all induction
    channels); NO runtime hook recomputes ocean EOS.

    Unlike :func:`build_tb_grid_cache`, this builder uses a **fixed regular
    grid with no transition refinement** — per-salinity refinement would place
    different Tb nodes at each ``w`` and could not be stacked row-major. The
    regular-grid choice is validated in Phase-1 pre-commit checks (max |ΔCMR²|
    vs the refined v2 cache < 0.25σ; ocean-onset & |Ae_synodic|=0.7 edges
    agree within one grid step) — see plans/europa-clipper-v3-salinity-plan.md.
    HP ices (III/V/VI) are structurally impossible on the Europa GSW-seawater
    path (basal ocean P ~156 MPa < the 200 MPa ``PminHPices_MPa`` cap), so the
    only discontinuity a regular grid must resolve is ocean-onset (~260 K),
    which coincides with the induction support edge.

    Parameters
    ----------
    planet_template_module
        See :func:`build_single_structure`.
    Tb_grid
        Iterable of basal-temperature values (K). Cast to 1-D float, sorted
        ascending. Use a regular spacing (e.g. 0.125 K).
    wOcean_ppt_grid
        Iterable of ocean salinities (g/kg). Cast to 1-D float, sorted
        ascending. Use log-spaced nodes (e.g. ``np.logspace(-1, 2, 16)`` for
        0.1–100 ppt) so the low-w physics (synodic de-saturation, the
        |Ae_synodic|=0.7 support contour) is resolved.
    output_path
        Pickle path. Atomic write via tempfile + ``os.replace``.
    progress
        Per-node INFO log lines (default True).
    ocean_overrides
        Optional ``Planet.Ocean.<key>`` overrides applied to EVERY node in
        addition to the per-node ``wOcean_ppt``. The one field that MUST come
        through here for salted-ocean caches is ``comp`` (e.g. ``'NaCl'``,
        ``'NH3'``, ``'MgSO4'``): the composition is baked into each cached
        structure at build time (it dispatches the ocean EOS in
        ``HydroEOS.py``), and there is NO runtime hook to recompute it. If
        ``comp`` is omitted the cache inherits the template module's
        ``Planet.Ocean.comp`` (so e.g. a ``PureH2O`` template yields a
        pure-water cache regardless of ``wOcean_ppt``). ``wOcean_ppt`` in this
        dict is ignored — the salinity axis is always set per-node from
        ``wOcean_ppt_grid``.

    Returns
    -------
    The saved-to-disk dict::

        {'Tb_K_grid': np.ndarray (n_Tb,),
         'wOcean_ppt_grid': np.ndarray (n_w,),
         'structures': [structure_dict | None, ...],  # row-major, len n_Tb*n_w
         'ocean_comp': str | None,  # baked-in composition (from ocean_overrides)
         'schema_version': 'v3.0'}

    ``structures`` is row-major: entry ``i_Tb * n_w + i_w`` is the structure at
    ``(Tb_K_grid[i_Tb], wOcean_ppt_grid[i_w])``. Failed/rejected builds are
    stored as explicit ``None`` (never silently dropped) so the 2D lookup can
    reject those (Tb, w) corners. Each structure carries both ``Tb_K`` and
    ``wOcean_ppt`` plus a ``has_ocean`` bool. The absence of ``wOcean_ppt_grid``
    in a cache signals the legacy 1D path, so v1/v2 artifacts stay servable.

    retry_frozen_as_no_ocean
        When True (used by the joint no-ocean+ocean Titan NH3 campaign), a node
        that raises :class:`NoIceLiquidTransitionError` (genuinely frozen: no
        ice-Ih->liquid transition for this Tb/w) is retried ONCE as a real
        no-ocean structure via ``do_overrides={'NO_OCEAN_EXCEPT_INNER_ICES':
        True}`` and stored (tagged ``has_ocean=False``) instead of ``None``.
        This is the ONLY failure that triggers the retry — thin-shell /
        PHydroMax / EOS-extrapolation failures (and any retry that itself
        fails) still store ``None``. Default False preserves the ocean-only
        behavior of every existing 2D cache (Europa v3/v5/v6/v7). The cache
        dict then also carries ``retry_frozen_as_no_ocean`` and
        ``n_no_ocean_nodes``.
    """
    # Deferred (defineStructs pulls in matplotlib/cmasher); only needed when the
    # frozen-node retry is active, but cheap to import unconditionally here.
    from PlanetProfile.Utilities.defineStructs import NoIceLiquidTransitionError

    Tb_arr = np.asarray(list(Tb_grid), dtype=np.float64)
    w_arr = np.asarray(list(wOcean_ppt_grid), dtype=np.float64)
    if Tb_arr.ndim != 1 or Tb_arr.size < 2:
        raise ValueError(
            f"Tb_grid must be a 1-D iterable with >= 2 points, got {Tb_arr!r}"
        )
    if w_arr.ndim != 1 or w_arr.size < 2:
        raise ValueError(
            f"wOcean_ppt_grid must be a 1-D iterable with >= 2 points, "
            f"got {w_arr!r}"
        )
    Tb_arr = Tb_arr[np.argsort(Tb_arr)]
    w_arr = w_arr[np.argsort(w_arr)]

    # Composition (and any other Ocean.* overrides) applied to every node.
    # ``wOcean_ppt`` is set per-node below, so drop it if present here to
    # avoid a stale scalar overriding the salinity axis.
    base_overrides = dict(ocean_overrides or {})
    base_overrides.pop("wOcean_ppt", None)
    ocean_comp = base_overrides.get("comp")
    if progress and ocean_comp is not None:
        log.info(f"Ocean composition baked into every node: comp={ocean_comp!r}")

    n_Tb, n_w = Tb_arr.size, w_arr.size
    total = n_Tb * n_w
    structures: List[Dict[str, Any] | None] = [None] * total
    n_ok = 0
    n_fail = 0
    for i_Tb, Tb_K in enumerate(Tb_arr):
        for i_w, w_ppt in enumerate(w_arr):
            flat = i_Tb * n_w + i_w
            if progress:
                log.info(
                    f"[{flat + 1}/{total}] Tb_K={Tb_K:.4f} "
                    f"wOcean_ppt={w_ppt:.4f}"
                )
            try:
                struct = build_single_structure(
                    planet_template_module, float(Tb_K),
                    ocean_overrides={**base_overrides,
                                     "wOcean_ppt": float(w_ppt)},
                    bulk_overrides=bulk_overrides,
                    planet_overrides=planet_overrides,
                    do_overrides=do_overrides,
                    extrap_ocean=extrap_ocean,
                )
                struct["has_ocean"] = True
            except NoIceLiquidTransitionError as frozen_exc:
                # Genuinely-frozen node: Tb is too cold (given this w's
                # freezing-point depression) for a liquid ocean to exist in the
                # Pfreeze window — GetPfreeze found no ice-Ih->liquid transition
                # (typed subclass; NOT a thin-shell / PHydroMax / EOS failure).
                # For the joint no-ocean+ocean posterior, retry this node ONCE
                # as a real no-ocean structure (NO_OCEAN_EXCEPT_INNER_ICES: the
                # "ocean" cells become high-pressure ice; k2/gravity computed on
                # the frozen column) so it enters the cache instead of being
                # stored as None (which the 2D interpolant would reject,
                # implicitly conditioning the posterior on "an ocean exists").
                # comp is kept (e.g. NH3): the NH3-depressed liquidus still sets
                # the HP-ice temperature profile via GetTfreeze, so these frozen
                # structures are NOT identical to a PureH2O no-ocean build.
                if not retry_frozen_as_no_ocean:
                    log.warning(
                        f"    × (Tb={Tb_K:.4f}, w={w_ppt:.4f}) → None — "
                        f"NoIceLiquidTransitionError: {frozen_exc}"
                    )
                    structures[flat] = None
                    n_fail += 1
                    continue
                try:
                    struct = build_single_structure(
                        planet_template_module, float(Tb_K),
                        ocean_overrides={**base_overrides,
                                         "wOcean_ppt": float(w_ppt)},
                        bulk_overrides=bulk_overrides,
                        planet_overrides=planet_overrides,
                        do_overrides={**(do_overrides or {}),
                                      "NO_OCEAN_EXCEPT_INNER_ICES": True},
                        extrap_ocean=extrap_ocean,
                    )
                    struct["has_ocean"] = False
                    if progress:
                        log.info(
                            f"    ↻ (Tb={Tb_K:.4f}, w={w_ppt:.4f}) frozen → "
                            f"rebuilt as no-ocean (NO_OCEAN_EXCEPT_INNER_ICES)"
                        )
                except Exception as retry_exc:
                    log.warning(
                        f"    × (Tb={Tb_K:.4f}, w={w_ppt:.4f}) frozen retry "
                        f"failed → None — "
                        f"{type(retry_exc).__name__}: {retry_exc}"
                    )
                    structures[flat] = None
                    n_fail += 1
                    continue
            except Exception as exc:
                log.warning(
                    f"    × (Tb={Tb_K:.4f}, w={w_ppt:.4f}) → None — "
                    f"{type(exc).__name__}: {exc}"
                )
                structures[flat] = None
                n_fail += 1
                continue
            structures[flat] = struct
            n_ok += 1
            if progress:
                log.info(
                    f"    → CMR²={struct.get('CMR2', float('nan')):.4f}, "
                    f"D_ocean={struct.get('D_ocean_km', 0.0):.1f} km, "
                    f"D_iceIh={struct.get('D_iceIh_km', 0.0):.1f} km, "
                    f"has_ocean={struct.get('has_ocean')}"
                )

    if n_ok == 0:
        raise RuntimeError(
            f"All {total} (Tb, w) nodes failed; no cache to write."
        )

    n_no_ocean = sum(
        1 for s in structures if s is not None and s.get("has_ocean") is False
    )
    cache = {
        "Tb_K_grid": Tb_arr,
        "wOcean_ppt_grid": w_arr,
        "structures": structures,
        "ocean_comp": ocean_comp,
        "schema_version": "v3.0",
        # Joint no-ocean+ocean metadata. When retry_frozen_as_no_ocean is on,
        # genuinely-frozen (low-w, cold-Tb) nodes are rebuilt as no-ocean
        # structures (tagged has_ocean=False) and enter the cache instead of
        # being stored as None — so the amortized posterior spans BOTH frozen
        # and ocean Titan interiors rather than being implicitly conditioned on
        # "an ocean exists". The frozen branch models NO residual brine (the
        # eutectic thin-brine regime is out of scope; ices are pure H2O with
        # NH3 rejected to the — here absent — liquid). Each structure carries a
        # per-node ``has_ocean`` bool.
        "retry_frozen_as_no_ocean": bool(retry_frozen_as_no_ocean),
        "n_no_ocean_nodes": int(n_no_ocean),
    }
    _save_cache_atomic(cache, output_path)
    if progress:
        log.info(
            f"2D cache written → {output_path} "
            f"({n_Tb} Tb × {n_w} w = {total} nodes: "
            f"{n_ok} built, {n_fail} None)"
        )
    return cache


# ---------------------------------------------------------------------------
# Frozen-branch (no-ocean) node builder -- frozen-branch design ruling, A3
# ---------------------------------------------------------------------------
#
# BUILD MECHANISM (established empirically, A1 completion; see the
# build_zbw_grid_cache docstring for the full argument):
#
#   Do.ConstantProps['Inner'] = True   -> uniform-density interior, so the
#     interior density is solved ANALYTICALLY from mass closure rather than
#     searched (LayerPropagators.FindInnerWithMoIAndConstantRho). This is
#     ruling F2's "analytically exact mass closure" and it is what makes the
#     mass residual land at machine precision instead of one radial step.
#   Do.Fe_CORE = False, Do.POROUS_ROCK = False -> single uniform rock sphere,
#     matching the config's declared uniform-interior assumption.
#   Do.SPECIFY_ICEI_BOTTOM_PRESSURE = True + Bulk.PbISet_MPa -> the ice-Ih
#     base pressure is SET, and Tb_K is derived from the melt curve. This is
#     the frozen branch's zb control: zb(PbI) is smooth and monotone
#     (~8.3 km/MPa at Enceladus), so the builder root-finds PbI to land on
#     each grid node. Driving zb through Do.ICEIh_THICKNESS instead was
#     measured to miss the target by ~2.1 km -- an order of magnitude past
#     invariant I-F4's 0.25 km -- so it is not used here.
#   Do.HYDROSPHERE_THICKNESS = True + Bulk.Dhsphere_m = 0.0 -> the seafloor
#     is chosen as the SHALLOWEST admissible candidate. Steps.iSilStart
#     pins that sweep to start at the ice-Ih base, so "shallowest" IS the ice
#     base: the whole hydrosphere is ice and no liquid layer exists.
#     Measured: 0 phase-0 layers at every probed node.
#
# WHY THIS MAKES THE MoI WINDOW NON-CONDITIONING (ruling 3.2). The window
# still filters the candidate list, but the selection among survivors is
# argmin|thickness - 0|, i.e. always the shallowest. So the window cannot
# move the chosen structure to a different admissible one -- it can only
# make the ice-base candidate inadmissible, in which case a deeper seafloor
# is selected, liquid appears, and invariant I-F2 HARD REJECTS the node.
# The window is therefore incapable of silently conditioning a SURVIVING
# node. See build_zbw_grid_cache's frozen_Cuncertainty argument.

# ---------------------------------------------------------------------------
# D1 -- the OCEAN-branch MoI acceptance window, read from the config.
# ---------------------------------------------------------------------------
# r5 build blocker D1: the ocean-branch MoI window was "unspecified and
# unenforceable" -- ``bulk_overrides`` was not a config field at all, so
# :func:`build_zbw_grid_cache` fell through to the Enceladus template's own
# ``Bulk.Cuncertainty = 0.001`` (PPEnceladus.py:23) no matter what the config
# asserted. This helper makes the window a DECLARED, READ, ENFORCED config
# quantity.
#
# WHAT THE WINDOW DOES ON THE OCEAN BRANCH (the structural half of the
# derivation; mirrors the frozen branch's ruled non-conditioning argument in
# the module comment above, but the mechanism is different):
# PlanetProfile builds the candidate list as the silicate-radius sweep,
# filters it to ``|C/MR^2 - Bulk.Cmeasured| <= half-width`` and then selects
# ``argmin |C/MR^2 - Cmeasured|`` over the survivors
# (LayerPropagators.FindInnerWithMoIAndEOS:1960-2032 for the Fe_CORE=False /
# POROUS_ROCK path Enceladus takes; FindInnerWithMoIAndConstantRho:1600-1659
# is the same predicate). The filtered sets are NESTED in the half-width, so:
#
#   * a node survives half-width H  <=>  min_i |CMR2_i - Cmeasured| <= H;
#   * whenever it survives, the SELECTED structure is IDENTICAL for every H
#     that admits it (the global argmin, once inside the window, stays the
#     argmin of every wider window).
#
# So on the ocean branch the half-width is a pure node-ACCEPTANCE threshold.
# It cannot bias WHICH structure a surviving node gets -- that is fixed by
# the argmin, which is the conditioning already DECLARED in config
# metadata.honesty_notes.cmr2_not_an_observable /
# ocean_branch_moi_selection_declared. Narrowing the window can only DELETE
# nodes (silently, to ``None``); widening it above the floor is inert. The
# non-conditioning direction is therefore WIDE, not narrow -- the opposite of
# the frozen branch, where the ruling had to pick the smallest workable
# window because there the selection is argmin|thickness| and a too-wide
# window admits liquid.
#
# The numerical half of the derivation (the measured floor) is recorded in
# the config block itself; see metadata.bulk_overrides_provenance.
_MOI_HALF_WIDTH_KEYS = ("Cuncertainty", "CuncertaintyUpper", "CuncertaintyLower")


def bulk_overrides_from_config(config: Any) -> Dict[str, Any]:
    """Resolve ``Planet.Bulk`` overrides declared by an inference config.

    Accepts an ``InferenceConfig``-like object (anything with a ``metadata``
    mapping), an already-parsed config ``dict``, or a path to the config
    JSON. Returns the ``metadata['bulk_overrides']`` block, with one
    expansion applied.

    THE EXPANSION. PlanetProfile's MoI acceptance test reads
    ``Bulk.CuncertaintyUpper`` / ``Bulk.CuncertaintyLower``, NOT
    ``Bulk.Cuncertainty``. Those two are only backfilled from
    ``Cuncertainty`` when the template left them ``None``
    (SetupInit.py:84-87), so a config that declared ``Cuncertainty`` alone
    would be a SILENT NO-OP against any template that sets them --
    the same defect class as the B7 dotted-key ``ocean_overrides`` no-op
    (config metadata.blockers_closed.B7). When the block declares
    ``Cuncertainty`` and does not declare an explicit upper/lower, both are
    set to the same value here so the declared window is enforced against
    every template.

    Raises ``ValueError`` when the config declares no MoI half-width at all.
    Falling through to the template default is exactly r5 blocker D1, so it
    is refused rather than defaulted.
    """
    if isinstance(config, str):
        import json as _json
        with open(config, "r") as fh:
            config = _json.load(fh)
    if isinstance(config, dict):
        meta = config.get("metadata") or {}
    else:
        meta = getattr(config, "metadata", None) or {}
    if not isinstance(meta, dict):
        raise ValueError(
            "bulk_overrides_from_config: config carries no metadata mapping; "
            "cannot resolve the MoI acceptance window (r5 blocker D1)."
        )
    raw = meta.get("bulk_overrides")
    if not isinstance(raw, dict) or not raw:
        raise ValueError(
            "bulk_overrides_from_config: config declares no "
            "metadata['bulk_overrides'] block. The ocean-branch MoI "
            "acceptance window MUST be declared in the config and read here "
            "-- silently inheriting the body template's own "
            "Bulk.Cuncertainty is r5 build blocker D1 (see "
            "validation_reports/enceladus_isostasy/r5_ADJUDICATION.md)."
        )
    out = {str(k): v for k, v in raw.items()}
    if not any(k in out for k in _MOI_HALF_WIDTH_KEYS):
        raise ValueError(
            "bulk_overrides_from_config: metadata['bulk_overrides'] declares "
            f"no MoI half-width (none of {_MOI_HALF_WIDTH_KEYS} present); "
            "got keys "
            f"{sorted(out)}. See r5 blocker D1."
        )
    if "Cuncertainty" in out:
        half = float(out["Cuncertainty"])
        out.setdefault("CuncertaintyUpper", half)
        out.setdefault("CuncertaintyLower", half)
    return out


def _template_cmeasured(planet_template_module: str) -> float | None:
    """``Bulk.Cmeasured`` of a body template, or ``None`` if unreadable.

    Diagnostic-only (the ``ocean_moi_window`` record); never gates a build,
    so an import failure must not take a cache down with it.
    """
    try:
        if planet_template_module in sys.modules:
            tmpl = sys.modules[planet_template_module].Planet
        else:
            tmpl = importlib.import_module(planet_template_module).Planet
        val = float(getattr(tmpl.Bulk, "Cmeasured", np.nan))
        return val if np.isfinite(val) else None
    except Exception:  # pragma: no cover - diagnostic only
        return None


def _moi_window_summary(
    structures: List[Dict[str, Any] | None],
    bulk_overrides: Dict[str, Any],
    Cmeasured: float | None,
) -> Dict[str, Any]:
    """Auditable record of the ocean-branch MoI window actually applied.

    Records the half-width the build ran with alongside the LARGEST achieved
    ``|C/MR^2 - Cmeasured|`` over the surviving nodes -- the measured floor
    below which a narrower window would have started deleting nodes. This is
    what makes the declared window checkable at production scale instead of
    resting on a probe (r5 blocker D1).
    """
    upper = bulk_overrides.get("CuncertaintyUpper",
                               bulk_overrides.get("Cuncertainty"))
    lower = bulk_overrides.get("CuncertaintyLower",
                               bulk_overrides.get("Cuncertainty"))
    devs = []
    if Cmeasured is not None and np.isfinite(Cmeasured):
        for s in structures:
            if s is None:
                continue
            c = s.get("CMR2")
            if c is None:
                continue
            c = float(c)
            if np.isfinite(c):
                devs.append(abs(c - float(Cmeasured)))
    max_dev = float(max(devs)) if devs else float("nan")
    half = float(upper) if upper is not None else float("nan")
    return {
        "Cmeasured": (float(Cmeasured) if Cmeasured is not None
                      and np.isfinite(Cmeasured) else None),
        "CuncertaintyUpper": (float(upper) if upper is not None else None),
        "CuncertaintyLower": (float(lower) if lower is not None else None),
        "source": bulk_overrides.get("_source", "explicit argument"),
        "n_nodes_measured": len(devs),
        "max_abs_cmr2_deviation": max_dev,
        "margin": (half - max_dev) if np.isfinite(half) and np.isfinite(max_dev)
                  else float("nan"),
        "reading": (
            "Ocean branch only. The half-width is a node-ACCEPTANCE "
            "threshold, not a selection knob: PP selects argmin|CMR2 - "
            "Cmeasured| over the in-window candidates, and those sets are "
            "nested in the half-width, so a surviving node's structure is "
            "identical for every half-width that admits it. "
            "max_abs_cmr2_deviation is the measured floor -- a half-width "
            "below it would have deleted at least one node to None."
        ),
    }


_FROZEN_G = 6.674e-11

# I-F2 audit (A8 item 2): the ratified window's upper edge is a CLOSED 0
# (mass_gate_blast_radius_audit.json's corrected signature: "EVERY residual
# ... NEGATIVE OR ZERO, max positive exactly 0.0000%"). This epsilon is pure
# floating-point slack on that boundary (outer-edge quadrature summing many
# shells), not a scientific tolerance -- it is many orders below the
# step-scaled lower edge it sits next to.
_FROZEN_IF2_POSITIVE_EPS = 1e-9


def _frozen_stored_mass_audit(
    struct: Dict[str, Any], M_target_kg: float
) -> Tuple[float, float]:
    """Re-derive achieved mass from a frozen node's STORED ``(r_m, rho)``
    arrays and return ``(residual_frac, step_scale_frac)`` (A8 item 2, the
    full ruling's I-F2 -- see :func:`build_zbw_grid_cache`'s docstring for
    the naming note).

    This is an AUDIT of the arrays that actually end up in the cache file,
    independent of ``mass_residual_frac`` (computed once, inside
    ``build_single_structure``, from the pre-storage arrays). Uses the
    SAME quadrature PlanetProfile's own layer convention and that field
    use: OUTER-EDGE -- ``rho[i]`` fills the shell ``[r[i-1], r[i]]`` with an
    implicit inner edge of 0 (``reducedPlanetModel`` / ``Geophysical.py:170``;
    the ``Mtot_achieved_kg`` note a few hundred lines above in this module).
    Re-deriving with the WRONG (midpoint) convention was exactly the error
    ``mass_gate_blast_radius_audit.json`` had to retract (a fabricated
    ~0.85% "systematic" that was pure quadrature-convention noise) --
    "an audit of existing data must adopt the convention the data was
    built with."

    ``step_scale_frac`` is the mass of the node's own OUTERMOST stored
    radial shell, as a fraction of ``M_target_kg``: the natural scale of
    the solver's documented -1-radial-step bias (``Silicates.py:118``,
    "first profile with ``Mtot_kg <= M_kg``"), and the scale
    :func:`build_zbw_grid_cache`'s ratified one-sided step-scaled window
    is defined against -- the acceptable range for ``residual_frac`` is
    ``(-step_scale_frac, 0]``, per the corrected
    ``mass_gate_blast_radius_audit.json`` signature ("EVERY residual ...
    NEGATIVE OR ZERO ... the one-signed ``(-delta_m_step, 0]`` signature").

    Returns ``(nan, nan)`` if the stored arrays are empty or ``M_target_kg``
    is non-positive/non-finite.
    """
    r = np.asarray(struct.get("r_m"), dtype=float)
    rho = np.asarray(struct.get("rho"), dtype=float)
    if r.size == 0 or rho.size != r.size:
        return float("nan"), float("nan")
    if not (np.isfinite(M_target_kg) and M_target_kg > 0):
        return float("nan"), float("nan")
    order = np.argsort(r)
    r_sorted = r[order]
    rho_sorted = rho[order]
    r_inner = np.concatenate(([0.0], r_sorted[:-1]))
    shell_mass = rho_sorted * (4.0 / 3.0) * np.pi * (r_sorted ** 3 - r_inner ** 3)
    Mtot_audit = float(np.sum(shell_mass))
    residual_frac = Mtot_audit / M_target_kg - 1.0
    step_scale_frac = float(shell_mass[-1]) / M_target_kg
    return residual_frac, step_scale_frac


def _volume_weighted(rho: np.ndarray, r_sorted: np.ndarray,
                     mask: np.ndarray) -> float:
    """Volume-weighted mean of ``rho`` over ``mask``, outer-edge convention."""
    r_in = np.concatenate(([0.0], r_sorted[:-1]))
    vol = (4.0 / 3.0) * np.pi * (r_sorted ** 3 - r_in ** 3)
    if not np.any(mask) or np.sum(vol[mask]) <= 0:
        return float("nan")
    return float(np.sum(rho[mask] * vol[mask]) / np.sum(vol[mask]))


def _frozen_node_fields(struct: Dict[str, Any]) -> Dict[str, Any]:
    """Derive the frozen branch's per-node diagnostics from a built structure.

    Returns ``rho_ice_mean_kgm3`` / ``rho_interior_mean_kgm3`` (volume-weighted
    over the ice-Ih and sub-hydrosphere layers respectively), ``n_liquid_layers``
    and ``zb_km_actual``.

    ``zb_km_actual`` is ``D_hsphere_km`` = ``(R_body - Sil.Rmean_m)/1e3``, the
    seafloor depth. On the frozen branch the whole hydrosphere is ice, so the
    seafloor IS the ice base and this is the shell thickness exactly. It is
    deliberately NOT ``D_iceIh_km`` (the summed Ih layer thicknesses), which
    undercounts by the inter-node radial gaps -- measured 0.86-1.18 km at
    Enceladus, i.e. 3-5x invariant I-F4's whole budget.
    """
    r = np.asarray(struct["r_m"], dtype=float)
    rho = np.asarray(struct["rho"], dtype=float)
    phases = np.asarray(struct["phases"])
    order = np.argsort(r)
    r, rho, phases = r[order], rho[order], phases[order]
    return {
        "n_liquid_layers": int(np.sum(phases == 0)),
        "rho_ice_mean_kgm3": _volume_weighted(rho, r, phases == 1),
        "rho_interior_mean_kgm3": _volume_weighted(rho, r, phases >= 50),
        "zb_km_actual": float(struct.get("D_hsphere_km", np.nan)),
    }


def _build_frozen_node(
    planet_template_module: str,
    zb_target_km: float,
    PbI_seed_MPa: float,
    ocean_overrides: Dict[str, Any],
    bulk_overrides: Dict[str, Any],
    planet_overrides: Dict[str, Any] | None,
    do_overrides: Dict[str, Any],
    extrap_ocean: bool,
    solve_tol_km: float,
    max_iter: int,
    progress: bool,
) -> Tuple[Dict[str, Any] | None, float, int]:
    """Root-find ``Bulk.PbISet_MPa`` so the frozen shell lands on ``zb_target_km``.

    ``zb(PbI)`` is smooth and monotone (hydrostatic to leading order), so a
    secant iteration converges in a handful of PlanetProfile runs. Seeded
    from the caller's hydrostatic estimate.

    Returns ``(struct | None, PbI_used_MPa, n_evals)``. ``struct`` carries the
    frozen diagnostics from :func:`_frozen_node_fields` merged in; ``None``
    means every attempt raised.
    """
    evals: List[Tuple[float, float, Dict[str, Any]]] = []  # (PbI, zb, struct)

    def evaluate(PbI: float):
        struct = build_single_structure(
            planet_template_module, 272.0,
            ocean_overrides=ocean_overrides,
            bulk_overrides={**bulk_overrides, "PbISet_MPa": float(PbI)},
            planet_overrides=planet_overrides,
            do_overrides=do_overrides,
            extrap_ocean=extrap_ocean,
        )
        fields = _frozen_node_fields(struct)
        struct.update(fields)
        zb = fields["zb_km_actual"]
        if not np.isfinite(zb):
            raise ValueError(f"non-finite achieved zb at PbISet_MPa={PbI}")
        evals.append((float(PbI), float(zb), struct))
        return float(zb)

    PbI = float(PbI_seed_MPa)
    try:
        zb = evaluate(PbI)
    except Exception as exc:
        log.warning(f"    frozen seed PbI={PbI:.4f} MPa failed — "
                    f"{type(exc).__name__}: {exc}")
        return None, PbI, len(evals)

    for _ in range(int(max_iter)):
        if abs(zb - zb_target_km) < solve_tol_km:
            break
        if len(evals) >= 2 and evals[-1][0] != evals[-2][0]:
            # Secant on the last two distinct evaluations.
            (p0, z0, _), (p1, z1, _) = evals[-2], evals[-1]
            if z1 == z0:
                break
            PbI_next = p1 + (zb_target_km - z1) * (p1 - p0) / (z1 - z0)
        else:
            # First correction: proportional, since zb ~ PbI to leading order.
            if zb <= 0:
                break
            PbI_next = PbI * (zb_target_km / zb)
        if not np.isfinite(PbI_next) or PbI_next <= 0:
            break
        PbI = float(PbI_next)
        try:
            zb = evaluate(PbI)
        except Exception as exc:
            log.warning(f"    frozen iterate PbI={PbI:.4f} MPa failed — "
                        f"{type(exc).__name__}: {exc}")
            break

    # Best achieved node, whether or not the loop converged; the I-F4
    # invariant in the caller decides whether it is good enough.
    best = min(evals, key=lambda e: abs(e[1] - zb_target_km))
    if progress:
        log.info(f"    frozen solve: {len(evals)} PP run(s), "
                 f"PbI={best[0]:.5f} MPa -> zb={best[1]:.4f} km "
                 f"(target {zb_target_km:.4f} km)")
    return best[2], best[0], len(evals)


def build_zbw_grid_cache(
    planet_template_module: str,
    zb_km_grid: Iterable[float],
    wOcean_ppt_grid: Iterable[float],
    output_path: str,
    progress: bool = True,
    ocean_overrides: Dict[str, Any] | None = None,
    bulk_overrides: Dict[str, Any] | None = None,
    planet_overrides: Dict[str, Any] | None = None,
    do_overrides: Dict[str, Any] | None = None,
    extrap_ocean: bool = False,
    tb_placeholder_K: float = 272.0,
    zb_tol_km: float | None = None,
    frozen_zb_km_grid: Iterable[float] | None = None,
    frozen_wOcean_ppt: float | None = None,
    frozen_Cuncertainty: float = 0.015,
    frozen_mass_tol: float = 1e-6,
    frozen_rho_closure_tol_kgm3: float = 12.0,
    frozen_zb_tol_km: float = 0.25,
    frozen_moi_nonconditioning_window: Tuple[float, float] | None = None,
    frozen_max_iter: int = 10,
    config: Any | None = None,
) -> Dict[str, Any]:
    """Build a (zb_km × wOcean_ppt) ocean grid plus an optional FROZEN zb axis.

    Schema ``v3.1-zbw`` (ocean only, unchanged) or ``v3.2-zbw-joint`` when
    ``frozen_zb_km_grid`` is supplied.

    THE TWO BRANCHES DO NOT SHARE AN AXIS (frozen-branch design ruling):
    the frozen arrays are SEPARATE (``frozen_zb_km_grid`` /
    ``frozen_structures``) and are NOT crossed with ``w``. Salinity is
    UNDEFINED without an ocean -- not merely unconstrained -- so the frozen
    segment is built once at a single nominal ``w`` (config
    ``structure_cache_spec.w_axis_scope``); multiplying 39 frozen zb nodes by
    40 salinities would be 39 copies of the same structure carrying a
    meaningless coordinate. A single ragged shared zb axis is likewise
    rejected: the ocean branch's support ends near 45 km and the frozen
    branch's begins near 46.7 km, and splicing them made the reported branch
    odds linearly rescalable by an arbitrary box edge (config
    ``branch_model.why_reparameterized``).

    OCEAN BRANCH (``zb_km_grid`` × ``wOcean_ppt_grid``, schema v3.1 rows):
    unchanged. A node whose achieved structure has no phase-0 (liquid) layer
    is rejected to ``None``, same as an out-of-tolerance zb miss -- a frozen
    result on the OCEAN axis is out of scope there, not a failure to retry.

    FROZEN BRANCH (``frozen_zb_km_grid``, one node per zb): built by the
    constant-density mass-closure route described in the module-level
    ``_build_frozen_node`` comment above. On this branch zb is NOT a free
    parameter -- at fixed ``rho_ice`` mass conservation makes zb <-> rho_rock
    a bijection (``isostasy.frozen_zb_from_mass``), so the grid indexes the
    sampled rock density through that map and the recorded support
    (zb in [46.74, 65.56] km for rho_rock ~ U[2200,2600] x rho_ice ~
    U[915,935]) is fixed by the density prior rather than by a box edge.

    BUILD-TIME INVARIANTS on every frozen node (ruling I-F1..I-F4).
    All are HARD REJECT TO ``None`` -- never a warning, never a nearest-node
    fallback:

    - **I-F1 mass closure**, ``abs(mass_residual_frac) <= frozen_mass_tol``
      (default 1e-6). The constant-rho path solves the interior density
      analytically, so this is exact to floating point (measured: 1e-16 at
      six probed nodes), NOT the one-radial-step residual the EOS path
      leaves. A node that cannot meet 1e-6 did not take the constant-rho
      path and must not be silently kept -- that is exactly the MAJOR-1
      defect, where the smoke cache stored the template ``Mtot_kg`` while
      being -22.16% / -9.21% off.
    - **I-F2 stored-array mass-quadrature audit** (A8 item 2). Re-derives
      the achieved mass from the node's own STORED ``(r_m, rho)`` arrays
      (outer-edge shell quadrature, PlanetProfile's own convention -- see
      :func:`_frozen_stored_mass_audit`) and gates it against the ratified
      ONE-SIDED, STEP-SCALED window ``(-step_scale_frac, 0]``, where
      ``step_scale_frac`` is the mass of the node's own outermost stored
      shell as a fraction of the template mass -- the scale of the
      solver's documented -1-radial-step bias
      (``frozen_node_mass_violation_rootcause.json``,
      ``mass_gate_blast_radius_audit.json``'s corrected convention). NO
      global constant is invented here; the window is per-node and derived
      from the node's own grid. Near-redundant with I-F1 under this
      builder's analytic constant-rho closure (both should read ~machine
      precision on a healthy node) but ORDERED separately because it
      audits what is actually WRITTEN to the cache, not the field computed
      inside ``build_single_structure`` before storage. Rejects are
      counted in ``n_frozen_if2_rejected`` (also on the returned cache
      dict's top level). NOTE ON LABELS: the pre-existing ``I-F2_liquid``
      rejects-dict key below (from A2-A7) is, per the manager ruling's own
      numbering, actually the full ruling's I-F3 ("branch flag from the
      phase array") -- its content is unchanged and confirmed correct, but
      a renumbering of that key was not ordered as part of this task, so
      the two "I-F2"-prefixed labels below name genuinely different checks.
    - **no liquid** (code label ``I-F2_liquid``, see the note just above),
      zero phase-0 layers, read from the PHASE ARRAY and never from
      ``D_ocean_km`` (frozen smoke nodes stored ``D_ocean_km`` 1.510 /
      1.151 km with no phase-0 layers at all). This is the invariant that
      catches a thin liquid film -- the failure mode of driving the
      seafloor to a target depth BELOW the ice base.
    - **I-F3 parameterization closure** (code label ``I-F3_closure``), the
      achieved interior density must agree with
      ``isostasy.frozen_rho_rock_from_zb`` evaluated at the node's own
      achieved zb and achieved mean ice density, to within
      ``frozen_rho_closure_tol_kgm3`` (default 12 kg/m^3). This is what
      makes the sampled coordinate ``rho_rock`` mean the same thing in the
      cache as it does in the prior.
    - **I-F4 zb placement**, ``abs(zb_km_actual - zb_km_node) <
      frozen_zb_tol_km`` (default 0.25 km).

    BUILD-LEVEL INVARIANT (ruling I-F6, non-conditioning). If EVERY
    surviving frozen node's derived ``CMR2`` lands inside the template's own
    MoI window (``Bulk.Cmeasured +/- Bulk.Cuncertainty``, Enceladus
    0.335 +/- 0.001), the BUILD FAILS with ``RuntimeError`` -- it is not a
    node-level reject. A frozen set that cannot leave a 0.002-wide window is
    a set that was SELECTED by the MoI, which would reintroduce ruling F1's
    undeclared conditioning under a new name. Measured on the intended
    Enceladus grid the derived C/MR^2 spans 0.322-0.344, i.e. 10x outside the
    window in both directions, so the invariant passes with wide margin.
    Override the window via ``frozen_moi_nonconditioning_window``.

    Mechanism (PlanetProfile's existing "two-of-three" zb/Tb pairing,
    ``Planet.Do.ICEIh_THICKNESS`` + ``Planet.Bulk.zb_approximate_km`` --
    Main.py's ``oceanComp`` inductogram path, ``LayerPropagators.
    GetIceShellTFreeze``): each node sets ``Bulk.zb_approximate_km`` to the
    grid's target shell thickness and ``Do.ICEIh_THICKNESS=True``. That
    clears the OTHER member of the pair -- SetupInit.py sets
    ``Bulk.zb_approximate_km = nan`` when ``ICEIh_THICKNESS`` is False and,
    symmetrically, ``GetIceShellTFreeze`` root-finds ``Tb_K`` and
    REASSIGNS ``Planet.Bulk.Tb_K`` to the solved value -- so the ``Tb_K``
    passed into :func:`build_single_structure` here is a placeholder only
    (never read on this path; verified empirically). The SOLVED Tb_K is
    read back via the ``"Tb_K"`` field of the returned structure dict
    (fixed to read ``Planet.Bulk.Tb_K`` post-run rather than echo the
    input -- see :func:`build_single_structure`).

    Root-finding is only tolerance-bounded in TEMPERATURE space
    (``GetIceShellTFreeze``'s internal ``xtol=Planet.TfreezeRes_K``), not
    in zb space. At Enceladus gravity d(zb)/d(Tb) is enormous (~0.27 K
    spans the WHOLE 5-45 km shell range per the config's own honesty
    note), so a node can converge in T while still missing its zb target
    by several km. This builder enforces the zb-placement invariant
    explicitly (module spec / RESUME_NOTE ``zb_placement``): a node is
    HARD REJECTED to ``None`` when
    ``abs(D_iceIh_km_achieved - zb_km_target) >= zb_tol_km``.
    ``D_iceIh_km`` (summed phase-Ih layer thickness, already computed by
    :func:`build_single_structure`) is used as the "achieved zb" reading --
    the existing, already-tested field for shell thickness in this
    package; it is NOT bit-identical to PlanetProfile's own internal
    ``Planet.zb_km`` (observed offset ~0.4 km on one probe structure,
    plausibly the Melosh conductive-layer correction), which is not
    exposed by :func:`build_single_structure`'s return contract.

    NOT implemented here (explicitly out of scope -- reviewer-blocked
    production physics, RESUME_NOTE.md "Blocked on reviewer sign-off"):
    the per-node MASS-CONSERVATION invariant (MAJOR-1) and its two-part
    ``iSilStart`` retry fix. This builder does not compute or gate on a
    mass residual; only the zb-placement check above and ordinary build
    failures reject nodes to ``None``.

    Parameters
    ----------
    zb_km_grid
        Iterable of target ice-shell thicknesses (km). Cast to 1-D float,
        sorted ascending.
    wOcean_ppt_grid
        See :func:`build_tbw_grid_cache`.
    tb_placeholder_K
        Value passed as the (unused, overwritten-by-solve) ``Tb_K``
        argument to :func:`build_single_structure`. Default 272.0 K is
        empirically verified to build successfully across the Enceladus
        template's zb range; override for a different body/template.
    zb_tol_km
        Half-width of the zb-placement acceptance window (km). Default
        (``None``) uses half the smallest spacing in ``zb_km_grid``.
    frozen_zb_km_grid
        Optional 1-D iterable of frozen-branch shell thicknesses (km). When
        ``None`` (default) no frozen branch is built and the output is
        byte-compatible schema ``v3.1-zbw``. The intended Enceladus grid is
        ``arange(46.5, 65.8 + eps, 0.5)`` (39 nodes), which brackets the
        true frozen support [46.74, 65.56] km on both sides.
    frozen_wOcean_ppt
        Nominal salinity for the frozen build (g/kg). Default ``None`` uses
        the template's own ``Ocean.wOcean_ppt``. It reaches the frozen
        structure only through the melt curve that converts the set ice-base
        pressure into ``Tb_K``, and is recorded per node as
        ``wOcean_ppt_nominal`` with ``w_is_defined = False``.
    frozen_Cuncertainty
        MoI acceptance half-width applied to the frozen build via
        ``Bulk.Cuncertainty`` (default 0.015, ruling 3.2's documented
        fallback). The ruling's PREFERRED route is no MoI gate at all; see
        this function's implementation notes and the A3 report for why
        PlanetProfile's no-gate path
        (``Do.SPECIFY_HYDROSPHERE_SEAFLOOR_PRESSURE``) cannot express a
        seafloor coincident with the ice base at Enceladus. 0.015 is chosen
        over the template's 0.001 because the frozen set's derived C/MR^2
        spans 0.322-0.344 (max deviation 0.0108), so 0.001 would reject
        every node while 0.015 admits all of them with >=0.004 margin; and
        because with ``Dhsphere_m = 0`` the window cannot bias WHICH
        admissible structure is selected, only whether the intended one
        survives (see the module-level note above).
    frozen_mass_tol, frozen_rho_closure_tol_kgm3, frozen_zb_tol_km
        Invariant I-F1 / I-F3 / I-F4 thresholds.
    frozen_moi_nonconditioning_window
        ``(C_centre, C_halfwidth)`` for invariant I-F6. Default ``None``
        reads the template's ``Bulk.Cmeasured`` / ``Bulk.Cuncertainty``.
    frozen_max_iter
        Maximum secant iterations in the per-node ``PbISet_MPa`` solve.
    config
        Optional inference config (``InferenceConfig``-like object, parsed
        dict, or path to the config JSON) supplying the OCEAN-branch
        ``Planet.Bulk`` overrides -- above all the MoI acceptance half-width
        -- via ``metadata['bulk_overrides']``
        (:func:`bulk_overrides_from_config`). Without this the builder
        inherits the body template's own ``Bulk.Cuncertainty``
        (0.001 at Enceladus), which is r5 build blocker D1: the campaign
        config's asserted window was unreadable and therefore
        unenforceable. Explicit ``bulk_overrides`` entries WIN over the
        config's, so a caller can still probe a different window; the
        resolved window is recorded in the returned cache under
        ``ocean_moi_window``. Raises ``ValueError`` when a config is passed
        that declares no MoI half-width -- defaulting silently is the
        defect.

    Returns
    -------
    The saved-to-disk dict::

        {'zb_km_grid': np.ndarray (n_zb,),
         'wOcean_ppt_grid': np.ndarray (n_w,),
         'structures': [structure_dict | None, ...],  # row-major, len n_zb*n_w
         'ocean_comp': str | None,
         'schema_version': 'v3.1-zbw' | 'v3.2-zbw-joint',
         'zb_tol_km': float,
         'n_zb_placement_rejected': int,
         'frozen_branch_supported': bool,
         # v3.2-zbw-joint only:
         'frozen_zb_km_grid': np.ndarray (n_fz,),
         'frozen_structures': [structure_dict | None, ...],  # len n_fz
         'frozen_build_spec': {...},
         # A8 item 2: count of nodes rejected by the stored-array I-F2
         # mass-quadrature audit (mirrors frozen_build_spec's own field).
         'n_frozen_if2_rejected': int}

    ``structures`` is row-major: entry ``i_zb * n_w + i_w`` is the
    structure at ``(zb_km_grid[i_zb], wOcean_ppt_grid[i_w])``. Each
    surviving structure additionally carries ``zb_km_node`` (the grid
    target) alongside its existing ``Tb_K`` (now the SOLVED value) and
    ``D_iceIh_km`` (the achieved zb) fields.

    ``frozen_structures`` is one node per ``frozen_zb_km_grid`` entry. Each
    surviving frozen structure carries, on top of the standard
    :func:`build_single_structure` fields:

    ``branch`` ('frozen'), ``has_ocean`` (False), ``zb_km_node``,
    ``zb_km_actual``, ``zb_residual_km``, ``n_liquid_layers`` (0),
    ``rho_ice_mean_kgm3``, ``rho_interior_mean_kgm3``,
    ``rho_rock_closure_kgm3`` (the analytic
    ``frozen_rho_rock_from_zb`` value at this node's own achieved zb and ice
    density), ``rho_rock_closure_residual_kgm3``, ``mass_residual_frac``,
    ``Mtot_achieved_kg``, ``cmr2_derived`` (DIAGNOSTIC -- never a gate),
    ``PbISet_MPa``, ``wOcean_ppt_nominal`` and ``w_is_defined`` (False).
    """
    zb_arr = np.asarray(list(zb_km_grid), dtype=np.float64)
    w_arr = np.asarray(list(wOcean_ppt_grid), dtype=np.float64)
    if zb_arr.ndim != 1 or zb_arr.size < 2:
        raise ValueError(
            f"zb_km_grid must be a 1-D iterable with >= 2 points, got {zb_arr!r}"
        )
    if w_arr.ndim != 1 or w_arr.size < 2:
        raise ValueError(
            f"wOcean_ppt_grid must be a 1-D iterable with >= 2 points, "
            f"got {w_arr!r}"
        )
    zb_arr = zb_arr[np.argsort(zb_arr)]
    w_arr = w_arr[np.argsort(w_arr)]

    if zb_tol_km is None:
        zb_tol_km = float(np.min(np.diff(zb_arr)) / 2.0)
    zb_tol_km = float(zb_tol_km)

    base_overrides = dict(ocean_overrides or {})
    base_overrides.pop("wOcean_ppt", None)
    ocean_comp = base_overrides.get("comp")
    if progress and ocean_comp is not None:
        log.info(f"Ocean composition baked into every node: comp={ocean_comp!r}")

    base_do_overrides = dict(do_overrides or {})
    base_do_overrides["ICEIh_THICKNESS"] = True
    # D1: the config's declared Bulk overrides (the MoI acceptance window)
    # underlie any explicitly-passed ones, which win.
    _window_source = "explicit argument"
    if config is not None:
        cfg_bulk = bulk_overrides_from_config(config)
        explicit = dict(bulk_overrides or {})
        base_bulk_overrides = {**cfg_bulk, **explicit}
        _window_source = (
            "config metadata['bulk_overrides']"
            if not any(k in explicit for k in _MOI_HALF_WIDTH_KEYS)
            else "explicit argument (overrides config "
                 "metadata['bulk_overrides'])")
        if progress:
            log.info(
                f"Bulk overrides from config: {cfg_bulk}; "
                f"resolved: {base_bulk_overrides}")
    else:
        base_bulk_overrides = dict(bulk_overrides or {})

    n_zb, n_w = zb_arr.size, w_arr.size
    total = n_zb * n_w
    structures: List[Dict[str, Any] | None] = [None] * total
    n_ok = 0
    n_fail = 0
    n_zb_reject = 0
    for i_zb, zb_km in enumerate(zb_arr):
        for i_w, w_ppt in enumerate(w_arr):
            flat = i_zb * n_w + i_w
            if progress:
                log.info(
                    f"[{flat + 1}/{total}] zb_km={zb_km:.4f} "
                    f"wOcean_ppt={w_ppt:.4f}"
                )
            try:
                struct = build_single_structure(
                    planet_template_module, float(tb_placeholder_K),
                    ocean_overrides={**base_overrides,
                                     "wOcean_ppt": float(w_ppt)},
                    bulk_overrides={**base_bulk_overrides,
                                    "zb_approximate_km": float(zb_km)},
                    planet_overrides=planet_overrides,
                    do_overrides=base_do_overrides,
                    extrap_ocean=extrap_ocean,
                )
            except Exception as exc:
                log.warning(
                    f"    × (zb={zb_km:.4f}, w={w_ppt:.4f}) → None — "
                    f"{type(exc).__name__}: {exc}"
                )
                structures[flat] = None
                n_fail += 1
                continue

            has_ocean = bool(np.any(np.asarray(struct["phases"]) == 0))
            if not has_ocean:
                # Ocean-branch-only builder (see docstring): a frozen
                # result at a requested zb is out of scope, not a failure
                # to retry.
                log.warning(
                    f"    × (zb={zb_km:.4f}, w={w_ppt:.4f}) → None — "
                    f"no phase-0 layer (frozen; ocean-branch-only builder)"
                )
                structures[flat] = None
                n_fail += 1
                continue

            zb_actual_km = float(struct.get("D_iceIh_km", np.nan))
            zb_residual_km = zb_actual_km - float(zb_km)
            if not (np.isfinite(zb_actual_km)
                    and abs(zb_residual_km) < zb_tol_km):
                log.warning(
                    f"    × (zb={zb_km:.4f}, w={w_ppt:.4f}) → None — "
                    f"zb_placement invariant failed: achieved "
                    f"{zb_actual_km:.4f} km, residual {zb_residual_km:+.4f} "
                    f"km >= tol {zb_tol_km:.4f} km"
                )
                structures[flat] = None
                n_zb_reject += 1
                n_fail += 1
                continue

            struct["has_ocean"] = True
            struct["zb_km_node"] = float(zb_km)
            struct["zb_km_actual"] = zb_actual_km
            struct["zb_residual_km"] = zb_residual_km
            structures[flat] = struct
            n_ok += 1
            if progress:
                log.info(
                    f"    → Tb_K_solved={struct['Tb_K']:.4f}, "
                    f"zb_actual={zb_actual_km:.4f} km "
                    f"(residual {zb_residual_km:+.4f} km), "
                    f"CMR²={struct.get('CMR2', float('nan')):.4f}"
                )

    if n_ok == 0:
        raise RuntimeError(
            f"All {total} (zb, w) nodes failed; no cache to write."
        )

    cache = {
        "zb_km_grid": zb_arr,
        "wOcean_ppt_grid": w_arr,
        "structures": structures,
        "ocean_comp": ocean_comp,
        "schema_version": "v3.1-zbw",
        "zb_tol_km": zb_tol_km,
        "n_zb_placement_rejected": int(n_zb_reject),
        # Ocean-branch-only unless a frozen axis was requested below.
        "frozen_branch_supported": False,
        # D1: the MoI acceptance window this build actually ran with, plus
        # the measured floor (largest achieved |CMR2 - Cmeasured|) that makes
        # the declared window auditable at production scale.
        "ocean_moi_window": _moi_window_summary(
            structures,
            {**base_bulk_overrides, "_source": _window_source},
            _template_cmeasured(planet_template_module)),
    }

    if frozen_zb_km_grid is not None:
        frozen_arr, frozen_structs, frozen_spec = _build_frozen_axis(
            planet_template_module=planet_template_module,
            frozen_zb_km_grid=frozen_zb_km_grid,
            base_ocean_overrides=base_overrides,
            base_bulk_overrides=base_bulk_overrides,
            planet_overrides=planet_overrides,
            base_do_overrides=do_overrides,
            extrap_ocean=extrap_ocean,
            frozen_wOcean_ppt=frozen_wOcean_ppt,
            frozen_Cuncertainty=frozen_Cuncertainty,
            frozen_mass_tol=frozen_mass_tol,
            frozen_rho_closure_tol_kgm3=frozen_rho_closure_tol_kgm3,
            frozen_zb_tol_km=frozen_zb_tol_km,
            frozen_moi_nonconditioning_window=(
                frozen_moi_nonconditioning_window),
            frozen_max_iter=frozen_max_iter,
            progress=progress,
        )
        cache["schema_version"] = "v3.2-zbw-joint"
        cache["frozen_branch_supported"] = True
        cache["frozen_zb_km_grid"] = frozen_arr
        cache["frozen_structures"] = frozen_structs
        cache["frozen_build_spec"] = frozen_spec
        # A8 item 2: top-level mirror of frozen_spec['n_frozen_if2_rejected']
        # (stored-array I-F2 mass-quadrature audit reject count).
        cache["n_frozen_if2_rejected"] = int(
            frozen_spec["n_frozen_if2_rejected"])

    _save_cache_atomic(cache, output_path)
    if progress:
        msg = (f"zb×w cache written → {output_path} "
               f"({n_zb} zb × {n_w} w = {total} ocean nodes: "
               f"{n_ok} built, {n_fail} None, "
               f"{n_zb_reject} rejected by zb_placement)")
        if frozen_zb_km_grid is not None:
            spec = cache["frozen_build_spec"]
            msg += (f"; frozen axis {spec['n_nodes']} nodes: "
                    f"{spec['n_built']} built, {spec['n_rejected']} None "
                    f"(schema {cache['schema_version']})")
        log.info(msg)
    return cache


def _build_frozen_axis(
    planet_template_module: str,
    frozen_zb_km_grid: Iterable[float],
    base_ocean_overrides: Dict[str, Any],
    base_bulk_overrides: Dict[str, Any],
    planet_overrides: Dict[str, Any] | None,
    base_do_overrides: Dict[str, Any] | None,
    extrap_ocean: bool,
    frozen_wOcean_ppt: float | None,
    frozen_Cuncertainty: float,
    frozen_mass_tol: float,
    frozen_rho_closure_tol_kgm3: float,
    frozen_zb_tol_km: float,
    frozen_moi_nonconditioning_window: Tuple[float, float] | None,
    frozen_max_iter: int,
    progress: bool,
) -> Tuple[np.ndarray, List[Dict[str, Any] | None], Dict[str, Any]]:
    """Build the frozen (no-ocean) zb axis of a v3.2-zbw-joint cache.

    See :func:`build_zbw_grid_cache` for the mechanism, the per-node field
    contract and invariants I-F1..I-F4 / I-F6. Split out only to keep the
    two branches' loops legible; it is not part of the public API.
    """
    from PlanetProfile.Gravity.isostasy import frozen_rho_rock_from_zb

    fz_arr = np.asarray(list(frozen_zb_km_grid), dtype=np.float64)
    if fz_arr.ndim != 1 or fz_arr.size < 1:
        raise ValueError(
            f"frozen_zb_km_grid must be a 1-D iterable with >= 1 point, "
            f"got {fz_arr!r}")
    fz_arr = fz_arr[np.argsort(fz_arr)]

    # Template scalars: body mass/radius for the analytic closure, and the
    # MoI window I-F6 tests non-conditioning against.
    if planet_template_module in sys.modules:
        importlib.reload(sys.modules[planet_template_module])
    else:
        importlib.import_module(planet_template_module)
    tmpl = sys.modules[planet_template_module].Planet
    M_body_kg = float(tmpl.Bulk.M_kg)
    R_body_m = float(tmpl.Bulk.R_m)
    if frozen_moi_nonconditioning_window is None:
        frozen_moi_nonconditioning_window = (
            float(getattr(tmpl.Bulk, "Cmeasured", np.nan)),
            float(getattr(tmpl.Bulk, "Cuncertainty", np.nan)))
    w_nominal = (float(frozen_wOcean_ppt) if frozen_wOcean_ppt is not None
                 else float(getattr(tmpl.Ocean, "wOcean_ppt", np.nan)))

    frozen_ocean_overrides = {**base_ocean_overrides,
                              "wOcean_ppt": w_nominal}
    frozen_bulk_overrides = {
        **base_bulk_overrides,
        # Seafloor = shallowest admissible candidate = the ice base, so the
        # hydrosphere is entirely ice (I-F2). See the module-level note.
        "Dhsphere_m": 0.0,
        "Cuncertainty": float(frozen_Cuncertainty),
        "CuncertaintyUpper": float(frozen_Cuncertainty),
        "CuncertaintyLower": float(frozen_Cuncertainty),
    }
    frozen_do_overrides = {
        **(base_do_overrides or {}),
        "ConstantProps": {"Inner": True, "Ocean": False, "Ice": False},
        "Fe_CORE": False,
        "POROUS_ROCK": False,
        "SPECIFY_ICEI_BOTTOM_PRESSURE": True,
        "HYDROSPHERE_THICKNESS": True,
        # The frozen branch drives zb through PbISet_MPa, NOT through the
        # ocean branch's ICEIh_THICKNESS root-find (measured 2.1 km miss).
        "ICEIh_THICKNESS": False,
    }

    # Hydrostatic seed for the PbI solve: P ~ rho_ice g zb, with g = GM/R^2.
    g_surf = _FROZEN_G * M_body_kg / R_body_m ** 2
    rho_ice_seed = 930.0

    structs: List[Dict[str, Any] | None] = [None] * fz_arr.size
    n_built = 0
    # NOTE (A8): "I-F2_liquid" here is the code label already in place from
    # A2-A7 and is left UNCHANGED by this task -- per the manager ruling
    # (frozen_branch_DESIGN_RULING.md, "Implementation record + manager
    # rulings", item 3) its content is confirmed correct but it is actually
    # the full ruling's I-F3 ("branch flag from the phase array"); a
    # renumbering was not ordered here, so both labels stand side by side.
    # "I-F2_mass_quad" below is the full ruling's ACTUAL I-F2 (stored-array
    # outer-edge mass-quadrature audit), added by this task.
    rejects = {"solve": 0, "I-F1_mass": 0, "I-F2_mass_quad": 0,
               "I-F2_liquid": 0, "I-F3_closure": 0, "I-F4_placement": 0}
    n_pp_runs = 0
    for i, zb_node in enumerate(fz_arr):
        if progress:
            log.info(f"[frozen {i + 1}/{fz_arr.size}] zb_km={zb_node:.4f}")
        PbI_seed = rho_ice_seed * g_surf * float(zb_node) * 1e3 / 1e6
        struct, PbI_used, n_evals = _build_frozen_node(
            planet_template_module, float(zb_node), PbI_seed,
            frozen_ocean_overrides, frozen_bulk_overrides, planet_overrides,
            frozen_do_overrides, extrap_ocean,
            solve_tol_km=0.4 * float(frozen_zb_tol_km),
            max_iter=int(frozen_max_iter), progress=progress)
        n_pp_runs += n_evals
        if struct is None:
            rejects["solve"] += 1
            continue

        # ---- invariants, in order; every failure is a HARD REJECT to None
        zb_actual = float(struct["zb_km_actual"])
        zb_residual = zb_actual - float(zb_node)
        mass_res = float(struct.get("mass_residual_frac", np.nan))
        rho_ice_mean = float(struct["rho_ice_mean_kgm3"])
        rho_int_mean = float(struct["rho_interior_mean_kgm3"])
        rho_closure = frozen_rho_rock_from_zb(
            M_body_kg, R_body_m, zb_actual * 1e3, rho_ice_mean)
        closure_residual = (rho_int_mean - rho_closure
                            if rho_closure is not None else np.nan)
        # I-F2 (A8 item 2): audit the STORED (r_m, rho) arrays, independent
        # of the mass_residual_frac field computed inside
        # build_single_structure. See _frozen_stored_mass_audit.
        if2_residual_frac, if2_step_scale_frac = _frozen_stored_mass_audit(
            struct, M_body_kg)

        failed = None
        if not (np.isfinite(mass_res) and abs(mass_res) <= frozen_mass_tol):
            failed = ("I-F1_mass",
                      f"mass_residual_frac {mass_res:+.6e} exceeds "
                      f"{frozen_mass_tol:.1e}")
        elif not (np.isfinite(if2_residual_frac)
                  and np.isfinite(if2_step_scale_frac)
                  and -if2_step_scale_frac < if2_residual_frac
                  <= _FROZEN_IF2_POSITIVE_EPS):
            failed = ("I-F2_mass_quad",
                      f"stored-array mass audit residual "
                      f"{if2_residual_frac:+.6e} outside ratified window "
                      f"(-{if2_step_scale_frac:.6e}, "
                      f"{_FROZEN_IF2_POSITIVE_EPS:.1e}]")
        elif int(struct["n_liquid_layers"]) != 0:
            failed = ("I-F2_liquid",
                      f"{int(struct['n_liquid_layers'])} phase-0 (liquid) "
                      "layers on a frozen node")
        elif not (np.isfinite(closure_residual)
                  and abs(closure_residual) <= frozen_rho_closure_tol_kgm3):
            failed = ("I-F3_closure",
                      f"achieved rho_interior {rho_int_mean:.3f} vs analytic "
                      f"{rho_closure if rho_closure is None else f'{rho_closure:.3f}'}"
                      f" kg/m^3 (residual {closure_residual:+.3f}) exceeds "
                      f"{frozen_rho_closure_tol_kgm3:.3f} kg/m^3")
        elif not abs(zb_residual) < frozen_zb_tol_km:
            failed = ("I-F4_placement",
                      f"achieved zb {zb_actual:.4f} km, residual "
                      f"{zb_residual:+.4f} km >= tol "
                      f"{frozen_zb_tol_km:.4f} km")
        if failed is not None:
            key, why = failed
            rejects[key] += 1
            log.warning(f"    × frozen zb={zb_node:.4f} → None — "
                        f"{key} invariant failed: {why}")
            continue

        struct.update({
            "branch": "frozen",
            "has_ocean": False,
            "zb_km_node": float(zb_node),
            "zb_residual_km": float(zb_residual),
            "rho_rock_closure_kgm3": float(rho_closure),
            "rho_rock_closure_residual_kgm3": float(closure_residual),
            # I-F2 stored-array mass-quadrature audit result (A8 item 2),
            # kept for downstream auditability -- never re-gates anything.
            "if2_mass_audit_residual_frac": float(if2_residual_frac),
            "if2_mass_audit_step_scale_frac": float(if2_step_scale_frac),
            # DIAGNOSTIC. Ruling 3.2: C/MR^2 is derived along the frozen
            # bijection and reported; it never selects a structure.
            "cmr2_derived": float(struct.get("CMR2", np.nan)),
            "PbISet_MPa": float(PbI_used),
            "wOcean_ppt_nominal": w_nominal,
            "w_is_defined": False,
        })
        structs[i] = struct
        n_built += 1
        if progress:
            log.info(
                f"    → zb={zb_actual:.4f} km ({zb_residual:+.4f}), "
                f"rho_ice={rho_ice_mean:.3f}, rho_rock={rho_int_mean:.3f} "
                f"(closure {closure_residual:+.4f}), "
                f"mass_res={mass_res:+.2e}, C/MR²={struct['cmr2_derived']:.5f}")

    if n_built == 0:
        raise RuntimeError(
            f"All {fz_arr.size} frozen zb nodes failed; rejects={rejects}. "
            "No frozen branch to write.")

    # ---- I-F6 (build-level, both-directions non-conditioning) ------------
    C_centre, C_half = frozen_moi_nonconditioning_window
    cmr2_built = np.array([s["cmr2_derived"] for s in structs
                           if s is not None], dtype=float)
    if np.isfinite(C_centre) and np.isfinite(C_half):
        inside = np.abs(cmr2_built - C_centre) <= C_half
        if bool(np.all(inside)):
            raise RuntimeError(
                "I-F6 non-conditioning invariant FAILED: all "
                f"{cmr2_built.size} surviving frozen nodes have derived "
                f"C/MR² inside {C_centre:.4f} ± {C_half:.4f} "
                f"(span {cmr2_built.min():.5f}-{cmr2_built.max():.5f}). "
                "A frozen set that cannot leave the MoI window was SELECTED "
                "by the MoI, which is the undeclared conditioning the "
                "frozen-branch design ruling (F1) exists to remove. This "
                "FAILS THE BUILD by design; it is not a node-level reject.")

    spec = {
        "mechanism": (
            "Do.ConstantProps['Inner'] + Fe_CORE=False + POROUS_ROCK=False "
            "(analytic constant-rho mass closure); "
            "Do.SPECIFY_ICEI_BOTTOM_PRESSURE + Bulk.PbISet_MPa (zb control, "
            "secant-solved per node); Do.HYDROSPHERE_THICKNESS + "
            "Bulk.Dhsphere_m=0 (seafloor = ice base => no liquid)"),
        "n_nodes": int(fz_arr.size),
        "n_built": int(n_built),
        "n_rejected": int(fz_arr.size - n_built),
        "rejects_by_invariant": dict(rejects),
        # A8 item 2: explicit named count of the stored-array I-F2 mass-
        # quadrature audit rejects, mirrored onto the top-level cache dict
        # by build_zbw_grid_cache as n_frozen_if2_rejected (same value as
        # rejects_by_invariant['I-F2_mass_quad'], surfaced under its own
        # name for direct access without reading through the spec).
        "n_frozen_if2_rejected": int(rejects["I-F2_mass_quad"]),
        "n_pp_runs": int(n_pp_runs),
        "wOcean_ppt_nominal": w_nominal,
        "w_axis_scope": ("salinity is UNDEFINED without an ocean; the frozen "
                         "axis is not crossed with w"),
        "Cuncertainty": float(frozen_Cuncertainty),
        "Cuncertainty_rationale": (
            "ruling 3.2 documented fallback. With Dhsphere_m=0 the window "
            "cannot bias WHICH admissible structure is selected (always the "
            "shallowest = the ice base), only whether the intended one "
            "survives; a node it excludes acquires liquid and is rejected by "
            "I-F2. C/MR² is derived and recorded, never matched."),
        "tolerances": {
            "I-F1_mass": float(frozen_mass_tol),
            "I-F2_mass_quad_window": (
                "one-sided, step-scaled: (-step_scale_frac, 0] per node "
                "(no single global constant; see "
                "_frozen_stored_mass_audit / mass_gate_blast_radius_audit"
                ".json's corrected signature)"),
            "I-F3_rho_closure_kgm3": float(frozen_rho_closure_tol_kgm3),
            "I-F4_zb_km": float(frozen_zb_tol_km),
        },
        "I-F6_window": [float(C_centre), float(C_half)],
        "I-F6_cmr2_span": [float(cmr2_built.min()), float(cmr2_built.max())],
        "M_body_kg": M_body_kg,
        "R_body_m": R_body_m,
    }
    return fz_arr, structs, spec
