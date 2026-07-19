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

    Returns
    -------
    A dict with keys: ``r_m``, ``rho``, ``K_Pa``, ``mu_Pa``, ``eta_Pa_base``,
    ``phases``, ``T_K``, ``P_MPa``, ``bulk_visc``, ``changeIndices``,
    ``n_layers``, ``layer_upper_radii``, ``layer_types``, ``region_phases``,
    ``omega``, ``eccentricity``, ``host_mass``, ``a_m``, ``R_body_m``,
    ``Mtot_kg``, ``CMR2``, ``Tb_K``, ``rhoSil_kgm3``, ``D_hsphere_km``,
    ``D_iceIh_km``, ``D_iceIII_km``, ``D_iceV_km``, ``D_iceVI_km``.
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

    Planet, Params = RunPP(Planet, configParams)
    Params.CALC_NEW_GRAVITY = True
    Planet, Params = SetupGravity(Planet, Params)

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
        "Mtot_kg": Planet.Bulk.M_kg,
        "CMR2": CMR2_pp,
        "J2": float(J2_pred),
        "C22": float(C22_pred),
        "Tb_K": float(Tb_K),
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
        # Ocean thickness = total hydrosphere − all ice phases. Required by
        # mcmc_plots.plot_layers_vs_docean and the per-sample D_ocean_km
        # extracted in the runner's recomputation pass. Without this key, the
        # plotter falls through to its `pt.get('D_ocean_km', 0.0)` default and
        # renders the (incorrect) "no-ocean: 100%" panel for Europa.
        "D_ocean_km": max(
            0.0,
            D_hsphere_km - D_iceIh_km - D_iceIII_km - D_iceV_km - D_iceVI_km,
        ),
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
            )
            above = _bisect_transition(
                planet_template_module, Tb_mid, Tb_hi,
                regions_mid, regions_hi, eps_T, max_iter,
                ocean_overrides=ocean_overrides,
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

    Returns
    -------
    The saved-to-disk dict::

        {'Tb_K_grid': np.ndarray (n_Tb,),
         'wOcean_ppt_grid': np.ndarray (n_w,),
         'structures': [structure_dict | None, ...],  # row-major, len n_Tb*n_w
         'schema_version': 'v3.0'}

    ``structures`` is row-major: entry ``i_Tb * n_w + i_w`` is the structure at
    ``(Tb_K_grid[i_Tb], wOcean_ppt_grid[i_w])``. Failed/rejected builds are
    stored as explicit ``None`` (never silently dropped) so the 2D lookup can
    reject those (Tb, w) corners. Each structure carries both ``Tb_K`` and
    ``wOcean_ppt``. The absence of ``wOcean_ppt_grid`` in a cache signals the
    legacy 1D path, so v1/v2 artifacts stay servable.
    """
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
                    ocean_overrides={"wOcean_ppt": float(w_ppt)},
                )
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
                    f"D_iceIh={struct.get('D_iceIh_km', 0.0):.1f} km"
                )

    if n_ok == 0:
        raise RuntimeError(
            f"All {total} (Tb, w) nodes failed; no cache to write."
        )

    cache = {
        "Tb_K_grid": Tb_arr,
        "wOcean_ppt_grid": w_arr,
        "structures": structures,
        "schema_version": "v3.0",
    }
    _save_cache_atomic(cache, output_path)
    if progress:
        log.info(
            f"2D cache written → {output_path} "
            f"({n_Tb} Tb × {n_w} w = {total} nodes: "
            f"{n_ok} built, {n_fail} None)"
        )
    return cache
