"""Per-sample radial interior profiles for the inference globe panel.

The sample picker selects one posterior draw; this module rebuilds that
draw's full radial structure by pushing its parameter vector through the
SAME bilinear (Tb x w) cache interpolation the likelihood used
(`forward_models.apply_parameters`), then renders standard
PlanetProfile-style property profiles (rho, T, P, sigma or sound
speeds), a wedge diagram, and a PP-format radial data table. All figures are matplotlib and are meant
to be displayed through Utilities.crisp_figs (vector SVG + cached
exports).

The structure cache is the file named by
``result.config.structure_cache_path`` (shipped with the app for the
amortized slots). When it is missing the profile/table views degrade
to the wedge figure, which needs only the per-sample derived arrays
every result already carries.
"""
from __future__ import annotations

import pickle
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import streamlit as st

from Utilities.globe_view import LAYER_COLORS, LAYER_LABELS

# PlanetProfile phase codes -> globe layer kinds (colors/labels shared
# with the 3D view). Codes: ocean=0, ice I..VI = 1..6, clathrate=30,
# silicate=50 (+ice-in-pore variants 5x), iron core=100.
_PHASE_KINDS = (
    (lambda p: p == 0, 'ocean'),
    (lambda p: p == 1, 'ice_Ih'),
    (lambda p: p in (2, 3), 'ice_III'),   # II folded into III color band
    (lambda p: p == 5, 'ice_V'),
    (lambda p: p == 6, 'ice_VI'),
    (lambda p: 30 <= p < 40, 'ice_Ih'),   # clathrate: show with Ih band
    (lambda p: 50 <= p < 100, 'silicate'),
    (lambda p: p >= 100, 'core'),
)


def phase_kind(p) -> str:
    p = int(p)
    for test, kind in _PHASE_KINDS:
        if test(p):
            return kind
    return 'silicate'


@st.cache_resource(show_spinner=False)
def _load_structure_cache(path_str: str):
    with open(path_str, 'rb') as f:
        return pickle.load(f)


def resolve_cache_path(result, parent_directory) -> Optional[Path]:
    p = getattr(result.config, 'structure_cache_path', '') or ''
    if not p:
        return None
    cand = Path(p)
    if not cand.is_absolute():
        cand = Path(parent_directory) / p
    return cand if cand.exists() else None


def _theta_for(result, sel_idx: Optional[int]) -> Dict[str, float]:
    samples = np.asarray(result.samples, float)
    names = list(result.param_names)
    row = (samples[sel_idx] if sel_idx is not None
           else np.median(samples, axis=0))
    return {n: float(row[j]) for j, n in enumerate(names)}


def _expanded_theta(result, sel_idx: Optional[int], structure_data
                    ) -> Tuple[Optional[Dict[str, float]], str]:
    """Full theta dict for one draw: mirror the runner's _expand_theta
    (group broadcast + fixed params) and derive Tb_K for the v5/v6
    thick-ice reparameterization by inverting the cached D_iceIh(Tb, w)
    field. Returns (theta, '') or (None, reason)."""
    theta = _theta_for(result, sel_idx)
    cfg = result.config
    for _gkey, _members in (getattr(cfg, 'param_groups', {}) or {}).items():
        if _gkey in theta:
            for _m in _members:
                theta[_m] = theta[_gkey]
    for _k, _v in (getattr(cfg, 'fixed_params', {}) or {}).items():
        try:
            theta[_k] = float(_v)
        except (TypeError, ValueError):
            pass
    if 'D_iceIh_km' in theta and 'Tb_K' not in theta:
        try:
            from PlanetProfile.Inference.grid_interp_2d import (
                invert_d_iceIh_to_Tb, wOcean_ppt_from_theta)
            d_flat = [(None if s is None else s.get('D_iceIh_km'))
                      for s in structure_data.get('structures', [])]
            Tb = invert_d_iceIh_to_Tb(
                np.asarray(structure_data['Tb_K_grid'], float),
                np.asarray(structure_data['wOcean_ppt_grid'], float),
                d_flat, float(theta['D_iceIh_km']),
                wOcean_ppt_from_theta(theta))
            if Tb is None:
                return None, ("selected draw's ice thickness is outside "
                              "the cached (Tb, w) band — no structure")
            theta['Tb_K'] = float(Tb)
        except Exception as e:
            return None, f"D_iceIh → Tb inversion unavailable: {e}"
    return theta, ''


def _arrhenius_of(result):
    cfg = result.config
    return (getattr(cfg, 'arrhenius_params', None)
            or (getattr(cfg, 'sampler_settings', {}) or {}).get(
                'arrhenius_params') or None)


def profile_for_sample(result, sel_idx: Optional[int], parent_directory
                       ) -> Tuple[Optional[Dict], str]:
    """Radial structure for one posterior draw (or the marginal median
    when sel_idx is None). Returns (profile dict, note) — profile None
    with the reason in note when the cache can't serve it."""
    cache_path = resolve_cache_path(result, parent_directory)
    if cache_path is None:
        return None, ("structure cache not found for this result "
                      "(config.structure_cache_path missing or not "
                      "shipped) — layer stack only")
    try:
        structure_data = _load_structure_cache(str(cache_path))
    except Exception as e:
        return None, f"structure cache unreadable: {e}"

    theta, _terr = _expanded_theta(result, sel_idx, structure_data)
    if theta is None:
        return None, _terr
    cfg = result.config

    try:
        from PlanetProfile.Inference.forward_models import (
            apply_parameters, UnservableSampleError)
        modified = apply_parameters(theta, structure_data)
    except UnservableSampleError:
        return None, ("selected draw falls in an unbuilt corner of the "
                      "(Tb, w) cache — no structure")
    except Exception as e:
        return None, f"structure interpolation failed: {e}"

    r_m = modified.get('r_m')
    if r_m is None:
        return None, ("cache did not yield a radial structure for this "
                      "parameterization")

    # Effective viscosity — the SAME eta the tidal-heating forward model
    # uses (forward_model_k2_flexible): the cache's base profile, replaced
    # per-phase by any sampled log10_eta_* (hook writes 'eta_Pa'), then
    # Arrhenius-scaled eta(T) = eta_ref exp(E/R (1/T - 1/T_ref)) with
    # T_ref = the draw's basal temperature when the run configures
    # arrhenius_params.
    try:
        _eta = modified.get('eta_Pa')
        if _eta is None and modified.get('eta_Pa_base') is not None:
            _eta = np.asarray(modified['eta_Pa_base'], float).copy()
        _ap = (getattr(cfg, 'arrhenius_params', None)
               or (getattr(cfg, 'sampler_settings', {}) or {}).get(
                   'arrhenius_params'))
        if (_eta is not None and _ap
                and all(k in modified for k in
                        ('phases', 'changeIndices', 'n_layers'))):
            from PlanetProfile.Inference.forward_models import (
                apply_arrhenius_viscosity)
            _ape = dict(_ap)
            if ('Tb_K_sampled' in modified
                    and 'reference_temp_K' not in _ap):
                _ape['reference_temp_K'] = float(modified['Tb_K_sampled'])
            _eta = apply_arrhenius_viscosity(
                np.asarray(_eta, float).copy(), modified['phases'],
                modified['changeIndices'], modified['n_layers'],
                modified.get('T_K'), _ape)
        if _eta is not None:
            modified = dict(modified)
            modified['eta_Pa'] = _eta
    except Exception:
        pass  # fall back to the base profile below

    order = np.argsort(np.asarray(r_m, float))  # plot outward
    prof = {'r_km': np.asarray(r_m, float)[order] / 1e3}
    for key in ('rho', 'T_K', 'P_MPa', 'sigma_Sm', 'eta_Pa_base',
                'eta_Pa', 'mu_Pa', 'K_Pa', 'phases'):
        v = modified.get(key)
        prof[key] = (np.asarray(v, float)[order]
                     if v is not None and np.size(v) == order.size else None)
    # Seismic velocities from the cached elastic moduli:
    # VP = sqrt((K + 4/3 mu)/rho), VS = sqrt(mu/rho) (liquid: mu=0 -> VS=0)
    rho, mu, K = prof.get('rho'), prof.get('mu_Pa'), prof.get('K_Pa')
    if rho is not None and mu is not None and K is not None:
        with np.errstate(invalid='ignore', divide='ignore'):
            prof['VP_kms'] = np.sqrt((K + 4.0 * mu / 3.0)
                                     / np.where(rho > 0, rho, np.nan)) / 1e3
            prof['VS_kms'] = np.sqrt(mu / np.where(rho > 0, rho, np.nan)) / 1e3
    else:
        prof['VP_kms'] = prof['VS_kms'] = None
    _augment_thermo(prof)
    prof['Tb_K'] = theta.get('Tb_K')
    return prof, ''


_G_GRAV = 6.674e-11
# phase code -> SeaFreeze material for Cp/alpha (H2O phases only; the
# ocean uses the pure-water spline — a few-percent approximation at
# seawater salinities. Silicate/core thermodynamics are not in the
# inference cache and stay nan.)
_SF_PHASE = {0: 'water1', 1: 'Ih', 2: 'II', 3: 'III', 5: 'V', 6: 'VI'}


def _augment_thermo(prof):
    """Fill the derivable PP-profile quantities the cache omits:
    g(r), shell volumes/masses (from the density profile), and Cp/alpha
    (SeaFreeze, H2O phases)."""
    r_m = np.asarray(prof['r_km'], float) * 1e3   # ascending
    n = r_m.size
    rho = prof.get('rho')
    if rho is not None:
        r3 = r_m ** 3
        dV = 4.0 / 3.0 * np.pi * np.diff(r3)
        V0 = 4.0 / 3.0 * np.pi * r3[0]
        rho_mid = 0.5 * (np.asarray(rho)[1:] + np.asarray(rho)[:-1])
        dM = np.concatenate([[V0 * rho[0]], dV * rho_mid])
        with np.errstate(invalid='ignore', divide='ignore'):
            prof['g_ms2'] = _G_GRAV * np.cumsum(dM) / r_m ** 2
        prof['VLayer_m3'] = np.concatenate([[V0], dV])
        prof['MLayer_kg'] = dM
    P, T, ph = prof.get('P_MPa'), prof.get('T_K'), prof.get('phases')
    _missing = [k for k, v in (('P_MPa', P), ('T_K', T), ('phases', ph))
                if v is None]
    if _missing:
        prof['thermo_notes'] = (
            'Cp/α skipped: cache lacks ' + ', '.join(_missing)
            + ' on the radial grid (absent or length-mismatched)')
    if P is not None and T is not None and ph is not None:
        Cp = np.full(n, np.nan)
        al = np.full(n, np.nan)
        notes = []
        try:
            import seafreeze.seafreeze as sf
            phi = np.asarray(ph, float).astype(int)
            for code, mat in _SF_PHASE.items():
                m = (phi == code) & np.isfinite(P) & np.isfinite(T)
                if not m.any():
                    continue
                pts = np.array(list(zip(np.asarray(P, float)[m],
                                        np.asarray(T, float)[m])),
                               dtype='f,f').astype(object)
                try:
                    out = sf.getProp(pts, mat)
                    Cp[m] = np.ravel(out.Cp)
                    al[m] = np.ravel(out.alpha)
                except Exception as e:
                    notes.append(f'SeaFreeze {mat}: {e}')
            # Clathrate layers (phase 30-39): PP's own parameterization
            # (ClathProps — Ning et al. 2015 Cp, thermal expansivity;
            # user 2026-07-25: Titan's surface clathrate lid left the
            # visible top of the table empty). ClathProps meshgrids
            # (P, T); the paired-point values are its diagonal.
            mcl = ((phi >= 30) & (phi < 40)
                   & np.isfinite(P) & np.isfinite(T))
            if mcl.any():
                try:
                    from PlanetProfile.Thermodynamics.Clathrates. \
                        ClathrateProps import ClathProps
                    _, Cp_cl, al_cl, _ = ClathProps(
                        np.asarray(P, float)[mcl],
                        np.asarray(T, float)[mcl])
                    Cp[mcl] = np.diagonal(Cp_cl)
                    al[mcl] = np.diagonal(al_cl)
                except Exception as e:
                    notes.append(f'clathrate Cp/α unavailable: {e}')
        except Exception as e:
            notes.append(f'SeaFreeze unavailable: {e}')
        if not np.any(np.isfinite(Cp)) and not notes:
            codes = sorted(set(np.asarray(ph)[np.isfinite(ph)].astype(int)
                               .tolist()))
            notes.append('no H2O-phase layers matched the SeaFreeze set '
                         f'(phase codes present: {codes})')
        prof['Cp_JkgK'] = Cp
        prof['alpha_pK'] = al
        if notes:
            prof['thermo_notes'] = '; '.join(notes)


def _phase_bands(r_km, phases) -> List[Tuple[float, float, str]]:
    """Contiguous (r_lo, r_hi, kind) runs of one material, outward."""
    bands = []
    if phases is None:
        return bands
    kinds = [phase_kind(p) for p in phases]
    lo = r_km[0]
    for i in range(1, len(kinds) + 1):
        if i == len(kinds) or kinds[i] != kinds[i - 1]:
            bands.append((lo, r_km[i - 1] if i == len(kinds) else r_km[i],
                          kinds[i - 1]))
            if i < len(kinds):
                lo = r_km[i]
    return bands


def build_profile_figure(prof: Dict, title: str):
    """Standard-property profiles with RADIUS ON THE Y-AXIS, surface at
    the top of the figure (user 2026-07-23): density, temperature,
    pressure, viscosity and shear modulus (the tidal-heating rheology
    inputs, user 2026-07-24), and conductivity — or sound speeds when
    the cache carries no conductivity data."""
    import matplotlib.pyplot as plt

    r = prof['r_km']
    sig = prof.get('sigma_Sm')
    # "No conductivity data": nothing above the numerical floor the cache
    # uses for non-conducting layers (1e-16 S/m)
    has_sigma = sig is not None and np.any(
        np.isfinite(sig) & (np.asarray(sig) > 1e-12))
    if has_sigma:
        last = [('sigma_Sm', r'Conductivity $\sigma$ (S m$^{-1}$)',
                 'log', 1.0)]
    else:
        last = [('VP_kms', r'$V_P$, $V_S$ (km s$^{-1}$)', 'linear', 1.0)]
    eta_key = 'eta_Pa' if prof.get('eta_Pa') is not None else 'eta_Pa_base'
    panels = [
        ('rho', r'Density $\rho$ (kg m$^{-3}$)', 'linear', 1.0),
        ('T_K', r'Temperature $T$ (K)', 'linear', 1.0),
        ('P_MPa', r'Pressure $P$ (MPa)', 'linear', 1.0),
        (eta_key, r'Viscosity $\eta$ (Pa s)', 'log', 1.0),
        ('mu_Pa', r'Shear modulus $\mu$ (GPa)', 'linear', 1e-9),
    ] + last
    fig, axes = plt.subplots(1, len(panels), figsize=(15.5, 5.6),
                             sharey=True)
    bands = _phase_bands(r, prof.get('phases'))
    seen_kinds: List[str] = []
    for (key, xlab, xscale, xfac), ax in zip(panels, axes):
        for lo, hi, kind in bands:
            ax.axhspan(lo, hi, color=_rgb01(LAYER_COLORS[kind]),
                       alpha=0.25, linewidth=0)
            if kind not in seen_kinds:
                seen_kinds.append(kind)
        v = prof.get(key)
        if v is None:
            ax.text(0.5, 0.5, 'not in cache', transform=ax.transAxes,
                    ha='center', color='0.5')
        else:
            vv = np.asarray(v, float) * xfac
            if xscale == 'log':
                vv = np.where(vv > 0, vv, np.nan)
            ax.plot(vv, r, color='k', linewidth=1.4, label=r'$V_P$'
                    if key == 'VP_kms' else None)
            if key == 'VP_kms' and prof.get('VS_kms') is not None:
                ax.plot(prof['VS_kms'], r, color='k', linewidth=1.1,
                        linestyle='--', label=r'$V_S$')
                ax.legend(fontsize=8, loc='lower right')
            ax.set_xscale(xscale)
        ax.set_xlabel(xlab)
    axes[0].set_ylabel('Radius $r$ (km)')
    axes[0].set_ylim(0, float(np.nanmax(r)) * 1.02)  # surface at top
    handles = [plt.Rectangle((0, 0), 1, 1,
                             color=_rgb01(LAYER_COLORS[k]), alpha=0.45)
               for k in seen_kinds]
    fig.legend(handles, [LAYER_LABELS.get(k, k) for k in seen_kinds],
               loc='lower center', ncol=min(4, len(seen_kinds) or 1),
               fontsize=8, frameon=False)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    return fig


def heating_for_sample(result, sel_idx: Optional[int], parent_directory,
                       qrad_override_Wkg: Optional[float] = None,
                       qrad_source: Optional[str] = None
                       ) -> Tuple[Optional[Dict], str]:
    """Radial heating budget for one draw.

    Tidal: TidalPy radial volumetric heating from the SAME per-draw solve
    the k2/h2 likelihood used (forward_models._solve_draw via
    forward_model_radial_heating) — structure, rheology, Arrhenius
    viscosity and orbital parameters (e, omega, a) all identical.

    Radiogenic: a modeling ASSUMPTION, not an inferred quantity (the
    inference carries no radiogenic observable). Q_rad is the body's
    PlanetProfile config constant Sil.Qrad_Wkg unless qrad_override_Wkg
    is given (e.g. a McDonough+ 2020 radionuclide-inventory value chosen
    in the GUI); qrad_source labels the choice in the output dict.
    """
    cache_path = resolve_cache_path(result, parent_directory)
    if cache_path is None:
        return None, "structure cache not found for this result"
    try:
        structure_data = _load_structure_cache(str(cache_path))
    except Exception as e:
        return None, f"structure cache unreadable: {e}"
    theta, _terr = _expanded_theta(result, sel_idx, structure_data)
    if theta is None:
        return None, _terr

    prof, note = profile_for_sample(result, sel_idx, parent_directory)
    if prof is None:
        return None, note

    try:
        from PlanetProfile.Inference.forward_models import (
            forward_model_radial_heating)
        rh = forward_model_radial_heating(
            theta, structure_data, arrhenius_params=_arrhenius_of(result))
    except Exception as e:
        return None, f"tidal heating solve unavailable: {e}"
    if rh is None:
        return None, ("tidal radial solve failed for this draw (or the "
                      "cache carries no eccentricity)")
    r_solver_m, H_solver = rh

    r_km = np.asarray(prof['r_km'], float)
    r_m = r_km * 1e3
    H_tid = np.interp(r_m, r_solver_m, H_solver, left=0.0, right=0.0)
    P_tidal_W = float(np.trapezoid(
        H_solver * 4.0 * np.pi * r_solver_m ** 2, r_solver_m))

    heat = {'r_km': r_km, 'H_tidal_Wm3': H_tid, 'P_tidal_W': P_tidal_W,
            'phases': prof.get('phases')}
    notes = []
    ph = prof.get('phases')
    rho = prof.get('rho')
    Qrad = None
    if qrad_override_Wkg is not None and np.isfinite(qrad_override_Wkg):
        Qrad = float(qrad_override_Wkg)
        heat['Qrad_source'] = qrad_source or 'override'
    else:
        try:
            from Utilities.PlanetLoader import load_planet_module
            body = getattr(result.config, 'bodyname', '') or ''
            Qrad = float(load_planet_module(
                parent_directory, body).Planet.Sil.Qrad_Wkg)
            heat['Qrad_source'] = 'body config constant Sil.Qrad_Wkg'
        except Exception as e:
            notes.append(f'radiogenic skipped: body config unavailable ({e})')
    if Qrad is not None and ph is not None and rho is not None:
        sil = (np.asarray(ph) >= 50) & (np.asarray(ph) < 100)
        rho_sil_used = np.asarray(rho, float)
        # Mass-conservation configs (sampled core + derived rho_sil):
        # the draw's silicate density is the UNIFORM derived value, not
        # the cache's Perple_X profile — use it so radiogenic heating is
        # consistent with the inferred mantle density (user 2026-07-25).
        rho_sil_draw = None
        _dp = (getattr(result.config, 'derived_params', {}) or {}).get(
            'rho_sil_kgm3', {})
        if (_dp.get('derivation') == 'mass_conservation'
                and 'R_core_km' in theta and 'rho_core_kgm3' in theta):
            try:
                from PlanetProfile.Inference.forward_models import (
                    apply_parameters)
                from PlanetProfile.Inference.structure_derivation import (
                    extract_hydrosphere_layers, derive_silicate_density)
                _struct = apply_parameters(dict(theta), structure_data)
                _hl, _R_ob_m, _M_hyd = extract_hydrosphere_layers(_struct)
                _R_c_m = float(theta['R_core_km']) * 1e3
                _rho_c = float(theta['rho_core_kgm3'])
                _rs, _ok = derive_silicate_density(
                    float(_struct['Mtot_kg']), _M_hyd, _R_ob_m,
                    _R_c_m, _rho_c,
                    bounds=tuple(_dp.get('bounds', (2200.0, 3500.0))))
                if _ok and np.isfinite(_rs):
                    rho_sil_draw = float(_rs)
                    rho_sil_used = np.where(sil, rho_sil_draw,
                                            rho_sil_used)
                    heat['P_rad_W'] = float(
                        (float(_struct['Mtot_kg']) - _M_hyd
                         - 4.0 / 3.0 * np.pi * _R_c_m ** 3 * _rho_c)
                        * Qrad)
                    heat['rho_sil_draw_kgm3'] = rho_sil_draw
            except Exception as e:
                notes.append(f'derived rho_sil unavailable for '
                             f'radiogenic scaling ({e}); using the '
                             'cache density profile instead')
        # Zero OUTSIDE the silicate shell by assumption: U/Th/K reside
        # in the rock; ices/ocean carry none and the Fe core is taken
        # radionuclide-free.
        H_rad = np.where(sil, rho_sil_used * Qrad, 0.0)
        heat['H_rad_Wm3'] = H_rad
        heat['Qrad_Wkg'] = Qrad
        if heat.get('P_rad_W') is None:
            ml = prof.get('MLayer_kg')
            if ml is not None:
                heat['P_rad_W'] = float(np.sum(np.asarray(ml)[sil])
                                        * Qrad)
        heat['H_total_Wm3'] = H_tid + H_rad
    else:
        heat['H_rad_Wm3'] = None
        heat['H_total_Wm3'] = H_tid
    heat['notes'] = notes
    return heat, ''


# Shipped Perple_X silicate tables (PlanetProfile user cache) and
# human-readable composition labels. The Fe-S core table is excluded:
# this check concerns the rocky mantle only.
PERPLEX_SILICATE_TABLES = (
    ('CV3hy1wt_678_1.tab',
     'CV3 chondrite, 1 wt% hydrated'),
    ('CI_undifferentiated_hhph_DEW17_nofluid_nomelt_685.tab',
     'CI chondrite, undifferentiated (hydrous HP phases, DEW17)'),
    ('CM_undifferentiated_hhph_DEW17_nofluid_nomelt_685.tab',
     'CM chondrite, undifferentiated (hydrous HP phases, DEW17)'),
    ('CV_undifferentiated_v4_687_DEW17_nofluid_nomelt_v2.tab',
     'CV chondrite, undifferentiated (DEW17)'),
    ('Comet_67P_CG_v7_excluding_fluid_properties.tab',
     'Comet 67P/C-G bulk composition'),
    ('CM_hydrous_differentiated_Ganymede_Core080Fe020S'
     '_excluding_fluid_properties.tab',
     'CM hydrous, differentiated (80/20 Fe/S core removed)'),
)


def mineralogy_for_sample(result, sel_idx: Optional[int], parent_directory,
                          phiTop_max_frac: float = 0.5,
                          Pclosure_MPa: float = 350.0,
                          tol_pct: float = 2.0
                          ) -> Tuple[Optional[Dict], str]:
    """Post-hoc Perple_X mineralogy consistency check for one draw.

    The inference constrains a BULK silicate density (uniform derived
    value for mass-conservation configs, or a sampled rho_sil). Each
    candidate Perple_X table gives a GRAIN density rho_grain(P, T) along
    this draw's silicate P, T profile. Porosity is the free parameter
    connecting them: with vacuum pores and the Han (2014) closure law
    phi(P) = phiTop * exp(-6.15 P / Pclosure), mass balance over the
    silicate shell is LINEAR in phiTop:

        phiTop = (sum(rho_g V) - rho_bulk sum(V))
                 / sum(rho_g exp(-6.15 P / Pclosure) V)

    Classification (scientific review 2026-07-28): the primary,
    well-conditioned quantity is mismatch_pct = (rho_grain_mean -
    rho_bulk)/rho_bulk. Verdicts use a symmetric ±tol_pct band (grain
    densities carry ~1-2% thermodynamic + 1-3% cold-edge uncertainty):
    'too light' below -tol_pct (no porosity can fix a deficit),
    'consistent' up to +tol_pct plus the POROSITY HEADROOM
    f_lever*phiTop_max (the most porosity can lighten the bulk, with
    f_lever the fraction of the grain-density integral shallow enough
    for pores to survive), 'too dense' beyond. phiTop itself has gain
    1/f_lever (~200x for deep mantles) and is reported only when
    identifiable (f_lever >= 0.08); phiTop_identifiable flags it.

    Display-only; feeds nothing back into the inference. Assumptions:
    VACUUM pores — the conservative/limiting case (maximum lightening
    per unit porosity, so 'too dense' exclusions are as generous as
    possible; water-filled pores would need ~1.4x more porosity), which
    also makes Peff = P exact (Ppore = 0). Grain density clamps to each
    table's grid edge outside its P, T range; layers within 25 K of a
    table's cold edge are counted in n_cold_edge (several tables show
    anomalous positive drho/dT there — assemblage-boundary/spline-edge
    artifacts).
    """
    prof, note = profile_for_sample(result, sel_idx, parent_directory)
    if prof is None:
        return None, note
    ph = prof.get('phases')
    P = prof.get('P_MPa')
    T = prof.get('T_K')
    V = prof.get('VLayer_m3')
    if any(v is None for v in (ph, P, T, V)):
        return None, 'draw profile lacks P, T, or layer volumes'
    ph = np.asarray(ph)
    sil = (ph >= 50) & (ph < 100)
    if not np.any(sil):
        return None, 'no silicate layers in this draw'
    Ps = np.asarray(P, float)[sil]
    Ts = np.asarray(T, float)[sil]
    Vs = np.asarray(V, float)[sil]
    ok = np.isfinite(Ps) & np.isfinite(Ts) & np.isfinite(Vs) & (Vs > 0)
    if not np.any(ok):
        return None, 'silicate P, T profile is not finite for this draw'
    Ps, Ts, Vs = Ps[ok], Ts[ok], Vs[ok]

    # Bulk silicate density this draw asserts (same precedence as the
    # heating tab): derived mass-conservation value > sampled parameter
    # > cache profile volume-weighted mean.
    rho_bulk = None
    rho_source = ''
    cache_path = resolve_cache_path(result, parent_directory)
    structure_data = None
    theta = None
    if cache_path is not None:
        try:
            structure_data = _load_structure_cache(str(cache_path))
            theta, _ = _expanded_theta(result, sel_idx, structure_data)
        except Exception:
            theta = None
    _dp = (getattr(result.config, 'derived_params', {}) or {}).get(
        'rho_sil_kgm3', {})
    if (theta is not None and structure_data is not None
            and _dp.get('derivation') == 'mass_conservation'
            and 'R_core_km' in theta and 'rho_core_kgm3' in theta):
        try:
            from PlanetProfile.Inference.forward_models import (
                apply_parameters)
            from PlanetProfile.Inference.structure_derivation import (
                extract_hydrosphere_layers, derive_silicate_density)
            _struct = apply_parameters(dict(theta), structure_data)
            _hl, _R_ob_m, _M_hyd = extract_hydrosphere_layers(_struct)
            _rs, _ok = derive_silicate_density(
                float(_struct['Mtot_kg']), _M_hyd, _R_ob_m,
                float(theta['R_core_km']) * 1e3,
                float(theta['rho_core_kgm3']),
                bounds=tuple(_dp.get('bounds', (2200.0, 3500.0))))
            if _ok and np.isfinite(_rs):
                rho_bulk = float(_rs)
                rho_source = ('mass-conservation derived rho_sil '
                              '(uniform by construction)')
        except Exception:
            pass
    if rho_bulk is None and theta is not None \
            and 'rho_sil_kgm3' in theta:
        rho_bulk = float(theta['rho_sil_kgm3'])
        rho_source = 'sampled rho_sil parameter'
    if rho_bulk is None:
        rho_prof = np.asarray(prof.get('rho'), float)[sil][ok]
        rho_bulk = float(np.sum(rho_prof * Vs) / np.sum(Vs))
        rho_source = ("volume-weighted mean of the cache's silicate "
                      "density profile")

    rows = []
    notes = []
    try:
        from PlanetProfile.Thermodynamics.InnerEOS import GetInnerEOS
    except Exception as e:
        return None, f'Perple_X EOS loader unavailable: {e}'
    for fname, label in PERPLEX_SILICATE_TABLES:
        try:
            eos = GetInnerEOS(fname)
            rho_g = np.asarray(
                eos.fn_rho_kgm3(Ps, Ts, grid=False), float)
        except Exception as e:
            notes.append(f'{label}: table unavailable ({e})')
            continue
        good = np.isfinite(rho_g) & (rho_g > 0)
        if not np.any(good):
            notes.append(f'{label}: no valid grain density on this '
                         'P, T profile')
            continue
        w = Vs[good]
        rg = rho_g[good]
        pp = Ps[good]
        num = float(np.sum(rg * w))
        rho_grain_mean = num / float(np.sum(w))
        # Porosity leverage: fraction of the grain-density integral
        # that sits shallow enough (P < ~Pclosure/6.15 e-folds) for
        # Han (2014) porosity to act on. Deep mantles have almost none,
        # so porosity can only lighten the bulk by f_lever * phiTop
        # ('porosity headroom').
        den = float(np.sum(
            rg * np.exp(-6.15 * pp / Pclosure_MPa) * w))
        f_lever = den / num if num > 0 else 0.0
        mismatch_pct = (rho_grain_mean - rho_bulk) / rho_bulk * 100.0
        headroom_pct = f_lever * phiTop_max_frac * 100.0
        phiTop = (num - rho_bulk * float(np.sum(w))) / den \
            if f_lever > 1e-6 else np.nan
        identifiable = f_lever >= 0.08
        # EOSwrapper only proxies Tmin for ice/ocean types; reach the
        # underlying PerplexEOSStruct for the table's cold edge
        try:
            from PlanetProfile.Utilities.defineStructs import EOSlist
            Tmin_tab = float(EOSlist.loaded[eos.key].Tmin)
        except Exception:
            Tmin_tab = np.nan
        n_cold = int(np.sum(Ts[good] < Tmin_tab + 25.0)) \
            if np.isfinite(Tmin_tab) else 0
        if mismatch_pct < -tol_pct:
            status = (f'too light by {-mismatch_pct:.1f}% '
                      f'(beyond the ±{tol_pct:.0f}% tolerance; no '
                      'porosity can fix a deficit)')
        elif mismatch_pct <= tol_pct + headroom_pct:
            status = f'consistent within ±{tol_pct:.0f}%'
            if identifiable and np.isfinite(phiTop) and phiTop > 0:
                status += (f' (phi_top ~ {min(phiTop, phiTop_max_frac):.0%}'
                           ', vacuum pores)')
        elif f_lever <= 1e-6:
            status = (f'too dense by {mismatch_pct:.1f}% — porosity '
                      'cannot act at these pressures')
        else:
            status = (f'too dense by {mismatch_pct:.1f}% — exceeds '
                      f'the ±{tol_pct:.0f}% tolerance plus the '
                      f'{headroom_pct:.2f}% porosity headroom')
        rows.append({'table': fname, 'label': label,
                     'rho_grain_mean_kgm3': rho_grain_mean,
                     'mismatch_pct': mismatch_pct,
                     'f_lever_frac': f_lever,
                     'headroom_pct': headroom_pct,
                     'phiTop_frac': float(phiTop)
                     if np.isfinite(phiTop) else np.nan,
                     'phiTop_identifiable': bool(identifiable),
                     'n_cold_edge': n_cold,
                     'status': status})
    if not rows:
        return None, 'no Perple_X table could be evaluated: ' \
            + '; '.join(notes)
    out = {'rho_bulk_kgm3': rho_bulk, 'rho_source': rho_source,
           'P_sil_range_MPa': (float(Ps.min()), float(Ps.max())),
           'T_sil_range_K': (float(Ts.min()), float(Ts.max())),
           'Pclosure_MPa': Pclosure_MPa,
           'phiTop_max_frac': phiTop_max_frac,
           'tol_pct': tol_pct,
           'rows': rows, 'notes': notes}
    return out, ''


def build_heating_figure(heat: Dict, title: str):
    """Radial heat production: volumetric rates (log) and the cumulative
    power they integrate to, radius on Y with the surface at the top
    (matching the property profiles)."""
    import matplotlib.pyplot as plt

    r = heat['r_km']
    r_m = np.asarray(r, float) * 1e3
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 5.6), sharey=True)
    bands = _phase_bands(r, heat.get('phases'))
    seen: List[str] = []
    for ax in (ax1, ax2):
        for lo, hi, kind in bands:
            ax.axhspan(lo, hi, color=_rgb01(LAYER_COLORS[kind]),
                       alpha=0.2, linewidth=0)
            if kind not in seen:
                seen.append(kind)

    def _cum_W(H):
        integrand = np.asarray(H, float) * 4.0 * np.pi * r_m ** 2
        out = np.zeros_like(r_m)
        out[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1])
                            * np.diff(r_m))
        return out

    series = [('H_tidal_Wm3', 'tidal', 'crimson', '-'),
              ('H_rad_Wm3', 'radiogenic', 'steelblue', '-'),
              ('H_total_Wm3', 'total', 'k', '--')]
    # Cumulative power in named units (GW/TW), not a 1e12 axis
    # multiplier (user 2026-07-25).
    _pmax = max((np.nanmax(_cum_W(heat[k])) for k, *_ in series
                 if heat.get(k) is not None), default=0.0)
    unit, fac = (('TW', 1e12) if _pmax >= 1e12
                 else ('GW', 1e9) if _pmax >= 1e9 else ('MW', 1e6))
    for key, lab, color, ls in series:
        H = heat.get(key)
        if H is None:
            continue
        Hp = np.where(np.asarray(H, float) > 0, H, np.nan)
        ax1.plot(Hp, r, color=color, linestyle=ls, linewidth=1.4,
                 label=lab)
        ax2.plot(_cum_W(H) / fac, r, color=color, linestyle=ls,
                 linewidth=1.4)
    ax1.set_xscale('log')
    ax1.set_xlabel(r'Volumetric heating $H$ (W m$^{-3}$)')
    ax1.set_ylabel('Radius $r$ (km)')
    ax1.legend(fontsize=8, loc='lower left')
    ax2.set_xlabel(r'Cumulative power $\int_0^r H\,4\pi r^2\,dr$'
                   f' ({unit})')
    ax1.set_ylim(0, float(np.nanmax(r)) * 1.02)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    return fig


def build_geotherm_figure(prof: Dict, title: str):
    """Geotherm T(P) from the draw's cached radial structure: pressure on
    Y increasing downward (surface at top), phase bands along the path.
    The mantle segment is where the mineralogy assumption (mantle EOS
    table, evaluated by Perple_X at cache-build time) lives."""
    import matplotlib.pyplot as plt

    T = prof.get('T_K')
    P = prof.get('P_MPa')
    if T is None or P is None:
        return None
    T = np.asarray(T, float)
    P = np.asarray(P, float)
    r = np.asarray(prof['r_km'], float)
    order = np.argsort(-r)  # surface first
    T, P, r_s = T[order], P[order], r[order]
    ph = prof.get('phases')
    ph_s = np.asarray(ph)[order] if ph is not None else None

    fig, ax = plt.subplots(figsize=(6.4, 5.8))
    # Phase bands as P-intervals along the path
    if ph_s is not None:
        kinds = [phase_kind(p) for p in ph_s]
        lo_i = 0
        seen: List[str] = []
        for i in range(1, len(kinds) + 1):
            if i == len(kinds) or kinds[i] != kinds[i - 1]:
                p_lo = P[lo_i]
                p_hi = P[i - 1] if i == len(kinds) else P[i]
                ax.axhspan(min(p_lo, p_hi), max(p_lo, p_hi),
                           color=_rgb01(LAYER_COLORS[kinds[i - 1]]),
                           alpha=0.2, linewidth=0)
                if kinds[i - 1] not in seen:
                    seen.append(kinds[i - 1])
                if i < len(kinds):
                    lo_i = i
        handles = [plt.Rectangle((0, 0), 1, 1,
                                 color=_rgb01(LAYER_COLORS[k]), alpha=0.45)
                   for k in seen]
        ax.legend(handles, [LAYER_LABELS.get(k, k) for k in seen],
                  fontsize=8, loc='lower left', frameon=False)
    ax.plot(T, P, color='k', linewidth=1.6)
    ax.set_xlabel('Temperature $T$ (K)')
    ax.set_ylabel('Pressure $P$ (MPa)')
    ax.set_yscale('log')
    Ppos = P[P > 0]
    if Ppos.size:
        ax.set_ylim(float(np.nanmax(P)) * 1.1, float(Ppos.min()) * 0.9)
    ax.set_title(title, fontsize=11)
    fig.tight_layout()
    return fig


def build_wedge_figure(R_km: float, thicknesses: List[Tuple[str, float]],
                       title: str, wedge_angle_deg: float = 25.0):
    """Standard PlanetProfile-style wedge diagram (interior structure as
    an angular wedge, surface arc at top) built from per-sample layer
    thicknesses, using the globe layer colors. Mirrors the CLI PlotWedge
    convention: wedge opens upward, angles 90 +/- wedge_angle_deg."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Wedge

    ang1, ang2 = 90.0 - wedge_angle_deg, 90.0 + wedge_angle_deg
    fig, ax = plt.subplots(figsize=(5.6, 6.0))
    # outer radii, marching inward; draw outermost first so inner layers
    # overplot the wedge center region
    r_out = float(R_km)
    bands = []
    for kind, d in thicknesses:
        if d <= 0:
            continue
        bands.append((kind, r_out, d))
        r_out -= d
    for kind, r_o, d in bands:
        ax.add_patch(Wedge((0.5, 0), r_o / R_km, ang1, ang2,
                           fc=_rgb01(LAYER_COLORS[kind]),
                           ec='white', linewidth=0.6))
    # innermost fill (below the last labeled layer, e.g. residual core gap)
    if bands and (bands[-1][1] - bands[-1][2]) > 0:
        ax.add_patch(Wedge((0.5, 0), (bands[-1][1] - bands[-1][2]) / R_km,
                           ang1, ang2,
                           fc=_rgb01(LAYER_COLORS[bands[-1][0]]),
                           ec='none'))
    # Labels along the wedge centerline; thin layers get side callouts
    side_y = 1.02
    for kind, r_o, d in bands:
        y_mid = (r_o - d / 2) / R_km
        label = f"{LAYER_LABELS.get(kind, kind)}  {d:.0f} km"
        if d / R_km > 0.05:
            ax.text(0.5, y_mid, label, ha='center', va='center',
                    fontsize=8)
        else:
            x_edge = 0.5 + (r_o / R_km) * np.cos(np.radians(ang1))
            ax.plot([x_edge, 0.98], [y_mid, side_y], color='0.4',
                    linewidth=0.6)
            ax.text(0.99, side_y, label, ha='left', va='center',
                    fontsize=7)
            side_y -= 0.06
    ax.set_xlim(-0.05, 1.45)
    ax.set_ylim(-0.02, 1.1)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    return fig


def pp_wedge_exports(geo: Dict, parent_directory, bodyname: str,
                     w_ppt: Optional[float] = None,
                     rho_sil_free: bool = False):
    """Render the interior wedge with PlanetProfile's OWN PlotWedge
    (Plotting/ProfilePlots.py), not a lookalike: deepcopy the body's PP
    config Planet, inject the selected draw's layer geometry, and let the
    CLI plotting code draw it. Returns
    (svg_str, pdf_bytes, png_bytes, notes) where notes lists any
    honesty adjustments made (e.g. suppressed fictitious core).

    Inference carries no ice-convection partition, so the shell is drawn
    fully conductive (eLid = zb, Dconv = 0); mantle/core composition
    labels come from the body's PP defaults, not the inference.
    """
    import os
    import tempfile
    from copy import deepcopy
    from types import SimpleNamespace

    from Utilities.PlanetLoader import load_planet_module
    mod = load_planet_module(parent_directory, bodyname)
    Planet = deepcopy(mod.Planet)
    # Whether this BODY's PlanetProfile model actually includes an iron
    # core — captured before any overrides. Some inference configs sample
    # R_core_km for coreless bodies purely as a silicate-interior density
    # degree of freedom (mass-conservation rho_sil); drawing an Fe core
    # for those would assert structure the model does not claim
    # (user 2026-07-24, Titan free-gravity).
    default_has_core = bool(getattr(mod.Planet.Do, 'Fe_CORE', False))
    notes: List[str] = []

    R_km = float(geo['R_km'])
    d_ih = float(geo.get('d_ih', 0.0))
    d_oc = float(geo.get('d_oc', 0.0))
    d3, d5, d6 = (float(geo.get(k, 0.0)) for k in ('d3', 'd5', 'd6'))
    d_hs = float(geo.get('d_hs', 0.0))
    rc = float(geo.get('r_core', 0.0) or 0.0)
    if rc > 1.0 and not default_has_core:
        notes.append(
            f"the sampled core radius acts as a silicate-density degree "
            f"of freedom for {bodyname} (no iron core in its "
            "PlanetProfile model) — no core is drawn")
        rc = 0.0

    Planet.Bulk.R_m = R_km * 1e3
    Planet.zb_km = d_ih
    Planet.zIceI_m = 0.0
    Planet.dzIceI_km = d_ih
    Planet.eLid_m = d_ih * 1e3      # fully conductive (no convection info)
    Planet.Dconv_m = 0.0
    Planet.deltaTBL_m = 0.0
    Planet.D_km = d_oc
    z_km = d_ih + d_oc
    Planet.dzIceII_km = 0.0
    Planet.zIceII_m = 0.0
    Planet.dzIceIIIund_km = 0.0
    Planet.zIceIIIund_m = 0.0
    Planet.dzIceVund_km = 0.0
    Planet.zIceVund_m = 0.0
    Planet.zIceIII_m = z_km * 1e3
    Planet.dzIceIII_km = d3
    z_km += d3
    Planet.zIceV_m = z_km * 1e3
    Planet.dzIceV_km = d5
    z_km += d5
    Planet.zIceVI_m = z_km * 1e3
    Planet.dzIceVI_km = d6
    Planet.Sil.Rmean_m = (R_km - d_hs) * 1e3
    Planet.Core.Rmean_m = rc * 1e3
    Planet.Do.Fe_CORE = rc > 1.0
    # Frozen-hydrosphere draws: PlotWedge prints 'no ocean' instead of
    # the body-default ocean composition (which the inference does not
    # assert).
    Planet.Do.NO_OCEAN = d_oc <= 0.0
    Planet.dzFeS_km = 0.0
    Planet.dzSilPorous_km = 0.0
    Planet.Do.POROUS_ROCK = False
    Planet.Do.CLATHRATE = False
    Planet.Bulk.asymIce = None
    Planet.Magnetic.ionosBounds_m = None
    if w_ppt is not None and np.isfinite(w_ppt):
        Planet.Ocean.wOcean_ppt = float(w_ppt)
    # State WHICH mantle EOS database the structure came from (user
    # 2026-07-25): the cache build ran full PlanetProfile with the body
    # config's Perple_X-derived table, so name it rather than hide it —
    # with the caveat that a rho_sil-sampling run rescales the density.
    _meos = str(getattr(Planet.Sil, 'mantleEOS', '') or '')
    if _meos:
        notes.append(
            f"mantle structure from Perple_X EOS table '{_meos}' "
            "(evaluated at structure-cache build)"
            + ("; this run's sampled/derived silicate density RESCALES "
               "that structure — a self-consistent Perple_X refit "
               "(requires porosity) is roadmap" if rho_sil_free else ""))

    from PlanetProfile.Plotting import ProfilePlots as _PPplots
    import matplotlib.pyplot as plt
    outs = {}
    tmpdir = tempfile.mkdtemp(prefix='ppwedge_')
    old_fmt = _PPplots.FigMisc.figFormat
    old_transparent = getattr(_PPplots.FigMisc, 'TRANSPARENT', False)
    # The Inference page forces text.usetex=False for its crisp figures,
    # but when TeX is installed GetConfig builds FigLbl labels with
    # siunitx macros (\si{km}) that mathtext cannot parse — render the
    # wedge exactly as the CLI does (usetex iff TEX_INSTALLED).
    try:
        plt.rcParams['text.usetex'] = bool(
            getattr(_PPplots.FigMisc, 'TEX_INSTALLED', False))
        _PPplots.FigMisc.TRANSPARENT = False  # white background in the app
        for fmt in ('svg', 'pdf', 'png'):
            path = os.path.join(tmpdir, f'wedge.{fmt}')
            params = SimpleNamespace(
                ALL_ONE_BODY=True, TITLES=False,
                FigureFiles=SimpleNamespace(vwedg=path))
            _PPplots.FigMisc.figFormat = fmt
            _PPplots.PlotWedge([Planet], params)
            with open(path, 'r' if fmt == 'svg' else 'rb') as f:
                outs[fmt] = f.read()
    finally:
        _PPplots.FigMisc.figFormat = old_fmt
        _PPplots.FigMisc.TRANSPARENT = old_transparent
        # NOT old_usetex: the snapshot above ran AFTER the PlanetProfile
        # imports, which set usetex=True when TeX is installed — restoring
        # it would leave LaTeX on for the app figures built later in the
        # same rerun (Europa v4 u-panel crash, user 2026-07-24). App
        # policy: matplotlib text is mathtext-only outside this function.
        plt.rcParams['text.usetex'] = False
    return outs['svg'], outs['pdf'], outs['png'], notes


# --- PlanetProfile-form radial data table -------------------------------
# Column set and row format follow Main.WriteProfile so the download can
# be read like a standard PlanetProfile output profile. Quantities the
# inference cache does not carry are written as nan.
_PP_COLUMNS = [
    ('P (MPa)', 'P_MPa'), ('T (K)', 'T_K'), ('r (m)', '_r_m'),
    ('phase ID', 'phases'), ('rho (kg/m3)', 'rho'),
    ('Cp (J/kg/K)', 'Cp_JkgK'), ('alpha (1/K)', 'alpha_pK'),
    ('g (m/s2)', 'g_ms2'),
    ('phi (void/solid frac)', None), ('sigma (S/m)', 'sigma_Sm'),
    ('k (W/m/K)', None), ('VP (km/s)', 'VP_kms'), ('VS (km/s)', 'VS_kms'),
    ('QS', None), ('KS (GPa)', '_KS_GPa'), ('GS (GPa)', '_GS_GPa'),
    ('Ppore (MPa)', None), ('rhoMatrix (kg/m3)', None),
    ('rhoPore (kg/m3)', None), ('MLayer (kg)', 'MLayer_kg'),
    ('VLayer (m3)', 'VLayer_m3'), ('Htidal (W/m3)', None),
    ('eta (Pa s)', '_eta_eff'),
]


def profile_table(prof: Dict) -> "list[dict]":
    """PP-ordered (surface inward, P ascending) rows of the radial data."""
    idx = np.argsort(-np.asarray(prof['r_km']))  # surface first
    n = idx.size
    derived = {
        '_r_m': np.asarray(prof['r_km']) * 1e3,
        '_KS_GPa': (np.asarray(prof['K_Pa']) / 1e9
                    if prof.get('K_Pa') is not None else None),
        '_GS_GPa': (np.asarray(prof['mu_Pa']) / 1e9
                    if prof.get('mu_Pa') is not None else None),
        # Effective viscosity (sampled overrides + Arrhenius) when the
        # profile carries it; the cache's base profile otherwise.
        '_eta_eff': (prof['eta_Pa'] if prof.get('eta_Pa') is not None
                     else prof.get('eta_Pa_base')),
    }
    rows = []
    for i in idx:
        row = {}
        for name, key in _PP_COLUMNS:
            if key is None:
                v = np.nan
            else:
                arr = derived.get(key, prof.get(key))
                v = (float(np.asarray(arr)[i])
                     if arr is not None and np.size(arr) == n else np.nan)
            row[name] = int(v) if name == 'phase ID' and np.isfinite(v) \
                else v
        rows.append(row)
    return rows


def profile_txt(prof: Dict, bodyname: str, sample_label: str) -> str:
    """Render the radial table in Main.WriteProfile's on-disk format
    (24-char right-aligned %.17e columns, 8-char phase ID) with a short
    provenance header. nHeadLines counts as in the CLI writer so generic
    skip-header readers work."""
    header = [
        f'PlanetProfile inference profile ({bodyname}, {sample_label})',
        'Source = PlanetProfileApp inference globe panel '
        '(structure-cache interpolation of one posterior draw)',
        'Note = columns not carried by the inference cache are nan',
    ]
    if prof.get('Tb_K') is not None:
        header.append(f'Tb_K = {prof["Tb_K"]:.3f}')
    n_head = len(header) + 3
    header = [f'  nHeadLines = {n_head:d}'] + header
    col_widths = [24] * len(_PP_COLUMNS)
    col_widths[3] = 8  # phase ID
    col_line = ' '.join(
        (' ' + name).ljust(w) if j == 0 else name.ljust(w)
        for j, ((name, _), w) in enumerate(zip(_PP_COLUMNS, col_widths)))
    lines = [f'{bodyname} inference sample profile',
             '\n  '.join(header), col_line]
    for row in profile_table(prof):
        parts = []
        for (name, _), w in zip(_PP_COLUMNS, col_widths):
            v = row[name]
            if name == 'phase ID':
                parts.append(f'{int(v) if np.isfinite(float(v)) else 0:8d}')
            else:
                parts.append(f'{float(v):24.17e}')
        lines.append(' '.join(parts))
    return '\n'.join(lines) + '\n'


def _rgb01(rgb_str: str):
    """'rgb(r, g, b)' (plotly convention in globe_view) -> matplotlib tuple."""
    vals = rgb_str[rgb_str.index('(') + 1:rgb_str.index(')')].split(',')
    return tuple(int(v) / 255.0 for v in vals)
