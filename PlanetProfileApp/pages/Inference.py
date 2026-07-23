"""
Inference Tab for PlanetProfileApp

Bayesian inference for constraining planetary interior parameters using
observed data (tidal Love numbers, magnetic induction, etc.).

Supports:
- MCMC (Markov Chain Monte Carlo) via pocoMC

Author: PlanetProfile Team
Date: 2026-04-29
"""
import os
import sys
import warnings
import streamlit as st
import numpy as np
from pathlib import Path
import time
import traceback

# nflows (sbi dependency) still calls the deprecated torch.triangular_solve
# internally; the warning is third-party noise, not actionable here.
warnings.filterwarnings(
    'ignore', message=r'torch\.triangular_solve is deprecated')


# ============================================================
# Path setup (critical for imports)
# ============================================================
this_file = os.path.abspath(__file__)
pages_directory = os.path.dirname(this_file)
app_directory = os.path.dirname(pages_directory)
parent_directory = os.path.dirname(app_directory)

if app_directory not in sys.path:
    sys.path.insert(0, app_directory)
if parent_directory not in sys.path:
    sys.path.insert(0, parent_directory)

# ============================================================
# Imports
# ============================================================
from Utilities.get_planet import get_planet
from Utilities.planet_sidebar import show_planet_status
# Public-serving detection (shared with the Exploreogram page; lives in
# Utilities because importing a pages/ module would execute the page).
from Utilities.app_mode import public_mode
from utils.mcmc_run_loader import (
    DEFAULT_RESULTS_ROOTS,
    scan_mcmc_results,
    load_mcmc_result,
    get_run_display_label,
)

# Lazy imports for PlanetProfile inference modules
def lazy_import_inference():
    """Lazy import to avoid startup cost if tab not used."""
    try:
        from PlanetProfile.Inference.parameter_registry import (
            PARAMETER_REGISTRY,
            PARAMETER_PRESETS,
            CATEGORY_LABELS,
            CATEGORY_ORDER,
            validate_parameter_combination,
            get_parameters_by_rheology
        )
        from PlanetProfile.Inference.inference_core import InferenceConfig
        from PlanetProfile.Inference.mcmc_runner import MCMCRunner
        return (PARAMETER_REGISTRY, PARAMETER_PRESETS, CATEGORY_LABELS,
                CATEGORY_ORDER, validate_parameter_combination,
                get_parameters_by_rheology, InferenceConfig, MCMCRunner)
    except ImportError as e:
        st.error(f"Failed to import PlanetProfile inference modules: {e}")
        return None

# ============================================================
# Session State Initialization
# ============================================================
def initialize_session_state():
    """Initialize session state variables for inference tab."""

    # (A stale-session migration used to wipe param_space whenever it
    # contained 'log10_zeta'. Removed 2026-07-18: modern configs — Test50
    # 8D, Europa 7D/Clipper — legitimately sample log10_zeta, and the
    # migration silently emptied any loaded config that used it.)

    # Control flags
    if 'inference_running' not in st.session_state:
        st.session_state.inference_running = False

    if 'inference_force_rerun' not in st.session_state:
        st.session_state.inference_force_rerun = False

    if 'inference_start_time' not in st.session_state:
        st.session_state.inference_start_time = None

    # Error tracking
    if 'inference_error' not in st.session_state:
        st.session_state.inference_error = None

    if 'inference_error_traceback' not in st.session_state:
        st.session_state.inference_error_traceback = None

    # Configuration
    if 'inference_mode' not in st.session_state:
        st.session_state.inference_mode = 'mcmc'

    if 'inference_preset' not in st.session_state:
        st.session_state.inference_preset = 'andrade_titan'

    if 'inference_custom_mode' not in st.session_state:
        st.session_state.inference_custom_mode = False

    if 'inference_selected_params' not in st.session_state:
        st.session_state.inference_selected_params = []

    if 'inference_param_space' not in st.session_state:
        st.session_state.inference_param_space = {}

    # Sync the MCMC observables to the GLOBALLY chosen planet: when the
    # user switches body in the sidebar, reseed the observable defaults
    # (and clear the instantiated widgets) for the new body. Triggered
    # only on a CHANGE of the global planet, so an explicitly applied
    # config for a different body is not clobbered on later reruns.
    _planet_now = str(st.session_state.get('Planet', '') or '')
    _planet_known = _planet_now.capitalize() in BODY_OBS_DEFAULTS
    if _planet_known and _planet_now != st.session_state.get(
            '_inference_planet_seen'):
        st.session_state['_inference_planet_seen'] = _planet_now
        _bd = BODY_OBS_DEFAULTS[_planet_now.capitalize()]
        st.session_state.inference_bodyname = _planet_now.capitalize()
        st.session_state.inference_observables = dict(_bd)
        for _wk in _OBS_WIDGET_KEYS:
            st.session_state.pop(_wk, None)

    if 'inference_observables' not in st.session_state:
        st.session_state.inference_observables = dict(
            BODY_OBS_DEFAULTS['Titan'])

    if 'inference_sampler_settings' not in st.session_state:
        st.session_state.inference_sampler_settings = {
            'n_effective': 500,
            'random_state': 42,
            'n_reeval': 500,
            'checkpoint_interval': 100,
            'max_iterations': 10000
        }

    if 'inference_structure_cache_path' not in st.session_state:
        st.session_state.inference_structure_cache_path = ''

    if 'inference_use_clathrate' not in st.session_state:
        st.session_state.inference_use_clathrate = True

    # Results caching
    if 'inference_results' not in st.session_state:
        st.session_state.inference_results = None

    if 'inference_cache_key' not in st.session_state:
        st.session_state.inference_cache_key = None

    # Progress tracking
    if 'inference_progress' not in st.session_state:
        st.session_state.inference_progress = {
            'iteration': 0,
            'n_total': 0,
            'n_samples': 0,
            'ess': 0,
            'acceptance_rate': None
        }


# ============================================================
# Config-file selector helpers
# ============================================================

def _load_available_configs():
    """Scan PlanetProfile/Inference/configs/*.json and return a dict.

    Returns
    -------
    dict
        ``{body_name: [(display_name, filepath), ...]}``
        sorted alphabetically by body and then by display name.
    """
    import json as _json
    import glob as _glob

    configs_dir = Path(parent_directory) / "PlanetProfile" / "Inference" / "configs"
    pattern = str(configs_dir / "*.json")
    files = sorted(_glob.glob(pattern))

    result = {}
    for fp in files:
        try:
            with open(fp) as fh:
                data = _json.load(fh)
        except Exception:
            continue

        body = data.get("bodyname") or data.get("body") or data.get("planet_name") or data.get("Body") or "Unknown"
        display = Path(fp).stem  # filename without .json
        result.setdefault(body, []).append((display, fp))

    # Sort bodies alphabetically; entries within each body already sorted by filename
    return dict(sorted(result.items()))


# Mission-appropriate observable defaults per body (fallbacks for the
# MCMC-mode observables panel; plans/mission-body-constraints.md). The old
# universal fallback was the Cassini-Titan pair — which silently baselined
# Titan k2 for every other body (user-reported 2026-07-20).
BODY_OBS_DEFAULTS = {
    'Titan': {'Re_k2': (0.608, 0.048), 'abs_Im_k2': (0.135, 0.035),
              'CMR2': (0.343, 0.001)},   # Cassini / Petricca et al. 2025
    'Europa': {'Re_k2': (0.23, 0.015), 'abs_Im_k2': (0.004, 0.015),
               'CMR2': (0.3547, 0.0024)},  # Mazarico proj. / GC21 MoI
    'Enceladus': {'Re_k2': (0.015, 0.02), 'abs_Im_k2': (0.005, 0.02),
                  'CMR2': (0.335, 0.003)},  # theory / Iess-Park context
    'Ganymede': {'Re_k2': (0.45, 0.10), 'abs_Im_k2': (0.008, 0.10),
                 'CMR2': (0.3115, 0.0028)},  # theory / Anderson 1996
    'Callisto': {'Re_k2': (0.3, 0.1), 'abs_Im_k2': (0.005, 0.1),
                 'CMR2': (0.3549, 0.0042)},  # theory / Anderson 2001
}

# Observable-panel widget keys: keyed-widget-state WINS over value= on
# rerun, so a config load must clear these for the loaded values to show
# (same reseed pattern as the amortized slot switch).
_OBS_WIDGET_KEYS = ('Re_k2_value', 'Re_k2_unc', 'Im_k2_value', 'Im_k2_unc',
                    'CMR2_value', 'CMR2_unc', 'use_k2_obs', 'use_cmr2_obs')


def render_config_file_selector():
    """Render a two-level body/config-variant selector above the run-config form.

    When the user picks a config the following session-state keys are populated
    from the JSON file so the existing form widgets pre-fill:
      - ``inference_preset``           (set to 'custom' so the form doesn't override)
      - ``inference_custom_mode``      (True — unlocks all form widgets)
      - ``inference_selected_params``  (keys of param_space in the config)
      - ``inference_param_space``      (full prior dicts)
      - ``inference_observables``      (from config 'observables', keys normalised)
      - ``inference_sampler_settings`` (merged with defaults)
      - ``inference_structure_cache_path``

    NOTE: the existing ``render_observables_config`` /
    ``render_sampler_settings`` / ``render_structure_cache_input`` functions read
    these same keys, so they will pre-fill after a config is loaded.  The
    ``render_prior_config`` widget reads ``inference_param_space`` for bounds, so
    that pre-fills too.  Fields that the form does *not* read from session-state
    (e.g. the body name used inside ``render_run_button`` which is derived from
    ``inference_preset``) are handled by setting ``inference_preset`` to a
    matching preset key when one exists, or to ``'custom'`` otherwise.
    """
    import json as _json

    configs_by_body = _load_available_configs()
    if not configs_by_body:
        st.info("No config files found in PlanetProfile/Inference/configs/.")
        return

    st.markdown("#### 📂 Load Config File")

    body_names = list(configs_by_body.keys())
    selected_body = st.selectbox(
        "Body:",
        options=body_names,
        key="cfg_selector_body",
    )

    variants = configs_by_body[selected_body]  # list of (display_name, filepath)
    variant_labels = [v[0] for v in variants]
    selected_variant_label = st.selectbox(
        "Config variant:",
        options=variant_labels,
        key="cfg_selector_variant",
    )

    # Resolve filepath
    selected_filepath = next(fp for (lbl, fp) in variants if lbl == selected_variant_label)

    # Load JSON for preview and button action
    try:
        with open(selected_filepath) as fh:
            cfg_data = _json.load(fh)
    except Exception as exc:
        st.error(f"Could not read config: {exc}")
        return

    # Config preview expander
    with st.expander("Config preview"):
        st.code(_json.dumps(cfg_data, indent=2), language="json")

    # Apply button
    if st.button("Apply config to form", key="cfg_selector_apply"):
        _apply_config_to_session_state(cfg_data)
        st.success(f"Loaded: {selected_variant_label}")
        st.rerun()


def _apply_config_to_session_state(cfg):
    """Populate inference session-state keys from a loaded config dict.

    Parameters
    ----------
    cfg : dict
        Parsed JSON config from PlanetProfile/Inference/configs/*.json
    """
    # --- param_space and selected_params ---
    param_space = cfg.get("param_space", {})
    if param_space:
        st.session_state.inference_param_space = dict(param_space)
        st.session_state.inference_selected_params = list(param_space.keys())

    # --- observables (normalise key names) ---
    raw_obs = cfg.get("observables", {})
    if raw_obs:
        normalised = {}
        for key, val in raw_obs.items():
            # Ensure value is a list/tuple of two floats
            if isinstance(val, (list, tuple)) and len(val) == 2:
                pair = (float(val[0]), float(val[1]))
            else:
                continue

            # Normalise Im_k2 -> abs_Im_k2
            if key == "Im_k2":
                normalised["abs_Im_k2"] = pair
            else:
                normalised[key] = pair
        if normalised:
            st.session_state.inference_observables = normalised
            # Clear the instantiated observable widgets so the loaded
            # config's values actually DISPLAY (keyed-widget-state wins
            # over value= — without this, e.g. Titan k2 defaults kept
            # showing after loading a Europa config; user 2026-07-20).
            for wk in _OBS_WIDGET_KEYS:
                st.session_state.pop(wk, None)

    # Record the body for mission-appropriate observable fallbacks.
    if cfg.get('bodyname'):
        st.session_state.inference_bodyname = str(cfg['bodyname'])

    # --- sampler settings (merge with existing defaults) ---
    raw_sampler = cfg.get("sampler_settings", {})
    if raw_sampler:
        current = dict(st.session_state.inference_sampler_settings)
        for field in ("n_effective", "n_reeval", "random_state",
                      "checkpoint_interval", "max_iterations"):
            if field in raw_sampler:
                current[field] = raw_sampler[field]
        st.session_state.inference_sampler_settings = current

    # Also honour top-level random_state
    if "random_state" in cfg and "random_state" not in cfg.get("sampler_settings", {}):
        s = dict(st.session_state.inference_sampler_settings)
        s["random_state"] = cfg["random_state"]
        st.session_state.inference_sampler_settings = s

    # --- structure cache path ---
    cache_path = cfg.get("structure_cache_path", "")
    if cache_path:
        st.session_state.inference_structure_cache_path = cache_path

    # --- preset / custom mode ---
    # Try to match config bodyname to a known preset so the run button can
    # resolve bodyname correctly.  Fall back to 'custom' if no match.
    bodyname = cfg.get("bodyname", "").lower()
    if "titan" in bodyname:
        # Pick the 8D noocean preset when the cache path suggests it, else andrade_titan
        cache = cfg.get("structure_cache_path", "")
        if "noocean" in cache or "allice" in cache:
            st.session_state.inference_preset = "andrade_titan_noocean_8D"
        else:
            st.session_state.inference_preset = "andrade_titan"
    elif "europa" in bodyname:
        st.session_state.inference_preset = "andrade_europa"
    else:
        # No matching preset — use custom mode so the form stays editable
        # and inference_preset is not silently wrong.
        st.session_state.inference_preset = "custom"

    # Always switch to custom mode so the full form is editable and no preset
    # silently overwrites the values we just loaded.
    st.session_state.inference_custom_mode = True

    # (The old preset radio and its widget-sync guard were removed
    # 2026-07-18: presets were redundant with this loader. inference_preset
    # (set above) still carries the matched preset key for bodyname
    # resolution in the run button.)

    # Same widget-state-wins rule applies to the keyed form inputs: drop
    # their stored widget state so each re-initializes from the values we
    # just loaded (value= is only honored when the key is absent).
    for wkey in ('Re_k2_value', 'Re_k2_unc', 'Im_k2_value', 'Im_k2_unc',
                 'CMR2_value', 'CMR2_unc', 'use_cmr2_obs', 'use_k2_obs',
                 'cache_path_input'):
        st.session_state.pop(wkey, None)


# ============================================================
# UI Render Functions
# ============================================================


def render_parameter_config(PARAMETER_REGISTRY, CATEGORY_LABELS, CATEGORY_ORDER, validate_parameter_combination):
    """Render dynamic parameter configuration UI."""

    if not st.session_state.inference_custom_mode:
        # Show read-only parameter list for preset
        st.markdown("#### 🔧 Selected Parameters")

        param_list = st.session_state.inference_selected_params
        for param_id in param_list:
            param_def = PARAMETER_REGISTRY[param_id]
            with st.expander(f"{param_def.label} ({param_def.latex_label})"):
                st.markdown(f"**Description:** {param_def.description}")
                st.markdown(f"**Category:** {param_def.category}")
                if param_def.units:
                    st.markdown(f"**Units:** {param_def.units}")

                # Show default prior
                st.markdown(f"**Default Prior:** {param_def.default_prior}")
                st.markdown(f"**Default Bounds:** {param_def.default_bounds}")

        return  # Exit early for preset mode

    # Custom mode: full parameter selection UI
    st.markdown("#### 🔧 Parameter Selection")

    # Multi-select for parameters
    param_options = {pid: f"{p.label} ({p.latex_label})"
                     for pid, p in PARAMETER_REGISTRY.items()}

    selected = st.multiselect(
        "Select parameters to infer:",
        options=list(param_options.keys()),
        default=st.session_state.inference_selected_params,
        format_func=lambda x: param_options[x],
        key='param_multiselect'
    )

    # Update session state
    st.session_state.inference_selected_params = selected

    # Validate parameter combination
    if selected:
        validation = validate_parameter_combination(selected)

        if not validation['valid']:
            st.error("❌ **Invalid parameter combination:**")
            for warning in validation['warnings']:
                st.error(f"- {warning}")
        else:
            # Show inferred rheology
            rheology = validation['rheology']
            if rheology:
                st.success(f"✅ Rheology: {rheology.capitalize()}")

            # Show warnings (non-blocking)
            for warning in validation['warnings']:
                st.warning(warning)

        # Show requires_rebuild warning if any param needs structure grid
        if validation['requires_rebuild']:
            st.warning("""
            ⚠️ **Structure Grid Required:** One or more selected parameters require
            a pre-computed structure grid. Generate grid with:
            ```bash
            python scripts/prepare_structure_variants.py
            ```
            """)


def render_prior_config(PARAMETER_REGISTRY):
    """Render prior configuration for selected parameters."""

    if not st.session_state.inference_selected_params:
        st.info("Select parameters above to configure priors.")
        return

    st.markdown("#### 📊 Prior Configuration")

    # Body-specific recommendations
    BODY_PRIOR_HINTS = {
        "Titan": {
            "log10_eta_Ih": "For Titan, Ice Ih viscosity is typically expected between $10^{13}$ and $10^{15}$ Pa·s.",
            "log10_eta_sil": "Titan's silicate mantle is likely partially hydrated, implying $10^{18}$ to $10^{20}$ Pa·s.",
            "alpha": "Andrade $\\alpha$ for ices is often taken as $1/3$ (0.33)."
        },
        "Europa": {
            "log10_eta_Ih": "Europa's shell is likely warmer/thinner, $10^{13}$ to $10^{14}$ Pa·s is common.",
            "log10_mu_Ih": "Maxwell shear modulus for Ice Ih is ~3.5 GPa ($\\log_{10}\\mu \\approx 9.5$)."
        }
    }

    # Infer current body from preset
    current_body = "Titan" if "titan" in st.session_state.inference_preset else \
                   ("Europa" if "europa" in st.session_state.inference_preset else None)

    # Initialize param_space if empty
    if not st.session_state.inference_param_space:
        st.session_state.inference_param_space = {}

    # Render prior inputs for each parameter
    for param_id in st.session_state.inference_selected_params:
        param_def = PARAMETER_REGISTRY[param_id]

        with st.expander(f"{param_def.label} ({param_def.latex_label})", expanded=False):
            st.markdown(f"*{param_def.description}*")

            # Show hint if available for this body/parameter
            if current_body and param_id in BODY_PRIOR_HINTS[current_body]:
                st.caption(f"💡 **Recommendation:** {BODY_PRIOR_HINTS[current_body][param_id]}")

            # Prior type selector
            prior_type = st.selectbox(
                "Prior type:",
                options=['uniform', 'normal', 'log-uniform'],
                index=0 if param_def.default_prior == 'uniform' else
                      (1 if param_def.default_prior == 'normal' else 2),
                key=f'prior_type_{param_id}'
            )

            # Bounds or mean/std inputs
            if prior_type in ['uniform', 'log-uniform']:
                col1, col2 = st.columns(2)
                with col1:
                    low = st.number_input(
                        "Lower bound:",
                        value=float(param_def.default_bounds[0]),
                        format="%.4f",
                        key=f'lower_{param_id}'
                    )
                with col2:
                    high = st.number_input(
                        "Upper bound:",
                        value=float(param_def.default_bounds[1]),
                        format="%.4f",
                        key=f'upper_{param_id}'
                    )

                # Store in param_space
                st.session_state.inference_param_space[param_id] = {
                    'prior_type': prior_type,
                    'bounds': [low, high]
                }

            elif prior_type == 'normal':
                col1, col2 = st.columns(2)
                with col1:
                    mean = st.number_input(
                        "Mean:",
                        value=float(param_def.default_mean) if param_def.default_mean else 0.0,
                        format="%.4f",
                        key=f'mean_{param_id}'
                    )
                with col2:
                    std = st.number_input(
                        "Std. deviation:",
                        value=float(param_def.default_std) if param_def.default_std else 1.0,
                        format="%.4f",
                        key=f'std_{param_id}'
                    )

                # Store in param_space
                st.session_state.inference_param_space[param_id] = {
                    'prior_type': prior_type,
                    'mean': mean,
                    'std': std
                }


def render_observables_config():
    """Render observable configuration (k2 values)."""
    st.markdown("#### 🎯 Observables")

    st.markdown("""
    Configure observed tidal Love number constraints. These are the data
    that the inference will try to match.
    """)

    # k2 inputs, checkbox-gated like CMR2. Not every loaded config carries a
    # k2 observable (e.g. callisto_mgso4 is CMR2-only); force-adding the
    # Titan display defaults to such a run would silently change its
    # likelihood, so absence of any k2 channel unchecks the box.
    _obs_now = st.session_state.inference_observables
    use_k2 = st.checkbox(
        "Include k₂ tidal-Love-number constraints",
        value=('Re_k2' in _obs_now or 'abs_Im_k2' in _obs_now
               or 'Im_k2' in _obs_now),
        key='use_k2_obs',
        help="Gaussian χ² terms for Re(k₂) and |Im(k₂)|."
    )
    Re_k2_value = Re_k2_uncertainty = None
    _body_defaults = BODY_OBS_DEFAULTS.get(
        st.session_state.get('inference_bodyname', 'Titan'),
        BODY_OBS_DEFAULTS['Titan'])
    re_k2_default = _obs_now.get('Re_k2', _body_defaults['Re_k2'])
    if use_k2:
        col1, col2 = st.columns(2)
        with col1:
            Re_k2_value = st.number_input(
                "Re(k₂) — Real part:",
                value=re_k2_default[0],
                format="%.4f",
                key='Re_k2_value'
            )
        with col2:
            Re_k2_uncertainty = st.number_input(
                "± Uncertainty:",
                value=re_k2_default[1],
                format="%.4f",
                key='Re_k2_unc'
            )

    # Im(k2) input - handle both old 'Im_k2' and new 'abs_Im_k2' keys for backward compatibility
    observables = st.session_state.inference_observables
    im_k2_default = observables.get(
        'abs_Im_k2', observables.get('Im_k2', _body_defaults['abs_Im_k2']))

    Im_k2_value = Im_k2_uncertainty = None
    if use_k2:
        col1, col2 = st.columns(2)
        with col1:
            Im_k2_value = st.number_input(
                "Im(k₂) — Imaginary part:",
                value=im_k2_default[0],
                format="%.4f",
                key='Im_k2_value'
            )
        with col2:
            Im_k2_uncertainty = st.number_input(
                "± Uncertainty:",
                value=im_k2_default[1],
                format="%.4f",
                key='Im_k2_unc'
            )

    # C/MR² (optional, checkbox-gated). No fallback tuple here: absence of the
    # key must leave the checkbox unchecked, so presets without a CMR2 term
    # (e.g. andrade_titan_noocean_8D) launch runs matching their config JSON.
    cmr2_default = st.session_state.inference_observables.get('CMR2')
    use_cmr2 = st.checkbox(
        "Include C/MR² moment-of-inertia constraint",
        value=cmr2_default is not None,
        key='use_cmr2_obs',
        help="Adds a Gaussian χ² term for the axial moment of inertia computed from the structure profile."
    )
    cmr2_obs = None
    if use_cmr2:
        col1, col2 = st.columns(2)
        with col1:
            cmr2_value = st.number_input(
                "C/MR² — Moment of inertia:",
                value=float(cmr2_default[0]) if cmr2_default else 0.343,
                format="%.4f",
                key='CMR2_value'
            )
        with col2:
            cmr2_uncertainty = st.number_input(
                "± Uncertainty:",
                value=float(cmr2_default[1]) if cmr2_default else 0.001,
                format="%.5f",
                key='CMR2_unc'
            )
        cmr2_obs = (cmr2_value, cmr2_uncertainty)

    # Update session state (use abs_Im_k2 as the canonical key). Start from
    # the existing dict so observables this form does not render (e.g. a
    # loaded config's Re_h2/Im_h2/Ae_* channels) survive the rerun instead
    # of being silently dropped from the run.
    new_observables = dict(st.session_state.inference_observables)
    new_observables.pop('Im_k2', None)  # superseded by canonical abs_Im_k2
    if use_k2:
        new_observables['Re_k2'] = (Re_k2_value, Re_k2_uncertainty)
        new_observables['abs_Im_k2'] = (Im_k2_value, Im_k2_uncertainty)
    else:
        new_observables.pop('Re_k2', None)
        new_observables.pop('abs_Im_k2', None)
    if cmr2_obs is not None:
        new_observables['CMR2'] = cmr2_obs
    else:
        new_observables.pop('CMR2', None)
    st.session_state.inference_observables = new_observables

    extra = [k for k in new_observables
             if k not in ('Re_k2', 'abs_Im_k2', 'CMR2')]
    if extra:
        st.caption(f"Additional observables from loaded config (kept in the "
                   f"run, edit via config JSON): {', '.join(extra)}")

    # Show reference
    st.caption("**Reference:** Petricca et al. (2025) *Nature* — Titan k₂ and C/MR² constraints")


def render_sampler_settings():
    """Render MCMC sampler settings."""
    st.markdown("#### ⚙️ Sampler Settings")

    # n_effective slider
    n_eff = st.slider(
        "Target effective sample size (ESS):",
        min_value=100,
        max_value=5000,
        value=st.session_state.inference_sampler_settings['n_effective'],
        step=100,
        key='n_eff_slider',
        help="Number of independent samples to collect. Higher = better convergence, longer runtime."
    )

    # n_reeval slider
    n_reeval = st.slider(
        "Heating re-evaluations:",
        min_value=100,
        max_value=1000,
        value=st.session_state.inference_sampler_settings['n_reeval'],
        step=100,
        key='n_reeval_slider',
        help="Number of posterior samples to re-compute full heating for."
    )

    # Advanced settings expander
    with st.expander("🔧 Advanced Settings"):
        random_state = st.number_input(
            "Random seed:",
            value=st.session_state.inference_sampler_settings['random_state'],
            format="%d",
            key='random_seed'
        )

        checkpoint_interval = st.number_input(
            "Checkpoint interval (samples):",
            value=st.session_state.inference_sampler_settings['checkpoint_interval'],
            format="%d",
            key='checkpoint_interval',
            help="Save progress every N samples (for long runs)."
        )

        max_iterations = st.number_input(
            "Max iterations:",
            value=st.session_state.inference_sampler_settings['max_iterations'],
            format="%d",
            key='max_iterations',
            help="Maximum MCMC steps before stopping (safety limit)."
        )

    # Update session state
    st.session_state.inference_sampler_settings = {
        'n_effective': n_eff,
        'n_reeval': n_reeval,
        'random_state': random_state,
        'checkpoint_interval': checkpoint_interval,
        'max_iterations': max_iterations
    }


def render_structure_cache_input():
    """Render structure cache path input with clathrate toggle."""
    st.markdown("#### 💾 Structure Cache")

    # Determine body-specific checkbox label
    preset = st.session_state.inference_preset
    if preset in ['andrade_titan', 'maxwell_titan']:
        clathrate_label = "Include clathrate cap (top layer)"
    elif preset == 'andrade_europa':
        clathrate_label = "Include clathrate underplate (bottom layer)"
    else:
        clathrate_label = "Include clathrate layer"

    # Clathrate toggle checkbox
    use_clathrate = st.checkbox(
        clathrate_label,
        value=st.session_state.inference_use_clathrate,
        key='clathrate_checkbox',
        help="Checked: *_clath.pkl | Unchecked: *_noclath.pkl"
    )

    # Update session state and regenerate path if toggled
    if use_clathrate != st.session_state.inference_use_clathrate:
        st.session_state.inference_use_clathrate = use_clathrate

        # Update cache path if it follows a known pattern
        current_path = st.session_state.inference_structure_cache_path
        if current_path:
            # Replace clath ↔ noclath in path
            if '_clath.pkl' in current_path:
                st.session_state.inference_structure_cache_path = current_path.replace(
                    '_clath.pkl', '_noclath.pkl'
                )
            elif '_noclath.pkl' in current_path:
                st.session_state.inference_structure_cache_path = current_path.replace(
                    '_noclath.pkl', '_clath.pkl'
                )

    # Show current path (read-only display for preset mode, editable for custom)
    if not st.session_state.inference_custom_mode:
        # Preset mode: show as info box (read-only)
        st.info(f"**Auto-configured path:** `{st.session_state.inference_structure_cache_path}`")

        # Validation
        if st.session_state.inference_structure_cache_path:
            full_path = Path(parent_directory) / st.session_state.inference_structure_cache_path
            if full_path.exists():
                st.success(f"✅ Cache file found ({full_path.stat().st_size / 1024:.1f} KB)")
            else:
                # NOTE: the andrade_europa entry (PPTest3 -> europa_cache/)
                # was removed 2026-07-12: it silently built a PureH2O
                # clathrate cache that mismatched the seawater Europa
                # campaign and lacked induction fields. Europa uses the
                # committed seawater grid via the config-file selector.
                _auto_gen_map = {
                    'andrade_titan': ('PlanetProfile.Test.PPTest41', 'titan_cache/', False, '1-3 minutes'),
                    'maxwell_titan': ('PlanetProfile.Test.PPTest42', 'titan_cache/', True, '15-30 minutes'),
                }
                preset = st.session_state.inference_preset
                gen_flag = f'_cache_gen_failed_{preset}'

                if preset in _auto_gen_map and not st.session_state.get(gen_flag):
                    test_module, output_dir, is_maxwell, est_time = _auto_gen_map[preset]
                    import subprocess
                    cmd = [
                        sys.executable, '-m',
                        'PlanetProfile.Inference.prepare_structure_variants',
                        '--test-module', test_module,
                        '--output-dir', output_dir,
                    ]
                    if is_maxwell:
                        cmd.append('--maxwell')
                    with st.spinner(f"Generating structure cache for {preset} — this takes {est_time}..."):
                        proc = subprocess.run(
                            cmd, capture_output=True, text=True,
                            cwd=parent_directory
                        )
                    if proc.returncode == 0:
                        st.rerun()
                    else:
                        st.session_state[gen_flag] = True
                        st.error("❌ Cache generation failed.")
                        st.code(proc.stderr or proc.stdout)
                else:
                    st.error(f"❌ Cache file not found: {full_path}")
                    st.markdown("""
                    **Generate structure cache manually:**
                    ```bash
                    python -m PlanetProfile.Inference.prepare_structure_variants \\
                        --test-module <module> --output-dir <dir>/
                    ```
                    """)
    else:
        # Custom mode: allow manual path input
        cache_path = st.text_input(
            "Path to structure cache (relative to project root):",
            value=st.session_state.inference_structure_cache_path,
            placeholder="titan_cache/titan_structure_clath.pkl",
            key='cache_path_input',
            help="Pre-computed planetary structure cache (see MCMC_INFERENCE_GUIDE.md)"
        )

        st.session_state.inference_structure_cache_path = cache_path

        # Validate path exists
        if cache_path:
            full_path = Path(parent_directory) / cache_path
            if full_path.exists():
                st.success(f"✅ Cache file found ({full_path.stat().st_size / 1024:.1f} KB)")
            else:
                st.error(f"❌ Cache file not found: {full_path}")
                st.markdown("""
                **To generate structure cache:**
                ```bash
                python scripts/prepare_structure_cache.py --body <Body> --cache-path <cache_dir>/
                ```
                """)


def validate_inference_config():
    """Validate configuration before running inference."""
    errors = []
    warnings = []

    # Check parameters selected
    if not st.session_state.inference_selected_params:
        errors.append("No parameters selected for inference.")

    # Check param_space configured
    if not st.session_state.inference_param_space:
        errors.append("Prior configuration incomplete.")

    # Check all selected params have prior config
    for param_id in st.session_state.inference_selected_params:
        if param_id not in st.session_state.inference_param_space:
            errors.append(f"Parameter '{param_id}' has no prior configuration.")

    # Check structure cache path
    if not st.session_state.inference_structure_cache_path:
        errors.append("Structure cache path not specified.")
    else:
        cache_path = Path(parent_directory) / st.session_state.inference_structure_cache_path
        if not cache_path.exists():
            errors.append(f"Structure cache file not found: {cache_path}")

    # Check observables
    if not st.session_state.inference_observables:
        errors.append("No observables configured.")

    # Warnings
    if st.session_state.inference_sampler_settings['n_effective'] < 500:
        warnings.append("n_effective < 500 may give poor convergence.")

    return {'valid': len(errors) == 0, 'errors': errors, 'warnings': warnings}


def render_run_button(InferenceConfig, MCMCRunner):
    """Render run button and execution logic."""
    st.markdown("---")

    # Validate configuration
    validation = validate_inference_config()

    # Show validation errors. Do NOT st.stop() here: that kills the whole
    # script, so an existing result in session state (results panel below)
    # silently vanished whenever the run-config form was incomplete — e.g.
    # a fresh Europa session with a loaded/completed run. Skip the run
    # button only; let the rest of the page (results) render.
    if not validation['valid']:
        st.error("❌ **Configuration errors:**")
        for error in validation['errors']:
            st.error(f"- {error}")
        return

    # Show validation warnings
    for warning in validation['warnings']:
        st.warning(f"⚠️ {warning}")

    # Run button
    if st.button("🚀 Run MCMC Inference", type="primary", use_container_width=True):
        st.session_state.inference_running = True
        st.session_state.inference_start_time = time.time()
        st.session_state.inference_error = None

        # Build InferenceConfig
        try:
            # Get body name from preset (not from session planet)
            preset = st.session_state.inference_preset
            if preset in ['andrade_titan', 'maxwell_titan', 'andrade_titan_noocean_8D',
                          'andrade_titan_noocean_diff_10D']:
                bodyname = 'Titan'
            elif preset == 'andrade_europa':
                bodyname = 'Europa'
            else:
                bodyname = 'Custom'

            config = InferenceConfig(
                mode='mcmc',
                bodyname=bodyname,
                param_space=st.session_state.inference_param_space,
                observables=st.session_state.inference_observables,
                sampler_settings=st.session_state.inference_sampler_settings,
                structure_cache_path=st.session_state.inference_structure_cache_path,
                random_state=st.session_state.inference_sampler_settings['random_state']
            )

            # Initialize runner
            runner = MCMCRunner(config)

            # Progress placeholder
            progress_placeholder = st.empty()
            status_placeholder = st.empty()

            # Progress callback
            def progress_callback(progress_dict):
                st.session_state.inference_progress = progress_dict

                # Update UI
                progress_placeholder.progress(
                    progress_dict['iteration'] / progress_dict['n_total'],
                    text=f"Iteration {progress_dict['iteration']} / {progress_dict['n_total']}"
                )

                ess = progress_dict.get('ess')
                acc = progress_dict.get('acceptance_rate')
                status_placeholder.markdown(f"""
                **Current Status:**
                - Samples: {progress_dict['n_samples']}
                - ESS: {f"{ess:.0f} / {config.sampler_settings.get('n_effective', '—')}" if ess is not None else "N/A"}
                - Acceptance: {f"{acc:.1%}" if acc is not None else "N/A"}
                """)

            # Run inference
            with st.spinner("Running MCMC inference... (typically 1-5 minutes for n_effective=500, longer for larger ESS targets)"):
                result = runner.run(progress_callback=progress_callback)

            # Store result
            st.session_state.inference_results = result
            st.session_state.inference_running = False

            # Show success
            elapsed = time.time() - st.session_state.inference_start_time
            st.success(f"✅ Inference complete in {elapsed/60:.1f} minutes!")
            st.balloons()

        except Exception as e:
            st.session_state.inference_error = str(e)
            st.session_state.inference_error_traceback = traceback.format_exc()
            st.session_state.inference_running = False

            st.error(f"❌ **Inference failed:** {e}")
            with st.expander("🔍 Full traceback"):
                st.code(st.session_state.inference_error_traceback)


# ============================================================
# Amortized (pretrained SBI) execution mode
# ============================================================
#
# The MCMC form's priors/sigmas are meaningful inputs to a sampler; for a
# pretrained NPE artifact they are FROZEN at training time. This mode
# therefore renders them read-only from the artifact's own metadata, keeps
# the observable VALUES live (that is what amortization buys), offers exact
# prior TRUNCATION (narrowing within the trained box = exact posterior
# conditioning for uniform priors), and routes through
# SBIRunner.infer_from_artifact — which never trains and reuses the MCMC
# forward model for k2/CMR2/thickness/heating recompute so every results
# panel below works identically.

# Slot registry: artifact file -> display + the structure cache its
# recompute needs (same cache the training config used).
_SBI_ARTIFACT_SLOTS = {
    'titan_andrade_noocean_posterior.pt': {
        'label': '1D · Cassini–Titan (Andrade, no ocean) — Test50 8D',
        'bodyname': 'Titan',
        'config_path': ('PlanetProfile/Inference/configs/'
                        'test50_titan_noocean_andrade_8D.json'),
        'cache_path': ('PlanetProfile/Test/mcmc_results/Titan/'
                       'Test50_andrade_noocean_yao2014/'
                       'titan_allice_yao2014_structure_grid.pkl'),
        'default_obs': {'Re_k2': 0.608, 'Im_k2': 0.135},
        # Hard validated-domain guard (user-ratified deployment condition,
        # 2026-07-11): the anchor gates validate conditioning only for
        # |Im k2| <= 0.20; above that the flow shows a directional low-eta
        # bias (onset ~0.18, W1 gate failure at 0.25). Runs outside these
        # limits are refused, not warned.
        'x_obs_limits': {'Im_k2': (0.0, 0.20)},
        'scope_note': ('Validated domain: |Im k2| <= 0.20 (SBC + crosscheck '
                       '+ 8/9 W1 anchors green). Known limitation: '
                       'directional low-viscosity bias in the bimodal '
                       'regime above Im k2 ~ 0.18; see sbi_artifacts/INDEX.md.'),
        # Cross-version sampling gate-validated on this machine (2026-07-13):
        # crosscheck vs the Test50 production MCMC PASSES all 8 params under
        # torch 2.8 (artifact trained on 2.11). See sbi_artifacts/INDEX.md.
        'validated_version_pairs': (('torch', '2.11.0', '2.8.0'),),
    },
    'europa_galileo_v1p1_8D_posterior_1m.pt': {
        'label': '1D · Galileo–Europa (Andrade, seawater) — v1.1 honest observables, 8D',
        'bodyname': 'Europa',
        'config_path': ('PlanetProfile/Inference/configs/'
                        'europa_galileo_v1p1_8D.json'),
        'cache_path': ('PlanetProfile/Test/mcmc_results/Europa/'
                       'Test51_seawater/europa_seawater_structure_grid.pkl'),
        # Honest Galileo-era data = CMR2 (GC21 MoI) + the synodic |Ae|>0.7
        # support cut. k2/h2 are LABELED HYPOTHETICAL channels at theory
        # widths (no Galileo measurement exists); zeta split ice/silicate.
        'default_obs': {'CMR2': 0.3547, 'Re_k2': 0.23, 'Im_k2': 0.004,
                        'Re_h2': 1.2, 'Im_h2': 0.0},
        # Conservative carryover from the v1 anchor sweep: no v1.1 anchor
        # walk was run; SBC (8/8, global) + fiducial crosscheck (8/8) are
        # the calibration evidence. Refuse far-off-fiducial Im k2.
        'x_obs_limits': {'Im_k2': (0.0, 0.15)},
        # Same synodic induction support edge as v1 (identical cut + cache):
        # no conductive ocean below ~261.5 K; restore the hard edge.
        'default_truncate': {'Tb_K': (261.5, None)},
        'scope_note': ('Galileo v1.1 (replaces v1): honest-data framing — '
                       'only C/MR² (Gomez Casajus et al. 2021) and the '
                       'synodic induction support cut are Galileo-era '
                       'measurements; the k2/h2 inputs are HYPOTHETICAL '
                       'exploration channels at theory widths (no Galileo '
                       'k2/h2 measurement exists). Independent ice/silicate '
                       'Andrade ζ. Gates: SBC 8/8, crosscheck 8/8 (no soft '
                       'fails); details sbi_artifacts/INDEX.md.'),
        # Trained torch 2.8.0 / sbi 0.26.1 — same runtime; no pair needed.
    },
    'europa_clipper_v4_geodesy_11D_posterior_1m.pt': {
        'label': ('1D · Clipper–Europa (Andrade, seawater) — v4 geodesy, '
                  '11D sampled salinity + non-hydrostatic gravity'),
        'bodyname': 'Europa',
        'config_path': ('PlanetProfile/Inference/configs/'
                        'europa_clipper_v4_geodesy_11D.json'),
        'cache_path': ('PlanetProfile/Test/mcmc_results/Europa/'
                       'Test52_seawater_v3/'
                       'europa_seawater_structure_grid_v3_2d.pkl'),
        # Mazarico-2023-projected k2 (sigma 0.015) + Park-style degree-2
        # gravity through the Clairaut hydrostatic map + sampled
        # non-hydrostatic offsets; CMR2 = the GC21 Galileo-derived MoI
        # prior; 14 Bind channels at 1.5 nT as v2/v3.
        'default_obs': {
            'CMR2': 0.3547,
            'C20': -4.579201e-4, 'C22': 1.377618e-4,
            'Re_k2': 0.23, 'Im_k2': 0.004,
            'Re_h2': 1.2, 'Im_h2': 0.0,
            'Bind_synodic_x_real': 91.8248, 'Bind_synodic_x_imag': -157.7708,
            'Bind_synodic_y_real': -59.4061, 'Bind_synodic_y_imag': -25.8874,
            'Bind_synodic_z_real': -5.6378, 'Bind_synodic_z_imag': -12.4064,
            'Bind_synodic 2nd_x_real': 14.7108, 'Bind_synodic 2nd_x_imag': 3.1779,
            'Bind_synodic 2nd_y_real': 1.8788, 'Bind_synodic 2nd_y_imag': -9.8686,
            'Bind_orbital_x_real': -0.3812, 'Bind_orbital_x_imag': 6.454,
            'Bind_orbital_y_real': -2.2124, 'Bind_orbital_y_imag': -0.2825,
        },
        'x_obs_limits': {'Im_k2': (0.0, 0.15)},
        # Inherited v2 grid-walk envelope (no v4 anchor walk; SBC-global +
        # fiducial gates are the v4 calibration evidence).
        'derived_ae_guards': [{
            'label': 'synodic', 'comp': 'x',
            'Be_comp': (128.466302323619, -170.084146795136),
            'ae_range': (0.75, 0.94),
            'warn_support_below': 0.7,
        }],
        # 2D (Tb, w) support — no 1D Tb truncation pre-applied.
        'scope_note': ('Clipper v4: 11 sampled parameters = v3 interior + '
                       'salinity PLUS two non-hydrostatic degree-2 gravity '
                       'offsets. REPORTABLE (calibrated): the identifiable '
                       'combination u = dC22_nh + dC20_nh/3.324 as an '
                       'upper limit, and the Tb–salinity joint. NOT '
                       'calibrated (do not cite): per-component '
                       'dC20_nh/dC22_nh marginals, the −C20/C22 ratio, and '
                       'the prior-dominated interior scalars ζ_Ih/ρ_core '
                       '(SBI conservatively wider). Re k2 = 0.23 sits a '
                       'mild 1.3σ below the joint posterior-predictive '
                       '(consistent; a v5 question, not a claim). Gate '
                       'details: sbi_artifacts/INDEX.md.'),
        # Trained torch 2.8.0 / sbi 0.26.1 — same runtime; no pair needed.
    },
    'europa_clipper_v5_geodesy_11D_posterior_1m.pt': {
        'label': ('1D · Clipper–Europa (Andrade, seawater) — v5 geodesy, '
                  '11D ice-thickness reparam + non-hydrostatic gravity'),
        'bodyname': 'Europa',
        'config_path': ('PlanetProfile/Inference/configs/'
                        'europa_clipper_v5_geodesy_11D.json'),
        'cache_path': ('PlanetProfile/Test/mcmc_results/Europa/'
                       'Test52_seawater_v5/'
                       'europa_seawater_structure_grid_v5_2d.pkl'),
        # v5 fiducial conditioning (config observables central values,
        # recomputed 2026-07-20 from the v5 cache). CMR2 retained here as the
        # GC21 Galileo MoI for the legacy display; the dual C/MR² readout
        # (actual structure integral vs RD-from-C22 hydrostatic reference) is
        # computed per sample regardless.
        'default_obs': {
            'CMR2': 0.3547,
            'C20': -4.578888e-4, 'C22': 1.377523e-4,
            'Re_k2': 0.23, 'Im_k2': 0.004,
            'Re_h2': 1.2, 'Im_h2': 0.0,
            'Bind_synodic_x_real': 91.07058126361886,
            'Bind_synodic_x_imag': -157.85676859240647,
            'Bind_synodic_y_real': -59.40533696629077,
            'Bind_synodic_y_imag': -25.61791852141546,
            'Bind_synodic_z_real': -5.675654206203697,
            'Bind_synodic_z_imag': -12.364182150469917,
            'Bind_synodic 2nd_x_real': 14.717630809871173,
            'Bind_synodic 2nd_x_imag': 3.1440689462941154,
            'Bind_synodic 2nd_y_real': 1.8561258037971171,
            'Bind_synodic 2nd_y_imag': -9.87259597254061,
            'Bind_orbital_x_real': -0.5352738290604214,
            'Bind_orbital_x_imag': 6.386436443696461,
            'Bind_orbital_y_real': -2.1855594672667733,
            'Bind_orbital_y_imag': -0.3339307437762303,
        },
        'x_obs_limits': {'Im_k2': (0.0, 0.15)},
        'derived_ae_guards': [{
            'label': 'synodic', 'comp': 'x',
            'Be_comp': (128.466302323619, -170.084146795136),
            'ae_range': (0.75, 0.94),
            'warn_support_below': 0.7,
        }],
        'scope_note': ('Clipper v5: v4 geodesy + D_iceIh ice-thickness '
                       'reparameterization (open uniform[5,80] km prior, '
                       'pivoted 2026-07-21). Dual C/MR² deliverable: ACTUAL '
                       '(structure moment integral) vs HYDROSTATIC REFERENCE '
                       '(Radau–Darwin from C22); the gap is the inferred '
                       'non-hydrostaticity, read against the ~0.0035 RD floor. '
                       'v5 ratification FAILED-this-pass (undertraining + '
                       'gate mis-spec, no model defect); do not cite as '
                       'ratified. Gate details: sbi_artifacts/INDEX.md.'),
    },
    # --- v5 channel-family siblings (reached via the channel selector, not
    # the dropdown). Same body/priors/cache as geodesy_11D; differ only in
    # which observable channels the flow conditions on. default_obs reuses
    # the geodesy central values — only the names present in each config's
    # observable set are read, so the extra keys are inert. ---
    'europa_clipper_v5_noinduction_7obs_posterior_1m.pt': {
        'label': ('1D · Clipper–Europa (Andrade, seawater) — v5 geodesy, '
                  '11D · channels: gravity + tidal k₂/h₂ (no induction)'),
        'bodyname': 'Europa',
        'config_path': ('PlanetProfile/Inference/configs/'
                        'europa_clipper_v5_noinduction_7obs.json'),
        'cache_path': ('PlanetProfile/Test/mcmc_results/Europa/'
                       'Test52_seawater_v5/'
                       'europa_seawater_structure_grid_v5_2d.pkl'),
        'default_obs': {
            'CMR2': 0.3547,
            'C20': -4.578888e-4, 'C22': 1.377523e-4,
            'Re_k2': 0.23, 'Im_k2': 0.004,
            'Re_h2': 1.2, 'Im_h2': 0.0,
        },
        'x_obs_limits': {'Im_k2': (0.0, 0.15)},
        'scope_note': ('Clipper v5 channel ablation: gravity + tidal k₂/h₂, '
                       'NO magnetic-induction channel (Galileo synodic support '
                       'still enters as a training-time cut). Same priors and '
                       '(Tb,w) cache as v5 geodesy. v5 ratification '
                       'FAILED-this-pass; do not cite as ratified.'),
    },
    'europa_clipper_v5_nok2_17obs_posterior_1m.pt': {
        'label': ('1D · Clipper–Europa (Andrade, seawater) — v5 geodesy, '
                  '11D · channels: gravity + induction (no tidal k₂/h₂)'),
        'bodyname': 'Europa',
        'config_path': ('PlanetProfile/Inference/configs/'
                        'europa_clipper_v5_nok2_17obs.json'),
        'cache_path': ('PlanetProfile/Test/mcmc_results/Europa/'
                       'Test52_seawater_v5/'
                       'europa_seawater_structure_grid_v5_2d.pkl'),
        'default_obs': {
            'CMR2': 0.3547,
            'C20': -4.578888e-4, 'C22': 1.377523e-4,
            'Bind_synodic_x_real': 91.07058126361886,
            'Bind_synodic_x_imag': -157.85676859240647,
            'Bind_synodic_y_real': -59.40533696629077,
            'Bind_synodic_y_imag': -25.61791852141546,
            'Bind_synodic_z_real': -5.675654206203697,
            'Bind_synodic_z_imag': -12.364182150469917,
            'Bind_synodic 2nd_x_real': 14.717630809871173,
            'Bind_synodic 2nd_x_imag': 3.1440689462941154,
            'Bind_synodic 2nd_y_real': 1.8561258037971171,
            'Bind_synodic 2nd_y_imag': -9.87259597254061,
            'Bind_orbital_x_real': -0.5352738290604214,
            'Bind_orbital_x_imag': 6.386436443696461,
            'Bind_orbital_y_real': -2.1855594672667733,
            'Bind_orbital_y_imag': -0.3339307437762303,
        },
        'derived_ae_guards': [{
            'label': 'synodic', 'comp': 'x',
            'Be_comp': (128.466302323619, -170.084146795136),
            'ae_range': (0.75, 0.94),
            'warn_support_below': 0.7,
        }],
        'scope_note': ('Clipper v5 channel ablation: gravity + magnetic '
                       'induction, NO tidal k₂/h₂ channel. Same priors and '
                       '(Tb,w) cache as v5 geodesy. v5 ratification '
                       'FAILED-this-pass; do not cite as ratified.'),
    },
}


# Channel-family groups: artifacts trained on the SAME body/priors/cache that
# differ ONLY in which observable channels they condition on. The amortized
# panel offers a mutually-exclusive channel selector that swaps between the
# members of a family (all three keep the static gravity channel CMR2/C20/C22;
# there is deliberately no "exclude gravity" member and no "exclude both"
# option — excluding both tides and induction is not a modelling case we
# support). Key = the 'all channels' (fullest) member; values map a channel
# choice to the sibling artifact filename.
_CHANNEL_FAMILIES = {
    'europa_clipper_v5_geodesy_11D_posterior_1m.pt': {
        'all': 'europa_clipper_v5_geodesy_11D_posterior_1m.pt',
        'no_induction': 'europa_clipper_v5_noinduction_7obs_posterior_1m.pt',
        'no_k2h2': 'europa_clipper_v5_nok2_17obs_posterior_1m.pt',
    },
}
# Reverse lookup: any family member -> its family key, so selecting a
# non-'all' member still resolves the full channel selector.
_CHANNEL_FAMILY_OF = {
    member: key
    for key, members in _CHANNEL_FAMILIES.items()
    for member in members.values()
}
_CHANNEL_CHOICE_LABELS = {
    'all': 'All channels',
    'no_induction': 'Exclude magnetic induction',
    'no_k2h2': 'Exclude tidal $k_2$/$h_2$',
}


@st.cache_data
def _europa_excitation_moments():
    """Complex excitation moments Be_{x,y,z}(label) [nT] from the MoonMag
    excitation table, keyed by the canonical config labels (the file's
    'adjusted orbital' row serves the 'orbital' label). Used to derive the
    covarying Bind vector components from a single Ae per frequency."""
    import csv
    path = (Path(parent_directory) / 'PlanetProfile' / 'MagneticInduction'
            / 'MoonMag' / 'excitation' / 'Be1xyz_Europa.txt')
    alias = {'adjusted orbital': 'orbital'}
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            name = alias.get(row['exc name'], row['exc name'])
            out[name] = {
                'x': complex(float(row['Bex_Re(nT)']), float(row['Bex_Im(nT)'])),
                'y': complex(float(row['Bey_Re(nT)']), float(row['Bey_Im(nT)'])),
                'z': complex(float(row['Bez_Re(nT)']), float(row['Bez_Im(nT)'])),
            }
    return out


def _sbi_artifacts_dir():
    return Path(parent_directory) / 'PlanetProfile' / 'Inference' / 'sbi_artifacts'


@st.cache_data(show_spinner=False)
def _body_gravity_constants(cache_path):
    """(omega, R_body_m, Mtot_kg) for the body of a structure cache.

    These are constant across the (Tb, w) grid, so we read them from the
    first structure node. Used for the live hydrostatic-C/MR² slider readout
    (Radau–Darwin needs the rotation parameter). Returns None on any failure
    so the readout degrades gracefully rather than crashing the config panel.
    """
    import pickle
    import numpy as _np
    try:
        p = Path(parent_directory) / cache_path
        with open(p, 'rb') as f:
            d = pickle.load(f)

        def _first(o, depth=0):
            if depth > 4:
                return None
            if isinstance(o, dict):
                if 'R_body_m' in o and 'omega' in o and 'Mtot_kg' in o:
                    return o
                for v in o.values():
                    r = _first(v, depth + 1)
                    if r:
                        return r
            if isinstance(o, (list, tuple)):
                for v in o:
                    if v is None:
                        continue
                    r = _first(v, depth + 1)
                    if r:
                        return r
            if isinstance(o, _np.ndarray):
                for v in o.tolist():
                    if v is None:
                        continue
                    r = _first(v, depth + 1)
                    if r:
                        return r
            return None

        s = _first(d)
        if not s:
            return None
        omega = float(s.get('omega'))
        R = float(s.get('R_body_m'))
        M = float(s.get('Mtot_kg'))
        if not all(map(_np.isfinite, (omega, R, M))):
            return None
        return (omega, R, M)
    except Exception:
        return None


# GC21 Table 2 SOL-A (unconstrained), UNNORMALIZED at 1565 km — the plausible
# Galileo-era band for the live readout shading. C22 = 138.62 ± 2.44 ×10⁻⁶.
_GC21_SOLA_C22 = (138.62e-6, 2.44e-6)


@st.fragment
def _render_live_hydrostatic_cmr2(x_obs, grav):
    """Live 'if hydrostatic' C/MR² from the current C22 slider value.

    Isolated in a fragment so it refreshes on C20/C22 edits without rerunning
    the whole page or touching the posterior (that stays behind the Generate
    button). ``grav`` = (omega, R_body_m, Mtot_kg) or None.
    """
    import numpy as _np
    from PlanetProfile.Inference.gravity_obs import cmr2_from_c22_rd
    c22 = x_obs.get('C22')
    if c22 is None:
        return
    if grav is None:
        st.caption("Live hydrostatic C/MR² readout unavailable "
                   "(body constants not found in the structure cache).")
        return
    omega, R, M = grav
    lo, hi = _GC21_SOLA_C22
    in_band = (lo - 2 * hi) <= c22 <= (lo + 2 * hi)
    try:
        cmr2_h = float(cmr2_from_c22_rd(float(c22), omega, R, M))
        if not _np.isfinite(cmr2_h):
            # RD returns NaN (not a raise) when k_f leaves its physical range.
            raise ValueError("non-finite RD C/MR²")
        c1, c2 = st.columns(2)
        c1.metric("Hydrostatic C/MR² (from current $C_{22}$)",
                  f"{cmr2_h:.5f}")
        band_lo = float(cmr2_from_c22_rd(lo - hi, omega, R, M))
        band_hi = float(cmr2_from_c22_rd(lo + hi, omega, R, M))
        c2.metric("GC21 SOL-A band (±1σ)",
                  f"{band_lo:.4f}–{band_hi:.4f}")
        if not in_band:
            st.caption(
                f"⚠️ $C_{{22}}$ = {c22:.3e} is **outside** the GC21 SOL-A "
                f"plausible band ({lo:.3e} ± {hi:.2e}); the readout still "
                "shows the corresponding hydrostatic C/MR² for exploration."
            )
        else:
            st.caption(
                "Hydrostatic reference: the C/MR² this $C_{22}$ implies **if** "
                "the body is in hydrostatic equilibrium (Radau–Darwin). The "
                "actual C/MR² is computed from the inferred interior after "
                "**Generate Posterior**; their gap is the non-hydrostaticity."
            )
    except (ValueError, ZeroDivisionError):
        st.caption(
            f"⚠️ $C_{{22}}$ = {c22:.3e} is Radau–Darwin-unphysical here "
            "(would require C/MR² > 2/3). Reduce $C_{22}$ toward the "
            "plausible band."
        )


@st.cache_resource(show_spinner="Loading SBI artifact...")
def _load_sbi_runner(artifact_path: str, _mtime: float, validated_pairs=()):
    """Load a pretrained artifact into a sampling-only SBIRunner, cached per
    (path, mtime) so a redeployed artifact invalidates the cache."""
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    return SBIRunner.load_artifact(
        artifact_path, validated_version_pairs=validated_pairs)


def _stale_library_banner():
    """Warn (loudly) if PlanetProfile library code changed on disk after this
    server imported it — Streamlit hot-reloads page files but NOT library
    modules, so results/signatures silently come from the old code."""
    import PlanetProfile.Inference.sbi_runner as _sr
    loaded_at = getattr(_sr, '_MODULE_LOADED_AT', None)
    # The stamp and this page ship in the same commit: a module without it
    # predates this page revision and is stale by construction.
    if loaded_at is None:
        st.error(
            "⚠️ **Stale server:** the in-memory PlanetProfile library "
            "predates this page's code. Library modules do not hot-reload — "
            "stop and restart the Streamlit server, then rerun."
        )
        return
    try:
        src_mtime = Path(_sr.__file__).stat().st_mtime
    except OSError:
        return
    if src_mtime > loaded_at:
        st.error(
            "⚠️ **Stale server:** `PlanetProfile/Inference/sbi_runner.py` "
            "changed on disk after this Streamlit server started. Library "
            "modules do not hot-reload — results below come from the OLD "
            "code. Stop and restart the Streamlit server, then rerun."
        )


# --- Amortized runtime estimate (rough figure of merit, no extra compute) ---
# Wall time is dominated by the forward-model recompute loop, which scales with
# the number of forward passes = 2·(derived-quantity subset) + (heating subset).
# The factor 2 is because each derived sample does two passes (k2/CMR2 + the
# true log-likelihood). The Posterior-samples count itself barely matters — the
# flow draw is ~instant and the recompute is capped by the derived subset — so
# the estimate keys off the two SUBSET sliders, which are the real cost levers.
#
# We do NOT run anything to estimate: seed a rough per-pass figure of merit for
# a typical laptop, then after the user's first real run replace it with the
# actually-observed per-pass time (result.metadata['elapsed_time_s']), stored in
# session_state per artifact. No EWMA, no probe runs — one plain overwrite.
_AMORT_OVERHEAD_S = 6.0          # rough fixed setup (flow load + packaging)
_AMORT_PERPASS_SEED_S = 0.008    # rough per-forward-pass, typical laptop


def _amort_n_derived_eff(n_post, n_reeval, n_derived):
    """Effective #samples the derived-quantity loop recomputes (mirrors
    sbi_runner.infer_from_artifact: superset of the heating subset, capped at
    n_post; None => all)."""
    if n_derived is None or int(n_derived) >= n_post:
        return n_post
    return min(max(int(n_derived), n_reeval), n_post)


def _amort_work_units(n_post, n_reeval, n_derived):
    """Forward-pass count driving wall time: 2·derived + heating."""
    n_der_eff = _amort_n_derived_eff(n_post, n_reeval, n_derived)
    return 2 * n_der_eff + n_reeval, n_der_eff


def _fmt_duration(seconds):
    if seconds < 90:
        return f"~{seconds:.0f} s"
    return f"~{seconds / 60:.1f} min"


def _estimate_amort_runtime(slot_key, n_post, n_reeval, n_derived):
    """Rough predicted wall time (s): overhead + per_pass·work_units. per_pass
    is the typical-laptop seed until the user's first run overwrites it with the
    observed value for this artifact."""
    work_units, n_der_eff = _amort_work_units(n_post, n_reeval, n_derived)
    per_pass = st.session_state.get(f'amort_perpass_{slot_key}',
                                    _AMORT_PERPASS_SEED_S)
    return _AMORT_OVERHEAD_S + per_pass * work_units, work_units, n_der_eff


def _record_amort_runtime(slot_key, result):
    """After a real run, overwrite the per-artifact per-pass figure with the
    observed value (no extra compute — just uses the run we already did)."""
    meta = result.metadata
    elapsed = meta.get('elapsed_time_s')
    if not elapsed:
        return
    # len() on numpy index arrays is safe; guard missing keys without truthiness
    # tests on arrays (an ndarray in a boolean context raises).
    di = meta.get('derived_indices')
    hi = meta.get('heating_indices')
    n_der_eff = len(di) if di is not None else 0
    n_heat = len(hi) if hi is not None else 0
    work_units = 2 * n_der_eff + n_heat
    if work_units <= 0:
        return
    measured = max(float(elapsed) - _AMORT_OVERHEAD_S, 0.0) / work_units
    st.session_state[f'amort_perpass_{slot_key}'] = measured


def render_amortized_config():
    """Amortized-mode configuration form. Returns a run-spec dict or None."""
    _stale_library_banner()
    available = []
    for fname, slot in _SBI_ARTIFACT_SLOTS.items():
        p = _sbi_artifacts_dir() / fname
        if not p.exists():
            continue
        # Collapse channel-family siblings: only the family key (the
        # 'all channels' member) appears in the dropdown; the reduced-channel
        # siblings are reached via the mutually-exclusive channel selector
        # rendered after selection, not as separate dropdown entries.
        fam = _CHANNEL_FAMILY_OF.get(fname)
        if fam is not None and fam != fname:
            continue
        available.append((slot['label'], fname))
    if not available:
        st.warning("No pretrained SBI artifacts found in "
                   "`PlanetProfile/Inference/sbi_artifacts/`. Train and "
                   "validate one via the SBI pipeline (see sbi_artifacts/"
                   "INDEX.md), then reload.")
        return None

    labels = [lbl for lbl, _ in available]
    choice = st.selectbox("Pretrained model:", labels, key='amort_artifact_choice')
    fname = dict(available)[choice]

    # Mutually-exclusive channel selector for artifacts that ship a
    # channel-family (same body/priors/cache, differ only in conditioned
    # channels). All members keep the static gravity channel; the choice
    # swaps between the full model and each single-channel-excluded sibling.
    # A member is only offered if its artifact actually exists on disk.
    fam = _CHANNEL_FAMILIES.get(fname)
    if fam:
        avail_choices = [c for c, mf in fam.items()
                         if (_sbi_artifacts_dir() / mf).exists()]
        if len(avail_choices) > 1:
            ch = st.segmented_control(
                "Conditioning channels:",
                options=avail_choices,
                format_func=lambda c: _CHANNEL_CHOICE_LABELS.get(c, c),
                default='all' if 'all' in avail_choices else avail_choices[0],
                key='amort_channel_choice',
                help="Swap between the sibling artifacts trained on the same "
                     "priors but different observable channels. All keep the "
                     "static gravity channel (CMR₂/C₂₀/C₂₂); this toggles the "
                     "tidal k₂/h₂ and magnetic-induction channels.")
            if ch is None:
                ch = 'all' if 'all' in avail_choices else avail_choices[0]
            fname = fam[ch]

    slot = _SBI_ARTIFACT_SLOTS[fname]
    artifact_path = _sbi_artifacts_dir() / fname

    # Namespace every observable/sigma/truncation widget key by the artifact
    # filename. Streamlit ignores a widget's value= once its key exists in
    # session_state, so slots that share an observable NAME (e.g. Re_k2 in
    # both the Titan and Clipper slots) would otherwise share the key
    # 'amort_obs_Re_k2' — whichever slot renders FIRST (Titan, 0.608) writes
    # that key and its value sticks when you switch slots. The per-transition
    # reseed below is not enough on its own: a hot-reload (e.g. after git pull)
    # keeps session_state while resetting the transition bookkeeping, so the
    # reseed never fires and the stale value leaks. Distinct keys per artifact
    # make the leak structurally impossible (user-reported 2026-07-20: Clipper
    # x_obs stuck at Titan's Re_k2=0.608 instead of Mazarico 0.23).
    slot_key = fname.replace('.', '_').replace(' ', '_')

    # Body-of-record banner. The dropdown has no body-aware default, so it
    # sits on the FIRST slot (Titan) on open — which silently seeds Titan's
    # k2 (Re_k2 = 0.608) into the observable defaults below. A user on
    # Europa then reads Titan's k2 as Europa's (user-reported 2026-07-20).
    # Ordering is intentionally unchanged; make the active body unmistakable
    # instead, and warn when it disagrees with the globally selected planet.
    _slot_body = slot['bodyname']
    # 'ChosenPlanet' is the canonical body-NAME string (path-validated, e.g.
    # "Europa"/"Titan"); 'Planet' is a PlanetStruct object, so str() on it
    # yields an object repr — never compare against that (it would make the
    # mismatch warning fire constantly with garbage text).
    _chosen = st.session_state.get('ChosenPlanet', None)
    _global_body = (str(_chosen).strip().capitalize()
                    if _chosen and not str(_chosen).startswith('--') else '')
    if _global_body and _global_body != _slot_body:
        st.warning(
            f"⚠️ You are viewing the **{_slot_body}** model, but your "
            f"globally selected body is **{_global_body}**. Every observable "
            f"default below — including k₂ — is a **{_slot_body}** value. "
            f"Switch the dropdown above to a {_global_body} model if that is "
            f"what you intend.")
    else:
        st.info(f"📍 Active model body: **{_slot_body}** — the observable "
                f"defaults below are {_slot_body} values.")

    # Reseed per-slot widget state when the artifact selection changes.
    # Keyed widget values win over value= defaults on rerun, so without
    # this the previous slot's inputs leak into the new slot wherever the
    # observable names overlap (user-visible 2026-07-13: switching from
    # Titan to Europa kept Titan's Re_k2=0.608 / Im_k2=0.135 and
    # conditioned the Europa flow far outside its physical k2 range,
    # producing a spurious silicate-dominated heating posterior).
    if st.session_state.get('amort_active_slot') != fname:
        st.session_state['amort_active_slot'] = fname
        for k in list(st.session_state.keys()):
            if k.startswith(('amort_obs_', 'amort_sigma_', 'amort_trunc_')):
                del st.session_state[k]

    # Slot-specific induction documentation (numbers sourced from each
    # training config; the generic model-assumptions expander stays
    # slot-agnostic because it renders before slot selection).
    with st.expander("🧲 How magnetic induction constrains this model"):
        if 'clipper' in fname:
            st.markdown(
                "This model conditions on **Europa Clipper-era induction "
                "measurements directly**: 14 channels "
                "`Bind_<frequency>_<component>_<re|im>` — the induced dipole "
                "coefficient expressed as equivalent surface field (nT), "
                "computed as the complex product `Ae(frequency) × Be_component` "
                "— at the synodic (11.23 h), 2nd synodic harmonic (5.62 h), "
                "and orbital (85.2 h) periods, each with σ = 1.5 nT (Europa "
                "Clipper magnetometer requirement, Kivelson et al. 2023). "
                "Imaginary parts are **signed** (not folded). In addition, "
                "the training set was cut to the Galileo-era support "
                "`|Ae_synodic| > 0.7`, `|Im Ae| < 0.4` — models without a "
                "conductive ocean are outside the flow's support. `Ae` "
                "depends on the ocean state only through the ice basal "
                "temperature T_b (fixed salinity in this artifact), so the "
                "induction channels constrain T_b, hence ocean and ice "
                "thickness.")
        else:
            st.markdown(
                "This model does **not** condition on an induction "
                "measurement channel. Instead, the Galileo-era synodic "
                "induction constraint enters as a hard **support cut at "
                "training time**: every training model satisfies "
                "`|Ae_synodic| > 0.7` and `|Im Ae| < 0.4` (a present, "
                "conductive ocean). Because `Ae` depends on the ocean state "
                "only through the ice basal temperature T_b in this model "
                "family, the cut removes T_b below ~261.5 K — which is why "
                "the T_b truncation below defaults to the induction-"
                "surviving range. The posterior Ae panel in the results "
                "shows where your posterior lands inside that support.")

    runner = _load_sbi_runner(str(artifact_path), artifact_path.stat().st_mtime,
                              slot.get('validated_version_pairs', ()))
    meta = runner.artifact_meta
    bounds = {n: (float(lo), float(hi)) for n, (lo, hi)
              in zip(runner.param_names, meta['param_bounds'])}
    noise_meta = (meta.get('rejection_stats') or {}).get('obs_noise') or {}
    train_sigma = noise_meta.get('sigma', {})

    st.caption(
        f"nsf/{meta.get('density_estimator')} · trained on "
        f"{meta.get('n_train_effective'):,} sims · config {meta.get('config_hash')} "
        f"· git {meta.get('git_sha')} · {meta.get('created_utc')}"
    )
    st.info("**Amortized mode:** priors and observation uncertainties are "
            "frozen at training time (shown read-only below). Observable "
            "*values* are free — conditioning is instant. Use MCMC mode to "
            "vary priors or σ freely. See `sbi_artifacts/INDEX.md` for the "
            "validated conditioning domain.")

    # --- Priors (read-only) + exact truncation ---
    st.markdown("#### 🔒 Trained priors")
    st.table({
        'parameter': list(bounds.keys()),
        'prior': [f"uniform [{lo:g}, {hi:g}]" for lo, hi in bounds.values()],
    })
    truncate_bounds = {}
    # Slot-specified default truncations (pre-applied ON) — e.g. restoring a
    # hard induction support edge the smooth flow cannot represent. None in a
    # bound means "use the trained edge" (no truncation on that side).
    default_trunc = slot.get('default_truncate') or {}
    has_default_trunc = any(name in default_trunc for name in bounds)
    with st.expander("✂️ Prior truncation (exact, within trained bounds)",
                     expanded=has_default_trunc):
        st.caption("Narrowing a uniform prior post hoc is exact (sample "
                   "filtering). Widening requires retraining.")
        if has_default_trunc:
            st.info("⚙️ This artifact ships with a **default truncation** "
                    "applied (validated-domain / support-edge correction). "
                    "Reset a slider to the full trained range to override.")
        for name, (lo, hi) in bounds.items():
            dlo, dhi = default_trunc.get(name, (None, None))
            init = (float(dlo) if dlo is not None else lo,
                    float(dhi) if dhi is not None else hi)
            sel = st.slider(name, min_value=lo, max_value=hi,
                            value=init, key=f'amort_trunc_{slot_key}_{name}')
            if sel != (lo, hi):
                truncate_bounds[name] = (float(sel[0]), float(sel[1]))

    # --- Observables: values live, sigmas locked ---
    st.markdown("#### 🔭 Observables (condition the posterior)")
    x_obs = {}
    obs_sigma = {}

    # Bind_ vector components COVARY: all components of one frequency are a
    # single complex Ae(frequency) times fixed excitation moments Be_comp,
    # so independent per-component inputs allow physically inconsistent
    # combinations. Default entry is therefore the PHYSICAL parameters —
    # complex Ae per frequency — and the components are derived on-manifold.
    # Raw per-component entry stays available for inverting actual
    # spacecraft measurement sets (where the components carry noise
    # independently).
    bind_names = [n for n in runner.obs_names if n.startswith('Bind_')]
    scalar_names = [n for n in runner.obs_names if not n.startswith('Bind_')]
    bind_raw = False
    if bind_names:
        bind_raw = st.checkbox(
            "Enter raw Bind vector components (only for inverting actual "
            "spacecraft measurements)", key='amort_bind_raw')

    # Physical (Radau–Darwin) clamp bounds for the static gravity coefficients.
    # C20/C22 are O(1e-4) — displayed and entered in scientific notation. Out-
    # of-range entries are auto-reverted by st.number_input to the nearest
    # bound (per user request). Bounds are the RD-physical extremes for THIS
    # body: C22 in (0, q_r·scale] keeps C/MR² in (0, 2/3]; C20 = −J2 spans the
    # matching negative range. Values inside this range but outside the GC21
    # plausible band are ALLOWED (and flagged by the live readout) — only truly
    # RD-unphysical values (C/MR² > 2/3) are refused.
    _grav_clamp = _body_gravity_constants(slot.get('cache_path'))
    _gc_bounds = {}
    if _grav_clamp is not None:
        from PlanetProfile.Inference.gravity_obs import (
            rotation_parameter as _rp, R_REF_GC21_M as _RREF,
            J2_OVER_C22 as _J2OC)
        _om, _R, _M = _grav_clamp
        _c22_hi = _rp(_om, _R, _M) * (_R / _RREF) ** 2   # C/MR²=2/3 edge
        _gc_bounds['C22'] = (0.0, float(_c22_hi))
        _gc_bounds['C20'] = (float(-_J2OC * _c22_hi), 0.0)

    def _clamp_to_bounds(key, lo, hi):
        """on_change callback: snap a widget's own session_state value into
        [lo, hi]. NOTE we deliberately do NOT pass min_value/max_value to the
        number_input — Streamlit's frontend enforces those by REJECTING an
        out-of-range typed entry and reverting to the previous value, so the
        entry never reaches session_state and cannot be snapped to the bound.
        With no min/max the typed value is accepted, this callback runs before
        the ensuing rerun, and rewriting the widget's key here makes the
        display snap to the nearest bound (per user request)."""
        v = st.session_state.get(key)
        if v is None:
            return
        try:
            st.session_state[key] = float(min(max(float(v), lo), hi))
        except (TypeError, ValueError):
            pass

    def _obs_input(name):
        c1, c2 = st.columns(2)
        default = slot.get('default_obs', {}).get(name, 0.0)
        # Scientific display + physical clamp for the tiny gravity coefficients.
        is_grav = name in _gc_bounds
        fmt = "%.4e" if is_grav else "%.4f"
        lo_hi = _gc_bounds.get(name)
        with c1:
            key = f'amort_obs_{slot_key}_{name}'
            kw = {'format': fmt, 'key': key}
            if lo_hi is not None:
                lo, hi = lo_hi
                # No min_value/max_value (see _clamp_to_bounds): they would
                # block the typed entry from reaching session_state, defeating
                # the snap. step still gives a sensible +/- increment.
                kw['step'] = (hi - lo) / 1000.0 or None
                kw['on_change'] = _clamp_to_bounds
                kw['args'] = (key, lo, hi)
                kw['help'] = (
                    f"Radau–Darwin-physical range for this body: "
                    f"[{lo:.3e}, {hi:.3e}]. Entries outside snap to the "
                    "nearest bound (C/MR² must stay ≤ 2/3). Values inside this "
                    "range but outside the GC21 plausible band are allowed.")
                default = min(max(float(default), lo), hi)
            x_obs[name] = st.number_input(
                f"{name}:", value=float(default), **kw)
        sig = float(train_sigma.get(name, float('nan')))
        obs_sigma[name] = sig
        with c2:
            st.number_input(
                "± σ (frozen at training):", value=sig,
                format=fmt, disabled=True,
                key=f'amort_sigma_{slot_key}_{name}',
                help="The flow was trained with this observation noise; "
                     "changing σ requires retraining (or MCMC mode).")

    for name in scalar_names:
        _obs_input(name)

    # Live hydrostatic-C/MR² readout: as the C22 (and C20) inputs change, show
    # the C/MR² the observed C22 would imply IF the body were hydrostatic
    # (Radau–Darwin). This is the "hydrostatic reference" — the ACTUAL C/MR²
    # comes from the inferred structure only after Generate. Works for values
    # ABOVE/BELOW the plausible band (per user request), flagging RD-unphysical
    # C22 (would need C/MR² > 2/3).
    if 'C22' in scalar_names:
        _render_live_hydrostatic_cmr2(x_obs, _grav_clamp)

    if bind_names and bind_raw:
        for name in bind_names:
            _obs_input(name)
    elif bind_names:
        # Physical entry: one complex Ae per excitation frequency.
        be = _europa_excitation_moments()
        # label order as they appear in obs_names
        labels, comps = [], {}
        for n in bind_names:
            lab, comp = n[len('Bind_'):].rsplit('_', 2)[0], n.rsplit('_', 2)[1]
            if lab not in labels:
                labels.append(lab)
            comps.setdefault(lab, [])
            if comp not in comps[lab]:
                comps[lab].append(comp)
        st.markdown("**Induction response Ae per frequency** (vector "
                    "components are derived as Ae × Be and conditioned "
                    "jointly):")
        defaults_obs = slot.get('default_obs', {})
        for lab in labels:
            ref_comp = comps[lab][0]
            be_ref = be.get(lab, {}).get(ref_comp)
            d_re = defaults_obs.get(f'Bind_{lab}_{ref_comp}_real', 0.0)
            d_im = defaults_obs.get(f'Bind_{lab}_{ref_comp}_imag', 0.0)
            ae_default = (complex(d_re, d_im) / be_ref) if be_ref else 0j
            c1, c2 = st.columns(2)
            with c1:
                ae_re = st.number_input(f"Re Ae ({lab}):",
                                        value=float(ae_default.real),
                                        format="%.4f",
                                        key=f'amort_obs_ae_{slot_key}_{lab}_re')
            with c2:
                ae_im = st.number_input(f"Im Ae ({lab}):",
                                        value=float(ae_default.imag),
                                        format="%.4f",
                                        key=f'amort_obs_ae_{slot_key}_{lab}_im')
            ae = complex(ae_re, ae_im)
            derived = []
            for comp in comps[lab]:
                be_c = be.get(lab, {}).get(comp)
                if be_c is None:
                    st.error(f"No excitation moment for {lab}/{comp} — "
                             f"cannot derive Bind components.")
                    return None
                bind = ae * be_c
                x_obs[f'Bind_{lab}_{comp}_real'] = float(bind.real)
                x_obs[f'Bind_{lab}_{comp}_imag'] = float(bind.imag)
                obs_sigma[f'Bind_{lab}_{comp}_real'] = float(
                    train_sigma.get(f'Bind_{lab}_{comp}_real', float('nan')))
                obs_sigma[f'Bind_{lab}_{comp}_imag'] = float(
                    train_sigma.get(f'Bind_{lab}_{comp}_imag', float('nan')))
                derived.append(f"{comp}: {bind.real:+.2f}{bind.imag:+.2f}i")
            st.caption(f"{lab} components (nT, σ = 1.5 each): "
                       + ";  ".join(derived))

    # --- Sampler settings ---
    st.markdown("#### ⚙️ Sampling")
    # Public serving caps: forward-model re-evaluation is the CPU cost
    # (~100 evals is tens of seconds on a shared cloud vCPU); flow sampling
    # is cheap but memory scales with draws.
    if public_mode():
        post_opts, post_default = [1000, 2000, 5000, 10000, 20000], 5000
        reeval_opts, reeval_default = [100, 250], 100
    else:
        post_opts, post_default = [1000, 2000, 5000, 10000, 20000, 50000], 10000
        reeval_opts, reeval_default = [100, 250, 500, 1000, 2000], 500
    n_post = st.select_slider("Posterior samples:",
                              options=post_opts,
                              value=post_default, key='amort_n_post')
    n_reeval = st.select_slider("Heating re-evaluation subset:",
                                options=reeval_opts,
                                value=reeval_default, key='amort_n_reeval')
    # Derived-quantity recompute subset. This is the dominant wall-time cost
    # of amortized conditioning (two forward-model passes per sample:
    # k2/CMR2/thicknesses + the true Gaussian log-likelihood). The flow draw
    # itself is instant; capping the recompute is what makes conditioning
    # interactive. Rows beyond the subset are NaN-padded — k2/CMR2 cloud
    # plots already filter NaN, and best-fit uses nanargmax. "All" recomputes
    # every sample (the artifact/validation-grade path).
    derived_opts = [500, 1000, 2000, 5000, 'All']
    n_derived_choice = st.select_slider(
        "Derived-quantity recompute subset (k₂/CMR²/likelihood):",
        options=derived_opts, value=1000, key='amort_n_derived',
        help="Caps the forward-model recompute that dominates conditioning "
             "time. The subset always includes the heating subset. Lower = "
             "faster conditioning; un-recomputed samples are omitted from "
             "the k₂/CMR² clouds. 'All' matches the validation-grade run.")
    n_derived = None if n_derived_choice == 'All' else int(n_derived_choice)
    seed = st.number_input("Seed:", value=42, step=1, key='amort_seed')

    # --- Rough estimated conditioning time (no extra compute) ------------
    # Driven by the two SUBSET sliders (derived recompute + heating), which are
    # the real cost levers: wall time ≈ overhead + per_pass·(2·derived + heating).
    # The Posterior-samples slider is intentionally NOT a big lever here — the
    # flow draw is ~instant and the recompute is capped by the derived subset —
    # so the estimate mostly moves when you change the recompute subsets. Seeded
    # with a typical-laptop figure; replaced by the observed time after your
    # first real run on this server.
    est_s, work_units, n_der_eff_est = _estimate_amort_runtime(
        slot_key, int(n_post), int(n_reeval), n_derived)
    calibrated = f'amort_perpass_{slot_key}' in st.session_state
    n_der_txt = ("all samples" if (n_derived is None or n_derived >= n_post)
                 else f"{n_der_eff_est:,} of {int(n_post):,}")
    st.info(
        f"⏱️ **Rough conditioning time: {_fmt_duration(est_s)}** "
        f"(≈ {work_units:,} forward-model passes: 2 × {n_der_txt} derived "
        f"+ {int(n_reeval):,} heating).  "
        + ("Based on your last run on this server."
           if calibrated else
           "Typical-laptop estimate; replaced by the measured time after your "
           "first run.")
        + "  Set mainly by the recompute-subset sliders, not posterior count.")

    # --- Validated-domain guard (hard refusal, per deployment conditions) ---
    if slot.get('scope_note'):
        st.caption(f"ℹ️ {slot['scope_note']}")
    domain_violations = []
    for name, (lo, hi) in (slot.get('x_obs_limits') or {}).items():
        val = x_obs.get(name)
        if val is not None and not (lo <= abs(val) <= hi):
            domain_violations.append(
                f"{name} = {val:g} outside validated domain [{lo:g}, {hi:g}]")
    # Derived-|Ae| guard (Clipper-era Bind channels): the W1 grid-walk
    # validates conditioning only over the |Ae(label)| envelope its anchors
    # exercised. |Ae| implied by the user's Bind inputs is |Bind| / |Be|
    # per component (Bind = Ae * Be); the highest-SNR component is checked.
    for g in (slot.get('derived_ae_guards') or []):
        label, comp = g['label'], g['comp']
        re_v = x_obs.get(f"Bind_{label}_{comp}_real")
        im_v = x_obs.get(f"Bind_{label}_{comp}_imag")
        if re_v is None or im_v is None:
            continue
        ae = abs(complex(re_v, im_v)) / abs(complex(*g['Be_comp']))
        lo, hi = g['ae_range']
        if not (lo <= ae <= hi):
            domain_violations.append(
                f"implied |Ae_{label}| = {ae:.3f} outside the "
                f"grid-walk-validated envelope [{lo:g}, {hi:g}]")
        warn_below = g.get('warn_support_below')
        if warn_below is not None and ae < warn_below:
            st.warning(
                f"⚠️ Your {label} Bind inputs imply |Ae_{label}| = {ae:.3f} "
                f"< {warn_below:g}, conflicting with the Galileo support cut "
                f"baked into training (|Ae_synodic| > {warn_below:g}) — the "
                f"flow cannot represent such models.")
    if domain_violations:
        st.error("🚫 **Outside validated conditioning domain** — this "
                 "artifact's gates only validate posteriors within its "
                 "domain; conditioning beyond it is refused (use MCMC mode "
                 "instead): " + "; ".join(domain_violations))
        return None

    # --- Structure cache (needed for the forward-model recompute) ---
    cache_path = slot['cache_path']
    full_cache = Path(parent_directory) / cache_path
    if full_cache.exists():
        st.success(f"✅ Structure cache found ({full_cache.stat().st_size/1024:.0f} KB)")
    else:
        st.error(f"❌ Structure cache not found: {full_cache} — rebuild it "
                 f"(see sbi_artifacts/INDEX.md) before running.")
        return None

    return {
        'artifact_path': str(artifact_path),
        'runner': runner,
        'bodyname': slot['bodyname'],
        'param_bounds': bounds,
        'x_obs': x_obs,
        'obs_sigma': obs_sigma,
        'truncate_bounds': truncate_bounds or None,
        'n_posterior_samples': int(n_post),
        'n_reeval': int(n_reeval),
        'n_derived': n_derived,
        'seed': int(seed),
        'cache_path': cache_path,
        'config_path': slot.get('config_path'),
        'validated_version_pairs': slot.get('validated_version_pairs', ()),
    }


def render_model_assumptions(exec_mode: str):
    """Model build-up + assumptions text, shared by both execution modes."""
    with st.expander("📖 Model build-up & assumptions"):
        st.markdown("""
**All current models are 1D** — spherically symmetric radial interior
structures (the "1D ·" tag on each pretrained model). Parameter counts
like 7D/8D refer to the number of *sampled parameters*, not structure
dimensionality. The roadmap targets one 1D model per mission–body pair
(Galileo/JUICE Ganymede & Callisto, Cassini Enceladus, JUICE Europa;
JUICE projections from Van Hoolst et al. 2024, SSRv 220:54). Laterally
varying 3D models, with spacecraft trajectory/measurement-geometry
visualization, are future work.

**Forward model chain** (identical physics for MCMC and amortized modes):

1. **Interior structure** — a precomputed PlanetProfile thermodynamic profile
   grid over the ice-shell basal temperature `T_b`: radial density, bulk and
   shear moduli, and phase assignments from the equation-of-state stack
   (SeaFreeze et al.), built once per body configuration and cached
   (the *structure cache*). Each posterior sample selects the nearest-`T_b`
   grid structure.
2. **Rheology** — sampled parameters modify the viscoelastic properties per
   ice/rock phase: Andrade (`α`, `ζ`, per-phase viscosities `η`) or Maxwell
   (viscosities only). Viscosities are sampled as log₁₀ quantities.
3. **Tidal response** — the TidalPy multilayer radial solver computes the
   complex Love numbers k₂ (and h₂) for the modified structure at the
   body's tidal forcing frequency; per-phase tidal heating comes from the
   same solution (volumetric heating integrated over each phase layer).
4. **Moment of inertia** — configurations that sample a core
   (`R_core`, `ρ_core`) derive the silicate density by mass conservation
   against the cached hydrosphere and integrate C/MR² analytically
   (plus a per-`T_b` discretization-offset anchor where a sidecar file
   exists). Core-blind configurations read the cached C/MR², which is then
   effectively constant.
5. **Magnetic induction** (when configured) — the cached conductivity
   profile yields the complex induction response `Ae` per excitation
   frequency.
6. **Likelihood** — a diagonal Gaussian χ² over the configured observables
   (|Im k₂| convention for the imaginary channel). No-ocean configurations
   apply a phase-stability guard that rejects samples whose ice-Ih shell
   would melt.
7. **Priors** — uniform boxes over the sampled parameters.

**Key assumptions:** spherical symmetry; linear viscoelastic response at
the tidal frequency; the chosen rheological law; hydrostatic interior
structure from the PlanetProfile profile; fixed ocean/rock composition per
configuration; independent Gaussian observational errors.
""")
        if exec_mode == 'amortized':
            st.markdown("""
**Amortized specifics:** a neural posterior estimator (normalizing flow,
NPE) was trained once on prior-predictive simulations of this forward
model, with observation noise injected to match the Gaussian likelihood.
Priors and observation σ are therefore **frozen at training time**;
conditioning on new observable values is instant. The artifact was
validated against the MCMC posterior via simulation-based calibration
(SBC), a full cross-check, and ground-truth anchor comparisons — see
`sbi_artifacts/INDEX.md` for the validated conditioning domain. The
plotted k₂ / C-MR² / heating values are recomputed with the *same forward
model* used by MCMC, not by the flow.
""")
        else:
            st.markdown("""
**MCMC specifics:** the posterior is sampled with pocoMC (preconditioned
Monte Carlo) against the full likelihood, so priors, observables, and
uncertainties are all freely configurable. Convergence is tracked via
effective sample size, acceptance rate, and R-hat.
""")


def _amortized_spec_fingerprint(spec):
    """Stable fingerprint of everything that changes the posterior, used to
    flag stale results when controls move after a run."""
    import json as _json
    return _json.dumps({
        'artifact': spec['artifact_path'],
        'x_obs': spec['x_obs'],
        'truncate': spec['truncate_bounds'],
        'n_post': spec['n_posterior_samples'],
        'n_reeval': spec['n_reeval'],
        'n_derived': spec.get('n_derived'),
        'seed': spec['seed'],
    }, sort_keys=True)


def render_amortized_run_button(spec, InferenceConfig):
    """Run button for amortized mode: instant conditioning + recompute."""
    st.markdown("---")
    if spec is None:
        st.button("⚡ Generate Posterior", disabled=True)
        return
    if not st.button("⚡ Generate Posterior", type="primary",
                     help="Conditions the pretrained flow (instant) and "
                          "recomputes k2/C-MR2/heating for the plots (~seconds)."):
        return

    try:
        import json as _json
        from PlanetProfile.Inference.sbi_runner import SBIRunner

        # Build the runner config from the TRAINING config JSON, not a
        # minimal reconstruction: derivation blocks (derived_params
        # mass-conservation) change how CMR2 is recomputed for the plots.
        # A bare config here silently fell back to the core-blind v1 CMR2
        # path and put Europa's posterior CMR2 ~3.6 sigma off. Only
        # observable VALUES/sigmas + sampler settings are GUI-overridable.
        config_path = spec.get('config_path')
        if config_path:
            cfg_dict = _json.loads(
                (Path(parent_directory) / config_path).read_text())
            cfg_dict['mode'] = 'sbi'
            cfg_dict['observables'] = {
                n: [spec['x_obs'][n], spec['obs_sigma'][n]]
                for n in spec['x_obs']}
            cfg_dict.setdefault('sampler_settings', {})['n_reeval'] = \
                spec['n_reeval']
            cfg_dict['structure_cache_path'] = spec['cache_path']
            cfg_dict['random_state'] = spec['seed']
            config = InferenceConfig.from_dict(cfg_dict)
        else:
            config = InferenceConfig(
                mode='sbi',
                bodyname=spec['bodyname'],
                param_space={n: {'prior_type': 'uniform', 'bounds': [lo, hi]}
                             for n, (lo, hi) in spec['param_bounds'].items()},
                observables={n: [spec['x_obs'][n], spec['obs_sigma'][n]]
                             for n in spec['x_obs']},
                sampler_settings={'n_reeval': spec['n_reeval']},
                structure_cache_path=spec['cache_path'],
                random_state=spec['seed'],
            )
        runner = SBIRunner(config)

        progress_placeholder = st.empty()

        def progress_callback(d):
            progress_placeholder.progress(
                d['iteration'] / d['n_total'],
                text=f"Stage {d['iteration']} / {d['n_total']}")

        with st.spinner("Conditioning pretrained flow + recomputing observables..."):
            result = runner.infer_from_artifact(
                spec['artifact_path'],
                x_obs=spec['x_obs'],
                n_posterior_samples=spec['n_posterior_samples'],
                seed=spec['seed'],
                n_reeval=spec['n_reeval'],
                n_derived=spec.get('n_derived'),
                truncate_bounds=spec['truncate_bounds'],
                progress_callback=progress_callback,
                validated_version_pairs=spec.get('validated_version_pairs', ()),
            )

        st.session_state.inference_results = result
        st.session_state.amort_last_run_fp = _amortized_spec_fingerprint(spec)
        # Self-calibrate the runtime estimator from this run's actual wall time.
        _slot_key = Path(spec['artifact_path']).name.replace('.', '_').replace(' ', '_')
        _record_amort_runtime(_slot_key, result)
        kept = result.metadata.get('kept_fraction')
        msg = "✅ Amortized posterior ready."
        if kept is not None:
            msg += f" Truncation kept {kept:.1%} of draws."
        st.success(msg)

    except Exception as e:
        st.session_state.inference_error = str(e)
        st.session_state.inference_error_traceback = traceback.format_exc()
        st.error(f"❌ **Amortized inference failed:** {e}")
        with st.expander("🔍 Full traceback"):
            st.code(st.session_state.inference_error_traceback)


# --- Crisp corner-plot rendering (vector inline + vector download) ---------
#
# st.pyplot() rasterizes a matplotlib figure to a screen-DPI PNG, then the
# browser downscales that big bitmap into the ~700 px column — so a 32" corner
# canvas with 10 pt titles ends up blurry and unreadable. Fix, three parts:
#   (A) render the figure INLINE as SVG (vector text, crisp at any zoom),
#   (B) offer PDF + PNG DOWNLOAD buttons (publication-quality vector copy),
#   (C) larger fonts, capped figure size, and corner's built-in Gaussian
#       smoothing on filled contours, RASTERIZED inside the vector frame.
#
# Why smooth + rasterize: a corner over tens of thousands of posterior samples,
# left as vector scatter/hist2d, produces a ~10 MB SVG that chokes the browser.
# corner's own Gaussian smoothing (smooth/smooth1d) gives clean filled contour
# bands near-instantly; rasterizing only those 2D collections (text/axes stay
# vector) drops the SVG to ~0.7 MB / PDF ~0.2 MB with no loss of legibility.
# (A scipy gaussian_kde density looks marginally smoother but adds ~10 s per
# plot across ~78 panels — not worth it; user chose Gaussian smoothing.)
_CORNER_MAX_FIGSIZE = 26.0   # inches; cap so fonts stay legible on many-param models


def _render_corner_figure(samples, labels, *, seed=0):
    """Build a crisp, Gaussian-smoothed corner figure (2D density rasterized).

    Returns a matplotlib Figure. Titles/labels/axes stay vector; only the heavy
    2D density collections are rasterized, so SVG/PDF export stays small.
    """
    import matplotlib.pyplot as plt
    import corner

    n_dim = samples.shape[1]
    fig_size = min(_CORNER_MAX_FIGSIZE, max(10.0, 2.5 * n_dim))
    fig = plt.figure(figsize=(fig_size, fig_size))

    corner.corner(
        samples,
        labels=labels,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_fmt='.2f',
        title_kwargs={'fontsize': 15},
        label_kwargs={'fontsize': 16},
        color='steelblue',
        smooth=1.0,             # Gaussian smoothing of the 2D histograms
        smooth1d=1.0,           # ... and the 1D marginals
        plot_datapoints=False,  # drop the per-sample scatter (the SVG-bloat source)
        fill_contours=True,
        hist_kwargs={'linewidth': 1.5},
        fig=fig,
    )

    for ax in fig.get_axes():
        ax.tick_params(labelsize=10)
    # rasterize only the heavy 2D density; text/lines stay vector
    _rasterize_heavy_artists(fig)

    return fig


from Utilities.crisp_figs import (
    rasterize_heavy_artists as _rasterize_heavy_artists,
    display_vector_fig as _crisp_display)


def _result_token():
    """Cache token for figure exports: stable while the same result object
    sits in session state; changes on any new run / loaded pickle."""
    res = st.session_state.get('inference_results')
    return None if res is None else (id(res), len(getattr(res, 'samples', ())))


def _display_vector_fig(fig, *, key, download_label='plot'):
    """Crisp SVG + cached PDF/PNG downloads (see Utilities/crisp_figs).
    Exports are cached per (key, result token): reruns skip all three
    savefig calls — the dominant per-rerun cost on the results page."""
    _crisp_display(fig, key=key, download_label=download_label,
                   token=_result_token())


def _display_corner(fig, *, key):
    """Show a corner figure inline as crisp SVG + PDF/PNG downloads."""
    _display_vector_fig(fig, key=key, download_label='corner plot')


def render_results():
    """Render inference results (corner plots, traces, etc.)."""

    if st.session_state.inference_results is None:
        st.info("📊 **Results will appear here after running inference.**")
        st.markdown("""
        Run inference above to see:
        - Corner plots (posterior marginals + covariances)
        - k₂ posterior scatter with 1σ/2σ observation ellipses
        - C/MR² moment-of-inertia posterior histogram
        - Per-phase heating distributions
        - Export to PKL
        """)
        return

    # Results available
    result = st.session_state.inference_results

    st.markdown("### 📊 Inference Results")

    def _phase_heating_total(h):
        """Total heating over PHYSICAL phases only. compute_heating returns
        per-phase keys ('Sil', 'Ih', 'VI', ...) PLUS *_W aggregate keys
        (Silicate_W, HP_Ice_W, Ice_Ih_W) that duplicate them — summing all
        values double-counts and halves every fraction (user-visible
        2026-07-13: k2 scatter showed silicate fraction 0.2-0.7 while the
        stacked-bar panel correctly showed ~1)."""
        return sum(v for k, v in h.items() if not k.endswith('_W'))

    # Convergence metrics (mode-aware: MCMC diagnostics don't apply to an
    # amortized flow — its calibration is established by the SBC/crosscheck
    # validation gates at training time, not per-run statistics).
    is_amortized = (result.metadata or {}).get('sampler') == 'sbi'
    with st.expander("📈 Sampling Diagnostics", expanded=True):
        metrics = result.convergence_metrics
        if is_amortized:
            col1, col2, col3 = st.columns(3)
            col1.metric("Posterior draws", f"{metrics['n_samples']}")
            kept = (result.metadata or {}).get('kept_fraction')
            col2.metric("Truncation kept",
                        f"{kept:.1%}" if kept is not None else "no truncation")
            col3.metric("Sampler", "Amortized NPE")
            st.caption(
                f"**Runtime:** {result.metadata['elapsed_time_s']/60:.1f} min · "
                f"artifact: `{Path(str(result.metadata.get('artifact_path', ''))).name}` "
                f"(trained on {result.metadata.get('n_train_effective', '?'):,} sims). "
                f"Acceptance rate and R-hat are MCMC convergence diagnostics and "
                f"do not apply here; this artifact's calibration was validated "
                f"via the SBC / cross-check / anchor gates (see "
                f"`sbi_artifacts/INDEX.md`)."
            )
        else:
            col1, col2, col3 = st.columns(3)
            col1.metric("ESS", f"{metrics['ess']:.0f}")
            acc = metrics['acceptance_rate']
            col2.metric("Acceptance Rate", f"{acc:.1%}" if acc is not None else "N/A")
            rhat = metrics.get('r_hat')
            col3.metric("R-hat", f"{rhat:.3f}" if rhat is not None else "N/A")
            st.caption(f"**Samples:** {metrics['n_samples']} | **Runtime:** {result.metadata['elapsed_time_s']/60:.1f} min")

    # Observables the posterior was actually conditioned on — makes any
    # widget/config mixup immediately visible (amortized runs record the
    # exact x_obs; older results fall back to the config's central values).
    with st.expander("🔭 Observables used in this run", expanded=False):
        used = (result.metadata or {}).get('x_obs')
        if used:
            st.table({'observable': list(used.keys()),
                      'conditioned value': [f"{v:.4g}" for v in used.values()]})
            st.caption("Values the pretrained flow was conditioned on for "
                       "this run (recorded at run time).")
        else:
            obs = getattr(result.config, 'observables', {}) or {}
            if obs:
                st.table({'observable': list(obs.keys()),
                          'value': [f"{v[0]:.4g}" for v in obs.values()],
                          'σ': [f"{v[1]:.4g}" for v in obs.values()]})
                st.caption("From this run's configuration (run predates "
                           "x_obs recording).")
            else:
                st.info("No observable record on this result.")

    # Summary statistics
    with st.expander("📋 Parameter Summary", expanded=True):
        st.markdown("**Posterior Statistics:**")

        # Build summary table
        summary_data = []
        for i, param_name in enumerate(result.param_names):
            samples_param = result.samples[:, i]
            summary_data.append({
                'Parameter': result.param_labels[i],
                'Mean': f"{np.mean(samples_param):.4f}",
                'Median': f"{np.median(samples_param):.4f}",
                'Std': f"{np.std(samples_param):.4f}",
                '16th %ile': f"{np.percentile(samples_param, 16):.4f}",
                '84th %ile': f"{np.percentile(samples_param, 84):.4f}",
            })

        st.table(summary_data)

    # Corner plot
    with st.expander("🔺 Corner Plot", expanded=True):
        try:
            import corner
            import matplotlib.pyplot as plt

            # Augment posterior samples with derived structure quantities when available.
            # D_ocean and D_iceIh are functions of Tb_K via the grid cache, so their
            # off-diagonal panels reveal which ocean/ice thickness range the posterior
            # prefers, and the diagonal panels give the marginal distributions.
            corner_samples = result.samples
            corner_labels = list(result.param_labels)

            D_ocean = getattr(result, 'D_ocean_results', None)
            D_iceIh = getattr(result, 'D_iceIh_results', None)

            if D_ocean is not None and np.any(np.isfinite(D_ocean)):
                corner_samples = np.column_stack([corner_samples, D_ocean])
                corner_labels.append(r'$D_{\rm ocean}$ (km)')
            if D_iceIh is not None and np.any(np.isfinite(D_iceIh)):
                corner_samples = np.column_stack([corner_samples, D_iceIh])
                corner_labels.append(r'$D_{\rm IceIh}$ (km)')

            # Filter out columns with no dynamic range (zero variance)
            # corner.corner fails if any column is constant.
            valid_cols = []
            for i in range(corner_samples.shape[1]):
                col_data = corner_samples[:, i]
                col_data = col_data[np.isfinite(col_data)]
                if len(col_data) > 0 and np.ptp(col_data) > 1e-12:
                    valid_cols.append(i)
            
            if not valid_cols:
                st.warning("No parameters with dynamic range found for corner plot.")
            else:
                corner_samples_plot = corner_samples[:, valid_cols]
                corner_labels_plot = [corner_labels[i] for i in valid_cols]
                # Drop rows with any non-finite value: derived columns
                # (D_ocean/D_iceIh) are NaN-padded outside the recompute
                # subset in amortized mode (n_derived), and corner.corner
                # errors on NaN. Parameter columns are always finite, so
                # this only prunes un-recomputed derived rows.
                row_ok = np.all(np.isfinite(corner_samples_plot), axis=1)
                if not np.all(row_ok):
                    corner_samples_plot = corner_samples_plot[row_ok]

                # builder path: corner() is the most expensive build on
                # the page (~seconds) — skip BOTH build and exports on
                # cache hit (token = current result).
                _crisp_display(
                    builder=lambda: _render_corner_figure(
                        corner_samples_plot, corner_labels_plot),
                    key='corner', download_label='corner plot',
                    token=_result_token())
        except ImportError:
            st.info("Install the `corner` library to view corner plots: `pip install corner`")
        except Exception as e:
            st.warning(f"Corner plot unavailable: {e}")

    # Interior layer-thickness distributions (posterior structural summary)
    with st.expander("🧅 Interior Layer Thicknesses", expanded=True):
        try:
            import matplotlib.pyplot as plt

            layers = []  # (label, values_km, color)
            D_ice = getattr(result, 'D_iceIh_results', None)
            D_oc = getattr(result, 'D_ocean_results', None)
            D_hs = getattr(result, 'D_hsphere_results', None)
            R_body = (result.metadata or {}).get('R_body_km')
            R_core = None
            if 'R_core_km' in result.param_names:
                R_core = result.samples[
                    :, result.param_names.index('R_core_km')].astype(float)

            if D_ice is not None and np.any(np.isfinite(D_ice)):
                layers.append(('Ice Ih shell thickness',
                               np.asarray(D_ice, float), '#AEE1F8'))
            if D_oc is not None and np.any(np.isfinite(D_oc)):
                layers.append(('Ocean thickness',
                               np.asarray(D_oc, float), '#1E90FF'))
            # Rock mantle thickness = R_body - hydrosphere - core radius.
            # Uses the packaged total hydrosphere D_hsphere (includes HP
            # ices) so the derivation is exact for HP-ice bodies too.
            if (R_body is not None and R_core is not None
                    and D_hs is not None and np.any(np.isfinite(D_hs))):
                layers.append(('Rock mantle thickness',
                               R_body - np.asarray(D_hs, float) - R_core,
                               '#C8A96E'))
            if R_core is not None:
                layers.append(('Core radius', R_core, '#8B5A2B'))

            if not layers:
                st.info(
                    "No layer-thickness data on this result. Thickness "
                    "arrays are packaged by runs made after 2026-07; older "
                    "saved results lack them. Mantle thickness additionally "
                    "needs a sampled core (R_core_km) and a result carrying "
                    "the hydrosphere thickness + body radius.")
            else:
                n_cols = 2
                n_rows = (len(layers) + n_cols - 1) // n_cols
                fig, axes = plt.subplots(n_rows, n_cols,
                                         figsize=(10, 3.2 * n_rows))
                axes = np.atleast_1d(axes).ravel()
                for ax, (label, vals, color) in zip(axes, layers):
                    v = vals[np.isfinite(vals)]
                    ax.hist(v, bins=40, color=color, edgecolor='none',
                            alpha=0.85)
                    lo, med, hi = np.percentile(v, [16, 50, 84])
                    ax.axvline(med, color='k', linewidth=1.2)
                    ax.axvline(lo, color='k', linewidth=0.8, linestyle=':')
                    ax.axvline(hi, color='k', linewidth=0.8, linestyle=':')
                    ax.set_title(f"{label}\n"
                                 f"{med:.1f} (+{hi - med:.1f} / "
                                 f"-{med - lo:.1f}) km", fontsize=10)
                    ax.set_xlabel('km')
                    ax.set_ylabel('samples')
                for ax in axes[len(layers):]:
                    ax.axis('off')
                fig.tight_layout()
                _display_vector_fig(fig, key='layer_thicknesses',
                                    download_label='layer thicknesses')
                plt.close(fig)
                bits = ["Median with 16–84% interval (dotted)."]
                if any(lbl.startswith('Rock mantle') for lbl, _, _ in layers):
                    bits.append(
                        "Mantle thickness is derived per sample as "
                        "R_body − hydrosphere thickness − core radius "
                        "(hydrosphere includes any HP-ice layers).")
                if R_core is None:
                    bits.append(
                        "No sampled core in this configuration, so core "
                        "radius and mantle thickness are not shown.")
                st.caption(" ".join(bits))
        except Exception as e:
            st.warning(f"Layer-thickness panel unavailable: {e}")

    # k2 scatter with error ellipse
    with st.expander("📡 k₂ Posterior Scatter", expanded=True):
        try:
            import matplotlib.pyplot as plt
            from matplotlib.patches import Ellipse

            # Observed-ellipse values come ONLY from this run's config —
            # never hardcoded fallbacks (the old defaults were Titan's k2
            # and drew a Titan ellipse on every non-Titan run whose config
            # keyed the imaginary channel 'Im_k2' instead of 'abs_Im_k2').
            def _obs_channel(*aliases):
                for a in aliases:
                    if a in result.config.observables:
                        v = result.config.observables[a]
                        return float(v[0]), float(v[1])
                return None, None
            Re_obs, Re_err = _obs_channel('Re_k2')
            Im_obs, Im_err = _obs_channel('Im_k2', 'abs_Im_k2')

            heating_results = result.heating_results or []
            eval_idx = result.metadata.get('heating_indices')

            # Slice k2_results to the heating-evaluated subset so lengths match
            if heating_results and eval_idx is not None and len(heating_results) == len(eval_idx):
                Re_arr = result.k2_results[eval_idx, 0]
                Im_arr = np.abs(result.k2_results[eval_idx, 1])
                eval_samples = result.samples[eval_idx]
            else:
                Re_arr = result.k2_results[:, 0]
                Im_arr = np.abs(result.k2_results[:, 1])
                eval_samples = result.samples

            f_sil = []
            for h in heating_results:
                total = _phase_heating_total(h) if h else 0.0
                sil = h.get('Sil', 0.0) + h.get('Fe', 0.0)
                f_sil.append(sil / total if total > 0 else 0.0)
            f_sil = np.array(f_sil) if f_sil else None

            # Scale point size by ocean thickness (Tb_K proxy) when available
            pt_size = 8
            if 'Tb_K' in result.param_names:
                tb_idx = result.param_names.index('Tb_K')
                tb_vals = eval_samples[:, tb_idx]
                tb_norm = (tb_vals - tb_vals.min()) / (np.ptp(tb_vals) + 1e-10)
                pt_size = 4 + 20 * tb_norm  # 4–24 px, larger = warmer = more ocean

            plt.rcParams['text.usetex'] = False
            fig, ax = plt.subplots(figsize=(8, 6))
            if f_sil is not None and len(f_sil) == len(Re_arr):
                sc = ax.scatter(Re_arr, Im_arr, c=f_sil, cmap='viridis',
                                s=pt_size, alpha=0.6, vmin=0, vmax=1)
                plt.colorbar(sc, ax=ax, label='Silicate heating fraction')
            else:
                ax.scatter(Re_arr, Im_arr, s=pt_size, alpha=0.6, color='steelblue')

            have_ellipse = Re_obs is not None and Im_obs is not None
            if have_ellipse:
                ax.add_patch(Ellipse((Re_obs, abs(Im_obs)), 2*Re_err, 2*Im_err,
                                     fill=False, color='red', linewidth=2,
                                     linestyle='--', label=r'1$\sigma$'))
                ax.add_patch(Ellipse((Re_obs, abs(Im_obs)), 4*Re_err, 4*Im_err,
                                     fill=False, color='red', linewidth=1,
                                     linestyle=':', label=r'2$\sigma$'))
            ax.set_xlabel(r'Re$(k_2)$')
            ax.set_ylabel(r'$|$Im$(k_2)|$')
            ax.set_title(r'$k_2$ Posterior')
            if have_ellipse:
                ax.legend()
            _rasterize_heavy_artists(fig)
            _display_vector_fig(fig, key='k2_scatter',
                                download_label='k₂ posterior')
            plt.close(fig)
            caption_bits = []
            if have_ellipse:
                caption_bits.append(
                    "Red ellipses: 1σ / 2σ observational constraint "
                    "(this run's k₂ conditioning values).")
            else:
                caption_bits.append(
                    "No k₂ observational constraint in this run's config — "
                    "ellipse omitted.")
            if isinstance(pt_size, np.ndarray):
                caption_bits.append(
                    "**Point size** scales with the sampled basal temperature "
                    "T_b (larger = warmer ice shell base; the structural "
                    "dimension of the model).")
            if f_sil is not None and len(f_sil) == len(Re_arr):
                caption_bits.append(
                    "**Color** is the fraction of tidal heating dissipated "
                    "in the silicate interior (viridis: dark = ice-dominated, "
                    "yellow = rock-dominated dissipation).")
            st.caption(" ".join(caption_bits))

            # Component panes: marginal posteriors of both k₂ components,
            # as (Re, Im) or equivalently (amplitude, phase lag).
            pane_mode = st.radio(
                'k₂ component view', ['Re / Im', 'Amplitude / Phase'],
                horizontal=True, key='k2_pane_mode')

            def _build_k2_panes():
                figp, (axl, axr) = plt.subplots(1, 2, figsize=(10, 3.8))
                if pane_mode == 'Re / Im':
                    panes = [
                        (axl, Re_arr, r'Re$(k_2)$', Re_obs, Re_err),
                        (axr, Im_arr, r'$|$Im$(k_2)|$',
                         abs(Im_obs) if Im_obs is not None else None, Im_err),
                    ]
                else:
                    amp = np.hypot(Re_arr, Im_arr)
                    phase = np.degrees(np.arctan2(Im_arr, Re_arr))
                    if Re_obs is not None and Im_obs is not None:
                        amp_obs = float(np.hypot(Re_obs, Im_obs))
                        ph_obs = float(np.degrees(
                            np.arctan2(abs(Im_obs), Re_obs)))
                        # first-order error propagation from (σ_Re, σ_Im)
                        amp_err = float(np.hypot(Re_obs * Re_err,
                                                 Im_obs * Im_err) / amp_obs
                                        ) if amp_obs > 0 else None
                        ph_err = float(np.degrees(np.hypot(
                            Im_obs * Re_err, Re_obs * Im_err)
                            / amp_obs ** 2)) if amp_obs > 0 else None
                    else:
                        amp_obs = ph_obs = amp_err = ph_err = None
                    panes = [
                        (axl, amp, r'$|k_2|$', amp_obs, amp_err),
                        (axr, phase, r'phase lag $\arg(k_2)$ (deg)',
                         ph_obs, ph_err),
                    ]
                for axi, arr, lab, obs, err in panes:
                    axi.hist(arr, bins=40, color='steelblue', alpha=0.8)
                    if obs is not None:
                        axi.axvline(obs, color='red', linestyle='--',
                                    linewidth=1.5, label='observed')
                        if err:
                            axi.axvspan(obs - err, obs + err, color='red',
                                        alpha=0.12, label=r'1$\sigma$')
                        axi.legend(fontsize=8)
                    axi.set_xlabel(lab)
                axl.set_ylabel('posterior samples')
                figp.suptitle(r'$k_2$ component posteriors', y=1.02)
                figp.tight_layout()
                return figp

            _crisp_display(builder=_build_k2_panes, key='k2_panes',
                           download_label='k₂ components',
                           token=(_result_token(), pane_mode))
        except Exception as e:
            st.warning(f"k₂ scatter unavailable: {e}")

    # Posterior induction response (complex plane) — only when the run
    # carried induction (metadata packaged by both runners since 2026-07).
    _ind_ae = (result.metadata or {}).get('induction_Ae')
    if _ind_ae:
        with st.expander("🧲 Induction Response (Ae) Posterior", expanded=True):
            try:
                import matplotlib.pyplot as plt

                _be = (result.metadata or {}).get('Be_nT') or {}
                unit_opts = ['Ae (dimensionless)']
                if _be:
                    unit_opts.append('Ae·|B| (nT, per component)')
                unit = st.radio("Units", unit_opts, horizontal=True,
                                key='ind_ae_units')
                if not _be:
                    st.caption("Ae·|B| (nT) view needs excitation moments, "
                               "which this run's config does not carry "
                               "(induction enters as a support bound only).")

                labels = list(_ind_ae.keys())
                if unit.startswith('Ae ('):
                    # One complex-plane panel per excitation label.
                    fig, axes = plt.subplots(1, len(labels),
                                             figsize=(5 * len(labels), 4.6),
                                             squeeze=False)
                    for ax, lab in zip(axes[0], labels):
                        re = np.asarray(_ind_ae[lab]['re'], float)
                        im = np.asarray(_ind_ae[lab]['im'], float)
                        ax.scatter(re, im, s=6, alpha=0.35, color='steelblue')
                        th = np.linspace(0, 2 * np.pi, 200)
                        ax.plot(np.cos(th), np.sin(th), color='0.6',
                                linewidth=0.8, linestyle=':')
                        ax.set_xlabel(r'Re$(A_e)$')
                        ax.set_ylabel(r'Im$(A_e)$')
                        ax.set_title(lab)
                        ax.set_aspect('equal', adjustable='datalim')
                    fig.suptitle('Posterior induction response '
                                 '(dotted: $|A_e| = 1$)', fontsize=10)
                    fig.tight_layout()
                    _rasterize_heavy_artists(fig)
                    _display_vector_fig(fig, key='ae_complex_plane',
                                        download_label='Aₑ complex plane')
                    plt.close(fig)
                else:
                    # Per-component induced surface field Bind = Ae * Be
                    # (plain complex product, nT) — one row per label.
                    obs = getattr(result.config, 'observables', {}) or {}
                    for lab in labels:
                        comps = _be.get(lab)
                        if not comps:
                            st.caption(f"{lab}: no excitation moments — skipped.")
                            continue
                        ae = (np.asarray(_ind_ae[lab]['re'], float)
                              + 1j * np.asarray(_ind_ae[lab]['im'], float))
                        cnames = list(comps.keys())
                        fig, axes = plt.subplots(1, len(cnames),
                                                 figsize=(4.4 * len(cnames), 4.2),
                                                 squeeze=False)
                        for ax, comp in zip(axes[0], cnames):
                            be_c = complex(*comps[comp])
                            bind = ae * be_c
                            ax.scatter(bind.real, bind.imag, s=6, alpha=0.35,
                                       color='steelblue')
                            # Observed conditioning value + 1.5 nT circle
                            ore = obs.get(f'Bind_{lab}_{comp}_real')
                            oim = obs.get(f'Bind_{lab}_{comp}_imag')
                            if ore is not None and oim is not None:
                                ax.plot(ore[0], oim[0], 'r*', markersize=10)
                                th = np.linspace(0, 2 * np.pi, 100)
                                ax.plot(ore[0] + ore[1] * np.cos(th),
                                        oim[0] + oim[1] * np.sin(th),
                                        'r--', linewidth=0.9)
                            ax.set_xlabel(fr'Re$(B_{{ind,{comp}}})$ (nT)')
                            ax.set_ylabel(fr'Im$(B_{{ind,{comp}}})$ (nT)')
                            ax.set_title(f'{lab} · {comp}')
                            ax.set_aspect('equal', adjustable='datalim')
                        fig.suptitle('Induced dipole coefficient as equivalent '
                                     'surface field (red: conditioned value '
                                     '± 1σ)', fontsize=10)
                        fig.tight_layout()
                        _rasterize_heavy_artists(fig)
                        _display_vector_fig(
                            fig, key=f'bind_{lab}',
                            download_label=f'Bind {lab}')
                        plt.close(fig)

                # k2 on the complex plane, organized by the discrete model
                # family: posterior samples quantize to the Ae/Tb grid
                # nodes ("models"), so per-node mean k2 traces a curve —
                # points connected by line segments, marker size scaling
                # with ocean thickness and color with salinity (constant-
                # color for fixed-salinity artifacts until the salinity-
                # sampled campaign lands). Replaces the earlier 3D view,
                # which was trivial: Ae quantizes to a few nodes spanning
                # a linear k2 range.
                if st.checkbox("Complex-plane signals by model (connected)",
                               key='ind_k2_plane', value=True):
                    # Every dimensionless complex signal of one interior
                    # model — k2, h2, and Ae at each excitation frequency —
                    # on a single Re-Im plane, connected per model node.
                    lab0 = labels[0]
                    ae0 = (np.asarray(_ind_ae[lab0]['re'], float)
                           + 1j * np.asarray(_ind_ae[lab0]['im'], float))
                    k2 = np.asarray(result.k2_results, float)
                    h2 = getattr(result, 'h2_results', None)
                    h2 = np.asarray(h2, float) if h2 is not None else None
                    n = min(len(ae0), len(k2))
                    tb = (result.samples[:n,
                          result.param_names.index('Tb_K')]
                          if 'Tb_K' in result.param_names else None)
                    d_oc = getattr(result, 'D_ocean_results', None)
                    d_oc = (np.asarray(d_oc, float)[:n]
                            if d_oc is not None else None)
                    sal_col = next((p for p in result.param_names
                                    if 'wOcean' in p), None)
                    sal = (result.samples[:n,
                           result.param_names.index(sal_col)]
                           if sal_col else None)

                    # Per-sample signal arrays shared by both branches:
                    # k₂, h₂, then Ae at each excitation frequency.
                    sig_list = [('k₂', k2[:n, 0] + 1j * k2[:n, 1])]
                    if h2 is not None and len(h2) >= n:
                        sig_list.append(('h₂', h2[:n, 0] + 1j * h2[:n, 1]))
                    for lab in labels:
                        sig_list.append((f'Ae {lab}',
                                         np.asarray(_ind_ae[lab]['re'],
                                                    float)[:n]
                                         + 1j * np.asarray(_ind_ae[lab]['im'],
                                                           float)[:n]))
                    sig_names = [nm for nm, _ in sig_list]
                    markers = ['o', 's', '^', 'D', 'v', 'P'][:len(sig_names)]

                    if sal is not None:
                        # Salinity-sampled artifact (v3+): one faint path
                        # per posterior SAMPLE — the full interior model
                        # (rheology + Tb + salinity) sets every signal, so
                        # the cloud shows the joint spread and the Tb–w
                        # degeneracy as a color gradient along the Ae arcs.
                        from matplotlib.collections import LineCollection
                        n_show = min(n, 1200)
                        idx = (np.random.default_rng(0).choice(
                                   n, n_show, replace=False)
                               if n > n_show else np.arange(n))
                        S = np.array([arr for _, arr in sig_list])
                        cmap = plt.cm.viridis
                        cnorm = plt.Normalize(np.nanmin(sal), np.nanmax(sal)
                                              or 1.0)
                        if d_oc is not None:
                            dspan = np.ptp(d_oc[np.isfinite(d_oc)])
                            dspan = dspan if dspan > 0 else 1.0
                            szs = 4 + 26 * (d_oc[idx]
                                            - np.nanmin(d_oc)) / dspan
                        else:
                            szs = np.full(n_show, 8.0)

                        figk, axk = plt.subplots(figsize=(7.2, 5.6))
                        paths = np.stack([S[:, idx].real, S[:, idx].imag],
                                         axis=-1).transpose(1, 0, 2)
                        axk.add_collection(LineCollection(
                            paths,
                            colors=cmap(cnorm(sal[idx])),
                            linewidths=0.5, alpha=0.12, zorder=1))
                        for (nm, arr), mk in zip(sig_list, markers):
                            axk.scatter(arr[idx].real, arr[idx].imag,
                                        s=szs, marker=mk, c=sal[idx],
                                        cmap=cmap, norm=cnorm, alpha=0.45,
                                        linewidths=0, zorder=2)
                        axk.autoscale_view()
                        for nm, mk in zip(sig_names, markers):
                            axk.scatter([], [], marker=mk, color='0.4',
                                        edgecolor='k', linewidth=0.4, s=60,
                                        label=nm)
                        axk.legend(fontsize=8, loc='best',
                                   title='signal', title_fontsize=8)
                        figk.colorbar(plt.cm.ScalarMappable(
                            norm=cnorm, cmap=cmap), ax=axk, label=sal_col)
                        axk.set_xlabel('Re(signal)')
                        axk.set_ylabel('Im(signal)')
                        axk.set_title('Complex-plane response signals per '
                                      'posterior sample')
                        _rasterize_heavy_artists(figk)
                        _display_vector_fig(figk, key='signals_salinity',
                                            download_label='complex-plane signals')
                        plt.close(figk)
                        st.caption(
                            f"One faint path per posterior sample "
                            f"({n_show} of {n} shown) linking its k₂, h₂, "
                            "and induction response Ae at each excitation "
                            "frequency. Marker shape = signal; marker "
                            f"size = ocean thickness; color = {sal_col}. "
                            "Ae varies with both T_b and ocean salinity "
                            "through the 2D structure grid; the color "
                            "gradient along the Ae arcs traces the "
                            "T_b–salinity degeneracy ridge.")
                    else:
                        # Fixed-salinity artifacts (v1/v2): samples
                        # quantize to the Ae/Tb grid nodes, so per-sample
                        # paths collapse onto a few curves — group by node
                        # (mean signals per node), ordered by Tb.
                        _, inv = np.unique(np.round(ae0[:n], 12),
                                           return_inverse=True)
                        node_rows = []
                        for g in range(inv.max() + 1):
                            m = inv == g
                            sigs = [(nm, complex(np.mean(arr[m])))
                                    for nm, arr in sig_list]
                            node_rows.append({
                                'tb': (float(np.mean(tb[m]))
                                       if tb is not None else g),
                                'd': (float(np.mean(d_oc[m]))
                                      if d_oc is not None else 20.0),
                                'sigs': sigs,
                            })
                        node_rows.sort(key=lambda r: r['tb'])

                        dvals = np.array([r['d'] for r in node_rows])
                        span = np.ptp(dvals) if np.ptp(dvals) > 0 else 1.0
                        sizes = 25 + 175 * (dvals - dvals.min()) / span
                        cmap_tb = plt.cm.plasma
                        tvals = np.array([r['tb'] for r in node_rows])
                        tspan = np.ptp(tvals) if np.ptp(tvals) > 0 else 1.0
                        colors_n = [cmap_tb(0.15 + 0.7 * (v - tvals.min())
                                            / tspan) for v in tvals]

                        figk, axk = plt.subplots(figsize=(7.2, 5.6))
                        # faint sample clouds for the rheology-spread signals
                        axk.scatter(k2[:n, 0], k2[:n, 1], s=3, alpha=0.10,
                                    color='0.6', zorder=1)
                        if h2 is not None and len(h2) >= n:
                            axk.scatter(h2[:n, 0], h2[:n, 1], s=3,
                                        alpha=0.10, color='0.75', zorder=1)
                        for row, col, sz in zip(node_rows, colors_n, sizes):
                            pts = np.array([[c.real, c.imag]
                                            for _, c in row['sigs']])
                            axk.plot(pts[:, 0], pts[:, 1], '-', color=col,
                                     linewidth=0.9, alpha=0.75, zorder=2)
                            for (nm, c), mk in zip(row['sigs'], markers):
                                axk.scatter(c.real, c.imag, s=sz, marker=mk,
                                            color=col, edgecolor='k',
                                            linewidth=0.4, zorder=3)
                        for nm, mk in zip(sig_names, markers):
                            axk.scatter([], [], marker=mk, color='0.4',
                                        edgecolor='k', linewidth=0.4, s=60,
                                        label=nm)
                        axk.legend(fontsize=8, loc='best',
                                   title='signal', title_fontsize=8)
                        axk.set_xlabel('Re(signal)')
                        axk.set_ylabel('Im(signal)')
                        axk.set_title('Complex-plane response signals per '
                                      'interior model')
                        _rasterize_heavy_artists(figk)
                        _display_vector_fig(figk, key='signals_fixed',
                                            download_label='complex-plane signals')
                        plt.close(figk)
                        st.caption(
                            "Each connected path is one interior model "
                            "(Ae/T_b grid node): its k₂, h₂, and induction "
                            "response(s) Ae on a single complex plane. "
                            "Marker shape = signal; marker size = ocean "
                            "thickness; line/marker color = T_b (per-sample "
                            "paths colored by salinity activate with a "
                            "salinity-sampled artifact). Grey clouds: "
                            "per-sample k₂/h₂ rheology spread.")

                if next((p for p in result.param_names
                         if 'wOcean' in p), None) is not None:
                    st.caption(
                        "Each point is one posterior sample's complex "
                        "induction response, from the same 2D (T_b, "
                        "salinity) response grid the likelihood/support "
                        "cut uses — Ae varies with both bottom temperature "
                        "and ocean salinity in this model family.")
                else:
                    st.caption(
                        "Each point is one posterior sample's complex "
                        "induction response, from the same per-T_b response "
                        "grid the likelihood/support cut uses. Ae depends "
                        "on the ocean state only through T_b in this "
                        "fixed-salinity model family.")
            except Exception as e:
                st.warning(f"Induction posterior panel unavailable: {e}")

    # Interactive 3D globe: textured body surface (semi-transparent) over
    # the posterior-median concentric interior layers; degree-2 (C20/C22)
    # figure deformation — even 1D inference carries a non-spherical shape
    # the C/MR^2 view never shows (user 2026-07-19).
    with st.expander("🌐 Interactive Globe", expanded=False):
        try:
            from Utilities.globe_view import (
                build_globe_figure, posterior_median_layers, texture_path)

            body = getattr(result.config, 'bodyname', '') or ''
            R_km, _ = posterior_median_layers(result)
            if R_km is None:
                st.info("Globe unavailable: this result predates the "
                        "R_body_km/layer packaging — rerun to enable.")
            else:
                # --- Sample picker: click a posterior sample to see the
                # interior it corresponds to (user 2026-07-19). Median is
                # the default until a point is selected.
                names = list(result.param_names)
                samples = np.asarray(result.samples, float)
                d_oc = np.asarray(getattr(result, 'D_ocean_results', []),
                                  float)
                d_hs = np.asarray(getattr(result, 'D_hsphere_results', []),
                                  float)
                d_ih = np.asarray(getattr(result, 'D_iceIh_results', []),
                                  float)
                d_hp = np.asarray(getattr(result, 'D_iceHP_results', []),
                                  float)
                cmr2s = np.asarray(getattr(result, 'cmr2_results', []),
                                   float)
                tb = (samples[:, names.index('Tb_K')]
                      if 'Tb_K' in names else None)
                sal_col = next((p_ for p_ in names if 'wOcean' in p_), None)
                sal = (samples[:, names.index(sal_col)]
                       if sal_col else None)

                sel_idx = None
                if tb is not None and d_oc.size == len(samples):
                    import plotly.graph_objects as go_
                    # --- X-axis selector (user 2026-07-20): plotting vs ICE
                    # THICKNESS (not Tb) reveals the full salinity distribution
                    # — a uniform-Tb view clusters high salinity at low Tb, an
                    # artifact of the Tb<->thickness<->salinity coupling. Each
                    # candidate is offered only when its array is present and
                    # has some finite values; default to ice thickness.
                    def _finite_ok(a):
                        return (a is not None and np.size(a) == len(samples)
                                and np.isfinite(a).any())
                    x_choices = []
                    if _finite_ok(d_ih):
                        x_choices.append(('Ice thickness', d_ih,
                                          'Ice thickness (km)',
                                          'ice %{x:.1f} km'))
                    if tb is not None:
                        x_choices.append(('T_b', tb, 'T_b (K)',
                                          'Tb %{x:.2f} K'))
                    if _finite_ok(d_oc):
                        x_choices.append(('Ocean thickness', d_oc,
                                          'Ocean thickness (km)',
                                          'ocean %{x:.1f} km'))
                    x_labels = [c[0] for c in x_choices]
                    x_pick = st.segmented_control(
                        'Sample-picker x-axis', x_labels,
                        default=x_labels[0], key='globe_xaxis',
                        label_visibility='collapsed')
                    if x_pick not in x_labels:
                        x_pick = x_labels[0]
                    _lab, x_arr, x_title, x_hover = x_choices[
                        x_labels.index(x_pick)]

                    n_show = min(len(samples), 1500)
                    rng = np.random.default_rng(0)
                    idx_show = (rng.choice(len(samples), n_show,
                                           replace=False)
                                if len(samples) > n_show
                                else np.arange(len(samples)))
                    marker = dict(size=5, opacity=0.6)
                    if sal is not None:
                        marker.update(color=sal[idx_show],
                                      colorscale='Viridis',
                                      colorbar=dict(title=sal_col,
                                                    thickness=12))
                    figsel = go_.Figure(go_.Scattergl(
                        x=np.asarray(x_arr, float)[idx_show],
                        y=d_oc[idx_show], mode='markers',
                        marker=marker, customdata=idx_show,
                        hovertemplate=(x_hover + ', ocean %{y:.1f} km'
                                       '<extra>click to view</extra>')))
                    figsel.update_layout(
                        height=260, margin=dict(l=0, r=0, t=24, b=0),
                        xaxis_title=x_title,
                        yaxis_title='Ocean thickness (km)',
                        title=dict(text='Click a posterior sample to view '
                                        'its interior', font=dict(size=12)))
                    ev = st.plotly_chart(
                        figsel, width='stretch', key='globe_picker',
                        on_select='rerun', selection_mode='points')
                    try:
                        pts = (ev.get('selection', {}) or {}).get('points') \
                            if isinstance(ev, dict) else ev.selection.points
                        if pts:
                            p0 = pts[0]
                            # customdata carries the ORIGINAL sample index;
                            # for a 1-D customdata array each point's value
                            # comes back as a list [idx]. Fall back to the
                            # plotted-position index mapped through idx_show.
                            cd = p0.get('customdata')
                            if isinstance(cd, (list, tuple)):
                                cd = cd[0] if cd else None
                            if cd is None:
                                pos = p0.get('point_index',
                                             p0.get('point_number'))
                                cd = (idx_show[int(pos)]
                                      if pos is not None else None)
                            sel_idx = int(cd) if cd is not None else None
                    except Exception:
                        sel_idx = None

                def _sample_layers(i):
                    lays = []
                    if 'R_core_km' in names and samples[i, names.index(
                            'R_core_km')] > 1.0:
                        lays.append({'name': 'core',
                                     'r_km': float(samples[i, names.index(
                                         'R_core_km')]),
                                     'kind': 'core'})
                    if d_hs.size > i and np.isfinite(d_hs[i]):
                        lays.append({'name': 'silicate (ocean floor)',
                                     'r_km': R_km - float(d_hs[i]),
                                     'kind': 'silicate'})
                    if (d_ih.size > i and np.isfinite(d_ih[i])
                            and d_oc.size > i and d_oc[i] > 0.5):
                        lays.append({'name': 'ocean (top)',
                                     'r_km': R_km - float(d_ih[i]),
                                     'kind': 'ocean'})
                    # HP-ice shell (III+V+VI) between ocean and silicate;
                    # absent for HP-free bodies (Europa -> d_hp ~ 0). No-ocean
                    # bodies have no d_oc[i] -> treat as 0 (HP under Ih).
                    if (d_hp.size > i and np.isfinite(d_hp[i])
                            and d_hp[i] > 0.5 and d_ih.size > i
                            and np.isfinite(d_ih[i])):
                        _doc = (float(d_oc[i]) if d_oc.size > i
                                and np.isfinite(d_oc[i]) else 0.0)
                        lays.append({'name': 'high-pressure ice (top)',
                                     'r_km': R_km - float(d_ih[i]) - _doc,
                                     'kind': 'hp_ice'})
                    return lays

                from PlanetProfile.Inference.gravity_obs import \
                    radau_darwin_kf

                obs = getattr(result.config, 'observables', {}) or {}
                c20s = getattr(result, 'c20_results', None)
                c22s = getattr(result, 'c22_results', None)
                if sel_idx is not None:
                    glayers = _sample_layers(sel_idx)
                    cmr2_i = (float(cmr2s[sel_idx])
                              if cmr2s.size > sel_idx else np.nan)
                    c20 = (float(c20s[sel_idx]) if c20s is not None
                           else obs.get('C20', [None])[0])
                    c22 = (float(c22s[sel_idx]) if c22s is not None
                           else obs.get('C22', [None])[0])
                    which = (f"sample #{sel_idx}: " + ", ".join(
                        f"{n} = {samples[sel_idx, j]:.3g}"
                        for j, n in enumerate(names)))
                else:
                    _, glayers = posterior_median_layers(result)
                    _c = cmr2s[np.isfinite(cmr2s)] if cmr2s.size else []
                    cmr2_i = float(np.median(_c)) if len(_c) else np.nan
                    c20 = obs.get('C20', [None])[0]
                    c22 = obs.get('C22', [None])[0]
                    which = "posterior median (click a point above to pick a sample)"
                kf = None
                if np.isfinite(cmr2_i):
                    try:
                        kf = radau_darwin_kf(cmr2_i)
                    except ValueError:
                        kf = None

                # Runs without C20/C22 observables (e.g. CMR2-only
                # configs) still get a figure: the hydrostatic pair
                # implied by the posterior k_f and the body's rotation —
                # otherwise the exaggeration slider has nothing to act
                # on (user-reported 2026-07-19).
                figure_src = "measured C20/C22 observables"
                if c20 is None and c22 is None:
                    from Utilities.globe_view import \
                        hydrostatic_display_pair
                    c20, c22 = hydrostatic_display_pair(body, R_km, kf)
                    figure_src = ("hydrostatic figure from the posterior "
                                  "k_f (no C20/C22 in this run)")

                dev = abs(c20 or 0.0) + 3 * abs(c22 or 0.0)
                default_ex = (min(500.0, max(1.0, 0.03 / dev))
                              if dev > 0 else 1.0)
                exagg = st.slider(
                    "Shape exaggeration (1 = true figure)", 1.0, 500.0,
                    float(round(default_ex)), 1.0, key='globe_exagg')

                fig3d = build_globe_figure(
                    body, R_km, glayers, c20=c20, c22=c22, kf=kf,
                    exaggeration=exagg)
                st.plotly_chart(fig3d, width='stretch', key='globe_chart')
                from Utilities.globe_view import triaxial_semiaxes
                _kf_amp = (1.0 + (kf or 1.0)) / (kf or 1.0)
                _a, _b, _c = triaxial_semiaxes(
                    R_km, (c20 or 0.0), (c22 or 0.0), _kf_amp, 1.0)
                st.caption(
                    f"Principal axes (true figure, no exaggeration): "
                    f"a = {_a:.2f} km (sub-parent, moment A), "
                    f"b = {_b:.2f} km (orbital, B), "
                    f"c = {_c:.2f} km (spin, C); "
                    f"a − c = {(_a - _c) * 1000:.0f} m. "
                    f"Figure source: {figure_src}.")
                tex_note = ("NASA/USGS global mosaic (public domain)"
                            if texture_path(body) else
                            "no shipped texture for this body yet — "
                            "shaded sphere")
                shape_note = (
                    "Surface deformed by the degree-2 figure implied by "
                    "C20/C22 (tidal bulge along the sub-parent x-axis, "
                    "polar flattening) with fluid amplification "
                    "(1+k_f)/k_f — the non-spherical shape the "
                    "spherically averaged C/MR² never captures."
                    if (c20 or c22) else
                    "No C20/C22 in this run — sphere shown; gravity-pair "
                    "configs deform the figure.")
                st.caption(
                    f"Showing {which}. Drag to rotate, scroll to zoom. "
                    f"Texture: {tex_note}. Cutaway wedge exposes the "
                    f"interior; the innermost body is solid. {shape_note} "
                    "Lateral (3D) structure becomes per-layer r(θ, φ) "
                    "fields on this same mesh — roadmap.")
        except Exception as e:
            st.warning(f"Globe unavailable: {e}")

    # v4 geodesy: the reviewer-binding REPORTABLE quantity — the
    # identifiable non-hydrostatic combination u = dC22_nh + dC20_nh/R.
    # Per-component marginals are NOT calibrated (do not cite).
    if ('dC20_nh' in result.param_names
            and 'dC22_nh' in result.param_names):
        with st.expander("🪨 Non-hydrostatic gravity: identifiable "
                         "combination u", expanded=True):
            try:
                import matplotlib.pyplot as plt
                md_cfg = (getattr(result.config, 'metadata', {}) or {})
                R_ratio = float(md_cfg.get('gravity_j2_over_c22', 3.324))
                _s = np.asarray(result.samples, float)
                u = (_s[:, result.param_names.index('dC22_nh')]
                     + _s[:, result.param_names.index('dC20_nh')] / R_ratio)
                obs = getattr(result.config, 'observables', {}) or {}
                s20 = float(obs.get('C20', [0, np.nan])[1])
                s22 = float(obs.get('C22', [0, np.nan])[1])
                sigma_u = float(np.sqrt(s22 ** 2 + (s20 / R_ratio) ** 2))

                med = float(np.median(u))
                lo68, hi68 = np.percentile(u, [16, 84])
                ul95 = float(np.percentile(np.abs(u), 95))

                c1, c2, c3 = st.columns(3)
                c1.metric("median u", f"{med:.2e}")
                c2.metric("95% upper limit |u|", f"{ul95:.2e}")
                c3.metric("|u|₉₅ / σ_u(obs)",
                          f"{ul95 / sigma_u:.1f}"
                          if np.isfinite(sigma_u) and sigma_u > 0 else "—")

                figu, axu = plt.subplots(figsize=(7, 3.2))
                axu.hist(u, bins=60, color='peru', alpha=0.8)
                axu.axvline(0, color='k', linewidth=0.8)
                axu.axvline(med, color='crimson', linewidth=1.2,
                            label=f'median {med:.2e}')
                axu.axvspan(lo68, hi68, color='crimson', alpha=0.12,
                            label='68%')
                if np.isfinite(sigma_u) and sigma_u > 0:
                    axu.axvspan(-sigma_u, sigma_u, color='steelblue',
                                alpha=0.10,
                                label=f'±σ_u(obs) = {sigma_u:.2e}')
                axu.set_xlabel(f'u = dC22_nh + dC20_nh / {R_ratio:g} '
                               '(unnormalized)')
                axu.set_ylabel('posterior samples')
                axu.legend(fontsize=8)
                _display_vector_fig(figu, key='nonhydro_u',
                                    download_label='non-hydrostatic u')
                plt.close(figu)
                st.caption(
                    "u is the RATIO-BREAKING combination of the sampled "
                    "non-hydrostatic Stokes offsets — the interior term "
                    "cancels exactly, making u interior-independent and "
                    "SBC-calibrated (the v4 deliverable). Report u as a "
                    "provisional upper limit in absolute and observation-σ "
                    "units. Do NOT cite the per-component dC20_nh/dC22_nh "
                    "marginals or the −C20/C22 ratio: the ratio-preserving "
                    "direction is pinned only through the weakly "
                    "identified interior. σ_u convention risk: if the "
                    "mission's σ(C20) is 4π-normalized, σ_u is ~1.4× "
                    "larger than shown.")
            except Exception as e:
                st.warning(f"u panel unavailable: {e}")

    # C/MR² posterior
    with st.expander("⚖️ C/MR² Moment-of-Inertia Posterior", expanded=True):
        with st.popover("📖 How degree-2 gravity constrains C/MR²"):
            st.markdown(
                "A body's degree-2 gravity field measures **differences** "
                "of its principal moments of inertia $A \\le B \\le C$ "
                "(for synchronous Europa: $A$ along the Jupiter-facing "
                "axis, $C$ along the spin axis). With unnormalized Stokes "
                "coefficients (Mazarico et al. 2023, Eq. 5):")
            st.latex(r"C_{20} = -J_2 = -\frac{C - \tfrac{1}{2}(A + B)}"
                     r"{M R^2}, \qquad C_{22} = \frac{B - A}{4 M R^2}")
            st.markdown(
                "So $C_{20}$ senses the polar flattening of the mass "
                "distribution and $C_{22}$ the equatorial (tidal) "
                "elongation — but neither gives the **mean** moment "
                "$C/MR^2$ (the differentiation measure this panel shows) "
                "directly. Getting $C/MR^2$ requires the **hydrostatic "
                "assumption**: for a fluid body in equilibrium, figure "
                "theory ties the fluid Love number $k_f$ — and hence the "
                "ratio $J_2/C_{22} = 10/3$ (3.324 with Europa's "
                "rapid-rotation correction, Tricarico 2014) — to "
                "$C/MR^2$ via the Radau–Darwin relation:")
            st.latex(r"\frac{C}{MR^2} = \frac{2}{3}\left[1 - \frac{2}{5}"
                     r"\sqrt{\frac{4 - k_f}{1 + k_f}}\,\right],\qquad "
                     r"C_{22} = \frac{k_f\, q_r}{4},\quad "
                     r"q_r = \frac{\omega^2 R^3}{G M}")
            st.markdown(
                "This panel shows **two** C/MR² distributions:\n\n"
                "- **Actual** (blue) — computed directly from the inferred "
                "interior density profile by exact moment-of-inertia "
                "integration; makes **no** hydrostatic assumption.\n"
                "- **Hydrostatic reference** (green) — the C/MR² the model's "
                "$C_{22}$ would imply *if* the body were hydrostatic, via the "
                "Radau–Darwin relation above.\n\n"
                "The **gap between them is the inferred non-hydrostaticity**. "
                "Because Radau–Darwin carries a ~1% $k_f$ systematic "
                "(≈0.0035 in C/MR²), gaps below that floor are the RD "
                "approximation, not a physical departure — the readout flags "
                "this. The Galileo-era C/MR² (0.3547 ± 0.0024, Gomez Casajus "
                "et al. 2021) is itself the hydrostatic reduction of the "
                "Galileo $C_{22}$; in the free-gravity design it is **not** "
                "imposed as an independent constraint (that would double-count "
                "$C_{22}$), but may be applied optionally as a literature "
                "reweight.")
        cmr2_obs = result.config.observables.get('CMR2')
        cmr2_results = getattr(result, 'cmr2_results', None)
        # Hydrostatic-REFERENCE C/MR² (Radau–Darwin from the model C22): the
        # value the body WOULD have if hydrostatic. Paired with the ACTUAL
        # (structure moment-integral) cmr2_results, the gap between the two is
        # the inferred non-hydrostaticity. None on older pkls / non-gravity runs.
        cmr2_hydro_results = getattr(result, 'cmr2_hydro_results', None)
        # ~1% Radau–Darwin k_f systematic in C/MR² units: an interpretation
        # FLOOR — gaps below this are the RD approximation, not physics.
        RD_FLOOR = 0.0035

        if cmr2_results is None or not np.any(np.isfinite(cmr2_results)):
            if cmr2_obs is None:
                st.info("C/MR² was not included as an observable in this run. "
                        "Enable it in **Observables** to add a moment-of-inertia constraint.")
            else:
                st.warning("C/MR² observable was specified but no finite values were computed. "
                           "Rebuild the structure cache so CMR2 is extracted from Planet.CMR2mean.")
        else:
            try:
                import matplotlib.pyplot as plt

                finite_mask = np.isfinite(cmr2_results)
                cmr2_vals = cmr2_results[finite_mask]

                # Detect effectively-fixed CMR2: either truly constant (no
                # structural parameters free) or varying only at the cache-
                # discretization level (core-blind models like Test50: the
                # 9-point Tb grid moves CMR2 by ~2.5e-5). Rendering that
                # numerical noise as a full histogram produces a misleading
                # offset-notation axis (ticks in 1e-5 above the mean), so
                # treat anything far below the observation's sigma (or an
                # absolute floor) as fixed and show the constraint-band view.
                cmr2_range = cmr2_vals.max() - cmr2_vals.min()
                _obs_err_for_scale = (cmr2_obs[1] if cmr2_obs is not None
                                      else None)
                fixed_tol = (0.1 * _obs_err_for_scale
                             if _obs_err_for_scale else 1e-4)
                is_fixed = cmr2_range < max(fixed_tol, 1e-8)

                fig, ax = plt.subplots(figsize=(7, 4))

                if cmr2_obs is not None:
                    obs_val, obs_err = cmr2_obs
                else:
                    obs_val, obs_err = None, None

                if is_fixed:
                    # Single value — show as a vertical line against the
                    # observed constraint bands so the tension is visible.
                    model_cmr2 = cmr2_vals[0]
                    if obs_val is not None:
                        x_lo = min(model_cmr2, obs_val) - 6 * obs_err
                        x_hi = max(model_cmr2, obs_val) + 6 * obs_err
                        x_range = np.linspace(x_lo, x_hi, 400)
                        gauss = np.exp(-0.5 * ((x_range - obs_val) / obs_err) ** 2)
                        ax.plot(x_range, gauss, 'r-', linewidth=2,
                                label=fr'Observed: {obs_val:.4f} ± {obs_err:.4f}')
                        ax.axvspan(obs_val - obs_err, obs_val + obs_err,
                                   alpha=0.15, color='red', label=r'1$\sigma$')
                        ax.axvspan(obs_val - 2 * obs_err, obs_val + 2 * obs_err,
                                   alpha=0.07, color='red', label=r'2$\sigma$')
                    ax.axvline(model_cmr2, color='steelblue', linewidth=2.5,
                               label=f'Model: {model_cmr2:.5f}')
                    ax.set_xlabel(r'$C/MR^2$')
                    ax.set_ylabel('Likelihood (observed)')
                    ax.set_title(r'C/MR² — Fixed Structure')
                    ax.legend(fontsize=9)

                    st.info(
                        f"C/MR² is **effectively constant** across the posterior for this "
                        f"model (range {cmr2_range:.2g} — cache-discretization level). "
                        f"This happens when no core/structural parameters are sampled "
                        f"(core-blind models like Test50): C/MR² acts as a uniform χ² "
                        f"offset and does not discriminate between rheological models. "
                        f"For a C/MR²-constrained posterior use a differentiated-core "
                        f"model (Test52-family)."
                    )
                else:
                    # Variable CMR2 — full histogram
                    ax.hist(cmr2_vals, bins=40, density=True, alpha=0.7,
                            color='steelblue',
                            label='Actual (structure integral)')

                    # Overlay the hydrostatic-reference C/MR² (RD-from-C22):
                    # what the body would be IF hydrostatic. The separation
                    # from the actual distribution IS the non-hydrostaticity.
                    cmr2_hydro_vals = None
                    if cmr2_hydro_results is not None:
                        _hm = np.isfinite(cmr2_hydro_results)
                        if np.any(_hm):
                            cmr2_hydro_vals = cmr2_hydro_results[_hm]
                            ax.hist(cmr2_hydro_vals, bins=40, density=True,
                                    alpha=0.45, color='seagreen',
                                    label='Hydrostatic reference (RD from $C_{22}$)')

                    if obs_val is not None:
                        x_range = np.linspace(
                            min(cmr2_vals.min(), obs_val - 4 * obs_err),
                            max(cmr2_vals.max(), obs_val + 4 * obs_err),
                            300,
                        )
                        gauss = np.exp(-0.5 * ((x_range - obs_val) / obs_err) ** 2)
                        hist_counts, _ = np.histogram(cmr2_vals, bins=40, density=True)
                        peak = hist_counts.max()
                        ax.plot(x_range, gauss * peak, 'r-', linewidth=2,
                                label=fr'Observed: {obs_val:.4f} ± {obs_err:.4f}')
                        ax.axvline(obs_val, color='red', linestyle='--',
                                   linewidth=1.5, alpha=0.7)
                        ax.axvspan(obs_val - obs_err, obs_val + obs_err,
                                   alpha=0.12, color='red', label=r'1$\sigma$')
                        ax.axvspan(obs_val - 2 * obs_err, obs_val + 2 * obs_err,
                                   alpha=0.06, color='red', label=r'2$\sigma$')

                    ax.set_xlabel(r'$C/MR^2$')
                    ax.set_ylabel('Probability density')
                    ax.set_title(r'Moment-of-Inertia Posterior')
                    ax.legend(fontsize=9)

                fig.tight_layout()
                _display_vector_fig(fig, key='cmr2_posterior',
                                    download_label='C/MR² posterior')
                plt.close(fig)

                # Summary metrics
                model_val = cmr2_vals[0] if is_fixed else np.median(cmr2_vals)
                col1, col2, col3 = st.columns(3)
                col1.metric("Model C/MR²" if is_fixed else "Median C/MR²",
                            f"{model_val:.5f}")
                if not is_fixed:
                    q16, q84 = np.percentile(cmr2_vals, [16, 84])
                    col2.metric("16th–84th %ile", f"{q16:.5f} – {q84:.5f}")
                if obs_val is not None:
                    tension = abs(model_val - obs_val) / obs_err
                    delta = model_val - obs_val
                    col3.metric("Tension (σ)", f"{tension:.1f}σ",
                                delta=f"{delta:+.5f}",
                                delta_color="inverse")

                if not is_fixed and obs_val is not None:
                    n_out = np.sum(np.abs(cmr2_vals - obs_val) > 2 * obs_err)
                    if n_out > 0:
                        st.caption(
                            f"{n_out} / {len(cmr2_vals)} samples "
                            f"({n_out/len(cmr2_vals):.1%}) fall outside the 2σ constraint."
                        )

                # Non-hydrostaticity: the gap between the ACTUAL C/MR² (from
                # the interior structure) and the HYDROSTATIC REFERENCE (what
                # C22 implies if hydrostatic). Reported per-sample (paired) so
                # the interior/gravity correlation is preserved. Read against
                # the ~1% Radau–Darwin floor — gaps below it are the RD
                # approximation, not physics.
                if (cmr2_hydro_results is not None
                        and np.any(np.isfinite(cmr2_hydro_results))):
                    _pair = np.isfinite(cmr2_results) & np.isfinite(cmr2_hydro_results)
                    if np.any(_pair):
                        gap = cmr2_results[_pair] - cmr2_hydro_results[_pair]
                        gap_med = float(np.median(gap))
                        gap_abs_med = float(np.median(np.abs(gap)))
                        st.markdown("#### Non-hydrostaticity (actual − hydrostatic reference)")
                        gc1, gc2, gc3 = st.columns(3)
                        gc1.metric("Median gap", f"{gap_med:+.5f}")
                        gc2.metric("Median |gap|", f"{gap_abs_med:.5f}")
                        gc3.metric("RD floor", f"±{RD_FLOOR:.4f}")
                        if gap_abs_med < RD_FLOOR:
                            st.caption(
                                f"Median |gap| {gap_abs_med:.5f} is **below** the "
                                f"~{RD_FLOOR:.4f} Radau–Darwin systematic floor — "
                                "consistent with hydrostatic equilibrium; the "
                                "residual is the RD closed-form approximation, "
                                "**not** a measured departure. Do not interpret "
                                "sub-floor gaps as physical non-hydrostaticity."
                            )
                        else:
                            st.caption(
                                f"Median |gap| {gap_abs_med:.5f} **exceeds** the "
                                f"~{RD_FLOOR:.4f} Radau–Darwin floor — a candidate "
                                "non-hydrostatic signal (the actual moment of "
                                "inertia departs from the hydrostatic C₂₂ "
                                "reduction beyond the RD approximation)."
                            )

            except Exception as e:
                st.warning(f"C/MR² posterior plot unavailable: {e}")

    # Heating distribution
    with st.expander("🔥 Heating Distribution", expanded=False):
        heating_results = list(result.heating_results) if result.heating_results is not None else []
        n_nonempty = sum(1 for h in heating_results if h)
        if not heating_results or n_nonempty == 0:
            st.info("Heating data unavailable. Set **n_reeval > 0** in sampler "
                    "settings (or the heating re-evaluation subset in amortized "
                    "mode) to compute per-phase heating.")
            st.caption(
                f"Diagnostic: result carries {len(heating_results)} heating "
                f"entries ({n_nonempty} non-empty); sampler = "
                f"{(result.metadata or {}).get('sampler')}, generated "
                f"{(result.metadata or {}).get('elapsed_time_s', 0):.0f}s run. "
                f"If this run was made with a re-evaluation subset > 0, the "
                f"app server is likely running stale code — restart it."
            )
        else:
            st.caption(f"{n_nonempty} forward-model heating evaluations "
                       f"(seeded subset of the posterior).")
            try:
                import matplotlib.pyplot as plt

                eval_idx = result.metadata.get('heating_indices',
                                               list(range(len(heating_results))))
                eval_samples = result.samples[eval_idx]
                n_params = eval_samples.shape[1]

                # Get heating for each reservoir/phase
                silicate_vals = []
                ice_ih_vals = []
                ice_iii_vals = []
                ice_v_vals = []
                ice_vi_vals = []

                for h in heating_results:
                    silicate_vals.append(h.get('Silicate_W', h.get('Sil', 0.0) + h.get('Fe', 0.0)))
                    ice_ih_vals.append(h.get('Ice_Ih_W', h.get('Ih', 0.0) + h.get('Clath', 0.0)))
                    ice_iii_vals.append(h.get('III', 0.0))
                    ice_v_vals.append(h.get('V', 0.0))
                    ice_vi_vals.append(h.get('VI', 0.0))

                # Per-model reservoir/phase fractions (Stacked Bar)
                st.markdown("#### 📉 Per-Model Partitioning")

                n_models = len(heating_results)
                sil_arr = np.array(silicate_vals)
                ih_arr = np.array(ice_ih_vals)
                iii_arr = np.array(ice_iii_vals)
                v_arr = np.array(ice_v_vals)
                vi_arr = np.array(ice_vi_vals)
                
                tot_arr = sil_arr + ih_arr + iii_arr + v_arr + vi_arr + 1e-30

                # Sort models by silicate fraction
                f_sil_arr = sil_arr / tot_arr
                sort_order = np.argsort(f_sil_arr)

                fig_stack, ax_stack = plt.subplots(figsize=(10, 4))
                x = np.arange(n_models)
                bottom = np.zeros(n_models)

                # Colors matching PlanetProfile conventions
                plot_colors = {
                    'Silicate': '#9c755f', 
                    'Ice VI': '#d62728', # Red
                    'Ice V': '#2ca02c',  # Green
                    'Ice III': '#ff7f0e', # Orange
                    'Ice Ih': '#4e79a7'   # Blue
                }

                # Silicate
                ax_stack.bar(x, sil_arr[sort_order]/tot_arr[sort_order], bottom=bottom, color=plot_colors['Silicate'], label='Silicate', width=1.0)
                bottom += sil_arr[sort_order]/tot_arr[sort_order]
                
                # HP Ices (Individual)
                ax_stack.bar(x, vi_arr[sort_order]/tot_arr[sort_order], bottom=bottom, color=plot_colors['Ice VI'], label='Ice VI', width=1.0)
                bottom += vi_arr[sort_order]/tot_arr[sort_order]
                
                ax_stack.bar(x, v_arr[sort_order]/tot_arr[sort_order], bottom=bottom, color=plot_colors['Ice V'], label='Ice V', width=1.0)
                bottom += v_arr[sort_order]/tot_arr[sort_order]
                
                ax_stack.bar(x, iii_arr[sort_order]/tot_arr[sort_order], bottom=bottom, color=plot_colors['Ice III'], label='Ice III', width=1.0)
                bottom += iii_arr[sort_order]/tot_arr[sort_order]
                
                # Ice Ih
                ax_stack.bar(x, ih_arr[sort_order]/tot_arr[sort_order], bottom=bottom, color=plot_colors['Ice Ih'], label='Ice Ih', width=1.0)

                ax_stack.set_xlim(0, n_models-1)
                ax_stack.set_ylim(0, 1)
                ax_stack.set_ylabel("Fraction")
                ax_stack.set_xlabel("Samples (sorted by silicate fraction)")
                ax_stack.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)

                _display_vector_fig(fig_stack, key='reservoir_stack',
                                    download_label='reservoir fractions')
                plt.close(fig_stack)

                # Aggregate Reservoir Visualization
                st.markdown("#### 🧊 Reservoir Partitioning (Median)")
                col_pie, col_stats = st.columns([1, 1])

                sil_med = np.median(silicate_vals)
                ih_med = np.median(ice_ih_vals)
                hp_med = np.median(iii_arr + v_arr + vi_arr)
                tot_med = sil_med + ih_med + hp_med

                with col_pie:
                    if tot_med > 0:
                        fig_pie, ax_pie = plt.subplots(figsize=(5, 5))
                        ax_pie.pie([sil_med, hp_med, ih_med],
                                   labels=['Silicate', 'HP Ice', 'Ice Ih'],
                                   autopct='%1.1f%%',
                                   colors=['#C8A96E', '#9B59B6', '#AEE1F8'],
                                   startangle=90)
                        ax_pie.axis('equal')
                        _display_vector_fig(fig_pie, key='heating_pie',
                                            download_label='heating pie')
                        plt.close(fig_pie)

                with col_stats:
                    st.metric("Total Power (Median)", f"{tot_med/1e9:.2f} GW")
                    st.write(f"- **Silicate:** {sil_med/1e9:.2f} GW")
                    st.write(f"- **HP Ices:** {hp_med/1e9:.2f} GW")
                    st.write(f"- **Ice Ih Shell:** {ih_med/1e9:.2f} GW")

                    if tot_med > 0:
                        st.info(f"Titan's heat budget is dominated by the **{['Silicate', 'HP Ice', 'Ice Ih'][np.argmax([sil_med, hp_med, ih_med])]}** reservoir.")

                st.markdown("---")

                # Consistent phase color scheme used across all plots
                phase_colors = {
                    'Clath': '#D4F1F9',
                    'Ih':    '#AEE1F8',
                    '0':     '#1E90FF',
                    'III':   '#C97BAE',
                    'II':    '#B0E0E6',
                    'V':     '#9B59B6',
                    'VI':    '#6C3483',
                    'Sil':   '#C8A96E',
                    'Rock':  '#C8A96E',
                    'Core':  '#8B5A2B',
                }

                # Use only phases actually present in the data.
                # (actual_keys was previously undefined here — a NameError
                # silently disabled this whole sub-panel via the outer except.)
                actual_keys = (set().union(*(h.keys() for h in heating_results if h))
                               if any(heating_results) else set())
                phases_to_show = [p for p in ['Ih', 'III', 'II', 'V', 'VI', 'Sil', 'Core', 'Clath']
                                  if p in actual_keys]

                heating_fracs = {}
                for phase in phases_to_show:
                    fracs = []
                    for h in heating_results:
                        total = _phase_heating_total(h) if h else 0.0
                        fracs.append(h.get(phase, 0) / total if total > 0 else 0)
                    heating_fracs[phase] = np.array(fracs)

                n_panels = n_params + 1
                n_cols = 3
                n_rows = (n_panels + n_cols - 1) // n_cols
                fig, axes = plt.subplots(n_rows, n_cols,
                                         figsize=(6 * n_cols, 5 * n_rows))
                axes = axes.flatten()
                for ax in axes[n_panels:]:
                    ax.set_visible(False)

                for ip, plabel in enumerate(result.param_labels):
                    ax = axes[ip]
                    x = eval_samples[:, ip]
                    for phase in phases_to_show:
                        ax.scatter(x, heating_fracs[phase], s=4, alpha=0.3,
                                   color=phase_colors.get(phase, 'gray'),
                                   label=phase)
                    ax.set_xlabel(plabel)
                    ax.set_ylabel('Heating fraction')
                    ax.set_ylim(-0.05, 1.05)
                    if ip == 0:
                        ax.legend(fontsize=8, loc='best')

                # Stacked bar sorted by silicate fraction
                ax = axes[n_params]
                f_sil_heat = heating_fracs.get('Sil', np.zeros(len(heating_results)))
                sort_order = np.argsort(f_sil_heat)
                bottom = np.zeros(len(heating_results))
                for phase in ['Core', 'Sil', 'VI', 'V', 'III', 'II', 'Ih', 'Clath']:
                    if phase in heating_fracs:
                        ax.bar(range(len(heating_results)),
                               heating_fracs[phase][sort_order],
                               bottom=bottom,
                               color=phase_colors.get(phase, 'gray'),
                               label=phase, width=1.0)
                        bottom += heating_fracs[phase][sort_order]
                ax.set_xlabel('Samples sorted by silicate fraction')
                ax.set_ylabel('Heating fraction')
                ax.legend(fontsize=8)
                ax.set_title('Per-phase heating across posterior')

                fig.suptitle('Heating Distribution', fontsize=14)
                fig.tight_layout()
                _rasterize_heavy_artists(fig)
                _display_vector_fig(fig, key='heating_distribution',
                                    download_label='heating distribution')
                plt.close(fig)
            except Exception as e:
                st.warning(f"Heating distribution unavailable: {e}")

    # Export button
    st.markdown("---")
    if st.button("💾 Export Results (PKL)", use_container_width=True):
        import pickle
        output_path = Path(parent_directory) / "inference_result.pkl"

        with open(output_path, 'wb') as f:
            pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)

        st.success(f"✅ Results saved to: {output_path}")


# ============================================================
# Browse-runs mode (T07b) — load saved *_result.pkl files and reuse render_results
# ============================================================
def render_browse_runs():
    """Browse and load previously saved MCMC result PKLs.

    Scans ``PlanetProfile/Test/mcmc_results`` (and any extra roots in
    DEFAULT_RESULTS_ROOTS) for ``*_result.pkl`` files, lets the user pick one
    via a selectbox, loads it into ``st.session_state.inference_results``, and
    surfaces summary metrics (log_Z, ESS, acceptance, n_samples).  After
    loading, the existing ``render_results()`` panel renders posterior plots.
    """
    st.markdown("### 📂 Browse Saved MCMC Runs")
    st.caption(
        "Loads serialized ``InferenceResult`` objects from disk and reuses the "
        "results-rendering panel to show summary metrics, corner plot, k₂ scatter, "
        "and C/MR² posterior."
    )

    # Scan all default roots and accumulate records
    all_records: list[dict] = []
    for root in DEFAULT_RESULTS_ROOTS:
        # Resolve relative to repo root (parent_directory)
        abs_root = root if os.path.isabs(root) else os.path.join(parent_directory, root)
        all_records.extend(scan_mcmc_results(abs_root))

    if not all_records:
        st.info(
            "No saved MCMC runs found. Run `Run new MCMC` mode (or copy "
            "``*_result.pkl`` files into ``PlanetProfile/Test/mcmc_results/``)."
        )
        return

    # Selectbox over all runs
    labels = [get_run_display_label(r) for r in all_records]
    chosen_idx = st.selectbox(
        "Select a run to load:",
        options=list(range(len(all_records))),
        format_func=lambda i: labels[i],
        key="browse_runs_selectbox",
    )
    record = all_records[chosen_idx]

    # Load button — loading is cached on (path, mtime)
    col_load, col_clear = st.columns([1, 1])
    with col_load:
        if st.button("Load run", key="browse_runs_load_btn", type="primary"):
            with st.spinner(f"Loading {record['name']}..."):
                result = load_mcmc_result(record["path"], record["mtime"])
            if result is None:
                st.error(f"Failed to load: {record['path']}")
            else:
                st.session_state.inference_results = result
                st.success(f"Loaded: {record['name']}")
    with col_clear:
        if st.button("Clear loaded run", key="browse_runs_clear_btn"):
            st.session_state.inference_results = None
            st.rerun()

    # Show summary metrics for the *currently loaded* result (if any)
    res = st.session_state.inference_results
    if res is not None:
        with st.expander("📈 Run summary metrics", expanded=True):
            metrics = getattr(res, "convergence_metrics", {}) or {}
            metadata = getattr(res, "metadata", {}) or {}

            col1, col2, col3 = st.columns(3)
            log_z = metadata.get("log_Z")
            log_z_err = metadata.get("log_Z_err")
            if log_z is not None:
                err_str = f" ± {log_z_err:.3f}" if log_z_err is not None else ""
                col1.metric("log Z", f"{log_z:.3f}{err_str}")
            else:
                col1.metric("log Z", "N/A")

            ess = metrics.get("ess")
            col2.metric("ESS", f"{ess:.0f}" if ess is not None else "N/A")

            acc = metrics.get("acceptance_rate")
            col3.metric("Acceptance",
                        f"{acc:.1%}" if acc is not None else "N/A")

            col4, col5, col6 = st.columns(3)
            n_samples = metrics.get("n_samples")
            col4.metric("Samples", f"{n_samples}" if n_samples is not None else "N/A")
            r_hat = metrics.get("r_hat")
            col5.metric("R-hat", f"{r_hat:.3f}" if r_hat is not None else "N/A")
            elapsed = metadata.get("elapsed_time_s")
            col6.metric("Runtime",
                        f"{elapsed/60:.1f} min" if elapsed is not None else "N/A")

            samples = getattr(res, "samples", None)
            param_names = getattr(res, "param_names", None)
            if samples is not None:
                st.caption(f"Posterior shape: {samples.shape}")
            if param_names:
                st.caption(f"Parameters: {', '.join(param_names)}")


# ============================================================
# Main Page
# ============================================================
def main():
    """Main Inference tab page."""

    # Initialize session state
    initialize_session_state()

    # Show planet status in sidebar
    show_planet_status()

    # Page title
    st.title("🧮 Bayesian Inference")

    # Note: Planet selection is optional for inference - body name comes from preset
    # (e.g., 'andrade_titan' preset uses Titan structure cache)

    # Mode toggle: run new MCMC vs browse saved runs (T07b).
    # Browse mode reuses the existing render_results() panel after loading a PKL
    # via PlanetProfileApp/utils/mcmc_run_loader.py.
    # Track previous mode so we can clear ghost-state when switching modes.
    _prev_mode = st.session_state.get("inference_page_mode_prev")
    inference_page_mode = st.radio(
        "Mode",
        options=["run", "browse"],
        format_func=lambda m: ("🚀 Run new MCMC" if m == "run"
                               else "📂 Browse saved runs"),
        index=0,
        horizontal=True,
        key="inference_page_mode",
        help="`Run new MCMC` configures and launches a fresh run.  "
             "`Browse saved runs` loads a previously saved *_result.pkl and "
             "renders its posterior summary, corner plot, and metrics."
    )
    # When switching browse → run, drop the loaded PKL so the run-mode results
    # column doesn't ghost-display the previously loaded posterior.
    if (_prev_mode == "browse" and inference_page_mode == "run"
            and st.session_state.get("inference_results") is not None):
        st.session_state.inference_results = None
    st.session_state["inference_page_mode_prev"] = inference_page_mode

    if inference_page_mode == "browse":
        # Browse-runs mode — load a saved PKL and reuse render_results
        render_browse_runs()
        st.markdown("---")
        st.subheader("📊 Results")
        render_results()
        return

    # --- Run mode (default) ---
    # Execution method: full MCMC sampling vs instant conditioning of a
    # pretrained (amortized) SBI artifact. The results pipeline is shared.
    if public_mode():
        # Public serving: amortized only. MCMC launches hour-scale compute
        # on a shared host; hide it rather than cap it.
        exec_mode = 'amortized'
        st.info("🌐 **Public demo** — amortized (pretrained SBI) inference "
                "only. Full MCMC sampling and free priors/σ are available "
                "by running the app from the "
                "[PlanetProfile repository](https://github.com/vancesteven/PlanetProfile).")
    else:
        exec_mode = st.radio(
            "Inference method",
            options=['mcmc', 'amortized'],
            format_func=lambda m: ("🧮 MCMC (pocomc — full sampling, free priors/σ)"
                                   if m == 'mcmc'
                                   else "⚡ Amortized (pretrained SBI — instant, frozen priors/σ)"),
            horizontal=True,
            key='inference_exec_mode',
            help="Amortized mode conditions a pretrained neural posterior "
                 "estimator: instant, but priors and observation σ are fixed at "
                 "training time. MCMC samples the full likelihood with whatever "
                 "priors/σ you configure."
        )

    # Config-file selector (body + variant) — sits above the run-configuration form.
    # Populates the same session-state keys that the form widgets read so fields
    # pre-fill after a config is applied. (MCMC mode only — amortized configs
    # come from the artifact itself.)
    if exec_mode == 'mcmc':
        render_config_file_selector()

    # Lazy import inference modules
    imports = lazy_import_inference()
    if imports is None:
        st.error("Failed to load inference modules. Check installation.")
        st.stop()

    (PARAMETER_REGISTRY, PARAMETER_PRESETS, CATEGORY_LABELS,
     CATEGORY_ORDER, validate_parameter_combination,
     get_parameters_by_rheology, InferenceConfig, MCMCRunner) = imports

    # Divider
    st.markdown("---")

    # Two-column layout: Config (left) + Results (right)
    col_config, col_results = st.columns([1, 2])

    if exec_mode == 'amortized':
        with col_config:
            st.subheader("⚡ Amortized Configuration")
            render_model_assumptions('amortized')
            spec = render_amortized_config()
            render_amortized_run_button(spec, InferenceConfig)
        with col_results:
            st.subheader("📊 Results")
            # Controls (sliders, observables, seed, ...) do NOT auto-rerun
            # the inference — flag when the displayed posterior no longer
            # matches the form so the required button click is obvious.
            if (spec is not None
                    and st.session_state.get('inference_results') is not None
                    and st.session_state.get('amort_last_run_fp') is not None
                    and _amortized_spec_fingerprint(spec)
                        != st.session_state.amort_last_run_fp):
                st.warning("⚠️ **Controls changed since this posterior was "
                           "generated** — the plots below show the *previous* "
                           "settings. Click **⚡ Generate Posterior** to update.")
            render_results()
        return

    with col_config:
        st.subheader("⚙️ Configuration")

        render_model_assumptions('mcmc')

        if st.button("🔄 Reset Configuration", help="Clear all inference settings and reload defaults"):
            for key in ['inference_param_space', 'inference_selected_params', 'inference_preset',
                        'inference_custom_mode', 'inference_observables', 'inference_sampler_settings',
                        'inference_structure_cache_path', 'inference_use_clathrate', 'inference_results',
                        'inference_cache_key', 'inference_error', 'inference_error_traceback',
                        '_default_cfg_applied']:
                if key in st.session_state:
                    del st.session_state[key]
            st.rerun()

        # The old configuration-preset radio was redundant with the config-
        # file loader above (and historically stomped loaded configs). A
        # fresh session instead auto-loads the Titan no-ocean 8D config as
        # the worked example of how an MCMC run is set up; every other
        # configuration comes from "Load config file". Only when the
        # globally chosen planet actually IS Titan (or unset) — auto-
        # loading a Titan config under a Europa session planted Titan k2
        # in the observables (user 2026-07-20).
        _global_planet = str(st.session_state.get('Planet', '') or '')
        if (not st.session_state.get('_default_cfg_applied')
                and _global_planet.lower() in ('', 'titan',
                                               '-- select a planet --')):
            try:
                import json as _json
                _default_cfg = _json.loads(
                    (Path(parent_directory) / 'PlanetProfile' / 'Inference'
                     / 'configs' / 'test50_titan_noocean_andrade_8D.json'
                     ).read_text())
                _apply_config_to_session_state(_default_cfg)
                st.session_state['_default_cfg_applied'] = True
                st.rerun()
            except Exception as _e:
                st.warning(f"Could not auto-load the default configuration: {_e}")
                st.session_state['_default_cfg_applied'] = True

        st.markdown("---")

        # Parameter configuration
        render_parameter_config(PARAMETER_REGISTRY, CATEGORY_LABELS, CATEGORY_ORDER, validate_parameter_combination)

        st.markdown("---")

        # Prior configuration
        render_prior_config(PARAMETER_REGISTRY)

        st.markdown("---")

        # Observables
        render_observables_config()

        st.markdown("---")

        # Sampler settings
        render_sampler_settings()

        st.markdown("---")

        # Structure cache
        render_structure_cache_input()

        # Run button
        render_run_button(InferenceConfig, MCMCRunner)

        # Help button
        with st.expander("📚 About Bayesian Inference"):
            st.markdown("""
            **What is Bayesian Inference?**

            Bayesian inference uses observed data to constrain the posterior probability
            distribution of model parameters. Given prior knowledge and a likelihood
            function, we can compute:

            ```
            P(θ | data) ∝ P(data | θ) × P(θ)
            ```

            Where:
            - `θ` = interior parameters (viscosities, rheology constants)
            - `P(θ | data)` = posterior (what we want)
            - `P(data | θ)` = likelihood (forward model)
            - `P(θ)` = prior (initial assumptions)

            **MCMC (Markov Chain Monte Carlo):**
            - Samples from posterior via random walk
            - Robust, proven method for 5-10D spaces
            - Takes ~1-3 hours for typical convergence

            **Reference:** Petricca et al. (2025) *Nature* - Titan's tidal Love numbers
            challenge the ocean-world paradigm.
            """)

    with col_results:
        st.subheader("📊 Results")

        # Render results (or placeholder)
        render_results()


# ============================================================
# Entry Point
# ============================================================
if __name__ == "__main__":
    main()
