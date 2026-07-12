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
import streamlit as st
import numpy as np
from pathlib import Path
import time
import traceback

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

    # Migrate stale session state: clear param_space if it contains old 'log10_zeta' key
    if ('inference_param_space' in st.session_state and
            'log10_zeta' in st.session_state.inference_param_space):
        st.session_state.inference_param_space = {}
        st.session_state.inference_selected_params = []

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

    if 'inference_observables' not in st.session_state:
        st.session_state.inference_observables = {
            'Re_k2': (0.608, 0.048),  # Petricca et al. 2025 defaults
            'abs_Im_k2': (0.135, 0.035),
            'CMR2': (0.343, 0.001),   # Petricca et al. 2025
        }

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

    NOTE: the existing ``render_preset_selector`` / ``render_observables_config`` /
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

    # Sync the preset radio WIDGET to 'custom'. The radio has key=
    # 'preset_radio', and Streamlit keyed-widget state wins over index= on
    # every rerun — without this, the radio keeps its old selection and its
    # non-custom branch immediately overwrites the config we just loaded
    # (params, observables, cache path). 'custom' is the one branch that
    # writes nothing back. inference_preset (set above) still carries the
    # matched preset for bodyname resolution in the run button.
    st.session_state.preset_radio = 'custom'

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

def render_preset_selector(PARAMETER_PRESETS):
    """Render preset configuration selector."""
    st.markdown("#### 📋 Configuration Preset")

    # NOTE (2026-07-12): the legacy 'andrade_europa' preset (5D PPTest46,
    # pure-water reference structure) is retired from this radio — its
    # parameter set cannot consume the seawater Tb-grid cache the Europa
    # campaign uses, and its auto-gen path built a wrong PureH2O cache.
    # Europa runs load configs/europa_seawater_andrade_7D.json via the
    # "Load Config File" selector above.
    preset_options = {
        'andrade_titan': f"🪐 {PARAMETER_PRESETS['andrade_titan']['name']}",
        'andrade_titan_noocean_8D': f"🪐 {PARAMETER_PRESETS['andrade_titan_noocean_8D']['name']}",
        'maxwell_titan': f"🪐 {PARAMETER_PRESETS['maxwell_titan']['name']}",
        'custom': "⚙️ Custom Parameter Selection"
    }

    _preset_index = {name: i for i, name in enumerate(preset_options)}
    preset_choice = st.radio(
        "Select configuration:",
        options=list(preset_options.keys()),
        format_func=lambda x: preset_options[x],
        index=_preset_index.get(st.session_state.inference_preset,
                                len(preset_options) - 1),
        key='preset_radio'
    )

    # Update session state
    if preset_choice != 'custom':
        st.session_state.inference_preset = preset_choice
        st.session_state.inference_custom_mode = False

        # Load preset parameters
        preset = PARAMETER_PRESETS[preset_choice]
        st.session_state.inference_selected_params = preset['parameters']

        # Filter param_space to only keep parameters in the preset
        # (prevents stale parameters from previous configurations)
        if st.session_state.inference_param_space:
            st.session_state.inference_param_space = {
                k: v for k, v in st.session_state.inference_param_space.items()
                if k in preset['parameters']
            }

        # Load preset observables
        st.session_state.inference_observables = preset['observables']

        # Auto-populate structure cache path based on preset
        if preset_choice == 'andrade_titan':
            clath_suffix = 'clath' if st.session_state.inference_use_clathrate else 'noclath'
            st.session_state.inference_structure_cache_path = f"titan_cache/titan_structure_{clath_suffix}.pkl"
        elif preset_choice == 'andrade_titan_noocean_8D':
            st.session_state.inference_structure_cache_path = 'PlanetProfile/Test/mcmc_results/Titan/Test50_andrade_noocean_yao2014/titan_allice_yao2014_structure_grid.pkl'
        elif preset_choice == 'maxwell_titan':
            st.session_state.inference_structure_cache_path = 'titan_cache/titan_maxwell_grid_cache.pkl'

        # Show preset description
        st.info(f"**Description:** {preset['description']}")
        st.caption(f"Rheology: {preset['rheology'].capitalize()} | Test: {preset['test_module']}")

    else:
        st.session_state.inference_custom_mode = True
        st.info("**Custom mode:** Select parameters manually below.")


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
    re_k2_default = _obs_now.get('Re_k2', (0.608, 0.048))
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
    im_k2_default = observables.get('abs_Im_k2', observables.get('Im_k2', (0.135, 0.035)))

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

    # Show validation errors
    if not validation['valid']:
        st.error("❌ **Configuration errors:**")
        for error in validation['errors']:
            st.error(f"- {error}")
        st.stop()

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
        'label': 'Titan (Andrade, no ocean) — Test50 8D',
        'bodyname': 'Titan',
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
    },
}


def _sbi_artifacts_dir():
    return Path(parent_directory) / 'PlanetProfile' / 'Inference' / 'sbi_artifacts'


@st.cache_resource(show_spinner="Loading SBI artifact...")
def _load_sbi_runner(artifact_path: str, _mtime: float):
    """Load a pretrained artifact into a sampling-only SBIRunner, cached per
    (path, mtime) so a redeployed artifact invalidates the cache."""
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    return SBIRunner.load_artifact(artifact_path)


def render_amortized_config():
    """Amortized-mode configuration form. Returns a run-spec dict or None."""
    available = []
    for fname, slot in _SBI_ARTIFACT_SLOTS.items():
        p = _sbi_artifacts_dir() / fname
        if p.exists():
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
    slot = _SBI_ARTIFACT_SLOTS[fname]
    artifact_path = _sbi_artifacts_dir() / fname

    runner = _load_sbi_runner(str(artifact_path), artifact_path.stat().st_mtime)
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
    with st.expander("✂️ Prior truncation (exact, within trained bounds)"):
        st.caption("Narrowing a uniform prior post hoc is exact (sample "
                   "filtering). Widening requires retraining.")
        for name, (lo, hi) in bounds.items():
            sel = st.slider(name, min_value=lo, max_value=hi,
                            value=(lo, hi), key=f'amort_trunc_{name}')
            if sel != (lo, hi):
                truncate_bounds[name] = (float(sel[0]), float(sel[1]))

    # --- Observables: values live, sigmas locked ---
    st.markdown("#### 🔭 Observables (condition the posterior)")
    x_obs = {}
    obs_sigma = {}
    for name in runner.obs_names:
        c1, c2 = st.columns(2)
        default = slot.get('default_obs', {}).get(name, 0.0)
        with c1:
            x_obs[name] = st.number_input(
                f"{name}:", value=float(default), format="%.4f",
                key=f'amort_obs_{name}')
        sig = float(train_sigma.get(name, float('nan')))
        obs_sigma[name] = sig
        with c2:
            st.number_input(
                "± σ (frozen at training):", value=sig, format="%.4f",
                key=f'amort_sigma_{name}', disabled=True,
                help="The flow was trained with this observation noise; "
                     "changing σ requires retraining (or MCMC mode).")

    # --- Sampler settings ---
    st.markdown("#### ⚙️ Sampling")
    n_post = st.select_slider("Posterior samples:",
                              options=[1000, 2000, 5000, 10000, 20000, 50000],
                              value=10000, key='amort_n_post')
    n_reeval = st.select_slider("Heating re-evaluation subset:",
                                options=[100, 250, 500, 1000, 2000],
                                value=500, key='amort_n_reeval')
    seed = st.number_input("Seed:", value=42, step=1, key='amort_seed')

    # --- Validated-domain guard (hard refusal, per deployment conditions) ---
    if slot.get('scope_note'):
        st.caption(f"ℹ️ {slot['scope_note']}")
    domain_violations = []
    for name, (lo, hi) in (slot.get('x_obs_limits') or {}).items():
        val = x_obs.get(name)
        if val is not None and not (lo <= abs(val) <= hi):
            domain_violations.append(
                f"{name} = {val:g} outside validated domain [{lo:g}, {hi:g}]")
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
        'seed': int(seed),
        'cache_path': cache_path,
    }


def render_model_assumptions(exec_mode: str):
    """Model build-up + assumptions text, shared by both execution modes."""
    with st.expander("📖 Model build-up & assumptions"):
        st.markdown("""
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
        from PlanetProfile.Inference.sbi_runner import SBIRunner

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
                truncate_bounds=spec['truncate_bounds'],
                progress_callback=progress_callback,
            )

        st.session_state.inference_results = result
        st.session_state.amort_last_run_fp = _amortized_spec_fingerprint(spec)
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

                n_dim = len(valid_cols)
                fig_size = max(10, 2.5 * n_dim)
                fig = plt.figure(figsize=(fig_size, fig_size))
                corner.corner(
                    corner_samples_plot,
                    labels=corner_labels_plot,
                    quantiles=[0.16, 0.5, 0.84],
                    show_titles=True,
                    title_fmt='.2f',
                    title_kwargs={'fontsize': 10},
                    color='steelblue',
                    fig=fig,
                )
                st.pyplot(fig)
                plt.close(fig)
        except ImportError:
            st.info("Install the `corner` library to view corner plots: `pip install corner`")
        except Exception as e:
            st.warning(f"Corner plot unavailable: {e}")

    # k2 scatter with error ellipse
    with st.expander("📡 k₂ Posterior Scatter", expanded=True):
        try:
            import matplotlib.pyplot as plt
            from matplotlib.patches import Ellipse

            Re_obs, Re_err = result.config.observables.get('Re_k2', (0.608, 0.048))
            Im_obs, Im_err = result.config.observables.get('abs_Im_k2', (0.135, 0.035))

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
                total = sum(h.values()) if h else 1e-30
                f_sil.append(h.get('Sil', 0) / total if total > 0 else 0)
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
                sc = ax.scatter(Re_arr, Im_arr, c=f_sil, cmap='RdYlBu_r',
                                s=pt_size, alpha=0.6, vmin=0, vmax=1)
                plt.colorbar(sc, ax=ax, label='Silicate heating fraction')
            else:
                ax.scatter(Re_arr, Im_arr, s=pt_size, alpha=0.6, color='steelblue')

            ax.add_patch(Ellipse((Re_obs, Im_obs), 2*Re_err, 2*Im_err,
                                 fill=False, color='red', linewidth=2,
                                 linestyle='--', label=r'1$\sigma$'))
            ax.add_patch(Ellipse((Re_obs, Im_obs), 4*Re_err, 4*Im_err,
                                 fill=False, color='red', linewidth=1,
                                 linestyle=':', label=r'2$\sigma$'))
            ax.set_xlabel(r'Re$(k_2)$')
            ax.set_ylabel(r'$|$Im$(k_2)|$')
            ax.set_title(r'$k_2$ Posterior')
            ax.legend()
            st.pyplot(fig)
            plt.close(fig)
            caption_bits = ["Red ellipses: 1σ / 2σ observational constraint."]
            if isinstance(pt_size, np.ndarray):
                caption_bits.append(
                    "**Point size** scales with the sampled basal temperature "
                    "T_b (larger = warmer ice shell base; the structural "
                    "dimension of the model).")
            if f_sil is not None and len(f_sil) == len(Re_arr):
                caption_bits.append(
                    "**Color** is the fraction of tidal heating dissipated "
                    "in the silicate interior (blue = ice-dominated, red = "
                    "rock-dominated dissipation).")
            st.caption(" ".join(caption_bits))
        except Exception as e:
            st.warning(f"k₂ scatter unavailable: {e}")

    # C/MR² posterior
    with st.expander("⚖️ C/MR² Moment-of-Inertia Posterior", expanded=True):
        cmr2_obs = result.config.observables.get('CMR2')
        cmr2_results = getattr(result, 'cmr2_results', None)

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
                            color='steelblue', label='Posterior')

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
                st.pyplot(fig)
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

                st.pyplot(fig_stack)
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
                        st.pyplot(fig_pie)
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
                        total = sum(h.values()) if h else 1e-30
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
                st.pyplot(fig)
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
                        'inference_cache_key', 'inference_error', 'inference_error_traceback']:
                if key in st.session_state:
                    del st.session_state[key]
            st.rerun()

        # Preset selector
        render_preset_selector(PARAMETER_PRESETS)

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
