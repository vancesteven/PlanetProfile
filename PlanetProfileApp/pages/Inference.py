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
# UI Render Functions
# ============================================================

def render_preset_selector(PARAMETER_PRESETS):
    """Render preset configuration selector."""
    st.markdown("#### 📋 Configuration Preset")

    preset_options = {
        'andrade_titan': f"🪐 {PARAMETER_PRESETS['andrade_titan']['name']}",
        'andrade_titan_noocean_8D': f"🪐 {PARAMETER_PRESETS['andrade_titan_noocean_8D']['name']}",
        'maxwell_titan': f"🪐 {PARAMETER_PRESETS['maxwell_titan']['name']}",
        'andrade_europa': f"🌊 {PARAMETER_PRESETS['andrade_europa']['name']}",
        'custom': "⚙️ Custom Parameter Selection"
    }

    preset_choice = st.radio(
        "Select configuration:",
        options=list(preset_options.keys()),
        format_func=lambda x: preset_options[x],
        index=0 if st.session_state.inference_preset == 'andrade_titan' else
              (1 if st.session_state.inference_preset == 'andrade_titan_noocean_8D' else
               (2 if st.session_state.inference_preset == 'maxwell_titan' else
                (3 if st.session_state.inference_preset == 'andrade_europa' else 4))),
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
            st.session_state.inference_structure_cache_path = 'PlanetProfile/Test/mcmc_results/titan_allice_yao2014_structure_grid.pkl'
        elif preset_choice == 'maxwell_titan':
            st.session_state.inference_structure_cache_path = 'titan_cache/titan_maxwell_grid_cache.pkl'
        elif preset_choice == 'andrade_europa':
            clath_suffix = 'clath' if st.session_state.inference_use_clathrate else 'noclath'
            st.session_state.inference_structure_cache_path = f"europa_cache/europa_structure_{clath_suffix}.pkl"

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

    # Re(k2) input
    col1, col2 = st.columns(2)
    with col1:
        Re_k2_value = st.number_input(
            "Re(k₂) — Real part:",
            value=st.session_state.inference_observables['Re_k2'][0],
            format="%.4f",
            key='Re_k2_value'
        )
    with col2:
        Re_k2_uncertainty = st.number_input(
            "± Uncertainty:",
            value=st.session_state.inference_observables['Re_k2'][1],
            format="%.4f",
            key='Re_k2_unc'
        )

    # Im(k2) input - handle both old 'Im_k2' and new 'abs_Im_k2' keys for backward compatibility
    observables = st.session_state.inference_observables
    im_k2_default = observables.get('abs_Im_k2', observables.get('Im_k2', (0.135, 0.035)))

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

    # C/MR² (optional, checkbox-gated)
    cmr2_default = st.session_state.inference_observables.get('CMR2', (0.343, 0.001))
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

    # Update session state (use abs_Im_k2 as the canonical key)
    new_observables = {
        'Re_k2': (Re_k2_value, Re_k2_uncertainty),
        'abs_Im_k2': (Im_k2_value, Im_k2_uncertainty),
    }
    if cmr2_obs is not None:
        new_observables['CMR2'] = cmr2_obs
    st.session_state.inference_observables = new_observables

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
                _auto_gen_map = {
                    'andrade_titan': ('PlanetProfile.Test.PPTest41', 'titan_cache/', False, '1-3 minutes'),
                    'maxwell_titan': ('PlanetProfile.Test.PPTest42', 'titan_cache/', True, '15-30 minutes'),
                    'andrade_europa': ('PlanetProfile.Test.PPTest3', 'europa_cache/', False, '1-3 minutes'),
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
            if preset in ['andrade_titan', 'maxwell_titan']:
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

                status_placeholder.markdown(f"""
                **Current Status:**
                - Samples: {progress_dict['n_samples']}
                - ESS: {progress_dict['ess']:.0f} / {config.sampler_settings['n_effective']}
                - Acceptance: {f"{progress_dict['acceptance_rate']:.1%}" if progress_dict['acceptance_rate'] is not None else "N/A"}
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

    # Convergence metrics
    with st.expander("📈 Convergence Metrics", expanded=True):
        col1, col2, col3 = st.columns(3)

        metrics = result.convergence_metrics
        col1.metric("ESS", f"{metrics['ess']:.0f}")
        acc = metrics['acceptance_rate']
        col2.metric("Acceptance Rate", f"{acc:.1%}" if acc is not None else "N/A")
        col3.metric("R-hat", f"{metrics['r_hat']:.3f}")

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
                tb_norm = (tb_vals - tb_vals.min()) / (tb_vals.ptp() + 1e-10)
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

                # Detect fixed-structure case: all samples share the same CMR2
                # (occurs when no structural parameters like Tb_K are free).
                cmr2_range = cmr2_vals.max() - cmr2_vals.min()
                is_fixed = cmr2_range < 1e-8

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
                        "The structure is **fixed** for this model (no structural parameters "
                        "such as Tb_K are free). C/MR² is constant across all posterior "
                        "samples and acts as a uniform chi² offset — it constrains the "
                        "structure but does not discriminate between rheological models."
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
        heating_results = result.heating_results or []
        if not heating_results:
            st.info("Heating data unavailable. Set **n_reeval > 0** in sampler settings to compute per-phase heating.")
        else:
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
                                   colors=['#9c755f', '#76b7b2', '#4e79a7'],
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

                phase_colors = {'Ih': 'C0', 'III': 'C1', 'V': 'C2', 'VI': 'C3',
                                'Sil': 'C4', 'Clath': 'C5'}

                # Use only phases actually present in the data
                phases_to_show = [p for p in ['Ih', 'III', 'V', 'VI', 'Sil']
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
                for phase in ['Sil', 'V', 'VI', 'III', 'Ih']:
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

    with col_config:
        st.subheader("⚙️ Configuration")

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
