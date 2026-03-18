"""
Utility helper functions for PlanetProfileApp
Provides validation, formatting, and common UI components
"""

import streamlit as st
import numpy as np
from datetime import datetime
import re

# Physical constants for validation
JUPITER_RADIUS_M = 69911e3
JUPITER_MASS_KG = 1.898e27
EARTH_RADIUS_M = 6.371e6
EARTH_MASS_KG = 5.972e24


def validate_radius(R_m, body_name=None):
    """
    Validate radius with helpful feedback

    Returns: (is_valid, message, severity)
    severity: 'error', 'warning', 'info', or None
    """
    if R_m < 100e3:  # Less than 100 km
        return False, "⚠️ Radius seems too small for a planetary body (< 100 km)", "error"
    elif R_m < 500e3:  # Less than 500 km
        return True, "💡 Small body - similar to Enceladus (252 km)", "info"
    elif R_m > JUPITER_RADIUS_M:
        return False, f"⚠️ Radius exceeds Jupiter ({JUPITER_RADIUS_M/1e3:.0f} km)", "warning"
    elif R_m > EARTH_RADIUS_M * 2:
        return True, f"⚠️ Large body (> 2× Earth radius)", "warning"
    else:
        return True, None, None


def validate_mass(M_kg, R_m=None, body_name=None):
    """
    Validate mass with helpful feedback

    If R_m provided, checks hydrostatic equilibrium
    """
    if M_kg < 1e20:  # Less than 100 billion billion kg
        return False, "⚠️ Mass seems too small for a planetary body", "error"
    elif M_kg > JUPITER_MASS_KG:
        return False, f"⚠️ Mass exceeds Jupiter ({JUPITER_MASS_KG:.2e} kg)", "warning"

    # Check density if radius provided
    if R_m is not None:
        volume = (4/3) * np.pi * R_m**3
        density = M_kg / volume

        if density < 500:  # kg/m³
            return True, f"⚠️ Very low density ({density:.0f} kg/m³) - mostly ice?", "warning"
        elif density > 8000:  # kg/m³
            return True, f"⚠️ Very high density ({density:.0f} kg/m³) - iron core?", "warning"
        elif density > 3000 and density < 4000:
            return True, f"✅ Density typical for rocky body ({density:.0f} kg/m³)", "info"
        else:
            return True, f"Density: {density:.0f} kg/m³", "info"

    return True, None, None


def suggest_mass_from_radius(R_m, body_type="icy"):
    """
    Suggest mass based on radius and body type

    Args:
        R_m: Radius in meters
        body_type: 'icy' (1500 kg/m³), 'mixed' (2500), or 'rocky' (3500)
    """
    densities = {
        'icy': 1500,
        'mixed': 2500,
        'rocky': 3500
    }

    density = densities.get(body_type, 2000)
    volume = (4/3) * np.pi * R_m**3
    return density * volume


def estimate_runtime(Planet_or_dict, do_settings=None):
    """
    Estimate runtime based on number of steps

    Args:
        Planet_or_dict: Planet object or dict with 'nIceI', 'nOcean', 'nSilicate', etc.
        do_settings: Planet.Do object with calculation flags (optional)

    Returns:
        String with estimated time (e.g., "~30s", "~2min")
    """
    base_time = 5  # Base overhead in seconds

    # Time per step varies by layer (seconds per step)
    step_times = {
        'nIceI': 0.01,
        'nOcean': 0.02,
        'nSilicate': 0.03,
        'nCore': 0.02
    }

    # Handle both Planet object and dict
    if hasattr(Planet_or_dict, 'Steps'):
        # It's a Planet object
        n_steps_dict = {
            'nIceI': getattr(Planet_or_dict.Steps, 'nIceI', 100) or 100,
            'nOcean': getattr(Planet_or_dict.Steps, 'nOcean', 100) or 100,
            'nSilicate': getattr(Planet_or_dict.Steps, 'nSilicate', 100) or 100,
            'nCore': getattr(Planet_or_dict.Steps, 'nCore', 50) or 50
        }
        if do_settings is None and hasattr(Planet_or_dict, 'Do'):
            do_settings = Planet_or_dict.Do
    else:
        # It's already a dict
        n_steps_dict = Planet_or_dict

    total_time = base_time
    for key, n_steps in n_steps_dict.items():
        if key in step_times and n_steps is not None:
            total_time += n_steps * step_times[key]

    # Additional time for optional calculations
    if do_settings:
        if getattr(do_settings, 'CALC_SEISMIC', False):
            total_time *= 1.5
        if getattr(do_settings, 'CALC_CONDUCT', False):
            total_time *= 1.3

    # Format output
    if total_time < 60:
        return f"~{int(total_time)}s"
    elif total_time < 3600:
        minutes = int(total_time / 60)
        return f"~{minutes}min"
    else:
        hours = int(total_time / 3600)
        return f"~{hours}h"


def format_scientific(value, precision=2):
    """Format large numbers in scientific notation"""
    if abs(value) >= 1e6 or abs(value) <= 1e-3:
        return f"{value:.{precision}e}"
    else:
        return f"{value:.{precision}f}"


def create_metric_card(label, value, delta=None, help_text=None):
    """
    Create a nice metric display card

    Args:
        label: Metric label
        value: Metric value (will be formatted)
        delta: Change indicator (e.g., "Custom", "+10%")
        help_text: Tooltip text
    """
    if isinstance(value, (int, float)):
        if abs(value) >= 1e6:
            value_str = format_scientific(value)
        else:
            value_str = f"{value:,.2f}"
    else:
        value_str = str(value)

    st.metric(label, value_str, delta=delta, help=help_text)


def show_parameter_info(param_name, default_val, current_val, units=""):
    """
    Show parameter with default and current values
    """
    col1, col2, col3 = st.columns([2, 1, 1])

    with col1:
        st.write(f"**{param_name}**")
    with col2:
        st.write(f"Default: {format_scientific(default_val)} {units}")
    with col3:
        changed = abs(current_val - default_val) > 1e-10
        if changed:
            st.write(f"✏️ {format_scientific(current_val)} {units}")
        else:
            st.write(f"{format_scientific(current_val)} {units}")


def create_progress_indicator(current_step=None, total_steps=None, step_names=None, steps_complete=None):
    """
    Create a progress bar showing completion status

    Args:
        current_step: Current step number (1-indexed), or
        total_steps: Total number of steps, or
        step_names: List of step names, or
        steps_complete: Dict with boolean values for each step (legacy format)
    """
    if steps_complete is not None:
        # Legacy format: dict with boolean values
        total_steps = len(steps_complete)
        completed = sum(steps_complete.values())
        progress = completed / total_steps if total_steps > 0 else 0

        st.progress(progress)

        # Show checklist
        col1, col2 = st.columns(2)
        for i, (step_name, is_complete) in enumerate(steps_complete.items()):
            icon = "✅" if is_complete else "⬜"
            with col1 if i % 2 == 0 else col2:
                st.write(f"{icon} {step_name}")

        return progress

    elif current_step is not None and total_steps is not None:
        # New format: current step out of total
        progress = (current_step - 1) / total_steps if total_steps > 0 else 0

        # Create progress bar
        st.progress(progress)

        # Show step indicator
        if step_names and len(step_names) == total_steps:
            # Show all steps with current highlighted
            cols = st.columns(total_steps)
            for i, name in enumerate(step_names):
                with cols[i]:
                    if i + 1 < current_step:
                        st.markdown(f"✅ **{name}**")
                    elif i + 1 == current_step:
                        st.markdown(f"🔵 **{name}**")
                    else:
                        st.markdown(f"⬜ {name}")
        else:
            st.caption(f"Step {current_step} of {total_steps}")

        return progress

    else:
        st.warning("⚠️ Invalid progress indicator arguments")
        return 0


def show_validation_message(is_valid, message, severity):
    """
    Display validation message with appropriate styling

    Args:
        is_valid: Boolean
        message: Message to display
        severity: 'error', 'warning', 'info', or None
    """
    if message is None:
        return

    if severity == 'error':
        st.error(message)
    elif severity == 'warning':
        st.warning(message)
    elif severity == 'info':
        st.info(message)


def create_collapsible_help(title, content):
    """
    Create a collapsible help section

    Args:
        title: Title of the help section
        content: Help text content (can be markdown)
    """
    with st.expander(f"❓ {title}"):
        st.markdown(content)


def create_comparison_columns(labels, show_diff=False):
    """
    Create columns for comparing configurations

    Args:
        labels: List of configuration names
        show_diff: Whether to show a difference column

    Returns:
        List of column objects
    """
    n_cols = len(labels) + (1 if show_diff else 0)
    cols = st.columns(n_cols)

    for i, label in enumerate(labels):
        with cols[i]:
            st.subheader(label)

    if show_diff:
        with cols[-1]:
            st.subheader("Δ Difference")

    return cols


def sanitize_filename(name):
    """
    Sanitize filename by replacing invalid characters

    Args:
        name: Original filename

    Returns:
        Sanitized filename safe for filesystem
    """
    # Replace invalid characters with underscores
    name = re.sub(r'[<>:"/\\|?*]', '_', name)
    # Replace spaces with underscores
    name = name.replace(' ', '_')
    # Remove leading/trailing underscores and dots
    name = name.strip('_.')
    # Ensure doesn't start with number
    if name and name[0].isdigit():
        name = 'run_' + name
    return name


def create_download_button(data, filename, label="📥 Download", file_type="text"):
    """
    Create a styled download button

    Args:
        data: Data to download (str, bytes, or file-like)
        filename: Filename for download
        label: Button label
        file_type: MIME type hint
    """
    mime_types = {
        "text": "text/plain",
        "json": "application/json",
        "csv": "text/csv",
        "pdf": "application/pdf"
    }

    mime = mime_types.get(file_type, "application/octet-stream")

    st.download_button(
        label=label,
        data=data,
        file_name=filename,
        mime=mime
    )


def format_timestamp(timestamp=None):
    """
    Format timestamp for display

    Args:
        timestamp: datetime object or Unix timestamp (None = now)

    Returns:
        Formatted string
    """
    if timestamp is None:
        dt = datetime.now()
    elif isinstance(timestamp, (int, float)):
        dt = datetime.fromtimestamp(timestamp)
    else:
        dt = timestamp

    return dt.strftime("%Y-%m-%d %H:%M:%S")


def create_sidebar_status(planet_name, changes_made=False, run_status=None):
    """
    Create a status indicator in the sidebar

    Args:
        planet_name: Name of selected planet
        changes_made: Whether user has made changes
        run_status: 'idle', 'running', 'complete', 'error'
    """
    with st.sidebar:
        st.markdown("---")
        st.subheader("Status")

        if planet_name and planet_name != "-- Select a Planet --":
            st.write(f"🌍 **Body:** {planet_name}")
        else:
            st.write("⚠️ No planet selected")

        if changes_made:
            st.write("✏️ **Modified**")
        else:
            st.write("✅ **Default config**")

        if run_status:
            status_icons = {
                'idle': '⏸️',
                'running': '🚀',
                'complete': '✅',
                'error': '❌'
            }
            icon = status_icons.get(run_status, '❓')
            st.write(f"{icon} **Run:** {run_status.title()}")


def show_tips_carousel(tips=None):
    """
    Show rotating tips for users

    Args:
        tips: Optional list of tip strings. If not provided, uses default tips.
    """
    if tips is None:
        tips = [
            "💡 **Tip:** Start with a similar body as your template",
            "💡 **Tip:** Check the validation messages for parameter guidance",
            "💡 **Tip:** Use presets for quick common configurations",
            "💡 **Tip:** Save your session to resume work later",
            "💡 **Tip:** Export your configuration to share with collaborators",
            "💡 **Tip:** Comparison mode lets you test multiple scenarios",
        ]

    # Rotate tips based on time
    import time
    tip_index = int(time.time() / 10) % len(tips)  # Change every 10 seconds

    st.info(tips[tip_index])
