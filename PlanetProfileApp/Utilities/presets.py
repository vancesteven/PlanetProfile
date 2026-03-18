"""
Preset Configurations for PlanetProfileApp
Common planetary body scenarios for quick setup
"""

# Preset configurations organized by body and scenario
PRESETS = {
    "Europa": {
        "Thin Ice Shell": {
            "description": "Europa with ~10-20 km ice shell, warm ocean",
            "settings": {
                "Planet.Bulk.Tb_K": 273.0,
                "Planet.Ocean.wOcean_ppt": 100.0,
                "Planet.Ocean.comp": "NaCl",
                "Planet.Steps.nIceI": 100,
                "Planet.Steps.nOcean": 150
            }
        },
        "Thick Ice Shell": {
            "description": "Europa with ~30-40 km ice shell, cold ocean",
            "settings": {
                "Planet.Bulk.Tb_K": 265.0,
                "Planet.Ocean.wOcean_ppt": 50.0,
                "Planet.Ocean.comp": "NaCl",
                "Planet.Steps.nIceI": 150,
                "Planet.Steps.nOcean": 120
            }
        },
        "MgSO4 Ocean": {
            "description": "Europa with magnesium sulfate ocean",
            "settings": {
                "Planet.Bulk.Tb_K": 270.0,
                "Planet.Ocean.wOcean_ppt": 100.0,
                "Planet.Ocean.comp": "MgSO4",
                "Planet.Steps.nIceI": 120,
                "Planet.Steps.nOcean": 140
            }
        }
    },
    "Ganymede": {
        "Standard": {
            "description": "Ganymede with multi-layer ice structure",
            "settings": {
                "Planet.Bulk.Tb_K": 255.0,
                "Planet.Ocean.wOcean_ppt": 35.0,
                "Planet.Ocean.comp": "Seawater",
                "Planet.Steps.nIceI": 100,
                "Planet.Steps.nOcean": 100,
                "Planet.Steps.nIceIII": 80,
                "Planet.Steps.nIceV": 50
            }
        },
        "Warm Deep Ocean": {
            "description": "Ganymede with warmer conditions at depth",
            "settings": {
                "Planet.Bulk.Tb_K": 265.0,
                "Planet.Ocean.wOcean_ppt": 10.0,
                "Planet.Ocean.comp": "PureH2O",
                "Planet.Steps.nIceI": 100,
                "Planet.Steps.nOcean": 150,
                "Planet.Steps.nIceIII": 100
            }
        }
    },
    "Callisto": {
        "Subsurface Ocean": {
            "description": "Callisto with deep subsurface ocean",
            "settings": {
                "Planet.Bulk.Tb_K": 260.0,
                "Planet.Ocean.wOcean_ppt": 35.0,
                "Planet.Ocean.comp": "Seawater",
                "Planet.Steps.nIceI": 150,
                "Planet.Steps.nOcean": 100
            }
        }
    },
    "Enceladus": {
        "Active Plumes": {
            "description": "Enceladus with thin ice, active ocean",
            "settings": {
                "Planet.Bulk.Tb_K": 273.0,
                "Planet.Ocean.wOcean_ppt": 20.0,
                "Planet.Ocean.comp": "NaCl",
                "Planet.Steps.nIceI": 80,
                "Planet.Steps.nOcean": 120
            }
        },
        "Regional Ocean": {
            "description": "Enceladus with localized ocean",
            "settings": {
                "Planet.Bulk.Tb_K": 268.0,
                "Planet.Ocean.wOcean_ppt": 10.0,
                "Planet.Ocean.comp": "NaCl",
                "Planet.Steps.nIceI": 60,
                "Planet.Steps.nOcean": 80
            }
        }
    },
    "Titan": {
        "Ammonia Ocean": {
            "description": "Titan with NH3-H2O ocean",
            "settings": {
                "Planet.Bulk.Tb_K": 250.0,
                "Planet.Ocean.wOcean_ppt": 100.0,
                "Planet.Ocean.comp": "NH3",
                "Planet.Steps.nIceI": 100,
                "Planet.Steps.nOcean": 100
            }
        },
        "Pure Water": {
            "description": "Titan with pure water ocean",
            "settings": {
                "Planet.Bulk.Tb_K": 260.0,
                "Planet.Ocean.wOcean_ppt": 0.0,
                "Planet.Ocean.comp": "PureH2O",
                "Planet.Steps.nIceI": 120,
                "Planet.Steps.nOcean": 120
            }
        }
    },
    "Triton": {
        "Standard": {
            "description": "Triton with subsurface ocean",
            "settings": {
                "Planet.Bulk.Tb_K": 255.0,
                "Planet.Ocean.wOcean_ppt": 100.0,
                "Planet.Ocean.comp": "NH3",
                "Planet.Steps.nIceI": 100,
                "Planet.Steps.nOcean": 100
            }
        }
    },
    "Pluto": {
        "Deep Ocean": {
            "description": "Pluto with deep subsurface ocean",
            "settings": {
                "Planet.Bulk.Tb_K": 250.0,
                "Planet.Ocean.wOcean_ppt": 150.0,
                "Planet.Ocean.comp": "NH3",
                "Planet.Steps.nIceI": 150,
                "Planet.Steps.nOcean": 100
            }
        }
    }
}


def get_presets_for_body(body_name):
    """
    Get all presets for a specific body

    Args:
        body_name: Name of the planetary body

    Returns:
        Dict of preset configurations or empty dict
    """
    return PRESETS.get(body_name, {})


def get_preset(body_name, preset_name):
    """
    Get a specific preset configuration

    Args:
        body_name: Name of the planetary body
        preset_name: Name of the preset

    Returns:
        Dict with 'description' and 'settings' or None
    """
    body_presets = PRESETS.get(body_name, {})
    return body_presets.get(preset_name, None)


def apply_preset(preset_config, session_state=None):
    """
    Apply preset settings to session state

    Args:
        preset_config: Preset config dict with 'settings' key, or direct settings dict
        session_state: Streamlit session_state object (uses st.session_state if not provided)

    Returns:
        List of changed keys
    """
    import streamlit as st

    if session_state is None:
        session_state = st.session_state

    # Handle both preset_config format and direct settings dict
    if isinstance(preset_config, dict) and 'settings' in preset_config:
        preset_settings = preset_config['settings']
    else:
        preset_settings = preset_config

    changed = []

    for key, value in preset_settings.items():
        session_state[key] = value
        changed.append(key)

        # Mark as changed for tracking
        if 'Bulk' in key:
            if 'changed_bulk_settings_flags' not in session_state:
                session_state['changed_bulk_settings_flags'] = {}
            if 'changed_bulk_settings' not in session_state:
                session_state['changed_bulk_settings'] = {}
            session_state['changed_bulk_settings_flags'][key] = True
            session_state['changed_bulk_settings'][key] = value
        elif 'Ocean' in key:
            session_state['custom_ocean_flag'] = True
        elif 'Steps' in key:
            if 'changed_step_settings_flags' not in session_state:
                session_state['changed_step_settings_flags'] = {}
            if 'changed_step_settings' not in session_state:
                session_state['changed_step_settings'] = {}
            session_state['changed_step_settings_flags'][key] = True
            session_state['changed_step_settings'][key] = value

    return changed


def list_all_bodies_with_presets():
    """
    Get list of all bodies that have presets

    Returns:
        List of body names
    """
    return list(PRESETS.keys())


def get_preset_summary(preset_config=None, body_name=None, preset_name=None):
    """
    Get a summary of a preset for display

    Args:
        preset_config: Preset configuration dict (direct), or
        body_name: Name of the planetary body (to lookup), and
        preset_name: Name of the preset (to lookup)

    Returns:
        Formatted string summary
    """
    # Allow passing preset directly or looking it up
    if preset_config is None:
        if body_name is not None and preset_name is not None:
            preset = get_preset(body_name, preset_name)
        else:
            return "Invalid arguments: provide either preset_config or both body_name and preset_name"
    else:
        preset = preset_config

    if not preset:
        return "Preset not found"

    lines = [preset['description'], ""]
    lines.append("**Settings:**")

    for key, value in preset['settings'].items():
        param_name = key.split('.')[-1]
        if isinstance(value, float):
            lines.append(f"  - {param_name}: {value:.2f}")
        else:
            lines.append(f"  - {param_name}: {value}")

    return "\n".join(lines)
