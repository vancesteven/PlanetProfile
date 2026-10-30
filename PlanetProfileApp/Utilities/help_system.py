"""
Help System for PlanetProfileApp
Contextual help and documentation for parameters
"""

# Help text for all parameters
HELP_TEXT = {
    # Bulk planetary properties
    "Planet.Bulk.R_m": {
        "short": "Mean radius of the planetary body",
        "long": """
        **Radius (m)**

        The mean radius of the planetary body in meters.

        **Typical ranges:**
        - Small moons (Enceladus): ~250 km
        - Medium moons (Europa): ~1,500 km
        - Large moons (Ganymede): ~2,600 km
        - Ice giants: ~25,000 km

        **What it affects:**
        - Total mass (with density)
        - Surface gravity
        - Pressure gradient
        - Layer thicknesses
        """,
        "units": "m"
    },

    "Planet.Bulk.M_kg": {
        "short": "Total mass of the planetary body",
        "long": """
        **Mass (kg)**

        The total mass of the planetary body in kilograms.

        **Typical ranges:**
        - Small moons: ~10^20 kg
        - Medium moons: ~10^22 kg
        - Large moons: ~10^23 kg

        **What it affects:**
        - Density (with radius)
        - Internal pressure
        - Gravitational compression
        - Moment of inertia

        **Tip:** Use the estimated density to check if mass and radius are consistent.
        """,
        "units": "kg"
    },

    "Planet.Bulk.Tsurf_K": {
        "short": "Temperature at the surface",
        "long": """
        **Surface Temperature (K)**

        The temperature at the planetary surface in Kelvin.

        **Typical ranges:**
        - Jupiter system: 90-130 K
        - Saturn system: 70-95 K
        - Uranus/Neptune: 50-75 K
        - Pluto/KBO: 40-60 K

        **What it affects:**
        - Ice shell thermal structure
        - Top boundary condition for heat flux
        - Phase of surface materials

        **Note:** This is usually well-constrained by observations.
        """,
        "units": "K"
    },

    "Planet.Bulk.Psurf_MPa": {
        "short": "Pressure at the surface",
        "long": """
        **Surface Pressure (MPa)**

        The pressure at the planetary surface in megapascals.

        **Typical values:**
        - Vacuum (most moons): ~0 MPa
        - Thin atmosphere (Triton): ~0.000014 MPa
        - Thick atmosphere (Titan): ~0.15 MPa

        **What it affects:**
        - Top boundary condition for pressure integration
        - Phase of surface volatiles
        - Very minor effect on interior structure
        """,
        "units": "MPa"
    },

    "Planet.Bulk.Cmeasured": {
        "short": "Normalized axial moment of inertia",
        "long": """
        **Moment of Inertia (C/MR²)**

        The normalized axial moment of inertia factor.

        **Physical meaning:**
        - Lower values (0.3-0.35) → denser core, differentiated
        - Higher values (0.35-0.4) → more uniform, less differentiated
        - Theoretical max for uniform sphere: 0.4

        **Typical values:**
        - Europa: ~0.346
        - Ganymede: ~0.311 (iron core)
        - Callisto: ~0.355 (less differentiated)

        **What it affects:**
        - Constraint on interior density distribution
        - Core size estimation
        - Degree of differentiation
        """,
        "units": "dimensionless"
    },

    "Planet.Bulk.Tb_K": {
        "short": "Temperature at ice-ocean interface",
        "long": """
        **Basal Temperature (K)**

        The temperature at the base of the ice shell / top of the ocean.

        **Typical ranges:**
        - Cold oceans: 250-265 K
        - Warm oceans: 265-273 K (melting point)
        - Depends on salinity and pressure

        **What it affects:**
        - Ice shell thickness
        - Ocean heat flux
        - Melting temperature at depth
        - Convection in ocean and ice

        **Tip:** Higher temperatures = thinner ice shells, more vigorous convection.
        """,
        "units": "K"
    },

    # Ocean properties
    "Planet.Ocean.wOcean_ppt": {
        "short": "Ocean salinity concentration",
        "long": """
        **Ocean Salinity (ppt)**

        Salt concentration in parts per thousand (g/kg).

        **Typical values:**
        - Pure water: 0 ppt
        - Low salinity: 1-10 ppt
        - Earth seawater: ~35 ppt
        - Concentrated: 50-150 ppt
        - Near saturation: 200-300 ppt (depends on salt type)

        **What it affects:**
        - Ocean density
        - Melting point depression
        - Electrical conductivity
        - Convection patterns
        - Magnetic induction signature

        **Tip:** Higher salinity → higher density and conductivity.
        """,
        "units": "ppt (g/kg)"
    },

    "Planet.Ocean.comp": {
        "short": "Ocean chemical composition",
        "long": """
        **Ocean Composition**

        The dominant dissolved salt in the ocean.

        **Options:**
        - **PureH2O**: Pure water, no dissolved salts
        - **NaCl**: Sodium chloride (table salt) - most common assumption
        - **MgSO4**: Magnesium sulfate (Epsom salt) - also common in icy moons
        - **NH3**: Ammonia - acts as antifreeze
        - **Seawater**: Earth-like seawater composition
        - **CustomSolution**: Define your own mix using Reaktoro

        **What it affects:**
        - Thermodynamic properties (density, heat capacity)
        - Electrical conductivity
        - Melting point depression
        - Phase diagram

        **Observations suggest:** NaCl or MgSO4 for most icy moons.
        """,
        "units": "composition type"
    },

    # Layer step settings
    "Planet.Steps.nIceI": {
        "short": "Number of grid points in ice I shell",
        "long": """
        **Ice I Steps**

        Number of computational grid points in the ice I (normal ice) shell.

        **Typical values:**
        - Fast calculation: 50-80
        - Standard: 100-150
        - High resolution: 200-300

        **What it affects:**
        - Resolution of temperature and pressure profiles
        - Calculation time (more steps = slower)
        - Accuracy of ice shell properties

        **Tip:** Start with 100, increase if you need more detail.
        """,
        "units": "grid points"
    },

    "Planet.Steps.nOcean": {
        "short": "Number of grid points in ocean layer",
        "long": """
        **Ocean Steps**

        Number of computational grid points in the liquid ocean.

        **Typical values:**
        - Fast: 80-100
        - Standard: 120-150
        - High resolution: 200-250

        **What it affects:**
        - Ocean temperature/pressure profile resolution
        - Convection calculations
        - Calculation time

        **Tip:** Ocean usually needs more points than ice for accuracy.
        """,
        "units": "grid points"
    },

    # General tips
    "simulation_tips": {
        "short": "Tips for running simulations",
        "long": """
        **Simulation Tips**

        1. **Start simple:** Begin with default values for your chosen body

        2. **Parameter sweep:** Change one parameter at a time to understand effects

        3. **Check convergence:** If results look strange, try more grid points

        4. **Physical constraints:** Make sure your parameters are physically reasonable
           - Mass and radius should give reasonable density
           - Temperature should be between surface and interior values
           - Salinity should be below saturation

        5. **Calculation time:** More grid points = more accuracy but slower
           - Development: 50-100 points per layer
           - Production: 150-200 points per layer
           - Publication: 200-300 points per layer

        6. **Compare to observations:** Check if your results match known constraints
           like moment of inertia, magnetic field measurements, etc.
        """
    }
}


def get_help(parameter_key, format='long'):
    """
    Get help text for a parameter

    Args:
        parameter_key: Key like "Planet.Bulk.R_m"
        format: 'short' or 'long'

    Returns:
        Help text string or None
    """
    help_dict = HELP_TEXT.get(parameter_key, None)
    if help_dict is None:
        return None

    return help_dict.get(format, help_dict.get('short', 'No help available'))


def get_units(parameter_key):
    """
    Get units for a parameter

    Args:
        parameter_key: Key like "Planet.Bulk.R_m"

    Returns:
        Units string or empty string
    """
    help_dict = HELP_TEXT.get(parameter_key, {})
    return help_dict.get('units', '')


def create_help_button(parameter_key, label="Help"):
    """
    Create a help button for a parameter (for use with Streamlit)

    Args:
        parameter_key: Key like "Planet.Bulk.R_m"
        label: Button label

    Returns:
        Help text to use with st.help parameter
    """
    return get_help(parameter_key, format='short')


def get_parameter_category(parameter_key):
    """
    Get the category of a parameter

    Args:
        parameter_key: Key like "Planet.Bulk.R_m"

    Returns:
        Category string ('bulk', 'ocean', 'core', 'steps', etc.)
    """
    if 'Bulk' in parameter_key:
        return 'bulk'
    elif 'Ocean' in parameter_key:
        return 'ocean'
    elif 'Core' in parameter_key or 'Sil' in parameter_key:
        return 'interior'
    elif 'Steps' in parameter_key:
        return 'computational'
    else:
        return 'other'


# FAQ content
FAQ = {
    "What affects ocean thickness?": """
    Ocean thickness is determined by several factors:

    1. **Heat flux from below** - More heat → thinner ice, thicker ocean
    2. **Surface temperature** - Colder surface → thicker ice
    3. **Ocean salinity** - Salt lowers freezing point → can prevent freezing
    4. **Internal heat sources** - Tidal heating, radioactive decay
    5. **Bulk composition** - Amount of water vs. rock

    The balance between heat loss and heat production determines where
    the ice-ocean boundary stabilizes.
    """,

    "Why does my simulation fail?": """
    Common reasons for simulation failures:

    1. **Unrealistic parameters** - Check that mass/radius give reasonable density
    2. **Too few grid points** - Try increasing nIceI, nOcean, etc.
    3. **Extreme conditions** - Very high/low temperatures or salinities
    4. **Numerical instability** - Try smaller steps or different parameters
    5. **Missing data** - Some ocean compositions need additional data files

    **Troubleshooting steps:**
    - Start with a default body configuration
    - Change one parameter at a time
    - Check the error message carefully
    - Increase grid resolution if results look jagged
    """,

    "How do I interpret the plots?": """
    **Main plot types:**

    1. **Hydrosphere plots** - Temperature, density, pressure vs depth in ice/ocean
       - Look for phase transitions (kinks in curves)
       - Check if ocean exists (liquid water layer)

    2. **Wedge diagram** - Cross-section showing all layers
       - Quick visual of interior structure
       - Layer thicknesses to scale

    3. **Gravity plots** - Internal gravity and pressure
       - Should be smooth curves
       - Discontinuities show layer boundaries

    4. **Seismic plots** - Sound speeds, elastic moduli
       - Used to compare with future seismic data

    5. **Magnetic induction** - Response to external fields
       - Depends on ocean conductivity and thickness
    """,

    "What's the difference between ocean compositions?": """
    **Common ocean types:**

    **NaCl (Sodium Chloride)**
    - Most commonly assumed for icy moons
    - Moderate conductivity
    - Similar to Earth's oceans
    - Evidence from spectroscopy on Europa

    **MgSO4 (Magnesium Sulfate)**
    - Also common assumption
    - Different thermal properties than NaCl
    - May match spectroscopic observations
    - Lower conductivity than NaCl

    **NH3 (Ammonia)**
    - Acts as antifreeze - lowers melting point significantly
    - Allows liquid at lower temperatures
    - Suggested for outer solar system bodies

    **PureH2O**
    - Simplest case
    - No dissolved salts
    - Higher freezing point

    **Seawater**
    - Earth-like composition
    - Mix of multiple salts
    - Good for comparative studies
    """
}


def get_faq_answer(question):
    """
    Get answer to FAQ question

    Args:
        question: FAQ question string

    Returns:
        Answer string or None
    """
    return FAQ.get(question, None)


def list_faq_questions():
    """
    Get list of all FAQ questions

    Returns:
        List of question strings
    """
    return list(FAQ.keys())
