Tutorials
=========

This section provides step-by-step tutorials for common PlanetProfile tasks.

Tutorial 1: Exploring Ocean Salinity Effects
---------------------------------------------

**Goal**: Understand how ocean salinity affects Europa's interior structure and magnetic induction.

**Time**: 15 minutes

**Prerequisites**: PlanetProfile installed, basic command line knowledge

Step 1: Create a Custom Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a new directory for this tutorial::

    cd ~/PlanetProfile_work
    mkdir Tutorial1_Salinity
    cd Tutorial1_Salinity

Copy the default Europa configuration::

    cp ../Europa/PPEuropa.py ./PPEuropa_Tutorial.py

Step 2: Modify the Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Edit ``PPEuropa_Tutorial.py`` to change the ocean salinity:

.. code-block:: python

    # Original
    Planet.Ocean.wOcean_ppt = 35.2  # Seawater salinity

    # Modified
    Planet.Ocean.wOcean_ppt = 100.0  # High salinity

Step 3: Run the Model
~~~~~~~~~~~~~~~~~~~~~~

Run PlanetProfile with your modified configuration::

    python -m PlanetProfile.Main ./PPEuropa_Tutorial.py

The model will:

1. Calculate interior structure with high-salinity ocean
2. Generate plots showing results
3. Save outputs to current directory

Step 4: Compare Results
~~~~~~~~~~~~~~~~~~~~~~~~

Run a second model with low salinity:

.. code-block:: python

    # Create PPEuropa_LowSalt.py
    Planet.Ocean.wOcean_ppt = 10.0  # Low salinity

Run it::

    python -m PlanetProfile.Main ./PPEuropa_LowSalt.py

Step 5: Analyze the Differences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Load both results in Python:

.. code-block:: python

    import pickle
    import matplotlib.pyplot as plt

    # Load high salinity model
    with open('PlanetProfile_Europa_Tutorial.pkl', 'rb') as f:
        High = pickle.load(f)

    # Load low salinity model
    with open('PlanetProfile_Europa_LowSalt.pkl', 'rb') as f:
        Low = pickle.load(f)

    # Compare ice shell thickness
    print(f"High salinity ice thickness: {High.zb_km:.2f} km")
    print(f"Low salinity ice thickness: {Low.zb_km:.2f} km")

    # Compare ocean conductivity
    print(f"High salinity conductivity: {High.Ocean.sigmaOcean_Sm:.2f} S/m")
    print(f"Low salinity conductivity: {Low.Ocean.sigmaOcean_Sm:.2f} S/m")

    # Compare magnetic induction amplitude
    if hasattr(High.Magnetic, 'Amp'):
        print(f"High salinity Amp: {High.Magnetic.Amp[0]:.3f}")
        print(f"Low salinity Amp: {Low.Magnetic.Amp[0]:.3f}")

**Expected Results**:

- Higher salinity → thinner ice shell (lower freezing point)
- Higher salinity → higher ocean conductivity
- Higher salinity → stronger induced magnetic field

Tutorial 2: Core Composition Study
-----------------------------------

**Goal**: Investigate how core composition (FeS fraction) affects Ganymede's interior.

**Time**: 20 minutes

Step 1: Setup
~~~~~~~~~~~~~

Create a working directory::

    mkdir Tutorial2_Core
    cd Tutorial2_Core
    cp ../Ganymede/PPGanymede.py ./PPGanymede_Core.py

Step 2: Vary Core FeS Fraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Modify ``PPGanymede_Core.py``:

.. code-block:: python

    # Enable iron core
    Planet.Do.Fe_CORE = True

    # Set FeS fraction (0 = pure Fe, higher values = more sulfur)
    Planet.Core.xFeS = 0.10  # 10% FeS by mole fraction

Step 3: Run Multiple Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a script to run multiple core compositions:

.. code-block:: python

    # run_core_suite.py
    import os
    from PlanetProfile.Main import RunPPfile

    xFeS_values = [0.0, 0.05, 0.10, 0.15, 0.20]

    for xFeS in xFeS_values:
        # Create configuration with this xFeS
        config_text = f"""
from PlanetProfile.Default.Ganymede.PPGanymede import *

Planet.Core.xFeS = {xFeS}
        """

        filename = f'PPGanymede_xFeS{xFeS:.2f}.py'
        with open(filename, 'w') as f:
            f.write(config_text)

        # Run model
        print(f"Running xFeS = {xFeS}")
        RunPPfile(filename)

Run the script::

    python run_core_suite.py

Step 4: Analyze Results
~~~~~~~~~~~~~~~~~~~~~~~~

Compare core properties:

.. code-block:: python

    import pickle
    import glob
    import numpy as np
    import matplotlib.pyplot as plt

    # Load all results
    results = []
    xFeS_values = []
    for pkl_file in sorted(glob.glob('PlanetProfile_Ganymede_xFeS*.pkl')):
        with open(pkl_file, 'rb') as f:
            Planet = pickle.load(f)
            results.append(Planet)
            xFeS_values.append(Planet.Core.xFeS)

    # Extract properties
    core_radii = [P.Core.Rcore_m/1e3 for P in results]  # km
    core_densities = [P.Core.rhoCore_kgm3 for P in results]  # kg/m³

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(xFeS_values, core_radii, 'o-')
    ax1.set_xlabel('FeS Fraction')
    ax1.set_ylabel('Core Radius (km)')
    ax1.set_title('Core Radius vs. Composition')
    ax1.grid(True)

    ax2.plot(xFeS_values, core_densities, 's-', color='orange')
    ax2.set_xlabel('FeS Fraction')
    ax2.set_ylabel('Core Density (kg/m³)')
    ax2.set_title('Core Density vs. Composition')
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig('core_comparison.png', dpi=300)
    plt.show()

**Expected Results**:

- Higher FeS → larger core radius (FeS is less dense than Fe)
- Higher FeS → lower core density
- Total core mass constrained by body mass and moment of inertia

Tutorial 3: Parameter Space Exploration (Exploreogram)
-------------------------------------------------------

**Goal**: Use exploreograms to visualize how two parameters affect model outputs.

**Time**: 30 minutes

What is an Exploreogram?
~~~~~~~~~~~~~~~~~~~~~~~~~

An exploreogram runs PlanetProfile many times across a 2D grid of parameter values, creating a heatmap of output properties. This is useful for:

- Sensitivity studies
- Parameter trade-off analysis
- Identifying valid model regions

Step 1: Using PlanetProfileApp (GUI)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The easiest way to run exploreograms is through the GUI:

1. Launch PlanetProfileApp::

    cd PlanetProfileApp
    streamlit run PlanetProfileApp.py

2. Navigate to the **Exploreogram** page

3. Configure your exploration:

   - **Planet**: Select Europa
   - **X-Axis**: Ocean Salinity (0-60 ppt)
   - **Y-Axis**: Ice Bottom Temperature (265-280 K)
   - **Grid**: 10×10 (100 models)
   - **Color Variable**: Ice Shell Thickness (zb_km)

4. Select plot type:

   - **Interactive (Plotly)**: For exploration, zooming, hovering
   - **Static (Matplotlib)**: For publication figures

5. Click **Run Exploreogram**

6. Results display as an interactive heatmap showing how ice thickness varies with salinity and temperature

Step 2: Using Command Line
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For more control, use the command line interface.

Create an exploration configuration ``PPEuropaExplore.py``:

.. code-block:: python

    from PlanetProfile.Default.Europa.PPEuropa import *

    # Exploration settings
    Params.Explore.xName = 'wOcean_ppt'  # X-axis parameter
    Params.Explore.yName = 'Tb_K'         # Y-axis parameter
    Params.Explore.zName = 'zb_km'        # Output to plot

    # Parameter ranges
    Params.Explore.xRange = [0, 60]  # Salinity from 0-60 ppt
    Params.Explore.yRange = [265, 280]  # Temperature from 265-280 K

    # Grid resolution
    Params.Explore.nx = 10  # 10 points on X-axis
    Params.Explore.ny = 10  # 10 points on Y-axis

    # Enable exploreogram mode
    Params.DO_EXPLOREOGRAM = True

Run it::

    python -m PlanetProfile.Main Europa exploreogram

Results are saved to ``output/exploreograms/Europa/``.

Step 3: Magnetic Induction Exploreogram
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Explore how magnetic field amplitude varies with ocean properties.

In the GUI:

1. **X-Axis**: Ocean Salinity (10-50 ppt)
2. **Y-Axis**: Ocean Depth (80-120 km)
3. **Color Variable**: Induced Magnetic Field Amplitude (Amp_nT)
4. **Magnetic Induction**:

   - Check "Enable magnetic induction"
   - Select frequencies: ☑ Orbital, ☑ Synodic

5. **Grid**: 8×8 (64 models, ~10 minutes)

The heatmap shows:

- Regions where induction is strongest (high conductivity + thick ocean)
- Trade-offs between salinity and ocean thickness
- Observable vs. unobservable magnetic signatures

Step 4: Analyzing Exploreogram Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Load and analyze saved results:

.. code-block:: python

    import pickle
    import numpy as np
    import matplotlib.pyplot as plt

    # Load exploreogram results
    with open('output/exploreograms/Europa/ExplorationResults_Europa.pkl', 'rb') as f:
        Exploration = pickle.load(f)

    # Extract data
    xData = Exploration.xData  # Salinity values
    yData = Exploration.yData  # Temperature values
    zData = Exploration.base.zb_km  # Ice thickness grid

    # Find minimum ice thickness
    valid_mask = Exploration.base.VALID
    zData_valid = zData.copy()
    zData_valid[~valid_mask] = np.nan

    min_idx = np.nanargmin(zData_valid)
    min_row, min_col = np.unravel_index(min_idx, zData.shape)

    print(f"Minimum ice thickness: {zData[min_row, min_col]:.2f} km")
    print(f"At salinity: {xData[min_col]:.1f} ppt")
    print(f"At temperature: {yData[min_row]:.1f} K")

    # Plot contours
    fig, ax = plt.subplots(figsize=(10, 8))

    contour = ax.contourf(xData, yData, zData, levels=20, cmap='viridis')
    ax.contour(xData, yData, zData, levels=10, colors='white', linewidths=0.5, alpha=0.3)

    ax.set_xlabel('Ocean Salinity (ppt)')
    ax.set_ylabel('Ice Bottom Temperature (K)')
    ax.set_title('Europa Ice Shell Thickness')

    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('Ice Thickness (km)')

    # Mark minimum
    ax.plot(xData[min_col], yData[min_row], 'r*', markersize=20, label='Minimum')
    ax.legend()

    plt.tight_layout()
    plt.savefig('exploreogram_analysis.png', dpi=300)
    plt.show()

Tutorial 4: Custom Ocean Composition
-------------------------------------

**Goal**: Model an ocean with custom salt composition using Reaktoro.

**Time**: 25 minutes

**Prerequisites**: Reaktoro installed

Step 1: Understanding Reaktoro
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reaktoro is a geochemical modeling tool that computes thermodynamic properties for arbitrary aqueous solutions. PlanetProfile uses it for custom ocean compositions beyond the standard presets (Seawater, MgSO4, NaCl).

Step 2: Define Custom Solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create configuration ``PPEuropa_CustomOcean.py``:

.. code-block:: python

    from PlanetProfile.Default.Europa.PPEuropa import *

    # Use Reaktoro for ocean composition
    Planet.Ocean.comp = 'CustomSolution'

    # Define ocean composition (mole fractions or molalities)
    Planet.Ocean.mixingRatioToH2O = {
        'Na+': 0.5,    # Sodium
        'Cl-': 0.5,    # Chloride
        'Mg+2': 0.1,   # Magnesium
        'SO4-2': 0.1,  # Sulfate
        'K+': 0.05,    # Potassium
    }

    # Reaktoro database
    Planet.Ocean.database = 'Supcrt98'

    # Enable custom solution calculations
    Params.CustomSolution.CALC_CUSTOM_SOLUTION = True

Step 3: Run the Model
~~~~~~~~~~~~~~~~~~~~~~

::

    python -m PlanetProfile.Main ./PPEuropa_CustomOcean.py

Reaktoro will compute:

- Freezing temperature
- Density
- Heat capacity
- Thermal conductivity
- Electrical conductivity

Step 4: Compare to Standard Seawater
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import pickle

    # Load custom ocean model
    with open('PlanetProfile_Europa_CustomOcean.pkl', 'rb') as f:
        Custom = pickle.load(f)

    # Load standard seawater model
    with open('../Europa/PlanetProfile_Europa_default.pkl', 'rb') as f:
        Seawater = pickle.load(f)

    # Compare properties
    print("\nOcean Properties Comparison:")
    print(f"{'Property':<30} {'Custom':<15} {'Seawater':<15}")
    print("-" * 60)
    print(f"{'Ice thickness (km)':<30} {Custom.zb_km:<15.2f} {Seawater.zb_km:<15.2f}")
    print(f"{'Ocean conductivity (S/m)':<30} {Custom.Ocean.sigmaOcean_Sm:<15.3f} {Seawater.Ocean.sigmaOcean_Sm:<15.3f}")
    print(f"{'Mean ocean T (K)':<30} {Custom.Ocean.Tmean_K:<15.2f} {Seawater.Ocean.Tmean_K:<15.2f}")

Tutorial 5: Tidal Love Numbers
-------------------------------

**Goal**: Calculate tidal Love numbers to compare with observations.

**Time**: 20 minutes

**Prerequisites**: PyALMA3 installed

Step 1: Enable Gravity Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Love numbers describe tidal deformation and are observable via spacecraft measurements.

Modify ``PPEuropa.py``:

.. code-block:: python

    # Enable gravity calculations
    Params.CALC_GRAVITY = True

    # Configure PyALMA3 settings
    Params.Gravity.harmonic_degrees = [2]  # Degree-2 Love numbers (k2, h2)
    Params.Gravity.time_log_kyrs = [-3, 3]  # Time range in log10(kyrs)

Step 2: Set Rheological Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Love numbers depend on material rheology:

.. code-block:: python

    # Ice rheology
    Planet.Ocean.VISCOUS_ICE = True
    Planet.Ocean.iceViscosity_kgms = 1e14  # Ice viscosity

    # Mantle rheology
    Planet.Sil.mantleViscosity_kgms = 1e19  # Mantle viscosity

Step 3: Run Model
~~~~~~~~~~~~~~~~~

::

    python -m PlanetProfile.Main Europa

Step 4: Analyze Love Numbers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import pickle
    import matplotlib.pyplot as plt
    import numpy as np

    with open('Europa/PlanetProfile_Europa_default.pkl', 'rb') as f:
        Planet = pickle.load(f)

    # Extract Love numbers
    k2 = Planet.Gravity.k  # k2 Love number
    h2 = Planet.Gravity.h  # h2 Love number
    times = Planet.Gravity.time_log_kyrs  # Time array

    # Plot complex Love number
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Amplitude
    ax1.semilogx(10**times, np.abs(k2[0, :]), label='k2')
    ax1.semilogx(10**times, np.abs(h2[0, :]), label='h2')
    ax1.set_xlabel('Period (kyrs)')
    ax1.set_ylabel('Amplitude')
    ax1.set_title('Love Number Amplitude')
    ax1.legend()
    ax1.grid(True)

    # Phase
    ax2.semilogx(10**times, np.angle(k2[0, :], deg=True), label='k2')
    ax2.semilogx(10**times, np.angle(h2[0, :], deg=True), label='h2')
    ax2.set_xlabel('Period (kyrs)')
    ax2.set_ylabel('Phase (degrees)')
    ax2.set_title('Love Number Phase')
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig('love_numbers.png', dpi=300)
    plt.show()

    # Compare to observations
    # Europa's observed k2 ~ 0.015 (from Galileo)
    k2_orbital = np.abs(k2[0, np.argmin(np.abs(times - np.log10(3.5)))])  # ~3.5 day period
    print(f"Model k2 at Europa's orbital period: {k2_orbital:.4f}")
    print(f"Observed k2: ~0.015")

Next Steps
----------

Explore more advanced topics:

- :doc:`configuration`: Complete parameter reference
- :doc:`exploreogram`: Advanced parameter space exploration
- :doc:`api_reference`: Python API for scripting
- :doc:`contributing`: Contribute to PlanetProfile development

Need Help?
----------

If you encounter issues with these tutorials:

1. Check the **Troubleshooting** section in :doc:`getting_started`
2. Review the **Examples** section for working code
3. Post an issue on GitHub: https://github.com/vancesteven/PlanetProfile/issues
4. Email: steven.d.vance@jpl.nasa.gov
