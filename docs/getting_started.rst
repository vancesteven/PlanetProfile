Getting Started with PlanetProfile
====================================

This guide will help you get started with PlanetProfile, from installation through running your first model.

Installation
------------

Prerequisites
~~~~~~~~~~~~~

PlanetProfile requires Python 3.8 or higher. Python 3.11 is recommended for users, and **required** for developers.

We strongly recommend using conda or mamba to manage your Python environment::

    conda create -n planetprofile python=3.11
    conda activate planetprofile

Required Dependencies
~~~~~~~~~~~~~~~~~~~~~

Install core dependencies with conda/mamba::

    conda install numpy=1.26.4 scipy matplotlib mpmath pandas
    conda install -c conda-forge gsw obspy spiceypy cmasher reaktoro

Install additional Python packages with pip::

    pip install SeaFreeze hdf5storage PyALMA3 TidalPy

Installing PlanetProfile
~~~~~~~~~~~~~~~~~~~~~~~~~

For Users
^^^^^^^^^

1. Install PlanetProfile via pip::

    pip install PlanetProfile

2. Create a working directory for your models::

    mkdir ~/PlanetProfile_work
    cd ~/PlanetProfile_work

3. Run the installation script to download configuration files and EOS data::

    python -m PlanetProfile.install

This will:

- Copy default configuration files to your working directory
- Download Perple_X equation of state data (~164 MB)
- Set up directory structure for each supported body

For Developers
^^^^^^^^^^^^^^

1. Clone the repository::

    git clone https://github.com/vancesteven/PlanetProfile.git
    cd PlanetProfile

2. Run the installation script::

    python -m PlanetProfile.install PPinstall

3. (Optional) Run the test suite to verify installation::

    python -m PlanetProfile.BuildTest

Your First Model
----------------

Running a Basic Model
~~~~~~~~~~~~~~~~~~~~~

Let's run a model for Europa, one of Jupiter's icy moons.

Command Line Interface
^^^^^^^^^^^^^^^^^^^^^^

The simplest way to run PlanetProfile::

    python -m PlanetProfile.Main Europa

This will:

1. Load the default configuration for Europa
2. Calculate the interior structure model
3. Generate plots showing the results
4. Save output files in the ``Europa/`` directory

Output files include:

- ``Europa/figures/``: PDF plots of the model
- ``Europa/PlanetProfile_Europa_*.txt``: Text file with model results
- ``Europa/PlanetProfile_Europa_*.pkl``: Python pickle file for reloading
- ``Europa/PlanetProfile_Europa_*.mat``: MATLAB-compatible output

As a Python Module
^^^^^^^^^^^^^^^^^^

You can also run PlanetProfile from within Python::

    from PlanetProfile.Main import RunPPfile

    # Run with default configuration
    RunPPfile('Europa', 'PPEuropa.py')

Using PlanetProfileApp (GUI)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PlanetProfile includes a graphical interface built with Streamlit, located in the ``PlanetProfileApp/`` directory of the repository.

1. Install Streamlit and additional dependencies::

    pip install streamlit pdf2image Pillow pandas
    conda install poppler  # macOS
    # For Windows, download Poppler separately

2. Launch the GUI::

    cd PlanetProfileApp
    streamlit run PlanetProfileApp.py

3. Your browser will open to the PlanetProfileApp interface

4. Follow the guided workflow:

   - **About**: Introduction and documentation
   - **Main Settings**: Select a planet or create custom body
   - **Bulk Planetary Settings**: Set mass, radius, etc.
   - **Ocean Settings**: Configure ocean composition
   - **Core and Silicate Settings**: Set interior properties
   - **Layer Step Settings**: Control resolution
   - **Run PlanetProfile**: Execute the model
   - **Outputs**: View results and figures
   - **Exploreogram**: Parameter space exploration

Understanding the Output
~~~~~~~~~~~~~~~~~~~~~~~~~

After running a model, you'll see several output files:

Text File (``.txt``)
^^^^^^^^^^^^^^^^^^^^

Contains tabulated model results::

    # Planet: Europa
    # Ocean composition: Seawater
    # Ocean salinity: 35.2 ppt
    # ...
    # Layer results:
    # r_m  P_MPa  T_K  rho_kgm3  phase  ...

Key sections:

- **Header**: Model configuration and parameters
- **Layer-by-layer**: Properties at each radial step
- **Summary**: Key results (ice thickness, ocean properties, etc.)

Figures
^^^^^^^

Default plots include:

1. **PT Plot**: Pressure-temperature profile showing phase transitions
2. **Profile Plot**: Radial profiles of properties (density, temperature, etc.)
3. **Seismic Plot**: P- and S-wave velocities
4. **Conductivity Plot**: Electrical conductivity profile
5. **Magnetic Plot**: Induced magnetic field (if calculated)

Pickle File (``.pkl``)
^^^^^^^^^^^^^^^^^^^^^^^

Python object containing the complete ``Planet`` structure. Load it with::

    import pickle
    with open('Europa/PlanetProfile_Europa_default.pkl', 'rb') as f:
        Planet = pickle.load(f)

    # Access model results
    print(Planet.zb_km)  # Ice shell thickness in km
    print(Planet.Ocean.comp)  # Ocean composition
    print(Planet.r_m)  # Radial array

Customizing Your Model
----------------------

Models are controlled by configuration files in the body directory.

Configuration File Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each body has a configuration file ``Body/PPBody.py``. For example, ``Europa/PPEuropa.py``::

    # Europa configuration
    from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

    Planet = PlanetStruct('Europa')

    # Bulk properties
    Planet.Bulk.R_m = 1560.8e3  # Radius in meters
    Planet.Bulk.M_kg = 4.7998e22  # Mass in kg
    Planet.Bulk.Tb_K = 269.8  # Basal temperature

    # Ocean composition
    Planet.Ocean.comp = 'Seawater'
    Planet.Ocean.wOcean_ppt = 35.2  # Salinity in ppt

    # Core properties
    Planet.Do.Fe_CORE = True
    Planet.Core.xFeS = 0.0  # FeS fraction

Key Parameters
~~~~~~~~~~~~~~

Bulk Properties
^^^^^^^^^^^^^^^

- ``Planet.Bulk.R_m``: Body radius (m)
- ``Planet.Bulk.M_kg``: Body mass (kg)
- ``Planet.Bulk.Tb_K``: Ice-ocean interface temperature (K)
- ``Planet.Bulk.Tsurf_K``: Surface temperature (K)
- ``Planet.Bulk.Psurf_MPa``: Surface pressure (MPa)

Ocean Properties
^^^^^^^^^^^^^^^^

- ``Planet.Ocean.comp``: Composition ('Seawater', 'MgSO4', 'NaCl', 'PureH2O', etc.)
- ``Planet.Ocean.wOcean_ppt``: Salinity (parts per thousand)
- ``Planet.Ocean.deltaP``: Pressure step size (MPa)

Ice Shell
^^^^^^^^^

- ``Planet.Ocean.deltaP``: Pressure step for ice layers (MPa)
- ``Planet.Sil.icePhi_frac``: Ice porosity (fraction)
- ``Planet.Sil.icePclosure_MPa``: Pore closure pressure (MPa)

Interior
^^^^^^^^

- ``Planet.Do.Fe_CORE``: Include iron core (True/False)
- ``Planet.Core.xFeS``: Iron sulfide fraction in core
- ``Planet.Sil.mantleEOS``: Mantle equation of state file
- ``Planet.Sil.silPhi_frac``: Silicate porosity
- ``Planet.Sil.Qrad_Wkg``: Radiogenic heating (W/kg)
- ``Planet.Sil.Htidal_Wm3``: Tidal heating (W/m³)

Computational Settings
^^^^^^^^^^^^^^^^^^^^^^

- ``Params.CALC_NEW``: Force recalculation (True) or reload (False)
- ``Params.DO_PARALLEL``: Use parallel processing
- ``Params.SKIP_INNER``: Skip mantle/core calculations
- ``Params.SKIP_INDUCTION``: Skip magnetic induction

Example: Changing Ocean Salinity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let's modify Europa's ocean salinity::

    # In Europa/PPEuropa.py
    Planet.Ocean.wOcean_ppt = 50.0  # Increase salinity to 50 ppt

Run the model::

    python -m PlanetProfile.Main Europa

The output will show how the increased salinity affects:

- Ocean conductivity (higher)
- Freezing temperature (lower)
- Ice shell thickness (thinner)
- Magnetic induction amplitude (higher)

Next Steps
----------

Now that you know the basics, explore:

- :doc:`tutorials`: Step-by-step examples for common tasks
- :doc:`configuration`: Complete configuration reference
- :doc:`exploreogram`: Parameter space exploration
- :doc:`gui_guide`: Comprehensive GUI usage guide
- :doc:`api_reference`: Python API documentation

Need Help?
----------

- **Documentation**: https://vancesteven.github.io/PlanetProfile
- **GitHub Issues**: https://github.com/vancesteven/PlanetProfile/issues
- **Email**: steven.d.vance@jpl.nasa.gov

Citing PlanetProfile
--------------------

If you use PlanetProfile in your research, please cite:

- Vance et al. (2018) Geophysical investigations of habitability in ice-covered ocean worlds. *Journal of Geophysical Research: Planets*, 10.1002/2017JE005341
- Styczinski, Vance, and Melwani Daswani (2023) PlanetProfile: Self-consistent interior structure modeling for ocean worlds and rocky dwarf planets in Python. *Earth and Space Science*, 10(8), 10.1029/2022ea002748
