# PlanetProfile v2.2
Software framework for constructing 1D interior structure models based on planetary properties. Self-consistent thermodynamics are used for fluid, rock, and mineral phases. Sound speeds, attenuation, and electrical conductivities are computed as outputs. The main code is called from an input file containing all the planetary data.

The main repository is mirrored at https://github.com/NASA-Planetary-Science/PlanetProfile; any pull requests should be submitted to https://github.com/vancesteven/PlanetProfile.

## Acknowledging PlanetProfile
We want to hear about your work with PlanetProfile! Please consider sending us a message alerting us to your work (svance@jpl.caltech.edu). Suggested acknowledgement in publications: "Data used in this work were generated using the open source PlanetProfile software hosted on GitHub."

## Getting started
PlanetProfile is available in Python and Matlab.

*For Python,* the recommended way to install is with with pip:
1. (Recommended) Install all dependencies listed in the next section before proceeding.
1. At a terminal:
`python -m pip install PlanetProfile`
Python 3.8 or higher is required. Pip will install dependencies, but a conda environment with the prerequisites listed below is recommended.
1. Create a directory where you'd like to store configurations and have folders for each body.
1. Navigate into the new directory.
1. At a terminal:
`python -m PlanetProfile.install PPinstall`
This will copy files from their defaults to the current directory.
1. Run the software with, for example:
`python -m PlanetProfile.Main Europa`
or
`python -m PlanetProfile.Main path/to/PPBody.py`
or in a Python script with
`from PlanetProfile.Main import RunPPfile
RunPPfile('Europa', 'PPEuropa.py')`

*For Matlab:*
  1. Download or clone this repository.
  1. Install prerequisites below.
  1. At a terminal: 
  `make install`
  Or, add everything in the top-level directory except the PlanetProfile sub-folder to the Matlab path.
  1. In Matlab, set the current directory to the top-level directory of the downloaded repository (top PlanetProfile folder).
  1. Run the software with `PPEuropa` in the Matlab command prompt, or by opening and running one of the files located at Body/PPBody.m (e.g. Titan/PPTitan.m).
  
## Prerequisites
A simple list with install commands for Python is in the next section.
* SeaFreeze -- see https://github.com/Bjournaux/SeaFreeze
  * Python: Installed with pip: `pip install SeaFreeze`
  * Matlab: Download the repository to Thermodynamics/SeaFreeze and add the contents to the Matlab path.
* MoonMag -- see https://github.com/itsmoosh/MoonMag
  * Python only -- Installed with pip: `pip install MoonMag`
  * In Matlab, magnetic induction calculations are performed numerically and cannot incorporate lateral asymmetry in conducting layers.
* Gibbs Seawater toolbox of TEOS-10 -- see https://www.teos-10.org/
  * Python: Installed with conda via conda-forge: `conda install -c conda-forge gsw` 
  * Matlab: Already packaged into the PlanetProfile repository along with the original license.
* Perple_X -- see http://www.perplex.ethz.ch/
  * For both Python and Matlab, Perple_X outputs are currently hosted as part of the installation, in Thermodynamics/Perple_X for Matlab and in PlanetProfile/Thermodynamics/EOSdata/Perple_X for Python. The files we use were generated with Perple_X v6.7.9.
* TauP/ObsPy (optional) -- see https://www.seis.sc.edu/taup/
  * Python: Installed with conda via conda-forge: `conda install -c conda-forge obspy`
  * Matlab: Download mMatTauP contents into Utilities/ and add-with-subfolders to the Matlab path.
* A working TeX/LaTeX distribution (such as TeXlive) is recommended for optimum plot labels. TeXlive is available at: https://tug.org/texlive/acquire-netinstall.html

## Installation of prerequisites
### Python 
1. Python version 3.8+ must be installed, preferably via Anaconda. Required modules can be installed in Miniconda with the following command:
  1. `conda install numpy scipy matplotlib mpmath`
1. Conda-forge modules can be installed in Anaconda or Miniconda with the following command:
  1. `conda install -c conda-forge gsw obspy spiceypy cmasher`
1. AFTER the above modules have been installed with conda, install SeaFreeze and MoonMag with the following command:
  1. `pip install SeaFreeze MoonMag`
  1. Finally, install PlanetProfile with `pip install PlanetProfile`
  
### Matlab
1. Download PlanetProfile repository.
1. Download SeaFreeze repository to PlanetProfile/Thermodynamics/SeaFreeze/ (NOT PlanetProfile/PlanetProfile/Thermodynamics).
1. Add SeaFreeze folder and sub-folders to Matlab path.
Some magnetic field features require use of the [SPICE toolkit through Mice](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/MATLAB/index.html). To install Mice:
1. Navigate to https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html
1. Follow the link for your operating system and download the .zip or .tar.Z file to PlanetProfile/Utilities/spice/
1. Unpack the archive (into PlanetProfile/spice/mice/)
1. Add PlanetProfile/Utilities/spice/mice/src/mice/ and PlanetProfile/Utilities/spice/mice/lib/ to your Matlab path.
1. Install necessary SPICE kernels by downloading them from https://naif.jpl.nasa.gov/pub/naif/generic_kernels/ and placing them in PlanetProfile/Utilities/spice/. The planetary constants kernel (PCK) and leap-seconds kernel (TLS) are saved in this repository, but the generic ephemeris kernels (SPK, .bsp files) are too large for us to save here. There is one for each planet's satellites, located at https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/. Currently in use are:
  1. jup365.bsp
  1. sat427.bsp
  1. ura111.bsp
  1. nep095.bsp

## Contributing
PlanetProfile is open source software. Please see the [LICENSE](https://github.com/vancesteven/PlanetProfile/blob/master/LICENSE) file and read the guidelines for contrbuting in [CONTRIBUTING.md](https://github.com/vancesteven/PlanetProfile/blob/master/CONTRIBUTING.md) if you are interested in joining the project. Also see our community guidelines in [CODE_OF_CONDUCT.md](https://github.com/vancesteven/PlanetProfile/blob/master/CODE_OF_CONDUCT.md).

## Notes
* With the PlanetProfile 2.0 release, both Python and Matlab are available. The two branches do not have the same functionality yet with this release--some features exist in the Python version that are not yet implemented in the Matlab. A later release will align their functionality as much as possible. For now, the Python version is recommended.
* As of 2020-09-28, PlanetProfile v1.1.0 was released along with code for making calculations regarding magnetic induction. The development (master) branch of PlanetProfile is set up to generate profiles from minimal inputs. Output profiles that may be used along with the induction calculations may be found in the v1.1.0 release.
* The default settings include a recalculation of all parameters. It is recommended to recalculate all parameters whenever PlanetProfile is updated and any time a change in input parameters may affect layer thicknesses or other intermediate variables.

Some calculations in Python use parallel computing with the multiprocessing builtin module. There are sometimes cross-platform compatibility issues that crop up. By default, multiprocessing is enabled; disable it by setting Params.DO_PARALLEL = False in configPP.py.
