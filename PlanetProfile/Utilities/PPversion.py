"""
PlanetProfile version checking module
This module checks the version of PlanetProfile and its dependencies.
It is intended to be imported by other modules in the package.

To update the version number, edit the version field in pyproject.toml (single source of truth).
Python automatically reads the version from installed package metadata.

To add packages to the compatibility check, add them to the compatNums dictionary, pkgNames dictionary, 
installInstruct dictionary, and upgradeInstruct dictionary.
Also add the package to the dependencies in pyproject.toml.
"""
from importlib.metadata import version, PackageNotFoundError

# Get PlanetProfile version from installed package metadata (pyproject.toml is source of truth)
try:
    ppVerNum = version('PlanetProfile')
except PackageNotFoundError:
    # Running in development mode (cloned repo, not installed via pip)
    # Read version from pyproject.toml to stay synchronized with source of truth
    import os
    from pathlib import Path

    current_file = Path(__file__)
    repo_root = current_file.parent.parent.parent
    pyproject_path = repo_root / 'pyproject.toml'

    ppVerNum = None
    if pyproject_path.exists():
        with open(pyproject_path, 'r') as f:
            for line in f:
                if line.strip().startswith('version'):
                    ppVerNum = line.split('=')[1].strip().strip('"').strip("'")
                    break

    if ppVerNum is None:
        ppVerNum = 'dev-unknown'

# Compatible version tag numbers
compatNums = {
    'seafreeze': '1.0.0',
    'gsw': '3.6.20',
    'obspy': '1.4.2',

    'reaktoro': '2.13.0',
    'pyalma3': '1.0.1'
}
# Printable package names
pkgNames = {
    'seafreeze': 'SeaFreeze',
    'gsw': 'GSW',
    'obspy': 'ObsPy.TauP',

    'reaktoro': 'reaktoro',
    'pyalma3': 'pyALMA3'
}
# Instructions for installation
installInstruct = {
    'seafreeze': 'pip install SeaFreeze',
    'gsw': 'conda install -c conda-forge gsw',
    'obspy': 'conda install -c conda-forge obspy',

    'reaktoro': 'pip install reaktoro',
    'pyalma3': 'pip install pyalma3'
}
# Instructions for upgrading
upgradeInstruct = {
    'seafreeze': 'pip install --upgrade SeaFreeze',
    'gsw': 'conda update gsw',
    'obspy': 'conda update obspy',

    'reaktoro': 'pip install --upgrade reaktoro',
    'pyalma3': 'pip install --upgrade pyalma3'
}


def CheckCompat(package):
    # Grab version number from package
    try:
        pkgVer = version(package)
    except ModuleNotFoundError:
        raise ModuleNotFoundError(f'{pkgNames[package]} is not installed. Install it with the command: ' +
                                  installInstruct[package])
    compatVer = compatNums[package]
    pkgCompatNums = [int(numStr) for numStr in compatVer.split('.')]
    pkgVerNums = [int(numStr) for numStr in pkgVer.split('.')[:3]]
    pkgVerWarning = f'WARNING: Installed {pkgNames[package]} version is {pkgVer} but this version of ' + \
                    f'PlanetProfile is marked compatible with v{compatVer}. Upgrade it with the command: ' + \
                    upgradeInstruct[package]
    # Check each version number tag hierarchically to see if we have tested with a newer version
    if((pkgCompatNums[0] > pkgVerNums[0]) or
      ((pkgCompatNums[0] == pkgVerNums[0]) and (pkgCompatNums[1] > pkgVerNums[1])) or
      ((pkgCompatNums[0] == pkgVerNums[0]) and (pkgCompatNums[1] == pkgVerNums[1]) and (pkgCompatNums[2] > pkgVerNums[2]))):
        print(pkgVerWarning)
