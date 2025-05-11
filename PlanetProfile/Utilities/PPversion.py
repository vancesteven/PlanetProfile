"""
PlanetProfile version checking module
This module checks the version of PlanetProfile and its dependencies.
It is intended to be imported by other modules in the package.

To update the version number, change the value of the file number in PPverNum.txt and update the setup.py file 
in the root directory.
The version number in PPverNum.txt should be the same as the version number in setup.py.

To add packages to the compatibility check, add them to the compatNums dictionary, pkgNames dictionary, 
installInstruct dictionary, and upgradeInstruct dictionary.
Also add the package to the dependencies in setup.py.
"""
from importlib.metadata import version
from pathlib import Path
from PlanetProfile import _PPverNumFile
# Current PlanetProfile version tag
ppVerNum = Path(_PPverNumFile).read_text().replace('\n', '')
# Compatible version tag numbers
compatNums = {
    'seafreeze': '0.9.3',
    'gsw': '3.4.0',
    'obspy': '1.2.2',
    'MoonMag': '1.5.1',
    'reaktoro': '2.12.5',
    'pyalma3': '1.0.1'
}
# Printable package names
pkgNames = {
    'seafreeze': 'SeaFreeze',
    'gsw': 'GSW',
    'obspy': 'ObsPy.TauP',
    'MoonMag': 'MoonMag',
    'reaktoro': 'reaktoro',
    'pyalma3': 'pyALMA3'
}
# Instructions for installation
installInstruct = {
    'seafreeze': 'pip install SeaFreeze',
    'gsw': 'conda install -c conda-forge gsw',
    'obspy': 'conda install -c conda-forge obspy',
    'MoonMag': 'pip install MoonMag',
    'reaktoro': 'pip install reaktoro',
    'pyalma3': 'pip install pyalma3'
}
# Instructions for upgrading
upgradeInstruct = {
    'seafreeze': 'pip install --upgrade SeaFreeze',
    'gsw': 'conda update gsw',
    'obspy': 'conda update obspy',
    'MoonMag': 'pip install --upgrade MoonMag',
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
