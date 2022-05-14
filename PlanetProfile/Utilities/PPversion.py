from importlib.metadata import version
# Current PlanetProfile version tag
ppVerNum = '2.0.0-dev'
# Compatible version tag numbers
compatNums = {
    'seafreeze': '0.9.2',
    'gsw': '3.4.0',
    'obspy': '1.2.2'
}
# Printable package names
pkgNames = {
    'seafreeze': 'SeaFreeze',
    'gsw': 'GSW',
    'obspy': 'ObsPy.TauP'
}
# Instructions for installation
installInstruct = {
    'seafreeze': 'pip3 install SeaFreeze',
    'gsw': 'conda install -c conda-forge gsw',
    'obspy': 'conda install -c conda-forge obspy'
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
                    f'PlanetProfile is marked compatible with v{compatVer}.'
    # Check each version number tag hierarchically to see if we have tested with a newer version
    if((pkgCompatNums[0] > pkgVerNums[0]) or
      ((pkgCompatNums[0] == pkgVerNums[0]) and (pkgCompatNums[1] > pkgVerNums[1])) or
      ((pkgCompatNums[0] == pkgVerNums[0]) and (pkgCompatNums[1] == pkgVerNums[1]) and (pkgCompatNums[2] > pkgVerNums[2]))):
        print(pkgVerWarning)