from importlib.metadata import version
# Current PlanetProfile version tag
verNum = '1.2.0-dev'
# Compatible SeaFreeze version tag
seaCompatVer = '0.9.2'
# Compatible GSW version tag
gswCompatVer = ''
# Compatible TauP version tag
taupCompatVer = ''


def CheckSeaFreeze(compatVer):
    # Grab version number from SeaFreeze package
    seaVer = version('seafreeze')
    seaCompatNums = [int(numStr) for numStr in compatVer.split('.')]
    seaVerNums = [int(numStr) for numStr in seaVer.split('.')[:3]]
    seaVerWarning = 'WARNING: Installed SeaFreeze version is ' + seaVer + ' but this version of ' + \
        'PlanetProfile is compatible with SeaFreeze v' + compatVer + '.'
    # Check each version number tag hierarchically to see if we have tested with a newer version
    if((seaCompatNums[0] > seaVerNums[0]) or
      ((seaCompatNums[0] == seaVerNums[0]) and (seaCompatNums[1] > seaVerNums[1])) or
      ((seaCompatNums[0] == seaVerNums[0]) and (seaCompatNums[1] == seaVerNums[1]) and (seaCompatNums[2] > seaVerNums[2]))):
        print(seaVerWarning)

def CheckGSW(compatVer):
    print('CheckGSW not implemented yet')


def CheckTauP(compatVer):
    print('CheckTauP not implemented yet')