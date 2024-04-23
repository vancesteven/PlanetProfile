import os, shutil
import matplotlib.pyplot as plt

def CopyCarefully(source, destination):
    try:
        if os.path.dirname(destination) != '':
            os.makedirs(os.path.dirname(destination), exist_ok=True)
        shutil.copy(source, destination)
    except OSError as err:
        raise OSError(f'Unable to copy from {source} to {destination}. ' +
                      f'Check that you have write permission for {os.getcwd()}. ' +
                      f'The error reported was:\n{err}')
    else:
        print(f'{destination} was copied from default at {source}.')

    return

def CopyOnlyIfNeeded(source, destination):
    if not os.path.isfile(destination):
        try:
            if os.path.dirname(destination) != '':
                os.makedirs(os.path.dirname(destination), exist_ok=True)
            shutil.copy(source, destination)
        except OSError as err:
            raise OSError(f'Unable to copy from {source} to {destination}. ' +
                          f'Check that you have write permission for {os.getcwd()}. ' +
                          f'The error reported was:\n{err}')
        else:
            print(f'{destination} was copied from default at {source}.')

    return

def RemoveCarefully(file):
    if os.path.isfile(file):
        os.remove(file)

# Set accessible file paths for installation directory and default configs
_ROOT = os.path.abspath(os.path.dirname(__file__))
_defaultConfig = os.path.join(_ROOT, 'defaultConfig.py')
_defaultConfigPlots = os.path.join(_ROOT, 'Plotting', 'defaultConfigPlots.py')
_defaultConfigInduct = os.path.join(_ROOT, 'MagneticInduction', 'defaultConfigInduct.py')
_defaultConfigTrajec = os.path.join(_ROOT, 'TrajecAnalysis', 'defaultConfigTrajec.py')
_Defaults = os.path.join(_ROOT, 'Default')
_DefaultList = next(os.walk(_Defaults))[1]
_Test = os.path.join(_ROOT, 'Test')
_TestImport = 'PlanetProfile.Test'
_PPverNumFile = os.path.join(_ROOT, 'Utilities', 'PPverNum.txt')
_healpixSphere = os.path.join(_ROOT, 'Utilities', 'healpixSphere.txt')
_SPICE = os.path.join(_ROOT, 'SPICE')

# Copy user config files to local dir if the user does not have them yet
_userConfig = 'configPP.py'
_userConfigPlots = 'configPPplots.py'
_userConfigInduct = 'configPPinduct.py'
_userConfigTrajec = 'configPPtrajec.py'
configTemplates = [_defaultConfig, _defaultConfigPlots, _defaultConfigInduct, _defaultConfigTrajec]
configLocals = [_userConfig, _userConfigPlots, _userConfigInduct, _userConfigTrajec]
if any([not os.path.isfile(cfg) for cfg in configLocals]):
    if input(f'configPP files not found in pwd: {os.getcwd()}. Copy from defaults to local dir? ' +
             f'[y]/n ') in ['', 'y', 'Y', 'yes', 'Yes']:
        for template, local in zip(configTemplates, configLocals):
            CopyOnlyIfNeeded(template, local)

_defaultCycler = plt.rcParams['axes.prop_cycle']  # Default matplotlib cycler for colors/styles
