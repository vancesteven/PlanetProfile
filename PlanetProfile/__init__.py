import os, shutil

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


# Set accessible file paths for installation directory and default configs
_ROOT = os.path.abspath(os.path.dirname(__file__))
_defaultConfig = os.path.join(_ROOT, 'defaultConfig.py')
_defaultConfigPlots = os.path.join(_ROOT, 'Plotting', 'defaultConfigPlots.py')
_defaultConfigInduct = os.path.join(_ROOT, 'MagneticInduction', 'defaultConfigInduct.py')
_Defaults = os.path.join(_ROOT, 'Default')
_Test = os.path.join(_ROOT, 'Test')
_TestImport = 'PlanetProfile.Test'
_SPICE = os.path.join(_ROOT, 'SPICE')

# Copy user config files to local dir if the user does not have them yet
_userConfig = 'configPP.py'
_userConfigPlots = 'configPPplots.py'
_userConfigInduct = 'configPPinduct.py'
configTemplates = [_defaultConfig, _defaultConfigPlots, _defaultConfigInduct]
configLocals = [_userConfig, _userConfigPlots, _userConfigInduct]
for template, local in zip(configTemplates, configLocals):
    CopyOnlyIfNeeded(template, local)

