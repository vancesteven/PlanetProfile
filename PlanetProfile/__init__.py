import os, shutil

def CopyCarefully(source, destination):
    try:
        if os.path.dirname(destination) != '':
            os.makedirs(os.path.dirname(destination), exist_ok=True)
        shutil.copy(source, destination)
    except OSError:
        raise OSError(f'Unable to copy from {source} to {destination}. ' +
                      f'Check that you have write permission for {os.getcwd()}.')
    else:
        print(f'{destination} was copied from default at {source}.')

    return

def CopyOnlyIfNeeded(source, destination):
    if not os.path.isfile(destination):
        try:
            if os.path.dirname(destination) != '':
                os.makedirs(os.path.dirname(destination), exist_ok=True)
            shutil.copy(source, destination)
        except OSError:
            raise OSError(f'Unable to copy from {source} to {destination}. ' +
                          f'Check that you have write permission for {os.getcwd()}.')
        else:
            print(f'{destination} was copied from default at {source}.')

    return


# Set accessible file paths for installation directory and default configs
_ROOT = os.path.abspath(os.path.dirname(__file__))
_defaultConfig = os.path.join(_ROOT, 'defaultConfig.py')
_userConfig = os.path.join(_ROOT, 'userConfig.py')
_configPlots = os.path.join(_ROOT, 'Plotting', 'defaultConfigPlots.py')
_userConfigPlots = os.path.join(_ROOT, 'Plotting', 'userConfigPlots.py')
_configInduct = os.path.join(_ROOT, 'MagneticInduction', 'defaultConfigInduct.py')
_userConfigInduct = os.path.join(_ROOT, 'MagneticInduction', 'userConfigInduct.py')
_Defaults = os.path.join(_ROOT, 'Default')
_Test = os.path.join(_ROOT, 'Test')
_TestImport = 'PlanetProfile.Test'

# Copy user config files to local dir if the user does not have them yet
_localConfig = 'config.py'
_localConfigPlots = 'configPlots.py'
_localConfigInduct = 'configInduct.py'
configTemplates = [_userConfig, _userConfigPlots, _userConfigInduct]
configLocals = [_localConfig, _localConfigPlots, _localConfigInduct]
for template, local in zip(configTemplates, configLocals):
    if not os.path.isfile(local):
        CopyCarefully(template, local)

