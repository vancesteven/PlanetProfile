import os, shutil
from glob import glob
from PlanetProfile import _ROOT, _Defaults, _SPICE, configTemplates, configLocals, CopyOnlyIfNeeded

def PPinstall():
    """ Copies default body files from PlanetProfile/Default/Body/ directories to Body/ directories found here. 
        The input files found in the Body/PPBody.py files here will override the defaults.
        Same for config files in this directory--they will override settings found in config files
        in the PlanetProfile/ directory.
    """

    print('Copying any Body/PPBody.py, Body/inductionData/*.txt, and config files that aren\'t ' +
             'in the right places from their defaults.')

    # Copy default profiles
    Defaults = glob(os.path.join(_Defaults, '*', 'PP*.py'))
    Copies = [file.split(f'Default{os.sep}')[-1] for file in Defaults]
    for Body, Destination in zip(Defaults, Copies):
        CopyOnlyIfNeeded(Body, Destination)

    # Copy inductionData contents
    inductionFiles = glob(os.path.join(_Defaults, '*', 'inductionData', '*.txt'))
    copies = [file.split(f'Default{os.sep}')[-1] for file in inductionFiles]
    for file, copy in zip(inductionFiles, copies):
        CopyOnlyIfNeeded(file, copy)

    # Copy config files
    for template, local in zip(configTemplates, configLocals):
        CopyOnlyIfNeeded(template, local)

    # Copy SPICE kernels and readme
    spiceFiles = glob(os.path.join(_SPICE, '*'))
    spiceCopies = [file.split(f'{_ROOT}{os.sep}')[-1] for file in spiceFiles]
    for repoSpice, localSpice in zip(spiceFiles, spiceCopies):
        CopyOnlyIfNeeded(repoSpice, localSpice)


def PPuninstall():
    """ Removes the same contents as what is added by PPinstall.
    """
    print('WARNING: This will remove all files within Body folders, the SPICE folder, the folders themselves, and configPP files.')
    answer = input('Continue? (y/n) ')

    if answer in ['y', 'Y', 'yes', 'Yes', 'YES']:
        print('\nRemoving all Body folders, SPICE kernels, and configPP files.')
        Bodies = [os.path.basename(bodyDir) for bodyDir in glob(os.path.join(_Defaults, '*'))]
        [shutil.rmtree(Body, ignore_errors=True) for Body in Bodies]
        shutil.rmtree(os.path.basename(_SPICE), ignore_errors=True)
        [os.remove(cfg) for cfg in configLocals]

    elif answer in ['n', 'N', 'no', 'No', 'NO']:
        print('Aborting.')
        exit(0)
    else:
        raise ValueError(f'Response "{answer}" not recognized.')


if __name__ == '__main__':
    PPinstall()
