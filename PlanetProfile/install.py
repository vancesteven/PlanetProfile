import os, shutil
from glob import glob

import numpy as np

from PlanetProfile import _ROOT, _Defaults, _SPICE, configTemplates, configLocals, CopyOnlyIfNeeded, RemoveCarefully


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


def PPuninstall(KEEP_NEW=None):
    """
    Removes the same contents as what is added by PPinstall.
    If any directory would then be empty, the directory is deleted.

    Args:
        KEEP_NEW: If True, retain any files the user has created atop the default install.
            If False, remove everything in the Body/ and SPICE/ directories. Defaults to None (prompt user).

    Returns:
        None
    """
    print(
        'WARNING: This will remove all files within Body folders, the SPICE folder, the folders themselves, and configPP files ' +
        'that have been copied from the PlanetProfile defaults.')
    answer = input('Continue? (y/N) ')

    if answer in ['y', 'Y', 'yes', 'Yes', 'YES']:
        print('\nRemoving all copied Body folder contents, SPICE kernels, and configPP files.')

        # Remove config files
        [RemoveCarefully(cfg) for cfg in configLocals]

        Bodies = [os.path.basename(bodyDir) for bodyDir in glob(os.path.join(_Defaults, '*'))]
        # Remove PPBody.py files
        Defaults = glob(os.path.join(_Defaults, '*', 'PP*.py'))
        PPfiles = [file.split(f'Default{os.sep}')[-1] for file in Defaults]
        [RemoveCarefully(file) for file in PPfiles]

        # Remove induction files
        inductionDefaults = glob(os.path.join(_Defaults, '*', 'inductionData', '*.txt'))
        inductionLocal = [file.split(f'Default{os.sep}')[-1] for file in inductionDefaults]
        inductionDirs = [os.path.join(Body, 'inductionData') for Body in Bodies]
        [RemoveCarefully(file) for file in inductionLocal]

        # Remove SPICE files
        spiceDefaults = glob(os.path.join(_SPICE, '*'))
        spiceLocal = [file.split(f'{_ROOT}{os.sep}')[-1] for file in spiceDefaults]
        spiceDir = os.path.basename(_SPICE)
        [RemoveCarefully(file) for file in spiceLocal]

        # Check if each Body folder (and the SPICE folder) is empty, and if so, delete it.
        # Otherwise, warn user and ask to delete
        if os.path.isdir(spiceDir):
            if not os.listdir(spiceDir):
                os.rmdir(spiceDir)
            else:
                if KEEP_NEW is None:
                    print(f'{spiceDir} contains non-default files.')
                    answer = input('Delete them anyway? (y/N) ')
                    if answer in ['y', 'Y', 'yes', 'Yes', 'YES']:
                        shutil.rmtree(os.path.basename(_SPICE), ignore_errors=True)
                elif not KEEP_NEW:
                    shutil.rmtree(os.path.basename(_SPICE), ignore_errors=True)

        for Body, inductDir in zip(Bodies, inductionDirs):
            if os.path.isdir(inductDir):
                if not os.listdir(inductDir):
                    os.rmdir(inductDir)

            if os.path.isdir(Body):
                if not os.listdir(Body):
                    os.rmdir(Body)
                else:
                    if KEEP_NEW is None:
                        print(f'{Body} dir contains non-default files.')
                        answer = input('Delete them anyway? (y/N) ')
                        if answer in ['y', 'Y', 'yes', 'Yes', 'YES']:
                            shutil.rmtree(Body, ignore_errors=True)
                    elif not KEEP_NEW:
                        shutil.rmtree(Body, ignore_errors=True)
    else:
        print('Aborting.')
        exit(0)


if __name__ == '__main__':
    PPinstall()
