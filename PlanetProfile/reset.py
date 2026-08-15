import os, shutil
from glob import glob
import logging as log
from PlanetProfile import _ROOT, _Defaults, _SPICE, configTemplates, configLocals, CopyCarefully

def PPreset():
    """ Copies default body files from PlanetProfile/Default/Body/ directories to Body/ directories found here.
        The input files found in the Body/PPBody.py files here will override the defaults.
        Same for config files in this directory--they will override settings found in config files
        in the PlanetProfile/ directory.
    """

    print('WARNING: This will reset all Body/PPBody.py, Body/inductionData/*.txt, and config files to their defaults.')
    answer = input('Continue? (y/n) ')

    if answer in ['y', 'Y', 'yes', 'Yes', 'YES']:
        # Copy default profiles
        Defaults = glob(os.path.join(_Defaults, '*', 'PP*.py'))
        Copies = [file.split(f'Default{os.sep}')[-1] for file in Defaults]
        for Body, Destination in zip(Defaults, Copies):
            CopyCarefully(Body, Destination)

        # Copy inductionData contents
        inductionFiles = glob(os.path.join(_Defaults, '*', 'inductionData', '*.txt'))
        copies = [file.split(f'Default{os.sep}')[-1] for file in inductionFiles]
        for file, copy in zip(inductionFiles, copies):
            CopyCarefully(file, copy)

        # Copy config files
        for template, local in zip(configTemplates, configLocals):
            CopyCarefully(template, local)

        # Copy SPICE kernels and readme
        spiceFiles = glob(os.path.join(_SPICE, '*'))
        spiceCopies = [file.split(f'{_ROOT}{os.sep}')[-1] for file in spiceFiles]
        for repoSpice, localSpice in zip(spiceFiles, spiceCopies):
            CopyCarefully(repoSpice, localSpice)

    elif answer in ['n', 'N', 'no', 'No', 'NO']:
        print('Aborting.')
        exit(0)
    else:
        raise ValueError(f'Response "{answer}" not recognized.')

if __name__ == '__main__':
    PPreset()
