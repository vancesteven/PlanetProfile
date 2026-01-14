import os, shutil
from glob import glob
from pathlib import Path
from urllib.request import urlopen, urlretrieve
import json

import numpy as np

from PlanetProfile import _ROOT, _Defaults, _SPICE, configTemplates, configLocals, CopyOnlyIfNeeded, RemoveCarefully


def DownloadPerplexFiles():
    """Download Perple_X folder from GitHub if not present."""
    perplex_dir = Path(_ROOT) / 'Thermodynamics' / 'EOStables' / 'Perple_X'
    # Check if folder exists and has files
    if perplex_dir.exists() and any(f.suffix == '.tab' for f in perplex_dir.iterdir()):
        print('\n Perple_X folder already populated.')
        return
    
    # Create folder and download
    perplex_dir.mkdir(parents=True, exist_ok=True)
    print('\nDownloading Perple_X folder from GitHub...')
    print('='*60)
    
    # Get file list from GitHub API
    api_url = "https://api.github.com/repos/vancesteven/PlanetProfile/contents/PlanetProfile/Thermodynamics/EOStables/Perple_X"
    try:
        with urlopen(api_url) as response:
            files = [item['name'] for item in json.loads(response.read()) if item['type'] == 'file']
    except:
        # Fallback list if API fails
        files = ['Fe-S_3D_EOS.mat', 'CM_undifferentiated_hhph_DEW17_nofluid_nomelt_685.tab',
                 'CI_undifferentiated_hhph_DEW17_nofluid_nomelt_685.tab', 'CV_undifferentiated_v4_687_DEW17_nofluid_nomelt_v2.tab',
                 'Comet_67P_CG_v7_excluding_fluid_properties.tab', 'CV3hy1wt_678_1.tab',
                 'CM_hydrous_differentiated_Ganymede_Core080Fe020S_excluding_fluid_properties.tab']
    
    # Download each file
    base_url = "https://raw.githubusercontent.com/vancesteven/PlanetProfile/main/PlanetProfile/Thermodynamics/EOStables/Perple_X/"
    for filename in files:
        print(f'  {filename}...', end='', flush=True)
        try:
            urlretrieve(base_url + filename, perplex_dir / filename)
            print(' ✓')
        except:
            print(' ✗')
    
    print('='*60)
    print('Perplex Download complete\n')


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
    
    # Download Perple_X folder if not present
    DownloadPerplexFiles()


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
