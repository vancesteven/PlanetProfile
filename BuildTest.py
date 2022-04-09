""" This is a test script meant to check that the current PlanetProfile
    build maintains all functionality across updates. Run this script
    as "python Testing.py" from the main PlanetProfile directory
    before committing any major updates. New major functionality added
    to PlanetProfile should be accompanied by a new PPTest#.py test
    body in the Test/ directory.
"""

from PlanetProfile import PlanetProfile
from config import Params
import logging as log
import importlib, os, fnmatch, sys, copy

def full():
    # Include timestamps in messages and force debug level logging
    log.basicConfig(level=log.DEBUG, format='[%(levelname)s] %(asctime)s - %(message)s')
    testBase = 'Test.PPTest'

    # Set general testing config atop standard config options
    Params.CALC_NEW = True
    Params.CALC_NEW_REF = True
    Params.CALC_NEW_INDUC = True
    Params.CALC_SEISMIC = True
    Params.CALC_CONDUCT = True
    Params.RUN_ALL_PROFILES = False
    Params.COMPARE = False

    # Get first test profile separately, as we will reuse it
    testPlanet1 = importlib.import_module(f'{testBase}{1}').Planet
    _ = PlanetProfile(copy.deepcopy(testPlanet1), Params)

    # Loop over remaining test profiles (2 onwards)
    nTests = len(fnmatch.filter(os.listdir('Test'), 'PPTest*'))
    for i in range(15, nTests+1):
        testPlanetN = importlib.import_module(f'{testBase}{i}').Planet
        log.info(f'Test case body: {testBase}{i}')
        _ = PlanetProfile(testPlanetN, Params)

        # Verify that we can reload things as needed in each case
        Params.CALC_NEW = False
        Params.CALC_NEW_REF = False
        Params.CALC_NEW_INDUC = False
        _ = PlanetProfile(testPlanetN, Params)

        Params.CALC_NEW = True
        Params.CALC_NEW_REF = True
        Params.CALC_NEW_INDUC = True

    testPlanet1.name = 'Test0'
    # Test that we can successfully run things not including parallelization options
    Params.DO_PARALLEL = False
    _ = PlanetProfile(copy.deepcopy(testPlanet1), Params)
    Params.DO_PARALLEL = True

    # Make sure our auxiliary calculation flags work correctly
    Params.CALC_SEISMIC = False
    Params.CALC_CONDUCT = False
    _ = PlanetProfile(copy.deepcopy(testPlanet1), Params)
    Params.CALC_SEISMIC = True
    Params.CALC_CONDUCT = True

    # Check that skipping layers/portions works correctly
    Params.SKIP_INNER = True
    _ = PlanetProfile(copy.deepcopy(testPlanet1), Params)
    Params.SKIP_INNER = False

    log.info('Testing complete!')
    return


def simple():
    # Include timestamps in messages and force debug level logging
    log.basicConfig(level=log.DEBUG, format='[%(levelname)s] %(asctime)s - %(message)s')
    testMod = 'Test.PPTest'

    # Set general testing config atop standard config options
    Params.CALC_NEW = True
    Params.CALC_NEW_REF = True
    Params.CALC_NEW_INDUC = True
    Params.CALC_SEISMIC = True
    Params.CALC_CONDUCT = True
    Params.RUN_ALL_PROFILES = False
    Params.COMPARE = False

    iTest = 9
    testPlanet = importlib.import_module(f'{testMod}{iTest}').Planet
    log.info(f'Test case body: {testMod}{iTest}')
    _ = PlanetProfile(testPlanet, Params)
    log.info('Simple test complete!')
    return


if __name__ == '__main__':
    if len(sys.argv) > 1:
        # Test type was passed as command line argument
        testType = sys.argv[1]
    else:
        testType = 'full'

    if testType == 'simple':
        simple()
    else:
        full()
