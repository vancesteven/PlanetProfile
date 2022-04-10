""" This is a test script meant to check that the current PlanetProfile
    build maintains all functionality across updates. Run this script
    as "python Testing.py" from the main PlanetProfile directory
    before committing any major updates. New major functionality added
    to PlanetProfile should be accompanied by a new PPTest#.py test
    body in the Test/ directory.
"""

import logging as log
import numpy as np
import importlib, os, fnmatch, sys, time
from copy import deepcopy
from PlanetProfile import PlanetProfile
from config import Params

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

    # Get total number of test files to run
    nTests = len(fnmatch.filter(os.listdir('Test'), 'PPTest*'))
    # Create list for tracking outputs
    TestPlanets = np.empty(0, dtype=object)
    # Record start time and create t record array
    tMarks = np.zeros(0)
    tMarks = np.append(tMarks, time.time())
    # Get first test profile separately, as we will reuse it
    testPlanet1 = importlib.import_module(f'{testBase}{1}').Planet
    TestPlanets = np.append(TestPlanets, PlanetProfile(deepcopy(testPlanet1), Params)[0])
    tMarks = np.append(tMarks, time.time())

    # Loop over remaining test profiles (2 onwards)
    for i in range(2, nTests+1):
        testPlanetN = importlib.import_module(f'{testBase}{i}').Planet
        log.info(f'Test case body: {testBase}{i}')
        TestPlanets = np.append(TestPlanets, PlanetProfile(testPlanetN, Params)[0])
        tMarks = np.append(tMarks, time.time())

        # Verify that we can reload things as needed in each case
        Params.CALC_NEW = False
        Params.CALC_NEW_REF = False
        Params.CALC_NEW_INDUC = False
        TestPlanets = np.append(TestPlanets, PlanetProfile(deepcopy(testPlanetN), Params)[0])
        TestPlanets[-1].saveLabel += ' RELOAD'
        tMarks = np.append(tMarks, time.time())

        Params.CALC_NEW = True
        Params.CALC_NEW_REF = True
        Params.CALC_NEW_INDUC = True

    testPlanet1.name = 'Test0'
    # Test that we can successfully run things not including parallelization options
    Params.DO_PARALLEL = False
    TestPlanets = np.append(TestPlanets, PlanetProfile(deepcopy(testPlanet1), Params)[0])
    TestPlanets[-1].saveLabel += ' NO_PARALLEL'
    tMarks = np.append(tMarks, time.time())
    Params.DO_PARALLEL = True

    # Make sure our auxiliary calculation flags work correctly
    Params.CALC_SEISMIC = False
    Params.CALC_CONDUCT = False
    TestPlanets = np.append(TestPlanets, PlanetProfile(deepcopy(testPlanet1), Params)[0])
    TestPlanets[-1].saveLabel += ' NO_SEISMIC_OR_CONDUCT'
    tMarks = np.append(tMarks, time.time())
    Params.CALC_SEISMIC = True
    Params.CALC_CONDUCT = True

    # Check that skipping layers/portions works correctly
    Params.SKIP_INNER = True
    TestPlanets = np.append(TestPlanets, PlanetProfile(deepcopy(testPlanet1), Params)[0])
    TestPlanets[-1].saveLabel += ' SKIP_INNER'
    tMarks = np.append(tMarks, time.time())
    Params.SKIP_INNER = False

    log.info('Testing complete!')

    dt = np.diff(tMarks)
    log.debug('Elapsed time:\n' + '\n'.join([f'    {dt[i]:.3f} s for {Planet.name} {Planet.saveLabel}' for i, Planet in enumerate(TestPlanets)]))
    log.debug(f'Total elapsed time: {tMarks[-1] - tMarks[0]:.1f} s')
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

    iTest = 1
    tStart = time.time()
    testPlanet = importlib.import_module(f'{testMod}{iTest}').Planet
    log.info(f'Test case body: {testMod}{iTest}')
    _ = PlanetProfile(testPlanet, Params)
    tEnd = time.time()
    log.info('Simple test complete!')

    log.debug(f'Elapsed time:\n    {tEnd - tStart:.1f} s')
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
