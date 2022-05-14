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
from PlanetProfile import _Test, _TestImport
from PlanetProfile.GetConfig import Params
from PlanetProfile.Main import PlanetProfile, InductOgram, ReloadInductOgram
from PlanetProfile.Plotting.ProfilePlots import PlotInductOgram
from PlanetProfile.Test.TestBayes import TestBayes

def full():
    testBase = f'{_TestImport}.PPTest'

    # Set general testing config atop standard config options
    Params.CALC_NEW = True
    Params.CALC_NEW_REF = True
    Params.CALC_NEW_INDUC = True
    Params.CALC_SEISMIC = True
    Params.CALC_CONDUCT = True
    Params.RUN_ALL_PROFILES = False
    Params.COMPARE = False
    Params.NO_SAVEFILE = False
    Params.DO_INDUCTOGRAM = False

    # Get total number of test files to run
    fList = fnmatch.filter(os.listdir(_Test), 'PPTest*')
    fList = [fName for fName in fList if 'Induct' not in fName and 'Bayes' not in fName]
    nTests = np.size(fList)
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

    # Run all types of inductogram on Test7, with Test11 for porosity
    Params.DO_INDUCTOGRAM = True
    Params.NO_SAVEFILE = True
    for inductOtype in ['sigma', 'Tb', 'rho']:
        Params.Induct.inductOtype = inductOtype
        _ = TestInductOgram(7, Params)
        Params.DO_PARALLEL = False
        InductionLbl = TestInductOgram(7, Params)
        Params.DO_PARALLEL = True
        TestPlanets = np.append(TestPlanets, deepcopy(InductionLbl))
        tMarks = np.append(tMarks, time.time())

        InductionLbl = TestInductOgram(7, Params, CALC_NEW=False)
        TestPlanets = np.append(TestPlanets, deepcopy(InductionLbl))
        tMarks = np.append(tMarks, time.time())

    for inductOtype in ['phi']:
        Params.Induct.inductOtype = inductOtype
        _ = TestInductOgram(11, Params)
        Params.DO_PARALLEL = False
        InductionLbl = TestInductOgram(11, Params)
        Params.DO_PARALLEL = True
        TestPlanets = np.append(TestPlanets, deepcopy(InductionLbl))
        tMarks = np.append(tMarks, time.time())

        InductionLbl = TestInductOgram(11, Params, CALC_NEW=False)
        TestPlanets = np.append(TestPlanets, deepcopy(InductionLbl))
        tMarks = np.append(tMarks, time.time())
    Params.DO_INDUCTOGRAM = False
    Params.NO_SAVEFILE = False
    Params.SKIP_INNER = False

    # Test Bayesian analysis UpdateRun capabilities
    PlanetBayes, _ = TestBayes('Test')
    PlanetBayes.saveLabel = 'Bayes'
    TestPlanets = np.append(TestPlanets, PlanetBayes)
    tMarks = np.append(tMarks, time.time())

    log.info('Testing complete!')

    dt = np.diff(tMarks)
    log.debug('Elapsed time:\n' + '\n'.join([f'    {dt[i]:.3f} s for {Planet.name} {Planet.saveLabel}' for i, Planet in enumerate(TestPlanets)]))
    log.debug(f'Total elapsed time: {tMarks[-1] - tMarks[0]:.1f} s')
    return


def TestInductOgram(testNum, Params, CALC_NEW=True):
    testName = f'Test{testNum}'

    # Set sizes low so things don't take ages to run
    Params.Induct.nwPts, Params.Induct.nTbPts, Params.Induct.nphiPts,  \
    Params.Induct.nrhoPts, Params.Induct.nSigmaPts, Params.Induct.nDpts \
        = (4 for _ in range(6))
    
    if CALC_NEW:
        InductionLbl, Params = InductOgram(testName, Params)
        end = ''
    else:
        InductionLbl, Params = ReloadInductOgram(testName, Params)
        end = ' RELOAD'
    PlotInductOgram(InductionLbl, Params)
    InductionLbl.name = testName
    InductionLbl.saveLabel = f'{Params.Induct.inductOtype} induct-o-gram{end}'
    return InductionLbl


def simple():
    # Include timestamps in messages and force debug level logging
    log.basicConfig(level=log.DEBUG, format='[%(levelname)s] %(asctime)s - %(message)s')
    testMod = f'{_TestImport}.PPTest'

    # Set general testing config atop standard config options
    Params.CALC_NEW = True
    Params.CALC_NEW_REF = True
    Params.CALC_NEW_INDUCT = True
    Params.CALC_SEISMIC = True
    Params.CALC_CONDUCT = True
    Params.RUN_ALL_PROFILES = False
    Params.COMPARE = False
    Params.DO_INDUCTOGRAM = False
    Params.DO_PARALLEL = False

    iTest = 3
    bodyname = f'{testMod}{iTest}'
    tStart = time.time()
    testPlanet = importlib.import_module(bodyname).Planet
    log.info(f'Test case body: {bodyname}')
    if Params.DO_INDUCTOGRAM:
        TestInductOgram(iTest, Params, CALC_NEW=Params.CALC_NEW_INDUCT)
    else:
        _ = PlanetProfile(testPlanet, Params)
    tEnd = time.time()
    log.info('Simple test complete!')

    log.debug(f'Elapsed time:\n    {tEnd - tStart:.1f} s')
    return


if __name__ == '__main__':
    # Include timestamps in messages and force debug level logging
    root = log.getLogger()
    if root.handlers:
        for handler in root.handlers:
            root.removeHandler(handler)
    log.basicConfig(level=log.DEBUG, format='[%(levelname)s] %(asctime)s - %(message)s')
    log.getLogger().setLevel(log.DEBUG)

    if len(sys.argv) > 1:
        # Test type was passed as command line argument
        testType = sys.argv[1]
    else:
        testType = 'full'

    if testType == 'simple':
        simple()
    elif testType == 'Bayes':
        _, _ = TestBayes('Test')
    else:
        full()
