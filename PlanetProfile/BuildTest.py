""" This is a test script meant to check that the current PlanetProfile
    build maintains all functionality across updates. Run this script
    as "python Testing.py" from the main PlanetProfile directory
    before committing any major updates. New major functionality added
    to PlanetProfile should be accompanied by a new PPTest#.py test
    body in the Test/ directory.
"""

import logging
import numpy as np
import importlib, os, fnmatch, sys, time
from copy import deepcopy
from PlanetProfile import _Test, _TestImport
from PlanetProfile.GetConfig import Params as configParams
from PlanetProfile.Main import PlanetProfile, InductOgram, ReloadInductOgram, ExploreOgram, ReloadExploreOgram
from PlanetProfile.Plotting.ProfilePlots import PlotExploreOgram, PlotExploreOgramDsigma
from PlanetProfile.Plotting.MagPlots import PlotInductOgram
from PlanetProfile.Test.TestBayes import TestBayes

# Include timestamps in messages and force debug level logging for all testing
log = logging.getLogger('PlanetProfile')
stream = logging.StreamHandler()
stream.setFormatter(logging.Formatter('[%(levelname)s] %(asctime)s - %(message)s'))
log.setLevel(logging.DEBUG)
for deftHandler in log.handlers:
    log.removeHandler(deftHandler)
log.addHandler(stream)
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('PIL').setLevel(logging.WARNING)
logging.getLogger('MoonMag').setLevel(logging.DEBUG)


def full(iTestStart=2, skipType=None):
    testBase = f'{_TestImport}.PPTest'

    # Set general testing config atop standard config options
    Params = configParams
    Params.CALC_NEW = True
    Params.CALC_NEW_REF = True
    Params.CALC_NEW_INDUCT = True
    Params.CALC_SEISMIC = True
    Params.CALC_CONDUCT = True
    Params.DO_PARALLEL = False
    Params.RUN_ALL_PROFILES = False
    Params.COMPARE = False
    Params.NO_SAVEFILE = False
    Params.DO_INDUCTOGRAM = False
    Params.SKIP_INDUCTION = False

    # Get total number of test files to run
    fList = fnmatch.filter(os.listdir(_Test), 'PPTest*')
    fList = [fName for fName in fList if 'Induct' not in fName and 'Bayes' not in fName and 'Explore' not in fName]
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
    if iTestStart is None:
        iTestStart = 2

    if skipType is None:
        for i in range(iTestStart, nTests+1):
            testPlanetN = importlib.import_module(f'{testBase}{i}').Planet
            log.info(f'Test case body: {testBase}{i}')
            TestPlanets = np.append(TestPlanets, PlanetProfile(testPlanetN, Params)[0])
            tMarks = np.append(tMarks, time.time())

            # Verify that we can reload things as needed in each case
            Params.CALC_NEW = False
            Params.CALC_NEW_REF = False
            Params.CALC_NEW_INDUCT = False
            TestPlanets = np.append(TestPlanets, PlanetProfile(deepcopy(testPlanetN), Params)[0])
            TestPlanets[-1].saveLabel += ' RELOAD'
            tMarks = np.append(tMarks, time.time())

            Params.CALC_NEW = True
            Params.CALC_NEW_REF = True
            Params.CALC_NEW_INDUCT = True

        testPlanet1.name = 'Test0'
        # Test that we can successfully run standard profiles with parallelization options
        Params.DO_PARALLEL = True
        TestPlanets = np.append(TestPlanets, PlanetProfile(deepcopy(testPlanet1), Params)[0])
        TestPlanets[-1].saveLabel += ' NO_PARALLEL'
        tMarks = np.append(tMarks, time.time())

        # Make sure our auxiliary calculation flags work correctly
        Params.CALC_SEISMIC = False
        Params.CALC_CONDUCT = False
        TestPlanets = np.append(TestPlanets, PlanetProfile(deepcopy(testPlanet1), Params)[0])
        TestPlanets[-1].saveLabel += ' NO_SEISMIC_OR_CONDUCT'
        tMarks = np.append(tMarks, time.time())
        Params.CALC_SEISMIC = True
        Params.CALC_CONDUCT = True

    if skipType is None or skipType.lower() == 'bayes':
        # Test Bayesian analysis UpdateRun capabilities
        PlanetBayes, _ = TestBayes('Test')
        PlanetBayes.saveLabel = 'Bayes'
        TestPlanets = np.append(TestPlanets, PlanetBayes)
        tMarks = np.append(tMarks, time.time())

    # Check that skipping layers/portions works correctly
    Params.SKIP_INNER = True
    TestPlanets = np.append(TestPlanets, PlanetProfile(deepcopy(testPlanet1), Params)[0])
    TestPlanets[-1].saveLabel += ' SKIP_INNER'
    tMarks = np.append(tMarks, time.time())
    Params.SKIP_INNER = False

    if skipType is None or skipType.lower() == 'induct':
        # Test out all the inductogram config options
        TestPlanets, Params, tMarks = TestAllInductOgrams(TestPlanets, Params, tMarks)

    if skipType is None or skipType.lower() == 'explore' or skipType.lower() == 'explorewaterless':
        # Test out all the exploreogram config options
        if skipType is not None and skipType.lower() == 'explorewaterless':
            SKIP_HYDRO = True
        else:
            SKIP_HYDRO = False
        TestPlanets, Params, tMarks = TestAllExploreOgrams(TestPlanets, Params, tMarks, SKIP_HYDRO=SKIP_HYDRO)

    log.info('Testing complete!')

    dt = np.diff(tMarks)
    log.debug('Elapsed time:\n' + '\n'.join([f'    {dt[i]:.3f} s for {Planet.name} {Planet.saveLabel}' for i, Planet in enumerate(TestPlanets)]))
    log.debug(f'Total elapsed time: {tMarks[-1] - tMarks[0]:.1f} s')
    return


def TestAllInductOgrams(TestPlanets, Params, tMarks):
    # Run all types of inductogram on Test7, with Test11 for porosity
    Params.DO_INDUCTOGRAM = True
    Params.NO_SAVEFILE = True
    Params.SKIP_INNER = True
    for inductOtype in ['sigma', 'Tb', 'rho']:
        Params.Induct.inductOtype = inductOtype
        Params.DO_PARALLEL = False
        _ = TestInductOgram(7, Params)
        Params.DO_PARALLEL = True
        Induction = TestInductOgram(7, Params)
        TestPlanets = np.append(TestPlanets, deepcopy(Induction))
        tMarks = np.append(tMarks, time.time())

        Induction = TestInductOgram(7, Params, CALC_NEW=False)
        TestPlanets = np.append(TestPlanets, deepcopy(Induction))
        tMarks = np.append(tMarks, time.time())

    for inductOtype in ['phi']:
        Params.Induct.inductOtype = inductOtype
        Params.DO_PARALLEL = False
        _ = TestInductOgram(11, Params)
        Params.DO_PARALLEL = True
        Induction = TestInductOgram(11, Params)
        TestPlanets = np.append(TestPlanets, deepcopy(Induction))
        tMarks = np.append(tMarks, time.time())

        Induction = TestInductOgram(11, Params, CALC_NEW=False)
        TestPlanets = np.append(TestPlanets, deepcopy(Induction))
        tMarks = np.append(tMarks, time.time())
    Params.DO_INDUCTOGRAM = False
    Params.NO_SAVEFILE = False
    Params.SKIP_INNER = False

    return TestPlanets, Params, tMarks


def TestAllExploreOgrams(TestPlanets, Params, tMarks, SKIP_HYDRO=False):
    # Run all types of exploreogram on Test7, with Test5 for waterless
    Params.DO_EXPLOREOGRAM = True
    Params.NO_SAVEFILE = True
    Params.SKIP_INNER = True
    hydroExploreBds = {
        'xFeS': [0, 1],
        'rhoSilInput_kgm3': [2000, 4500],
        'wOcean_ppt': [0, 100],
        'Tb_K': [255, 273],
        'ionosTop_km': [5, 50],
        'sigmaIonos_Sm': [1e-6, 1e-2],
        'silPhi_frac': [0, 0.8],
        'silPclosure_MPa': [200, 750],
        'icePhi_frac': [0, 0.8],
        'icePclosure_MPa': [5, 40],
        'Htidal_Wm3': [1e-18, 1e-11],
        'Qrad_Wkg': [1e-20, 1e-14]
    }
    waterlessExploreBds = {
        'xFeS': [0, 1],
        'rhoSilInput_kgm3': [2000, 4500],
        'ionosTop_km': [5, 50],
        'sigmaIonos_Sm': [1e-6, 1e-2],
        'silPhi_frac': [0, 0.8],
        'silPclosure_MPa': [200, 750],
        'Htidal_Wm3': [1e-18, 1e-11],
        'Qrad_Wkg': [1e-20, 1e-14],
        'qSurf_Wm2': [50e-3, 400e-3]
    }

    if not SKIP_HYDRO:
        log.info('Running exploreOgrams for icy bodies for all input types.')
        for xName in hydroExploreBds.keys():
            for yName in hydroExploreBds.keys():
                if xName != yName:
                    Params.Explore.xName = xName
                    Params.Explore.yName = yName
                    Params.Explore.xRange = hydroExploreBds[xName]
                    Params.Explore.yRange = hydroExploreBds[yName]
                    Params.DO_PARALLEL = False
                    _ = TestExploreOgram(7, Params)
                    Params.DO_PARALLEL = True
                    Exploration = TestExploreOgram(7, Params)
                    TestPlanets = np.append(TestPlanets, deepcopy(Exploration))
                    tMarks = np.append(tMarks, time.time())

                    Exploration = TestExploreOgram(7, Params, CALC_NEW=False)
                    TestPlanets = np.append(TestPlanets, deepcopy(Exploration))
                    tMarks = np.append(tMarks, time.time())

    log.info('Running exploreOgrams for waterless bodies for all input types.')
    for xName in waterlessExploreBds.keys():
        for yName in waterlessExploreBds.keys():
            if xName != yName:
                Params.Explore.xName = xName
                Params.Explore.yName = yName
                Params.Explore.xRange = waterlessExploreBds[xName]
                Params.Explore.yRange = waterlessExploreBds[yName]
                Params.DO_PARALLEL = False
                _ = TestExploreOgram(5, Params)
                Params.DO_PARALLEL = True
                Exploration = TestExploreOgram(5, Params)
                TestPlanets = np.append(TestPlanets, deepcopy(Exploration))
                tMarks = np.append(tMarks, time.time())

                Exploration = TestExploreOgram(5, Params, CALC_NEW=False)
                TestPlanets = np.append(TestPlanets, deepcopy(Exploration))
                tMarks = np.append(tMarks, time.time())

    Params.DO_EXPLOREOGRAM = False
    Params.NO_SAVEFILE = False
    Params.SKIP_INNER = False

    return TestPlanets, Params, tMarks


def TestInductOgram(testNum, Params, CALC_NEW=True):
    testName = f'Test{testNum}'

    # Set sizes low so things don't take ages to run
    Params.Induct.nwPts, Params.Induct.nTbPts, Params.Induct.nphiPts,  \
    Params.Induct.nrhoPts, Params.Induct.nSigmaPts, Params.Induct.nDpts \
        = (4 for _ in range(6))
    
    if CALC_NEW:
        Induction, Params = InductOgram(testName, Params)
        end = ''
    else:
        Induction, _, Params = ReloadInductOgram(testName, Params)
        end = ' RELOAD'
    PlotInductOgram(Induction, Params)
    Induction.name = testName
    Induction.saveLabel = f'{Params.Induct.inductOtype} induct-o-gram{end}'

    return Induction


def TestExploreOgram(testNum, Params, CALC_NEW=True):
    testName = f'Test{testNum}'

    # Set sizes low so things don't take ages to run
    Params.Explore.nx, Params.Explore.ny \
        = (4 for _ in range(2))

    if CALC_NEW:
        Exploration, Params = ExploreOgram(testName, Params)
        end = ''
    else:
        Exploration, Params = ReloadExploreOgram(testName, Params)
        end = ' RELOAD'

    if isinstance(Params.Explore.zName, list):
        figNames = Params.FigureFiles.explore + []
        for zName, figName in zip(Params.Explore.zName, figNames):
            Exploration.zName = zName
            Params.FigureFiles.explore = figName
            PlotExploreOgram([Exploration], Params)
        Params.FigureFiles.explore = figNames
    PlotExploreOgramDsigma([Exploration], Params)
    Exploration.name = testName
    Exploration.saveLabel = f'{Params.Explore.xName} x {Params.Explore.yName} explore-o-gram{end}'

    return Exploration


def simple(iTests=None):
    testMod = f'{_TestImport}.PPTest'

    # Set general testing config atop standard config options
    Params = configParams
    Params.CALC_NEW = True
    Params.CALC_NEW_REF = True
    Params.CALC_NEW_INDUCT = True
    Params.CALC_SEISMIC = True
    Params.CALC_CONDUCT = True
    Params.RUN_ALL_PROFILES = False
    Params.COMPARE = False
    Params.DO_INDUCTOGRAM = False
    Params.DO_EXPLOREOGRAM = False
    Params.DO_PARALLEL = False
    Params.SKIP_INDUCTION = False

    tStart = time.time()
    for iTest in iTests:
        bodyname = f'{testMod}{iTest}'
        testPlanet = importlib.import_module(bodyname).Planet
        log.info(f'Test case body: {bodyname}')
        if Params.DO_INDUCTOGRAM:
            TestInductOgram(iTest, Params, CALC_NEW=Params.CALC_NEW_INDUCT)
        elif Params.DO_EXPLOREOGRAM:
            TestExploreOgram(iTest, Params, CALC_NEW=Params.CALC_NEW_INDUCT)
        else:
            _ = PlanetProfile(testPlanet, Params)
    tEnd = time.time()
    log.info('Simple test complete!')

    log.debug(f'Elapsed time:\n    {tEnd - tStart:.1f} s')
    return


if __name__ == '__main__':
    skipType = None
    if len(sys.argv) > 1:
        # Test type was passed as command line argument
        testType = sys.argv[1]
        if len(sys.argv) > 2:
            if sys.argv[2].isdigit():
                # Test profile number was passed as command line argument
                iTest = int(sys.argv[2])
            else:
                # Skip to specific testing section was passed as command line arg
                iTest = None
                skipType = sys.argv[2]
        else:
            iTest = None
    else:
        testType = 'full'
        iTest = None

    if testType == 'simple':
        simple([iTest])
    elif testType == 'Bayes':
        _, _ = TestBayes('Test')
    else:
        full(iTest, skipType=skipType)
