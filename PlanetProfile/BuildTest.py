""" This is a test script meant to check that the current PlanetProfile
    build maintains all functionality across updates. Run this script
    as "python -m PlanetProfile.BuildTest" from the main PlanetProfile directory
    before committing any major updates. New major functionality added
    to PlanetProfile should be accompanied by a new PPTest#.py test
    body in the Test/ directory.
"""

import logging
import numpy as np
import importlib, os, fnmatch, sys, time
from copy import deepcopy
from PlanetProfile import _Test, _TestImport
from PlanetProfile.GetConfig import Params as configParams, FigMisc
from PlanetProfile.Main import PlanetProfile, InductOgram, ReloadInductOgram, ExploreOgram, ReloadExploreOgram
from PlanetProfile.Plotting.ExplorationPlots import GenerateExplorationPlots, PlotExploreOgramMultiSubplot
from PlanetProfile.Plotting.MagPlots import GenerateMagPlots, GenerateExplorationMagPlots
from PlanetProfile.Test.TestBayes import TestBayes
from PlanetProfile.Utilities.defineStructs import EOSlist, EOSlistStruct

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
    Params = setFullSettings(Params)
    

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
            try:
                TestPlanets = np.append(TestPlanets, PlanetProfile(testPlanetN, Params)[0])
                tMarks = np.append(tMarks, time.time())
            except Exception as e:
                log.error(f"Error running PlanetProfile for test {i}: {e}")
                print(f"Error in test {i}: {e}")
                raise  # Stop execution on error

            # Verify that we can reload things as needed in each case
            Params.CALC_NEW = False
            Params.CALC_NEW_REF = False
            Params.CALC_NEW_INDUCT = False
            Params.CALC_NEW_GRAVITY = False
            Params.CALC_NEW_ASYM = False
            TestPlanets = np.append(TestPlanets, PlanetProfile(deepcopy(testPlanetN), Params)[0])
            TestPlanets[-1].saveLabel += ' RELOAD'
            tMarks = np.append(tMarks, time.time())

            Params.CALC_NEW = True
            Params.CALC_NEW_REF = True
            Params.CALC_NEW_INDUCT = True
            Params.CALC_NEW_GRAVITY = True
            Params.CALC_NEW_ASYM = False

        testPlanet1.name = 'Test0'
        # Test that we can successfully run standard profiles with parallelization options
        Params.DO_PARALLEL = True
        TestPlanets = np.append(TestPlanets, PlanetProfile(deepcopy(testPlanet1), Params)[0])
        TestPlanets[-1].saveLabel += '  PARALLEL'
        tMarks = np.append(tMarks, time.time())

        # Make sure our auxiliary calculation flags work correctly
        Params.CALC_SEISMIC = False
        Params.CALC_CONDUCT = False
        Params.CALC_VISCOSITY = False
        Params.CALC_OCEAN_PROPS = False
        TestPlanets = np.append(TestPlanets, PlanetProfile(deepcopy(testPlanet1), Params)[0])
        TestPlanets[-1].saveLabel += ' NO_SEISMIC_OR_CONDUCT_OR_VISCOSITY'
        tMarks = np.append(tMarks, time.time())
        Params.CALC_SEISMIC = True
        Params.CALC_CONDUCT = True
        Params.CALC_VISCOSITY = True
        Params.CALC_OCEAN_PROPS = True
    """ Bayes testing disabled for now
    if skipType is None or skipType.lower() == 'bayes':
        # Test Bayesian analysis UpdateRun capabilities
        PlanetBayes, _ = TestBayes('Test')
        PlanetBayes.saveLabel = 'Bayes'
        TestPlanets = np.append(TestPlanets, PlanetBayes)
        tMarks = np.append(tMarks, time.time())
    """

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
    Params.CALC_NEW_INDUCT = True
    Params.CALC_NEW = True
    Params.NO_SAVEFILE = True
    Params.SKIP_INNER = True
    Params.Induct.oceanCompList = ['Seawater', 'MgSO4', 'NaCl', 'PureH2O']
    for inductOtype in ['sigma', 'Tb', 'rho', 'oceanComp']:
        Params.Induct.inductOtype = inductOtype
        # Test non-parallel and parallel processing
        Params.DO_PARALLEL = False
        Params.PRELOAD_EOS = False
        _ = TestInductOgram(7, Params)
        Params.DO_PARALLEL = True
        Params.PRELOAD_EOS = True
        Induction = TestInductOgram(7, Params)
        TestPlanets = np.append(TestPlanets, deepcopy(Induction))
        tMarks = np.append(tMarks, time.time())
        # Test reloading
        Induction = TestInductOgram(7, Params, CALC_NEW=False)
        TestPlanets = np.append(TestPlanets, deepcopy(Induction))
        tMarks = np.append(tMarks, time.time())

    for inductOtype in ['phi']:
        Params.Induct.inductOtype = inductOtype
        Params.DO_PARALLEL = False
        Params.PRELOAD_EOS = False
        _ = TestInductOgram(11, Params)
        Params.DO_PARALLEL = True
        Params.PRELOAD_EOS = True
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
    Params.CALC_NEW = True
    Params.NO_SAVEFILE = True
    Params.Explore.zName = 'CMR2mean'
    Params.Explore.oceanCompRangeList = ['Seawater', 'MgSO4', 'NaCl', 'PureH2O']
    allZNames = ["CMR2mean", "D_km", "Dconv_m", "dzIceI_km", "dzClath_km", "dzIceIII_km", "dzIceIIIund_km",
    "dzIceV_km", "dzIceVund_km", "dzIceVI_km", "dzWetHPs_km", "eLid_km", "phiSeafloor_frac", "Rcore_km", "rhoSilMean_kgm3", "rhoCoreMean_kgm3",
    "sigmaMean_Sm", "silPhiCalc_frac", "zb_km", "zSeafloor_km",  
    "hLoveAmp", "kLoveAmp", "lLoveAmp", "deltaLoveAmp", "hLovePhase", "kLovePhase", "lLovePhase", "deltaLovePhase",
    'InductionAmp', 'InductionPhase', 'InductionrBi1Tot_nT', 'InductioniBi1Tot_nT', 'InductionrBi1x_nT', 'InductionrBi1y_nT', 
    'InductionrBi1z_nT', 'InductioniBi1x_nT', 'InductioniBi1y_nT', 'InductioniBi1z_nT']
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
        'Qrad_Wkg': [1e-20, 1e-14],
        'oceanComp': [-12, -3], # Placeholder - range is determined by input range file data,
        'zb_approximate_km': [10, 90]
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
            global EOSlist
            EOSlist = EOSlistStruct()
            for yName in hydroExploreBds.keys():
                if xName != yName:
                    Params.Explore.xName = xName
                    Params.Explore.yName = yName
                    Params.Explore.xRange = hydroExploreBds[xName]
                    Params.Explore.yRange = hydroExploreBds[yName]
                    Params.DO_PARALLEL = False
                    Params.PRELOAD_EOS = False
                    _ = TestExploreOgram(7, Params)
                    Params.DO_PARALLEL = True
                    Params.PRELOAD_EOS = True
                    Exploration = TestExploreOgram(7, Params)
                    TestPlanets = np.append(TestPlanets, deepcopy(Exploration))
                    tMarks = np.append(tMarks, time.time())

                    Exploration = TestExploreOgram(7, Params, CALC_NEW=False)
                    TestPlanets = np.append(TestPlanets, deepcopy(Exploration))
                    tMarks = np.append(tMarks, time.time())
        
        log.info('Testing all zName options for exploreOgrams.')
        Exploration = TestExploreOgram(7, Params, CALC_NEW=False)
        # Now test the multi-subplot function if not already tested (dividing into groups of 10 to prevent errors in too big of file size)
        for i in range(0, len(allZNames), 10):
            Params.Explore.zName = allZNames[i:i+10]
            Exploration = TestExploreOgram(7, Params, CALC_NEW=False)
    
    log.info('Running exploreOgrams for waterless bodies for all input types.')
    Params.Explore.zName = 'CMR2mean'
    for xName in waterlessExploreBds.keys():
        for yName in waterlessExploreBds.keys():
            if xName != yName:
                Params.Explore.xName = xName
                Params.Explore.yName = yName
                Params.Explore.xRange = waterlessExploreBds[xName]
                Params.Explore.yRange = waterlessExploreBds[yName]
                Params.DO_PARALLEL = False
                Params.PRELOAD_EOS = False
                _ = TestExploreOgram(5, Params)
                Params.DO_PARALLEL = True
                Params.PRELOAD_EOS = True
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
    Params.Induct.nrhoPts, Params.Induct.nSigmaPts, Params.Induct.nDpts, Params.Induct.nOceanCompPts, Params.Induct.nZbPts \
        = (4 for _ in range(8))
    
    if CALC_NEW:
        Induction, Params = InductOgram(testName, Params)
        end = ''
    else:
        Induction, _, Params = ReloadInductOgram(testName, Params)
        end = ' RELOAD'
    GenerateExplorationMagPlots([Induction], [Params.FigureFiles], Params)
    Induction.name = testName
    Induction.saveLabel = f'{Params.Induct.inductOtype} induct-o-gram{end}'

    return Induction


def TestExploreOgram(testNum, Params, CALC_NEW=True):
    testName = f'Test{testNum}'

    # Set sizes low so things don't take ages to run
    Params.Explore.nx, Params.Explore.ny \
        = (4 for _ in range(2))

    if CALC_NEW:
        Params.CALC_NEW = True
        Exploration, Params = ExploreOgram(testName, Params)
        end = ''
    else:
        Params.CALC_NEW = False
        Exploration, Params = ReloadExploreOgram(testName, Params)
        end = ' RELOAD'

    GenerateExplorationPlots([Exploration], [Params.FigureFiles], Params)
    Exploration.name = testName
    Exploration.saveLabel = f'{Params.Explore.xName} x {Params.Explore.yName} explore-o-gram{end}'

    return Exploration


def simple(iTests=None):
    testMod = f'{_TestImport}.PPTest'

    # Set general testing config atop standard config options
    Params = configParams
    Params = setFullSettings(Params)

    tStart = time.time()
    # Normalize iTests into a list
    if isinstance(iTests, int):
        iTests = [iTests]
    else:
        iTests = list(iTests) 
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
            # Verify that we can reload things as needed in each case
            Params.CALC_NEW = False
            Params.CALC_NEW_REF = False
            Params.CALC_NEW_INDUCT = False
            _ = PlanetProfile(testPlanet, Params)
            Params.CALC_NEW = True
            Params.CALC_NEW_REF = True
            Params.CALC_NEW_INDUCT = True
    tEnd = time.time()
    log.info('Simple test complete!')

    log.debug(f'Elapsed time:\n    {tEnd - tStart:.1f} s')
    return


def setFullSettings(Params):
    Params.CALC_NEW = True
    Params.CALC_NEW_REF = True
    Params.CALC_NEW_INDUCT = True
    Params.CALC_NEW_GRAVITY = True
    Params.CALC_NEW_ASYM = False
    Params.CALC_SEISMIC = True
    Params.CALC_CONDUCT = True
    Params.CALC_VISCOSITY = True
    Params.CALC_OCEAN_PROPS = True
    Params.CALC_ASYM = True
    Params.DO_PARALLEL = False
    Params.RUN_ALL_PROFILES = False
    Params.COMPARE = False
    Params.NO_SAVEFILE = False
    Params.FORCE_EOS_RECALC = False
    Params.SKIP_INNER = False
    Params.DISP_LAYERS = True
    Params.DISP_TABLE  = True
    Params.ALLOW_BROKEN_MODELS = False
    Params.DEPRECATED = False
    Params.TIME_AND_DATE_LABEL = False
    # Set plotting options
    Params.SKIP_PLOTS = False
    Params.PLOT_GRAVITY = True
    Params.PLOT_HYDROSPHERE = True
    Params.PLOT_HYDROSPHERE_THERMODYNAMICS = True
    Params.PLOT_MELTING_CURVES = True
    Params.PLOT_SPECIES_HYDROSPHERE = True
    Params.PLOT_REF = True
    Params.PLOT_SIGS = True
    Params.PLOT_SOUNDS = True
    Params.PLOT_REF = True
    Params.PLOT_SIGS = True
    Params.PLOT_SOUNDS = True
    Params.PLOT_TRADEOFF = True
    Params.PLOT_POROSITY = True
    Params.PLOT_SEISMIC = True
    Params.PLOT_PRESSURE_DEPTH = True
    Params.PLOT_VISCOSITY = True
    Params.PLOT_WEDGE = True
    Params.PLOT_HYDRO_PHASE = True
    Params.PLOT_PVT_HYDRO = True
    Params.PLOT_PVT_ISOTHERMAL_HYDRO = True
    Params.PLOT_PVT_INNER = True
    Params.PLOT_BDIP = True
    Params.PLOT_BSURF = True
    Params.PLOT_ASYM = True
    Params.PLOT_TRAJECS = True
    Params.PLOT_BINVERSION = True
    Params.LEGEND = True
    Params.TITLES = True
    Params.PLOT_COMBO_EXPLORATIONS = True
    FigMisc.propsToPlot = ['rho', 'Cp', 'alpha', 'VP', 'KS', 'sig', 'VS', 'GS']
    FigMisc.PmaxHydro_MPa = None
    FigMisc.TmaxHydro_K = None
    FigMisc.TminHydro_K = None
    FigMisc.SHOW_GEOTHERM = True
    
    # Disable other features not relevant to single model runs
    Params.DO_INDUCTOGRAM = False
    Params.DO_EXPLOREOGRAM = False
    Params.DO_MONTECARLO = False
    Params.SKIP_INDUCTION = False
    Params.SKIP_GRAVITY = False
    Params.PRELOAD_EOS = False
    
    return Params
    
if __name__ == '__main__':
    skipType = None
    if len(sys.argv) > 1:
        # Test type was passed as command line argument
        testType = sys.argv[1]
        if len(sys.argv) == 3:
            if sys.argv[2].isdigit():
                # Test profile number was passed as command line argument
                iTest = int(sys.argv[2])
            else:
                # Skip to specific testing section was passed as command line arg
                iTest = None
                skipType = sys.argv[2]
        elif len(sys.argv) > 3:
            iTest = []
            for arg in sys.argv[2:]:
                if arg.isdigit():
                    iTest.append(int(arg))
        else:
            iTest = None
    else:
        testType = 'full'
        iTest = None

    if testType == 'simple':
        simple(iTest)
    elif testType == 'Bayes':
        _, _ = TestBayes('Test')
    else:
        full(iTest, skipType=skipType)
