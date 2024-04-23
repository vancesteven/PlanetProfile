""" Default trajectory analysis settings """
from PlanetProfile.Utilities.defineStructs import TrajecParamsStruct

configTrajecVersion = 1  # Integer number for config file version. Increment when new settings are added to the default config file.

def trajecAssign():
    Trajec = TrajecParamsStruct()

    Trajec.targetBody = 'Europa'  # Body for which to analyze spacecraft data to invert interior structure
    Trajec.scSelect = ['Galileo', 'Juno']  # List of spacecraft to include in analysis
    Trajec.plasmaType = 'Alfven'  # A string describing the type of large-scale plasma model to apply. Options are 'Alfven' and 'none'.
    Trajec.trajecAppend = ''  # Custom string to use to save/reload specific settings
    Trajec.MAGdir = 'SpacecraftMAGdata'  # Directory where spacecraft magnetic data is stored
    Trajec.FORCE_MAG_RECALC = False  # Whether to read in MAG data from disk and regenerate reformatted HDF5 version.
    Trajec.EXPANDED_RANGE = False  # Whether to plot an expanded set of B measurements farther from the encounter CA, with range set by etExpandRange_s
    Trajec.PLANETMAG_MODEL = False  # Whether to load in evaluated magnetic field models printed to disk from PlanetMag instead of directly evaluating excitation moments

    Trajec.fbInclude = {  # Set to 'all' or a list of strings of encounter ID numbers
        'Cassini': 'all',
        'Clipper': 'all',
        'Galileo': 'all',
        'Juno': 'all',
        'JUICE': 'all'
    }
    Trajec.fbRange_Rp = 0.5  # Maximum distance in planetary radii (of parent planet, for moons) within which to mark flyby encounters
    Trajec.etPredRange_s = 5400  # Range in seconds for the span across closest approach to use for predicted spacecraft trajectories, i.e. those for which we do not yet have data
    Trajec.etExpandRange_s = 18000  # Range in seconds for the span across closest approach to use for expanded spacecraft trajectories, i.e. beyond the main PDS files near CA
    Trajec.etStep_s = 1.0  # Step size to use in spanning the above
    Trajec.fbDescrip = {  # String to prepend to flyby ID number based on how each mission labels orbits
        'Cassini': 'Rev ',
        'Clipper': Trajec.targetBody[0],
        'Galileo': Trajec.targetBody[0],
        'Juno': 'PJ',
        'JUICE': Trajec.targetBody[0]
    }
    Trajec.spiceSPK = {
        'Cassini': [
            'cas_2004_v26.tm',
            'cas_2005_v27.tm',
            'cas_2006_v26.tm',
            'cas_2007_v24.tm',
            'cas_2008_v22.tm',
            'cas_2009_v19.tm',
            'cas_2010_v18.tm',
            'cas_2011_v17.tm',
            'cas_2012_v14.tm',
            'cas_2013_v12.tm',
            'cas_2014_v08.tm',
            'cas_2015_v08.tm',
            'cas_2016_v06.tm',
            'cas_2017_v03.tm'
        ],
        'Clipper': ['Rnd7_T3_pad_scpse.bsp'],
        'Galileo': ['s980326a.bsp', 's000131a.bsp', 's030916a.bsp'],
        'JUICE': [],
        'Juno': ['juno_pred_orbit.bsp', 'juno_rec_orbit.bsp'],
        'Voyager 1': ['vgr1_jup230.bsp', 'vgr1_sat337.bsp'],
        'Voyager 2': ['vgr2_jup230.bsp', 'vgr2_sat337.bsp', 'vgr2.ura111.bsp', 'vgr2_nep097.bsp']
    }
    Trajec.nFitParamsGlobal = 3  # Number of fit parameters that persist across flybys, e.g. ocean salinity or ice melting temp
    Trajec.nFitParamsFlybys = 3  # Number of fit parameters that vary with each flyby, e.g. Alfven wave characteristics

    Trajec.SCera = None  # Spacecraft era to use for excitation moments. If None, the sc name will be used to select the era.
    Trajec.BextModel = None  # Planetary magnetic field from which to use excitation moments derived using PlanetMag for induction calculations. If None, use the model from the defaults listed in the class definition.

    return Trajec
