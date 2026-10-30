KPL/FK


Custom Dynamic Frames Kernel
========================================================================

    This FK contains dynamic frames definitions for modeling and analysis
    of MAG data from multiple missions to the outer planets.


Version and Date
========================================================================

    Version 1.4 -- November 27, 2023 -- Marshall J. Styczinski, JPL
    
      Added Earth-Sun-Magnetic and Earth-Sun-Orbital frames akin to
      those for the giant planets.

    Version 1.3 -- August 28, 2022 -- Marshall J. Styczinski, JPL
    
      Added NLS frame for Neptune, used in magnetic field modeling.

    Version 1.2 -- August 27, 2022 -- Marshall J. Styczinski, JPL
    
      Added ULS frame for Uranus, used in magnetic field modeling.
    
    Version 1.1 -- July 24, 2022 -- Marshall J. Styczinski, JPL
    
      Added US3 frame, identical to IAU_URANUS but with the z axis
      along the angular momentum vector instead of opposite to it.
      Added Solar-Magnetic-Planet frames for each giant planet.

    Version 1.0 -- July 20, 2022 -- Marshall J. Styczinski, JPL
    
      Added Planet-Dipole-Solar-Zenith (PDSZ) frame for each giant
      planet for use in evaluating magnetopause current and field
      models.
    
    Version 0.9 -- July 20, 2022 -- Marshall J. Styczinski, JPL
   
      Initial version; includes MPhiO frames for the Galilean moons,
      as well as Planet-Sun-Orbital and Planet-Sun-Magnetic (PSO and
      PSM) frames for each giant planet.


References
========================================================================

    1. ``Frames Required Reading''

    2. ``Kernel Pool Required Reading''

    3. https://lasp.colorado.edu/home/mop/files/2015/02/CoOrd_systems12.pdf
   
    4. https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=GO-J-POS-6-SC-TRAJ-JUP-COORDS-V1.0
   
    5. https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=GO-J-POS-6-SC-TRAJ-MOON-COORDS-V1.0
    
    6. O4 model for Jupiter:
       Acuna, M.H. and Ness, N.F., 1976. The main magnetic field of Jupiter.
       Journal of Geophysical Research, 81(16), pp.2917-2922.
       https://doi.org/10.1029/JA081i016p02917
    
    7. Cassini 11 model for Saturn:
       Dougherty, M.K., Cao, H., Khurana, K.K., Hunt, G.J., Provan, G.,
       Kellock, S., Burton, M.E., Burk, T.A., Bunce, E.J., Cowley, S.W. and
       Kivelson, M.G., 2018. Saturn's magnetic field revealed by the Cassini
       Grand Finale. Science, 362(6410), p.eaat5434.
       https://doi.org/10.1126/science.aat5434
    
    8. OTD model for Uranus:
       Ness, N.F., Acuna, M.H., Behannon, K.W., Burlaga, L.F., Connerney, J.E.,
       Lepping, R.P. and Neubauer, F.M., 1986. Magnetic fields at Uranus.
       Science, 233(4759), pp.85-89.
       https://doi.org/10.1126/science.233.4759.85
    
    9. O8 model for Neptune:
       Connerney, J.E.P., Acuna, M.H. and Ness, N.F., 1991. The magnetic
       field of Neptune. Journal of Geophysical Research: Space Physics,
       96(S01), pp.19023-19042.
       https://doi.org/10.1029/91JA01165
    
    10. Definition of what is now called Jupiter-Solar-Magnetospheric
        frame:
        Connerney, J.E.P., Acuna, M.H. and Ness, N.F., 1981. Modeling the
        Jovian current sheet and inner magnetosphere. Journal of
        Geophysical Research: Space Physics, 86(A10), pp.8370-8384.
        https://doi.org/10.1029/JA086iA10p08370
        
    11. Justification of 0.01 degree off-pole, 180 degrees longitude
        orientation of KSM model:
        cas_dyn_v03.tf
    
    12. Alexeev, I.I. and Belenkaya, E.S., 2005, March. Modeling of the
        Jovian magnetosphere. In Annales Geophysicae (Vol. 23, No. 3,
        pp. 809-826). Copernicus GmbH.
        https://doi.org/10.5194/angeo-23-809-2005


Frames Definitions
========================================================================

    This FK currently defines the following dynamic frames of use:
    
        Name                  Relative to            NAIF ID
    ======================  =====================  ============

      US3                     J2000                  1850007
      ULS                     J2000                  1850017
      NLS                     J2000                  1850018
      NLS_RADEC               J2000                  1850028
            
      US3 frame:
        -- +Z is aligned with the Uranus spin pole (primary axis)
        -- +X is in the direction of the IAU prime meridian
        -- centered on the planet
      
      ULS frame:
        -- +Z is aligned with the Uranus spin pole (primary axis)
        -- 302 degrees W longitude is at the location of Voyager 2
           closest approach, as defined in [8]. The Voyager 2 CA
           is located at 167.3113 degrees E longitude in the
           IAU_URANUS frame, at the time the reconstructed trajectory
           minimizes the distance to the Uranus center of mass.
           This results in a shift of -134.6887 degrees E longitude,
           i.e. +X of ULS points toward IAU E longitude 225.3113.
        -- centered on the planet
      
      NLS frame:
        -- +Z is aligned with the Neptune spin pole (primary axis),
           defined in [9] to be RA 298.90 degrees, DEC 42.84 degrees.
           These values are assumed to be in reference to J2000 or ICRF.
        -- 167.7 degrees W longitude is at the location of Voyager 2
           near closest approach, at 0356 SCET day 237 as defined in [9].
           Voyager 2 was at -155.6860 degrees E longitude in the
           IAU_NEPTUNE frame at this time, according to the VG2 SPK
           generated using the nep097.bsp generic kernel, and evaluated
           using that same kernel.
           This results in a shift of 12.0140 degrees E longitude,
           i.e. +X of NLS points toward IAU E longitude 12.0140.
        -- centered on the planet
      
      NLS_RADEC frame:
        -- Exactly as NLS, but using the RA/DEC values defined in [9].
           The spin pole is displaced about 125 km along the 1-bar
           surface of the planet.
    
        Name                  Relative to            NAIF ID
    ======================  =====================  ============

      IO_PHI_O                IAU_JUPITER            1859501
      EUROPA_PHI_O            IAU_JUPITER            1859502
      GANYMEDE_PHI_O          IAU_JUPITER            1859503
      CALLISTO_PHI_O          IAU_JUPITER            1859504
      
      MOON_PHI_O frames:
        -- +Z is aligned with the Jupiter spin pole (primary axis)
        -- +Y is in the direction of the component of the Moon-Jupiter
           vector perpendicular to Z (secondary axis directed at Jupiter)
        -- centered on the moon
         
        Name                  Relative to            NAIF ID
    ======================  =====================  ============
                                    
      ESO                     J2000                  1850063
      JSO                     J2000                  1850065
      KSO                     J2000                  1850066
      USO                     J2000                  1850067
      NSO                     J2000                  1850068
      
      Planet-Sun-Orbit frames:
        -- +X axis is along the geometric position of the Sun as seen
           from Jupiter (primary axis)
        -- +Y axis is in the direction of the inertial geometric velocity
           of the Sun as seen from Jupiter
        -- centered on the planet
         
        Name                  Relative to            NAIF ID
    ======================  =====================  ============
    
      ESM                     J2000                  1850073
      JSM                     J2000                  1850075
      KSM                     J2000                  1850076
      USM                     J2000                  1850077
      NSM                     J2000                  1850078
      
      Planet-Sun-Magnetic frames:
        -- as defined in [3], with dipole orientations as defined in [6-9]
        -- originally referred to as Magnetic Equatorial coordinates
           in [10], with the O4 model for Jupiter.
        -- +X axis is along the geometric position of the Sun as seen
           from Jupiter
        -- +Y axis is along M x X, where M is the direction of the
           magnetic dipole moment as defined in the models referenced
           above
        -- centered on the planet
         
        Name                  Relative to            NAIF ID
    ======================  =====================  ============
            
      JDSZ                    J2000                  1850085
      KDSZ                    J2000                  1850086
      UDSZ                    J2000                  1850087
      NDSZ                    J2000                  1850088
      
      Planet-Dipole-Solar-Zenith frames:
        -- Cartesian form of the spherical coordinates defined in [12]
           (but not with this name)
        -- +Z axis is along the geometric position of the Sun as seen
           from Jupiter
        -- +phi direction is positive counterclockwise (right-handed)
           around +Z direction, starting from 0 at the XZ-plane of the
           PSM frame, i.e. opposite the direction of the magnetic dipole
           moment as defined in the models referenced in [6-9]
        -- The above makes the +Y axis along M x Z, where M is the
           direction of the magnetic dipole moment as defined in the
           models referenced in [6-9]
        -- centered on the planet
         
        Name                  Relative to            NAIF ID
    ======================  =====================  ============
            
      SMJ                     J2000                  1850095
      SMK                     J2000                  1850096
      SMU                     J2000                  1850097
      SMN                     J2000                  1850098
      
      Solar-Magnetic-Planet frames:
        -- Not to be confused with Planet-Solar-Magnetospheric frames
        -- +Z axis is along the dipole axis as defined in the models
           referenced in [6-9]
        -- +phi direction is positive (right-handed) around +Z axis,
           starting from 0 at the planet-Sun direction
        -- centered on the planet
            

\begindata

    FRAME_US3                      = 1850007
    FRAME_1850007_NAME             = 'US3'
    FRAME_1850007_CLASS            = 5
    FRAME_1850007_CLASS_ID         = 1850007
    FRAME_1850007_CENTER           = 799
    FRAME_1850007_RELATIVE         = 'J2000'
    FRAME_1850007_DEF_STYLE        = 'PARAMETERIZED'
    FRAME_1850007_FAMILY           = 'TWO-VECTOR'
    FRAME_1850007_PRI_AXIS         = 'Z'
    FRAME_1850007_PRI_VECTOR_DEF   = 'CONSTANT'
    FRAME_1850007_PRI_FRAME        = 'IAU_URANUS'
    FRAME_1850007_PRI_SPEC         = 'RECTANGULAR'
    FRAME_1850007_PRI_VECTOR       = ( 0, 0, -1 )
    FRAME_1850007_SEC_AXIS         = 'X'
    FRAME_1850007_SEC_VECTOR_DEF   = 'CONSTANT'
    FRAME_1850007_SEC_FRAME        = 'IAU_URANUS'
    FRAME_1850007_SEC_SPEC         = 'RECTANGULAR'
    FRAME_1850007_SEC_VECTOR       = ( 1, 0, 0 )

    FRAME_ULS                      = 1850017
    FRAME_1850017_NAME             = 'ULS'
    FRAME_1850017_CLASS            = 5
    FRAME_1850017_CLASS_ID         = 1850017
    FRAME_1850017_CENTER           = 799
    FRAME_1850017_RELATIVE         = 'J2000'
    FRAME_1850017_DEF_STYLE        = 'PARAMETERIZED'
    FRAME_1850017_FAMILY           = 'TWO-VECTOR'
    FRAME_1850017_PRI_AXIS         = 'Z'
    FRAME_1850017_PRI_VECTOR_DEF   = 'CONSTANT'
    FRAME_1850017_PRI_FRAME        = 'IAU_URANUS'
    FRAME_1850017_PRI_SPEC         = 'RECTANGULAR'
    FRAME_1850017_PRI_VECTOR       = ( 0, 0, -1 )
    FRAME_1850017_SEC_AXIS         = 'X'
    FRAME_1850017_SEC_VECTOR_DEF   = 'CONSTANT'
    FRAME_1850017_SEC_FRAME        = 'IAU_URANUS'
    FRAME_1850017_SEC_SPEC         = 'LATITUDINAL'
    FRAME_1850017_SEC_UNITS        = 'DEGREES'
    FRAME_1850017_SEC_LONGITUDE    = 225.3113
    FRAME_1850017_SEC_LATITUDE     =  0.0

    FRAME_NLS                      = 1850018
    FRAME_1850018_NAME             = 'NLS'
    FRAME_1850018_CLASS            = 5
    FRAME_1850018_CLASS_ID         = 1850018
    FRAME_1850018_CENTER           = 899
    FRAME_1850018_RELATIVE         = 'J2000'
    FRAME_1850018_DEF_STYLE        = 'PARAMETERIZED'
    FRAME_1850018_FAMILY           = 'TWO-VECTOR'
    FRAME_1850018_PRI_AXIS         = 'Z'
    FRAME_1850018_PRI_VECTOR_DEF   = 'CONSTANT'
    FRAME_1850018_PRI_FRAME        = 'IAU_NEPTUNE'
    FRAME_1850018_PRI_SPEC         = 'RECTANGULAR'
    FRAME_1850018_PRI_VECTOR       = ( 0, 0, 1 )
    FRAME_1850018_SEC_AXIS         = 'X'
    FRAME_1850018_SEC_VECTOR_DEF   = 'CONSTANT'
    FRAME_1850018_SEC_FRAME        = 'IAU_NEPTUNE'
    FRAME_1850018_SEC_SPEC         = 'LATITUDINAL'
    FRAME_1850018_SEC_UNITS        = 'DEGREES'
    FRAME_1850018_SEC_LONGITUDE    = 12.0140
    FRAME_1850018_SEC_LATITUDE     =  0.0

    FRAME_NLS_RADEC                = 1850028
    FRAME_1850028_NAME             = 'NLS_RADEC'
    FRAME_1850028_CLASS            = 5
    FRAME_1850028_CLASS_ID         = 1850028
    FRAME_1850028_CENTER           = 899
    FRAME_1850028_RELATIVE         = 'J2000'
    FRAME_1850028_DEF_STYLE        = 'PARAMETERIZED'
    FRAME_1850028_FAMILY           = 'TWO-VECTOR'
    FRAME_1850028_PRI_AXIS         = 'Z'
    FRAME_1850028_PRI_VECTOR_DEF   = 'CONSTANT'
    FRAME_1850028_PRI_SPEC         = 'RA/DEC'
    FRAME_1850028_PRI_FRAME        = 'J2000'
    FRAME_1850028_PRI_UNITS        = 'DEGREES'
    FRAME_1850028_PRI_RA           = 298.90
    FRAME_1850028_PRI_DEC          = 42.84
    FRAME_1850028_SEC_AXIS         = 'X'
    FRAME_1850028_SEC_VECTOR_DEF   = 'CONSTANT'
    FRAME_1850028_SEC_FRAME        = 'IAU_NEPTUNE'
    FRAME_1850028_SEC_SPEC         = 'LATITUDINAL'
    FRAME_1850028_SEC_UNITS        = 'DEGREES'
    FRAME_1850028_SEC_LONGITUDE    = 12.0152
    FRAME_1850028_SEC_LATITUDE     =  0.0
        
    FRAME_IO_PHI_O                 = 1859501
    FRAME_1859501_NAME             = 'IO_PHI_O'
    FRAME_1859501_CLASS            = 5
    FRAME_1859501_CLASS_ID         = 1859501
    FRAME_1859501_CENTER           = 501
    FRAME_1859501_RELATIVE         = 'J2000'
    FRAME_1859501_DEF_STYLE        = 'PARAMETERIZED'
    FRAME_1859501_FAMILY           = 'TWO-VECTOR'
    FRAME_1859501_PRI_AXIS         = 'Z'
    FRAME_1859501_PRI_VECTOR_DEF   = 'CONSTANT'
    FRAME_1859501_PRI_FRAME        = 'IAU_JUPITER'
    FRAME_1859501_PRI_SPEC         = 'RECTANGULAR'
    FRAME_1859501_PRI_VECTOR       = ( 0, 0, 1 )
    FRAME_1859501_SEC_AXIS         = 'Y'
    FRAME_1859501_SEC_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
    FRAME_1859501_SEC_OBSERVER     = 'IO'
    FRAME_1859501_SEC_TARGET       = 'JUPITER'
    FRAME_1859501_SEC_ABCORR       = 'NONE'

    FRAME_EUROPA_PHI_O             = 1859502
    FRAME_1859502_NAME             = 'EUROPA_PHI_O'
    FRAME_1859502_CLASS            = 5
    FRAME_1859502_CLASS_ID         = 1859502
    FRAME_1859502_CENTER           = 502
    FRAME_1859502_RELATIVE         = 'J2000'
    FRAME_1859502_DEF_STYLE        = 'PARAMETERIZED'
    FRAME_1859502_FAMILY           = 'TWO-VECTOR'
    FRAME_1859502_PRI_AXIS         = 'Z'
    FRAME_1859502_PRI_VECTOR_DEF   = 'CONSTANT'
    FRAME_1859502_PRI_FRAME        = 'IAU_JUPITER'
    FRAME_1859502_PRI_SPEC         = 'RECTANGULAR'
    FRAME_1859502_PRI_VECTOR       = ( 0, 0, 1 )
    FRAME_1859502_SEC_AXIS         = 'Y'
    FRAME_1859502_SEC_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
    FRAME_1859502_SEC_OBSERVER     = 'EUROPA'
    FRAME_1859502_SEC_TARGET       = 'JUPITER'
    FRAME_1859502_SEC_ABCORR       = 'NONE'

    FRAME_GANYMEDE_PHI_O           = 1859503
    FRAME_1859503_NAME             = 'GANYMEDE_PHI_O'
    FRAME_1859503_CLASS            = 5
    FRAME_1859503_CLASS_ID         = 1859503
    FRAME_1859503_CENTER           = 503
    FRAME_1859503_RELATIVE         = 'J2000'
    FRAME_1859503_DEF_STYLE        = 'PARAMETERIZED'
    FRAME_1859503_FAMILY           = 'TWO-VECTOR'
    FRAME_1859503_PRI_AXIS         = 'Z'
    FRAME_1859503_PRI_VECTOR_DEF   = 'CONSTANT'
    FRAME_1859503_PRI_FRAME        = 'IAU_JUPITER'
    FRAME_1859503_PRI_SPEC         = 'RECTANGULAR'
    FRAME_1859503_PRI_VECTOR       = ( 0, 0, 1 )
    FRAME_1859503_SEC_AXIS         = 'Y'
    FRAME_1859503_SEC_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
    FRAME_1859503_SEC_OBSERVER     = 'GANYMEDE'
    FRAME_1859503_SEC_TARGET       = 'JUPITER'
    FRAME_1859503_SEC_ABCORR       = 'NONE'

    FRAME_CALLISTO_PHI_O           = 1859504
    FRAME_1859504_NAME             = 'CALLISTO_PHI_O'
    FRAME_1859504_CLASS            = 5
    FRAME_1859504_CLASS_ID         = 1859504
    FRAME_1859504_CENTER           = 504
    FRAME_1859504_RELATIVE         = 'J2000'
    FRAME_1859504_DEF_STYLE        = 'PARAMETERIZED'
    FRAME_1859504_FAMILY           = 'TWO-VECTOR'
    FRAME_1859504_PRI_AXIS         = 'Z'
    FRAME_1859504_PRI_VECTOR_DEF   = 'CONSTANT'
    FRAME_1859504_PRI_FRAME        = 'IAU_JUPITER'
    FRAME_1859504_PRI_SPEC         = 'RECTANGULAR'
    FRAME_1859504_PRI_VECTOR       = ( 0, 0, 1 )
    FRAME_1859504_SEC_AXIS         = 'Y'
    FRAME_1859504_SEC_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
    FRAME_1859504_SEC_OBSERVER     = 'CALLISTO'
    FRAME_1859504_SEC_TARGET       = 'JUPITER'
    FRAME_1859504_SEC_ABCORR       = 'NONE'
      
    FRAME_ESO                     = 1850063
    FRAME_1850063_NAME            = 'ESO'
    FRAME_1850063_CLASS           = 5
    FRAME_1850063_CLASS_ID        = 1850063
    FRAME_1850063_CENTER          = 399
    FRAME_1850063_RELATIVE        = 'J2000'
    FRAME_1850063_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850063_FAMILY          = 'TWO-VECTOR'
    FRAME_1850063_PRI_AXIS        = 'X'
    FRAME_1850063_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850063_PRI_OBSERVER    = 'EARTH'
    FRAME_1850063_PRI_TARGET      = 'SUN'
    FRAME_1850063_PRI_ABCORR      = 'NONE'
    FRAME_1850063_SEC_AXIS        = 'Y'
    FRAME_1850063_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_VELOCITY'
    FRAME_1850063_SEC_OBSERVER    = 'EARTH'
    FRAME_1850063_SEC_TARGET      = 'SUN'
    FRAME_1850063_SEC_ABCORR      = 'NONE'
    FRAME_1850063_SEC_FRAME       = 'J2000'
      
    FRAME_JSO                     = 1850065
    FRAME_1850065_NAME            = 'JSO'
    FRAME_1850065_CLASS           = 5
    FRAME_1850065_CLASS_ID        = 1850065
    FRAME_1850065_CENTER          = 599
    FRAME_1850065_RELATIVE        = 'J2000'
    FRAME_1850065_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850065_FAMILY          = 'TWO-VECTOR'
    FRAME_1850065_PRI_AXIS        = 'X'
    FRAME_1850065_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850065_PRI_OBSERVER    = 'JUPITER'
    FRAME_1850065_PRI_TARGET      = 'SUN'
    FRAME_1850065_PRI_ABCORR      = 'NONE'
    FRAME_1850065_SEC_AXIS        = 'Y'
    FRAME_1850065_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_VELOCITY'
    FRAME_1850065_SEC_OBSERVER    = 'JUPITER'
    FRAME_1850065_SEC_TARGET      = 'SUN'
    FRAME_1850065_SEC_ABCORR      = 'NONE'
    FRAME_1850065_SEC_FRAME       = 'J2000'
      
    FRAME_KSO                     = 1850066
    FRAME_1850066_NAME            = 'KSO'
    FRAME_1850066_CLASS           = 5
    FRAME_1850066_CLASS_ID        = 1850066
    FRAME_1850066_CENTER          = 699
    FRAME_1850066_RELATIVE        = 'J2000'
    FRAME_1850066_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850066_FAMILY          = 'TWO-VECTOR'
    FRAME_1850066_PRI_AXIS        = 'X'
    FRAME_1850066_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850066_PRI_OBSERVER    = 'SATURN'
    FRAME_1850066_PRI_TARGET      = 'SUN'
    FRAME_1850066_PRI_ABCORR      = 'NONE'
    FRAME_1850066_SEC_AXIS        = 'Y'
    FRAME_1850066_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_VELOCITY'
    FRAME_1850066_SEC_OBSERVER    = 'SATURN'
    FRAME_1850066_SEC_TARGET      = 'SUN'
    FRAME_1850066_SEC_ABCORR      = 'NONE'
    FRAME_1850066_SEC_FRAME       = 'J2000'
      
    FRAME_USO                     = 1850067
    FRAME_1850067_NAME            = 'USO'
    FRAME_1850067_CLASS           = 5
    FRAME_1850067_CLASS_ID        = 1850067
    FRAME_1850067_CENTER          = 799
    FRAME_1850067_RELATIVE        = 'J2000'
    FRAME_1850067_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850067_FAMILY          = 'TWO-VECTOR'
    FRAME_1850067_PRI_AXIS        = 'X'
    FRAME_1850067_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850067_PRI_OBSERVER    = 'URANUS'
    FRAME_1850067_PRI_TARGET      = 'SUN'
    FRAME_1850067_PRI_ABCORR      = 'NONE'
    FRAME_1850067_SEC_AXIS        = 'Y'
    FRAME_1850067_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_VELOCITY'
    FRAME_1850067_SEC_OBSERVER    = 'URANUS'
    FRAME_1850067_SEC_TARGET      = 'SUN'
    FRAME_1850067_SEC_ABCORR      = 'NONE'
    FRAME_1850067_SEC_FRAME       = 'J2000'
      
    FRAME_NSO                     = 1850068
    FRAME_1850068_NAME            = 'NSO'
    FRAME_1850068_CLASS           = 5
    FRAME_1850068_CLASS_ID        = 1850068
    FRAME_1850068_CENTER          = 899
    FRAME_1850068_RELATIVE        = 'J2000'
    FRAME_1850068_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850068_FAMILY          = 'TWO-VECTOR'
    FRAME_1850068_PRI_AXIS        = 'X'
    FRAME_1850068_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850068_PRI_OBSERVER    = 'NEPTUNE'
    FRAME_1850068_PRI_TARGET      = 'SUN'
    FRAME_1850068_PRI_ABCORR      = 'NONE'
    FRAME_1850068_SEC_AXIS        = 'Y'
    FRAME_1850068_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_VELOCITY'
    FRAME_1850068_SEC_OBSERVER    = 'NEPTUNE'
    FRAME_1850068_SEC_TARGET      = 'SUN'
    FRAME_1850068_SEC_ABCORR      = 'NONE'
    FRAME_1850068_SEC_FRAME       = 'J2000'
      
    FRAME_ESM                     = 1850073
    FRAME_1850073_NAME            = 'ESM'
    FRAME_1850073_CLASS           = 5
    FRAME_1850073_CLASS_ID        = 1850073
    FRAME_1850073_CENTER          = 399
    FRAME_1850073_RELATIVE        = 'J2000'
    FRAME_1850073_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850073_FAMILY          = 'TWO-VECTOR'
    FRAME_1850073_PRI_AXIS        = 'X'
    FRAME_1850073_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850073_PRI_OBSERVER    = 'EARTH'
    FRAME_1850073_PRI_TARGET      = 'SUN'
    FRAME_1850073_PRI_ABCORR      = 'LT+S'
    FRAME_1850073_SEC_AXIS        = 'Z'
    FRAME_1850073_SEC_VECTOR_DEF  = 'CONSTANT'
    FRAME_1850073_SEC_FRAME       = 'IAU_EARTH'
    FRAME_1850073_SEC_SPEC        = 'LATITUDINAL'
    FRAME_1850073_SEC_UNITS       = 'DEGREES'
    FRAME_1850073_SEC_LONGITUDE   = 107.32
    FRAME_1850073_SEC_LATITUDE    = -80.59
      
    FRAME_JSM                     = 1850075
    FRAME_1850075_NAME            = 'JSM'
    FRAME_1850075_CLASS           = 5
    FRAME_1850075_CLASS_ID        = 1850075
    FRAME_1850075_CENTER          = 599
    FRAME_1850075_RELATIVE        = 'J2000'
    FRAME_1850075_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850075_FAMILY          = 'TWO-VECTOR'
    FRAME_1850075_PRI_AXIS        = 'X'
    FRAME_1850075_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850075_PRI_OBSERVER    = 'JUPITER'
    FRAME_1850075_PRI_TARGET      = 'SUN'
    FRAME_1850075_PRI_ABCORR      = 'LT+S'
    FRAME_1850075_SEC_AXIS        = 'Z'
    FRAME_1850075_SEC_VECTOR_DEF  = 'CONSTANT'
    FRAME_1850075_SEC_FRAME       = 'IAU_JUPITER'
    FRAME_1850075_SEC_SPEC        = 'LATITUDINAL'
    FRAME_1850075_SEC_UNITS       = 'DEGREES'
    FRAME_1850075_SEC_LONGITUDE   = 158.0
    FRAME_1850075_SEC_LATITUDE    =  80.4
    
    FRAME_KSM                     = 1850076
    FRAME_1850076_NAME            = 'KSM'
    FRAME_1850076_CLASS           = 5
    FRAME_1850076_CLASS_ID        = 1850076
    FRAME_1850076_CENTER          = 699
    FRAME_1850076_RELATIVE        = 'J2000'
    FRAME_1850076_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850076_FAMILY          = 'TWO-VECTOR'
    FRAME_1850076_PRI_AXIS        = 'X'
    FRAME_1850076_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850076_PRI_OBSERVER    = 'SATURN'
    FRAME_1850076_PRI_TARGET      = 'SUN'
    FRAME_1850076_PRI_ABCORR      = 'LT+S'
    FRAME_1850076_SEC_AXIS        = 'Z'
    FRAME_1850076_SEC_VECTOR_DEF  = 'CONSTANT'
    FRAME_1850076_SEC_FRAME       = 'IAU_SATURN'
    FRAME_1850076_SEC_SPEC        = 'LATITUDINAL'
    FRAME_1850076_SEC_UNITS       = 'DEGREES'
    FRAME_1850076_SEC_LONGITUDE   = 180.00
    FRAME_1850076_SEC_LATITUDE    =  89.99
    
    FRAME_USM                     = 1850077
    FRAME_1850077_NAME            = 'USM'
    FRAME_1850077_CLASS           = 5
    FRAME_1850077_CLASS_ID        = 1850077
    FRAME_1850077_CENTER          = 799
    FRAME_1850077_RELATIVE        = 'J2000'
    FRAME_1850077_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850077_FAMILY          = 'TWO-VECTOR'
    FRAME_1850077_PRI_AXIS        = 'X'
    FRAME_1850077_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850077_PRI_OBSERVER    = 'URANUS'
    FRAME_1850077_PRI_TARGET      = 'SUN'
    FRAME_1850077_PRI_ABCORR      = 'LT+S'
    FRAME_1850077_SEC_AXIS        = 'Z'
    FRAME_1850077_SEC_VECTOR_DEF  = 'CONSTANT'
    FRAME_1850077_SEC_FRAME       = 'US3'
    FRAME_1850077_SEC_SPEC        = 'LATITUDINAL'
    FRAME_1850077_SEC_UNITS       = 'DEGREES'
    FRAME_1850077_SEC_LONGITUDE   = -48.0
    FRAME_1850077_SEC_LATITUDE    =  30.0
    
    FRAME_NSM                     = 1850078
    FRAME_1850078_NAME            = 'NSM'
    FRAME_1850078_CLASS           = 5
    FRAME_1850078_CLASS_ID        = 1850078
    FRAME_1850078_CENTER          = 899
    FRAME_1850078_RELATIVE        = 'J2000'
    FRAME_1850078_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850078_FAMILY          = 'TWO-VECTOR'
    FRAME_1850078_PRI_AXIS        = 'X'
    FRAME_1850078_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850078_PRI_OBSERVER    = 'NEPTUNE'
    FRAME_1850078_PRI_TARGET      = 'SUN'
    FRAME_1850078_PRI_ABCORR      = 'LT+S'
    FRAME_1850078_SEC_AXIS        = 'Z'
    FRAME_1850078_SEC_VECTOR_DEF  = 'CONSTANT'
    FRAME_1850078_SEC_FRAME       = 'IAU_NEPTUNE'
    FRAME_1850078_SEC_SPEC        = 'LATITUDINAL'
    FRAME_1850078_SEC_UNITS       = 'DEGREES'
    FRAME_1850078_SEC_LONGITUDE   = -71.96
    FRAME_1850078_SEC_LATITUDE    =  43.10
      
    FRAME_JDSZ                     = 1850085
    FRAME_1850085_NAME            = 'JDSZ'
    FRAME_1850085_CLASS           = 5
    FRAME_1850085_CLASS_ID        = 1850085
    FRAME_1850085_CENTER          = 599
    FRAME_1850085_RELATIVE        = 'J2000'
    FRAME_1850085_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850085_FAMILY          = 'TWO-VECTOR'
    FRAME_1850085_PRI_AXIS        = 'Z'
    FRAME_1850085_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850085_PRI_OBSERVER    = 'JUPITER'
    FRAME_1850085_PRI_TARGET      = 'SUN'
    FRAME_1850085_PRI_ABCORR      = 'LT+S'
    FRAME_1850085_SEC_AXIS         = 'X'
    FRAME_1850085_SEC_VECTOR_DEF   = 'CONSTANT'
    FRAME_1850085_SEC_FRAME        = 'JSM'
    FRAME_1850085_SEC_SPEC         = 'RECTANGULAR'
    FRAME_1850085_SEC_VECTOR       = ( 0, 0, -1 )
    
    FRAME_KDSZ                     = 1850086
    FRAME_1850086_NAME            = 'KDSZ'
    FRAME_1850086_CLASS           = 5
    FRAME_1850086_CLASS_ID        = 1850086
    FRAME_1850086_CENTER          = 699
    FRAME_1850086_RELATIVE        = 'J2000'
    FRAME_1850086_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850086_FAMILY          = 'TWO-VECTOR'
    FRAME_1850086_PRI_AXIS        = 'Z'
    FRAME_1850086_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850086_PRI_OBSERVER    = 'SATURN'
    FRAME_1850086_PRI_TARGET      = 'SUN'
    FRAME_1850086_PRI_ABCORR      = 'LT+S'
    FRAME_1850086_SEC_AXIS        = 'X'
    FRAME_1850086_SEC_VECTOR_DEF  = 'CONSTANT'
    FRAME_1850086_SEC_FRAME       = 'KSM'
    FRAME_1850086_SEC_SPEC        = 'RECTANGULAR'
    FRAME_1850086_SEC_VECTOR      = ( 0, 0, -1 )
    
    FRAME_UDSZ                     = 1850087
    FRAME_1850087_NAME            = 'UDSZ'
    FRAME_1850087_CLASS           = 5
    FRAME_1850087_CLASS_ID        = 1850087
    FRAME_1850087_CENTER          = 799
    FRAME_1850087_RELATIVE        = 'J2000'
    FRAME_1850087_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850087_FAMILY          = 'TWO-VECTOR'
    FRAME_1850087_PRI_AXIS        = 'Z'
    FRAME_1850087_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850087_PRI_OBSERVER    = 'URANUS'
    FRAME_1850087_PRI_TARGET      = 'SUN'
    FRAME_1850087_PRI_ABCORR      = 'LT+S'
    FRAME_1850087_SEC_AXIS        = 'X'
    FRAME_1850087_SEC_VECTOR_DEF  = 'CONSTANT'
    FRAME_1850087_SEC_FRAME       = 'USM'
    FRAME_1850087_SEC_SPEC        = 'RECTANGULAR'
    FRAME_1850087_SEC_VECTOR      = ( 0, 0, -1 )
    
    FRAME_NDSZ                    = 1850088
    FRAME_1850088_NAME            = 'NDSZ'
    FRAME_1850088_CLASS           = 5
    FRAME_1850088_CLASS_ID        = 1850088
    FRAME_1850088_CENTER          = 899
    FRAME_1850088_RELATIVE        = 'J2000'
    FRAME_1850088_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850088_FAMILY          = 'TWO-VECTOR'
    FRAME_1850088_PRI_AXIS        = 'Z'
    FRAME_1850088_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850088_PRI_OBSERVER    = 'NEPTUNE'
    FRAME_1850088_PRI_TARGET      = 'SUN'
    FRAME_1850088_PRI_ABCORR      = 'LT+S'
    FRAME_1850088_SEC_AXIS        = 'X'
    FRAME_1850088_SEC_VECTOR_DEF  = 'CONSTANT'
    FRAME_1850088_SEC_FRAME       = 'NSM'
    FRAME_1850088_SEC_SPEC        = 'RECTANGULAR'
    FRAME_1850088_SEC_VECTOR      = ( 0, 0, -1 )
      
    FRAME_SMJ                     = 1850095
    FRAME_1850095_NAME            = 'SMJ'
    FRAME_1850095_CLASS           = 5
    FRAME_1850095_CLASS_ID        = 1850095
    FRAME_1850095_CENTER          = 599
    FRAME_1850095_RELATIVE        = 'J2000'
    FRAME_1850095_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850095_FAMILY          = 'TWO-VECTOR'
    FRAME_1850095_PRI_AXIS        = 'Z'
    FRAME_1850095_PRI_VECTOR_DEF  = 'CONSTANT'
    FRAME_1850095_PRI_FRAME       = 'IAU_JUPITER'
    FRAME_1850095_PRI_SPEC        = 'LATITUDINAL'
    FRAME_1850095_PRI_UNITS       = 'DEGREES'
    FRAME_1850095_PRI_LONGITUDE   = 158.0
    FRAME_1850095_PRI_LATITUDE    =  80.4
    FRAME_1850095_SEC_AXIS        = 'X'
    FRAME_1850095_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850095_SEC_OBSERVER    = 'JUPITER'
    FRAME_1850095_SEC_TARGET      = 'SUN'
    FRAME_1850095_SEC_ABCORR      = 'LT+S'
    
    FRAME_SMK                     = 1850096
    FRAME_1850096_NAME            = 'SMK'
    FRAME_1850096_CLASS           = 5
    FRAME_1850096_CLASS_ID        = 1850096
    FRAME_1850096_CENTER          = 699
    FRAME_1850096_RELATIVE        = 'J2000'
    FRAME_1850096_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850096_FAMILY          = 'TWO-VECTOR'
    FRAME_1850096_PRI_AXIS        = 'Z'
    FRAME_1850096_PRI_VECTOR_DEF  = 'CONSTANT'
    FRAME_1850096_PRI_FRAME       = 'IAU_SATURN'
    FRAME_1850096_PRI_SPEC        = 'LATITUDINAL'
    FRAME_1850096_PRI_UNITS       = 'DEGREES'
    FRAME_1850096_PRI_LONGITUDE   = 180.00
    FRAME_1850096_PRI_LATITUDE    =  89.99
    FRAME_1850096_SEC_AXIS        = 'X'
    FRAME_1850096_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850096_SEC_OBSERVER    = 'SATURN'
    FRAME_1850096_SEC_TARGET      = 'SUN'
    FRAME_1850096_SEC_ABCORR      = 'LT+S'
    
    FRAME_SMU                     = 1850097
    FRAME_1850097_NAME            = 'SMU'
    FRAME_1850097_CLASS           = 5
    FRAME_1850097_CLASS_ID        = 1850097
    FRAME_1850097_CENTER          = 799
    FRAME_1850097_RELATIVE        = 'J2000'
    FRAME_1850097_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850097_FAMILY          = 'TWO-VECTOR'
    FRAME_1850097_PRI_AXIS        = 'Z'
    FRAME_1850097_PRI_VECTOR_DEF  = 'CONSTANT'
    FRAME_1850097_PRI_FRAME       = 'IAU_URANUS'
    FRAME_1850097_PRI_SPEC        = 'LATITUDINAL'
    FRAME_1850097_PRI_UNITS       = 'DEGREES'
    FRAME_1850097_PRI_LONGITUDE   =  48.0
    FRAME_1850097_PRI_LATITUDE    = -30.0
    FRAME_1850097_SEC_AXIS        = 'X'
    FRAME_1850097_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850097_SEC_OBSERVER    = 'URANUS'
    FRAME_1850097_SEC_TARGET      = 'SUN'
    FRAME_1850097_SEC_ABCORR      = 'LT+S'
    
    FRAME_SMN                     = 1850098
    FRAME_1850098_NAME            = 'SMN'
    FRAME_1850098_CLASS           = 5
    FRAME_1850098_CLASS_ID        = 1850098
    FRAME_1850098_CENTER          = 899
    FRAME_1850098_RELATIVE        = 'J2000'
    FRAME_1850098_DEF_STYLE       = 'PARAMETERIZED'
    FRAME_1850098_FAMILY          = 'TWO-VECTOR'
    FRAME_1850098_PRI_AXIS        = 'Z'
    FRAME_1850098_PRI_VECTOR_DEF  = 'CONSTANT'
    FRAME_1850098_PRI_FRAME       = 'IAU_NEPTUNE'
    FRAME_1850098_PRI_SPEC        = 'LATITUDINAL'
    FRAME_1850098_PRI_UNITS       = 'DEGREES'
    FRAME_1850098_PRI_LONGITUDE   = -71.96
    FRAME_1850098_PRI_LATITUDE    =  43.10
    FRAME_1850098_SEC_AXIS        = 'X'
    FRAME_1850098_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
    FRAME_1850098_SEC_OBSERVER    = 'NEPTUNE'
    FRAME_1850098_SEC_TARGET      = 'SUN'
    FRAME_1850098_SEC_ABCORR      = 'LT+S'

\begintext

End of FK.
