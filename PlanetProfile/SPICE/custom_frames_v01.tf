KPL/FK


Custom Dynamic Frames Kernel
========================================================================

    This FK contains dynamic frames definitions for modeling and analysis
    of MAG data from multiple missions to the outer planets.


Version and Date
========================================================================

    Version 1.0 -- July 20, 2022 -- Marshall J. Styczinski, JPL
    
      Added Planet-Dipole-Solar-Zenith (PDSZ) frame for each giant
      planet for use in evaluating magnetopause current and field
      models.
    
    Version 0.9 -- July 20, 2022 -- Marshall J. Styczinski, JPL
   
      Initial version; includes MPhiO frames for the Galilean moons,
      as well as Planet-Solar-Orbital and Planet-Solar-Magnetospheric
      (PSO and PSM) frames for each giant planet.


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
       Kivelson, M.G., 2018. Saturn’s magnetic field revealed by the Cassini
       Grand Finale. Science, 362(6410), p.eaat5434.
       https://doi.org/10.1126/science.aat5434
    
    8. AH5 model for Uranus:
       Herbert, F., 2009. Aurora and magnetic field of Uranus. Journal of
       Geophysical Research: Space Physics, 114(A11).
       https://doi.org/10.1029/2009JA014394
    
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
                                    
      JSO                     J2000                  1850065
      KSO                     J2000                  1850066
      USO                     J2000                  1850067
      NSO                     J2000                  1850068
      
      Planet-Solar-Orbital frames:
        -- +X axis is along the geometric position of the Sun as seen
           from Jupiter (primary axis)
        -- +Y axis is in the direction of the inertial geometric velocity
           of the Sun as seen from Jupiter
        -- centered on the planet
         
        Name                  Relative to            NAIF ID
    ======================  =====================  ============
    
      JSM                     J2000                  1850075
      KSM                     J2000                  1850076
      USM                     J2000                  1850077
      NSM                     J2000                  1850078
      
      Planet-Solar-Magnetospheric frames:
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
            

\begindata

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
    FRAME_1850075_PRI_ABCORR      = 'NONE'
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
    FRAME_1850076_PRI_ABCORR      = 'NONE'
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
    FRAME_1850077_PRI_ABCORR      = 'NONE'
    FRAME_1850077_SEC_AXIS        = 'Z'
    FRAME_1850077_SEC_VECTOR_DEF  = 'CONSTANT'
    FRAME_1850077_SEC_FRAME       = 'IAU_URANUS'
    FRAME_1850077_SEC_SPEC        = 'LATITUDINAL'
    FRAME_1850077_SEC_UNITS       = 'DEGREES'
    FRAME_1850077_SEC_LONGITUDE   = -55.7
    FRAME_1850077_SEC_LATITUDE    =  30.2
    
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
    FRAME_1850078_PRI_ABCORR      = 'NONE'
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
    FRAME_1850085_PRI_ABCORR      = 'NONE'
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
    FRAME_1850086_PRI_ABCORR      = 'NONE'
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
    FRAME_1850087_PRI_ABCORR      = 'NONE'
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
    FRAME_1850088_PRI_ABCORR      = 'NONE'
    FRAME_1850088_SEC_AXIS        = 'X'
    FRAME_1850088_SEC_VECTOR_DEF  = 'CONSTANT'
    FRAME_1850088_SEC_FRAME       = 'NSM'
    FRAME_1850088_SEC_SPEC        = 'RECTANGULAR'
    FRAME_1850088_SEC_VECTOR      = ( 0, 0, -1 )

\begintext

End of FK.
