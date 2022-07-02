KPL/FK

\beginlabel
PDS_VERSION_ID               = PDS3
RECORD_TYPE                  = STREAM
RECORD_BYTES                 = "N/A"
^SPICE_KERNEL                = "cas_dyn_v03.tf"
MISSION_NAME                 = "CASSINI-HUYGENS"
SPACECRAFT_NAME              = "CASSINI ORBITER"
DATA_SET_ID                  = "CO-S/J/E/V-SPICE-6-V1.0"
KERNEL_TYPE_ID               = FK
PRODUCT_ID                   = "cas_dyn_v03.tf"
PRODUCT_CREATION_TIME        = 2018-09-27T10:35:23
PRODUCER_ID                  = "NAIF/JPL"
MISSION_PHASE_NAME           = "N/A"
PRODUCT_VERSION_TYPE         = ACTUAL
PLATFORM_OR_MOUNTING_NAME    = "N/A"
START_TIME                   = "N/A"
STOP_TIME                    = "N/A"
SPACECRAFT_CLOCK_START_COUNT = "N/A"
SPACECRAFT_CLOCK_STOP_COUNT  = "N/A"
TARGET_NAME                  = "N/A"
INSTRUMENT_NAME              = "N/A"
NAIF_INSTRUMENT_ID           = "N/A"
SOURCE_PRODUCT_ID            = "N/A"
NOTE                         = "See comments in the file for details"
OBJECT                       = SPICE_KERNEL
  INTERCHANGE_FORMAT         = ASCII
  KERNEL_TYPE                = FRAMES
  DESCRIPTION                = "SPICE Frames Kernel file containing
definitions of a few dynamic frames used by various project science and
engineering teams. "
END_OBJECT                   = SPICE_KERNEL
\endlabel


CASSINI Dynamic Frame Definitions
===========================================================================

   This kernel contains dynamic frame definitions to support CASSINI
   mission.


Version and Date
---------------------------------------------------------------

   Version 3.0 -- September 20, 2018 -- Boris Semenov

      Added CASSINI_RINGS_PLANNING, and CASSINI_RINGS_SHA frames.

   Version 2.0 -- July 10, 2013 -- Boris Semenov

      Added CASSINI_KRTP, CASSINI_KSM, and CASSINI_KSO frames.

   Version 1.0 -- October 27, 2010 -- Boris Semenov

      Initial release. Includes only CASSINI_RTN frame.


References
---------------------------------------------------------------

      1. "Kernel Pool Required Reading"

      2. "Frames Required Reading"

      3. http://mapsview.engin.umich.edu/data_descriptions/
         coordinate_systems.php

      4. CO_MAG_CAL_1SEC_DS.CAT

      5. "Rings Into SPICE, Phase 1", Jeff Boyer, 2017-09-13,
         SPICEyRings-PSG73-ws.docx
 
      6. "Solar Hour Angle (SHA) Code Sets", Jeff Boyer: MSAUL-278, A.I.,
         2017-10-20, CAS_SHA_CodeSets-JSB_MSAUL-278-AI.pdf


Contact Information
---------------------------------------------------------------

   Direct questions, comments, or concerns about the contents of this 
   kernel to:

      Boris Semenov, NAIF/JPL, Boris.Semenov@jpl.nasa.gov


Implementation Notes
---------------------------------------------------------------

   This file is used by the SPICE system as follows: programs that make
   use of this frame kernel must "load" the kernel using SPICE's FURNSH
   routine, normally during program initialization.  Loading the kernel
   associates the data items with their names in a data structure
   called the "kernel pool". This file was created and may be updated
   with a text editor or word processor.


Frames Defined in This Kernel
---------------------------------------------------------------

   The following frames are defined in this kernel file:
 
      Frame Name                Relative To              Type     NAIF ID
      =======================   ===================      =======  =======
 
      Dynamic Frames:
      ------------------
      CASSINI_RTN               J2000                    DYNAMIC   -82902
      CASSINI_KRTP              J2000                    DYNAMIC   -82903
      CASSINI_KSM               J2000                    DYNAMIC   -82904
      CASSINI_KSO               J2000                    DYNAMIC   -82905
      CASSINI_RINGS_PLANNING    J2000                    DYNAMIC   -82906
      CASSINI_RINGS_SHA         J2000                    DYNAMIC   -82907


CASSINI RTN Frame
---------------------------------------------------------------

   CASSINI RTN Frame (CASSINI_RTN) is defined as follows:

      R is positive from the Sun to the spacecraft.
      T is omega cross R, where omega is the sun spin axis.
      N is R cross T, which completes the right-handed system.

   This frame assumes the instantaneous center is located at the Sun,
   and not the CASSINI spacecraft.  Further, the axes are associated
   with the normal X, Y, and Z in the following manner:

      R -> X
      T -> Y
      N -> Z

   This set of keywords defines the CASSINI_RTN frame.

      \begindata

         FRAME_CASSINI_RTN            = -82902
         FRAME_-82902_NAME            = 'CASSINI_RTN'
         FRAME_-82902_CLASS           = 5
         FRAME_-82902_CLASS_ID        = -82902
         FRAME_-82902_CENTER          = 10
         FRAME_-82902_RELATIVE        = 'J2000'
         FRAME_-82902_DEF_STYLE       = 'PARAMETERIZED'
         FRAME_-82902_FAMILY          = 'TWO-VECTOR'
         FRAME_-82902_PRI_AXIS        = 'X'
         FRAME_-82902_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
         FRAME_-82902_PRI_OBSERVER    = 'SUN'
         FRAME_-82902_PRI_TARGET      = 'CASSINI'
         FRAME_-82902_PRI_ABCORR      = 'NONE'
         FRAME_-82902_SEC_AXIS        = 'Z'
         FRAME_-82902_SEC_VECTOR_DEF  = 'CONSTANT'
         FRAME_-82902_SEC_SPEC        = 'RECTANGULAR'
         FRAME_-82902_SEC_VECTOR      = ( 0.0, 0.0, 1.0 )
         FRAME_-82902_SEC_FRAME       = 'IAU_SUN'

      \begintext 


CASSINI_KRTP
---------------------------------------------------------------

   The CASSINI KRTP frame is defined in [4] as follows:

      Kronocentric body-fixed, J2000 spherical Coordinates (KRTP)
      -----------------------------------------------------------
      KRTP magnetic field vector components form the standard right-handed
      spherical triad (R, Theta, Phi) for a planet-centered system. 
      Namely, R is radial (along the line from the center of Saturn to the
      center of the spacecraft), and positive away from Saturn. Phi, the 
      azimuthal component, is parallel to the Kronographic equator (Omega
      x R) and positive in the direction of corotation. Theta, the 
      'southward' component, completes the right-handed set.
   
   This frame assumes the instantaneous center is located at Saturn,
   and not the CASSINI spacecraft.  Further, the axes are associated
   with the normal X, Y, and Z in the following manner:

      R     -> X
      Theta -> Y
      Phi   -> Z

   This set of keywords defines the CASSINI_KRTP frame.

      \begindata

         FRAME_CASSINI_KRTP           = -82903
         FRAME_-82903_NAME            = 'CASSINI_KRTP'
         FRAME_-82903_CLASS           = 5
         FRAME_-82903_CLASS_ID        = -82903
         FRAME_-82903_CENTER          = 699
         FRAME_-82903_RELATIVE        = 'J2000'
         FRAME_-82903_DEF_STYLE       = 'PARAMETERIZED'
         FRAME_-82903_FAMILY          = 'TWO-VECTOR'
         FRAME_-82903_PRI_AXIS        = 'X'
         FRAME_-82903_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
         FRAME_-82903_PRI_OBSERVER    = 'SATURN'
         FRAME_-82903_PRI_TARGET      = 'CASSINI'
         FRAME_-82903_PRI_ABCORR      = 'NONE'
         FRAME_-82903_SEC_AXIS        = 'Y'
         FRAME_-82903_SEC_VECTOR_DEF  = 'CONSTANT'
         FRAME_-82903_SEC_SPEC        = 'RECTANGULAR'
         FRAME_-82903_SEC_VECTOR      = ( 0.0, 0.0, -1.0 )
         FRAME_-82903_SEC_FRAME       = 'IAU_SATURN'

      \begintext 


CASSINI_KSM
---------------------------------------------------------------

   The CASSINI KSM frame is defined in [4] as follows:

      Kronocentric Solar Magnetospheric Coordinates (KSM)
      ---------------------------------------------------
      KSM is a cartesian Saturn-center coordinate system where X points 
      from Saturn to the Sun, the X-Z plane contains Saturn's centered 
      magnetic dipole axis, M, and Y completes right handed set.

   For consistency with the custom software at the PPI Node of PDS
   the Saturn's centered magnetic dipole axis is set to 180 degrees 
   longitude, 89.99 degrees latitude in the IAU_SATURN frame. Note that
   other source (e.g. [3]) make assume that the dipole axis is 
   parallel to the spin axis.
   
   This set of keywords defines the CASSINI_KSM frame.

      \begindata

         FRAME_CASSINI_KSM            = -82904
         FRAME_-82904_NAME            = 'CASSINI_KSM'
         FRAME_-82904_CLASS           = 5
         FRAME_-82904_CLASS_ID        = -82904
         FRAME_-82904_CENTER          = 699
         FRAME_-82904_RELATIVE        = 'J2000'
         FRAME_-82904_DEF_STYLE       = 'PARAMETERIZED'
         FRAME_-82904_FAMILY          = 'TWO-VECTOR'
         FRAME_-82904_PRI_AXIS        = 'X'
         FRAME_-82904_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
         FRAME_-82904_PRI_OBSERVER    = 'SATURN'
         FRAME_-82904_PRI_TARGET      = 'SUN'
         FRAME_-82904_PRI_ABCORR      = 'NONE'
         FRAME_-82904_SEC_AXIS        = 'Z'
         FRAME_-82904_SEC_VECTOR_DEF  = 'CONSTANT'
         FRAME_-82904_SEC_SPEC        = 'LATITUDINAL'
         FRAME_-82904_SEC_UNITS       = 'DEGREES'
         FRAME_-82904_SEC_LONGITUDE   = 180.00
         FRAME_-82904_SEC_LATITUDE    =  89.99
         FRAME_-82904_SEC_FRAME       = 'IAU_SATURN'

      \begintext 


CASSINI_KSO
---------------------------------------------------------------

   The CASSINI KSO frame is defined in [4] as follows:

       Kronocentric Solar Orbital Coordinates (KSO)
       --------------------------------------------
       KSO is a cartesian Saturn-centered coordinate system where X points
       from Saturn to the Sun, Z is parallel to Saturn's orbital plane 
       upward normal, and Y completes the right handed set.

   This set of keywords defines the CASSINI_KSO frame.

      \begindata

         FRAME_CASSINI_KSO            = -82905
         FRAME_-82905_NAME            = 'CASSINI_KSO'
         FRAME_-82905_CLASS           = 5
         FRAME_-82905_CLASS_ID        = -82905
         FRAME_-82905_CENTER          = 699
         FRAME_-82905_RELATIVE        = 'J2000'
         FRAME_-82905_DEF_STYLE       = 'PARAMETERIZED'
         FRAME_-82905_FAMILY          = 'TWO-VECTOR'
         FRAME_-82905_PRI_AXIS        = 'X'
         FRAME_-82905_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
         FRAME_-82905_PRI_OBSERVER    = 'SATURN' 
         FRAME_-82905_PRI_TARGET      = 'SUN'
         FRAME_-82905_PRI_ABCORR      = 'NONE'
         FRAME_-82905_SEC_AXIS        = 'Y'
         FRAME_-82905_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_VELOCITY'
         FRAME_-82905_SEC_OBSERVER    = 'SATURN' 
         FRAME_-82905_SEC_TARGET      = 'SUN'
         FRAME_-82905_SEC_ABCORR      = 'NONE'
         FRAME_-82905_SEC_FRAME       = 'J2000'

      \begintext 


CASSINI_RINGS_PLANNING
---------------------------------------------------------------

   The CASSINI_RINGS_PLANNING frame facilitates calculation of the ring
   longitude used in the CASSINI project planning/targeting tools.

   The CASSINI_RINGS_PLANNING frame is defined as follows (per section
   1 in [5]):

      +Z     is along the Saturn pole (+Z of IAU_SATURN.)

      +Y     is in the direction of the J2000 pole (+Z of J2000.) 

      +X     completes the right-handed frame (is along the ascending
             node of the Saturn equator on the J2000 equator..)

      origin is at the center of Saturn.

   The ring longitude of a point is its planetographic (west-positive)
   longitude in this frame, computed by RECPGR for BODY='SATURN'.

   This set of keywords defines the CASSINI_RINGS_PLANNING frame.

      \begindata

         FRAME_CASSINI_RINGS_PLANNING = -82906
         FRAME_-82906_NAME            = 'CASSINI_RINGS_PLANNING'
         FRAME_-82906_CLASS           = 5
         FRAME_-82906_CLASS_ID        = -82906
         FRAME_-82906_CENTER          = 699
         FRAME_-82906_RELATIVE        = 'J2000'
         FRAME_-82906_DEF_STYLE       = 'PARAMETERIZED'
         FRAME_-82906_FAMILY          = 'TWO-VECTOR'
         FRAME_-82906_PRI_AXIS        = 'Z'
         FRAME_-82906_PRI_VECTOR_DEF  = 'CONSTANT'
         FRAME_-82906_PRI_SPEC        = 'RECTANGULAR'
         FRAME_-82906_PRI_VECTOR      = ( 0.0, 0.0, 1.0 )
         FRAME_-82906_PRI_FRAME       = 'IAU_SATURN'
         FRAME_-82906_SEC_AXIS        = 'Y'
         FRAME_-82906_SEC_VECTOR_DEF  = 'CONSTANT'
         FRAME_-82906_SEC_SPEC        = 'RECTANGULAR'
         FRAME_-82906_SEC_VECTOR      = ( 0.0, 0.0, 1.0 )
         FRAME_-82906_SEC_FRAME       = 'J2000'

      \begintext 


CASSINI_RINGS_SHA
---------------------------------------------------------------

   The CASSINI_RINGS_SHA frame facilitates calculation of the solar
   hour angle (SHA) used in the CASSINI project planning/targeting
   tools.

   The CASSINI_RINGS_SHA frame is defined as follows (based on PDT
   Calculate_Solar_Hour_Angle function from [6]):

      +Z     is along the Saturn pole (+Z of IAU_SATURN.)

      -X     is in the direction of the apparent position of the Sun. 

      +Y     completes the right-handed frame.

      origin is at the center of Saturn.

   The solar hour angle of a point is its planetographic
   (west-positive) longitude in this frame, computed by RECPGR for
   BODY='SATURN'.

   This set of keywords defines the CASSINI_RINGS_SHA frame.

      \begindata

         FRAME_CASSINI_RINGS_SHA      = -82907
         FRAME_-82907_NAME            = 'CASSINI_RINGS_SHA'
         FRAME_-82907_CLASS           = 5
         FRAME_-82907_CLASS_ID        = -82907
         FRAME_-82907_CENTER          = 699
         FRAME_-82907_RELATIVE        = 'J2000'
         FRAME_-82907_DEF_STYLE       = 'PARAMETERIZED'
         FRAME_-82907_FAMILY          = 'TWO-VECTOR'
         FRAME_-82907_PRI_AXIS        = 'Z'
         FRAME_-82907_PRI_VECTOR_DEF  = 'CONSTANT'
         FRAME_-82907_PRI_SPEC        = 'RECTANGULAR'
         FRAME_-82907_PRI_VECTOR      = ( 0.0, 0.0, 1.0 )
         FRAME_-82907_PRI_FRAME       = 'IAU_SATURN'
         FRAME_-82907_SEC_AXIS        = '-X'
         FRAME_-82907_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
         FRAME_-82907_SEC_OBSERVER    = 'SATURN'
         FRAME_-82907_SEC_TARGET      = 'SUN'
         FRAME_-82907_SEC_ABCORR      = 'LT+S'

      \begintext 

End of FK file.
