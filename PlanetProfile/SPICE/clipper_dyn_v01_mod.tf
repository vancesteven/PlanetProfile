KPL/FK


EUROPAM Dynamic Frames Kernel
========================================================================

   This is a placeholder FK for Europa Clipper dynamic frames.


Version and Date
========================================================================

   Version 0.1_mod -- February 20, 2022 -- Marshall J. Styczinski, JPL
   
      Edited the EPhiO definition to be named as EUROPA_PHI_O and
      added analogous versions for Io, Ganymede, and Callisto.
   
   Version 0.1 -- December 18, 2019 -- Boris Semenov, NAIF

      Added the EUROPAM_EUROPA_E_PHI_O frame.

   Version 0.0 -- November 15, 2017 -- Boris Semenov, NAIF

      Initial version; includes the EUROPAM_JUP_POLE_EUROPA and 
      EUROPAM_JUP_EUROPA_POLE frames.


References
========================================================================

   1. ``Frames Required Reading''

   2. ``Kernel Pool Required Reading''

   3. E-mail ``RE: spacecraft frame kernel'' from Corey Cochrane, JPL,
      from November 15, 2017

   4. https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=GO-J-POS-6-SC-TRAJ-MOON-COORDS-V1.0


Frames Definitions
========================================================================

   This FK currently defines six dynamic frames of use:

      EUROPAM_JUP_POLE_EUROPA (-159951), per [3]:

         -- +Z is aligned with the Jupiter pole (primary axis) 
         -- +X is in the direction Jupiter-Europa vector (secondary axis)
         -- centered on Jupiter

      EUROPAM_JUP_EUROPA_POLE (-159952), per [3]:

         -- +X is aligned with the Jupiter-Europa vector (primary axis)
         -- +Z is in the direction of Jupiter pole (secondary axis)
         -- centered on Jupiter

      IO_PHI_O (159501), per [4]:

         -- +Z is aligned with the Jupiter pole (primary axis) 
         -- +Y is in the direction of the component of the Io-Jupiter
            vector perpendicular to Z (secondary axis directed at Jupiter)
         -- centered on Io
         
      EUROPA_PHI_O (159502), as above:

         -- +Z is aligned with the Jupiter pole (primary axis)
         -- +Y approximately toward Jupiter (secondary axis toward
            Jupiter)
         -- centered on Europa
         
      GANYMEDE_PHI_O (159503), as above:

         -- +Z is aligned with the Jupiter pole (primary axis)
         -- +Y approximately toward Jupiter (secondary axis toward
            Jupiter)
         -- centered on Ganymede
         
      CALLISTO_PHI_O (159504), as above:

         -- +Z is aligned with the Jupiter pole (primary axis)
         -- +Y approximately toward Jupiter (secondary axis toward
            Jupiter)
         -- centered on Callisto

   \begindata

      FRAME_EUROPAM_JUP_POLE_EUROPA =  -159951
      FRAME_-159951_NAME             = 'EUROPAM_JUP_POLE_EUROPA'
      FRAME_-159951_CLASS            =  5
      FRAME_-159951_CLASS_ID         =  -159951
      FRAME_-159951_CENTER           =  599
      FRAME_-159951_RELATIVE         = 'J2000'
      FRAME_-159951_DEF_STYLE        = 'PARAMETERIZED'
      FRAME_-159951_FAMILY           = 'TWO-VECTOR'
      FRAME_-159951_PRI_AXIS         = 'Z'
      FRAME_-159951_PRI_VECTOR_DEF   = 'CONSTANT'
      FRAME_-159951_PRI_FRAME        = 'IAU_JUPITER'
      FRAME_-159951_PRI_SPEC         = 'RECTANGULAR'
      FRAME_-159951_PRI_VECTOR       = ( 0, 0, 1 )
      FRAME_-159951_SEC_AXIS         = 'X'
      FRAME_-159951_SEC_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
      FRAME_-159951_SEC_OBSERVER     = 'JUPITER'
      FRAME_-159951_SEC_TARGET       = 'EUROPA'
      FRAME_-159951_SEC_ABCORR       = 'NONE'

      FRAME_EUROPAM_JUP_EUROPA_POLE =  -159952
      FRAME_-159952_NAME             = 'EUROPAM_JUP_EUROPA_POLE'
      FRAME_-159952_CLASS            =  5
      FRAME_-159952_CLASS_ID         =  -159952
      FRAME_-159952_CENTER           =  599
      FRAME_-159952_RELATIVE         = 'J2000'
      FRAME_-159952_DEF_STYLE        = 'PARAMETERIZED'
      FRAME_-159952_FAMILY           = 'TWO-VECTOR'
      FRAME_-159952_PRI_AXIS         = 'X'
      FRAME_-159952_PRI_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
      FRAME_-159952_PRI_OBSERVER     = 'JUPITER'
      FRAME_-159952_PRI_TARGET       = 'EUROPA'
      FRAME_-159952_PRI_ABCORR       = 'NONE'
      FRAME_-159952_SEC_AXIS         = 'Z'
      FRAME_-159952_SEC_VECTOR_DEF   = 'CONSTANT'
      FRAME_-159952_SEC_FRAME        = 'IAU_JUPITER'
      FRAME_-159952_SEC_SPEC         = 'RECTANGULAR'
      FRAME_-159952_SEC_VECTOR       = ( 0, 0, 1 )

      FRAME_IO_PHI_O   =  159501
      FRAME_159501_NAME             = 'IO_PHI_O'
      FRAME_159501_CLASS            =  5
      FRAME_159501_CLASS_ID         =  159501
      FRAME_159501_CENTER           =  501
      FRAME_159501_RELATIVE         = 'J2000'
      FRAME_159501_DEF_STYLE        = 'PARAMETERIZED'
      FRAME_159501_FAMILY           = 'TWO-VECTOR'
      FRAME_159501_PRI_AXIS         = 'Z'
      FRAME_159501_PRI_VECTOR_DEF   = 'CONSTANT'
      FRAME_159501_PRI_FRAME        = 'IAU_JUPITER'
      FRAME_159501_PRI_SPEC         = 'RECTANGULAR'
      FRAME_159501_PRI_VECTOR       = ( 0, 0, 1 )
      FRAME_159501_SEC_AXIS         = 'Y'
      FRAME_159501_SEC_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
      FRAME_159501_SEC_OBSERVER     = 'IO'
      FRAME_159501_SEC_TARGET       = 'JUPITER'
      FRAME_159501_SEC_ABCORR       = 'NONE'

      FRAME_EUROPA_PHI_O   =  159502
      FRAME_159502_NAME             = 'EUROPA_PHI_O'
      FRAME_159502_CLASS            =  5
      FRAME_159502_CLASS_ID         =  159502
      FRAME_159502_CENTER           =  502
      FRAME_159502_RELATIVE         = 'J2000'
      FRAME_159502_DEF_STYLE        = 'PARAMETERIZED'
      FRAME_159502_FAMILY           = 'TWO-VECTOR'
      FRAME_159502_PRI_AXIS         = 'Z'
      FRAME_159502_PRI_VECTOR_DEF   = 'CONSTANT'
      FRAME_159502_PRI_FRAME        = 'IAU_JUPITER'
      FRAME_159502_PRI_SPEC         = 'RECTANGULAR'
      FRAME_159502_PRI_VECTOR       = ( 0, 0, 1 )
      FRAME_159502_SEC_AXIS         = 'Y'
      FRAME_159502_SEC_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
      FRAME_159502_SEC_OBSERVER     = 'EUROPA'
      FRAME_159502_SEC_TARGET       = 'JUPITER'
      FRAME_159502_SEC_ABCORR       = 'NONE'

      FRAME_GANYMEDE_PHI_O   =  159503
      FRAME_159503_NAME             = 'GANYMEDE_PHI_O'
      FRAME_159503_CLASS            =  5
      FRAME_159503_CLASS_ID         =  159503
      FRAME_159503_CENTER           =  503
      FRAME_159503_RELATIVE         = 'J2000'
      FRAME_159503_DEF_STYLE        = 'PARAMETERIZED'
      FRAME_159503_FAMILY           = 'TWO-VECTOR'
      FRAME_159503_PRI_AXIS         = 'Z'
      FRAME_159503_PRI_VECTOR_DEF   = 'CONSTANT'
      FRAME_159503_PRI_FRAME        = 'IAU_JUPITER'
      FRAME_159503_PRI_SPEC         = 'RECTANGULAR'
      FRAME_159503_PRI_VECTOR       = ( 0, 0, 1 )
      FRAME_159503_SEC_AXIS         = 'Y'
      FRAME_159503_SEC_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
      FRAME_159503_SEC_OBSERVER     = 'GANYMEDE'
      FRAME_159503_SEC_TARGET       = 'JUPITER'
      FRAME_159503_SEC_ABCORR       = 'NONE'

      FRAME_CALLISTO_PHI_O   =  159504
      FRAME_159504_NAME             = 'CALLISTO_PHI_O'
      FRAME_159504_CLASS            =  5
      FRAME_159504_CLASS_ID         =  159504
      FRAME_159504_CENTER           =  504
      FRAME_159504_RELATIVE         = 'J2000'
      FRAME_159504_DEF_STYLE        = 'PARAMETERIZED'
      FRAME_159504_FAMILY           = 'TWO-VECTOR'
      FRAME_159504_PRI_AXIS         = 'Z'
      FRAME_159504_PRI_VECTOR_DEF   = 'CONSTANT'
      FRAME_159504_PRI_FRAME        = 'IAU_JUPITER'
      FRAME_159504_PRI_SPEC         = 'RECTANGULAR'
      FRAME_159504_PRI_VECTOR       = ( 0, 0, 1 )
      FRAME_159504_SEC_AXIS         = 'Y'
      FRAME_159504_SEC_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
      FRAME_159504_SEC_OBSERVER     = 'CALLISTO'
      FRAME_159504_SEC_TARGET       = 'JUPITER'
      FRAME_159504_SEC_ABCORR       = 'NONE'

   \begintext

End of FK.
