KPL/FK

Juno Frames Kernel
===============================================================================

   This frame kernel contains complete set of frame definitions for the
   Juno (JUNO) spacecraft, its structures and science instruments. This
   frame kernel also contains name - to - NAIF ID mappings for JUNO
   science instruments and s/c structures (see the last section of the
   file.)


Version and Date
-------------------------------------------------------------------------------

   Version 1.2_mod -- February 11, 2022 -- Marshall Styczinski, JPL
   
      Added JUNO_JMAG_O4 frame based on the description in [15]
      and JUNO_JSM frame based on description in [16].
   
   Version 1.2 -- June 24, 2017 -- Boris Semenov, NAIF

      Updated the JUNO_JUNOCAM frame alignment based on [14].

   Version 1.1 -- June 8, 2017 -- Boris Semenov, NAIF

      Updated JUNO_JADE_* frame alignments based on [13].

   Version 1.0 -- June 1, 2017 -- Boris Semenov, NAIF

      Updated JUNO_UVS_BASE alignment angles per [12].

      Added a warning about time-dependency of the JUNO_CHU*_TO_SC,
      JUNO_MOB_*B_TO_SC, and JUNO_SC_TO_MOB_*B alignments to the 
      comments in their sections.

   Version 0.9 -- September 8, 2016 -- Boris Semenov, NAIF

      Added JUNO_MAG_VIP4, JUNO_JSS, JUNO_JSO, JUNO_JH, and
      JUNO_SUN_EQU_RTN frames for magnetospheric studies at Jupiter
      defined in [9].

      Updated JUNO_CHU*_TO_FGM_*B and JUNO_MOB_*B_TO_FGM_*B alignments
      based on [10].

      Updated JUNO_CHU*_TO_SC, JUNO_MOB_*B_TO_SC, and JUNO_SC_TO_MOB_*B
      alignments based on a quick comparison of CHU, MOB and S/C CKs
      for January-March 2016 time period, then again based on
      quick comparison for August 6-19, 2016 period.

      Updated HGA alignment based on [11].

   Version 0.8 -- November 5, 2014 -- Boris Semenov, NAIF

      Redefined SRU frames according to [6]. Note that the updated 
      SRU frames are very different from those defined in the
      earlier version of this kernel. 

   Version 0.7 -- September 16, 2014 -- Boris Semenov, NAIF

      Updated JIRAM alignments based on the information provided by
      Alessandro Mura, INAF-IAPS on September 9 and 16, 2014.

   Version 0.6 -- September 17, 2013 -- Boris Semenov, NAIF

      Corrected the rotation order in the JEDI_090 and JEDI_270 frame
      definitions.

      Added the JUNOCAM_CUBE frame; incorporated the JUNOCAM pre-launch
      DCM and CCD twist into the JUNOCAM_CUBE and JUNOCAM frame
      definitions.

      Incorporated actual JIRAM alignments provided by Alberto Adriani,
      INAF-IAPS on September 16, 2013.
 
   Version 0.5 -- November 15, 2011 -- Boris Semenov, NAIF

      Changed MAG frames:

         -- added layers of "connection" fixed-offset frames
            incorporating CHU? and MOB_?? to FGM_?? alignments
            (CHU?_TO_FGM_??, MOB_??_TO_FGM_??, 6 frames), CHU? and
            MOB_?? to spacecraft alignments (CHU?_TO_SC, MOB_??_TO_SC,
            6 frames), and spacecraft to MOB_?? alignments
            (SC_TO_MOB_??, 2 frames)

   Version 0.4 -- November 4, 2011 -- Boris Semenov, NAIF

      Changed MAG frames and instrument names based on the calibration
      reports and discussions with the MAG team as follows:

         -- renamed ASC_CHU1 -> ASC_CHUD, ASC_CHU2 -> ASC_CHUC,
            ASC_CHU3 -> ASC_CHUB, ASC_CHU4 -> ASC_CHUA

         -- redefined MOB frames to have +Z between the directions of
            CHU +Zs (s/c -Z)

         -- redefined MOB_OB and CHUA/B frames to have +X along s/c -X and
            MOB_IB and CHUC/D frames to have +X along s/c +X

         -- redefined FGM_IB and FGM_OB frames to be CK-based with with
            orientation provided by CKs relative to CHU*, MOB_*B or 
            J2000 frames

   Version 0.3 -- September 1, 2010 -- Boris Semenov, NAIF

      Changed MAG frames and instrument names based on review comments
      by Jack Connerney as follows:

         -- terminology change: in-bound -> inboard, out-bound ->
            outboard

         -- frame and instrument name change: ASC1 -> ASC_CHU1, ASC2 ->
            ASC_CHU2, ASC3 -> ASC_CHU3, ASC4 -> ASC_CHU4

         -- re-defined MOB_OB and MOB_IB frames to be co-aligned with 
            the FGM_OB and FGM_IB frames correspondingly

   Version 0.2 -- July 13, 2009 -- Boris Semenov, NAIF

      Re-defined JUNO_FGM_IB and JUNO_FGM_IB frames based on 
      the input from Jack Connerney.

   Version 0.1 -- June 12, 2009 -- Boris Semenov, NAIF

      Defined preliminary frames for all structures and instruments
      based on CDR materials and draft MICDs.

   Version 0.0 -- February 3, 2009 -- Boris Semenov, NAIF

      Initial Release: bare-bones version with just two frames needed
      to access JUNO nominal attitude CK file.


References
-------------------------------------------------------------------------------

   1. ``Frames Required Reading''

   2. ``Kernel Pool Required Reading''

   3. ``C-Kernel Required Reading''

   4. JUNO Instrument ICDs and MICDs

   5. JUNO CDR Materials.

   6. JUNO-ME-11-5159_Pointing_Alignment Verification.xlsx

   7. Updated JIRAM alignment angles provided in the updated FK
      ``juno_v02.corretto.new.2.tf'' by Alberto Adriani, INAF-IAPS on
      September 16, 2013.

   8. E-mail ``Modifications to JUNO/JIRAM *tf and *ti kernel files.''
      from Alessandro Mura, INAF-IAPS, September 9 and 16, 2014.

   9. ``Jupiter Coordinate Systems'', F. Bagenal & R. J. Wilson, LASP -
      U of Colorado, TBD -- 07/15/16

   10. JN-DTU-TN-3105_v1_1_Juno_MAG_to_mASC_relative_orientation.pdf

   11. E-mail stating the measured HGA electrical boresight direction,
       from Kristen Francis, LMCO to William Folkner, JPL, 09/07/16

   12. juno_v08rev_tkg.tf; modified FK used by the UVS team, SWRI; 
       provided by Brad Trantham on 05/31/17

   13. ANODE_LOOK_ELC_DEFL_NONE_V02.CSV and ANODE_LOOK_ION_DEFL_NONE_V02.CSV
       from JNO-J/SW-JAD-3-CALIBRATED-V1.0, provided by Rob Wilson, LASP
       on 05/25/17.

   14. E-mail from Mike Caplinger, MSSS re updated Junocam camera model, 
       06/20/17
       
   15. https://lasp.colorado.edu/home/mop/files/2015/02/CoOrd_systems12.pdf
   
   16. https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=GO-J-POS-6-SC-TRAJ-JUP-COORDS-V1.0


Contact Information
-------------------------------------------------------------------------------

   Boris V. Semenov, NAIF/JPL, (818)-354-8136, bsemenov@spice.jpl.nasa.gov


Implementation Notes
-------------------------------------------------------------------------------

   This file is used by the SPICE system as follows: programs that make
   use of this frame kernel must ``load'' the kernel, normally during
   program initialization. The SPICE routine FURNSH loads a kernel file
   into the pool as shown below.

      CALL FURNSH ( 'frame_kernel_name; )    -- FORTRAN
      furnsh_c ( "frame_kernel_name" );      -- C
      cspice_furnsh, frame_kernel_name       -- IDL
      cspice_furnsh( 'frame_kernel_name' )   -- MATLAB

   This file was created and may be updated with a text editor or word
   processor.


JUNO Frames
-------------------------------------------------------------------------------

   The following JUNO frames are defined in this kernel file:

           Name                  Relative to           Type       NAIF ID
      ======================  =====================  ============   =======

   Frames for Magnetospheric Studies at Jupiter:
   ---------------------------------------------
      JUNO_MAG_VIP4           IAU_JUPITER            FIXED          -61952
      JUNO_JSS                J2000                  DYNAMIC        -61953
      JUNO_JSO                J2000                  DYNAMIC        -61954
      JUNO_JH                 J2000                  DYNAMIC        -61955
      JUNO_SUN_EQU_RTN        J2000                  DYNAMIC        -61956
      JUNO_JMAG_O4            IAU_JUPITER            FIXED         1661957
      JUNO_JSM                J2000                  DYNAMIC       1661958

   Spacecraft frame:
   -----------------
      JUNO_SPACECRAFT         varies                 CK             -61000
      JUNO_SPIN_AXIS          J2000                  CK             -61900

   Magnetometer and ASC CHU frames:
   --------------------------------
      JUNO_ASC_CHUD           J2000, SPACECRAFT      CK             -61111
      JUNO_ASC_CHUC           J2000, SPACECRAFT      CK             -61112
      JUNO_MOB_IB             J2000                  CK             -61113
      JUNO_FGM_IB             MOB_IB, CHUD, CHUC     CK             -61114

      JUNO_ASC_CHUB           J2000, SPACECRAFT      CK             -61121
      JUNO_ASC_CHUA           J2000, SPACECRAFT      CK             -61122
      JUNO_MOB_OB             J2000                  CK             -61123
      JUNO_FGM_OB             MOB_OB, CHUB, CHUA     CK             -61124

      JUNO_CHUA_TO_FGM_OB     CHUA                   FIXED          -61131 (a)
      JUNO_CHUB_TO_FGM_OB     CHUB                   FIXED          -61132 (b)
      JUNO_CHUC_TO_FGM_IB     CHUC                   FIXED          -61133 (c)
      JUNO_CHUD_TO_FGM_IB     CHUD                   FIXED          -61134 (d)
      JUNO_MOB_OB_TO_FGM_OB   MOB_OB                 FIXED          -61135 (e)
      JUNO_MOB_IB_TO_FGM_IB   MOB_IB                 FIXED          -61136 (f)

      JUNO_CHUA_TO_SC         CHUA                   FIXED          -61141 (g)
      JUNO_CHUB_TO_SC         CHUB                   FIXED          -61142 (h)
      JUNO_CHUC_TO_SC         CHUC                   FIXED          -61143 (i)
      JUNO_CHUD_TO_SC         CHUD                   FIXED          -61144 (j)
      JUNO_MOB_OB_TO_SC       MOB_OB                 FIXED          -61145 (k)
      JUNO_MOB_IB_TO_SC       MOB_IB                 FIXED          -61146 (l)

      JUNO_SC_TO_MOB_OB       SPACECRAFT             FIXED          -61155 (m)
      JUNO_SC_TO_MOB_IB       SPACECRAFT             FIXED          -61156 (n)

   JADE Frames:
   ------------
      JUNO_JADE_E060          SPACECRAFT             FIXED          -61201
      JUNO_JADE_E180          SPACECRAFT             FIXED          -61202
      JUNO_JADE_E300          SPACECRAFT             FIXED          -61203
      JUNO_JADE_I             SPACECRAFT             FIXED          -61204

   JEDI Frames:
   ------------
      JUNO_JEDI_090           SPACECRAFT             FIXED          -61301
      JUNO_JEDI_A180          SPACECRAFT             FIXED          -61302
      JUNO_JEDI_270           SPACECRAFT             FIXED          -61303

   JIRAM Frames:
   -------------
      JUNO_JIRAM_URF          SPACECRAFT             FIXED          -61401
      JUNO_JIRAM_I            JIRAM_URF              FIXED          -61410
      JUNO_JIRAM_I_LBAND      JIRAM_IMAGER           FIXED          -61411
      JUNO_JIRAM_I_MBAND      JIRAM_IMAGER           FIXED          -61412
      JUNO_JIRAM_S            JIRAM_IMAGER           FIXED          -61420

   JUNOCAM Frames:
   ---------------
      JUNO_JUNOCAM_CUBE       SPACECRAFT             FIXED          -61505
      JUNO_JUNOCAM            JUNOCAM_CUBE           FIXED          -61500

   MWR Frames:
   -----------
      JUNO_MWR_A1             SPACECRAFT             FIXED          -61601
      JUNO_MWR_A2             SPACECRAFT             FIXED          -61602
      JUNO_MWR_A3             SPACECRAFT             FIXED          -61603
      JUNO_MWR_A4             SPACECRAFT             FIXED          -61604
      JUNO_MWR_A5             SPACECRAFT             FIXED          -61605
      JUNO_MWR_A6             SPACECRAFT             FIXED          -61606

   UVS Frames:
   -----------
      JUNO_UVS_BASE           SPACECRAFT             FIXED          -61700
      JUNO_UVS                UVS_BASE               CK             -61701

   WAVES Frames:
   -----------
      JUNO_WAVES_MSC          SPACECRAFT             FIXED          -61810
      JUNO_WAVES_ANTENNA      SPACECRAFT             FIXED          -61820

   Solar Array frames:
   -------------------
      JUNO_SA1_HINGE          SPACECRAFT             FIXED          -61010
      JUNO_SA1                SA1_HINGE              CK             -61011
      JUNO_SA2_HINGE          SPACECRAFT             FIXED          -61020
      JUNO_SA2                SA2_HINGE              CK             -61021
      JUNO_SA3_HINGE          SPACECRAFT             FIXED          -61030
      JUNO_SA3                SA3_HINGE              CK             -61031

   Antenna frames:
   ---------------
      JUNO_HGA                SPACECRAFT             FIXED          -61040
      JUNO_MGA                SPACECRAFT             FIXED          -61050
      JUNO_LGA_FORWARD        SPACECRAFT             FIXED          -61061
      JUNO_LGA_AFT            SPACECRAFT             FIXED          -61062
      JUNO_LGA_TOROID         SPACECRAFT             FIXED          -61063

   ACS Sensor frames:
   ------------------
      JUNO_SRU1               SPACECRAFT             FIXED          -61071
      JUNO_SRU2               SPACECRAFT             FIXED          -61072
      JUNO_SSS1               SPACECRAFT             FIXED          -61073
      JUNO_SSS2               SPACECRAFT             FIXED          -61074

   Truster frames:
   ---------------
      JUNO_REM_FL1            SPACECRAFT             FIXED          -61081
      JUNO_REM_FL2            SPACECRAFT             FIXED          -61082
      JUNO_REM_FL3            SPACECRAFT             FIXED          -61083
      JUNO_REM_FL4            SPACECRAFT             FIXED          -61084
      JUNO_REM_FA1            SPACECRAFT             FIXED          -61085
      JUNO_REM_FA2            SPACECRAFT             FIXED          -61086

      JUNO_REM_RL1            SPACECRAFT             FIXED          -61091
      JUNO_REM_RL2            SPACECRAFT             FIXED          -61092
      JUNO_REM_RL3            SPACECRAFT             FIXED          -61093
      JUNO_REM_RL4            SPACECRAFT             FIXED          -61094
      JUNO_REM_RA1            SPACECRAFT             FIXED          -61095
      JUNO_REM_RA2            SPACECRAFT             FIXED          -61096


JUNO Frames Hierarchy
-------------------------------------------------------------------------------

   The diagram below shows the JUNO frames hierarchy:

    JUNO_MAG_VIP4      JUNO_JSO                             JUNO_SUN_EQU_RTN 
    --------------     --------                             ----------------
        ^                 ^                                              ^
        |<-fxd            |<-dyn                                    dyn->|
        |                 |                                              |
        |     JUNO_JSS    |                                  JUNO_JH     |
        |     --------    |                                  -------     |
        |        ^        |                                     ^        |
        |        |<-dyn   |                                     |<-dyn   |
        |        |        |                                     |        |
        |        |        |      "J2000" INERTIAL               |        |
        +----------------------------------------------------------------+
        |        | | |                 |  |                 | | |        |
        |<-pck   | | |             ck->|  |<-ck             | | |   pck->|
        |        | | |                 |  |                 | | |        |
        V        | | |                 V  |                 | | |        V
   "IAU_JUPITER" | | |   "JUNO_SPIN_AXIS" |                 | | |  "IAU_EARTH"
   ------------- | | |   ---------------- |                 | | |  -----------
                 | | |                 |  |                 | | |
             ck->| | |<-ck         ck->|  |             ck->| | |<-ck 
                 | | |                 |  |                 | | |
                 V | V                 |  |                 V | V
   "JUNO_ASC_CHUD" | "JUNO_ASC_CHUC"   |  |   "JUNO_ASC_CHUB" | "JUNO_ASC_CHUA"
   --------------- | ---------------   |  |   --------------- | ---------------
        |  |       |       |  |        |  |        |  |       |       |  |
        |  |       |<-ck   |  |        |  |        |  |       |<-ck   |  |
        |  |       |       |  |        |  |        |  |       |       |  |
        |  |       |       |  |        |  |        |  |       |       |  |
        |  |       |       |  |   fxd  |  | fxd    |  |       |       |  |
        |  | "(n)"<-----------------.  |  |  .----------------->"(m)" |  |
        |  | ----- |       |  |     |  |  |  |     |  |       | ----- |  |
        |  |   |   |       |  | ck  |  |  |  | ck  |  |       |   |   |  |
        |  |   |   | "(l)"-------.  |  |  |  |  .-------"(k)" |   |   |  |
        |  |   |   | ----- |  |  |  |  |  |  |  |  |  | ----- |   |   |  |
        |  |   |   |   ^   |  |  |  |  |  |  |  |  |  |   ^   |   |   |  |
        |  |ck>|   |fx>|   |  |  |  |  |  |  |  |  |  |   |<fx|   |<ck|  |
        |  |   |   |   |   |  |  |  |  |  |  |  |  |  |   |   |   |   |  |
        |  |   V   V       |  |  |  |  |  |  |  |  |  |       V   V   |  |
        |  | "JUNO_MOB_IB" |  |  |  |  |  |  |  |  |  | "JUNO_MOB_OB" |  |
        |  | ------------- |  |  |  |  |  |  |  |  |  | ------------- |  |
        |  |       |       |  |  |  |  |  |  |  |  |  |       |       |  |
        |  |<-fxd  |  fxd->|  |  |  |  |  |  |  |  |  |<-fxd  |  fxd->|  |
        |  |       |       |  |  |  |  |  |  |  |  |  |       |       |  |
        |  |  fxd->|       |  |  |  |  |  |  |  |  |  |  fxd->|       |  |
        |  |       |       |  |  |  |  |  |  |  |  |  |       |       |  |
        |  V       V       V  |  |  |  |  |  |  |  |  V       V       v  |
        | "(d)"  "(f)"  "(c)" |  |  |  |  |  |  |  | "(b)"  "(e)"  "(a)" |
        | -----  -----  ----- |  |  |  |  |  |  |  | -----  -----  ----- |
        |  |       |       |  |  |  |  |  |  |  |  |  |       |       |  |
        |  |<-ck   |   ck->|  |  |  |  |  |  |  |  |  |<-ck   |   ck->|  |
        |  |       |       |  |  |  |  |  |  |  |  |  |       |       |  |
        |  |       |<-ck   |  |  |  |  |  |  |  |  |  |       |<-ck   |  |
        |  |       |       |  |  |  |  |  |  |  |  |  |       |       |  |
        |  V       V       V  |  |  |  |  |  |  |  |  V       V       V  |
        |    "JUNO_FGM_IB"    |  |  |  |  |  |  |  |    "JUNO_FGM_OB"    |
        |    -------------    |  |  |  |  |  |  |  |    -------------    |
        |                     |  |  |  |  |  |  |  |                     |
        |                     |  |  |  |  |  |  |  |                     |
        |<-fxd           fxd->|  |  |  |  |  |  |  |<-fxd           fxd->|
        |                     |  |  |  |  |  |  |  |                     |
        V                     V  |  |  |  |  |  |  V                     V
      "(j)"                "(i)" |  |  |  |  |  | "(h)"                "(g)"
      -----                ----- |  |  |  |  |  | -----                -----   
        |                     |  |  |  |  |  |  |  |                     |
        |<-ck             ck->|  |  |  |  |  |  |  |<-ck                 |<-ck
        |                     |  |  |  |  |  |  |  |                     |
        V                     V  V  |  V  V  |  V  V                     V
                                "JUNO_SPACECRAFT"
        +-----------------------------------------------------------------+
            |    |    |    |    |    |    |            |  |  |  |  |  |  |
       fxd->|    |    |    |    |    |    |       fxd->|  |  |  |  |  |  |
            |    |    |    |    |    |    |            |  |  |  |  |  |  |
            V    |    |    |    |    |    |            V  |  |  |  |  |  |
   "JUNO_JADE_*" |    |    |    |    |    |    "JUNO_HGA" |  |  |  |  |  |
   ------------- |    |    |    |    |    |    ---------- |  |  |  |  |  |
                 |    |    |    |    |    |               |  |  |  |  |  |
            fxd->|    |    |    |    |    |          fxd->|  |  |  |  |  |
                 |    |    |    |    |    |               |  |  |  |  |  |
                 V    |    |    |    |    |               V  |  |  |  |  |
      "JUNO_JEDI_*"   |    |    |    |    |       "JUNO_MGA" |  |  |  |  |
      -------------   |    |    |    |    |       ---------- |  |  |  |  |
                      |    |    |    |    |                  |  |  |  |  |
                 fxd->|    |    |    |    |             fxd->|  |  |  |  |
                      |    |    |    |    |                  |  |  |  |  |
                      V    |    |    |    |                  V  |  |  |  |
        "JUNO_JIRAM_URF"   |    |    |    |        "JUNO_LGA_*" |  |  |  |
        ----------------   |    |    |    |        ------------ |  |  |  |
                 |         |    |    |    |                     |  |  |  |
            fxd->|         |    |    |    |                fxd->|  |  |  |
                 |         |    |    |    |                     |  |  |  |
                 V         |    |    |    |                     V  |  |  |
          "JUNO_JIRAM_S"   |    |    |    |           "JUNO_REM_*" |  |  |
          --------------   |    |    |    |           ------------ |  |  |
                 |         |    |    |    |                        |  |  |
            fxd->|         |    |    |    |<-fxd              fxd->|  |  |
                 |         |    |    |    |                        |  |  |
                 V         |    |    |    V                        V  |  |
          "JUNO_JIRAM_I"   |    |    |   "JUNO_WAVES_*"   "JUNO_SRU*" |  |
          --------------   |    |    |   --------------   ----------- |  |
                 |         |    |    |                                |  |
            fxd->|         |    |    |<-fxd                      fxd->|  |
                 |         |    |    |                                |  |
                 V         |    |    V                                V  |
    "JUNO_JIRAM_I_*BAND"   |    |   "JUNO_UVS_BASE"          "JUNO_SSS*" |
    --------------------   |    |   ---------------          ----------- |
                           |    |    |                                   |
                      fxd->|    |    |<-ck                          fxd->|
                           |    |    |                                   |
                           V    |    V                                   V
          "JUNO_JUNOCAM_CUBE"   |   "JUNO_UVS"             "JUNO_SA*_HINGE"
          -------------------   |   ----------             ----------------
                           |    |                                        |
                      fxd->|    |                                        |
                           |    |                                        |
                           V    |                                        |
               "JUNO_JUNOCAM"   |                                        |
               --------------   |                                        |
                                |                                        |
                           fxd->|                                    ck->|
                                |                                        |
                                V                                        V
                     "JUNO_MWR_A*"                                "JUNO_SA*"
                     -------------                                ----------


Spacecraft Orientation Frame Chains

   Possible "J2000" -> "JUNO_SPACECRAFT" frame chains are:

      -- using data from a nominal spacecraft CK:

            "J2000" 
            -ck->  "JUNO_SPIN_AXIS"
            -ck->  "JUNO_SPACECRAFT"

      -- using data from a reconstructed spacecraft CK:

            "J2000"
            -ck->  "JUNO_SPACECRAFT"

      -- using data from a reconstructed MOB_## CK, the MOB_##-to-S/C
         alignment fixed-offset frame, and an alignment frame-to-s/c
         "coverage" CK:

            "J2000"  
            -ck->  "JUNO_MOB_##"
            -fxd-> "JUNO_MOB_##_TO_SC"  
            -ck->  "JUNO_SPACECRAFT"

      -- using data from a reconstructed CHU# CK, the CHU#-to-S/C
         alignment fixed-offset frame, and an alignment frame-to-s/c
         "coverage" CK:

            "J2000"
            -ck->  "JUNO_ASC_CHU#"
            -fxd-> "JUNO_CHU#_TO_SC"
            -ck->  "JUNO_SPACECRAFT"


Magnetometer Sensor Frame Chains

   Possible "J2000" -> "JUNO_FGM_IB" frame chains are (same for OB):

      -- using data from a nominal spacecraft CK, the s/c-to-MOB_IB
         alignment fixed-offset frame, an alignment frame-to-MOB_IB
         "coverage" CK, the MOB_IB-to-FGM_IB alignment fixed-offset
         frame, and an alignment-to-FGM_IB "coverage" CK:

            "J2000" 
            -ck->  "JUNO_SPIN_AXIS"
            -ck->  "JUNO_SPACECRAFT"
            -fxd-> "JUNO_SC_TO_MOB_IB"
            -ck->  "JUNO_MOB_IB"
            -fxd-> "JUNO_MOB_IB_TO_FGM_IB"
            -ck->  "JUNO_FGM_IB" 

      -- using data from a reconstructed spacecraft CK, the s/c-to-MOB_IB
         alignment fixed-offset frame, an alignment frame-to-MOB_IB
         "coverage" CK, the MOB_IB-to-FGM_IB alignment fixed-offset
         frame, and an alignment-to-FGM_IB "coverage" CK:

            "J2000" 
            -ck->  "JUNO_SPACECRAFT"
            -fxd-> "JUNO_SC_TO_MOB_IB"
            -ck->  "JUNO_MOB_IB"
            -fxd-> "JUNO_MOB_IB_TO_FGM_IB"
            -ck->  "JUNO_FGM_IB" 

      -- using data from a reconstructed MOB_IB CK, the MOB_IB-to-FGM_IB
         alignment fixed-offset frame, and an alignment-to-FGM_IB
         "coverage" CK:

            "J2000" 
            -ck->  "JUNO_MOB_IB"
            -fxd-> "JUNO_MOB_IB_TO_FGM_IB"
            -ck->  "JUNO_FGM_IB" 

      -- using data from a reconstructed CHU# CK, the CHU#-to-FGM_IB
         alignment fixed-offset frame, and an alignment-to-FGM_IB
         "coverage" CK:

            "J2000" 
            -ck->  "JUNO_ASC_CHU#"
            -fxd-> "JUNO_CHU#_TO_FGM_IB"
            -ck->  "JUNO_FGM_IB" 


Frames for Magnetospheric Studies at Jupiter
-------------------------------------------------------------------------------

   This section defines a few frames for magnetospheric studies at 
   Jupiter described in [9].


Jupiter Magnetic VIP4 Frame

   The JUNO_MAG_VIP4 frame implements the JUNO mission Jupiter
   Magnetic VIP4 (MAG_VIP4) reference frame described in [9].
 
   It is defined as a fixed offset frame with respect to the
   IAU_JUPITER frame (which is equivalent to the System III right
   handed frame when used with pck00010.tpc) as follows:

      -  +Z axis is along planetocentric LON=159.2 deg/LAT=80.5 deg in
         in the IAU_JUPITER frame
        
      -  +Y axis is along planetocentric LON=249.2 deg/LAT=0 deg in
         the IAU_JUPITER frame
        
      -  +X completes the right handed frame

      -  the center is at the center of Jupiter.

   Two rotations are needed to align the IAU_JUPITER frame with the
   JUNO_MAG_VIP4 frame: first by +159.2 degrees about Z, second by +9.5
   degrees about Y.
 
   The keywords below implement the JUNO_MAG_VIP4 frame. Since the
   frame definition below contains the reverse transformation -- from
   the JUNO_MAG_VIP4 frame to the IAU_JUPITER frame -- the order of
   rotations is reversed and the signs of rotation angles are changed
   to the opposite ones compared to the description above.

   \begindata

      FRAME_JUNO_MAG_VIP4          = -61952
      FRAME_-61952_NAME            = 'JUNO_MAG_VIP4'
      FRAME_-61952_CLASS           = 4
      FRAME_-61952_CLASS_ID        = -61952
      FRAME_-61952_CENTER          = 599
      TKFRAME_-61952_SPEC          = 'ANGLES'
      TKFRAME_-61952_RELATIVE      = 'IAU_JUPITER'
      TKFRAME_-61952_ANGLES        = ( 0, -159.2, -9.5 )
      TKFRAME_-61952_AXES          = ( 2,    3,    2   )
      TKFRAME_-61952_UNITS         = 'DEGREES'

   \begintext


Jupiter-De-Spun-Sun Frame

   The JUNO_JSS frame implements the JUNO mission Jupiter-De-Spun-Sun
   (JSS) reference frame described in [9].

   It is defined as a dynamic frame as follows:

      -  +Z axis is along the Jupiter north pole (+Z of the IAU_JUPITER
         frame)
 
      -  +X axis is in the direction of the geometric position of the
         Sun as seen from Jupiter
 
      -  +Y completes the right handed frame
 
      -  the center is at the center of Jupiter.

   The keywords below implement the JUNO_JSS frame as a dynamic frame.

   \begindata

      FRAME_JUNO_JSS               = -61953
      FRAME_-61953_NAME            = 'JUNO_JSS'
      FRAME_-61953_CLASS           = 5
      FRAME_-61953_CLASS_ID        = -61953
      FRAME_-61953_CENTER          = 599
      FRAME_-61953_RELATIVE        = 'J2000'
      FRAME_-61953_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_-61953_FAMILY          = 'TWO-VECTOR'
      FRAME_-61953_PRI_AXIS        = 'Z'
      FRAME_-61953_PRI_VECTOR_DEF  = 'CONSTANT'
      FRAME_-61953_PRI_FRAME       = 'IAU_JUPITER'
      FRAME_-61953_PRI_SPEC        = 'RECTANGULAR'
      FRAME_-61953_PRI_VECTOR      =  ( 0, 0, 1 )
      FRAME_-61953_SEC_AXIS        = 'X'
      FRAME_-61953_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_-61953_SEC_OBSERVER    = 'JUPITER'
      FRAME_-61953_SEC_TARGET      = 'SUN'
      FRAME_-61953_SEC_ABCORR      = 'NONE'

   \begintext


Jupiter-Sun-Orbit Frame

   The JUNO_JSO frame implements the JUNO mission Jupiter-Sun-Orbit
   (JSO) reference frame described in [9].

   It is defined as a dynamic frame as follows:

      -  +X axis is along the geometric position of the Sun as seen
         from Jupiter
 
      -  +Y axis is in the direction of the inertial geometric velocity 
         of the Sun as seen from Jupiter
 
      -  +Z completes the right handed frame
 
      -  the center is at the center of Jupiter.

   The keywords below implement the JUNO_JSO frame as a dynamic frame.

   \begindata

      FRAME_JUNO_JSO               = -61954
      FRAME_-61954_NAME            = 'JUNO_JSO'
      FRAME_-61954_CLASS           = 5
      FRAME_-61954_CLASS_ID        = -61954
      FRAME_-61954_CENTER          = 599
      FRAME_-61954_RELATIVE        = 'J2000'
      FRAME_-61954_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_-61954_FAMILY          = 'TWO-VECTOR'
      FRAME_-61954_PRI_AXIS        = 'X'
      FRAME_-61954_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_-61954_PRI_OBSERVER    = 'JUPITER'
      FRAME_-61954_PRI_TARGET      = 'SUN'
      FRAME_-61954_PRI_ABCORR      = 'NONE'
      FRAME_-61954_SEC_AXIS        = 'Y'
      FRAME_-61954_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_VELOCITY'
      FRAME_-61954_SEC_OBSERVER    = 'JUPITER'
      FRAME_-61954_SEC_TARGET      = 'SUN'
      FRAME_-61954_SEC_ABCORR      = 'NONE'
      FRAME_-61954_SEC_FRAME       = 'J2000'

   \begintext


Jupiter Heliospheric Frame

   The JUNO_JH frame implements the JUNO mission Jupiter Heliospheric
   (JH) reference frame described in [9].

   It is defined as a dynamic frame as follows:

      -  +X axis is along the geometric position of the Sun as seen
         from Jupiter
 
      -  +Z axis is in the direction of the Sun north pole.
 
      -  +Y completes the right handed frame
 
      -  the center is at the center of Jupiter.

   The keywords below implement the JUNO_JH frame as a dynamic frame.

   \begindata

      FRAME_JUNO_JH               = -61955
      FRAME_-61955_NAME            = 'JUNO_JH'
      FRAME_-61955_CLASS           = 5
      FRAME_-61955_CLASS_ID        = -61955
      FRAME_-61955_CENTER          = 599
      FRAME_-61955_RELATIVE        = 'J2000'
      FRAME_-61955_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_-61955_FAMILY          = 'TWO-VECTOR'
      FRAME_-61955_PRI_AXIS        = 'X'
      FRAME_-61955_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_-61955_PRI_OBSERVER    = 'JUPITER'
      FRAME_-61955_PRI_TARGET      = 'SUN'
      FRAME_-61955_PRI_ABCORR      = 'NONE'
      FRAME_-61955_SEC_AXIS        = 'Z'
      FRAME_-61955_SEC_VECTOR_DEF  = 'CONSTANT'
      FRAME_-61955_SEC_FRAME       = 'IAU_SUN'
      FRAME_-61955_SEC_SPEC        = 'RECTANGULAR'
      FRAME_-61955_SEC_VECTOR      =  ( 0, 0, 1 )

   \begintext


Juno Solar Equatorial RTN Frame

   The JUNO_SUN_EQU_RTN frame implements the JUNO mission JUNO RTN
   reference frame described in [9].

   It is defined as a dynamic frame as follows:

      -  +X axis is along the geometric position of Juno as seen
         from the Sun
 
      -  +Z axis is in the direction of the Sun north pole.
 
      -  +Y completes the right handed frame
 
      -  the center is at the JUNO spacecraft.

   The keywords below implement the JUNO_SUN_EQU_RTN frame as a dynamic
   frame.

   \begindata

      FRAME_JUNO_SUN_EQU_RTN       = -61956
      FRAME_-61956_NAME            = 'JUNO_SUN_EQU_RTN'
      FRAME_-61956_CLASS           = 5
      FRAME_-61956_CLASS_ID        = -61956
      FRAME_-61956_CENTER          = -61
      FRAME_-61956_RELATIVE        = 'J2000'
      FRAME_-61956_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_-61956_FAMILY          = 'TWO-VECTOR'
      FRAME_-61956_PRI_AXIS        = 'X'
      FRAME_-61956_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_-61956_PRI_OBSERVER    = 'SUN'
      FRAME_-61956_PRI_TARGET      = 'JUNO'
      FRAME_-61956_PRI_ABCORR      = 'NONE'
      FRAME_-61956_SEC_AXIS        = 'Z'
      FRAME_-61956_SEC_VECTOR_DEF  = 'CONSTANT'
      FRAME_-61956_SEC_FRAME       = 'IAU_SUN'
      FRAME_-61956_SEC_SPEC        = 'RECTANGULAR'
      FRAME_-61956_SEC_VECTOR      =  ( 0, 0, 1 )

   \begintext
   

Jupiter O4 Magnetic Frame

   The JUNO_JMAG_O4 frame is identical to JUNO_MAG_VIP4, but with
   the dipole vector defined using the O4 model instead of the VIP4
   model, as described in [15].
 
   It is defined as a fixed offset frame with respect to the
   IAU_JUPITER frame (which is equivalent to the System III right
   handed frame when used with pck00010.tpc) as follows:

      -  +Z axis is along planetocentric LON=158.0 deg/LAT=80.4 deg in
         in the IAU_JUPITER frame
        
      -  +Y axis is along planetocentric LON=248.0 deg/LAT=0 deg in
         the IAU_JUPITER frame
        
      -  +X completes the right handed frame

      -  the center is at the center of Jupiter.

   Two rotations are needed to align the IAU_JUPITER frame with the
   JUNO_JMAG_O4 frame: first by +158.0 degrees about Z, second by +9.6
   degrees about Y.
 
   The keywords below implement the JUNO_JMAG_O4 frame. Since the
   frame definition below contains the reverse transformation -- from
   the JUNO_JMAG_O4 frame to the IAU_JUPITER frame -- the order of
   rotations is reversed and the signs of rotation angles are changed
   to the opposite ones compared to the description above.

   \begindata

      FRAME_JUNO_JMAG_O4            = 1661957
      FRAME_1661957_NAME            = 'JUNO_JMAG_O4'
      FRAME_1661957_CLASS           = 4
      FRAME_1661957_CLASS_ID        = 1661957
      FRAME_1661957_CENTER          = 599
      TKFRAME_1661957_SPEC          = 'ANGLES'
      TKFRAME_1661957_RELATIVE      = 'IAU_JUPITER'
      TKFRAME_1661957_ANGLES        = ( 0, -158.0, -9.6 )
      TKFRAME_1661957_AXES          = ( 2,    3,    2   )
      TKFRAME_1661957_UNITS         = 'DEGREES'

   \begintext


Jupiter Solar Magnetospheric frame

   The JUNO_JSM frame implements the Jupiter Solar Magnetospheric
   frame described in [16].

   It is defined as a dynamic frame as follows:

      -  +X axis is along the geometric position of the Sun as seen
         from Jupiter
  
      -  +Y axis is along M x X, where M is the direction of the
         magnetic dipole moment as defined in the O4 model, along
         LON=158.0 deg/LAT=80.4 deg in the IAU_JUPITER frame
      
      -  +Z axis completes the right-handed set, along X x Y.
 
      -  the center is at the center of Jupiter.

   The keywords below implement the JUNO_JSM frame as a dynamic frame.
   The +Z axis is specified as a secondary vector, so the actual
   direction will be the component of the dipole axis perpendicular
   to +X according to [1], with +Y completing the right-handed set.

   \begindata

      FRAME_JUNO_JSM                = 1661958
      FRAME_1661958_NAME            = 'JUNO_JSM'
      FRAME_1661958_CLASS           = 5
      FRAME_1661958_CLASS_ID        = 1661958
      FRAME_1661958_CENTER          = 599
      FRAME_1661958_RELATIVE        = 'J2000'
      FRAME_1661958_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_1661958_FAMILY          = 'TWO-VECTOR'
      FRAME_1661958_PRI_AXIS        = 'X'
      FRAME_1661958_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_1661958_PRI_OBSERVER    = 'JUPITER'
      FRAME_1661958_PRI_TARGET      = 'SUN'
      FRAME_1661958_PRI_ABCORR      = 'NONE'
      FRAME_1661958_SEC_AXIS        = 'Z'
      FRAME_1661958_SEC_VECTOR_DEF  = 'CONSTANT'
      FRAME_1661958_SEC_FRAME       = 'IAU_JUPITER'
      FRAME_1661958_SEC_SPEC        = 'LATITUDINAL'
      FRAME_1661958_SEC_UNITS       = 'DEGREES'
      FRAME_1661958_SEC_LONGITUDE   = 158.0
      FRAME_1661958_SEC_LATITUDE    = 80.4

   \begintext


Spacecraft Bus Frame
-------------------------------------------------------------------------------

   The spacecraft frame (or AACS control frame) is defined by the s/c design
   as follows [from 5]:

      -  +Z axis is along the nominal spin axis and points in the
         direction of the nominal HGA boresight

      -  +X axis is along the solar array 1 symmetry axis and points
         towards the magnetometer boom

      -  +Y axis completes the right-handed frame

      -  the origin of the frame is centered on the launch vehicle
         separation plane.

   These diagrams illustrates the s/c frame:

      Spacecraft +Z side view:
      ------------------------

                .-\
             .-'   \
          .-'       \ Solar Array 3
       .-'           \
       \           .-'\
        \       .-'    \
         \   ,-'        \
          \-'         .-'\
           \       .-'    \
            \   .-'        \
             \-'         .-'\             <-. Spin Direction
              \       .-'    \               `.
               \   .-'        \.
                \-'         .-' `-.
                 \       .-'       `-.
                  \   .-'    ----- .  `-.      Solar Array 1
                   .-'    +Ysc ^    `.   `-.--------------------.
                   |    /      |      \    |      |      |      |
                   |   /       |       \   |      |      |      |``--..
                   |  '        |    +Xsc   |      |      |      |      `--..
                   |  |        o------> |  |      |      |      |           |
                   |  .      +Zsc       ,  |      |      |      |     ..--''
                   |   \               /   |      |      |      |..--'
                   |    \ HGA         /    |      |      |      | Magnetometer
                   `-.   `.         .'   .-'--------------------'     Boom
                  /   `-.  ` ----- '  .-'
                 /       `-.       .-'
                /           `-. .-'
               /`-.           /'
              /    `-.       /
             /        `-.   /
            /`-.         `-/
           /    `-.       /
          /        `-.   /
         /`-.         `-/
        /    `-.       /
       /        `-.   /
      /.           `-/
        `-.         / Solar Array 2                  +Zsc is out of
           `-.     /                                    the page.
              `-. /
                 `


      Spacecraft -Y side view:
      ------------------------
                               o
                              / \
                      .-----------------, HGA
                       \               /
                        `.           .'
                           `-------'
                           |       |                            Magnetometer
                           |       |                                Boom
      ================o==o==o==o-----------o======o======o======o===========*
          Solar    |           |           |   Solar Array 1
         Array 2   |           |           |
                   |                       |
                   |                       |
                   |      +Zsc ^           |
                   |           |           |
                   |           |           |
                   |           |           |            +Ysc is into
                   `----- +Ysc x------>   -'              the page.
                            /_____\  +Xsc


   Since the S/C bus attitude with respect to an inertial frame is provided
   by a C kernel (see [3] for more information), this frame is defined as
   a CK-based frame.

   \begindata

      FRAME_JUNO_SPACECRAFT        = -61000
      FRAME_-61000_NAME            = 'JUNO_SPACECRAFT'
      FRAME_-61000_CLASS           = 3
      FRAME_-61000_CLASS_ID        = -61000
      FRAME_-61000_CENTER          = -61
      CK_-61000_SCLK               = -61
      CK_-61000_SPK                = -61

   \begintext


Spin Axis Frame
-------------------------------------------------------------------------------


   The JUNO_SPIN_AXIS frame is a special frame used in the nominal
   orientation CK files. In these files the JUNO_SPACECRAFT frame
   orientation is not stored relative to the J2000 frame. Instead it is
   "decomposed" into two orientations: the nominal spin axis
   orientation captured in the segments providing the orientation of
   the JUNO_SPIN_AXIS frame relative to the J2000 frame and the nominal
   rotation about the spin axis captured in the segments providing the
   orientation of the JUNO_SPACECRAFT frame relative to the
   JUNO_SPIN_AXIS frame.

   JUNO_SPIN_AXIS is defined as a CK-based frame.

   \begindata

      FRAME_JUNO_SPIN_AXIS         = -61900
      FRAME_-61900_NAME            = 'JUNO_SPIN_AXIS'
      FRAME_-61900_CLASS           = 3
      FRAME_-61900_CLASS_ID        = -61900
      FRAME_-61900_CENTER          = -61
      CK_-61900_SCLK               = -61
      CK_-61900_SPK                = -61

   \begintext


Magnetometer Frames
-------------------------------------------------------------------------------

   The set of frames for the Magnetometer experiment includes 

      -  two fluxgate magnetometer (FGM) frames (JUNO_FGM_IB and
         JUNO_FGM_OB)

      -  two magnetometer optical bench (MOB) frames (JUNO_MOB_IB and
         JUNO_MOB_OB)
 
      -  four Advanced Stellar Compass (ASC) Camera Head Unit (CHU)
         frames (JUNO_ASC_CHUA, JUNO_ASC_CHUB, JUNO_ASC_CHUC, and
         JUNO_ASC_CHUD.

   The set also includes an additional set of fixed-offset frames
   incorporating relative alignment data for 

      -  FGM_## relative to CHU# (JUNO_CHU#_TO_FGM_##)

      -  FGM_## relative to MOB_## (JUNO_MOB_##_TO_FGM_##)

      -  spacecraft relative to CHU# (JUNO_CHU#_TO_SC)
 
      -  spacecraft relative to MOB_##  (JUNO_MOB_##_TO_SC)

      -  MOB_## relative to spacecraft (JUNO_SC_TO_MOB_##)

   The "alignment" frames facilitate the ability to easily change the
   values of alignments and to switch different alignment sets by
   storing alignments directly in the FK (instead of CKs) while using
   "coverage" CKs to store zero-offset rotations from these alignment
   frames to the desired structure frames by simply constraining
   coverages available to SPICE.


FGM Frames
 
   The FGM frames -- JUNO_FGM_IB and JUNO_FGM_OB -- are the frames with
   respect to which the FGM sensors were calibrated. These frames are
   defined by the orthogonal faces of the alignment cubes mounted on
   the optical benches as follows:

      -  +Z axis is the normal to the top side of the alignment cube
         (the side parallel to the optical bench plane) and pointing
         away from the bench; nominally points in the same direction as
         the s/c +Z axis.

      -  +X axis is the normal to the side of the cube facing the FGM
         sensor mounted on the bench; for OB FGM nominally points in
         the same direction as the s/c +X axis, for IB FGM nominally
         points in the same direction as the s/c -X axis
         
      -  +Y axis completes the right-handed frame

      -  the origin of the frame is in the geometric center of the FGM
         enclosures

   The OB FGM frame is nominally co-aligned with the spacecraft frame.

   The IB FGM frame is nominally rotated from the spacecraft frame by
   180 degrees about Z.

   This diagram illustrates the FGM frames:

      Spacecraft -Z side view:
      ------------------------


                         Magnetometer Boom

                 /```---...          
                 |         ```---.. +Yfgm_ib
                 |                        ^ .            
                 |                        |  ```---...       
                 |                        |   CHUD   CHUB...
    +Zsc    +Xsc |                      .-|--@.         .@----.|    +Xfgm_ob
      x------>   |                 <------x @ |         | @ x------>
      |          |          +Xfgm_ib    `----@'         `@--|-'|
      |          |                            CHUC   CHUA'''|
      |          |                           ...---'''      |
      v          |                       /'''               V
       +Ysc      |         ...---'''\   /                  +Yfgm_ob 
                 \...---'''          ---


                                           +Zsc, +Zfgm_ib, +Zfgm_ob, 
                                               are into the page.


   The FGM frames are defined as CK-based frames because their
   orientation can be provided with respect to either one of the CHU
   frames (of the two CHUs mounted on the bench) or with respect to the
   bench's MOB frame.

   \begindata

      FRAME_JUNO_FGM_IB            = -61114
      FRAME_-61114_NAME            = 'JUNO_FGM_IB'
      FRAME_-61114_CLASS           = 3
      FRAME_-61114_CLASS_ID        = -61114
      FRAME_-61114_CENTER          = -61
      CK_-61114_SCLK               = -61
      CK_-61114_SPK                = -61

      FRAME_JUNO_FGM_OB            = -61124
      FRAME_-61124_NAME            = 'JUNO_FGM_OB'
      FRAME_-61124_CLASS           = 3
      FRAME_-61124_CLASS_ID        = -61124
      FRAME_-61124_CENTER          = -61
      CK_-61124_SCLK               = -61
      CK_-61124_SPK                = -61

   \begintext


MOB Frames

   The MOB frames -- JUNO_MOB_IB and JUNO_MOB_OB -- are the principal
   optical bench orientation frames. These frames are not tied to a
   particular physical feature but are instead defined by two vectors
   that are the sum and the difference of the boresight directions of
   the two CHUs mounted on each bench:

      -  +Z axis is along the vector that is the sum of the CHU
         boresight directions (CHUA + CHUB for OB, CHUD + CHUC for IB)

      -  +Y axis is along the vector that is the difference of the
         boresight directions (CHUA - CHUB for OB, CHUD - CHUC for IB)

      -  +X axis completes the right-handed frame

      -  the origin of the frame is between the focal points of the two
         CHUs whose boresights define the frame.

   The OB MOB frame is nominally rotated from the s/c frame by 180
   degrees about X, then by 180 degrees about Z. 

   The IB MOB nominally rotated from the s/c frame by 180 degrees about
   X.

   This diagram illustrates the MOB frames:

      Spacecraft -Z side view:
      ------------------------


                         Magnetometer Boom

                 /```---...          ___     
                 |         ```---.../   \  +Ymob_ib
                 |                       \..^              
                 |                          |```---...       
                 |                          | CHUD   CHUB...
    +Zsc    +Xsc |                      .---|@.  +X     .@----.. 
      x------>   |               FGM_IB | @ o-----> <-----o @ ||FGM_OB
      |          |                      `----@'     +X  `@|---'' 
      |          |                            CHUC   CHUA'|'
      |          |                           ...---'''    |    
      v          |                       /'''             V    
       +Ysc      |         ...---'''\   /               +Ymob_ob
                 \...---'''          ---
 
                                           +Zsc is into the page.

                                           +Zmob_ib, and +Zmob_ob
                                            are out of the page.


   The MOB frames are defined as CK-based frames because their
   orientation is provided in CKs with respect to the J2000 frame.

   \begindata

      FRAME_JUNO_MOB_IB            = -61113
      FRAME_-61113_NAME            = 'JUNO_MOB_IB'
      FRAME_-61113_CLASS           = 3
      FRAME_-61113_CLASS_ID        = -61113
      FRAME_-61113_CENTER          = -61
      CK_-61113_SCLK               = -61
      CK_-61113_SPK                = -61

      FRAME_JUNO_MOB_OB            = -61123
      FRAME_-61123_NAME            = 'JUNO_MOB_OB'
      FRAME_-61123_CLASS           = 3
      FRAME_-61123_CLASS_ID        = -61123
      FRAME_-61123_CENTER          = -61
      CK_-61123_SCLK               = -61
      CK_-61123_SPK                = -61

   \begintext


CHU frames

   The ASC CHU frames -- JUNO_ASC_CHUA, JUNO_ASC_CHUB, JUNO_ASC_CHUC,
   and JUNO_ASC_CHUD -- are defined as follows:

      -  +Z axis is along the boresight

      -  +X axis is along the sensor side roughly aligned with the cable
         direction, pointing away from the cable side.

      -  +Y axis completes the right-handed frame

      -  the origin of the frame is at the CHU's focal point

   The mapping between CHUs and their telemetry designations (0-3) is 
   as follows:

      CHUA (-61122) -> TLM "0"
      CHUB (-61121) -> TLM "1"
      CHUC (-61112) -> TLM "2"
      CHUD (-61111) -> TLM "3"

   These diagrams illustrate the ASC CHU frames:

      Spacecraft -Z side view:
      ------------------------


                        Magnetometer Boom
                 /```---...          ___
                 |         ```---.../   \    ^ +YchuD
                 |                       \...|
                 |                           ^ +YchuC 
                 |                      CHUD |        ```- CHUB
    +Zsc    +Xsc |                      .--- o----> <----o----..
      x------>   |                FGM_IB| @ @|| +X   +X ||@ @ ||FGM_OB
      |          |                      `--- o----> <----o----''
      |          |                      CHUC          ...| CHUA
      |          |                           ...- +YchuB V 
      v          |                       /'''            |
       +Ysc      |         ...---'''\   /         +YchuA V 
                 \...---'''          ---

                                               +Zsc is into the page.

                                        +Zchu1, +Zchu2, +Zchu3, and +Zchu4
                                         are out of the page, inclined at
                                                   ~13 degrees.

   This diagram illustrates all magnetometer frames together.

      Spacecraft +X side view:
      ------------------------


            Outboard Bench                          Inboard Bench
            --------------                          -------------

                    +Z            +Z                   +Z
                     ^             ^                    ^
                     |             |                    |
                     | FGM_OB      | S/C         FGM_IB |
      +X     +Y      |             |                    |     +Y     +X
       x------>      o------>      o------>      <------x      <------o
       |            +X     +Y     +X     +Y     +Y     +X             |  
       | MOB_OB                                                MOB_IB |  
       |                                                              |
       v                                                              v
      +Z                                                           +Z
                           .> +Y               +Y <.
                    +X  .-'                         `-. +X
            +X x      x'  CHUA                         `o      o +X
              / `-.    \                               /    .-' \
             /     `>   \                       CHUD  /   <'     \
            / CHUB  +Y   \                           /  +Y  CHUC  \
           v              v                         v              v
         +Z              +Z                       +Z               +Z

         /<--------|-------->\                    /<--------|-------->\
            13 deg | 13 deg                          13 deg | 13 deg



       +Xs of s/c, FGM_OB, MOB_IB, CHUC, and CHUD are out of the page.
       +Xs of      FGM_IB, MOB_OB, CHUA, and CHUB are into   the page.


   The ASC CHU frames are defined as CK-based frames because their
   orientation is provided either with respect to the J2000 frame,
   based on the CHU's attitude solution, or with respect to the
   spacecraft frame, based on the in-flight relative alignment
   calibration.

   \begindata

      FRAME_JUNO_ASC_CHUD          = -61111
      FRAME_-61111_NAME            = 'JUNO_ASC_CHUD'
      FRAME_-61111_CLASS           = 3
      FRAME_-61111_CLASS_ID        = -61111
      FRAME_-61111_CENTER          = -61
      CK_-61111_SCLK               = -61
      CK_-61111_SPK                = -61

      FRAME_JUNO_ASC_CHUC          = -61112
      FRAME_-61112_NAME            = 'JUNO_ASC_CHUC'
      FRAME_-61112_CLASS           = 3
      FRAME_-61112_CLASS_ID        = -61112
      FRAME_-61112_CENTER          = -61
      CK_-61112_SCLK               = -61
      CK_-61112_SPK                = -61

      FRAME_JUNO_ASC_CHUB          = -61121
      FRAME_-61121_NAME            = 'JUNO_ASC_CHUB'
      FRAME_-61121_CLASS           = 3
      FRAME_-61121_CLASS_ID        = -61121
      FRAME_-61121_CENTER          = -61
      CK_-61121_SCLK               = -61
      CK_-61121_SPK                = -61

      FRAME_JUNO_ASC_CHUA          = -61122
      FRAME_-61122_NAME            = 'JUNO_ASC_CHUA'
      FRAME_-61122_CLASS           = 3
      FRAME_-61122_CLASS_ID        = -61122
      FRAME_-61122_CENTER          = -61
      CK_-61122_SCLK               = -61
      CK_-61122_SPK                = -61

   \begintext


Alignment Frames

   In the FK versions 0.5 (November 2011) to 0.8 (up to April 2016) 
   the JUNO_CHU#_TO_FGM_## frames below incorporated the following CHU#
   to FGM_## alignments based on the "Micro Advanced Stellar Compass
   JUNO Reference Cube Calibration Test Report"
   (JN-DTU-TR-3008_v1_0_ref_cube_calibration_report-2.pdf):
 
      --------------------------------------------------------
                          Rot3 X       Rot2 Y       Rot1 Z
      --------------------------------------------------------
      FGM_OB -> CHUA   166.98004591  0.83914744  -179.77479685
      FGM_OB -> CHUB  -167.65012849  0.22864645  -179.36168926
      FGM_IB -> CHUC  -166.89008185  0.16658654   179.10988258
      FGM_IB -> CHUD   166.92526364  0.41005034   178.56549108
      --------------------------------------------------------

   These alignments were derived by creating the transformation
   matrices from the CHU# frames to the FGM frames that have their +X
   axes (primary axes) along the vectors that are the negative of the
   measured optical cube X side normals given in the CHU# frames and
   their +Y axes (plane defining axes) along the vectors that are the
   measured optical cube Y side normals given in the CHU# frames.

   These angles were provided in the frame definitions using the
   following keywords:

      tkframe_-61131_angles      = ( 166.98004591, 0.83914744, -179.77479685 )
      tkframe_-61132_angles      = ( -167.65012849, 0.22864645, -179.36168926 )
      tkframe_-61133_angles      = ( -166.89008185, 0.16658654, 179.10988258 )
      tkframe_-61134_angles      = ( 166.92526364, 0.41005034, 178.56549108 )


   Starting with the FK version 0.9 (April 2016)
   the JUNO_CHU#_TO_FGM_## frames below incorporate the following CHU#
   to FGM_## alignments from [10] table 4:

      --------------------------------------------------------
                          Rot3 X       Rot2 Y       Rot1 Z
      --------------------------------------------------------
      FGM_OB -> CHUA   166.981554    0.867687    -179.789384
      FGM_OB -> CHUB  -167.680175    0.231891    -179.384575
      FGM_IB -> CHUC  -166.964941    0.211558     179.078816
      FGM_IB -> CHUD   166.887958    0.457271     178.543535
      --------------------------------------------------------

   \begindata

      FRAME_JUNO_CHUA_TO_FGM_OB      = -61131
      FRAME_-61131_NAME            = 'JUNO_CHUA_TO_FGM_OB'
      FRAME_-61131_CLASS           = 4
      FRAME_-61131_CLASS_ID        = -61131
      FRAME_-61131_CENTER          = -61
      TKFRAME_-61131_SPEC          = 'ANGLES'
      TKFRAME_-61131_RELATIVE      = 'JUNO_ASC_CHUA'
      TKFRAME_-61131_ANGLES      = (  166.981554, 0.867687, -179.789384 )
      TKFRAME_-61131_AXES          = ( 1, 2, 3 )
      TKFRAME_-61131_UNITS         = 'DEGREES'

      FRAME_JUNO_CHUB_TO_FGM_OB      = -61132
      FRAME_-61132_NAME            = 'JUNO_CHUB_TO_FGM_OB'
      FRAME_-61132_CLASS           = 4
      FRAME_-61132_CLASS_ID        = -61132
      FRAME_-61132_CENTER          = -61
      TKFRAME_-61132_SPEC          = 'ANGLES'
      TKFRAME_-61132_RELATIVE      = 'JUNO_ASC_CHUB'
      TKFRAME_-61132_ANGLES      = ( -167.680175, 0.231891, -179.384575 )
      TKFRAME_-61132_AXES          = ( 1, 2, 3 )
      TKFRAME_-61132_UNITS         = 'DEGREES'

      FRAME_JUNO_CHUC_TO_FGM_IB      = -61133
      FRAME_-61133_NAME            = 'JUNO_CHUC_TO_FGM_IB'
      FRAME_-61133_CLASS           = 4
      FRAME_-61133_CLASS_ID        = -61133
      FRAME_-61133_CENTER          = -61
      TKFRAME_-61133_SPEC          = 'ANGLES'
      TKFRAME_-61133_RELATIVE      = 'JUNO_ASC_CHUC'
      TKFRAME_-61133_ANGLES      = ( -166.964941, 0.211558, 179.078816 )
      TKFRAME_-61133_AXES          = ( 1, 2, 3 )
      TKFRAME_-61133_UNITS         = 'DEGREES'

      FRAME_JUNO_CHUD_TO_FGM_IB      = -61134
      FRAME_-61134_NAME            = 'JUNO_CHUD_TO_FGM_IB'
      FRAME_-61134_CLASS           = 4
      FRAME_-61134_CLASS_ID        = -61134
      FRAME_-61134_CENTER          = -61
      TKFRAME_-61134_SPEC          = 'ANGLES'
      TKFRAME_-61134_RELATIVE      = 'JUNO_ASC_CHUD'
      TKFRAME_-61134_ANGLES      = (  166.887958, 0.457271, 178.543535 )
      TKFRAME_-61134_AXES          = ( 1, 2, 3 )
      TKFRAME_-61134_UNITS         = 'DEGREES'

   \begintext

   In the FK versions 0.5 (November 2011) to 0.8 (up to April 2016) 
   the JUNO_MOB_##_TO_FGM_## frames below incorporated the following
   MOB_## to FGM_## alignments based on the "Micro Advanced Stellar
   Compass JUNO Reference Cube Calibration Test Report"
   (JN-DTU-TR-3008_v1_0_ref_cube_calibration_report-2.pdf):
 
      -----------------------------------------------------------
                            Rot3 X        Rot2 Y        Rot1 Z
      -----------------------------------------------------------
      FGM_OB -> MOB_OB   179.65236785   0.57188198   179.06999179
      FGM_IB -> MOB_IB  -179.98504194   0.35177733   178.31458807
      -----------------------------------------------------------

   These angles were provided in the frame definitions using the
   following keywords:

      tkframe_-61135_angles      = ( 179.65236785, 0.57188198, 179.06999179 )
      tkframe_-61136_angles      = ( -179.98504194, 0.35177733, 178.31458807 )


   Starting with the FK version 0.9 (April 2016)
   the JUNO_CHU#_TO_FGM_## frames below incorporate the following CHU#
   to FGM_## alignments based on [10]:

      -----------------------------------------------------------
                            Rot3 X        Rot2 Y        Rot1 Z
      -----------------------------------------------------------
      FGM_OB -> MOB_OB   179.63717677   0.58603721   178.99314287
      FGM_IB -> MOB_IB   179.95841684   0.39617415   178.28125011
      -----------------------------------------------------------

   These alignments were derived by creating the transformation
   matrices from the FGM_## frames to the MOB_## frames that have their
   +Z axes (primary axes) along the sums of a particular MOB's CHU#
   boresight vectors and their +Y axes (plane defining axes) along the
   differences of a particular MOB's CHU# boresight vectors.

   \begindata

      FRAME_JUNO_MOB_OB_TO_FGM_OB      = -61135
      FRAME_-61135_NAME            = 'JUNO_MOB_OB_TO_FGM_OB'
      FRAME_-61135_CLASS           = 4
      FRAME_-61135_CLASS_ID        = -61135
      FRAME_-61135_CENTER          = -61
      TKFRAME_-61135_SPEC          = 'ANGLES'
      TKFRAME_-61135_RELATIVE      = 'JUNO_MOB_OB'
      TKFRAME_-61135_ANGLES      = ( 179.63717677, 0.58603721, 178.99314287 )
      TKFRAME_-61135_AXES          = ( 1, 2, 3 )
      TKFRAME_-61135_UNITS         = 'DEGREES'

      FRAME_JUNO_MOB_IB_TO_FGM_IB      = -61136
      FRAME_-61136_NAME            = 'JUNO_MOB_IB_TO_FGM_IB'
      FRAME_-61136_CLASS           = 4
      FRAME_-61136_CLASS_ID        = -61136
      FRAME_-61136_CENTER          = -61
      TKFRAME_-61136_SPEC          = 'ANGLES'
      TKFRAME_-61136_RELATIVE      = 'JUNO_MOB_IB'
      TKFRAME_-61136_ANGLES      = ( 179.95841684, 0.39617415, 178.28125011 )
      TKFRAME_-61136_AXES          = ( 1, 2, 3 )
      TKFRAME_-61136_UNITS         = 'DEGREES'

   \begintext

   In the FK versions 0.5 (November 2011) to 0.8 (up to April 2016) 
   the JUNO_CHU#_TO_SC frames below incorporated the following CHU# to
   SC alignments based on a quick comparison of the CHU and s/c
   orientations during the ASC turn-on activity on August 25, 2011:

      -------------------------------------------------------
                      Rot3 Z          Rot2 Y          Rot1 X
      -------------------------------------------------------
      S/C->CHUA       179.215          0.207         -167.634
      S/C->CHUB       178.862         -0.010          167.032
      S/C->CHUC        -0.724         -0.633         -166.482
      S/C->CHUD        -0.460         -0.608          167.377
      -------------------------------------------------------

   These angles were provided in the frame definitions using the
   following keywords:

      tkframe_-61141_angles        = ( 179.215, 0.207, -167.634 )
      tkframe_-61142_angles        = ( 178.862, -0.010, 167.032 )
      tkframe_-61143_angles        = ( -0.724, -0.633, -166.482 )
      tkframe_-61144_angles        = ( -0.460, -0.608, 167.377 )


   Starting with the FK version 0.9, draft 2 (April 2016)
   the JUNO_CHU#_TO_SC frames below incorporate the following CHU# to
   SC alignments based on a quick comparison of the CHU and s/c
   orientations during January-March 2016:

      -------------------------------------------------------
                      Rot3 Z          Rot2 Y          Rot1 X
      -------------------------------------------------------
      S/C->CHUA       178.970          1.290         -167.630
      S/C->CHUB       179.100          1.065          167.035
      S/C->CHUC        -0.980          0.375         -166.475
      S/C->CHUD        -0.240          0.400          167.385
      -------------------------------------------------------

   These angles were provided in the frame definitions using the
   following keywords:

      tkframe_-61141_angles        = ( 178.970, 1.290, -167.630 )
      tkframe_-61142_angles        = ( 179.100, 1.065, 167.035 )
      tkframe_-61143_angles        = ( -0.980, 0.375, -166.475 )
      tkframe_-61144_angles        = ( -0.240, 0.400, 167.385 )

   Starting with the FK version 0.9, draft 5 (August 2016)
   the JUNO_CHU#_TO_SC frames below incorporate the following CHU# to
   SC alignments based on a quick comparison of the CHU and s/c
   orientations for August 6-19, 2016:

      -------------------------------------------------------
                      Rot3 Z          Rot2 Y          Rot1 X
      -------------------------------------------------------
      S/C->CHUA       178.950          1.370         -167.635
      S/C->CHUB       179.125          1.150          167.035
      S/C->CHUC        -1.000          0.480         -166.480
      S/C->CHUD        -0.220          0.510          167.380
      -------------------------------------------------------

   WARNING: the fixed alignments below are just a snapshot of the actual
   alignments at a particular point in time. The actual alignments proved
   to be significantly varying due a temperature-dependent solar array 
   flexing, both during cruise and during orbital operations, with the 
   Y rotation changing back and forth by as much as 0.08 deg during each 
   peri-jove pass. Because of this variability the frames below and any 
   frame chains containing them should NOT be used in science data 
   analysis.

   \begindata

      FRAME_JUNO_CHUA_TO_SC        = -61141
      FRAME_-61141_NAME            = 'JUNO_CHUA_TO_SC'
      FRAME_-61141_CLASS           = 4
      FRAME_-61141_CLASS_ID        = -61141
      FRAME_-61141_CENTER          = -61
      TKFRAME_-61141_SPEC          = 'ANGLES'
      TKFRAME_-61141_RELATIVE      = 'JUNO_ASC_CHUA'
      TKFRAME_-61141_ANGLES        = ( 178.950, 1.370, -167.635 )
      TKFRAME_-61141_AXES          = ( 3, 2, 1 )
      TKFRAME_-61141_UNITS         = 'DEGREES'

      FRAME_JUNO_CHUB_TO_SC        = -61142
      FRAME_-61142_NAME            = 'JUNO_CHUB_TO_SC'
      FRAME_-61142_CLASS           = 4
      FRAME_-61142_CLASS_ID        = -61142
      FRAME_-61142_CENTER          = -61
      TKFRAME_-61142_SPEC          = 'ANGLES'
      TKFRAME_-61142_RELATIVE      = 'JUNO_ASC_CHUB'
      TKFRAME_-61142_ANGLES        = ( 179.125, 1.150, 167.035 )
      TKFRAME_-61142_AXES          = ( 3, 2, 1 )
      TKFRAME_-61142_UNITS         = 'DEGREES'

      FRAME_JUNO_CHUC_TO_SC        = -61143
      FRAME_-61143_NAME            = 'JUNO_CHUC_TO_SC'
      FRAME_-61143_CLASS           = 4
      FRAME_-61143_CLASS_ID        = -61143
      FRAME_-61143_CENTER          = -61
      TKFRAME_-61143_SPEC          = 'ANGLES'
      TKFRAME_-61143_RELATIVE      = 'JUNO_ASC_CHUC'
      TKFRAME_-61143_ANGLES        = ( -1.000, 0.480, -166.480 )
      TKFRAME_-61143_AXES          = ( 3, 2, 1 )
      TKFRAME_-61143_UNITS         = 'DEGREES'

      FRAME_JUNO_CHUD_TO_SC        = -61144
      FRAME_-61144_NAME            = 'JUNO_CHUD_TO_SC'
      FRAME_-61144_CLASS           = 4
      FRAME_-61144_CLASS_ID        = -61144
      FRAME_-61144_CENTER          = -61
      TKFRAME_-61144_SPEC          = 'ANGLES'
      TKFRAME_-61144_RELATIVE      = 'JUNO_ASC_CHUD'
      TKFRAME_-61144_ANGLES        = ( -0.220, 0.510, 167.380 )
      TKFRAME_-61144_AXES          = ( 3, 2, 1 )
      TKFRAME_-61144_UNITS         = 'DEGREES'

   \begintext

   In the FK versions 0.5 (November 2011) to 0.8 (up to April 2016) 
   the JUNO_MOB_##_TO_SC and JUNO_SC_TO_MOB_## frames below incorporated
   the following MOB_## to SC alignments based on a quick comparison of
   the CHU and s/c orientations during the ASC turn-on activity on
   August 25, 2011:

      -------------------------------------------------------
                        Rot3 X        Rot2 Y        Rot1 Z
      -------------------------------------------------------
      SC -> MOB_OB  -179.69984062   0.10355166   179.50575594
      SC -> MOB_IB  -179.55188944   0.63655035     0.06025323
      -------------------------------------------------------

      -------------------------------------------------------
                        Rot3 Z        Rot2 Y        Rot1 X
      -------------------------------------------------------
      MOB_OB -> SC  -179.50575594  -0.10355166   179.69984062
      MOB_IB -> SC    -0.06025323  -0.63655035   179.55188944
      -------------------------------------------------------

   These angles were provided in the frame definitions using the
   following keywords:

      tkframe_-61145_angles      = ( -179.69984062, 0.10355166, 179.50575594 )
      tkframe_-61146_angles      = ( -179.55188944, 0.63655035, 0.06025323 )
      tkframe_-61155_angles      = ( -179.50575594, -0.10355166, 179.69984062 )
      tkframe_-61156_angles      = ( -0.06025323, -0.63655035, 179.55188944 )


   Starting with the FK version 0.9, draft 2 (April 2016)
   the JUNO_MOB_##_TO_SC and JUNO_SC_TO_MOB_## frames below incorporate
   the following MOB_## to SC alignments based on a quick comparison of
   the CHU and s/c orientations during January-March 2016:

      -------------------------------------------------------
                        Rot3 X        Rot2 Y        Rot1 Z
      -------------------------------------------------------
      SC -> MOB_OB  -179.71273540   1.20947262   179.49305746
      SC -> MOB_IB  -179.54535339  -0.39823128     0.05211686
      -------------------------------------------------------

      -------------------------------------------------------
                        Rot3 Z        Rot2 Y        Rot1 X
      -------------------------------------------------------
      MOB_OB -> SC  -179.49305746  -1.20947262   179.71273540
      MOB_IB -> SC    -0.05211686   0.39823128   179.54535339
      -------------------------------------------------------

   These angles were provided in the frame definitions using the
   following keywords:

      tkframe_-61145_angles      = ( -179.71273540, 1.20947262, 179.49305746 )
      tkframe_-61146_angles      = ( -179.54535339, -0.39823128, 0.05211686 )
      tkframe_-61155_angles      = ( -179.49305746, -1.20947262, 179.71273540 )
      tkframe_-61156_angles      = ( -0.05211686,  0.39823128, 179.54535339 )


   Starting with the FK version 0.9, draft 5 (August 2016)
   the JUNO_MOB_##_TO_SC and JUNO_SC_TO_MOB_## frames below incorporate
   the following MOB_## to SC alignments based on a quick comparison of
   the CHU and s/c orientations for August 6-19, 2016:

      -------------------------------------------------------
                        Rot3 X        Rot2 Y        Rot1 Z
      -------------------------------------------------------
      SC -> MOB_OB  -179.71070379   1.29397715   179.50482902
      SC -> MOB_IB  -179.55054081  -0.50866857     0.06234192
      -------------------------------------------------------

      -------------------------------------------------------
                        Rot3 Z        Rot2 Y        Rot1 X
      -------------------------------------------------------
      MOB_OB -> SC  -179.50482902  -1.29397715   179.71070379
      MOB_IB -> SC    -0.06234192   0.50866857   179.55054081
      -------------------------------------------------------

   These alignments were derived by creating the transformation
   matrices from the SC frame to the MOB_## frames that have their +Z
   axes (primary axes) along the sums of a particular MOB's CHU#
   boresight vectors and their +Y axes (plane defining axes) along the
   differences of a particular MOB's CHU# boresight vectors.

   WARNING: the fixed alignments below are just a snapshot of the actual
   alignments at a particular point in time. The actual alignments proved
   to be significantly varying due a temperature-dependent solar array 
   flexing, both during cruise and during orbital operations, with the 
   Y rotation changing back and forth by as much as 0.08 deg during each 
   peri-jove pass. Because of this variability the frames below and any 
   frame chains containing them should NOT be used in science data 
   analysis.

   \begindata

      FRAME_JUNO_MOB_OB_TO_SC      = -61145
      FRAME_-61145_NAME            = 'JUNO_MOB_OB_TO_SC'
      FRAME_-61145_CLASS           = 4
      FRAME_-61145_CLASS_ID        = -61145
      FRAME_-61145_CENTER          = -61
      TKFRAME_-61145_SPEC          = 'ANGLES'
      TKFRAME_-61145_RELATIVE      = 'JUNO_MOB_OB'
      TKFRAME_-61145_ANGLES      = ( -179.71070379, 1.29397715, 179.50482902 )
      TKFRAME_-61145_AXES          = ( 1, 2, 3 )
      TKFRAME_-61145_UNITS         = 'DEGREES'

      FRAME_JUNO_MOB_IB_TO_SC      = -61146
      FRAME_-61146_NAME            = 'JUNO_MOB_IB_TO_SC'
      FRAME_-61146_CLASS           = 4
      FRAME_-61146_CLASS_ID        = -61146
      FRAME_-61146_CENTER          = -61
      TKFRAME_-61146_SPEC          = 'ANGLES'
      TKFRAME_-61146_RELATIVE      = 'JUNO_MOB_IB'
      TKFRAME_-61146_ANGLES      = ( -179.55054081, -0.50866857, 0.06234192 )
      TKFRAME_-61146_AXES          = ( 1, 2, 3 )
      TKFRAME_-61146_UNITS         = 'DEGREES'

      FRAME_JUNO_SC_TO_MOB_OB      = -61155
      FRAME_-61155_NAME            = 'JUNO_SC_TO_MOB_OB'
      FRAME_-61155_CLASS           = 4
      FRAME_-61155_CLASS_ID        = -61155
      FRAME_-61155_CENTER          = -61
      TKFRAME_-61155_SPEC          = 'ANGLES'
      TKFRAME_-61155_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61155_ANGLES      = ( -179.50482902, -1.29397715, 179.71070379 )
      TKFRAME_-61155_AXES          = ( 3, 2, 1 )
      TKFRAME_-61155_UNITS         = 'DEGREES'

      FRAME_JUNO_SC_TO_MOB_IB      = -61156
      FRAME_-61156_NAME            = 'JUNO_SC_TO_MOB_IB'
      FRAME_-61156_CLASS           = 4
      FRAME_-61156_CLASS_ID        = -61156
      FRAME_-61156_CENTER          = -61
      TKFRAME_-61156_SPEC          = 'ANGLES'
      TKFRAME_-61156_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61156_ANGLES      = ( -0.06234192, 0.50866857, 179.55054081 )
      TKFRAME_-61156_AXES          = ( 3, 2, 1 )
      TKFRAME_-61156_UNITS         = 'DEGREES'

   \begintext


JADE Frames
-------------------------------------------------------------------------------

   The set of frames for the JADE experiment includes three Electron
   Sensor frames -- JUNO_JADE_E060, JUNO_JADE_E180, and JUNO_JADE_E300
   -- and one Ion Sensor frame -- JUNO_JADE_I.

   The JADE Elector Sensor frames -- JUNO_JADE_E060, JUNO_JADE_E180, and
   JUNO_JADE_E300 -- are fixed w.r.t. to the spacecraft and defined as
   follows:

      -  +Z axis is normal to the plane between upper and lower
         entrance deflectors ("principal" plane) and points away from
         the sensor base

      -  +X axis is in the "principal" plane nominally points towards the
         center of the sensor FOV

      -  +Y axis completes the right-handed frame

      -  the origin of the frame is at the intersection of the "principal"
         plane and the sensor symmetry axis.

   This diagram illustrates the JADE Electron sensor frames:

      Spacecraft +Z side view:
      ------------------------

            Solar Array 3
           ~ ~ ~ ~ ~ ~ ~ ~
            \   .-'        \      +Xe060
             \-'         .-'\           ^    <-. Spin Direction
              \       .-'              /        `.
               \   .-'        <.      /
                \-'     +Ye060  `-.  /
                 \       .-'       `o.+Ze060
                  \   .-'    ----- .  `-.          Solar Array 1
                   .-'    +Ysc ^    `.   `-.--------------------.
                   |    /      |      \    |      |      |      |
                               |       \   |      |      |      |``--..
         +Xe180     +Ze180     |    +Xsc   |      |      |      |      `--..
            <------o           o------> |  |      |      |      |           |
                   |  .      +Zsc       ,  |      |      |      |     ..--''
                   |   \               /   |      |      |      |..--'
                   |    \ HGA         /           |      |      | Magnetometer
                   V .   `.         .'   .>  -------------------'     Boom
             +Ye180   `-.  ` ----- '  .-' +Ye300
                 /       `-  +Ze300 o'
                /           `-. .-'  \
               /`-.           /'      \
              /    `-.       /         \
             /        `-.   /           v
            /`-.         `-/      +Xe300
           /~ ~ ~ ~ ~ ~ ~ ~                       +Zsc, +Zje060, +Zje180,
            Solar Array 2                          +Zje300 are out of
                                                         the page.

   As seen on the diagram the JUNO_JADE_E060, JUNO_JADE_E180, and
   JUNO_JADE_E300 frames are nominally rotated from the spacecraft
   frame by +60, 180 and -60 degrees about +Z axis correspondingly.
   These nominal rotations were provided in the FK versions 0.0-1.0
   as

      tkframe_-61201_angles        = ( -60.0, 0.0, 0.0 )
      tkframe_-61201_axes          = (   3,   2,   1   )

      tkframe_-61202_angles        = ( 180.0, 0.0, 0.0 )
      tkframe_-61202_axes          = (   3,   2,   1   )

      tkframe_-61203_angles        = (  60.0, 0.0, 0.0 )
      tkframe_-61203_axes          = (   3,   2,   1   )

   Based on the JADE_E* anode look directions from [13], the
   JUNO_JADE_E060 frame is rotated from the spacecraft frame first by
   +60 deg about Z then by +0.35 deg about Y, JUNO_JADE_E180 -- first
   by 180 deg about Z then by +0.48 deg about Y, and JUNO_JADE_E300 --
   first by -60 degrees about Z then by +0.34 deg about Y.

   Since the frame definitions below contains the reverse
   transformations -- from the JADE Electron Sensor frames to the
   spacecraft frame -- the order of rotations is reversed and the signs
   of rotation angles are changed to the opposite ones compared to the
   description above.

   \begindata

      FRAME_JUNO_JADE_E060         = -61201
      FRAME_-61201_NAME            = 'JUNO_JADE_E060'
      FRAME_-61201_CLASS           = 4
      FRAME_-61201_CLASS_ID        = -61201
      FRAME_-61201_CENTER          = -61
      TKFRAME_-61201_SPEC          = 'ANGLES'
      TKFRAME_-61201_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61201_ANGLES        = ( -60.0, -0.35, 0.0 )
      TKFRAME_-61201_AXES          = (   3,    2,    1   )
      TKFRAME_-61201_UNITS         = 'DEGREES'

      FRAME_JUNO_JADE_E180         = -61202
      FRAME_-61202_NAME            = 'JUNO_JADE_E180'
      FRAME_-61202_CLASS           = 4
      FRAME_-61202_CLASS_ID        = -61202
      FRAME_-61202_CENTER          = -61
      TKFRAME_-61202_SPEC          = 'ANGLES'
      TKFRAME_-61202_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61202_ANGLES        = ( 180.0, -0.48, 0.0 )
      TKFRAME_-61202_AXES          = (   3,    2,    1   )
      TKFRAME_-61202_UNITS         = 'DEGREES'

      FRAME_JUNO_JADE_E300         = -61203
      FRAME_-61203_NAME            = 'JUNO_JADE_E300'
      FRAME_-61203_CLASS           = 4
      FRAME_-61203_CLASS_ID        = -61203
      FRAME_-61203_CENTER          = -61
      TKFRAME_-61203_SPEC          = 'ANGLES'
      TKFRAME_-61203_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61203_ANGLES        = (  60.0, -0.34, 0.0 )
      TKFRAME_-61203_AXES          = (   3,    2,    1   )
      TKFRAME_-61203_UNITS         = 'DEGREES'

   \begintext

   The JADE Ion Sensor frame -- JUNO_JADE_I -- is fixed w.r.t. to the
   spacecraft and defined as follows:

      -  +Z axis is normal to the plane between upper and lower
         entrance deflectors ("principal" plane) and points away from
         the anode

      -  +X axis is in the "principal" plane nominally points towards the
         center of the sensor FOV

      -  +Y axis completes the right-handed frame

      -  the origin of the frame is at the intersection of the "principal"
         plane and the sensor symmetry axis.

   This diagram illustrates the JADE Ion sensor frame:


      Spacecraft +Z side view:
      ------------------------

                Solar Array 3             <-. Spin Direction
               ~ ~ ~ ~ ~ ~ ~ ~                `.
               \   .-'        \.
                \-'         .-' `-.
                 \       .-'       `-.
    (45 deg above \   .-'    ----- .  `-.      Solar Array 1
      the page)    .-'    +Ysc ^    `.   `-.--------------------.
          +Xi    ..@    /      |      \    |      |      |      |
             <-''  |\          |       \   |      |      |      |``--..
          +Yi      | \         |    +Xsc   |      |      |      |      `--..
    (45 deg below  |15\        o------> |  |      |      |      |           |
      the page)    |   v     +Zsc       ,  |      |      |      |     ..--''
                   |  +Zi              /   |      |      |      |..--'
                   |             HGA  /    |      |      |      | Magnetometer
                   `-.   `.         .'   .-'--------------------'     Boom
                  /   `-.  ` ----- '  .-'
                 /       `-.       .-'
                /           `-. .-'
               /`-.           /'
               ~ ~ ~ ~ ~ ~ ~ ~                         +Zsc is out of
               Solar Array 2                              the page.


   As seen on the diagram the JUNO_JADE_I frame is nominally rotated
   from the spacecraft first by -75 degrees about +Z, then by +90
   degrees about +Y, then by -135 degrees about +Z.

   Based on the JADE_I* anode look directions from [13], the nominal
   orientation described above has been used for JADE_I observation
   geometry computations.

   Since the frame definitions below contains the reverse
   transformation -- from the JADE Ion Sensor frame to the spacecraft
   frame -- the order of rotations is reversed and the signs of
   rotation angles are changed to the opposite ones compared to the
   description above.

   \begindata

      FRAME_JUNO_JADE_I            = -61204
      FRAME_-61204_NAME            = 'JUNO_JADE_I'
      FRAME_-61204_CLASS           = 4
      FRAME_-61204_CLASS_ID        = -61204
      FRAME_-61204_CENTER          = -61
      TKFRAME_-61204_SPEC          = 'ANGLES'
      TKFRAME_-61204_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61204_ANGLES        = ( 75.0, -90.0, 135.0 )
      TKFRAME_-61204_AXES          = (  3,     2,     3   )
      TKFRAME_-61204_UNITS         = 'DEGREES'

   \begintext


JEDI Frames
-------------------------------------------------------------------------------

   The set of frames for the JEDI experiment includes three JEDI Sensor
   frames -- JUNO_JEDI_090, JUNO_JEDI_A180, and JUNO_JEDI_270, -- fixed
   w.r.t. to the spacecraft and defined as follows:

      -  +Z axis is along the sensor head symmetry axis, pointing away
         from the sensor electronics box

      -  +X axis is in the FOV center plane, nominally pointing towards
         the center of the full (160 x 12) sensor FOV

      -  +Y axis completes the right-handed frame

      -  the origin of the frame is at the intersection of the FOV
         center plane and the sensor head symmetry axis.

   Note that for the JEDI/090 and JEDI/270 sensors the frame is as
   defined in the JEDI MICD. For the JEDI/A180 sensor the frame is not
   as defined in the JEDI MICD: it is rotated by -90 degrees about +X
   to align the +Z axis with the sensor's symmetry axis.

   This diagram illustrates the JEDI sensor frames:

      Spacecraft +Z side view:
      ------------------------

           Solar Array 3
          ~ ~ ~ ~ ~ ~ ~ ~
           \       .-'    \       +X090
            \   .-'        \    ^
             \-'         .-'\   |         <-. Spin Direction
              \              \  |            `.
               \   .  +Y090   \ |
                \-'      <------o +Z090
                 \       .-        `-.
                  \   .-'    ----- .  `-.      Solar Array 1
                   .-'    +Ysc ^    `.   `-.--------------------.
                   |    /      |      \    |      |      |      |
                   |   /       |       \   |      |      |      |``--..
                   |  '        |    +Xsc   |      |      |      |      `--..
                               o------> |  |      |      |      |           |
         +Xa180   +Ya180     +Zsc       ,  |      |      |      |     ..--''
            <------x                   /   |      |      |      |..--'
                   |    \ HGA         /    |      |      |      | Magnetometer
                   | .   `.         .'   .-'--------------------'     Boom
                   |  `-.  ` ----- '  .-'
                   V       +Z270   .-'
             +Za180          -. o------>         +X090 and +X270 are out
                              / |     +Y270      of the page, ~10 degrees
              /    `-.       /  |                above the s/c XY plane.
             /        `-.   /   |
            /`-.         `-/    V                +Z090 and +Z270 are out
           /    `-.       /      +X270          of the page, tilted by ~13
           ~ ~ ~ ~ ~ ~ ~ ~                    degrees towards the s/c +Z axis.
            Solar Array 2
                                                 +Y090 is out of the page,
                                               ~8 deg above the s/c XY plane.

                                                 +Y270 is into the page,
                                               ~8 deg below the s/c XY plane.

                                                 +Yja180 is into the page.

                                                 +Zsc is out of the page.

   The JUNO_JEDI_090 frame is nominally rotated from the spacecraft
   frame first by +90 degrees about +Z, then by +8 degrees about +X,
   then by -10 degrees about +Y.
 
   The JUNO_JEDI_A180 frame is nominally rotated from the spacecraft
   frame first by 180 degrees about +Z, then by -90 degrees about +X.
 
   The JUNO_JEDI_270 frame is nominally rotated from the spacecraft
   frame first by -90 degrees about +Z, then by -8 degrees about +X,
   then by -10 degrees about +Y.

   Since the frame definitions below contains the reverse
   transformations -- from the JEDI Sensor frames to the spacecraft
   frame -- the order of rotations is reversed and the signs of
   rotation angles are changed to the opposite ones compared to the
   description above.

   \begindata

      FRAME_JUNO_JEDI_090          = -61301
      FRAME_-61301_NAME            = 'JUNO_JEDI_090'
      FRAME_-61301_CLASS           = 4
      FRAME_-61301_CLASS_ID        = -61301
      FRAME_-61301_CENTER          = -61
      TKFRAME_-61301_SPEC          = 'ANGLES'
      TKFRAME_-61301_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61301_ANGLES        = ( -90.0, -8.0, 10.0 )
      TKFRAME_-61301_AXES          = (   3,    1,    2   )
      TKFRAME_-61301_UNITS         = 'DEGREES'

      FRAME_JUNO_JEDI_A180         = -61302
      FRAME_-61302_NAME            = 'JUNO_JEDI_A180'
      FRAME_-61302_CLASS           = 4
      FRAME_-61302_CLASS_ID        = -61302
      FRAME_-61302_CENTER          = -61
      TKFRAME_-61302_SPEC          = 'ANGLES'
      TKFRAME_-61302_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61302_ANGLES        = ( 180.0, 0.0, 90.0 )
      TKFRAME_-61302_AXES          = (   3,   2,    1   )
      TKFRAME_-61302_UNITS         = 'DEGREES'

      FRAME_JUNO_JEDI_270          = -61303
      FRAME_-61303_NAME            = 'JUNO_JEDI_270'
      FRAME_-61303_CLASS           = 4
      FRAME_-61303_CLASS_ID        = -61303
      FRAME_-61303_CENTER          = -61
      TKFRAME_-61303_SPEC          = 'ANGLES'
      TKFRAME_-61303_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61303_ANGLES        = ( 90.0, 8.0, 10.0 )
      TKFRAME_-61303_AXES          = (  3,    1,   2   )
      TKFRAME_-61303_UNITS         = 'DEGREES'

   \begintext


JIRAM Frames:
-------------------------------------------------------------------------------

   The set of frames for the JIRAM experiment includes the unit
   reference frame (URF) -- JUNO_JIRAM_URF, -- the imager and
   spectrometer frames -- JUNO_JIRAM_I and JUNO_JIRAM_S, -- and the
   imager L-band and M-band channel frames -- JUNO_JIRAM_I_LBAND and
   JUNO_JIRAM_I_MBAND. This set of frames does not attempt to account
   for the de-spinning mirror motion and assumes that the de-spinning
   compensation results in effectively inertially fixed instrument
   pointing for a given observation.

   The JIRAM unit reference frame -- JUNO_JIRAM_URF -- is the JIRAM
   "hardware" frame. It is fixed w.r.t. to the spacecraft and is
   defined as follows:

      -  +Z axis is normal to the JIRAM optical head base, points from
         the base towards the spacecraft, and is nominally co-aligned
         with the spacecraft +Z axis;

      -  +X axis is along the nominal FOV center axis

      -  +Y axis completes the right-handed frame

      -  the origin of the frame is at center of the JIRAM optical
         head reference mounting hole.

   This diagram illustrates the JIRAM unit reference frame:

      Spacecraft +Z side view:
      ------------------------

                Solar Array 3             <-. Spin Direction
               ~ ~ ~ ~ ~ ~ ~ ~               `.
               \   .-'        \.
      -.        \-'         .-' `-.
        ` +Xurf  \       .-'       `-.
      22    <-.   \   .-'    ----- .  `-.      Solar Array 1
      deg      `-. .-'    +Ysc ^    `.   `-.--------------------.
      ---         `o    /      |      \    |      |      |      |
                  /|+Zurf      |       \   |      |      |      |``--..
                 / |  '        |    +Xsc   |      |      |      |      `--..
                /  |  |        o------> |  |      |      |      |           |
               V   |  .      +Zsc       ,  |      |      |      |     ..--''
          +Yurf    |   \               /   |      |      |      |..--'
                   |    \ HGA         /    |      |      |      | Magnetometer
                   `-.   `.         .'   .-'--------------------'     Boom
                  /   `-.  ` ----- '  .-'
                 /       `-.       .-'
                /           `-. .-'
               /`-.           /'
               ~ ~ ~ ~ ~ ~ ~ ~
               Solar Array 2                  +Zurf and +Zsc are out of
                                                     the page.

   As seen on the diagram the JUNO_JIRAM_URF frame is nominally rotated
   from the spacecraft frame by +158 degrees about +Z.

   Since the frame definition below contain the reverse transformation
   -- from the JIRAM unit reference frame to the spacecraft frame --
   the order of rotations is reversed and the signs of rotation angles
   are changed to the opposite ones compared to the description above.

   \begindata

      FRAME_JUNO_JIRAM_URF         = -61401
      FRAME_-61401_NAME            = 'JUNO_JIRAM_URF'
      FRAME_-61401_CLASS           = 4
      FRAME_-61401_CLASS_ID        = -61401
      FRAME_-61401_CENTER          = -61
      TKFRAME_-61401_SPEC          = 'ANGLES'
      TKFRAME_-61401_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61401_ANGLES        = ( -158.0, 0.0, 0.0 )
      TKFRAME_-61401_AXES          = (    3,   2,   1   )
      TKFRAME_-61401_UNITS         = 'DEGREES'

   \begintext

   The JIRAM imager and spectrometer frames -- JUNO_JIRAM_I and
   JUNO_JIRAM_S, -- and the imager L-band and M-band channel frames --
   JUNO_JIRAM_I_LBAND and JUNO_JIRAM_I_MBAND -- are the "image" frames
   defined as follows:

      -  +Z axis is along the boresight (view direction of the center
         pixel of the spectrometer line, combined imager CCD or individual
         L-band and M-band CCDs);

      -  +X axis is along the CCD columns (spectral direction of the
         spectrometer), pointing in the "along-track" direction resulting
         from the spacecraft rotation

      -  +Y axis completes the right-handed frame, pointing along the
         spatial spectrometer direction.

      -  the origin of the frame is the instrument's focal point.

   This diagram illustrates the spectrometer and imager frames:

      Spacecraft +Z side view:
      ------------------------

                Solar Array 3             <-. Spin Direction
               ~ ~ ~ ~ ~ ~ ~ ~               `.
               \   .-'        \.
                  '         .-' `-.
          +Zs/i/m/l      .-'       `-.
            <-.       .-'    ----- .  `-.      Solar Array 1
               `-. .-'    +Ysc ^    `.   `-.--------------------.
                  `o    /      |      \    |      |      |      |
                  /|+Ys/i/m/l  |       \   |      |      |      |``--..
                 / |  '        |    +Xsc   |      |      |      |      `--..
                /  |  |        o------> |  |      |      |      |           |
               V   |  .      +Zsc       ,  |      |      |      |     ..--''
      +Xs/i/m/l    |   \               /   |      |      |      |..--'
                   |    \ HGA         /    |      |      |      | Magnetometer
                   `-.   `.         .'   .-'--------------------'     Boom
                  /   `-.  ` ----- '  .-'
                 /       `-.       .-'
                /           `-. .-'
               /`-.           /'
               ~ ~ ~ ~ ~ ~ ~ ~
               Solar Array 2                  +Ys/i/m/l and +Zsc are out of
                                                      the page.

   The spectrometer frame is defined w.r.t. to the unit reference
   frame; the imager frame is defined w.r.t. to the spectrometer frame;
   the imager L-band and M-band channel frames are defined w.r.t. to
   the imager frame.

   Based on [7] the JIRAM spectrometer frame is rotated from the unit
   reference frame first by +90 degrees about +Y, then by +90 degrees
   about +Z, then by +6.125 degrees about +Y.

   The JIRAM imager frame is rotated from the spectrometer frame by
   +0.22689 degrees about +Y.

   The JIRAM imager L-band frame is rotated from the imager frame by
   +0.96257 degrees (1/2 of L-band CCD) about +Y.

   The JIRAM imager M-band frame is rotated from the imager frame by
   -0.93507 degrees (1/2 of M-band CCD) about +Y.

   The rotation angles above were used in the FK version 0.6 and were
   incorporated into the frame definitions using the following 
   keywords:

      TKFRAME_-61420_ANGLES        = ( -90.0, -90.0, -6.125 ) 
      TKFRAME_-61410_ANGLES        = ( 0.0, -0.22689, 0.0 ) 
      TKFRAME_-61411_ANGLES        = ( 0.0, -0.96257, 0.0 )  
      TKFRAME_-61412_ANGLES        = ( 0.0, +0.93507, 0.0 ) 


   Based on [8], which provided a summary of alignment changes that
   resulted from star and Moon calibrations and from correcting pixel
   offsets between detector and band centers, the JIRAM spectrometer
   frame is rotated from the unit reference frame first by +90.063
   degrees about +Y, then by +90 degrees about +Z, then by +6.2765
   degrees about +Y.

   The JIRAM imager frame is rotated from the spectrometer frame by
   +0.23844 degrees about +Y due to 17.5 pixel offset between the
   center of spectrometer and imager.

   The JIRAM imager L-band frame is rotated from the imager frame by
   +0.94012 degrees (1/2 of L-band CCD) about +Y due to 69 pixel offset
   between L center and I center.

   The JIRAM imager M-band frame is rotated from the imager frame by
   -0.94012 degrees (1/2 of M-band CCD) about +Y due to 69 pixel offset
   between M center and I center.

   Since the frame definitions below contain the reverse
   transformations, the order of rotations is reversed and the signs of
   rotation angles are changed to the opposite ones compared to the
   description above.

   \begindata

      FRAME_JUNO_JIRAM_S           = -61420
      FRAME_-61420_NAME            = 'JUNO_JIRAM_S'
      FRAME_-61420_CLASS           = 4
      FRAME_-61420_CLASS_ID        = -61420
      FRAME_-61420_CENTER          = -61
      TKFRAME_-61420_SPEC          = 'ANGLES'
      TKFRAME_-61420_RELATIVE      = 'JUNO_JIRAM_URF'
      TKFRAME_-61420_ANGLES        = ( -90.063, -90.0,  -6.2765 )
      TKFRAME_-61420_AXES          = (   2,       3,     2      )
      TKFRAME_-61420_UNITS         = 'DEGREES'

      FRAME_JUNO_JIRAM_I           = -61410
      FRAME_-61410_NAME            = 'JUNO_JIRAM_I'
      FRAME_-61410_CLASS           = 4
      FRAME_-61410_CLASS_ID        = -61410
      FRAME_-61410_CENTER          = -61
      TKFRAME_-61410_SPEC          = 'ANGLES'
      TKFRAME_-61410_RELATIVE      = 'JUNO_JIRAM_S'
      TKFRAME_-61410_ANGLES        = ( 0.0, -0.23844, 0.0 )
      TKFRAME_-61410_AXES          = ( 3,    2,       1   )
      TKFRAME_-61410_UNITS         = 'DEGREES'

      FRAME_JUNO_JIRAM_I_LBAND     = -61411
      FRAME_-61411_NAME            = 'JUNO_JIRAM_I_LBAND'
      FRAME_-61411_CLASS           = 4
      FRAME_-61411_CLASS_ID        = -61411
      FRAME_-61411_CENTER          = -61
      TKFRAME_-61411_SPEC          = 'ANGLES'
      TKFRAME_-61411_RELATIVE      = 'JUNO_JIRAM_I'
      TKFRAME_-61411_ANGLES        = ( 0.0, -0.94012, 0.0 )
      TKFRAME_-61411_AXES          = ( 3,    2,       1   )
      TKFRAME_-61411_UNITS         = 'DEGREES'

      FRAME_JUNO_JIRAM_I_MBAND     = -61412
      FRAME_-61412_NAME            = 'JUNO_JIRAM_I_MBAND'
      FRAME_-61412_CLASS           = 4
      FRAME_-61412_CLASS_ID        = -61412
      FRAME_-61412_CENTER          = -61
      TKFRAME_-61412_SPEC          = 'ANGLES'
      TKFRAME_-61412_RELATIVE      = 'JUNO_JIRAM_I'
      TKFRAME_-61412_ANGLES        = ( 0.0, +0.94012, 0.0 )
      TKFRAME_-61412_AXES          = ( 3,    2,       1   )
      TKFRAME_-61412_UNITS         = 'DEGREES'

   \begintext


JUNOCAM Frames
-------------------------------------------------------------------------------

   The JUNOCAM frame -- JUNO_JUNOCAM -- is fixed w.r.t. to the
   spacecraft and defined as follows:

      -  +Z axis is along the camera boresight;

      -  +X axis is along the CCD lines, pointing in the increasing
         pixel direction

      -  +Y axis completes the right-handed frame

      -  the origin of the frame is at the camera focal point.

   This diagram illustrates the JUNOCAM frame:

      Spacecraft +Z side view:
      ------------------------

                Solar Array 3             <-. Spin Direction
               ~ ~ ~ ~ ~ ~ ~ ~               `.
               \   .-'        \.
                \-'         .-' `-.
                 \       .-'       `-.
                  \   .-'    ----- .  `-.      Solar Array 1
                   .-'    +Ysc ^    `.   `-.--------------------.
                   |    /      |      \    |      |      |      |
                   |   /       |       \   |      |      |      |``--..
                      '        |    +Xsc   |      |      |      |      `--..
                 +Xjc |        o------> |  |      |      |      |           |
            <------x  .      +Zsc       ,  |      |      |      |     ..--''
          +Zjc     |   \               /   |      |      |      |..--'
                   |    \ HGA         /    |      |      |      | Magnetometer
                   | .   `.         .'   .-'--------------------'     Boom
                   V  `-.  ` ----- '  .-'
                 /  +Yjc `-.       .-'
                /           `-. .-'                +Xjc is into
               /`-.           /'                     the page.
               ~ ~ ~ ~ ~ ~ ~ ~
               Solar Array 2                       +Zsc is out of
                                                     the page.

   As seen on the diagram the JUNO_JUNOCAM frame is nominally rotated
   from the spacecraft frame first by 180 degrees about +Z, then by +90
   degrees about +Y.

   In order to incorporate the pre-launch calibrated misalignment
   between the spacecraft frame and JUNOCAM optical cube, specified in
   [6] as

      DCM - SMRF to JUNOCam (*)

         -0.0059163 -0.0142817 -0.9998805
          0.0023828 -0.9998954  0.0142678
         -0.9999797 -0.0022981  0.0059497

      (*) Measured from the instrument alignment cube.
          No correction from cube axes to instrument axes

   and the pre-launch calibrated misalignment between optical cube and
   the CCD specified in [6] as CW rotation by 0.69 degrees about +Z
   axis of the camera frame as two separate rotations, an additional
   frame -- JUNO_JUNOCAM_CUBE -- was added between the JUNO_SPACECRAFT
   frame and the JUNO_JUNOCAM frame in the FK ver. 0.6. The rotation
   between these two frames was implemented in FKs up to ver. 1.1 using
   these keywords:

      tkframe_-61500_axes          = (  3,     2,    1   )
      tkframe_-61500_angles        = (  0.69,  0.0,  0.0 )

   Based on the initial reanalysis of Junocam cruise star imaging (per
   [14]) the angles aligning the JUNO_JUNOCAM_CUBE frame with the
   JUNO_JUNOCAM frame were changed to: the first rotation by -0.69
   degrees about Z, the second rotation by +0.469 degrees about Y, and
   the third rotation by -0.583 degrees about X.
 
   The JUNO_JUNOCAM_CUBE misalignment from [6] and the JUNO_JUNOCAM
   from [14] are incorporated in the two frame definitions below.

   \begindata

      FRAME_JUNO_JUNOCAM_CUBE      = -61505
      FRAME_-61505_NAME            = 'JUNO_JUNOCAM_CUBE'
      FRAME_-61505_CLASS           = 4
      FRAME_-61505_CLASS_ID        = -61505
      FRAME_-61505_CENTER          = -61
      TKFRAME_-61505_SPEC          = 'MATRIX'
      TKFRAME_-61505_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61505_MATRIX        = ( 
                                       -0.0059163 -0.0142817 -0.9998805
                                        0.0023828 -0.9998954  0.0142678
                                       -0.9999797 -0.0022981  0.0059497
                                     )

      FRAME_JUNO_JUNOCAM           = -61500
      FRAME_-61500_NAME            = 'JUNO_JUNOCAM'
      FRAME_-61500_CLASS           = 4
      FRAME_-61500_CLASS_ID        = -61500
      FRAME_-61500_CENTER          = -61
      TKFRAME_-61500_SPEC          = 'ANGLES'
      TKFRAME_-61500_RELATIVE      = 'JUNO_JUNOCAM_CUBE'
      TKFRAME_-61500_ANGLES        = (  0.69,  -0.469,  0.583 )
      TKFRAME_-61500_AXES          = (  3,      2,      1     )
      TKFRAME_-61500_UNITS         = 'DEGREES'

   \begintext


MWR Frames:
-------------------------------------------------------------------------------

   The set of frames for the MWR experiment includes six MWR antenna
   frames -- JUNO_MWR_A1, JUNO_MWR_A2, JUNO_MWR_A3, JUNO_MWR_A4,
   JUNO_MWR_A5, and JUNO_MWR_A6, -- fixed w.r.t. to the spacecraft and
   defined as follows:

      -  +Z axis is along the antenna boresight

      -  +Y axis is nominally along the spacecraft +Z axis

      -  +X axis completes the right-handed frame

      -  the origin of the frame is at the geometric center of the
         antenna patch or outer rim.

   This diagram illustrates the MWR antenna frames:

      Spacecraft +Z side view:
      ------------------------

              Solar Array 3              ^        <-. Spin Direction
               ~ ~ ~ ~ ~ ~ ~   +Xa1     / +Za1       `.
               \   .-'        <.       /
                \-'         .   `-.   /
                 \       .-'       `-o +Ya1
                  \   .-'    ----- .   -.      Solar Array 1
                   .-'    +Ysc ^    `.   `-.--------------------.
                   |    /      |      \    |      |      |      |
                   |   /       |       \   |      |      |      |``--..
                   |  '        |    +Xsc   |      |      |      |      `--..
                   |  |        o------> |  |      |      |      |           |
                   |  .      +Zsc       ,  |      |      |      |     ..--''
                   |   \               /                 |      |..--'
                   |    \ HGA         /     +Xa2/3/4/5/6 |      | Magnetometer
                   `-.   `.         .'   .->   -----------------'     Boom
                  /   `-.  ` ----- '  .-'
                 /       `-.         o +Ya2/3/4/5/6
                /           `-. .-'   \
               /`-.           /'       \
               ~ ~ ~ ~ ~ ~ ~ ~          \ +Za2/3/4/5/6
               Solar Array 2             V
                                                   +Ya* and +Zsc are
                                                    out of the page.


   The JUNO_MWR_A1 frame is nominally rotated from the spacecraft frame
   first by +150 degrees about +Z, then by +90 degrees about +X.

   The JUNO_MWR_A2..A6 frames are nominally rotated from the spacecraft
   frame first by +30 degrees about +Z, then by +90 degrees about +X.

   Since the frame definitions below contain the reverse
   transformations, the order of rotations is reversed and the signs of
   rotation angles are changed to the opposite ones compared to the
   description above.

   \begindata

      FRAME_JUNO_MWR_A1            = -61601
      FRAME_-61601_NAME            = 'JUNO_MWR_A1'
      FRAME_-61601_CLASS           = 4
      FRAME_-61601_CLASS_ID        = -61601
      FRAME_-61601_CENTER          = -61
      TKFRAME_-61601_SPEC          = 'ANGLES'
      TKFRAME_-61601_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61601_ANGLES        = ( -150.0, 0.0, -90.0 )
      TKFRAME_-61601_AXES          = (    3,   2,     1   )
      TKFRAME_-61601_UNITS         = 'DEGREES'

      FRAME_JUNO_MWR_A2            = -61602
      FRAME_-61602_NAME            = 'JUNO_MWR_A2'
      FRAME_-61602_CLASS           = 4
      FRAME_-61602_CLASS_ID        = -61602
      FRAME_-61602_CENTER          = -61
      TKFRAME_-61602_SPEC          = 'ANGLES'
      TKFRAME_-61602_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61602_ANGLES        = ( -30.0, 0.0, -90.0 )
      TKFRAME_-61602_AXES          = (   3,   2,     1   )
      TKFRAME_-61602_UNITS         = 'DEGREES'

      FRAME_JUNO_MWR_A3            = -61603
      FRAME_-61603_NAME            = 'JUNO_MWR_A3'
      FRAME_-61603_CLASS           = 4
      FRAME_-61603_CLASS_ID        = -61603
      FRAME_-61603_CENTER          = -61
      TKFRAME_-61603_SPEC          = 'ANGLES'
      TKFRAME_-61603_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61603_ANGLES        = ( -30.0, 0.0, -90.0 )
      TKFRAME_-61603_AXES          = (   3,   2,     1   )
      TKFRAME_-61603_UNITS         = 'DEGREES'

      FRAME_JUNO_MWR_A4            = -61604
      FRAME_-61604_NAME            = 'JUNO_MWR_A4'
      FRAME_-61604_CLASS           = 4
      FRAME_-61604_CLASS_ID        = -61604
      FRAME_-61604_CENTER          = -61
      TKFRAME_-61604_SPEC          = 'ANGLES'
      TKFRAME_-61604_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61604_ANGLES        = ( -30.0, 0.0, -90.0 )
      TKFRAME_-61604_AXES          =   ( 3,   2,     1   )
      TKFRAME_-61604_UNITS         = 'DEGREES'

      FRAME_JUNO_MWR_A5            = -61605
      FRAME_-61605_NAME            = 'JUNO_MWR_A5'
      FRAME_-61605_CLASS           = 4
      FRAME_-61605_CLASS_ID        = -61605
      FRAME_-61605_CENTER          = -61
      TKFRAME_-61605_SPEC          = 'ANGLES'
      TKFRAME_-61605_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61605_ANGLES        = ( -30.0, 0.0, -90.0 )
      TKFRAME_-61605_AXES          = (   3,   2,     1   )
      TKFRAME_-61605_UNITS         = 'DEGREES'

      FRAME_JUNO_MWR_A6            = -61606
      FRAME_-61606_NAME            = 'JUNO_MWR_A6'
      FRAME_-61606_CLASS           = 4
      FRAME_-61606_CLASS_ID        = -61606
      FRAME_-61606_CENTER          = -61
      TKFRAME_-61606_SPEC          = 'ANGLES'
      TKFRAME_-61606_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61606_ANGLES        = ( -30.0, 0.0, -90.0 )
      TKFRAME_-61606_AXES          = (   3,   2,     1   )
      TKFRAME_-61606_UNITS         = 'DEGREES'

   \begintext



UVS Frames
-------------------------------------------------------------------------------

   The set of frames for the UVS experiment includes two frames
   -- JUNO_UVS_BASE and JUNO_UVS.

   The UVS "hardware" frame, JUNO_UVS_BASE, is fixed w.r.t. to the
   spacecraft and defined as follows:

      -  +Z axis is normal to the sensor mounting plane and points
         into the sensor; it is nominally along the slit

      -  +X axis is parallel to the scan mirror rotation axis and points
         from the sensor electronics towards the entrance baffle

      -  +Y axis completes the right-handed frame

      -  the origin of the frame is at the center of the reference
         mounting hole.

   The UVS "observation" frame, JUNO_UVS, is a CK-based frame defined
   as follows:

      -  +X axis is the instrument boresight

      -  +Z axis is along the slit and for the "zero" scan mirror position
         points in the same direction as the +Z axis of the JUNO_UVS_BASE
         frame

      -  +Y axis completes the right-handed frame

      -  the origin of the frame is at the intersection of the
         reflected boresight and the scan mirror axis.

   This diagram illustrates the UVS frames (the JUNO_UVS frame is shown
   in the "zero" scan mirror position):

      Spacecraft +Z side view:
      ------------------------

                Solar Array 3             <-. Spin Direction
               ~ ~ ~ ~ ~ ~ ~ ~               `.
               \   .-'        \.
                \-'         .-' `-.
                 \       .-'       `-.
                      .-'    ----- .  `-.      Solar Array 1
                +Yu  '    +Ysc ^    `.   `-.--------------------.
                   ^    /      |      \    |      |      |      |
                   | +Xub      |       \   |      |      |      |``--..
                   |^          |    +Xsc   |      |      |      |      `--..
         +Xu       || |        o------> |  |      |      |      |           |
            <------x| .      +Zsc       ,  |      |      |      |     ..--''
                +Zu |  \               /   |      |      |      |..--'
                    x------>          /    |      |      |      | Magnetometer
                 +Zub    +Yub   HGA .'   .-'--------------------'     Boom
                  /   `-.    ----- '  .-'
                 /       `-.       .-'
                /           `-. .-'               +Zub and +Zu are
               /`-.           /'                    into the page.
               ~ ~ ~ ~ ~ ~ ~ ~
               Solar Array 2                       +Zsc is out of
                                                     the page.


      Scan Mirror Plane and Axis (Spacecraft +Z side view:):
      ------------------------------------------------------

                     |
                     | Scan Mirror Axis
                     |

                     ^ +Yu
                .    |
                 `.  |
            +Xu    `.|
              <------x.
    Reflected      +Zu `.    Scan Mirror plane
    boresight            `  in "zero" position

                     ^ +Xub                             ^ +Ysc
                     |                                  |
                     |                                  |
                     |                                  |
                     x------> Yub                       o------>
                   +Zub                              +Zsc       +Xsc


   As seen on the diagram the JUNO_UVS_BASE frame is nominally rotated
   from the spacecraft frame first by +90 degrees about +Z, then by 180
   degrees about +X. This nominal alignment was provided in FK versions
   0.0-0.9 as 

      tkframe_-61700_angles        = ( -90.0, 0.0, 180.0 )

   The calibrated rotations that align the spacecraft frame with the
   JUNO_UVS_BASE frame, provided in [12], are first by +89.855 deg
   about Z, then by -0.625 deg about Y, then by 180 degrees about X.
   These angles are incorporated into the JUNO_UVS_BASE frame
   definition below.

   For a perfect nominal alignment of the boresight and scan mirror
   axis and a perfect nominal 45 degrees position of the mirror plane,
   at the "zero" mirror position the JUNO_UVS frame is simply rotated
   by -90 degrees about +Z of the JUNO_UVS_BASE frame. For any other
   mirror position within the +30/-30 degrees range the two additional
   rotations are needed: first by the "scan angle" about +Y, then by
   the "scan angle" about +X. These rotations will be stored in CK files.

   Since the frame JUNO_UVS_BASE definition below contain the reverse
   transformation -- from the UVS_BASE frame to the spacecraft frame --
   the order of rotations is reversed and the signs of rotation angles
   are changed to the opposite ones compared to the description above.

   \begindata

      FRAME_JUNO_UVS_BASE          = -61700
      FRAME_-61700_NAME            = 'JUNO_UVS_BASE'
      FRAME_-61700_CLASS           = 4
      FRAME_-61700_CLASS_ID        = -61700
      FRAME_-61700_CENTER          = -61
      TKFRAME_-61700_SPEC          = 'ANGLES'
      TKFRAME_-61700_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61700_ANGLES        = ( -89.855, 0.625, 180.0 )
      TKFRAME_-61700_AXES          = (   3,     2,     1     )
      TKFRAME_-61700_UNITS         = 'DEGREES'

      FRAME_JUNO_UVS               = -61701
      FRAME_-61701_NAME            = 'JUNO_UVS'
      FRAME_-61701_CLASS           = 3
      FRAME_-61701_CLASS_ID        = -61701
      FRAME_-61701_CENTER          = -61
      CK_-61701_SCLK               = -61
      CK_-61701_SPK                = -61

   \begintext


WAVES Frames
-------------------------------------------------------------------------------

   The set of frames for the WAVES experiment includes two frames
   -- JUNO_WAVES_MSC and  JUNO_WAVES_ANTENNA.

   The WAVES MSC frame, JUNO_WAVES_MSC, is fixed w.r.t. to the
   spacecraft and defined as follows:

      -  +Z axis is along the MSC center axis and points in the same
         direction as the spacecraft +Z axis

      -  +Y axis is in the mounting plate plane and points away from
         the spacecraft

      -  +X axis completes the right-handed frame

      -  the origin of the frame is at the geometric center of the  MSC
         center axis.

   The WAVES electrical antenna frame, JUNO_WAVES_ANTENNA, is fixed
   w.r.t. to the spacecraft and defined as follows:

      -  +Z axis is normal to the antenna mounting plane and points
         away from the spacecraft; it is nominally co-aligned with the
         spacecraft -Z axis

      -  +Y axis is in the mounting plate plane and points between the
         two antennas is stowed position

      -  +X axis completes the right-handed frame

      -  the origin of the frame is at the geometric center of the
         antenna mounting plate.

   This diagram illustrates the WAVES frames:

      Spacecraft +Z side view:
      ------------------------

                Solar Array 3             <-. Spin Direction
               ~ ~ ~ ~ ~ ~ ~ ~               `.
               \   .-'        \.
                \-'         .-' `-.
                 \       .-'       `-.
                  \   .-'    ----- .  `-.      Solar Array 1
                   .-'    +Ysc ^    `.   `-.--------------------.
                   |    /      |      \    |      |      |      |
                   |   /       |           |      |      |      |``--..
                   |  '        |     +Yant               |      |      `--..
                   |  |        o----> <----x +Zant       |      |           |
                   |  .    +Zsc   +Xsc  ,  |             |      |     ..--''
                   |   \               /   |      |      |      |..--'
                   |    \ HGA         /    |      |      |      | Magnetometer
                   `-.   `.         .'   . V -------------------'     Boom
                  /   `-.  ` ----- '  .-'   +Xant
        Solar    /       `-.         '
       Array 2  /           `  .o +Zmsc
               /`-.         .-'  \
               ~ ~ ~ ~   <-'      \
                       +Xmsc       \                +Zsc is out of
                                    V +Ymsc            the page.

   The JUNO_WAVES_MSC frame is nominally rotated from the spacecraft
   frame by -150 degrees about +Z.

   The JUNO_WAVES_ANTENNA frame is nominally rotated from the spacecraft
   frame first by -90 degrees about +Z, then by 180 degrees about +X.

   Since the frame definitions below contain the reverse
   transformations, the order of rotations is reversed and the signs of
   rotation angles are changed to the opposite ones compared to the
   description above.

   \begindata

      FRAME_JUNO_WAVES_MSC         = -61810
      FRAME_-61810_NAME            = 'JUNO_WAVES_MSC'
      FRAME_-61810_CLASS           = 4
      FRAME_-61810_CLASS_ID        = -61810
      FRAME_-61810_CENTER          = -61
      TKFRAME_-61810_SPEC          = 'ANGLES'
      TKFRAME_-61810_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61810_ANGLES        = ( 150.0, 0.0, 0.0 )
      TKFRAME_-61810_AXES          = (   3,   2,   1   )
      TKFRAME_-61810_UNITS         = 'DEGREES'

      FRAME_JUNO_WAVES_ANTENNA     = -61820
      FRAME_-61820_NAME            = 'JUNO_WAVES_ANTENNA'
      FRAME_-61820_CLASS           = 4
      FRAME_-61820_CLASS_ID        = -61820
      FRAME_-61820_CENTER          = -61
      TKFRAME_-61820_SPEC          = 'ANGLES'
      TKFRAME_-61820_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61820_ANGLES        = ( 90.0, 0.0, 180.0 )
      TKFRAME_-61820_AXES          = (  3,   2,     1   )
      TKFRAME_-61820_UNITS         = 'DEGREES'

   \begintext


Solar Array Frames
-------------------------------------------------------------------------------

   Two frames are defined for each of the three solar arrays.

   The first frame -- a fixed-offset frame named JUNO_SA#_HINGE (where
   # is 1, 2, or 3) -- is fixed w.r.t. to the spacecraft and is defined
   as follows:

      -  +Z axis is co-aligned with the s/c +Z axis

      -  +X axis is along the solar array symmetry axis and points
         from the hinge to the outer side of the array

      -  +Y axis is along the hinge axis and points to complete the
         right-handed frame

      -  the origin of the frame is in the middle of the hinge axis

   The second frame -- a CK-based frame named JUNO_SA# (where # is 1,
   2, or 3) -- rotates about the hinge (+Y) with respect to the
   corresponding JUNO_SA#_HINGE frame and is defined as follows:

      -  +Z axis is along the normal to the array surface on the solar
         cell side

      -  +X axis is along the solar array symmetry axis and points
         from the hinge towards the outer side of the array

      -  +Y axis is along the hinge axis and points to complete the
         right-handed frame

      -  the origin of the frame is in the middle of between two center
         sections of the array

   In the non-articulated ("zero") position the two frames for each of the
   arrays are co-aligned, as shown on this diagram:

      Spacecraft +Z side view:
      ------------------------

                .-\
             .-'   \
          .-'       \ Solar Array 3
       .-'           \
       \           .-'\
        \      ^ +Xsa3
         \   ,  \       \
          \-'    \    .-'\
           \      o +Zsa3 \
        +Ysa3   .-'        \
             <-'         .-'\             <-. Spin Direction
                     ^ +Xsa3h                `.
               \   .  \       \.
                \-'    \    .-' `-.
                 \      o +Zsa3h   `-.
             +Ysa3h   .-'    ----- .  `-. +Ysa1h       +Ysa1    Solar Array 1
                  <.-'    +Ysc ^    `.   ` ^  ---------  ^  ----.
                        /      |      \    |      |      |      |
                   |   /       |       \   |      |      |      |``--..
                   |  '        |    +Xsc   |   +Xsa1h    |    +Xsa1    `--..
                   |  |        o------> |  o----->       o----->            |
                   |  .      +Zsc       , +Zsa1h       +Zsa1    |     ..--''
                   |   \               /   |      |      |      |..--'
                   |    \ HGA         /    |      |      |      | Magnetometer
                   `-.   `.         .'   .-'--------------------'     Boom
                  /   +Zsa2h ----- '  .-'
                 /      o`-.       .-'
                /      /    `->
               +Xsa2h /        +Ysa2h
              /      V       /
             /         -.   /
            /`-. +Zsa2   `-/
           /      o       /
          /      / `-.
         +Xsa2  /     `-> +Ysa2
        /      V
       /         -.   /
      /.           `-/
        `-.         / Solar Array 2              All +Z axes are out of
           `-.     /                                    the page.
              `-. /
                 `


   The keywords below define the solar array frames.

   \begindata

      FRAME_JUNO_SA1_HINGE         = -61010
      FRAME_-61010_NAME            = 'JUNO_SA1_HINGE'
      FRAME_-61010_CLASS           = 4
      FRAME_-61010_CLASS_ID        = -61010
      FRAME_-61010_CENTER          = -61
      TKFRAME_-61010_SPEC          = 'ANGLES'
      TKFRAME_-61010_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61010_ANGLES        = ( 0.0, 0.0, 0.0 )
      TKFRAME_-61010_AXES          = ( 1,   2,   3   )
      TKFRAME_-61010_UNITS         = 'DEGREES'

      FRAME_JUNO_SA1               = -61011
      FRAME_-61011_NAME            = 'JUNO_SA1'
      FRAME_-61011_CLASS           = 3
      FRAME_-61011_CLASS_ID        = -61011
      FRAME_-61011_CENTER          = -61
      CK_-61011_SCLK               = -61
      CK_-61011_SPK                = -61

      FRAME_JUNO_SA2_HINGE         = -61020
      FRAME_-61020_NAME            = 'JUNO_SA2_HINGE'
      FRAME_-61020_CLASS           = 4
      FRAME_-61020_CLASS_ID        = -61020
      FRAME_-61020_CENTER          = -61
      TKFRAME_-61020_SPEC          = 'ANGLES'
      TKFRAME_-61020_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61020_ANGLES        = ( 0.0, 0.0, 0.0 )
      TKFRAME_-61020_AXES          = ( 1,   2,   3   )
      TKFRAME_-61020_UNITS         = 'DEGREES'

      FRAME_JUNO_SA2               = -61021
      FRAME_-61021_NAME            = 'JUNO_SA2'
      FRAME_-61021_CLASS           = 3
      FRAME_-61021_CLASS_ID        = -61021
      FRAME_-61021_CENTER          = -61
      CK_-61021_SCLK               = -61
      CK_-61021_SPK                = -61

      FRAME_JUNO_SA3_HINGE         = -61030
      FRAME_-61030_NAME            = 'JUNO_SA3_HINGE'
      FRAME_-61030_CLASS           = 4
      FRAME_-61030_CLASS_ID        = -61030
      FRAME_-61030_CENTER          = -61
      TKFRAME_-61030_SPEC          = 'ANGLES'
      TKFRAME_-61030_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61030_ANGLES        = ( 0.0, 0.0, 0.0 )
      TKFRAME_-61030_AXES          = ( 1,   2,   3   )
      TKFRAME_-61030_UNITS         = 'DEGREES'

      FRAME_JUNO_SA3               = -61031
      FRAME_-61031_NAME            = 'JUNO_SA3'
      FRAME_-61031_CLASS           = 3
      FRAME_-61031_CLASS_ID        = -61031
      FRAME_-61031_CENTER          = -61
      CK_-61031_SCLK               = -61
      CK_-61031_SPK                = -61

   \begintext


Antenna Frames
-------------------------------------------------------------------------------

   The JUNO HGA, MGA, and "forward" and "aft" LGA antenna frames --
   JUNO_HGA, JUNO_MGA, JUNO_LGA_FORWARD, and JUNO_LGA_AFT -- are
   defined as follows:

      -  +Z axis is along the antenna boresight

      -  +X axis is along the clock reference direction of the antenna
         pattern

      -  +Y axis completes the right-handed frame

      -  the origin of the frame is in the geometric center of the
         antenna's outer rim or patch

   The JUNO "toroid" LGA frame -- JUNO_LGA_TOROID -- is defined to be
   co-aligned with the spacecraft frame.

   All antenna frames are defined as fixed-offset frames with respect
   to the spacecraft frame.

   Nominally all antenna frames (except for JUNO_LGA_AFT) are
   co-aligned with the spacecraft frame. The JUNO_LGA_AFT frame is
   rotated from the spacecraft frame by 180 degrees about +X.

   This diagram illustrates the antenna frames:

      Spacecraft -Y side view:
      ------------------------

                         +Zhga
                               ^
                               | +Zmga,+Zlgaf
                               |    ^ ^
                               |    | |
                         +Yhga x------> +Xhga
                              / \   | |
                  HGA .-------------x-x---->-> +Xmga,+Xlgaf
                       \      +Ymga,+Ylgaf
                        `.           .'
                           `-------'
                           |       |                            Magnetometer
                           |       |                                Boom
      ================o==o==o==o-----------o======o======o======o===========*
          Solar    |           |           |   Solar Array 1
         Array 2   |           |           |
                   |                       |
                   |                       |
                   |      +Zsc ^           |
                   |    +Zlgat |           |
                   |           |           |
                   |           |           |          +Ysc, +Yhga, Ylgaf,
                   `----- +Ysc x------>   -'          and +Ylgat are into
            +Ylgaa o------>  /_____\  +Xsc                  the page.
                   |    +Xlgaa       +Xlgat
                   |                                   +Ylgaa is out of
                   |                                      the page.
                   V +Zlgaa


   According to [11] the measured HGA electrical boresight [in the 
   spacecraft frame, BVS] is:

      [-3.022155372263570e-04 5.638899817258940e-04 0.999999795346908]

   Three rotations are needed to align the HGA +Z axis with this
   direction while keeping the HGA +X axis close to the spacecraft +X
   axis: first by +118.1889804 degrees about Z, second by +0.03665615
   degrees about Y and third by -118.1889804 degrees about Z. These
   rotations are incorporated in reverse order with opposite signs in
   the HGA frame definition below.
 
   All other antenna frame definitions are nominal as described above.
 
   The keywords below define the antenna frames.

   \begindata

      FRAME_JUNO_HGA               = -61040
      FRAME_-61040_NAME            = 'JUNO_HGA'
      FRAME_-61040_CLASS           = 4
      FRAME_-61040_CLASS_ID        = -61040
      FRAME_-61040_CENTER          = -61
      TKFRAME_-61040_SPEC          = 'ANGLES'
      TKFRAME_-61040_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61040_ANGLES        = ( -118.1889804, -0.03665615, 118.1889804 )
      TKFRAME_-61040_AXES          = ( 3,   2,   3   )
      TKFRAME_-61040_UNITS         = 'DEGREES'

      FRAME_JUNO_MGA               = -61050
      FRAME_-61050_NAME            = 'JUNO_MGA'
      FRAME_-61050_CLASS           = 4
      FRAME_-61050_CLASS_ID        = -61050
      FRAME_-61050_CENTER          = -61
      TKFRAME_-61050_SPEC          = 'ANGLES'
      TKFRAME_-61050_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61050_ANGLES        = ( 0.0, 0.0, 0.0 )
      TKFRAME_-61050_AXES          = ( 1,   2,   3   )
      TKFRAME_-61050_UNITS         = 'DEGREES'

      FRAME_JUNO_LGA_FORWARD       = -61061
      FRAME_-61061_NAME            = 'JUNO_LGA_FORWARD'
      FRAME_-61061_CLASS           = 4
      FRAME_-61061_CLASS_ID        = -61061
      FRAME_-61061_CENTER          = -61
      TKFRAME_-61061_SPEC          = 'ANGLES'
      TKFRAME_-61061_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61061_ANGLES        = ( 0.0, 0.0, 0.0 )
      TKFRAME_-61061_AXES          = ( 1,   2,   3   )
      TKFRAME_-61061_UNITS         = 'DEGREES'

      FRAME_JUNO_LGA_AFT           = -61062
      FRAME_-61062_NAME            = 'JUNO_LGA_AFT'
      FRAME_-61062_CLASS           = 4
      FRAME_-61062_CLASS_ID        = -61062
      FRAME_-61062_CENTER          = -61
      TKFRAME_-61062_SPEC          = 'ANGLES'
      TKFRAME_-61062_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61062_ANGLES        = ( 180.0, 0.0, 0.0 )
      TKFRAME_-61062_AXES          = (   1,   2,   3   )
      TKFRAME_-61062_UNITS         = 'DEGREES'

      FRAME_JUNO_LGA_TOROID        = -61063
      FRAME_-61063_NAME            = 'JUNO_LGA_TOROID'
      FRAME_-61063_CLASS           = 4
      FRAME_-61063_CLASS_ID        = -61063
      FRAME_-61063_CENTER          = -61
      TKFRAME_-61063_SPEC          = 'ANGLES'
      TKFRAME_-61063_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61063_ANGLES        = ( 0.0, 0.0, 0.0 )
      TKFRAME_-61063_AXES          = ( 1,   2,   3   )
      TKFRAME_-61063_UNITS         = 'DEGREES'

   \begintext


ACS Sensor Frames
-------------------------------------------------------------------------------

   This section defined frames for ACS sensors -- SRUs and SSSs.


SRU Frames

   The JUNO SRU frames -- JUNO_SRU1 and JUNO_SRU2 -- are defined as
   follows:

      -  +X axis is along the SRU boresight

      -  +Y axis is nominally along the spacecraft +Z axis

      -  +Z axis completes the right-handed frame

      -  the origin of the frame is at the SRU focal point

   Nominally the SRU1 frame is rotated from the spacecraft frame first
   by +90 degrees about X, then by +65 degrees about Y.

   Nominally the SRU2 frame is rotated from the spacecraft frame first
   by +90 degrees about X, then by +55 degrees about X.

   This diagram illustrates the SRU frames:

      Spacecraft +Z side view:
      ------------------------
                                        +Xsru1
             Solar Array 3                ^            <-. Spin Direction
                                         /  10 deg        `.
               ~ ~ ~ ~ ~ ~ ~ ~\         /
                \-'         .-'        /      .> +Xsru2
                 \       .-'   +Ysru1 o.   .-'
                  \   .-'    ---  +Ysru2`o'              Solar Array 1
                   .-'    +Ysc ^          `.`> +Zsru1  ---------.
                   |    /      |      \     `.           |      |
                   |   /       |       \   |  > +Zsru2   |      |``--..
                   |  '        |    +Xsc   |             |      |      `--..
                   |  |        o------> |  |      |      |      |           |
                   |  .      +Zsc       ,  |      |      |      |     ..--''
                   |   \               /   |      |      |      |..--'
                   |    \ HGA         /    |      |      |      | Magnetometer
                   `-.   `.         .'   .-'--------------------'     Boom
                  /   `-.  ` ----- '  .-'
                 /       `-.       .-'
                /           `-. .-'              +Xsru1 and +Xsru2 are
               ~ ~ ~ ~ ~ ~ ~ ~ '                    into the page.

            Solar Array 2                            +Zsc is out of
                                                        the page.

   The following pre-launch calibrated misalignments between the
   spacecraft frame and the SRU frames were provided in [6]:

      DCM - SMRF to SRU1 (*)

          0.4204194  0.9073288  0.0013891
         -0.0007953 -0.0011625  0.9999990
          0.9073296 -0.4204200  0.0002328

      DCM - SMRF to SRU2 (*)

          0.5720852  0.8201927 -0.0015696
          0.0015127  0.0008586  0.9999985
          0.8201928 -0.5720867 -0.0007496

      (*) Measured from the SRU alignment cube. Correction from 
          cube axes (ARF) to instrument axes (BRF) applied via 
          vendor data

   The misalignments from [6] are incorporated in the two frame
   definitions below.

   \begindata

      FRAME_JUNO_SRU1              = -61071
      FRAME_-61071_NAME            = 'JUNO_SRU1'
      FRAME_-61071_CLASS           = 4
      FRAME_-61071_CLASS_ID        = -61071
      FRAME_-61071_CENTER          = -61
      TKFRAME_-61071_SPEC          = 'MATRIX'
      TKFRAME_-61071_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61071_MATRIX        = (  
                                        0.4204194  0.9073288  0.0013891
                                       -0.0007953 -0.0011625  0.9999990
                                        0.9073296 -0.4204200  0.0002328
                                     )

      FRAME_JUNO_SRU2              = -61072
      FRAME_-61072_NAME            = 'JUNO_SRU2'
      FRAME_-61072_CLASS           = 4
      FRAME_-61072_CLASS_ID        = -61072
      FRAME_-61072_CENTER          = -61
      TKFRAME_-61072_SPEC          = 'MATRIX'
      TKFRAME_-61072_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61072_MATRIX        = (  
                                        0.5720852  0.8201927 -0.0015696
                                        0.0015127  0.0008586  0.9999985
                                        0.8201928 -0.5720867 -0.0007496
                                     )

   \begintext


SSS Frames

   The JUNO SSS frames -- JUNO_SSS1 and JUNO_SSS2 -- are defined as
   follows:

      -  +Z axis is along the SSS boresight and is nominally along the
         spacecraft -X axis.

      -  +X axis is nominally along the spacecraft -Z axis

      -  +Y axis completes the right-handed frame

      -  the origin of the frame is at the SSS geometric center

   Nominally the SSS frames are rotated from the spacecraft frame first
   by +90 degrees about Y, then by 180 degrees about X.

   This diagram illustrates the SSS frames:

      Spacecraft +Z side view:
      ------------------------

             Solar Array 3                <-. Spin Direction
                                             `.
               ~ ~ ~ ~ ~ ~ ~ ~\.
                \-'         .-' `-.
                 \       .-'       `-.
                  \   .-'    ----- .  `-.      Solar Array 1
                   .-'    +Ysc ^    `.   `-.--------------------.
     +Zss1                     |      \    |      |      |      |
           <------x +Xsss1     |       \   |      |      |      |``--..
           <------x +Xsss2     |    +Xsc   |      |      |      |      `--..
     +Zsss2       |            o------> |  |      |      |      |           |
                  ||  .      +Zsc       ,  |      |      |      |     ..--''
           +Ysss1 V|   \               /   |      |      |      |..--'
           +Ysss2 V|    \ HGA         /    |      |      |      | Magnetometer
                   `-.   `.         .'   .-'--------------------'     Boom
                  /   `-.  ` ----- '  .-'
                 /       `-.       .-'
                /           `-. .-'             +Xsss1 and +Xsss2 are
               ~ ~ ~ ~ ~ ~ ~ ~/'                    into the page.

             Solar Array 2                          +Zsc is out of
                                                      the page.

   Since the frame definitions below contains the reverse
   transformations -- from the SSS frames to the spacecraft frame --
   the order of rotations is reversed and the signs of rotation angles
   are changed to the opposite ones compared to the description above.

   \begindata

      FRAME_JUNO_SSS1              = -61073
      FRAME_-61073_NAME            = 'JUNO_SSS1'
      FRAME_-61073_CLASS           = 4
      FRAME_-61073_CLASS_ID        = -61073
      FRAME_-61073_CENTER          = -61
      TKFRAME_-61073_SPEC          = 'ANGLES'
      TKFRAME_-61073_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61073_ANGLES        = ( 0.0, -90.0, 180.0 )
      TKFRAME_-61073_AXES          = ( 3,     2,     1   )
      TKFRAME_-61073_UNITS         = 'DEGREES'

      FRAME_JUNO_SSS2              = -61074
      FRAME_-61074_NAME            = 'JUNO_SSS2'
      FRAME_-61074_CLASS           = 4
      FRAME_-61074_CLASS_ID        = -61074
      FRAME_-61074_CENTER          = -61
      TKFRAME_-61074_SPEC          = 'ANGLES'
      TKFRAME_-61074_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61074_ANGLES        = ( 0.0, -90.0, 180.0 )
      TKFRAME_-61074_AXES          = ( 3,     2,     1   )
      TKFRAME_-61074_UNITS         = 'DEGREES'

   \begintext


Truster Frames
-------------------------------------------------------------------------------

   The JUNO REM thruster frames -- JUNO_REM_[FL1,2,3,4] [FA1,2]
   [RL1,2,3,4] [RA1,2] -- are defined as follows:

      -  +Z axis is along the thrust vector.

      -  +X axis is in the plane containing the thrust vector and the
         spacecraft +Z axis and points in the direction of the
         spacecraft +Z axis

      -  +Y axis completes the right-handed frame

      -  the origin of the frame is center of the nozzle outer edge.

   These diagrams approximately illustrate the REM thruster vector
   directions -- +Z axes of the REM frames, -- other axes are not
   shown:

      Spacecraft +Z side view:
      ------------------------

             Solar Array 3                <-. Spin Direction
                               ^ +Zfa1        `.
               ~ ~ ~ ~ ~ ~ ~ ~\|
                 +Zfl4 <------@@@------> +Zfl1
                 \        FL4 FA1 FL1
                  \   .-'    ----- .  `-.      Solar Array 1
                   .-'    +Ysc ^    `.   `-.--------------------.
                   |    /      |      \    |      |      |      |
                   |   /       |       \   |      |      |      |``--..
                   |  '        |    +Xsc   |      |      |      |      `--..
                   |  |        o------> |  |      |      |      |           |
                   |  .      +Zsc       ,  |      |      |      |     ..--''
                   |   \               /   |      |      |      |..--'
                   |    \ HGA         /    |      |      |      | Magnetometer
                   `-.   `.         .'   .-'--------------------'     Boom
                  /   `-.  ` ----- '  .-'
                 /        FL3 FA2 FL2
                /+Zfl3 <------@@@------> +Zfl2
               ~ ~ ~ ~ ~ ~ ~ ~/|
                               v +Zfa2
             Solar Array 2                           +Zsc is out of
                                                      the page.

      Spacecraft -Z side view:
      ------------------------

             Solar Array 2                 -. Spin Direction
                               ^ +Zra2       `.
               ~ ~ ~ ~ ~ ~ ~ ~\|               v
                 +Zrl3 <------@@@------> +Zrl2
                 \        RL3 RA2 RL2
                  \   .-'             `-.      Solar Array 1
                   .-'                   `-.--------------------.
                   |                       |      |      |      |
                   |                       |      |      |      |``--..
                   |        +Zsc    +Xsc   |      |      |      |      `--..
                   |           o------>    |      |      |      |           |
                   |           |           |      |      |      |     ..--''
                   |           |           |      |      |      |..--'
                   |           |           |      |      |      | Magnetometer
                   `-.    +Ysc V         .-'--------------------'     Boom
                  /   `-.             .-'
                 /        RL4 RA1 RL1
                 +Zrl4 <------@@@------> +Zrl1
               ~ ~ ~ ~ ~ ~ ~ ~/|
                               v +Zra1
             Solar Array 3                           +Zsc is into
                                                      the page.

   The angles in the REM frame definitions below are based the
   information from [4].

   \begindata

      FRAME_JUNO_REM_FL1           = -61081
      FRAME_-61081_NAME            = 'JUNO_REM_FL1'
      FRAME_-61081_CLASS           = 4
      FRAME_-61081_CLASS_ID        = -61081
      FRAME_-61081_CENTER          = -61
      TKFRAME_-61081_SPEC          = 'ANGLES'
      TKFRAME_-61081_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61081_ANGLES        = (   22.009651,  -76.551858,  157.431346 )
      TKFRAME_-61081_AXES          = ( 1,   2,   3   )
      TKFRAME_-61081_UNITS         = 'DEGREES'

      FRAME_JUNO_REM_FL2           = -61082
      FRAME_-61082_NAME            = 'JUNO_REM_FL2'
      FRAME_-61082_CLASS           = 4
      FRAME_-61082_CLASS_ID        = -61082
      FRAME_-61082_CENTER          = -61
      TKFRAME_-61082_SPEC          = 'ANGLES'
      TKFRAME_-61082_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61082_ANGLES        = (  -22.009640,  -76.551859, -157.431357 )
      TKFRAME_-61082_AXES          = ( 1,   2,   3   )
      TKFRAME_-61082_UNITS         = 'DEGREES'

      FRAME_JUNO_REM_FL3           = -61083
      FRAME_-61083_NAME            = 'JUNO_REM_FL3'
      FRAME_-61083_CLASS           = 4
      FRAME_-61083_CLASS_ID        = -61083
      FRAME_-61083_CENTER          = -61
      TKFRAME_-61083_SPEC          = 'ANGLES'
      TKFRAME_-61083_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61083_ANGLES        = (  -22.009651,   76.551858,  -22.568654 )
      TKFRAME_-61083_AXES          = ( 1,   2,   3   )
      TKFRAME_-61083_UNITS         = 'DEGREES'

      FRAME_JUNO_REM_FL4           = -61084
      FRAME_-61084_NAME            = 'JUNO_REM_FL4'
      FRAME_-61084_CLASS           = 4
      FRAME_-61084_CLASS_ID        = -61084
      FRAME_-61084_CENTER          = -61
      TKFRAME_-61084_SPEC          = 'ANGLES'
      TKFRAME_-61084_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61084_ANGLES        = (   22.009640,   76.551859,   22.568643 )
      TKFRAME_-61084_AXES          = ( 1,   2,   3   )
      TKFRAME_-61084_UNITS         = 'DEGREES'

      FRAME_JUNO_REM_FA1           = -61085
      FRAME_-61085_NAME            = 'JUNO_REM_FA1'
      FRAME_-61085_CLASS           = 4
      FRAME_-61085_CLASS_ID        = -61085
      FRAME_-61085_CENTER          = -61
      TKFRAME_-61085_SPEC          = 'ANGLES'
      TKFRAME_-61085_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61085_ANGLES        = (   10.000059,    0.000000,   90.000000 )
      TKFRAME_-61085_AXES          = ( 1,   2,   3   )
      TKFRAME_-61085_UNITS         = 'DEGREES'

      FRAME_JUNO_REM_FA2           = -61086
      FRAME_-61086_NAME            = 'JUNO_REM_FA2'
      FRAME_-61086_CLASS           = 4
      FRAME_-61086_CLASS_ID        = -61086
      FRAME_-61086_CENTER          = -61
      TKFRAME_-61086_SPEC          = 'ANGLES'
      TKFRAME_-61086_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61086_ANGLES        = (  -10.000059,    0.000000,  -90.000000 )
      TKFRAME_-61086_AXES          = ( 1,   2,   3   )
      TKFRAME_-61086_UNITS         = 'DEGREES'

      FRAME_JUNO_REM_RL1           = -61091
      FRAME_-61091_NAME            = 'JUNO_REM_RL1'
      FRAME_-61091_CLASS           = 4
      FRAME_-61091_CLASS_ID        = -61091
      FRAME_-61091_CENTER          = -61
      TKFRAME_-61091_SPEC          = 'ANGLES'
      TKFRAME_-61091_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61091_ANGLES        = (  157.990360,  -76.551859,   22.568643 )
      TKFRAME_-61091_AXES          = ( 1,   2,   3   )
      TKFRAME_-61091_UNITS         = 'DEGREES'

      FRAME_JUNO_REM_RL2           = -61092
      FRAME_-61092_NAME            = 'JUNO_REM_RL2'
      FRAME_-61092_CLASS           = 4
      FRAME_-61092_CLASS_ID        = -61092
      FRAME_-61092_CENTER          = -61
      TKFRAME_-61092_SPEC          = 'ANGLES'
      TKFRAME_-61092_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61092_ANGLES        = ( -157.990349,  -76.551858,  -22.568654 )
      TKFRAME_-61092_AXES          = ( 1,   2,   3   )
      TKFRAME_-61092_UNITS         = 'DEGREES'

      FRAME_JUNO_REM_RL3           = -61093
      FRAME_-61093_NAME            = 'JUNO_REM_RL3'
      FRAME_-61093_CLASS           = 4
      FRAME_-61093_CLASS_ID        = -61093
      FRAME_-61093_CENTER          = -61
      TKFRAME_-61093_SPEC          = 'ANGLES'
      TKFRAME_-61093_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61093_ANGLES        = ( -157.990360,   76.551859, -157.431357 )
      TKFRAME_-61093_AXES          = ( 1,   2,   3   )
      TKFRAME_-61093_UNITS         = 'DEGREES'

      FRAME_JUNO_REM_RL4           = -61094
      FRAME_-61094_NAME            = 'JUNO_REM_RL4'
      FRAME_-61094_CLASS           = 4
      FRAME_-61094_CLASS_ID        = -61094
      FRAME_-61094_CENTER          = -61
      TKFRAME_-61094_SPEC          = 'ANGLES'
      TKFRAME_-61094_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61094_ANGLES        = (  157.990349,   76.551858,  157.431346 )
      TKFRAME_-61094_AXES          = ( 1,   2,   3   )
      TKFRAME_-61094_UNITS         = 'DEGREES'

      FRAME_JUNO_REM_RA1           = -61095
      FRAME_-61095_NAME            = 'JUNO_REM_RA1'
      FRAME_-61095_CLASS           = 4
      FRAME_-61095_CLASS_ID        = -61095
      FRAME_-61095_CENTER          = -61
      TKFRAME_-61095_SPEC          = 'ANGLES'
      TKFRAME_-61095_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61095_ANGLES        = (  169.999941,    0.000000,   90.000000 )
      TKFRAME_-61095_AXES          = ( 1,   2,   3   )
      TKFRAME_-61095_UNITS         = 'DEGREES'

      FRAME_JUNO_REM_RA2           = -61096
      FRAME_-61096_NAME            = 'JUNO_REM_RA2'
      FRAME_-61096_CLASS           = 4
      FRAME_-61096_CLASS_ID        = -61096
      FRAME_-61096_CENTER          = -61
      TKFRAME_-61096_SPEC          = 'ANGLES'
      TKFRAME_-61096_RELATIVE      = 'JUNO_SPACECRAFT'
      TKFRAME_-61096_ANGLES        = ( -169.999941,    0.000000,  -90.000000 )
      TKFRAME_-61096_AXES          = ( 1,   2,   3   )
      TKFRAME_-61096_UNITS         = 'DEGREES'

   \begintext


Juno NAIF ID Codes -- Definitions
========================================================================

   This section contains name to NAIF ID mappings for the JUNO mission.
   Once the contents of this file is loaded into the KERNEL POOL, these
   mappings become available within SPICE, making it possible to use
   names instead of ID code in the high level SPICE routine calls.

   Spacecraft:
   -----------

      JUNO                           -61
      JUNO                           -61
      JUNO_SPACECRAFT                -61000
      JUNO_SPACECRAFT_BUS            -61000
      JUNO_SC_BUS                    -61000

   Magnetometers and CHUs:
   -----------------------

      JUNO_ASC_CHUD                  -61111
      JUNO_ASC_CHUC                  -61112
      JUNO_FGM_IB                    -61114
      JUNO_ASC_CHUB                  -61121
      JUNO_ASC_CHUA                  -61122
      JUNO_FGM_OB                    -61124

   JADE:
   -----

      JUNO_JADE_E060                 -61201
      JUNO_JADE_E180                 -61202
      JUNO_JADE_E300                 -61203
      JUNO_JADE_I                    -61204

   JEDI:
   -----

      JUNO_JEDI_090                  -61301
      JUNO_JEDI_A180                 -61302
      JUNO_JEDI_270                  -61303

   JIRAM:
   ------

      JUNO_JIRAM_I                   -61410
      JUNO_JIRAM_I_LBAND             -61411
      JUNO_JIRAM_I_MBAND             -61412
      JUNO_JIRAM_S                   -61420

   JUNOCAM:
   --------

      JUNO_JUNOCAM                   -61500
      JUNO_JUNOCAM_BLUE              -61501
      JUNO_JUNOCAM_GREEN             -61502
      JUNO_JUNOCAM_RED               -61503
      JUNO_JUNOCAM_METHANE           -61504

   MWR:
   ----

      JUNO_MWR_A1                    -61601
      JUNO_MWR_A2                    -61602
      JUNO_MWR_A3                    -61603
      JUNO_MWR_A4                    -61604
      JUNO_MWR_A5                    -61605
      JUNO_MWR_A6                    -61606

   UVS:
   ----

      JUNO_UVS_BASE                  -61700
      JUNO_UVS                       -61701

   WAVES:
   ------

      JUNO_WAVES_MSC                 -61810
      JUNO_WAVES_ANTENNA             -61820

   Solar Arrays:
   -------------

      JUNO_SA1_HINGE                 -61010
      JUNO_SA1                       -61011
      JUNO_SA2_HINGE                 -61020
      JUNO_SA2                       -61021
      JUNO_SA3_HINGE                 -61030
      JUNO_SA3                       -61031

   Antennas:
   ---------

      JUNO_HGA                       -61040
      JUNO_MGA                       -61050
      JUNO_LGA_FORWARD               -61061
      JUNO_LGA_AFT                   -61062
      JUNO_LGA_TOROID                -61063

   ACS Sensors:
   ------------

      JUNO_SRU1                      -61071
      JUNO_SRU2                      -61072

      JUNO_SSS1                      -61073
      JUNO_SSS2                      -61074

   Thrusters:
   ----------

      JUNO_REM_FL1                   -61081
      JUNO_REM_FL2                   -61082
      JUNO_REM_FL3                   -61083
      JUNO_REM_FL4                   -61084
      JUNO_REM_FA1                   -61085
      JUNO_REM_FA2                   -61086

      JUNO_REM_RL1                   -61091
      JUNO_REM_RL2                   -61092
      JUNO_REM_RL3                   -61093
      JUNO_REM_RL4                   -61094
      JUNO_REM_RA1                   -61095
      JUNO_REM_RA2                   -61096

   The mappings summarized in this table are implemented by the keywords
   below.

   \begindata

      NAIF_BODY_NAME += ( 'JUNO'                        )
      NAIF_BODY_CODE += ( -61                           )

      NAIF_BODY_NAME += ( 'JUNO'                        )
      NAIF_BODY_CODE += ( -61                           )

      NAIF_BODY_NAME += ( 'JUNO_SPACECRAFT'             )
      NAIF_BODY_CODE += ( -61000                        )

      NAIF_BODY_NAME += ( 'JUNO_SPACECRAFT_BUS'         )
      NAIF_BODY_CODE += ( -61000                        )

      NAIF_BODY_NAME += ( 'JUNO_SC_BUS'                 )
      NAIF_BODY_CODE += ( -61000                        )

      NAIF_BODY_NAME += ( 'JUNO_ASC_CHUD'               )
      NAIF_BODY_CODE += ( -61111                        )

      NAIF_BODY_NAME += ( 'JUNO_ASC_CHUC'               )
      NAIF_BODY_CODE += ( -61112                        )

      NAIF_BODY_NAME += ( 'JUNO_FGM_IB'                 )
      NAIF_BODY_CODE += ( -61114                        )

      NAIF_BODY_NAME += ( 'JUNO_ASC_CHUB'               )
      NAIF_BODY_CODE += ( -61121                        )

      NAIF_BODY_NAME += ( 'JUNO_ASC_CHUA'               )
      NAIF_BODY_CODE += ( -61122                        )

      NAIF_BODY_NAME += ( 'JUNO_FGM_OB'                 )
      NAIF_BODY_CODE += ( -61124                        )

      NAIF_BODY_NAME += ( 'JUNO_JADE_E060'              )
      NAIF_BODY_CODE += ( -61201                        )

      NAIF_BODY_NAME += ( 'JUNO_JADE_E180'              )
      NAIF_BODY_CODE += ( -61202                        )

      NAIF_BODY_NAME += ( 'JUNO_JADE_E300'              )
      NAIF_BODY_CODE += ( -61203                        )

      NAIF_BODY_NAME += ( 'JUNO_JADE_I'                 )
      NAIF_BODY_CODE += ( -61204                        )

      NAIF_BODY_NAME += ( 'JUNO_JEDI_090'               )
      NAIF_BODY_CODE += ( -61301                        )

      NAIF_BODY_NAME += ( 'JUNO_JEDI_A180'              )
      NAIF_BODY_CODE += ( -61302                        )

      NAIF_BODY_NAME += ( 'JUNO_JEDI_270'               )
      NAIF_BODY_CODE += ( -61303                        )

      NAIF_BODY_NAME += ( 'JUNO_JIRAM_I'                )
      NAIF_BODY_CODE += ( -61410                        )

      NAIF_BODY_NAME += ( 'JUNO_JIRAM_I_LBAND'          )
      NAIF_BODY_CODE += ( -61411                        )

      NAIF_BODY_NAME += ( 'JUNO_JIRAM_I_MBAND'          )
      NAIF_BODY_CODE += ( -61412                        )

      NAIF_BODY_NAME += ( 'JUNO_JIRAM_S'                )
      NAIF_BODY_CODE += ( -61420                        )

      NAIF_BODY_NAME += ( 'JUNO_JUNOCAM'                )
      NAIF_BODY_CODE += ( -61500                        )

      NAIF_BODY_NAME += ( 'JUNO_JUNOCAM_BLUE'           )
      NAIF_BODY_CODE += ( -61501                        )

      NAIF_BODY_NAME += ( 'JUNO_JUNOCAM_GREEN'          )
      NAIF_BODY_CODE += ( -61502                        )

      NAIF_BODY_NAME += ( 'JUNO_JUNOCAM_RED'            )
      NAIF_BODY_CODE += ( -61503                        )

      NAIF_BODY_NAME += ( 'JUNO_JUNOCAM_METHANE'        )
      NAIF_BODY_CODE += ( -61504                        )

      NAIF_BODY_NAME += ( 'JUNO_MWR_A1'                 )
      NAIF_BODY_CODE += ( -61601                        )

      NAIF_BODY_NAME += ( 'JUNO_MWR_A2'                 )
      NAIF_BODY_CODE += ( -61602                        )

      NAIF_BODY_NAME += ( 'JUNO_MWR_A3'                 )
      NAIF_BODY_CODE += ( -61603                        )

      NAIF_BODY_NAME += ( 'JUNO_MWR_A4'                 )
      NAIF_BODY_CODE += ( -61604                        )

      NAIF_BODY_NAME += ( 'JUNO_MWR_A5'                 )
      NAIF_BODY_CODE += ( -61605                        )

      NAIF_BODY_NAME += ( 'JUNO_MWR_A6'                 )
      NAIF_BODY_CODE += ( -61606                        )

      NAIF_BODY_NAME += ( 'JUNO_UVS_BASE'               )
      NAIF_BODY_CODE += ( -61700                        )

      NAIF_BODY_NAME += ( 'JUNO_UVS'                    )
      NAIF_BODY_CODE += ( -61701                        )

      NAIF_BODY_NAME += ( 'JUNO_WAVES_MSC'              )
      NAIF_BODY_CODE += ( -61810                        )

      NAIF_BODY_NAME += ( 'JUNO_WAVES_ANTENNA'          )
      NAIF_BODY_CODE += ( -61820                        )

      NAIF_BODY_NAME += ( 'JUNO_SA1_HINGE'              )
      NAIF_BODY_CODE += ( -61010                        )

      NAIF_BODY_NAME += ( 'JUNO_SA1'                    )
      NAIF_BODY_CODE += ( -61011                        )

      NAIF_BODY_NAME += ( 'JUNO_SA2_HINGE'              )
      NAIF_BODY_CODE += ( -61020                        )

      NAIF_BODY_NAME += ( 'JUNO_SA2'                    )
      NAIF_BODY_CODE += ( -61021                        )

      NAIF_BODY_NAME += ( 'JUNO_SA3_HINGE'              )
      NAIF_BODY_CODE += ( -61030                        )

      NAIF_BODY_NAME += ( 'JUNO_SA3'                    )
      NAIF_BODY_CODE += ( -61031                        )

      NAIF_BODY_NAME += ( 'JUNO_HGA'                    )
      NAIF_BODY_CODE += ( -61040                        )

      NAIF_BODY_NAME += ( 'JUNO_MGA'                    )
      NAIF_BODY_CODE += ( -61050                        )

      NAIF_BODY_NAME += ( 'JUNO_LGA_FORWARD'            )
      NAIF_BODY_CODE += ( -61061                        )

      NAIF_BODY_NAME += ( 'JUNO_LGA_AFT'                )
      NAIF_BODY_CODE += ( -61062                        )

      NAIF_BODY_NAME += ( 'JUNO_LGA_TOROID'             )
      NAIF_BODY_CODE += ( -61063                        )

      NAIF_BODY_NAME += ( 'JUNO_SRU1'                   )
      NAIF_BODY_CODE += ( -61071                        )

      NAIF_BODY_NAME += ( 'JUNO_SRU2'                   )
      NAIF_BODY_CODE += ( -61072                        )

      NAIF_BODY_NAME += ( 'JUNO_SSS1'                   )
      NAIF_BODY_CODE += ( -61073                        )

      NAIF_BODY_NAME += ( 'JUNO_SSS2'                   )
      NAIF_BODY_CODE += ( -61074                        )

      NAIF_BODY_NAME += ( 'JUNO_REM_FL1'                )
      NAIF_BODY_CODE += ( -61081                        )

      NAIF_BODY_NAME += ( 'JUNO_REM_FL2'                )
      NAIF_BODY_CODE += ( -61082                        )

      NAIF_BODY_NAME += ( 'JUNO_REM_FL3'                )
      NAIF_BODY_CODE += ( -61083                        )

      NAIF_BODY_NAME += ( 'JUNO_REM_FL4'                )
      NAIF_BODY_CODE += ( -61084                        )

      NAIF_BODY_NAME += ( 'JUNO_REM_FA1'                )
      NAIF_BODY_CODE += ( -61085                        )

      NAIF_BODY_NAME += ( 'JUNO_REM_FA2'                )
      NAIF_BODY_CODE += ( -61086                        )

      NAIF_BODY_NAME += ( 'JUNO_REM_RL1'                )
      NAIF_BODY_CODE += ( -61091                        )

      NAIF_BODY_NAME += ( 'JUNO_REM_RL2'                )
      NAIF_BODY_CODE += ( -61092                        )

      NAIF_BODY_NAME += ( 'JUNO_REM_RL3'                )
      NAIF_BODY_CODE += ( -61093                        )

      NAIF_BODY_NAME += ( 'JUNO_REM_RL4'                )
      NAIF_BODY_CODE += ( -61094                        )

      NAIF_BODY_NAME += ( 'JUNO_REM_RA1'                )
      NAIF_BODY_CODE += ( -61095                        )

      NAIF_BODY_NAME += ( 'JUNO_REM_RA2'                )
      NAIF_BODY_CODE += ( -61096                        )

   \begintext
