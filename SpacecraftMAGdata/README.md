# SpacecraftMAGdata

This folder should contain symlinks or locally copied folders for each spacecraft, then folders for each target body, then .tab files with the names as copied from the Planetary Data System (PDS). Like this:

SpacecraftMAGdata
    -> Galileo
        -> Jupiter
            -> ORB01_SYS3.TAB
            -> ORB02_SYS3.TAB
            -> etc.
        -> Europa
            -> ORB04_EUR_SYS3.TAB
            -> etc.
        -> etc.
    -> Cassini
        -> Saturn
            -> 04135_04240_00_FGM_KRTP_1S.TAB
            -> etc.
            
This is the file structure expected by PlanetProfile by default, and can be edited in configPP.py after installation.
