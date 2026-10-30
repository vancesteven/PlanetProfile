# PlanetProfile SPICE support

PlanetProfile supports SPICE through the python module SpiceyPy. In this folder should be symlinks to (or files of) generic kernels for the outer planets (.bsp), a leapseconds kernel (.tls), and a planetary constants kernel (.tpc). Folders containing kernels with spacecraft information should also be placed here. Like this:

SPICE
    -> jup365.bsp
    -> naif0012.tls
    -> pck00010.tpc
    -> Galileo
        -> s980326a.bsp
        -> s000131a.bsp
        -> etc.
    -> Juno
        -> juno_rec_210513_210630_210707.bsp
        -> etc.
        
This is the directory structure expected by PlanetProfile by default, and can be changed in configPP.py after PlanetProfile has been installed. All .bsp kernels within the SPICE/Spacecraft directories will be loaded.
    
