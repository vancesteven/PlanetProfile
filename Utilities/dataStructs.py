"""
dataStructs: Create empty classes and subclasses for holding body-specific data

Example usage:
Planet = PlanetStruct("nameOfBody")
Planet.R_m = 1560e3
Planet.Ocean.comp = "MgSO4"
Planet.Silicate.mantleEOS = "CV3hy1wt_678_1.tab"
Planet.Core.Fe_core = False
"""
class PlanetStruct:
    def __init__(self, name):
        self.name = name

    class Ocean:
        def fnTfreeze_K(self, PPg, wwg, TT):
            # Somehow make an interpolator a la:
            # load L_Ice_MgSO4.mat
            # Planet.Ocean.fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');
            pass

    class Silicate: # Quantities from the old "Seismic" struct should go in here.
        pass

    class Core:
        pass

    class Seismic:
        pass

    class Magnetic:
        pass

class ParamsStruct:

    class lbls:
        pass