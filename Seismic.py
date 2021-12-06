import numpy as np

def SeismicCalcs(Planet, Params):
    """ Calculation of seismic properties, including wave speeds

        Assigns Planet attributes:
            Seismic.VP_kms, Seismic.VS_kms, Seismic.QS, Seismic.KS_GPa, Seismic.GS_GPa
    """
    VP_kms, VS_kms, QS, KS_GPa, GS_GPa = (np.zeros(Planet.Steps.nTotal) for _ in range(5))

    Planet.Seismic.VP_kms = VP_kms
    Planet.Seismic.VS_kms = VS_kms
    Planet.Seismic.QS = QS
    Planet.Seismic.KS_GPa = KS_GPa
    Planet.Seismic.GS_GPa = GS_GPa

    return Planet