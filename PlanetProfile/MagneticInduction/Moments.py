import numpy as np

""" Classes for keeping track of excitation and induced magnetic moments """

class InductionStruct:
    def __init__(self):
        self.bodyname = None  # Name of body modeled.
        self.yName = None  # Name of variable along y axis. Options are "Tb", "phi", "rho", "sigma", where the first 3 are vs. salinity, and sigma is vs. thickness.
        self.Texc_hr = None  # Dict of excitation periods modeled.
        self.Amp = None  # Amplitude of dipole response (modulus of complex dipole response).
        self.phase = None  # (Positive) phase delay in degrees.
        self.Bix_nT = None  # Induced Bx dipole moments relative to body surface in nT for each excitation.
        self.Biy_nT = None  # Induced By dipole moments relative to body surface in nT for each excitation.
        self.Biz_nT = None  # Induced Bz dipole moments relative to body surface in nT for each excitation.
        self.wOcean_ppt = None  # Values of salinity used.
        self.oceanComp = None  # Ocean composition used.
        self.Tb_K = None  # Values of Bulk.Tb_K used.
        self.rhoSilMean_kgm3 = None  # Values of Sil.rhoMean_kgm3 resulted (also equal to those set for all but phi inductOtype).
        self.phiSilMax_frac = None  # Values of Sil.phiRockMax_frac set.
        self.Tmean_K = None  # Ocean mean temperature result in K.
        self.sigmaMean_Sm = None  # Mean ocean conductivity. Used to map plots vs. salinity onto D/σ plots.
        self.sigmaTop_Sm = None  # Ocean top conductivity. Used to map plots vs. salinity onto D/σ plots.
        self.D_km = None  # Ocean layer thickness in km. Used to map plots vs. salinity onto D/σ plots.
        self.zb_km = None  # Upper ice shell thickness in km.
        self.R_m = None  # Body radius in m, used to scale amplitudes.
        self.rBds_m = None  # Conducting layer upper boundaries in m.
        self.sigmaLayers_Sm = None  # Conductivities below each boundary in S/m.

        self.x = None  # Variable to plot on x axis of inductogram plots
        self.y = None  # Variable to plot on y axis of inductogram plots
        self.compsList = None  # Linear list of compositions for each model point
        self.comps = None  # Minimal list of compositions, with 1 entry per comp
        self.SINGLE_COMP = None  # Boolean flag for tracking if all of the models have the same composition

    def SetAxes(self, inductOtype):
        # Set the x and y variables to plot in inductograms based on inductOtype
        if inductOtype == 'sigma':
            self.x = self.sigmaMean_Sm
            self.y = self.D_km
        else:
            self.x = self.wOcean_ppt
            if inductOtype == 'Tb':
                self.y = self.Tb_K
            elif inductOtype == 'rho':
                self.y = self.rhoSilMean_kgm3
            elif inductOtype == 'phi':
                self.y = self.phiSilMax_frac
            else:
                raise ValueError(f'inductOtype {inductOtype} not recognized.')

    def SetComps(self, inductOtype):
        # Set some attributes pertaining to handling multiple ocean compositions in plots
        self.compsList = self.oceanComp.flatten()
        if np.all(self.compsList == self.compsList[0]) and inductOtype != 'sigma':
            self.SINGLE_COMP = True
            self.comps = [self.compsList[0]]
        else:
            self.SINGLE_COMP = False
            self.comps = np.unique(self.compsList)


class ExcitationsList:
    def __init__(self):
        self.nprmMax = 1
        Texc_hr = {}
        # Approximate (.2f) periods to select from excitation spectrum for each body in hr
        Texc_hr['Europa'] =    {'synodic': 11.23, 'orbital': 85.15, 'true anomaly': 84.63, 'synodic harmonic':  5.62}
        Texc_hr['Ganymede'] =  {'synodic': 10.53, 'orbital':171.71, 'true anomaly':  None, 'synodic harmonic':  5.27}
        Texc_hr['Enceladus'] = {'synodic':  None, 'orbital':  None, 'true anomaly': 32.93, 'synodic harmonic':  None}
        self.Texc_hr = Texc_hr

    def __call__(self, bodyname):
        if bodyname[:4] == 'Test':
            name = 'Test'
        else:
            name = bodyname
        return self.Texc_hr[name]


Excitations = ExcitationsList()
InductionResults = InductionStruct()

it = 'Europa'
Excitations.Texc_hr['Test'] = Excitations.Texc_hr[it]
