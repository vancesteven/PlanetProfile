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
        self.phiRockMax_frac = None  # Values of Sil.phiRockMax_frac set.
        self.Tmean_K = None  # Ocean mean temperature result in K.
        self.sigmaMean_Sm = None  # Mean ocean conductivity. Used to map plots vs. salinity onto D/sigma plots.
        self.sigmaTop_Sm = None  # Ocean top conductivity. Used to map plots vs. salinity onto D/sigma plots.
        self.D_km = None  # Ocean layer thickness in km. Used to map plots vs. salinity onto D/sigma plots.
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
                self.y = self.phiRockMax_frac
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
        Texc_hr['Io'] = {
            'synodic 4th': 3.24,
            'synodic 3rd': 4.32,
            'synodic 2nd': 6.48,
            'synodic 2nd-TA beat': 7.65,
            'synodic+TA beat': 9.92,
            'synodic': 12.95,
            'synodic-TA beat': 18.67,
            'true anomaly': 42.31,
            'orbital+year beat': 42.44,
        }
        Texc_hr['Europa'] = {
            'synodic 3rd': 3.74,
            'synodic 2nd': 5.62,
            'synodic 2nd-TA beat': 6.02,
            'synodic+TA beat': 9.92,
            'synodic': 11.23,
            'synodic-TA beat': 12.95,
            'TA+year beat': 84.57,
            'true anomaly': 84.64,
            'TA-year beat': 84.71,
            'orbital': 85.21,
            'orbital-year beat': 85.28
        }
        Texc_hr['Ganymede'] = {
            'synodic 2nd': 5.27,
            'synodic': 10.53,
            'orbital': 171.70
        }
        Texc_hr['Callisto'] = {
            'synodic 5th': 2.04,
            'synodic 3rd': 3.39,
            'synodic 2nd': 5.09,
            'synodic': 10.18,
            'orbital': 400.54
        }
        Texc_hr['Mimas'] = {
            'orbital 2nd': 11.31,
            'orbital 4th-TA 2nd-half year beat': 11.34,
            'TA 2nd-orbital beat': 22.50,
            'true anomaly': 22.56,
            'orbital': 22.62,
            'orbital 2nd-TA beat': 22.677,
            'orbital 2nd-TA-year beat': 22.679,
            'orbital 2nd-TA-half year beat': 22.681
        }
        Texc_hr['Enceladus'] = {'true anomaly': 32.93}
        Texc_hr['Tethys'] = {'true anomaly': 45.26, 'TA-year beat': 45.27}
        Texc_hr['Dione'] = {'orbital-half year beat': 65.72}
        Texc_hr['Rhea'] = {'orbital': 108.42}
        Texc_hr['Titan'] = {'orbital': 382.69}
        Texc_hr['Iapetus'] = {'orbital': 1903.94}
        Texc_hr['Miranda'] = {
            'synodic 4th': 8.76,
            'synodic 2nd+orbital beat': 11.56,
            'synodic 3rd': 11.69,
            'synodic+orbital beat': 17.24,
            'synodic 2nd': 17.53,
            'orbital': 33.92,
            'synodic': 35.06,
            'orbital-synodic beat': 1050.93
        }
        Texc_hr['Ariel'] = {
            'synodic 3rd': 8.04,
            'synodic 2nd': 12.06,
            'synodic': 24.11
        }
        Texc_hr['Umbriel'] = {
            'synodic 2nd': 10.43,
            'synodic': 20.85,
            'synodic-orbital beat': 26.39
        }
        Texc_hr['Titania'] = {'synodic 2nd': 9.40, 'synodic': 18.79}
        Texc_hr['Oberon'] = {'synodic': 18.21}
        Texc_hr['Triton'] = {
            'synodic': 14.46,
            'synodic-orbital beat': 16.11,
            'orbital': 141.04
        }
        self.Texc_hr = Texc_hr

    def __call__(self, bodyname):
        if bodyname[:4] == 'Test':
            name = 'Test'
        else:
            name = bodyname
        if name in self.Texc_hr.keys():
            return self.Texc_hr[name]
        else:
            return None


Excitations = ExcitationsList()
InductionResults = InductionStruct()

it = 'Europa'
Excitations.Texc_hr['Test'] = Excitations.Texc_hr[it]
