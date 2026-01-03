import numpy as np

""" Classes for keeping track of excitation moments """
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

it = 'Europa'
Excitations.Texc_hr['Test'] = Excitations.Texc_hr[it]
