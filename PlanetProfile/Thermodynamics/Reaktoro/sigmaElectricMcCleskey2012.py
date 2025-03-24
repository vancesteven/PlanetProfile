import numpy as np
import logging
from PlanetProfile.Utilities.defineStructs import Constants

# Assign logger
log = logging.getLogger('PlanetProfile')

def elecCondMcCleskey2012(P_MPa, T_C, ions):
    # Full data structure including all ions from the MATLAB data matrix
    data = {
        'K_p1': {'z': 1, 'lam2': 0.003046, 'lam1': 1.261, 'lam0': 40.7, 'A2': 0.00535, 'A1': 0.9316, 'A0': 22.59, 'B': 1.5},
        'Na_p1': {'z': 1, 'lam2': 0.003763, 'lam1': 0.877, 'lam0': 26.23, 'A2': 0.00027, 'A1': 1.141, 'A0': 32.07, 'B': 1.7},
        'H_p1': {'z': 1, 'lam2': 0.01414, 'lam1': 5.355, 'lam0': 224.2, 'A2': 0.00918, 'A1': 1.842, 'A0': 39.23, 'B': 0.3},
        'Li_p1': {'z': 1, 'lam2': 0.002628, 'lam1': 0.7079, 'lam0': 19.2, 'A2': 0.00412, 'A1': 0.4632, 'A0': 13.71,'B': 0.2},
        'Cs_p1': {'z': 1, 'lam2': 0.003453, 'lam1': 1.249, 'lam0': 43.94, 'A2': 0.00646, 'A1': 0.7023, 'A0': 21.79, 'B': 1.3},
        'NH4_p1': {'z': 1, 'lam2': 0.003341, 'lam1': 1.285, 'lam0': 39.04, 'A2': 0.00132, 'A1': 0.607, 'A0': 11.19, 'B': 0.3},
        'Ca_p2': {'z': 2, 'lam2': 0.009645, 'lam1': 1.984, 'lam0': 62.28, 'A2': 0.03174, 'A1': 2.334, 'A0': 132.3, 'B': 2.8},
        'Mg_p2': {'z': 2, 'lam2': 0.01068, 'lam1': 1.695, 'lam0': 57.16, 'A2': 0.02453, 'A1': 1.915, 'A0': 80.5, 'B': 2.1},
        'Ba_p2': {'z': 2, 'lam2': 0.01059, 'lam1': 2.09, 'lam0': 68.1, 'A2': 0.03127, 'A1': 2.248, 'A0': 93.91, 'B': 1.9},
        'Sr_p2': {'z': 2, 'lam2': 0.006649, 'lam1': 2.069, 'lam0': 61.63, 'A2': 0.00702, 'A1': 0.9009, 'A0': 33.41, 'B': 0.1},
        'SO4_m2': {'z': -2, 'lam2': 0.01037, 'lam1': 2.838, 'lam0': 82.37, 'A2': 0.03324, 'A1': 5.889, 'A0': 193.5, 'B': 2.6},
        'Cl_m1': {'z': -1, 'lam2': 0.003817, 'lam1': 1.337, 'lam0': 40.99, 'A2': 0.00613, 'A1': 0.9469, 'A0': 22.01, 'B': 1.5},
        'F_m1': {'z': -1, 'lam2': 0.002764, 'lam1': 1.087, 'lam0': 26.66, 'A2': 0.00178, 'A1': 0.6202, 'A0': 19.34, 'B': 0.5},
        'Br_m1': {'z': -1, 'lam2': 0.000709, 'lam1': 1.477, 'lam0': 40.91, 'A2': 0.00251, 'A1': 0.5398, 'A0': 12.01, 'B': 0.1},
        'CO3_m2': {'z': -2, 'lam2': 0.00326, 'lam1': 2.998, 'lam0': 64.03, 'A2': 0.00181, 'A1': 5.542, 'A0': 120.2, 'B': 2.3},
        'HCO3_m1': {'z': -1, 'lam2': 0.000614, 'lam1': 0.9048, 'lam0': 21.14, 'A2': 0.00503, 'A1': 0.8957, 'A0': 10.97, 'B': 0.1},
        'NO3_m1': {'z': -1, 'lam2': 0.001925, 'lam1': 1.214, 'lam0': 39.9, 'A2': 0.00118, 'A1': 0.5045, 'A0': 23.31, 'B': 0.1},
        'OH_m1': {'z': -1, 'lam2': 0.003396, 'lam1': 2.925, 'lam0': 121.3, 'A2': 0.00933, 'A1': 0.1086, 'A0': 35.9, 'B': 0.01},
        'Al_p3': {'z': 3, 'lam2': 0.02376, 'lam1': 3.227, 'lam0': 90.24, 'A2': 0.06484, 'A1': 5.149, 'A0': 76.79, 'B': 3},
        'Cu_p2': {'z': 2, 'lam2': 0.00818, 'lam1': 1.939, 'lam0': 53.26, 'A2': 0.0292, 'A1': 6.745, 'A0': 151.5, 'B': 8},
        'Fe_p2': {'z': 2, 'lam2': 0.009939, 'lam1': 1.878, 'lam0': 54.8, 'A2': 0.03997, 'A1': 3.217, 'A0': 164.5, 'B': 4},
        'Fe_p3': {'z': 3, 'lam2': 0.02077, 'lam1': 4.39, 'lam0': 82.42, 'A2': 0.09676, 'A1': 20.76, 'A0': 22.18, 'B': 4},
        'Mn_p2': {'z': 2, 'lam2': 0.01275, 'lam1': 2.109, 'lam0': 46.19, 'A2': 0.1071, 'A1': 9.023, 'A0': 135.4, 'B': 7.6},
        'Zn_p2': {'z': 2, 'lam2': 0.01249, 'lam1': 1.912, 'lam0': 48.2, 'A2': 0.08284, 'A1': 5.188, 'A0': 75.73, 'B': 7},
        'KSO4_m1': {'z': -1, 'lam2': 0.002439, 'lam1': 4.253, 'lam0': 129.7, 'A2': 0.01576, 'A1': 6.21, 'A0': 146.8, 'B': 1.3},
        'NaSO4_m1': {'z': -1, 'lam2': 0.002309, 'lam1': 5.459, 'lam0': 219.2, 'A2': 0.01454, 'A1': 5.193, 'A0': 253.6, 'B': 0.5},
        'HSO4_m1': {'z': -1, 'lam2': 0.000927, 'lam1': 0.8337, 'lam0': 29.56, 'A2': 0.02887, 'A1': 0.873, 'A0': 36.25, 'B': 7},
        'NaCO3_m1': {'z': -1, 'lam2': 0.00336, 'lam1': 3.845, 'lam0': 89.51, 'A2': 0.00061, 'A1': 6.387, 'A0': 141.7, 'B': 2}
    }
    # Ensure there is at least one ion in the ions list, otherwise we will just set the conductivity to pure water
    # Make deep copy of ions
    if ions:
        sigma_mS_cm = np.zeros_like(T_C)
    else:
        log.error("No ion data available. Assuming pure water")
        return np.zeros_like(T_C) + Constants.sigmaH2O_Sm
    for ion_name, ion_data in ions.items():
        if ion_name in data:
            # Get relevant data
            mols = np.array(ion_data['mols'])
            z = data[ion_name]['z']
            B = data[ion_name]['B']
            Ad = np.array([data[ion_name]['A2'], data[ion_name]['A1'], data[ion_name]['A0']])
            l210 = np.array([data[ion_name]['lam2'], data[ion_name]['lam1'], data[ion_name]['lam0']])

            # Ionic strength calculations
            I = 0.5 * mols * z**2
            I2 = np.sqrt(I)

            # Evaluate polynomial coefficients for temperature dependency
            lam0 = np.polyval(l210, T_C)
            At = np.polyval(Ad, T_C)

            # Evaluate (I^1/2)/(1+B*I1/2)
            xI = I2/(1+B*I2)

            # Evaluate molal conductivity
            lamda = lam0 - At * xI  # Use broadcasting

            # IMPLEMENT PRESSURE ADJUSTMENT HERE

            # Update sigma_Sm with conductivity contributions from each ion
            sigma_mS_cm += lamda * mols  # Broadcasting mols across num_TC

    return (sigma_mS_cm * 100 / 1000)  # convert from mS/cm to S/m

