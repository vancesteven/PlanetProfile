from numpy import (ndarray, array, round, loadtxt,
                   vstack, power, gradient, arange,
                   logspace, linspace, geomspace,
                   flipud, where, ones, insert)

from toml import load as tomlload
from alma import (read_model_pp, infer_rheology_pp,
                  build_model, love_numbers)
import os

def GravityParameters(Planet, Params):
    """ Calculate induced magnetic moments for the body and prints them to disk.
    """
    SKIP = False
    if Planet.Do.VALID and Params.CALC_NEW_GRAVITY:

        alma_params = tomlload(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'paramsPyALMA3.toml'))
        # Set Magnetic struct layer arrays as we need for induction calculations
        model = SetupGravity(Planet, Params)


def SetupGravity(Planet, Params):
    """ Setup config and parameters as necessary for compatibility with PyAlma3.
        Optionally also identify asymmetric shape information from gravity.

        Requires Planet attributes:
        #TODO
    """
    # Combine data into columns to reuse existing infrastructure
    r_m = Planet.r_m[:Planet.Steps.nTotal] #TODO FIX THIS
    MODEL = vstack(
        [Planet.P_MPa, Planet.T_K, r_m, Planet.phase, Planet.rho_kgm3, Planet.Cp_JkgK, Planet.alpha_pK,
         Planet.g_ms2, Planet.phi_frac, Planet.sigma_Sm, Planet.kTherm_WmK, Planet.Seismic.VP_kms, Planet.Seismic.VS_kms, Planet.Seismic.QS,
         Planet.Seismic.KS_GPa, Planet.Seismic.GS_GPa, Planet.Ppore_MPa, Planet.rhoMatrix_kgm3,
         Planet.rhoPore_kgm3, Planet.MLayer_kg, Planet.VLayer_m3, Planet.Htidal_Wm3, Planet.eta_Pas]).T

    # List of variable names
    COLUMNS = ['P', 'T', 'r', 'phase', 'rho', 'Cp', 'alpha',
               'g', 'phi', 'sigma', 'kTherm', 'VP', 'VS',
               'QS', 'KS', 'GS', 'Ppore', 'rhoMatrix',
               'rhoPore', 'MLayer', 'VLayer', 'Htidal', 'eta']

    # List of data units
    UNITS = ['MPa', 'K', 'm', '', 'kg m-3', 'J kg-1 K-1', 'K-1',
             'm s-2', '-', 'S m-1', 'W m-1 K-1', 'km s-1', 'km s-1',
             '', 'GPa', 'GPa', 'MPa', 'kg m-3',
             'kg m-3', 'kg', 'm3', 'W m-3', 'Pa s']

    # Check GPA and MPA to Pa
    changeheader = ['P', 'VP', 'VS', 'r', 'rho', 'KS', 'GS']
    for index, (header, unit) in enumerate(zip(COLUMNS, UNITS)):
        if header in changeheader:
            if unit.lower() == 'mpa':
                UNITS[index] = 'Pa'
                MODEL[:, index] = MODEL[:, index] * 1e6

            elif unit.lower() == 'gpa':
                UNITS[index] = 'Pa'
                MODEL[:, index] = MODEL[:, index] * 1e9

            elif unit.lower() == 'km':
                UNITS[index] = 'm'
                MODEL[:, index] = MODEL[:, index] * 1e3

            elif unit.lower() == 'km s-1':
                UNITS[index] = 'm s-1'
                MODEL[:, index] = MODEL[:, index] * 1e3

    # Round radius
    MODEL[:, COLUMNS.index('r')] = round(MODEL[:, COLUMNS.index('r')], 0)

    # Flip model: core at top and surface at bottom
    if MODEL[0, COLUMNS.index('r')] > MODEL[-1, COLUMNS.index('r')]:
        MODEL = flipud(MODEL)

    ## Calculate elastic parameters
    # LAMBDA = rho (Vp^2 - 2 * Vs^2)
    # 1st Lame parameter
    LAMBDA = MODEL[:, COLUMNS.index('rho')] * (
            power(MODEL[:, COLUMNS.index('VP')], 2) - 2. * power(MODEL[:, COLUMNS.index('VS')], 2))

    # shear modulus G or MU = rho Vs^2
    # MU = MODEL[:, COLUMNS.index('rho')] * power(MODEL[:, COLUMNS.index('VS')], 2)
    MU = MODEL[:, COLUMNS.index('rho')] * power(MODEL[:, COLUMNS.index('VP')], 2)

    # Bulk Modulus K = lambda + 2/3 mu
    K = LAMBDA + 2 * MU / 3

    # Poissons ratio sigma = lambda / 2*(lambda + mu)
    SIGMA = LAMBDA / (2 * LAMBDA + 2 * MU)

    # Youngs modulus Y = 2 * MU * (1 + sigma)
    Y = 2. * MU * (1 + SIGMA)

    # Rigidity RIG = 2/3 mu
    RIG = 2 * MU / 3

    # Viscosity VIS = RIG / GRAD(Vs)
    gradvs = abs(gradient(MODEL[:, COLUMNS.index('VS')], MODEL[:, COLUMNS.index('r')]))
    gradvs[gradvs == 0.] = 1e-30

    VIS = RIG / gradvs

    # JUST FOR TESTING PP MODELS
    # SET MU and VIS to min value
    # VIS[VIS == 0.] = min(VIS[VIS>0]) / 10000
    # MU[MU == 0.] = min(MU[MU>0]) / 10000
    # VIS[VIS == 0.] = (0. * where(VIS==0.)[0] + 1.) * linspace(min(VIS[VIS>0]) - 20000, min(VIS[VIS>0]) - 150000, num=len(where(VIS==0.)[0]))
    # MU[MU == 0.] = (0. * where(MU==0.)[0] + 1.) * linspace(min(MU[MU>0]) - 20000, min(MU[MU>0]) - 150000, num=len(where(MU==0.)[0]))

    # Return output variables
    return {'columns': COLUMNS, 'units': UNITS, 'model': MODEL, 'lambda': LAMBDA, 'mu': MU, 'k': K, 'sigma': SIGMA,
            'y': Y, 'rig': RIG, 'vis': VIS}
