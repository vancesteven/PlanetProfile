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
        Planet, Params = SetupGravity(Planet, Params)


def SetupGravity(Planet, Params):
    """ Setup config and parameters as necessary for compatibility with PyAlma3.

        NOTE: Function models the read_model_pp of PyALMA3 so we ensure compatibility with it.

        Requires Planet attributes:
        #TODO
    """
    # Combine data into columns to reuse existing infrastructure
    # Note we have to use r_m[:-1] since r_m has one extra value than other arrays
    Planet.Gravity.model = vstack(
        [Planet.P_MPa, Planet.T_K, Planet.r_m[:-1], Planet.phase, Planet.rho_kgm3, Planet.Cp_JkgK, Planet.alpha_pK,
         Planet.g_ms2, Planet.phi_frac, Planet.sigma_Sm, Planet.kTherm_WmK, Planet.Seismic.VP_kms, Planet.Seismic.VS_kms, Planet.Seismic.QS,
         Planet.Seismic.KS_GPa, Planet.Seismic.GS_GPa, Planet.Ppore_MPa, Planet.rhoMatrix_kgm3,
         Planet.rhoPore_kgm3, Planet.MLayer_kg, Planet.VLayer_m3, Planet.Htidal_Wm3, Planet.eta_Pas]).T

    # Convert parameter units to Pa and meters
    for index, (header, unit) in enumerate(zip(Planet.Gravity.columns, Planet.Gravity.units_PyALMA3)):
        if header in Planet.Gravity.parameters_to_convert:
            Planet.Gravity.model[:, index] = Planet.Gravity.model[:, index] * Planet.Gravity.parameters_to_convert[header]

    # Round radius - Is similar PyALMA3 function, it rounds radius so will do so here as well
    Planet.Gravity.model[:, 2] = round(Planet.Gravity.model[:, 2], 0)

    # Flip model: core at top and surface at bottom
    Planet.Gravity.model = flipud(Planet.Gravity.model)

    ## Calculate elastic parameters
    # LAMBDA = rho (Vp^2 - 2 * Vs^2)
    # 1st Lame parameter
    Planet.Gravity.LAMBDA_Pa = Planet.Gravity.model[:, 4] * (
            power(Planet.Gravity.model[:, 11], 2) -
            2. * power(Planet.Gravity.model[:, 12], 2))

    # shear modulus G or MU = rho Vs^2
    # MU = MODEL[:, COLUMNS.index('rho')] * power(MODEL[:, COLUMNS.index('VS')], 2)
    Planet.Gravity.MU_Pa = Planet.Gravity.model[:, 4] * power(Planet.Gravity.model[:, 11], 2)

    # Bulk Modulus K = lambda + 2/3 mu
    Planet.Gravity.K_Pa = Planet.Gravity.LAMBDA_Pa + 2 * Planet.Gravity.MU_Pa / 3

    # Poissons ratio sigma = lambda / 2*(lambda + mu)
    Planet.Gravity.SIGMA = Planet.Gravity.LAMBDA_Pa / (2 * Planet.Gravity.LAMBDA_Pa + 2 * Planet.Gravity.MU_Pa)

    # Youngs modulus Y = 2 * MU * (1 + sigma)
    Planet.Gravity.Y_Pa = 2. * Planet.Gravity.MU_Pa * (1 + Planet.Gravity.SIGMA)

    # Rigidity RIG = 2/3 mu
    Planet.Gravity.RIGIDITY_Pa = 2 * Planet.Gravity.MU_Pa / 3

    # Viscosity VIS = RIG / GRAD(Vs) where Grad(Vs) uses Vs and r parameters
    Planet.Gravity.grad_Vs_s = abs(gradient(Planet.Gravity.model[:, 12], Planet.Gravity.model[:, 2]))
    Planet.Gravity.grad_Vs_s[Planet.Gravity.grad_Vs_s == 0.] = 1e-30
    Planet.Gravity.VISCOSITY_kgms = Planet.Gravity.RIGIDITY_Pa / Planet.Gravity.grad_Vs_s

    # Return Planet and Params
    return Planet, Params
