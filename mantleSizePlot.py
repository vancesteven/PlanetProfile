import matplotlib.pyplot as plt
import numpy as np
#from Utilities.dataStructs import *

""" function MantleSizePlot()
Create a plot of density vs. radius in the mantle.
Returns:
    None
Inputs:
    Planet: PlanetStruct. Contains body configuration info from PPBody.py file. 
    rho_sil_kgm3 : 2d float numpy array of size (n,m,)
        densities of mantle [kg/m^3]
    R_sil_m : 2d float numpy array of size (n,m,)
        radii of silicate layer [m]
"""
def mantleSizePlot(Planet, rho_sil_kgm3, R_sil_m, figbase):

    Tb_K = [thisPlanet["Tb_K"] for thisPlanet in Planet]
    Cmeasured = [thisPlanet["Cmeasured"] for thisPlanet in Planet]
    Cuncertainty = [thisPlanet["Cuncertainty"] for thisPlanet in Planet]
    wo = [thisPlanet["ocean_wpct"] for thisPlanet in Planet]

    lstr_3 = []
    for iT in range(nTbs):
        plt.plot(rho_sil_kgm3 , R_sil_m*1e-3 , linewidth = Params.lw)
        lstr_3.append( f"$T_{{b}}$: {Tb_K[iT]:0.1f} K" )

    plt.legend(lstr_3)

    plt.xlabel("$\\rho_{\\mathrm{sil}} \\, (\\mathrm{kg} \\, \\mathrm{m}^{-3}$)")
    plt.ylabel("$R_{\\mathrm{sil}} \\, (\\mathrm{km})$")
    plt.title(f"No Fe core ; $C/MR^2 = {Cmeasured} \\pm {Cuncertainty}$ ; $ W = {wo} \\, wt \\%$ ")

    plt.savefig(fpath+cfg.xtn)
    print("Printed mantle size plot to file: " + fpath+cfg.xtn)

    if show:
        plt.ion() # turns on interactive mode, allowing execution to continue while plot is shown
        plt.show()