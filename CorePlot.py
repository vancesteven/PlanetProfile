import matplotlib.pyplot as plt
import numpy as np

def CorePlot( R_Fe_m , R_sil_m , C2inds , Planet:dict , rho_Fe_kgm3 , nTbs , wo , saveStr , lw = 1 , show = True ):
    """
        Shows and saves a plot of radius of silicon to iron
        as implemented in PlanetProfile.m line ~1100

        Parameters:
        -----------
        R_Fe_m : 2d float numpy array of size (n,m,)
            radii of iron core [m]
        R_sil_m : 2d float numpy array of size (n,m,)
            radii of silicate layer [m]
        C2inds : 2d int numpy array of size (n,p,)
            chooses indices for display in plot (determined by constraints on C/MR2)
        Planet : dictionary with keys Tb_K,Cmeasured,Cuncertainty,rho_sil_withcore_kgm3
        rho_Fe_kgm3 : float
            density of iron core [kg/m^3]
        nTbs : int
            number of temperature profiles = n (length of Tb_K)
        wo : float
            percent concentration of solute in ocean
        saveStr : string
            save location of plot file
        lw : float (optional)
            width of lines in plot
        show : boolean (optional)
            determines if plot should be shown upon execution
    """

    Tb_K = Planet["Tb_K"]
    Cmeasured = Planet["Cmeasured"]
    Cuncertainty = Planet["Cuncertainty"]
    rho_sil_withcore_kgm3 = Planet["rho_sil_withcore_kgm3"]

    lstr_3 = []
    for iT in range(0,nTbs):
        plt.plot(R_Fe_m[iT][C2inds[iT]]*1e-3 , R_sil_m[iT][C2inds[iT]]*1e-3 , linewidth = lw)
        lstr_3.append( f"$T_{{b}}$: {Tb_K[iT]:0.1f} K" )

    plt.legend(lstr_3)

    plt.xlabel("$R_{\\mathrm{Fe}} \\, (\\mathrm{km})$")
    plt.ylabel("$R_{\\mathrm{Si}} \\, (\\mathrm{km})$")
    plt.title(f"Fe core ; $C/MR^2 = {Cmeasured} \\pm {Cuncertainty}$ ; $W = {wo} \\,$wt%; $\\rho_{{\\mathrm{{sil}}}}$ : {rho_sil_withcore_kgm3:0.0f}; $\\rho_{{\\mathrm{{Fe}}}}$ : {rho_Fe_kgm3}")

    plt.savefig(saveStr)

    if show:
        plt.show()