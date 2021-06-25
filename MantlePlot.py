import matplotlib.pyplot as plt
import numpy as np

def MantlePlot( rho_sil_kgm3 , R_sil_m , C2inds , Planet:dict , nTbs , wo , lw = 1):
    Tb_K = Planet["Tb_K"]
    Cmeasured = Planet["Cmeasured"]
    Cuncertainty = Planet["Cuncertainty"]

    # technically, C2inds should have values subtracted by 1
    # due to difference in MATLAB and python indexing
    lstr_3 = []
    for iT in range(0,nTbs):
        plt.plot(rho_sil_kgm3[iT][C2inds[iT]] , R_sil_m[iT][C2inds[iT]]*1e-3 , linewidth = lw)
        lstr_3.append( f"$T_{{b}}$: {Tb_K[iT]:0.1f} K" )

    plt.legend(lstr_3)

    plt.xlabel("$\\rho_{\\mathrm{sil}} \\, (\\mathrm{kg} \\, \\mathrm{m}^{-3}$)")
    plt.ylabel("$R_{\\mathrm{sil}} \\, (\\mathrm{km})$")
    plt.title(f"No Fe core ; $C/MR^2 = {Cmeasured} \\pm {Cuncertainty}$ ; $ W = {wo} \\, wt \\%$ ")

    plt.show()