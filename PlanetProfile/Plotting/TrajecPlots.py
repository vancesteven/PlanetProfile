import os
import numpy as np
import logging
import matplotlib.pyplot as plt
import spiceypy as spice
from matplotlib.gridspec import GridSpec
from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc
from PlanetProfile.Utilities.defineStructs import Constants
from PlanetProfile.MagneticInduction.Moments import Excitations
from PlanetProfile.TrajecAnalysis.FlybyEvents import FlybyCA
from PlanetProfile.TrajecAnalysis.SpiceFuncs import BiTrajec, parentGM, spiceCode
from MoonMag.asymmetry_funcs import getMagSurf as GetMagSurf

# Assign logger
log = logging.getLogger('PlanetProfile')

def PlotMagCA(PlanetList, Params, scName):
    for Planet in PlanetList:
        Bx, By, Bz = BiTrajec(Planet, Params, scName, np.fromiter(FlybyCA[scName].etCA[Planet.bodyname].values(), dtype=np.float_))
        # Sum over all excitation periods
        Bx = np.sum(Bx, axis=0)
        By = np.sum(By, axis=0)
        Bz = np.sum(Bz, axis=0)
        Bmag_nT = np.sqrt(Bx**2 + By**2 + Bz**2)
        hCA_km = np.fromiter(FlybyCA[scName].rCA_km[Planet.bodyname].values(), dtype=np.float_) - Planet.Bulk.R_m / 1e3
        CAlbls = list(FlybyCA[scName].rCA_km[Planet.bodyname].keys())

        fig = plt.figure(figsize=FigSize.MagCA)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        if Style.GRIDS:
            ax.grid()
            ax.set_axisbelow(True)

        ax.set_title(f'{Planet.bodyname}: {FigLbl.MagCAtitle}')
        ax.set_xlabel(FigLbl.rCA)
        ax.set_ylabel(FigLbl.BCA)

        if FigMisc.SHOW_MAG_THRESH:
            ax.hlines(xmin=0, xmax=FigMisc.hCAmax_km, y=FigMisc.thresh_nT, linestyle=Style.LS_thresh, linewidth=Style.LW_thresh,
                                           color=Color.thresh)
            ax.annotate(FigLbl.thresh, (FigMisc.threshCenter, FigMisc.thresh_nT),
                        fontsize=FigMisc.CAlblSize)
        ax.scatter(hCA_km, Bmag_nT, s=Style.MW_CA, marker=Style.MS_CA, c=Color.CAdot)
        for CAlbl, hCA, Bmag in zip(CAlbls, hCA_km, Bmag_nT):
            ax.annotate(CAlbl, (hCA, Bmag), fontsize=FigMisc.CAlblSize)

        ax.set_xlim([0, FigMisc.hCAmax_km])
        ax.set_ylim(bottom=0)
        plt.tight_layout()
        fig.savefig(Params.FigureFiles.MagCA, format=FigMisc.figFormat, dpi=FigMisc.dpi)
        log.debug(f'Plot saved to file: {Params.FigureFiles.MagCA}')
        plt.close()

    return


def PlotApsidalPrec(bodynames, Params, tStart_UTC, tEnd_UTC, tRes_d=1):

    fig = plt.figure(figsize=(6, 6))
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(FigLbl.tPastJ2000)
    ax.set_ylabel(FigLbl.argPeri)
    ax.set_title(f'{FigLbl.apsidalTitle} from ' +
                 f'{tStart_UTC.replace(" ", "T").split("T")[0]} to ' +
                 f'{tEnd_UTC.replace(" ", "T").split("T")[0]}')

    for body in np.unique(bodynames):
        _, _, parent = spiceCode(body)
        GM = parentGM[parent]
        ets = np.arange(spice.str2et(tStart_UTC), spice.str2et(tEnd_UTC), tRes_d*86400)
        states, _ = spice.spkezr(body.upper(), ets, 'J2000', 'NONE', parent.upper())
        argPeri_deg = np.degrees(np.array([spice.oscltx(state, et, GM)[4] for state, et in zip(states, ets)]))
        ax.plot(ets * FigLbl.tJ2000mult, argPeri_deg, label=body)

    if Params.LEGEND:
        ax.legend()
    plt.tight_layout()
    fig.savefig(Params.FigureFiles.apsidal, format=FigMisc.figFormat, dpi=FigMisc.dpi)
    log.debug(f'Plot saved to file: {Params.FigureFiles.apsidal}')
    plt.close()

    return
