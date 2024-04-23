import os
import numpy as np
import logging
import matplotlib.pyplot as plt
import spiceypy as spice
import matplotlib.dates as mdt
from matplotlib.gridspec import GridSpec
from matplotlib import patches
from PlanetProfile import _healpixSphere
from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc
from PlanetProfile.TrajecAnalysis.MagneticFields import BiTrajecSingle
from PlanetProfile.TrajecAnalysis.FlybyEvents import GetFlybyCA
from PlanetProfile.TrajecAnalysis.SpiceFuncs import parentGM, spiceCode, spiceSCname, BodyVel_kms
from PlanetProfile.Utilities.defineStructs import xyzComps, Constants

# Assign logger
log = logging.getLogger('PlanetProfile')

#plt.rcParams['axes.prop_cycle'] = PPcycler
def PlotFlybys(Params, magData, modelData):

    FlybyCA = GetFlybyCA()
    FigLbl.SetTrajec(Params.Trajec.targetBody, modelData.allFlybys)
    if Params.PLOT_BINVERSION:
        PlotModelBxyz(Params, magData, modelData, FlybyCA)
    if Params.PLOT_TRAJECS:
        PlotTrajecs(Params, modelData, FlybyCA)

    return


def PlotModelBxyz(Params, magData, modelData, FlybyCA):
    """ Plot a comparison of measurements and inverted model for magnetic fields of selected flybys.
    """

    for scName, data in magData.items():
        for fbID, thisFname in Params.Trajec.FigureFiles.Btrajec[scName].items():

            # Create figure
            fig = plt.figure(figsize=FigSize.BtrajecCombo)
            grid = GridSpec(3, 1)
            axes = np.array([fig.add_subplot(grid[i, 0]) for i in range(3)])
            if Style.GRIDS:
                [ax.grid() for ax in axes]
                [ax.set_axisbelow(True) for ax in axes]

            # Labels and titles
            tCA = FlybyCA[scName].etCA[Params.Trajec.targetBody][fbID]
            if FigLbl.tCA_RELATIVE:
                t = (data.ets[fbID] - tCA) * FigLbl.tCArelMult
                tCA = 0
                axes[2].set_xlabel(FigLbl.tCArelLabel)
            else:
                t = np.array(data.t_UTC[fbID], dtype='datetime64')
                datefmt = mdt.ConciseDateFormatter(mdt.AutoDateLocator())
                axes[2].xaxis.set_major_formatter(datefmt)

            [ax.set_xlim([np.min(t), np.max(t)]) for ax in axes]
            [ax.set_ylabel(FigLbl.yLabelsBtrajec[comp]) for ax, comp in zip(axes, xyzComps)]
            if Params.TITLES:
                fig.suptitle(FigLbl.BtrajecTitle[scName][fbID])

            # Plot the data
            MAGdata = [data.BxIAU_nT[fbID], data.ByIAU_nT[fbID], data.BzIAU_nT[fbID]]
            [ax.plot(t, MAG, label=FigLbl.MAGdataLabel, color=Color.MAGdata,
                     linestyle=Style.LS_MAGdata, linewidth=Style.LW_MAGdata)
                for ax, MAG in zip(axes, MAGdata)]
            modelNet = [modelData.BxIAU_nT[scName][fbID], modelData.ByIAU_nT[scName][fbID],
                        modelData.BzIAU_nT[scName][fbID]]
            [ax.plot(t, net, label=FigLbl.modelNetLabel, color=Color.BcompsModelNet,
                     linestyle=Style.LS_modelNet, linewidth=Style.LW_modelNet)
                for ax, net in zip(axes, modelNet)]

            if FigMisc.SHOW_EXCITATION:
                modelExc = [modelData.BxIAUexc_nT[scName][fbID], modelData.ByIAUexc_nT[scName][fbID],
                            modelData.BzIAUexc_nT[scName][fbID]]
                [ax.plot(t, exc, label=FigLbl.modelExcLabel, color=Color.BcompsModelExc,
                         linestyle=Style.LS_modelExc, linewidth=Style.LW_modelExc)
                 for ax, exc in zip(axes, modelExc)]

            if FigMisc.SHOW_INDUCED:
                modelInd = [
                    modelData.BxIAUexc_nT[scName][fbID] + modelData.BxIAUind_nT[scName][fbID],
                    modelData.ByIAUexc_nT[scName][fbID] + modelData.ByIAUind_nT[scName][fbID],
                    modelData.BzIAUexc_nT[scName][fbID] + modelData.BzIAUind_nT[scName][fbID]]
                [ax.plot(t, ind, label=FigLbl.modelIndLabel, color=Color.BcompsModelInd,
                         linestyle=Style.LS_modelInd, linewidth=Style.LW_modelInd)
                 for ax, ind in zip(axes, modelInd)]

            if FigMisc.SHOW_PLASMA:
                modelPls = [
                    modelData.BxIAUexc_nT[scName][fbID] + modelData.BxIAUpls_nT[scName][fbID],
                    modelData.ByIAUexc_nT[scName][fbID] + modelData.ByIAUpls_nT[scName][fbID],
                    modelData.BzIAUexc_nT[scName][fbID] + modelData.BzIAUpls_nT[scName][fbID]]
                [ax.plot(t, pls, label=FigLbl.modelPlsLabel, color=Color.BcompsModelPls,
                         linestyle=Style.LS_modelPls, linewidth=Style.LW_modelPls)
                 for ax, pls in zip(axes, modelPls)]

            if FigMisc.MARK_CA_B:
                for ax in axes:
                    ax.axvline(x=tCA, color=Color.CAline, linestyle=Style.LS_CA,
                               linewidth=Style.LW_CA)
                ymin, ymax = axes[0].get_ylim()
                axes[0].text(tCA + FigLbl.CAoffset[0], ymax + (ymax - ymin) * FigLbl.CAoffset[1],
                        FigLbl.CAtxt, ha='center')

            if Params.LEGEND:
                axes[2].legend()

            plt.tight_layout()
            fig.savefig(thisFname, format=FigMisc.figFormat, dpi=FigMisc.dpi, metadata=FigLbl.meta)
            log.debug(f'Magnetic field trajectory comparison plot saved to file: {thisFname}')
            plt.close()

    return


def PlotTrajecs(Params, modelData, FlybyCA):
    """ Plot trajectories of selected flybys in the frame of the body in two planar projections.
    """

    # Create figures
    fig = plt.figure(figsize=FigSize.SCtrajecCombo)
    fig3D = plt.figure(figsize=FigSize.SCtrajec3D)
    grid = GridSpec(1, 3)
    grid3D = GridSpec(1, 1)
    axes = np.array([fig.add_subplot(grid[0, i]) for i in range(3)])
    ax3D = fig3D.add_subplot(grid3D[0, 0], projection='3d')
    if Style.GRIDS:
        [ax.grid() for ax in axes]
        [ax.set_axisbelow(True) for ax in axes]

    # Labels and titles
    axes[0].set_xlabel(FigLbl.xLabelsTrajec)
    axes[0].set_ylabel(FigLbl.yLabelsTrajec)
    axes[1].set_xlabel(FigLbl.xLabelsTrajec)
    axes[1].set_ylabel(FigLbl.zLabelsTrajec)
    axes[2].set_xlabel(FigLbl.yLabelsTrajec)
    axes[2].set_ylabel(FigLbl.zLabelsTrajec)
    ax3D.set_xlabel(FigLbl.xLabelsTrajec)
    ax3D.set_ylabel(FigLbl.yLabelsTrajec)
    if FigLbl.AXES_INFO:
        ax3D.zaxis.set_rotate_label(False)
    ax3D.set_zlabel(FigLbl.zLabelsTrajec, rotation=90)

    if Params.TITLES:
        fig.suptitle(FigLbl.trajecTitle)
        fig3D.suptitle(FigLbl.trajecTitle)
        axes[0].set_title('$xy$ plane')
        axes[1].set_title('$xz$ plane')
        axes[2].set_title('$yz$ plane')

    # Initialize axis max size tracker
    axLim = 0
    [ax.set_aspect(1) for ax in axes]
    ax3D.set_aspect('equal')

    # Show body surface
    [ax.add_artist(patches.Circle((0, 0), 1.0, color=Color.bodySurface, ec='none')) for ax in axes]
    # Load 1D lists of healpix surface pixels
    thHealpix, phiHealpix = np.loadtxt(_healpixSphere, delimiter=',', skiprows=2, unpack=True)
    # Interpolate healpix locations for compatibility with matplotlib 3D plotting functions
    th = np.unique(np.concatenate((np.array([0]), thHealpix, np.array([np.pi]))))
    phi = np.unique(np.concatenate((np.array([0]), phiHealpix, np.array([2*np.pi]))))
    thSurf, phiSurf = np.meshgrid(th, phi, indexing='ij')
    xSurf = np.sin(thSurf) * np.cos(phiSurf)
    ySurf = np.sin(thSurf) * np.sin(phiSurf)
    zSurf = np.cos(thSurf)
    ax3D.plot_surface(xSurf, ySurf, zSurf, color=Color.bodySurface)

    # Plot the trajectories
    for scName, fbList in modelData.allFlybys.items():
        for fbID, fbName in fbList.items():

            # Copy to locals for convenience
            x = modelData.x_Rp[scName][fbID]
            y = modelData.y_Rp[scName][fbID]
            z = modelData.z_Rp[scName][fbID]
            if FigMisc.trajLims is not None:
                inBounds = np.logical_and(modelData.r_Rp[scName][fbID] > -FigMisc.trajLims,
                                          modelData.r_Rp[scName][fbID] < FigMisc.trajLims)
                x = x[inBounds]
                y = y[inBounds]
                z = z[inBounds]
            else:
                inBounds = None

            # Skip this flyby if it doesn't fit our plotting criteria
            if np.size(x) == 0:
                continue
            axLim = np.maximum(axLim, np.max([np.abs(x), np.abs(y), np.abs(z)]))

            # Get CA xyz
            tCA = FlybyCA[scName].etCA[Params.Trajec.targetBody][fbID]
            iCA = np.argmin(np.abs(modelData.ets[scName][fbID][inBounds] - tCA))
            xCA = x[iCA]
            yCA = y[iCA]
            zCA = z[iCA]

            xSign = np.sign(x)
            ySign = np.sign(y)
            zSign = np.sign(z)
            # Protect from edge cases where a coordinate is exactly 0
            xSign[xSign==0] = 1
            ySign[ySign==0] = 1
            zSign[zSign==0] = 1

            if np.all(zSign == zSign[0]):
                xy = axes[0].plot(x, y, label=FigLbl.FBlabel[scName][fbID], zorder=np.mean(z),
                    linestyle=Style.LS_SCtrajec[scName], linewidth=Style.LW_SCtrajec)
                mcolor = xy[0].get_color()
            else:
                iTrans = np.where([zsm != zsp for zsm, zsp in zip(zSign[:-1], zSign[1:])])[0]
                xy = axes[0].plot(x[:iTrans[0]], y[:iTrans[0]], label=FigLbl.FBlabel[scName][fbID],
                    zorder=np.mean(z[:iTrans[0]]), linestyle=Style.LS_SCtrajec[scName],
                    linewidth=Style.LW_SCtrajec)
                mcolor = xy[0].get_color()

                # The sign can only change once or twice for flyby trajectories
                if np.size(iTrans) > 1:
                    axes[0].plot(x[iTrans[0]:iTrans[1]], y[iTrans[0]:iTrans[1]],
                                 zorder=np.mean(z[iTrans[0]:iTrans[1]]),
                                 linestyle=Style.LS_SCtrajec[scName], linewidth=Style.LW_SCtrajec,
                                 color=mcolor)
                    axes[0].plot(x[iTrans[1]:], y[iTrans[1]:], zorder=np.mean(z[iTrans[1]:]),
                                 linestyle=Style.LS_SCtrajec[scName], linewidth=Style.LW_SCtrajec,
                                 color=mcolor)
                else:
                    axes[0].plot(x[iTrans[0]:], y[iTrans[0]:], zorder=np.mean(z[iTrans[0]:]),
                                 linestyle=Style.LS_SCtrajec[scName], linewidth=Style.LW_SCtrajec,
                                 color=mcolor)

            if np.all(ySign == ySign[0]):
                axes[1].plot(x, z, label=FigLbl.FBlabel[scName][fbID], zorder=-np.mean(y),
                    linestyle=Style.LS_SCtrajec[scName], linewidth=Style.LW_SCtrajec,
                    color=mcolor)
            else:
                iTrans = np.where([ysm != ysp for ysm, ysp in zip(ySign[:-1], ySign[1:])])[0]
                axes[1].plot(x[:iTrans[0]], z[:iTrans[0]], label=FigLbl.FBlabel[scName][fbID],
                    zorder=-np.mean(y[:iTrans[0]]), linestyle=Style.LS_SCtrajec[scName],
                    linewidth=Style.LW_SCtrajec, color=mcolor)
                # The sign can only change once or twice for flyby trajectories
                if np.size(iTrans) > 1:
                    axes[1].plot(x[iTrans[0]:iTrans[1]], z[iTrans[0]:iTrans[1]],
                                 zorder=-np.mean(y[iTrans[0]:iTrans[1]]),
                                 linestyle=Style.LS_SCtrajec[scName], linewidth=Style.LW_SCtrajec,
                                 color=mcolor)
                    axes[1].plot(x[iTrans[1]:], z[iTrans[1]:], zorder=-np.mean(y[iTrans[1]:]),
                                 linestyle=Style.LS_SCtrajec[scName], linewidth=Style.LW_SCtrajec,
                                 color=mcolor)
                else:
                    axes[1].plot(x[iTrans[0]:], z[iTrans[0]:], zorder=-np.mean(y[iTrans[0]:]),
                                 linestyle=Style.LS_SCtrajec[scName], linewidth=Style.LW_SCtrajec,
                                 color=mcolor)

            if np.all(xSign == xSign[0]):
                axes[2].plot(y, z, label=FigLbl.FBlabel[scName][fbID], zorder=np.mean(x),
                    linestyle=Style.LS_SCtrajec[scName], linewidth=Style.LW_SCtrajec,
                    color=mcolor)
            else:
                iTrans = np.where([xsm != xsp for xsm, xsp in zip(xSign[:-1], xSign[1:])])[0]
                axes[2].plot(y[:iTrans[0]], z[:iTrans[0]], label=FigLbl.FBlabel[scName][fbID],
                    zorder=np.mean(x[:iTrans[0]]), linestyle=Style.LS_SCtrajec[scName],
                    linewidth=Style.LW_SCtrajec, color=mcolor)
                # The sign can only change once or twice for flyby trajectories
                if np.size(iTrans) > 1:
                    axes[2].plot(y[iTrans[0]:iTrans[1]], z[iTrans[0]:iTrans[1]],
                                 zorder=np.mean(x[iTrans[0]:iTrans[1]]),
                                 linestyle=Style.LS_SCtrajec[scName], linewidth=Style.LW_SCtrajec,
                                 color=mcolor)
                    axes[2].plot(y[iTrans[1]:], z[iTrans[1]:], zorder=np.mean(x[iTrans[1]:]),
                                 linestyle=Style.LS_SCtrajec[scName], linewidth=Style.LW_SCtrajec,
                                 color=mcolor)
                else:
                    axes[2].plot(y[iTrans[0]:], z[iTrans[0]:], zorder=np.mean(x[iTrans[0]:]),
                                 linestyle=Style.LS_SCtrajec[scName], linewidth=Style.LW_SCtrajec,
                                 color=mcolor)

            ax3D.plot(x, y, z, label=FigLbl.FBlabel[scName][fbID],
                linestyle=Style.LS_SCtrajec[scName], linewidth=Style.LW_SCtrajec, color=mcolor)

            # Mark exit points
            xExit, yExit, zExit = x[-1], y[-1], z[-1]
            if FigMisc.EXIT_ARROWS:
                vx, vy, vz, _ = BodyVel_kms(spiceSCname[scName], Params.Trajec.targetBody,
                                            modelData.ets[scName][fbID][inBounds][-1])
                axes[0].quiver(xExit, yExit, vx, vy, zorder=zExit, clip_on=False,
                               angles='xy', units='inches',
                               width=Style.LW_SCtrajec/72, color=mcolor)
                axes[1].quiver(xExit, zExit, vx, vz, zorder=-yExit, clip_on=False,
                               angles='xy', units='inches',
                               width=Style.LW_SCtrajec/72, color=mcolor)
                axes[2].quiver(yExit, zExit, vy, vz, zorder=xExit, clip_on=False,
                               angles='xy', units='inches',
                               width=Style.LW_SCtrajec/72, color=mcolor)
                ax3D.quiver(xExit, yExit, zExit, vx, vy, vz, length=1, normalize=True,
                            linewidth=Style.LW_SCtrajec, color=mcolor)
            else:
                axes[0].scatter(xExit, yExit, zorder=zExit, marker=Style.MS_exit, s=Style.MW_exit,
                                edgecolors=mcolor, c='none')
                axes[1].scatter(xExit, zExit, zorder=-yExit, marker=Style.MS_exit, s=Style.MW_exit,
                                edgecolors=mcolor, c='none')
                axes[2].scatter(yExit, zExit, zorder=xExit, marker=Style.MS_exit, s=Style.MW_exit,
                                edgecolors=mcolor, c='none')

            if FigMisc.MARK_CA_POS:
                axes[0].scatter(xCA, yCA, zorder=zCA, marker=Style.MS_CA, s=Style.MW_CA,
                                color=mcolor)
                axes[1].scatter(xCA, zCA, zorder=-yCA, marker=Style.MS_CA, s=Style.MW_CA,
                                color=mcolor)
                axes[2].scatter(yCA, zCA, zorder=xCA, marker=Style.MS_CA, s=Style.MW_CA,
                                color=mcolor)
                ax3D.scatter(xCA, yCA, zCA, marker=Style.MS_CA, s=Style.MW_CA,
                             color=mcolor)

    # Set axis limits
    if FigMisc.EXIT_ARROWS:
        # Expand plot limits to make room for arrow tips
        axLim *= 1.10
    [ax.set_xlim([-axLim, axLim]) for ax in axes]
    [ax.set_ylim([-axLim, axLim]) for ax in axes]
    ax3D.set_xlim([-axLim, axLim])
    ax3D.set_ylim([-axLim, axLim])
    ax3D.set_zlim([-axLim, axLim])

    if Params.LEGEND:
        axes[2].legend()
        ax3D.legend()

    fig.tight_layout()
    fig.savefig(Params.Trajec.FigureFiles.SCtrajec, format=FigMisc.figFormat,
                dpi=FigMisc.dpi, metadata=FigLbl.meta)
    log.debug(f'Spacecraft trajectory plot saved to file: {Params.Trajec.FigureFiles.SCtrajec}')

    fig3D.tight_layout()
    fig3D.savefig(Params.Trajec.FigureFiles.SCtrajec3D, format=FigMisc.figFormat,
                  dpi=FigMisc.dpi, metadata=FigLbl.meta)
    log.debug(f'3D trajectory plot saved to file: {Params.Trajec.FigureFiles.SCtrajec3D}')
    plt.close()

    return


def PlotMagCA(PlanetList, Params, scName):

    FlybyCA = GetFlybyCA()
    for Planet in PlanetList:
        Bx, By, Bz = BiTrajecSingle(Planet, Params, scName, np.fromiter(FlybyCA[scName].etCA[Planet.bodyname].values(), dtype=np.float_))
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
