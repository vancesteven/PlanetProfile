"""
Plotting tools for 3D laterally-varying interior structure.

Provides lat/lon surface maps for ice thickness, tidal heating,
basal temperature, ocean conductivity, and clathrate fraction.
Follows the PlotMagSurface/PlotAsym patterns from MagPlots.py.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import logging

from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc

log = logging.getLogger('PlanetProfile')


def _InterpolateToLatLon(field, theta_rad, phi_rad, latMap_deg, lonMap_deg):
    """ Interpolate a field from pixel grid to regular lat/lon grid.

        Args:
            field: 1D array of values at pixel centers (nPix,).
            theta_rad: Colatitude of pixel centers in radians (nPix,).
            phi_rad: Longitude of pixel centers in radians (nPix,).
            latMap_deg: 1D array of latitudes in degrees for output grid.
            lonMap_deg: 1D array of longitudes in degrees for output grid.

        Returns:
            data2D: 2D array (nLat, nLon) on the output grid.
    """
    from scipy.interpolate import griddata

    # Convert pixel coordinates to lat/lon in degrees
    lat_pix = 90.0 - np.degrees(theta_rad)
    lon_pix = np.degrees(phi_rad)

    # Wrap longitudes to match output range
    lonMin = lonMap_deg[0]
    lonMax = lonMap_deg[-1]
    lon_pix = np.where(lon_pix > lonMax, lon_pix - 360, lon_pix)
    lon_pix = np.where(lon_pix < lonMin, lon_pix + 360, lon_pix)

    # Build output meshgrid
    lonGrid, latGrid = np.meshgrid(lonMap_deg, latMap_deg)

    # Interpolate
    data2D = griddata(
        (lon_pix, lat_pix), field,
        (lonGrid, latGrid),
        method='linear', fill_value=np.nanmean(field)
    )

    return data2D


def _SetMap(ax):
    """ Configure lat/lon tick formatting on a map axis. """
    ax.set_xticks(FigMisc.lonMapTicks_deg)
    ax.set_yticks(FigMisc.latMapTicks_deg)
    ax.tick_params(axis='both', which='major', labelsize=FigMisc.latlonSize)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(FigMisc.LonMapFormatter))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(FigMisc.LatMapFormatter))


def _PlotSurfaceMap(data2D, latMap_deg, lonMap_deg, cmap, cbarLabel,
                    fName, title=None, contourLevels=None, diverging=False):
    """ Generic lat/lon surface map with pcolormesh and optional contours.

        Args:
            data2D: 2D array (nLat, nLon).
            latMap_deg: 1D latitude array in degrees.
            lonMap_deg: 1D longitude array in degrees.
            cmap: Matplotlib colormap name or object.
            cbarLabel: Colorbar label string.
            fName: Output file path.
            title: Optional figure title.
            contourLevels: Optional contour levels (int or array).
            diverging: If True, center colormap at zero.
    """
    fig, ax = plt.subplots(1, 1, figsize=FigSize.asym)

    vmin, vmax = np.nanmin(data2D), np.nanmax(data2D)
    if diverging:
        vlim = max(abs(vmin), abs(vmax))
        vmin, vmax = -vlim, vlim

    mesh = ax.pcolormesh(lonMap_deg, latMap_deg, data2D,
                         shading='auto', cmap=cmap,
                         vmin=vmin, vmax=vmax,
                         rasterized=FigMisc.PT_RASTER)

    if contourLevels is not None:
        try:
            contours = ax.contour(lonMap_deg, latMap_deg, data2D,
                                  levels=contourLevels, colors='black',
                                  linewidths=0.5)
            ax.clabel(contours, fmt=ticker.FuncFormatter(FigMisc.Cformat),
                      fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
        except Exception as e:
            log.debug(f'Could not add contours: {e}')

    _SetMap(ax)
    ax.set_aspect(1)
    fig.colorbar(mesh, ax=ax, label=cbarLabel)

    if title is not None:
        if not FigMisc.TEX_INSTALLED:
            title = FigLbl.StripLatexFromString(title)
        ax.set_title(title)

    fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat,
                dpi=FigMisc.dpi, metadata=FigLbl.meta,
                transparent=FigMisc.TRANSPARENT)
    log.info(f'Lateral plot saved to {fName}')
    plt.close(fig)


def PlotIceThickness(Planet, Params):
    """ Plot lat/lon map of ice shell thickness variation. """
    Lateral = Planet.Lateral
    if Lateral.dIce_m is None:
        return

    data2D = _InterpolateToLatLon(
        Lateral.dIce_m / 1e3,
        Lateral.theta_rad, Lateral.phi_rad,
        FigMisc.latMap_deg, FigMisc.lonMap_deg
    )

    fName = os.path.join(Params.FigureFiles.path,
                         f'{Planet.name}_dIce{FigMisc.xtn}')

    nLevels = 8
    _PlotSurfaceMap(
        data2D, FigMisc.latMap_deg, FigMisc.lonMap_deg,
        cmap='viridis', cbarLabel='Ice thickness (km)',
        fName=fName, title=f'{Planet.name} ice shell thickness',
        contourLevels=nLevels
    )


def PlotTidalHeating(Planet, Params):
    """ Plot lat/lon map of tidal surface heat flux in ice shell. """
    Lateral = Planet.Lateral
    if Lateral.HtidalIce_Wm3 is None or Lateral.dIce_m is None:
        return

    # Surface heat flux: integrate volumetric heating over ice thickness
    qTidal_mWm2 = Lateral.HtidalIce_Wm3 * Lateral.dIce_m * 1e3  # mW/m^2

    data2D = _InterpolateToLatLon(
        qTidal_mWm2,
        Lateral.theta_rad, Lateral.phi_rad,
        FigMisc.latMap_deg, FigMisc.lonMap_deg
    )

    fName = os.path.join(Params.FigureFiles.path,
                         f'{Planet.name}_Htidal{FigMisc.xtn}')

    _PlotSurfaceMap(
        data2D, FigMisc.latMap_deg, FigMisc.lonMap_deg,
        cmap='magma', cbarLabel=r'Tidal surface heat flux ($\mathrm{mW\,m^{-2}}$)',
        fName=fName, title=f'{Planet.name} tidal surface heat loss'
    )


def PlotBasalTemperature(Planet, Params):
    """ Plot lat/lon map of basal ice temperature. """
    Lateral = Planet.Lateral
    if Lateral.Tb_K is None or np.all(np.isnan(Lateral.Tb_K)):
        return

    Tb_mean = np.nanmean(Lateral.Tb_K)
    dTb = Lateral.Tb_K - Tb_mean

    data2D = _InterpolateToLatLon(
        dTb, Lateral.theta_rad, Lateral.phi_rad,
        FigMisc.latMap_deg, FigMisc.lonMap_deg
    )

    fName = os.path.join(Params.FigureFiles.path,
                         f'{Planet.name}_Tb{FigMisc.xtn}')

    _PlotSurfaceMap(
        data2D, FigMisc.latMap_deg, FigMisc.lonMap_deg,
        cmap='seismic', cbarLabel=f'$T_b - \\overline{{T_b}}$ (K), mean={Tb_mean:.1f} K',
        fName=fName, title=f'{Planet.name} basal temperature deviation',
        diverging=True, contourLevels=8
    )


def PlotOceanConductivity(Planet, Params):
    """ Plot lat/lon map of mean ocean electrical conductivity. """
    Lateral = Planet.Lateral
    if Lateral.sigma_mean_Sm is None or np.all(np.isnan(Lateral.sigma_mean_Sm)):
        return

    data2D = _InterpolateToLatLon(
        Lateral.sigma_mean_Sm,
        Lateral.theta_rad, Lateral.phi_rad,
        FigMisc.latMap_deg, FigMisc.lonMap_deg
    )

    fName = os.path.join(Params.FigureFiles.path,
                         f'{Planet.name}_sigmaOcean{FigMisc.xtn}')

    _PlotSurfaceMap(
        data2D, FigMisc.latMap_deg, FigMisc.lonMap_deg,
        cmap='cividis',
        cbarLabel=r'$\overline{\sigma}_\mathrm{ocean}$ (S/m)',
        fName=fName, title=f'{Planet.name} mean ocean conductivity'
    )


def PlotClathrateFraction(Planet, Params):
    """ Plot lat/lon map of clathrate volume fraction. """
    Lateral = Planet.Lateral
    if not Lateral.DO_CLATH_LATERAL or Lateral.fClath is None:
        return

    data2D = _InterpolateToLatLon(
        Lateral.fClath,
        Lateral.theta_rad, Lateral.phi_rad,
        FigMisc.latMap_deg, FigMisc.lonMap_deg
    )

    fName = os.path.join(Params.FigureFiles.path,
                         f'{Planet.name}_fClath{FigMisc.xtn}')

    _PlotSurfaceMap(
        data2D, FigMisc.latMap_deg, FigMisc.lonMap_deg,
        cmap='YlOrBr', cbarLabel='Clathrate fraction',
        fName=fName, title=f'{Planet.name} clathrate volume fraction',
        contourLevels=6
    )


def PlotEffectiveConductivity(Planet, Params):
    """ Plot lat/lon map of effective thermal conductivity. """
    Lateral = Planet.Lateral
    if Lateral.kThermEff_WmK is None:
        return

    data2D = _InterpolateToLatLon(
        Lateral.kThermEff_WmK,
        Lateral.theta_rad, Lateral.phi_rad,
        FigMisc.latMap_deg, FigMisc.lonMap_deg
    )

    fName = os.path.join(Params.FigureFiles.path,
                         f'{Planet.name}_kTherm{FigMisc.xtn}')

    _PlotSurfaceMap(
        data2D, FigMisc.latMap_deg, FigMisc.lonMap_deg,
        cmap='coolwarm',
        cbarLabel=r'$k_\mathrm{eff}$ (W m$^{-1}$ K$^{-1}$)',
        fName=fName, title=f'{Planet.name} effective thermal conductivity'
    )


def PlotLateralSummary(Planet, Params):
    """ Multi-panel figure combining key lateral structure fields. """
    Lateral = Planet.Lateral
    if Lateral.dIce_m is None:
        return

    latMap = FigMisc.latMap_deg
    lonMap = FigMisc.lonMap_deg

    # Determine which panels to show
    panels = []
    panels.append(('Ice thickness (km)', Lateral.dIce_m / 1e3, 'viridis', False))

    if Lateral.HtidalIce_Wm3 is not None and Lateral.dIce_m is not None:
        panels.append(('Tidal heat flux (mW/m$^2$)', Lateral.HtidalIce_Wm3 * Lateral.dIce_m * 1e3, 'magma', False))

    if Lateral.Tb_K is not None and not np.all(np.isnan(Lateral.Tb_K)):
        panels.append(('$\\Delta T_b$ (K)', Lateral.Tb_K - np.nanmean(Lateral.Tb_K), 'seismic', True))

    if Lateral.sigma_mean_Sm is not None and not np.all(np.isnan(Lateral.sigma_mean_Sm)):
        panels.append(('$\\overline{\\sigma}$ (S/m)', Lateral.sigma_mean_Sm, 'cividis', False))

    if Lateral.DO_CLATH_LATERAL and Lateral.fClath is not None:
        panels.append(('Clathrate fraction', Lateral.fClath, 'YlOrBr', False))

    nPanels = len(panels)
    if nPanels == 0:
        return

    nCols = min(nPanels, 3)
    nRows = int(np.ceil(nPanels / nCols))
    fig, axes = plt.subplots(nRows, nCols, figsize=(5 * nCols, 4 * nRows))
    if nPanels == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for i, (label, field, cmap, diverging) in enumerate(panels):
        ax = axes[i]
        data2D = _InterpolateToLatLon(
            field, Lateral.theta_rad, Lateral.phi_rad, latMap, lonMap
        )

        vmin, vmax = np.nanmin(data2D), np.nanmax(data2D)
        if diverging:
            vlim = max(abs(vmin), abs(vmax))
            vmin, vmax = -vlim, vlim

        mesh = ax.pcolormesh(lonMap, latMap, data2D, shading='auto',
                             cmap=cmap, vmin=vmin, vmax=vmax,
                             rasterized=FigMisc.PT_RASTER)
        _SetMap(ax)
        ax.set_aspect(1)
        fig.colorbar(mesh, ax=ax, label=label, shrink=0.8)

    # Hide unused axes
    for j in range(nPanels, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(f'{Planet.name} lateral structure summary', fontsize=12)
    fig.tight_layout()

    fName = os.path.join(Params.FigureFiles.path,
                         f'{Planet.name}_lateralSummary{FigMisc.xtn}')
    fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat,
                dpi=FigMisc.dpi, metadata=FigLbl.meta,
                transparent=FigMisc.TRANSPARENT)
    log.info(f'Lateral summary plot saved to {fName}')
    plt.close(fig)


def PlotTidalHeatingByLayer(Planet, Params):
    """ Multi-panel plot of tidal dissipation at top and bottom of each
        ice layer type, following Tobie et al. (2005) Figure 10.

        Shows % deviation from the mean at each position. Panels are
        arranged as rows (layer type) x columns (bottom, top).
    """
    Lateral = Planet.Lateral
    latMap = FigMisc.latMap_deg
    lonMap = FigMisc.lonMap_deg

    # Collect panels: (label, field, mean_value)
    panels = []
    if Lateral.HtidalIceI_bot_Wm3 is not None and np.any(Lateral.HtidalIceI_bot_Wm3 > 0):
        mean_bot = np.mean(Lateral.HtidalIceI_bot_Wm3[Lateral.HtidalIceI_bot_Wm3 > 0])
        mean_top = np.mean(Lateral.HtidalIceI_top_Wm3[Lateral.HtidalIceI_top_Wm3 > 0])
        panels.append(('Bottom of ice I', Lateral.HtidalIceI_bot_Wm3, mean_bot))
        panels.append(('Top of ice I', Lateral.HtidalIceI_top_Wm3, mean_top))
    if Lateral.HtidalHP_bot_Wm3 is not None and np.any(Lateral.HtidalHP_bot_Wm3 > 0):
        hp_bot_valid = Lateral.HtidalHP_bot_Wm3[Lateral.HtidalHP_bot_Wm3 > 0]
        hp_top_valid = Lateral.HtidalHP_top_Wm3[Lateral.HtidalHP_top_Wm3 > 0]
        if len(hp_bot_valid) > 0:
            mean_hp_bot = np.mean(hp_bot_valid)
            mean_hp_top = np.mean(hp_top_valid) if len(hp_top_valid) > 0 else 0
            panels.append(('Bottom of HP ice', Lateral.HtidalHP_bot_Wm3, mean_hp_bot))
            panels.append(('Top of HP ice', Lateral.HtidalHP_top_Wm3, mean_hp_top))

    if len(panels) == 0:
        return

    nPanels = len(panels)
    nRows = nPanels // 2
    nCols = 2
    fig, axes = plt.subplots(nRows, nCols, figsize=(5 * nCols, 4 * nRows))
    if nRows == 1:
        axes = axes.reshape(1, -1)

    for idx, (label, field, mean_val) in enumerate(panels):
        row = idx // 2
        col = idx % 2
        ax = axes[row, col]

        # Compute % deviation from mean
        pct_dev = np.where(field > 0, 100.0 * (field - mean_val) / mean_val, 0.0)

        data2D = _InterpolateToLatLon(
            pct_dev, Lateral.theta_rad, Lateral.phi_rad, latMap, lonMap
        )

        vlim = max(abs(np.nanmin(data2D)), abs(np.nanmax(data2D)))
        if vlim == 0:
            vlim = 1.0

        mesh = ax.pcolormesh(lonMap, latMap, data2D, shading='auto',
                             cmap='RdBu_r', vmin=-vlim, vmax=vlim,
                             rasterized=FigMisc.PT_RASTER)
        _SetMap(ax)
        ax.set_aspect(1)
        fig.colorbar(mesh, ax=ax, label='% deviation from mean', shrink=0.8)

        title = f'{label}, mean = {mean_val:.1e} W/m$^3$'
        if not FigMisc.TEX_INSTALLED:
            title = FigLbl.StripLatexFromString(title)
        ax.set_title(title, fontsize=10)

    fig.suptitle(f'{Planet.name} tidal dissipation by layer (cf. Tobie et al. 2005, Fig. 10)',
                 fontsize=12)
    fig.tight_layout()

    fName = os.path.join(Params.FigureFiles.path,
                         f'{Planet.name}_HtidalByLayer{FigMisc.xtn}')
    fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat,
                dpi=FigMisc.dpi, metadata=FigLbl.meta,
                transparent=FigMisc.TRANSPARENT)
    log.info(f'Tidal heating by-layer plot saved to {fName}')
    plt.close(fig)


def PlotTidalHeatingVsDepth(Planet, Params):
    """ Plot tidal dissipation rate vs depth for selected columns.

        Shows H(r) profiles for a representative set of grid points
        (e.g., sub-planetary, anti-planetary, pole, equatorial) to
        illustrate how dissipation varies through the ice shell.
    """
    Lateral = Planet.Lateral
    if Lateral.HtidalIceI_profile_Wm3 is None:
        return

    # Select representative columns: find indices closest to key locations
    # Sub-planetary point (0, 0), anti-planetary (0, 180), north pole (90, 0), equator-flank (0, 90)
    target_locations = {
        'Sub-planetary (0N, 0E)': (np.pi / 2, 0),
        'Anti-planetary (0N, 180E)': (np.pi / 2, np.pi),
        'North pole (90N)': (0.01, 0),
        'Equator flank (0N, 90E)': (np.pi / 2, np.pi / 2),
    }

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    ax_iceI = axes[0]
    ax_hp = axes[1]
    colors = plt.cm.tab10(np.linspace(0, 1, len(target_locations)))

    has_iceI = False
    has_hp = False
    for idx, (label, (theta_t, phi_t)) in enumerate(target_locations.items()):
        # Find nearest pixel
        dist = np.sqrt((Lateral.theta_rad - theta_t)**2 + (Lateral.phi_rad - phi_t)**2)
        iPix = np.argmin(dist)

        # Ice I profile
        prof = Lateral.HtidalIceI_profile_Wm3[iPix]
        radii = Lateral.rIceI_profile_m[iPix]
        if prof is not None and radii is not None:
            depth_km = (Planet.Bulk.R_m - radii) / 1e3
            ax_iceI.plot(prof * 1e9, depth_km, color=colors[idx], label=label, linewidth=1.5)
            has_iceI = True

        # HP ice profile
        hp_prof = Lateral.HtidalHP_profile_Wm3[iPix] if Lateral.HtidalHP_profile_Wm3 is not None else None
        hp_radii = Lateral.rHP_profile_m[iPix] if Lateral.rHP_profile_m is not None else None
        if hp_prof is not None and hp_radii is not None:
            depth_km = (Planet.Bulk.R_m - hp_radii) / 1e3
            ax_hp.plot(hp_prof * 1e9, depth_km, color=colors[idx], label=label, linewidth=1.5)
            has_hp = True

    for ax, title_str, has_data in [(ax_iceI, 'Ice I', has_iceI), (ax_hp, 'HP ice', has_hp)]:
        if has_data:
            ax.set_xlabel(r'Tidal heating ($10^{-9}$ W/m$^3$)')
            ax.set_ylabel('Depth (km)')
            ax.invert_yaxis()
            ax.legend(fontsize=8)
            ax.set_title(title_str)
        else:
            ax.set_visible(False)

    fig.suptitle(f'{Planet.name} tidal dissipation vs depth', fontsize=12)
    fig.tight_layout()

    fName = os.path.join(Params.FigureFiles.path,
                         f'{Planet.name}_HtidalVsDepth{FigMisc.xtn}')
    fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat,
                dpi=FigMisc.dpi, metadata=FigLbl.meta,
                transparent=FigMisc.TRANSPARENT)
    log.info(f'Tidal heating vs depth plot saved to {fName}')
    plt.close(fig)


def GenerateLateralPlots(Planet, Params):
    """ Generate all lateral structure plots.

        Called from RunLateral3D after saving results.

        Args:
            Planet: PlanetStruct with completed Lateral substruct.
            Params: ParamsStruct.
    """
    if not Planet.Lateral.DO_3D or Planet.Lateral.dIce_m is None:
        return

    if not Params.PLOT_LATERAL:
        log.debug('Lateral plotting disabled (Params.PLOT_LATERAL = False)')
        return

    log.info('Generating lateral structure plots...')

    PlotIceThickness(Planet, Params)
    PlotTidalHeating(Planet, Params)
    PlotTidalHeatingByLayer(Planet, Params)
    PlotTidalHeatingVsDepth(Planet, Params)
    PlotBasalTemperature(Planet, Params)
    PlotOceanConductivity(Planet, Params)
    PlotClathrateFraction(Planet, Params)
    PlotEffectiveConductivity(Planet, Params)
    PlotLateralSummary(Planet, Params)

    log.info('Lateral structure plots complete.')
