"""
Diagnostic plots for ice convection and tidal dissipation.

Provides visual comparison of PlanetProfile convection parameters against
published scaling laws (Kalousova & Sotin 2018, Deschamps & Sotin 2001)
and tidal dissipation profiles (Tobie et al. 2003, 2005).
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import logging

from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc
from PlanetProfile.Utilities.defineStructs import Constants

log = logging.getLogger('PlanetProfile')


def PlotConvectionDiagnostics(PlanetList, Params):
    """ Multi-panel diagnostic figure showing HP ice convection parameters.

        Panel layout (2 rows x 3 columns):
        - (a) Ra* vs Ra*_c per phase with convection regime
        - (b) Layer structure: eLid, Dconv, deltaTBL stacked bar per phase
        - (c) Convection velocity per phase
        - (d) Temperature profile through hydrosphere with convection zones
        - (e) Viscosity profile through hydrosphere
        - (f) Tidal dissipation H(z) through ice (computed from local mu, eta)

        Only generated when at least one Planet has HP ice convection.
    """
    for Planet in PlanetList:
        if not Planet.Do.VALID:
            continue
        _PlotConvectionDiagnosticsOne(Planet, Params)


def _PlotConvectionDiagnosticsOne(Planet, Params):
    """ Generate convection diagnostic figure for a single Planet. """

    # Collect per-phase convection data
    phases = []
    phase_labels = []
    Ra_star = []
    Ra_crit = []
    eLid_km = []
    Dconv_km = []
    deltaTBL_km = []
    vConv_myr = []
    meltFrac = []
    layerThick_km = []

    # Ice Ih
    if Planet.RaConvect is not None and np.isfinite(Planet.RaConvect) and Planet.RaConvect > 0:
        phases.append('Ih')
        phase_labels.append('Ice Ih')
        Ra_star.append(Planet.RaConvect)
        Ra_crit.append(Planet.RaCrit if Planet.RaCrit is not None else np.nan)
        eLid_km.append((Planet.eLid_m / 1e3) if Planet.eLid_m is not None else 0)
        Dconv_km.append((Planet.Dconv_m / 1e3) if Planet.Dconv_m is not None else 0)
        deltaTBL_km.append((Planet.deltaTBL_m / 1e3) if Planet.deltaTBL_m is not None else 0)
        vConv_myr.append(Planet.vConv_ms * 3.156e7 if np.isfinite(Planet.vConv_ms) else 0)
        meltFrac.append(0)
        layerThick_km.append(Planet.zb_km if hasattr(Planet, 'zb_km') and Planet.zb_km is not None else 0)

    # Ice III
    if Planet.RaConvectIII is not None and np.isfinite(Planet.RaConvectIII) and Planet.RaConvectIII > 0:
        phases.append('III')
        phase_labels.append('Ice III')
        Ra_star.append(Planet.RaConvectIII)
        Ra_crit.append(Planet.RaCritIII if Planet.RaCritIII is not None else np.nan)
        eLid_km.append((Planet.eLidIII_m / 1e3) if Planet.eLidIII_m is not None else 0)
        Dconv_km.append((Planet.DconvIII_m / 1e3) if Planet.DconvIII_m is not None else 0)
        deltaTBL_km.append((Planet.deltaTBLIII_m / 1e3) if Planet.deltaTBLIII_m is not None else 0)
        vConv_myr.append(Planet.vConvIII_ms * 3.156e7 if np.isfinite(Planet.vConvIII_ms) else 0)
        meltFrac.append(Planet.meltFractionIII)
        layerThick_km.append(Planet.dzIceIII_km if hasattr(Planet, 'dzIceIII_km') and Planet.dzIceIII_km is not None else 0)

    # Ice V
    if Planet.RaConvectV is not None and np.isfinite(Planet.RaConvectV) and Planet.RaConvectV > 0:
        phases.append('V')
        phase_labels.append('Ice V')
        Ra_star.append(Planet.RaConvectV)
        Ra_crit.append(Planet.RaCritV if Planet.RaCritV is not None else np.nan)
        eLid_km.append((Planet.eLidV_m / 1e3) if Planet.eLidV_m is not None else 0)
        Dconv_km.append((Planet.DconvV_m / 1e3) if Planet.DconvV_m is not None else 0)
        deltaTBL_km.append((Planet.deltaTBLV_m / 1e3) if Planet.deltaTBLV_m is not None else 0)
        vConv_myr.append(Planet.vConvV_ms * 3.156e7 if np.isfinite(Planet.vConvV_ms) else 0)
        meltFrac.append(Planet.meltFractionV)
        layerThick_km.append(Planet.dzIceV_km if hasattr(Planet, 'dzIceV_km') and Planet.dzIceV_km is not None else 0)

    # Ice VI
    if Planet.RaConvectVI is not None and np.isfinite(Planet.RaConvectVI) and Planet.RaConvectVI > 0:
        phases.append('VI')
        phase_labels.append('Ice VI')
        Ra_star.append(Planet.RaConvectVI)
        Ra_crit.append(Planet.RaCritVI if Planet.RaCritVI is not None else np.nan)
        eLid_km.append((Planet.eLidVI_m / 1e3) if Planet.eLidVI_m is not None else 0)
        Dconv_km.append((Planet.DconvVI_m / 1e3) if Planet.DconvVI_m is not None else 0)
        deltaTBL_km.append((Planet.deltaTBLVI_m / 1e3) if Planet.deltaTBLVI_m is not None else 0)
        vConv_myr.append(Planet.vConvVI_ms * 3.156e7 if np.isfinite(Planet.vConvVI_ms) else 0)
        meltFrac.append(Planet.meltFractionVI)
        layerThick_km.append(Planet.dzIceVI_km if hasattr(Planet, 'dzIceVI_km') and Planet.dzIceVI_km is not None else 0)

    if len(phases) == 0:
        return

    Ra_star = np.array(Ra_star)
    Ra_crit = np.array(Ra_crit)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # --- Panel (a): Ra* vs Ra*_c ---
    ax = axes[0, 0]
    x = np.arange(len(phases))
    width = 0.35
    bars1 = ax.bar(x - width/2, Ra_star, width, label=r'Ra$^*$', color='steelblue')
    bars2 = ax.bar(x + width/2, Ra_crit, width, label=r'Ra$^*_c$', color='salmon')
    ax.set_yscale('log')
    ax.set_xticks(x)
    ax.set_xticklabels(phase_labels)
    ax.set_ylabel('Rayleigh number')
    ax.set_title('(a) Convection regime')
    ax.legend()
    # Mark supercritical with a checkmark
    for i in range(len(phases)):
        if Ra_star[i] > Ra_crit[i]:
            ax.annotate('supercritical', (x[i], Ra_star[i]), textcoords="offset points",
                        xytext=(0, 10), ha='center', fontsize=8, color='green')
        else:
            ax.annotate('subcritical', (x[i], Ra_star[i]), textcoords="offset points",
                        xytext=(0, 10), ha='center', fontsize=8, color='gray')

    # --- Panel (b): Layer structure stacked bar ---
    ax = axes[0, 1]
    bottom_pos = np.zeros(len(phases))
    ax.bar(x, eLid_km, width=0.6, label='Elastic lid', color='lightcoral', bottom=bottom_pos)
    bottom_pos += np.array(eLid_km)
    ax.bar(x, Dconv_km, width=0.6, label='Convecting interior', color='lightskyblue', bottom=bottom_pos)
    bottom_pos += np.array(Dconv_km)
    ax.bar(x, deltaTBL_km, width=0.6, label='Basal TBL', color='sandybrown', bottom=bottom_pos)
    # Show total layer thickness as reference line
    for i in range(len(phases)):
        if layerThick_km[i] > 0:
            ax.plot([x[i] - 0.35, x[i] + 0.35], [layerThick_km[i], layerThick_km[i]],
                    'k--', linewidth=1.5)
    ax.set_xticks(x)
    ax.set_xticklabels(phase_labels)
    ax.set_ylabel('Thickness (km)')
    ax.set_title('(b) Layer structure')
    ax.legend(fontsize=8)

    # --- Panel (c): Convection velocity ---
    ax = axes[0, 2]
    colors = ['steelblue' if v > 0 else 'lightgray' for v in vConv_myr]
    ax.bar(x, vConv_myr, color=colors, width=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(phase_labels)
    ax.set_ylabel('Convection velocity (m/yr)')
    ax.set_title('(c) Convection velocity')
    # Add melt fraction annotations
    for i in range(len(phases)):
        if meltFrac[i] > 0:
            ax.annotate(f'melt: {meltFrac[i]:.1%}', (x[i], vConv_myr[i]),
                        textcoords="offset points", xytext=(0, 10),
                        ha='center', fontsize=8, color='red')

    # --- Panel (d): Temperature profile ---
    ax = axes[1, 0]
    if hasattr(Planet, 'T_K') and Planet.T_K is not None and hasattr(Planet, 'z_m') and Planet.z_m is not None:
        nHydro = Planet.Steps.nHydro if Planet.Steps.nHydro is not None else len(Planet.T_K)
        z_km = Planet.z_m[:nHydro] / 1e3
        T_K = Planet.T_K[:nHydro]
        phase = Planet.phase[:nHydro] if Planet.phase is not None else np.zeros(nHydro)

        # Color by phase
        phase_colors = {0: 'royalblue', 1: 'lightskyblue', 2: 'cyan',
                        3: 'mediumpurple', 5: 'orchid', 6: 'hotpink',
                        -1: 'lightskyblue', -3: 'mediumpurple', -5: 'orchid'}
        for j in range(len(z_km) - 1):
            p = int(phase[j]) if np.isfinite(phase[j]) else 0
            c = phase_colors.get(p, 'gray')
            ax.plot(T_K[j:j+2], z_km[j:j+2], color=c, linewidth=1.5)

        # Add phase labels
        phase_names = {0: 'Ocean', 1: 'Ice Ih', 3: 'Ice III', 5: 'Ice V', 6: 'Ice VI',
                       -1: 'Ice Ih', -3: 'Ice III', -5: 'Ice V'}
        prev_phase = None
        for j in range(len(phase)):
            p = int(phase[j]) if np.isfinite(phase[j]) else 0
            if p != prev_phase and abs(p) in phase_names:
                ax.annotate(phase_names[abs(p)], (T_K[j], z_km[j]),
                            fontsize=7, color='black', alpha=0.7)
                prev_phase = p

    ax.invert_yaxis()
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Depth (km)')
    ax.set_title('(d) Temperature profile')

    # --- Panel (e): Viscosity profile ---
    ax = axes[1, 1]
    if hasattr(Planet, 'eta_Pas') and Planet.eta_Pas is not None:
        nHydro = Planet.Steps.nHydro if Planet.Steps.nHydro is not None else len(Planet.eta_Pas)
        z_km = Planet.z_m[:nHydro] / 1e3
        eta = Planet.eta_Pas[:nHydro]
        # Only plot ice layers (eta is defined)
        ice_mask = eta > 0
        if np.any(ice_mask):
            ax.plot(eta[ice_mask], z_km[ice_mask], 'k-', linewidth=1.5)
    ax.set_xscale('log')
    ax.invert_yaxis()
    ax.set_xlabel('Viscosity (Pa s)')
    ax.set_ylabel('Depth (km)')
    ax.set_title('(e) Viscosity profile')

    # --- Panel (f): Tidal dissipation H(z) ---
    ax = axes[1, 2]
    _PlotTidalDissipationProfile1D(Planet, ax)
    ax.set_title('(f) Tidal dissipation rate')

    fig.suptitle(f'{Planet.name} {Planet.saveLabel} — Convection diagnostics', fontsize=14)
    fig.tight_layout()

    fName = Params.FigureFiles.fName + f'ConvDiag{FigMisc.xtn}'
    fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat,
                dpi=FigMisc.dpi, metadata=FigLbl.meta,
                transparent=FigMisc.TRANSPARENT)
    log.info(f'Convection diagnostics plot saved to {fName}')
    plt.close(fig)


def _PlotTidalDissipationProfile1D(Planet, ax):
    """ Compute and plot tidal dissipation H(r) from local mu, eta for a 1D model.

        Uses the Maxwell dissipation formula:
            H = D * eps0^2  where D = omega * mu / (mu^2 + omega^2 * eta^2)

        Shows dissipation vs depth for all ice layers.
    """
    if not hasattr(Planet, 'eta_Pas') or Planet.eta_Pas is None:
        ax.text(0.5, 0.5, 'No viscosity data', transform=ax.transAxes,
                ha='center', va='center')
        return
    if not hasattr(Planet.Seismic, 'GS_GPa') or Planet.Seismic.GS_GPa is None:
        ax.text(0.5, 0.5, 'No shear modulus data', transform=ax.transAxes,
                ha='center', va='center')
        return

    nHydro = Planet.Steps.nHydro if Planet.Steps.nHydro is not None else 0
    if nHydro == 0:
        return

    # Tidal frequency
    if Planet.Magnetic.Texc_hr is not None and len(Planet.Magnetic.Texc_hr) > 0:
        T_hr = Planet.Magnetic.Texc_hr[0]
    elif hasattr(Planet.Bulk, 'meanMotion_radps') and Planet.Bulk.meanMotion_radps is not None:
        T_hr = 2 * np.pi / Planet.Bulk.meanMotion_radps / 3600
    else:
        T_hr = 85.2  # Default Europa
    omega = 2 * np.pi / (T_hr * 3600)

    # Eccentricity and tidal strain amplitude (Tobie et al. 2003)
    e = Planet.Bulk.eccentricity if hasattr(Planet.Bulk, 'eccentricity') and Planet.Bulk.eccentricity is not None else 0.01
    n = Planet.Bulk.meanMotion_radps if hasattr(Planet.Bulk, 'meanMotion_radps') and Planet.Bulk.meanMotion_radps is not None else omega
    from PlanetProfile.Utilities.defineStructs import Constants
    g_surf = Constants.G * Planet.Bulk.M_kg / Planet.Bulk.R_m**2
    eps0 = 1.5 * e * n**2 * Planet.Bulk.R_m / g_surf

    z_km = Planet.z_m[:nHydro] / 1e3
    mu_Pa = Planet.Seismic.GS_GPa[:nHydro] * 1e9
    eta = Planet.eta_Pas[:nHydro]
    phase = Planet.phase[:nHydro] if Planet.phase is not None else np.zeros(nHydro)

    # Compute Maxwell dissipation for ice layers only
    ice_mask = np.abs(phase) > 0  # Non-ocean
    H_Wm3 = np.full(nHydro, np.nan)
    valid = ice_mask & (mu_Pa > 0) & (eta > 0) & np.isfinite(mu_Pa) & np.isfinite(eta)
    if np.any(valid):
        D = omega**2 * eta[valid] * mu_Pa[valid]**2 / (mu_Pa[valid]**2 + omega**2 * eta[valid]**2)
        H_Wm3[valid] = D * eps0**2

    if not np.any(np.isfinite(H_Wm3)):
        ax.text(0.5, 0.5, 'No ice layers with\nfinite dissipation', transform=ax.transAxes,
                ha='center', va='center')
        return

    # Plot by phase
    phase_colors = {1: 'lightskyblue', 3: 'mediumpurple', 5: 'orchid', 6: 'hotpink',
                    -1: 'lightskyblue', -3: 'mediumpurple', -5: 'orchid'}
    phase_names = {1: 'Ice Ih', 3: 'Ice III', 5: 'Ice V', 6: 'Ice VI'}
    plotted_phases = set()
    for j in range(nHydro - 1):
        if np.isfinite(H_Wm3[j]) and np.isfinite(H_Wm3[j+1]):
            p = abs(int(phase[j]))
            c = phase_colors.get(p, phase_colors.get(int(phase[j]), 'gray'))
            lbl = phase_names.get(p, None) if p not in plotted_phases else None
            ax.plot(H_Wm3[j:j+2] * 1e9, z_km[j:j+2], color=c, linewidth=1.5, label=lbl)
            if lbl:
                plotted_phases.add(p)

    # Add user-specified H_tidal as reference line
    if hasattr(Planet.Ocean, 'HtidalIce_Wm3') and Planet.Ocean.HtidalIce_Wm3 > 0:
        ax.axvline(Planet.Ocean.HtidalIce_Wm3 * 1e9, color='red', linestyle='--',
                   linewidth=1, alpha=0.7, label=f'User H_tidal = {Planet.Ocean.HtidalIce_Wm3:.1e}')

    ax.set_xscale('log')
    ax.invert_yaxis()
    ax.set_xlabel(r'$H_\mathrm{tidal}$ ($10^{-9}$ W/m$^3$)')
    ax.set_ylabel('Depth (km)')
    ax.legend(fontsize=7)


def PlotTidalDissipationProfile(PlanetList, Params):
    """ Standalone tidal dissipation profile plot for 1D models.

        Computes H(r) from local shear modulus and viscosity using Maxwell
        rheology and plots dissipation rate vs depth through all ice layers.
        Shows how dissipation varies with temperature/pressure through the
        shell, comparable to Tobie et al. (2003) Figure 3.
    """
    for Planet in PlanetList:
        if not Planet.Do.VALID:
            continue
        if not hasattr(Planet, 'eta_Pas') or Planet.eta_Pas is None:
            continue
        if not hasattr(Planet.Seismic, 'GS_GPa') or Planet.Seismic.GS_GPa is None:
            continue

        fig, ax = plt.subplots(1, 1, figsize=(6, 8))
        _PlotTidalDissipationProfile1D(Planet, ax)
        ax.set_title(f'{Planet.name} {Planet.saveLabel}')
        fig.tight_layout()

        fName = Params.FigureFiles.fName + f'HtidalProfile{FigMisc.xtn}'
        fig.savefig(fName, bbox_inches='tight', format=FigMisc.figFormat,
                    dpi=FigMisc.dpi, metadata=FigLbl.meta,
                    transparent=FigMisc.TRANSPARENT)
        log.info(f'Tidal dissipation profile saved to {fName}')
        plt.close(fig)
