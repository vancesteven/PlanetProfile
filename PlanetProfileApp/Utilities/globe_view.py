"""Interactive 3D globe of the modeled body for the inference results view.

Renders a plotly figure the user can grab and rotate: the body's surface
with the best available global mosaic projected onto it (NASA/USGS public
domain, shipped in PlanetProfileApp/assets/globes/<body>.jpg), made
semi-transparent so the concentric interior layers (ice shell, ocean,
high-pressure ices, silicate, core) shine through at their posterior
relative depths.

Degree-2 figure: even 1D inference carries a NON-SPHERICAL shape — the
measured C20/C22 define a triaxial equipotential figure that the
spherically averaged C/MR^2 view never shows. The surface mesh is
deformed by the degree-2 shape implied by (C20, C22) with a fluid-body
amplification h_f ~ 1 + k_f (surface follows the equipotential), under a
user exaggeration factor so the (tiny, for most bodies) flattening is
visible. For future 3D (laterally varying) models the same mesh accepts
per-layer r(theta, phi) fields.

Conventions: colatitude theta from +z (spin axis), longitude phi from +x
(the sub-parent axis for synchronous rotators — the tidal bulge lies
along +/-x). Unnormalized C20/C22; shape radius
r(theta, phi) = R * [1 - h_amp * (C20 * P20(cos theta)
                                  + 3 * C22 * cos(2 phi) * sin^2 theta)]
with P20 = (3 cos^2 theta - 1)/2 and the minus sign because C20 < 0
(oblate) must FLATTEN the poles.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

_ASSETS = Path(__file__).resolve().parents[1] / 'assets' / 'globes'

# Spin rate + mass for the hydrostatic display figure when a run has no
# C20/C22 observables (values match the inference cache_builder table).
BODY_ROTATION = {
    'europa':    {'omega': 2.0479e-5, 'M_kg': 4.7998e22},
    'ganymede':  {'omega': 1.0163e-5, 'M_kg': 1.4819e23},
    'callisto':  {'omega': 4.358e-6,  'M_kg': 1.0759e23},
    'titan':     {'omega': 4.560e-6,  'M_kg': 1.3452e23},
    'enceladus': {'omega': 5.307e-5,  'M_kg': 1.08022e20},
}
_G = 6.674e-11


def hydrostatic_display_pair(bodyname, R_km, kf):
    """Hydrostatic (C20, C22) for DISPLAY when a run carries no gravity
    observables: C22 = kf*q_r/4 at the physical radius, C20 = -(10/3)*C22
    (classical ratio suffices for visualization)."""
    rot = BODY_ROTATION.get(str(bodyname).strip().lower())
    if rot is None or kf is None or not np.isfinite(kf):
        return None, None
    R_m = float(R_km) * 1e3
    q_r = rot['omega'] ** 2 * R_m ** 3 / (_G * rot['M_kg'])
    c22 = 0.25 * float(kf) * q_r
    return -(10.0 / 3.0) * c22, c22

# Layer colors. Ice phases run pale (Ih) -> deep cyan (VI) with rising
# pressure; the ocean is a distinct saturated navy so it never reads as
# an ice band.
LAYER_COLORS = {
    'core': 'rgb(120, 110, 100)',
    'silicate': 'rgb(150, 96, 52)',
    # Ocean: saturated royal blue (liquid). Ice phases walk across hue
    # (pale white -> cyan -> teal -> violet) so adjacent HP shells are
    # clearly distinct, not one blue mass.
    'ocean': 'rgb(30, 90, 220)',
    'ice_Ih': 'rgb(222, 238, 250)',   # pale white-blue skin
    'ice_III': 'rgb(120, 224, 224)',  # cyan
    'ice_V': 'rgb(70, 170, 150)',     # teal-green
    'ice_VI': 'rgb(150, 120, 210)',   # violet
    # legacy aliases
    'hp_ice': 'rgb(120, 224, 224)',
    'ice': 'rgb(222, 238, 250)',
}

# Friendly legend/hover names per layer kind.
LAYER_LABELS = {
    'core': 'Metallic core',
    'silicate': 'Silicate mantle',
    'ocean': 'Liquid ocean',
    'ice_Ih': 'Ice Ih',
    'ice_III': 'Ice III',
    'ice_V': 'Ice V',
    'ice_VI': 'Ice VI',
}


def texture_path(bodyname: str) -> Optional[Path]:
    p = _ASSETS / f'{bodyname.strip().lower()}.jpg'
    return p if p.exists() else None


def _sphere_mesh(n_lat: int = 60, n_lon: int = 120):
    theta = np.linspace(1e-3, np.pi - 1e-3, n_lat)   # colatitude
    phi = np.linspace(0.0, 2.0 * np.pi, n_lon)
    T, P = np.meshgrid(theta, phi, indexing='ij')
    return T, P


def _degree2_radius(R: float, T: np.ndarray, P: np.ndarray,
                    c20: float = 0.0, c22: float = 0.0,
                    h_amp: float = 1.0, exaggeration: float = 1.0):
    """Degree-2 equipotential shape r(theta, phi).

    Geoid of the measured field: r/R = 1 + C20*P20(cos theta)
    + 3*C22*cos(2 phi)*sin^2(theta) — C20 < 0 flattens the poles,
    C22 > 0 bulges the sub-/anti-parent (+/-x) axis. h_amp carries the
    fluid-body amplification (1 + k_f)/k_f: the measured coefficients
    are the INDUCED response k_f times the direct tidal/centrifugal
    potential, and the fluid surface displacement follows
    h_f = 1 + k_f times the direct potential.
    """
    p20 = 0.5 * (3.0 * np.cos(T) ** 2 - 1.0)
    dev = h_amp * (c20 * p20 + 3.0 * c22 * np.cos(2.0 * P) * np.sin(T) ** 2)
    return R * (1.0 + exaggeration * dev)


def _surface_trace(R_km: float, T, P, texture_img=None,
                   c20=0.0, c22=0.0, h_amp=1.0, exaggeration=1.0,
                   opacity=0.65, name='surface'):
    import plotly.graph_objects as go

    r = _degree2_radius(R_km, T, P, c20, c22, h_amp, exaggeration)
    x = r * np.sin(T) * np.cos(P)
    y = r * np.sin(T) * np.sin(P)
    z = r * np.cos(T)

    if texture_img is not None:
        # Grid-exact texture: LANCZOS-resample the equirectangular map to
        # the mesh grid, then map longitudes (mesh phi range may be a
        # cutaway sub-interval). Grayscale + 'gray' colorscale (robust in
        # plotly Surface; icy-moon mosaics read well in monochrome).
        from PIL import Image
        h, w = texture_img.shape[:2]
        img = Image.fromarray(texture_img).convert('L')
        n_lat_m, n_lon_m = T.shape
        full = np.asarray(img.resize((2 * n_lon_m, n_lat_m),
                                     Image.LANCZOS), dtype=float)
        lon_frac = (P / (2.0 * np.pi)) % 1.0
        col = np.clip((lon_frac * (full.shape[1] - 1)).astype(int),
                      0, full.shape[1] - 1)
        rows = np.arange(n_lat_m)[:, None] * np.ones((1, n_lon_m), int)
        shade = full[rows, col]
        surfargs = dict(surfacecolor=shade, colorscale='gray',
                        cmin=0, cmax=255)
    else:
        surfargs = dict(surfacecolor=np.cos(T), colorscale='Greys',
                        cmin=-2, cmax=2)

    return go.Surface(
        x=x, y=y, z=z, opacity=opacity, showscale=False, showlegend=False,
        name=name, hovertemplate=f'{name}<extra></extra>', **surfargs)


def _layer_sphere_trace(r_km: float, T, P, color: str, name: str,
                        opacity: float = 0.35, c20: float = 0.0,
                        c22: float = 0.0, h_amp: float = 1.0,
                        exaggeration: float = 1.0):
    import plotly.graph_objects as go
    # Interior interfaces follow (to first order for visualization) the
    # same relative degree-2 figure as the surface — without this, at
    # high exaggeration the spherical shells poke through the deformed
    # surface.
    r = _degree2_radius(r_km, T, P, c20, c22, h_amp, exaggeration)
    x = r * np.sin(T) * np.cos(P)
    y = r * np.sin(T) * np.sin(P)
    z = r * np.cos(T)
    return go.Surface(
        x=x, y=y, z=z, opacity=opacity, showscale=False, showlegend=False,
        name=name, surfacecolor=np.zeros_like(x),
        colorscale=[[0, color], [1, color]],
        hovertemplate=f'{name}: r = {r_km:.0f} km<extra></extra>')


def _cut_face_trace(r_in, r_out, phi, color, name,
                    c20=0.0, c22=0.0, h_amp=1.0, exaggeration=1.0):
    """A FILLED radial cross-section face on the meridional half-plane at
    longitude phi, spanning colatitude [0, pi] and radius [r_in, r_out].
    These caps on the two wedge boundaries give each layer real, filled
    thickness in the cutaway (a hollow spherical shell would show zero
    thickness / see-through)."""
    import plotly.graph_objects as go
    th = np.linspace(0.0, np.pi, 60)
    rr = np.linspace(float(r_in), float(r_out), 6)
    TH, RR = np.meshgrid(th, rr, indexing='ij')
    P = np.full_like(TH, float(phi))
    Rdef = _degree2_radius(RR, TH, P, c20, c22, h_amp, exaggeration)
    x = Rdef * np.sin(TH) * np.cos(P)
    y = Rdef * np.sin(TH) * np.sin(P)
    z = Rdef * np.cos(TH)
    return go.Surface(
        x=x, y=y, z=z, opacity=1.0, showscale=False, showlegend=False,
        name=name, surfacecolor=np.zeros_like(x),
        colorscale=[[0, color], [1, color]],
        hovertemplate=(f'{name}: {r_in:.0f}–{r_out:.0f} km'
                       '<extra></extra>'))


def triaxial_semiaxes(R_km, c20, c22, h_amp, exaggeration=1.0):
    """Semi-axes (a, b, c) of the degree-2 figure [km]: a along the
    sub-parent x-axis, b along the orbital y-axis, c along the spin
    z-axis."""
    a = _degree2_radius(R_km, np.array(np.pi / 2), np.array(0.0),
                        c20, c22, h_amp, exaggeration)
    b = _degree2_radius(R_km, np.array(np.pi / 2), np.array(np.pi / 2),
                        c20, c22, h_amp, exaggeration)
    c = _degree2_radius(R_km, np.array(0.0), np.array(0.0),
                        c20, c22, h_amp, exaggeration)
    return float(a), float(b), float(c)


def _axis_traces(a, b, c):
    """Principal-axis lines with labels: the principal moments satisfy
    A < B < C about the a (sub-parent), b (orbital), and c (spin) axes."""
    import plotly.graph_objects as go
    ext = 1.22
    specs = [
        ('A', (a * ext, 0, 0), 'rgb(255, 120, 90)',
         f'a-axis (sub-parent): {a:.1f} km — smallest principal moment A'),
        ('B', (0, b * ext, 0), 'rgb(120, 230, 120)',
         f'b-axis (orbital): {b:.1f} km — intermediate principal moment B'),
        ('C', (0, 0, c * ext), 'rgb(120, 190, 255)',
         f'c-axis (spin): {c:.1f} km — largest principal moment C'),
    ]
    traces = []
    for label, tip, color, hover in specs:
        traces.append(go.Scatter3d(
            x=[-tip[0], tip[0]], y=[-tip[1], tip[1]],
            z=[-tip[2], tip[2]],
            mode='lines+text', text=['', label],
            textposition='top center',
            textfont=dict(color=color, size=15),
            line=dict(color=color, width=5), showlegend=False,
            name=label, hovertemplate=hover + '<extra></extra>'))
    return traces


def build_globe_figure(
    bodyname: str,
    R_body_km: float,
    layers: List[Dict],
    c20: Optional[float] = None,
    c22: Optional[float] = None,
    kf: Optional[float] = None,
    exaggeration: float = 1.0,
    show_axes: bool = True,
    n_lat: int = 160,
    n_lon: int = 320,
):
    """Assemble the interactive globe.

    layers: list of {'name': str, 'r_km': float, 'kind': key into
    LAYER_COLORS}, ordered inner->outer; each is drawn as a
    semi-transparent concentric sphere. The outer surface uses the body
    texture when shipped, deformed by (c20, c22) with amplification
    h_f ~ 1 + k_f and the given exaggeration.
    """
    import plotly.graph_objects as go

    # Cutaway design (user 2026-07-20): a 90-degree longitude wedge
    # [0, pi/2] is removed. The intact 3/4 shows the textured exterior;
    # the wedge exposes a FILLED radial cross-section — each layer painted
    # as a solid annular face on the two meridional cut planes, so every
    # phase (incl. Ice Ih) has real proportional thickness rather than a
    # see-through hollow shell.
    theta_s = np.linspace(1e-3, np.pi - 1e-3, n_lat)
    phi_s = np.linspace(0.5 * np.pi, 2.0 * np.pi, n_lon)
    T_surf, P_surf = np.meshgrid(theta_s, phi_s, indexing='ij')

    texture_img = None
    tp = texture_path(bodyname)
    if tp is not None:
        from PIL import Image
        texture_img = np.asarray(Image.open(tp).convert('RGB'))

    _kf = kf if (kf is not None and np.isfinite(kf) and kf > 0.05) else 1.0
    h_amp = (1.0 + _kf) / _kf
    _dk = dict(c20=(c20 or 0.0), c22=(c22 or 0.0), h_amp=h_amp,
               exaggeration=exaggeration)

    fig = go.Figure()
    valid_layers = [l for l in layers
                    if l.get('r_km') and np.isfinite(l['r_km'])
                    and l['r_km'] > 0]
    ordered = sorted(valid_layers, key=lambda l: -l['r_km'])

    # Textured exterior over the intact 3/4 (opaque, drawn first).
    fig.add_trace(_surface_trace(
        float(R_body_km), T_surf, P_surf, texture_img=texture_img,
        opacity=1.0, name=f'{bodyname} surface', **_dk))

    # Filled cross-section: consecutive outer radii give each layer's
    # [r_in, r_out]; innermost fills to the center. Two faces per layer
    # (the phi=0 and phi=pi/2 wedge boundaries).
    radii = [float(l['r_km']) for l in ordered] + [0.0]
    for k, layer in enumerate(ordered):
        r_out, r_in = radii[k], radii[k + 1]
        if r_out - r_in < 1e-3:
            continue
        color = LAYER_COLORS.get(layer.get('kind', 'silicate'),
                                 LAYER_COLORS['silicate'])
        name = layer.get('name', layer.get('kind', 'layer'))
        for phi in (0.0, 0.5 * np.pi):
            fig.add_trace(_cut_face_trace(r_in, r_out, phi, color, name,
                                          **_dk))

    # Legend: one color swatch per layer present (radial order, incl. the
    # Ih surface skin), via invisible Scatter3d proxies — go.Surface has no
    # reliable legend entry.
    import plotly.graph_objects as _go
    for layer in ordered:
        color = LAYER_COLORS.get(layer.get('kind', 'silicate'),
                                 LAYER_COLORS['silicate'])
        fig.add_trace(_go.Scatter3d(
            x=[None], y=[None], z=[None], mode='markers',
            marker=dict(size=9, color=color, symbol='square'),
            name=layer.get('name', layer.get('kind', 'layer')),
            showlegend=True, hoverinfo='skip'))

    a, b, c = triaxial_semiaxes(float(R_body_km), (c20 or 0.0),
                                (c22 or 0.0), h_amp, exaggeration)
    if show_axes:
        for tr in _axis_traces(a, b, c):
            fig.add_trace(tr)

    lim = 1.32 * max(R_body_km, a, b, c)
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False, range=[-lim, lim]),
            yaxis=dict(visible=False, range=[-lim, lim]),
            zaxis=dict(visible=False, range=[-lim, lim]),
            aspectmode='cube',
            bgcolor='rgb(3, 6, 18)',
            # Angle showing both the cut face (wedge bisector phi = pi/4)
            # and a good stretch of the textured surface
            camera=dict(eye=dict(x=1.75, y=0.75, z=0.65)),
        ),
        margin=dict(l=0, r=0, t=0, b=0),
        height=560,
        showlegend=True,
        legend=dict(x=0.01, y=0.99, bgcolor='rgba(3,6,18,0.6)',
                    font=dict(color='rgb(220,225,235)', size=11),
                    itemsizing='constant'),
        paper_bgcolor='rgb(3, 6, 18)',
    )
    return fig


def posterior_median_layers(result) -> Tuple[float, List[Dict]]:
    """Extract (R_body_km, inner->outer layer list) from an
    InferenceResult using posterior medians. Works for both MCMC and
    amortized results (same packaging)."""
    md = getattr(result, 'metadata', {}) or {}
    R_km = md.get('R_body_km')
    if isinstance(R_km, dict):
        R_km = list(R_km.values())[0]
    if R_km is None or not np.isfinite(R_km):
        return None, []
    R_km = float(R_km)

    def _med(attr):
        v = getattr(result, attr, None)
        if v is None:
            return None
        v = np.asarray(v, float)
        v = v[np.isfinite(v)]
        return float(np.median(v)) if v.size else None

    d_hs = _med('D_hsphere_results') or 0.0
    d_oc = _med('D_ocean_results') or 0.0
    d_ih = _med('D_iceIh_results') or 0.0

    phases = getattr(result, 'D_icePhase_results', None) or {}

    def _phase_med(key):
        v = phases.get(key)
        if v is None:
            return 0.0
        v = np.asarray(v, float)
        v = v[np.isfinite(v)]
        return float(np.median(v)) if v.size else 0.0

    d_iii, d_v, d_vi = _phase_med('III'), _phase_med('V'), _phase_med('VI')

    r_core = None
    names = list(getattr(result, 'param_names', []))
    if 'R_core_km' in names:
        r_core = float(np.median(
            np.asarray(result.samples, float)[:, names.index('R_core_km')]))

    return R_km, _build_radial_layers(
        R_km, r_core, d_hs, d_ih, d_oc, d_iii, d_v, d_vi)


def _build_radial_layers(R_km, r_core, d_hs, d_ih, d_oc,
                         d_iii, d_v, d_vi, min_km=0.5):
    """Assemble the concentric layer list (each drawn at its OUTER radius)
    from the surface inward, one labeled shell per material with finite
    proportional thickness: Ice Ih (the textured surface skin) -> ocean ->
    Ice III -> Ice V -> Ice VI -> silicate mantle -> core. Layers thinner
    than min_km are omitted (a body simply lacks that phase)."""
    layers: List[Dict] = []
    r = float(R_km)  # running OUTER radius, marching inward

    # Ice Ih is the outer skin (the textured surface at R). Emit it as a
    # legend/label proxy only — its band is the gap from R down to the
    # next interface. (r moves to its base.)
    if d_ih > min_km:
        layers.append({'name': LAYER_LABELS['ice_Ih'], 'r_km': r,
                       'kind': 'ice_Ih', 'surface_skin': True})
    r -= d_ih

    for thick, kind in ((d_oc, 'ocean'), (d_iii, 'ice_III'),
                        (d_v, 'ice_V'), (d_vi, 'ice_VI')):
        if thick > min_km:
            layers.append({'name': LAYER_LABELS[kind], 'r_km': r,
                           'kind': kind})
        r -= thick

    # Silicate mantle: outer radius = hydrosphere base = R - D_hsphere.
    layers.append({'name': LAYER_LABELS['silicate'],
                   'r_km': R_km - d_hs, 'kind': 'silicate'})
    if r_core and r_core > 1.0:
        layers.append({'name': LAYER_LABELS['core'], 'r_km': r_core,
                       'kind': 'core'})
    return layers
