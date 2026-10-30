"""Crisp vector figure display with per-session export caching.

Streamlit reruns the whole script on every widget interaction; without
caching, each matplotlib figure pays artist creation + THREE savefig
exports (SVG inline + PDF + PNG downloads) per rerun — ~1 s per heavy
figure, ~13 s per rerun on the Inference results page. This module keys
the exported bytes on a caller-supplied token (typically id(result)):
reruns render the cached SVG string and serve cached download bytes
without touching matplotlib. Callers that pass a `builder` callable
also skip figure CONSTRUCTION on cache hits (use for expensive builds
like corner plots).

Pattern:
    # cheap build, cached exports (no closure needed):
    fig = ...build...
    display_vector_fig(fig, key='k2', download_label='k2', token=tok)

    # expensive build, cached exports AND skipped build:
    display_vector_fig(builder=lambda: build_corner(...),
                       key='corner', download_label='corner', token=tok)

Gotchas baked in (from the corner pilot): st.image needs the SVG as a
STRING (bytes go to PIL and raise); rasterize collections/images only,
never the whole figure, or the text blurs.
"""
from __future__ import annotations

import io
import re
import subprocess
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path

import streamlit as st

_CACHE_KEY = '_crisp_fig_cache'


@lru_cache(maxsize=1)
def _app_git_sha() -> str:
    """Return the app checkout SHA when the deployed snapshot includes git."""
    try:
        proc = subprocess.run(
            ['git', 'rev-parse', '--short', 'HEAD'],
            cwd=Path(__file__).resolve().parents[2],
            capture_output=True, text=True, timeout=1, check=False)
    except (OSError, subprocess.SubprocessError):
        return 'unavailable'
    return proc.stdout.strip() or 'unavailable'


def _slug(value: str, fallback: str) -> str:
    value = re.sub(r'[^A-Za-z0-9._-]+', '-', str(value)).strip('-._')
    return value or fallback


def _normalise_provenance(provenance=None) -> dict[str, str]:
    """Make backend-safe, complete export provenance from caller context."""
    provenance = provenance or {}
    artifact = str(provenance.get('artifact_filename') or 'not-applicable')
    slot_short = provenance.get('slot_short')
    if not slot_short:
        slot_short = Path(artifact).stem if artifact != 'not-applicable' else 'custom'
    return {
        'model_label': str(provenance.get('model_label') or 'Custom inference run'),
        'artifact_filename': artifact,
        'slot_short': _slug(slot_short, 'custom'),
        'conditioning': str(provenance.get('conditioning') or 'not recorded'),
        'app_git_sha': str(provenance.get('app_git_sha') or _app_git_sha()),
        'date': str(provenance.get('date') or
                    datetime.now(timezone.utc).date().isoformat()),
    }


def _metadata_text(figure_name: str, provenance: dict[str, str]):
    title = f"PlanetProfileApp {figure_name} — {provenance['model_label']}"
    subject = (
        f"Model slot: {provenance['model_label']}; "
        f"artifact: {provenance['artifact_filename']}; "
        f"conditioning: {provenance['conditioning']}; "
        f"app git SHA: {provenance['app_git_sha']}")
    return title, subject


def _filename_stem(figure_name: str, provenance: dict[str, str]) -> str:
    return (f"{_slug(figure_name, 'figure')}_"
            f"{provenance['slot_short']}_{provenance['date']}")


def rasterize_heavy_artists(fig):
    """Rasterize scatter/line-collections + images so SVG/PDF stay small
    while text/axes/legends/colorbars remain vector. Contour sets are
    skipped: matplotlib ignores set_rasterized on them (with a warning
    per call) — corner() contours stay vector, which is fine since their
    path count is modest."""
    from matplotlib.contour import ContourSet
    for ax in fig.get_axes():
        for coll in ax.collections:
            if isinstance(coll, ContourSet):
                continue
            coll.set_rasterized(True)
        for im in ax.images:
            im.set_rasterized(True)


def _cache() -> dict:
    return st.session_state.setdefault(_CACHE_KEY, {})


def _export(fig, *, figure_name='figure', provenance=None):
    # App figures use plain-text/mathtext labels with unicode (±, σ) and
    # raw underscores; PlanetProfile's config import flips the GLOBAL
    # text.usetex to True when TeX is installed (first wedge render), and
    # any app figure exported after that would be fed to real LaTeX and
    # crash (user 2026-07-24: Europa v4 "u panel unavailable"). Text is
    # rendered at savefig time, so pin usetex off around the exports;
    # PlanetProfile's own plot exports manage usetex themselves.
    import matplotlib.pyplot as plt
    provenance = _normalise_provenance(provenance)
    title, subject = _metadata_text(figure_name, provenance)
    pdf_metadata = {
        'Title': title,
        'Subject': subject,
        'Author': 'PlanetProfileApp',
        'Creator': 'PlanetProfileApp',
        'Keywords': 'PlanetProfile inference export',
    }
    svg_metadata = {
        'Title': title,
        'Description': subject,
        'Creator': 'PlanetProfileApp',
        'Keywords': ['PlanetProfile', 'inference', 'provenance'],
        'Date': provenance['date'],
    }
    png_metadata = {
        'Title': title,
        'Description': subject,
        'Software': 'PlanetProfileApp',
    }
    old_usetex = plt.rcParams.get('text.usetex', False)
    try:
        plt.rcParams['text.usetex'] = False
        fig.set_dpi(150)
        svg_buf = io.BytesIO()
        fig.savefig(svg_buf, format='svg', bbox_inches='tight',
                    metadata=svg_metadata)
        pdf_buf = io.BytesIO()
        fig.savefig(pdf_buf, format='pdf', bbox_inches='tight',
                    metadata=pdf_metadata)
        png_buf = io.BytesIO()
        fig.savefig(png_buf, format='png', dpi=200, bbox_inches='tight',
                    metadata=png_metadata)
    finally:
        plt.rcParams['text.usetex'] = old_usetex
    return (svg_buf.getvalue().decode('utf-8'), pdf_buf.getvalue(),
            png_buf.getvalue())


def _render(entry, *, key, download_label, provenance=None):
    svg, pdf, png = entry['svg'], entry['pdf'], entry['png']
    provenance = entry.get('provenance') or _normalise_provenance(provenance)
    filename_stem = entry.get('filename_stem') or \
        _filename_stem(key, provenance)
    st.image(svg, width='stretch')
    c1, c2 = st.columns(2)
    c1.download_button(
        f'Download {download_label} (PDF)', pdf,
        file_name=f'{filename_stem}.pdf', mime='application/pdf',
        icon=':material/download:', width='stretch', key=f'{key}_pdf')
    c2.download_button(
        f'Download {download_label} (PNG)', png,
        file_name=f'{filename_stem}.png', mime='image/png',
        icon=':material/download:', width='stretch', key=f'{key}_png')


def display_vector_fig(fig=None, *, key, download_label='plot',
                       token=None, builder=None, heavy=False,
                       close=True, provenance=None):
    """Render a figure as crisp inline SVG + PDF/PNG downloads, caching
    the exported bytes per (key, token).

    fig: an already-built matplotlib figure (exports cached; build cost
        still paid by the caller each rerun).
    builder: zero-arg callable returning a figure — called ONLY on cache
        miss, so both build and export are skipped on reruns. Provide
        exactly one of fig/builder.
    token: cache token; None disables caching (always rebuild/export).
    heavy: rasterize collections/images before export (builder path;
        for the fig path call rasterize_heavy_artists yourself, as
        before).
    close: plt.close() the figure after export (builder path only).
    provenance: model/conditioning context embedded in PDF/SVG/PNG metadata
        and used for the provenance-aware download filename.
    """
    cache = _cache()
    provenance = _normalise_provenance(provenance)
    entry = cache.get(key)
    if (token is not None and entry is not None and entry['token'] == token
            and entry.get('provenance') == provenance):
        _render(entry, key=key, download_label=download_label)
        return

    if fig is None:
        # Text objects capture usetex at creation — pin it off during the
        # build so a PlanetProfile import earlier in the rerun (which
        # flips the global rcParam when TeX is installed) can't route
        # app-figure labels through real LaTeX.
        import matplotlib.pyplot as plt
        old_usetex = plt.rcParams.get('text.usetex', False)
        try:
            plt.rcParams['text.usetex'] = False
            fig = builder()
        finally:
            plt.rcParams['text.usetex'] = old_usetex
        if heavy:
            rasterize_heavy_artists(fig)
    svg, pdf, png = _export(fig, figure_name=key, provenance=provenance)
    if builder is not None and close:
        import matplotlib.pyplot as plt
        plt.close(fig)
    entry = {
        'token': token, 'svg': svg, 'pdf': pdf, 'png': png,
        'provenance': provenance,
        'filename_stem': _filename_stem(key, provenance),
    }
    if token is not None:
        cache[key] = entry
    _render(entry, key=key, download_label=download_label)
