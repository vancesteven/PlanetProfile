"""
CLI producer for CMR2 discretization-offset sidecar files.

Writes ``<cache_stem>_offsets.json`` next to a Tb-grid structure cache,
recording, per grid point, the offset between the cache's native
PlanetProfile CMR2 and the value reconstructed by
:func:`PlanetProfile.Inference.structure_derivation.compute_cmr2` for a
no-core (R_core = 0) uniform-silicate reconstruction:

    CMR2_offset = CMR2_pp_native (cached) - recon_nocore

``MCMCRunner._derive_cmr2_via_mass_conservation`` consumes the sidecar
automatically (np.interp in Tb) whenever the file exists next to the
config's ``structure_cache_path``.

VALIDITY GUARD (scientific, do not remove): the definition above isolates
PlanetProfile's radial-integration *discretization* error ONLY when the
cached silicate interior has (near-)uniform density, as in the Test52
Titan grid (built from a uniform-rho_sil template). For a cache whose
silicate interior is radially structured (e.g. the pre-existing Callisto
C2_andrade grids, rho_sil spanning ~1300-3200 kg/m3), native - recon is
dominated by the PHYSICAL structured-vs-uniform difference (measured
-0.035 = -8.4 sigma_obs for Callisto) and is NOT a discretization offset;
applying it would corrupt the CMR2 likelihood. This tool therefore refuses
to write a sidecar when the silicate density is non-uniform or the offset
is implausibly large, and directs the user to rebuild the grid from a
uniform-silicate template first (opus design review, 2026-07-10).

Usage::

    python -m PlanetProfile.Inference.build_cmr2_offset_sidecar \\
        --config PlanetProfile/Inference/configs/test52_titan_noocean_andrade_10D.json \\
        [--cache path/to/grid.pkl] [--out path/to/sidecar.json] [--force]

Exit codes: 0 written; 1 setup/IO error; 3 validity guard refused.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import pickle
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

log = logging.getLogger('PlanetProfile')

# Guard thresholds (see module docstring). A genuine discretization offset
# is ~1e-3 in CMR2 (Test52: 9.5e-4); anything beyond MAX_ABS_OFFSET signals
# a physical artifact. Silicate-density spread beyond RHO_SIL_SPAN_MAX
# means the uniform-silicate reconstruction is not comparing like with like.
MAX_ABS_OFFSET = 0.005
RHO_SIL_SPAN_MAX_KGM3 = 1.0

DEFINITION = ('CMR2_offset = CMR2_pp_native (cached) - recon_nocore '
              '(structure_derivation.compute_cmr2 on R_core=0 reconstruction)')


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'],
            cwd=Path(__file__).resolve().parent,
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return ''


def compute_offsets(grid: dict):
    """Per-grid-point (offset, diagnostics) for a Test50-format Tb grid.

    Returns (Tb_K_grid list, offsets list, diagnostics list). Raises
    ValueError with a scientific explanation when the validity guard trips.
    """
    from .structure_derivation import (
        compute_cmr2,
        derive_silicate_density,
        extract_hydrosphere_layers,
    )

    if 'Tb_K_grid' not in grid or 'structures' not in grid:
        raise KeyError(
            "Cache is not a Test50-format Tb grid ({'Tb_K_grid', 'structures'})")

    tb_grid = [float(t) for t in np.asarray(grid['Tb_K_grid'])]
    offsets = []
    diags = []
    for tb, struct in zip(tb_grid, grid['structures']):
        native = float(struct['CMR2'])
        R_body_m = float(struct['R_body_m'])
        M_total_kg = float(struct['Mtot_kg'])
        hydro_layers, R_oceanbot_m, M_hydro_kg = extract_hydrosphere_layers(struct)
        if not hydro_layers:
            raise ValueError(f'Tb={tb}: empty hydrosphere in cache structure')

        # Silicate-uniformity guard: cells strictly below the hydrosphere
        # bottom must share (near-)one density for native-recon to isolate
        # discretization error.
        r_m = np.asarray(struct['r_m'], float)
        rho = np.asarray(struct['rho'], float)
        sil_mask = r_m <= R_oceanbot_m + 1e-6
        rho_sil_cells = rho[sil_mask]
        span = (float(np.nanmax(rho_sil_cells) - np.nanmin(rho_sil_cells))
                if rho_sil_cells.size else 0.0)
        if span > RHO_SIL_SPAN_MAX_KGM3:
            raise ValueError(
                f'Tb={tb}: cached silicate density is radially structured '
                f'(span {span:.1f} kg/m3 > {RHO_SIL_SPAN_MAX_KGM3}). '
                f'native - recon would measure the physical structured-vs-'
                f'uniform silicate difference, not discretization error. '
                f'Rebuild the grid from a uniform-silicate template before '
                f'generating a sidecar.')

        # R_core=0 -> rho_core irrelevant; wide bounds so derivation never
        # rejects (mirrors the sidecar-definition reconstruction, not the
        # likelihood's config-bounded call).
        rho_sil, _ = derive_silicate_density(
            M_total_kg=M_total_kg, M_hydrosphere_kg=M_hydro_kg,
            R_oceanbot_m=R_oceanbot_m, R_core_m=0.0, rho_core_kgm3=0.0,
            bounds=(0.0, 1.0e9))
        recon = compute_cmr2(
            [(0.0, R_oceanbot_m, rho_sil)] + list(hydro_layers),
            R_body_m=R_body_m, M_body_kg=M_total_kg)
        offset = native - recon
        if abs(offset) > MAX_ABS_OFFSET:
            raise ValueError(
                f'Tb={tb}: |offset| = {abs(offset):.4g} > {MAX_ABS_OFFSET} '
                f'(a genuine discretization offset is ~1e-3). This cache '
                f'cannot yield a valid sidecar; rebuild from a uniform-'
                f'silicate template.')
        offsets.append(float(offset))
        diags.append({'Tb_K': tb, 'native': native, 'recon_nocore': float(recon),
                      'rho_sil_span_kgm3': span})
    return tb_grid, offsets, diags


def _sha256(path, chunk=1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


# Hydrosphere-match thresholds for measurement mode (opus review 2026-07-10):
# the measured offset is the hydrosphere cells' discretization term, so it
# transfers to the production cache only if both grids share (essentially)
# the same hydrosphere. Refuse above REFUSE; record measured values always.
HYDRO_MATCH_REFUSE_FRAC = 1e-4


def check_hydrosphere_match(measure_grid: dict, production_grid: dict):
    """Per-Tb hydrosphere comparison between measurement and production grids.

    Returns a dict with the measured max fractional deviations; raises
    ValueError when grids are not comparable or deviate beyond the refuse
    threshold.
    """
    from .structure_derivation import extract_hydrosphere_layers

    tb_m = np.asarray(measure_grid['Tb_K_grid'], float)
    tb_p = np.asarray(production_grid['Tb_K_grid'], float)
    if tb_m.shape != tb_p.shape or not np.allclose(tb_m, tb_p, rtol=0, atol=1e-9):
        raise ValueError(
            f'Tb grids differ between measurement ({tb_m}) and production '
            f'({tb_p}) caches — rebuild the measurement grid on the '
            f'production Tb grid.')

    max_dr = 0.0
    max_drho = 0.0
    n_cells_equal = True
    for sm, sp in zip(measure_grid['structures'], production_grid['structures']):
        hm, _, _ = extract_hydrosphere_layers(sm)
        hp, _, _ = extract_hydrosphere_layers(sp)
        if len(hm) != len(hp):
            n_cells_equal = False
            continue
        am, ap = np.asarray(hm, float), np.asarray(hp, float)
        r_scale = np.maximum(np.abs(ap[:, 1]), 1.0)
        rho_scale = np.maximum(np.abs(ap[:, 2]), 1.0)
        max_dr = max(max_dr, float(np.max(np.abs(am[:, :2] - ap[:, :2])
                                          / r_scale[:, None])))
        max_drho = max(max_drho, float(np.max(np.abs(am[:, 2] - ap[:, 2])
                                              / rho_scale)))
    result = {'n_cells_equal': bool(n_cells_equal),
              'max_frac_dr': max_dr, 'max_frac_drho': max_drho}
    if not n_cells_equal or max_dr > HYDRO_MATCH_REFUSE_FRAC \
            or max_drho > HYDRO_MATCH_REFUSE_FRAC:
        raise ValueError(
            f'Hydrosphere mismatch between measurement and production grids '
            f'({result}); transferability of the discretization offset is '
            f'unproven — refusing to write a sidecar.')
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description='Build a CMR2 discretization-offset sidecar for a Tb-grid cache.')
    parser.add_argument('--config', help='InferenceConfig JSON (cache path read from it).')
    parser.add_argument('--cache', help='Explicit cache .pkl (overrides --config).')
    parser.add_argument('--out', help='Explicit sidecar output path '
                                      '(default: <cache_stem>_offsets.json next to cache).')
    parser.add_argument('--force', action='store_true',
                        help='Overwrite an existing sidecar.')
    parser.add_argument('--production-cache',
                        help='MEASUREMENT MODE: --cache is a uniform-silicate '
                             'metrology grid; the sidecar is written for THIS '
                             'production cache (default output: its '
                             '<stem>_offsets.json) after a per-Tb hydrosphere-'
                             'match check between the two grids.')
    args = parser.parse_args(argv)

    from .structure_cache import resolve_cache_path

    if args.cache:
        cache_path = resolve_cache_path(args.cache)
    elif args.config:
        with open(args.config) as f:
            cfg = json.load(f)
        cache_path = resolve_cache_path(cfg['structure_cache_path'])
    else:
        parser.error('one of --config or --cache is required')

    if not Path(cache_path).exists():
        print(f'ERROR: cache not found: {cache_path}', file=sys.stderr)
        return 1

    production_path = None
    if args.production_cache:
        production_path = resolve_cache_path(args.production_cache)
        if not Path(production_path).exists():
            print(f'ERROR: production cache not found: {production_path}',
                  file=sys.stderr)
            return 1

    # The sidecar sits next to (and is named for) the cache it corrects:
    # the production cache in measurement mode, else the measured cache.
    sidecar_target = Path(production_path) if production_path else Path(cache_path)
    out_path = Path(args.out) if args.out else (
        sidecar_target.with_name(sidecar_target.stem + '_offsets.json'))
    if out_path.exists() and not args.force:
        print(f'ERROR: sidecar exists: {out_path} (use --force to overwrite)',
              file=sys.stderr)
        return 1

    with open(cache_path, 'rb') as f:
        grid = pickle.load(f)

    measurement_meta = None
    if production_path:
        with open(production_path, 'rb') as f:
            production_grid = pickle.load(f)
        try:
            hydro_match = check_hydrosphere_match(grid, production_grid)
        except ValueError as e:
            print(f'REFUSED: {e}', file=sys.stderr)
            return 3
        measurement_meta = {
            'measurement_mode': 'uniform_silicate_measurement_grid',
            'measured_on_cache': str(cache_path),
            'measured_on_cache_sha256': _sha256(cache_path),
            'production_cache': str(production_path),
            'production_cache_sha256': _sha256(production_path),
            'hydrosphere_match': hydro_match,
            'transferability_note': (
                'Offset is the hydrosphere cells\' discretization term, '
                'measured on a uniform-silicate metrology grid (where '
                'native - reconstruction isolates integration-scheme error) '
                'and applied to the production cache\'s reconstruction path; '
                'validity rests on the hydrosphere_match evidence above '
                '(opus review 2026-07-10).'),
        }

    try:
        tb_grid, offsets, diags = compute_offsets(grid)
    except (KeyError, ValueError) as e:
        print(f'REFUSED: {e}', file=sys.stderr)
        return 3

    sidecar = {
        'cache_path': str(sidecar_target),
        'Tb_K_grid': tb_grid,
        'CMR2_offset_per_point': offsets,
        'offset_mean': float(np.mean(offsets)),
        'offset_std': float(np.std(offsets)),
        'offset_min': float(np.min(offsets)),
        'offset_max': float(np.max(offsets)),
        'definition': DEFINITION,
        'git_sha': _git_sha(),
        'build_timestamp_utc': datetime.now(timezone.utc).isoformat().replace(
            '+00:00', 'Z'),
        'note': ('Built by PlanetProfile.Inference.build_cmr2_offset_sidecar; '
                 'consumed by MCMCRunner._derive_cmr2_via_mass_conservation '
                 '(np.interp in Tb) whenever this file sits next to the cache.'),
        'diagnostics': diags,
    }
    if measurement_meta:
        sidecar.update(measurement_meta)
    with open(out_path, 'w') as f:
        json.dump(sidecar, f, indent=2)
    print(f'WROTE {out_path} (mean {sidecar["offset_mean"]:.6g} '
          f'+/- {sidecar["offset_std"]:.2g})')
    return 0


if __name__ == '__main__':
    sys.exit(main())
