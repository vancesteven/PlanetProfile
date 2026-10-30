#!/usr/bin/env python3
"""
Audit HP-ice Rayleigh numbers (Ra*) in cached MCMC structures.

Detects which structure caches built with KALOUSOVA_CONVECTION=False
might have triggered the buggy Deschamps & Sotin 2001 HP-ice dispatcher.

Usage:
    python audit_hp_ice_ra.py [--dir DIR] [--format json]
"""

import argparse
import json
import pickle
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# Constants
G_SI = 6.67430e-11  # gravitational constant, m^3 kg^-1 s^-2
ALPHA_ICE = 1.0e-4  # thermal expansivity, K^-1
KAPPA_ICE = 1.0e-6  # thermal diffusivity, m^2/s
RA_STAR_OK_THRESHOLD = 1e4
RA_STAR_AMBIGUOUS_THRESHOLD = 1e6


def compute_gravity_profile(r_m: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """
    Compute local gravity g(r) from cumulative mass below each radius.

    REQUIRES r_m to be ascending (center -> surface). Caller must normalize.
    Returns gravity in the same ordering as the (ascending) input.
    """
    # Vectorized cumulative-mass: dM_j = 4 pi r_mid^2 rho_j dr_j
    r_mid = 0.5 * (r_m[:-1] + r_m[1:])
    dr = np.diff(r_m)
    dM = 4.0 * np.pi * r_mid**2 * rho[:-1] * dr
    M_below = np.concatenate(([0.0], np.cumsum(dM)))

    g = np.zeros_like(r_m, dtype=float)
    nonzero = r_m > 0
    g[nonzero] = G_SI * M_below[nonzero] / (r_m[nonzero] ** 2)
    return g


def find_contiguous_regions(phases: np.ndarray, phase_id: int) -> List[Tuple[int, int]]:
    """
    Find contiguous runs of cells with a given phase.

    Args:
        phases: phase array
        phase_id: target phase (3, 5, or 6)

    Returns:
        List of (start_index, end_index) tuples (inclusive)
    """
    mask = phases == phase_id
    indices = np.where(mask)[0]

    if len(indices) == 0:
        return []

    regions = []
    start = indices[0]
    for i in range(1, len(indices)):
        if indices[i] != indices[i - 1] + 1:
            regions.append((start, indices[i - 1]))
            start = indices[i]
    regions.append((start, indices[-1]))

    return regions


def compute_ra_star(
    r_m: np.ndarray,
    rho: np.ndarray,
    T_K: np.ndarray,
    eta_Pa_base: np.ndarray,
    i_top: int,
    i_bot: int,
) -> float:
    """
    Compute modified Rayleigh number Ra* for an HP-ice layer.

    Args:
        r_m: radius array (m)
        rho: density array (kg/m^3)
        T_K: temperature array (K)
        eta_Pa_base: base viscosity array (Pa·s)
        i_top: index of shallower (upper) cell
        i_bot: index of deeper (lower) cell

    Returns:
        Ra* (dimensionless)
    """
    with np.errstate(all='ignore'):
        # Extract layer parameters
        r_top = r_m[i_top]
        r_bot = r_m[i_bot]
        D = abs(r_top - r_bot)

        T_top = T_K[i_top]
        T_m = T_K[i_bot]
        Delta_T = T_m - T_top

        # Mean values over the layer
        layer_indices = np.arange(min(i_top, i_bot), max(i_top, i_bot) + 1)
        rho_mean = np.mean(rho[layer_indices])
        eta_m = np.mean(eta_Pa_base[layer_indices])

        # Local gravity at layer midpoint
        r_mid = 0.5 * (r_top + r_bot)
        g_profile = compute_gravity_profile(r_m, rho)
        # Interpolate g to midpoint
        idx_mid = np.searchsorted(r_m, r_mid)
        if idx_mid == 0:
            g_layer = g_profile[0]
        elif idx_mid >= len(g_profile):
            g_layer = g_profile[-1]
        else:
            w = (r_mid - r_m[idx_mid - 1]) / (r_m[idx_mid] - r_m[idx_mid - 1])
            g_layer = (1 - w) * g_profile[idx_mid - 1] + w * g_profile[idx_mid]

        # Ra* = rho * g * alpha * Delta_T * D^3 / (kappa * eta_m)
        if eta_m > 0:
            ra_star = (rho_mean * g_layer * ALPHA_ICE * Delta_T * D**3) / (KAPPA_ICE * eta_m)
        else:
            ra_star = np.nan

    return float(ra_star)


def audit_structure(structure: Dict[str, Any], cache_name: str) -> Optional[Dict[str, Any]]:
    """
    Audit a single structure dict for HP-ice Ra* values.

    Returns:
        Dict with phase results or None on error/skip
    """
    try:
        # Validate required fields
        required_fields = ['r_m', 'rho', 'T_K', 'phases', 'eta_Pa_base', 'Tb_K']
        for field in required_fields:
            if field not in structure:
                print(f"[SKIP-PARTIAL] {cache_name}: missing field '{field}'")
                return None

        r_m = np.array(structure['r_m'], dtype=float)
        rho = np.array(structure['rho'], dtype=float)
        T_K = np.array(structure['T_K'], dtype=float)
        phases = np.array(structure['phases'], dtype=int)
        eta_Pa_base = np.array(structure['eta_Pa_base'], dtype=float)
        Tb_K = float(structure['Tb_K'])

        # Normalize to ascending r_m (center -> surface) so all downstream
        # logic (gravity profile, np.searchsorted, top/bot orientation) is
        # consistent regardless of how this cache was written.
        if r_m.size > 1 and r_m[0] > r_m[-1]:
            r_m = r_m[::-1]
            rho = rho[::-1]
            T_K = T_K[::-1]
            phases = phases[::-1]
            eta_Pa_base = eta_Pa_base[::-1]

        results = {
            'cache_name': cache_name,
            'Tb_K': Tb_K,
            'phase_results': [],
        }

        # Audit each HP-ice phase
        for phase_id, phase_name in [(3, 'III'), (5, 'V'), (6, 'VI')]:
            regions = find_contiguous_regions(phases, phase_id)

            for start_idx, end_idx in regions:
                if end_idx - start_idx < 1:  # Skip single-cell regions
                    continue

                # Determine shallower (upper) and deeper (lower) indices
                if r_m[start_idx] < r_m[end_idx]:
                    # r_m increasing: end_idx is shallower
                    i_top = end_idx
                    i_bot = start_idx
                else:
                    # r_m decreasing: start_idx is shallower
                    i_top = start_idx
                    i_bot = end_idx

                D_km = abs(r_m[i_top] - r_m[i_bot]) / 1000.0
                T_top = T_K[i_top]
                T_m = T_K[i_bot]
                eta_m = np.mean(eta_Pa_base[min(i_top, i_bot):max(i_top, i_bot) + 1])

                ra_star = compute_ra_star(r_m, rho, T_K, eta_Pa_base, i_top, i_bot)

                # Determine verdict
                if np.isnan(ra_star) or np.isinf(ra_star):
                    verdict = 'ERROR'
                elif ra_star < RA_STAR_OK_THRESHOLD:
                    verdict = 'OK'
                elif ra_star < RA_STAR_AMBIGUOUS_THRESHOLD:
                    verdict = 'AMBIGUOUS'
                else:
                    verdict = 'SUSPECT'

                results['phase_results'].append({
                    'phase': phase_name,
                    'D_km': D_km,
                    'T_top_K': T_top,
                    'T_m_K': T_m,
                    'eta_m': eta_m,
                    'Ra_star': ra_star,
                    'verdict': verdict,
                })

        return results

    except Exception as e:
        print(f"[ERROR] {cache_name}: {type(e).__name__}: {e}")
        return None


def load_cache(pkl_path: Path) -> Optional[Dict[str, Any]]:
    """Load a pickle cache file."""
    try:
        with open(pkl_path, 'rb') as f:
            return pickle.load(f)
    except Exception as e:
        print(f"[ERROR] {pkl_path.name}: Failed to load pickle: {e}")
        return None


def audit_directory(cache_dir: Path) -> Tuple[List[Dict], int, int, int, int]:
    """
    Audit all MCMC structure caches in a directory.

    Returns:
        (results_list, total_scanned, ok_count, ambiguous_count, suspect_count, error_count)
    """
    # Find cache files matching pattern
    import re
    pattern = re.compile(r'^(titan|callisto|europa|ganymede).*\.pkl$')
    exclude_pattern = re.compile(r'_results\.pkl$|_mcmc\.pkl$|_corner\.png$')

    cache_files = sorted([
        f for f in cache_dir.rglob('*.pkl')
        if pattern.match(f.name) and not exclude_pattern.search(f.name)
    ])

    all_results = []
    total_scanned = 0
    ok_caches = 0
    ambiguous_caches = 0
    suspect_caches = 0
    error_caches = 0

    for cache_file in cache_files:
        cache_data = load_cache(cache_file)
        if cache_data is None:
            error_caches += 1
            continue

        # Determine cache format
        structures = []
        if 'Tb_K_grid' in cache_data and 'structures' in cache_data:
            # Multi-structure grid format
            structures = cache_data['structures']
        elif all(k in cache_data for k in ['r_m', 'phases', 'T_K']):
            # Single structure format
            structures = [cache_data]
        else:
            print(f"[SKIP] {cache_file.name} — unrecognized cache format "
                  f"(top-level keys: {list(cache_data.keys())})")
            error_caches += 1
            continue

        # Audit each structure
        for struct_idx, structure in enumerate(structures):
            total_scanned += 1
            struct_name = f"{cache_file.name}[{struct_idx}]" if len(structures) > 1 else cache_file.name

            result = audit_structure(structure, struct_name)
            if result is None:
                continue

            all_results.append(result)

            # Tally verdicts for this cache
            has_ok = any(r['verdict'] == 'OK' for r in result['phase_results'])
            has_ambiguous = any(r['verdict'] == 'AMBIGUOUS' for r in result['phase_results'])
            has_suspect = any(r['verdict'] == 'SUSPECT' for r in result['phase_results'])
            has_error = any(r['verdict'] == 'ERROR' for r in result['phase_results'])

            if has_error:
                pass  # Don't count in ok/ambiguous/suspect yet
            elif has_suspect:
                suspect_caches += 1
            elif has_ambiguous:
                ambiguous_caches += 1
            elif has_ok or len(result['phase_results']) == 0:
                ok_caches += 1

    return all_results, total_scanned, ok_caches, ambiguous_caches, suspect_caches


def print_table(results: List[Dict]) -> None:
    """Print audit results as a formatted table."""
    print("\n" + "=" * 120)
    print(f"{'Cache':<48} {'Tb_K':>8} {'Phase':>6} {'D_km':>8} {'T_top':>7} {'T_m':>7} "
          f"{'eta_m':>12} {'Ra*':>12} {'Verdict':<12}")
    print("=" * 120)

    for result in results:
        cache_name = result['cache_name'][:47]
        Tb_K = result['Tb_K']

        if not result['phase_results']:
            print(f"{cache_name:<48} {Tb_K:>8.2f}  (no HP-ice phases)")
            continue

        for idx, phase_result in enumerate(result['phase_results']):
            if idx == 0:
                print(f"{cache_name:<48} {Tb_K:>8.2f}  {phase_result['phase']:>6}  "
                      f"{phase_result['D_km']:>8.1f}  {phase_result['T_top_K']:>7.1f}  "
                      f"{phase_result['T_m_K']:>7.1f}  {phase_result['eta_m']:>12.2e}  "
                      f"{phase_result['Ra_star']:>12.2e}  {phase_result['verdict']:<12}")
            else:
                print(f"{'':48}  {'':8}  {phase_result['phase']:>6}  "
                      f"{phase_result['D_km']:>8.1f}  {phase_result['T_top_K']:>7.1f}  "
                      f"{phase_result['T_m_K']:>7.1f}  {phase_result['eta_m']:>12.2e}  "
                      f"{phase_result['Ra_star']:>12.2e}  {phase_result['verdict']:<12}")

    print("=" * 120 + "\n")


def print_summary(total: int, ok: int, ambiguous: int, suspect: int) -> None:
    """Print audit summary."""
    print("\n" + "=" * 60)
    print("AUDIT SUMMARY")
    print("=" * 60)
    print(f"Total caches scanned:           {total}")
    print(f"Caches with all-OK verdicts:    {ok}")
    print(f"Caches with ambiguous results:  {ambiguous}")
    print(f"Caches with suspect results:    {suspect}")
    print("=" * 60 + "\n")


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Audit HP-ice Rayleigh numbers in MCMC structure caches.'
    )
    parser.add_argument(
        '--dir',
        type=Path,
        default=Path('/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai/'
                    'PlanetProfile/Test/mcmc_results/'),
        help='Directory containing .pkl structure caches'
    )
    parser.add_argument(
        '--format',
        choices=['text', 'json'],
        default='text',
        help='Output format'
    )

    args = parser.parse_args()
    cache_dir = args.dir

    if not cache_dir.is_dir():
        print(f"ERROR: Cache directory not found: {cache_dir}", file=sys.stderr)
        return 1

    print(f"Auditing HP-ice Ra* in: {cache_dir}\n")

    results, total, ok, ambiguous, suspect = audit_directory(cache_dir)

    if args.format == 'text':
        print_table(results)
        print_summary(total, ok, ambiguous, suspect)
    elif args.format == 'json':
        summary = {
            'total_scanned': total,
            'ok_caches': ok,
            'ambiguous_caches': ambiguous,
            'suspect_caches': suspect,
            'results': results,
        }
        json_path = cache_dir / 'audit_ra_star.json'
        with open(json_path, 'w') as f:
            json.dump(summary, f, indent=2)
        print(f"JSON summary written to: {json_path}\n")
        print_summary(total, ok, ambiguous, suspect)

    return 0


if __name__ == '__main__':
    sys.exit(main())
