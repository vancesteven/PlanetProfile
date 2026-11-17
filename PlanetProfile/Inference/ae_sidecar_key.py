"""Stdlib-only helpers defining the induction-Ae sidecar key.

Shared by ``MCMCRunner._precompute_ae_grid`` (which writes/reads the sidecar)
and ``plans/scripts/extract_sbi_cache_paths.py`` (which validates, at deploy
time, that a shipped sidecar's key equals what the slot's runner computes).
Kept free of numpy/scipy so the deploy script can load it by file path.
"""
from __future__ import annotations

from typing import Any, Dict, Iterable, Optional, Set


def parse_bind_channel(name: str):
    """Parse ``Bind_<label>_<comp>_<part>`` into (label, comp, part) or None."""
    if not name.startswith('Bind_'):
        return None
    body = name[len('Bind_'):]
    for part in ('real', 'imag'):
        suffix = '_' + part
        if not body.endswith(suffix):
            continue
        rest = body[:-len(suffix)]
        for comp in ('x', 'y', 'z'):
            comp_suffix = '_' + comp
            if rest.endswith(comp_suffix):
                label = rest[:-len(comp_suffix)]
                if label:
                    return label, comp, part
    return None


def has_induction_channels(observables: Iterable[str],
                           bound_labels: Iterable[str]) -> bool:
    obs = [k for k in observables
           if k.startswith(('Ae_', 'BiAmp_', 'BiPhase_', 'Bind_'))]
    return bool(obs) or bool(list(bound_labels))


def requested_ae_labels(observables: Iterable[str],
                        bound_labels: Iterable[str]) -> Set[str]:
    """Excitation labels the Ae precompute must cover."""
    requested: Set[str] = set(bound_labels)
    for k in observables:
        if k.startswith('Ae_') and k.endswith('_real'):
            requested.add(k[len('Ae_'):-len('_real')])
        elif k.startswith('Ae_') and k.endswith('_imag'):
            requested.add(k[len('Ae_'):-len('_imag')])
        elif k.startswith('BiAmp_'):
            requested.add(k[len('BiAmp_'):])
        elif k.startswith('BiPhase_') and k.endswith('_deg'):
            requested.add(k[len('BiPhase_'):-len('_deg')])
        else:
            parsed = parse_bind_channel(k)
            if parsed is not None:
                requested.add(parsed[0])
    return requested


def sidecar_key(labels: Iterable[str], n: int, cache_size: int) -> Dict[str, Any]:
    return {'labels': tuple(sorted(labels)), 'n': n, 'cache_size': cache_size}


def is_stale_same_labels(stored: Any, ours: Dict[str, Any]) -> bool:
    """True when ``stored`` differs from ``ours`` ONLY in cache identity
    (``n`` / ``cache_size``) with identical labels: a sidecar left over from
    a cache that has since been rebuilt."""
    if not isinstance(stored, dict) or stored == ours:
        return False
    if set(stored) != set(ours):
        return False
    if stored.get('labels') != ours.get('labels'):
        return False
    return True
