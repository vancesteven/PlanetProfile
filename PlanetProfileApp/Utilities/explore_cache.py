"""
Exploration cache utilities for the Exploreogram GUI.

Enables instant z-variable switching by caching the full ExplorationResultsStruct
(which contains ALL extracted fields, not just the current z-variable).
Also provides disk persistence for cross-session reuse and future MC/NN analyses.
"""
import os
import json
import logging
import pickle

log = logging.getLogger('PlanetProfile')


def generate_cache_key(bodyname, xName, xRange, yName, yRange, nx, ny,
                       skip_induction=True, ocean_comp='default',
                       exc_names=None):
    """
    Generate a unique cache key from computational parameters.

    z-variable is deliberately EXCLUDED — changing it is a visualization-only
    operation that should not trigger recomputation.

    Args:
        bodyname: Planet name (e.g. 'Europa')
        xName: X-axis exploration parameter name
        xRange: [min, max] for X parameter
        yName: Y-axis exploration parameter name
        yRange: [min, max] for Y parameter
        nx, ny: Grid resolution
        skip_induction: Whether induction was skipped
        ocean_comp: Ocean composition identifier
        exc_names: List of excitation frequency names (e.g. ['orbital', 'synodic', 'synodic 2nd'])

    Returns:
        String cache key suitable for filenames
    """
    key = (
        f"{bodyname}_{xName}_{xRange[0]:.4g}-{xRange[1]:.4g}_"
        f"{yName}_{yRange[0]:.4g}-{yRange[1]:.4g}_{int(nx)}x{int(ny)}_"
        f"induct{int(not skip_induction)}"
    )
    # Include excitation selection so changing frequencies triggers recompute
    if exc_names:
        exc_str = '+'.join(sorted(exc_names))
        key += f"_{exc_str}"
    # Sanitize for filesystem
    return key.replace('/', '_').replace(' ', '_')


def get_cache_path(cache_dir, cache_key):
    """Return the full path to a cache file for the given key."""
    return os.path.join(cache_dir, f"{cache_key}.pkl")


def save_to_cache(Exploration, cache_dir, cache_key, meta=None):
    """
    Save an ExplorationResultsStruct to disk cache.

    Uses the same pickle format as PlanetProfile's WriteResults().

    meta: optional JSON-serializable dict describing the run (widget
    parameters: bodyname, xName/yName/zName, ranges, nx/ny, drivers,
    label). Written as a <key>.meta.json sidecar so browsable libraries
    (e.g. the public demo loader) can seed the GUI to match this cache
    entry without parsing the key string.
    """
    os.makedirs(cache_dir, exist_ok=True)
    path = get_cache_path(cache_dir, cache_key)
    with open(path, 'wb') as f:
        pickle.dump(Exploration, f)
    if meta is not None:
        with open(os.path.join(cache_dir, f"{cache_key}.meta.json"), 'w') as f:
            json.dump(meta, f, indent=1)
    log.info(f"Exploration cached to {path}")
    return path


def load_cache_meta(cache_dir, cache_key):
    """Load the .meta.json sidecar for a cache entry (None if absent)."""
    path = os.path.join(cache_dir, f"{cache_key}.meta.json")
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except Exception as e:
        log.warning(f"Cache meta unreadable ({path}): {e}")
        return None


def load_from_cache(cache_dir, cache_key):
    """
    Load an ExplorationResultsStruct from disk cache.

    Returns None if cache file does not exist or is corrupt.
    """
    path = get_cache_path(cache_dir, cache_key)
    if not os.path.isfile(path):
        return None
    try:
        with open(path, 'rb') as f:
            Exploration = pickle.load(f)
        log.info(f"Exploration loaded from cache: {path}")
        return Exploration
    except Exception as e:
        log.warning(f"Cache file corrupt or incompatible ({path}): {e}")
        return None


def list_cached_explorations(cache_dir):
    """
    List all cached exploration keys in the cache directory.

    Returns list of (cache_key, file_size_mb, mtime) tuples.
    """
    if not os.path.isdir(cache_dir):
        return []
    entries = []
    for fname in os.listdir(cache_dir):
        if fname.endswith('.pkl'):
            fpath = os.path.join(cache_dir, fname)
            size_mb = os.path.getsize(fpath) / (1024 * 1024)
            mtime = os.path.getmtime(fpath)
            key = fname[:-4]  # strip .pkl
            entries.append((key, size_mb, mtime))
    entries.sort(key=lambda x: x[2], reverse=True)  # newest first
    return entries
