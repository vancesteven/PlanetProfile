"""
Inference cache utilities for the Inference GUI.

Enables instant result loading by caching InferenceResult objects.
Provides three-tier caching:
1. Session state (fastest, current session only)
2. Disk persistence (cross-session reuse)
3. Recomputation (if cache misses)

Pattern adapted from explore_cache.py
"""
import os
import logging
import pickle
import json
import hashlib
from typing import Optional, Dict, Any

log = logging.getLogger('PlanetProfile')


def generate_inference_cache_key(
    bodyname: str,
    param_space: Dict[str, Dict[str, Any]],
    observables: Dict[str, tuple],
    sampler_type: str = 'mcmc',
    sampler_settings: Optional[Dict] = None,
    random_state: int = 42
) -> str:
    """
    Generate a unique cache key from computational parameters.

    Visualization-only parameters (plot settings, colors) are excluded
    to avoid unnecessary recomputation.

    Args:
        bodyname: Planet name (e.g. 'Titan')
        param_space: Dict of {param_name: {prior_type, bounds/mean/std}}
        observables: Dict of {observable_name: (value, uncertainty)}
        sampler_type: 'mcmc' or 'sbi'
        sampler_settings: Sampler hyperparameters (n_effective, etc.)
        random_state: Random seed for reproducibility

    Returns:
        String cache key suitable for filenames
    """
    # Build deterministic representation
    key_data = {
        'bodyname': bodyname,
        'param_space': sorted(param_space.items()),  # Sort for consistency
        'observables': sorted(observables.items()),
        'sampler_type': sampler_type,
        'sampler_settings': sampler_settings or {},
        'random_state': random_state
    }

    # Hash to fixed-length string
    key_str = json.dumps(key_data, sort_keys=True)
    key_hash = hashlib.md5(key_str.encode()).hexdigest()[:16]

    # Human-readable prefix
    param_names = '_'.join(sorted(param_space.keys())[:3])  # First 3 params
    prefix = f"{bodyname}_{sampler_type}_{param_names}"

    # Sanitize for filesystem
    safe_prefix = prefix.replace('/', '_').replace(' ', '_')

    return f"{safe_prefix}_{key_hash}"


def get_cache_path(cache_dir: str, cache_key: str) -> str:
    """
    Resolve cache file path.

    Args:
        cache_dir: Base cache directory
        cache_key: Cache key from generate_inference_cache_key()

    Returns:
        Full path to cache file (.pkl extension)
    """
    os.makedirs(cache_dir, exist_ok=True)
    return os.path.join(cache_dir, f"{cache_key}.pkl")


def save_inference_to_cache(
    result: Any,
    cache_dir: str,
    cache_key: str
) -> None:
    """
    Save InferenceResult to disk cache.

    Args:
        result: InferenceResult object to cache
        cache_dir: Base cache directory
        cache_key: Cache key from generate_inference_cache_key()
    """
    cache_path = get_cache_path(cache_dir, cache_key)

    try:
        with open(cache_path, 'wb') as f:
            pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)
        log.info(f"Saved inference result to cache: {cache_path}")
    except Exception as e:
        log.warning(f"Failed to save inference cache: {e}")


def load_inference_from_cache(
    cache_dir: str,
    cache_key: str
) -> Optional[Any]:
    """
    Load InferenceResult from disk cache.

    Args:
        cache_dir: Base cache directory
        cache_key: Cache key from generate_inference_cache_key()

    Returns:
        InferenceResult object if found, None otherwise
    """
    cache_path = get_cache_path(cache_dir, cache_key)

    if not os.path.exists(cache_path):
        return None

    try:
        with open(cache_path, 'rb') as f:
            result = pickle.load(f)
        log.info(f"Loaded inference result from cache: {cache_path}")
        return result
    except Exception as e:
        log.warning(f"Failed to load inference cache: {e}")
        return None


def list_cached_inferences(cache_dir: str) -> list:
    """
    List all cached inference runs in directory.

    Args:
        cache_dir: Base cache directory

    Returns:
        List of tuples (cache_key, filepath, file_size_mb, modified_time)
    """
    if not os.path.exists(cache_dir):
        return []

    cached_runs = []
    for filename in os.listdir(cache_dir):
        if filename.endswith('.pkl'):
            filepath = os.path.join(cache_dir, filename)
            cache_key = filename[:-4]  # Remove .pkl extension
            file_size = os.path.getsize(filepath) / (1024 * 1024)  # MB
            modified_time = os.path.getmtime(filepath)

            cached_runs.append((cache_key, filepath, file_size, modified_time))

    # Sort by modified time (newest first)
    cached_runs.sort(key=lambda x: x[3], reverse=True)

    return cached_runs


def delete_cache(cache_dir: str, cache_key: Optional[str] = None) -> None:
    """
    Delete cached inference results.

    Args:
        cache_dir: Base cache directory
        cache_key: Specific cache key to delete, or None to delete all
    """
    if cache_key is not None:
        # Delete specific cache
        cache_path = get_cache_path(cache_dir, cache_key)
        if os.path.exists(cache_path):
            os.remove(cache_path)
            log.info(f"Deleted inference cache: {cache_path}")
    else:
        # Delete all caches in directory
        if os.path.exists(cache_dir):
            for filename in os.listdir(cache_dir):
                if filename.endswith('.pkl'):
                    filepath = os.path.join(cache_dir, filename)
                    os.remove(filepath)
            log.info(f"Deleted all inference caches in: {cache_dir}")


def get_cache_metadata(cache_dir: str, cache_key: str) -> Optional[Dict]:
    """
    Extract metadata from cached inference result without loading full object.

    Args:
        cache_dir: Base cache directory
        cache_key: Cache key from generate_inference_cache_key()

    Returns:
        Dict with metadata (bodyname, n_samples, param_names, etc.) or None
    """
    cache_path = get_cache_path(cache_dir, cache_key)

    if not os.path.exists(cache_path):
        return None

    try:
        with open(cache_path, 'rb') as f:
            result = pickle.load(f)

        # Extract lightweight metadata
        metadata = {
            'bodyname': result.config.bodyname if hasattr(result, 'config') else 'Unknown',
            'mode': result.config.mode if hasattr(result, 'config') else 'Unknown',
            'n_samples': len(result.samples) if hasattr(result, 'samples') else 0,
            'param_names': result.param_names if hasattr(result, 'param_names') else [],
            'file_size_mb': os.path.getsize(cache_path) / (1024 * 1024),
            'modified_time': os.path.getmtime(cache_path)
        }

        return metadata
    except Exception as e:
        log.warning(f"Failed to extract cache metadata: {e}")
        return None
