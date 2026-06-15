"""
mcmc_run_loader.py — Utility for scanning and loading saved MCMC result files.

Provides a lightweight index of available *_result.pkl files and a
Streamlit-cached loader for InferenceResult objects.  Intended for use by
the "Browse runs" mode in PlanetProfileApp.

No heavy scientific imports at module level — InferenceResult is imported
inside load_mcmc_result() to avoid circular-import issues during app startup.
"""
import os
import pickle
import logging
import streamlit as st
from pathlib import Path
from typing import Optional

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default directories to scan (relative to repo root).
# Users can extend this list or pass an explicit results_root to
# scan_mcmc_results() to include additional locations.
# ---------------------------------------------------------------------------
DEFAULT_RESULTS_ROOTS: list[str] = [
    "PlanetProfile/Test/mcmc_results",
    # Add more default scan roots here, e.g.:
    # "PlanetProfile/Inference/results",
]


def _derive_bodyname(pkl_path: Path) -> str:
    """
    Heuristically extract the body name from the file path.

    Strategy (in order):
      1. If a parent directory name looks like a known body (capitalised word),
         use it.
      2. Strip common result suffixes from the filename stem.
      3. Fall back to the immediate parent directory name.
    """
    _KNOWN_BODIES = {
        "europa", "titan", "ganymede", "callisto", "enceladus",
        "io", "triton", "miranda", "ariel", "umbriel", "titania", "oberon",
        "pluto", "ceres", "dione", "rhea",
    }
    for part in pkl_path.parts:
        if part.lower() in _KNOWN_BODIES:
            return part.capitalize()
    # Fallback: stem minus trailing _result / _results
    stem = pkl_path.stem
    for suffix in ("_result", "_results"):
        if stem.lower().endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return stem


def _derive_name(pkl_path: Path, bodyname: str) -> str:
    """
    Build a human-readable run name from the file path.

    Format:  "<BodyName> / <run_dir>"  where run_dir is the immediate parent
    of the pkl.  If the parent directory IS the body directory, use the file
    stem instead.
    """
    parent = pkl_path.parent
    # If parent folder matches the body name (case-insensitive), use file stem
    if parent.name.lower() == bodyname.lower():
        return f"{bodyname} / {pkl_path.stem}"
    return f"{bodyname} / {parent.name}"


def scan_mcmc_results(results_root: str) -> list[dict]:
    """
    Recursively scan *results_root* for ``*_result.pkl`` files and return a
    structured index — newest files first.

    Each entry contains:
      - ``path``     : absolute path string to the pkl
      - ``name``     : human-readable label, e.g. ``"Europa / T25_seawater"``
      - ``bodyname`` : body name derived from directory/filename heuristic
      - ``mtime``    : last-modified timestamp (float, seconds since epoch)
      - ``size_mb``  : file size in megabytes (float)

    The pkl is NOT loaded; only ``os.stat`` metadata is collected.

    Parameters
    ----------
    results_root:
        Directory to search.  May be absolute or relative; relative paths are
        resolved from the current working directory.

    Returns
    -------
    list[dict]
        Run records sorted by ``mtime`` descending (newest first).
        Returns an empty list if the directory does not exist.
    """
    root = Path(results_root).resolve()
    if not root.is_dir():
        log.warning("scan_mcmc_results: directory not found: %s", root)
        return []

    records: list[dict] = []
    for pkl_path in root.rglob("*_result.pkl"):
        try:
            stat = pkl_path.stat()
        except OSError as exc:
            log.warning("scan_mcmc_results: cannot stat %s: %s", pkl_path, exc)
            continue

        bodyname = _derive_bodyname(pkl_path)
        records.append(
            {
                "path": str(pkl_path),
                "name": _derive_name(pkl_path, bodyname),
                "bodyname": bodyname,
                "mtime": stat.st_mtime,
                "size_mb": stat.st_size / (1024 * 1024),
            }
        )

    records.sort(key=lambda r: r["mtime"], reverse=True)
    return records


@st.cache_data(show_spinner=False)
def load_mcmc_result(path: str, mtime: float):  # -> Optional[InferenceResult]
    """
    Load an ``InferenceResult`` from *path* with Streamlit result caching.

    The *mtime* parameter is included in the cache key so that the cache is
    automatically invalidated whenever the file on disk changes.

    Parameters
    ----------
    path:
        Absolute path to the ``*_result.pkl`` file.
    mtime:
        File modification time (``os.stat().st_mtime``).  Pass the value
        from the corresponding ``scan_mcmc_results`` record.

    Returns
    -------
    InferenceResult or None
        The loaded result object, or ``None`` if the file is missing or
        cannot be unpickled (a warning is logged in that case).
    """
    # Deferred import to avoid circular imports at app startup
    from PlanetProfile.Inference.inference_core import InferenceResult  # noqa: F401

    pkl = Path(path)
    if not pkl.is_file():
        log.warning("load_mcmc_result: file not found: %s", path)
        return None
    try:
        with open(pkl, "rb") as fh:
            result = pickle.load(fh)
        return result
    except Exception as exc:  # broad catch — corrupt pkl, version mismatch, etc.
        log.warning("load_mcmc_result: failed to load %s: %s", path, exc)
        return None


def get_run_display_label(run_record: dict) -> str:
    """
    Return a short display string suitable for ``st.selectbox``.

    Example output::

        "Europa — T25_seawater  (2026-06-13, 45.2 MB)"

    Parameters
    ----------
    run_record:
        A dict as returned by :func:`scan_mcmc_results`.
    """
    import datetime

    date_str = datetime.datetime.fromtimestamp(run_record["mtime"]).strftime("%Y-%m-%d")
    size_str = f"{run_record['size_mb']:.1f} MB"
    return f"{run_record['name']}  ({date_str}, {size_str})"
