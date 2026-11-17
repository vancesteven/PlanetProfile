"""Pure (Streamlit-free) helpers behind the Inference page's MCMC config loader.

The page keeps the FULL loaded config dict as the single source of truth and
edits only a working copy of the user-editable blocks; ``build_run_config``
folds those edits back into the loaded dict so every non-edited block
(fixed_params, derived_params, param_groups, metadata, induction_bounds,
gravity_forward_model, arrhenius_params, ocean_overrides, bodyname, ...)
reaches the runner unchanged.
"""

from __future__ import annotations

import copy
import functools
import json
import math
import pickle
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from .inference_core import InferenceConfig
from .structure_cache import resolve_cache_path

CONFIGS_DIR = Path(__file__).resolve().parent / "configs"
REPO_ROOT = Path(__file__).resolve().parents[2]

GUI_STATUSES = ("production", "sensitivity", "smoke", "candidate", "retired")
PRODUCTION_STATUS = "production"

SAMPLER_EDIT_FIELDS = ("n_effective", "n_reeval", "random_state",
                       "checkpoint_interval", "max_iterations")
SAMPLER_FORM_DEFAULTS = {
    "n_effective": 500, "random_state": 42, "n_reeval": 500,
    "checkpoint_interval": 100, "max_iterations": 10000,
}

# prior_type -> the numeric fields the form edits for it. Semantics are the
# runner's (mcmc_runner._build_prior): log-uniform and uniform-in-linear
# bounds are log10 EXPONENTS; truncated_gaussian bounds are in linear units.
PRIOR_FIELDS = {
    "uniform": ("bounds",),
    "log-uniform": ("bounds",),
    "uniform-in-linear": ("bounds",),
    "normal": ("mean", "std"),
    "truncated_gaussian": ("mean", "std", "bounds"),
}

# Cache axis -> (axis key in the pickle, sampled-parameter name, transform
# from axis units to the parameter's units).
_GRID_AXES = (
    ("Tb_K_grid", "Tb_K", lambda v: v),
    ("wOcean_ppt_grid", "log10_wOcean_ppt", math.log10),
    ("zb_km_grid", "zb_km", lambda v: v),
)

BUILD_COMMAND = ("python -m PlanetProfile.Inference.build_phase_c1_cache "
                 "--config {config} --n-grid <N_TB>{wflag}")


# --------------------------------------------------------------------------
# Config catalogue / status
# --------------------------------------------------------------------------

def gui_info(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """``metadata.gui`` as {status, note, runnable}; unlabeled -> candidate."""
    gui = (cfg.get("metadata") or {}).get("gui") or {}
    status = gui.get("status")
    if status not in GUI_STATUSES:
        status = "candidate"
        note = gui.get("note") or "No metadata.gui status recorded."
    else:
        note = gui.get("note", "")
    return {"status": status, "note": note,
            "runnable": gui.get("runnable", True) is not False}


def config_catalogue(configs_dir: Optional[Path] = None) -> List[Dict[str, Any]]:
    """One entry per config JSON: path, stem, body, status, note, runnable."""
    entries = []
    for fp in sorted(Path(configs_dir or CONFIGS_DIR).glob("*.json")):
        try:
            data = json.loads(fp.read_text(encoding="utf-8"))
        except Exception:
            continue
        body = (data.get("bodyname") or data.get("body")
                or data.get("planet_name") or data.get("Body") or "Unknown")
        entries.append({"path": str(fp), "stem": fp.stem, "body": body,
                        **gui_info(data)})
    return entries


# --------------------------------------------------------------------------
# Run-config assembly
# --------------------------------------------------------------------------

def _pair(v) -> tuple:
    return (float(v[0]), float(v[1]))


def merge_observables(loaded_obs: Dict[str, Any],
                      form_obs: Dict[str, Any]) -> Dict[str, Any]:
    """Loaded observables with the form's edits applied.

    Entries whose (value, sigma) match the loaded ones keep the loaded
    object verbatim (so the loaded key spelling, e.g. ``Im_k2``, and list
    type are preserved); edited ones keep the loaded key; keys the form
    removed are dropped; keys the form added are appended as [value, sigma].
    """
    out: Dict[str, Any] = {}
    for key, val in form_obs.items():
        orig = key
        if key not in loaded_obs and key == "abs_Im_k2" and "Im_k2" in loaded_obs:
            orig = "Im_k2"
        if orig in loaded_obs:
            try:
                same = _pair(loaded_obs[orig]) == _pair(val)
            except (TypeError, ValueError, IndexError):
                same = False
            out[orig] = loaded_obs[orig] if same else [val[0], val[1]]
        else:
            out[key] = [val[0], val[1]]
    # preserve the loaded ordering
    ordered = {k: out[k] for k in loaded_obs if k in out}
    ordered.update({k: v for k, v in out.items() if k not in ordered})
    return ordered


def build_run_config(loaded: Optional[Dict[str, Any]],
                     *,
                     param_space: Dict[str, Any],
                     observables: Dict[str, Any],
                     sampler_settings: Dict[str, Any],
                     structure_cache_path: Optional[str],
                     selected_params: Optional[Sequence[str]] = None,
                     fallback_bodyname: str = "Custom") -> InferenceConfig:
    """InferenceConfig for a run: the loaded dict plus the form's edits.

    With no edits the result equals ``InferenceConfig.from_dict(loaded)``.
    With ``loaded=None`` (no config file loaded) a bare config is built from
    the form fields alone (legacy preset/custom path).
    """
    if selected_params is not None:
        ps = {k: copy.deepcopy(param_space[k]) for k in selected_params
              if k in param_space}
    else:
        ps = copy.deepcopy(param_space)

    if loaded is None:
        return InferenceConfig(
            mode="mcmc", bodyname=fallback_bodyname, param_space=ps,
            observables=observables, sampler_settings=sampler_settings,
            structure_cache_path=structure_cache_path,
            random_state=sampler_settings["random_state"])

    cfg = copy.deepcopy(loaded)
    cfg["param_space"] = ps
    cfg["observables"] = merge_observables(loaded.get("observables", {}),
                                           observables)

    ss = dict(cfg.get("sampler_settings") or {})
    loaded_ss = loaded.get("sampler_settings") or {}
    for field in SAMPLER_EDIT_FIELDS:
        if field not in sampler_settings:
            continue
        if field == "random_state":
            base = loaded_ss.get("random_state", loaded.get("random_state",
                                 SAMPLER_FORM_DEFAULTS["random_state"]))
        else:
            base = loaded_ss.get(field, SAMPLER_FORM_DEFAULTS[field])
        if sampler_settings[field] != base:
            ss[field] = sampler_settings[field]
            if field == "random_state":
                cfg["random_state"] = sampler_settings[field]
    cfg["sampler_settings"] = ss

    if structure_cache_path and structure_cache_path != loaded.get(
            "structure_cache_path"):
        cfg["structure_cache_path"] = structure_cache_path
    return InferenceConfig.from_dict(cfg)


def run_config_from_state(state, fallback_bodyname: str = "Custom"
                          ) -> InferenceConfig:
    """``build_run_config`` fed from the page's session-state keys."""
    return build_run_config(
        state.get("inference_loaded_config"),
        param_space=state["inference_param_space"],
        observables=state["inference_observables"],
        sampler_settings=state["inference_sampler_settings"],
        structure_cache_path=state["inference_structure_cache_path"],
        selected_params=state["inference_selected_params"],
        fallback_bodyname=fallback_bodyname)


# --------------------------------------------------------------------------
# Prior editing helpers
# --------------------------------------------------------------------------

def default_prior_entry(param_def, prior_type: str,
                        from_entry: Optional[Dict[str, Any]] = None
                        ) -> Dict[str, Any]:
    """Fresh prior dict of ``prior_type`` for a newly added/retyped parameter."""
    entry: Dict[str, Any] = {"prior_type": prior_type}
    fields = PRIOR_FIELDS[prior_type]
    bounds = list(from_entry["bounds"]) if (from_entry and "bounds" in from_entry
                                            ) else list(param_def.default_bounds)
    if prior_type == "log-uniform" and not (from_entry and from_entry.get(
            "prior_type") in ("log-uniform", "uniform-in-linear")):
        if bounds[0] > 0 and bounds[1] > 0:
            bounds = [math.log10(bounds[0]), math.log10(bounds[1])]
    if "bounds" in fields:
        entry["bounds"] = [float(bounds[0]), float(bounds[1])]
    if "mean" in fields:
        mean = param_def.default_mean
        if mean is None:
            mean = 0.5 * (bounds[0] + bounds[1])
        std = param_def.default_std
        if not std:
            std = 0.25 * abs(bounds[1] - bounds[0]) or 1.0
        entry["mean"] = float(mean)
        entry["std"] = float(std)
    return entry


# --------------------------------------------------------------------------
# Structure-grid check
# --------------------------------------------------------------------------

@functools.lru_cache(maxsize=8)
def _read_axes(path: str, mtime: float, size: int) -> Dict[str, Any]:
    with open(path, "rb") as fh:
        data = pickle.load(fh)
    axes = {}
    if isinstance(data, dict):
        for axis_key, param, fn in _GRID_AXES:
            grid = data.get(axis_key)
            if grid is None or len(grid) == 0:
                continue
            vals = [float(v) for v in grid]
            lo, hi = min(vals), max(vals)
            if axis_key == "zb_km_grid" and data.get("frozen_zb_km_grid") is not None \
                    and len(data["frozen_zb_km_grid"]):
                hi = max(hi, float(max(data["frozen_zb_km_grid"])))
            axes[param] = {"axis": axis_key, "range": (fn(lo), fn(hi)),
                           "n": len(vals)}
    return axes


def cache_axes(path: Path) -> Dict[str, Any]:
    st = path.stat()
    return _read_axes(str(path), st.st_mtime, st.st_size)


def _prior_extent(spec: Dict[str, Any]) -> Optional[tuple]:
    ptype = spec.get("prior_type")
    if ptype in ("uniform", "uniform-in-linear", "truncated_gaussian"):
        b = spec.get("bounds")
        return (float(b[0]), float(b[1])) if b else None
    if ptype == "normal":
        m, s = spec.get("mean"), spec.get("std")
        return (m - 3 * s, m + 3 * s) if m is not None and s else None
    return None


def check_structure_cache(param_space: Dict[str, Any],
                          fixed_params: Optional[Dict[str, Any]],
                          cache_path: Optional[str],
                          config_path: Optional[str] = None
                          ) -> Dict[str, Any]:
    """Real structure-grid validation.

    Returns ``{'status': 'unset'|'missing'|'unreadable'|'ok', 'resolved',
    'axes', 'issues', 'command'}``. ``issues`` lists each sampled/fixed
    parameter whose prior (or fixed value) extends outside the cached axis,
    as ``{'param', 'prior_range', 'cache_range'}``.
    """
    if not cache_path:
        return {"status": "unset", "resolved": None, "axes": {}, "issues": [],
                "command": None}
    resolved = resolve_cache_path(cache_path)
    if not resolved.is_absolute() and not resolved.exists():
        resolved = REPO_ROOT / resolved
    needs_w = "log10_wOcean_ppt" in (param_space or {})
    cfg_ref = config_path or "<config.json>"
    command = BUILD_COMMAND.format(
        config=cfg_ref, wflag=" --n-wgrid <N_W>" if needs_w else "")
    out = {"status": "missing", "resolved": str(resolved), "axes": {},
           "issues": [], "command": command}
    if not resolved.exists():
        return out
    try:
        axes = cache_axes(resolved)
    except Exception as exc:  # corrupt / unreadable pickle
        out.update(status="unreadable", error=str(exc))
        return out
    out["status"] = "ok"
    out["axes"] = axes
    extents = {k: _prior_extent(v) for k, v in (param_space or {}).items()}
    for k, v in (fixed_params or {}).items():
        try:
            extents[k] = (float(v), float(v))
        except (TypeError, ValueError):
            pass
    for param, info in axes.items():
        ext = extents.get(param)
        if ext is None:
            continue
        lo, hi = info["range"]
        tol = 1e-6 * max(abs(hi - lo), 1e-12)
        if ext[0] < lo - tol or ext[1] > hi + tol:
            out["issues"].append({"param": param, "prior_range": ext,
                                  "cache_range": (lo, hi)})
    return out
