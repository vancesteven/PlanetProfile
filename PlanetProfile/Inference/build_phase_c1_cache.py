"""
CLI driver for Phase C1 v2 structure-cache builds.

Reads the Tb_K bounds and ``structure_cache_path`` from a v2 JSON config
and dispatches a Tb-grid build via
:mod:`PlanetProfile.Inference.cache_builder`.

Usage
-----

::

    # Build Europa cache (default 9-point Tb grid, default body module
    # PlanetProfile.Default.Europa.PPEuropa)
    python -m PlanetProfile.Inference.build_phase_c1_cache \\
        --config PlanetProfile/Inference/configs/europa_seawater_andrade_7D.json

    # Override grid resolution
    python -m PlanetProfile.Inference.build_phase_c1_cache \\
        --config PlanetProfile/Inference/configs/ganymede_pureh2o_andrade_8D.json \\
        --n-grid 13

    # Use a custom planet-template module (rather than the
    # convention-derived default)
    python -m PlanetProfile.Inference.build_phase_c1_cache \\
        --config PlanetProfile/Inference/configs/callisto_mgso4_andrade_8D.json \\
        --template PlanetProfile.Default.Callisto.PPCallisto

The grid is uniformly spaced between ``param_space.Tb_K.bounds[0]`` and
``param_space.Tb_K.bounds[1]`` from the JSON. The output pickle is saved to
``structure_cache_path`` (also from the JSON), making the cache immediately
runnable through ``MCMCRunner(InferenceConfig.from_json(...))``.

Phase C1 Stage 3.
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path

import numpy as np


def _default_template_module(bodyname: str) -> str:
    """Convention: ``PlanetProfile.Default.<Body>.PP<Body>``."""
    return f"PlanetProfile.Default.{bodyname}.PP{bodyname}"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build a Tb-grid structure cache for a Phase C1 v2 "
            "BodyConfig JSON."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--config", required=True,
        help=(
            "Path to a v2 BodyConfig JSON (e.g. "
            "PlanetProfile/Inference/configs/europa_seawater_andrade_7D.json)"
        ),
    )
    parser.add_argument(
        "--n-grid", type=int, default=9,
        help="Number of Tb_K grid points (default: 9)",
    )
    parser.add_argument(
        "--n-wgrid", type=int, default=0,
        help=(
            "Number of log-spaced ocean-salinity grid points. When > 0, build "
            "a 2D (Tb × wOcean_ppt) cache (schema v3.0) using the config's "
            "'log10_wOcean_ppt' param_space bounds (0.1–100 ppt for [-1, 2]). "
            "Default 0 → the legacy 1D Tb-only cache."
        ),
    )
    parser.add_argument(
        "--template", default=None,
        help=(
            "PP planet-template module (e.g. "
            "PlanetProfile.Default.Europa.PPEuropa). Defaults to the "
            "convention 'PlanetProfile.Default.<bodyname>.PP<bodyname>' "
            "based on the JSON's bodyname."
        ),
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Overwrite the cache file if it already exists.",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help=(
            "Print the resolved Tb grid + template + output path and exit "
            "without running PP."
        ),
    )
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(name)s — %(message)s")
    log = logging.getLogger("PlanetProfile.Inference.build_phase_c1_cache")

    # --- Resolve config + bounds ---
    config_path = Path(args.config).resolve()
    if not config_path.is_file():
        log.error(f"Config not found: {config_path}")
        return 2
    with open(config_path) as f:
        cfg = json.load(f)

    try:
        Tb_bounds = cfg["param_space"]["Tb_K"]["bounds"]
    except KeyError:
        log.error(
            "Config does not sample Tb_K. This builder is intended for "
            "v2 configs with 'Tb_K' in param_space."
        )
        return 2
    if len(Tb_bounds) != 2 or Tb_bounds[1] <= Tb_bounds[0]:
        log.error(f"Bad Tb_K bounds in config: {Tb_bounds!r}")
        return 2

    bodyname = cfg.get("bodyname", "")
    if not bodyname:
        log.error("Config missing 'bodyname'.")
        return 2

    template = args.template or _default_template_module(bodyname)

    # Optional Ocean.* overrides (e.g. composition switch / concentration sweep)
    ocean_overrides = cfg.get("ocean_overrides") or None

    out_rel = cfg.get("structure_cache_path")
    if not out_rel:
        log.error("Config missing 'structure_cache_path'.")
        return 2
    # Resolve relative paths against the repo root (config files use repo-
    # relative paths so they're portable).
    repo_root = Path(__file__).resolve().parents[2]
    out_path = (repo_root / out_rel).resolve()

    Tb_grid = np.linspace(float(Tb_bounds[0]), float(Tb_bounds[1]), int(args.n_grid))

    # --- 2D salinity axis (v3): read log10_wOcean_ppt bounds if requested ---
    build_2d = args.n_wgrid and args.n_wgrid > 0
    w_grid = None
    if build_2d:
        try:
            w_bounds_log10 = cfg["param_space"]["log10_wOcean_ppt"]["bounds"]
        except KeyError:
            log.error(
                "--n-wgrid > 0 requires 'log10_wOcean_ppt' in param_space "
                "(v3 salinity config)."
            )
            return 2
        if len(w_bounds_log10) != 2 or w_bounds_log10[1] <= w_bounds_log10[0]:
            log.error(f"Bad log10_wOcean_ppt bounds: {w_bounds_log10!r}")
            return 2
        # Log-spaced salinity nodes (a linear grid in log10 w = geometric in w),
        # so the low-w physics is resolved (reviewer: Ae moves fastest at low w).
        w_grid = np.logspace(
            float(w_bounds_log10[0]), float(w_bounds_log10[1]), int(args.n_wgrid)
        )

    log.info(f"Body                : {bodyname}")
    log.info(f"Template module     : {template}")
    log.info(f"Tb_K bounds         : [{Tb_bounds[0]:.4f}, {Tb_bounds[1]:.4f}] K")
    log.info(f"Grid points         : {args.n_grid}")
    log.info(f"Tb_K grid           : {Tb_grid.tolist()}")
    if build_2d:
        log.info(f"wOcean_ppt bounds   : "
                 f"[{10**w_bounds_log10[0]:.4f}, {10**w_bounds_log10[1]:.4f}] ppt "
                 f"(log10 [{w_bounds_log10[0]}, {w_bounds_log10[1]}])")
        log.info(f"w grid points       : {args.n_wgrid}")
        log.info(f"wOcean_ppt grid     : {[round(w, 4) for w in w_grid.tolist()]}")
    log.info(f"Output cache path   : {out_path}")
    if ocean_overrides:
        log.info(f"Ocean overrides     : {ocean_overrides}")
        if build_2d:
            log.warning(
                "ocean_overrides is ignored for 2D builds — the salinity axis "
                "is set per-node from wOcean_ppt_grid."
            )

    if out_path.exists() and not args.force:
        log.error(
            f"Output already exists: {out_path}. Pass --force to overwrite."
        )
        return 3

    if args.dry_run:
        log.info("Dry run — exiting without invoking PlanetProfile.")
        return 0

    # --- Build ---
    t0 = time.time()
    if build_2d:
        from PlanetProfile.Inference.cache_builder import build_tbw_grid_cache
        build_tbw_grid_cache(
            template, Tb_grid, w_grid, str(out_path), progress=True,
        )
    else:
        from PlanetProfile.Inference.cache_builder import build_tb_grid_cache
        build_tb_grid_cache(
            template, Tb_grid, str(out_path),
            progress=True, ocean_overrides=ocean_overrides,
        )
    dt = time.time() - t0
    log.info(f"Done in {dt:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
