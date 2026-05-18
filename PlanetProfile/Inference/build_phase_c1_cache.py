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

    log.info(f"Body                : {bodyname}")
    log.info(f"Template module     : {template}")
    log.info(f"Tb_K bounds         : [{Tb_bounds[0]:.4f}, {Tb_bounds[1]:.4f}] K")
    log.info(f"Grid points         : {args.n_grid}")
    log.info(f"Tb_K grid           : {Tb_grid.tolist()}")
    log.info(f"Output cache path   : {out_path}")
    if ocean_overrides:
        log.info(f"Ocean overrides     : {ocean_overrides}")

    if out_path.exists() and not args.force:
        log.error(
            f"Output already exists: {out_path}. Pass --force to overwrite."
        )
        return 3

    if args.dry_run:
        log.info("Dry run — exiting without invoking PlanetProfile.")
        return 0

    # --- Build ---
    from PlanetProfile.Inference.cache_builder import build_tb_grid_cache

    t0 = time.time()
    build_tb_grid_cache(
        template, Tb_grid, str(out_path),
        progress=True, ocean_overrides=ocean_overrides,
    )
    dt = time.time() - t0
    log.info(f"Done in {dt:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
