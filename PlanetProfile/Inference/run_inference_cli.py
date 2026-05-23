#!/usr/bin/env python3
"""
Command-line interface for running inference on HPC or in subprocess.

This script is the entry point for headless inference execution. It:
1. Loads InferenceConfig from JSON file
2. Runs inference (MCMC or SBI)
3. Saves InferenceResult to pickle file
4. Optionally writes progress updates to JSON file for GUI monitoring

Usage:
    # Basic usage
    python -m PlanetProfile.Inference.run_inference_cli \\
        --config inference_config.json \\
        --output results.pkl

    # With progress tracking (for GUI)
    python -m PlanetProfile.Inference.run_inference_cli \\
        --config inference_config.json \\
        --output results.pkl \\
        --progress progress.json

    # HPC batch submission example
    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --time=4:00:00
    #SBATCH --mem=16GB

    module load python/3.11
    python -m PlanetProfile.Inference.run_inference_cli \\
        --config $SLURM_SUBMIT_DIR/config.json \\
        --output $SLURM_SUBMIT_DIR/result_$SLURM_JOB_ID.pkl

Author: PlanetProfile Team
Date: 2026-04-29
"""
import argparse
import json
import logging
import os
import sys
import time
import tempfile
from pathlib import Path
from typing import Dict, Optional

# Ensure numba has a writable cache dir BEFORE TidalPy gets imported
# anywhere downstream. TidalPy.rheology.partial_melt.melting_models uses
# @njit(cacheable=True), which raises "no locator available" when numba
# can't establish a cache location — seen under sandboxed runtimes where
# the default cache path is unwritable. Set this only if not already set.
if not os.environ.get('NUMBA_CACHE_DIR'):
    _default_numba_cache = os.path.join(tempfile.gettempdir(), 'pp_numba_cache')
    try:
        os.makedirs(_default_numba_cache, exist_ok=True)
        os.environ['NUMBA_CACHE_DIR'] = _default_numba_cache
    except OSError:
        # If creation fails, leave NUMBA_CACHE_DIR unset so numba uses its
        # own default and surfaces any error itself.
        pass

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
log = logging.getLogger('PlanetProfile.Inference')


def setup_progress_callback(progress_file: Optional[str], update_interval: int = 50):
    """
    Create progress callback that writes to JSON file.

    Args:
        progress_file: Path to JSON file for progress updates
        update_interval: Update file every N iterations

    Returns:
        Callback function or None
    """
    if progress_file is None:
        return None

    progress_path = Path(progress_file)
    progress_path.parent.mkdir(parents=True, exist_ok=True)

    last_update = {'iteration': -update_interval}  # Force first update

    def progress_callback(data: Dict):
        """Write progress data to JSON file."""
        current_iter = data.get('iteration', 0)

        # Only update every N iterations
        if current_iter - last_update['iteration'] >= update_interval:
            try:
                with open(progress_path, 'w') as f:
                    json.dump(data, f, indent=2)
                last_update['iteration'] = current_iter
            except Exception as e:
                log.warning(f"Failed to write progress file: {e}")

    return progress_callback


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Run Bayesian inference for PlanetProfile (MCMC/SBI)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run MCMC with default settings
  python -m PlanetProfile.Inference.run_inference_cli \\
      --config titan_inference.json \\
      --output titan_result.pkl

  # Run with progress tracking (for GUI)
  python -m PlanetProfile.Inference.run_inference_cli \\
      --config config.json \\
      --output result.pkl \\
      --progress progress.json \\
      --update-interval 100

  # Validate config without running
  python -m PlanetProfile.Inference.run_inference_cli \\
      --config config.json \\
      --validate-only
        """
    )

    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Path to InferenceConfig JSON file'
    )
    parser.add_argument(
        '--output',
        type=str,
        help='Path to save InferenceResult pickle file (required unless --validate-only)'
    )
    parser.add_argument(
        '--progress',
        type=str,
        default=None,
        help='Path to JSON file for progress updates (optional)'
    )
    parser.add_argument(
        '--update-interval',
        type=int,
        default=50,
        help='Write progress every N iterations (default: 50)'
    )
    parser.add_argument(
        '--validate-only',
        action='store_true',
        help='Validate config and exit without running inference'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    # Set log level
    if args.verbose:
        logging.getLogger('PlanetProfile').setLevel(logging.DEBUG)

    # Validate arguments
    if not args.validate_only and args.output is None:
        parser.error("--output is required unless --validate-only is set")

    # Load configuration
    log.info(f"Loading inference configuration from {args.config}")
    try:
        from PlanetProfile.Inference.inference_core import InferenceConfig, validate_config

        config = InferenceConfig.from_json(args.config)
        log.info(f"  Mode: {config.mode}")
        log.info(f"  Body: {config.bodyname}")
        log.info(f"  Parameters: {list(config.param_space.keys())}")
        log.info(f"  Observables: {list(config.observables.keys())}")

    except Exception as e:
        log.error(f"Failed to load config: {e}")
        return 1

    # Validate configuration
    log.info("Validating configuration...")
    is_valid, errors = validate_config(config)

    if not is_valid:
        log.error("Configuration validation failed:")
        for err in errors:
            log.error(f"  - {err}")
        return 1

    log.info("Configuration valid")

    if args.validate_only:
        log.info("Validation complete (--validate-only), exiting")
        return 0

    # Setup progress callback
    progress_callback = setup_progress_callback(args.progress, args.update_interval)

    # Run inference
    log.info("=" * 70)
    log.info(f"Starting {config.mode.upper()} inference")
    log.info("=" * 70)

    t0 = time.time()

    try:
        from PlanetProfile.Inference.inference_core import run_inference

        result = run_inference(config, progress_callback=progress_callback)

        elapsed = time.time() - t0
        log.info("=" * 70)
        log.info(f"Inference completed in {elapsed/60:.1f} minutes")
        log.info("=" * 70)

        # Log summary statistics
        log.info("Posterior summary:")
        summary = result.get_summary_stats()
        for i, param_name in enumerate(result.param_names):
            log.info(f"  {param_name}:")
            log.info(f"    median = {summary['median'][i]:.4g}")
            log.info(f"    mean   = {summary['mean'][i]:.4g}")
            log.info(f"    std    = {summary['std'][i]:.4g}")

        # Log convergence metrics
        if result.convergence_metrics:
            log.info("Convergence metrics:")
            for key, val in result.convergence_metrics.items():
                if isinstance(val, (int, float)):
                    log.info(f"  {key}: {val:.4g}")
                else:
                    log.info(f"  {key}: {val}")

        # Save result
        log.info(f"Saving result to {args.output}")
        result.save(args.output)

        file_size_mb = Path(args.output).stat().st_size / (1024 * 1024)
        log.info(f"Result saved ({file_size_mb:.1f} MB)")

        # Write final progress update
        if args.progress:
            final_progress = {
                'iteration': -1,  # Sentinel for "completed"
                'status': 'completed',
                'n_samples': len(result.samples),
                'elapsed_seconds': elapsed,
                'convergence_metrics': result.convergence_metrics
            }
            with open(args.progress, 'w') as f:
                json.dump(final_progress, f, indent=2)

        return 0

    except KeyboardInterrupt:
        log.warning("Inference interrupted by user (Ctrl+C)")
        if args.progress:
            with open(args.progress, 'w') as f:
                json.dump({'status': 'interrupted'}, f)
        return 130  # Standard exit code for SIGINT

    except Exception as e:
        log.error(f"Inference failed: {e}", exc_info=True)
        if args.progress:
            with open(args.progress, 'w') as f:
                json.dump({'status': 'failed', 'error': str(e)}, f)
        return 1


if __name__ == '__main__':
    sys.exit(main())
