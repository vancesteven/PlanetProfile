"""
Helper CLI to generate both clathrate and no-clathrate structure cache variants.

Usage:
    python -m PlanetProfile.Inference.prepare_structure_variants \
        --test-module PlanetProfile.Test.PPTest41 \
        --output-dir cache/

Generates:
    - cache/{bodyname}_structure_clath.pkl
    - cache/{bodyname}_structure_noclath.pkl

Author: PlanetProfile Team
Date: 2026-04-29
"""
import argparse
import logging
import sys
from pathlib import Path

log = logging.getLogger('PlanetProfile')


def main():
    """Generate both clathrate variants of structure cache."""
    parser = argparse.ArgumentParser(
        description='Generate clathrate and no-clathrate structure cache variants',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Titan variants (PPTest41 + PPTest41NoClath)
  python -m PlanetProfile.Inference.prepare_structure_variants \\
      --test-module PlanetProfile.Test.PPTest41 \\
      --output-dir titan_cache/

  # Europa variants (PPTest3 + PPTest3NoClath)
  python -m PlanetProfile.Inference.prepare_structure_variants \\
      --test-module PlanetProfile.Test.PPTest3 \\
      --output-dir europa_cache/
        """
    )
    parser.add_argument('--test-module', required=True,
                       help='PPTest module with clathrates (e.g., PlanetProfile.Test.PPTest41)')
    parser.add_argument('--output-dir', default='.',
                       help='Output directory for cache files (default: current directory)')
    parser.add_argument('--rheology', default='andrade',
                       choices=['andrade', 'maxwell'],
                       help='Rheology model for TidalPy (default: andrade)')
    parser.add_argument('--force', action='store_true',
                       help='Force rebuild even if cache exists')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose logging')
    parser.add_argument('--maxwell', action='store_true',
                       help='Build Tb_K grid cache for Maxwell rheology instead of clathrate variants')
    parser.add_argument('--tb-min', type=float, default=251.2,
                       help='Minimum Tb_K for Maxwell grid (default: 251.2 K, above triple point)')
    parser.add_argument('--tb-max', type=float, default=270.0,
                       help='Maximum Tb_K for Maxwell grid (default: 270.0 K; '
                            'Titan clathrate models fail above ~269 K)')
    parser.add_argument('--tb-step', type=float, default=2.0,
                       help='Tb_K step size in K (default: 2.0)')

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    from .structure_cache import build_structure_from_pptest, save_structure_cache

    # Infer no-clathrate module name
    test_module = args.test_module
    if 'PPTest41' in test_module:
        test_module_noclath = test_module.replace('PPTest41', 'PPTest41NoClath')
        bodyname = 'Titan'
    elif 'PPTest3' in test_module:
        test_module_noclath = test_module.replace('PPTest3', 'PPTest3NoClath')
        bodyname = 'Europa'
    elif 'PPTest42' in test_module:
        test_module_noclath = test_module.replace('PPTest42', 'PPTest42NoClath')
        bodyname = 'Titan'
    elif 'PPTest43' in test_module:
        test_module_noclath = test_module.replace('PPTest43', 'PPTest43NoClath')
        bodyname = 'Titan'
    elif 'PPTest44' in test_module:
        test_module_noclath = test_module.replace('PPTest44', 'PPTest44NoClath')
        bodyname = 'Titan'
    else:
        log.error(f"Unknown test module: {test_module}")
        log.error("Supported: PPTest3, PPTest41, PPTest42, PPTest43, PPTest44")
        sys.exit(1)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    bodyname_lower = bodyname.lower()

    # Maxwell grid mode: build Tb_K-keyed structure dict for grid cache
    if args.maxwell:
        from .structure_cache import build_structure_grid
        import numpy as np
        tb_grid = list(np.arange(args.tb_min, args.tb_max + args.tb_step / 2, args.tb_step))
        grid_cache_path = output_dir / f'{bodyname_lower}_maxwell_grid_cache.pkl'
        log.info("=" * 70)
        log.info(f"Building Maxwell Tb_K grid cache: {len(tb_grid)} points "
                 f"[{tb_grid[0]:.1f} – {tb_grid[-1]:.1f} K]")
        log.info(f"Test module: {test_module}")
        log.info("=" * 70)
        try:
            build_structure_grid(
                test_module, 'Tb_K', tb_grid, str(grid_cache_path),
                rheology='maxwell', force_rebuild=args.force
            )
            log.info(f"✓ Maxwell grid cache saved: {grid_cache_path}")
        except Exception as e:
            log.error(f"✗ Failed to build Maxwell grid cache: {e}")
            if args.verbose:
                import traceback
                traceback.print_exc()
            sys.exit(1)
        sys.exit(0)

    # Output file paths (andrade / clathrate variants)
    cache_clath = output_dir / f'{bodyname_lower}_structure_clath.pkl'
    cache_noclath = output_dir / f'{bodyname_lower}_structure_noclath.pkl'

    # Check if caches exist
    if not args.force:
        if cache_clath.exists() and cache_noclath.exists():
            log.info("Both structure caches already exist. Use --force to rebuild.")
            log.info(f"  {cache_clath}")
            log.info(f"  {cache_noclath}")
            sys.exit(0)

    # Generate clathrate variant
    log.info("=" * 70)
    log.info(f"Generating structure WITH clathrates: {test_module}")
    log.info("=" * 70)

    try:
        structure_clath = build_structure_from_pptest(
            test_module,
            rheology=args.rheology,
            force_rebuild=args.force
        )
        save_structure_cache(structure_clath, str(cache_clath))
        log.info(f"✓ Clathrate structure saved: {cache_clath}")
    except Exception as e:
        log.error(f"✗ Failed to build clathrate structure: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

    # Generate no-clathrate variant
    log.info("")
    log.info("=" * 70)
    log.info(f"Generating structure WITHOUT clathrates: {test_module_noclath}")
    log.info("=" * 70)

    try:
        structure_noclath = build_structure_from_pptest(
            test_module_noclath,
            rheology=args.rheology,
            force_rebuild=args.force
        )
        save_structure_cache(structure_noclath, str(cache_noclath))
        log.info(f"✓ No-clathrate structure saved: {cache_noclath}")
    except Exception as e:
        log.error(f"✗ Failed to build no-clathrate structure: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

    # Summary
    log.info("")
    log.info("=" * 70)
    log.info("Structure variant generation complete")
    log.info("=" * 70)
    log.info(f"Body: {bodyname}")
    log.info(f"Rheology: {args.rheology}")
    log.info(f"Output directory: {output_dir.absolute()}")
    log.info("")
    log.info("Generated files:")
    log.info(f"  1. {cache_clath.name} ({cache_clath.stat().st_size / (1024**2):.1f} MB)")
    log.info(f"  2. {cache_noclath.name} ({cache_noclath.stat().st_size / (1024**2):.1f} MB)")
    log.info("")
    log.info("Next steps:")
    log.info("  1. Run validation:")
    log.info(f"     python -m PlanetProfile.Inference.validate_framework \\")
    log.info(f"         --test-module {test_module} \\")
    log.info(f"         --save-cache {cache_clath}")
    log.info("")
    log.info("  2. Create InferenceConfig with one of the cache files:")
    log.info("     config.structure_cache_path = 'path/to/{bodyname_lower}_structure_clath.pkl'")


if __name__ == '__main__':
    main()
