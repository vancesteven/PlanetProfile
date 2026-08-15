#!/usr/bin/env python3
"""
Validation script for inference framework.

Tests all infrastructure components without running full MCMC:
1. Structure caching (load/save)
2. Forward model execution (k2 computation)
3. Log-likelihood evaluation
4. Config serialization (JSON round-trip)
5. Result serialization (pickle round-trip)

Run this before deploying to HPC to ensure framework is working correctly.

Usage:
    python -m PlanetProfile.Inference.validate_framework

    # With specific test module
    python -m PlanetProfile.Inference.validate_framework --test-module PlanetProfile.Test.PPTest41

    # Save structure cache for later use
    python -m PlanetProfile.Inference.validate_framework --save-cache structure.pkl

Author: PlanetProfile Team
Date: 2026-04-29
"""
import argparse
import logging
import tempfile
import numpy as np
from pathlib import Path
import sys

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s'
)
log = logging.getLogger('PlanetProfile.Inference')


def test_structure_cache(test_module: str, save_path: str = None):
    """Test structure caching functionality."""
    log.info("=" * 70)
    log.info("TEST 1: Structure Cache")
    log.info("=" * 70)

    from PlanetProfile.Inference import (
        build_structure_from_pptest,
        save_structure_cache,
        load_structure_cache,
        validate_structure_cache
    )

    # Build structure from test module
    log.info(f"Building structure from {test_module}...")
    try:
        structure_data = build_structure_from_pptest(
            test_module,
            rheology='andrade',
            force_rebuild=True
        )
        log.info("✓ Structure build successful")
    except Exception as e:
        log.error(f"✗ Structure build failed: {e}")
        return False

    # Save to temp file
    with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as tmp:
        tmp_path = tmp.name

    try:
        save_structure_cache(structure_data, tmp_path)
        log.info(f"✓ Structure saved to {tmp_path}")

        # Load back
        loaded_data = load_structure_cache(tmp_path)
        log.info("✓ Structure loaded successfully")

        # Validate
        is_valid, warnings = validate_structure_cache(
            loaded_data,
            expected_bodyname=structure_data['bodyname']
        )
        if is_valid:
            log.info("✓ Structure validation passed")
        else:
            log.warning(f"Structure validation warnings: {warnings}")

        # Save to user-specified path if requested
        if save_path:
            save_structure_cache(structure_data, save_path)
            log.info(f"✓ Structure saved to {save_path} for reuse")

        return True

    except Exception as e:
        log.error(f"✗ Structure cache test failed: {e}")
        return False
    finally:
        # Clean up temp file
        Path(tmp_path).unlink(missing_ok=True)


def test_forward_model(structure_data):
    """Test forward model execution."""
    log.info("")
    log.info("=" * 70)
    log.info("TEST 2: Forward Model")
    log.info("=" * 70)

    from PlanetProfile.Inference import forward_model_k2

    # Test parameters (Andrade, 5D)
    theta = np.array([
        0.3,      # alpha
        0.0,      # log10_zeta
        14.0,     # log10_eta_Ih
        14.0,     # log10_eta_HP
        20.0      # log10_eta_sil
    ])

    log.info("Testing forward model with Andrade rheology...")
    try:
        Re_k2, Im_k2, _ = forward_model_k2(
            theta,
            structure_data,
            rheology='andrade',
            return_heating=False
        )

        if np.isnan(Re_k2):
            log.error("✗ Forward model returned NaN")
            return False

        log.info(f"✓ Forward model executed: Re(k2)={Re_k2:.4f}, Im(k2)={Im_k2:.4f}")

        # Test with heating computation
        log.info("Testing with heating computation...")
        Re_k2, Im_k2, heating = forward_model_k2(
            theta,
            structure_data,
            rheology='andrade',
            return_heating=True
        )

        if heating is not None:
            log.info(f"✓ Heating computed: {list(heating.keys())}")
            for phase, power_W in heating.items():
                log.info(f"    {phase}: {power_W/1e12:.3f} TW")
        else:
            log.warning("  Heating returned None (eccentricity may be zero)")

        return True

    except Exception as e:
        log.error(f"✗ Forward model test failed: {e}", exc_info=True)
        return False


def test_log_likelihood(structure_data):
    """Test log-likelihood function."""
    log.info("")
    log.info("=" * 70)
    log.info("TEST 3: Log-Likelihood")
    log.info("=" * 70)

    from PlanetProfile.Inference import create_log_likelihood

    # Petricca et al. 2025 Titan observables
    observables = {
        'Re_k2': (0.608, 0.048),
        'abs_Im_k2': (0.135, 0.035)
    }

    log.info("Creating log-likelihood function...")
    try:
        log_likelihood = create_log_likelihood(
            observables,
            structure_data,
            rheology='andrade'
        )

        # Test evaluation
        theta = np.array([0.3, 0.0, 14.0, 14.0, 20.0])
        log_like = log_likelihood(theta)

        if np.isfinite(log_like):
            log.info(f"✓ Log-likelihood evaluated: {log_like:.2f}")
        else:
            log.error(f"✗ Log-likelihood returned non-finite value: {log_like}")
            return False

        return True

    except Exception as e:
        log.error(f"✗ Log-likelihood test failed: {e}")
        return False


def test_config_serialization():
    """Test InferenceConfig JSON serialization."""
    log.info("")
    log.info("=" * 70)
    log.info("TEST 4: Config Serialization")
    log.info("=" * 70)

    from PlanetProfile.Inference import InferenceConfig

    # Create test config
    config = InferenceConfig(
        mode='mcmc',
        bodyname='Titan',
        param_space={
            'alpha': {'prior_type': 'uniform', 'bounds': [0.2, 0.4]},
            'log10_eta_Ih': {'prior_type': 'uniform', 'bounds': [12.0, 16.0]},
        },
        observables={
            'Re_k2': (0.608, 0.048),
            'Im_k2': (0.135, 0.035)
        },
        sampler_settings={'n_effective': 500},
        random_state=42
    )

    # Test JSON round-trip
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as tmp:
        tmp_path = tmp.name

    try:
        config.to_json(tmp_path)
        log.info(f"✓ Config saved to JSON: {tmp_path}")

        loaded_config = InferenceConfig.from_json(tmp_path)
        log.info("✓ Config loaded from JSON")

        # Verify fields match
        assert loaded_config.mode == config.mode
        assert loaded_config.bodyname == config.bodyname
        assert loaded_config.param_space == config.param_space
        assert loaded_config.random_state == config.random_state

        log.info("✓ Config round-trip successful")

        # Test hash generation
        hash1 = config.generate_hash()
        hash2 = loaded_config.generate_hash()
        assert hash1 == hash2, "Config hashes don't match"
        log.info(f"✓ Config hash: {hash1}")

        return True

    except Exception as e:
        log.error(f"✗ Config serialization test failed: {e}")
        return False
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def test_result_serialization():
    """Test InferenceResult pickle serialization."""
    log.info("")
    log.info("=" * 70)
    log.info("TEST 5: Result Serialization")
    log.info("=" * 70)

    from PlanetProfile.Inference import InferenceConfig, InferenceResult

    # Create dummy result
    config = InferenceConfig(
        mode='mcmc',
        bodyname='Titan',
        param_space={'alpha': {'prior_type': 'uniform', 'bounds': [0.2, 0.4]}},
        observables={'Re_k2': (0.608, 0.048)},
        sampler_settings={'n_effective': 100}
    )

    samples = np.random.randn(1000, 2)  # 1000 samples, 2 params
    log_likes = -np.random.rand(1000) * 10

    result = InferenceResult(
        config=config,
        samples=samples,
        log_likelihoods=log_likes,
        param_names=['alpha', 'log10_eta_Ih'],
        param_labels=[r'$\alpha$', r'$\log_{10}(\eta_{\rm Ih})$'],
        convergence_metrics={'R_hat': 1.01, 'ESS': 950}
    )

    # Test pickle round-trip
    with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as tmp:
        tmp_path = tmp.name

    try:
        result.save(tmp_path)
        log.info(f"✓ Result saved to pickle: {tmp_path}")

        loaded_result = InferenceResult.load(tmp_path)
        log.info("✓ Result loaded from pickle")

        # Verify data matches. equal_nan=True so a NaN-padded log-likelihood
        # (amortized SBI n_derived subset) round-trips cleanly.
        assert np.allclose(loaded_result.samples, result.samples, equal_nan=True)
        assert np.allclose(loaded_result.log_likelihoods,
                           result.log_likelihoods, equal_nan=True)
        assert loaded_result.param_names == result.param_names

        log.info("✓ Result round-trip successful")

        # Test summary stats
        summary = loaded_result.get_summary_stats()
        log.info(f"✓ Summary stats computed: {list(summary.keys())}")

        # Test best fit
        best_fit = loaded_result.get_best_fit()
        log.info(f"✓ Best fit: {best_fit}")

        return True

    except Exception as e:
        log.error(f"✗ Result serialization test failed: {e}")
        return False
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def main():
    """Run all validation tests."""
    parser = argparse.ArgumentParser(
        description='Validate PlanetProfile inference framework'
    )
    parser.add_argument(
        '--test-module',
        type=str,
        default='PlanetProfile.Test.PPTest41',
        help='Test module to use for structure building (default: PPTest41)'
    )
    parser.add_argument(
        '--save-cache',
        type=str,
        default=None,
        help='Save structure cache to this path for reuse'
    )
    parser.add_argument(
        '--skip-structure',
        action='store_true',
        help='Skip structure building (faster, but skips forward model tests)'
    )

    args = parser.parse_args()

    log.info("PlanetProfile Inference Framework Validation")
    log.info("=" * 70)

    all_passed = True
    structure_data = None

    # Test 1: Structure cache
    if not args.skip_structure:
        if not test_structure_cache(args.test_module, args.save_cache):
            all_passed = False
        else:
            # Load structure for remaining tests
            from PlanetProfile.Inference import build_structure_from_pptest
            structure_data = build_structure_from_pptest(args.test_module, force_rebuild=False)

    # Test 2: Forward model (requires structure)
    if structure_data is not None:
        if not test_forward_model(structure_data):
            all_passed = False
    else:
        log.warning("Skipping forward model test (no structure data)")

    # Test 3: Log-likelihood (requires structure)
    if structure_data is not None:
        if not test_log_likelihood(structure_data):
            all_passed = False
    else:
        log.warning("Skipping log-likelihood test (no structure data)")

    # Test 4: Config serialization
    if not test_config_serialization():
        all_passed = False

    # Test 5: Result serialization
    if not test_result_serialization():
        all_passed = False

    # Summary
    log.info("")
    log.info("=" * 70)
    if all_passed:
        log.info("✓ ALL TESTS PASSED")
        log.info("=" * 70)
        log.info("Framework is ready for HPC deployment")
        return 0
    else:
        log.error("✗ SOME TESTS FAILED")
        log.info("=" * 70)
        log.error("Fix errors before deploying to HPC")
        return 1


if __name__ == '__main__':
    sys.exit(main())
